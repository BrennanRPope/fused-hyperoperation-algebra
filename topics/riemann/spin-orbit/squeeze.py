"""Can the zeros-as-spins picture force Lambda <= 0?  Two sides of t = 0.

Forward (capture): an off-line zero rho = beta + i*gamma sits at height y = 2(beta - 1/2) off the line in z units.
  Mean-field capture time t_cap = ln(1 + 2 S y^2)/(4S), S = sum 1/(x - z_k)^2 over the other zeros.
  S grows with height (S ~ n^2 pi^2/3, n = log(z/4pi)/(4pi) zeros per unit z), so t_cap -> 0 as height -> infinity,
  but t_cap > 0 at every finite height for every y > 0. Rigorous counterpart: Ki-Kim-Lee 2009 (for every t > 0 only
  finitely many non-real zeros of H_t; Lambda < 1/2), de Bruijn 1950 (t_cap <= y^2/2).
Backward (break): rewinding from t = 0 to -tau, a pair collides when t_c > -tau. Small GUE gaps (density ~ (pi^2/3)u^2)
  make the broken fraction ~ (pi^2/9)(n sqrt(8 tau))^3, which grows with height for fixed tau: every rewind breaks RH
  somewhere high up. Rigorous counterpart: Rodgers-Tao 2018 (Lambda >= 0), Newman's conjecture.
"""
import json, math
import numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt

K = np.load("../known-tools/rewind_snaps.npz")
g = np.array(json.load(open("../round9/zeros.json"))["gamma"])
n_of = lambda z: np.log(z / (4 * math.pi)) / (4 * math.pi)
out = {}

# ---- forward: measured S between our zeros vs mean field ---------------------------------
mid = K["mid"]; S = K["S"]
Smf = n_of(mid)**2 * math.pi**2 / 3 * 2      # sum over both sides of k^-2 / n^-2, two sides: 2 * zeta(2) n^2
# (S at a midpoint of a gap: nearest zeros at +-1/(2n): sum 2 * sum_k 1/((k-1/2)/n)^2 = 2 n^2 * (pi^2/2) = pi^2 n^2)
Smid = math.pi**2 * n_of(mid)**2
ratio = np.median(S / Smid)
out["S_over_meanfield_at_gap_midpoints_median"] = float(ratio)
H = np.logspace(3, 22, 200)                 # height gamma
zH = 2 * H
Sbg = math.pi**2 * n_of(zH)**2 * 0.5         # an off-line zero sits at a typical spot: use pi^2 n^2/2 (between
# the midpoint value pi^2 n^2 and the on-zero value pi^2 n^2/3); the factor-of-2 choice is stated in the README
def tcap(y, S): return np.log1p(2 * S * y**2) / (4 * S)
caps = {}
for beta in (1.0, 0.75, 0.6, 0.51):
    y = 2 * (beta - 0.5)
    caps[str(beta)] = tcap(y, Sbg)
out["tcap_at_3e12_beta1"] = float(np.interp(3e12, H, caps["1.0"]))
out["tcap_at_3e12_beta1_S_range"] = [float(tcap(1.0, math.pi**2 * n_of(6e12)**2 / 3)),
                                     float(tcap(1.0, math.pi**2 * n_of(6e12)**2))]
out["debruijn_beta1"] = 0.5

# ---- backward: broken fraction vs rewind depth --------------------------------------------
tc = K["tc"]
taus = np.logspace(-4, 0, 40)
lo, hi = 1000, 4500               # middle zeros, to avoid edge effects
n_mid = float(np.median(n_of(mid[lo:hi])))
frac = [float(np.mean(np.isfinite(tc[lo:hi]) & (tc[lo:hi] > -t))) for t in taus]
pred = [(math.pi**2 / 9) * (n_mid * math.sqrt(8 * t))**3 for t in taus]
out["backward"] = {"taus": taus.tolist(), "broken_fraction": frac, "gue_small_gap_prediction": pred, "n_mid": n_mid}
# predicted broken fraction at fixed tau = 1e-3 vs height
Hb = np.logspace(3, 22, 200)
fb = (math.pi**2 / 9) * (n_of(2 * Hb) * math.sqrt(8e-3))**3
out["broken_fraction_tau1e-3_at_1e3_1e22"] = [float(fb[0]), float(fb[-1])]
json.dump(out, open("squeeze.json", "w"), indent=1)

fig, ax = plt.subplots(1, 3, figsize=(18, 5.6))
a = ax[0]
for b, c in zip(caps, ("C3", "C1", "C2", "C0")):
    a.loglog(H, caps[b], color=c, label=f"off-line zero at β = {b}")
a.axhline(0.5, color="k", ls="--", lw=0.8); a.text(2e3, 0.55, "de Bruijn 1950: Λ ≤ ½ (β < 1)", fontsize=8)
a.axhline(0.2, color="k", ls=":", lw=0.8); a.text(2e3, 0.13, "Platt–Trudgian 2021: Λ ≤ 0.2", fontsize=8)
a.axvline(3e12, color="0.5", lw=0.8); a.text(4e12, 2e-3, "RH checked\nbelow 3·10¹²", fontsize=8, color="0.4")
a.set_xlabel("height γ"); a.set_ylabel("heat time to pull it onto the line")
a.set_ylim(1e-4, 1)
a.set_title("F. Forward: an off-line zero is pulled in faster the higher it sits\n"
            "(the crowd S grows), but at every finite height it takes t > 0", fontsize=10)
a.legend(fontsize=8, loc="lower left")
a = ax[1]
mk = np.asarray(frac) > 0
a.loglog(np.asarray(taus)[mk], np.asarray(frac)[mk], "o", color="C3", ms=4, label="our zeros 1000–4500: pairs broken by rewind depth τ (none below τ ≈ 0.004)")
a.loglog(taus, pred, "k-", lw=1, label="GUE small gaps: (π²/9)(n√(8τ))³")
a.set_xlabel("rewind depth τ (run heat back to t = −τ)"); a.set_ylabel("fraction of pairs that left the line")
a.set_title("B. Backward: rewinding by any τ breaks some pairs.\n"
            "The broken share grows like the GUE small-gap law (τ^1.5), at about half its size", fontsize=10)
a.legend(fontsize=8, loc="lower right"); a.set_ylim(1e-5, 1.5)
a = ax[2]
a.loglog(Hb, fb, color="C3", label="share of pairs broken by τ = 0.001")
a.loglog(H, caps["1.0"], color="C0", label="time to capture a β = 1 zero")

a.set_xlabel("height γ")
a.set_title("The squeeze: going higher, rewinding breaks more pairs (Λ ≥ 0, Rodgers–Tao)\n"
            "and capture gets faster (Ki–Kim–Lee). They meet at t = 0 only at infinite height", fontsize=10)
a.legend(fontsize=8)
plt.tight_layout(); plt.savefig("figures/fig_squeeze.png", dpi=105)
print(json.dumps({k: v for k, v in out.items() if k != "backward"}, indent=1))
print("backward sample", list(zip(taus[::6].round(5), np.round(frac[::6], 4), np.round(pred[::6], 4))))
