import json, sys, math
import numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
sys.path.insert(0, "../unified-hypothesis")
import kern
R = json.load(open("rewind.json")); D = np.load("rewind_snaps.npz")
g = np.array(json.load(open("../round9/zeros.json"))["gamma"]); z = 2 * g
fig, ax = plt.subplots(1, 3, figsize=(18, 5.6))
# 1. low-height check
K = kern.riemann(500)
a = ax[0]
for t in np.arange(0, -4.6, -0.05):
    rz = kern.real_zeros(K, 74.0, dx=0.01, t=float(t))
    a.plot(rz, np.full(len(rz), t), "k.", ms=1.5)
tf = R["low_formula"]["pairs"][3]["t_c"]
a.axhline(tf, color="C3", ls="--", lw=1, label=f"pair formula: collision at t = {tf:.2f}")
a.axhline(R["low_direct"]["t_first_nonreal"], color="C0", ls=":", lw=1.2,
          label=f"direct rewind of H_t: first non-real zero by t = {R['low_direct']['t_first_nonreal']:.2f}")
a.set_xlabel("z = 2γ"); a.set_ylabel("heat time t (rewind is down)")
a.set_title("A1. Rewinding H_t: real zeros (dots) drift together; the pair at γ = 30.42, 32.94\n"
            "meets and leaves the line where the two-zero formula says", fontsize=10)
a.legend(fontsize=8, loc="lower left"); a.set_xlim(20, 74)
# 2. collision time per pair
a = ax[1]
tc = D["tc"]; mid = D["mid"] / 2
ok = np.isfinite(tc)
a.semilogy(mid[ok], -tc[ok], ".", ms=1.5, color="C7", label="each close pair: −t_c")
run = np.minimum.accumulate(np.where(ok, -tc, np.inf))
a.semilogy(mid, run, color="C3", lw=1.2, label="closest to 0 so far")
b = R["pairs"]["top10"][0]
a.annotate(f"γ ≈ {b['gamma'][0]:.2f}: t_c = {b['t_c']:.2e}\n(N-body rewind: {R['nbody_closest']['t_collision']:.2e})",
           (b["gamma"][0], -b["t_c"]), (1800, 1.5e-2), fontsize=8, arrowprops=dict(arrowstyle="->"))
a.set_xlabel("height γ"); a.set_ylabel("−t_c (how far back the pair survives)")
a.set_title("A2. Each close pair of zeros gives a lower bound Λ ≥ t_c. The bound creeps toward 0\n"
            "with height (Lehmer pairs; Csordas–Smith–Varga 1994); Rodgers–Tao 2018: Λ ≥ 0", fontsize=10)
a.legend(fontsize=8, loc="lower left")
# 3. forward flow
a = ax[2]
dens = lambda w: np.log(w / (4 * math.pi)) / (4 * math.pi)
u = np.linspace(0, 3, 300)
for k, c in (("0.0", "C0"), ("0.05", "C1"), ("0.2", "C2")):
    w = D[k][1000:4500]; m = (w[:-1] + w[1:]) / 2
    ng = np.diff(w) * dens(m)
    a.hist(ng, bins=60, range=(0, 3), density=True, histtype="step", color=c, lw=1.4,
           label=f"t = {k}: var {R['forward'][k]['var']:.3f}, smallest gap {R['forward'][k]['min']:.2f}")
a.plot(u, 32 / math.pi**2 * u**2 * np.exp(-4 * u**2 / math.pi), "k--", lw=0.8, label="GUE (Wigner surmise)")
a.set_xlabel("gap / mean gap (zeros 1000 to 4500)")
a.set_title("A3. Heat forward (t > 0) evens the gaps out, toward a lattice\n"
            "(Rodgers–Tao). RH puts the zeros exactly at the border, t = 0", fontsize=10)
a.legend(fontsize=8)
plt.tight_layout(); plt.savefig("figures/fig_rewind.png", dpi=105)
