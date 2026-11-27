import json, math
import numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
R = json.load(open("spinorbit.json")); TR = np.load("spinorbit_traj.npz")
g = np.array(json.load(open("../round9/zeros.json"))["gamma"]); T = g[-1]
K = np.load("../known-tools/rewind_snaps.npz")

# ---- figure 1: primes as spins ----
fig, ax = plt.subplots(1, 3, figsize=(18, 5.4))
a = ax[0]
ps = sorted(int(p) for p in R["lean"])
pp = np.linspace(2, 32, 300)
Lbar = np.mean(1 / np.log(g / (2 * math.pi)))
a.plot(pp, -np.log(pp) / np.sqrt(pp) * Lbar, "k-", lw=1, label="Landau's formula: −(log p/√p)·⟨1/log(γ/2π)⟩")
a.errorbar(ps, [R["lean"][str(p)]["mean_cos"] for p in ps], yerr=[R["lean"][str(p)]["se"] for p in ps],
           fmt="o", color="C3", capsize=2, label="each prime p: mean cos(γ log p) over the 5500 zeros")
rat = list(R["pairs"])
xr = [eval(r) * 8 + 4 for r in rat]
a.errorbar([2 + 4 * i for i in range(len(rat))], [R["pairs"][r]["mean_cos"] for r in rat],
           yerr=[R["pairs"][r]["se"] for r in rat], fmt="s", color="C0", capsize=2,
           label="two primes against each other: mean cos(γ log p/q)")
for i, r in enumerate(rat):
    a.text(2 + 4 * i, 0.02, r, ha="center", fontsize=8, color="C0")
a.axhline(0, color="k", lw=0.5)
a.set_xlabel("p (red) / position of the pair p/q (blue)"); a.set_ylabel("mean cosine of the phase at the zeros")
a.set_title("I-a. Primes as spins: at the zeros every prime leans toward 180°,\n"
            "matched to within 0.001; two primes never lock to each other (0 to within 0.0005)", fontsize=10)
a.legend(fontsize=8, loc="lower right"); a.set_ylim(-0.17, 0.05)
a = ax[1]
w = R["lean2_windows"]
a.errorbar([v["mid"] for v in w], [v["mean"] for v in w], yerr=[v["se"] for v in w], fmt="o", color="C3",
           capsize=2, label="p = 2, zeros in each height window")
hh = np.linspace(200, T, 300)
a.plot(hh, -(math.log(2) / math.sqrt(2)) / np.log(hh / (2 * math.pi)), "k-", lw=1, label="−(log 2/√2)/log(γ/2π)")
a.axhline(0, color="k", lw=0.5)
a.set_xlabel("height γ"); a.set_ylabel("mean cos(γ log 2)")
a.set_title("I-b. The lean fades with height like 1/log γ:\n"
            "the opposite of a resonance lock, which tightens as the orbit shrinks", fontsize=10)
a.legend(fontsize=8)
a = ax[2]
th2 = np.mod(g * math.log(2), 2 * math.pi); th32 = np.mod(g * math.log(1.5), 2 * math.pi)
bins = np.linspace(0, 2 * math.pi, 25)
a.hist(th2, bins=bins, density=True, histtype="step", lw=1.5, color="C3", label="phase of 2^{iγ}")
a.hist(th32, bins=bins, density=True, histtype="step", lw=1.5, color="C0", label="phase of (3/2)^{iγ}: 2 against 3")
c2 = R["lean"]["2"]["mean_cos"]
tt = np.linspace(0, 2 * math.pi, 200)
Lb = np.mean(1 / np.log(g / (2 * math.pi)))
dd = 1 + 2 * sum(-(math.log(2) / 2**(k / 2)) * Lb * np.cos(k * tt) for k in range(1, 40))
a.plot(tt, dd / (2 * math.pi), "k--", lw=1, label="Landau over all powers of 2: 1 + 2Σ_k ⟨cos kθ⟩ cos kθ")
a.axhline(1 / (2 * math.pi), color="k", lw=0.5)
a.set_xticks([0, math.pi, 2 * math.pi]); a.set_xticklabels(["0°", "180°", "360°"])
a.set_xlabel("phase at the zeros")
a.set_title("I-c. Where the phases sit: prime 2 avoids 0° (its powers 4, 8, ... add up there);\n"
            "the angle between 2 and 3 is flat (no resonance between primes)", fontsize=10)
a.legend(fontsize=8, loc="lower center")
plt.tight_layout(); plt.savefig("figures/fig_primes_as_spins.png", dpi=105)

# ---- figure 2: zeros as spins ----
fig, ax = plt.subplots(1, 2, figsize=(14, 5.4))
a = ax[0]
for i, cp in enumerate(R["capture"]):
    ys = TR[f"{cp['j']}_{cp['y0']}"]
    col = f"C{i}"
    a.plot(ys[:, 0], ys[:, 1], color=col, lw=1.5, label=f"inserted pair at γ ≈ {cp['x0']/2:.0f}, y₀ = {cp['y0']}: "
           f"captured at t = {cp['t_capture_nbody']:.3f} (formula {cp['t_formula']:.3f})")
    S = cp["S"]; y0 = cp["y0"]; t = np.linspace(0, cp["t_formula"], 200)
    a.plot(t, np.sqrt(np.maximum((y0**2 + 1 / (2 * S)) * np.exp(-4 * S * t) - 1 / (2 * S), 0)), ":", color=col, lw=1)
    tb = np.linspace(0, y0**2 / 2, 200)
    a.plot(tb, np.sqrt(np.maximum(y0**2 - 2 * tb, 0)), "--", color="0.6", lw=0.8)
a.plot([], [], "--", color="0.6", label="de Bruijn's bound √(y₀² − 2t) (1950)")
a.set_xlabel("forward heat time t"); a.set_ylabel("height of the pair off the line (z units)")
a.set_title("II-a. Zeros as spins: an off-line pair (a toy, inserted among the real zeros)\n"
            "is pulled onto the line by forward heat flow, faster than de Bruijn's bound", fontsize=10)
a.legend(fontsize=7.5)
a = ax[1]
x = K["S"] * K["d0"]**2 / 2
a.hist(np.log10(x), bins=80, color="C7")
a.axvline(0, color="C3", lw=1.5)
f = float((x < 1).mean())
a.text(-2.75, a.get_ylim()[1] * 0.55, f"{100*f:.0f}% of pairs: S d₀² < 2\nrewind makes them a\nconjugate pair (one up,\none down): the 180° family",
       fontsize=9, color="C3")
a.text(0.22, a.get_ylim()[1] * 0.55, f"{100*(1-f):.0f}%: S d₀² ≥ 2\nthey stay on the line\nuntil a neighbour joins:\nthe 0° family", fontsize=9)
a.set_xlabel("log₁₀(S d₀²/2) for each neighbouring pair of the 5500 zeros")
a.set_title("II-b. Two families, decided by the starting state (as in Schnittman 2004,\n"
            "where the spin nearer the orbit's axis decides): here the gap against the crowd", fontsize=10)
plt.tight_layout(); plt.savefig("figures/fig_zeros_as_spins.png", dpi=105)
print("frac x<1", f)
