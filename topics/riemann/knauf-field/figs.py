import numpy as np, json, mpmath as mp, matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
plt.rcParams.update({"font.size": 10})
gam = [float(mp.zetazero(k).imag) for k in range(1, 6)]
# ---------- Fig A: temperature is roundness p; turning p imaginary makes it a field
fig, ax = plt.subplots(1, 3, figsize=(15, 4.6))
n = np.arange(1, 41); E = np.log(n)
for p, col in ((0.5, "tab:red"), (1.0, "tab:orange"), (2.0, "tab:blue")):
    ax[0].vlines(E + {0.5: -0.02, 1.0: 0, 2.0: 0.02}[p], 0, n**(-p), color=col, lw=1.6, label=f"p = {p}: weight n^(-p)")
ax[0].set_xlabel("height s = ln n  (energy of level n)"); ax[0].set_ylabel("weight of level")
ax[0].set_title("A1. Real p is temperature: it sets each level's weight")
ax[0].legend(frameon=False)
ax[0].text(1.6, 0.62, "ζ(p) = Σ n^(−p) = (p-norm of 1/n)^p\nzeta is the p-th power of a p-norm", fontsize=9)
for a, (edge1, edge2, label, zero_re) in zip(ax[1:], ((1.0, 0.5, "primon gas: ζ(p)", 0.5), (2.0, 1.5, "Knauf chain: ζ(p−1)/ζ(p)", 1.5))):
    a.axvspan(edge1, edge1 + 1.6, color="tab:green", alpha=0.10)
    a.axvline(edge1, color="tab:green", lw=2); a.axvline(edge2, color="tab:blue", lw=2, ls="--")
    a.plot([zero_re]*5, gam, "o", color="tab:red", ms=6, label="zeros (RH: on the dashed line)")
    a.annotate("", xy=(edge1 + 1.4, 9), xytext=(edge1 + 0.3, 9), arrowprops=dict(arrowstyle="->", lw=1.5))
    a.text(edge1 + 0.32, 10.2, "real p = temperature", fontsize=9)
    a.annotate("", xy=(edge1 + 1.05, 33), xytext=(edge1 + 1.05, 21), arrowprops=dict(arrowstyle="->", lw=1.5))
    a.text(edge1 + 1.12, 25, "imaginary p\n= field on height", fontsize=9)
    a.text(edge1 + 0.04, 37.3, "ℓ¹ edge = phase transition", color="tab:green", fontsize=9)
    a.text(edge2 - 0.62, 39.4, "ℓ² edge = RH line", color="tab:blue", fontsize=9)
    a.set_xlim(edge2 - 0.7, edge1 + 1.6); a.set_ylim(0, 41); a.set_xlabel("Re p (temperature)"); a.set_ylabel("Im p (field)")
    a.legend(loc="lower right", frameon=False, fontsize=8)
ax[1].set_title("A2. Wick-turning p: the primon gas"); ax[2].set_title("A3. Same turn on Knauf's chain")
fig.tight_layout(); fig.savefig("figures/figA_temperature_is_p.png", dpi=150); plt.close(fig)
# ---------- Fig B: l1 -> l2 edge, three models
d = json.load(open("l2edge.json"))
fig, ax = plt.subplots(1, 2, figsize=(12, 4.4))
for key, lab, col in (("zeta", "primon gas: Σ 1/n at line ½", "tab:red"), ("knauf", "Knauf chain: Σ φ(q)²/q³ at line 3/2", "tab:blue"), ("tau", "Ramanujan τ: Σ τ(n)²/n¹² at line 6", "tab:purple")):
    X = np.array(d[key]["X"]); S = np.array(d[key]["S"])
    ax[0].plot(np.log(X), S, "o-", ms=3, color=col, label=f"{lab}  (slope {d[key]['log_slope']:.3f})")
ax[0].set_xlabel("ln X"); ax[0].set_ylabel("sum of squared level weights up to X")
ax[0].set_title("B1. At the RH line the squared weights grow like ln X:\nexactly on the ℓ² edge (marginal), in all three models"); ax[0].legend(frameon=False, fontsize=8)
p = np.linspace(0.6, 3.2, 300)
ax[1].plot(p, 1/p, color="k"); ax[1].axhline(1, color="tab:green", lw=1); ax[1].axhline(0.5, color="tab:blue", lw=1, ls="--")
ax[1].plot([1, 2], [1, 0.5], "o", color="tab:red", ms=8)
ax[1].text(1.05, 1.06, "p = 1: absolute convergence = phase transition", fontsize=9, color="tab:green")
ax[1].text(2.05, 0.55, "p = 2 (round): mean square = RH line", fontsize=9, color="tab:blue")
ax[1].set_xlabel("roundness p of the norm"); ax[1].set_ylabel("edge temperature σ = 1/p\n(weights n^(−σ) are in ℓ^p above it)")
ax[1].set_title("B2. The RH line sits ½ below the transition because\nroundness goes from 1 to 2 (shown for ζ)")
fig.tight_layout(); fig.savefig("figures/figB_l1_to_l2.png", dpi=150); plt.close(fig)
# ---------- Fig C: finite chain zeros
f = json.load(open("finite_zeros.json")); r = json.load(open("refined.json")); e = json.load(open("errscan.json"))
fig, ax = plt.subplots(1, 3, figsize=(16, 4.8))
cols = {"12": "tab:orange", "14": "tab:green", "16": "tab:purple", "18": "tab:blue"}
for nn, v in f.items():
    z = np.array([[a, b] for a, b, w in v["zeros"]])
    ax[0].plot(z[:, 0], z[:, 1], "o", ms=5, color=cols[nn], label=f"chain length n = {nn} (2^{nn} states)")
ax[0].axvline(1.5, color="tab:blue", ls="--"); ax[0].axvline(2.0, color="tab:green")
ax[0].plot([1.5]*5, gam, "x", color="k", ms=9, mew=2, label="ζ(β−1) zeros 3/2 + iγ")
ax[0].text(1.53, 8.5, "zeros sit near 3/2 + iγ\nalready at n = 12", fontsize=9)
ax[0].text(2.05, 2.0, "finite-size zeros creep\ntoward the transition at 2", fontsize=8)
ax[0].set_xlabel("Re β (temperature)"); ax[0].set_ylabel("Im β (field)"); ax[0].set_xlim(1.0, 3.0); ax[0].set_ylim(0, 36)
ax[0].set_title("C1. Zeros of finite Knauf chains"); ax[0].legend(frameon=False, fontsize=7, loc="upper right")
ns = sorted(int(k) for k in r["rows"])
for j in range(5):
    ax[1].plot(ns, [r["rows"][str(k)][j][0] - 1.5 for k in ns], "o-", ms=3, label=f"zero near γ = {gam[j]:.2f}")
ax[1].axhline(0, color="tab:blue", ls="--"); ax[1].set_xlabel("chain length n"); ax[1].set_ylabel("Re(zero) − 3/2")
ax[1].set_title("C2. They wobble ±0.03 around the line\nand do not settle as n grows (n = 10…24)"); ax[1].legend(frameon=False, fontsize=7)
for sg, col in (("1.2", "tab:red"), ("1.5", "tab:blue"), ("1.8", "tab:green"), ("2.2", "tab:gray")):
    ax[2].semilogy(range(10, 25), e[sg], "o-", ms=3, color=col, label=f"Re β = {sg}")
ax[2].set_xlabel("chain length n"); ax[2].set_ylabel("RMS |Z_n − ζ(β−1)/ζ(β)| over Im β ∈ [5, 35]")
ax[2].set_title("C3. The finite chain approaches ζ(β−1)/ζ(β) slowly\non the RH line, and not visibly at Re β = 1.2"); ax[2].legend(frameon=False, fontsize=8)
fig.tight_layout(); fig.savefig("figures/figC_finite_chain_zeros.png", dpi=150); plt.close(fig)
# ---------- Fig D: magnetic Lee-Yang circle test
c = json.load(open("circle.json"))
fig, ax = plt.subplots(1, 2, figsize=(10, 4.6))
th = np.linspace(0, 2*np.pi, 300)
for s, col in (("1.5", "tab:blue"), ("3.0", "tab:red")):
    z = np.array(c["12"][s]["roots"]); ax[0].plot(z[:, 0], z[:, 1], "o", color=col, label=f"temperature β = {s}")
ax[0].plot(np.cos(th), np.sin(th), "k--", lw=1, label="unit circle (Lee–Yang)")
ax[0].set_aspect("equal"); ax[0].set_xlabel("Re w"); ax[0].set_ylabel("Im w")
ax[0].set_title("D1. Magnetic field zeros, n = 12:\noff the unit circle, in a ring"); ax[0].legend(frameon=False, fontsize=7, loc="center")
for nn, col in (("8", "tab:orange"), ("10", "tab:green"), ("12", "tab:blue")):
    S = sorted(float(s) for s in c[nn]); lo = [c[nn][str(s) if str(s) in c[nn] else s]["absmin"] for s in c[nn]]
    keys = list(c[nn].keys()); sv = [float(k) for k in keys]
    ax[1].fill_between(sv, [c[nn][k]["absmin"] for k in keys], [c[nn][k]["absmax"] for k in keys], color=col, alpha=0.35, label=f"n = {nn}")
ax[1].axhline(1, color="k", ls="--", lw=1); ax[1].set_xlabel("temperature β"); ax[1].set_ylabel("|w| of zeros (band)")
ax[1].set_title("D2. Zeros spread over a band, not a circle,\nand the band drifts off |w| = 1 as β grows"); ax[1].legend(frameon=False, fontsize=8)
fig.tight_layout(); fig.savefig("figures/figD_magnetic_circle.png", dpi=150); plt.close(fig)
print("ok")
