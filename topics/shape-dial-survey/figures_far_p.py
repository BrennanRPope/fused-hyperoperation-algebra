"""Figure 23: the Lee-Yang far zone as formulas in p (run after far_bounds_p.py)."""
import json, math
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.ticker
import far_bounds_p as F
from mpmath import iv
D = json.load(open("far_bounds_p.json")); N = json.load(open("general_p_far.json"))["scan"]
BLUE, ORANGE, GREY, RED, GREEN, PURPLE = "#2b6cb0", "#dd6b20", "#888888", "#c53030", "#2f855a", "#6b46c1"
plt.rcParams.update({"font.size": 11, "axes.spines.top": False, "axes.spines.right": False})
fig = plt.figure(figsize=(15.5, 11.5))
gs = fig.add_gridspec(2, 2, wspace=0.24, hspace=0.36)
# A: certified takeover vs p
ax = fig.add_subplot(gs[0, :])
for k, st in enumerate(D["F2_strips"]):
    ax.plot(st["p"], [st["T"]]*2, color=BLUE, lw=2.4, label="certified takeover T*(p), boxes covering 1.1 ≤ p ≤ 8" if k == 0 else None)
f3 = {float(p): max(Ts) for p, Ts in D["F3"].items() if all(Ts)}
ax.plot(list(f3), list(f3.values()), "o", color=BLUE, mfc="white", ms=8, mew=2, label="certified at sample roundness (step boxes)")
ps = sorted(float(p) for p in N)
ax.plot(ps, [max(N[str(p)]["takeover"], 2.0) for p in ps], "s", color=ORANGE, ms=7, label="where it really takes over (exact pieces, numerical)")
pp = np.linspace(1.0, 17, 400)
ax.plot(pp, 4*pp**2, "--", color=GREY, lw=1.2, label="4 p²")
ax.plot([3], [21], "D", color=PURPLE, ms=8, label="p = 3 certificate's own bound (21)")
ax.set_xscale("log"); ax.set_yscale("log"); ax.set_xlim(0.98, 17.5); ax.set_ylim(1.5, 2500)
ax.set_xticks([1, 1.5, 2, 3, 4, 5, 7, 10, 16]); ax.set_xticklabels(["1", "1.5", "2", "3", "4", "5", "7", "10", "16"])
ax.set_xlabel("roundness p"); ax.set_ylabel("s past which all three far-zone tests are proved")
ax.text(1.02, 140, "near p = 1 the kink margin\nshrinks like 2.2 (p − 1)²,\nso the takeover climbs", fontsize=9.5, color=GREY)
ax.legend(fontsize=9.5, loc="upper left", bbox_to_anchor=(0.2, 1.0))
ax.set_title("Where the far-out argument is proved to take over, at every roundness: about p², with a constant near 3")
# B: margins rising with s at p = 3, top step
ax = fig.add_subplot(gs[1, 0])
ss = np.geomspace(18, 400, 60); P, B = iv.mpf(3), iv.mpf(2)
best = []
for s in ss:
    ms = [F.margins(P, B, iv.mpf(float(s)), k, l) for k, l in F.PARAMS if s >= F.monotone_from(P, l)]
    ms = [m for m in ms if m is not None]
    best.append(max(ms, key=min) if ms else (np.nan, np.nan, np.nan))
best = np.array(best)
for j, (nm, col) in enumerate([("A: |ε| < 1", GREEN), ("B: the phase moves faster than ε", ORANGE), ("C: the sign test (tight one)", BLUE)]):
    ax.plot(ss, best[:, j], color=col, lw=2, label=nm)
km = 1 - abs(math.sin(math.pi*1.5))/3**2
ax.axhline(km, color=BLUE, ls=":", lw=1.2); ax.text(120, km + 0.02, f"kink margin 8/9 ≈ {km:.3f}", color=BLUE, fontsize=9.5)
ax.axhline(0, color="k", lw=1)
T3 = next(st["T"] for st in D["F2_strips"] if st["p"][0] <= 3 <= st["p"][1])
ax.axvline(T3, color=GREY, ls="--"); ax.text(T3*1.03, 0.05, f"T* ≈ {T3:.0f}", color=GREY)
ax.set_xscale("log"); ax.set_ylim(0, 1.1); ax.set_xlim(18, 420)
ax.xaxis.set_minor_formatter(matplotlib.ticker.NullFormatter()); ax.set_xticks([20, 30, 50, 100, 200, 400]); ax.set_xticklabels(["20", "30", "50", "100", "200", "400"])
ax.text(27.5, 0.2, "below about 24 the formulas\ngive no positive bound yet", fontsize=9, color=GREY)
ax.set_xlabel("s"); ax.set_ylabel("proved lower bound on the margin")
ax.legend(fontsize=9, loc="lower right")
ax.set_title("p = 3, top step: the formulas' margins rise to the kink margin")
# C: why p^2
ax = fig.add_subplot(gs[1, 1])
p = 8.0; c = p; y = np.linspace(1e-4, 1.0, 800)
y1 = 1/(p + 2); yz = 2*math.sin(math.pi/p)
for s, col in [(40, RED), (190, GREEN)]:
    w = np.exp(c*np.log(y) - s*y); w /= w.max()
    ax.plot(y, w, color=col, lw=2, label=f"weight yᶜ e^(−sy), s = {s}")
ax.axvspan(0, y1, color=BLUE, alpha=0.12)
ax.text(y1*0.5, 1.04, "window\n1/(p+2)", ha="center", fontsize=9, color=BLUE)
ax.axvline(yz, color=GREY, ls=":"); ax.text(yz + 0.01, 0.75, "first zero of w\nat 2 sin(π/p)", fontsize=9, color=GREY)
ax.set_xlabel("y (the turned path at x = 1)"); ax.set_ylabel("relative weight"); ax.set_ylim(0, 1.18); ax.set_xlim(0, 1.0)
ax.legend(fontsize=9, loc="upper right")
ax.set_title("Why p²: the weight peaks at y ≈ p/s and must sit\ninside a window of width about 1/p (p = 8 shown)")
fig.suptitle("Figure 23. Lee–Yang far zone as formulas in p: proved takeover for every roundness (DRAFT)", fontsize=14, y=0.995)
fig.savefig("figures/fig23-lee-yang-far-p.png", dpi=110, bbox_inches="tight"); print("saved")
