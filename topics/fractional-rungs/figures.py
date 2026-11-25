import sys, json, math
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
sys.path.insert(0, "/mnt/project-files/math-theory/fractional-rungs")
from regular_iteration import Es, B

OUT = "/mnt/project-files/math-theory/fractional-rungs/figures/"
plt.rcParams.update({"font.family": "serif", "font.size": 11, "axes.spines.top": False, "axes.spines.right": False})
RED, BLUE, INK, SOFT = "#c1121f", "#0b4f8a", "#10131a", "#5a6070"

def col(s):
    if s < 0: return RED
    if s == 0: return INK
    return BLUE

# ---- F1: the bend is the rung (level lines near and far) ----
S = [-1, -0.5, 0, 0.5, 1, 2]
names = {-1: "s = −1\nlog(eᵃ+eᵇ−1)", -0.5: "s = −½", 0: "s = 0\naddition", 0.5: "s = ½\nhalf rung", 1: "s = 1\ncircle operation", 2: "s = 2"}
fig, axes = plt.subplots(2, 6, figsize=(16, 6.2))
for row, (lim, levels) in enumerate([(3, np.arange(0.5, 6.01, 0.5)), (40, np.arange(5, 60.1, 5))]):
    g = np.linspace(lim/400, lim, 170)
    A, Bg = np.meshgrid(g, g)
    for j, s in enumerate(S):
        ax = axes[row, j]
        Z = B(s, A.ravel(), Bg.ravel()).reshape(A.shape)
        ax.contour(A, Bg, Z, levels=levels, colors=col(s), linewidths=1.1)
        ax.set_xlim(0, lim); ax.set_ylim(0, lim); ax.set_aspect("equal")
        ax.set_xticks([0, lim]); ax.set_yticks([0, lim])
        if row == 0: ax.set_title(names[s], fontsize=11, color=col(s))
        if j == 0: ax.set_ylabel("near 0 (zoom in)" if row == 0 else "far out (zoom out)")
fig.suptitle("Level lines a ∘ b = const at real rung heights s. Flat at addition; bowing out below it, in above it. "
             "Zoomed out, the lower rungs become max (corners), the upper ones multiplication (hyperbolas).", fontsize=11.5)
fig.tight_layout(rect=[0, 0, 1, 0.94])
fig.savefig(OUT + "f1_level_lines.png", dpi=130); plt.close(fig)

# ---- F2: one smooth curve through every rung ----
ss = np.linspace(-2.5, 2.4, 99)
v11 = np.array([float(B(s, 1.0, 1.0)) for s in ss])
e1 = np.array([float(Es(s, np.array([1.0]))[0]) for s in ss])
fig, axes = plt.subplots(1, 3, figsize=(15, 4.6))
ax = axes[0]
ax.plot(ss, v11, color=BLUE, lw=2, label="1 ∘ 1 at height s")
ax.plot(ss, e1, color=SOFT, lw=1.6, ls="--", label="the rung's own e: Eˢ(1)")
for s in [-2, -1, 0, 1, 2]:
    y = float(B(s, 1.0, 1.0)); ax.plot([s], [y], "o", color=INK, ms=5); ax.annotate(f"{y:.3f}", (s, y), xytext=(4, -12), textcoords="offset points", fontsize=9)
ax.set_xlabel("rung height s"); ax.set_yscale("log"); ax.legend(fontsize=9, loc="upper left")
ax.set_title("A smooth curve through every rung")
ax = axes[1]
ax.plot(ss, ss, color=BLUE, lw=2, label="corner weight = s (a straight line)")
ax.plot(ss, ss*(ss - 1)/4, color=RED, lw=2, label="bend = s(s−1)/4 (a parabola)")
ax.axhline(0, color=SOFT, lw=.8)
for s in [0, 1]: ax.plot([s], [0], "o", color=INK, ms=5)
ax.annotate("addition and the circle operation:\nthe only flat ones", (0.5, 0), xytext=(-0.6, 1.4), fontsize=9, arrowprops=dict(arrowstyle="->", color=SOFT))
ax.set_xlabel("rung height s"); ax.legend(fontsize=9, loc="upper left"); ax.set_ylim(-1.5, 3.0)
ax.set_title("Its shape near 0: a ∘ b = a + b + s·ab + bend·(a²b + ab²) + …", fontsize=10.5)
ax = axes[2]
consts = json.load(open("/mnt/project-files/math-theory/fused-constants/constants.json"))
ks = [1, 2, 3] + list(range(4, 13))
vals = [-1.0, 0.0, float(consts["rung3"]["c"]) - 1] + [float(consts["linear"][str(k)]["c"]) - 1 for k in range(4, 13)]
ax.plot(ks, vals, "o", color=INK, ms=5, label="tower constants (fused), from the constants work")
from scipy.interpolate import PchipInterpolator
kk = np.linspace(1, 12, 300)
ax.plot(kk, PchipInterpolator(ks, vals)(kk), color=SOFT, lw=1.2, ls="--", label="one smooth curve through them")
cs = np.polyfit(ks[2:], np.log(1 - np.array(vals[2:])), 1)
ax.plot(kk, 1 - np.exp(np.polyval(cs, kk)), color=RED, lw=1, ls=":", label=f"another: 1 − C·{math.exp(cs[0]):.2f}ᵏ")
ax.axhline(1, color=SOFT, lw=.8)
ax.set_xlabel("rung"); ax.legend(fontsize=8.5, loc="lower right"); ax.set_ylim(-1.1, 1.1)
ax.set_title("Tower constants: the dots are fixed, the in-between is a choice", fontsize=10.5)
fig.tight_layout(); fig.savefig(OUT + "f2_smooth_curve.png", dpi=130); plt.close(fig)

# ---- F3: the corner dial and KdV ----
fig, axes = plt.subplots(1, 2, figsize=(12, 4.2))
ax = axes[0]
r = np.logspace(-2, 2, 400); Akdv = ((r - 1)/(r + 1))**2
ax.semilogx(r, Akdv, color=BLUE, lw=2)
ax.axhline(1, color=SOFT, lw=.8, ls="--"); ax.axhline(0, color=SOFT, lw=.8, ls="--")
ax.text(0.012, 0.92, "corner full: the circle operation", fontsize=9, color=SOFT); ax.text(0.6, 0.04, "corner empty: addition", fontsize=9, color=SOFT)
ax.set_xlabel("speed ratio of the two solitons, k₂/k₁"); ax.set_ylabel("corner weight A")
ax.set_title("KdV collisions live between addition and the circle operation")
ax = axes[1]
g = np.linspace(0, 1, 200); A_, B_ = np.meshgrid(g, g)
for kd, c_, ls in [(-1, RED, "-"), (0, INK, "-"), (0.5, BLUE, "--"), (1, BLUE, "-")]:
    ax.contour(A_, B_, A_ + B_ + kd*A_*B_, levels=[0.8], colors=c_, linestyles=ls, linewidths=1.6)
ax.set_aspect("equal"); ax.set_xticks([0, 1]); ax.set_yticks([0, 1])
ax.text(0.62, 0.62, "k = −1 (OR)", color=RED, fontsize=9); ax.text(0.42, 0.44, "k = 0", color=INK, fontsize=9); ax.text(0.47, 0.17, "k = ½", color=BLUE, fontsize=9)
ax.text(0.03, 0.42, "k = 1", color=BLUE, fontsize=9)
ax.set_title("The dial a + b + k·ab = 0.8: the corner sets the bend", fontsize=10.5)
fig.tight_layout(); fig.savefig(OUT + "f3_dial_kdv.png", dpi=130); plt.close(fig)
print("figures written")
