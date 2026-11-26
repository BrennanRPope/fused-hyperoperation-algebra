"""Figure 6: the other four dials (odds t, similarity ratio r, growth degree beta, holonomy). Run: python3 figures_dials6.py"""
import os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon, Rectangle
from scipy.optimize import brentq

OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "figures")
NAVY, MAG, BLUE, ORANGE, GREY, GREEN = "#0b1f4b", "#e0115f", "#1f6fd1", "#e8892b", "#8a93a8", "#2a9d6f"
plt.rcParams.update({"font.size": 11, "axes.titleweight": "bold", "axes.titlesize": 13,
                     "axes.edgecolor": "#c9cfdd", "axes.labelcolor": NAVY, "text.color": NAVY,
                     "xtick.color": NAVY, "ytick.color": NAVY})
BOX = dict(fc="white", ec="none", pad=0.8, alpha=0.92)
fig, axs = plt.subplots(1, 4, figsize=(19, 5.2), gridspec_kw=dict(wspace=0.32))

# A. odds t = tan(theta): a slope; complement mirrors in the diagonal
ax = axs[0]
ax.plot([0, 1.05], [0, 1.05], color=GREY, ls=":", lw=1.2)
ax.text(0.8, 0.62, "t = 1:\nthe full set", color=GREY, fontsize=9, ha="left", bbox=BOX)
t0 = np.sqrt(2) - 1
ax.plot([0, 1], [0, t0], color=BLUE, lw=2.6)
ax.plot([0, t0], [0, 1], color=MAG, lw=2.2, ls="--")
th = np.linspace(0, np.arctan(t0), 40); ax.plot(0.3*np.cos(th), 0.3*np.sin(th), color=BLUE, lw=1.2)
ax.text(0.33, 0.05, "22.5°", color=BLUE, fontsize=9)
ax.text(1.0, t0 + 0.04, "odds t = tan θ\nt = √2 − 1: the union half", color=BLUE, fontsize=9.5, ha="right", va="bottom", bbox=BOX)
ax.text(t0 + 0.03, 0.97, "complement 1/t:\nthe mirror in the diagonal", color=MAG, fontsize=9.5, va="top", bbox=BOX)
ax.set_xlim(0, 1.08); ax.set_ylim(0, 1.08)
ax.set_title("t · odds (a slope)")
ax.set_xlabel("reading: the chance t/(1 + t)")

# B. similarity ratio r: inputs per piece; dimension D comes out of Moran's equation
ax = axs[1]
s = np.linspace(0, 2.2, 300)
sets = [([0.25, 0.25], BLUE, "half box: r = 1/4, 1/4"), ([1/3, 1/3], GREEN, "Cantor: r = 1/3, 1/3"),
        ([0.5]*3, ORANGE, "gasket: r = 1/2, 1/2, 1/2")]
ax.axhline(1, color=MAG, lw=1.4)
ax.text(2.18, 1.07, "wall: Σ rᵢᴰ = 1", color=MAG, fontsize=9.5, ha="right")
for rs, col, lab in sets:
    ax.plot(s, [sum(r**x for r in rs) for x in s], color=col, lw=2.3, label=lab)
    D = brentq(lambda x: sum(r**x for r in rs) - 1, 0.01, 3)
    ax.plot([D], [1], "o", ms=8, color=col, mec="white", mew=1.8, zorder=4)
    ax.annotate(f"D = {D:.3f}", (D, 1), xytext={BLUE: (0.12, 0.3), GREEN: (0.85, 0.55), ORANGE: (1.62, 1.6)}[col], color=col, fontsize=9,
                arrowprops=dict(arrowstyle="-", color=col, lw=0.8), bbox=BOX)
ax.set_ylim(0, 3.2); ax.set_xlim(0, 2.2)
ax.legend(fontsize=8.5, loc="upper right", frameon=False)
ax.set_title("r · similarity ratio (per piece)")
ax.set_xlabel("s;  the dimension D is a reading")

# C. growth degree beta: a setting of the size measure
ax = axs[2]
X = np.logspace(0, 2, 50)
for bb, col in zip((0, 1, 2, 3), (GREY, BLUE, GREEN, ORANGE)):
    ax.loglog(X, X**bb if bb else np.ones_like(X), color=col, lw=2.3)
    ax.text(105, 100.0**bb if bb else 1.0, f"β = {bb}: doubling weight {['1', '1/2', '1/4', '1/8'][bb]}", color=col, fontsize=9, va="center")
ax.set_xlim(1, 8000); ax.set_ylim(0.5, 3e6)
ax.set_title("β · growth degree (a setting)")
ax.set_xlabel("size X;  count of elements up to X ~ X^β")

# D. holonomy (twist): a reading; the oriented area between two moves
ax = axs[3]
A_, B_ = np.array([1.0, 0.35]), np.array([0.3, 1.0])
ax.add_patch(Polygon([[0, 0], A_, A_ + B_, B_], closed=True, fc=MAG, alpha=0.18, ec=MAG, lw=1.4))
for v, col, lab, off in [(A_, BLUE, "move 1: (a, c)", (0.04, -0.13)), (B_, GREEN, "move 2: (b, e)", (-0.42, 0.06))]:
    ax.annotate("", v, (0, 0), arrowprops=dict(arrowstyle="-|>", color=col, lw=2.2))
    ax.text(v[0] + off[0], v[1] + off[1], lab, color=col, fontsize=9.5)
ax.text(0.65, 0.68, "area a e − b c\n= the order gap", color=MAG, fontsize=10, ha="center", bbox=BOX)
ax.text(0.0, 1.62, "holonomy ω = curvature × area:\nthe turn picked up round a loop", fontsize=9.5, color=NAVY)
ax.set_xlim(-0.15, 1.75); ax.set_ylim(-0.2, 1.9)
ax.set_title("ω · holonomy (twist), a reading")
ax.set_xlabel("read off the curvature dial and two directions")
for a in axs:
    a.set_box_aspect(1)
for a in axs[[0, 3]]:
    a.set_xticks([]); a.set_yticks([])
fig.savefig(os.path.join(OUT, "fig6-other-dials.png"), dpi=150, bbox_inches="tight")
print("wrote fig6-other-dials.png")
