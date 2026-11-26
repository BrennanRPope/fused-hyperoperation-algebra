"""DRAFT figure O9: gravity's 6 through the double copy (Kerr-Schild form). Run after checks_dcopy.py."""
import os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
HERE = os.path.dirname(os.path.abspath(__file__)); OUT = os.path.join(HERE, "..", "figures")
NEWT, CROSS, SPACE, SR, EN = "#cfcfcf", "#2b6cb0", "#888888", "#d9822b", "#7aa6d6"
plt.rcParams.update({"font.size": 11, "axes.titlesize": 12})
fig, axs = plt.subplots(1, 3, figsize=(18, 5.6), gridspec_kw={"width_ratios": [1.15, 1, 1.1]})

# A: the binomial square
ax = axs[0]; ax.set_aspect("equal"); ax.axis("off")
e, p = 2.0, 1.0; x0, y0 = 0.3, 0.2
ax.add_patch(Rectangle((x0, y0 + p), e, e, color=NEWT)); ax.text(x0 + e / 2, y0 + p + e / 2, "E²\nNewton's pull\nadds 0\n(orbit closes,\nat any strength)", ha="center", va="center", fontsize=9.5)
ax.add_patch(Rectangle((x0 + e, y0 + p), p, e, color=CROSS)); ax.add_patch(Rectangle((x0, y0), e, p, color=CROSS))
ax.text(x0 + e + p / 2, y0 + p + e / 2, "E·p", ha="center", va="center", color="white", fontsize=10, fontweight="bold")
ax.text(x0 + e / 2, y0 + p / 2, "E·p", ha="center", va="center", color="white", fontsize=10, fontweight="bold")
ax.add_patch(Rectangle((x0 + e, y0), p, p, color=SPACE)); ax.text(x0 + e + p / 2, y0 + p / 2, "p²\nspace", ha="center", va="center", color="white", fontsize=9.5)
ax.text(x0 + e + p + 0.15, y0 + p + e / 2, "cross term 2E·p:\nstrength 2M, so\nadds (2)² × 1 = 4", color=CROSS, fontsize=10, va="center")
ax.text(x0 + e + p + 0.15, y0 + p / 2, "with Newton:\nadds 2", color=SPACE, fontsize=10, va="center")
ax.text(x0, y0 + e + p + 0.25, "gravity:  (2M/r)(k·p)²,   k·p = E + p", fontsize=11, fontweight="bold")
ax.text(x0, -0.55, "single copy (electric):  (2Q/r)(k·p) = (2Q/r)(E + p)\nE adds 0, p adds 1 (Sommerfeld). One power, no square.", fontsize=9.5, color=SR)
ax.set_xlim(0, 5.2); ax.set_ylim(-0.9, 3.9)
ax.set_title("A. Gravity's coupling is the electric one squared,\nand the square has three pieces")

# B: two readings of the same 6
ax = axs[1]
cases = ["electric orbit\n(single copy)", "GR orbit:\nclock-warp\nreading", "GR orbit:\ndouble-copy\nreading"]
ax.bar(0, 1, color=CROSS); ax.text(0, 0.5, "cross\nterm 1", ha="center", va="center", color="white", fontsize=9, fontweight="bold")
ax.bar(1, 2, color=SPACE); ax.bar(1, 1, bottom=2, color=SR); ax.bar(1, 3, bottom=3, color=EN)
for y, t in [(1, "space 2"), (2.5, "SR 1"), (4.5, "energy 3")]: ax.text(1, y, t, ha="center", color="white", fontsize=9.5, fontweight="bold")
ax.bar(2, 2, color=SPACE); ax.bar(2, 4, bottom=2, color=CROSS)
ax.text(2, 1, "space 2", ha="center", color="white", fontsize=9.5, fontweight="bold"); ax.text(2, 4, "cross term\n2² = 4", ha="center", color="white", fontsize=9.5, fontweight="bold")
for i, t in enumerate([1, 6, 6]): ax.text(i, t + 0.12, str(t), ha="center", fontsize=13, fontweight="bold")
ax.annotate("", xy=(1.62, 4.0), xytext=(1.38, 4.0), arrowprops=dict(arrowstyle="->", color="black"))
ax.text(1.5, 6.6, "1 + 3 = 2²", ha="center", fontsize=11, fontweight="bold")
ax.set_xticks(range(3)); ax.set_xticklabels(cases, fontsize=9.5); ax.set_ylim(0, 7.3)
ax.set_ylabel("turn per lap ÷ (π M²/L²), first order")
ax.set_title("B. Same 6, two coordinate pictures:\nboth give space 2; the double copy makes 1 + 3 one square")

# C: exact factors
ax = axs[2]
x = np.linspace(0, 0.2, 300); v2 = x / (1 - 2 * x)
ax.plot(x, np.ones_like(x), color="black", lw=2, label="E² only: T = 1 exactly")
ax.plot(x, 1 - 2 * x, color=SPACE, lw=2.2, label="E² + p²: T = 1 − 2M/r exactly")
ax.plot(x, 1 / (1 + 4 * x), color=CROSS, lw=2, ls="-.", label="E² + 2E·p: T = 1/(1 + 4M/r) exactly")
ax.plot(x, 1 - 6 * x, color=CROSS, lw=2.4, label="all three: T = 1 − 6M/r = (1 − 2M/r)(1 − (2v)²)")
ax.axvline(1 / 6, color=CROSS, lw=0.8, ls=":"); ax.text(1 / 6 + 0.002, 0.85, "r = 6M", color=CROSS, fontsize=9)
ax.axhline(0, color="black", lw=0.5); ax.set_ylim(-0.25, 1.1)
ax.set_xlabel("M/r"); ax.set_ylabel("T = (radial ÷ orbital frequency)²  (1 = closed orbit)")
ax.legend(fontsize=8.5, loc="lower left")
ax.set_title("C. Switch the pieces on one at a time\n(exact; simulated orbits agree)")
fig.suptitle("Figure O9. Gravity's 6 through string theory's double copy: (k·p)² = E² + 2E·p + p², and the 2 in 2v is the 2 of the square", y=1.03)
fig.savefig(os.path.join(OUT, "figO9-double-copy-six.png"), dpi=140, bbox_inches="tight"); print("wrote figO9")
