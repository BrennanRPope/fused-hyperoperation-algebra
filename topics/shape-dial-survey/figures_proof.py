"""Figure 16: what is proved about the Lee-Yang rule (run after proof_checks.py)."""
import json
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from matplotlib.lines import Line2D
PC = json.load(open("proof_checks.json"))
BLUE, ORANGE, GREY, RED, GREEN = "#2b6cb0", "#dd6b20", "#888888", "#c53030", "#2f855a"
plt.rcParams.update({"font.size": 11, "axes.spines.top": False, "axes.spines.right": False})
fig, (ax, bx) = plt.subplots(1, 2, figsize=(14.5, 5.4), gridspec_kw=dict(width_ratios=[1.25, 1], wspace=0.25))
P = np.linspace(1, 7, 601); amax = 9.6
# break region a > p (proved), except even p
ax.fill_between(P, P, amax, color=RED, alpha=.18, lw=0)
# gap: 1 < a <= p
ax.fill_between(P, 1, np.minimum(P, amax), color=GREY, alpha=.22, lw=0)
# proved holding: a <= 1
ax.fill_between(P, 0, 1, color=BLUE, alpha=.35, lw=0)
ax.plot(P, P, color="k", lw=1.2); ax.text(6.3, 6.65, "a = p", rotation=38, fontsize=10)
for pe in [2, 4, 6]:
    ax.plot([pe, pe], [0, amax], color=GREEN if pe == 2 else "#9ae6b4", lw=4 if pe == 2 else 3, solid_capstyle="butt", zorder=3)
ax.scatter([4]*9, np.arange(1, 10), s=46, color=GREEN, zorder=4, edgecolor="white")
# numerical checks from proof_checks
for k, (tot, real) in PC["disc"].items():
    p = float(k.split(",")[0][2:]); a = float(k.split("a=")[1].split(" ")[0])
    ax.scatter([p], [a], marker="x", s=55, color="k", zorder=5)
ax.set_xlim(1, 7); ax.set_ylim(0, 9.6); ax.set_xlabel("roundness p   (1 = cross-polytope, 2 = ball, large = cube)")
ax.set_ylabel("a = (n − 1)/p   (more spin components →)")
ax.set_title("Lee–Yang rule: what is proved")
ax.legend(handles=[Patch(color=BLUE, alpha=.35, label="holds, proved (Pólya 1918): a ≤ 1"),
                   Line2D([], [], color=GREEN, lw=4, label="holds, proved: p = 2 (Bessel), p = 4 whole a (dots)"),
                   Line2D([], [], color="#9ae6b4", lw=3, label="even p, other a: rule says holds, not proved"),
                   Patch(color=GREY, alpha=.22, label="gap: 1 < a ≤ p, holds numerically, not proved"),
                   Patch(color=RED, alpha=.18, label="breaks, proved (p not even)"),
                   Line2D([], [], marker="x", ls="", color="k", label="checked: every zero in |s| < 60 real")],
          loc="upper center", bbox_to_anchor=(0.5, -0.16), ncol=2, fontsize=8.8, frameon=False)
# right: g(x) = a p x^(p-1)(1-x^p)^(a-1) at p = 3
x = np.linspace(0, 0.999, 800); p = 3.0
for a, c, ls in [(0.5, BLUE, "-"), (1.0, BLUE, "--"), (2.0, GREY, "-"), (3.5, RED, "-")]:
    g = a*p*x**(p - 1)*(1 - x**p)**(a - 1); g = g/np.max(g[x < 0.98])
    bx.plot(x, g, color=c, ls=ls, lw=2.2, label=f"a = {a:g}" + ("  rises to the rim: Pólya applies" if a <= 1 else "  turns down: Pólya silent"))
bx.set_ylim(0, 1.6); bx.set_xlabel("x (centre → rim)"); bx.set_ylabel("g(x), scaled")
bx.set_title("p = 3: the weight in s·G(s) = ∫ g(x) sin(sx) dx")
bx.legend(fontsize=9, loc="upper left")
fig.savefig("figures/fig16-lee-yang-proof-map.png", dpi=130, bbox_inches="tight"); print("wrote fig16")
