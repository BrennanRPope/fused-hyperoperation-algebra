"""Figure 26: what the Lean file LeeYangNearZone.lean proves, and what it takes as named inputs."""
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch
BLUE, ORANGE, GREY = "#2b6cb0", "#dd6b20", "#888888"
plt.rcParams.update({"font.size": 11})
fig = plt.figure(figsize=(15.5, 10))
gs = fig.add_gridspec(1, 2, width_ratios=[1.55, 1], wspace=0.08)

ax = fig.add_subplot(gs[0, 0]); ax.set_xlim(0, 10); ax.set_ylim(0, 10); ax.axis("off")
def box(x, y, w, h, title, body, kind):
    col = {"lean": BLUE, "cite": ORANGE, "core": "#2f855a"}[kind]
    ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.08", fc=col, alpha=0.13, ec=col, lw=2))
    ax.text(x + w/2, y + h - 0.22, title, ha="center", va="top", fontsize=11.5, fontweight="bold", color=col)
    ax.text(x + w/2, y + h - 0.62, body, ha="center", va="top", fontsize=9.6)
def arrow(x0, y0, x1, y1):
    ax.annotate("", (x1, y1), (x0, y0), arrowprops=dict(arrowstyle="-|>", lw=1.6, color="k"))
box(0.1, 7.6, 3.0, 2.2, "NearCert (cited)", "near_zone_p.py: for each box,\nthe runs of cells, their end signs,\nslopes and the next step's signs", "cite")
box(3.5, 7.6, 3.0, 2.2, "Taylor cell bound", "taylor_cell: a bound on the\n(n+1)-st derivative bounds G\non a whole cell", "lean")
box(6.9, 7.6, 3.0, 2.2, "FarCert (cited)", "far_bounds_p.py: past T*(p),\nG at b + 1 alternates at the\nzeros of G at b", "cite")
box(0.6, 4.6, 4.2, 2.3, "Run certificate → alternation", "RunCert.nearAlt: one zero per run\n(intermediate value + strict monotone),\nno zeros in the gaps, so the next\nstep's signs flip from zero to zero", "lean")
box(5.2, 4.6, 4.2, 2.3, "Near + far → Alternates p b", "alternates_of_near_far: a zero of G\npast T*(p) inside the near zone glues\nthe two halves, so every positive\nzero is covered", "lean")
box(0.1, 1.5, 3.0, 2.4, "Toolkit (cited)", "Pólya: base a ≤ 1\nKrein: the climb b → b + 1\nbreak past a = p\n(from the shared core)", "cite")
box(3.5, 1.5, 3.0, 2.4, "Core theorem", "holds_of_alternatesBelow\n(LeeYangCore, already built):\nbase + climb on every step", "core")
box(6.9, 1.5, 3.0, 2.4, "rule_range", "Holds p a ⇔ a ≤ p\nfor every 1.1 ≤ p ≤ 5\nnot an even whole number;\nholds_range: Lee–Yang\nfor every a ≤ p", "lean")
arrow(1.6, 7.6, 2.2, 6.95); arrow(5.0, 7.6, 3.6, 6.95); arrow(8.4, 7.6, 7.8, 6.95)
arrow(4.8, 5.75, 5.2, 5.75); arrow(7.3, 4.6, 6.0, 3.95); arrow(3.1, 2.7, 3.5, 2.7); arrow(6.5, 2.7, 6.9, 2.7)
ax.text(0.1, 0.75, "Blue: proved in Lean against Mathlib, no sorry (LeeYangNearZone.lean, required in CI).  Green: proved earlier.\n"
        "Orange: named inputs. The two certificates are the Python runs; the toolkit is the cited literature.", fontsize=10, color=GREY)
ax.set_title("What the Lean file proves, and what it takes as inputs", fontsize=13)

# right: the Taylor cell bound on a real cell
ax = fig.add_subplot(gs[0, 1])
xg, wg = np.polynomial.legendre.leggauss(400); xg = (xg + 1)/2; wg = wg/2
def G(p, c, s): return np.array([np.sum(wg*(1 - xg**p)**c*np.cos(si*xg)) for si in np.atleast_1d(s)])
p, c, x0, h, n = 3.0, 2.0, 18.0, 0.125, 3
t = np.linspace(x0, x0 + h, 200)
g = G(p, c, t); g0 = G(p, c, x0)[0]
eps = 1e-6
d = [(G(p, c, x0 + k*eps)[0]) for k in range(6)]
# derivatives at x0 from the integral directly
der = [np.sum(wg*(1 - xg**p)**c*(xg**m)*np.cos(x0*xg + m*np.pi/2)) for m in range(1, n + 1)]
M = 1/(n + 2)          # |G^(n+1)| <= 1/(n+2) since int x^(n+1) dx = 1/(n+2) and the weight is at most 1
u = t - x0
import math
bound = sum(abs(der[m])*u**(m + 1)/math.factorial(m + 1) for m in range(n)) + M*u**(n + 1)/math.factorial(n + 1)
ax.fill_between(t, g0 - bound, g0 + bound, color=BLUE, alpha=0.18, label="Lean's bound: |G(t) − G(x₀)| ≤ Taylor terms + M u⁴/4!")
ax.plot(t, g, color="k", lw=2, label="G at p = 3, b = 2")
ax.axhline(0, color="#c53030", lw=1.2, ls="--")
ax.text(x0 + h*0.98, 0.00004, "G = 0", ha="right", va="bottom", color="#c53030", fontsize=10)
ax.set_ylim(-0.0003, 0.0021)
ax.set_xlabel("s, across one cell of width 1/8"); ax.set_ylabel("G")
ax.legend(fontsize=9.5, loc="upper right")
ax.set_title("The cell bound in action: the curve stays inside\nthe band, which stays above zero, so G > 0\non the whole cell without sampling it", fontsize=12)
for sp in ("top", "right"): ax.spines[sp].set_visible(False)
fig.suptitle("Figure 26. The near-zone reasoning, proved in Lean (DRAFT)", fontsize=14, y=0.99)
fig.savefig("figures/fig26-lee-yang-near-lean.png", dpi=110, bbox_inches="tight"); print("saved")
