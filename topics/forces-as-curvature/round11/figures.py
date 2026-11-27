"""DRAFT picture for round 11. Run after checks.py. Blue = the law, orange = glued flat pieces (Regge), grey = reference."""
import json, math
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon, FancyArrowPatch

OR, BL, GR = "#e07b00", "#1f6fc2", "#777777"
R = json.load(open("results.json"))
plt.rcParams.update({"font.size": 10.5, "figure.dpi": 130})
fig, ax = plt.subplots(1, 3, figsize=(16, 5))

# A: five flat triangles round a corner of an icosahedron, laid flat: a 60 degree gap is left (the deficit)
a = ax[0]
for i in range(5):
    t0, t1 = math.radians(60 * i), math.radians(60 * (i + 1))
    a.add_patch(Polygon([(0, 0), (math.cos(t0), math.sin(t0)), (math.cos(t1), math.sin(t1))], closed=True,
                        fc="#fbe3c7", ec=OR, lw=1.5))
a.plot([0, 1], [0, 0], color=OR, lw=1.5); a.plot([0, math.cos(math.radians(300))], [0, math.sin(math.radians(300))], color=OR, lw=1.5, ls="--")
arc = np.radians(np.linspace(300, 360, 30)); a.plot(0.35 * np.cos(arc), 0.35 * np.sin(arc), color=GR)
a.text(0.42, -0.28, "gap = deficit\n(60° here)", color=GR, fontsize=9)
# laid flat, carrying an arrow across the pieces never turns it; gluing the gap shut turns the last piece by 60 degrees
for ang, lab, dy in ((30, "arrow starts here", 0.13), (270, "same arrow, carried round;\nclosing the gap turns it 60°", -0.3)):
    t = math.radians(ang); p = np.array([0.55 * math.cos(t), 0.55 * math.sin(t)]) - np.array([0.12, 0])
    a.add_patch(FancyArrowPatch(p, p + np.array([0.25, 0.0]), arrowstyle="-|>", mutation_scale=14, color=BL, lw=2))
    a.text(p[0] - 0.15, p[1] + dy, lab, color=BL, fontsize=8.5)
a.set_xlim(-1.15, 1.25); a.set_ylim(-1.15, 1.15); a.set_aspect("equal"); a.axis("off")
a.set_title("A. Flat pieces round a corner, laid flat:\nthe curvature is the gap a loop comes back turned by")

# B: 2D sphere: deficits vs the law's twist
X = R["sphere2d"]; lev = [r[0] for r in X]
b = ax[1]
b.semilogy(lev, [r[5] for r in X], "o-", color=BL, label="any patch: deficits vs law's twist (falls 4× per level)")
b.semilogy(lev, [r[7] for r in X], "s-", color=OR, label="worst six-fold corner (levels off near 0.3%)")
b.semilogy(lev, [r[6] for r in X], "^-", color=GR, label="the 12 five-fold corners (stay about 14.6% off)")
b.set_xticks(lev); b.set_xlabel("refinement level of the glued sphere (each level halves the pieces)")
b.set_ylabel("relative mismatch"); b.legend(fontsize=8, loc="lower left")
b.set_title("B. 2-sphere: the corners' deficits match the law's\ntwist over any patch, though not corner by corner")

# C: 3-sphere: Regge's action and the cosmological-constant size
Y = R["sphere3d"]; lv = [r[0] for r in Y]
target = 6 * math.pi**2
c = ax[2]
c.semilogy(lv, [abs(r[3] / target - 1) for r in Y], "o-", color=OR, label="Regge's action vs ½∫R = 6π²")
c.semilogy(lv, [abs(x - 1) for x in R["lambda_rho2"]], "s-", color=BL, label="with Λ: best size gives Λρ² → 1")
c.semilogy(lv, [abs(Y[0][3] / target - 1) / 4**l for l in lv], ":", color=GR, label="second order (÷4 per level)")
c.set_xticks(lv); c.set_xticklabels([f"{l}\n({r[1]:,} tets)" for l, r in zip(lv, Y)])
c.set_xlabel("refinement of the 600-cell"); c.set_ylabel("relative error"); c.legend(fontsize=8)
c.set_title("C. 3-sphere glued from flat tetrahedra:\nRegge's action tends to Einstein's")
fig.tight_layout(); fig.savefig("figures/f18_spacetime_from_flat_pieces.png", bbox_inches="tight")
print("f18 written")
