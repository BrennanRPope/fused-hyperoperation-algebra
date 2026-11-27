"""DRAFT picture for round 14. Run after checks.py. Blue = the law, orange = Regge read-back, grey = reference."""
import json, math
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon

OR, BL, GR = "#e07b00", "#1f6fc2", "#777777"
R = json.load(open("results.json"))
plt.rcParams.update({"font.size": 10.5, "figure.dpi": 130})
fig, ax = plt.subplots(1, 3, figsize=(16.5, 5.2))

# A. slabs of spacetime, drawn in one space direction and time
a = ax[0]
sizes = [1.0, 1.12, 1.32, 1.62]; ts = [0, 1, 2, 3]
xs = np.array([-1.5, -0.5, 0.5, 1.5]) / 1.5          # corners of the space pieces at size 1
for k in range(3):
    for j in range(3):
        quad = [(sizes[k] * xs[j], ts[k]), (sizes[k] * xs[j + 1], ts[k]), (sizes[k + 1] * xs[j + 1], ts[k + 1]), (sizes[k + 1] * xs[j], ts[k + 1])]
        a.add_patch(Polygon(quad, closed=True, fc=("#fde3c2" if (j + k) % 2 else "#fbd29f"), ec=OR, lw=1.2))
for k in range(4):
    a.plot(sizes[k] * xs, [ts[k]] * 4, "o-", color="k", ms=4, lw=2)
a.annotate("space piece at one moment\n(size a)", (-0.33, 0), (-0.8, -0.6), ha="center", fontsize=8.5, arrowprops=dict(arrowstyle="-", color=GR))
a.annotate("strut: lasts proper\ntime τ", (sizes[1] * xs[3] * 0.5 + sizes[2] * xs[3] * 0.5, 1.5), (1.9, 1.75), fontsize=8.5,
           arrowprops=dict(arrowstyle="->", color=GR))
a.annotate("angle hinge: a space\nedge swept through time", (sizes[2] * xs[1] * 0.5 + sizes[3] * xs[1] * 0.5, 2.5), (-2.75, 1.2), fontsize=8.5,
           color=OR, arrowprops=dict(arrowstyle="->", color=OR))
a.annotate("boost hinge: a space\ntriangle at one moment", (sizes[1] * xs[2], 1), (1.3, 0.25), fontsize=8.5, color=OR,
           arrowprops=dict(arrowstyle="->", color=OR))
a.text(0, 3.45, "each slab is flat Minkowski spacetime;\ncurvature sits only on the hinges", ha="center", fontsize=8.5)
a.set_xlim(-2.9, 3.0); a.set_ylim(-0.8, 3.95); a.set_xlabel("space"); a.set_ylabel("time")
a.set_xticks([]); a.set_yticks([])
a.set_title("A. Glue flat slabs of spacetime:\na space piece now, the same piece a moment later")

# B. the two terms of Einstein's action, read back
b = ax[1]
lv = np.arange(4)
st = np.abs(np.array(R["static_ratio"]) - 1); kn = np.abs(np.array([r[1] for r in R["kinetic_ratio"]]) - 1)
b.semilogy(lv, st, "o-", color=OR, ms=8, label="slab at rest: space curvature  6π² a τ")
b.semilogy(lv, kn, "s-", color=OR, ms=8, mfc="white", label="growing slab: kinetic term  −6π² a (Δa)²/τ")
b.semilogy(lv, [abs(R["throat"][l] - 1) for l in lv], "^-", color=BL, ms=7, label="with Λ: the throat a² = 3/Λ")
b.semilogy(lv, kn[0] / 4.0**lv, ":", color=GR, label="second order (÷4 per level)")
b.set_xticks(lv); b.set_xticklabels([f"{l}\n({n:,} pieces)" for l, n in zip(lv, [600, 4800, 38400, 307200])])
b.set_xlabel("halvings of the space pieces (law's half kicks)"); b.set_ylabel("relative error vs Einstein")
b.legend(fontsize=8, loc="lower left")
b.set_title("B. Read back: both parts of Einstein's action,\nthe kinetic one from the boosts, with its minus sign")

# C. a de Sitter universe grown slab by slab
c = ax[2]
tt = np.linspace(0, 2.4, 300)
c.plot(tt, np.cosh(tt), color=GR, lw=2, label="de Sitter: a = cosh t  (Λ = 3)")
r3 = R["run_level3_tau0.025"]; c.plot(r3["t"][::2], r3["a"][::2], "o", color=OR, ms=5, label="glued slabs, level 3 (τ = 0.025)")
h1 = R["horizon_run_level1"]
c.plot(h1["t"], h1["a"], "^", color=OR, mfc="white", ms=6, label="glued slabs, level 1 (coarse, τ = 0.1)")
po = np.array(h1["piece_over_horizon"]); i = int(np.argmax(po >= 1.0))
c.annotate("level 1 pieces reach the\nhorizon size here, and\nthe slabs stop following", (h1["t"][i], h1["a"][i]), (0.55, 4.3), fontsize=8.5,
           arrowprops=dict(arrowstyle="->", color=GR))
c.set_xlabel("proper time t (in units of √(3/Λ))"); c.set_ylabel("size of the universe a")
c.legend(fontsize=8, loc="upper left"); c.set_ylim(0.8, 6)
c.set_title("C. Grow the universe slab by slab:\nit follows de Sitter while pieces are inside the horizon")
fig.tight_layout(); fig.savefig("figures/f21_spacetime_slabs_with_time.png", bbox_inches="tight")
print("f21 written")
