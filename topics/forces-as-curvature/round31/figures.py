"""DRAFT round 31 picture: f38_hidden_sphere_flux.png, from results.json."""
import json, math, os
import numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch
HERE = os.path.dirname(os.path.abspath(__file__))
R = json.load(open(os.path.join(HERE, "results.json")))
GREY, BLUE, ORANGE, RED, PURPLE = "#777777", "#1f6fc5", "#e08a1e", "#c8372d", "#7b4fa8"
fig = plt.figure(figsize=(20, 13.5))
gs = fig.add_gridspec(2, 2, hspace=0.25, wspace=0.12)
fig.suptitle("Round 31 (DRAFT): does the hidden 2-sphere, and its flux, come out of the law's hidden 3-sphere?", fontsize=16)

# A: Hopf circles
ax = fig.add_subplot(gs[0, 0], projection="3d")
psi = np.linspace(0, 4 * np.pi, 400)
cmap = plt.cm.viridis
for th in [0.7, 1.4, 2.1]:
    for k, ph in enumerate(np.linspace(0, 2 * np.pi, 6)[:-1]):
        z1 = np.cos(th / 2) * np.exp(1j * (psi + ph) / 2); z2 = np.sin(th / 2) * np.exp(1j * (psi - ph) / 2)
        x1, x2, x3, x4 = z1.real, z1.imag, z2.real, z2.imag
        d = 1 - x4
        if np.min(d) < 0.05: continue
        ax.plot(x1 / d, x2 / d, x3 / d, color=cmap(th / 3.0), lw=1.6)
ax.set_box_aspect((1, 1, 1)); ax.set_xlim(-2.6, 2.6); ax.set_ylim(-2.6, 2.6); ax.set_zlim(-2.6, 2.6)
ax.set_xticks([]); ax.set_yticks([]); ax.set_zticks([])
ax.set_title("A. The law's hidden 3-sphere is made of linked circles (H1)\n(drawn by stereographic projection; colour = which point of the 2-sphere below)")

# B: base sphere with flux
ax = fig.add_subplot(gs[0, 1], projection="3d")
u, v = np.mgrid[0:2 * np.pi:40j, 0:np.pi:20j]
ax.plot_surface(np.cos(u) * np.sin(v), np.sin(u) * np.sin(v), np.cos(v), color=BLUE, alpha=0.25, linewidth=0)
for th in np.linspace(0.3, np.pi - 0.3, 5):
    for ph in np.linspace(0, 2 * np.pi, 9)[:-1]:
        n = np.array([np.sin(th) * np.cos(ph), np.sin(th) * np.sin(ph), np.cos(th)])
        ax.quiver(*n, *(0.45 * n), color=ORANGE, lw=1.4, arrow_length_ratio=0.3)
for p in ([0, 0, 1], [0, 0, -1]):
    ax.quiver(*p, *(0.45 * np.array(p)), color=ORANGE, lw=1.4, arrow_length_ratio=0.3)
ax.set_box_aspect((1, 1, 1)); ax.set_xticks([]); ax.set_yticks([]); ax.set_zticks([])
ax.set_xlim(-1.5, 1.5); ax.set_ylim(-1.5, 1.5); ax.set_zlim(-1.5, 1.5)
ax.set_title("B. Each circle is one point of a 2-sphere (radius 1/(2√q), the law's sphere at dial 4q);\nthe way the circles twist round it is exactly one unit of magnetic flux (H1, H2)")

# C: zero modes per sector
ax = fig.add_subplot(gs[1, 0]); rows = R["H3"]
fl = np.array([r["flux"] for r in rows]); L = np.array([r["left"] for r in rows]); Rr = np.array([r["right"] for r in rows])
ax.bar(fl, L, color=BLUE, width=0.6, label="left-handed massless modes")
ax.bar(fl, -Rr, color=ORANGE, width=0.6, label="right-handed massless modes")
for r in rows:
    if r["flux"] == 0: continue
    sp_ = (r["spin_left"] or r["spin_right"])[0]
    name = {0.0: "single", 0.5: "doublet", 1.0: "triplet", 1.5: "quartet"}[sp_]
    y = r["left"] + 0.25 if r["left"] else -r["right"] - 0.25
    ax.text(r["flux"], y, name, ha="center", va="bottom" if r["left"] else "top", fontsize=10)
ax.axhline(0, color="k", lw=1)
for x, lab, xy in [(1, "round 10 (charge ½)", (-1.2, 2.4)), (2, "round 9 (charge 1)", (-0.6, 3.5))]:
    ax.annotate(lab, (x - 0.25, x + 0.05), xy, fontsize=10, color=BLUE, ha="center",
                arrowprops=dict(arrowstyle="->", color=BLUE))
ax.set_xticks(fl); ax.set_ylim(-5.2, 5.4)
ax.set_xlabel("flux seen by the matter, n = 2 × its weight along the circles")
ax.set_ylabel("massless modes (left up, right down)")
ax.set_title("C. Every sector of the hidden 3-sphere: hands come in mirror pairs (H3, H4)")
ax.legend(loc="lower right", fontsize=10)
ax.text(-4.3, 4.6, f"sum over all sectors: left − right = {R['H4']['net_all']}", fontsize=11, color=RED)

# D: chain
ax = fig.add_subplot(gs[1, 1]); ax.axis("off"); ax.set_xlim(0, 10); ax.set_ylim(0, 10)
ax.set_title("D. What follows from the law and what is still a choice")
boxes = [(5, 9.0, "the weak force's hidden 3-sphere (rounds 4 and 7) = the law's 3-sphere", GREY),
         (5, 7.0, "read along its circles: a 2-sphere with one unit of flux\nround 9's assumption (a) now follows (H1, H2)", BLUE),
         (5, 4.9, "geometry carries every sector with its mirror:\nnet handedness 0 (Witten 1981) (H4)", RED),
         (5, 2.8, "one hand only if matter keeps one sign of weight:\nround 10's bookkeeping, suggested by the law, not forced (assumption (b) stays)", ORANGE),
         (5, 0.9, "left over: SU(2) × U(1), the weak and hypercharge pattern; flux 2 gives doublets", PURPLE)]
for xx, yy, txt, col in boxes:
    ax.add_patch(FancyBboxPatch((xx - 4.6, yy - 0.75), 9.2, 1.5, boxstyle="round,pad=0.08", fc="white", ec=col, lw=2))
    ax.text(xx, yy, txt, ha="center", va="center", fontsize=11, color=col)
for y0, y1 in [(8.25, 7.75), (6.25, 5.65), (4.15, 3.55), (2.05, 1.65)]:
    ax.annotate("", (5, y1), (5, y0), arrowprops=dict(arrowstyle="->", color=GREY, lw=2))
os.makedirs(os.path.join(HERE, "figures"), exist_ok=True)
fig.savefig(os.path.join(HERE, "figures", "f38_hidden_sphere_flux.png"), dpi=120, bbox_inches="tight")
