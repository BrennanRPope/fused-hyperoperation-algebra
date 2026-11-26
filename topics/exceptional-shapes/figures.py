"""Figures for the exceptional-shapes pass. Run: python3 figures.py"""
import itertools
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch
from shapes import *

OR, BL, RD, GR = "#e07b00", "#1f6fb5", "#c0392b", "#9a9a9a"
plt.rcParams.update({"font.size": 10, "axes.titlesize": 11, "figure.dpi": 130})
out = "figures/"
T2 = binary_tetrahedral(); O2 = binary_octahedral(); I2 = binary_icosahedral(); Q8 = cross_polytope4()
C120 = onetwenty_cell(); R8 = e8_roots()


def f4_roots():
    return np.vstack([np.sqrt(2) * dual24(), Q8, tesseract_half()])


def edges(V, d=None, tol=1e-6):
    D = np.linalg.norm(V[:, None] - V[None], axis=2)
    if d is None:
        d = D[D > tol].min()
    return [(i, j) for i in range(len(V)) for j in range(i + 1, len(V)) if abs(D[i, j] - d) < tol]


def draw(ax, V, P, col, lw=0.4, ms=8, alpha=0.5, E=None):
    X = V @ P.T
    for i, j in (E if E is not None else edges(V)):
        ax.plot(*X[[i, j]].T, color=col, lw=lw, alpha=alpha)
    ax.scatter(*X.T, s=ms, color=col, zorder=3)
    ax.set_aspect("equal"); ax.axis("off")


PF4, hF4, _ = coxeter_plane(f4_roots())
PH4, hH4, _ = coxeter_plane(I2)
PE8, hE8, _ = coxeter_plane(R8)

# ---- fig2: the shapes ----
fig, axs = plt.subplots(1, 4, figsize=(15, 4.2))
draw(axs[0], T2, PF4, BL, lw=0.8, ms=18, alpha=0.7)
axs[0].set_title(f"24-cell: 24 corners, 24 octahedra\n(Coxeter plane, 12-fold)")
draw(axs[1], I2, PH4, BL, lw=0.35, ms=8)
axs[1].set_title("600-cell: 120 corners, 600 tetrahedra\n(Coxeter plane, 30-fold)")
draw(axs[2], C120, PH4, OR, lw=0.25, ms=3)
axs[2].set_title("120-cell: 600 corners, 120 dodecahedra\n(dual of the 600-cell)")
draw(axs[3], R8, PE8, "#555555", lw=0.08, ms=4, alpha=0.25)
axs[3].set_title("E8: 240 roots in 8D\n(Petrie picture, 8 rings of 30)")
fig.suptitle("Figure 2. The four missing shapes, each drawn in its own most symmetric flat view", y=1.02)
fig.savefig(out + "fig2-the-four-missing-shapes.png", bbox_inches="tight"); plt.close(fig)

# ---- fig3: nesting of corner sets under the core law ----
fig, ax = plt.subplots(figsize=(12, 6.2)); ax.axis("off"); ax.set_xlim(0, 12); ax.set_ylim(0, 6.6)
def box(x, y, w, h, title, sub, col):
    ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.05", fc="white", ec=col, lw=2))
    ax.text(x + w / 2, y + h - 0.28, title, ha="center", va="top", fontsize=11, weight="bold", color=col)
    ax.text(x + w / 2, y + 0.18, sub, ha="center", va="bottom", fontsize=8.5)
def mini(x, y, V, P, col, s=0.55, lw=0.3):
    X = V @ P.T; X = X / np.abs(X).max() * s + [x, y]
    for i, j in edges(V):
        ax.plot(*X[[i, j]].T, color=col, lw=lw, alpha=0.6)
    ax.scatter(*X.T, s=3, color=col, zorder=3)
box(0.2, 0.3, 2.4, 2.5, "cross-polytope", "8 corners: Paulis\n±1, ±i, ±j, ±k", OR); mini(1.4, 1.75, Q8, PF4, OR)
box(0.2, 3.5, 2.4, 2.5, "tesseract", "16 corners\n(±1±i±j±k)/2", OR); mini(1.4, 4.95, tesseract_half(), PF4, OR)
box(3.4, 1.9, 2.6, 2.6, "24-cell", "24 corners: a group\ncross-polytope + tesseract", BL); mini(4.7, 3.35, T2, PF4, BL)
box(6.8, 3.6, 2.4, 2.6, "Clifford 48", "24-cell + dual 24-cell\na group (Q# page)", OR); mini(8.0, 5.05, O2, PF4, OR, lw=0.15)
box(6.8, 0.2, 2.4, 2.9, "600-cell", "120 corners: a group\nholds five 24-cells", BL); mini(8.0, 1.85, I2, PH4, BL, lw=0.15)
box(9.7, 0.2, 2.2, 2.9, "120-cell", "600 corners: NOT a group\n= five 600-cells", RD); mini(10.8, 1.85, C120, PH4, RD, s=0.5, lw=0.1)
for (a, b, txt) in (((2.6, 1.6), (3.4, 2.8), "+"), ((2.6, 4.7), (3.4, 3.8), "+"), ((6.0, 3.7), (6.8, 4.6), "+ dual"), ((6.0, 2.8), (6.8, 1.9), "×5"), ((9.2, 1.6), (9.7, 1.6), "×5")):
    ax.annotate("", xy=b, xytext=a, arrowprops=dict(arrowstyle="->", color=GR, lw=1.5))
    ax.text((a[0] + b[0]) / 2, (a[1] + b[1]) / 2 + 0.15, txt, color=GR, fontsize=9, ha="center")
ax.text(6, 6.45, "Corner sets that are closed under our law x ∘ₖ y = x + y + kxy (lift u = 1 + kx, every k ≠ 0)", ha="center", fontsize=11)
ax.text(0.2, -0.15, "Blue = closed and new to the catalog; orange = already in the catalog; red = not closed (only the 600-cell's corners act on it).", fontsize=8.5, color="#333")
fig.suptitle("Figure 3. How the corner sets nest", y=0.98)
fig.savefig(out + "fig3-corner-sets-nest.png", bbox_inches="tight"); plt.close(fig)

# ---- fig4: swap ladder ----
fig, axs = plt.subplots(1, 2, figsize=(13, 5), gridspec_kw={"width_ratios": [3, 1.3]})
ax = axs[0]; ax.axis("off"); ax.set_xlim(-0.5, 9.5); ax.set_ylim(-1.8, 2.2)
lad = [("Clifford 48", O2, PF4, OR, "octahedral\n24 turns (S4)"), ("24-cell", T2, PF4, BL, "tetrahedral\n12 turns (A4)"),
       ("cross-polytope", Q8, PF4, OR, "Klein four\n4 turns"), ("±1", np.array([[1, 0, 0, 0], [-1, 0, 0, 0]]), PF4, GR, "no turn\n1")]
for n, (nm, V, P, col, sub) in enumerate(lad):
    x = n * 2.6
    X = V @ P.T if len(V) > 2 else np.array([[-0.5, 0], [0.5, 0]])
    X = X / max(np.abs(X).max(), 1e-9) * 0.8 + [x, 0.5]
    if len(V) > 2:
        for i, j in edges(V): ax.plot(*X[[i, j]].T, color=col, lw=0.4, alpha=0.6)
    ax.scatter(*X.T, s=10, color=col, zorder=3)
    ax.text(x, 1.6, f"{nm}\n{len(V)} corners", ha="center", fontsize=10, weight="bold", color=col)
    ax.text(x, -0.75, "3D shadow:\n" + sub, ha="center", fontsize=8.5)
    if n < 3:
        ax.annotate("", xy=(x + 1.75, 0.5), xytext=(x + 0.85, 0.5), arrowprops=dict(arrowstyle="->", lw=2, color="k"))
        ax.text(x + 1.3, 0.75, "swap", ha="center", fontsize=9)
ax.text(4, -1.7, "Swapping two corners (a b a⁻¹ b⁻¹) of each shape lands on the next shape; after three swaps nothing is left.\nKnown: this is why the quartic equation can be solved by roots (S4 is solvable).", ha="center", fontsize=9)
ax = axs[1]; ax.axis("off"); ax.set_xlim(-1.3, 1.3); ax.set_ylim(-1.8, 2.2)
X = I2 @ PH4.T; X = X / np.abs(X).max() * 0.8 + [0, 0.5]
for i, j in edges(I2): ax.plot(*X[[i, j]].T, color=BL, lw=0.2, alpha=0.6)
ax.scatter(*X.T, s=4, color=BL, zorder=3)
ax.text(0, 1.6, "600-cell\n120 corners", ha="center", weight="bold", color=BL)
ax.annotate("", xy=(0.95, 0.75), xytext=(0.95, 0.25), arrowprops=dict(arrowstyle="->", lw=2, connectionstyle="arc3,rad=-1.4"))
ax.text(1.25, 0.5, "swap", fontsize=9, rotation=90, va="center")
ax.text(0, -1.0, "Every corner is itself one swap:\nswapping never shrinks it.\nKnown: icosahedral A5 is simple, which is\nwhy the quintic can't be solved by roots (Klein).", ha="center", fontsize=8.5)
fig.suptitle("Figure 4. The swap rule on the new corner sets: the Clifford shapes run down to nothing, the 600-cell stays whole", y=1.0)
fig.savefig(out + "fig4-swap-ladder.png", bbox_inches="tight"); plt.close(fig)

# ---- fig5: Hopf shadows ----
def hopf(G, deg):
    for g in G:
        if abs(rot_angle_deg(g) - deg) < 1e-6:
            n = np.array([0, *(g[1:] / np.linalg.norm(g[1:]))]); break
    imgs = uniq([qmul(qmul(g, n), qconj(g))[1:] for g in G], nd=5)
    # rotate so that n -> z
    z = n[1:]; a = np.cross(z, [0, 0, 1]); s = np.linalg.norm(a); c = z[2]
    if s > 1e-9:
        a /= s; K = np.array([[0, -a[2], a[1]], [a[2], 0, -a[0]], [-a[1], a[0], 0]]); Rm = np.eye(3) + s * K + (1 - c) * K @ K
        imgs = imgs @ Rm.T
    return imgs
fig = plt.figure(figsize=(15, 4.4))
cases = [(I2, 72, "600-cell, 10 corners per fibre\n→ 12 points: icosahedron", BL), (I2, 120, "600-cell, 6 per fibre\n→ 20 points: dodecahedron", BL),
         (I2, 180, "600-cell, 4 per fibre\n→ 30 points: icosidodecahedron", BL), (T2, 180, "24-cell, 4 per fibre\n→ 6 points: octahedron", OR), (T2, 120, "24-cell, 6 per fibre\n→ 4 points: tetrahedron", OR)]
for n, (G, deg, ttl, col) in enumerate(cases):
    ax = fig.add_subplot(1, 5, n + 1, projection="3d")
    u, v = np.mgrid[0:2 * np.pi:30j, 0:np.pi:15j]
    ax.plot_wireframe(np.cos(u) * np.sin(v), np.sin(u) * np.sin(v), np.cos(v), color="#dddddd", lw=0.3)
    P = hopf(G, deg); ax.scatter(*P.T, s=30, color=col, depthshade=False)
    D = P @ P.T; m = np.max(D - 2 * np.eye(len(P)))
    for i, j in itertools.combinations(range(len(P)), 2):
        if abs(D[i, j] - m) < 1e-6: ax.plot(*P[[i, j]].T, color=col, lw=0.8)
    ax.set_title(ttl, fontsize=9); ax.set_axis_off(); ax.set_box_aspect((1, 1, 1))
fig.suptitle("Figure 5. Polar coordinates on S³ (a point of the sphere S² plus a phase): each corner set's shadow on S² is a 3D solid", y=1.0)
fig.savefig(out + "fig5-hopf-shadows.png", bbox_inches="tight"); plt.close(fig)

# ---- fig6: p dial ----
fig, ax = plt.subplots(figsize=(10, 4.6))
p = np.concatenate([np.linspace(1, 2, 200), np.linspace(2, 60, 400)])
ax.plot(p, 4.0 ** (-np.abs(0.5 - 1 / p)), color=GR, lw=2, label="ℓᵖ balls in 4D (the p dial)")
ax.set_xscale("log"); ax.set_xticks([1, 4 / 3, 1.8, 2.25, 4, 10, 60]); ax.set_xticklabels(["1", "4/3", "1.80", "2.25", "4", "10", "∞"])
pts = [(1, 0.5, "cross-polytope", OR, "left"), (60, 0.5, "tesseract", OR, "right"), (4 / 3, 2 ** -0.5, "24-cell", BL, "left"), (4, 2 ** -0.5, "24-cell", BL, "left"),
       (1.7993, PHI ** 2 / 8 ** .5, "600-cell", BL, "right"), (2.2515, PHI ** 2 / 8 ** .5, "120-cell", BL, "left"), (2, 1, "sphere (p = 2)", OR, "center")]
for x, y, t, c, ha in pts:
    ax.scatter([x], [y], s=70, color=c, zorder=3)
    ax.annotate(t, (x, y), xytext=(0, 9 if t != "sphere" else 8), textcoords="offset points", ha=ha if ha != "left" else "center", color=c, fontsize=9)
ax.plot([4 / 3, 4], [2 ** -.5] * 2, ls=":", color=BL); ax.text(2.3, 0.655, "self-dual: one ratio, two mirror spots (1/p + 1/q = 1)", fontsize=8, color=BL, ha="center")
ax.annotate("simplex: ratio 1/4,\nbelow every ℓᵖ ball\n(off this dial; it has its own)", (1.05, 0.27), fontsize=8.5, color=RD)
ax.scatter([1.0], [0.25], color=RD, marker="v", s=60)
ax.set_ylim(0.2, 1.1); ax.set_ylabel("inner radius ÷ outer radius"); ax.set_xlabel("p")
ax.axvline(2, color="#cccccc", lw=0.8); ax.text(1.12, 0.40, "← cornered side\n(corners poke out)", fontsize=8); ax.text(12, 0.40, "boxy side →\n(faces poke out)", fontsize=8)
ax.set_title("Figure 6. Where each shape sits on the p dial, matched by one number (how round it is)")
fig.savefig(out + "fig6-p-dial-positions.png", bbox_inches="tight"); plt.close(fig)

# ---- fig7: gates ----
fig, axs = plt.subplots(1, 2, figsize=(13, 5))
ax = axs[0]; ax.set_aspect("equal"); ax.axis("off")
th = np.linspace(0, np.pi, 200); ax.plot(np.cos(th), np.sin(th), color="#cccccc")
for a, lab in ((0, "0"), (45, "45°"), (72, "72°"), (90, "90°"), (120, "120°"), (144, "144°"), (180, "180°")):
    t = np.radians(180 - a); ax.text(1.12 * np.cos(t), 1.12 * np.sin(t), lab, ha="center", va="center", fontsize=8)
for a in (0, 72, 120, 144, 180):
    t = np.radians(180 - a); ax.plot([0, np.cos(t)], [0, np.sin(t)], color=BL, lw=3, alpha=0.6)
for a in (0, 90, 120, 180):
    t = np.radians(180 - a); ax.plot([0, 0.75 * np.cos(t)], [0, 0.75 * np.sin(t)], color=OR, lw=3, alpha=0.8)
for a, lab, r in ((180, "X, Y, Z, H", 1.0), (90, "S, √X", 1.3), (45, "T", 1.3)):
    t = np.radians(180 - a); ax.text(r * np.cos(t), r * np.sin(t) + (0.12 if a == 180 else 0), lab, ha="center", fontsize=10, weight="bold", color=RD if a == 45 else "k")
ax.text(0, -0.25, "Turn angles a corner can have.\nBlue (long): 600-cell.  Orange (short): Clifford 48.\nT turns 45°: neither shape has it, in any orientation.\nS and √X (90°) are Clifford only; H (180°) fits both\nbut never together with X, Y, Z in one 600-cell.", ha="center", va="top", fontsize=9)
ax.set_ylim(-1.3, 1.5)
ax = axs[1]
counts = [60, 1848, 50244]
ax.bar(["600-cell\nalone", "with one T", "with two T's"], counts, color=[BL, GR, GR]); ax.set_yscale("log")
for i, c in enumerate(counts): ax.text(i, c * 1.15, f"{c:,}", ha="center")
ax.set_ylabel("different gates (up to phase)"); ax.set_title("A finite corner set plus one off-corner gate keeps growing\n(known: 'golden gates', Parzanchevski–Sarnak 2018)", fontsize=10)
fig.suptitle("Figure 7. Quantum gates on the 600-cell", y=1.0)
fig.savefig(out + "fig7-gates-on-the-600-cell.png", bbox_inches="tight"); plt.close(fig)

# ---- fig8: regular polytopes by dimension ----
fig, ax = plt.subplots(figsize=(11, 4.8))
dims = np.arange(1, 9); fam = [1, 0, 3, 3, 3, 3, 3, 3]; extra = [0, 0, 2, 3, 0, 0, 0, 0]
ax.bar(dims[2:], fam[2:], color=GR, zorder=2, label="simplex, cube, cross-polytope (every dimension)")
ax.bar(dims[2:], extra[2:], bottom=fam[2:], color=BL, zorder=2, label="extra regular shapes")
ax.bar([2], [6.5], color=BL, hatch="//", alpha=0.5, zorder=2); ax.text(2, 6.6, "∞\n(every polygon)", ha="center", fontsize=8)
ax.text(3, 5.15, "icosahedron\ndodecahedron", ha="center", fontsize=7.5); ax.text(4, 6.15, "24-cell\n600-cell, 120-cell", ha="center", fontsize=7.5)
for d in (1, 2, 4, 8):
    ax.axvspan(d - 0.45, d + 0.45, color=OR, alpha=0.12, zorder=0)
ax.text(8, 3.4, "E8 lattice\n(not a regular\npolytope)", ha="center", fontsize=8, color=OR)
ax.text(4.65, 4.9, "Shaded: dimensions where our law splits\ninto radius × sphere point (1, 2, 4, 8).\n3D's extras are 4D's shadows (the 600-cell's\nHopf shadow is the icosahedron).", ha="left", fontsize=8.5)
ax.set_xticks(dims); ax.set_ylim(0, 8.4); ax.set_xlabel("dimension"); ax.set_ylabel("number of regular convex shapes"); ax.legend(loc="upper right", fontsize=8.5)
ax.set_title("Figure 8. Regular shapes by dimension (counted by Schläfli's test in checks G1)")
fig.savefig(out + "fig8-regular-shapes-by-dimension.png", bbox_inches="tight"); plt.close(fig)

# ---- fig9: E8 ----
E, X, al, be = e8_from_icosians()
fig, axs = plt.subplots(1, 3, figsize=(15, 4.8))
ax = axs[0]; Y = R8 @ PE8.T
for i, j in edges(R8, d=np.sqrt(2)):
    ax.plot(*Y[[i, j]].T, color="#666", lw=0.05, alpha=0.3)
ax.scatter(*Y.T, s=5, color=OR); ax.set_aspect("equal"); ax.axis("off"); ax.set_title("E8's 240 roots, Petrie picture\n(known; Coxeter number 30)")
ax = axs[1]
r = np.linalg.norm(E[:, :4], axis=1); outer = r > 0.7
for sel, col, lab in ((outer, BL, "600-cell of radius 1"), (~outer, OR, "600-cell of radius 1/φ")):
    V = E[sel, :4]; Yp = V @ PH4.T
    for i, j in edges(V): ax.plot(*Yp[[i, j]].T, color=col, lw=0.2, alpha=0.5)
    ax.scatter(*Yp.T, s=5, color=col, label=lab)
ax.set_aspect("equal"); ax.axis("off"); ax.legend(loc="lower center", fontsize=8, bbox_to_anchor=(0.5, -0.12))
ax.set_title("E8 seen from 4D: two 600-cells, sizes 1 : 1/φ\n(the other 4 coordinates are the golden twin)")
ax = axs[2]
vals = [0.0974, 1 - 0.0974]
ax.pie(vals, labels=["(ab)c = a(bc)\n9.7%", "grouping changes\nthe answer\n90.3%"], colors=[BL, "#dddddd"], startangle=90, textprops={"fontsize": 9})
ax.set_title("E8's corners as octonion units:\nclosed under our law, but grouping matters\n(random triples of the 240 units)")
fig.suptitle("Figure 9. E8: the 8D exceptional shape", y=1.02)
fig.savefig(out + "fig9-e8.png", bbox_inches="tight"); plt.close(fig)

# ---- fig10: 24-cell extras ----
fig, axs = plt.subplots(1, 2, figsize=(13, 4.6))
ax = axs[0]
n = np.arange(1, 21); sig = [sum(d for d in range(1, k + 1) if k % d == 0 and d % 2) for k in n]
ax.bar(n, [24 * s for s in sig], color=BL, label="Hurwitz integers of size n (24-cell corners scaled)")
ax.bar(n, [8 * sum(d for d in range(1, k + 1) if k % d == 0 and d % 4) for k in n], color=OR, width=0.45, label="ways to write n as four squares (Jacobi)")
ax.set_xlabel("n"); ax.set_xticks(n); ax.legend(fontsize=8.5); ax.set_title("The 24 corners make four-square counting exact:\nsize-n Hurwitz integers = 24 × (sum of odd divisors of n)", fontsize=10)
ax = axs[1]; ax.axis("off"); ax.set_aspect("equal")
Pc = PF4
V = T2 @ Pc.T
for i, j in edges(T2): ax.plot(*V[[i, j]].T, color=BL, lw=0.8)
roots = np.array([v for v in itertools.product((-1, 0, 1), repeat=4) if sum(abs(x) for x in v) == 2], float)
Rr = roots @ Pc.T
ax.scatter(*Rr.T, s=40, color=OR, zorder=4, label="the 24 nearest lattice points (kissing number 24)")
ax.scatter([0], [0], s=60, color="k", zorder=5)
ax.legend(loc="lower center", fontsize=8.5, bbox_to_anchor=(0.5, -0.1))
ax.set_title("The 24-cell is the tile of the D4 lattice:\ncopies of it fill 4D space (above 4D only the cube does that).\nLattice points sit one beyond each of its 24 faces", fontsize=10)
fig.suptitle("Figure 10. What the 24-cell brings with it (all known)", y=1.02)
fig.savefig(out + "fig10-24-cell-extras.png", bbox_inches="tight"); plt.close(fig)

# ---- fig1: sorting board ----
rows = [("24-cell", "closed group (24)", "runs down: 48→24→8→2", "p = 4/3 and 4 (self-dual)", "X,Y,Z; H,S in Clifford 48", "Bridge (four squares, D4 tiling)", OR),
        ("600-cell", "closed group (120)", "never shrinks (perfect)", "p = 1.80, cornered side", "X,Y,Z yes; S, T never", "Adds: first corner set the swap\nrule can't wear down; brings φ", BL),
        ("120-cell", "NOT closed; = 5 × 600-cell", "(acted on by 600-cell)", "p = 2.25, boxy side", "—", "Nothing new beyond the 600-cell", RD),
        ("E8 (8D)", "closed loop (240);\ngrouping matters", "swaps plus a grouping rule", "—", "—", "Bridge: 600-cell + golden twin;\noctonions; densest packing", OR),
        ("Only 3 families\nabove 4D", "—", "—", "—", "—", "Confirmed (Schläfli test, dims 5–8)", OR)]
cols = ["shape", "under our law ∘ₖ", "swap rule", "roundness (p dial)", "quantum gates", "verdict"]
fig, ax = plt.subplots(figsize=(15, 4.4)); ax.axis("off")
xs = [0, 1.4, 3.6, 5.8, 7.9, 10.0]; ax.set_xlim(0, 13.3); ax.set_ylim(-0.3, len(rows) + 0.6)
for x, c in zip(xs, cols): ax.text(x, len(rows) + 0.2, c, weight="bold", fontsize=10)
for r, row in enumerate(rows):
    y = len(rows) - r - 0.5
    ax.add_patch(plt.Rectangle((-0.1, y - 0.45), 13.4, 0.9, color=row[-1], alpha=0.10))
    for x, t in zip(xs, row[:-1]): ax.text(x, y, t, va="center", fontsize=9)
ax.text(0, -0.25, "Row colour: orange = known math in our language (bridge), blue = adds something to the theory, red = adds nothing.", fontsize=8.5)
ax.set_title("Figure 1. What each missing shape does in Hypercube Composition Theory")
fig.savefig(out + "fig1-sorting-board.png", bbox_inches="tight"); plt.close(fig)
print("figures done")
