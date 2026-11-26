"""Figures for beyond-the-cube/README.md. Run after checks.py (needs assoc6.json)."""
import json, math, itertools
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d.art3d import Poly3DCollection

OUT = "figures/"
BLUE, ORANGE, GREEN, GRAY, RED, PURPLE = "#2b6cb0", "#dd6b20", "#2f855a", "#a0aec0", "#c53030", "#6b46c1"
plt.rcParams.update({"font.size": 11, "axes.spines.top": False, "axes.spines.right": False})


def cube_edges(ax, lo=-1, hi=1, **kw):
    for a, b in itertools.combinations(itertools.product([lo, hi], repeat=3), 2):
        if sum(x != y for x, y in zip(a, b)) == 1:
            ax.plot(*zip(a, b), **kw)


# ---------------------------------------------------------------- fig0 verdict board
gaps = ["Directions that\ninteract", "Why quantum stops\nat 2√2", "North Star:\nquantum + gravity",
        "G past\none loop", "Frame dragging\n(by hand now)", "Navier–Stokes\nstretching corner", "A new physics\nprediction"]
shapes = ["Associahedron\n(ABHY, u-coords)", "Elliptope\n(+ cut polytope)", "Amplituhedron\n(total positivity)",
          "Permutohedron /\nzonotope", "Simplex /\ncross-polytope", "Birkhoff\npolytope", "Cyclic\npolytope"]
# 2 = carries it further, 1 = partial / testbed, 0 = restatement or nothing
V = np.array([
    [2, 0, 1, 1, 0, 0, 0],  # interact
    [0, 2, 0, 0, 0, 1, 0],  # 2√2
    [1, 0, 1, 0, 0, 0, 0],  # north star
    [1, 0, 0, 0, 0, 0, 0],  # G
    [0, 0, 0, 0, 0, 0, 0],  # frame dragging
    [0, 0, 0, 0, 0, 0, 0],  # NS
    [0, 0, 0, 0, 0, 0, 0],  # prediction
])
txt = [
    ["union rule\nglues it", "", "positivity\nonly", "orders,\nno coupling", "", "", ""],
    ["", "the answer\n(Tsirelson)", "", "", "", "quantum\nsubset", ""],
    ["door via\ndouble copy", "", "planar N=4\nonly", "", "", "", ""],
    ["exact 2D flows\n(Y-system)", "", "", "", "", "", ""],
    ["", "", "", "", "", "", ""],
    ["", "", "", "", "", "", ""],
    ["", "", "", "", "", "", ""],
]
cmap = {0: "#edf2f7", 1: "#fbd38d", 2: "#9ae6b4"}
fig, ax = plt.subplots(figsize=(13, 7.4))
for i in range(len(gaps)):
    for j in range(len(shapes)):
        ax.add_patch(plt.Rectangle((j, -i), 0.96, 0.92, color=cmap[V[i, j]]))
        ax.text(j + 0.48, -i + 0.46, txt[i][j], ha="center", va="center", fontsize=9)
for j, sname in enumerate(shapes):
    ax.text(j + 0.48, 1.25, sname, ha="center", va="center", fontsize=10, weight="bold")
for i, g in enumerate(gaps):
    ax.text(-0.1, -i + 0.46, g, ha="right", va="center", fontsize=10, weight="bold")
ax.set_xlim(-2.2, len(shapes)); ax.set_ylim(-len(gaps) + 0.8, 1.7); ax.axis("off")
for c, lab in [(cmap[2], "carries the theory further on this gap"), (cmap[1], "partial: a door or a testbed, not a result"),
               (cmap[0], "restatement, or nothing to say")]:
    pass
hand = [plt.Rectangle((0, 0), 1, 1, color=cmap[k]) for k in (2, 1, 0)]
ax.legend(hand, ["carries the theory further on this gap", "partial: a door or a testbed, not a result",
                 "restatement, or nothing to say"], loc="lower center", bbox_to_anchor=(0.45, -0.08), ncol=3, frameon=False)
ax.set_title("Seven open gaps (rows) against seven shapes (columns)", fontsize=14, weight="bold", pad=18)
fig.tight_layout(); fig.savefig(OUT + "fig0-verdict-board.png", dpi=140); plt.close(fig)

# ---------------------------------------------------------------- fig1 elliptope
fig = plt.figure(figsize=(12, 5.6))
ax = fig.add_subplot(121, projection="3d")
cube_edges(ax, color=GRAY, lw=1)
tet = [(1, 1, 1), (1, -1, -1), (-1, 1, -1), (-1, -1, 1)]
faces = [[tet[a], tet[b], tet[c]] for a, b, c in itertools.combinations(range(4), 3)]
ax.add_collection3d(Poly3DCollection(faces, facecolor=BLUE, alpha=0.35, edgecolor=BLUE))
th = np.linspace(0, 2*np.pi, 70); ph = np.linspace(0, np.pi, 40)
X, Y = np.meshgrid(np.linspace(-1, 1, 60), np.linspace(-1, 1, 60))
R = np.sqrt(np.clip((1 - X**2)*(1 - Y**2), 0, None))
for sgn in (1, -1):
    ax.plot_wireframe(X, Y, X*Y + sgn*R, color=ORANGE, lw=0.35, alpha=0.6, rstride=3, cstride=3)
for p in tet:
    ax.scatter(*p, color=BLUE, s=40)
ax.text(1.05, 1.05, 1.15, "corner (1,1,1)", fontsize=9)
ax.set_xlabel("x = corr(A,B)"); ax.set_ylabel("y = corr(A,C)"); ax.set_zlabel("z = corr(B,C)")
ax.set_title("Gray: the cube of all ±1 correlations\nBlue: classical (cut polytope, 4 cube corners)\nOrange: quantum (elliptope)", fontsize=10)
ax.view_init(18, 35)
ax2 = fig.add_subplot(122)
# slice z = 0: elliptope x^2 + y^2 <= 1, tetrahedron |x| + |y| <= 1
t = np.linspace(0, 2*np.pi, 300)
ax2.fill(np.cos(t), np.sin(t), color=ORANGE, alpha=0.35, label="quantum (elliptope): disc x² + y² ≤ 1")
ax2.fill([1, 0, -1, 0], [0, 1, 0, -1], color=BLUE, alpha=0.5, label="classical (cut polytope): |x| + |y| ≤ 1")
ax2.plot([-1, 1, 1, -1, -1], [-1, -1, 1, 1, -1], color=GRAY, label="cube face (no rule at all)")
ax2.scatter([1/math.sqrt(2)], [1/math.sqrt(2)], color=RED, zorder=5)
ax2.annotate("quantum reaches here:\nthe 2√2 point", (0.707, 0.707), (0.15, 1.12), arrowprops=dict(arrowstyle="->"), fontsize=9)
ax2.set_aspect("equal"); ax2.set_xlim(-1.25, 1.35); ax2.set_ylim(-1.15, 1.45)
ax2.set_xlabel("x = corr(A,B)"); ax2.set_ylabel("y = corr(A,C)")
ax2.set_title("Slice at z = 0: the disc pokes out past the diamond", fontsize=10)
ax2.legend(loc="upper center", bbox_to_anchor=(0.5, -0.13), fontsize=8.5, frameon=False)
fig.suptitle("Fig 1. The cube's corners make the classical shape; the quantum shape is the elliptope around it", weight="bold")
fig.tight_layout(); fig.savefig(OUT + "fig1-elliptope.png", dpi=140); plt.close(fig)

# ---------------------------------------------------------------- fig2 CHSH bars
fig, ax = plt.subplots(figsize=(8, 4.2))
vals = [2, 2*math.sqrt(2), 4]
labs = ["cube corners\n(classical, any hidden loop)", "elliptope\n(quantum, unit vectors)", "no-signalling box\n(Popescu–Rohrlich)"]
ax.bar(labs, vals, color=[BLUE, ORANGE, GRAY])
for i, v in enumerate(vals):
    ax.text(i, v + 0.08, ["2", "2√2 ≈ 2.83", "4"][i], ha="center", weight="bold")
ax.set_ylabel("best CHSH score"); ax.set_ylim(0, 4.6)
ax.set_title("Fig 2. Three shapes, three ceilings: the 2√2 comes from the elliptope", weight="bold", fontsize=11)
fig.tight_layout(); fig.savefig(OUT + "fig2-chsh-ceilings.png", dpi=140); plt.close(fig)

# ---------------------------------------------------------------- fig3 independence surface inside COR(2)
fig = plt.figure(figsize=(7.5, 6))
ax = fig.add_subplot(111, projection="3d")
T = np.array([[0, 0, 0], [1, 0, 0], [0, 1, 0], [1, 1, 1]], float)
faces = [[T[a], T[b], T[c]] for a, b, c in itertools.combinations(range(4), 3)]
ax.add_collection3d(Poly3DCollection(faces, facecolor=BLUE, alpha=0.12, edgecolor=BLUE))
P, Q = np.meshgrid(np.linspace(0, 1, 25), np.linspace(0, 1, 25))
ax.plot_surface(P, Q, P*Q, color=GREEN, alpha=0.75)
labels = ["never A, never B", "A only", "B only", "always both"]
for p, l in zip(T, labels):
    ax.scatter(*p, color=BLUE, s=30); ax.text(*(p + np.array([0.03, 0.03, 0.05])), l, fontsize=9)
ax.set_xlabel("P(A) = p"); ax.set_ylabel("P(B) = q"); ax.set_zlabel("P(A and B)")
ax.set_title("Fig 3. Blue: every classical pair of events (Pitowsky's tetrahedron).\nGreen: the independent ones, P(A and B) = pq,"
             " where our union\nrule p ∘₋₁ q = p + q − pq is exact.", fontsize=10)
ax.view_init(22, -60)
fig.tight_layout(); fig.savefig(OUT + "fig3-independence-surface.png", dpi=140); plt.close(fig)

# ---------------------------------------------------------------- fig4 pentagon: ABHY + u-equation as union
D = json.load(open("assoc6.json"))
pent = np.array(D["pent"])
fig, axs = plt.subplots(1, 2, figsize=(12.5, 5.6))
ax = axs[0]
ax.fill(*np.vstack([pent, pent[:1]]).T, color=GREEN, alpha=0.3)
ax.plot(*np.vstack([pent, pent[:1]]).T, color=GREEN, lw=2)
names = ["X₁₃ = 0", "X₁₄ = 0", "X₂₄ = 0", "X₂₅ = 0", "X₃₅ = 0"]
# vertex i is between facet order[i] and order[i+1]; edge between vertex i-1 and i is facet order[i]
for i in range(5):
    a, b = pent[i-1], pent[i]
    m = (a + b)/2; c = pent.mean(0); d = (m - c)/np.linalg.norm(m - c)
    ax.text(*(m + 0.22*d), names[i], ha="center", va="center", fontsize=10, color=GREEN, weight="bold")
ax.set_aspect("equal"); ax.set_xlabel("X₁₃ (a planar channel)"); ax.set_ylabel("X₁₄")
ax.set_title("ABHY kinematic associahedron, 5 particles:\nits canonical form is the φ³ amplitude\n1/(X₁₃X₁₄) + 1/(X₁₄X₂₄) + 1/(X₂₄X₂₅) + 1/(X₂₅X₃₅) + 1/(X₃₅X₁₃)", fontsize=10)
ax.margins(0.25)
ax = axs[1]
ang = np.pi/2 - 2*np.pi*np.arange(5)/5
pts = np.c_[np.cos(ang), np.sin(ang)]
ax.plot(*np.vstack([pts, pts[:1]]).T, color="k", lw=1)
for i, p in enumerate(pts):
    ax.text(*(1.13*p), str(i+1), ha="center", va="center", fontsize=12, weight="bold")
def chord(i, j, **kw):
    ax.plot(*zip(pts[i-1], pts[j-1]), **kw)
chord(1, 3, color=BLUE, lw=3)
chord(2, 4, color=ORANGE, lw=2.5, ls="--"); chord(2, 5, color=ORANGE, lw=2.5, ls="--")
ax.text(0, -1.45, "u₁₃ = p₂₄ ∘₋₁ p₂₅    (p = 1 − u)\n"
        "the blue chord's coordinate is the union over the\norange chords that cross it", ha="center", fontsize=11)
ax.text(-0.55, 0.62, "chord 13", color=BLUE, fontsize=10, weight="bold")
ax.text(0.1, -0.68, "chords 24, 25 (cross 13)", color=ORANGE, fontsize=10, weight="bold")
ax.set_aspect("equal"); ax.axis("off"); ax.set_ylim(-1.8, 1.3)
ax.set_title("Its coordinates: u_X + ∏ u_Y (over chords Y crossing X) = 1", fontsize=10)
fig.suptitle("Fig 4. The associahedron has our union rule built into its coordinates", weight="bold")
fig.tight_layout(); fig.savefig(OUT + "fig4-pentagon-union.png", dpi=140); plt.close(fig)

# ---------------------------------------------------------------- fig5 cube to 3D associahedron
fig = plt.figure(figsize=(12, 5.4))
ax = fig.add_subplot(121, projection="3d")
cube_edges(ax, 0, 1, color=BLUE, lw=2)
for p in itertools.product([0, 1], repeat=3):
    ax.scatter(*p, color=BLUE, s=25)
ax.set_title("No coupling: the 3-cube\n8 corners, 12 edges (three independent intervals)", fontsize=10)
ax.set_axis_off()
ax = fig.add_subplot(122, projection="3d")
V6 = np.array(D["V6"]); V6 = (V6 - V6.min(0))/(V6.max(0) - V6.min(0))  # affine rescale, same shape type
for a, b in D["edges"]:
    ax.plot(*zip(V6[a], V6[b]), color=GREEN, lw=2)
ax.scatter(*V6.T, color=GREEN, s=25)
ax.set_title("Coupled (six particles): the 3D associahedron\n14 corners, 21 edges, 9 faces (one per chord of a hexagon)", fontsize=10)
ax.set_axis_off(); ax.view_init(22, -55)
fig.suptitle("Fig 5. Same three directions; coupling them bends the cube into the associahedron (checked: B7)", weight="bold")
fig.tight_layout(); fig.savefig(OUT + "fig5-cube-to-associahedron.png", dpi=140); plt.close(fig)

# ---------------------------------------------------------------- fig6 gluing on a line (AdS3 link)
fig, ax = plt.subplots(figsize=(11, 3.6))
z = [0.0, 1.0, 2.6, 3.7, 5.2]
ax.axhline(0, color="k", lw=1)
spans = [(z[0], z[1], "A", BLUE), (z[2], z[3], "B₁", ORANGE), (z[3], z[4], "B₂", PURPLE)]
for a, b, l, c in spans:
    ax.plot([a, b], [0, 0], color=c, lw=8, solid_capstyle="butt"); ax.text((a+b)/2, -0.22, l, ha="center", color=c, fontsize=13, weight="bold")
for a, b, c, h in [(0.5, 3.15, ORANGE, 1.0), (0.5, 4.45, PURPLE, 1.45)]:
    t = np.linspace(0, np.pi, 80); r = (b - a)/2
    ax.plot((a+b)/2 + r*np.cos(t), h*np.sin(t), color=c, lw=1.5)
eta = lambda a1, a2, b1, b2: (a2 - a1)*(b2 - b1)/((b1 - a1)*(b2 - a2))
x1 = eta(z[0], z[1], z[2], z[3]); x2 = eta(z[0], z[1], z[3], z[4]); x12 = eta(z[0], z[1], z[2], z[4])
ax.text(5.6, 0.9, f"x(A,B₁) = {x1:.3f}\nx(A,B₂) = {x2:.3f}\nx(A,B₁B₂) = {x12:.3f}\n= x(A,B₁) ∘₋₁ x(A,B₂)\n   = {x1 + x2 - x1*x2:.3f}",
        fontsize=11, va="center")
ax.set_xlim(-0.3, 8.6); ax.set_ylim(-0.3, 1.7); ax.axis("off")
ax.set_title("Fig 6. One number, two fields: x is the AdS₃ linking probability on our probability page and 1 − u in the associahedron;\n"
             "gluing neighbouring intervals is our union rule in both", fontsize=10.5, weight="bold")
fig.tight_layout(); fig.savefig(OUT + "fig6-gluing-one-number.png", dpi=140); plt.close(fig)

# ---------------------------------------------------------------- fig7 budget region
fig, ax = plt.subplots(figsize=(6.6, 7.0))
xs = np.linspace(0, 1, 400)
ks = [(-0.999, "k = −1: the whole square (independent)", BLUE), (-0.6, "k = −0.6", "#63b3ed"),
      (1e-9, "k = 0: the triangle (simplex, exclusive)", "k"), (1.0, "k = 1", "#f6ad55"), (5.0, "k = 5", ORANGE), (60.0, "k = 60: shrinking to the axes", RED)]
for k, lab, c in ks:
    yv = (1 - xs)/(1 + k*xs)
    ax.plot(xs, np.minimum(yv, 1), color=c, lw=2.2, label=lab)
ax.plot([0, 1, 1, 0, 0], [0, 0, 1, 1, 0], color=GRAY, lw=1)
ax.set_aspect("equal"); ax.set_xlabel("x"); ax.set_ylabel("y")
ax.legend(fontsize=8.5, loc="upper center", bbox_to_anchor=(0.5, -0.12), ncol=2, frameon=False)
ax.set_title("Fig 7. What the square becomes as the dial turns:\nthe region x ∘ₖ y ≤ 1 (edge = Sugeno's λ-measure rule)", fontsize=10.5, weight="bold")
fig.tight_layout(); fig.savefig(OUT + "fig7-budget-region.png", dpi=140); plt.close(fig)

# ---------------------------------------------------------------- fig8 Y-system pentagon and Lee-Yang
fig, axs = plt.subplots(1, 2, figsize=(11.5, 4.6))
y = [0.37, 2.9]
for _ in range(10):
    y.append((1 + y[-1]) / y[-2])
ax = axs[0]
ax.plot(range(12), y[:12], "o-", color=GREEN)
for n in (0, 5, 10):
    ax.axvline(n, color=GRAY, ls=":")
ax.set_xlabel("step n"); ax.set_ylabel("y(n)")
ax.set_title("y(n+1)·y(n−1) = 1 + y(n) comes back every 5 steps\n(the pentagon's 5 corners; Zamolodchikov 1991)", fontsize=10)
ax = axs[1]
phi = (1 + 5**0.5)/2
ax.bar(["fixed point y", "Rogers L(1/(1+y))\n÷ (π²/6)", "c_eff of the\nLee–Yang flow"], [phi, 6/math.pi**2*(math.pi**2/15), 0.4],
       color=[GREEN, GREEN, ORANGE])
ax.text(0, phi + 0.05, "golden ratio\n1.618", ha="center", fontsize=9)
ax.text(1, 0.45, "0.4", ha="center", fontsize=9); ax.text(2, 0.45, "2/5 exactly", ha="center", fontsize=9)
ax.set_ylim(0, 2.1)
ax.set_title("At the fixed point the dilogarithm gives the\nexact central charge, past any loop order", fontsize=10)
fig.suptitle("Fig 8. The pentagon as an exactly solved flow: a testbed for 'past one loop'", weight="bold")
fig.tight_layout(); fig.savefig(OUT + "fig8-ysystem-lee-yang.png", dpi=140); plt.close(fig)
print("figures written")
