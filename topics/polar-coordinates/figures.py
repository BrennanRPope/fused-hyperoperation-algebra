"""Figures for the polar-coordinates DRAFT. Run: python3 figures.py -> figures/*.png"""
import itertools
import os

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import Circle, Ellipse, Polygon, Rectangle
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
from scipy.spatial import ConvexHull

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "figures")
os.makedirs(OUT, exist_ok=True)
plt.rcParams.update({"font.size": 11, "axes.titlesize": 12, "figure.dpi": 130})
BLUE, ORANGE, GREEN, RED, PURPLE, GRAY = "#3366cc", "#e67e22", "#2e9e5b", "#c0392b", "#8e44ad", "#777777"


def save(fig, name):
    fig.savefig(os.path.join(OUT, name), bbox_inches="tight")
    plt.close(fig)
    print("wrote", name)


def hull_poly(pts):
    pts = np.asarray(pts, float)
    h = ConvexHull(pts)
    return pts[h.vertices]


def to_plane(pts):
    """Points on the plane t0+t1+t2 = const -> 2D drawing coordinates."""
    pts = np.asarray(pts, float)
    e1 = np.array([1, -1, 0]) / np.sqrt(2)
    e2 = np.array([1, 1, -2]) / np.sqrt(6)
    return np.c_[pts @ e1, pts @ e2]


def to_space(pts):
    """Points on t0+..+t3 = const in R^4 -> 3D drawing coordinates."""
    pts = np.asarray(pts, float)
    B = np.array([[1, -1, 0, 0], [1, 1, -2, 0], [1, 1, 1, -3]], float)
    B /= np.linalg.norm(B, axis=1, keepdims=True)
    return pts @ B.T


def minkowski(sets):
    pts = np.zeros((1, sets[0].shape[1]))
    for s in sets:
        pts = np.unique(np.round((pts[:, None] + s[None]).reshape(-1, s.shape[1]), 9), axis=0)
    return pts


def simplex_on(I, n, w=1.0):
    return w * np.array([np.eye(n)[i] for i in I])


def shape(y, n):
    """Sum over corners I of y_I * simplex on I (Postnikov's generalized permutohedron)."""
    sets = []
    for mask, w in y.items():
        if w > 0:
            sets.append(simplex_on([i for i in range(n) if mask >> i & 1], n, w))
    return minkowski(sets)


# ------------------------------------------------------------- fig 1
def fig1():
    fig, ax = plt.subplots(1, 3, figsize=(15, 4.8))
    a0, b0 = 0.9, 0.5
    A = ax[0]
    A.add_patch(Rectangle((0, 0), 1, 1, fc="#dddddd", ec="k"))
    A.add_patch(Rectangle((1, 0), a0, 1, fc="#cfe0fa", ec="k"))
    A.add_patch(Rectangle((0, 1), 1, b0, fc="#fbe0c4", ec="k"))
    A.add_patch(Rectangle((1, 1), a0, b0, fc="#d5f0dd", ec="k"))
    for (x, y, t) in [(0.5, 0.5, "corner ∅\npiece 1"), (1 + a0 / 2, 0.5, "corner {a}\npiece a"),
                      (0.5, 1 + b0 / 2, "corner {b}\npiece b"), (1 + a0 / 2, 1 + b0 / 2, "corner {a,b}\npiece ab")]:
        A.text(x, y, t, ha="center", va="center", fontsize=10)
    A.set_xlim(-0.1, 2.1), A.set_ylim(-0.1, 1.8), A.set_aspect("equal")
    A.set_title("Our tesseract coordinates (Cartesian):\nbox with sides 1+a, 1+b cut at 1")
    A.set_xlabel("1 + a"), A.set_ylabel("1 + b")

    c = 0.35
    seq = [(1 + c) ** i - 1 for i in range(7)]
    B = ax[1]
    for i, x in enumerate(seq):
        for j, y in enumerate(seq):
            B.plot(x, y, "o", color=BLUE, ms=4)
    for i in range(6):
        B.annotate("", (seq[i + 1], 0), (seq[i], 0), arrowprops=dict(arrowstyle="->", color=ORANGE))
    B.set_title(f"Composing a fixed step c = {c} again and again\n(a ∘ c ∘ c ∘ …) in Cartesian a, b: steps grow")
    B.set_xlabel("a"), B.set_ylabel("b")
    C = ax[2]
    u = np.log1p(np.array(seq))
    for x in u:
        for y in u:
            C.plot(x, y, "o", color=BLUE, ms=4)
    for i in range(6):
        C.annotate("", (u[i + 1], 0), (u[i], 0), arrowprops=dict(arrowstyle="->", color=ORANGE))
    C.set_title("Same points in log-radius coordinates\nu = log(1+a), v = log(1+b): ∘ becomes +, an even grid")
    C.set_xlabel("u = log(1 + a)   (radial coordinate of axis a)")
    C.set_ylabel("v = log(1 + b)")
    C.set_aspect("equal")
    save(fig, "fig1-cartesian-vs-log-radius.png")


# ------------------------------------------------------------- fig 2
def fig2():
    fig, ax = plt.subplots(1, 3, figsize=(15, 5))
    cases = [(0, 1, "k² − 4d < 0  (here k=0, d=1)\nruler is an ANGLE: comes back around", BLUE),
             (0, 0, "k² − 4d = 0  (here k=0, d=0)\nparabolic: plain adding of slopes", GREEN),
             (0, -1, "k² − 4d > 0  (here k=0, d=−1)\nruler is a RADIUS / rapidity: never returns", RED)]
    t = np.linspace(-3, 3, 600)
    for A, (k, d, title, col) in zip(ax, cases):
        # unit norm curve x^2 + kxy + d y^2 = 1
        X, Y = np.meshgrid(np.linspace(-2.2, 2.2, 400), np.linspace(-2.2, 2.2, 400))
        A.contour(X, Y, X ** 2 + k * X * Y + d * Y ** 2, levels=[1], colors=[col], linewidths=2)
        A.axhline(0, color="#bbb", lw=0.8), A.axvline(0, color="#bbb", lw=0.8)
        # orbit of repeated composition with a fixed step a
        a = 0.35
        z = np.array([1.0, 0.0])
        pts = [z]
        for _ in range(14):
            x0, x1 = z
            z = np.array([x0 - d * x1 * a, x0 * a + x1 + k * x1 * a])
            nz = np.sqrt(abs(z[0] ** 2 + k * z[0] * z[1] + d * z[1] ** 2))
            z = z / nz
            pts.append(z)
        pts = np.array(pts)
        pts = pts[np.all(np.abs(pts) < 2.2, axis=1)]
        A.plot(pts[:, 0], pts[:, 1], "o-", color=ORANGE, ms=5, lw=1)
        A.annotate("start 1", pts[0], (pts[0][0] + 0.15, pts[0][1] - 0.35), fontsize=10)
        if d == -1:
            A.plot(t, t, "--", color=GRAY, lw=1), A.plot(t, -t, "--", color=GRAY, lw=1)
            A.text(-2.05, 1.45, "dashed: light cone,\nthe two fixed points\n(speed of light)", fontsize=9, color=GRAY)
        A.set_xlim(-2.2, 2.2), A.set_ylim(-2.2, 2.2), A.set_aspect("equal")
        A.set_title(title)
        A.set_xlabel("x  (the 1 part)"), A.set_ylabel("y  (the ε part)")
    fig.suptitle("The plane dial d decides whether an axis is an angle or a radius.\n"
                 "Each point x + yε with ε² = kε − d; composing with a (orange dots) multiplies by 1 + εa; "
                 "the slope y/x is the number, the curve is radius 1", y=1.04)
    save(fig, "fig2-plane-dial-angle-or-radius.png")


# ------------------------------------------------------------- fig 3
def fig3():
    fig, ax = plt.subplots(1, 3, figsize=(19, 5.5))
    fig.subplots_adjust(wspace=0.35)
    A = ax[0]
    th = np.linspace(0, 2 * np.pi, 200)
    for r in (0.5, 1, 1.5, 2):
        A.plot(r * np.cos(th), r * np.sin(th), color="#bbb", lw=0.8)
    for a in np.linspace(0, 2 * np.pi, 13)[:-1]:
        A.plot([0, 2.1 * np.cos(a)], [0, 2.1 * np.sin(a)], color="#bbb", lw=0.8)
    z = 1.5 * np.exp(1j * 0.8)
    A.annotate("", (z.real, z.imag), (0, 0), arrowprops=dict(arrowstyle="->", color=BLUE, lw=2))
    A.text(-1.6, 1.25, "radius r = |1 + ka|", color=BLUE)
    A.add_patch(matplotlib.patches.Arc((0, 0), 0.8, 0.8, theta1=0, theta2=np.degrees(0.8), color=ORANGE, lw=2))
    A.text(0.5, 0.12, "angle θ", color=ORANGE)
    A.set_aspect("equal"), A.set_xlim(-2.3, 2.3), A.set_ylim(-2.3, 2.3)
    A.set_title("One complex axis: the lift 1 + ka in polar form\nlog r = radial ruler, θ = the one period of ∘ₖ")
    A.axis("off")

    B = ax[1]
    for u in np.linspace(-1.5, 1.5, 4):
        for v in np.linspace(-1.5, 1.5, 4):
            B.add_patch(Circle((u, v), 0.22, fill=False, color=BLUE, lw=1.2))
            B.add_patch(Circle((u + 0.09, v), 0.12, fill=False, color=ORANGE, lw=1.2))
    B.set_xlim(-2.1, 2.1), B.set_ylim(-2.1, 2.1), B.set_aspect("equal")
    B.set_xlabel("u₁ = log r₁"), B.set_ylabel("u₂ = log r₂")
    B.set_title("Two axes: each point of the log-radius plane\ncarries a torus of angles (θ₁, θ₂)\n[torus drawn as two circles]")

    C = ax[2]
    xs = np.linspace(0, 1, 7)
    C.plot([0, 1], [0, 0], color="k", lw=2)
    for x in xs:
        w1 = 0.35 * np.sqrt(x)  # circle of θ1 shrinks at x=0
        w2 = 0.35 * np.sqrt(1 - x)
        if w1 > 1e-6:
            C.add_patch(Ellipse((x, 0.55), 0.09, 2 * w1, fill=False, color=BLUE, lw=1.3))
        else:
            C.plot(x, 0.55, "o", color=BLUE)
        if w2 > 1e-6:
            C.add_patch(Ellipse((x, -0.55), 0.09, 2 * w2, fill=False, color=ORANGE, lw=1.3))
        else:
            C.plot(x, -0.55, "o", color=ORANGE)
    C.text(0.5, 1.05, "θ₁ circle (size |z₁|)", ha="center", color=BLUE)
    C.text(0.5, -1.1, "θ₂ circle (size |z₂|)", ha="center", color=ORANGE)
    C.text(-0.02, 0.1, "|z₁|² = 0", ha="right"), C.text(1.02, 0.1, "|z₁|² = 1")
    C.set_xlim(-0.35, 1.35), C.set_ylim(-1.3, 1.3)
    C.axis("off")
    C.set_title("Hyperspherical (Hopf) coordinates on the\n3-sphere |z₁|²+|z₂|² = 1: squared radii fill a\nSIMPLEX (here a segment), angles ride on top")
    save(fig, "fig3-polar-in-many-dimensions.png")


# ------------------------------------------------------------- fig 4
def moment(u, verts):
    w = np.exp(2 * (verts @ u.T).T)
    return (w @ verts) / w.sum(1, keepdims=True)


def fig4():
    sq = np.array([[0, 0], [1, 0], [1, 1], [0, 1]], float)
    tri = np.array([[0, 0], [1, 0], [0, 1]], float)
    hexa = to_plane(np.array(list(itertools.permutations([0, 1, 2])), float))
    pent3 = shape({0b011: 1, 0b110: 1, 0b111: 1}, 3)
    pent = to_plane(pent3)
    shapes = [("Square = (P¹)²\n4 charts = our 4 tesseract corners", sq, ["∅", "{a}", "{a,b}", "{b}"]),
              ("Triangle = CP²\n(hyperspherical / Hopf)", tri, None),
              ("Hexagon = permutohedron\n(braid-arrangement surface)", hexa, None),
              ("Pentagon = associahedron\n(Loday)", pent, None)]
    fig, ax = plt.subplots(1, 4, figsize=(18, 5))
    g = np.linspace(-2.5, 2.5, 11)
    s = np.linspace(-6, 6, 300)
    for A, (title, V, labels) in zip(ax, shapes):
        # work in the polygon's own 2D coordinates; log-radius grid u = (s, const) etc.
        Vc = V - V.mean(0)
        P = hull_poly(V)
        A.add_patch(Polygon(P, closed=True, fc="#f3f3f3", ec="k", lw=2))
        for c in g:
            for line, col in ((np.c_[s, np.full_like(s, c)], BLUE), (np.c_[np.full_like(s, c), s], ORANGE)):
                m = moment(line, Vc) + V.mean(0)
                A.plot(m[:, 0], m[:, 1], color=col, lw=0.8)
        for i, p in enumerate(P):
            A.plot(*p, "ko", ms=6)
        if labels:
            off = [(-0.13, -0.1), (0.03, -0.1), (0.03, 0.04), (-0.13, 0.04)]
            for p, lab, o in zip(sq, labels, off):
                A.text(p[0] + o[0], p[1] + o[1], lab, fontsize=10)
        A.set_aspect("equal"), A.axis("off")
        A.set_title(title)
    fig.suptitle("Same log-radius grid (blue: u₁ lines, orange: u₂ lines), same angles, four different shapes.\n"
                 "The moment map squeezes the whole log-radius plane into the shape; each black corner is one chart "
                 "(one fixed point)", y=1.03)
    save(fig, "fig4-same-grid-four-shapes.png")


# ------------------------------------------------------------- fig 5
def draw_poly3d(A, pts3, color, title, sub):
    h = ConvexHull(pts3)
    # merge coplanar triangles by drawing all faces with light edges, then true edges
    faces = [pts3[s] for s in h.simplices]
    A.add_collection3d(Poly3DCollection(faces, facecolor=color, alpha=0.35, edgecolor="none"))
    eq = h.equations
    V = h.vertices
    for i, j in itertools.combinations(V, 2):
        mid = (pts3[i] + pts3[j]) / 2
        on = np.sum(np.abs(eq[:, :-1] @ mid + eq[:, -1]) < 1e-7)
        normals = eq[np.abs(eq[:, :-1] @ mid + eq[:, -1]) < 1e-7, :-1]
        if on >= 2 and np.linalg.matrix_rank(np.round(normals, 6)) >= 2:
            A.plot(*zip(pts3[i], pts3[j]), color="k", lw=1.2)
    A.scatter(*pts3[V].T, color="k", s=10)
    lo, hi = pts3.min(0), pts3.max(0)
    c, r = (lo + hi) / 2, (hi - lo).max() / 2
    A.set_xlim(c[0] - r, c[0] + r), A.set_ylim(c[1] - r, c[1] + r), A.set_zlim(c[2] - r, c[2] + r)
    A.set_box_aspect((1, 1, 1)), A.set_axis_off()
    A.set_title(f"{title}\n{sub}", fontsize=10.5)


def fig5():
    n = 4
    cube = shape({0b0011: 1, 0b0101: 1, 0b1001: 1}, n)
    simp = shape({0b1111: 1}, n)
    perm = shape({m: 1 for m in (0b0011, 0b0101, 0b1001, 0b0110, 0b1010, 0b1100)}, n)
    assoc = shape({m: 1 for m in (0b0011, 0b0110, 0b1100, 0b0111, 0b1110, 0b1111)}, n)
    items = [(cube[:, 1:], "#3366cc", "CUBE: 8 vertices", "Newt((1+x)(1+y)(1+z)) = Newt(1 + x∘y∘z)\nsum of 3 segments [0, eᵢ]"),
             (to_space(simp), "#2e9e5b", "SIMPLEX: 4 vertices", "Newt(1 + x + y + z)\none simplex on the full corner"),
             (to_space(perm), "#e67e22", "PERMUTOHEDRON: 24 vertices", "Newt(∏ (xᵢ + xⱼ))\nsum of all 6 segments [eᵢ, eⱼ]"),
             (to_space(assoc), "#8e44ad", "ASSOCIAHEDRON: 14 vertices", "Newt(∏ (xᵢ + … + xⱼ))\nsum of interval simplices (Loday)")]
    fig = plt.figure(figsize=(18, 5.2))
    for i, (p, col, t, s) in enumerate(items):
        A = fig.add_subplot(1, 4, i + 1, projection="3d")
        A.view_init(elev=22, azim=35 + 10 * i)
        draw_poly3d(A, p, col, t, s)
    fig.suptitle("In log-radius coordinates every shape is the Newton polytope of a polynomial.\n"
                 "Multiplying polynomials adds shapes (Minkowski sum); the cube is the shape of our own box coordinate",
                 y=1.04)
    save(fig, "fig5-four-shapes-as-polynomials.png")


# ------------------------------------------------------------- fig 6
CUBE_POS = {0b000: (0, 0), 0b001: (1, 0), 0b010: (0.45, 0.35), 0b100: (0, 1),
            0b011: (1.45, 0.35), 0b101: (1, 1), 0b110: (0.45, 1.35), 0b111: (1.45, 1.35)}
CUBE_LAB = {0b000: "∅", 0b001: "0", 0b010: "1", 0b100: "2", 0b011: "01", 0b101: "02", 0b110: "12", 0b111: "012"}


def draw_weight_cube(A, y, col):
    for m1, m2 in itertools.combinations(range(8), 2):
        if bin(m1 ^ m2).count("1") == 1:
            A.plot(*zip(CUBE_POS[m1], CUBE_POS[m2]), color="#bbb", lw=1)
    for m, p in CUBE_POS.items():
        w = y.get(m, 0)
        if w > 0:
            A.add_patch(Circle(p, 0.09 + 0.07 * w, color=col, alpha=0.85))
        else:
            A.add_patch(Circle(p, 0.04, color="#999"))
        A.text(p[0] + 0.12, p[1] - 0.12, CUBE_LAB[m], fontsize=9)
    A.set_xlim(-0.3, 1.8), A.set_ylim(-0.3, 1.65), A.set_aspect("equal"), A.axis("off")


def fig6():
    shapes = [("Simplex (triangle)", {0b111: 1}, GREEN),
              ("Cube (square)", {0b011: 1, 0b101: 1}, BLUE),
              ("Associahedron (pentagon)", {0b011: 1, 0b110: 1, 0b111: 1}, PURPLE),
              ("Permutohedron (hexagon)", {0b011: 1, 0b101: 1, 0b110: 1}, ORANGE)]
    fig, ax = plt.subplots(2, 4, figsize=(17, 8.5), gridspec_kw={"height_ratios": [1, 1.1]})
    R = max(np.abs(hull_poly(to_plane(shape(y, 3))) - hull_poly(to_plane(shape(y, 3))).mean(0)).max()
            for _, y, _ in shapes) * 1.1
    for j, (t, y, col) in enumerate(shapes):
        draw_weight_cube(ax[0, j], y, col)
        ax[0, j].set_title(f"{t}\nweights on the corners of the 3-cube")
        P = hull_poly(to_plane(shape(y, 3)))
        P = P - P.mean(0)
        A = ax[1, j]
        A.add_patch(Polygon(P, closed=True, fc=col, alpha=0.3, ec="k", lw=2))
        A.plot(*np.vstack([P, P[:1]]).T, "k-")
        A.set_xlim(-R, R), A.set_ylim(-R, R), A.set_aspect("equal"), A.axis("off")
        A.set_title(f"{len(P)} vertices")
    fig.suptitle("A shape is a point in tesseract coordinates: one weight per corner of the cube.\n"
                 "Shape = Σ (weight on corner I) × (simplex on I).  Facet heights = our zeta transform of the weights "
                 "(Postnikov's z_I).  Axis 0 is the constant 1.", y=1.02)
    save(fig, "fig6-shape-is-a-point-in-tesseract-coordinates.png")


# ------------------------------------------------------------- fig 7
def fig7():
    simp = {0b111: 1.0}
    assoc = {0b011: 1.0, 0b110: 1.0, 0b111: 1.0}
    perm = {0b011: 1.0, 0b101: 1.0, 0b110: 1.0}
    cube = {0b011: 1.0, 0b101: 1.0}

    def mix(y1, y2, s):
        keys = set(y1) | set(y2)
        return {k: (1 - s) * y1.get(k, 0) + s * y2.get(k, 0) for k in keys}

    path = [(simp, assoc, 0), (simp, assoc, 0.25), (simp, assoc, 0.5), (simp, assoc, 1),
            (assoc, perm, 0.5), (assoc, perm, 1), (perm, cube, 0.5), (perm, cube, 1)]
    names = ["simplex", "", "", "associahedron", "", "permutohedron", "", "cube"]
    fig, ax = plt.subplots(2, 8, figsize=(22, 6), gridspec_kw={"height_ratios": [1, 0.75]})
    cmap = plt.get_cmap("viridis")
    for j, (y1, y2, s) in enumerate(path):
        y = mix(y1, y2, s)
        y = {k: v for k, v in y.items() if v > 1e-9}
        P = hull_poly(to_plane(shape(y, 3)))
        P = P - P.mean(0)
        A = ax[0, j]
        col = cmap(j / 7)
        A.add_patch(Polygon(P, closed=True, fc=col, alpha=0.45, ec="k", lw=1.5))
        A.set_xlim(-1.5, 1.5), A.set_ylim(-1.5, 1.5), A.set_aspect("equal"), A.axis("off")
        A.set_title(f"{names[j]}\n{len(P)} vertices" if names[j] else f"\n{len(P)} vertices")
        B = ax[1, j]
        keys = [0b001, 0b010, 0b100, 0b011, 0b101, 0b110, 0b111]
        B.bar(range(7), [y.get(k, 0) for k in keys], color=col)
        B.set_xticks(range(7)), B.set_xticklabels([CUBE_LAB[k] for k in keys], fontsize=8)
        B.set_ylim(0, 1.1)
        if j == 0:
            B.set_ylabel("corner weight")
        else:
            B.set_yticks([])
    fig.suptitle("Switching shapes = sliding corner weights (bottom bars). Every in-between shape is real; "
                 "vertex count jumps only when a weight starts or stops being zero", y=1.02)
    save(fig, "fig7-sliding-weights-switches-shapes.png")


# ------------------------------------------------------------- fig 8
def fig8():
    fig = plt.figure(figsize=(18, 5.5))
    A = fig.add_subplot(1, 4, 1)
    th = np.linspace(0, 2 * np.pi, 800)
    ps = [1, 1.5, 2, 4, 12]
    cols = plt.get_cmap("plasma")(np.linspace(0, 0.85, len(ps) + 1))
    for p, c in zip(ps, cols):
        x, y = np.cos(th), np.sin(th)
        r = (np.abs(x) ** p + np.abs(y) ** p) ** (-1 / p)
        A.plot(r * x, r * y, color=c, lw=2, label=f"p = {p}")
    A.plot([-1, 1, 1, -1, -1], [-1, -1, 1, 1, -1], color=cols[-1], lw=2, label="p = ∞")
    A.set_aspect("equal"), A.legend(fontsize=9, loc="center")
    A.set_title("One dial p: |u₁|ᵖ + |u₂|ᵖ = 1\ncross-polytope → circle → square")
    A.set_xlabel("u₁ (log-radius)"), A.set_ylabel("u₂")
    oct_v = np.vstack([np.eye(3), -np.eye(3)])
    cube_v = np.array(list(itertools.product([-1, 1], repeat=3)), float)
    for k, (p, name) in enumerate([(1, "p = 1: CROSS-POLYTOPE\n(octahedron)"), (2, "p = 2: SPHERE"),
                                   (None, "p = ∞: CUBE")]):
        B = fig.add_subplot(1, 4, k + 2, projection="3d")
        col = cols[[0, 2, 5][k]]
        if p == 2:
            uu, vv = np.meshgrid(np.linspace(0, 2 * np.pi, 60), np.linspace(0, np.pi, 30))
            B.plot_surface(np.cos(uu) * np.sin(vv), np.sin(uu) * np.sin(vv), np.cos(vv), color=col,
                           alpha=0.75, linewidth=0.2, edgecolor="k")
        else:
            V = oct_v if p == 1 else cube_v
            h = ConvexHull(V)
            B.add_collection3d(Poly3DCollection([V[t] for t in h.simplices], facecolor=col, alpha=0.6,
                                                edgecolor="none"))
            for i, j in itertools.combinations(range(len(V)), 2):
                d = np.abs(V[i] - V[j]).sum()
                if (p == 1 and abs(d - 2) < 1e-9 and not np.allclose(V[i], -V[j])) or (p is None and abs(d - 2) < 1e-9 and np.count_nonzero(V[i] - V[j]) == 1):
                    B.plot(*zip(V[i], V[j]), color="k", lw=1.2)
        B.set_xlim(-1, 1), B.set_ylim(-1, 1), B.set_zlim(-1, 1)
        B.set_box_aspect((1, 1, 1)), B.set_axis_off()
        B.set_title(name)
    fig.suptitle("The p dial is a temperature: (|u₁|ᵖ+…+|uₙ|ᵖ)^(1/p) = exp( soft-max of log|uᵢ| at temperature 1/p ).\n"
                 "p → ∞ is the same cooling limit (sum → max) as Maslov dequantization", y=1.04)
    save(fig, "fig8-p-dial-cross-polytope-sphere-cube.png")


# ------------------------------------------------------------- fig 9
def amoeba_points(polys, N=4000):
    """Points (log|x|, log|y|) on the zero set of a product of linear-ish factors."""
    out = []
    r = np.exp(np.random.uniform(-6, 6, N))
    al = np.random.uniform(0, 2 * np.pi, N)
    x = r * np.exp(1j * al)
    for kind in polys:
        if kind == "1+x":       # x = -1, y free
            y = np.exp(np.random.uniform(-6, 6, N)) * np.exp(1j * al)
            out.append(np.c_[np.zeros(N), np.log(np.abs(y))])
        elif kind == "1+y":
            out.append(np.c_[np.log(np.abs(x)), np.zeros(N)])
        elif kind == "x+y":     # y = -x
            out.append(np.c_[np.log(np.abs(x)), np.log(np.abs(x))])
        elif kind == "1+x+y":
            y = -1 - x
            out.append(np.c_[np.log(np.abs(x)), np.log(np.abs(y))])
            yy = r * np.exp(1j * al)  # also sample by y, so the downward ray shows
            xx = -1 - yy
            out.append(np.c_[np.log(np.abs(xx)), np.log(np.abs(yy))])
    return np.vstack(out)


def fig9():
    np.random.seed(3)
    rows = [("Cube (square)", ["1+x", "1+y"], "(1+x)(1+y) = 1 + x∘y", np.array([[0, 0], [1, 0], [1, 1], [0, 1]]),
             "Cartesian axes\nuᵢ = 0"),
            ("Simplex (triangle)", ["1+x+y"], "1 + x + y", np.array([[0, 0], [1, 0], [0, 1]]),
             "tropical line\n(3 rays)"),
            ("Permutohedron (hexagon)", ["1+x", "1+y", "x+y"], "(1+x)(1+y)(x+y)",
             np.array([[0, 1], [1, 0], [2, 0], [2, 1], [1, 2], [0, 2]]), "braid arrangement\nuᵢ = 0, u₁ = u₂")]
    fig, ax = plt.subplots(2, 3, figsize=(15, 10.5))
    fig.subplots_adjust(hspace=0.45)
    for j, (t, polys, f, newt, grid) in enumerate(rows):
        A = ax[0, j]
        P = hull_poly(newt)
        A.add_patch(Polygon(P, closed=True, fc="#e8eefc", ec="k", lw=2))
        A.scatter(*newt.T, color="k", zorder=3)
        A.set_xlim(-0.5, 2.5), A.set_ylim(-0.5, 2.5), A.set_aspect("equal")
        A.set_title(f"{t}\nNewton polygon of {f}")
        B = ax[1, j]
        pts = amoeba_points(polys)
        B.scatter(pts[:, 0], pts[:, 1], s=1, color=BLUE, alpha=0.4)
        B.set_xlim(-6, 6), B.set_ylim(-6, 6), B.set_aspect("equal")
        B.axhline(0, color="#ddd", lw=0.5, zorder=0), B.axvline(0, color="#ddd", lw=0.5, zorder=0)
        B.set_title(f"its zero set in log-radius coordinates:\n{grid}")
        B.set_xlabel("log|x|"), B.set_ylabel("log|y|")
    fig.suptitle("Why the hypercube ideas looked Cartesian: in log-radius coordinates the cube's own zero set IS the "
                 "Cartesian axes.\nOther shapes come with other grids. Switching coordinate grid = switching shape "
                 "(each ray is perpendicular to an edge of the shape above it)", y=1.02)
    save(fig, "fig9-each-shape-has-its-own-grid.png")


# ------------------------------------------------------------- fig 10
def fig10():
    fig, ax = plt.subplots(1, 2, figsize=(13, 6))
    ks = [-0.95, -0.6, 0.0001, 1.0, 5.0]
    cols = plt.get_cmap("viridis")(np.linspace(0, 0.9, len(ks)))
    xs = np.linspace(0, 1, 400)
    for k, c in zip(ks, cols):
        # boundary x o_k y = 1  ->  y = (1 - x)/(1 + k x)
        y = (1 - xs) / (1 + k * xs)
        lab = f"k = {k:g}" if abs(k) > 1e-3 else "k = 0"
        ax[0].fill_between(xs, 0, y, color=c, alpha=0.12)
        ax[0].plot(xs, y, color=c, lw=2, label=lab)
        size = np.log1p(k) / k
        ax[1].fill([0, size, 0], [0, 0, size], color=c, alpha=0.12)
        ax[1].plot([size, 0], [0, size], color=c, lw=2, label=f"{lab}: size {size:.2f}")
    ax[0].plot([0, 1, 1], [1, 1, 0], "--", color=GRAY)
    ax[0].text(0.62, 1.03, "k = −1: the whole square", color=GRAY)
    ax[0].set_xlim(0, 1.1), ax[0].set_ylim(0, 1.15), ax[0].set_aspect("equal")
    ax[0].set_title("Cartesian x, y: the budget  x ∘ₖ y ≤ 1\nbends from square (k = −1) to triangle (k = 0) and thinner")
    ax[0].set_xlabel("x"), ax[0].set_ylabel("y"), ax[0].legend(fontsize=9)
    ax[1].set_xlim(0, 3.2), ax[1].set_ylim(0, 3.2), ax[1].set_aspect("equal")
    ax[1].set_title("Log-radius u = log(1+kx)/k: the SAME budget is always a triangle;\nonly its size moves "
                    "(→ ∞ as k → −1, which is the square)")
    ax[1].set_xlabel("u"), ax[1].set_ylabel("v"), ax[1].legend(fontsize=9)
    fig.suptitle("A shape switch in one coordinate system can be a pure size change in another "
                 "(budget region from the beyond-the-cube pass)", y=1.02)
    save(fig, "fig10-budget-region-one-simplex.png")


if __name__ == "__main__":
    for f in (fig1, fig2, fig3, fig4, fig5, fig6, fig7, fig8, fig9, fig10):
        f()
