"""Pictures for section 7 of the shape-dials page: topological shapes (DRAFT). Run: python3 topo_figures.py"""
import os
import sys
import json
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch, Circle, Arc
from skimage import measure

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import topo_shapes as T

FIG = os.path.join(HERE, "figures")
TEAL, ORANGE, PURPLE, GREY, BLUE = "#2a9d8f", "#e76f51", "#8d6cab", "#6c757d", "#1f77b4"
LIGHT = np.array([0.5, -0.6, 0.8]); LIGHT /= np.linalg.norm(LIGHT)
plt.rcParams.update({"font.size": 10, "axes.titlesize": 10})
RES = json.load(open(os.path.join(HERE, "topo_results.json")))


def shade(nrm, color, two_sided=True):
    sh = np.abs(nrm @ LIGHT) if two_sided else np.clip(nrm @ LIGHT, 0, 1)
    sh = np.clip(sh, 0.15, 1)
    base = np.array(matplotlib.colors.to_rgb(color))
    return np.clip(base * (0.35 + 0.75 * sh[..., None]), 0, 1)


def psurf(ax, X, Y, Z, color, m=None, view=(28, 35), alpha=1.0):
    P = np.stack([X, Y, Z], -1)
    nrm = np.cross(np.gradient(P, axis=1), np.gradient(P, axis=0))
    nrm /= np.linalg.norm(nrm, axis=-1, keepdims=True) + 1e-12
    fc = shade(nrm, color)
    if alpha < 1:
        fc = np.concatenate([fc, np.full(fc.shape[:-1] + (1,), alpha)], -1)
    ax.plot_surface(X, Y, Z, facecolors=fc, linewidth=0, antialiased=False, shade=False,
                    rcount=X.shape[0], ccount=X.shape[1])
    finish(ax, m if m is not None else np.nanmax(np.abs(P)) * 0.62, view)


def msurf(ax, D, L, color, view=(32, 35), m=None, scale=(1, 1, 1)):
    sp = tuple((L[1] - L[0]) * s for s in scale)
    v, f, n, _ = measure.marching_cubes(D, 0.0, spacing=sp)
    v = v + np.array([L[0] * s for s in scale])
    fn = np.cross(v[f[:, 1]] - v[f[:, 0]], v[f[:, 2]] - v[f[:, 0]])
    fn /= np.linalg.norm(fn, axis=1, keepdims=True) + 1e-12
    tri = ax.plot_trisurf(v[:, 0], v[:, 1], v[:, 2], triangles=f, linewidth=0, antialiased=False, shade=False)
    tri.set_facecolor(shade(fn, color))
    finish(ax, m if m is not None else np.abs(v).max() * 0.62, view)


def finish(ax, m, view):
    ax.set_xlim(-m, m); ax.set_ylim(-m, m); ax.set_zlim(-m, m)
    ax.set_box_aspect((1, 1, 1)); ax.set_axis_off(); ax.view_init(*view)


def sphere(n=80):
    u, v = np.meshgrid(np.linspace(0, 2 * np.pi, 2 * n), np.linspace(0.02, np.pi - 0.02, n))
    return np.sin(v) * np.cos(u), np.sin(v) * np.sin(u), np.cos(v)


# ---------------------------------------------------------------- fig15: the gallery
def fig15():
    fig = plt.figure(figsize=(15, 7.8))
    L = np.linspace(-4.2, 4.2, 170)
    X, Y, Z = np.meshgrid(L, L, L, indexing="ij")
    top = [("Every convex shape\n(any p, any depth)\none type: a ball, χ = 2, no holes", None, TEAL),
           ("Torus\n1 hole, χ = 0", 1, BLUE), ("Genus-2 surface\n2 holes, χ = −2", 2, BLUE),
           ("Genus-3 surface\n3 holes, χ = −4", 3, BLUE)]
    for i, (title, g, col) in enumerate(top):
        ax = fig.add_subplot(2, 4, i + 1, projection="3d")
        if g is None:
            psurf(ax, *sphere(), col, m=0.75)
        else:
            msurf(ax, T.genus_sdf(g, X, Y, Z), L, col, m=0.62 * max(1.5, 1.05 * g + 0.4))
        ax.set_title(title, fontsize=10.5)
    bot = [("Möbius strip\nthe twist itself (a strip with one edge)", T.mobius(), ORANGE, (30, 30)),
           ("Projective plane: 1 twist, χ = 1\n(Steiner's Roman surface; it must\ncross itself in 3D)", T.roman_surface(),
            ORANGE, (24, 40)),
           ("Klein bottle: 2 twists, χ = 0\n(figure-8 picture; it must\ncross itself in 3D)", T.klein_figure8(), ORANGE,
            (30, 30))]
    for i, (title, S, col, view) in enumerate(bot):
        ax = fig.add_subplot(2, 4, 5 + i, projection="3d")
        psurf(ax, *S, col, view=view)
        ax.set_title(title, fontsize=10.5)
    ax = fig.add_subplot(2, 4, 8); ax.set_axis_off()
    ax.text(0.0, 0.92, "Every closed surface is one of:", fontsize=11, weight="bold", transform=ax.transAxes)
    ax.text(0.0, 0.80, "• a sphere with g holes (handles)", fontsize=10.5, transform=ax.transAxes)
    ax.text(0.0, 0.70, "• a sphere with c twists (crosscaps)", fontsize=10.5, transform=ax.transAxes)
    ax.text(0.0, 0.57, "Holes: χ = 2 − 2g.   Twists: χ = 2 − c.", fontsize=10.5, transform=ax.transAxes)
    ax.text(0.0, 0.36, "A hole and a twist together make\nthree twists (Dyck 1888).", fontsize=10.5,
            transform=ax.transAxes)
    ax.text(0.0, 0.18, "g and c are whole numbers you pick,\nnot dials: no turning of p or depth\nchanges them (checks T5, T6).",
            fontsize=10.5, transform=ax.transAxes, color=GREY)
    fig.suptitle("Figure 15. Shape types. All the convex shapes so far are one type; holes and twists are the others "
                 "(checks T1, T4, T5)", fontsize=12.5, y=1.0)
    fig.subplots_adjust(wspace=0.02, hspace=0.25)
    fig.savefig(os.path.join(FIG, "fig15-shape-types-holes-and-twists.png"), dpi=110, bbox_inches="tight"); plt.close(fig)


# ---------------------------------------------------------------- fig16: gluing polygons
LCOL = {"a": BLUE, "b": ORANGE, "c": PURPLE, "a1": BLUE, "b1": ORANGE, "a2": TEAL, "b2": "#c9184a",
        "c1": BLUE, "c2": ORANGE, "c3": TEAL}


def draw_word(ax, word, title, rot=0.0):
    S = T.Surface(word)
    n = S.n
    ang = np.pi / 2 + np.pi / n + rot + 2 * np.pi * np.arange(n + 1) / n
    P = np.stack([np.cos(ang), np.sin(ang)], 1)
    if n == 4:
        P = np.array([[-1, 1], [-1, -1], [1, -1], [1, 1], [-1, 1]]) * 0.8
    ax.fill(P[:, 0], P[:, 1], color="#f1f3f5", zorder=0)
    for i, (l, s) in enumerate(S.edges):
        a, b = P[i], P[i + 1]
        col = LCOL.get(l, GREY)
        ax.plot(*zip(a, b), color=col, lw=2.5)
        st, en = (a, b) if s == 1 else (b, a)
        mid = (a + b) / 2
        d = (en - st) * 0.18
        ax.add_patch(FancyArrowPatch(mid - d, mid + d, arrowstyle="-|>", mutation_scale=16, color=col, lw=2.5))
        out = mid * 1.18 if np.linalg.norm(mid) > 1e-9 else mid
        ax.text(*out, l, color=col, fontsize=12, ha="center", va="center", weight="bold")
    ax.set_xlim(-1.35, 1.35); ax.set_ylim(-1.35, 1.35); ax.set_aspect("equal"); ax.set_axis_off()
    ax.set_title(title, fontsize=10.5)


def fig16():
    fig, axs = plt.subplots(1, 5, figsize=(16, 4.2))
    items = [("a b a' b'", "Torus  (word a b a⁻¹ b⁻¹)\nglue arrows to matching arrows\nχ = 1 − 2 + 1 = 0, 1 hole"),
             ("a b a b'", "Klein bottle  (a b a b⁻¹)\none pair glued with a flip\nχ = 0, 2 twists"),
             ("a b a b", "Projective plane  (a b a b)\nopposite sides glued with a flip\nχ = 2 − 2 + 1 = 1, 1 twist"),
             (T.genus_word(2), "Genus-2 surface\n(a₁b₁a₁⁻¹b₁⁻¹ a₂b₂a₂⁻¹b₂⁻¹)\nχ = 1 − 4 + 1 = −2, 2 holes"),
             (T.crosscap_word(3), "Dyck's surface  (c₁c₁c₂c₂c₃c₃)\nχ = 1 − 3 + 1 = −1, 3 twists\n(= torus with one twist)")]
    for ax, (w, t) in zip(axs, items):
        draw_word(ax, w, t)
    fig.suptitle("Figure 16. Each surface is a polygon we already have with its edges glued in pairs. The word round the "
                 "edge is the whole recipe;\nχ = corners after gluing − edge pairs + 1. Writing two words one after the "
                 "other joins the surfaces, and 2 − χ adds (checks T1, T2)", fontsize=11.5, y=1.08)
    fig.savefig(os.path.join(FIG, "fig16-surfaces-from-glued-polygons.png"), dpi=120, bbox_inches="tight"); plt.close(fig)


# ---------------------------------------------------------------- fig17: curvature is forced
def poincare_polygon(n, alpha):
    rho = T.regular_polygon_radius(n, alpha, -1)
    rE = np.tanh(rho / 2)                       # Euclidean radius of the corners in the disk
    ang = np.pi / 2 + np.pi / n + 2 * np.pi * np.arange(n) / n
    V = rE * np.stack([np.cos(ang), np.sin(ang)], 1)
    pts = []
    for i in range(n):
        a, b = V[i], V[(i + 1) % n]
        # geodesic through a, b: circle orthogonal to the unit circle
        A = np.array([[2 * a[0], 2 * a[1]], [2 * b[0], 2 * b[1]]])
        rhs = np.array([a @ a + 1, b @ b + 1])
        c = np.linalg.solve(A, rhs)
        R = np.sqrt(c @ c - 1)
        t0, t1 = np.arctan2(*(a - c)[::-1]), np.arctan2(*(b - c)[::-1])
        dt = (t1 - t0 + np.pi) % (2 * np.pi) - np.pi
        t = t0 + np.linspace(0, dt, 60)
        pts.append(c + R * np.stack([np.cos(t), np.sin(t)], 1))
    return np.concatenate(pts)


def fig17():
    fig, axs = plt.subplots(1, 4, figsize=(16, 4.6))
    # positive: projective plane square with 180-degree corners = a hemisphere seen from above
    ax = axs[0]
    ax.add_patch(Circle((0, 0), 1, color="#fde2d8", zorder=0)); ax.add_patch(Circle((0, 0), 1, fill=False, color=ORANGE, lw=2.5))
    for k in range(4):
        ax.plot(np.cos(k * np.pi / 2 + np.pi / 4), np.sin(k * np.pi / 2 + np.pi / 4), "o", color="k")
    ax.text(0, 0, "hemisphere", ha="center", va="center", fontsize=10, color=GREY)
    ax.set_title("Projective plane, χ = 1\n4 corners glue to 2 points: each corner 180°\n"
                 "fits only on a sphere (curvature > 0)", fontsize=10)
    # flat
    ax = axs[1]
    sq = np.array([[-1, -1], [1, -1], [1, 1], [-1, 1], [-1, -1]]) * 0.75
    ax.fill(sq[:, 0], sq[:, 1], color="#d8ecf7"); ax.plot(sq[:, 0], sq[:, 1], color=BLUE, lw=2.5)
    ax.set_title("Torus and Klein bottle, χ = 0\n4 corners glue to 1 point: each 90°\nfits only on the flat plane", fontsize=10)
    # hyperbolic octagon and hexagon
    for ax, (n, nm, chi) in zip(axs[2:], [(8, "Genus-2 surface", -2), (6, "Dyck's surface", -1)]):
        ax.add_patch(Circle((0, 0), 1, fill=False, color=GREY, lw=1))
        P = poincare_polygon(n, 2 * np.pi / n)
        ax.fill(P[:, 0], P[:, 1], color="#e6dcef"); ax.plot(P[:, 0], P[:, 1], color=PURPLE, lw=2.5)
        ax.set_title(f"{nm}, χ = {str(chi).replace('-', '−')}\n{n} corners glue to 1 point: each {360 // n}°\n"
                     "fits only in the hyperbolic plane (curvature < 0)", fontsize=10)
    for ax in axs:
        ax.set_xlim(-1.1, 1.1); ax.set_ylim(-1.25, 1.1); ax.set_aspect("equal"); ax.set_axis_off()
    axs[2].text(0, -1.2, "drawn in the Poincaré disk", ha="center", fontsize=9, color=GREY)
    fig.suptitle("Figure 17. Curvature is not a free choice once you pick the holes and twists. The glued corners must "
                 "add to 360°, and Gauss–Bonnet\n(curvature × area = 2π χ) then fixes the sign: positive for χ > 0, flat "
                 "for χ = 0, negative for χ < 0. So the curvature reading d − k²/4 is forced (check T3)",
                 fontsize=11.5, y=1.07)
    fig.savefig(os.path.join(FIG, "fig17-curvature-forced-by-holes.png"), dpi=120, bbox_inches="tight"); plt.close(fig)


# ---------------------------------------------------------------- fig18: p rounds the torus, aspect ratio changes the type
def fig18():
    fig = plt.figure(figsize=(16, 8))
    for i, p in enumerate([1.3, 2.0, 4.0, np.inf]):
        ax = fig.add_subplot(2, 5, i + 1, projection="3d")
        psurf(ax, *T.torus_of_revolution(1.1, 0.42, 120, 60, p=p), BLUE, m=1.05)
        lab = "p = ∞ (square)" if np.isinf(p) else ("p = 2 (round)" if p == 2 else f"p = {p}")
        ax.set_title(f"{lab}\n1 hole, χ = 0", fontsize=10.5)
    ax = fig.add_subplot(2, 5, 5)
    u = np.linspace(0, 2 * np.pi, 200)
    ax.add_patch(plt.Rectangle((-1, -1), 2, 2, color="#d8ecf7"))
    for y in (-1, 1):
        ax.annotate("", (0.15, y), (-0.15, y), arrowprops=dict(arrowstyle="-|>", color=BLUE, lw=2))
    for x in (-1, 1):
        ax.annotate("", (x, 0.15), (x, -0.15), arrowprops=dict(arrowstyle="-|>", color=ORANGE, lw=2))
    ax.set_xlim(-1.3, 1.3); ax.set_ylim(-1.3, 1.3); ax.set_aspect("equal"); ax.set_axis_off()
    ax.set_title("The flat torus is the cube family's\nsquare (p = ∞) with opposite sides\nglued. At p = 2 in 4D it is the\n"
                 "Clifford torus (check T9)", fontsize=10)
    # aspect ratio row: cut-open tori
    rows = RES["T7"]["rows"]
    for i, (q, lab) in enumerate([(2.0, "ring"), (1.3, "ring"), (1.0, "horn: R = r"), (0.55, "spindle"), (0.0, "sphere")]):
        ax = fig.add_subplot(2, 5, 6 + i, projection="3d")
        r = 0.5; R = q * r
        u, v = np.meshgrid(np.linspace(0.25 * np.pi, 2 * np.pi, 120), np.linspace(0, 2 * np.pi, 80))
        rho = R + r * np.cos(v)
        rho = np.where(rho < 0, np.nan, rho)          # keep the outer surface only
        psurf(ax, rho * np.cos(u), rho * np.sin(u), r * np.sin(v), TEAL if q >= 1 or q == 0 else ORANGE, m=0.85,
              view=(25, 20))
        hole = "1 hole, χ = 0" if q > 1 else ("the hole just closes" if q == 1 else "no hole, χ = 2")
        ax.set_title(f"R/r = {q:g} ({lab})\n{hole}", fontsize=10.5)
    fig.suptitle("Figure 18. Top: the roundness dial p reshapes a torus's tube and core from diamond-like to round to square, "
                 "and its hole count never changes (check T6).\nBottom: shrinking the core radius R against the tube "
                 "radius r (cut open to show inside) changes the type only at one point, the horn torus R = r (check T7)",
                 fontsize=11.5, y=1.0)
    fig.subplots_adjust(wspace=0.02, hspace=0.25)
    fig.savefig(os.path.join(FIG, "fig18-p-on-the-torus-and-the-horn.png"), dpi=110, bbox_inches="tight"); plt.close(fig)


# ---------------------------------------------------------------- fig19: quantum gravity tie, genus expansion
def fig19():
    fig = plt.figure(figsize=(15, 4.6))
    L = np.linspace(-3.4, 3.4, 140)
    X, Y, Z = np.meshgrid(L, L, L, indexing="ij")
    items = [(None, r"sphere: weight $g_s^{-2}$"), (1, r"torus: weight $g_s^{0} = 1$"), (2, r"genus 2: weight $g_s^{2}$"),
             (3, r"genus 3: weight $g_s^{4}$")]
    for i, (g, t) in enumerate(items):
        ax = fig.add_subplot(1, 5, i + 1, projection="3d")
        if g is None:
            psurf(ax, *sphere(50), BLUE, m=0.75)
        else:
            msurf(ax, T.genus_sdf(g, X, Y, Z), L, BLUE, m=0.62 * max(1.5, 1.05 * g + 0.4))
        ax.set_title(t, fontsize=11)
    ax = fig.add_subplot(1, 5, 5); ax.set_axis_off()
    ax.text(0, 0.85, "In string theory each surface\nshape adds to the answer with\n" + r"weight $g_s^{-\chi} = g_s^{2g-2}$.",
            fontsize=10.5, transform=ax.transAxes, va="top")
    ax.text(0, 0.40, r"Each extra hole multiplies by $g_s^2$:" + "\nthe same odds, once per hole. Holes\nadd, so their odds multiply.",
            fontsize=10.5, transform=ax.transAxes, va="top")
    ax.text(0, 0.02, "Known: 't Hooft 1974; Polyakov 1981.", fontsize=9.5, color=GREY, transform=ax.transAxes)
    fig.suptitle("Figure 19. Quantum gravity test: in 2D gravity and string theory the hole count is summed over, and "
                 "each hole costs a fixed odds " + r"$g_s^2$." + "\nThat matches holes being whole-number picks with 2 − χ adding, "
                 "not a continuous dial", fontsize=11.5, y=1.04)
    fig.subplots_adjust(wspace=0.02)
    axs3 = fig.axes[:4]
    for a, b in zip(axs3[:-1], axs3[1:]):
        pa, pb = a.get_position(), b.get_position()
        fig.text((pa.x1 + pb.x0) / 2, (pa.y0 + pa.y1) / 2 - 0.04, "+", fontsize=24, ha="center", va="center")
    fig.savefig(os.path.join(FIG, "fig19-quantum-gravity-genus-expansion.png"), dpi=110, bbox_inches="tight"); plt.close(fig)


if __name__ == "__main__":
    for f in (fig15, fig16, fig17, fig18, fig19):
        f(); print("done", f.__name__)
