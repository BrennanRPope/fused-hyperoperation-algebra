"""Pictures for section 6 of the shape-dials page: is depth a new dial? (DRAFT). Run: python3 depth_figures.py"""
import sys, os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon as MPoly
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import figures as F                     # reuse the 3D surface helper and depth shapes
from facedial import FaceDial, regular_polygon, circle_dirs, lse_radial, cube, cross_polytope, dodecahedron
from wythoff import roots_for, chamber_corners, orbit, wythoff_polytope

FIG = os.path.join(HERE, "figures")
DIRS = F.DIRS
TEAL, ORANGE, PURPLE, GREY = "#2a9d8f", "#e76f51", "#8d6cab", "#6c757d"


def rV(V, U=DIRS):
    """Radial function of conv(V) about its centroid along directions U."""
    from scipy.spatial import ConvexHull
    V = np.asarray(V, float) - np.mean(V, 0)
    hull = ConvexHull(V)
    A, b = hull.equations[:, :-1], -hull.equations[:, -1]
    return 1.0 / (U @ A.T / b).max(1)


def draw(ax, V, color, title):
    F.surf(ax, rV(V), color)
    ax.set_title(title, fontsize=10)


def sym(ax, s):
    ax.text(0.5, 0.5, s, fontsize=30, ha="center", va="center", transform=ax.transAxes); ax.set_axis_off()


# ---------------------------------------------------------------- fig10: depth is a Minkowski sum
def fig10():
    fig = plt.figure(figsize=(14, 7.6))
    rows = [("B3", "cube", "cuboctahedron", "truncated cube", TEAL),
            ("A3", "tetrahedron", "octahedron", "truncated tetrahedron", ORANGE)]
    for i, (fam, a, b, c, col) in enumerate(rows):
        R = roots_for(fam); C = chamber_corners(R)
        A_, B_ = orbit(R, C[0]), orbit(R, C[1])
        S = orbit(R, 0.5 * C[0] + 0.5 * C[1])
        axs = [fig.add_subplot(2, 5, i * 5 + k + 1, projection="3d") if k % 2 == 0 else fig.add_subplot(2, 5, i * 5 + k + 1)
               for k in range(5)]
        draw(axs[0], A_, col, f"½ × {a}\n(corner 0, odds 1)")
        sym(axs[1], "+")
        draw(axs[2], B_, col, f"½ × {b}\n(corner 1, odds 1)")
        sym(axs[3], "=")
        draw(axs[4], S, col, f"{c}\n(depth ½ on the dial)")
    fig.suptitle("Figure 10. Depth is not a new ingredient: the shape at depth ½ is the Minkowski sum of the two corner "
                 "shapes, each weighted ½.\nAny depth k + f is (1 − f)·corner_k + f·corner_(k+1): weights are chances, "
                 "f = t/(1 + t) for odds t (checks M1, M2, O1)", fontsize=11.5, y=1.0)
    fig.subplots_adjust(wspace=0.0, hspace=0.18)
    fig.savefig(os.path.join(FIG, "fig10-depth-is-a-minkowski-sum.png"), dpi=110, bbox_inches="tight"); plt.close(fig)


# ---------------------------------------------------------------- fig11: two readings of one slider (2D, square family)
def fig11():
    sq = np.array([[1, 1], [-1, 1], [-1, -1], [1, -1]], float) / np.sqrt(2)           # circumradius 1
    di = np.array([[1, 0], [0, 1], [-1, 0], [0, -1]], float)
    def hull2(P):
        from scipy.spatial import ConvexHull
        h = ConvexHull(P); return P[h.vertices]
    def msum(A, B, f):
        return hull2(np.array([(1 - f) * a + f * b for a in A for b in B]))
    def inter(f):
        # face side: square's faces (normals ±x, ±y) at offset s0, diamond's faces (normals (±1,±1)/√2) at s1
        R = roots_for("I2(n)", 4); C = chamber_corners(R)
        x = (1 - f) * C[0] + f * C[1]
        from scipy.spatial import HalfspaceIntersection
        Ns = []
        for j in range(2):
            for nvec in orbit(R, C[j]):
                Ns.append(np.append(nvec, -(x @ C[j])))
        hs = HalfspaceIntersection(np.array(Ns), np.zeros(2))
        return hull2(hs.intersections), x, C
    R = roots_for("I2(n)", 4); C = chamber_corners(R)
    A_, B_ = orbit(R, C[0]), orbit(R, C[1])
    fs = [0, 0.25, 0.5, 0.75, 1]
    fig, axs = plt.subplots(2, 5, figsize=(14, 6.2))
    for j, f in enumerate(fs):
        ax = axs[0, j]
        P = msum(A_, B_, f)
        ax.add_patch(MPoly(hull2(A_) * (1 - f), fill=False, ec=TEAL, ls="--", lw=1))
        ax.add_patch(MPoly(hull2(B_) * f, fill=False, ec=PURPLE, ls="--", lw=1))
        ax.add_patch(MPoly(P, fc=ORANGE, alpha=0.35, ec=ORANGE, lw=2))
        ax.set_title(f"{1 - f:g} × square + {f:g} × diamond", fontsize=9.5)
        ax = axs[1, j]
        Q, x, Cc = inter(f)
        s0, s1 = x @ Cc[0], x @ Cc[1]
        from scipy.spatial import HalfspaceIntersection
        for jj, colj in ((0, TEAL), (1, PURPLE)):
            Nj = np.array([np.append(nv, -(x @ Cc[jj])) for nv in orbit(R, Cc[jj])])
            Pj = hull2(HalfspaceIntersection(Nj, np.zeros(2)).intersections)
            ax.add_patch(MPoly(Pj, fill=False, ec=colj, ls="--", lw=1))
        ax.add_patch(MPoly(Q, fc=TEAL, alpha=0.35, ec=TEAL, lw=2))
        ax.set_title(f"face side: offsets {s0:.2f} and {s1:.2f}", fontsize=9.5)
        for a in axs[:, j]:
            a.set_xlim(-1.15, 1.15); a.set_ylim(-1.15, 1.15); a.set_aspect("equal"); a.set_xticks([]); a.set_yticks([])
    axs[0, 0].set_ylabel("corner side:\nadd shapes (Minkowski)", fontsize=10)
    axs[1, 0].set_ylabel("face side:\ncut by both face sets", fontsize=10)
    fig.suptitle("Figure 11. One slider, two existing readings (square family). Top: the corner shapes added with weights "
                 "1 − f and f.\nBottom: the same shapes as the intersection of the two face sets, each pushed out by a "
                 "ratio (checks M2, H1). Middle column: the regular octagon", fontsize=11.5, y=1.02)
    fig.savefig(os.path.join(FIG, "fig11-depth-two-readings.png"), dpi=115, bbox_inches="tight"); plt.close(fig)


# ---------------------------------------------------------------- fig12: depth is a path on the odds triangle
def fig12():
    R = roots_for("B3"); C = chamber_corners(R)
    tri = np.array([[0, 0], [1, 0], [0.5, np.sqrt(3) / 2]])          # barycentric picture of the chamber
    pts = {
        (1, 0, 0): "cube", (0, 1, 0): "cuboctahedron", (0, 0, 1): "octahedron",
        (.5, .5, 0): "truncated cube", (0, .5, .5): "truncated\noctahedron", (.5, 0, .5): "rhombicuboctahedron",
        (1 / 3, 1 / 3, 1 / 3): "truncated\ncuboctahedron",
    }
    fig = plt.figure(figsize=(11.5, 10))
    ax = fig.add_axes([0.08, 0.08, 0.84, 0.82])
    ax.add_patch(MPoly(tri, fc="#f1f3f5", ec="#adb5bd", lw=1.5))
    ax.plot(*tri[[0, 1]].T, color=ORANGE, lw=6, alpha=0.6, solid_capstyle="round")
    ax.plot(*tri[[1, 2]].T, color=ORANGE, lw=6, alpha=0.6, solid_capstyle="round", label="the depth path (depth 0 → 1 → 2)")
    ax.set_xlim(-0.25, 1.25); ax.set_ylim(-0.2, 1.08); ax.set_aspect("equal"); ax.set_axis_off()
    ax.legend(loc="upper left", frameon=False, fontsize=11)
    ax.text(0.5, -0.17, "Each point is a set of odds on the three corners of the mirror chamber (weights summing to 1).\n"
            "Weighted Minkowski sums of the three corner shapes reach every point; the depth slider only walks two edges.",
            ha="center", fontsize=10.5)
    for lam, nm in pts.items():
        lam = np.array(lam, float)
        xy = lam @ tri
        on_path = lam[2] == 0 or lam[0] == 0
        ax.plot(*xy, "o", color=ORANGE if on_path else PURPLE, ms=9, zorder=5)
        fx, fy = ax.transData.transform(xy)
        fx, fy = fig.transFigure.inverted().transform((fx, fy))
        sub = fig.add_axes([fx - 0.06, fy + 0.012, 0.12, 0.12], projection="3d")
        sub.set_facecolor("none"); sub.patch.set_alpha(0)
        F.surf(sub, rV(orbit(R, lam @ C)), TEAL if on_path else PURPLE)
        sub.set_title(nm + ("" if on_path else "\n(off the depth path)"), fontsize=9, y=0.95)
    fig.suptitle("Figure 12. The cube family on its odds triangle. Depth is one route across it (orange);\nthe "
                 "rhombicuboctahedron and the truncated cuboctahedron need odds on all three corners (check I5)",
                 fontsize=12, y=0.99)
    fig.savefig(os.path.join(FIG, "fig12-depth-on-the-odds-triangle.png"), dpi=110, bbox_inches="tight"); plt.close(fig)


# ---------------------------------------------------------------- fig13: independence of p and depth
def fig13():
    fig = plt.figure(figsize=(14, 10.5))
    cubeFD = FaceDial(cube(3), "cube")
    row1 = [(1, "p = 1"), (2, "p = 2 (round)"), (4, "p = 4"), (np.inf, "p = ∞")]
    for k, (p, lab) in enumerate(row1):
        ax = fig.add_subplot(3, 4, k + 1, projection="3d")
        F.surf(ax, F.rad(cubeFD, DIRS, p), "#1f77b4")
        ax.set_title(f"cube faces, {lab}", fontsize=10)
    for k, t in enumerate([0, 0.5, 1, 1.5]):
        ax = fig.add_subplot(3, 4, 5 + k, projection="3d")
        F.surf(ax, rV(wythoff_polytope("B3", t)), TEAL)
        ax.set_title(f"depth {t}: {F.WNAMES[('B3', t)]}", fontsize=10)
    crossings = [("cube", cube(3), "octahedron (B3 depth 2)"), ("octahedron", cross_polytope(3), "cuboctahedron (B3 depth 1)"),
                 ("dodecahedron", dodecahedron(), "icosidodecahedron (H3 depth 1)")]
    for k, (nm, V, res) in enumerate(crossings):
        ax = fig.add_subplot(3, 4, 9 + k, projection="3d")
        fd = FaceDial(V, nm)
        F.surf(ax, fd.radial(DIRS @ np.eye(3), 1), PURPLE)
        ax.set_title(f"{nm} faces at p = 1\n= {res}", fontsize=10)
    ax = fig.add_subplot(3, 4, 12); ax.set_axis_off()
    ax.text(0.0, 0.5, "Row 1: for any p between 1 and ∞, p never\nmakes a flat face (strictly convex, check I1).\n\n"
                      "Row 2: depth never makes a round shape;\nevery depth shape is a polytope (I4).\n\n"
                      "Row 3: p = 1 gives the polar of a zonotope,\nso it lands on a few whole depths, but never on\n"
                      "a truncation, whose polar has triangle faces (I2, I3).", fontsize=10, va="center")
    fig.text(0.02, 0.83, "p alone", rotation=90, fontsize=12, va="center")
    fig.text(0.02, 0.52, "depth alone", rotation=90, fontsize=12, va="center")
    fig.text(0.02, 0.2, "where they cross", rotation=90, fontsize=12, va="center")
    fig.suptitle("Figure 13. p and depth are independent of each other: each reaches shapes the other cannot", fontsize=12,
                 y=0.995)
    fig.subplots_adjust(wspace=0.0, hspace=0.2)
    fig.savefig(os.path.join(FIG, "fig13-p-and-depth-independent.png"), dpi=110, bbox_inches="tight"); plt.close(fig)


# ---------------------------------------------------------------- fig14: the simplex gap, settled
def fig14():
    tri = FaceDial(regular_polygon(3), "triangle")
    th, U2 = circle_dirs(2048)

    def modes(rr):
        rr = rr / rr.mean()
        return [np.mean(rr * np.exp(-1j * n * th)) for n in (3, 6)]
    ph = modes(tri.radial_polytope(U2))[0]; ph /= abs(ph)
    fr = lambda rr: (lambda c: ((c[0] / ph).real, (c[1] / ph ** 2).real))(modes(rr))
    ps = np.geomspace(2.0001, 400, 60); bs = np.geomspace(0.02, 400, 80)
    fp = np.array([fr(tri.radial(U2, p)) for p in ps])
    te = np.array([fr(lse_radial(tri, U2, b)) for b in bs])
    fig, axs = plt.subplots(1, 3, figsize=(15.5, 4.9))
    ax = axs[0]
    ax.plot(fp[:, 0], fp[:, 1], color=TEAL, lw=2.2, label="face power (offset ½), p from 2 to ∞")
    ax.plot(te[:, 0], te[:, 1], color=ORANGE, lw=2.2, ls="--", label="temperature, β from 0 to ∞")
    for p in [2.0001, 3, 5, 10, 30]:
        c = fr(tri.radial(U2, p)); ax.plot(*c, "o", color=TEAL, ms=5)
        ax.annotate(f"p={p:.0f}", c, textcoords="offset points", xytext=(4, -11), fontsize=8, color=TEAL)
    tc = fr(tri.radial_polytope(U2)); ax.plot(*tc, "k*", ms=12); ax.annotate("triangle", tc, xytext=(-40, 8),
                                                                              textcoords="offset points")
    ax.annotate("circle", (0, 0), xytext=(6, 6), textcoords="offset points")
    ax.set_xlabel("3-fold content of the outline"); ax.set_ylabel("6-fold content of the outline")
    ax.set_title("A. Two different curves from circle to triangle:\none family would give one curve (check S4)")
    ax.legend(frameon=False, fontsize=8.5); ax.grid(alpha=0.3)
    ax = axs[1]
    o = np.argsort(te[:, 0])
    gap = fp[:, 1] - np.interp(fp[:, 0], te[o, 0], te[o, 1])
    ok = (fp[:, 0] > te[:, 0].min()) & (fp[:, 0] < te[:, 0].max())
    ax.semilogx(ps[ok], np.abs(gap[ok]), color=PURPLE, lw=2.2)
    ax.set_xlabel("p"); ax.set_ylabel("6-fold gap at equal 3-fold content")
    ax.set_title("B. The gap is zero only at the two ends\n(round, p = 2; triangle, p = ∞)")
    ax.grid(alpha=0.3, which="both")
    ax = axs[2]
    As = np.array([1, 2, 4, 8, 16, 32, 64, 128, 256]); beta = 3.0
    r_l = lse_radial(tri, U2, beta)
    errs = [np.abs(FaceDial(regular_polygon(3), "t", a=a).radial(U2, beta * a) - r_l).max() for a in As]
    ax.loglog(As, errs, "o-", color=TEAL, lw=2, label="face power with p = 3a vs temperature β = 3")
    ax.loglog(As, errs[0] * As[0] / As, "k:", label="1/a")
    ax.axvline(0.5, color=GREY, lw=0.8)
    ax.set_xlabel("offset a"); ax.set_ylabel("largest radius gap"); ax.legend(frameon=False, fontsize=8.5)
    ax.set_title("C. Exact in one limit only: a → ∞ with p = β·a\n((1 + y/a)^(βa) → e^(βy), Lean: face_to_temperature)")
    ax.grid(alpha=0.3, which="both")
    fig.suptitle("Figure 14. The simplex gap, settled: the round-to-triangle (temperature) route is a different curve "
                 "from face power at the triangle's own offset; the two meet only as a → ∞", y=1.03, fontsize=12)
    fig.savefig(os.path.join(FIG, "fig14-simplex-gap-settled.png"), dpi=125, bbox_inches="tight"); plt.close(fig)


if __name__ == "__main__":
    which = sys.argv[1:] or ["10", "11", "12", "13", "14"]
    for w in which:
        globals()["fig" + w]()
        print("fig", w, "done", flush=True)
