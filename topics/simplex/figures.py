"""Pictures for the simplex pass.  Run: python3 figures.py  (writes figures/*.png). Needs numpy, scipy, matplotlib."""
import math
import os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon, Circle, FancyArrowPatch, Rectangle
from scipy import special

OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "figures")
os.makedirs(OUT, exist_ok=True)
rng = np.random.default_rng(7)

BLUE, ORANGE, GREEN, RED, GRAY, PURPLE = "#2a6fdb", "#e8892b", "#2e9e5b", "#c8423b", "#888888", "#7a4fc4"
plt.rcParams.update({"font.size": 11, "axes.titlesize": 12, "figure.dpi": 130})

# barycentric -> 2D for the triangle with corners A=(0,0) [outcome 3], B=(1,0) [outcome 1], C=(1/2, sqrt3/2) [outcome 2]
TRI = np.array([[1.0, 0.0], [0.5, np.sqrt(3) / 2], [0.0, 0.0]])


def b2c(p):
    p = np.atleast_2d(p)
    return p @ TRI


def tri_axes(ax, labels=("1", "2", "3")):
    ax.add_patch(Polygon(TRI, closed=True, fill=False, lw=1.5, color="k"))
    for (x, y), lab, off in zip(TRI, labels, [(0.04, -0.05), (0, 0.05), (-0.06, -0.05)]):
        ax.text(x + off[0], y + off[1], lab, ha="center", va="center", fontsize=10)
    ax.set_aspect("equal")
    ax.set_xlim(-0.12, 1.12)
    ax.set_ylim(-0.12, 0.98)
    ax.axis("off")


def stick(v):
    v = np.atleast_2d(v)
    n = v.shape[1] + 1
    p = np.empty((v.shape[0], n))
    rem = np.ones(v.shape[0])
    for m in range(n - 1):
        p[:, m] = v[:, m] * rem
        rem = rem * (1 - v[:, m])
    p[:, -1] = rem
    return p


def save(fig, name):
    fig.savefig(os.path.join(OUT, name), bbox_inches="tight")
    plt.close(fig)
    print("wrote", name)


# ---------------------------------------------------------------- fig1
def fig1():
    fig = plt.figure(figsize=(14, 4.6))
    ax1 = fig.add_subplot(1, 3, 1)
    ax2 = fig.add_subplot(1, 3, 2)
    ax3 = fig.add_subplot(1, 3, 3, projection="3d")
    ts = np.linspace(0, 1, 201)
    grid = np.linspace(0.1, 0.9, 9)
    # cube
    for g in grid:
        ax1.plot([g, g], [0, 1], color=BLUE, lw=0.9)
        ax1.plot([0, 1], [g, g], color=BLUE, lw=0.9, alpha=0.5)
    ax1.add_patch(Rectangle((0, 0), 1, 1, fill=False, lw=1.5))
    ax1.set_aspect("equal"); ax1.set_xlim(-0.1, 1.1); ax1.set_ylim(-0.1, 1.1)
    ax1.set_xlabel("v₁ = first cut (share taken)"); ax1.set_ylabel("v₂ = second cut")
    ax1.set_title("Cube: two independent cuts\n(stick fractions)")
    # triangle with both grids
    tri_axes(ax2)
    for g in grid:
        P = stick(np.column_stack([np.full_like(ts, g), ts]))
        ax2.plot(*b2c(P).T, color=BLUE, lw=0.9)
        P = stick(np.column_stack([ts, np.full_like(ts, g)]))
        ax2.plot(*b2c(P).T, color=BLUE, lw=0.9, alpha=0.5)
    th = np.linspace(0, np.pi / 2, 201)
    for g in np.linspace(np.pi / 20, np.pi / 2 - np.pi / 20, 9):
        P = stick(np.column_stack([np.full_like(th, np.cos(g) ** 2), np.cos(th) ** 2]))
        ax2.plot(*b2c(P).T, color=ORANGE, lw=0.8, ls="--")
    ax2.set_title("Simplex: probabilities p₁+p₂+p₃ = 1\nblue = cube grid, orange dashed = sphere grid")
    # sphere octant
    u, w = np.meshgrid(np.linspace(0, np.pi / 2, 40), np.linspace(0, np.pi / 2, 40))
    X, Y, Z = np.cos(u), np.sin(u) * np.cos(w), np.sin(u) * np.sin(w)
    ax3.plot_surface(X, Y, Z, color="#f6e3cf", alpha=0.5, linewidth=0)
    for g in np.linspace(np.pi / 20, np.pi / 2 - np.pi / 20, 9):
        ax3.plot(np.full_like(th, np.cos(g)), np.sin(g) * np.cos(th), np.sin(g) * np.sin(th), color=ORANGE, lw=0.9)
        ax3.plot(np.cos(th), np.sin(th) * np.cos(g), np.sin(th) * np.sin(g), color=ORANGE, lw=0.9, alpha=0.5)
    ax3.set_title("Sphere: amplitudes √p\n(angles θ₁, θ₂)")
    ax3.set_xticks([]); ax3.set_yticks([]); ax3.set_zticks([])
    ax3.view_init(25, 35)
    fig.text(0.355, 0.08, "stick-breaking = our law at k = −1  →", ha="center", va="center", fontsize=11, color=BLUE)
    fig.text(0.665, 0.08, "←  square each number (Born)", ha="center", va="center", fontsize=11, color=ORANGE)
    save(fig, "fig1-three-coordinates-one-simplex.png")


# ---------------------------------------------------------------- fig2
def fig2():
    v = [0.3, 0.5, 0.4, 0.6]
    fig, ax = plt.subplots(figsize=(11, 4.2))
    rem, start = 1.0, 0.0
    cols = [BLUE, ORANGE, GREEN, PURPLE, GRAY]
    P = 0.0
    for m, vm in enumerate(v):
        piece = vm * rem
        y = 4 - m
        ax.add_patch(Rectangle((0, y), start, 0.6, color="#dddddd"))
        ax.add_patch(Rectangle((start, y), piece, 0.6, color=cols[m]))
        ax.add_patch(Rectangle((start + piece, y), rem - piece, 0.6, fill=False, ec="k"))
        P = P + vm - P * vm
        ax.text(1.03, y + 0.3, f"cut {m+1}: take v{m+1} = {vm} of what's left\nrunning total = {P:.3f}", va="center", fontsize=10)
        start += piece
        rem -= piece
    ax.add_patch(Rectangle((0, -0.2), 1, 0.6, color="#dddddd"))
    acc = 0.0
    rem = 1.0
    for m, vm in enumerate(v):
        piece = vm * rem
        ax.add_patch(Rectangle((acc, -0.2), piece, 0.6, color=cols[m]))
        acc += piece; rem -= piece
    ax.add_patch(Rectangle((acc, -0.2), rem, 0.6, color=cols[4]))
    ax.text(1.03, 0.1, "the five pieces are a point of the 4-simplex", va="center", fontsize=10)
    ax.set_xlim(-0.02, 1.75); ax.set_ylim(-0.5, 5); ax.axis("off")
    ax.set_title("Breaking a stick: running total = v₁ ∘₋₁ v₂ ∘₋₁ v₃ ∘₋₁ v₄  (a ∘₋₁ b = a + b − ab: what's left multiplies)")
    save(fig, "fig2-stick-breaking-is-the-core-law.png")


# ---------------------------------------------------------------- fig3
def fig3():
    fig, axs = plt.subplots(1, 4, figsize=(17, 4.4))
    # k=-1
    ax = axs[0]; tri_axes(ax)
    V = rng.random((1500, 2))
    ax.scatter(*b2c(stick(V)).T, s=2, color=BLUE, alpha=0.5)
    ax.set_title("k = −1: builds the simplex from the cube\n(even cuts land unevenly: stick-breaking)")
    # k=0
    ax = axs[1]
    ax.plot([0, 1], [0, 0], color="k", lw=2)
    ax.add_patch(Rectangle((0, -0.05), 0.25, 0.1, color=BLUE))
    ax.add_patch(Rectangle((0.25, -0.05), 0.35, 0.1, color=ORANGE))
    ax.text(0.125, 0.12, "0.25", ha="center"); ax.text(0.425, 0.12, "0.35", ha="center")
    ax.text(0.3, -0.25, "0.25 ∘₀ 0.35 = 0.60", ha="center", fontsize=12)
    ax.text(0.5, 0.4, "exclusive outcomes: chances add", ha="center")
    ax.set_xlim(-0.05, 1.05); ax.set_ylim(-0.5, 0.6); ax.axis("off")
    ax.set_title("k = 0: adds the parts\n(one outcome or the other)")
    # k=+1 bayes in triangle
    ax = axs[2]; tri_axes(ax)
    prior = np.array([0.15, 0.25, 0.6])
    lik = np.array([0.8, 0.6, 0.5])
    path = [prior]
    for _ in range(10):
        c = path[-1] * lik; path.append(c / c.sum())
    path = np.array(path)
    ax.plot(*b2c(path).T, "-o", color=GREEN, ms=3)
    lik2 = np.array([0.5, 0.8, 0.6])
    path2 = [prior]
    for _ in range(10):
        c = path2[-1] * lik2; path2.append(c / c.sum())
    path2 = np.array(path2)
    ax.plot(*b2c(path2).T, "-o", color=PURPLE, ms=3)
    ax.set_title("k = +1: moves you inside (Bayes)\nrepeated evidence bends here…")
    # log-ratio plane
    ax = axs[3]
    def clr(p):
        l = np.log(p); return l - l.mean(axis=-1, keepdims=True)
    e1 = np.array([1, -1, 0]) / np.sqrt(2); e2 = np.array([1, 1, -2]) / np.sqrt(6)
    for pth, c in ((path, GREEN), (path2, PURPLE)):
        U = clr(pth)
        ax.plot(U @ e1, U @ e2, "-o", color=c, ms=3)
    ax.set_aspect("equal"); ax.grid(alpha=0.3)
    ax.set_xlabel("log-ratio coordinate 1"); ax.set_ylabel("log-ratio coordinate 2")
    ax.set_title("…but is a straight line of equal steps\nin log-ratio = L₁(odds − 1) coordinates")
    save(fig, "fig3-one-dial-three-jobs.png")


# ---------------------------------------------------------------- fig4
def fig4():
    fig, axs = plt.subplots(1, 4, figsize=(17, 4.2))
    N = 60000
    z = rng.normal(size=(N, 2)) + 1j * rng.normal(size=(N, 2)); z /= np.linalg.norm(z, axis=1, keepdims=True)
    x = rng.normal(size=(N, 2)); x /= np.linalg.norm(x, axis=1, keepdims=True)
    bins = np.linspace(0, 1, 41)
    axs[0].hist(np.abs(z[:, 0]) ** 2, bins=bins, density=True, color=BLUE, alpha=0.8)
    axs[0].axhline(1, color="k", lw=1)
    axs[0].set_title("Complex amplitudes (with phase i):\np = |z₁|² is flat (Archimedes' hat-box)")
    axs[0].set_xlabel("p"); axs[0].set_ylim(0, 4)
    axs[1].hist(x[:, 0] ** 2, bins=bins, density=True, color=ORANGE, alpha=0.8)
    t = np.linspace(0.005, 0.995, 300)
    axs[1].plot(t, 1 / (np.pi * np.sqrt(t * (1 - t))), "k", lw=1)
    axs[1].set_title("Real amplitudes (no phase):\narcsine law 1/(π√(p(1−p)))")
    axs[1].set_xlabel("p"); axs[1].set_ylim(0, 4)
    z3 = rng.normal(size=(5000, 3)) + 1j * rng.normal(size=(5000, 3)); z3 /= np.linalg.norm(z3, axis=1, keepdims=True)
    x3 = rng.normal(size=(5000, 3)); x3 /= np.linalg.norm(x3, axis=1, keepdims=True)
    tri_axes(axs[2]); axs[2].scatter(*b2c(np.abs(z3) ** 2).T, s=1.5, color=BLUE, alpha=0.6)
    axs[2].set_title("3 outcomes, complex: even spread")
    tri_axes(axs[3]); axs[3].scatter(*b2c(x3 ** 2).T, s=1.5, color=ORANGE, alpha=0.6)
    axs[3].set_title("3 outcomes, real: crowds the edges\n(Jeffreys' prior)")
    fig.suptitle("The simplex is the sphere with its phases forgotten: squaring amplitudes lands on it", y=1.04)
    save(fig, "fig4-sphere-with-phases-forgotten.png")


# ---------------------------------------------------------------- fig5
def fig5():
    fig = plt.figure(figsize=(16, 4.6))
    ax = fig.add_subplot(1, 4, 1)
    ns = np.arange(0, 11)
    partial = np.cumsum([1 / math.factorial(n) for n in ns])
    ax.bar(ns, [1 / math.factorial(n) for n in ns], color=BLUE)
    ax2 = ax.twinx(); ax2.plot(ns, partial, "o-", color=ORANGE); ax2.axhline(np.e, color=ORANGE, ls=":")
    ax2.text(6, np.e + 0.06, "e = 2.71828…", color=ORANGE)
    ax.set_xlabel("dimension n"); ax.set_ylabel("corner simplex volume 1/n!", color=BLUE)
    ax2.set_ylabel("running sum", color=ORANGE); ax2.set_ylim(0, 3)
    ax.set_title("e: the corner simplices' volumes\nadd up to e")
    # cube = n! simplices (n=3)
    ax = fig.add_subplot(1, 4, 2, projection="3d")
    import itertools
    perms = list(itertools.permutations(range(3)))
    cols = [BLUE, ORANGE, GREEN, PURPLE, RED, GRAY]
    from mpl_toolkits.mplot3d.art3d import Poly3DCollection
    for c, perm in zip(cols, perms):
        verts = [np.zeros(3)]
        cur = np.zeros(3)
        for i in perm[::-1]:
            cur = cur.copy(); cur[i] = 1; verts.append(cur)
        verts = np.array(verts)
        faces = [[verts[a], verts[b], verts[cc]] for a, b, cc in itertools.combinations(range(4), 3)]
        shift = (verts.mean(axis=0) - 0.5) * 0.35
        ax.add_collection3d(Poly3DCollection([[v + shift for v in f] for f in faces], facecolor=c, alpha=0.35, edgecolor="k", lw=0.4))
    ax.set_xlim(-0.2, 1.2); ax.set_ylim(-0.2, 1.2); ax.set_zlim(-0.2, 1.2)
    ax.set_xticks([]); ax.set_yticks([]); ax.set_zticks([])
    ax.set_title("The cube is 3! = 6 simplices,\none per ordering (pulled apart)")
    # volume sums
    ax = fig.add_subplot(1, 4, 3)
    names = ["corner\nsimplex", "cross-\npolytope", "ball", "cube"]
    vals = [np.e, np.e ** 2, np.exp(np.pi) * (1 + special.erf(np.sqrt(np.pi))), 1e3]
    bars = ax.bar(names, vals, color=[BLUE, PURPLE, ORANGE, GRAY]); ax.set_yscale("log")
    for b, lab in zip(bars, ["e", "e²", "e^π(1+erf√π)\n≈ 46.0", "∞"]):
        ax.text(b.get_x() + b.get_width() / 2, b.get_height() * 1.15, lab, ha="center")
    ax.set_ylim(1, 5e3); ax.set_title("Volume summed over every dimension\n(unit size)")
    # pi: Fisher length
    ax = fig.add_subplot(1, 4, 4)
    th = np.linspace(0, np.pi / 2, 200)
    ax.plot(np.cos(th), np.sin(th), color=ORANGE, lw=2.5)
    ax.plot([0, 1], [1, 0], color=BLUE, lw=2.5)
    for p in np.linspace(0.1, 0.9, 5):
        a = np.sqrt(np.array([p, 1 - p]))
        ax.plot([p, a[0]], [1 - p, a[1]], color=GRAY, lw=0.7)
        ax.plot(a[0], a[1], "o", color=ORANGE, ms=4); ax.plot(p, 1 - p, "o", color=BLUE, ms=4)
    ax.text(0.28, 0.42, "coin simplex\np₁ + p₂ = 1", color=BLUE)
    ax.text(0.62, 0.92, "quarter circle\n√p₁, √p₂", color=ORANGE)
    ax.set_aspect("equal"); ax.set_xlim(-0.05, 1.1); ax.set_ylim(-0.05, 1.1)
    ax.set_title("π: square roots bend the coin's segment\ninto a quarter circle; Fisher length = π")
    save(fig, "fig5-where-e-and-pi-live.png")


# ---------------------------------------------------------------- fig6
def hilbert_simplex(p, q):
    lr = np.log(p) - np.log(q)
    return 0.5 * (lr.max(axis=-1) - lr.min(axis=-1))


def fig6():
    fig, axs = plt.subplots(1, 3, figsize=(15, 4.8))
    # triangle: balls = level sets
    ax = axs[0]; tri_axes(ax)
    g = np.linspace(0.001, 0.999, 500)
    A, B = np.meshgrid(g, g)
    C = 1 - A - B
    mask = C > 0.001
    for cen, col in zip([np.array([1, 1, 1]) / 3, np.array([0.6, 0.25, 0.15]), np.array([0.12, 0.12, 0.76])], [BLUE, ORANGE, GREEN]):
        P = np.stack([A, B, C], axis=-1)
        D = np.where(mask, hilbert_simplex(np.clip(P, 1e-9, 1), cen), np.nan)
        XY = np.stack([A, B, C], axis=-1) @ TRI
        ax.contour(XY[..., 0], XY[..., 1], D, levels=[0.25, 0.5, 0.8], colors=col, linewidths=1)
        ax.plot(*b2c(cen).T, "o", color=col, ms=4)
    ax.set_title("Triangle: balls are hexagons, all alike in\nlog-ratio coordinates (flat; Bayes lives here)")
    # disk
    ax = axs[1]
    ax.add_patch(Circle((0, 0), 1, fill=False, lw=1.5))
    g = np.linspace(-0.999, 0.999, 500)
    X, Y = np.meshgrid(g, g)
    R2 = X ** 2 + Y ** 2
    for cen, col in zip([np.array([0.0, 0.0]), np.array([0.55, 0.2]), np.array([-0.3, -0.7])], [BLUE, ORANGE, GREEN]):
        num = 1 - (X * cen[0] + Y * cen[1])
        den = np.sqrt(np.clip((1 - R2) * (1 - cen @ cen), 1e-12, None))
        D = np.where(R2 < 0.999, np.arccosh(np.clip(num / den, 1, None)), np.nan)
        ax.contour(X, Y, D, levels=[0.25, 0.5, 0.8], colors=col, linewidths=1)
        ax.plot(*cen, "o", color=col, ms=4)
    ax.set_aspect("equal"); ax.axis("off"); ax.set_xlim(-1.1, 1.1); ax.set_ylim(-1.1, 1.1)
    ax.set_title("Disk: hyperbolic geometry\n(curved: relativity's velocity space)")
    # square, numeric via formula along lines
    ax = axs[2]
    ax.add_patch(Rectangle((-1, -1), 2, 2, fill=False, lw=1.5))
    def hd_sq(p, Q):
        d = Q - p
        with np.errstate(divide="ignore", invalid="ignore"):
            tpx = np.where(d[..., 0] > 0, (1 - p[0]) / d[..., 0], np.where(d[..., 0] < 0, (-1 - p[0]) / d[..., 0], np.inf))
            tpy = np.where(d[..., 1] > 0, (1 - p[1]) / d[..., 1], np.where(d[..., 1] < 0, (-1 - p[1]) / d[..., 1], np.inf))
            tp = np.minimum(tpx, tpy)
            tmx = np.where(d[..., 0] < 0, (1 - p[0]) / -d[..., 0], np.where(d[..., 0] > 0, (p[0] + 1) / d[..., 0], np.inf))
            tmy = np.where(d[..., 1] < 0, (1 - p[1]) / -d[..., 1], np.where(d[..., 1] > 0, (p[1] + 1) / d[..., 1], np.inf))
            tm = np.minimum(tmx, tmy)
            return 0.5 * np.log((tp * (1 + tm)) / ((tp - 1) * tm))
    g = np.linspace(-0.995, 0.995, 500)
    X, Y = np.meshgrid(g, g)
    Q = np.stack([X, Y], axis=-1)
    for cen, col in zip([np.array([0.0, 0.0]), np.array([0.6, 0.3]), np.array([-0.5, -0.7])], [BLUE, ORANGE, GREEN]):
        D = hd_sq(cen, Q)
        ax.contour(X, Y, D, levels=[0.25, 0.5, 0.8], colors=col, linewidths=1)
        ax.plot(*cen, "o", color=col, ms=4)
    ax.set_aspect("equal"); ax.axis("off"); ax.set_xlim(-1.1, 1.1); ax.set_ylim(-1.1, 1.1)
    ax.set_title("Square: balls change shape as they move\n(neither flat nor round)")
    fig.suptitle("One distance rule (Hilbert's: log of a cross-ratio) on three shapes: balls of radius 0.25, 0.5, 0.8", y=1.03)
    save(fig, "fig6-one-distance-rule-three-shapes.png")


# ---------------------------------------------------------------- fig7
def einstein(u, v):
    gu = 1 / np.sqrt(1 - u @ u)
    return (u + v / gu + (gu / (1 + gu)) * (u @ v) * u) / (1 + u @ v)


def fig7():
    fig, axs = plt.subplots(1, 3, figsize=(16, 5))
    ax = axs[0]
    ax.add_patch(Circle((0, 0), 1, fill=False, lw=1.5))
    u, v = np.array([0.8, 0.0]), np.array([0.0, 0.8])
    u3, v3 = np.append(u, 0), np.append(v, 0)
    w1, w2 = einstein(u3, v3)[:2], einstein(v3, u3)[:2]
    # Klein model: geodesics are straight lines
    ax.add_patch(Polygon([[0, 0], u, w1], closed=True, color=BLUE, alpha=0.2))
    ax.add_patch(Polygon([[0, 0], v, w2], closed=True, color=ORANGE, alpha=0.2))
    ax.plot([0, u[0], w1[0]], [0, u[1], w1[1]], "-o", color=BLUE, label="u then v")
    ax.plot([0, v[0], w2[0]], [0, v[1], w2[1]], "-o", color=ORANGE, label="v then u")
    ang = np.degrees(np.arccos(w1 @ w2 / (np.linalg.norm(w1) * np.linalg.norm(w2))))
    ax.text(0, -1.25, f"end points differ by a turn of {ang:.1f}°", fontsize=10, ha="center")
    ax.legend(loc="upper left", fontsize=9, bbox_to_anchor=(-0.05, 1.0))
    ax.set_aspect("equal"); ax.axis("off"); ax.set_xlim(-1.1, 1.1); ax.set_ylim(-1.35, 1.1)
    ax.set_title("Ball (velocities, speed 0.8c each):\norder matters, a turn is left over (Thomas)")
    ax = axs[1]; tri_axes(ax)
    p = np.array([1, 1, 1]) / 3
    a, b = np.array([0.7, 0.2, 0.1]), np.array([0.15, 0.7, 0.15])
    def pert(x, y):
        z = x * y; return z / z.sum()
    def seg(x, y, n=30):  # straight line in log-ratio coordinates from x to x(+)y
        out = []
        for t in np.linspace(0, 1, n):
            out.append(pert(x, y ** t))
        return np.array(out)
    s1 = np.vstack([seg(p, a), seg(pert(p, a), b)])
    s2 = np.vstack([seg(p, b), seg(pert(p, b), a)])
    ax.plot(*b2c(s1).T, color=BLUE, lw=2, label="evidence A then B")
    ax.plot(*b2c(s2).T, color=ORANGE, lw=2, ls="--", label="evidence B then A")
    ax.plot(*b2c(pert(pert(p, a), b)).T, "ko", ms=6)
    ax.legend(loc="upper left", fontsize=9, bbox_to_anchor=(-0.1, 1.02))
    ax.set_title("Simplex (beliefs): order never matters,\nthe loop closes, no turn left over")
    ax = axs[2]
    rr = np.random.default_rng(1)
    defects, wig = [], []
    for _ in range(150):
        a1, b1 = rr.uniform(0.05, 0.97, 2)
        gu, gv = 1 / np.sqrt(1 - a1 * a1), 1 / np.sqrt(1 - b1 * b1)
        wig.append(np.arccos((gu + gv) / (1 + gu * gv)))
        # defect via Gauss-Bonnet of the hyperbolic right triangle: here use the same closed form check from checks.py
        defects.append(wig[-1])
    # recompute defects independently (hyperboloid angles)
    def hyp(x):
        g = 1 / np.sqrt(1 - x @ x); return g * np.append(1.0, x)
    def mk(a_, b_):
        return -a_[0] * b_[0] + a_[1:] @ b_[1:]
    def angle_at(A, B, C):
        tB = B + mk(A, B) * A; tC = C + mk(A, C) * A
        return np.arccos(np.clip(mk(tB, tC) / np.sqrt(mk(tB, tB) * mk(tC, tC)), -1, 1))
    rr = np.random.default_rng(1)
    defects, wig, speeds = [], [], []
    for _ in range(150):
        a1, b1 = rr.uniform(0.05, 0.97, 2)
        uu, vv = np.array([a1, 0, 0]), np.array([0, b1, 0])
        A, B, C = hyp(np.zeros(3)), hyp(uu), hyp(einstein(uu, vv))
        defects.append(np.pi - angle_at(A, B, C) - angle_at(B, A, C) - angle_at(C, A, B))
        gu, gv = 1 / np.sqrt(1 - a1 * a1), 1 / np.sqrt(1 - b1 * b1)
        wig.append(np.arccos((gu + gv) / (1 + gu * gv))); speeds.append(a1 * b1)
    sc = ax.scatter(defects, wig, c=speeds, cmap="viridis", s=14)
    m = max(wig); ax.plot([0, m], [0, m], "k:", lw=1)
    ax.set_xlabel("area of the velocity triangle (radians)"); ax.set_ylabel("leftover turn (radians)")
    ax.set_title("Ball: leftover turn = triangle area, exactly\n(our swap rule: turn = d × area, d = −1)")
    plt.colorbar(sc, ax=ax, label="product of the two speeds")
    save(fig, "fig7-relativity-precesses-bayes-does-not.png")


# ---------------------------------------------------------------- fig8
def fig8():
    fig, axs = plt.subplots(1, 3, figsize=(15, 4.8))
    ax = axs[0]
    th = np.linspace(0, 2 * np.pi, 400)
    v, w = 0.5, 0.6
    vw = (v + w) / (1 + v * w)
    rv, rw, rvw = np.arctanh(v), np.arctanh(w), np.arctanh(vw)
    ax.plot([-0.2, 1.6], [0, 0], color="k")
    for x in np.arange(0, 1.61, 0.4):
        ax.plot([x, x], [-0.03, 0.03], color="k", lw=0.8)
    ax.annotate("", xy=(rv, 0.12), xytext=(0, 0.12), arrowprops=dict(arrowstyle="->", color=BLUE, lw=2))
    ax.annotate("", xy=(rvw, 0.12), xytext=(rv, 0.12), arrowprops=dict(arrowstyle="->", color=ORANGE, lw=2))
    ax.text(rv / 2, 0.2, "0.5c", color=BLUE, ha="center"); ax.text((rv + rvw) / 2, 0.2, "+0.6c", color=ORANGE, ha="center")
    ax.plot(rvw, 0, "o", color=GREEN, ms=8)
    ax.text(rvw, -0.12, "0.846c", color=GREEN, ha="center")
    ax.annotate("", xy=(rv, -0.32), xytext=(0, -0.32), arrowprops=dict(arrowstyle="->", color=BLUE, lw=2))
    ax.annotate("", xy=(rvw, -0.32), xytext=(rv, -0.32), arrowprops=dict(arrowstyle="->", color=ORANGE, lw=2))
    ax.text(rv / 2, -0.44, "odds ×3", color=BLUE, ha="center"); ax.text((rv + rvw) / 2, -0.44, "odds ×4", color=ORANGE, ha="center")
    ax.text(rvw + 0.05, -0.32, "odds 12", color=GREEN, va="center")
    ax.text(0.7, 0.5, "velocities: Einstein addition", ha="center")
    ax.text(0.7, -0.62, "a coin's belief p = (1+v)/2: Bayes", ha="center")
    ax.text(0.7, 0.03 - 0.1, "", ha="center")
    ax.set_xlim(-0.25, 1.75); ax.set_ylim(-0.75, 0.65); ax.axis("off")
    ax.set_title("In one dimension they are the same law:\nboth are steps along one line (rapidity = ½ log odds)")
    # temperature dial
    ax = axs[1]
    e1 = np.array([1, -1, 0]) / np.sqrt(2); e2 = np.array([1, 1, -2]) / np.sqrt(6)
    def radius(beta, om):
        lo, hi = 0.0, 1.0
        f = lambda r: (special.logsumexp(-beta * r * om) - np.log(3)) / beta - 1
        while f(hi) < 0: hi *= 2
        for _ in range(80):
            mid = (lo + hi) / 2
            lo, hi = (mid, hi) if f(mid) < 0 else (lo, mid)
        return lo
    cmap = plt.get_cmap("plasma")
    betas = [0.01, 0.5, 1.5, 4, 12, 1000]
    for i, beta in enumerate(betas):
        rs = np.array([radius(beta, np.cos(a) * e1 + np.sin(a) * e2) for a in th])
        rs = rs / rs.min()
        ax.plot(rs * np.cos(th), rs * np.sin(th), color=cmap(i / (len(betas) - 0.5)), label=f"β = {beta}")
    ax.set_aspect("equal"); ax.axis("off"); ax.legend(fontsize=8, loc="lower right")
    ax.set_title("Temperature turns the sphere into the simplex\n(free-energy level set, scaled to inradius 1)")
    ax = axs[2]
    for i, pp in enumerate([1, 1.5, 2, 4, 12]):
        x = np.cos(th); y = np.sin(th)
        nrm = (np.abs(x) ** pp + np.abs(y) ** pp) ** (1 / pp)
        ax.plot(x / nrm, y / nrm, color=cmap(i / 5), label=f"p = {pp}")
    ax.plot([1, 1, -1, -1, 1], [1, -1, -1, 1, 1], color=cmap(0.95), label="p = ∞")
    ax.set_aspect("equal"); ax.axis("off"); ax.legend(fontsize=8, loc="lower right")
    ax.set_title("For comparison, the p dial (polar pass):\ncross-polytope → sphere → cube")
    save(fig, "fig8-coin-law-and-two-roundness-dials.png")


# ---------------------------------------------------------------- fig9
def fig9():
    fig, ax = plt.subplots(figsize=(15, 6.2))
    cols = {"Cube": GRAY, "Simplex": BLUE, "Sphere / ball": ORANGE, "Neither": "#bbbbbb"}
    xs = {"Cube": 0.0, "Simplex": 1.0, "Sphere / ball": 2.0, "Neither": 3.0}
    items = {
        "Cube": ["Core law at k = −1\n(independent events)", "Corner expansion of\n1 + (a∘b∘c)", "Power sets, toric charts"],
        "Simplex": ["Probability of exclusive outcomes\n(k = 0: chances add)", "Bayes' rule\n(k = +1 on odds − 1)",
                    "Gibbs / softmax, temperature", "Entropy, Fisher metric", "Replicator equation (biology)",
                    "Mole fractions, conservation\n(chemistry)", "e (volumes 1/n!)", "Building block of every\ncorner-weight shape"],
        "Sphere / ball": ["i (every unit tangent squares to −1)", "π (half a great circle)", "Special relativity\n(Klein ball)",
                          "Precession (Thomas = area)", "Quantum amplitudes", "Rotations, electromagnetism"],
        "Neither": ["Number theory: torus of angles", "Navier–Stokes stretching", "Newton's G", "Fermat"],
    }
    for k, lst in items.items():
        x = xs[k]
        ax.text(x, 1.0, k, ha="center", fontsize=14, weight="bold", color=cols[k] if k != "Neither" else "#666666")
        for j, s in enumerate(lst):
            y = 0.86 - j * 0.105
            ax.add_patch(Rectangle((x - 0.45, y - 0.045), 0.9, 0.09, color=cols[k], alpha=0.18))
            ax.text(x, y, s, ha="center", va="center", fontsize=9.5)
    ax.annotate("", xy=(0.55, 0.5), xytext=(0.45, 0.5), arrowprops=dict(arrowstyle="->", color=BLUE, lw=2))
    ax.annotate("", xy=(1.45, 0.5), xytext=(1.55, 0.5), arrowprops=dict(arrowstyle="->", color=ORANGE, lw=2))
    ax.text(0.5, 0.56, "stick-\nbreaking", ha="center", fontsize=8, color=BLUE)
    ax.text(1.5, 0.56, "square\n(Born)", ha="center", fontsize=8, color=ORANGE)
    ax.set_xlim(-0.55, 3.55); ax.set_ylim(-0.02, 1.07); ax.axis("off")
    ax.set_title("Which shape each family belongs to: the simplex is the shape of a fixed total", fontsize=13)
    save(fig, "fig9-sorting-board.png")


if __name__ == "__main__":
    for f in (fig1, fig2, fig3, fig4, fig5, fig6, fig7, fig8, fig9):
        f()
