"""DRAFT figures for hypersphere/README.md. Run: python3 figures.py  (numpy, scipy, matplotlib)."""
import math
import os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch, Polygon
from scipy import special

OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "figures")
os.makedirs(OUT, exist_ok=True)
rng = np.random.default_rng(7)
KNOWN, OURS, GREY, RED = "#d9822b", "#2b6cb0", "#888888", "#c53030"
plt.rcParams.update({"font.size": 11, "axes.titlesize": 12})


def save(fig, name):
    fig.savefig(os.path.join(OUT, name), dpi=140, bbox_inches="tight")
    plt.close(fig)
    print("wrote", name)


def sphere_wire(ax, R=1.0, alpha=0.12):
    u, v = np.mgrid[0:2 * np.pi:40j, 0:np.pi:20j]
    ax.plot_wireframe(R * np.cos(u) * np.sin(v), R * np.sin(u) * np.sin(v), R * np.cos(v),
                      color=GREY, alpha=alpha, linewidth=0.6)


# ---------------------------------------------------------------- fig1 i and pi
fig = plt.figure(figsize=(15, 5))
ax = fig.add_subplot(1, 3, 1)
t = np.linspace(0, 2 * np.pi, 400)
ax.plot(np.cos(t), np.sin(t), color=GREY)
th = np.linspace(0, np.pi, 200)
ax.plot(np.cos(th), np.sin(th), color=OURS, lw=4, alpha=0.7)
ax.annotate("", xy=(1, 0.75), xytext=(1, 0), arrowprops=dict(arrowstyle="-|>", color=RED, lw=2.5))
ax.text(1.06, 0.4, "i = the one step\nthat stays on\nthe circle", color=RED)
ax.plot([1, -1], [0, 0], "ko")
ax.text(1.05, -0.15, "1 (lift of 0)")
ax.text(-1.45, -0.18, "−1")
ax.text(-0.45, 1.1, "walk length π reaches −1", color=OURS)
ax.set_aspect("equal"); ax.set_xlim(-1.6, 2.1); ax.set_ylim(-1.3, 1.4); ax.axis("off")
ax.set_title("A. Circle S¹: i is the tangent at 1,\nπ is the distance from 1 to −1")

ax = fig.add_subplot(1, 3, 2, projection="3d")
sphere_wire(ax)
for lab, p in (("i", (1, 0, 0)), ("j", (0, 1, 0)), ("k", (0, 0, 1))):
    ax.quiver(0, 0, 0, *p, color=RED, lw=2, arrow_length_ratio=0.12)
    ax.text(*(1.15 * np.array(p)), lab, color=RED, fontsize=14)
pts = rng.normal(size=(60, 3)); pts /= np.linalg.norm(pts, axis=1)[:, None]
ax.scatter(*pts.T, s=8, color=OURS)
ax.set_title("B. 3-sphere S³: the steps that stay on it\nform a whole sphere of i's (each squares to −1)")
ax.text2D(0.05, 0.02, "i·j = k but j·i = −k: order matters", transform=ax.transAxes, color=RED)
ax.set_axis_off(); ax.set_box_aspect((1, 1, 1))

ax = fig.add_subplot(1, 3, 3)
for a in np.linspace(0, 2 * np.pi, 16, endpoint=False):
    ax.annotate("", xy=(2.6 * np.cos(a), 2.6 * np.sin(a)), xytext=(0, 0),
                arrowprops=dict(arrowstyle="-|>", color=OURS, lw=1.2))
for frac, lab in ((0.25, "π/4"), (0.5, "π/2"), (0.75, "3π/4")):
    r = math.tan(frac * math.pi / 2)
    ax.plot(r * np.cos(t), r * np.sin(t), color=GREY, lw=0.6, ls=":")
    ax.text(r * np.cos(0.3) + 0.05, r * np.sin(0.3), lab, fontsize=9, color=GREY)
ax.plot(0, 0, "ko"); ax.text(0.08, 0.08, "1")
ax.text(-2.1, -2.95, "each ray = one i; dotted rings = distance from 1 on S³;\nwalking length π along a ray reaches −1, at infinity here", fontsize=10)
ax.set_aspect("equal"); ax.set_xlim(-2.8, 2.8); ax.set_ylim(-3.2, 2.8); ax.axis("off")
ax.set_title("C. S³ flattened (seen from −1): every\ndirection out of 1 is an i, every one reaches −1 at π")
save(fig, "fig1-i-and-pi-on-the-sphere.png")

# ---------------------------------------------------------------- fig2 k route vs sphere route
fig, ax = plt.subplots(figsize=(13, 4.6))
ax.axis("off")
rows = [("", "k route (Nov 25 page)", "Sphere route (this page)", "Verdict"),
        ("i", "chain ∘ₖ steps at homothety k; size is\n1 + Re(k)·ε; keeping size forces Re k = 0",
         "the sphere is the set of size-1 lifts;\nits tangent at 1 is i, and i² = −1", "same fact, one line shorter"),
        ("π", "the k = i flow first lands on the\nwall's mirror point at x = π",
         "distance from 1 to the opposite\npoint −1, in every dimension", "same fact, shorter"),
        ("new", "one i", "S³: a sphere of i's, ij = k ≠ ji;\na 2π turn in space = −1 (spin)", "new content"),
        ("√π", "not reached", "Gaussian: the circle's 2π gives\n∫e^{−x²} = √π (Poisson)", "new, known")]
colw = [0.06, 0.32, 0.34, 0.2]
y = 0.95
for r_i, row in enumerate(rows):
    x = 0.0
    for c_i, cell in enumerate(row):
        col = "black"
        if c_i == 3 and r_i > 0:
            col = OURS if "new" in cell else KNOWN
        ax.text(x, y, cell, va="top", fontsize=11, weight="bold" if r_i == 0 or c_i == 0 else "normal", color=col)
        x += colw[c_i]
    y -= 0.22 if r_i else 0.12
ax.set_title("Figure 2. i and π: the sphere states the k route's facts more simply; the new content comes from S³")
save(fig, "fig2-k-route-vs-sphere-route.png")

# ---------------------------------------------------------------- fig3 e from cube and simplex
fig = plt.figure(figsize=(15, 4.8))
ax = fig.add_subplot(1, 3, 1)
ax.add_patch(Polygon([[0, 0], [1, 0], [1, 1], [0, 1]], fill=False, ec="black"))
ax.add_patch(Polygon([[0, 0], [1, 0], [0, 1]], color=OURS, alpha=0.5))
ax.add_patch(Polygon([[0, 0], [1, 1], [0, 1]], color=KNOWN, alpha=0.25))
ax.text(0.15, 0.2, "x+y ≤ 1\narea 1/2!", color="white", weight="bold")
ax.text(0.15, 0.75, "x < y\n(sorted)", color=KNOWN)
ax.set_aspect("equal"); ax.axis("off"); ax.set_xlim(-0.1, 1.1); ax.set_ylim(-0.1, 1.1)
ax.set_title("A. The square holds 2! triangles;\nthe corner one has area 1/2!")
ax = fig.add_subplot(1, 3, 2, projection="3d")
for s, e in [((0, 0, 0), (1, 0, 0)), ((0, 0, 0), (0, 1, 0)), ((0, 0, 0), (0, 0, 1)), ((1, 0, 0), (1, 1, 0)),
             ((1, 0, 0), (1, 0, 1)), ((0, 1, 0), (1, 1, 0)), ((0, 1, 0), (0, 1, 1)), ((0, 0, 1), (1, 0, 1)),
             ((0, 0, 1), (0, 1, 1)), ((1, 1, 0), (1, 1, 1)), ((1, 0, 1), (1, 1, 1)), ((0, 1, 1), (1, 1, 1))]:
    ax.plot(*zip(s, e), color="black", lw=0.8)
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
V = np.array([[0, 0, 0], [1, 0, 0], [0, 1, 0], [0, 0, 1]])
ax.add_collection3d(Poly3DCollection([V[[0, 1, 2]], V[[0, 1, 3]], V[[0, 2, 3]], V[[1, 2, 3]]], color=OURS, alpha=0.45))
ax.set_title("B. The cube holds 3! = 6 corner-to-corner\ntetrahedra; the corner one has volume 1/3!")
ax.set_axis_off(); ax.set_box_aspect((1, 1, 1))
ax = fig.add_subplot(1, 3, 3)
vals = [1 / math.factorial(n) for n in range(0, 9)]
cum = np.cumsum(vals)
ax.bar(range(9), vals, color=OURS, label="corner simplex volume 1/n!")
ax.plot(range(9), cum, "o-", color=KNOWN, label="running total")
ax.axhline(math.e, color=RED, ls="--"); ax.text(4.5, math.e + 0.07, "e = 2.71828…", color=RED)
ax.set_xlabel("dimension n"); ax.set_ylim(0, 3.1); ax.legend(loc="center right")
ax.set_title("C. Add the corner simplex of every\ndimension's cube: the total is e")
save(fig, "fig3-e-from-the-cube-corner.png")

# ---------------------------------------------------------------- fig4 random-sum puzzle
M = 400000
total = np.zeros(M); cnt = np.zeros(M, int); alive = np.ones(M, bool)
while alive.any():
    total[alive] += rng.random(alive.sum()); cnt[alive] += 1; alive = total <= 1
fig, axs = plt.subplots(1, 2, figsize=(12, 4.3))
ks = np.arange(2, 9)
emp = [np.mean(cnt == k) for k in ks]
exact = [(k - 1) / math.factorial(k) for k in ks]
axs[0].bar(ks - 0.2, emp, 0.4, color=OURS, label="simulated (400,000 runs)")
axs[0].bar(ks + 0.2, exact, 0.4, color=KNOWN, label="exact (n−1)/n!")
axs[0].set_xlabel("tries needed to pass 1"); axs[0].legend()
axs[0].set_title(f"A. Keep adding random numbers in [0,1]\nuntil the total passes 1: average = {cnt.mean():.4f} ≈ e")
N = np.arange(1, 7)
axs[1].plot(N, [1 / math.factorial(n) for n in N], "o-", color=OURS)
axs[1].set_yscale("log"); axs[1].set_xlabel("n")
axs[1].set_ylabel("chance the first n numbers sum below 1")
axs[1].set_title("B. That chance is exactly the corner\nsimplex's volume 1/n!, so the average is Σ 1/n! = e")
save(fig, "fig4-random-sum-puzzle.png")

# ---------------------------------------------------------------- fig5 e from the sphere: e^pi and the bell curve
fig, axs = plt.subplots(1, 2, figsize=(13, 4.5))
Vb = lambda n: math.pi ** (n / 2) / special.gamma(n / 2 + 1)
ns = np.arange(0, 21)
axs[0].bar(ns, [Vb(n) for n in ns], color=[OURS if n % 2 == 0 else GREY for n in ns])
axs[0].set_xlabel("dimension n"); axs[0].set_ylabel("volume of the unit n-ball")
axs[0].set_title("A. Unit-ball volumes peak at n = 5. Even ones (blue)\nare π^m/m! = angles × simplex, and they sum to e^π, not e")
axs[0].text(9, 4.2, f"Σ even = e^π = {math.exp(math.pi):.3f}\nto get e: shrink to disc area 1", color=OURS)
t = np.linspace(-4, 4, 800)
for n, c in ((5, "#9ec5e8"), (20, "#5a9bd5"), (200, OURS)):
    d = np.clip(1 - t * t / n, 0, None) ** ((n - 3) / 2); d /= np.trapezoid(d, t)
    axs[1].plot(t, d, color=c, label=f"one coordinate on the sphere of radius √{n}")
axs[1].plot(t, np.exp(-t * t / 2) / math.sqrt(2 * math.pi), "k--", label="bell curve e^{−t²/2}/√(2π)")
axs[1].legend(fontsize=9); axs[1].set_xlabel("t")
axs[1].set_title("B. In high dimension the sphere's (1 − t²/n)^{n/2}\nbecomes e^{−t²/2}: the ∘ₖ chain again, at k = −½")
save(fig, "fig5-e-on-the-sphere.png")

# ---------------------------------------------------------------- fig6 ball <-> simplex, polydisc <-> cube
fig = plt.figure(figsize=(14, 5))
Z = rng.normal(size=(6000, 6)); Z /= np.linalg.norm(Z, axis=1)[:, None]
s_sph = Z[:, 0::2] ** 2 + Z[:, 1::2] ** 2
rad = rng.random(6000) ** (1 / 6)
s_ball = s_sph * rad[:, None] ** 2
s_poly = rng.random((6000, 3))
for i, (S, ttl, col) in enumerate([(s_sph, "A. Sphere S⁵ in C³: squared radii\nfill the triangle (the 2-simplex)", OURS),
                                   (s_ball, "B. Ball in C³: squared radii fill\nthe solid corner simplex", OURS),
                                   (s_poly, "C. Polydisc (|zᵢ| ≤ 1 each): squared\nradii fill the cube", KNOWN)]):
    ax = fig.add_subplot(1, 3, i + 1, projection="3d")
    ax.scatter(*S.T, s=1.5, color=col, alpha=0.5)
    ax.set_xlabel("|z₁|²"); ax.set_ylabel("|z₂|²"); ax.set_zlabel("|z₃|²")
    ax.set_xlim(0, 1); ax.set_ylim(0, 1); ax.set_zlim(0, 1); ax.set_title(ttl); ax.set_box_aspect((1, 1, 1))
fig.suptitle("Figure 6. In polar coordinates the sphere and the simplex are one object (with an angle on every axis); "
             "the cube is the polydisc", y=1.0)
save(fig, "fig6-sphere-is-simplex-cube-is-polydisc.png")

# ---------------------------------------------------------------- fig7 polar split only in 1,2,4,8
def conj(a):
    if len(a) == 1: return a.copy()
    h = len(a) // 2; return np.concatenate([conj(a[:h]), -a[h:]])
def mul(a, b):
    if len(a) == 1: return a * b
    h = len(a) // 2; p, q, r, s = a[:h], a[h:], b[:h], b[h:]
    return np.concatenate([mul(p, r) - mul(conj(s), q), mul(s, p) + mul(q, conj(r))])
fig, ax = plt.subplots(figsize=(11, 4.4))
labels = ["1\nreal", "2\ncomplex", "4\nquaternion", "8\noctonion", "16\nsedenion"]
errs = []
for n in (1, 2, 4, 8, 16):
    e = 0
    for _ in range(500):
        u, v = rng.normal(size=n), rng.normal(size=n)
        e = max(e, abs(np.linalg.norm(mul(u, v)) - np.linalg.norm(u) * np.linalg.norm(v)) / (np.linalg.norm(u) * np.linalg.norm(v)))
    errs.append(max(e, 1e-17))
ax.bar(range(5), errs, color=[OURS] * 4 + [RED]); ax.set_yscale("log")
ax.set_xticks(range(5)); ax.set_xticklabels(labels)
ax.set_ylabel("worst | |uv| − |u||v| | (relative)")
for i, txt in enumerate(["sphere S⁰\ngroup", "circle S¹\ngroup,\ncommutative", "S³\ngroup,\norder matters",
                         "S⁷\nnot even\nassociative", "split\nbreaks"]):
    ax.text(i, errs[i] * 4, txt, ha="center", fontsize=9)
ax.set_ylim(1e-18, 10)
ax.set_title("Figure 7. 'Radius adds, sphere multiplies' needs a multiplication that keeps size.\n"
             "It exists only in 1, 2, 4 and 8 dimensions (Hurwitz 1898); in 3, 5, 6, 7 there is none to test")
save(fig, "fig7-polar-split-only-in-1-2-4-8.png")

# ---------------------------------------------------------------- fig8 curvature = 1/R^2
fig, axs = plt.subplots(1, 4, figsize=(19, 4.6), gridspec_kw={"wspace": 0.35})
def tri_sphere(ax):
    t = np.linspace(0, 2 * np.pi, 300); ax.plot(np.cos(t), np.sin(t), color=GREY)
    ax.fill([-0.55, 0.6, 0.05], [-0.35, -0.3, 0.6], color=OURS, alpha=0.25)
    for a, b in (((-0.55, -0.35), (0.6, -0.3)), ((0.6, -0.3), (0.05, 0.6)), ((0.05, 0.6), (-0.55, -0.35))):
        m = ((a[0] + b[0]) / 2 * 1.12, (a[1] + b[1]) / 2 * 1.12)
        xs = np.linspace(0, 1, 50)
        ax.plot((1 - xs) ** 2 * a[0] + 2 * xs * (1 - xs) * m[0] + xs ** 2 * b[0],
                (1 - xs) ** 2 * a[1] + 2 * xs * (1 - xs) * m[1] + xs ** 2 * b[1], color=OURS, lw=2)
    ax.set_title("sphere, radius R\nangles add to more than π\ncurvature +1/R²")
def tri_flat(ax):
    ax.fill([-0.55, 0.6, 0.05], [-0.35, -0.3, 0.6], color=GREY, alpha=0.25, ec="black")
    ax.set_title("plane, R = ∞\nangles add to π\ncurvature 0")
def tri_saddle(ax):
    t = np.linspace(0, 2 * np.pi, 300); ax.plot(np.cos(t), np.sin(t), color=GREY, ls=":")
    pts = [(-0.55, -0.35), (0.6, -0.3), (0.05, 0.6)]
    for a, b in ((pts[0], pts[1]), (pts[1], pts[2]), (pts[2], pts[0])):
        m = ((a[0] + b[0]) / 2 * 0.6, (a[1] + b[1]) / 2 * 0.6)
        xs = np.linspace(0, 1, 50)
        ax.plot((1 - xs) ** 2 * a[0] + 2 * xs * (1 - xs) * m[0] + xs ** 2 * b[0],
                (1 - xs) ** 2 * a[1] + 2 * xs * (1 - xs) * m[1] + xs ** 2 * b[1], color=KNOWN, lw=2)
    ax.set_title("'sphere of imaginary radius' iR (hyperbolic)\nangles add to less than π\ncurvature −1/R²")
for ax, f in zip(axs[:3], (tri_sphere, tri_flat, tri_saddle)):
    f(ax); ax.set_aspect("equal"); ax.axis("off"); ax.set_xlim(-1.1, 1.1); ax.set_ylim(-1.1, 1.1)
d = np.linspace(-1, 1, 5)
axs[3].plot(d, d, color=OURS, label="excess / area = curvature")
axs[3].plot([0.25, -1], [0.33013 / (4 * 0.33080), -0.39874 / 0.39872], "o", color=KNOWN, label="measured triangles")
axs[3].legend(fontsize=9)
axs[3].set_xlabel("curvature = 1/R²  (the curvature dial d, at k = 0)"); axs[3].set_ylabel("(angle sum − π) ÷ area")
axs[3].set_title("one line: angle excess = curvature × area\n(dots: sphere R = 2 and hyperbolic, P1, P2)")
fig.suptitle("Figure 8. Curvature is the sphere's 1/radius² (d at k = 0). Flat is R = ∞; hyperbolic is an imaginary radius (Lambert 1766)", y=1.04)
save(fig, "fig8-curvature-is-one-over-radius-squared.png")

# ---------------------------------------------------------------- fig9 swap rule on the sphere
fig = plt.figure(figsize=(14, 5))
ax = fig.add_subplot(1, 2, 1, projection="3d")
sphere_wire(ax, alpha=0.1)
th0 = 0.9
phis = np.linspace(0, 2 * np.pi, 13)
ax.plot(np.sin(th0) * np.cos(np.linspace(0, 2 * np.pi, 200)), np.sin(th0) * np.sin(np.linspace(0, 2 * np.pi, 200)),
        np.cos(th0) * np.ones(200), color=OURS, lw=2)
uu, vv = np.mgrid[0:2 * np.pi:30j, 0:th0:10j]
ax.plot_surface(np.cos(uu) * np.sin(vv), np.sin(uu) * np.sin(vv), np.cos(vv), color=OURS, alpha=0.15)
vec = np.array([np.cos(th0), 0, -np.sin(th0)])
N = 1200
pp = np.linspace(0, 2 * np.pi, N + 1)
step = N // 12
for i, p in enumerate(pp):
    nrm = np.array([np.sin(th0) * np.cos(p), np.sin(th0) * np.sin(p), np.cos(th0)])
    vec = vec - np.dot(vec, nrm) * nrm; vec /= np.linalg.norm(vec)
    if i % step == 0:
        ax.quiver(*nrm, *(0.3 * vec), color=RED, lw=1.5, arrow_length_ratio=0.25)
ax.set_title("A. Carry an arrow round the blue loop without twisting it:\n"
             "it comes back turned by (cap area) / R²  = d × area")
ax.set_axis_off(); ax.set_box_aspect((1, 1, 1))
ax = fig.add_subplot(1, 2, 2)
th = np.linspace(0.01, np.pi, 100)
solid = 2 * np.pi * (1 - np.cos(th))
ax.plot(solid, solid, color=OURS, label="arrow on the sphere: turn = area (d = 1)")
ax.plot(solid, 0.5 * solid, color=KNOWN, label="spin ½ (Hopf fibration): phase = ½ × area")
ax.plot([0.1 * 0.07], [0.1 * 0.07], "ko", label="S³ swap of turns a, b: leftover turn = a·b")
ax.set_xlabel("enclosed area on the unit sphere (solid angle)"); ax.set_ylabel("turn or phase")
ax.legend(fontsize=9)
ax.set_title("B. The swap rule is the sphere's curvature; spin's\nBerry phase is the same rule at half weight")
save(fig, "fig9-swap-rule-is-the-spheres-curvature.png")

# ---------------------------------------------------------------- fig10 celestial sphere
fig, axs = plt.subplots(1, 2, figsize=(13, 5.5))
for ax, (beta, ttl) in zip(axs, ((0.0, "A. At rest"), (0.8, "B. Moving at 0.8c toward the pole"))):
    D = math.sqrt((1 + beta) / (1 - beta))
    for th in np.radians([15, 30, 45, 60, 75, 90, 105, 120, 135, 150, 165]):
        r = math.tan(th / 2) / D
        t = np.linspace(0, 2 * np.pi, 200)
        ax.plot(r * np.cos(t), r * np.sin(t), color=OURS if th < np.pi / 2 else GREY, lw=1)
    for a in np.linspace(0, 2 * np.pi, 12, endpoint=False):
        ax.plot([0, 4 * np.cos(a)], [0, 4 * np.sin(a)], color=KNOWN, lw=0.6)
    # a star field circle: directions on a small circle, boosted
    c0 = np.array([0.5, 0.2, 0.0]); c0[2] = math.sqrt(1 - 0.29)
    e1 = np.cross(c0, [0, 0, 1.0]); e1 /= np.linalg.norm(e1); e2 = np.cross(c0, e1)
    xs, ys = [], []
    for t_ in np.linspace(0, 2 * np.pi, 100):
        dd = math.cos(0.3) * c0 + math.sin(0.3) * (math.cos(t_) * e1 + math.sin(t_) * e2)
        czp = (dd[2] + beta) / (1 + beta * dd[2]); thp = math.acos(czp)
        ph = math.atan2(dd[1], dd[0]); r = math.tan(thp / 2)
        xs.append(r * math.cos(ph)); ys.append(r * math.sin(ph))
    ax.plot(xs, ys, color=RED, lw=2)
    ax.set_aspect("equal"); ax.set_xlim(-3, 3); ax.set_ylim(-3, 3); ax.axis("off"); ax.set_title(ttl)
axs[1].text(-3, -3.4, "boost = every radius shrinks by the Doppler factor (here 3); turn = add to the angle.\n"
            "The red circle of stars stays a circle (Penrose–Terrell).", fontsize=10)
fig.suptitle("Figure 10. The sky in polar coordinates (stereographic): radius = boost, angle = rotation. "
             "Two boosts along one line multiply their Doppler factors", y=1.0)
save(fig, "fig10-celestial-sphere-boost-is-radius.png")

# ---------------------------------------------------------------- fig11 Kepler to S^3
from scipy import integrate
def rhs_k(t, y):
    r = y[:2]; return np.concatenate([y[2:], -r / np.linalg.norm(r) ** 3])
fig, axs = plt.subplots(1, 3, figsize=(15, 4.8))
for v0, col in ((1.2, OURS), (0.8, KNOWN), (1.05, RED)):
    sol = integrate.solve_ivp(rhs_k, (0, 40), [1, 0, 0, v0], rtol=1e-10, atol=1e-12, max_step=0.01)
    axs[0].plot(sol.y[0], sol.y[1], color=col)
    axs[1].plot(sol.y[2], sol.y[3], color=col)
    E = 0.5 * v0 ** 2 - 1; p0 = math.sqrt(-2 * E)
    P = sol.y[2:].T; p2 = np.sum(P ** 2, axis=1)
    W = np.c_[(p2 - p0 ** 2) / (p2 + p0 ** 2), 2 * p0 * P / (p2 + p0 ** 2)[:, None]]
    axs[2].plot(W[:, 0], W[:, 2], color=col)
axs[0].plot(0, 0, "k*", ms=12); axs[0].set_aspect("equal"); axs[0].set_title("A. Kepler orbits (ellipses) around the sun")
axs[1].plot(0, 0, "k+"); axs[1].set_aspect("equal"); axs[1].set_title("B. Their velocities trace circles\n(Hamilton 1846)")
t = np.linspace(0, 2 * np.pi, 200); axs[2].plot(np.cos(t), np.sin(t), color=GREY, ls=":")
axs[2].set_aspect("equal"); axs[2].set_title("C. Mapped onto the 3-sphere (Fock 1935):\neach is a great circle (shown edge-on slice)")
for a in axs: a.set_xticks([]); a.set_yticks([])
fig.suptitle("Figure 11. Newton's gravity, bound orbits only: every Kepler orbit becomes a great circle on S³ "
             "(Moser 1970). The same S³ gives hydrogen's n² levels", y=1.02)
save(fig, "fig11-kepler-orbits-are-great-circles.png")

# ---------------------------------------------------------------- fig12 sorting board
fig, ax = plt.subplots(figsize=(14, 6.6))
ax.axis("off")
rows = [
    ("i", "simpler to state", "same content as the k route on S¹; S³ adds a sphere of i's"),
    ("π", "simpler to state", "distance from 1 to −1 in every dimension; S³: 2π turn = −1 (spin)"),
    ("e", "not from the sphere", "right shape: the simplex in the cube's corner (1/n!); sphere gives e^π"),
    ("circle operation ∘ₖ", "comes out, limited", "radius adds, sphere multiplies: only dims 1, 2, 4, 8"),
    ("rungs", "does not come out", "the sphere has radius 0; rung 3 survives only as angle × number"),
    ("curvature", "simpler, unifies", "1/R² (d at k = 0): sphere, plane, imaginary sphere on one line"),
    ("swap rule", "comes out directly", "turn = curvature × area is the sphere's curvature; spin = half weight"),
    ("probability", "partly", "Born rule = cap area; union rule only for orthogonal hemispheres (cube)"),
    ("electromagnetism", "simpler", "4π in Coulomb = sphere's area; cube gives the same flux with 6 integrals"),
    ("special relativity", "simpler, unifies", "sky as sphere: boost = radius, turn = angle; Doppler factors multiply"),
    ("gravity", "shape yes, G no", "Kepler orbits = great circles on S³; de Sitter = S⁴; G not fixed"),
]
colors = {"simpler": OURS, "comes out": OURS, "unifies": OURS, "partly": KNOWN, "shape": KNOWN,
          "not": RED, "does not": RED}
y = 0.97
ax.text(0.0, y, "Result", weight="bold"); ax.text(0.17, y, "From the hypersphere?", weight="bold"); ax.text(0.42, y, "Why", weight="bold")
for name, verdict, why in rows:
    y -= 0.085
    col = OURS if any(w in verdict for w in ("simpler", "comes out directly", "unifies")) else \
        (RED if verdict.startswith(("not", "does not")) else KNOWN)
    ax.text(0.0, y, name, fontsize=11); ax.text(0.17, y, verdict, color=col, fontsize=11, weight="bold")
    ax.text(0.42, y, why, fontsize=10)
ax.set_title("Figure 12. Each result, started from the hypersphere (blue = comes out or simpler, orange = partly, red = no)")
save(fig, "fig12-sorting-board.png")
