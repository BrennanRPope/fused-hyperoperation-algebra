"""Figures for shape-set-operations/README.md. Run: python3 figures.py (writes figures/*.png)."""
import math
import itertools
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle, Circle, Polygon, Wedge
from mpl_toolkits.mplot3d.art3d import Poly3DCollection

OR, BL, GR, RD, PU = "#e07b39", "#3b75af", "#5a9e5a", "#c0392b", "#8e5ea2"
plt.rcParams.update({"font.size": 10, "axes.titlesize": 11})


def save(fig, name):
    fig.tight_layout()
    fig.savefig(f"figures/{name}.png", dpi=130, bbox_inches="tight")
    plt.close(fig)


# Fig 1: the operation is a union -------------------------------------------------
fig, ax = plt.subplots(1, 4, figsize=(15, 4))
# k = -1: box of side 1
a, b = 0.35, 0.45
A = ax[0]
A.add_patch(Rectangle((0, 0), 1, 1, fill=False, lw=1.5))
A.add_patch(Rectangle((0, 0), a, 1, color=OR, alpha=0.5, label="A x X (size a)"))
A.add_patch(Rectangle((0, 0), 1, b, color=BL, alpha=0.5, label="X x B (size b)"))
A.set_title("k = -1: ambient box of side 1\nunion = a + b - ab = a ∘(-1) b")
A.text(a / 2, 0.75, "a", ha="center"); A.text(0.75, b / 2, "b", ha="center"); A.text(a / 2, b / 2, "overlap\nab", ha="center", fontsize=8)
A.legend(loc="upper right", fontsize=7)
# k = -1/M, M large
A = ax[1]
M = 4
A.add_patch(Rectangle((0, 0), M, M, fill=False, lw=1.5))
A.add_patch(Rectangle((0, 0), a, M, color=OR, alpha=0.5))
A.add_patch(Rectangle((0, 0), M, b, color=BL, alpha=0.5))
A.set_title(f"k = -1/M (here M = {M}): bigger box\noverlap ab is a smaller share; k -> 0 = plain sum")
# k = +1: grown box minus unit square
A = ax[2]
A.add_patch(Rectangle((0, 0), 1 + a, 1 + b, color=GR, alpha=0.25))
A.add_patch(Rectangle((0, 0), 1, 1, color="white", ec="k", lw=1.5))
A.add_patch(Rectangle((1, 0), a, 1, color=OR, alpha=0.5)); A.add_patch(Rectangle((0, 1), 1, b, color=BL, alpha=0.5))
A.add_patch(Rectangle((1, 1), a, b, color=PU, alpha=0.6))
A.text(0.5, 0.5, "reference\nbox", ha="center", va="center")
A.text(1 + a / 2, 0.5, "a", ha="center"); A.text(0.5, 1 + b / 2, "b", ha="center"); A.text(1 + a / 2, 1 + b / 2, "ab", ha="center", fontsize=8)
A.set_title("k = +1, reading 1: grown box minus the reference\n= a + b + ab = a ∘(+1) b")
# k = +1 via Euler characteristic
A = ax[3]
A.set_xlim(-0.2, 3.2); A.set_ylim(-0.2, 3.2)
A.plot([0, 3], [0, 0], color="k", lw=2); A.plot([0, 3], [0, 0], "o", color="white", mec="k")
A.text(1.5, -0.15, "open line: Euler characteristic -1", ha="center", fontsize=8)
A.plot([0, 0], [0.3, 3], color="k", lw=2); A.plot([0, 0], [0.3, 3], "o", color="white", mec="k")
A.add_patch(Rectangle((0.6, 0.3), 0.8, 2.7, color=OR, alpha=0.5)); A.add_patch(Rectangle((0.05, 1.2), 2.95, 0.9, color=BL, alpha=0.5))
A.set_title("k = +1, reading 2: same union, counted by\nEuler characteristic; ambient line counts -1")
A.text(1.6, 2.6, "chi(union) / chi(line)\n= a ∘(+1) b", fontsize=9)
for A in ax:
    A.set_aspect("equal"); A.set_xticks([]); A.set_yticks([])
ax[0].set_xlim(-0.05, 1.05); ax[0].set_ylim(-0.05, 1.05); ax[1].set_xlim(-0.1, M + .1); ax[1].set_ylim(-0.1, M + .1)
ax[2].set_xlim(-0.05, 1.5); ax[2].set_ylim(-0.05, 1.5)
fig.suptitle("Fig 1. One operation, one picture: the size of a union of two independent slabs. The dial k is -1 / (size of the ambient).", y=1.02)
save(fig, "fig1-operation-is-a-union")

# Fig 2: the five set operations on two slabs -----------------------------------
fig, ax = plt.subplots(1, 5, figsize=(15, 3.4))
a, b = 0.4, 0.5
xs, ys = np.meshgrid(np.linspace(0, 1, 400), np.linspace(0, 1, 400))
inA, inB = xs < a, ys < b
ops = [("union  A ∪ B", inA | inB, "a ∘(-1) b"), ("intersection  A ∩ B", inA & inB, "ab (the product)"),
       ("symmetric difference", inA ^ inB, "a ∘(-2) b  (XOR)"), ("difference  A \\ B", inA & ~inB, "a(1 - b)"),
       ("complement of A ∪ B", ~(inA | inB), "(1-a)(1-b)  (De Morgan)")]
for A, (t, m, f) in zip(ax, ops):
    A.imshow(m, origin="lower", extent=(0, 1, 0, 1), cmap=matplotlib.colors.ListedColormap(["white", PU]))
    A.axvline(a, color=OR, lw=1.5); A.axhline(b, color=BL, lw=1.5)
    A.set_title(f"{t}\nsize = {f}", fontsize=9); A.set_xticks([]); A.set_yticks([])
fig.suptitle("Fig 2. Every Boolean operation on two slabs lands on the dial family or on the product (unit box; a = 0.4, b = 0.5)", y=1.04)
save(fig, "fig2-five-set-operations")

# Fig 3: pair dial ---------------------------------------------------------------
def lens(r, d):
    if d >= 2 * r:
        return 0.0
    return 2 * r * r * math.acos(d / (2 * r)) - d / 2 * math.sqrt(4 * r * r - d * d)

fig, ax = plt.subplots(1, 3, figsize=(15, 4.2))
A = ax[0]
for i, d in enumerate([0.0, 0.6, 1.3, 2.1]):
    cx = i * 4.2
    A.add_patch(Circle((cx, 0), 1, color=OR, alpha=0.45)); A.add_patch(Circle((cx + d, 0), 1, color=BL, alpha=0.45))
    kk = -lens(1, d) / math.pi ** 2
    A.text(cx + d / 2, -1.5, f"d = {d}\nk*area = {kk*math.pi:.2f}", ha="center", fontsize=8)
A.set_xlim(-1.3, 16.5); A.set_ylim(-2.1, 1.3); A.set_aspect("equal"); A.axis("off")
A.set_title("Two disks: union = a ∘(k) b with k = -overlap/(a b).\nOn top of each other: k*area = -1, the wall")
A = ax[1]
ds = np.linspace(0, 2.2, 200)
A.plot(ds, [-lens(1, d) / math.pi for d in ds], color=PU, lw=2)
A.axhline(-1, color="gray", ls=":"); A.axhline(0, color="gray", ls=":")
A.text(0.05, -0.95, "wall: shape joined with itself", fontsize=8); A.text(1.3, 0.04, "apart: plain sum (k = 0)", fontsize=8)
A.set_xlabel("distance between centres"); A.set_ylabel("k x area")
A.set_title("Pair dial of two unit disks vs distance")
A = ax[2]
th = np.linspace(0, math.pi, 200)
A.plot(th, -2 * (1 - th / math.pi), color=OR, lw=2, label="two hemispheres (sets): -2(1 - theta/pi)")
A.plot(th, -(1 + np.cos(th)), color=BL, lw=2, ls="--", label="quantum spin, same pair: -(1 + cos theta)")
A.plot([math.pi / 2], [-1], "ko"); A.text(math.pi / 2 + 0.05, -1.08, "right angle: k = -1, independent", fontsize=8)
A.set_xlabel("angle between the two hemispheres (theta)"); A.set_ylabel("pair dial k (ambient = whole sphere)")
A.legend(fontsize=8, loc="lower right")
A.set_title("Sets give a straight line; quantum gives a cosine\n(Bell's 1964 picture, read as a dial)")
fig.suptitle("Fig 3. Correlated shapes: every pair gets its own dial, k(A, B) = -size(A ∩ B) / (size A x size B)", y=1.03)
save(fig, "fig3-pair-dial")

# Fig 4: three disks, defect --------------------------------------------------------
fig, ax = plt.subplots(1, 2, figsize=(11, 4.6))
A = ax[0]
cs = [(0, 0.5), (-0.45, -0.3), (0.45, -0.3)]
for c, col in zip(cs, [OR, BL, GR]):
    A.add_patch(Circle(c, 1, color=col, alpha=0.35))
A.text(0, -0.05, "triple\noverlap", ha="center", fontsize=9, weight="bold")
A.set_xlim(-1.6, 1.6); A.set_ylim(-1.5, 1.7); A.set_aspect("equal"); A.axis("off")
A.set_title("Three disks: the centre is where associativity\ncan fail (Yeung's I-measure centre)")
A = ax[1]
A.add_patch(Rectangle((0, 0), 1, 1, fill=False))
for i, (l, col) in enumerate(zip([0.3, 0.5, 0.7], [OR, BL, GR])):
    if i == 0: A.add_patch(Rectangle((0, 0), l, 1, color=col, alpha=0.35))
    if i == 1: A.add_patch(Rectangle((0, 0), 1, l, color=col, alpha=0.35))
    if i == 2: A.add_patch(Rectangle((0, 0), 1, 1, color=col, alpha=0.12))
A.text(0.5, 1.05, "three independent slabs (third axis into the page)", ha="center", fontsize=8)
A.set_xlim(-0.05, 1.05); A.set_ylim(-0.05, 1.15); A.set_aspect("equal"); A.axis("off")
A.set_title("Three independent slabs: centre = k^2 abc exactly,\nso the union is a o b o c (defect 0)")
fig.suptitle("Fig 4. Two shapes always fit the law with some dial. Three shapes fit one dial only if the triple overlap is k^2 abc. "
             "Disks: defect 0.24", y=1.02, fontsize=10)
save(fig, "fig4-three-shapes-defect")

# Fig 5: e on the l_p exponent ------------------------------------------------------
from scipy.special import gamma
fig, ax = plt.subplots(1, 2, figsize=(12, 4.6))
A = ax[0]
A.add_patch(Rectangle((0, 0), 1, 1, fill=False, lw=1.5))
t = np.linspace(0, 1, 300)
for p, col in [(1, OR), (2, BL), (3, GR), (6, PU)]:
    A.plot(t, (1 - t ** p) ** (1 / p), color=col, lw=2, label=f"p = {p}")
A.fill_between(t, 0, 1 - t, color=OR, alpha=0.15)
A.set_aspect("equal"); A.legend(fontsize=8); A.set_title("Unit square cut by the corner of the l_p ball\np = 1 corner triangle (area 1/2) ... p = 2 quarter disc (pi/4)")
A = ax[1]
ps = np.linspace(1, 6, 200)
def Ep(p):
    n = np.arange(0, 400)
    from scipy.special import gammaln
    return np.exp(n * gammaln(1 + 1 / p) - gammaln(1 + n / p)).sum()
A.plot(ps, [Ep(p) for p in ps], color=PU, lw=2)
for p, lab in [(1, "e = 2.718"), (2, "e^(pi/4)(1 + erf(sqrt(pi)/2)) = 3.926"), (3, "5.161")]:
    A.plot([p], [Ep(p)], "ko"); A.text(p + 0.08, Ep(p) - 0.25, lab, fontsize=8)
A.set_xlabel("l_p exponent p"); A.set_ylabel("sum over dimensions of the cut-off volume")
A.set_title("e generalizes along p to a Mittag-Leffler value:\nE_(1/p)(Gamma(1 + 1/p)) = average draws until sum x^p > 1")
fig.suptitle("Fig 5. The e result is an intersection: cube ∩ (corner of the l_p ball). The l_p exponent moves it continuously.", y=1.02)
save(fig, "fig5-e-on-the-lp-exponent")

# Fig 6: e on the k dial ---------------------------------------------------------------
fig, ax = plt.subplots(1, 2, figsize=(12, 4.6))
A = ax[0]
A.add_patch(Rectangle((0, 0), 1, 1, fill=False, lw=1.5))
for kv, col in [(-0.7, RD), (-0.4, OR), (1e-9, "k"), (1, BL), (4, GR)]:
    y = ((1 + kv) / (1 + kv * t) - 1) / kv
    y = np.clip(y, 0, 1)
    A.plot(t, y, color=col, lw=2, label=f"k = {round(kv, 2)}")
A.set_aspect("equal"); A.legend(fontsize=8)
A.set_title("Unit square cut by x ∘(k) y <= 1\nk = 0 straight (triangle); k > 0 sags; k < 0 bulges")
A = ax[1]
ks = np.concatenate([np.linspace(-0.9, -0.01, 120), np.linspace(0.01, 12, 300)])
E = (1 + ks) ** (1 / ks) + ks / (1 + ks)
A.plot(ks, E, color=PU, lw=2)
for kv, lab in [(-0.5, "k = -1/2: exactly 3"), (0, "k = 0: e"), (1, "k = 1: exactly 5/2"), (math.e - 1, "k = e-1: Vandervelde 2008")]:
    val = math.e if kv == 0 else (1 + kv) ** (1 / kv) + kv / (1 + kv)
    A.plot([kv], [val], "ko"); A.text(kv + 0.25, val + 0.05, lab, fontsize=8)
A.axhline(2, color="gray", ls=":"); A.text(8, 2.05, "k -> infinity: 2", fontsize=8)
A.set_ylim(1.8, 4.2); A.set_xlabel("dial k"); A.set_ylabel("average draws until x1 ∘k x2 ∘k ... passes 1")
A.set_title("E(k) = (1+k)^(1/k) + k/(1+k)\nat k = 1/m: the staircase (1+1/m)^m plus 1/(m+1)")
fig.suptitle("Fig 6. The same intersection along the homothety dial k: cube ∩ (the ∘ₖ ball of radius 1)", y=1.02)
save(fig, "fig6-e-on-the-k-dial")

# Fig 7: Steiner ------------------------------------------------------------------------
fig, ax = plt.subplots(1, 2, figsize=(11, 4.6))
A = ax[0]
tt = 0.35
A.add_patch(Rectangle((0, 0), 1, 1, color=OR, alpha=0.6))
for (x0, y0, w, h) in [(0, -tt, 1, tt), (0, 1, 1, tt), (-tt, 0, tt, 1), (1, 0, tt, 1)]:
    A.add_patch(Rectangle((x0, y0), w, h, color=BL, alpha=0.5))
for (cx, cy, th1) in [(0, 0, 180), (1, 0, 270), (1, 1, 0), (0, 1, 90)]:
    A.add_patch(Wedge((cx, cy), tt, th1, th1 + 90, color=GR, alpha=0.6))
A.text(0.5, 0.5, "layer 0:\nC(2,2) = 1 square", ha="center", fontsize=8)
A.text(0.5, 1.15, "layer 1: C(2,1) = ... 4 edge strips, t x perimeter", ha="center", fontsize=7)
A.text(1.2, 1.25, "corners: one disc", fontsize=7)
A.set_xlim(-0.5, 1.6); A.set_ylim(-0.5, 1.5); A.set_aspect("equal"); A.axis("off")
A.set_title("Square + disc of radius t: square + strips + one whole disc")
A = ax[1]
n = 4
labels = [f"j = {j}" for j in range(n + 1)]
A.bar(range(n + 1), [math.comb(n, j) for j in range(n + 1)], color=OR, alpha=0.7)
A.set_xticks(range(n + 1)); A.set_xticklabels(labels)
A.set_title("The 4-cube's intrinsic volumes V_j = C(4, j):\nthe same Pascal row as the o-corners of the 4-cube")
fig.suptitle("Fig 7. Minkowski sum (a set operation that adds shapes): vol(cube + tB) = sum_j C(n,j) kappa_(n-j) t^(n-j)  (Steiner 1840)", y=1.02)
save(fig, "fig7-steiner-cube-plus-ball")

# Fig 8: cube ∩ cross-polytope ---------------------------------------------------------
fig = plt.figure(figsize=(15, 4.6))
A = fig.add_subplot(1, 3, 1)
for s, col in [(1, OR), (1.5, BL), (2, GR)]:
    A.add_patch(Polygon([(s, 0), (0, s), (-s, 0), (0, -s)], fill=False, ec=col, lw=2, label=f"|x| + |y| <= {s}"))
A.add_patch(Rectangle((-1, -1), 2, 2, fill=False, lw=2, ec="k", label="square"))
A.set_xlim(-2.1, 2.1); A.set_ylim(-2.1, 2.1); A.set_aspect("equal"); A.legend(fontsize=7)
A.set_title("Square ∩ diamond, three sizes\ns = 1: four corner triangles (e's piece)")
def poly3d(ax, verts, col):
    from scipy.spatial import ConvexHull
    h = ConvexHull(verts)
    ax.add_collection3d(Poly3DCollection([verts[s] for s in h.simplices], alpha=0.35, facecolor=col, edgecolor="k", lw=0.3))
    ax.set_xlim(-1.2, 1.2); ax.set_ylim(-1.2, 1.2); ax.set_zlim(-1.2, 1.2); ax.set_box_aspect((1, 1, 1))
A = fig.add_subplot(1, 3, 2, projection="3d")
V = np.array([v for v in itertools.product([-1, 0, 1], repeat=3) if sum(map(abs, v)) == 2], float)
poly3d(A, V, BL); A.set_title("3D: cube ∩ octahedron (s = 2)\n= cuboctahedron, 12 corners")
A = fig.add_subplot(1, 3, 3)
n = np.arange(1, 9)
from math import factorial
def ih(n, s):
    return sum((-1) ** j * math.comb(n, j) * max(s - j, 0) ** n for j in range(n + 1)) / factorial(n)
A.plot(n, [ih(k, 1) for k in n], "o-", color=OR, label="s = 1: 1/n!  (sum = e)")
A.plot(n, [ih(k, 2) for k in n], "o-", color=GR, label="s = 2  (n = 4: the 24-cell, volume 16 x 1/2 = 8)")
A.plot(n, [ih(k, k / 2) for k in n], "o-", color=BL, label="s = n/2: exactly 1/2")
A.set_yscale("log"); A.set_xlabel("dimension n"); A.set_ylabel("volume of cube corner ∩ {sum <= s}")
A.legend(fontsize=7); A.set_title("Volumes are inclusion-exclusion over the cube's corners\n(Irwin-Hall); slabs give Eulerian numbers")
fig.suptitle("Fig 8. One intersection family, cube ∩ cross-polytope, passes through e, the cuboctahedron and the 24-cell", y=1.02)
save(fig, "fig8-cube-cap-cross-polytope")

# Fig 9: CHSH -----------------------------------------------------------------------------
fig, ax = plt.subplots(1, 2, figsize=(12, 4.4))
A = ax[0]
A.bar([0, 1, 2], [2, 2 * math.sqrt(2), 4], color=[OR, BL, PU], alpha=0.7)
A.set_xticks([0, 1, 2]); A.set_xticklabels(["sets with ordinary sizes\n(cube corners)", "quantum\n(Tsirelson)", "no-signalling\nmaximum (PR box)"])
A.set_ylabel("best CHSH value"); A.set_title("How far correlations can go")
for i, v in enumerate([2, 2 * math.sqrt(2), 4]):
    A.text(i, v + 0.05, f"{v:.3f}", ha="center")
A = ax[1]
S = np.linspace(2, 4, 50)
A.plot(S, (S - 2) / 4, color=RD, lw=2)
A.plot([2 * math.sqrt(2)], [(math.sqrt(2) - 1) / 2], "ko")
A.text(2.85, 0.19, "Tsirelson: least negative size\n0.2071 = (sqrt 2 - 1)/2 (linear program)", fontsize=8)
A.set_xlabel("CHSH value"); A.set_ylabel("negative size needed (total)")
A.set_title("Past 2, the sets need negative sizes\n(like Euler characteristic, where an open line counts -1)")
fig.suptitle("Fig 9. The limit of 'union at the core': ordinary sizes give only classical correlations. Signed sizes reach quantum ones.", y=1.02)
save(fig, "fig9-sets-vs-quantum")
print("ok")


# Fig 10: all 16 operations, granular -------------------------------------------------------
names = [((0, 0, 0, 0), "false", "0"), ((0, 0, 0, 1), "AND", "-k ab"), ((0, 0, 1, 0), "A minus B", "a + k ab"), ((0, 0, 1, 1), "first", "a"),
         ((0, 1, 0, 0), "B minus A", "b + k ab"), ((0, 1, 0, 1), "second", "b"), ((0, 1, 1, 0), "XOR", "a ∘(2k) b"), ((0, 1, 1, 1), "OR", "a ∘(k) b"),
         ((1, 0, 0, 0), "NOR", "-1/k - a∘(k)b"), ((1, 0, 0, 1), "XNOR", "mirror of XOR"), ((1, 0, 1, 0), "not second", "-1/k - b"),
         ((1, 0, 1, 1), "B implies A", "-1/k - b - k ab"), ((1, 1, 0, 0), "not first", "-1/k - a"), ((1, 1, 0, 1), "A implies B", "-1/k - a - k ab"),
         ((1, 1, 1, 0), "NAND", "-1/k + k ab"), ((1, 1, 1, 1), "true", "-1/k (the wall)")]
monoid = {"OR", "XOR", "AND", "XNOR"}
fig, axs = plt.subplots(4, 4, figsize=(12, 12.5))
a, b = 0.4, 0.5
for A, (tt, nm, fm) in zip(axs.flat, names):
    cells = [(0, 0, a, b, tt[3]), (a, 0, 1 - a, b, tt[1]), (0, b, a, 1 - b, tt[2]), (a, b, 1 - a, 1 - b, tt[0])]
    # (x0, y0, w, h, on): A = x < a, B = y < b; tt = f(00), f(01), f(10), f(11) with first = A
    for x0, y0, w, h, on in cells:
        A.add_patch(Rectangle((x0, y0), w, h, color=PU if on else "white", alpha=0.75 if on else 1, ec="0.6", lw=0.5))
    A.set_xlim(0, 1); A.set_ylim(0, 1); A.set_aspect("equal"); A.set_xticks([]); A.set_yticks([])
    col = "#2e7d32" if nm in monoid else "k"
    A.set_title(f"{nm}\nsize = {fm}", fontsize=9, color=col, weight="bold" if nm in monoid else "normal")
    for sp_ in A.spines.values():
        sp_.set_edgecolor(col); sp_.set_linewidth(2.5 if nm in monoid else 0.8)
fig.suptitle("Fig 10. All 16 operations on two independent sets (A = left strip, B = bottom strip), sizes in units where the ambient is -1/k.\n"
             "Green: the only four that are commutative and associative: OR = ∘k, XOR = ∘2k, and their mirrors AND, XNOR", y=1.0)
save(fig, "fig10-sixteen-operations")

# Fig 11: what each op does to the cube's corners ---------------------------------------------
fig, axs = plt.subplots(1, 4, figsize=(16, 4.4))
proj = {}
for v in itertools.product([0, 1], repeat=3):
    proj[v] = (v[0] + 0.45 * v[2], v[1] + 0.35 * v[2])
edges = [(u, v) for u in proj for v in proj if sum(abs(x - y) for x, y in zip(u, v)) == 1 and u < v]
def weight(op, v):
    s = sum(v)
    if s == 0:
        return None
    if op == "OR":
        return f"k^{s-1}" if s > 1 else "1"
    if op == "XOR":
        return f"(2k)^{s-1}" if s > 1 else "1"
    if op == "AND":
        return "k^2" if s == 3 else None
    if op == "diff":
        return None if v[0] == 0 else ("1" if s == 1 else f"k^{s-1}")
titles = {"OR": "union a∘b∘c: every corner,\nweight k^(|S|-1)  (inclusion-exclusion)",
          "XOR": "XOR: every corner,\nweight (2k)^(|S|-1)  (piling-up lemma)",
          "AND": "intersection: the top corner only\n(abc, times k^2 in these units)",
          "diff": "a minus b minus c: only the face\nthat contains a  (a times survival of b, c)"}
for A, op in zip(axs, ["OR", "XOR", "AND", "diff"]):
    for u, v in edges:
        A.plot(*zip(proj[u], proj[v]), color="0.7", lw=1)
    for v, (x, y) in proj.items():
        w = weight(op, v)
        A.plot(x, y, "o", ms=14, color=PU if w else "white", mec="k")
        lab = "".join(ch for ch, on in zip("abc", v) if on) or "∅"
        A.text(x, y - 0.17, lab, ha="center", fontsize=8)
        if w:
            A.text(x + 0.1, y + 0.1, w, fontsize=8, color=PU)
    A.set_xlim(-0.4, 1.8); A.set_ylim(-0.4, 1.6); A.set_aspect("equal"); A.axis("off")
    A.set_title(titles[op], fontsize=9)
fig.suptitle("Fig 11. On the cube's corners: union and XOR use every corner (one dial, k or 2k); intersection keeps one corner; "
             "difference keeps one face", y=1.03)
save(fig, "fig11-ops-on-cube-corners")

# Fig 12: the tetrahedron of Venn cells; distributivity with the shared-ambient dial ---------------
fig = plt.figure(figsize=(13, 5))
A = fig.add_subplot(1, 2, 1, projection="3d")
T = np.array([[1, 1, 1], [1, -1, -1], [-1, 1, -1], [-1, -1, 1]], float)   # corners: neither, B only, A only, both
labels = ["neither", "B only", "A only", "both"]
for i, j in itertools.combinations(range(4), 2):
    A.plot(*zip(T[i], T[j]), color="0.5", lw=1)
for i in range(4):
    A.text(*(T[i] * 1.12), labels[i], fontsize=9)
uu, vv = np.meshgrid(np.linspace(0, 1, 25), np.linspace(0, 1, 25))
W = [(1 - uu) * (1 - vv), (1 - uu) * vv, uu * (1 - vv), uu * vv]
X = sum(W[i][..., None] * T[i] for i in range(4))
A.plot_surface(X[..., 0], X[..., 1], X[..., 2], color=BL, alpha=0.5, lw=0)
A.set_title("Two sets = a point in the tetrahedron of the four Venn cells.\nIndependent pairs = the blue saddle (odds ratio 1, pair dial -1).\n"
            "Above it: overlap more (k < -1); below: overlap less", fontsize=9)
A.set_axis_off()
A = fig.add_subplot(1, 2, 2)
A.add_patch(Rectangle((0, 0), 1, 1, fill=False))
pa_, pb_, pc_ = 0.55, 0.45, 0.6
A.add_patch(Rectangle((0, 0), pa_, 1, color=OR, alpha=0.25))
A.add_patch(Rectangle((0, 0), pa_, pb_, color=PU, alpha=0.55))
A.add_patch(Rectangle((0, 0), pa_ * pc_, 1, color=GR, alpha=0.0, hatch="//", ec=GR))
A.text(pa_ / 2, 0.93, "A (size a) = the new ambient", ha="center", fontsize=8)
A.text(pa_ / 2, pb_ / 2, "A ∩ B", ha="center", fontsize=9)
A.text(0.78, 0.5, "outside A:\nnever counted", ha="center", fontsize=8)
A.set_xlim(0, 1); A.set_ylim(0, 1); A.set_aspect("equal"); A.set_xticks([]); A.set_yticks([])
A.set_title("a ∩ (b ∘k c) = (a∩b) ∘(-1/a) (a∩c): the two pieces are independent\ninside A, so their dial is -1/(size of A). "
            "Restricting = turning the dial", fontsize=9)
fig.suptitle("Fig 12. Correlated sets and the distributive law: the dial is -1/(size of whatever ambient the two sets are independent in)", y=1.02)
save(fig, "fig12-tetrahedron-and-distributivity")
print("ok2")
