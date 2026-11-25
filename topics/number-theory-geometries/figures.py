"""Figures for the number-theory geometries survey (DRAFT). Run: python3 figures.py  -> figures/*.png"""
import itertools, math
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d.art3d import Poly3DCollection, Line3DCollection
from scipy.spatial import ConvexHull
import importlib.util, os, sys, io, contextlib

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, 'figures'); os.makedirs(OUT, exist_ok=True)
# reuse helpers from checks.py without printing its output
spec = importlib.util.spec_from_file_location('checks', os.path.join(HERE, 'checks.py'))
C = importlib.util.module_from_spec(spec)
with contextlib.redirect_stdout(io.StringIO()):
    spec.loader.exec_module(C)

BLUE, ORANGE, GREEN, RED, PURPLE, GRAY = '#2a6fdb', '#e8871e', '#2e9e5b', '#d64545', '#7b5cc4', '#888888'
plt.rcParams.update({'font.size': 11, 'axes.titlesize': 12})

def cube_edges(n=3):
    V = list(itertools.product([0, 1], repeat=n))
    return [(a, b) for a, b in itertools.combinations(V, 2) if sum(abs(i - j) for i, j in zip(a, b)) == 1]

# ---------------------------------------------------------------- fig 0: the map
fig, ax = plt.subplots(figsize=(11, 6.2)); ax.axis('off'); ax.set_xlim(-0.1, 11.2); ax.set_ylim(-0.4, 6.3)
ax.add_patch(plt.Rectangle((4.2, 2.5), 2.6, 1.3, fc='#eef3fc', ec=BLUE, lw=2))
ax.text(5.5, 3.35, 'the hypercube', ha='center', fontsize=13, weight='bold')
ax.text(5.5, 2.85, 'corners = sets, paths = orderings\naggregate a ∘ₖ b = a+b+k·ab', ha='center', fontsize=9.5)
boxes = [
    (0.2, 5.0, 'Tropical geometry', 'bends: corner HEIGHTS\n→ cuts the cube into cells', GREEN, 'lead'),
    (3.9, 4.95, 'Matroid / permutohedron', 'bends: the PATH polytope\nk sets its size; sign of k = shape type', GREEN, 'lead'),
    (7.7, 5.0, 'Cluster algebras', 'bends: the cube ITSELF\ncube → associahedra', GREEN, 'lead'),
    (8.1, 2.5, 'p-adic / formal groups', 'bends: the OPERATION\nheight 1 → height 2', ORANGE, 'smaller lead'),
    (7.7, 0.2, 'Field with one element', 'bends: the poset\ncube = q→1 of subspaces', ORANGE, 'placement'),
    (3.9, 0.0, 'Arakelov / adelic', 'size only (heights)\nno bending', GRAY, 'not a lead'),
    (0.2, 0.2, 'Modular / hyperbolic', 'same as groups-lie page\n+ class numbers', GRAY, 'not a lead'),
    (0.0, 2.5, 'q-deformations', 'already on cube-weights page\n(q-orbit dials)', GRAY, 'already have'),
]
for x, y, title, body, col, tag in boxes:
    ax.add_patch(plt.Rectangle((x, y), 2.9, 1.15, fc='white', ec=col, lw=2))
    ax.text(x + 1.45, y + 0.85, title, ha='center', weight='bold', fontsize=10, color=col)
    ax.text(x + 1.45, y + 0.18, body, ha='center', fontsize=8.5)
    ax.text(x + 2.85, y - 0.05, tag, ha='right', va='top', fontsize=8, color=col, style='italic')
    cx, cy = x + 1.45, y + 0.575
    ax.annotate('', xy=(5.5 + (cx - 5.5)*0.28, 3.15 + (cy - 3.15)*0.28), xytext=(5.5 + (cx - 5.5)*0.72, 3.15 + (cy - 3.15)*0.72),
                arrowprops=dict(arrowstyle='->', color=col, lw=1.5))
ax.set_title('Fig 0. Which part of the hypercube each geometry bends (green = worth following)', loc='left')
fig.savefig(os.path.join(OUT, 'fig0-map.png'), dpi=140, bbox_inches='tight'); plt.close(fig)

# ---------------------------------------------------------------- fig 1: square lifted by heights
fig, axs = plt.subplots(1, 3, figsize=(12, 3.8))
sq = [(0, 0), (1, 0), (0, 1), (1, 1)]
cases = [('plain dial: heights |S|·v(k)\n(flat, no cut)', [0, 1, 1, 2]),
         ('heights lift corner {1,2}\n→ cut along the other diagonal', [0, 1, 1, 2.6]),
         ('heights sink corner {1,2}\n→ cut along the main diagonal', [0, 1, 1, 1.4])]
for ax, (title, h) in zip(axs, cases):
    cells = C.lower_cells(sq, h)
    for cell, col in zip(cells, [BLUE, ORANGE, GREEN]):
        pts = np.array([sq[i] for i in cell]); hull = ConvexHull(pts) if len(pts) > 2 else None
        poly = pts[hull.vertices]
        ax.fill(poly[:, 0], poly[:, 1], color=col, alpha=0.35, ec='k', lw=1.5)
    for (x, y), hh, lab in zip(sq, h, ['∅', '{1}', '{2}', '{1,2}']):
        ax.plot(x, y, 'ko'); ax.text(x + (0.04 if x == 0 else -0.04), y + (0.06 if y == 0 else -0.1), f'{lab}\nh={hh:g}',
                                     ha='left' if x == 0 else 'right', fontsize=9)
    ax.set_title(title, fontsize=10); ax.set_aspect('equal'); ax.set_xlim(-0.15, 1.15); ax.set_ylim(-0.15, 1.15); ax.axis('off')
fig.suptitle('Fig 1. Tropical geometry on the square: give each corner a height, keep the lower shell, read off the cut', x=0.02, y=1.05, ha='left')
fig.savefig(os.path.join(OUT, 'fig1-square-heights.png'), dpi=140, bbox_inches='tight'); plt.close(fig)

# ---------------------------------------------------------------- fig 2: 3-cube subdivisions
C3 = C.corners(3); rank = [sum(c) for c in C3]
rng = np.random.default_rng(7)
while True:
    hr = rng.random(8); cr = C.lower_cells(C3, hr)
    if len(cr) == 6: break
cases = [('plain dial (heights |S|·v(k)):\none cell, nothing cut', C.lower_cells(C3, [float(r) for r in rank])),
         ('q-binomial heights |S|(|S|−1)/2:\nthree level slabs', C.lower_cells(C3, [r*(r - 1)/2 for r in rank])),
         ('random heights: one of the\n74 triangulations (6 tetrahedra)', cr)]
fig = plt.figure(figsize=(13, 4.6))
cols = [BLUE, ORANGE, GREEN, RED, PURPLE, '#c9a227']
for idx, (title, cells) in enumerate(cases):
    ax = fig.add_subplot(1, 3, idx + 1, projection='3d')
    ax.add_collection3d(Line3DCollection([[a, b] for a, b in cube_edges()], colors='k', lw=1))
    for ci, cell in enumerate(sorted(cells, key=lambda c: sorted(c))):
        pts = np.array([C3[i] for i in cell], float)
        ctr = pts.mean(0); shrink = 0.0 if len(cells) == 1 else 0.12
        pts = ctr + (pts - ctr)*(1 - shrink) + (ctr - 0.5)*0.35*(len(cells) > 1)
        hull = ConvexHull(pts)
        faces = [pts[s] for s in hull.simplices]
        ax.add_collection3d(Poly3DCollection(faces, fc=cols[ci % 6], alpha=0.35, ec=cols[ci % 6], lw=0.3))
    for c in C3:
        ax.scatter(*c, color='k', s=12)
    ax.text(0, 0, -0.18, '∅', fontsize=9); ax.text(1, 1, 1.12, '{1,2,3}', fontsize=9)
    ax.set_title(title, fontsize=10); ax.set_axis_off(); ax.view_init(18, -58)
    ax.set_xlim(-0.2, 1.2); ax.set_ylim(-0.2, 1.2); ax.set_zlim(-0.2, 1.2)
fig.suptitle('Fig 2. The 3-cube cut by corner heights (pieces pulled apart to show them). Our dial sits at the uncut point.', x=0.02, y=1.05, ha='left')
fig.savefig(os.path.join(OUT, 'fig2-cube-subdivisions.png'), dpi=140, bbox_inches='tight'); plt.close(fig)

# ---------------------------------------------------------------- fig 3: hexagons by k
x0 = [0.3, 0.5, 0.8]
def to_plane(V):
    e1 = np.array([1, -1, 0])/math.sqrt(2); e2 = np.array([1, 1, -2])/math.sqrt(6)
    return np.c_[V @ e1, V @ e2]
fig, axs = plt.subplots(1, 2, figsize=(12.5, 5.4))
ax = axs[0]
for kk, col in [(-0.9, BLUE), (-0.5, '#6f9be8'), (0.0, 'k'), (0.5, '#f0ae62'), (0.9, ORANGE)]:
    V = C.path_vertices(x0, kk); P = to_plane(V - V.mean(0))
    if kk == 0:
        ax.plot(0, 0, 'ko', ms=7); ax.text(0.01, -0.03, 'k = 0: a single point', fontsize=9)
        continue
    hull = ConvexHull(P); poly = P[np.r_[hull.vertices, hull.vertices[0]]]
    ax.plot(poly[:, 0], poly[:, 1], '-', color=col, lw=2, label=f'k = {kk:+.1f}')
    ax.plot(P[:, 0], P[:, 1], 'o', color=col, ms=4)
V = C.path_vertices(x0, -0.9); P = to_plane(V - V.mean(0))
for perm, p in zip(itertools.permutations('123'), P):
    ax.text(p[0]*1.12, p[1]*1.12, ''.join(perm), fontsize=8.5, ha='center', va='center', color=BLUE)
ax.set_aspect('equal'); ax.legend(loc='lower right', fontsize=9)
ax.set_title('Each path (ordering, e.g. 213) gives a corner.\nk < 0 (blue, submodular) and k > 0 (orange, supermodular)\ngrow a hexagon out of a point', fontsize=10)
ax.set_xticks([]); ax.set_yticks([])
ax = axs[1]
ks = np.linspace(-1, 1, 201)
def second_diff(kk, x=(0.3, 0.5, 0.8)):
    f = lambda S: C.fold([x[i] for i in S], kk)
    return f((0, 1)) - f((0,)) - f((1,)) + f(())
ax.axhline(0, color='k', lw=0.8)
d = [second_diff(kk) for kk in ks]
ax.fill_between(ks, d, 0, where=np.array(d) < 0, color=BLUE, alpha=0.3, label='submodular: diminishing returns (union side, k = −1)')
ax.fill_between(ks, d, 0, where=np.array(d) > 0, color=ORANGE, alpha=0.3, label='supermodular: synergy (convex game)')
ax.plot(ks, d, 'k-')
ax.set_xlabel('dial k'); ax.set_ylabel('f({1,2}) − f({1}) − f({2}) + f(∅)  =  k·x₁·x₂')
ax.set_title('The sign of the dial is the type of the corner function', fontsize=10); ax.legend(fontsize=8.5, loc='upper left')
fig.suptitle('Fig 3. The paths through the 3-cube make a hexagon (a permutohedron) whose size is set by k', x=0.02, y=1.06, ha='left')
fig.savefig(os.path.join(OUT, 'fig3-path-hexagons.png'), dpi=140, bbox_inches='tight'); plt.close(fig)

# ---------------------------------------------------------------- fig 4: 3-D permutohedron n=4
x4 = [0.2, 0.35, 0.6, 0.75]
fig = plt.figure(figsize=(11, 4.8))
for idx, kk in enumerate([-0.7, 0.7]):
    V = C.path_vertices(x4, kk); ctr = V.mean(0)
    B = np.linalg.svd(V - ctr)[2][:3]; Y = (V - ctr) @ B.T
    hull = ConvexHull(Y)
    ax = fig.add_subplot(1, 2, idx + 1, projection='3d')
    ax.add_collection3d(Poly3DCollection([Y[s] for s in hull.simplices], fc=BLUE if kk < 0 else ORANGE, alpha=0.18, ec='none'))
    E = set()
    for s in hull.simplices:   # draw only edges between path vertices one swap apart
        pass
    perms = list(itertools.permutations(range(4)))
    for i, j in itertools.combinations(range(24), 2):
        a, b = perms[i], perms[j]
        diff = [t for t in range(4) if a[t] != b[t]]
        if len(diff) == 2 and diff[1] == diff[0] + 1:
            ax.plot(*zip(Y[i], Y[j]), color='k', lw=0.8)
    ax.scatter(*Y.T, color='k', s=10)
    ax.set_title(f'k = {kk:+.1f}: {"base polytope (submodular)" if kk < 0 else "core of the convex game"}\n24 orderings = 24 corners, swaps = edges', fontsize=10)
    ax.set_axis_off()
fig.suptitle('Fig 4. Four directions: the 24 paths through the tesseract make a 3-D permutohedron', x=0.02, y=1.05, ha='left')
fig.savefig(os.path.join(OUT, 'fig4-permutohedron-4d.png'), dpi=140, bbox_inches='tight'); plt.close(fig)

# ---------------------------------------------------------------- fig 5: cluster polygons
fig, axs = plt.subplots(1, 5, figsize=(15, 3.6))
data = [(4, 'no arrows (A1×A1)\nthe square = 2-cube', 'bc = 0', 'x′x = 2'),
        (5, 'A2', 'bc = 1', 'x′x = y + 1'), (6, 'B2', 'bc = 2', 'x′x = y² + 1'), (8, 'G2', 'bc = 3', 'x′x = y³ + 1'),
        (None, 'bc ≥ 4', 'infinite', 'never closes')]
for ax, (m, name, bc, rel) in zip(axs, data):
    ax.set_aspect('equal'); ax.axis('off')
    if m:
        th = np.linspace(0, 2*np.pi, m + 1) + np.pi/m + (np.pi/4 - np.pi/m if m == 4 else 0)
        ax.plot(np.cos(th), np.sin(th), '-o', color=GREEN if m > 4 else BLUE, lw=2)
        ax.text(0, 0, f'{m} clusters\nperiod {m}', ha='center', va='center', fontsize=9)
    else:
        th = np.linspace(0, 6*np.pi, 400); rr = 0.25 + th/(6*np.pi)*0.8
        ax.plot(rr*np.cos(th), rr*np.sin(th), color=GRAY)
    ax.set_xlim(-1.3, 1.3); ax.set_ylim(-1.3, 1.3)
    ax.set_title(f'{name}\n{bc}:  {rel}', fontsize=9.5)
fig.suptitle('Fig 5. Cluster algebras in rank 2: with no coupling the flips go round a square; turning the coupling bc up bends it', x=0.02, y=1.05, ha='left')
fig.savefig(os.path.join(OUT, 'fig5-cluster-polygons.png'), dpi=140, bbox_inches='tight'); plt.close(fig)

# ---------------------------------------------------------------- fig 6: cube vs associahedron (exchange graphs)
def eg_layout(B, seed):
    import random as R
    R.seed(seed)
    n = len(B)
    from fractions import Fraction as Fr
    X0 = [Fr(R.randint(2, 97), R.randint(2, 97)) for _ in range(n)]
    start = frozenset(X0); seen = {start: (B, X0)}; edges = set(); todo = [start]
    while todo:
        key = todo.pop(); Bc, Xc = seen[key]
        for kk in range(n):
            B2, X2 = C.mutate(Bc, Xc, kk); k2 = frozenset(X2); edges.add(frozenset([key, k2]))
            if k2 not in seen: seen[k2] = (B2, X2); todo.append(k2)
    nodes = list(seen); idx = {v: i for i, v in enumerate(nodes)}
    Aj = np.zeros((len(nodes), len(nodes)))
    for e in edges:
        a, b = tuple(e); Aj[idx[a], idx[b]] = Aj[idx[b], idx[a]] = 1
    Lp = np.diag(Aj.sum(1)) - Aj
    w, v = np.linalg.eigh(Lp)
    pos = v[:, 1:4]
    return pos, [(idx[tuple(e)[0]], idx[tuple(e)[1]]) for e in edges]
fig = plt.figure(figsize=(11, 4.8))
for idx, (B, title) in enumerate([([[0, 0, 0], [0, 0, 0], [0, 0, 0]], 'no arrows: 8 clusters, 12 flips\n= the 3-cube'),
                                   ([[0, 1, 0], [-1, 0, 1], [0, -1, 0]], 'arrows 1→2→3: 14 clusters, 21 flips\n= the A3 associahedron')]):
    pos, E = eg_layout(B, 3)
    ax = fig.add_subplot(1, 2, idx + 1, projection='3d')
    for a, b in E: ax.plot(*zip(pos[a], pos[b]), color=BLUE if idx == 0 else GREEN, lw=1.5)
    ax.scatter(*pos.T, color='k', s=14); ax.set_axis_off(); ax.set_title(title, fontsize=10)
fig.suptitle('Fig 6. Same three directions, with and without coupling (layout drawn from the flip graph itself)', x=0.02, y=1.05, ha='left')
fig.savefig(os.path.join(OUT, 'fig6-cube-to-associahedron.png'), dpi=140, bbox_inches='tight'); plt.close(fig)

# ---------------------------------------------------------------- fig 7: Boolean vs subspace lattice
fig, axs = plt.subplots(1, 2, figsize=(12, 4.8))
def hasse(ax, levels, edges_fn, title, col):
    pos = {}
    for r, lev in enumerate(levels):
        for i, v in enumerate(lev):
            pos[v] = ((i - (len(lev) - 1)/2)*(1.0 if len(lev) < 5 else 0.55), r)
    for a in pos:
        for b in pos:
            if edges_fn(a, b): ax.plot(*zip(pos[a], pos[b]), color=col, lw=0.8, alpha=0.7)
    for v, p in pos.items(): ax.plot(*p, 'o', color='k', ms=5)
    ax.set_title(title, fontsize=10); ax.axis('off')
B3 = [[s for s in itertools.combinations(range(3), r)] for r in range(4)]
hasse(axs[0], [[frozenset(s) for s in lev] for lev in B3], lambda a, b: a < b and len(b) == len(a) + 1,
      'q = 1: subsets of {1,2,3}: 1, 3, 3, 1 corners\n(the 3-cube; Möbius bottom-to-top = −1)', BLUE)
L3 = C.subspaces_F2(3)
levels = [sorted([s for s in L3 if len(s) == 2**d], key=lambda s: sorted(s)) for d in range(4)]
hasse(axs[1], levels, lambda a, b: a < b and len(b) == 2*len(a),
      'q = 2: subspaces of (F₂)³: 1, 7, 7, 1\n(Gaussian binomials; Möbius = −8)', ORANGE)
fig.suptitle('Fig 7. Field with one element: the cube is the q → 1 shadow of the lattice of subspaces (Tits 1957)', x=0.02, y=1.05, ha='left')
fig.savefig(os.path.join(OUT, 'fig7-cube-vs-subspaces.png'), dpi=140, bbox_inches='tight'); plt.close(fig)

# ---------------------------------------------------------------- fig 8: heights of formal groups + Newton polygon
fig, axs = plt.subplots(1, 2, figsize=(12.5, 4.6))
ax = axs[0]
p = 7
pts = [(j - 1, C.vp(math.comb(p, j), p)) for j in range(1, p + 1)]
ax.plot(*zip(*pts), 'ko'); ax.plot([0, p - 1], [1, 0], color=GREEN, lw=2)
for (xx, yy) in pts: ax.text(xx, yy + 0.05, f'{yy}', ha='center', fontsize=9)
ax.set_xlabel('power of x in ((1+x)⁷ − 1)/x'); ax.set_ylabel('how many times 7 divides the coefficient')
ax.set_title('Newton polygon at p = 7: one slope, −1/6\nso every 7-clock point of ∘₁ has 7-adic size 7^(−1/6)', fontsize=10)
ax.set_ylim(-0.3, 1.4)
ax = axs[1]
labels = ['∘₀  (plain +)', '∘₁  (multiplicative)', 'elliptic y²=x³−x\np = 5', 'elliptic y²=x³−x\np = 3, 7']
vals = [3.2, 1, 1, 2]; cols_ = [GRAY, BLUE, PURPLE, RED]
ax.bar(range(4), vals, color=cols_)
ax.set_xticks(range(4)); ax.set_xticklabels(labels, fontsize=9)
ax.set_yticks([1, 2, 3.2]); ax.set_yticklabels(['1', '2', '∞']); ax.set_ylim(0, 3.5)
ax.set_ylabel('height (z^(p^h) is the first term of [p](z) mod p)')
ax.set_title('Height: the one dial p-adics see and the real line cannot\n(every law on the landscape map has height 1 or ∞)', fontsize=10)
fig.suptitle('Fig 8. p-adic geometry: the size of the clock points, and the height of the law', x=0.02, y=1.05, ha='left')
fig.savefig(os.path.join(OUT, 'fig8-padic.png'), dpi=140, bbox_inches='tight'); plt.close(fig)

# ---------------------------------------------------------------- fig 9: product formula bars
from fractions import Fraction as Fr
import sympy as sp
u = 1 + Fr(3, 4)*Fr(10, 9)   # k = 3/4, x = 10/9 -> u = 1 + 5/6 = 11/6
places = ['real'] + [f'p = {p}' for p in sp.primefactors(u.numerator*u.denominator)]
vals = [math.log(abs(u))] + [math.log(C.abs_p(u, p)) for p in sp.primefactors(u.numerator*u.denominator)]
fig, ax = plt.subplots(figsize=(7.5, 3.8))
ax.bar(places, vals, color=[BLUE] + [ORANGE]*(len(vals) - 1)); ax.axhline(0, color='k', lw=0.8)
for i, v in enumerate(vals): ax.text(i, v + (0.05 if v > 0 else -0.15), f'{v:+.3f}', ha='center', fontsize=9)
ax.set_ylabel('ruler log|1 + kx| at that place')
ax.set_title(f'Fig 9. Arakelov: k = 3/4, x = 10/9, so 1 + kx = 11/6. The rulers at all places add to 0 (to rounding)\n(the product formula). True, tidy, and it bends nothing.', fontsize=10, loc='left')
fig.savefig(os.path.join(OUT, 'fig9-product-formula.png'), dpi=140, bbox_inches='tight'); plt.close(fig)
print('figures written to', OUT)
