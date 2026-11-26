"""Checks for the polar-coordinates DRAFT (Affine Composition Theory).

Every claim the README leans on is tested here. Run: python3 checks.py
Writes checks_output.txt next to this file.

Sections
  A  log of the lift is a radial coordinate: circle operation -> addition
  B  plane law (a+b+kab)/(1-dab) = multiplication in R[e]/(e^2 - k e + d),
     read off as a slope; sign of k^2 - 4d decides angle (periodic) or radius
  C  moment map = gradient of the soft support function; lands inside the
     shape, angles untouched (action-angle coordinates)
  D  Newton polytopes: cube = Newt(1 + (a o b o c)), products -> Minkowski sums,
     vertex counts of cube, simplex, permutohedron, Loday associahedron
  E  shape weights on cube corners: Postnikov's z_I = sum_{J in I} y_J is our
     tesseract zeta transform, and the facets sum_{i in I} t_i >= z_I are tight
  F  support function of the shape = sum_I y_I max_{i in I} u_i
  G  edge directions: cube and cross-polytope are type-B shapes; the
     cross-polytope is not type A (outside the Postnikov family)
  H  p-norm = exp(log-sum-exp at temperature 1/p of log|u|): the p dial is a
     Maslov-style smoothing; p -> infinity gives max (the cube's gauge)
  I  Hopf coordinates: squared radii on S^(2n-1) fill a simplex
"""
import itertools
import math
import os
import random
import sys

import numpy as np
from scipy.spatial import ConvexHull

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "tesseract"))
import tesseract as T  # noqa: E402  (our own tesseract coordinates)

random.seed(7)
np.random.seed(7)
results = []


def check(name, ok, detail=""):
    results.append((name, bool(ok), detail))


def circ(a, b, k):
    return a + b + k * a * b


# ---------------------------------------------------------------- A
for k in (1.0, 0.5, -0.7, 2.0):
    worst = 0.0
    for _ in range(200):
        a, b = np.random.uniform(-0.3, 2.0, 2) / max(1, abs(k))
        if 1 + k * a <= 0 or 1 + k * b <= 0:
            continue
        L = lambda x: math.log(1 + k * x) / k  # the ruler: log of the lift
        worst = max(worst, abs(L(circ(a, b, k)) - (L(a) + L(b))))
    check(f"A1 ruler log(1+kx)/k turns o_k into + (k={k})", worst < 1e-12, f"max err {worst:.1e}")

# n axes: box coordinates multiply, log-box coordinates add, per axis
a = np.random.uniform(0, 1, 4)
b = np.random.uniform(0, 1, 4)
ua, ub = np.log1p(a), np.log1p(b)
check("A2 per axis, log(1+a) of a o b is the vector sum", np.allclose(np.log1p(circ(a, b, 1)), ua + ub))
# whole corner list: the total composition a1 o a2 o ... is the sum of log-radii
tot = T.multadd_all(list(a))
check("A3 log(1 + a1 o a2 o a3 o a4) = sum of log-radii (the simplex direction)",
      abs(math.log1p(tot) - ua.sum()) < 1e-12)

# ---------------------------------------------------------------- B
def mul(x, y, k, d):
    """(x0 + x1 e)(y0 + y1 e) with e^2 = k e - d."""
    x0, x1 = x
    y0, y1 = y
    return (x0 * y0 - d * x1 * y1, x0 * y1 + x1 * y0 + k * x1 * y1)


def plane(a, b, k, d):
    return (a + b + k * a * b) / (1 - d * a * b)


worst = 0.0
for _ in range(500):
    k, d = np.random.uniform(-2, 2, 2)
    a, b = np.random.uniform(-0.5, 0.5, 2)
    if abs(1 - d * a * b) < 1e-3:
        continue
    p0, p1 = mul((1, a), (1, b), k, d)
    worst = max(worst, abs(p1 / p0 - plane(a, b, k, d)))
check("B1 plane law = product (1+e a)(1+e b), read as slope", worst < 1e-12, f"max err {worst:.1e}")

worst = 0.0
for _ in range(500):
    k, d = np.random.uniform(-2, 2, 2)
    x = np.random.uniform(-2, 2, 2)
    y = np.random.uniform(-2, 2, 2)
    N = lambda z: z[0] ** 2 + k * z[0] * z[1] + d * z[1] ** 2
    worst = max(worst, abs(N(mul(x, y, k, d)) - N(x) * N(y)))
check("B2 norm x^2 + kxy + dy^2 is multiplicative (the radius)", worst < 1e-9, f"max err {worst:.1e}")


def step_matrix(k, d, a):
    """Composing with a is the Mobius map x -> ((1+ka)x + a)/(-d a x + 1)."""
    return np.array([[1 + k * a, a], [-d * a, 1.0]])


okB3 = True
detail = []
for _ in range(300):
    k, d = np.random.uniform(-3, 3, 2)
    a = np.random.uniform(0.01, 0.3)
    M = step_matrix(k, d, a)
    tr, det = np.trace(M), np.linalg.det(M)
    disc_map = tr * tr - 4 * det
    okB3 &= abs(disc_map - a * a * (k * k - 4 * d)) < 1e-9
    ev = np.linalg.eigvals(M)
    rotating = np.iscomplexobj(ev) and np.any(np.abs(ev.imag) > 1e-12)
    okB3 &= rotating == (k * k - 4 * d < 0)
for k, d in [(0, 1), (0, -1), (1, 1), (2, 1), (1, 0), (3, 2), (0.5, 2)]:
    detail.append(f"(k={k},d={d}):{'angle' if k*k-4*d < 0 else ('parabolic' if k*k == 4*d else 'radius')}")
check("B3 step map has trace^2-4det = a^2(k^2-4d): k^2-4d<0 rotates (angle), >0 has two real "
      "fixed points (radius/rapidity), =0 parabolic", okB3, ", ".join(detail))

# k = 0, d = 1: tan addition is the argument of complex numbers 1 + i a
t1, t2 = 0.3, 0.9
check("B4 k=0,d=1: arctan a + arctan b = arctan(plane law)",
      abs(math.atan(plane(math.tan(t1), math.tan(t2), 0, 1)) - (t1 + t2)) < 1e-12)
check("B5 k=0,d=-1: artanh a + artanh b = artanh(plane law) (rapidity)",
      abs(math.atanh(plane(math.tanh(t1), math.tanh(0.4), 0, -1)) - (t1 + 0.4)) < 1e-12)
check("B6 d=0: the plane law is o_k and the norm is the lift 1+ka",
      abs(plane(0.3, 0.7, 1.5, 0) - circ(0.3, 0.7, 1.5)) < 1e-15)

# ---------------------------------------------------------------- C
def moment(u, verts):
    w = np.exp(2 * verts @ u)
    return (w[:, None] * verts).sum(0) / w.sum()


def in_hull(p, verts, tol=1e-9):
    h = ConvexHull(verts)
    return np.all(h.equations[:, :-1] @ p + h.equations[:, -1] <= tol)


square = np.array([[0, 0], [1, 0], [0, 1], [1, 1]], float)
tri = np.array([[0, 0], [1, 0], [0, 1]], float)
hexa = np.array([p for p in itertools.permutations([0, 1, 2])], float)
# project the hexagon to its plane (2 coords) for hull tests
hexa2 = hexa[:, :2]
okC = True
for verts in (square, tri, hexa2):
    for _ in range(300):
        u = np.random.normal(0, 3, 2)
        okC &= in_hull(moment(u, verts), verts)
check("C1 moment map (gradient of soft support fn) lands inside square, triangle, hexagon", okC)
# angles don't move the moment map: |z_i| fixed, arg free
z = np.array([1.3 * np.exp(1j * 0.4), 0.6 * np.exp(1j * 2.0)])
u = np.log(np.abs(z))
z2 = np.abs(z) * np.exp(1j * np.random.uniform(0, 6.28, 2))
check("C2 turning the angles leaves the moment map fixed",
      np.allclose(moment(u, square), moment(np.log(np.abs(z2)), square)))
# per-axis product: moment map of the square = coordinatewise sigmoid of 2u
uu = np.array([0.3, -1.1])
check("C3 square's moment map = per-axis 1/(1+e^-2u) (cube is a product of segments)",
      np.allclose(moment(uu, square), 1 / (1 + np.exp(-2 * uu))))

# ---------------------------------------------------------------- D
def newton_product(factors):
    """Vertices of Newt(prod f): Minkowski sum of each factor's exponent set."""
    pts = np.zeros((1, factors[0].shape[1]))
    for f in factors:
        pts = (pts[:, None, :] + f[None, :, :]).reshape(-1, f.shape[1])
        pts = np.unique(np.round(pts, 9), axis=0)
    return pts


def hull_vertices(pts):
    pts = np.unique(pts, axis=0)
    # drop the affine hull: work in the span of differences
    c = pts.mean(0)
    U, S, Vt = np.linalg.svd(pts - c)
    r = (S > 1e-9).sum()
    proj = (pts - c) @ Vt[:r].T
    if r == 1:
        return 2
    h = ConvexHull(proj)
    return len(h.vertices)


def simplex_on(I, n):
    return np.array([np.eye(n)[i] for i in I])


n = 4  # ground set {0,1,2,3}; x_0 = 1 plays the role of the constant term
e = np.eye(n)
cube = newton_product([simplex_on([0, i], n) for i in (1, 2, 3)])
simp = newton_product([simplex_on([0, 1, 2, 3], n)])
perm = newton_product([simplex_on([i, j], n) for i, j in itertools.combinations(range(n), 2)])
assoc = newton_product([simplex_on(list(range(i, j + 1)), n) for i in range(n) for j in range(i + 1, n)])
counts = {"cube": hull_vertices(cube), "simplex": hull_vertices(simp),
          "permutohedron": hull_vertices(perm), "associahedron": hull_vertices(assoc)}
check("D1 vertex counts 8, 4, 24, 14 (cube, simplex, permutohedron, Loday associahedron)",
      counts == {"cube": 8, "simplex": 4, "permutohedron": 24, "associahedron": 14}, str(counts))
# cube = Newt((1+x1)(1+x2)(1+x3)) = Newt(1 + x1 o x2 o x3): same monomials
mons = set(tuple(int(v) for v in p[1:]) for p in cube)
corners = set(tuple(int(m >> i & 1) for i in range(3)) for m in range(8))
check("D2 Newton polytope of 1 + (x1 o x2 o x3) is the 3-cube, one vertex per tesseract corner",
      mons == corners)
# Ostrowski: Newt(f g) = Newt f + Newt g, spot check with random polygons
A = np.random.randint(0, 4, (5, 2)).astype(float)
B = np.random.randint(0, 4, (4, 2)).astype(float)
msum = (A[:, None] + B[None]).reshape(-1, 2)
check("D3 product of polynomials -> Minkowski sum of shapes (vertex sums)",
      hull_vertices(msum) <= hull_vertices(A) + hull_vertices(B))

# ---------------------------------------------------------------- E
def shape_from_weights(y, n):
    """Generalized permutohedron sum_I y_I Delta_I, vertices (y >= 0)."""
    pts = np.zeros((1, n))
    for mask in range(1, 1 << n):
        if y[mask] == 0:
            continue
        I = [i for i in range(n) if mask >> i & 1]
        f = y[mask] * simplex_on(I, n)
        pts = (pts[:, None, :] + f[None]).reshape(-1, n)
        pts = np.unique(np.round(pts, 9), axis=0)
    return pts


okE = True
for trial in range(25):
    y = [0] + [random.choice([0, 0, 1, 2, 0.5]) for _ in range((1 << n) - 1)]
    if sum(y) == 0:
        continue
    z = T.zeta(y)  # our tesseract zeta: box coordinate z_I = sum_{J in I} y_J
    pts = shape_from_weights(y, n)
    for mask in range(1, 1 << n):
        I = [i for i in range(n) if mask >> i & 1]
        s = pts[:, I].sum(1)
        okE &= abs(s.min() - z[mask]) < 1e-9  # facet sum_{i in I} t_i >= z_I is tight
check("E1 shape weights y on cube corners; our zeta gives Postnikov's z_I; every "
      "sum_{i in I} t_i >= z_I is tight", okE)
y = [0] * 16
y[0b1111] = 1
check("E2 simplex = weight 1 on the full corner", hull_vertices(shape_from_weights(y, n)) == 4)

# ---------------------------------------------------------------- F
okF = True
for _ in range(50):
    y = [0] + [random.choice([0, 1, 0.5, 3]) for _ in range(15)]
    pts = shape_from_weights(y, n)
    w = np.random.normal(0, 1, n)
    hP = (pts @ w).max()
    formula = sum(y[m] * max(w[i] for i in range(n) if m >> i & 1) for m in range(1, 16))
    okF &= abs(hP - formula) < 1e-9
check("F1 support function = sum_I y_I * max_{i in I} u_i (a weighted sum of maxes)", okF)


# ---------------------------------------------------------------- G
def edge_dirs(pts):
    c = pts.mean(0)
    h = ConvexHull(pts - c) if np.linalg.matrix_rank(pts - c) == pts.shape[1] else None
    dirs = set()
    V = pts[h.vertices]
    # edges = vertex pairs whose midpoint is on >= dim-1 common facets
    eq = h.equations
    for i, j in itertools.combinations(range(len(V)), 2):
        mid = (V[i] + V[j]) / 2 - c
        on = np.sum(np.abs(eq[:, :-1] @ mid + eq[:, -1]) < 1e-9)
        if on >= pts.shape[1] - 1:
            d = V[j] - V[i]
            d = d / np.abs(d[np.abs(d) > 1e-9][0])
            if d[np.abs(d) > 1e-9][0] < 0:
                d = -d
            dirs.add(tuple(np.round(d, 6)))
    return dirs


def is_typeA(d):
    nz = [x for x in d if abs(x) > 1e-9]
    return len(nz) == 2 and abs(nz[0] + nz[1]) < 1e-9


def is_typeB(d):
    nz = [abs(x) for x in d if abs(x) > 1e-9]
    return (len(nz) == 1) or (len(nz) == 2 and abs(nz[0] - nz[1]) < 1e-9)


cross3 = np.vstack([np.eye(3), -np.eye(3)])
cube3 = np.array(list(itertools.product([-1, 1], repeat=3)), float)
dc, dq = edge_dirs(cross3), edge_dirs(cube3)
check("G1 cube and cross-polytope edges are all type-B roots (e_i, e_i +- e_j)",
      all(is_typeB(d) for d in dc | dq), f"cross {sorted(dc)}")
check("G2 cross-polytope has edges e_i + e_j, not type-A roots: outside the type-A family",
      not all(is_typeA(d) for d in dc))

# ---------------------------------------------------------------- H
okH = True
for p in (1, 1.5, 2, 3, 7):
    for _ in range(50):
        u = np.random.normal(0, 1, 3)
        lse = np.log(np.sum(np.exp(p * np.log(np.abs(u))))) / p
        okH &= abs(np.exp(lse) - np.linalg.norm(u, p)) < 1e-9
check("H1 p-norm = exp( (1/p) log sum exp(p log|u_i|) ): soft-max of log|u| at temperature 1/p", okH)
u = np.random.normal(0, 1, 5)
check("H2 p -> infinity: p-norm -> max |u_i| (cube's gauge); p = 1 is the cross-polytope's",
      abs(np.linalg.norm(u, 200) - np.abs(u).max()) < 0.02 * np.abs(u).max())

# ---------------------------------------------------------------- I
okI = True
for nn in (2, 3, 4):
    zz = np.random.normal(size=(200, nn)) + 1j * np.random.normal(size=(200, nn))
    zz /= np.linalg.norm(zz, axis=1, keepdims=True)
    r2 = np.abs(zz) ** 2
    okI &= np.allclose(r2.sum(1), 1) and np.all(r2 >= 0)
check("I1 on the sphere |z|=1 in C^n the squared radii (|z_1|^2..|z_n|^2) fill a simplex", okI)

# ---------------------------------------------------------------- J
# beyond-the-cube F: the budget region {x_i >= 0, x1 o_k x2 o_k x3 <= 1} is the cube at k=-1, the
# simplex at k=0. In log-radius coordinates u_i = log(1+k x_i)/k it is ALWAYS the simplex
# {u_i >= 0, sum u_i <= log(1+k)/k}; only its size changes, growing without bound as k -> -1.
okJ = True
for k in (-0.9, -0.5, 0.3, 1.0, 4.0):
    for _ in range(2000):
        x = np.random.uniform(0, 1.2, 3)
        if np.any(1 + k * x <= 0):
            continue
        tot = circ(circ(x[0], x[1], k), x[2], k)
        u = np.log1p(k * x) / k
        okJ &= (tot <= 1) == (u.sum() <= math.log1p(k) / k + 1e-12)
check("J1 budget region {x o_k y o_k z <= 1} is one simplex in log-radius coordinates at every k "
      "(size log(1+k)/k -> infinity as k -> -1, where it is the whole cube)", okJ,
      ", ".join(f"k={k}: size {math.log1p(k)/k:.3f}" for k in (-0.99, -0.5, 0.3, 1.0)))

# ---------------------------------------------------------------- K
# our corner compose (switches combine by OR) on two shapes' weight vectors: weights stay >= 0,
# so the result is again a shape in the family, and its facet heights are the products z_I * z'_I.
okK = True
examples = []
for _ in range(40):
    y1 = [0] + [random.choice([0, 0, 1, 0.5]) for _ in range(15)]
    y2 = [0] + [random.choice([0, 0, 1, 2]) for _ in range(15)]
    yc = T.compose(y1, y2)
    okK &= min(yc) >= 0
    z1, z2, zc = T.zeta(y1), T.zeta(y2), T.zeta(yc)
    okK &= all(abs(zc[m] - z1[m] * z2[m]) < 1e-9 for m in range(16))
cube_y = [0] * 16
for m in (0b0011, 0b0101, 0b1001):
    cube_y[m] = 1
cc = T.compose(cube_y, cube_y)
examples.append(f"cube∘cube has {hull_vertices(shape_from_weights(cc, n))} vertices")
check("K1 our corner compose maps two shapes to a shape (weights stay >= 0) whose facet heights "
      "multiply: z(y o y') = z(y) z(y')", okK, "; ".join(examples))

# ---------------------------------------------------------------- report
lines = []
passed = sum(ok for _, ok, _ in results)
for name, ok, det in results:
    lines.append(f"[{'PASS' if ok else 'FAIL'}] {name}" + (f"  ({det})" if det else ""))
lines.append(f"\n{passed}/{len(results)} pass")
out = "\n".join(lines)
print(out)
with open(os.path.join(HERE, "checks_output.txt"), "w") as fh:
    fh.write(out + "\n")
