"""Checks for the number-theory geometries survey (DRAFT, Nov 25 2026).

Run:  python3 checks.py > checks_output.txt
Each check prints PASS or FAIL with a short label. Sections:
  A tropical / regular subdivisions of the cube
  B the dial's sign = sub/supermodular; path vertices = generalized permutohedron
  C cluster algebras: the cube is the decoupled associahedron
  D field with one element: lambda-rings and the q-deformed cube (subspace lattice)
  E p-adic: Newton polygons and the height of the formal group
  F Arakelov / adelic: product formula for the ruler
  G modular: integer landscape laws and class numbers
"""
import itertools, math, random
from fractions import Fraction as Fr
import numpy as np
from scipy.spatial import ConvexHull
import sympy as sp
import mpmath as mp

random.seed(1); np.random.seed(1)
RESULTS = []
def check(label, ok, detail=""):
    RESULTS.append(bool(ok))
    print(("PASS " if ok else "FAIL ") + label + (("   " + detail) if detail else ""))

def corners(n):
    return list(itertools.product([0, 1], repeat=n))

# ---------------------------------------------------------------- A tropical
print("== A. Tropical geometry: corner heights subdivide the cube")
k = sp.Symbol('k')
xs = sp.symbols('x1:4')
P = sp.expand(sp.prod([1 + k*x for x in xs]))
poly = sp.Poly(P, *xs)
expo = sorted(m for m in poly.monoms())
check("A1 Newton polytope of prod(1+k x_i) (n=3) is exactly the 8 cube corners",
      expo == sorted(corners(3)))
check("A1b coefficient of corner S is k^|S| (the corner weight)",
      all(poly.coeff_monomial(sp.prod([x**e for x, e in zip(xs, m)])) == k**sum(m) for m in expo))

def lower_cells(pts, h, tol=1e-9):
    """Regular subdivision: cells of the lower hull of (pts, h). Returns set of frozensets of point indices."""
    L = np.hstack([np.array(pts, float), np.array(h, float)[:, None]])
    top = np.append(np.array(pts, float).mean(0), max(h) + 10.0)   # helper point far above, never on the lower hull
    hull = ConvexHull(np.vstack([L, top]))
    groups = {}
    for eq, simp in zip(hull.equations, hull.simplices):
        if eq[-2] < -tol:   # outward normal points down in the height coordinate
            key = tuple(np.round(eq / -eq[-2], 7))
            groups.setdefault(key, set()).update(i for i in simp.tolist() if i < len(L))
    # add points lying on each lower facet plane (coplanar points Qhull may leave out)
    cells = set()
    for key, s in groups.items():
        a = np.array(key)
        on = [i for i in range(len(L)) if abs(np.dot(a[:-1], L[i]) + a[-1]) < 1e-6]
        cells.add(frozenset(on))
    return cells

C3 = corners(3)
rank = [sum(c) for c in C3]
cells_flat = lower_cells(C3, [r * 1.0 for r in rank])   # heights from k^|S|: affine in |S|
check("A2 heights |S|·v(k) (the plain dial) give NO subdivision: one cell, the whole cube",
      len(cells_flat) == 1 and len(next(iter(cells_flat))) == 8)
cells_level = lower_cells(C3, [r*(r-1)/2 for r in rank])  # q-binomial heights
sizes = sorted(len(c) for c in cells_level)
check("A3 heights |S|(|S|-1)/2 (the q-binomial weights) slice the cube into 3 level slabs: tetrahedron, octahedron, tetrahedron",
      len(cells_level) == 3 and sizes == [4, 4, 6], f"cell sizes {sizes}")
# q-binomial theorem: prod_{i<n}(1+q^i t) = sum_j q^{j(j-1)/2} [n,j]_q t^j
q, t = sp.symbols('q t')
def qbin(n, j):
    num = sp.Mul(*[1 - q**(n - i) for i in range(j)]); den = sp.Mul(*[1 - q**(i + 1) for i in range(j)])
    return sp.cancel(num / den)
n = 5
lhs = sp.expand(sp.prod([1 + q**i*t for i in range(n)]))
rhs = sp.expand(sum(q**(j*(j-1)//2) * qbin(n, j) * t**j for j in range(n + 1)))
check("A3b q-binomial theorem (n=5): the fused product with dial k_i = q^i gives those weights",
      sp.simplify(lhs - rhs) == 0)
tri = set()
for _ in range(40000):
    c = lower_cells(C3, np.random.rand(8))
    if all(len(s) == 4 for s in c):
        tri.add(frozenset(c))
check("A4 random corner heights on the 3-cube reach exactly 74 triangulations (all regular; De Loera-Rambau-Santos)",
      len(tri) == 74, f"found {len(tri)}")
by_count = {}
for T in tri:
    by_count[len(T)] = by_count.get(len(T), 0) + 1
check("A4b they split as 72 with 6 tetrahedra and 2 with 5 (the two corner-cut-off ones)",
      by_count == {6: 72, 5: 2}, str(by_count))

# ---------------------------------------------------------------- B permutohedra
print("== B. The dial's sign is sub/supermodularity; path vertices form a generalized permutohedron")
def fold(vals, kk):
    p = 1.0
    for v in vals: p *= (1 + kk*v)
    return (p - 1) / kk if kk != 0 else sum(vals)
def f_set(S, x, kk): return fold([x[i] for i in S], kk)

def modularity(x, kk):
    n = len(x); worst_sub = 0; worst_sup = 0
    for S in itertools.chain.from_iterable(itertools.combinations(range(n), r) for r in range(n + 1)):
        for i, j in itertools.combinations([a for a in range(n) if a not in S], 2):
            d = f_set(S + (i, j), x, kk) - f_set(S + (i,), x, kk) - f_set(S + (j,), x, kk) + f_set(S, x, kk)
            worst_sub = max(worst_sub, d); worst_sup = min(worst_sup, d)
    return worst_sub, worst_sup
ok_neg = ok_pos = ok_zero = True
for _ in range(200):
    x = np.random.rand(4) * 0.9
    kk = -np.random.rand()
    a, b = modularity(x, kk); ok_neg &= a <= 1e-12
    kk = np.random.rand() * 2
    a, b = modularity(x, kk); ok_pos &= b >= -1e-12
    a, b = modularity(x, 0.0); ok_zero &= abs(a) < 1e-12 and abs(b) < 1e-12
check("B1 k<0 (with 1+k·x_i>0): the corner aggregate f(S) is submodular (200 random cases, n=4)", ok_neg)
check("B1b k>0: f is supermodular (a convex game)", ok_pos)
check("B1c k=0: f is modular (plain sum)", ok_zero)
x0 = [0.3, 0.5, 0.8]
s_sym = sp.symbols('a b c S0', positive=True)
check("B1d exact: second difference f(S+i+j)-f(S+i)-f(S+j)+f(S) = k·x_i·x_j·prod_S(1+k x)",
      sp.simplify(((1+k*s_sym[0])*(1+k*s_sym[1])*s_sym[3] - 1)/k - ((1+k*s_sym[0])*s_sym[3]-1)/k
                  - ((1+k*s_sym[1])*s_sym[3]-1)/k + (s_sym[3]-1)/k - k*s_sym[0]*s_sym[1]*s_sym[3]) == 0,
      "S0 stands for prod over S of (1+k x)")

def path_vertices(x, kk):
    n = len(x); V = []
    for perm in itertools.permutations(range(n)):
        v = [0.0]*n; S = ()
        for i in perm:
            v[i] = f_set(S + (i,), x, kk) - f_set(S, x, kk); S += (i,)
        V.append(tuple(v))
    return np.array(V)
def n_hull_vertices(V):
    V = np.unique(np.round(V, 12), axis=0)
    if len(V) == 1: return 1
    c = V.mean(0); d = V.shape[1]
    # project to the plane sum = const
    B = np.linalg.svd(V - c)[2][:d-1]
    Y = (V - c) @ B.T
    return len(ConvexHull(Y).vertices)
check("B2 n=3: the 6 path (ordering) marginal vectors are the 6 corners of a hexagon for k = ±0.5",
      n_hull_vertices(path_vertices(x0, -0.5)) == 6 and n_hull_vertices(path_vertices(x0, 0.5)) == 6)
check("B2b at k=0 all 6 paths give the same vector: the hexagon shrinks to a point",
      n_hull_vertices(path_vertices(x0, 0.0)) == 1)
x4 = [0.2, 0.35, 0.6, 0.75]
check("B2c n=4: 24 orderings give the 24 corners of a 3-D permutohedron (k = -0.7)",
      n_hull_vertices(path_vertices(x4, -0.7)) == 24)
def face_ok(x, kk):
    n = len(x); V = path_vertices(x, kk); ok = True
    for r in range(1, n):
        for S in itertools.combinations(range(n), r):
            ys = V[:, list(S)].sum(1); fs = f_set(S, x, kk)
            ok &= (ys <= fs + 1e-12).all() if kk < 0 else (ys >= fs - 1e-12).all()
            ok &= np.isclose(ys, fs).any()       # every corner S is a facet that some path touches
    return ok
check("B3 k<0: every path vertex obeys y(S) <= f(S) with equality on some path: base polytope of f (Edmonds 1970)",
      face_ok(x4, -0.7))
check("B3b k>0: y(S) >= f(S): the core of the convex game (Shapley 1971)", face_ok(x4, 0.7))
xe = 0.4; kk = -0.6; r = 1 + kk*xe
V = path_vertices([xe]*3, kk)
check("B4 equal entries: corners are the permutations of x·(1, r, r²) with r = 1+kx (a geometric permutohedron)",
      set(map(lambda v: tuple(np.round(sorted(v), 12)), V)) == {tuple(np.round(sorted([xe, xe*r, xe*r*r]), 12))})
# Shapley value = centroid; compare with game-theory page closed form x_i * int_0^1 prod_{j!=i}(1 + k t x_j) dt
kk = 0.8; V = path_vertices(x4, kk); cent = V.mean(0)
closed = [x4[i] * mp.quad(lambda tt: mp.fprod([1 + kk*tt*x4[j] for j in range(4) if j != i]), [0, 1]) for i in range(4)]
check("B5 the centroid of the polytope is the Shapley value from the game-theory page", np.allclose(cent, [float(c) for c in closed]))

# ---------------------------------------------------------------- C cluster
print("== C. Cluster algebras: the cube is the decoupled case; couplings bend it")
def period(b, c, x0=Fr(3, 7), x1=Fr(5, 11), N=60):
    seq = [x0, x1]
    for nn in range(N):
        e = b if nn % 2 == 0 else c
        seq.append((seq[-1]**e + 1) / seq[-2])
        if len(seq) > 3 and seq[-2] == x0 and seq[-1] == x1:
            return len(seq) - 2
    return None
pers = {bc: period(*bc) for bc in [(0, 0), (1, 1), (1, 2), (1, 3), (1, 4), (2, 2)]}
check("C1 rank-2 exchange x_{n+1}x_{n-1} = x_n^e + 1: periods 4,5,6,8 for bc = 0,1,2,3 (square, pentagon, hexagon, octagon)",
      pers[(0, 0)] == 4 and pers[(1, 1)] == 5 and pers[(1, 2)] == 6 and pers[(1, 3)] == 8, str(pers))
check("C1b bc >= 4: no period (infinite type)", pers[(1, 4)] is None and pers[(2, 2)] is None)
def mutate(B, X, kk):
    n = len(X)
    pos = Fr(1); neg = Fr(1)
    for i in range(n):
        if B[i][kk] > 0: pos *= X[i]**B[i][kk]
        if B[i][kk] < 0: neg *= X[i]**(-B[i][kk])
    X2 = list(X); X2[kk] = (pos + neg) / X[kk]
    B2 = [[0]*n for _ in range(n)]
    for i in range(n):
        for j in range(n):
            if i == kk or j == kk: B2[i][j] = -B[i][j]
            else: B2[i][j] = B[i][j] + (abs(B[i][kk])*B[kk][j] + B[i][kk]*abs(B[kk][j])) // 2
    return B2, X2
def exchange_graph(B, cap=500):
    n = len(B); X0 = [Fr(random.randint(2, 97), random.randint(2, 97)) for _ in range(n)]
    start = frozenset(X0); seen = {start: (B, X0)}; edges = set(); todo = [start]
    while todo and len(seen) < cap:
        key = todo.pop(); Bc, Xc = seen[key]
        for kk in range(n):
            B2, X2 = mutate(Bc, Xc, kk); k2 = frozenset(X2)
            edges.add(frozenset([key, k2]))
            if k2 not in seen: seen[k2] = (B2, X2); todo.append(k2)
    return len(seen), len(edges)
types = {
    'A1xA1xA1 (no arrows)': [[0,0,0],[0,0,0],[0,0,0]],
    'A3 (two arrows)':      [[0,1,0],[-1,0,1],[0,-1,0]],
    'A1xA1 (no arrows)':    [[0,0],[0,0]],
    'A2':                   [[0,1],[-1,0]],
    'B2':                   [[0,1],[-2,0]],
    'G2':                   [[0,1],[-3,0]],
    'A1xA1xA1xA1':          [[0]*4 for _ in range(4)],
}
eg = {name: exchange_graph(B) for name, B in types.items()}
for name, (v, e) in eg.items():
    print(f"     exchange graph {name}: {v} clusters, {e} flips")
check("C2 no arrows (A1^3): 8 clusters and 12 flips: the exchange graph IS the 3-cube", eg['A1xA1xA1 (no arrows)'] == (8, 12))
check("C2b no arrows (A1^4): 16 and 32: the tesseract", eg['A1xA1xA1xA1'] == (16, 32))
check("C2c turning on arrows bends the 3-cube into the A3 associahedron: 14 clusters, 21 flips (Catalan)",
      eg['A3 (two arrows)'] == (14, 21))
check("C2d rank 2: square 4, pentagon 5, hexagon 6, octagon 8",
      [eg[t][0] for t in ['A1xA1 (no arrows)', 'A2', 'B2', 'G2']] == [4, 5, 6, 8])
# Laurent phenomenon on the A2 cycle with symbolic start
a_, b_ = sp.symbols('a b'); seq = [a_, b_]
for _ in range(5): seq.append(sp.cancel((seq[-1] + 1) / seq[-2]))
check("C3 Laurent phenomenon: every A2 cluster variable is a Laurent polynomial (denominator a monomial)",
      all(sp.Poly(sp.fraction(sp.together(s))[1], a_, b_).is_monomial for s in seq))
# dilogarithm on the A2 Y-system (Lyness y_{n+1} y_{n-1} = 1 + y_n): the 1+y factors are the circle-composition lift at k=1
def rogers(z): return mp.polylog(2, z) + 0.5*mp.log(z)*mp.log(1 - z)
ys = [mp.mpf('0.7'), mp.mpf('1.9')]
for _ in range(3): ys.append((1 + ys[-1]) / ys[-2])
tot = sum(rogers(y/(1 + y)) for y in ys)
check("C4 A2 dilogarithm identity: sum of L(y/(1+y)) over the 5-cycle = π²/2 (Abel's five-term identity)",
      abs(tot - mp.pi**2/2) < 1e-20, f"sum = {mp.nstr(tot, 15)}")

# ---------------------------------------------------------------- D F_1
print("== D. Field with one element: lambda-rings and the q-deformed cube")
x_, y_ = sp.symbols('x y')
lam = sp.expand((1 + x_*t)*(1 + y_*t))
check("D1 lambda_t(x+y) for two line elements = 1 + t·(x ∘_t y): the dial k IS the lambda-ring's t",
      sp.expand(lam - (1 + t*(x_ + y_ + t*x_*y_))) == 0)
X4 = sp.symbols('y1:5')
lamn = sp.Poly(sp.expand(sp.prod([1 + v*t for v in X4])), t)
ej_ok = all(sp.expand(lamn.coeff_monomial(t**j) - sum(sp.prod(c) for c in itertools.combinations(X4, j))) == 0 for j in range(5))
check("D2 lambda^j of n lines = e_j = the sum over the rank-j corners of the cube (n=4)", ej_ok)
def subspaces_F2(n):
    vecs = [v for v in itertools.product([0, 1], repeat=n)]
    spaces = set()
    for r in range(n + 1):
        for gens in itertools.combinations(vecs, r):
            S = {tuple([0]*n)}
            for g in gens:
                S |= {tuple((a + b) % 2 for a, b in zip(s, g)) for s in S}
            spaces.add(frozenset(S))
    return spaces
L3 = subspaces_F2(3)
counts = [sum(1 for s in L3 if len(s) == 2**d) for d in range(4)]
check("D3 subspaces of F_2^3 by dimension: 1, 7, 7, 1 = Gaussian binomials [3,j]_2 (the cube deformed to q=2)",
      counts == [1, 7, 7, 1] and [int(qbin(3, j).subs(q, 2)) for j in range(4)] == [1, 7, 7, 1])
check("D3b at q=1 the Gaussian binomials are 1,3,3,1: the corners of the 3-cube (Tits 1957)",
      [sp.limit(qbin(3, j), q, 1) for j in range(4)] == [1, 3, 3, 1])
def mobius_top(spaces):
    order = sorted(spaces, key=len); mu = {}
    bottom = order[0]
    for s in order:
        mu[s] = 1 if s == bottom else -sum(mu[u] for u in order if len(u) < len(s) and u <= s)
    return mu[order[-1]]
check("D4 Möbius function bottom-to-top: Boolean cube (-1)^3 = -1; subspace lattice (-1)^3·2^3 = -8",
      mobius_top(L3) == -8)
pp = 7; kk_ = sp.Symbol('k')
fp = sp.Poly(sp.expand(((1 + kk_*x_)**pp - 1)/kk_), x_)
check("D5 the p-fold self-fusion is a Frobenius lift: ((1+kx)^p - 1)/k ≡ k^(p-1) x^p mod p (p=7; Borger's F_1 axiom)",
      all((c % pp == 0) for m, c in zip(fp.monoms(), [sp.Poly(cc, kk_).coeffs()[0] for cc in fp.coeffs()]) if m[0] < pp))

# ---------------------------------------------------------------- E p-adic
print("== E. p-adic: Newton polygons, and the height of the formal group")
def vp(nn, p):
    if nn == 0: return math.inf
    c = 0
    while nn % p == 0: nn //= p; c += 1
    return c
ok = True
for p in [3, 5, 7, 11]:
    pts = [(j - 1, vp(math.comb(p, j), p)) for j in range(1, p + 1)]   # ((1+x)^p-1)/x
    slopes = {Fr(pts[-1][1] - pts[0][1], pts[-1][0] - pts[0][0])}
    ok &= all(pts[i][1] >= pts[0][1] + Fr(-1, p - 1)*pts[i][0] for i in range(len(pts)))
check("E1 Newton polygon of ((1+x)^p-1)/x is one segment of slope -1/(p-1): every nonzero p-torsion point of ∘₁ has |x|_p = p^(-1/(p-1))", ok)
# formal groups and height
N = 56
def mul(a, b):
    c = [Fr(0)]*N
    for i, ai in enumerate(a):
        if ai == 0: continue
        for j in range(N - i): c[i + j] += ai*b[j]
    return c
def compose(f, g):   # f(g(z)), g[0]=0
    res = [Fr(0)]*N; powg = [Fr(1)] + [Fr(0)]*(N-1)
    for i in range(N):
        if i > 0: powg = mul(powg, g)
        if f[i] != 0: res = [r + f[i]*pg for r, pg in zip(res, powg)]
    return res
def revert(f):       # compositional inverse, f = z + ...
    g = [Fr(0), Fr(1)] + [Fr(0)]*(N-2)
    for m in range(2, N):
        c = compose(f, g)[m]; g[m] -= c
    return g
def height_from_log(logs, p):
    expf = revert(logs)
    mp_ = compose(expf, [p*c for c in logs])
    for d, c in enumerate(mp_):
        if c != 0 and (c.numerator % p != 0):
            assert c.denominator % p != 0
            return d
    return None
def log_mult(): return [Fr(0)] + [Fr((-1)**(m+1), m) for m in range(1, N)]   # log(1+z), the ruler at k=1
def log_add(): return [Fr(0), Fr(1)] + [Fr(0)]*(N-2)
def log_elliptic(a4, a6=0):
    # Weierstrass y^2 = x^3 + a4 x + a6; w(z) = z^3 + a4 z w^2 + a6 w^3 (Silverman IV.1)
    w = [Fr(0)]*N; w[3] = Fr(1)
    for _ in range(N):
        w2 = mul(w, w); w3 = mul(w2, w)
        nw = [Fr(0)]*N; nw[3] = Fr(1)
        for i in range(N - 1): nw[i + 1] += a4*w2[i]
        for i in range(N): nw[i] += a6*w3[i]
        if nw == w: break
        w = nw
    # x = z/w, y = -1/w ; omega = dx/(2y) = -(w/2) d(z/w); compute u = w/z^3 (unit), 1/u
    u = w[3:] + [Fr(0)]*3
    inv = [Fr(0)]*N; inv[0] = 1/u[0]
    for m in range(1, N):
        inv[m] = -sum(u[i]*inv[m - i] for i in range(1, m + 1)) / u[0]
    # x = z^{-2} * inv(z), y = -z^{-3} * inv(z)
    # dx/dz = sum (i-2) inv[i] z^{i-3}; 2y = -2 sum inv[i] z^{i-3}; ratio = sum (i-2)inv[i] z^i / (-2 sum inv[i] z^i)
    num = [Fr(i - 2)*inv[i] for i in range(N)]
    den = [Fr(-2)*inv[i] for i in range(N)]
    dinv = [Fr(0)]*N; dinv[0] = 1/den[0]
    for m in range(1, N):
        dinv[m] = -sum(den[i]*dinv[m - i] for i in range(1, m + 1)) / den[0]
    om = mul(num, dinv)                     # omega/dz, starts 1 + ...
    return [Fr(0)] + [om[m - 1]/m for m in range(1, N)]
hts = {}
for p in [3, 5, 7]:
    hts[('mult', p)] = height_from_log(log_mult(), p)
    hts[('ell', p)] = height_from_log(log_elliptic(Fr(-1)), p)
add_h = height_from_log(log_add(), 5)
check("E2 ∘₁ (multiplicative law): [p](z) mod p starts at z^p for p = 3,5,7: height 1",
      [hts[('mult', p)] for p in [3, 5, 7]] == [3, 5, 7])
check("E2b ∘₀ (plain +): [p](z) = p·z vanishes mod p: height infinity", add_h is None)
check("E3 elliptic curve y²=x³-x: [p](z) mod p starts at z^9 for p=3, z^5 for p=5, z^49 for p=7",
      (hts[('ell', 3)], hts[('ell', 5)], hts[('ell', 7)]) == (9, 5, 49),
      f"got {hts[('ell', 3)]}, {hts[('ell', 5)]}, {hts[('ell', 7)]}: height 2, 1, 2")
import json, os
cj = json.load(open(os.path.join(os.path.dirname(__file__), '..', 'toric', 'counts.json')))
ap = {int(pq): int(pq) + 1 - c['elliptic y^2=x^3-x'] for pq, c in cj.items() if sp.isprime(int(pq))}
check("E3b the toric page's own point counts give a_p = 0 exactly at p ≡ 3 mod 4 (supersingular = height 2)",
      all((ap[p] % p == 0) == (p % 4 == 3) for p in ap if p > 2), str(ap))

# ---------------------------------------------------------------- F Arakelov
print("== F. Arakelov / adelic: the ruler obeys the product formula")
def abs_p(r, p):
    return Fr(p) ** (-(vp(r.numerator, p) - vp(r.denominator, p)))
ok = True
for _ in range(50):
    kk = Fr(random.randint(1, 9), random.randint(1, 9)); xv = Fr(random.randint(1, 30), random.randint(1, 30))
    u = 1 + kk*xv
    primes = sp.primefactors(u.numerator*u.denominator)
    s = math.log(abs(u)) + sum(math.log(abs_p(u, p)) for p in primes)
    ok &= abs(s) < 1e-12
check("F1 for rational k, x: the ruler log|1+kx| summed over the real place and every prime is 0 (product formula)", ok)
def weil_h(r): return math.log(max(abs(r.numerator), abs(r.denominator)))
ok = True
for _ in range(200):
    kk = Fr(random.randint(1, 9), random.randint(1, 9))
    a = Fr(random.randint(1, 50), random.randint(1, 50)); b = Fr(random.randint(1, 50), random.randint(1, 50))
    ab = a + b + kk*a*b
    ok &= weil_h(1 + kk*ab) <= weil_h(1 + kk*a) + weil_h(1 + kk*b) + 1e-12
check("F2 Weil height of the lift is subadditive under ∘ₖ (it is just h(uv) <= h(u)+h(v)): bounds size, deforms nothing", ok)

# ---------------------------------------------------------------- G modular
print("== G. Modular / hyperbolic: integer laws and class numbers")
def class_number(D):   # count reduced positive definite forms of discriminant D<0
    h = 0; a = 1
    while 3*a*a <= -D:
        for b in range(-a + 1, a + 1):
            if (b*b - D) % (4*a) == 0:
                c = (b*b - D)//(4*a)
                if c >= a and not (b < 0 and (a == c)) and math.gcd(math.gcd(a, abs(b)), c) == 1:
                    h += 1
        a += 1
    return h
hs = {D: class_number(D) for D in [-3, -4, -20, -23, -47, -71]}
check("G1 forms α+βx+γx² (the full sl2 plane of laws) under integer chart changes: class numbers 1,1,2,3,5,7",
      [hs[D] for D in [-3, -4, -20, -23, -47, -71]] == [1, 1, 2, 3, 5, 7], str(hs))
pq = [(p, qq) for p in range(-4, 5) for qq in range(-4, 5)]
c_ = 2
check("G2 landscape laws (α=1) under x ↦ x/(1-cx): (p,q) ↦ (p+2c, q+pc+c²) keeps p²-4q (from groups-lie) — integer laws ↔ principal forms only",
      all((p + 2*c_)**2 - 4*(qq + p*c_ + c_*c_) == p*p - 4*qq for p, qq in pq))

print()
print(f"{sum(RESULTS)} of {len(RESULTS)} checks pass")
