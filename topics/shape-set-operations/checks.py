"""DRAFT checks for shape-set-operations/README.md (Affine Composition Theory).

Run: python3 checks.py > checks_output.txt
Needs numpy, scipy, sympy, mpmath. Every check prints PASS or FAIL with the numbers it compared.
Monte Carlo checks use a fixed seed and a tolerance of about 5 standard errors.
"""
import math
import itertools
import numpy as np
from scipy import special
import sympy as sp
import mpmath as mp

rng = np.random.default_rng(20261126)
results = []


def check(name, ok, detail):
    results.append((name, bool(ok)))
    print(f"[{'PASS' if ok else 'FAIL'}] {name}: {detail}")


def circ(a, b, k):
    return a + b + k * a * b


# ---------------------------------------------------------------------------
# A. The fused operation is a union of independent sets, at every k != 0
# ---------------------------------------------------------------------------
# A1: P = -k a turns a o_k b into the union rule 1 - (1-P)(1-Q).
a, b, k = sp.symbols('a b k')
lhs = -k * circ(a, b, k)
rhs = 1 - (1 - (-k * a)) * (1 - (-k * b))
check("A1 union rule at every k", sp.expand(lhs - rhs) == 0, "-k(a o_k b) - [1-(1+ka)(1+kb)] simplifies to 0")

# A2: cylinders in a box of side M. phi(A x X  u  X x B)/M = a o_{-1/M} b, with A, B finite unions of intervals.
def interval_union_length(iv):
    iv = sorted(iv)
    tot, cur_s, cur_e = 0.0, None, None
    for s, e in iv:
        if cur_e is None or s > cur_e:
            if cur_e is not None:
                tot += cur_e - cur_s
            cur_s, cur_e = s, e
        else:
            cur_e = max(cur_e, e)
    return tot + (cur_e - cur_s if cur_e is not None else 0.0)


worst = 0.0
for trial in range(200):
    M = rng.uniform(0.5, 4.0)
    A = [tuple(sorted(rng.uniform(0, M, 2))) for _ in range(rng.integers(1, 4))]
    B = [tuple(sorted(rng.uniform(0, M, 2))) for _ in range(rng.integers(1, 4))]
    la, lb = interval_union_length(A), interval_union_length(B)
    # exact area of the union of the two cylinders on the M x M box
    union_area = la * M + M * lb - la * lb
    worst = max(worst, abs(union_area / M - circ(la, lb, -1.0 / M)))
check("A2 cylinders: union area / M = a o_(-1/M) b", worst < 1e-12, f"200 random boxes, worst gap {worst:.1e}")

# A3: the same by Monte Carlo on one picture (M = 2, A = [0.3,1.1], B = [0.5,1.9])
M, A, B = 2.0, (0.3, 1.1), (0.5, 1.9)
pts = rng.uniform(0, M, (400000, 2))
inA = (pts[:, 0] > A[0]) & (pts[:, 0] < A[1])
inB = (pts[:, 1] > B[0]) & (pts[:, 1] < B[1])
la, lb = A[1] - A[0], B[1] - B[0]
est = {"union": (inA | inB).mean() * M * M, "intersection": (inA & inB).mean() * M * M,
       "symmetric difference": (inA ^ inB).mean() * M * M, "difference A minus B": (inA & ~inB).mean() * M * M}
pred = {"union": M * circ(la, lb, -1 / M), "intersection": la * lb,
        "symmetric difference": M * circ(la, lb, -2 / M), "difference A minus B": la * (M - lb)}
se = M * M * math.sqrt(0.25 / len(pts))
for key in est:
    check(f"A3 {key} (Monte Carlo)", abs(est[key] - pred[key]) < 5 * se,
          f"estimate {est[key]:.4f}, formula {pred[key]:.4f}")

# A4: three cylinders in the unit cube: union volume = a o b o c at k = -1, and the corner sums are inclusion-exclusion.
la3 = [0.3, 0.5, 0.7]
pts = rng.uniform(0, 1, (400000, 3))
inside = np.stack([pts[:, i] < la3[i] for i in range(3)], axis=1)
vol_union = inside.any(axis=1).mean()
pred = circ(circ(la3[0], la3[1], -1), la3[2], -1)
check("A4 three cylinders = a o b o c at k = -1", abs(vol_union - pred) < 5 * math.sqrt(0.25 / len(pts)),
      f"estimate {vol_union:.4f}, a o b o c = {pred:.4f}")
corner = 0.0
for r in range(1, 4):
    for S in itertools.combinations(range(3), r):
        corner += (-1) ** (r + 1) * np.prod([la3[i] for i in S])
check("A4b corner sum of the 3-cube at k = -1 = inclusion-exclusion", abs(corner - pred) < 1e-15,
      f"corner sum {corner:.6f}")

# A5: idempotence <=> wall. a o a = a only at a = 0 or a = -1/k (a shape joined with itself sits on the wall).
sol = sp.solve(sp.Eq(circ(a, a, k), a), a)
check("A5 idempotent elements are 0 and the wall -1/k", set(sol) == {0, -1 / k}, f"solutions {sol}")

# A6: complement turns union into a product (De Morgan): M - (union) = (M-a)(M-b)/M on the box.
cA, cB = M - la, M - lb
comp = M * M - M * circ(la, lb, -1 / M)
check("A6 De Morgan on the box: what the union misses is the product of the two complements", abs(comp - cA * cB) < 1e-12,
      f"{comp:.6f} vs {cA * cB:.6f}")

# ---------------------------------------------------------------------------
# B. Positive k: the same union with Euler characteristic, where an open line counts -1
# ---------------------------------------------------------------------------
# A set on the line is a union of cells of a fixed decomposition: points (chi = +1) and open intervals (chi = -1).
# Product cells have chi = product. chi of a union of product cells = sum over the cells.
def chi_cells(cells):
    return sum(1 if c[0] == 'pt' else -1 for c in cells)


def decomposition(nbreaks):
    # breakpoints 1..nbreaks; cells: (-inf,1), {1}, (1,2), {2}, ..., {n}, (n,inf)
    cells = [('iv', 0)]
    for i in range(1, nbreaks + 1):
        cells += [('pt', i), ('iv', i)]
    return cells


def chi(c):
    return 1 if c[0] == 'pt' else -1


worst = 0
for trial in range(300):
    cells = decomposition(int(rng.integers(1, 5)))
    A = [c for c in cells if rng.random() < 0.5]
    B = [c for c in cells if rng.random() < 0.5]
    a_, b_ = chi_cells(A), chi_cells(B)
    X = chi_cells(cells)                 # the whole open line: -1
    union = 0
    for c1 in cells:
        for c2 in cells:
            if c1 in A or c2 in B:
                union += chi(c1) * chi(c2)
    gap = abs(union / X - circ(a_, b_, -1 / X))
    worst = max(worst, gap)
check("B1 chi(A x R  u  R x B) / chi(R) = a o_(+1) b, since chi(R) = -1", worst == 0 and X == -1,
      f"300 random pairs of sets on the open line, exact integers, worst gap {worst}")

# B2: three cylinders in R^3: chi(union)/chi(R)^2 = a o_1 b o_1 c (exact).
worst = 0
cells = decomposition(3)
for trial in range(100):
    S = [[c for c in cells if rng.random() < 0.5] for _ in range(3)]
    vals = [chi_cells(s) for s in S]
    union = 0
    for c in itertools.product(cells, repeat=3):
        if any(c[i] in S[i] for i in range(3)):
            union += chi(c[0]) * chi(c[1]) * chi(c[2])
    worst = max(worst, abs(union - circ(circ(vals[0], vals[1], 1), vals[2], 1)))
check("B2 three cylinders in R^3 give a o_1 b o_1 c", worst == 0, f"100 random triples, worst gap {worst}")

# B3: the sign of k is set by the ambient's Euler characteristic: closed segment (chi 1) -> k = -1,
# open line (chi -1) -> k = +1; any cell count: k = -1/chi(X).
table = []
for npts, niv in [(1, 0), (2, 1), (1, 1), (0, 1), (1, 2), (2, 5)]:
    X = npts - niv
    table.append((npts, niv, X, None if X == 0 else -1 / X))
check("B3 k = -1/chi(ambient)", table[1][3] == -1 and table[3][3] == 1,
      "; ".join(f"{p} pts + {i} open = chi {x}, k = {kk}" for p, i, x, kk in table))

# ---------------------------------------------------------------------------
# C. Positive k, second reading: growth beyond a reference box (Minkowski shell)
# ---------------------------------------------------------------------------
# C1: area of [0,1+a] x [0,1+b] minus the unit square = a o_1 b.
check("C1 grown box minus unit square = a o_1 b", sp.expand((1 + a) * (1 + b) - 1 - circ(a, b, 1)) == 0, "exact")
# C2: vol(K + xK) - vol(K) = vol(K)((1+x)^n - 1) for a ball (Monte Carlo, n = 3, x = 0.3): homothety, any convex K.
x, n = 0.3, 3
R = 1 + x
pts = rng.uniform(-R, R, (400000, n))
r = np.linalg.norm(pts, axis=1)
shell = ((r < R) & (r > 1)).mean() * (2 * R) ** n
ball = 4 / 3 * math.pi
check("C2 Minkowski shell of the ball = vol * n-fold o_1", abs(shell - ball * ((1 + x) ** n - 1)) < 5 * (2 * R) ** n * 0.5 / math.sqrt(len(pts)),
      f"shell {shell:.4f}, vol*((1+x)^3-1) = {ball * ((1 + x) ** n - 1):.4f}")

# ---------------------------------------------------------------------------
# D. The pair dial: any two shapes, k(A,B) = -phi(A n B)/(phi(A) phi(B))
# ---------------------------------------------------------------------------
def lens_area(r, d):
    if d >= 2 * r:
        return 0.0
    return 2 * r * r * math.acos(d / (2 * r)) - d / 2 * math.sqrt(4 * r * r - d * d)


worst = 0.0
for d in np.linspace(0, 2.2, 12):
    A = math.pi
    I = lens_area(1.0, d)
    kd = -I / (A * A)
    worst = max(worst, abs((2 * A - I) - circ(A, A, kd)))
check("D1 two unit disks: union area = a o_k(d) a with k(d) = -lens/area^2", worst < 1e-12,
      "identity holds at every separation (it is inclusion-exclusion with the overlap renamed)")
check("D2 coincident disks sit on the wall", abs(-1 / (-(lens_area(1, 0)) / math.pi ** 2) - math.pi) < 1e-12,
      "k(0) = -1/pi, so the wall -1/k = pi = the disk itself")

# D3: two hemispheres at angle theta on the sphere: pair dial k = -2(1 - theta/pi) (linear in theta).
ths = [0.0, math.pi / 4, math.pi / 2, 3 * math.pi / 4]
v = rng.normal(size=(400000, 3))
v /= np.linalg.norm(v, axis=1)[:, None]
ok = True
det = []
for th in ths:
    n1 = np.array([0, 0, 1.0])
    n2 = np.array([math.sin(th), 0, math.cos(th)])
    pa, pb = (v @ n1 > 0), (v @ n2 > 0)
    kk = -(pa & pb).mean() / (pa.mean() * pb.mean())
    pred = -2 * (1 - th / math.pi)
    ok &= abs(kk - pred) < 0.01
    det.append(f"{th:.2f}: {kk:.3f} vs {pred:.3f}")
check("D3 hemispheres: k(theta) = -2(1 - theta/pi), k = -1 at right angles", ok, "; ".join(det))

# ---------------------------------------------------------------------------
# E. Three shapes: where associativity can fail. Defect = triple overlap - k^2 abc
# ---------------------------------------------------------------------------
# For three cylinders the defect is 0 (A4). For three disks (a Venn diagram) it is not.
centers = np.array([[0, 0.5], [-0.45, -0.3], [0.45, -0.3]])
pts = rng.uniform(-2, 2, (1200000, 2))
ins = np.stack([np.linalg.norm(pts - c, axis=1) < 1 for c in centers], axis=1)
box = 16.0
area = ins.mean(axis=0) * box
pair = {(i, j): (ins[:, i] & ins[:, j]).mean() * box for i, j in itertools.combinations(range(3), 2)}
triple = ins.all(axis=1).mean() * box
union = ins.any(axis=1).mean() * box
ie = area.sum() - sum(pair.values()) + triple
check("E1 Venn of three disks: union = inclusion-exclusion", abs(union - ie) < 1e-9, f"union {union:.4f}")
kk = -pair[(0, 1)] / (area[0] * area[1])
fused = circ(circ(area[0], area[1], kk), area[2], kk)
defect = union - fused
check("E2 three disks are not a o b o c: the defect is the triple-overlap mismatch",
      abs(defect) > 0.05, f"union {union:.4f}, a o b o c {fused:.4f}, defect {defect:.4f} "
      f"(= triple {triple:.4f} - k^2 abc {kk**2*area.prod():.4f} up to the unequal pair overlaps)")

# ---------------------------------------------------------------------------
# F. e from the cube and a corner shape, along the ℓp exponent p
# ---------------------------------------------------------------------------
def orthant_vol(n, p):
    return math.exp(n * math.lgamma(1 + 1 / p) - math.lgamma(1 + n / p))


check("F1 p = 1: cube n corner simplex = 1/n!", all(abs(orthant_vol(n, 1) - 1 / math.factorial(n)) < 1e-15 for n in range(10)),
      "n = 0..9")
for p in [1, 2, 3]:
    E = sum(orthant_vol(n, p) for n in range(200))
    ml = float(mp.nsum(lambda j: mp.gamma(1 + mp.mpf(1) / p) ** j / mp.gamma(1 + j / mp.mpf(p)), [0, mp.inf]))
    # Monte Carlo: draws uniform on [0,1] until the sum of p-th powers passes 1
    N = 300000
    s = np.zeros(N)
    cnt = np.zeros(N)
    alive = np.ones(N, bool)
    while alive.any():
        s[alive] += rng.uniform(0, 1, alive.sum()) ** p
        cnt[alive] += 1
        alive &= s <= 1
    mc, sem = cnt.mean(), cnt.std() / math.sqrt(N)
    check(f"F2 p = {p}: sum of orthant volumes = Mittag-Leffler E_(1/p)(Gamma(1+1/p)) = average draws",
          abs(E - ml) < 1e-12 and abs(mc - E) < 5 * sem, f"series {E:.6f}, Mittag-Leffler {ml:.6f}, draws {mc:.4f} +- {sem:.4f}")
e2 = math.exp(math.pi / 4) * (1 + math.erf(math.sqrt(math.pi) / 2))
check("F3 p = 2 closed form e^(pi/4)(1 + erf(sqrt(pi)/2))", abs(e2 - sum(orthant_vol(n, 2) for n in range(200))) < 1e-12,
      f"{e2:.10f}")

# ---------------------------------------------------------------------------
# G. e along the homothety dial k: the cube cut by the o_k-ball  {x1 o_k ... o_k xn <= 1}
# ---------------------------------------------------------------------------
def E_k(kv):
    return (1 + kv) ** (1 / kv) + kv / (1 + kv)


def vol_k_series(kv, nmax=12):
    z = sp.symbols('z')
    kk_ = sp.nsimplify(kv) if isinstance(kv, int) else sp.Float(kv, 30)
    G = 1 + z * ((1 + kk_) * sp.exp(z * sp.log(1 + kk_) / kk_) - 1) / (kk_ + z)
    ser = sp.series(G, z, 0, nmax + 1).removeO()
    return [float(ser.coeff(z, j)) for j in range(nmax + 1)]


for kv in [-0.5, 0.5, 1.0, math.e - 1, 3.0]:
    N = 300000
    prod = np.ones(N)
    cnt = np.zeros(N)
    alive = np.ones(N, bool)
    thr = 1 + kv
    while alive.any():
        prod[alive] *= 1 + kv * rng.uniform(0, 1, alive.sum())
        cnt[alive] += 1
        alive &= (prod <= thr) if kv > 0 else (prod >= thr)
    mc, sem = cnt.mean(), cnt.std() / math.sqrt(N)
    vols = vol_k_series(kv, 25)
    check(f"G1 k = {kv:.4f}: average draws until the o_k-composite passes 1 = (1+k)^(1/k) + k/(1+k)",
          abs(mc - E_k(kv)) < 5 * sem and abs(sum(vols) - E_k(kv)) < 1e-6,
          f"draws {mc:.4f} +- {sem:.4f}, formula {E_k(kv):.6f}, sum of volumes {sum(vols):.6f}")
check("G2 k -> 0 gives e", abs(E_k(1e-8) - math.e) < 1e-6, f"{E_k(1e-8):.8f}")
check("G3 k = e - 1 gives Vandervelde's (e-1)/e + e^(1/(e-1)) (Math. Intelligencer 2008)",
      abs(E_k(math.e - 1) - ((math.e - 1) / math.e + math.exp(1 / (math.e - 1)))) < 1e-12, f"{E_k(math.e - 1):.10f}")
v1 = vol_k_series(1, 3)
check("G4 k = 1, n = 2: area under the hyperbola (1+x)(1+y) <= 2 is 2 ln 2 - 1", abs(v1[2] - (2 * math.log(2) - 1)) < 1e-12,
      f"{v1[2]:.10f}")
check("G5 limits: k -> -1 gives infinity, k -> infinity gives 2", E_k(-1 + 1e-9) > 15 and abs(E_k(1e9) - 2) < 1e-3,
      f"E(-1+1e-9) = {E_k(-1 + 1e-9):.2f}, E(1e9) = {E_k(1e9):.5f}")


ok = all(abs(E_k(1 / m) - ((1 + 1 / m) ** m + 1 / (m + 1))) < 1e-12 for m in range(1, 40))
check("G6 at k = 1/m the count is the compound-interest staircase (1+1/m)^m plus 1/(m+1)", ok,
      f"m = 1..39; k = 1 gives 5/2 = {E_k(1):.6f}, k = -1/2 gives 3 = {E_k(-0.5):.6f}")

# ---------------------------------------------------------------------------
# H. Minkowski sum of the cube and a ball: Steiner's formula with the cube's layers C(n,j)
# ---------------------------------------------------------------------------
def kappa(m):
    return math.pi ** (m / 2) / math.gamma(m / 2 + 1)


for n, t in [(2, 0.4), (3, 0.3)]:
    pts = rng.uniform(-t, 1 + t, (600000, n))
    dist = np.linalg.norm(pts - np.clip(pts, 0, 1), axis=1)
    est = (dist <= t).mean() * (1 + 2 * t) ** n
    pred = sum(math.comb(n, j) * kappa(n - j) * t ** (n - j) for j in range(n + 1))
    se = (1 + 2 * t) ** n * 0.5 / math.sqrt(len(pts))
    check(f"H1 n = {n}: vol(cube + tB) = sum C(n,j) kappa_(n-j) t^(n-j)", abs(est - pred) < 5 * se,
          f"estimate {est:.4f}, Steiner {pred:.4f}")

# ---------------------------------------------------------------------------
# I. Cube n cross-polytope: one family through e, the cuboctahedron and the 24-cell
# ---------------------------------------------------------------------------
def irwin_hall_cdf(n, s):
    # vol([0,1]^n n {sum <= s}): inclusion-exclusion over the cube's corners
    return sum((-1) ** j * math.comb(n, j) * max(s - j, 0) ** n for j in range(n + 1)) / math.factorial(n)


pts = rng.uniform(0, 1, (400000, 4))
est = (pts.sum(axis=1) <= 2.0).mean()
check("I1 vol(cube n {sum <= s}) by corner inclusion-exclusion (n = 4, s = 2)", abs(est - irwin_hall_cdf(4, 2)) < 5 * 0.5 / math.sqrt(4e5),
      f"Monte Carlo {est:.4f}, formula {irwin_hall_cdf(4, 2):.4f}")
check("I2 {|x_i| <= 1, sum|x_i| <= 2} in 4D has volume 8 (the 24-cell with vertices (+-1,+-1,0,0))",
      abs(16 * irwin_hall_cdf(4, 2) - 8) < 1e-12, f"16 x {irwin_hall_cdf(4, 2):.4f} = {16 * irwin_hall_cdf(4, 2):.4f}")
# vertices of {|x_i|<=1, sum|x_i|<=2}: count by brute force over candidate points with coordinates in {-1,0,1}
def vertices(n):
    vs = [v for v in itertools.product([-1, 0, 1], repeat=n) if sum(map(abs, v)) == 2]
    return len(vs)
check("I3 vertex counts: n = 3 gives 12 (cuboctahedron), n = 4 gives 24 (24-cell)", vertices(3) == 12 and vertices(4) == 24,
      f"{vertices(3)}, {vertices(4)}")
# Eulerian slabs: vol(cube n {j-1 <= sum <= j}) = A(n, j)/n!
def eulerian(n, m):
    return sum((-1) ** i * math.comb(n + 1, i) * (m + 1 - i) ** n for i in range(m + 2))
ok = all(abs((irwin_hall_cdf(n, j) - irwin_hall_cdf(n, j - 1)) - eulerian(n, j - 1) / math.factorial(n)) < 1e-12
         for n in range(1, 8) for j in range(1, n + 1))
check("I4 slabs of the cube between sum = j-1 and j are Eulerian numbers / n! (Laplace)", ok, "n = 1..7")

# ---------------------------------------------------------------------------
# J. Correlated sets and their limit: sets with ordinary sizes cap CHSH at 2; signed sizes reach 2 sqrt 2
# ---------------------------------------------------------------------------
from scipy.optimize import linprog
# global assignments: (a0, a1, b0, b1) in {+1,-1}^4, 16 of them. Target: Tsirelson correlations E(ai, bj).
G = list(itertools.product([1, -1], repeat=4))
target = {(0, 0): 1 / math.sqrt(2), (0, 1): 1 / math.sqrt(2), (1, 0): 1 / math.sqrt(2), (1, 1): -1 / math.sqrt(2)}
# joint outcome probabilities p(x, y | i, j) = (1 + x y E_ij)/4 with fair marginals
rows, rhs = [], []
for (i, jj), E in target.items():
    for xv in [1, -1]:
        for yv in [1, -1]:
            rows.append([1.0 if (g[i] == xv and g[2 + jj] == yv) else 0.0 for g in G])
            rhs.append((1 + xv * yv * E) / 4)
Aeq = np.array(rows)
# split w = u - v with u, v >= 0; minimize total mass sum(u + v)
c = np.ones(32)
res = linprog(c, A_eq=np.hstack([Aeq, -Aeq]), b_eq=np.array(rhs), bounds=[(0, None)] * 32, method="highs")
neg = res.x[16:].sum()
chsh = sum(target[(0, 0)] + target[(0, 1)] + target[(1, 0)] - target[(1, 1)] for _ in [0])
res_pos = linprog(np.zeros(16), A_eq=Aeq, b_eq=np.array(rhs), bounds=[(0, None)] * 16, method="highs")
check("J1 Tsirelson's CHSH = 2 sqrt 2 has no realization by sets with ordinary sizes", res_pos.status == 2,
      f"CHSH {chsh:.4f}; linear program with sizes >= 0 is infeasible (Bell, Pitowsky's correlation polytope)")
check("J2 it does have a realization with signed sizes (Abramsky-Brandenburger 2011, Thm 5.5)", res.status == 0 and neg > 1e-6,
      f"least total negative size {neg:.4f} (total size {res.x[:16].sum() - neg:.4f})")
# classical bound: every deterministic assignment (a corner of the 4-cube) gives |CHSH| <= 2
vals = [g[0] * g[2] + g[0] * g[3] + g[1] * g[2] - g[1] * g[3] for g in G]
check("J3 corners of the 4-cube give CHSH in {-2, 2}", set(vals) == {-2, 2}, f"{sorted(set(vals))}")

# ---------------------------------------------------------------------------
# K. The other set operations, at the most granular level (two independent sets, sizes a, b, dial k)
# ---------------------------------------------------------------------------
a, b, c, k, q = sp.symbols('a b c k q')
P = lambda x: -k * x          # share of the ambient
size = lambda Pf: sp.simplify(-Pf / k)
union = size(1 - (1 - P(a)) * (1 - P(b)))
xor = size(P(a) + P(b) - 2 * P(a) * P(b))
inter = size(P(a) * P(b))
diff = size(P(a) * (1 - P(b)))
comp = size(1 - P(a))
check("K1 coordinator's algebra: union a o_k b, XOR a o_2k b, AND -k ab, A minus B a + k ab, complement -1/k - a",
      sp.simplify(union - circ(a, b, k)) == 0 and sp.simplify(xor - circ(a, b, 2 * k)) == 0 and sp.simplify(inter + k * a * b) == 0
      and sp.simplify(diff - (a + k * a * b)) == 0 and sp.simplify(comp - (-1 / k - a)) == 0, "all five exact")

# K2: all 16 Boolean operations on two independent sets; which are commutative and associative (and not constant)?
names = {(0, 0, 0, 1): "AND", (0, 1, 1, 1): "OR", (0, 1, 1, 0): "XOR", (1, 0, 0, 1): "XNOR", (1, 1, 1, 0): "NAND", (1, 0, 0, 0): "NOR",
         (0, 0, 1, 1): "first", (0, 1, 0, 1): "second", (0, 0, 1, 0): "A minus B", (0, 1, 0, 0): "B minus A",
         (1, 1, 0, 1): "B implies A", (1, 0, 1, 1): "A implies B", (1, 1, 0, 0): "not first", (1, 0, 1, 0): "not second",
         (0, 0, 0, 0): "false", (1, 1, 1, 1): "true"}
def op_size(tt, x, y):
    # tt lists f(0,0), f(0,1), f(1,0), f(1,1) with first argument = A
    Px, Py = P(x), P(y)
    w = [(1 - Px) * (1 - Py), (1 - Px) * Py, Px * (1 - Py), Px * Py]
    return sp.expand(size(sum(t * wi for t, wi in zip(tt, w))))
monoids = []
for tt, nm in names.items():
    f = lambda x, y: op_size(tt, x, y)
    comm = sp.simplify(f(a, b) - f(b, a)) == 0
    assoc = sp.simplify(f(f(a, b), c) - f(a, f(b, c))) == 0
    const = sp.simplify(sp.diff(f(a, b), a)) == 0 and sp.simplify(sp.diff(f(a, b), b)) == 0
    if comm and assoc and not const:
        monoids.append(nm)
check("K2 of the 16 Boolean operations, exactly four are commutative, associative and non-constant",
      sorted(monoids) == ["AND", "OR", "XNOR", "XOR"], f"{sorted(monoids)}")
cmap = lambda x: -1 / k - x
check("K3 AND and XNOR are OR and XOR seen in the mirror (the complement, a reflection about -1/(2k))",
      sp.simplify(cmap(circ(cmap(a), cmap(b), k)) - inter) == 0
      and sp.simplify(cmap(circ(cmap(a), cmap(b), 2 * k)) - op_size((1, 0, 0, 1), a, b)) == 0, "conjugation by the complement")
check("K4 De Morgan in sizes: complement(a o_k b) = complement(a) AND complement(b)",
      sp.simplify(cmap(circ(a, b, k)) - (-k * cmap(a) * cmap(b))) == 0, "exact")

# K5: distributivity of AND over the union, on numbers and with the right dial
AND = lambda x, y: -k * x * y
lhs = AND(a, circ(b, c, k))
naive = circ(AND(a, b), AND(a, c), k)
gap = sp.factor(sp.expand(lhs - naive))
right = circ(AND(a, b), AND(a, c), -1 / a)
check("K5 distributivity fails on numbers by -k^2 abc(1 + ka): zero only for crisp a (0 or the wall)",
      sp.simplify(gap - (-k ** 2 * a * b * c * (1 + k * a))) == 0, f"gap = {gap}")
check("K5b it holds exactly when the two pieces use the dial of their shared ambient A: k' = -1/a",
      sp.simplify(lhs - right) == 0, "a AND (b o_k c) = (a AND b) o_(-1/a) (a AND c)")
pts = rng.uniform(0, 1, (400000, 3))
pa, pb, pc = 0.6, 0.3, 0.45
A_, B_, C_ = pts[:, 0] < pa, pts[:, 1] < pb, pts[:, 2] < pc
mc = (A_ & (B_ | C_)).mean()
x_, y_ = pa * pb, pa * pc
check("K5c Monte Carlo: P(A and (B or C)) = (AB) o_(-1/P_A) (AC), not (AB) o_(-1) (AC)",
      abs(mc - (x_ + y_ - x_ * y_ / pa)) < 5 * 0.5 / math.sqrt(4e5) and abs(mc - (x_ + y_ - x_ * y_)) > 0.01,
      f"MC {mc:.4f}; shared-ambient dial {x_ + y_ - x_ * y_ / pa:.4f}; plain k = -1 {x_ + y_ - x_ * y_:.4f}")

# K6: difference is the o_k monoid acting on sizes; a one-dial 'second slot weight' q is associative only at q = 1
act = lambda x, y: x * (1 + k * y)
check("K6 A minus B is an action: a minus (b o_k c) = (a minus b) minus c", sp.simplify(act(a, circ(b, c, k)) - act(act(a, b), c)) == 0,
      "(A minus B) minus C = A minus (B or C)")
Fq = lambda x, y: x + q * y + k * x * y
cond = sp.Poly(sp.expand(Fq(Fq(a, b), c) - Fq(a, Fq(b, c))), a, b, c).coeffs()
solq = sp.solve(cond, [q], dict=True)
check("K6b the order-style family a + q b + k ab is associative only at q = 1 (the union)", solq == [{q: 1}], f"solutions {solq}")
check("K7 the asymmetry of difference is the odds ratio: (A minus B)/(B minus A) = odds(A)/odds(B)",
      sp.simplify(act(a, b) / act(b, a) - (a / (1 + k * a)) / (b / (1 + k * b))) == 0, "odds in units: a/(1 + ka)")
# K8: n-fold XOR = piling-up lemma = o_2k fold
Ps = [0.1, 0.35, 0.2, 0.4]
pile = (1 - np.prod([1 - 2 * p for p in Ps])) / 2
fold = 0.0
for p in Ps:
    fold = circ(fold, p, -2)
check("K8 n-fold XOR is the o_(-2) fold = Matsui's piling-up lemma", abs(pile - fold) < 1e-15, f"{pile:.6f} vs {fold:.6f}")
# K9: pair dial and the odds ratio both mark the independence saddle inside the tetrahedron of the four Venn cells
okk = True
for _ in range(2000):
    w = rng.dirichlet([1, 1, 1, 1])     # p00, p01, p10, p11
    PA, PB = w[2] + w[3], w[1] + w[3]
    kp = -w[3] / (PA * PB)
    orat = w[3] * w[0] / (w[1] * w[2])
    okk &= (abs(kp + 1) < 1e-9) == (abs(orat - 1) < 1e-9)
    okk &= np.sign(-kp - 1) == np.sign(orat - 1)
check("K9 pair dial below -1 <=> odds ratio above 1 (positive overlap); = -1 <=> independent (the saddle surface)",
      okk, "2000 random correlated pairs")

print()
npass = sum(ok for _, ok in results)
print(f"{npass}/{len(results)} checks pass")
