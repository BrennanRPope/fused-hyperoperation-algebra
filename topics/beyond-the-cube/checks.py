"""Checks for beyond-the-cube/README.md (DRAFT, Nov 26 2026).

Run: python3 checks.py > checks_output.txt
Needs numpy, scipy, sympy, mpmath.
"""
import itertools, math, json
import numpy as np
import sympy as sp
import mpmath as mp
from scipy.spatial import ConvexHull

rng = np.random.default_rng(20261126)
results = {}
passed = 0
total = 0


def check(name, ok, detail=""):
    global passed, total
    total += 1
    passed += bool(ok)
    print(("PASS " if ok else "FAIL ") + name + ("  | " + detail if detail else ""))
    results[name] = {"ok": bool(ok), "detail": detail}


def circ(a, b, k):
    return a + b + k * a * b


# ---------------------------------------------------------------- A. cube corners, cut polytope, elliptope
print("\n== A. Classical versus quantum correlations: cube corners, cut polytope, elliptope")
signs = list(itertools.product([1, -1], repeat=3))
cut_pts = sorted({(s[0]*s[1], s[0]*s[2], s[1]*s[2]) for s in signs})
check("A1 three-variable cut polytope has 4 vertices, the cube corners with x*y*z = +1",
      len(cut_pts) == 4 and all(p[0]*p[1]*p[2] == 1 for p in cut_pts), str(cut_pts))


def corr(x, y, z):
    return np.array([[1, x, y], [x, 1, z], [y, z, 1]], float)


check("A2a every cut vertex is a valid correlation matrix (inside the elliptope)",
      all(np.linalg.eigvalsh(corr(*p)).min() > -1e-12 for p in cut_pts))
ev = np.linalg.eigvalsh(corr(-.5, -.5, -.5))
check("A2b (-1/2,-1/2,-1/2) is in the elliptope (eigenvalues 0, 1.5, 1.5) but outside the tetrahedron (x+y+z = -1.5 < -1)",
      ev.min() > -1e-12 and -1.5 < -1, str(np.round(ev, 6)))
# elliptope volume vs tetrahedron vs cube, Monte Carlo
P = rng.uniform(-1, 1, (400000, 3))
x, y, z = P.T
inE = 1 + 2*x*y*z - x*x - y*y - z*z >= 0
inT = (x + y + z >= -1) & (x - y - z >= -1) & (-x + y - z >= -1) & (-x - y + z >= -1)
vol_E = 8 * inE.mean(); vol_T = 8 * inT.mean()
check("A2c volumes: cube 8, elliptope ~ pi^2/2 = 4.935, tetrahedron 8/3",
      abs(vol_E - math.pi**2/2) < 0.05 and abs(vol_T - 8/3) < 0.05 and np.all(inE[inT]),
      f"elliptope {vol_E:.3f}, tetrahedron {vol_T:.3f}; tetrahedron inside elliptope")


def chsh(a1, a2, b1, b2):
    return a1@b1 + a1@b2 + a2@b1 - a2@b2


best_det = max(s[0]*s[2] + s[0]*s[3] + s[1]*s[2] - s[1]*s[3] for s in itertools.product([1, -1], repeat=4))
check("A3a CHSH over the cube's corners (deterministic +-1 answers) peaks at 2", best_det == 2)
v = lambda t: np.array([math.cos(t), math.sin(t)])
q = chsh(v(0), v(math.pi/2), v(math.pi/4), v(-math.pi/4))
check("A3b CHSH with unit vectors (elliptope points) reaches 2*sqrt(2) = 2.828", abs(q - 2*math.sqrt(2)) < 1e-12, f"{q:.6f}")
best = 0
for _ in range(20000):
    V = rng.normal(size=(4, 4)); V /= np.linalg.norm(V, axis=1, keepdims=True)
    best = max(best, chsh(*V))
check("A3c 20,000 random unit-vector strategies in 4D never beat 2*sqrt(2) (Tsirelson 1980)",
      best <= 2*math.sqrt(2) + 1e-9, f"best {best:.4f}")
check("A3d no-signalling box (Popescu-Rohrlich) scores 4", 3*1 - (-1) == 4)
# independence surface inside the correlation polytope COR(2)
T = np.array([[0, 0, 0], [1, 0, 0], [0, 1, 0], [1, 1, 1]], float)
ok = True
for _ in range(5000):
    p, qq = rng.uniform(0, 1, 2)
    pt = np.array([p, qq, p*qq])
    # inside tetrahedron iff barycentric coords >= 0
    M = np.vstack([(T[1:] - T[0]).T])
    lam = np.linalg.solve(M, pt - T[0])
    ok &= lam.min() >= -1e-12 and lam.sum() <= 1 + 1e-12
    ok &= abs(circ(p, qq, -1) - (p + qq - p*qq)) < 1e-15
check("A4 our union-rule world (p, q, pq) is a saddle surface inside Pitowsky's correlation tetrahedron COR(2)", ok)

# ---------------------------------------------------------------- B. associahedron coordinates obey the union rule
print("\n== B. The associahedron's u-coordinates (binary geometry) and the union rule")


def u_vars(z):
    n = len(z)
    zz = lambda a, b: z[a % n] - z[b % n]
    U = {}
    for i in range(n):
        for j in range(i+2, n):
            if (j + 1) % n == i:
                continue
            U[(i, j)] = (zz(i, j+1)*zz(i+1, j)) / (zz(i, j)*zz(i+1, j+1))
    return U


def crosses(X, Y):
    i, j = X; a, b = Y
    return (i < a < j < b) or (a < i < b < j)


worst = 0; worst_union = 0; inside = True
for n in (5, 6, 7, 8):
    for _ in range(200):
        z = np.sort(rng.uniform(0, 1, n))
        U = u_vars(z)
        for X, uX in U.items():
            prodY = np.prod([uY for Y, uY in U.items() if crosses(X, Y)])
            worst = max(worst, abs(uX + prodY - 1))
            un = 0.0
            for Y, uY in U.items():
                if crosses(X, Y):
                    un = circ(un, 1 - uY, -1)
            worst_union = max(worst_union, abs(uX - un))
            inside &= 0 < uX < 1
check("B1 u-equations u_X + prod(u_Y over chords Y crossing X) = 1 hold for n = 5..8, all chords",
      worst < 1e-12 and inside, f"max error {worst:.1e}; every u in (0,1)")
check("B2 same equations in our notation: u_X = union (circ_-1) of p_Y = 1 - u_Y over crossing chords",
      worst_union < 1e-12, f"max error {worst_union:.1e}")


def eta(a1, a2, b1, b2):
    return (a2 - a1)*(b2 - b1)/((b1 - a1)*(b2 - a2))


z = np.sort(rng.uniform(0, 1, 6)); U = u_vars(z)
e = max(abs(U[(i, j)] - (1 - eta(z[i], z[i+1], z[j], z[(j+1)] if j+1 < 6 else z[0]+1))) for (i, j) in U if j+1 < 6)
check("B3 each u_ij is 1 - (AdS3 cross-ratio of polygon sides i and j), the probability-entanglement page's x",
      e < 1e-12, f"max error {e:.1e}")
werr = 0
for _ in range(2000):
    a1, a2, b1, m, b2 = np.sort(rng.uniform(0, 1, 5))
    werr = max(werr, abs(eta(a1, a2, b1, b2) - circ(eta(a1, a2, b1, m), eta(a1, a2, m, b2), -1)))
check("B4 gluing: x(A, B1B2) = x(A,B1) circ_-1 x(A,B2) for adjacent B1, B2 (same as that page's H4)",
      werr < 1e-12, f"max error {werr:.1e}")
# B5: the n=4 equation alone, then gluing, gives every u-equation
def u_arcs(x1, x2, y1, y2):
    return (x1 - y2)*(x2 - y1)/((x1 - y1)*(x2 - y2))


derr = 0; gerr = 0
for _ in range(500):
    n = 7; z = np.sort(rng.uniform(0, 1, n)); U = u_vars(z)
    for (i, j) in U:
        p1, p2, p3, p4 = z[i], z[i+1], z[j], z[(j+1) % n]
        u_ac = u_arcs(p1, p2, p3, p4)          # side i with side j
        u_bd = u_arcs(p2, p3, p4, p1)          # block i+1..j with block j+1..i
        derr = max(derr, abs(u_ac + u_bd - 1))
        prodY = np.prod([uY for Y, uY in U.items() if crosses((i, j), Y)])
        gerr = max(gerr, abs(u_bd - prodY))
check("B5 every u-equation = the four-point identity u(a,c) + u(b,d) = 1, with u(b,d) glued from the crossing chords (derivation, ours as written)",
      derr < 1e-12 and gerr < 1e-12, f"identity error {derr:.1e}, gluing error {gerr:.1e}")

# B6 ABHY pentagon canonical function = phi^3 five-point amplitude
c13, c14, c24 = rng.uniform(0.5, 2, 3)
# coordinates (X13, X14); facets as linear forms L(X) = a.X + b >= 0
facets = {
    "X13": (np.array([1., 0]), 0.),
    "X14": (np.array([0., 1]), 0.),
    "X24": (np.array([-1., 1]), c13),          # X24 = c13 + X14 - X13
    "X25": (np.array([-1., 0]), c13 + c14),    # X25 = c13 + c14 - X13
    "X35": (np.array([0., -1]), c14 + c24),    # X35 = c14 + c24 - X14
}
order = ["X13", "X14", "X24", "X25", "X35"]  # cyclic order around the pentagon (adjacent = compatible)


def L(name, X):
    a, b = facets[name]; return a @ X + b


def vertex(f1, f2):
    a1, b1 = facets[f1]; a2, b2 = facets[f2]
    return np.linalg.solve(np.vstack([a1, a2]), -np.array([b1, b2]))


verts = [vertex(order[i], order[(i+1) % 5]) for i in range(5)]
feasible = all(min(L(f, V) for f in facets) > -1e-12 for V in verts)
hull = ConvexHull(np.array(verts))
check("B6a ABHY region in the (X13, X14) plane is a pentagon: 5 feasible vertices, all on the hull",
      feasible and len(hull.vertices) == 5)


def tri_canon(P0, P1, P2, X):
    # canonical function of triangle = |det| (area form) / product of the 3 edge lines, normalised
    Y = np.array([1, *X]); Ws = []
    pts = [P0, P1, P2]
    for k in range(3):
        A = np.array([1, *pts[(k+1) % 3]]); B = np.array([1, *pts[(k+2) % 3]])
        W = np.cross(A, B)
        if W @ np.array([1, *pts[k]]) < 0:
            W = -W
        Ws.append(W)
    num = abs(np.linalg.det(np.vstack(Ws)))
    return num / np.prod([W @ Y for W in Ws])


def m5(X):
    v = {f: L(f, X) for f in facets}
    return sum(1/(v[order[i]]*v[order[(i+1) % 5]]) for i in range(5))


err = 0
cen = np.mean(verts, axis=0)
for _ in range(200):
    w = rng.dirichlet(np.ones(5)); X = w @ np.array(verts)
    tri = sum(tri_canon(verts[0], verts[i], verts[i+1], X) for i in (1, 2, 3))
    err = max(err, abs(tri / m5(X) - 1))
check("B6b canonical function of the ABHY pentagon = 1/(X13X14) + 1/(X14X24)... = five-point phi^3 tree amplitude (ABHY 2017)",
      err < 1e-10, f"max relative error {err:.1e}")

# B7 ABHY n = 6: 3D associahedron
chords6 = [(i, j) for i in range(1, 7) for j in range(i+2, 7) if not (i == 1 and j == 6)]
Xs = sp.symbols("x13 x14 x15")
cvals = {}
for (i, j) in [(1, 3), (1, 4), (1, 5), (2, 4), (2, 5), (3, 5)]:
    cvals[(i, j)] = sp.Rational(int(rng.integers(2, 6)), 2)
Xexpr = {(1, 3): Xs[0], (1, 4): Xs[1], (1, 5): Xs[2]}


def getX(i, j):
    if j - i == 1 or (i == 1 and j == 6):
        return sp.Integer(0)
    return Xexpr[(i, j)]


# X_{i+1,j+1} = c_ij + X_{i,j+1} + X_{i+1,j} - X_ij
for (i, j) in [(1, 3), (1, 4), (1, 5), (2, 4), (2, 5), (3, 5)]:
    Xexpr[(i+1, j+1)] = sp.expand(cvals[(i, j)] + getX(i, j+1) + getX(i+1, j) - getX(i, j))
lin = {ch: sp.Poly(Xexpr[ch], *Xs) for ch in chords6}
A_ = np.array([[float(lin[ch].coeff_monomial(s)) for s in Xs] for ch in chords6])
b_ = np.array([float(lin[ch].coeff_monomial(1)) for ch in chords6])
V6 = []
for tri3 in itertools.combinations(range(9), 3):
    Mx = A_[list(tri3)]
    if abs(np.linalg.det(Mx)) < 1e-9:
        continue
    pt = np.linalg.solve(Mx, -b_[list(tri3)])
    if (A_ @ pt + b_).min() > -1e-9:
        V6.append(tuple(np.round(pt, 9)))
V6 = sorted(set(V6))
h6 = ConvexHull(np.array(V6))
edges6 = set()
for simp, eq in zip(h6.simplices, h6.equations):
    pass
# count edges: pairs of vertices sharing two tight facets
tight = [set(np.where(abs(A_ @ np.array(v) + b_) < 1e-9)[0]) for v in V6]
for a_, b2 in itertools.combinations(range(len(V6)), 2):
    if len(tight[a_] & tight[b2]) >= 2:
        edges6.add((a_, b2))
check("B7 ABHY region for n = 6 is the 3D associahedron: 14 vertices, 21 edges, 9 facets",
      len(V6) == 14 and len(edges6) == 21, f"{len(V6)} vertices, {len(edges6)} edges")
json.dump({"V6": V6, "edges": sorted(edges6), "A": A_.tolist(), "b": b_.tolist(),
           "chords": chords6, "pent": [v.tolist() for v in verts]},
          open("assoc6.json", "w"))

# ---------------------------------------------------------------- C. Y-system pentagon and Lee-Yang
print("\n== C. The pentagon as an exact RG flow (Y-system, Lee-Yang)")
y = [0.37, 2.9]
for _ in range(10):
    y.append((1 + y[-1]) / y[-2])
check("C1 Y-system y(n+1) y(n-1) = 1 + y(n) repeats every 5 steps (the pentagon; Zamolodchikov periodicity)",
      abs(y[5] - y[0]) < 1e-12 and abs(y[6] - y[1]) < 1e-12, str(np.round(y[:6], 4)))
phi = (1 + math.sqrt(5))/2
Lr = lambda t: mp.polylog(2, t) + 0.5*mp.log(t)*mp.log(1 - t)
c_eff = 6/mp.pi**2 * Lr(1/(1 + phi))
check("C2 fixed point y = golden ratio; Rogers L(1/(1+y)) = pi^2/15 gives c_eff = 2/5 (scaling Lee-Yang model)",
      abs(c_eff - mp.mpf(2)/5) < 1e-15, f"c_eff = {mp.nstr(c_eff, 15)}")
ya, yb = rng.uniform(0, 3, 2)
check("C3 the Y-system's factors 1 + y are our lift: (1+ya)(1+yb) = 1 + (ya circ_1 yb)",
      abs((1+ya)*(1+yb) - (1 + circ(ya, yb, 1))) < 1e-12)

# ---------------------------------------------------------------- D. double copy at 4 points (scalar part)
print("\n== D. Double copy: two ordered pieces make a gravity piece")
s, t = sp.symbols("s t"); u = -s - t
m_same = 1/s + 1/t          # m(1234|1234): channels s, t
m_other = -1/s              # m(1234|1243): shared channel s
S_klt = sp.simplify(1/m_other)
check("D1 KLT kernel at 4 points is the inverse of the bi-adjoint (associahedron) matrix: S = -s",
      sp.simplify(S_klt + s) == 0)
A1234 = 1/(s*t); A1243 = 1/(s*u)
M = sp.simplify(A1234 * S_klt * A1243)
sym = all(sp.simplify(M - M.subs({s: a, t: b}, simultaneous=True)) == 0 for a, b in [(t, s), (s, u), (u, t)])
check("D2 gauge-stripped double copy M = A(1234) S A(1243) = -1/(s t u)-type: fully symmetric in s, t, u",
      sp.simplify(M - 1/(s*t*(s+t))) == 0 and sym, str(M))

# ---------------------------------------------------------------- E. total positivity (amplituhedron family)
print("\n== E. Total positivity: where the amplituhedron lives")


def toeplitz_minors(seq, size=5, maxk=3):
    T_ = np.zeros((size, size))
    for i in range(size):
        for j in range(size):
            if 0 <= i - j < len(seq):
                T_[i, j] = seq[i - j]
    mins = []
    for r in range(1, maxk+1):
        for R in itertools.combinations(range(size), r):
            for C in itertools.combinations(range(size), r):
                mins.append(np.linalg.det(T_[np.ix_(R, C)]))
    return min(mins)


def esym(xs):
    co = np.array([1.0])
    for xi in xs:
        co = np.convolve(co, [1, xi])
    return co


okp = True; neg_found = False
for _ in range(60):
    xs = rng.uniform(0.1, 2, 4); kk = rng.uniform(0.1, 2)
    w = esym(xs) * kk**np.arange(5)
    okp &= toeplitz_minors(w) > -1e-9
    w2 = esym(xs) * (-kk)**np.arange(5)
    neg_found |= toeplitz_minors(w2) < -1e-6
check("E1 for k > 0 and positive inputs, the corner weights k^m e_m form a totally nonnegative (Polya frequency) sequence",
      okp, "all minors up to 3x3 of the 5x5 Toeplitz matrix >= 0 (Aissen-Schoenberg-Whitney, Edrei 1952)")
check("E2 for k < 0 (union side) total positivity fails unless inputs are flipped x -> -x", neg_found)

# ---------------------------------------------------------------- F. what the cube becomes as k turns: budget region
print("\n== F. The budget region {x circ_k y <= 1}: cube at k = -1, simplex at k = 0")


def area(k):
    if abs(k) < 1e-12:
        return 0.5
    return ((1 + k)*math.log(1 + k) - k)/k**2 if k > -1 else 1.0


g = np.linspace(0, 1, 2001)
okF = True
for k in (-0.999999, -0.6, 0.0, 1.0, 5.0, 50.0):
    X_, Y_ = np.meshgrid(g, g)
    num = np.mean(circ(X_, Y_, k) <= 1 + 1e-12)
    okF &= abs(num - area(k)) < 2e-3
check("F1 area of the budget square region = ((1+k)ln(1+k) - k)/k^2: 1 at k = -1, 1/2 at k = 0, -> 0 as k grows",
      okF, f"areas k=-0.6: {area(-0.6):.4f}, k=1: {area(1):.4f}, k=50: {area(50):.4f}")
Pc = rng.uniform(0, 1, (300000, 3))
vol = lambda k: np.mean([1]) if False else np.mean((np.prod(1 + k*Pc, axis=1) - 1)/k <= 1) if k != 0 else np.mean(Pc.sum(1) <= 1)
check("F2 in 3D the same region is the whole cube at k = -1 and the simplex (volume 1/6) at k = 0",
      abs(vol(-1) - 1) < 1e-12 and abs(vol(0) - 1/6) < 3e-3, f"{vol(-1):.4f}, {vol(0):.4f}")
kk = 0.7; gx = rng.uniform(0, 1); gy = ((1 + kk)/(1 + kk*gx) - 1)/kk
check("F3 its boundary is Sugeno's lambda-measure normalisation (1+k gx)(1+k gy) = 1 + k (known, Sugeno 1974)",
      abs(circ(gx, gy, kk) - 1) < 1e-12)

# ---------------------------------------------------------------- G. Birkhoff, permutohedron, zonotope
print("\n== G. Birkhoff polytope, permutohedron, zonotope")
perms = list(itertools.permutations(range(3)))
img = {tuple(np.eye(3)[list(p)] @ np.array([1., 2, 3])) for p in perms}
check("G1 the permutohedron is the Birkhoff polytope's shadow under M -> M(1,2,3): 6 distinct corners", len(img) == 6)
segs = [np.array([1., -1, 0]), np.array([1., 0, -1]), np.array([0., 1, -1])]
sums = [sum(c*sg for c, sg in zip(cs, segs)) for cs in itertools.product([0, 1], repeat=3)]
B2 = np.array([[1/math.sqrt(2), -1/math.sqrt(2), 0], [1/math.sqrt(6), 1/math.sqrt(6), -2/math.sqrt(6)]])
hz = ConvexHull(np.array(sums) @ B2.T)
check("G2 the permutohedron is a zonotope: the Minkowski sum of the 3 segments e_i - e_j is a hexagon", len(hz.vertices) == 6)
Bm = 0.5*np.array([[0, 1, 1], [1, 0, 1], [1, 1, 0]])
r1 = np.sqrt(Bm[0]); r2 = np.sqrt(Bm[1])
# |<r1, r2 e^{i theta}>| can never vanish: only one overlapping entry
check("G3 a doubly stochastic matrix that no unitary produces (not unistochastic): the quantum subset of Birkhoff is smaller",
      abs((r1*r2).sum() - 0.5) < 1e-12 and np.count_nonzero(r1*r2) == 1)

print(f"\n{passed}/{total} checks passed")
json.dump({"passed": passed, "total": total, "results": results}, open("results.json", "w"), indent=1)
