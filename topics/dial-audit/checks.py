"""Dial audit (DRAFT, Nov 26 2026): is each of the seven dials fundamental or derived?

Three tests, applied to every dial:
  T1 reach / rank (fixed ruler): does turning the dial reach laws or shapes the others can't,
     and is the Jacobian from the dials to the object full rank?
  T2 units: does the dial change when the ruler's unit is changed (x -> lam*x)? A dial that only
     moves with the unit is a choice of unit, not a number the law carries.
  T3 known reductions: exact identities that write the dial through the others.

Every law here is a one-dimensional composition law F(a, b) with F(a, 0) = a. Its generator is the
vector field P(x) = dF/db (x, 0) = 1 + c1 x + c2 x^2 + c3 x^3 + ...  (Lie; Aczel 1948). The law is
recovered from P by L' = 1/P and F(a, b) = L^-1(L(a) + L(b)), so the coefficients c1, c2, c3, ... are
coordinates on the space of laws in a fixed ruler.

Run: python3 checks.py  (sympy, numpy, mpmath). Writes results.json.
"""
import os, json, math, random, itertools
import numpy as np
import sympy as sp
import mpmath as mp

random.seed(7); np.random.seed(7); mp.mp.dps = 30
RESULTS, PASS = {}, []


def check(name, ok, detail=""):
    PASS.append(bool(ok))
    print(("PASS " if ok else "FAIL ") + name + ("  [" + detail + "]" if detail else ""))


x, a, b, k, d, lam, s = sp.symbols('x a b k d lambda s')
N = 6  # series order for generator coefficients


def trunc(expr, n=N):
    return sp.series(expr, x, 0, n + 1).removeO()


def compose(f, g, n=N):
    """f(g(x)) truncated, f and g polynomials in x with no constant term."""
    return sp.expand(trunc(sp.expand(f.subs(x, g)), n))


def iterate_coeffs(G, n=N, kmax=None):
    """Coefficients of the regular iterate G^s for G = x + O(x^2): the coefficients of the integer
    iterates G^m are polynomials in m of degree <= n-1, so interpolate them in m and put s for m."""
    pts = list(range(0, n + 1))
    its = [x]
    for m in range(1, n + 1):
        its.append(compose(G, its[-1], n))
    coeffs = {}
    for j in range(2, n + 1):
        vals = [sp.expand(it).coeff(x, j) for it in its]
        coeffs[j] = sp.expand(sp.interpolate(list(zip(pts, vals)), s))
    return x + sum(coeffs[j] * x**j for j in range(2, n + 1))


def inverse_series(G, n=N):
    """Compositional inverse of G = x + ..., truncated."""
    H = x
    for _ in range(n):
        H = sp.expand(H - (compose(G, H, n) - x))
        H = trunc(H, n)
    return sp.expand(H)


def generator_of_rung(G, n=N):
    """Rung s of the law with generator G (rung 0 = +, rung 1 = the law): F_s(a,b) = G^s(G^-s a + G^-s b).
    Its own generator is P_s(x) = (G^s)'(G^-s(x)). Returns the series in x with s symbolic."""
    Gs = iterate_coeffs(G, n)
    Gms = Gs.subs(s, -s)
    dGs = sp.diff(Gs, x)
    return sp.expand(trunc(sp.expand(dGs.subs(x, Gms)), n - 1))


# ---------------------------------------------------------------------------
print("== A. k and d: the first two coordinates of the generator")
F_kd = (a + b + k * a * b) / (1 - d * a * b)
P_kd = sp.simplify(sp.diff(F_kd, b).subs(b, 0).subs(a, x))
check("A1 the plane law's generator is P(x) = 1 + k x + d x^2", sp.expand(P_kd - (1 + k * x + d * x**2)) == 0, str(P_kd))

# the plane law solves dF/db = P(F)/P(b) (it is the flow of P), so P determines F
lhs = sp.diff(F_kd, b)
rhs = (1 + k * F_kd + d * F_kd**2) / (1 + k * b + d * b**2)
check("A2 dF/db = P(F)/P(b) holds exactly, so (k, d) fix the law and the law fixes (k, d)", sp.simplify(lhs - rhs) == 0)

Jac = sp.Matrix([[sp.diff(P_kd.coeff(x, i), v) for v in (k, d)] for i in (1, 2)])
check("A3 Jacobian of (c1, c2) in (k, d) is the identity: rank 2, neither is built from the other", Jac == sp.eye(2))

# units: the same law read in a ruler with unit lam
G_units = sp.simplify(F_kd.subs({a: lam * a, b: lam * b}) / lam)
F_scaled = (a + b + (lam * k) * a * b) / (1 - (lam**2 * d) * a * b)
check("A4 changing the unit x -> lam x sends k -> lam k and d -> lam^2 d (same law, new ruler)", sp.simplify(G_units - F_scaled) == 0)

# general generator: c_n -> lam^n c_n
cs = sp.symbols('c1:7')
P_gen = 1 + sum(cs[i] * x**(i + 1) for i in range(6))
P_resc = sp.expand(P_gen.subs(x, lam * x))
check("A5 for any law, a unit change sends c_n -> lam^n c_n", all(sp.expand(P_resc.coeff(x, n) - lam**n * cs[n - 1]) == 0 for n in range(1, 7)))

inv = sp.simplify((lam**2 * d) / (lam * k)**2 - d / k**2)
check("A6 d/k^2 is unchanged by every unit change; it is the plane's one unit-free number (k != 0)", inv == 0)
Delta_over = sp.simplify(((lam * k)**2 - 4 * lam**2 * d) / (lam * k)**2 - (k**2 - 4 * d) / k**2)
check("A7 so is Delta/k^2 = 1 - 4 d/k^2; the curvature -Delta/4 itself scales like 1/unit^2 (a reading, not unit-free)", Delta_over == 0)

cvel = sp.symbols('c', positive=True)
u, v = sp.symbols('u v')
vel = (u + v) / (1 + u * v / cvel**2)
check("A8 Einstein's velocity addition is k = 0, d = -1/c^2; in units with c = 1, d = -1 (only the sign of d is unit-free at k = 0)",
      sp.simplify(vel - F_kd.subs({k: 0, d: -1 / cvel**2, a: u, b: v})) == 0)
RESULTS['A'] = {'generator': str(P_kd), 'unit_free': 'd/k^2 (k != 0); sign of d (k = 0); sign of k'}

# ---------------------------------------------------------------------------
print("== B. s, the height: rung s of the fused law, in generator coordinates")
E_k = sp.expand(trunc((sp.exp(k * x) - 1) / k))
P_rung = generator_of_rung(E_k)
cr = [sp.factor(P_rung.coeff(x, n)) for n in range(0, N)]
print("   rung-s generator coefficients of the k law:")
for n, c in enumerate(cr):
    print("    c%d = %s" % (n, c))
check("B1 c1 = s*k: the rung-s law's own homothety reading is s times the base k", sp.simplify(cr[1] - s * k) == 0)
homog = all(sp.simplify(cr[n].subs(k, lam * k) - lam**n * cr[n]) == 0 for n in range(1, N))
check("B2 c_n = k^n f_n(s) for every n: k is the unit, s is unit-free", homog)

f = [sp.factor(sp.simplify(cr[n] / k**n)) for n in range(N)]
f2, f3, f4 = f[2], f[3], f[4]
check("B3 every f_n (n >= 2) carries the factor s(s - 1): the rung laws meet the k line only at s = 0 (addition) and s = 1 (the law itself)",
      all(sp.simplify(f[nn].subs(s, 0)) == 0 and sp.simplify(f[nn].subs(s, 1)) == 0 for nn in range(2, N)), "f2 = " + str(f2))
r3 = sp.solve(sp.numer(sp.together(f3)), s)
r4 = sp.solve(sp.numer(sp.together(f4)), s)
common = [r for r in r3 if any(sp.simplify(r - q) == 0 for q in r4)]
check("B4 f3 and f4 share only the roots s = 0 and s = 1, so for every other s the rung law has a cubic or quartic term: it is not a (k, d) law",
      sorted([sp.nsimplify(c) for c in common]) == [0, 1], "f3 roots %s, f4 roots %s" % (r3, r4))
RESULTS['B'] = {'f': [str(t) for t in f]}

# s = 2 closed form: P_2(x) = (1 + x)(1 + ln(1 + x)) at k = 1
P2 = sp.expand(trunc((1 + x) * (1 + sp.log(1 + x)), N - 1))
check("B5 at k = 1, s = 2 the series matches the closed form (1 + x)(1 + ln(1 + x))", sp.expand(P_rung.subs({k: 1, s: 2}) - P2) == 0)

# rank test: rung s of the plane law (k, d): Jacobian of (c1, c2, c3, c4) in (k, d, s)
# exponential of the plane law: G with G' = P(G), G(0) = 0, as a series
g = sp.symbols('g1:8')
Gser = x + sum(g[i] * x**(i + 2) for i in range(N - 1))
eqs = sp.expand(sp.diff(Gser, x) - (1 + k * Gser + d * Gser**2))
sol = sp.solve([eqs.coeff(x, n) for n in range(1, N)], g[:N - 1], dict=True)[0]
G_kd = sp.expand(Gser.subs(sol))
P_rung_kd = generator_of_rung(G_kd)
coefs = [P_rung_kd.coeff(x, n) for n in range(1, 5)]
Jr = sp.Matrix([[sp.diff(c, v) for v in (k, d, s)] for c in coefs])
ranks = []
for _ in range(5):
    pt = {k: sp.Rational(random.randint(1, 9), random.randint(2, 7)), d: sp.Rational(random.randint(-9, 9), random.randint(2, 7)),
          s: sp.Rational(random.randint(1, 30), 7)}
    ranks.append(Jr.subs(pt).rank())
check("B6 rank test: (c1..c4) of rung s of the plane law have rank 3 in (k, d, s) at 5 random points: s is a third independent number",
      all(r == 3 for r in ranks), str(ranks))
check("B7 rung 1 of the plane law is the plane law itself (generator 1 + kx + dx^2)", sp.expand(P_rung_kd.subs(s, 1) - (1 + k * x + d * x**2)) == 0)

# polish: f3 at high precision on a grid, never zero away from 0 and 1
f3n = sp.lambdify(s, f3, 'mpmath')
f4n = sp.lambdify(s, f4, 'mpmath')
grid = [mp.mpf(i) / 40 for i in range(-120, 121)]
worst = min(max(abs(f3n(t)), abs(f4n(t))) for t in grid if abs(t) > 0.02 and abs(t - 1) > 0.02)
check("B8 polish (30 digits): max(|f3|, |f4|) stays away from 0 on s in [-3, 3] outside s = 0, 1", worst > mp.mpf('1e-6'), mp.nstr(worst, 6))

# ---------------------------------------------------------------------------
print("== C. p, the roundness: the same height, applied across the axes")
ok = True
for _ in range(200):
    A, B = np.random.rand(2) * 5 + 0.01
    sh = np.random.uniform(-2, 3)
    # generator of the Pythagorean sum: G(y) = y^(1/2); G^sh(y) = y^(2^-sh); rung sh = G^sh(G^-sh A + G^-sh B)
    rung = (A**(2**sh) + B**(2**sh))**(2**-sh)
    lp = (A**(2**sh) + B**(2**sh))**(1 / 2**sh)
    ok &= abs(rung - lp) < 1e-9 * max(1, lp)
check("C1 rung s of the Pythagorean sum (a^2 + b^2)^(1/2) is the l_p sum with p = 2^s (200 random cases)", ok)


def conic_residual(pts):
    M = np.array([[X * X, X * Y, Y * Y, X, Y, 1.0] for X, Y in pts])
    return np.linalg.svd(M, compute_uv=False)[-1] / np.linalg.svd(M, compute_uv=False)[0]


th = np.linspace(0.1, 6.2, 40)
res = {}
for pp in (1.5, 2.0, 3.0, 4.0):
    c, sn = np.cos(th), np.sin(th)
    nrm = (np.abs(c)**pp + np.abs(sn)**pp)**(1 / pp)
    res[pp] = conic_residual(list(zip(c / nrm, sn / nrm)))
check("C2 the l_p unit curve lies on a conic only at p = 2 (smallest singular value of the conic fit)",
      res[2.0] < 1e-12 and all(res[q] > 1e-4 for q in (1.5, 3.0, 4.0)), str({q: float('%.2e' % r) for q, r in res.items()}))
KD = [(0.7, 0.3), (-1.0, 0.0), (0.0, 1.0), (1.2, -0.5)]
okc = True
for kk, dd in KD:
    # unit curve of the law's norm a^2 + k ab + d b^2 = 1 (Coleman-McGuinness norm), a conic by construction
    pts = []
    for t_ in np.linspace(0, 2 * np.pi, 30, endpoint=False):
        cA, sA = np.cos(t_), np.sin(t_)
        q = cA * cA + kk * cA * sA + dd * sA * sA
        if q > 1e-6:
            pts.append((cA / math.sqrt(q), sA / math.sqrt(q)))
    okc &= conic_residual(pts) < 1e-10
check("C3 every (k, d) law's unit curve is a conic, so p != 2 shapes are out of the (k, d) family's reach", okc)
ok = all(abs(((3 * 2.0)**pp + (3 * 5.0)**pp)**(1 / pp) - 3 * (2.0**pp + 5.0**pp)**(1 / pp)) < 1e-9 for pp in (0.5, 1, 2, 3.7))
check("C4 the l_p norm is homogeneous, so no unit change moves p: p is unit-free", ok)


# rank test: object with height s along each axis (rung s of the k = 1 law) and roundness p across axes
Gs_sym = iterate_coeffs(sp.expand(trunc(sp.exp(x) - 1, 7)), 7)
Gs_f = sp.lambdify((x, s), Gs_sym, 'math')
Gms_f = sp.lambdify((x, s), Gs_sym.subs(s, -s), 'math')


def rung_num(A, B, sv):
    return Gs_f(Gms_f(A, sv) + Gms_f(B, sv), sv)


def phi2(sv, pv, a1, a2, b1, b2):
    u1, u2 = rung_num(a1, b1, sv), rung_num(a2, b2, sv)
    return (abs(u1)**pv + abs(u2)**pv)**(1 / pv)


rk = []
for _ in range(5):
    sv, pv = random.uniform(0.2, 1.8), random.uniform(1.2, 4)
    rows = []
    for _ in range(4):
        args = [random.uniform(0.02, 0.12) for _ in range(4)]
        h = 1e-5
        rows.append([(phi2(sv + h, pv, *args) - phi2(sv - h, pv, *args)) / (2 * h),
                     (phi2(sv, pv + h, *args) - phi2(sv, pv - h, *args)) / (2 * h)])
    sv_ = np.linalg.svd(np.array(rows), compute_uv=False)
    rk.append(sv_[-1] / sv_[0])
check("C5 rank test: a shape with height s along each axis and roundness p across them has rank 2 in (s, p): both are needed",
      all(r > 1e-4 for r in rk), str([float('%.1e' % r) for r in rk]))

# ---------------------------------------------------------------------------
print("== D. t, the odds: the input of the k = 1 law, one per element")
t1, t2, t3 = sp.symbols('t1 t2 t3', positive=True)


def circk(kk, A, B):
    return A + B + kk * A * B


lhs = circk(k, circk(k, t1, t2), t3)
rhs = circk(1, circk(1, k * t1, k * t2), k * t3) / k
check("D1 three odds fused at dial k = (1/k) x the k = 1 fusion of the odds k*t (k is a common odds scale)", sp.simplify(lhs - rhs) == 0)

tt = sp.symbols('t1:5', positive=True)
n = 4
subsets = [S for r_ in range(0, n + 1) for S in itertools.combinations(range(n), r_)]
w = [k**(len(S) - 1) * sp.prod([tt[i] for i in S]) if S else sp.Integer(1) / k for S in subsets]
Z = sum(w)
Pd = [wi / Z for wi in w]
vars_ = [k] + list(tt)
pt = {k: sp.Rational(3, 2), tt[0]: sp.Rational(1, 3), tt[1]: 2, tt[2]: sp.Rational(5, 7), tt[3]: sp.Rational(4, 3)}
Jw = sp.Matrix([[sp.diff(sp.log(wi), vv) for vv in vars_] for wi in w[1:]]).subs(pt)
JP = sp.Matrix([[sp.diff(pi, vv) for vv in vars_] for pi in Pd]).subs(pt)
check("D2 rank test: the corner weights k^(|S|-1) prod t see k (rank 5 = 4 odds + k), but the chance of each subset does not (rank 4): there k is absorbed into the odds",
      Jw.rank() == 5 and JP.rank() == 4, "weights %d, chances %d" % (Jw.rank(), JP.rank()))

tz = sp.symbols('z', positive=True)
logZ = -sp.log(1 - k * tz) / k
check("D3 per-element factor (1 - k z)^(-1/k): log Z = (1/k) f(k z), so only the sign of k survives a unit change (as dial-survey found)",
      sp.simplify(logZ.subs({k: lam * k, tz: tz / lam}) - logZ / lam) == 0)
m = sp.symbols('m', positive=True, integer=True)
check("D4 k = -1/m is literally m one-copy slots, each at odds z/m: (1 + z/m)^m",
      sp.simplify((1 - (-1 / m) * tz)**(-1 / (-1 / m)) - (1 + tz / m)**m) == 0)
P1 = t1 / (1 + t1); P2_ = t2 / (1 + t2)
U = P1 + P2_ - P1 * P2_
check("D5 'in A or in B' for independent elements has odds t1 o1 t2 and the complement is 1/t (repo's odds_union)",
      sp.simplify(U / (1 - U) - circk(1, t1, t2)) == 0 and sp.simplify((1 - P1) / P1 - 1 / t1) == 0)

# ---------------------------------------------------------------------------
print("== E. r, the similarity ratio: odds across pieces plus one zoom")


def moran_D(rs):
    f_ = lambda D: sum(ri**D for ri in rs) - 1
    lo, hi = mp.mpf('1e-9'), mp.mpf(60)
    return mp.findroot(f_, (lo, hi), solver='bisect') if f_(lo) * f_(hi) < 0 else None


ok = True; okJ = True
for _ in range(30):
    nn = random.randint(2, 6)
    rs = [mp.mpf(random.uniform(0.05, 0.7)) for _ in range(nn)]
    D = moran_D(rs)
    if D is None:
        continue
    wv = [ri**D for ri in rs]
    back = [wi**(1 / D) for wi in wv]
    ok &= abs(sum(wv) - 1) < mp.mpf('1e-25') and max(abs(p_ - q_) for p_, q_ in zip(back, rs)) < mp.mpf('1e-25')
check("E1 polish (30 digits): n ratios <-> (dimension D, piece weights w_i = r_i^D on the simplex) is a bijection (30 random fractals)", ok)

# Jacobian of r -> (D, w_1..w_{n-1}) is full rank
for _ in range(10):
    nn = random.randint(2, 5)
    rs = np.random.uniform(0.1, 0.6, nn)

    def mapping(rv):
        D = float(moran_D([mp.mpf(z) for z in rv]))
        return np.array([D] + [z**D for z in rv[:-1]])
    h = 1e-7
    Jn = np.array([(mapping(rs + h * np.eye(nn)[i]) - mapping(rs - h * np.eye(nn)[i])) / (2 * h) for i in range(nn)]).T
    sv_ = np.linalg.svd(Jn, compute_uv=False)
    okJ &= sv_[-1] / sv_[0] > 1e-6
check("E2 rank test: r -> (D, n-1 free weights) has full rank n (10 random cases): r carries n numbers, D plus n-1 relative odds", okJ)

rs = [mp.mpf('0.3')] * 4
D = moran_D(rs)
check("E3 N equal pieces weigh exactly 1/N however r is turned, so the uniform part of r is D alone (ratio-dial-survey C3)",
      all(abs(ri**D - mp.mpf(1) / 4) < mp.mpf('1e-25') for ri in rs))
rs2 = [mp.mpf('0.3'), mp.mpf('0.3'), mp.mpf('0.3'), mp.mpf('0.45')]
D2 = moran_D(rs2)
wv2 = [ri**D2 for ri in rs2]
check("E4 turning one piece moves the weights around the simplex (relative odds change) and moves D", abs(wv2[3] / wv2[0] - 1) > 0.1 and abs(D2 - D) > 0.05,
      "D %s -> %s, odds of piece 4 vs 1: %s" % (mp.nstr(D, 5), mp.nstr(D2, 5), mp.nstr(wv2[3] / wv2[0], 5)))
check("E5 ratios multiply as the k law's multipliers: (1 + k a)(1 + k b) = 1 + k (a o_k b) (repo's circK lift)",
      sp.expand((1 + k * a) * (1 + k * b) - (1 + k * circk(k, a, b))) == 0)

# ---------------------------------------------------------------------------
print("== F. q, the order dial")
qv = sp.symbols('q')


def circ_power(kk, A, nn):
    out = 0
    for _ in range(nn):
        out = circk(kk, out, A)
    return sp.expand(out)


check("F1 counting shadow: the q-integer [n]_q is the n-fold o_k power of 1 at k = q - 1 (n = 1..8)",
      all(sp.simplify(circ_power(qv - 1, 1, nn) - (qv**nn - 1) / (qv - 1)) == 0 for nn in range(1, 9)))


def gauss_from_words(nn, j):
    tot = 0
    for pos in itertools.combinations(range(nn), j):
        word = ['y' if i in pos else 'x' for i in range(nn)]
        inv = sum(1 for i1 in range(nn) for i2 in range(i1 + 1, nn) if word[i1] == 'y' and word[i2] == 'x')
        tot += qv**inv
    return sp.expand(tot)


def gauss_formula(nn, j):
    num = sp.prod([circ_power(qv - 1, 1, nn - i) for i in range(j)])
    den = sp.prod([circ_power(qv - 1, 1, i + 1) for i in range(j)])
    return sp.expand(sp.cancel(num / den))


check("F2 (x + y)^n with yx = q xy has Gaussian binomial coefficients, built from k = q - 1 powers (n <= 7)",
      all(sp.expand(gauss_from_words(nn, j) - gauss_formula(nn, j)) == 0 for nn in range(1, 8) for j in range(nn + 1)))


# moves built from the (k = 0, d) law: translations on the surface of curvature d
def gens(dv):
    Kx = np.array([[0, 0, 1], [0, 0, 0], [-dv, 0, 0]], float)
    Ky = np.array([[0, 0, 0], [0, 0, 1], [0, -dv, 0]], float)
    Jz = np.array([[0, -1, 0], [1, 0, 0], [0, 0, 0]], float)
    return Kx, Ky, Jz


def expm(M):
    from scipy.linalg import expm as _e
    return _e(M)


okF = True; ratios = {}
for dv in (-1.0, -0.4, 0.0, 0.3, 1.0):
    Kx, Ky, Jz = gens(dv)
    br = Kx @ Ky - Ky @ Kx
    okF &= np.allclose(br, dv * Jz)
    e = 1e-3
    C = expm(e * Kx) @ expm(e * Ky) @ expm(-e * Kx) @ expm(-e * Ky)
    ang = math.atan2(C[1, 0], C[0, 0])
    ratios[dv] = ang / (e * e)
    okF &= abs(ang / (e * e) - dv) < 1e-3
check("F3 moves built from the plane law: [Kx, Ky] = d Jz, and two small moves fail to commute by a turn = d x area (d = -1, -0.4, 0, 0.3, 1)",
      okF, str({kk: round(vv, 5) for kk, vv in ratios.items()}))

# finite triangle: holonomy = d x area (Gauss-Bonnet), sphere d = 1
A_ = np.array([0, 0, 1.0]); B_ = np.array([math.sin(0.7), 0, math.cos(0.7)]); C_ = np.array([0, math.sin(0.9), math.cos(0.9)])


def angle_at(P, Q, R):
    t1_ = Q - np.dot(Q, P) * P; t2_ = R - np.dot(R, P) * P
    return math.acos(np.dot(t1_, t2_) / np.linalg.norm(t1_) / np.linalg.norm(t2_))


excess = angle_at(A_, B_, C_) + angle_at(B_, C_, A_) + angle_at(C_, A_, B_) - math.pi
# area by L'Huilier
a_ = math.acos(np.dot(B_, C_)); b_ = math.acos(np.dot(A_, C_)); c_ = math.acos(np.dot(A_, B_)); s_ = (a_ + b_ + c_) / 2
area = 4 * math.atan(math.sqrt(math.tan(s_ / 2) * math.tan((s_ - a_) / 2) * math.tan((s_ - b_) / 2) * math.tan((s_ - c_) / 2)))
check("F4 a finite loop: the turn equals curvature x area exactly (Gauss-Bonnet, sphere d = 1)", abs(excess - area) < 1e-12,
      "excess %.12f, area %.12f" % (excess, area))

# Weyl/Moyal: hbar is a unit for continuous moves; clock-shift q is not removable by rescaling
nn = 5; qq = np.exp(2j * np.pi * 2 / nn)
Uc = np.diag([qq**j for j in range(nn)]); Vc = np.roll(np.eye(nn), 1, axis=0)
qrel = (Vc @ Uc @ np.linalg.inv(Uc @ Vc))[0, 0]
okW = np.allclose(Vc @ Uc, qrel * Uc @ Vc) and abs(abs(qrel) - 1) < 1e-12 and abs(qrel - 1) > 0.1
al, be = 1.7 * np.exp(0.3j), 0.6 * np.exp(-1.1j)
okW &= np.allclose((be * Vc) @ (al * Uc), qrel * (al * Uc) @ (be * Vc))
check("F5 clock and shift matrices obey VU = qUV, and rescaling U, V leaves q unchanged: q is unit-free", okW)
hb, sx, sy = sp.symbols('hbar s_x s_y', positive=True)
check("F6 for continuous moves the Weyl phase e^(i hbar s t) is unchanged by (hbar, s) -> (1, hbar s): hbar's size is a unit (Stone-von Neumann)",
      sp.simplify(sp.exp(sp.I * hb * sx * sy) - sp.exp(sp.I * 1 * (hb * sx) * sy)) == 0)

# ---------------------------------------------------------------------------
print("== G. Summary")
verdict = {
    'k': 'fundamental (first generator coefficient); its size is a unit choice, only sign survives a unit change',
    'd': 'fundamental (second generator coefficient); with k: d/k^2 is the unit-free number; curvature -Delta/4 is a reading',
    's': 'fundamental number (rank 3 with k, d); same operation as p, applied along an axis',
    'p': 'fundamental number (rank 2 with s); it is the height across axes, p = 2^s on the Pythagorean sum',
    't': 'fundamental as an input (one per element), not a law dial: odds compose by the k = 1 law; k is absorbed into t in subset chances',
    'r': 'partly derived: n ratios = D plus n-1 relative odds; only the overall zoom (equivalently D) is its own',
    'q': 'partly derived: for moves built from the plane law the order turn is d x area; counts see q only as k = q - 1; new only when put in by hand (hbar is then a unit for continuous moves; on a lattice q is a real invariant)'}
RESULTS['verdict'] = verdict
RESULTS['passed'] = '%d of %d' % (sum(PASS), len(PASS))
print("%d of %d pass" % (sum(PASS), len(PASS)))
with open(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'results.json'), 'w') as fh:
    json.dump(RESULTS, fh, indent=1, default=str)
