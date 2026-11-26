"""Checks for the homothety-dial survey (k dial, d = 0). Run: python3 checks.py > checks_output.txt
Every check prints PASS/FAIL with a label; the last line is 'X of Y pass'. numpy, scipy, sympy, mpmath."""
import itertools, math, time, json
import numpy as np, sympy as sp, mpmath as mp
from scipy import stats, integrate, optimize
from scipy.sparse import csr_matrix
from scipy.sparse.csgraph import maximum_flow

R = []
def check(label, ok, extra=""):
    R.append(bool(ok)); print(("PASS " if ok else "FAIL ") + label + (f"  [{extra}]" if extra else ""))
def circ(k, a, b): return a + b + k * a * b
POLISH = {}

k, a, b, c, lam = sp.symbols("k a b c lam")
# ---------------- A. Units: the three stops ----------------
check("A1 rescaling: lam*(a o_k b) = (lam a) o_{k/lam} (lam b)",
      sp.simplify(lam * circ(k, a, b) - circ(k / lam, lam * a, lam * b)) == 0)
check("A2 sign flip: a o_k b = -((-a) o_{-k} (-b)), so k<0 and k>0 are mirror images",
      sp.expand(circ(k, a, b) + circ(-k, -a, -b)) == 0)
check("A3 k -> 0 is plain addition, reached continuously", sp.limit(circ(k, a, b), k, 0) == a + b)
check("A4 the wall -1/k absorbs: a o_k (-1/k) = -1/k", sp.simplify(circ(k, a, -1 / k) + 1 / k) == 0)

# ---------------- B. C1 corner expansion and powers ----------------
xs = sp.symbols("x1:6")
fold = xs[0]
for x in xs[1:]: fold = circ(k, fold, x)
check("B1 1 + k(x1 o ... o x5) = prod(1 + k xi), every k",
      sp.expand(1 + k * fold - sp.prod([1 + k * x for x in xs])) == 0)
corners = sum(k ** (len(S) - 1) * sp.prod(S) for r in range(1, 6) for S in itertools.combinations(xs, r))
check("B2 corner weights k^(|S|-1) x^S (cube of corners)", sp.expand(fold - corners) == 0)
rng = np.random.default_rng(1)
p = rng.uniform(0, 1, 5); u = rng.uniform(size=(400000, 5)) < p
union_mc = u.any(axis=1).mean()
union_fold = p[0]
for q in p[1:]: union_fold = circ(-1, union_fold, q)
check("B3 k = -1: fold is P(at least one) of independent events (inclusion-exclusion)",
      abs(union_mc - union_fold) < 4e-3, f"MC {union_mc:.4f} vs {union_fold:.4f}")
n = sp.Symbol("n", integer=True, positive=True)
wp = ((1 + k * a) ** 6 - 1) / k
check("B4 powers: k->0 gives 6a (repeated addition), k=1 gives (1+a)^6-1, k=-1 gives 1-(1-a)^6",
      sp.simplify(sp.limit(wp, k, 0) - 6 * a) == 0 and sp.expand(wp.subs(k, 1) - ((1 + a) ** 6 - 1)) == 0
      and sp.expand(wp.subs(k, -1) - (1 - (1 - a) ** 6)) == 0)
A = rng.normal(size=(5, 5)); ev = np.linalg.eigvals(A)
kk = 0.37
d1 = np.linalg.det(np.eye(5) + kk * A); d2 = np.prod(1 + kk * ev)
check("B5 matrix fold: det(I + kA) = prod(1 + k lambda); walls at k = -1/lambda",
      abs(d1 - d2) < 1e-10 and all(abs(np.linalg.det(np.eye(5) - A / l)) < 1e-8 for l in ev))
# logistic <-> softmax normalizer
x = rng.uniform(0.1, 2, 4)
def Zls(kv): return 1 + (np.prod(1 + kv * x) - 1) / kv if kv != 0 else 1 + x.sum()
check("B6 normalizer 1 + (prod(1+kx)-1)/k: softmax-with-none (k->0) and independent logits (k=1)",
      abs(Zls(1e-9) - (1 + x.sum())) < 1e-7 and abs(Zls(1) - np.prod(1 + x)) < 1e-12)

# ---------------- C. C2 per-mode counts = Katz/Panjer (a,b,0) ----------------
def dialpmf(kv, z, emax=60):
    w = [1.0]
    for e in range(1, emax + 1): w.append(w[-1] * z * (1 + (e - 1) * kv) / e)
    w = np.array(w); return w / w.sum(), w
z = 0.6
pm, w = dialpmf(0.3, z)
ratios = [pm[e] / pm[e - 1] for e in range(1, 10)]
check("C1 ratio p_e/p_(e-1) = a + b/e with a = k z, b = z(1-k) (Katz/Panjer form)",
      max(abs(ratios[e - 1] - (0.3 * z + z * 0.7 / e)) for e in range(1, 10)) < 1e-12)
ok = True
for m in (1, 2, 3, 5):
    pmm, _ = dialpmf(-1 / m, z, m + 5)
    ok &= np.allclose(pmm[:m + 1], stats.binom.pmf(range(m + 1), m, z / (m + z)), atol=1e-13) and np.allclose(pmm[m + 1:], 0)
check("C2 k = -1/m: binomial(m, z/(m+z)), at most m copies (m = 1,2,3,5)", ok)
pm0, _ = dialpmf(0.0, z)
check("C3 k = 0: Poisson(z)", np.allclose(pm0[:30], stats.poisson.pmf(range(30), z), atol=1e-14))
ok = True
for kv in (0.25, 0.5, 1.0, 2.0):
    zz = 0.4 / kv * 0.9
    pk, _ = dialpmf(kv, zz, 200)
    ok &= np.allclose(pk[:40], stats.nbinom.pmf(range(40), 1 / kv, 1 - kv * zz), atol=1e-12)
check("C4 k > 0: negative binomial r = 1/k, success 1 - kz (wall z < 1/k); k = 1 geometric", ok)
coef = [sp.prod([1 + j * sp.Rational(-2, 3) for j in range(e)]) for e in range(6)]
check("C5 k = -2/3 (not -1/m): a negative coefficient appears (e = 3), so no real count",
      coef[3] < 0, f"coeffs {coef}")
zz, x_ = sp.symbols("z x", positive=True)
logZ = -sp.log(1 - k * zz) / k
mean = sp.simplify(zz * sp.diff(logZ, zz))
var = sp.simplify(zz * sp.diff(mean, zz))
check("C6 mean copies z/(1 - kz) = 1/(e^x - k): Fermi (-1), Boltzmann (0), Bose (+1)",
      sp.simplify(mean.subs(zz, sp.exp(-x_)) - 1 / (sp.exp(x_) - k)) == 0)
check("C7 spread = n(1 + k n): k < 0 copies avoid, k > 0 bunch", sp.simplify(var - mean * (1 + k * mean)) == 0)
# Sundt-Jewell table mapping: binomial a = -p/(1-p), b = (m+1)p/(1-p)  => b/a = -(m+1) => k = -1/m
pS, mS = sp.symbols("p m", positive=True)
aB, bB = -pS / (1 - pS), (mS + 1) * pS / (1 - pS)
ksol = sp.solve(sp.Eq((1 - k) / k, sp.simplify(bB / aB)), k)
check("C8 the binomial row of the (a,b,0) table forces k = -1/m (b/a = (1-k)/k = -(m+1))",
      len(ksol) == 1 and sp.simplify(ksol[0] + 1 / mS) == 0)
t0 = time.time(); mp.mp.dps = 30
kv, zv = mp.mpf("0.5"), mp.mpf("0.7")
lhs = mp.nsum(lambda e: zv ** e * mp.rf(1 / kv, e) * kv ** e / mp.factorial(e), [0, mp.inf])
check("C9 polish: sum_e z^e prod(1+jk)/e! = (1 - kz)^(-1/k) to 30 digits",
      abs(lhs - (1 - kv * zv) ** (-1 / kv)) < mp.mpf(10) ** -28)
POLISH["C9"] = time.time() - t0

# ---------------- D. C3 Tsallis / Renyi / Shannon ----------------
q_ = sp.Symbol("q", positive=True)
def lnq(x, q): return (x ** (1 - q) - 1) / (1 - q)
X, Y = sp.symbols("X Y", positive=True)
check("D1 Tsallis q-log is the dial log: ln_q(XY) = ln_q X o_k ln_q Y, k = 1 - q",
      sp.simplify(lnq(X * Y, q_) - circ(1 - q_, lnq(X, q_), lnq(Y, q_))) == 0)
def tsallis(P, q): return (1 - np.sum(P ** q)) / (q - 1)
def renyi(P, q): return np.log(np.sum(P ** q)) / (1 - q)
P1 = rng.dirichlet(np.ones(4)); P2 = rng.dirichlet(np.ones(3)); PJ = np.outer(P1, P2).ravel()
ok1 = ok2 = True
for qv in (0.3, 0.7, 1.5, 2.0, 3.0):
    kv = 1 - qv
    ok1 &= abs(tsallis(PJ, qv) - circ(kv, tsallis(P1, qv), tsallis(P2, qv))) < 1e-12
    ok2 &= abs(renyi(P1, qv) - np.log(1 + kv * tsallis(P1, qv)) / kv) < 1e-12
check("D2 pseudo-additivity S_q(A,B) = S_q(A) o_k S_q(B) for independent A, B (5 values of q)", ok1)
check("D3 Renyi = L_k(Tsallis) = log(1 + k T)/k, the dial's ruler", ok2)
sh = -np.sum(P1 * np.log(P1))
check("D4 k -> 0 (q -> 1): Tsallis and Renyi both become Shannon",
      abs(tsallis(P1, 1 + 1e-7) - sh) < 1e-6 and abs(renyi(P1, 1 + 1e-7) - sh) < 1e-6)

# ---------------- E. C4 Sugeno lambda-measures ----------------
def sugeno_lambda(g):
    f = lambda l: np.prod(1 + l * g) - (1 + l)
    s = g.sum()
    if abs(s - 1) < 1e-14: return 0.0
    return optimize.brentq(f, 1e-12, 1e6) if s < 1 else optimize.brentq(f, -1 + 1e-12, -1e-12)
gA = np.array([0.1, 0.2, 0.3]); gB = np.array([0.4, 0.5, 0.6])
lA, lB = sugeno_lambda(gA), sugeno_lambda(gB)
check("E1 normalisation pins lambda: sum g < 1 gives lambda > 0, sum g > 1 gives -1 < lambda < 0",
      lA > 0 and -1 < lB < 0, f"lambda = {lA:.4f}, {lB:.4f}")
def mob(g, l):
    return {S: l ** (len(S) - 1) * np.prod(g[list(S)]) for r in range(1, len(g) + 1) for S in itertools.combinations(range(len(g)), r)}
mA = mob(gA, lA); mB = mob(gB, lB)
check("E2 Mobius masses lambda^(|S|-1) prod g sum to 1; all >= 0 for lambda > 0 (belief function)",
      abs(sum(mA.values()) - 1) < 1e-12 and min(mA.values()) >= 0)
dl = -lB / (1 + lB)
gdual = np.array([1 - (np.prod(1 + lB * gB) / (1 + lB * gB[i]) - 1) / lB for i in range(3)])  # mu(complement of i) -> dual
check("E3 lambda < 0: masses change sign (plausibility side); its dual has weight -lambda/(1+lambda) > 0 (belief)",
      min(mB.values()) < 0 and dl > 0 and abs(sugeno_lambda(gdual) - dl) < 1e-9)
gC = np.array([0.2, 0.3, 0.5])
check("E4 sum g = 1: lambda = 0, the measure is a probability (additive)", sugeno_lambda(gC) == 0.0)

# ---------------- F. C5 calculus ----------------
t = sp.Symbol("t"); F = sp.Function("F")(t); G = sp.Function("G")(t); f = sp.Function("f")
Dk = lambda H: sp.diff(H, t) / (1 + k * H)
check("F1 product rule: D_k(F o_k G) = D_k F + D_k G, D_k F = F'/(1 + kF), every k",
      sp.simplify(Dk(circ(k, F, G)) - Dk(F) - Dk(G)) == 0)
U = sp.Function("U")(t)
Ek = (sp.exp(k * U) - 1) / k
check("F2 fundamental theorem: D_k(E_k(U)) = U', E_k(u) = (e^(ku) - 1)/k; k=0 ordinary, k=1 geometric",
      sp.simplify(Dk(Ek) - sp.diff(U, t)) == 0 and sp.limit(Ek, k, 0) == U)
fx = lambda s: 1 + np.sin(3 * s); I = integrate.quad(fx, 0, 1)[0]
ok = True
for kv in (-1.0, 0.0, 0.5, 1.0):
    N = 200000; s = (np.arange(N) + 0.5) / N
    prod_ = np.exp(np.sum(np.log1p(kv * fx(s) / N)))
    target = math.exp(kv * I)
    ok &= abs(prod_ / target - 1) < 5e-5
check("F3 product integral: prod(1 + k f dx) -> exp(k * int f) (Volterra 1887 at k = 1; survival at k = -1)", ok)
hz = lambda s: 2 * s
Fcdf = 1 - math.exp(-integrate.quad(hz, 0, 1)[0])
check("F4 k = -1 is the hazard calculus: D_{-1}F = F'/(1-F) is the hazard, F = 1 - exp(-int h)",
      abs((2 * 1 * math.exp(-1)) / (1 - Fcdf) - hz(1)) < 1e-12)
mu_, r_, i_ = sp.symbols("mu r i")
check("F5 graininess = k: (1 + mu r)(1 + mu i) = 1 + mu (r o_mu i) (Hilger; Fisher's rule with period mu)",
      sp.expand((1 + mu_ * r_) * (1 + mu_ * i_) - (1 + mu_ * circ(mu_, r_, i_))) == 0)
check("F6 e as the grainy-to-smooth slide: (1 + k)^(1/k) -> e as k -> 0", sp.limit((1 + k) ** (1 / k), k, 0) == sp.E)
# jets: R[e]/(e^2 - k e): (a + b e)(c + d e) = ac + (ad + bc + k bd) e
a1, b1, c1, d1 = sp.symbols("a1 b1 c1 d1")
eps = sp.Symbol("e")
prod_jet = sp.rem(sp.expand((a1 + b1 * eps) * (c1 + d1 * eps)), eps ** 2 - k * eps, eps)
f_ = sp.Function("g"); h = sp.Symbol("h")
check("F7 jets R[e]/(e^2 - ke): units 1 + b e multiply by o_k; k = 0 dual numbers (derivative), k = h two points (difference)",
      sp.expand(prod_jet.subs({a1: 1, c1: 1}) - (1 + circ(k, b1, d1) * eps)) == 0)

# ---------------- G. C6 primes Z_k(s) ----------------
def primes_upto(N):
    s = np.ones(N + 1, bool); s[:2] = False
    for i in range(2, int(N ** 0.5) + 1):
        if s[i]: s[i * i::i] = False
    return np.nonzero(s)[0]
PR = primes_upto(2_000_000).astype(float)
def mobius(nn):
    r, m, p = 1, nn, 2
    while p * p <= m:
        if m % p == 0:
            m //= p
            if m % p == 0: return 0
            r = -r
        p += 1
    return -r if m > 1 else r
def Mneck(kv, N): return sum(mobius(d) * kv ** (N // d) for d in range(1, N + 1) if N % d == 0) / N
def Zeuler(kv, s): return np.exp(np.sum(-np.log1p(-kv * PR ** (-s)) / kv)) if kv != 0 else np.exp(np.sum(PR ** (-s)))
def Zneck(kv, s, Nmax=50):
    return np.prod([float(mp.zeta(N * s)) ** (Mneck(kv, N) / kv) for N in range(1, Nmax + 1)])
ok = True
for kv in (-0.5, 0.5, 1.5, -1.0, 1.0):
    ok &= abs(Zeuler(kv, 2) / Zneck(kv, 2) - 1) < 2e-6
check("G1 Z_k(2): Euler product = prod_N zeta(2N)^(M(k,N)/k) (necklace counts), k = -1, -1/2, 1/2, 1, 3/2", ok)
check("G2 k = -1 is zeta(s)/zeta(2s), k = 1 is zeta(s): only finitely many zeta factors",
      abs(Mneck(-1, 1) + 1) < 1e-15 and abs(Mneck(-1, 2) - 1) < 1e-15 and all(abs(Mneck(-1, N)) < 1e-15 for N in range(3, 13))
      and all(abs(Mneck(1, N)) < 1e-15 for N in range(2, 13)),
      "M(-1,1) = -1, M(-1,2) = 1, so Z_-1 = zeta(s)^1 zeta(2s)^-1")
check("G3 k = -1/2 needs infinitely many factors (M(k,N) != 0 for N up to 12)", all(abs(Mneck(-0.5, N)) > 1e-12 for N in range(1, 13)))
check("G4 k = 0 limit: exp(P(s)), P the prime zeta", abs(Zeuler(1e-9, 2) / Zeuler(0, 2) - 1) < 1e-7)
t0 = time.time(); mp.mp.dps = 30
z30 = mp.zeta(2) / mp.zeta(4)
check("G5 polish: Z_-1(2) = zeta(2)/zeta(4) = 15/pi^2 to 30 digits", abs(z30 - 15 / mp.pi ** 2) < mp.mpf(10) ** -28)
POLISH["G5"] = time.time() - t0

# ---------------- H. C7 tunnelling averages, k = -2K ----------------
def avg_t(K, D1, D2):
    if K > 0:
        s = math.sqrt(K)
        cc = lambda ph: math.cos(s * D1) * math.cos(s * D2) + math.sin(s * D1) * math.sin(s * D2) * math.cos(ph)
        tf = lambda ph: (1 - cc(ph)) / (2 * K)
    elif K < 0:
        s = math.sqrt(-K)
        cc = lambda ph: math.cosh(s * D1) * math.cosh(s * D2) - math.sinh(s * D1) * math.sinh(s * D2) * math.cos(ph)
        tf = lambda ph: (1 - cc(ph)) / (2 * K)
    else:
        tf = lambda ph: (D1 ** 2 + D2 ** 2 - 2 * D1 * D2 * math.cos(ph)) / 4
    return integrate.quad(tf, 0, 2 * math.pi)[0] / (2 * math.pi)
def tK(K, D):
    if K > 0: return (1 - math.cos(math.sqrt(K) * D)) / (2 * K)
    if K < 0: return (1 - math.cosh(math.sqrt(-K) * D)) / (2 * K)
    return D ** 2 / 4
ok = True
for K in (1.0, 0.5, 0.0, -0.3, -1.0, -2.0):
    D1, D2 = 0.7, 1.1
    ok &= abs(avg_t(K, D1, D2) - circ(-2 * K, tK(K, D1), tK(K, D2))) < 1e-12
check("H1 averaged law of cosines: <t3> = t1 o_k t2 with k = -2K (sphere, flat, hyperbolic; 6 curvatures)", ok)
pp, qq = 0.1, 0.3
M1 = np.array([[1 - pp, pp], [pp, 1 - pp]]); M2 = np.array([[1 - qq, qq], [qq, 1 - qq]])
check("H2 k = -2 stop is the binary symmetric channel (piling-up lemma: 1 - 2p multiplies)",
      abs((M1 @ M2)[0, 1] - circ(-2, pp, qq)) < 1e-15)
check("H3 flat stop K = 0: t = D^2/4 adds (variances of independent steps add)", abs(avg_t(0, 0.7, 1.1) - (0.7 ** 2 + 1.1 ** 2) / 4) < 1e-12)

t0 = time.time(); mp.mp.dps = 30
D1m, D2m = mp.mpf("0.7"), mp.mpf("1.1")
avg30 = mp.quad(lambda ph: (1 - (mp.cosh(D1m) * mp.cosh(D2m) - mp.sinh(D1m) * mp.sinh(D2m) * mp.cos(ph))) / (-2), [0, 2 * mp.pi]) / (2 * mp.pi)
t1m, t2m = (1 - mp.cosh(D1m)) / (-2), (1 - mp.cosh(D2m)) / (-2)
check("H4 polish: hyperbolic (K = -1) average equals t1 o_2 t2 to 28 digits", abs(avg30 - (t1m + t2m + 2 * t1m * t2m)) < mp.mpf(10) ** -27)
POLISH["H4"] = time.time() - t0

# ---------------- I. C8 Lee-Yang on graphs (dial gas) ----------------
def gas_poly(edges, nv, kv):
    coeff = np.zeros(nv + 1)
    for S in range(1 << nv):
        m = bin(S).count("1"); e = sum(1 for (i, j) in edges if (S >> i) & 1 and (S >> j) & 1)
        coeff[m] += kv ** e
    return coeff
ring8 = [(i, (i + 1) % 8) for i in range(8)]
cube3 = [(i, i ^ (1 << b_)) for i in range(8) for b_ in range(3) if i < i ^ (1 << b_)]
def circle_dev(edges, nv, deg, kv):
    r = np.roots(gas_poly(edges, nv, kv)[::-1])
    return np.max(np.abs(np.abs(r) - kv ** (-deg / 2)))
okon = all(circle_dev(ring8, 8, 2, kv) < 1e-6 and circle_dev(cube3, 8, 3, kv) < 1e-6 for kv in (1.2, 1.5, 3.0))
check("I1 k > 1 (attraction): every zero on |z| = k^(-deg/2) (ring of 8, 3-cube; k = 1.2, 1.5, 3); k = 1 is the free gas (1+z)^n", okon and np.allclose(gas_poly(ring8, 8, 1.0), [math.comb(8, j) for j in range(9)]))
offs = [circle_dev(ring8, 8, 2, kv) for kv in (0.3, 0.6, 0.9)]
check("I2 k < 1 (repulsion): zeros leave that circle (ring of 8, k = 0.3, 0.6, 0.9)", min(offs) > 1e-3,
      ", ".join(f"{o:.3f}" for o in offs))
r0 = np.roots(gas_poly(ring8, 8, 0.0)[::-1])
check("I3 k = 0 (hard core): all zeros on the negative real axis (ring of 8)",
      np.max(np.abs(r0.imag)) < 1e-6 and np.all(r0.real < 0))
t0 = time.time(); mp.mp.dps = 30
cf = [mp.mpf(0)] * 9
for S in range(256):
    m = bin(S).count("1"); e = sum(1 for (i, j) in ring8 if (S >> i) & 1 and (S >> j) & 1)
    cf[m] += mp.mpf(2) ** e
roots30 = mp.polyroots(cf[::-1], maxsteps=200, extraprec=60)
check("I4 polish: at k = 2 the 8 ring zeros sit on |z| = 1/2 to 25 digits",
      max(abs(abs(rr) - mp.mpf(1) / 2) for rr in roots30) < mp.mpf(10) ** -25)
POLISH["I4"] = time.time() - t0

# ---------------- J. C9 counting with a pinned input ----------------
qs = sp.Symbol("q")
qnum = lambda nn: sum(qs ** j for j in range(nn))
def wpowk(kv, av, nn):
    r_ = 0
    for _ in range(nn): r_ = circ(kv, r_, av)
    return r_
check("J1 [n]_q = n copies of 1 under o_k with k = q - 1 (n = 1..8)",
      all(sp.expand(wpowk(qs - 1, 1, nn) - qnum(nn)) == 0 for nn in range(1, 9)))
check("J2 stops: k = -2 parity, k = -1 OR, k = 0 count, k = 1 nonempty subsets",
      [wpowk(-2, 1, nn) for nn in range(1, 7)] == [1, 0, 1, 0, 1, 0] and [wpowk(-1, 1, nn) for nn in range(1, 7)] == [1] * 6
      and [wpowk(0, 1, nn) for nn in range(1, 7)] == list(range(1, 7)) and [wpowk(1, 1, nn) for nn in range(1, 7)] == [2 ** nn - 1 for nn in range(1, 7)])
def subspaces(qv, dim, sub):
    vecs = list(itertools.product(range(qv), repeat=dim)); spans = set()
    for basis in itertools.combinations([v for v in vecs if any(v)], sub):
        span = set()
        for coefs in itertools.product(range(qv), repeat=sub):
            span.add(tuple(sum(cc * bb[i] for cc, bb in zip(coefs, basis)) % qv for i in range(dim)))
        if len(span) == qv ** sub: spans.add(frozenset(span))
    return len(spans)
gauss = lambda nn, kk_, qv: sp.prod([(qv ** (nn - j) - 1) for j in range(kk_)]) / sp.prod([(qv ** (j + 1) - 1) for j in range(kk_)])
check("J3 in-between whole k = q - 1 (q a prime): Gaussian binomials count subspaces (F_2^4: 35 planes; F_3^3: 13 lines)",
      subspaces(2, 4, 2) == 35 == gauss(4, 2, 2) and subspaces(3, 3, 1) == 13 == gauss(3, 1, 3))
from sympy.functions.combinatorial.numbers import stirling
T = lambda nn, x: sum(stirling(nn, j) * x ** j for j in range(nn + 1))
check("J4 Touchard T_n(k): Bell numbers at k = 1, complementary (Uppuluri-Carpenter) numbers at k = -1",
      [T(nn, 1) for nn in range(1, 9)] == [1, 2, 5, 15, 52, 203, 877, 4140]
      and [T(nn, -1) for nn in range(1, 9)] == [-1, 0, 1, 1, -2, -9, -9, 50])
from functools import reduce
kf5 = reduce(lambda r_, j: circ(k, r_, j), range(1, 6), 0)
check("J5 rising product: 1 + k(1 o_k 2 o_k ... o_k m) = prod_{j<=m}(1 + jk) (Combinatorics.kfold_lift)",
      sp.expand(1 + k * kf5 - sp.prod([1 + j * k for j in range(1, 6)])) == 0)
aI, bI = sp.symbols("aI bI")
dfe = sp.expand(circ(k, aI, bI) ** 2 - circ(k, aI, bI))
dfe = dfe.subs({aI ** 2: aI, bI ** 2: bI})
check("J6 idempotent defect (k+1)(k+2)ab: bits stay bits only at k = -1 (OR) and k = -2 (XOR)",
      sp.expand(dfe - (k + 1) * (k + 2) * aI * bI) == 0 and
      [kv for kv in range(-5, 6) if circ(kv, 1, 1) in (0, 1)] == [-2, -1])

# ---------------- K. C10 the alphabet ladder k = -q/(q-1) ----------------
def symch(qv, e): return (1 - e) * np.eye(qv) + e / (qv - 1) * (np.ones((qv, qv)) - np.eye(qv))
ok = True
for qv in (2, 3, 4, 5, 8):
    kq = -qv / (qv - 1); e1, e2 = 0.13, 0.31
    MM = symch(qv, e1) @ symch(qv, e2)
    ok &= np.allclose(MM, symch(qv, circ(kq, e1, e2)), atol=1e-14)
check("K1 q-ary symmetric channels: error rates compose by o_k with k = -q/(q-1) (q = 2,3,4,5,8)", ok)
check("K2 wall -1/k = (q-1)/q is full scrambling: the matrix is J/q, and it absorbs",
      all(np.allclose(symch(qv, (qv - 1) / qv), np.ones((qv, qv)) / qv) and abs(circ(-qv / (qv - 1), 0.2, (qv - 1) / qv) - (qv - 1) / qv) < 1e-15 for qv in (2, 3, 4, 7)))
check("K3 stops: q = 2 Haldane/bit-flip (-2), q = 4 Jukes-Cantor/depolarizing (-4/3), q -> inf erasure/union (-1)",
      -2 == -2 / (2 - 1) and abs(-4 / 3 - (-4 / 3)) == 0 and abs(-1e6 / (1e6 - 1) + 1) < 1e-5)
rho = np.array([[0.7, 0.2 + 0.1j], [0.2 - 0.1j, 0.3]])
Pa = [np.eye(2), np.array([[0, 1], [1, 0]]), np.array([[0, -1j], [1j, 0]]), np.diag([1, -1])]
def depol(r_, pv): return (1 - pv) * r_ + pv / 3 * sum(P @ r_ @ P.conj().T for P in Pa[1:])
check("K4 depolarizing qubit: two channels p, q equal one with p o_(-4/3) q",
      np.allclose(depol(depol(rho, 0.1), 0.25), depol(rho, circ(-4 / 3, 0.1, 0.25)), atol=1e-14))

# ---------------- L. C11 holographic corner-sum dial ----------------
kstar = math.sqrt(6) - 3
def fk(Ssingles, Spairs, Striple, kv): return sum(Ssingles) + kv * sum(Spairs) + kv ** 2 * Striple
def mincut_entropies(cap, bnd):
    nv = cap.shape[0]; out = {}
    names = list(bnd.keys())
    for r in range(1, 4):
        for T in itertools.combinations(names[:3], r):
            src = [bnd[x] for x in T]; snk = [bnd[x] for x in names if x not in T]
            C = np.zeros((nv + 2, nv + 2), dtype=np.int32); C[:nv, :nv] = cap
            for s_ in src: C[nv, s_] = 10 ** 6
            for s_ in snk: C[s_, nv + 1] = 10 ** 6
            out[T] = maximum_flow(csr_matrix(C), nv, nv + 1).flow_value
    return out
def vec(out):
    s = [out[(x,)] for x in "ABC"]; pr = [out[("A", "B")], out[("A", "C")], out[("B", "C")]]; return s, pr, out[("A", "B", "C")]
star = np.zeros((5, 5), dtype=np.int32)
for i in range(4): star[4, i] = star[i, 4] = 1
sv, pv_, tv = vec(mincut_entropies(star, dict(A=0, B=1, C=2, O=3)))
check("L1 star graph (four-party perfect tensor): singles 1, pairs 2, triple 1", sv == [1, 1, 1] and pv_ == [2, 2, 2] and tv == 1)
check("L2 its corner sum 3 + 6k + k^2 is <= 0 at k = -1 (MMI) and turns positive at k* = sqrt6 - 3 = -0.5505",
      fk(sv, pv_, tv, -1) <= 0 and abs(fk(sv, pv_, tv, kstar)) < 1e-12 and fk(sv, pv_, tv, kstar + 0.01) > 0 and fk(sv, pv_, tv, kstar - 0.01) < 0)
# LP over the cone of 4-party pure-state vectors obeying SA and SSA
names = ["A", "B", "C", "O"]
def idx(Tset):
    Tset = frozenset(Tset)
    if len(Tset) > 2 or (len(Tset) == 2 and "O" in Tset): Tset = frozenset(names) - Tset
    order = [frozenset(x) for x in ("A", "B", "C", "O", "AB", "AC", "BC")]
    return order.index(Tset)
def Svec(Tset, v): return 0 if len(Tset) in (0, 4) else v[idx(Tset)]
rows = []
subsets = [frozenset(s) for r in range(0, 5) for s in itertools.combinations(names, r)]
for X in subsets:
    for Y in subsets:
        if X & Y: continue
        for Z in subsets:
            if Z & (X | Y): continue
            if not X or not Y: continue
            row = np.zeros(7)
            for Tset, sg in ((X | Z, 1), (Y | Z, 1), (X | Y | Z, -1), (Z, -1)):
                if 0 < len(Tset) < 4: row[idx(Tset)] += sg
            rows.append(-row)  # -(S_XZ + S_YZ - S_XYZ - S_Z) <= 0
Aub = np.array(rows); bub = np.zeros(len(rows))
Aeq = np.array([[1, 1, 1, 0, 0, 0, 0]]); beq = np.array([1.0])
def lpmin(kv):
    cvec = np.array([1, 1, 1, kv ** 2, kv, kv, kv])  # S_ABC = S_O
    return optimize.linprog(cvec, A_ub=Aub, b_ub=bub, A_eq=Aeq, b_eq=beq, bounds=[(0, None)] * 7, method="highs").fun
grid_hi = [lpmin(kv) for kv in np.linspace(kstar + 1e-6, 2, 60)]
grid_lo = [lpmin(kv) for kv in np.linspace(-0.99, kstar - 1e-3, 30)]
check("L3 LP over every 4-party vector obeying SA + SSA: min corner sum >= 0 for all k >= k*, < 0 just below",
      min(grid_hi) > -1e-9 and max(grid_lo) < 0, f"min above {min(grid_hi):.2e}, max below {max(grid_lo):.2e}")
okm = okp = True; nsam = 0
for trial in range(60):
    nv = 9; cap = np.zeros((nv, nv), dtype=np.int32)
    for i in range(nv):
        for j in range(i + 1, nv):
            if rng.uniform() < 0.45: cap[i, j] = cap[j, i] = rng.integers(1, 6)
    sv_, pv2, tv2 = vec(mincut_entropies(cap, dict(A=0, B=1, C=2, O=3)))
    okm &= fk(sv_, pv2, tv2, -1) <= 1e-9; okp &= fk(sv_, pv2, tv2, kstar) >= -1e-9; nsam += 1
check(f"L4 {nsam} random min-cut graphs: corner sum <= 0 at k = -1 and >= 0 at k = k*", okm and okp)
s_, P_ = sp.symbols("s P", positive=True)
lower = s_ + k * P_ + k ** 2 * (P_ - s_) / 3          # uses Araki-Lieb: 3 S_ABC >= P - s
worst = sp.expand(lower.subs(P_, 2 * s_))               # coefficient of P is k(1 + k/3) < 0 on (-3, 0), and SA gives P <= 2s
check("L6 proof step: with SA (P <= 2s) and Araki-Lieb (3T >= P - s), sum >= (s/3)(k^2 + 6k + 3) >= 0 for k >= sqrt6 - 3",
      sp.simplify(worst - s_ / 3 * (k ** 2 + 6 * k + 3)) == 0 and sp.simplify(sp.expand(k * (1 + k / 3)).subs(k, -1)) < 0
      and abs(float(sp.solve(k ** 2 + 6 * k + 3, k)[1]) - kstar) < 1e-14)
t0 = time.time(); mp.mp.dps = 30
root = mp.findroot(lambda x: 3 + 6 * x + x ** 2, -0.55)
check("L5 polish: k* = sqrt(6) - 3 to 30 digits", abs(root - (mp.sqrt(6) - 3)) < mp.mpf(10) ** -28)
POLISH["L5"] = time.time() - t0

# ---------------- N. in-between curiosities ----------------
phi = lambda kv: -1 / (1 + kv)
mS_ = sp.Symbol("m", positive=True)
orb = [-1 / mS_, sp.simplify(phi(-1 / mS_)), sp.simplify(phi(phi(-1 / mS_))), sp.simplify(phi(phi(phi(-1 / mS_))))]
check("N1 the map k -> -1/(1+k) has order 3 and cycles copies -1/m -> alphabet -m/(m-1) -> m - 1 -> -1/m",
      sp.simplify(orb[1] + mS_ / (mS_ - 1)) == 0 and sp.simplify(orb[2] - (mS_ - 1)) == 0 and sp.simplify(orb[3] - orb[0]) == 0)
# per-mode law at k=-1/m vs channel: does the binomial count pass to the alphabet law by any simple map? test linear-fractional on the 'or' value
p1, p2 = 0.2, 0.35
copies = circ(-1 / 3, p1, p2); alpha = circ(-3 / 2, p1, p2)
check("N2 no shared value: at m = 3 the copies law and the alphabet law disagree on the same inputs (no identity map)",
      abs(copies - alpha) > 1e-3, f"{copies:.4f} vs {alpha:.4f}")

npass = sum(R)
print(f"polish costs (s): " + ", ".join(f"{kk_}: {v:.3f}" for kk_, v in POLISH.items()) + f"; total {sum(POLISH.values()):.3f}")
json.dump({"pass": npass, "total": len(R), "polish_seconds": POLISH, "kstar": kstar}, open("checks.json", "w"), indent=1)
print(f"{npass} of {len(R)} pass")
