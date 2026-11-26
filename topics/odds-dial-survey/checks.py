"""Checks for the odds-dial survey. Each check prints PASS/FAIL with a label; last line 'X of Y pass'.
Run from this folder: python3 checks.py > checks_output.txt   (numpy, sympy, mpmath; about 20 s)."""
import itertools, math, time, json
import numpy as np, sympy as sp, mpmath as mp

R = []
def check(label, ok, info=""):
    R.append(bool(ok)); print(("PASS " if ok else "FAIL ") + label + (f"  [{info}]" if info else ""))

def subsets(n):
    return itertools.product([0, 1], repeat=n)

def circ(a, b): return a + b + a * b
rng = np.random.default_rng(7)
out = {}

# ---------------- C1. Power set: empty, plain, full, and the exponent ----------------
n = 6
def total(ts):
    return sum(np.prod([t if s else 1.0 for t, s in zip(ts, S)]) for S in subsets(len(ts)))
check("C1a all t = 0 counts 1 (empty set only)", abs(total([0.0]*n) - 1) < 1e-12)
check("C1a all t = 1 counts 2^n (plain power set)", abs(total([1.0]*n) - 2**n) < 1e-9)
errs = [abs(total([2**a - 1]*n) - 2**(a*n)) / 2**(a*n) for a in np.linspace(0, 2, 9)]
check("C1b uniform t = 2^a - 1 counts 2^(a n), a in [0, 2]", max(errs) < 1e-12, f"max rel err {max(errs):.1e}")
a, b = sp.symbols('a b', positive=True)
check("C1c exponent ruler (2^a-1) o_1 (2^b-1) = 2^(a+b)-1", sp.simplify(circ(2**a-1, 2**b-1) - (2**(a+b)-1)) == 0)
for T in [10.0, 100.0, 1e4]:
    pfull = (T/(1+T))**n
pf = [(T/(1+T))**n for T in [1, 10, 100, 1e4]]
check("C1d forcing: P(S = everything) -> 1 as t -> infinity", pf[-1] > 0.999 and all(np.diff(pf) > 0), str([round(x, 4) for x in pf]))
ts = rng.uniform(0.1, 3, n)
ok = True
for S in subsets(n):
    w = np.prod([t if s else 1 for t, s in zip(ts, S)])
    wc = np.prod([1/t if (1 - s) else 1 for t, s in zip(ts, S)])
    ok &= abs(w - np.prod(ts) * wc) < 1e-9
check("C1e complement: weight(S; t) = (prod t) * weight(S^c; 1/t)", ok)
# one element alone: total is linear in t1 between 'element left out' and 'element forced'
rest = ts[1:]
lin = [total([t1] + list(rest)) for t1 in [0, 0.5, 1, 2]]
check("C1f one element alone: total = (1 + t1) * total(rest), straight line in t1",
      np.allclose(lin, [(1 + t1) * total(list(rest)) for t1 in [0, 0.5, 1, 2]]))

# ---------------- C2. Bayes: multiplying the odds ----------------
p, lH, lN = sp.symbols('p l_H l_N', positive=True)
odds = lambda q: q / (1 - q)
post = p*lH / (p*lH + (1-p)*lN)
check("C2a Bayes: posterior odds = prior odds x likelihood ratio", sp.simplify(odds(post) - odds(p)*lH/lN) == 0)
# log-odds add over many independent looks
L = rng.uniform(0.3, 3, 40); prior = 0.2
q = prior
for l in L: q = q*l / (q*l + (1-q))
check("C2b 40 looks: log posterior odds = log prior odds + sum log L",
      abs(math.log(q/(1-q)) - (math.log(prior/(1-prior)) + np.log(L).sum())) < 1e-9)
O1, O2 = sp.symbols('O1 O2', positive=True)
v = lambda O: (O - 1)/(O + 1)
check("C2c multiplying odds is velocity addition: v(O1 O2) = (v1 + v2)/(1 + v1 v2)",
      sp.simplify(v(O1*O2) - (v(O1)+v(O2))/(1+v(O1)*v(O2))) == 0)
# haploid selection: allele freq recurrence vs odds x (1+s)^g
s_, x0, g = 0.07, 0.01, 60
x = x0
for _ in range(g): x = x*(1+s_)/(x*(1+s_) + (1-x))
check("C2d selection: allele odds after g generations = odds0 * (1+s)^g",
      abs(x/(1-x) - x0/(1-x0)*(1+s_)**g) < 1e-12*(1+s_)**g)

# ---------------- C3. Logistic reading, cube vs simplex, binding, entropy ----------------
tt = np.logspace(-3, 3, 13)
chance = tt/(1+tt)
check("C3a chance t/(1+t) = logistic of log t; 1/2 exactly at t = 1",
      np.allclose(chance, 1/(1+np.exp(-np.log(tt)))) and abs(1/(1+1) - 0.5) == 0)
T = np.array([0.5, 1.0, 2.0, 0.3]); T2 = T.copy(); T2[0] = 5.0
cube = lambda T: T/(1+T); simp = lambda T: T/(1+T.sum())
c_ok = np.allclose(cube(T)[1:], cube(T2)[1:])
s_drop = np.all(simp(T2)[1:] < simp(T)[1:]) and np.allclose(simp(T2)[1:]/simp(T2)[1], simp(T)[1:]/simp(T)[1])
check("C3b turning one odds alone: cube chances of the others stay; simplex chances drop, ratios kept (Luce)", c_ok and s_drop)
small = np.array([1e-4, 2e-4, 3e-4])
check("C3c small odds: cube and simplex readings agree to first order",
      np.max(np.abs(cube(small) - simp(small))/small) < 1e-3)
w = sp.symbols('omega', positive=True); xx = sp.symbols('x', real=True)
disc = sp.discriminant(w*xx**2 + 2*xx + 1, xx)
check("C3d two binding sites 1 + 2x + w x^2: no real zero exactly when w > 1 (discriminant 4(1 - w))",
      sp.simplify(disc - 4*(1 - w)) == 0)
h = lambda q: -(q*np.log2(q) + (1-q)*np.log2(1-q))
grid = np.logspace(-3, 3, 2001); H = n*h(grid/(1+grid))
check("C3e entropy of the dialled power set peaks at t = 1 with n bits",
      abs(grid[np.argmax(H)] - 1) < 0.01 and abs(H.max() - n) < 1e-6)

# ---------------- C4. Small odds: the classical limit ----------------
def occ(k, t, cut=4000):
    # occupation of one element with k-cells, repetition-dial sign: weight (1 + k t)^(1/k)
    if k == 0: return t
    if k > 0:
        m = round(1/k); wts = [math.comb(m, j)*(t/m)**j for j in range(m+1)]
    else:
        wts = [(t)**j for j in range(cut)]   # k = -1: any number of copies (bosons)
    Z = sum(wts); return sum(j*w_ for j, w_ in enumerate(wts))/Z
okc = all(abs(occ(1, t) - t/(1+t)) < 1e-12 and abs(occ(-1, t) - t/(1-t)) < 1e-9 and abs(occ(0.5, t) - t/(1+0.5*t)) < 1e-12
          for t in [0.01, 0.1, 0.5, 0.9])
check("C4a occupation n = t/(1 + k t) at k = 1 (Fermi), 1/2 (two copies), 0, -1 (Bose)", okc)
rel = [abs(t/(1+t) - t)/t for t in [1e-1, 1e-2, 1e-3]]
check("C4b Fermi minus Boltzmann is relative size t: the gap closes linearly", np.allclose(rel, [0.1/1.1, 0.01/1.01, 0.001/1.001]))
from scipy.stats import binom, poisson
tvs = []
for N in [10, 100, 1000]:
    kk = np.arange(0, 40); lam = 2.0
    tvs.append(0.5*np.abs(binom.pmf(kk, N, (lam/N)/(1+lam/N)) - poisson.pmf(kk, lam)).sum())
check("C4c N elements at odds 2/N: count law -> Poisson(2), distance falls like 1/N", tvs[0] > tvs[1] > tvs[2] and tvs[2] < 3e-3,
      ", ".join(f"{x:.1e}" for x in tvs))
check("C4d (1 + 1/n)^n climbs to e (E.stair)", all(np.diff([(1+1/m)**m for m in [1, 10, 100, 1000]]) > 0) and abs((1+1e-6)**1e6 - math.e) < 2e-6)

# ---------------- C5. Lee-Yang: uniform odds is the fugacity ----------------
def gas_poly(n_sites, edges, k):
    c = np.zeros(n_sites+1)
    for S in subsets(n_sites):
        e = sum(S[i] and S[j] for i, j in edges)
        c[sum(S)] += k**e
    return c
ring = [(i, (i+1) % 12) for i in range(12)]
k = 2.0
roots = np.roots(gas_poly(12, ring, k)[::-1])
check("C5a ring of 12, pair weight k = 2: every zero on |k t| = 1", np.max(np.abs(np.abs(k*roots) - 1)) < 1e-6,
      f"max dev {np.max(np.abs(np.abs(k*roots)-1)):.1e}")
L4 = 4; tor = []
for i in range(L4):
    for j in range(L4):
        s0 = i*L4+j; tor += [(s0, i*L4+(j+1) % L4), (s0, ((i+1) % L4)*L4+j)]
ctor = gas_poly(16, tor, k)
zt = np.roots(ctor[::-1]); zeta = zt * k**2   # normalised odds: t k^(deg/2), deg = 4
check("C5b 4x4 torus, k = 2: every zero on |t k^2| = 1 (Lee-Yang circle)", np.max(np.abs(np.abs(zeta) - 1)) < 1e-6,
      f"max dev {np.max(np.abs(np.abs(zeta)-1)):.1e}")
cz = ctor * (k**2) ** (-np.arange(17))   # coefficients in normalised odds zeta
check("C5c complement t -> 1/t: in normalised odds the polynomial is palindromic", np.allclose(cz, cz[::-1]))
# multivariate: all |zeta_i| < 1 -> Z != 0
def Zmulti(zs, edges, k, deg):
    tot = 0
    for S in subsets(len(zs)):
        e = sum(S[i] and S[j] for i, j in edges)
        tot += np.prod([zs[i]/k**(deg/2) if S[i] else 1 for i in range(len(zs))]) * k**e
    return tot
ring6 = [(i, (i+1) % 6) for i in range(6)]
mins = min(abs(Zmulti(rng.uniform(0, 0.999, 6)*np.exp(2j*np.pi*rng.uniform(size=6)), ring6, 2.0, 2)) for _ in range(1500))
check("C5d one element at a time: Z != 0 whenever every normalised odds has |zeta_i| < 1 (1500 samples)", mins > 1e-6, f"min |Z| {mins:.2e}")
# real positive odds: density 1/2 at zeta = 1
zz = 1.0
dens = sum(i*c*zz**i for i, c in enumerate(cz))/sum(c*zz**i for i, c in enumerate(cz))/16
dens5 = dens
check("C5e on the real axis the density is exactly 1/2 at normalised odds 1", abs(dens - 0.5) < 1e-12)
hc = np.roots(gas_poly(12, ring, 0.0)[::-1])
check("C5f hard core (k = 0): zeros on the negative axis, none in (-1/4, 0)", np.max(np.abs(hc.imag)) < 1e-6 and np.max(hc.real) < -0.25)
angs = []
for kk_ in [1.5, 2.5, 4.0]:
    zr = np.roots(gas_poly(16, tor, kk_)[::-1]) * kk_**2
    angs.append(np.min(np.abs(np.angle(zr))))
check("C5g stronger pairs pull the nearest zero toward normalised odds 1 (angle shrinks)", angs[0] > angs[1] > angs[2], str([round(x, 3) for x in angs]))
out["C5"] = dict(ring=[[z.real, z.imag] for z in k*roots], torus=[[z.real, z.imag] for z in zeta], angles=angs)

# ---------------- C6. Primes: t_p = p^(-s) ----------------
from sympy import primerange
P = list(primerange(2, 200000))
for s in [2.0, 3.0]:
    prod = np.prod([1 + p_**-s for p_ in P])
    target = float(mp.zeta(s)/mp.zeta(2*s))
    check(f"C6a product over primes < 2e5 of (1 + p^-{int(s)}) = zeta({int(s)})/zeta({int(2*s)})", abs(prod - target) < 1e-5, f"{prod:.8f} vs {target:.8f}")
X = 10**6
mu2 = np.ones(X+1, dtype=bool); mu2[0] = False
for p_ in primerange(2, 1001): mu2[p_*p_::p_*p_] = False
dens = mu2.sum()/X
check("C6b squarefree numbers up to 1e6: density 6/pi^2", abs(dens - 6/math.pi**2) < 1e-4, f"{dens:.6f} vs {6/math.pi**2:.6f}")
xs = 0.37
tower = np.prod([1 + xs**(2**j) for j in range(8)])
check("C6c a boson is a tower of fermions: prod (1 + x^(2^j)) = 1/(1 - x) (PowersetDials.fermi_dirac)", abs(tower - 1/(1-xs)) < 1e-12)
odd_sf = sum(nn**-2.0 for nn in range(1, 200001) if mu2[nn] and nn % 2 == 1)
check("C6d one prime alone: t_2 -> 0 leaves the odd squarefree numbers, zeta(2)/zeta(4)/(1 + 1/4)",
      abs(odd_sf - (15/math.pi**2)/(1.25)) < 1e-5, f"{odd_sf:.7f}")

# ---------------- C7. Holographic gluing: odds of the two RT phases ----------------
def X_of(A, B):
    (a1, a2), (b1, b2) = A, B
    return (a2-a1)*(b2-b1)/((b1-a2)*(b2-a1))
def I_rt(A, B):  # in units of c/3
    return max(0.0, math.log(X_of(A, B)))
ok_g, ok_m = True, True
for _ in range(500):
    a1, a2, m0, d1, d2 = np.cumsum(rng.uniform(0.05, 2, 5))
    A, B1, B2, B12 = (a1, a2), (m0, d1), (d1, d2), (m0, d2)
    X1, X2, X12 = X_of(A, B1), X_of(A, B2), X_of(A, B12)
    ok_g &= abs(X12 - circ(X1, X2)) < 1e-9*max(1, X12)
    I3 = I_rt(A, B1) + I_rt(A, B2) - I_rt(A, B12)
    ok_m &= I3 <= 1e-12
check("C7a 500 interval triples: odds glue X(A, B1B2) = X1 o_1 X2", ok_g)
check("C7b monogamy I3 <= 0 in every configuration (HoloGluing.monogamy_glued)", ok_m)
lp = lambda X: math.log(X) if X > 1 else 0.0
I3f = lambda X1, X2: lp(X1) + lp(X2) - lp(circ(X1, X2))
check("C7c mutual information turns on at odds X = 1 exactly", lp(1.0) == 0 and lp(0.999) == 0 and lp(1.0001) > 0)

# ---------------- C8. k copies: tilting multiplies the odds ----------------
def pgf(k, mu, z):
    return math.exp(-mu*(1-z)) if k == 0 else (1 + k*mu*(1-z))**(-1/k)
ok_t = True
for kk_ in [-1.0, 0.0, 1.0, 0.5]:
    for mu in [0.2, 0.7]:
        for t in [0.3, 2.0]:
            mu_t = t*mu/(1 + kk_*mu*(1-t))
            for z in [0.1, 0.6]:
                ok_t &= abs(pgf(kk_, mu, t*z)/pgf(kk_, mu, t) - pgf(kk_, mu_t, z)) < 1e-12
            ok_t &= abs(mu_t/(1+kk_*mu_t) - t*mu/(1+kk_*mu)) < 1e-12
check("C8a tilting one mode by t^N stays in its family; the odds mu/(1 + k mu) is multiplied by t (k = -1, 0, 1/2, 1)", ok_t)
Lm = rng.normal(size=(5, 5)); Lm = Lm @ Lm.T / 5
Km = Lm @ np.linalg.inv(np.eye(5) + Lm)
t = 2.5
Kt = t*Km @ np.linalg.inv(np.eye(5) - Km + t*Km)
wts = {S: t**sum(S)*np.linalg.det(Lm[np.ix_([i for i in range(5) if S[i]], [i for i in range(5) if S[i]])]) if sum(S) else 1.0 for S in subsets(5)}
Zt = sum(wts.values()); incl = [sum(w_ for S, w_ in wts.items() if S[i])/Zt for i in range(5)]
check("C8b L-ensemble: tilting by t^|S| is the DPP with kernel tK(I - K + tK)^-1 (inclusion chances)", np.allclose(incl, np.diag(Kt)))
check("C8c t = 1: one fermion level half full; one boson level's mean diverges as t -> 1",
      abs(occ(1, 1.0) - 0.5) < 1e-12 and occ(-1, 0.999, cut=200000) > 900)

# ---------------- C9. Geometric odds: the q-binomial theorem ----------------
qq, T_ = sp.symbols('q t')
nq = 6
lhs = sp.expand(sp.prod([1 + qq**i*T_ for i in range(nq)]))
def gauss(n_, k_):
    num = sp.prod([1 - qq**(n_-j) for j in range(k_)]); den = sp.prod([1 - qq**(j+1) for j in range(k_)])
    return sp.cancel(num/den)
rhs = sp.expand(sum(qq**(kk*(kk-1)//2)*gauss(nq, kk)*T_**kk for kk in range(nq+1)))
check("C9a odds t_i = q^i t: prod (1 + q^i t) = sum q^(k(k-1)/2) [n,k]_q t^k (n = 6)", sp.expand(lhs - rhs) == 0)
check("C9b q = 1 brings back the binomial count (1 + t)^n", sp.expand(lhs.subs(qq, 1) - (1+T_)**nq) == 0)
Nn = 30
dist = sp.Poly(sp.expand(sp.prod([1 + qq**i for i in range(1, Nn+1)])), qq).all_coeffs()[::-1][:Nn+1]
def distinct_parts(m):
    c = [0]*(m+1); c[0] = 1
    for part in range(1, m+1):
        for s0 in range(m, part-1, -1): c[s0] += c[s0-part]
    return c
pent = sp.Poly(sp.expand(sp.prod([1 - qq**i for i in range(1, Nn+1)])), qq).all_coeffs()[::-1][:Nn+1]
pent_ok = all(int(pent[m]) == sum((1 if j % 2 == 0 else -1) for j in range(-10, 11) if j*(3*j-1)//2 == m) for m in range(Nn+1))
check("C9c odds q^i: distinct-part partitions; odds -q^i: Euler's pentagonal numbers (to q^30)",
      [int(c) for c in dist] == distinct_parts(Nn) and pent_ok)

# ---------------- C10. One bond alone: deletion <-> contraction ----------------
# Wheatstone bridge: nodes s,a,b,t; edges sa, sb, at, bt, ab(bridge)
E = [(0, 1), (0, 2), (1, 3), (2, 3), (1, 2)]
def connected(open_edges):
    par = list(range(4))
    def f(x):
        while par[x] != x: x = par[x]
        return x
    for (u, v_), o in zip(E, open_edges):
        if o: par[f(u)] = f(v_)
    return f(0) == f(3)
def rel(ps):
    return sum(np.prod([pp if o else 1-pp for pp, o in zip(ps, S)]) for S in subsets(5) if connected(S))
pp = 0.3
vals = [rel([pp]*4 + [pe]) for pe in [0, 0.25, 0.5, 1]]
lin_ok = np.allclose(vals, [vals[0] + pe*(vals[3]-vals[0]) for pe in [0, 0.25, 0.5, 1]])
par_ = lambda x, y: x + y - x*y
check("C10a bridge bond alone: reliability is a straight line in its chance t/(1+t)", lin_ok)
check("C10b ends: t = 0 is two series pairs in parallel; t = infinity two parallel pairs in series (GraphTheory.wheatstone)",
      abs(vals[0] - par_(pp*pp, pp*pp)) < 1e-12 and abs(vals[3] - par_(pp, pp)**2) < 1e-12)
check("C10c whole curve matches 2p^2 + 2p^3 - 5p^4 + 2p^5 at bond chance p",
      abs(rel([pp]*5) - (2*pp**2 + 2*pp**3 - 5*pp**4 + 2*pp**5)) < 1e-12)

# ---------------- N1. t = 1, the self-dual stop ----------------
fix = sp.solve(sp.Eq(1/sp.Symbol('t', positive=True), sp.Symbol('t', positive=True)))
check("N1a the complement t -> 1/t has one positive fixed point, t = 1 (PowersetDials.compl_fixed_dial)", fix == [1])
check("N1b at t = 1: chance 1/2, Bayes neutral (L = 1), fermion half full, MI onset, entropy max, Lee-Yang crossing",
      abs(1/(1+1) - 0.5) == 0 and abs(occ(1, 1.0) - 0.5) < 1e-12 and lp(1.0) == 0 and abs(dens5 - 0.5) < 1e-12 and abs(grid[np.argmax(H)] - 1) < 0.01)

# ---------------- N2. t = sqrt2 - 1, the union half ----------------
r2 = math.sqrt(2) - 1
partner = lambda t: (1 - t)/(1 + t)
check("N2a union half: t o_1 t = 1 and partner(t) = t (PowersetDials.half_dial, partner_fixed_iff, Thermodynamics.lens_fixed_point)",
      abs(circ(r2, r2) - 1) < 1e-15 and abs(partner(r2) - r2) < 1e-15)
on = [I3f(x_, x_) < 0 for x_ in [0.40, 0.41, 0.4142, 0.4143, 0.42, 0.9]]
check("N2b two equal neighbours, neither sharing MI with A: I3 < 0 exactly when their odds pass sqrt2 - 1", on == [False, False, False, True, True, True])
okb = True
for _ in range(400):
    X1 = rng.uniform(0, 1); X2 = rng.uniform(0, 1)
    okb &= (I3f(X1, X2) < 0) == (X2 > partner(X1))
check("N2c the edge of that region is the partner curve X2 = (1 - X1)/(1 + X1) (400 random pairs)", okb)
# real intervals: A = [-1, 0], B1 = [g, g+L], B2 = [g+L, g+2L]; find g with X1 = X2? not possible for adjacent -> use X(A,B12)=1 edge
def I3_geom(g, L):
    A, B1, B2, B12 = (-1.0, 0.0), (g, g+L), (g+L, g+2*L), (g, g+2*L)
    return I_rt(A, B1) + I_rt(A, B2) - I_rt(A, B12), X_of(A, B1), X_of(A, B2), X_of(A, B12)
gs = np.linspace(0.05, 3, 300)
agree, edge = True, True
for g_ in gs:
    I3, X1, X2, X12 = I3_geom(g_, 1.0)
    agree &= abs(I3 - I3f(X1, X2)) < 1e-9
    if X1 < 1 and X2 < 1:
        edge &= (I3 < -1e-12) == (X12 > 1)
edge_ok = agree and edge
check("N2d real intervals on a line: RT lengths give the odds formula, and I3 turns negative exactly where X1 o_1 X2 passes 1", edge_ok)
check("N2e CHSH's pi/8: the union half is tan(pi/8), but the Born odds there is tan^2(pi/8) = 3 - 2 sqrt2: no shared law",
      abs(math.tan(math.pi/8) - r2) < 1e-15 and abs(math.tan(math.pi/8)**2 - (3 - 2*math.sqrt(2))) < 1e-15)

# ---------------- N3. golden ratio ----------------
phi = (1 + math.sqrt(5))/2; tg = 1/phi; pg = tg/(1+tg)
check("N3a t = 1/phi: t o_1 t = 1/t, i.e. 'both miss' (1-p)^2 equals 'one hits' p = 1/phi^2",
      abs(circ(tg, tg) - 1/tg) < 1e-14 and abs((1-pg)**2 - pg) < 1e-15 and abs(pg - 1/phi**2) < 1e-15)
fib = gas_poly(12, [(i, i+1) for i in range(11)], 0.0).sum()
check("N3b hard-core path at t = 1 counts Fibonacci F(14) = 377: growth phi per site, not a special t", fib == 377)

# ---------------- Polish (mpmath, 30 digits) ----------------
mp.mp.dps = 30
t0 = time.time()
r = mp.findroot(lambda X: X*(2+X) - 1, 0.4)
check("P1 polish: onset of I3 < 0 is sqrt2 - 1 to 30 digits", abs(r - (mp.sqrt(2)-1)) < mp.mpf(10)**-29, mp.nstr(r, 30))
zr = mp.polyroots([mp.mpf(int(c)) for c in gas_poly(12, ring, 2.0)[::-1]], maxsteps=200, extraprec=200)
check("P2 polish: ring zeros refined at 30 digits sit on |2t| = 1", max(abs(abs(2*z) - 1) for z in zr) < mp.mpf(10)**-25)
check("P3 polish: zeta(2)/zeta(4) = 15/pi^2 to 30 digits", abs(mp.zeta(2)/mp.zeta(4) - 15/mp.pi**2) < mp.mpf(10)**-29)
check("P4 polish: golden point (1-p)^2 = p at 30 digits", abs((1 - 1/mp.phi**2)**2 - 1/mp.phi**2) < mp.mpf(10)**-29)
pol = time.time() - t0
print(f"polish cost {pol:.2f} s")
out["polish_seconds"] = pol
json.dump(out, open("results.json", "w"), indent=1, default=float)
print(f"{sum(R)} of {len(R)} pass")
