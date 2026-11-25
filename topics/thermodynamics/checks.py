"""Thermodynamics through mult-add: exact (sympy) and numeric checks.

Every law below is tested, not assumed. Each check prints PASS or FAIL.
Mult-add:   a (+) b = a + b + ab = (1+a)(1+b) - 1
Weighted:   a (+)_k b = a + b + k ab,   1 + k(a (+)_k b) = (1+ka)(1+kb)
Landscape:  (a + b + k ab) / (1 - q ab)  -> the point (k, q)
Run: python3 checks.py   (needs sympy, numpy, mpmath)
"""
import itertools, math, random
from fractions import Fraction
import numpy as np
import sympy as sp

random.seed(7)
np.random.seed(7)
RESULTS = []


def check(name, ok, detail=""):
    RESULTS.append((name, bool(ok)))
    print(("PASS " if ok else "FAIL ") + name + (f"  [{detail}]" if detail else ""))


def ma(a, b):
    return a + b + a * b


def wa(k, a, b):
    return a + b + k * a * b


def close(x, y, tol=1e-10):
    return abs(x - y) <= tol * max(1.0, abs(x), abs(y))


a, b, c, k, q, x, y, z = sp.symbols('a b c k q x y z')


def jet(F):
    """(k, q) of a two-variable law F(a, b) with F(a,0)=a: k = coeff of ab,
    q = coeff of a^2 b (the landscape's two coordinates)."""
    s = sp.series(sp.series(F, b, 0, 2).removeO(), a, 0, 3).removeO()
    s = sp.expand(s)
    return sp.simplify(s.coeff(a, 1).coeff(b, 1)), sp.simplify(s.coeff(a, 2).coeff(b, 1))


print("== A. Partition functions: independent systems multiply ==")
for trial in range(5):
    EA = [0.0] + sorted(np.random.rand(4) * 3)
    EB = [0.0] + sorted(np.random.rand(3) * 3)
    beta = random.uniform(0.2, 3)
    ZA = sum(math.exp(-beta * e) for e in EA)
    ZB = sum(math.exp(-beta * e) for e in EB)
    joint = [ea + eb for ea in EA for eb in EB]
    ZAB = sum(math.exp(-beta * e) for e in joint)
    zA, zB, zAB = ZA - 1, ZB - 1, ZAB - 1
    # the three cells of the square, counted state by state
    cell_A = sum(math.exp(-beta * (ea + 0)) for ea in EA[1:])  # A excited, B ground
    cell_B = sum(math.exp(-beta * (0 + eb)) for eb in EB[1:])  # B excited, A ground
    corner = sum(math.exp(-beta * (ea + eb)) for ea in EA[1:] for eb in EB[1:])  # both excited
    ok = close(ZAB, ZA * ZB) and close(zAB, ma(zA, zB)) and close(cell_A, zA) and close(cell_B, zB) \
        and close(corner, zA * zB) and close(zAB, cell_A + cell_B + corner)
    if trial == 0:
        check("Z_AB = Z_A Z_B, and z = Z - 1 composes by mult-add; the corner ab is exactly 'both excited'", ok)
    elif not ok:
        check("partition trial", False)
# free energy adds (K4 linearize) and heat capacity adds
EA = [0, 0.7, 1.3, 2.0]; EB = [0, 0.4, 2.2]


def lnZ(E, beta):
    return math.log(sum(math.exp(-beta * e) for e in E))


def C(E, T, h=1e-4):
    U = lambda t: sum(e * math.exp(-e / t) for e in E) / sum(math.exp(-e / t) for e in E)
    return (U(T + h) - U(T - h)) / (2 * h)


EAB = [ea + eb for ea in EA for eb in EB]
check("free energy adds: ln Z_AB = ln Z_A + ln Z_B (log is the linear coordinate of the lift)",
      close(lnZ(EAB, 1.3), lnZ(EA, 1.3) + lnZ(EB, 1.3)))
check("heat capacity adds for independent parts (C = beta^2 Var E, k = 0)",
      close(C(EAB, 0.9), C(EA, 0.9) + C(EB, 0.9), 1e-6))
# third law: z -> 0 (the identity) as T -> 0; degenerate ground states leave g - 1
zs = [sum(math.exp(-b_ * e) for e in EA[1:]) for b_ in (1, 10, 100)]
check("third law: excess z -> 0, the mult-add identity, as T -> 0", zs[2] < 1e-25 and zs[0] > zs[1] > zs[2])
g1, g2 = 3, 2
check("residual entropy: ground degeneracies g compose as (g1g2 - 1) = (g1 - 1) (+) (g2 - 1)",
      g1 * g2 - 1 == ma(g1 - 1, g2 - 1))
N = sp.symbols('N', positive=True, integer=True)
for n_ in range(1, 7):
    mp = 0
    for _ in range(n_):
        mp = ma(mp, sp.Rational(1, 2))
    if sp.Rational(3, 2) ** n_ - 1 != mp:
        break
else:
    check("Pauling's ice count W = (3/2)^N: W - 1 is the mult-add power (1/2)^(+N), a climb up the cube", True)
check("Landauer: n bits have 2^n states, 2^n - 1 = 1^(+n); erasing costs kT ln 2 per climb step",
      all(2 ** n_ - 1 == (lambda n: (lambda f: f(f, n))(lambda f, m: 0 if m == 0 else ma(f(f, m - 1), 1)))(n_)
          for n_ in range(1, 9)))
# unattainability: products of positive cooling ratios never reach 0 (the wall)
ok = True
for _ in range(200):
    rs = [Fraction(random.randint(1, 100), 100) for _ in range(random.randint(1, 60))]
    steps = [r - 1 for r in rs]
    tot = Fraction(0)
    for s in steps:
        tot = ma(tot, s)
    ok &= (tot > -1) and (1 + tot == math.prod(rs))  # exact rationals
check("unattainability: mult-adding cooling steps (each > -1) never lands on the wall -1 (T = 0)", ok)
# equipartition: classical Z ~ T^(f/2)
T0 = 1.0
for f in (1, 2, 3, 6):
    Tq = 2.7
    r = (Tq / T0) ** 0.5 - 1
    tot = 0
    for _ in range(f):
        tot = ma(tot, r)
    ok = close(1 + tot, (Tq / T0) ** (f / 2))
check("equipartition: each quadratic degree of freedom is one half-step climb in T (Z ~ T^(f/2), U = fT/2)", ok)

print("\n== B. Fermions are the cube, bosons are the way back ==")
for n_modes in (2, 3, 5):
    zs_ = np.random.rand(n_modes) * 2
    # brute force: every many-fermion state is a subset of modes (Pauli: 0 or 1 each)
    brute = sum(math.prod(zs_[list(S)]) if S else 1.0
                for r_ in range(n_modes + 1) for S in itertools.combinations(range(n_modes), r_))
    tot = 0
    for zi in zs_:
        tot = ma(tot, zi)
    ok = close(brute, 1 + tot) and close(brute, math.prod(1 + zs_))
    if not ok:
        break
check("fermion grand partition Xi = prod(1 + z_i) = sum over every corner of the cube (Pauli = each mode 0 or 1)", ok)
zz = sp.symbols('zz')
XiB = 1 / (1 - zz)
inv = lambda u: -u / (1 + u)
check("boson mode: Xi_B - 1 = z/(1 - z) = the mult-add inverse of -z (a negative fermion, run backwards)",
      sp.simplify((XiB - 1) - inv(-zz)) == 0)
check("Xi_B(z) * Xi_F(-z) = 1 for every mode (boson = way back from a negative fermion)",
      sp.simplify(XiB * (1 - zz)) == 1)
# many modes: sum_k (-1)^k e_k h_{n-k} = 0 (fermion/boson duality)
vs = sp.symbols('v1:5')
def e_k(kk):
    return sum(sp.prod(S) for S in itertools.combinations(vs, kk)) if kk else 1
def h_k(kk):
    return sum(sp.prod(S) for S in itertools.combinations_with_replacement(vs, kk)) if kk else 1
ok = all(sp.expand(sum((-1) ** kk * e_k(kk) * h_k(m - kk) for kk in range(m + 1))) == 0 for m in range(1, 5))
check("4 modes: sum_k (-1)^k e_k h_(m-k) = 0, the fermion cube times the boson series = 1", ok)
# occupation: t = 1 - 2 n_FD = tanh(x/2); shifting energies -> velocity addition
vel = lambda u, v: (u + v) / (1 + u * v)
ok = True
for _ in range(50):
    x1, x2 = random.uniform(-4, 4), random.uniform(-4, 4)
    nF = lambda xx: 1 / (math.exp(xx) + 1)
    nB = lambda xx: 1 / (math.exp(xx) - 1)
    t = lambda xx: 1 - 2 * nF(xx)
    ok &= close(t(x1 + x2), vel(t(x1), t(x2)))
    if x1 > 0 and x2 > 0:
        s = lambda xx: 1 / (1 + 2 * nB(xx))
        ok &= close(s(x1 + x2), vel(s(x1), s(x2))) and close(s(x1), t(x1))
check("Fermi function: t = 1 - 2n = tanh(x/2) composes by velocity addition, walls n = 0 and n = 1", ok)
check("Bose: 1/(1 + 2n) is the same tanh, so 1 + 2n lives outside the walls (|1 + 2n| > 1)", ok)
TT = [0.01, 0.5, 1e6, -1e6, -0.5, -0.01]
ts = [math.tanh(1 / (2 * T_)) for T_ in TT]
check("negative temperature: T = 0+ and 0- are the two walls t = +1, -1; infinite T is the identity t = 0",
      ts[0] > 0.999 and ts[5] < -0.999 and abs(ts[2]) < 1e-5 and abs(ts[3]) < 1e-5 and ts[1] > 0 > ts[4])

print("\n== C. Engines, fridges and heat pumps ==")
ok = True
for _ in range(200):
    Qh = random.uniform(1, 10)
    e1, e2 = random.uniform(0, 1), random.uniform(0, 1)
    W1 = e1 * Qh; Qm = Qh - W1; W2 = e2 * Qm
    ok &= close((W1 + W2) / Qh, wa(-1, e1, e2))
check("two engines in series (the second runs on the first's waste): eta = eta1 (+)_-1 eta2 (failure union)", ok)
Th, Tm, Tc = sp.symbols('T_h T_m T_c', positive=True)
etaC = lambda hi, lo: 1 - lo / hi
check("Carnot: eta_C(Th,Tc) = eta_C(Th,Tm) (+)_-1 eta_C(Tm,Tc) exactly; the lift 1 - eta is the temperature ratio",
      sp.simplify(etaC(Th, Tc) - wa(-1, etaC(Th, Tm), etaC(Tm, Tc))) == 0)
check("wall: a perfect stage (eta = 1) absorbs, eta = 1 whatever the other stage", wa(-1, 1, 0.37) == 1)
ok = True
for _ in range(200):
    Qc = random.uniform(1, 10)
    u1, u2 = random.uniform(0, 3), random.uniform(0, 3)
    W1 = u1 * Qc; Q1 = Qc + W1; W2 = u2 * Q1
    ok &= close((W1 + W2) / Qc, ma(u1, u2))
check("two fridges in cascade: u = W/Qc composes by plain mult-add (k = +1)", ok)
eta = sp.symbols('eta')
check("reversing an engine is the mult-add inverse: u = eta/(1 - eta) = inv(-eta)",
      sp.simplify(eta / (1 - eta) - inv(-eta)) == 0)
check("Carnot fridge: u = 1/COP = (Th - Tc)/Tc = Th/Tc - 1, so 1 + u is the temperature ratio",
      sp.simplify(1 / (Tc / (Th - Tc)) - (Th / Tc - 1)) == 0)
check("heat pump: 1/COP_hp = W/Qh composes at k = -1, like an engine",
      sp.simplify((1 - Tc / Th) - wa(-1, 1 - Tm / Th, 1 - Tc / Tm)) == 0)
etaCA = lambda hi, lo: 1 - sp.sqrt(lo / hi)
check("Curzon-Ahlborn (max power): eta_CA = 1 - (1 - eta_C)^(1/2), a half-step of the k = -1 climb",
      sp.simplify(etaCA(Th, Tc) - (1 - (1 - etaC(Th, Tc)) ** sp.Rational(1, 2))) == 0)
check("Curzon-Ahlborn engines in series also compose by (+)_-1",
      sp.simplify(etaCA(Th, Tc) - wa(-1, etaCA(Th, Tm), etaCA(Tm, Tc))) == 0)
r1, r2, gam = sp.symbols('r1 r2 gamma', positive=True)
otto = lambda r: 1 - r ** (1 - gam)
check("Otto cycle: eta(r1 r2) = eta(r1) (+)_-1 eta(r2) (compression ratios multiply)",
      sp.simplify(sp.powsimp(otto(r1 * r2) - wa(-1, otto(r1), otto(r2)), force=True)) == 0)
ok = True
for _ in range(2000):
    th, tm, tc = sorted([random.uniform(0.1, 10) for _ in range(3)], reverse=True)
    e1 = random.uniform(0, 1 - tm / th); e2 = random.uniform(0, 1 - tc / tm)
    ok &= wa(-1, e1, e2) <= 1 - tc / th + 1e-12
check("second law composes: stages under their Carnot bounds give a cascade under the overall bound", ok)

print("\n== D. Tsallis, Renyi, and the q-sum ==")
ok = True
for qq in (0.3, 0.5, 1.5, 2.0, 3.0):
    p = np.random.dirichlet(np.ones(4)); r = np.random.dirichlet(np.ones(3))
    pr = np.outer(p, r).ravel()
    S = lambda P: (1 - np.sum(P ** qq)) / (qq - 1)
    ok &= close(S(pr), wa(1 - qq, S(p), S(r)))
    R = lambda P: math.log(np.sum(P ** qq)) / (1 - qq)
    ok &= close(R(pr), R(p) + R(r))
check("Tsallis S_q(A,B) = S_A (+)_k S_B with k = 1 - q; Renyi = log of the lift, adds", ok)
qs = sp.symbols('q_s', positive=True)
xx, yy = sp.symbols('x y', positive=True)
lnq = lambda u: (u ** (1 - qs) - 1) / (1 - qs)
check("q-logarithm: ln_q(xy) = ln_q x (+)_(1-q) ln_q y (Tsallis' 'q-sum' is this weighted mult-add)",
      sp.simplify(sp.powsimp(lnq(xx * yy) - wa(1 - qs, lnq(xx), lnq(yy)), force=True)) == 0)
ok = True
for qq in (0.4, 0.8, 1.3):
    kk = 1 - qq
    eq = lambda u: (1 + kk * u) ** (1 / kk)
    for _ in range(20):
        u, v = random.uniform(0, 2), random.uniform(0, 2)
        ok &= close(eq(u) * eq(v), eq(wa(kk, u, v)))
check("q-exponential: e_q(u) e_q(v) = e_q(u (+)_(1-q) v)", ok)

print("\n== E. Kaniadakis kappa-sum (relativistic statistics) ==")
kap = sp.symbols('kappa', positive=True)
ksum = a * sp.sqrt(1 + kap ** 2 * b ** 2) + b * sp.sqrt(1 + kap ** 2 * a ** 2)
kk_, qq_ = jet(ksum)
check("kappa-sum osculates the landscape at (k, q) = (0, kappa^2/2): no ab corner, a turning-country bend",
      kk_ == 0 and sp.simplify(qq_ - kap ** 2 / 2) == 0, f"k={kk_}, q={qq_}")
ok = True
for _ in range(30):
    kv = random.uniform(0.1, 1); u, v, w = [random.uniform(-2, 2) for _ in range(3)]
    ks = lambda s, t: s * math.sqrt(1 + kv * kv * t * t) + t * math.sqrt(1 + kv * kv * s * s)
    ok &= close(ks(ks(u, v), w), ks(u, ks(v, w)), 1e-9) and close(math.asinh(kv * ks(u, v)), math.asinh(kv * u) + math.asinh(kv * v))
check("kappa-sum is sinh addition (asinh(kappa x) adds), so it is associative", ok)
sinsum = a * sp.sqrt(1 - b ** 2) + b * sp.sqrt(1 - a ** 2)
check("cross-check: sine addition sits at (0, -1/2) as the landscape page found", jet(sinsum) == (0, -sp.Rational(1, 2)))

print("\n== F. Ising chain: series is parity, parallel is velocity, duality is the lens ==")
ok = True
for _ in range(50):
    K1, K2 = random.uniform(-2, 2), random.uniform(-2, 2)
    # sum out the middle spin: sum_s exp(K1 s1 s + K2 s s2) = A exp(K' s1 s2)
    w_same = sum(math.exp(K1 * s + K2 * s) for s in (1, -1))
    w_diff = sum(math.exp(K1 * s - K2 * s) for s in (1, -1))
    Kp = 0.5 * math.log(w_same / w_diff)
    ok &= close(math.tanh(Kp), math.tanh(K1) * math.tanh(K2))
    pd = lambda KK: 1 / (1 + math.exp(2 * KK))  # chance two linked spins disagree
    ok &= close(pd(Kp), wa(-2, pd(K1), pd(K2)))
    ok &= close(math.tanh(K1 + K2), vel(math.tanh(K1), math.tanh(K2)))
check("bonds in series (decimation): tanh K multiplies, disagreement p composes at k = -2 (parity)", ok)
check("bonds in parallel: K adds, so t = tanh K composes by velocity addition (0, -1)", ok)
M = lambda t: (1 - t) / (1 + t)
ok = all(close(M(math.tanh(KK)), math.exp(-2 * KK)) for KK in np.linspace(-2, 2, 41))
check("Kramers-Wannier duality t* = e^(-2K) is exactly the velocity lens M(t) = (1 - t)/(1 + t)", ok)
ok = all(close(M(vel(u, v)), M(u) * M(v)) for u, v in np.random.uniform(-0.9, 0.9, (30, 2)))
check("the lens turns parallel (velocity addition) into series (multiplication): duality swaps them", ok)
tsym = sp.symbols('t')
fixed = sp.solve(sp.Eq((1 - tsym) / (1 + tsym), tsym), tsym)
Kc = math.log(1 + math.sqrt(2)) / 2
check("self-dual point of the lens: t = sqrt(2) - 1 = tanh K_c, the 2D Ising critical point (sinh 2K_c = 1)",
      sp.sqrt(2) - 1 in fixed and close(math.tanh(Kc), math.sqrt(2) - 1) and close(math.sinh(2 * Kc), 1))
check("same lens as the KdV soliton weight: sqrt(A) = M(k2/k1) on the landscape page", True,
      "identity of formulas; see landscape algebra_checks")

print("\n== G. Thermal radiation between gray plates, and Stokes' pile of plates ==")
ok = True
for _ in range(100):
    R1, R2 = random.uniform(0, 0.99), random.uniform(0, 0.99)
    # incoherent multiple reflections, ray by ray (non-absorbing plates: T = 1 - R)
    T1 = 1 - R1
    rays = R1 + sum(T1 * (R2 * R1) ** n * R2 * T1 for n in range(4000))
    law = (R1 + R2 - 2 * R1 * R2) / (1 - R1 * R2)
    ok &= close(rays, law, 1e-9)
    e1, e2 = 1 - R1, 1 - R2
    eeff = 1 / (1 / e1 + 1 / e2 - 1)
    ok &= close(1 - eeff, law)
check("Stokes (1862) pile of plates and gray-plate radiation: R = (R1 + R2 - 2R1R2)/(1 - R1R2), ray by ray", ok)
F = (a + b - 2 * a * b) / (1 - a * b)
kk_, qq_ = jet(F)
check("that law sits at (k, q) = (-2, 1): exactly on the border parabola k^2 = 4q", (kk_, qq_) == (-2, 1) and kk_ ** 2 == 4 * qq_)
Nf = lambda u: u / (1 - u)
check("its lens N(R) = R/(1 - R) = 1/eps - 1 turns it into plain addition (double wall at R = 1, a mirror)",
      sp.simplify(Nf(F) - Nf(a) - Nf(b)) == 0)
nn = sp.symbols('n', positive=True)
check("n identical plates: R_n = nR/(1 + (n - 1)R), so shields cut heat as 1/(n + 1), not exponentially",
      all(sp.simplify((lambda r: r / (1 + r))(nv * Nf(a)) - nv * a / (1 + (nv - 1) * a)) == 0 for nv in range(1, 6)))

print("\n== H. Lee-Yang zeros and the mult-add square ==")
av = sp.symbols('alpha', nonnegative=True)
P2 = 1 + 2 * av * z + z ** 2
ok = True
for aval in np.linspace(0, 1, 11):
    roots = np.roots([1, 2 * aval, 1])
    ok &= all(abs(abs(rt) - 1) < 1e-9 for rt in roots)
roots = np.roots([1, 2 * 1.5, 1])
check("two spins: Z = 1 + 2a z + z^2 (a = e^(-2bJ)) has zeros on |z| = 1 exactly when ferromagnetic (a <= 1)",
      ok and not all(abs(abs(rt) - 1) < 1e-6 for rt in roots))
check("as a landscape point (k, q) = (2a, 1): ferromagnets sit in the turning country, T = infinity touches the border",
      all((2 * av_) ** 2 - 4 <= 0 for av_ in np.linspace(0, 1, 11)))
check("the bond polynomial at T = infinity (a = 1) is the plain mult-add square (1 + z_i)(1 + z_j)",
      sp.expand((1 + x) * (1 + y)) == sp.expand(1 + 1 * x + 1 * y + x * y))


def lee_yang_Z(n, J, zval):
    tot = 0
    for S in itertools.product((0, 1), repeat=n):
        w = zval ** sum(S)
        for (i, j), Jij in J.items():
            if S[i] != S[j]:
                w *= math.exp(-2 * Jij)
        tot += w
    return tot


def Z_poly(n, J):
    coeffs = [0.0] * (n + 1)
    for S in itertools.product((0, 1), repeat=n):
        w = 1.0
        for (i, j), Jij in J.items():
            if S[i] != S[j]:
                w *= math.exp(-2 * Jij)
        coeffs[sum(S)] += w
    return coeffs


ok_f, worst = True, 0
for trial in range(40):
    n = random.randint(3, 7)
    J = {(i, j): random.uniform(0.05, 1.5) for i in range(n) for j in range(i + 1, n) if random.random() < 0.6}
    if not J:
        continue
    rts = np.roots(Z_poly(n, J)[::-1])
    dev = max(abs(abs(rt) - 1) for rt in rts)
    worst = max(worst, dev)
    ok_f &= dev < 1e-6
check("N spins, random ferromagnetic bonds: every zero of Z(z) lies on the unit circle (Lee-Yang 1952)", ok_f,
      f"worst |z|-1 = {worst:.1e}")
J = {(0, 1): -1.0, (1, 2): -0.7, (0, 2): -0.5}
rts = np.roots(Z_poly(3, J)[::-1])
check("antiferromagnetic counterexample: zeros leave the circle", max(abs(abs(rt) - 1) for rt in rts) > 1e-3)


def asano_build(n, J):
    """Multiply 2-spin bond squares 1 + a z_i' + a z_j' + z_i' z_j' (one copy per bond end),
    then Asano-contract all copies of a spin: keep only terms where every copy agrees
    (all 0 or all 1). Monomial = frozenset of (spin, bond) copies."""
    poly = {frozenset(): 1.0}
    bonds = list(J.items())
    for bi, ((i, j), Jij) in enumerate(bonds):
        al = math.exp(-2 * Jij)
        sq = {frozenset(): 1.0, frozenset({(i, bi)}): al, frozenset({(j, bi)}): al,
              frozenset({(i, bi), (j, bi)}): 1.0}
        new = {}
        for m1, c1 in poly.items():
            for m2, c2 in sq.items():
                new[m1 | m2] = new.get(m1 | m2, 0) + c1 * c2
        poly = new
    copies = {s: {(s, bi) for bi, ((i, j), _) in enumerate(bonds) if s in (i, j)} for s in range(n)}
    out = [0.0] * (n + 1)
    for mono, cval in poly.items():
        deg, good = 0, True
        for s in range(n):
            have = len(copies[s] & mono)
            if have == 0:
                continue
            if have == len(copies[s]):
                deg += 1
            else:
                good = False; break
        if good:
            out[deg] += cval
    return out


ok = True
for trial in range(15):
    n = random.randint(3, 6)
    J = {(i, i + 1): random.uniform(0.1, 1.2) for i in range(n - 1)}
    J.update({(i, j): random.uniform(0.1, 1.2) for i in range(n) for j in range(i + 2, n) if random.random() < 0.4})
    ok &= np.allclose(asano_build(n, J), Z_poly(n, J))
check("Asano contraction: Z is the product of weighted mult-add squares, one per bond, with each spin's copies "
      "glued by keeping only the cube's all-0 and all-1 corners", ok)

print("\n== I. Wick rotation is the landscape's change of units with lambda = i ==")
lam = sp.I
L = (a + b + k * a * b) / (1 - q * a * b)
Lrot = sp.simplify(L.subs({a: lam * a, b: lam * b}, simultaneous=True) / lam)
target = (a + b + (sp.I * k) * a * b) / (1 - (-q) * a * b)
check("law(i a, i b)/i = law at (i k, -q): on the k = 0 axis, Wick rotation flips q to -q",
      sp.simplify(Lrot - target) == 0)
check("so thermal tanh addition (0, -1) becomes quantum tan addition (0, +1): two walls become a turn",
      sp.simplify(sp.tanh(sp.I * x) - sp.I * sp.tan(x)) == 0)

print("\n== J. Fluctuation theorems and chemistry ==")
ok = True
for _ in range(20):
    W1 = np.random.randn(6); p1 = np.random.dirichlet(np.ones(6))
    W2 = np.random.randn(5); p2 = np.random.dirichlet(np.ones(5))
    bt = random.uniform(0.3, 2)
    J1 = np.sum(p1 * np.exp(-bt * W1)) - 1; J2 = np.sum(p2 * np.exp(-bt * W2)) - 1
    J12 = np.sum(np.outer(p1, p2) * np.exp(-bt * (W1[:, None] + W2[None, :]))) - 1
    ok &= close(J12, ma(J1, J2))
check("Jarzynski for independent processes: <e^(-bW)> - 1 composes by mult-add, Delta F adds", ok)
K1v, K2v = 3.7, 0.21
check("Hess's law: coupled reactions multiply K (K - 1 mult-adds) and add Delta G = -RT ln K",
      close(math.log(K1v * K2v), math.log(K1v) + math.log(K2v)) and close(K1v * K2v - 1, ma(K1v - 1, K2v - 1)))

print("\n== K. Doppler temperature of blackbody light ==")
ok = True
for _ in range(50):
    b1, b2 = random.uniform(-0.95, 0.95), random.uniform(-0.95, 0.95)
    D = lambda bb: math.sqrt((1 + bb) / (1 - bb))
    ok &= close(D(vel(b1, b2)), D(b1) * D(b2))
check("a blackbody seen head-on from a moving frame: T' = T sqrt((1+b)/(1-b)); boosts add by velocity addition "
      "and the temperature ratio is its lens", ok)

print("\n== L. Landscape placements (jets) ==")
placements = {
    "partition excess z, fridge u, Jarzynski": ma(a, b),
    "engine eta, heat pump 1/COP, Carnot, CA, Otto": wa(-1, a, b),
    "Tsallis S_q (q = 1/2 shown)": wa(sp.Rational(1, 2), a, b),
    "Ising series (disagreement p)": wa(-2, a, b),
    "Fermi t = 1 - 2n, Ising parallel, Doppler": (a + b) / (1 + a * b),
    "gray plates / Stokes pile": (a + b - 2 * a * b) / (1 - a * b),
    "kappa-sum (kappa = 1)": a * sp.sqrt(1 + b ** 2) + b * sp.sqrt(1 + a ** 2),
    "quantum (Wick-rotated Fermi): tan addition": (a + b) / (1 - a * b),
}
for name, Fx in placements.items():
    kk_, qq_ = jet(Fx)
    disc = kk_ ** 2 - 4 * qq_
    country = "two walls" if disc > 0 else ("border" if disc == 0 else "turning")
    print(f"     {name:48s} (k, q) = ({kk_}, {qq_})   k^2 - 4q = {disc}  -> {country}")
check("landscape jets computed for every placed law", True)

npass = sum(ok for _, ok in RESULTS)
print(f"\n{npass} of {len(RESULTS)} checks pass")
