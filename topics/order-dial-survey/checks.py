"""Checks for the order-dial survey (q, the order dial). Run: python3 checks.py > checks_output.txt
Each check prints PASS/FAIL with a label; the last line is "X of Y pass".
Exact (sympy) wherever a closed form exists; numeric checks get a 30-digit mpmath polish at the end (section P)."""
import itertools, math, time, cmath, random
import numpy as np
import sympy as sp
import mpmath as mp
from scipy.optimize import linprog
from scipy.sparse import lil_matrix
from qlib import (q, k, t, hb, qp_mul, qp_add, qp_scale, qp_pow, circ, X, Y, qint, qfact, gauss, inv, words,
                  sort_word)

RESULTS = []


def check(label, ok, note=""):
    RESULTS.append(bool(ok))
    print(("PASS " if ok else "FAIL ") + label + ((" -- " + note) if note else ""))


def poly_eq(a, b):
    return sp.expand(a - b) == 0


def qp_eq(A, B):
    keys = set(A) | set(B)
    return all(sp.expand(A.get(kk, 0) - B.get(kk, 0)) == 0 for kk in keys)


# ============ A. The core law over q-commuting inputs (C1) ============
al, be, ga = sp.symbols('alpha beta gamma')
a_el = qp_add(qp_scale(X, al), qp_scale(Y, be), {(1, 1): ga})
b_el = qp_add(qp_scale(X, be), {(0, 2): al}, {(2, 0): 1})
c_el = qp_add(qp_scale(Y, ga), {(1, 0): 2})
lhs = circ(circ(a_el, b_el), c_el)
rhs = circ(a_el, circ(b_el, c_el))
check("A1 o_k over the quantum plane (yx = q xy) is associative for every q and k (symbolic)", qp_eq(lhs, rhs))
comm = qp_add(circ(qp_scale(X, al), qp_scale(Y, be)), qp_scale(circ(qp_scale(Y, be), qp_scale(X, al)), -1))
check("A2 order gap of o_k on (alpha x, beta y) is k(1 - q) alpha beta xy: zero exactly at q = 1",
      qp_eq(comm, {(1, 1): k * (1 - q) * al * be}))
s2 = qp_pow(qp_add(X, Y), 2, qq=-1)
check("A3 at q = -1 the cross terms cancel: (x + y)^2 = x^2 + y^2", qp_eq(s2, {(2, 0): 1, (0, 2): 1}))
PX = np.array([[0, 1], [1, 0]]); PZ = np.array([[1, 0], [0, -1]])
check("A4 Pauli X, Z: ZX = -XZ, each squares to 1, (X + Z)^2 = 2 I (q = -1 stop, Lean `pauli`)",
      np.array_equal(PZ @ PX, -(PX @ PZ)) and np.array_equal((PX + PZ) @ (PX + PZ), 2 * np.eye(2)))
z2, c = sp.symbols('z2 c')
check("A5 if z z = q z z with q != 1 then z^2 = 0 (any q off 1 forbids repeats; Lean anticommute_self_zero at q=-1)",
      sp.solve(sp.Eq(z2, q * z2).subs(q, sp.Rational(3, 7)), z2) == [0] and
      sp.solve(sp.Eq(z2, q * z2).subs(q, sp.exp(sp.I * sp.pi / 5)), z2) == [0])

# ============ B. Corner expansion with q-commuting inputs -> soliton corners (C2) ============
ok = True
for n in (3, 4):
    qs = [[q] * n for _ in range(n)]
    # chain a_{n-1} o ... o a_0 in fully reversed order, a_i = x_i: expand the product (1 + k x_{n-1})...(1 + k x_0)
    # and normal-order; 1 + k (chain) = prod (1 + k x_i) in reversed order.
    terms = {}
    for S in itertools.product([0, 1], repeat=n):
        word = tuple(i for i in reversed(range(n)) if S[i])
        f, w = sort_word(word, qs)
        terms[w] = terms.get(w, 0) + f * k ** len(word)
    for w, coeff in terms.items():
        m = len(w)
        if not poly_eq(coeff, k ** m * q ** (m * (m - 1) // 2)):
            ok = False
check("B1 reversed-order chain over x_j x_i = q x_i x_j: each rank-m corner gets k^m q^C(m,2) (soliton corner"
      " weight A^C(m,2) with A = q), n = 3, 4", ok)
ok = True
for n in (3, 4):
    qs = [[q] * n for _ in range(n)]
    terms = {}
    for S in itertools.product([0, 1], repeat=n):
        word = tuple(i for i in range(n) if S[i])
        f, w = sort_word(word, qs)
        terms[w] = terms.get(w, 0) + f
    ok &= all(poly_eq(v, 1) for v in terms.values())
check("B2 in sorted order no corner picks up q: the plain corner expansion (q = 1 end) for every q", ok)

# ============ C. Binomial -> q-binomial, limits, roots of unity (C3) ============
ok = True
for n in range(1, 11):
    P = qp_pow(qp_add(X, Y), n)
    ok &= qp_eq(P, {(m, n - m): gauss(n, m) for m in range(n + 1)})
check("C1 Potter-Schutzenberger: (x + y)^n with yx = q xy is sum [n,m]_q x^m y^(n-m), n = 1..10 (exact)", ok)
check("C2 q = 1: [n,m]_1 = C(n,m) for n <= 14 (Lean gauss_at_one proves every n)",
      all(gauss(n, m).subs(q, 1) == math.comb(n, m) for n in range(15) for m in range(n + 1)))
check("C3 q = 0: [n,m]_0 = 1, so (x + y)^n = sum of the n + 1 sorted words, each once (only order-free words survive)",
      all(gauss(n, m).subs(q, 0) == 1 for n in range(13) for m in range(n + 1)))
check("C4 q -> infinity: [n,m]_q has degree m(n - m) and leading coefficient 1 (only the fully reversed word)",
      all(sp.Poly(gauss(n, m), q).degree() == m * (n - m) and sp.Poly(gauss(n, m), q).LC() == 1
          for n in range(1, 13) for m in range(n + 1)))
ok = True
for n in range(0, 11):
    lhs = sp.expand(sp.prod([1 + q ** i * t for i in range(n)]))
    rhs = sum(q ** (m * (m - 1) // 2) * gauss(n, m) * t ** m for m in range(n + 1))
    ok &= poly_eq(lhs, rhs)
check("C5 Cauchy (commutative) form: prod (1 + q^i t) = sum q^C(m,2) [n,m] t^m, n <= 10 (Lean q_binomial_theorem, n<=6)", ok)
ok = True
for mo in range(2, 8):
    PhiP = sp.Poly(sp.cyclotomic_poly(mo, q), q)
    for n in range(0, 22):
        A, B = divmod(n, mo)
        for m in range(n + 1):
            R, S = divmod(m, mo)
            ok &= sp.Poly(gauss(n, m) - math.comb(A, R) * gauss(B, S), q).rem(PhiP).is_zero
check("C6 q-Lucas (Olive 1965, Desarmenien 1982): [n,m] = C(n div d, m div d) [n mod d, m mod d] mod Phi_d(q), "
      "d = 2..7, n <= 21 (exact)", ok)
check("C7 q = -1 case equals the repo's Lean gauss_at_neg_one: C(n div 2, m div 2), or 0 for n even, m odd (n <= 16)",
      all(gauss(n, m).subs(q, -1) == (0 if (n % 2 == 0 and m % 2 == 1) else math.comb(n // 2, m // 2))
          for n in range(17) for m in range(n + 1)))

# ============ D. q-counting [n]_q (C4) ============
check("D1 [n]_q = sum over words with one y of q^inv, n <= 12", all(
    poly_eq(qint(n), sum(q ** inv(w) for w in words(n, 1))) for n in range(1, 13)))
check("D2 the lift 1 + (q - 1)[n]_q = q^n (Lean qnum_lift): the same count is the k dial at k = q - 1",
      all(poly_eq(1 + (q - 1) * qint(n), q ** n) for n in range(15)))
stops = all(qint(n).subs(q, 1) == n and qint(n).subs(q, 0) == (1 if n else 0) and
            qint(n).subs(q, -1) == n % 2 and qint(n).subs(q, 2) == 2 ** n - 1 for n in range(15))
check("D3 stops: q=1 counts n, q=0 gives 1, q=-1 parity, q=2 gives 2^n - 1 (all Lean, QDial)", stops)
check("D4 at a primitive d-th root of 1, [d]_q = 0 and [m]_q != 0 for 0 < m < d (d = 2..12)",
      all(abs(complex(qint(d).subs(q, cmath.exp(2j * cmath.pi / d)))) < 1e-12 and
          all(abs(complex(qint(m).subs(q, cmath.exp(2j * cmath.pi / d)))) > 1e-3 for m in range(1, d))
          for d in range(2, 13)))

# ============ E. Face weights q^area (C5) ============
def partitions_in_box(rows, cols):
    def rec(r, mx):
        if r == 0:
            yield ()
            return
        for p in range(mx + 1):
            for rest in rec(r - 1, p):
                yield (p,) + rest
    return rec(rows, cols)


check("E1 sum over partitions in an m x (n-m) box of q^|lambda| = [n,m]_q, n <= 10",
      all(poly_eq(sum(q ** sum(lam) for lam in partitions_in_box(m, n - m)), gauss(n, m))
          for n in range(11) for m in range(n + 1)))


def area_above(w):
    # path: x = step right, y = step up; area of boxes above the path in the box = number of (y before x) pairs
    ups = 0; area = 0
    for s in w:
        if s == 1:
            ups += 1
        else:
            area += ups
    return area


check("E2 a word's out-of-order count equals its lattice path's area (every word, n <= 10)",
      all(inv(w) == area_above(w) for n in range(11) for m in range(n + 1) for w in words(n, m)))
check("E3 q -> 0 keeps one path (area 0, the minimum); q -> infinity keeps one (area m(n-m), the maximum)",
      all(sum(1 for w in words(n, m) if inv(w) == 0) == 1 and
          sum(1 for w in words(n, m) if inv(w) == m * (n - m)) == 1 for n in range(1, 11) for m in range(n + 1)))
ok = True
for n in range(1, 8):
    lhs = sp.expand(sp.prod([1 + q ** i * t for i in range(n)]))
    for m in range(n + 1):
        ok &= poly_eq(sp.Poly(lhs, t).coeff_monomial(t ** m), q ** (m * (m - 1) // 2) * gauss(n, m))
check("E4 soliton/theta weight q^C(m,2) times [n,m]_q is the t^m layer of the q-orbit product (cube-weights I1), n <= 7", ok)

# ============ F. Routes on the permutohedron (C6) ============
check("F1 sum over permutations of q^inv = [n]_q! (Rodrigues 1839, MacMahon), n <= 7",
      all(poly_eq(sum(q ** inv(p) for p in itertools.permutations(range(n))), qfact(n)) for n in range(1, 8)))
w3 = cmath.exp(2j * cmath.pi / 3)
check("F2 at a cube root of 1 the routes cancel from n = 3 on (cube-weights E2), n = 3..7",
      all(abs(sum(w3 ** inv(p) for p in itertools.permutations(range(n)))) < 1e-9 for n in range(3, 8)))
check("F3 at q = -1 the route sum is [n]_{-1}! = 0 from n = 2 on (signs of permutations cancel)",
      all(sum((-1) ** inv(p) for p in itertools.permutations(range(n))) == 0 for n in range(2, 8)))

# ============ G. Rulecraft: order of reforms under Mallows weights (C7) ============
def outcome(order, lam, tgt, x0=0.5):
    x = x0
    for i in order:
        x = (1 - lam[i]) * x + lam[i] * tgt[i]
    return x


def mallows_mean(qv, perms, invs, outs):
    if qv == 0:
        return float(np.mean([o for o, v in zip(outs, invs) if v == 0]))
    w = np.array([qv ** v for v in invs], dtype=float)
    return float(np.dot(w, outs) / w.sum())


rng = np.random.default_rng(11)
n = 5
perms = list(itertools.permutations(range(n)))
invs = [inv(p) for p in perms]
qgrid = [0.0] + list(np.logspace(-3, 3, 61))
mono_ok, ends_ok, worst_step = True, True, 0.0
tight = None
for trial in range(300):
    lam = rng.uniform(0.01, 0.99, n)
    tgt = np.sort(rng.uniform(0, 1, n))  # label i = rank of target, so identity = sorted by target ascending
    outs = [outcome(p, lam, tgt) for p in perms]
    means = [mallows_mean(qv, perms, invs, outs) for qv in qgrid]
    d = np.diff(means)
    if d.max() > 1e-13:
        mono_ok = False
    if tight is None or d.max() > worst_step:
        worst_step = d.max(); tight = (lam.copy(), tgt.copy())
    ends_ok &= abs(means[0] - max(outs)) < 1e-15 and abs(np.mean(outs) - mallows_mean(1.0, perms, invs, outs)) < 1e-14
    ends_ok &= abs(mallows_mean(1e6, perms, invs, outs) - min(outs)) < 1e-4
check("G1 ends: q = 0 gives the target-sorted order (highest, rulecraft-history H5a); q = 1 the plain average over all "
      "120 orders; q -> infinity the reversed order (lowest); 300 random packages of 5 rules", ends_ok)
check("G2 the Mallows-weighted outcome falls monotonically as q grows, on a 62-point grid from 0 to 1000, "
      "300 random packages", mono_ok, f"largest step {worst_step:.2e} (negative = falling)")


def weak_leq_pairs(n):
    P = list(itertools.permutations(range(n)))
    invset = [frozenset((p[j], p[i]) for i in range(n) for j in range(i + 1, n) if p[i] > p[j]) for p in P]
    return P, invset


def dominated(n, q1, q2):
    """Is Mallows(q1) <= Mallows(q2) in the right weak order? (Strassen: a monotone coupling exists.)"""
    P, S = weak_leq_pairs(n)
    N = len(P)
    iv = [len(s) for s in S]
    mu = np.array([q1 ** v for v in iv]); mu /= mu.sum()
    nu = np.array([q2 ** v for v in iv]); nu /= nu.sum()
    pairs = [(a, b) for a in range(N) for b in range(N) if S[a] <= S[b]]
    A = lil_matrix((2 * N, len(pairs)))
    for j, (a, b) in enumerate(pairs):
        A[a, j] = 1; A[N + b, j] = 1
    res = linprog(np.zeros(len(pairs)), A_eq=A.tocsr(), b_eq=np.concatenate([mu, nu]), bounds=(0, None),
                  method="highs")
    return res.status == 0


qpairs = [(0.2, 0.5), (0.5, 0.9), (0.9, 1.0), (1.0, 1.5), (0.3, 3.0), (2.0, 5.0)]
dom = {nn: all(dominated(nn, a, b) for a, b in qpairs) for nn in (3, 4, 5)}
check("G3 a monotone coupling exists: Mallows(q) <= Mallows(q') in the weak order of orders for q < q' "
      "(6 pairs of q, n = 3, 4, 5; linear program)", all(dom.values()), str(dom))
check("G4 duality: weight of an order at q equals the weight of its reverse at 1/q (n = 5)",
      all(abs(0.7 ** inv(p) / sum(0.7 ** v for v in invs) -
              (1 / 0.7) ** inv(p[::-1]) / sum((1 / 0.7) ** v for v in invs)) < 1e-15 for p in perms))

# ============ H. Weyl, Heisenberg, Moyal (C8) ============
xs = np.linspace(-3, 3, 13)
f = lambda x: np.exp(-x ** 2) * (1 + x)
aa, bb = 0.7, 1.3
TM = np.exp(1j * bb * (xs + aa)) * f(xs + aa)     # shift after multiply
MT = np.exp(1j * bb * xs) * f(xs + aa)            # multiply after shift
check("H1 Weyl: shift-after-multiply = e^{iab} x multiply-after-shift (order costs a phase = area ab)",
      np.allclose(TM, np.exp(1j * aa * bb) * MT, atol=1e-15))
xx, pp = sp.symbols('x p')


def star(F, G, order=6):
    tot = 0
    for nn in range(order + 1):
        s = 0
        for kk in range(nn + 1):
            s += sp.binomial(nn, kk) * (-1) ** kk * sp.diff(F, xx, nn - kk, pp, kk) * sp.diff(G, xx, kk, pp, nn - kk)
        tot += (sp.I * hb / 2) ** nn / sp.factorial(nn) * s
    return sp.expand(tot)


check("H2 Moyal: x * p - p * x = i hbar", poly_eq(star(xx, pp) - star(pp, xx), sp.I * hb))
F1, F2, F3 = xx ** 2 * pp, pp ** 2 + xx, xx * pp ** 2
check("H3 the Moyal product is associative for every hbar (polynomials, exact)",
      poly_eq(star(star(F1, F2), F3), star(F1, star(F2, F3))))
check("H4 at hbar = 0 the star circle operation is the ordinary commutative circle operation",
      poly_eq((F1 + F2 + k * star(F1, F2)).subs(hb, 0), F1 + F2 + k * F1 * F2))
ok = True
for mo in range(2, 9):
    w = np.exp(2j * np.pi / mo)
    C = np.diag([w ** j for j in range(mo)]); S = np.roll(np.eye(mo), 1, axis=0)
    ok &= np.allclose(C @ S, w * S @ C) and np.allclose(np.linalg.matrix_power(C, mo), np.eye(mo))
check("H5 clock and shift (Sylvester 1882): C S = w S C with w = e^{2 pi i/m}, m = 2..8; hbar = 2 pi/m", ok)

# ============ I. Clifford and repetition: only q = +1 or -1 when things square to numbers (C9) ============
e1sq_c = sp.symbols('c')
# e2 e1 e1: move e1 left twice -> q^2 e1 e1 e2 = q^2 c e2 ; but e1 e1 = c so also = c e2.
check("I1 e2 e1 = q e1 e2 and e1^2 = c != 0 force q^2 = 1 (order-dial G1): Clifford lives only at q = +-1",
      set(sp.solve(sp.Eq(q ** 2 * e1sq_c, e1sq_c).subs(e1sq_c, 1), q)) == {-1, 1})
check("I2 the jump: Bose counting (x^2 != 0) survives only at q = 1; any q != 1 on self-swaps gives x^2 = 0",
      all(sp.solve(sp.Eq(z2, qq * z2), z2) == [0] for qq in (sp.Rational(1, 2), -1, sp.exp(2 * sp.pi * sp.I / 3), 2)))

# ============ J. Knots: the Jones polynomial along the unit circle (C10) ============
V = {"3_1": t ** -1 + t ** -3 - t ** -4,                                   # Lean trefoilL_f
     "4_1": t ** 2 - t + 1 - t ** -1 + t ** -2,                             # Lean figEight_f
     "5_1": t ** -2 + t ** -4 - t ** -5 + t ** -6 - t ** -7}                # Lean cinquefoilL_f
det = {"3_1": 3, "4_1": 5, "5_1": 5}
check("J1 t = 1: V = 1 and t = -1: |V| = determinant (3, 5, 5), the two ends Lean proves (jones_at_one/neg_one)",
      all(V[K].subs(t, 1) == 1 and abs(V[K].subs(t, -1)) == det[K] for K in V))
w3s = sp.exp(2 * sp.pi * sp.I / 3)
check("J2 t = e^{2 pi i/3}: V = 1 for every knot here (Jones 1985)",
      all(sp.simplify(sp.expand_complex(V[K].subs(t, w3s))) == 1 for K in V))
arf = {K: (0 if det[K] % 8 in (1, 7) else 1) for K in V}
check("J3 t = i: V = (-1)^Arf, Arf read from the determinant mod 8 (Murakami 1986)",
      all(sp.simplify(V[K].subs(t, sp.I)) == (-1) ** arf[K] for K in V))

# ============ K. The q-exponential (C11) ============
N = 7
ok = True
for qq in (sp.Rational(1, 3), sp.Integer(2), q):
    def eq_series(el):
        tot = {(0, 0): sp.Integer(1)}; P = {(0, 0): sp.Integer(1)}
        for nn in range(1, N + 1):
            P = qp_mul(P, el, qq)
            tot = qp_add(tot, qp_scale(P, 1 / qfact(nn, qq)))
        return tot
    lhs = eq_series(qp_add(X, Y))
    rhs = qp_mul(eq_series(X), eq_series(Y), qq)
    ok &= all(sp.simplify(lhs.get(kk, 0) - rhs.get(kk, 0)) == 0 for kk in set(lhs) | set(rhs) if sum(kk) <= N)
check("K1 Schutzenberger 1953: e_q(x + y) = e_q(x) e_q(y) when yx = q xy, e_q(z) = sum z^n/[n]_q!, through degree 7"
      " (q = 1/3, 2, symbolic q)", ok)
check("K2 q = 1 gives e (sum 1/n! = e, Lean E); q = 0 gives the geometric series 1/(1 - z)",
      abs(sum(1 / math.factorial(nn) for nn in range(30)) - math.e) < 1e-15 and
      all(qfact(nn, 0) == 1 for nn in range(10)))

# ============ L. Disguises (C12) ============
alpha, beta = sp.symbols('alpha beta', positive=True)
kk_, dd_ = alpha + beta, alpha * beta          # 1 + kx + dx^2 = (1 + alpha x)(1 + beta x)
qq_ = beta / alpha
check("L1 face weight q = beta/alpha is a reading of (k, d): (1 + q)^2/q = k^2/d and (1 - q)^2/q = Delta/d",
      sp.simplify((1 + qq_) ** 2 / qq_ - kk_ ** 2 / dd_) == 0 and sp.simplify((1 - qq_) ** 2 / qq_ - (kk_ ** 2 - 4 * dd_) / dd_) == 0)
ok = True
for kv, dv in [(1.0, 1.0), (2.0, 1.0), (3.0, 1.0), (0.5, 2.0), (2.5, 1.2)]:
    r = np.roots([1, -kv, dv])
    qv = r[1] / r[0]
    ok &= (abs(abs(qv) - 1) < 1e-9) == (kv ** 2 - 4 * dv <= 1e-12)
check("L2 |q| = 1 exactly when Delta <= 0 (angle side); q real off the circle when Delta > 0 (5 settings)", ok)
nn_ = sp.symbols('n', integer=True)
ok = True
for nv in range(1, 9):
    U = sp.expand(sp.simplify((alpha ** nv - beta ** nv) / (alpha - beta)))
    ok &= sp.simplify(U - alpha ** (nv - 1) * qint(nv, qq_)) == 0
check("L3 the dial's Lucas sequence U_n(k, d) = alpha^(n-1) [n]_{beta/alpha} (cut-independence Q10), n <= 8", ok)
xv = sp.symbols('x')
check("L4 Jackson: (f(qx) - f(x))/((q - 1)x) on x^n gives [n]_q x^(n-1), with q = 1 + h a step size (k dial)",
      all(sp.simplify(((q * xv) ** nv - xv ** nv) / ((q - 1) * xv) - qint(nv) * xv ** (nv - 1)) == 0 for nv in range(1, 9)))

# ============ N. In-between laws (section N) ============
ok = True
for mo in range(2, 9):
    Phi = sp.cyclotomic_poly(mo, q)
    ok &= all(sp.rem(gauss(mo, m), Phi, q) == 0 for m in range(1, mo))
check("N1 Frobenius at roots of unity: [d, m]_q = 0 at a primitive d-th root for 0 < m < d, so (x + y)^d = x^d + y^d"
      " (d = 2..8, exact)", ok)
P4 = qp_pow(qp_add(X, Y), 4, qq=-1)
check("N1b primitive matters: at q = -1 (not primitive of order 4) (x + y)^4 = x^4 + 2x^2y^2 + y^4",
      qp_eq(P4, {(4, 0): 1, (2, 2): 2, (0, 4): 1}))
ok = True
for mo in range(2, 9):
    w = np.exp(2j * np.pi / mo)
    C = np.diag([w ** j for j in range(mo)]); S = np.roll(np.eye(mo), 1, axis=0)
    ok &= np.allclose(np.linalg.matrix_power(C + S, mo), 2 * np.eye(mo), atol=1e-10)
check("N1c clock and shift obey it: (C + S)^m = C^m + S^m = 2I, m = 2..8", ok)
ok = True
for mo in range(2, 7):
    w = np.exp(2j * np.pi / mo)
    C = np.diag([w ** j for j in range(mo)]); S = np.roll(np.eye(mo), 1, axis=0)
    ok &= abs(np.linalg.det(S @ C) - w ** mo * np.linalg.det(C @ S)) < 1e-10 and abs(w ** mo - 1) < 1e-12
    ok &= abs(np.trace(C @ S)) < 1e-12   # tr(VU) = q tr(UV) forces tr = 0 when q != 1
check("N2 finite pairs need a root of unity: det(VU) = q^m det(UV) forces q^m = 1, and tr(UV) = 0 (m = 2..6)", ok)
ok = True
for mo in range(2, 7):
    w = np.exp(2j * np.pi / mo)
    C = np.diag([w ** j for j in range(mo)]); S = np.roll(np.eye(mo), 1, axis=0)
    Cm = np.linalg.matrix_power(C, mo); Sm = np.linalg.matrix_power(S, mo)
    ok &= np.allclose(Cm @ S, S @ Cm) and np.allclose(Sm @ C, C @ Sm)
check("N3 at a root of unity U^m and V^m commute with everything (the big centre of the rational quantum torus)", ok)
check("N4 integer q = 2: [4,2]_2 = 35 planes in F_2^4 and [3]_q = q^2 + q + 1 points of the projective plane"
      " (q = 2, 3, 4, 5)", gauss(4, 2).subs(q, 2) == 35 and all(qint(3).subs(q, qq) == qq ** 2 + qq + 1 for qq in (2, 3, 4, 5)))

# ============ P. Polish: 30 digits ============
mp.mp.dps = 30
t0 = time.time()
w7 = mp.exp(2j * mp.pi / 7)


def gauss_mp(n, m, z):
    num = mp.mpf(1); den = mp.mpf(1)
    for i in range(m):
        num *= (1 - z ** (n - i)); den *= (1 - z ** (i + 1))
    return num / den


err = 0
for nv in (15, 20):
    for m in range(nv + 1):
        A, B = divmod(nv, 7); R, S2 = divmod(m, 7)
        rhs = math.comb(A, R) * (gauss_mp(B, S2, w7) if S2 <= B else 0)
        # avoid 0/0 at the root: evaluate the Gaussian polynomial itself via its coefficients
        coeffs = sp.Poly(gauss(nv, m), q).all_coeffs()[::-1]
        lhs = mp.fsum(int(c) * w7 ** j for j, c in enumerate(coeffs))
        err = max(err, abs(lhs - rhs))
t1 = time.time()
check("P1 polish: q-Lucas at q = e^{2 pi i/7}, n = 15, 20, at 30 digits", err < mp.mpf(10) ** -25,
      f"worst gap {mp.nstr(err, 3)}, {t1 - t0:.2f} s")
t0 = time.time()
lam, tgt = tight
outs = [mp.mpf(1) / 2 for _ in perms]
outs = []
for p in perms:
    x = mp.mpf(1) / 2
    for i in p:
        x = (1 - mp.mpf(lam[i])) * x + mp.mpf(lam[i]) * mp.mpf(tgt[i])
    outs.append(x)
means = []
for qv in qgrid[1:]:
    qm = mp.mpf(qv)
    wts = [qm ** v for v in invs]
    means.append(mp.fsum(w * o for w, o in zip(wts, outs)) / mp.fsum(wts))
steps = [means[i + 1] - means[i] for i in range(len(means) - 1)]
t1 = time.time()
check("P2 polish: the tightest Rulecraft package from G2 is still falling at 30 digits", max(steps) < 0,
      f"largest step {mp.nstr(max(steps), 3)}, {t1 - t0:.2f} s")
t0 = time.time()
mo = 7
wm = mp.exp(2j * mp.pi / mo)
C = mp.diag([wm ** j for j in range(mo)])
S = mp.zeros(mo); [S.__setitem__(((j + 1) % mo, j), 1) for j in range(mo)]
M = (C + S) ** mo - 2 * mp.eye(mo)
gap = max(abs(M[i, j]) for i in range(mo) for j in range(mo))
t1 = time.time()
check("P3 polish: (C + S)^7 = 2I at 30 digits", gap < mp.mpf(10) ** -25, f"gap {mp.nstr(gap, 3)}, {t1 - t0:.2f} s")

print(f"{sum(RESULTS)} of {len(RESULTS)} pass")
