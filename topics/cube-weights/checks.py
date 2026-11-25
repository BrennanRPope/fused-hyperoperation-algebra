"""Cube weights: which parts of the n-cube can carry a weight (DRAFT, Nov 25 2026).

Run: python3 checks.py  (needs sympy, numpy). Prints PASS/FAIL per check.
Sections: A corners, B edges, C squares, D higher faces, E paths,
F orientations (Koszul), G Mobius with the dial, H dial plane d, I q-orbit.
"""
import itertools, math, random
from fractions import Fraction
import numpy as np
import sympy as sp

random.seed(1); np.random.seed(1)
RESULTS = []
def check(name, ok, info=""):
    RESULTS.append(bool(ok))
    print(("PASS " if ok else "FAIL ") + name + (("   " + info) if info else ""))

def subsets(n):
    for r in range(n + 1):
        for S in itertools.combinations(range(n), r):
            yield S

k, d, q, A = sp.symbols('k d q A')

# ---------------------------------------------------------------- A corners
print("\n== A. Corners")
def circ(a, b, kk): return a + b + kk * a * b
for n in range(2, 6):
    xs = sp.symbols(f'x0:{n}')
    c = xs[0]
    for x in xs[1:]: c = circ(c, x, k)
    cube = sum(k**(len(S) - 1) * sp.Mul(*[xs[i] for i in S]) for S in subsets(n) if S)
    check(f"A1 n={n}: chain of a∘k b = sum over non-empty corners of k^(m-1)·product",
          sp.expand(c - cube) == 0)

# A2: which rank-only corner weights W_m (W_1 = 1) come from ONE associative binary law?
W = sp.symbols('W2:7')
def chain_from_weights(xs, Ws):
    n = len(xs)
    tot = 0
    for S in subsets(n):
        if not S: continue
        m = len(S)
        w = 1 if m == 1 else Ws[m - 2]
        tot += w * sp.Mul(*[xs[i] for i in S])
    return tot
x, y, z, w_ = sp.symbols('x y z w')
c2 = lambda a, b: a + b + W[0] * a * b
lhs3 = sp.expand(c2(c2(x, y), z)); rhs3 = sp.expand(chain_from_weights([x, y, z], W))
sol3 = sp.solve(sp.Poly(lhs3 - rhs3, x, y, z).coeffs(), W[1], dict=True)
lhs4 = sp.expand(c2(c2(c2(x, y), z), w_)); rhs4 = sp.expand(chain_from_weights([x, y, z, w_], W))
sol4 = sp.solve(sp.Poly((lhs4 - rhs4).subs(W[1], W[0]**2), x, y, z, w_).coeffs(), W[2], dict=True)
check("A2 associativity forces W3 = W2^2 and W4 = W2^3: only k^(m-1) layer weights compose",
      sol3 == [{W[1]: W[0]**2}] and sol4 == [{W[2]: W[0]**3}], f"{sol3} {sol4}")
# constant pair weight A: W_m = A^C(m,2)
defect = sp.factor(sp.expand(c2(c2(x, y), z).subs(W[0], A) - chain_from_weights([x, y, z], [A, A**3])))
check("A3 a constant pair weight A (corner weight A^C(m,2)) is NOT one associative law on numbers",
      sp.simplify(defect - (A**2 - A**3) * x * y * z) == 0, f"defect = {defect}  (zero only at A = 0, 1)")

# ---------------------------------------------------------------- B edges
print("\n== B. Edges")
n = 3; xs = sp.symbols('x0:3')
lift = sp.expand(sp.Mul(*[1 + k * xi for xi in xs]))
ok = all(lift.coeff(xs[0], int(0 in S)).coeff(xs[1], int(1 in S)).coeff(xs[2], int(2 in S)) == k**len(S)
         for S in subsets(3))
check("B1 the dial k is a weight on edges: every step up the cube multiplies by k·a_i", ok)

def holonomies(edge, n):
    """edge(S, i) for i not in S (S a frozenset). Ratio of the two routes round each square."""
    out = []
    for S in subsets(n):
        S = frozenset(S)
        for i, j in itertools.combinations([t for t in range(n) if t not in S], 2):
            r1 = edge(S, i) * edge(S | {i}, j)
            r2 = edge(S, j) * edge(S | {j}, i)
            out.append(r1 / r2)
    return out
def path_products(edge, n):
    vals = []
    for perm in itertools.permutations(range(n)):
        S = frozenset(); p = 1.0
        for i in perm:
            p *= edge(S, i); S = S | {i}
        vals.append(p)
    return vals
n = 4
Wc = {frozenset(S): random.uniform(0.5, 2) for S in subsets(n)}
flat = lambda S, i: Wc[S | {i}] / Wc[S]
hol = holonomies(flat, n); pp = path_products(flat, n)
check("B2 edge weights made from corner weights are flat: every square closes, all 24 routes agree",
      max(abs(h - 1) for h in hol) < 1e-12 and max(pp) - min(pp) < 1e-12)
Er = {(frozenset(S), i): random.uniform(0.5, 2) for S in subsets(n) for i in range(n) if i not in S}
curv = lambda S, i: Er[(S, i)]
hol = holonomies(curv, n); pp = path_products(curv, n)
check("B3 random edge weights are curved: squares don't close and routes disagree",
      max(abs(h - 1) for h in hol) > 0.1 and max(pp) - min(pp) > 0.1,
      f"{len(hol)} squares, route products {min(pp):.2f}..{max(pp):.2f}")
# Flat edge weights carry no more than corner weights: rebuild W from any flat edge field
W0 = Wc[frozenset()]
rebuilt = {}
for S in subsets(n):
    T = frozenset(); p = W0
    for i in S:
        p *= flat(T, i); T = T | {i}
    rebuilt[frozenset(S)] = p
check("B4 a flat edge field is just corner weights in disguise (rebuilt exactly from one route)",
      max(abs(rebuilt[S] - Wc[S]) for S in Wc) < 1e-12)
# matrix edges
kk, dd = 1.3, -0.7
N = np.array([[kk, 1.0], [-dd, 0.0]])
a, b = 0.4, -0.9
M1, M2 = np.eye(2) + a * N, np.eye(2) + b * N
check("B5 matrix edges I + a·N with one N (the plane's Lucas weights) commute: flat",
      np.allclose(M1 @ M2, M2 @ M1))
P, Q = np.random.randn(3, 3), np.random.randn(3, 3)
g1, g2 = np.eye(3) + a * P, np.eye(3) + b * Q
check("B6 different matrices: going round a square misses by a·b·[P, Q] (the order dial's commutator)",
      np.allclose(g1 @ g2 - g2 @ g1, a * b * (P @ Q - Q @ P)))

# ---------------------------------------------------------------- C squares
print("\n== C. Squares (2-faces)")
def interactions(Wf, n):
    out = {}
    for S in subsets(n):
        S = frozenset(S)
        for i, j in itertools.combinations([t for t in range(n) if t not in S], 2):
            out[(S, i, j)] = sp.simplify(Wf(S | {i, j}) * Wf(S) / (Wf(S | {i}) * Wf(S | {j})))
    return out
n = 4; xs = sp.symbols('x0:4')
Wk = lambda S: k**len(S) * sp.Mul(*[xs[i] for i in S])
check("C1 lifted k-dial: the square interaction W(S+i+j)W(S)/(W(S+i)W(S+j)) is 1 everywhere (no pair term)",
      all(v == 1 for v in interactions(Wk, n).values()))
Aij = {(i, j): sp.Symbol(f'A{i}{j}') for i in range(n) for j in range(i + 1, n)}
Wsol = lambda S: sp.Mul(*[Aij[(i, j)] for i, j in itertools.combinations(sorted(S), 2)]) * sp.Mul(*[xs[i] for i in S])
I = interactions(Wsol, n)
check("C2 soliton corners: the square interaction is exactly the pair weight A_ij, on every square",
      all(v == Aij[(i, j)] for (S, i, j), v in I.items()))
U = [0, 1]
for m in range(2, 8): U.append(k * U[-1] - d * U[-2])
Wl = lambda S: U[len(S)] if S else 1
vals = {sp.simplify(Wl(frozenset(range(r + 2))) * Wl(frozenset(range(r))) / (Wl(frozenset(range(r + 1)))**2)) for r in range(1, 4)}
check("C3 Lucas numerator weights are not a pair weight: their square interaction changes with rank",
      len(vals) > 1, "U_(m+2)U_m/U_(m+1)^2 for m=1..3 differ")

# C4: multiparameter quantum plane. Element = dict {exponent tuple: coeff}. X_j X_i = q_ij X_i X_j (i<j).
def qp_mul(u, v, qs, n):
    out = {}
    for a_, ca in u.items():
        for b_, cb in v.items():
            f = 1
            for i in range(n):
                for j in range(i + 1, n):
                    if a_[j] and b_[i]:
                        f *= qs[(i, j)]**(a_[j] * b_[i])
            key = tuple(a_[t] + b_[t] for t in range(n))
            out[key] = sp.expand(out.get(key, 0) + f * ca * cb)
    return {kk_: vv for kk_, vv in out.items() if vv != 0}
def gen(i, n, c=1):
    e = [0] * n; e[i] = 1
    return {tuple(e): c}
one = lambda n: {tuple([0] * n): 1}
def add(u, v):
    out = dict(u)
    for kk_, vv in v.items(): out[kk_] = sp.expand(out.get(kk_, 0) + vv)
    return {kk_: vv for kk_, vv in out.items() if vv != 0}
n = 4
qs = {(i, j): sp.Symbol(f'q{i}{j}') for i in range(n) for j in range(i + 1, n)}
cs = sp.symbols('c0:4')
def ordered_lift(order):
    prod = one(n)
    for i in order:
        prod = qp_mul(prod, add(one(n), gen(i, n, k * cs[i])), qs, n)
    return prod
def expected(order):
    pos = {i: p for p, i in enumerate(order)}
    out = {}
    for S in subsets(n):
        e = tuple(1 if t in S else 0 for t in range(n))
        f = k**len(S) * sp.Mul(*[cs[i] for i in S])
        for i, j in itertools.combinations(S, 2):
            if pos[i] > pos[j]: f *= qs[(i, j)]
        out[e] = sp.expand(f)
    return out
ok = all(ordered_lift(list(o)) == expected(list(o)) for o in itertools.permutations(range(n)))
check("C4 in a quantum plane (X_j X_i = q_ij X_i X_j), a chain in any order, written in normal order, "
      "puts q_ij on each corner for each pair the order reverses (all 24 orders, n=4)", ok)
rev = ordered_lift([3, 2, 1, 0])
full = all(rev[tuple(1 if t in S else 0 for t in range(n))] ==
           sp.expand(k**len(S) * sp.Mul(*[cs[i] for i in S]) * sp.Mul(*[qs[(i, j)] for i, j in itertools.combinations(S, 2)]))
           for S in subsets(n))
check("C5 the fully reversed chain has corner weight k^m · product of ALL its pair weights: the soliton corner form", full)

# C6: real matrices with prescribed commutation (generalized Pauli, qutrits), check C4/C5 numerically
w3 = np.exp(2j * np.pi / 3)
Zm = np.diag([1, w3, w3**2]); Xm = np.roll(np.eye(3), 1, axis=0)
I3 = np.eye(3)
mats = [np.kron(Zm, I3), np.kron(Xm, Zm), np.kron(I3, Xm), np.kron(Zm @ Xm, Xm @ Zm)]
nm = len(mats); qnum = {}; okc = True
for i in range(nm):
    for j in range(i + 1, nm):
        L = mats[j] @ mats[i]; R = mats[i] @ mats[j]
        ratio = L[np.abs(R) > 1e-9] / R[np.abs(R) > 1e-9]
        okc &= np.allclose(ratio, ratio[0]); qnum[(i, j)] = ratio[0]
check("C6a four 9x9 clock-and-shift matrices commute up to cube roots of 1 (a real quantum plane)", okc,
      "q_ij in {1, w, w^2}")
kv = 0.7; cv = [0.3, -1.1, 0.8, 0.5]
Id = np.eye(9)
lhs = Id.copy()
for i in reversed(range(nm)): lhs = lhs @ (Id + kv * cv[i] * mats[i])
rhs = Id.copy() * 0
for S in subsets(nm):
    term = Id.copy()
    for i in S: term = term @ mats[i]
    f = kv**len(S) * np.prod([cv[i] for i in S]) * np.prod([qnum[(i, j)] for i, j in itertools.combinations(S, 2)])
    rhs = rhs + f * term
check("C6b with real matrices, the reversed product equals the corner sum with pair weights", np.allclose(lhs, rhs))
cm = lambda a_, b_: a_ + b_ + kv * a_ @ b_
aa, bb, cc = [kv * 0 + c_ * m_ for c_, m_ in zip(cv[:3], mats[:3])]
check("C6c and mult-add on these matrices is still associative (it is the circle operation of a ring)",
      np.allclose(cm(cm(aa, bb), cc), cm(aa, cm(bb, cc))))
chain = cv[3] * mats[3]
for i in reversed(range(3)): chain = cm(chain, cv[i] * mats[i])
check("C6d so the pair-weighted cube is (1 + k·(a4 ∘k a3 ∘k a2 ∘k a1)): one associative law, non-commuting inputs",
      np.allclose(Id + kv * chain, lhs))

# C7: KdV 3-soliton from the reversed chain with q_ij -> A_ij, X_i -> exp(eta_i)
X_, T_ = sp.symbols('x t', real=True)
ps = [sp.Rational(1, 2), sp.Rational(1), sp.Rational(7, 4)]
eta = [p * X_ - p**3 * T_ for p in ps]
Ak = {(i, j): ((ps[i] - ps[j]) / (ps[i] + ps[j]))**2 for i in range(3) for j in range(i + 1, 3)}
tau = sum(sp.Mul(*[Ak[(i, j)] for i, j in itertools.combinations(S, 2)]) * sp.exp(sum(eta[i] for i in S)) for S in subsets(3))
Dx = lambda f, m: sp.diff(f, X_, m)
bil = tau * sp.diff(tau, X_, T_) - Dx(tau, 1) * sp.diff(tau, T_) + tau * Dx(tau, 4) - 4 * Dx(tau, 1) * Dx(tau, 3) + 3 * Dx(tau, 2)**2
res = max(abs(sp.N(bil.subs({X_: xv, T_: tv}), 30)) / abs(sp.N(tau.subs({X_: xv, T_: tv})**2, 30))
          for xv, tv in [(sp.Rational(3, 10), sp.Rational(1, 10)), (sp.Rational(-6, 5), sp.Rational(2, 5)), (sp.Integer(2), sp.Rational(-3, 10))])
check("C7 put the KdV pair weights in as q_ij: the reversed chain is the 3-soliton tau and solves Hirota's KdV",
      res < 1e-20, f"relative residual {float(res):.1e}")

# ---------------------------------------------------------------- D higher faces
print("\n== D. Higher faces")
n = 4; avs = sp.symbols('a0:4'); bvs = sp.symbols('b0:4'); t = sp.Symbol('t')
facesum = 0
for f in itertools.product([0, 1, 2], repeat=n):   # 0 out, 1 bottom corner, 2 free direction
    term = 1
    for i, s in enumerate(f):
        term *= 1 if s == 0 else (avs[i] if s == 1 else t * bvs[i])
    facesum += term
check("D1 any face weight built element by element just factors: sum over all 81 faces = prod(1 + a_i + t·b_i)",
      sp.expand(facesum - sp.Mul(*[1 + avs[i] + t * bvs[i] for i in range(n)])) == 0)
fpoly = sp.expand(sum(t**f.count(2) for f in itertools.product([0, 1, 2], repeat=5)))
check("D2 face count of the 5-cube by dimension is (2 + t)^5 (classical f-polynomial)",
      sp.expand(fpoly - (2 + t)**5) == 0)
def mobius_log(Wf, n):
    """Mobius transform of log|W|: the |T|-body interaction living on face type T."""
    phi = {}
    for T in subsets(n):
        phi[T] = sum((-1)**(len(T) - len(U)) * math.log(abs(Wf(U))) for r in range(len(T) + 1) for U in itertools.combinations(T, r))
    return phi
def minors(M):
    return lambda S: 1.0 if not S else float(np.linalg.det(M[np.ix_(S, S)]))
pv = np.array([0.5, 1.0, 1.7, 2.6]); qv = np.array([3.1, 0.2, 1.3, 0.8])
Cauchy = np.array([[1.0 / (pv[i] + qv[j]) for j in range(4)] for i in range(4)])
phiC = mobius_log(minors(Cauchy), 4)
G = np.random.uniform(0.5, 1.5, (4, 4)) + 3 * np.eye(4)
phiG = mobius_log(minors(G), 4)
check("D3 principal minors of a Cauchy matrix: interactions stop at pairs (3- and 4-body terms vanish)",
      max(abs(v) for T, v in phiC.items() if len(T) >= 3) < 1e-12)
check("D4 a generic matrix's minors have 3- and 4-body interactions (that is the sibling corner-weights stream's ground)",
      max(abs(v) for T, v in phiG.items() if len(T) >= 3) > 1e-3,
      f"max 3+-body |phi| = {max(abs(v) for T, v in phiG.items() if len(T) >= 3):.3f}")

# ---------------------------------------------------------------- E paths
print("\n== E. Paths (maximal chains = orderings)")
def inv(p): return sum(1 for i in range(len(p)) for j in range(i + 1, len(p)) if p[i] > p[j])
def qint(m, qq): return sum(qq**i for i in range(m))
for n in range(1, 7):
    lhs = sp.expand(sum(q**inv(p) for p in itertools.permutations(range(n))))
    rhs = sp.expand(sp.Mul(*[qint(m, q) for m in range(1, n + 1)]))
    check(f"E1 n={n}: weighting each route through the cube by q^(inversions) sums to the q-factorial [n]_q!",
          lhs == rhs)
w = np.exp(2j * np.pi / 3)
vals = [abs(sum(w**inv(p) for p in itertools.permutations(range(n)))) for n in range(1, 6)]
check("E2 at q a cube root of 1 the routes cancel from n = 3 on", vals[0] > .5 and vals[1] > .5 and max(vals[2:]) < 1e-9,
      ", ".join(f"{v:.3g}" for v in vals))

# ---------------------------------------------------------------- F orientations / Koszul
print("\n== F. Orientations: Koszul complex and d^N = 0")
def koszul(n, qq, xv):
    idx = {S: r for r, S in enumerate(subsets(n))}
    D = np.zeros((2**n, 2**n), dtype=complex)
    for S in subsets(n):
        for i in S:
            T = tuple(s for s in S if s != i)
            D[idx[T], idx[S]] += qq**sum(1 for j in S if j < i) * xv[i]
    return D
n = 6; xv = np.random.uniform(0.5, 1.5, n)
D = koszul(n, -1, xv)
check("F1 orientation signs on edges (holonomy -1 round every square) give the Koszul complex: d^2 = 0",
      np.abs(D @ D).max() < 1e-12)
okF = True; info = []
for Nn in range(2, 6):
    Dq = koszul(n, np.exp(2j * np.pi / Nn), xv)
    powN = np.linalg.matrix_power(Dq, Nn); powN1 = np.linalg.matrix_power(Dq, Nn - 1)
    okF &= np.abs(powN).max() < 1e-9 and np.abs(powN1).max() > 1e-3
    info.append(f"N={Nn}")
check("F2 holonomy q = e^(2πi/N) round every square: d^N = 0 but d^(N-1) ≠ 0 (an N-complex)", okF, ", ".join(info))
Dg = koszul(n, 0.37 + 0.2j, xv)
check("F3 for a generic q no power of d below 2^n vanishes early (d^6 ≠ 0 on the 6-cube)",
      np.abs(np.linalg.matrix_power(Dg, 6)).max() > 1e-6)

# ---------------------------------------------------------------- G Mobius with the dial
print("\n== G. The Mobius function with the dial")
def stirling2(n_, m_):
    return sum((-1)**j * math.comb(m_, j) * (m_ - j)**n_ for j in range(m_ + 1)) // math.factorial(m_)
def mu_k_formula(r, kk_):
    return sum((-kk_)**m * math.factorial(m) * stirling2(r, m) for m in range(r + 1))
def zeta_k_inverse(n, kk_):
    """Invert zeta_k = delta + k*eta on the subsets of n (exact fractions, back-substitution)."""
    L = list(subsets(n))
    # solve sum_{S<=U<=T} zeta_k(S,U) mu(U,T) = delta(S,T), from the top down
    res = {}
    for T in L:
        Ts = set(T)
        for S in sorted([S for S in L if set(S) <= Ts], key=len, reverse=True):
            if S == T: res[(S, T)] = Fraction(1); continue
            sm = sum(kk_ * res[(U, T)] for U in L if set(S) < set(U) <= Ts)
            res[(S, T)] = -sm
    return res
okG = True
for kk_ in [Fraction(1), Fraction(-1), Fraction(3, 7), Fraction(2)]:
    mu = zeta_k_inverse(5, kk_)
    for (S, T), v in mu.items():
        okG &= v == mu_k_formula(len(T) - len(S), kk_)
check("G1 inverting zeta_k = 1 + k·(strict steps) on the 5-cube gives sum_m (-k)^m · m!·S(r, m) (Fubini numbers at -k)", okG)
def perm_faces(r):
    return [math.factorial(r - j) * stirling2(r, r - j) for j in range(r)]  # f_j for j = 0..r-1
okE = all(sum((-1)**(r - j) * f for j, f in enumerate(perm_faces(r))) == (-1)**r for r in range(1, 8))
okK = all(mu_k_formula(r, Fraction(3, 5)) == sum((-Fraction(3, 5))**(r - j) * f for j, f in enumerate(perm_faces(r))) for r in range(1, 8))
check("G2 μ_k(empty, full) = the permutohedron's face counts read at -k; at k = 1 Euler's formula gives (-1)^n",
      okE and okK, "faces of the hexagon: 6 corners, 6 sides, 1 face")
xx = sp.Symbol('x')
ser = sp.series(1 / (1 + k * (sp.exp(xx) - 1)), xx, 0, 7).removeO()
okS = all(sp.simplify(ser.coeff(xx, r) * math.factorial(r) - sum((-k)**m * math.factorial(m) * stirling2(r, m) for m in range(r + 1))) == 0 for r in range(7))
check("G3 their generating function is 1/(1 + k(e^x - 1)): the lift 1 + k·(...) of e^x - 1, inverted", okS)

# ---------------------------------------------------------------- H dial plane d
print("\n== H. The plane dial d as a square fold-back")
Ns = sp.Matrix([[k, 1], [-d, 0]])
check("H1 Cayley-Hamilton: N^2 = k·N - d·I, so each square's top corner folds back to an edge (k) and its bottom (-d)",
      sp.simplify(Ns**2 - (k * Ns - d * sp.eye(2))) == sp.zeros(2))

# ---------------------------------------------------------------- I q-orbit
print("\n== I. The q-orbit: geometric direction dials")
def qbinom(n_, m_, qq):
    num = sp.Mul(*[1 - qq**(n_ - i) for i in range(m_)]); den = sp.Mul(*[1 - qq**(i + 1) for i in range(m_)])
    return sp.cancel(num / den)
okI = True
for n_ in range(1, 8):
    for m_ in range(n_ + 1):
        lhs = sp.expand(sum(q**sum(S) for S in itertools.combinations(range(n_), m_)))
        okI &= sp.expand(lhs - q**(m_ * (m_ - 1) // 2) * qbinom(n_, m_, q)) == 0
check("I1 direction dials 1, q, q^2, ... (a q-orbit) make each layer q^C(m,2)·[n choose m]_q: a constant pair weight q "
      "appears (Cauchy's q-binomial theorem)", okI)

print(f"\n{sum(RESULTS)} of {len(RESULTS)} checks pass")
