"""Checks for the k^(m-1) cube-weighting naming page (DRAFT).

Exact rational arithmetic throughout (fractions.Fraction). Run: python3 checks.py
Sections
  A  the weighting itself: fold = corners weighted k^(m-1); rank-only weights that compose
  B  Sugeno lambda-measure: Mobius transform (Harsanyi dividends) = lambda^(m-1) prod g
  C  Grothendieck: single-box stable Grothendieck polynomial (set-valued tableaux)
  D  quasi-shuffle of weight k (Hoffman; Guo-Keigher mixable shuffle; Rota-Baxter weight k)
  E  Hoffman's exp/log = Psi_f with f = (e^{kt}-1)/k and L_k = log(1+kt)/k
"""
import itertools, json, random
from fractions import Fraction as Fr
from math import factorial

random.seed(20261127)
results = []
def check(name, ok, note=""):
    results.append({"id": name, "pass": bool(ok), "note": note})
    print(("PASS " if ok else "FAIL ") + name + ("  " + note if note else ""))

def rnd():
    return Fr(random.randint(-9, 9), random.randint(1, 9))

def fold(xs, k):
    acc = Fr(0)
    for x in xs:
        acc = acc + x + k * acc * x
    return acc

def corners(xs, k):
    n = len(xs); tot = Fr(0)
    for m in range(1, n + 1):
        for S in itertools.combinations(range(n), m):
            p = Fr(1)
            for i in S: p *= xs[i]
            tot += k ** (m - 1) * p
    return tot

# ---------------- A ----------------
okA1 = all(fold(xs, k) == corners(xs, k)
           for n in range(1, 7) for _ in range(5)
           for xs in [[rnd() for _ in range(n)]] for k in [rnd()])
check("A1 fold under x+y+kxy = sum over corners S of k^(|S|-1) prod x_S, n = 1..6", okA1)
okA1b = all(1 + k * fold(xs, k) == __import__("functools").reduce(lambda a, b: a * b, [1 + k * x for x in xs], Fr(1))
            for n in range(1, 7) for xs in [[rnd() for _ in range(n)]] for k in [rnd()] if k != 0)
check("A1b same thing as (prod(1 + k x_i) - 1)/k", okA1b)

# rank-only weights W1=1, W2, W3: associativity of 3 inputs forces W3 = W2^2 (symbolic by sampling)
def rank_fold3(x, y, z, W2, W3):
    return x + y + z + W2 * (x*y + y*z + x*z) + W3 * x*y*z
def two(x, y, W2):
    return x + y + W2 * x * y
okA2 = True
for _ in range(20):
    W2 = rnd(); x, y, z = rnd(), rnd(), rnd()
    # grouping (x y) z with the two-input law gives W3 = W2^2 exactly
    lhs = two(two(x, y, W2), z, W2)
    okA2 &= lhs == rank_fold3(x, y, z, W2, W2 ** 2)
    okA2 &= (x*y*z == 0) or lhs != rank_fold3(x, y, z, W2, W2 ** 2 + 1)
check("A2 rank-only weights from one two-input law must be W_m = W_2^(m-1) (here W3 = W2^2)", okA2,
      "re-check of cube-weights A2")

# ---------------- B Sugeno ----------------
def sugeno_from_singletons(g, lam, order):
    n = len(g); val = {frozenset(): Fr(0)}
    # build each set by adding its elements in the given order with the lambda rule
    for m in range(1, n + 1):
        for S in itertools.combinations(range(n), m):
            elems = sorted(S, key=lambda i: order.index(i))
            acc = Fr(0)
            for i in elems:
                acc = acc + g[i] + lam * acc * g[i]
            val[frozenset(S)] = acc
    return val
okB1 = okB2 = True
for n in range(2, 6):
    g = [rnd() for _ in range(n)]; lam = rnd()
    v1 = sugeno_from_singletons(g, lam, list(range(n)))
    v2 = sugeno_from_singletons(g, lam, random.sample(range(n), n))
    okB1 &= v1 == v2
    # lambda rule on every disjoint pair
    for A in v1:
        for B in v1:
            if not (A & B):
                okB1 &= v1[A | B] == v1[A] + v1[B] + lam * v1[A] * v1[B]
    # Mobius inversion
    for S in v1:
        if not S: continue
        m = sum((-1) ** (len(S) - len(T)) * v1[frozenset(T)]
                for r in range(len(S) + 1) for T in itertools.combinations(sorted(S), r))
        p = Fr(1)
        for i in S: p *= g[i]
        okB2 &= m == lam ** (len(S) - 1) * p
check("B1 Sugeno lambda-measure from singletons: order-free, and A,B disjoint -> g(AuB)=gA+gB+lam gA gB", okB1)
check("B2 its Mobius transform (Harsanyi dividends) is lam^(|S|-1) prod g_i, n = 2..5", okB2)

# ---------------- C Grothendieck ----------------
# set-valued tableaux of the one-box shape = nonempty subsets S of [n]; Buch Thm 3.1 sign (-1)^(|T|-|lambda|)
def G1(xs):
    n = len(xs); tot = Fr(0)
    for m in range(1, n + 1):
        for S in itertools.combinations(range(n), m):
            p = Fr(1)
            for i in S: p *= xs[i]
            tot += (-1) ** (m - 1) * p
    return tot
okC1 = okC2 = True
for n in range(1, 7):
    xs = [rnd() for _ in range(n)]
    pr = Fr(1)
    for x in xs: pr *= 1 - x
    okC1 &= G1(xs) == 1 - pr == fold(xs, Fr(-1))
    k = rnd()
    if k != 0:
        okC2 &= fold(xs, k) == -G1([-k * x for x in xs]) / k
check("C1 single-box stable Grothendieck G_1 = 1 - prod(1 - x_i) = the union rule (k = -1 fold)", okC1)
check("C2 every k: fold_k(x) = -G_1(-k x)/k, i.e. the beta-Grothendieck single box with beta = k", okC2)

# ---------------- D quasi-shuffle ----------------
# letters: sorted tuples of generator ids (monomials); bracket = merge of tuples; words: tuples of letters
def br(a, b): return tuple(sorted(a + b))
def add(d, w, c):
    if c:
        d[w] = d.get(w, 0) + c
        if d[w] == 0: del d[w]
from functools import lru_cache
def make_qsh(lam):
    @lru_cache(maxsize=None)
    def q(u, v):
        if not u: return {v: Fr(1)}
        if not v: return {u: Fr(1)}
        out = {}
        for w, c in q(u[1:], v).items(): add(out, (u[0],) + w, c)
        for w, c in q(u, v[1:]).items(): add(out, (v[0],) + w, c)
        if lam:
            for w, c in q(u[1:], v[1:]).items(): add(out, (br(u[0], v[0]),) + w, lam * c)
        return out
    return q
def prod_elems(X, Y, q):
    out = {}
    for u, a in X.items():
        for v, b in Y.items():
            for w, c in q(u, v).items(): add(out, w, a * b * c)
    return out
def rword(maxlen=3):
    return tuple((random.randint(1, 4),) for _ in range(random.randint(1, maxlen)))

okD1 = True
for lam in [Fr(1), Fr(-1), Fr(2, 3), Fr(-5, 2)]:
    q = make_qsh(lam)
    for _ in range(6):
        u, v, w = {rword(): Fr(1)}, {rword(): Fr(1)}, {rword(): Fr(1)}
        okD1 &= prod_elems(prod_elems(u, v, q), w, q) == prod_elems(u, prod_elems(v, w, q), q)
        okD1 &= prod_elems(u, v, q) == prod_elems(v, u, q)
check("D1 weight-k quasi-shuffle (a u)*(b v) = a(u*bv) + b(au*v) + k[ab](u*v) is associative and commutative", okD1)

def compositions(n):
    if n == 0: yield (); return
    for first in range(1, n + 1):
        for rest in compositions(n - first): yield (first,) + rest
def Psi(X, c):
    """Hoffman-Ihara Psi_f: c[m] is the t^m coefficient of f."""
    out = {}
    for w, a in X.items():
        if not w: add(out, w, a); continue
        for I in compositions(len(w)):
            coef = Fr(1); pos = 0; neww = []
            for m in I:
                coef *= c(m)
                letter = ()
                for L in w[pos:pos + m]: letter = br(letter, L)
                neww.append(letter); pos += m
            add(out, tuple(neww), a * coef)
    return out

okD2 = True
for lam in [Fr(2), Fr(-1, 3), Fr(7, 4)]:
    q1, ql = make_qsh(Fr(1)), make_qsh(lam)
    lin = lambda m, lam=lam: lam if m == 1 else Fr(0)
    inv = lambda m, lam=lam: 1 / lam if m == 1 else Fr(0)
    for _ in range(6):
        u, v = {rword(): Fr(1)}, {rword(): Fr(1)}
        lhs = Psi(prod_elems(Psi(u, lin), Psi(v, lin), q1), inv)
        okD2 &= lhs == prod_elems(u, v, ql)
check("D2 weight k = the Diehl-Ebrahimi-Fard-Tapia deformation with the linear f(t) = k t (a rescaling of weight 1)", okD2)

okD3 = True; okD4 = True
for k in [Fr(1), Fr(-1), Fr(3, 2)]:
    q = make_qsh(k)
    for n in range(1, 6):
        xs = [(i + 1,) for i in range(n)]
        P = {(): Fr(1)}
        for x in xs: P = prod_elems(P, {(x,): Fr(1)}, q)
        # every term is an ordered set partition; coefficient k^(n - #blocks)
        for w, c in P.items():
            okD3 &= c == k ** (n - len(w)) and sum(len(L) for L in w) == n
        # one-letter part of prod (eps + [x_i]) equals the fold, letters evaluated at numbers
        E = {(): Fr(1)}
        for x in xs: E = prod_elems(E, {(): Fr(1), (x,): Fr(1)}, q)
        vals = {i + 1: rnd() for i in range(n)}
        one = Fr(0)
        for w, c in E.items():
            if len(w) == 1:
                p = Fr(1)
                for gen in w[0]: p *= vals[gen]
                one += c * p
        okD4 &= one == fold([vals[i + 1] for i in range(n)], k)
check("D3 [x1]*...*[xn] at weight k: each merged block of m letters carries k^(m-1) (n = 1..5)", okD3)
check("D4 one-letter part of the product of (1 + [x_i]) at weight k = the fold under x+y+kxy", okD4)

# ---------------- E exp / log ----------------
def cE(k):   # (e^{kt}-1)/k = sum k^(m-1) t^m / m!
    return lambda m: k ** (m - 1) / factorial(m)
def cL(k):   # log(1+kt)/k = sum (-k)^(m-1) t^m / m
    return lambda m: (-k) ** (m - 1) / Fr(m)
okE1 = okE2 = okE3 = True
q0 = make_qsh(Fr(0))
for k in [Fr(1), Fr(2), Fr(-1, 2)]:
    qk = make_qsh(k)
    for _ in range(5):
        u, v = {rword(): Fr(1)}, {rword(): Fr(1)}
        okE1 &= Psi(prod_elems(u, v, qk), cL(k)) == prod_elems(Psi(u, cL(k)), Psi(v, cL(k)), q0)
        okE2 &= Psi(prod_elems(u, v, q0), cE(k)) == prod_elems(Psi(u, cE(k)), Psi(v, cE(k)), qk)
        okE3 &= Psi(Psi(u, cL(k)), cE(k)) == u
check("E1 Psi_{L_k}, L_k = log(1+kt)/k, turns the weight-k quasi-shuffle into the shuffle (Hoffman's log)", okE1)
check("E2 Psi_{(e^{kt}-1)/k} turns the shuffle into the weight-k quasi-shuffle (Hoffman's exp)", okE2)
check("E3 they are inverse (Psi_f Psi_g = Psi_(f o g), Hoffman-Ihara Thm 3.1)", okE3)
# same function on numbers: L_k(x + y + kxy) = L_k x + L_k y, coefficients match
okE4 = all(cL(k)(m) == (-k) ** (m - 1) / Fr(m) for k in [Fr(2), Fr(-3)] for m in range(1, 8))
import math
okE4 &= all(abs(math.log1p(k*(x+y+k*x*y))/k - math.log1p(k*x)/k - math.log1p(k*y)/k) < 1e-12
            for k, x, y in [(0.7, 0.2, 0.3), (-0.5, 0.4, 0.1), (2.0, 0.05, 0.11)])
check("E4 the series in Hoffman's log is our flat coordinate L_k: L_k(x o_k y) = L_k x + L_k y", okE4)

# Yamamoto's interpolation map: Psi_f with f(x) = x/(1 - kx), coefficients exactly k^(m-1)
okE5 = True
for k in [Fr(1), Fr(-2), Fr(3, 5)]:
    cY = lambda m, k=k: k ** (m - 1)
    cYi = lambda m, k=k: (-k) ** (m - 1)          # x/(1 + kx), the inverse series
    for _ in range(5):
        w = rword(4); u = {w: Fr(1)}
        Y = Psi(u, cY)
        okE5 &= all(c == k ** (len(w) - len(v)) for v, c in Y.items())
        okE5 &= Psi(Y, cYi) == u
check("E5 Yamamoto's t-interpolation map is Psi_f with f = x/(1 - kx): a block of m letters gets k^(m-1); inverse f = x/(1 + kx)", okE5)

# ---------------- F negative controls (each must FAIL the identity, so the checks above can fail) ----------------
u, v = {((1,), (2,)): Fr(1)}, {((3,),): Fr(1)}
k = Fr(2)
bad1 = Psi(prod_elems(u, v, make_qsh(k)), cL(Fr(1))) != prod_elems(Psi(u, cL(Fr(1))), Psi(v, cL(Fr(1))), q0)
bad2 = Psi(prod_elems(u, v, make_qsh(Fr(1))), lambda m: Fr(1)) != prod_elems(Psi(u, lambda m: Fr(1)), Psi(v, lambda m: Fr(1)), q0)
bad3 = fold([Fr(1,2), Fr(1,3), Fr(1,5)], Fr(2)) != corners([Fr(1,2), Fr(1,3), Fr(1,5)], Fr(3))
check("F1 controls: wrong k in the log, a wrong series, and a wrong corner weight each break the identity", bad1 and bad2 and bad3)

npass = sum(r["pass"] for r in results)
print(f"\n{npass}/{len(results)} pass")
json.dump({"pass": npass, "total": len(results), "checks": results}, open("checks.json", "w"), indent=1)
