"""Boolean lattice checks for the circle operation (Nov 25, 2026).

Every check is exact (fractions or integers). Run: python3 checks.py
"""
from fractions import Fraction as Fr
from itertools import combinations, permutations
from math import factorial, prod
import random

random.seed(1)
ok = 0
def check(name, cond):
    global ok
    assert cond, name
    ok += 1
    print("PASS", name)

def circ(k, a, b):            # the circle operation with dial k
    return a + b + k * a * b

def plane(k, d, a, b):        # the two-dial law (Coleman-McGuinness), a, b fractions
    return (a + b + k * a * b) / (1 - d * a * b)

def subsets(n):
    for m in range(n + 1):
        for S in combinations(range(n), m):
            yield S

def chain(op, xs):
    t = xs[0]
    for x in xs[1:]:
        t = op(t, x)
    return t

def rnd():
    return Fr(random.randint(-9, 9), random.randint(1, 7))

# 1. A chain is a sum over the non-empty corners of the cube, weighted k^(rank-1)
for n in range(2, 6):
    for k in [Fr(1), Fr(0), Fr(-1), Fr(1, 2), Fr(3)]:
        a = [rnd() for _ in range(n)]
        lhs = chain(lambda x, y: circ(k, x, y), a)
        rhs = sum(k ** (len(S) - 1) * prod(a[i] for i in S) for S in subsets(n) if S)
        check(f"chain of {n} = sum over corners S of k^(|S|-1) prod a_S, k={k}", lhs == rhs)

# 2. Inclusion-exclusion is the k = -1 chain (independent events)
p = [Fr(1, 10), Fr(1, 5), Fr(1, 4), Fr(1, 3)]
union = 1 - prod(1 - x for x in p)
check("P(at least one) = chain at k = -1 (inclusion-exclusion)", chain(lambda x, y: circ(-1, x, y), p) == union)

# 3. Noncommutative: ordered corners. 2x2 integer matrices.
def mm(A, B):
    return [[sum(A[i][t] * B[t][j] for t in range(2)) for j in range(2)] for i in range(2)]
def madd(A, B):
    return [[A[i][j] + B[i][j] for j in range(2)] for i in range(2)]
def msc(c, A):
    return [[c * A[i][j] for j in range(2)] for i in range(2)]
Z = [[0, 0], [0, 0]]
def mcirc(k, A, B):
    return madd(madd(A, B), msc(k, mm(A, B)))
for n in [3, 4]:
    for k in [1, 2, -1]:
        As = [[[random.randint(-3, 3) for _ in range(2)] for _ in range(2)] for _ in range(n)]
        lhs = chain(lambda X, Y: mcirc(k, X, Y), As)
        rhs = Z
        for S in subsets(n):
            if not S:
                continue
            P = As[S[0]]
            for i in S[1:]:
                P = mm(P, As[i])          # factors kept in chain order
            rhs = madd(rhs, msc(k ** (len(S) - 1), P))
        check(f"matrix chain of {n} = sum over corners, factors in chain order, k={k}", lhs == rhs)
A1 = [[0, 1], [0, 0]]; A2 = [[0, 0], [1, 0]]
check("order matters for matrices: A1 o A2 != A2 o A1", mcirc(1, A1, A2) != mcirc(1, A2, A1))

# 4. The two-dial plane weights corners by Lucas numbers U_m(k, d)
def lucas(k, d, m):           # U_0 = 0, U_1 = 1, U_m = k U_{m-1} - d U_{m-2}
    u0, u1 = Fr(0), Fr(1)
    if m == 0:
        return u0
    for _ in range(m - 1):
        u0, u1 = u1, k * u1 - d * u0
    return u1
for k, d in [(0, -1), (0, 1), (1, -1), (1, 0), (2, 3), (Fr(1, 2), Fr(-2, 3))]:
    k, d = Fr(k), Fr(d)
    for n in range(2, 6):
        a = [Fr(random.randint(-3, 3), random.randint(5, 11)) for _ in range(n)]
        lhs = chain(lambda x, y: plane(k, d, x, y), a)
        num = sum(lucas(k, d, len(S)) * prod(a[i] for i in S) for S in subsets(n) if S)
        den = 1 - d * sum(lucas(k, d, len(S) - 1) * prod(a[i] for i in S) for S in subsets(n) if len(S) >= 2)
        check(f"plane chain of {n} = Lucas-weighted corners, (k, d) = ({k}, {d})", lhs == num / den)
check("(k,d) = (0,-1) weights 0,1,0,1: velocity addition = odd corners / even corners",
      [lucas(0, -1, m) for m in range(6)] == [0, 1, 0, 1, 0, 1])
check("(k,d) = (0,1) weights 0,1,0,-1: the tangent of a sum of angles",
      [lucas(0, 1, m) for m in range(6)] == [0, 1, 0, -1, 0, 1])
check("(k,d) = (1,-1) weights are the Fibonacci numbers",
      [lucas(1, -1, m) for m in range(8)] == [0, 1, 1, 2, 3, 5, 8, 13])
check("(k,d) = (k,0) weights are k^(m-1), the single-dial case",
      all(lucas(Fr(3), 0, m) == 3 ** (m - 1) for m in range(1, 7)))

# 5. Incidence algebra of B_3: Moebius = circle inverse of the strict zeta; Hall's chains
n = 3
els = list(subsets(n)); idx = {S: i for i, S in enumerate(els)}; N = len(els)
def le(S, T): return set(S) <= set(T)
eta = [[1 if (le(S, T) and S != T) else 0 for T in els] for S in els]
I = [[1 if i == j else 0 for j in range(N)] for i in range(N)]
def MM(A, B): return [[sum(A[i][t] * B[t][j] for t in range(N)) for j in range(N)] for i in range(N)]
def MA(A, B): return [[A[i][j] + B[i][j] for j in range(N)] for i in range(N)]
def MS(c, A): return [[c * A[i][j] for j in range(N)] for i in range(N)]
powers = [I]
for _ in range(n + 1):
    powers.append(MM(powers[-1], eta))
check("eta^4 = 0 on B_3 (no chain longer than 3)", all(v == 0 for r in powers[4] for v in r))
for k in [1, 2, -1]:
    inv = I
    for m in range(1, n + 1):
        inv = MA(inv, MS((-k) ** m, powers[m]))          # (1 + k eta)^(-1) = sum (-k)^m eta^m
    check(f"(1 + k eta)(sum (-k)^m eta^m) = 1 on B_3, k={k}", MM(MA(I, MS(k, eta)), inv) == I)
    etap = MS(Fr(1, k), MA(inv, MS(-1, I)))              # circle inverse of eta with dial k
    zero = MA(MA(eta, etap), MS(k, MM(eta, etap)))
    check(f"eta o_k eta' = 0: the way back of the strict zeta, k={k}", all(v == 0 for r in zero for v in r))
mu = [[sum((-1) ** m * powers[m][i][j] for m in range(n + 1)) for j in range(N)] for i in range(N)]
check("Hall: mu = alternating count of chains equals (-1)^|T-S| (Rota's Boolean Moebius)",
      all(mu[idx[S]][idx[T]] == ((-1) ** (len(T) - len(S)) if le(S, T) else 0) for S in els for T in els))

# 6. Chains of B_n from bottom to top are ordered set partitions; maximal ones are orderings
def chains_bottom_top(n):
    full = tuple(range(n)); out = []
    def go(cur, path):
        if cur == full:
            out.append(path); return
        rest = [i for i in range(n) if i not in cur]
        for m in range(1, len(rest) + 1):
            for add in combinations(rest, m):
                nxt = tuple(sorted(cur + add)); go(nxt, path + [add])
    go((), []); return out
for n in range(1, 6):
    ch = chains_bottom_top(n)
    maximal = [c for c in ch if len(c) == n]
    check(f"B_{n}: maximal chains = n! = vertices of the permutohedron", len(maximal) == factorial(n))
    fubini = [1, 1, 3, 13, 75, 541][n]
    check(f"B_{n}: all chains bottom-to-top = ordered set partitions (Fubini {fubini}) = faces of the permutohedron", len(ch) == fubini)

# 7. Rank weights at equal entries: binomial layers
for n in range(1, 7):
    for k in [Fr(1), Fr(2), Fr(-1, 2)]:
        a = Fr(2, 3)
        lhs = chain(lambda x, y: circ(k, x, y), [a] * n)
        layers = sum(Fr(factorial(n), factorial(m) * factorial(n - m)) * k ** (m - 1) * a ** m for m in range(1, n + 1))
        check(f"equal entries, n={n}, k={k}: layers C(n,m) k^(m-1) a^m = ((1+ka)^n - 1)/k", lhs == layers == ((1 + k * a) ** n - 1) / k)

print(f"\n{ok} checks passed.")
