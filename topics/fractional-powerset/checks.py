"""Checks for the fractional power set write-up (Nov 25, 2026).

Pure Python (fractions only). Run:  python3 checks.py
Each check prints PASS/FAIL with a short note.
"""
from fractions import Fraction as Fr
from math import comb, factorial, log2
import itertools, random

results = []
def check(name, ok, note=""):
    results.append(ok)
    print(("PASS " if ok else "FAIL ") + name + (("   " + note) if note else ""))

N = 9  # series order

# ---------- truncated power series helpers (lists of Fractions, index = power) ----------
def mul(a, b):
    c = [Fr(0)] * (N + 1)
    for i, x in enumerate(a):
        if x == 0: continue
        for j, y in enumerate(b):
            if i + j > N: break
            c[i + j] += x * y
    return c

def compose(f, g):  # f(g(x)), g(0)=0
    res = [Fr(0)] * (N + 1)
    p = [Fr(1)] + [Fr(0)] * N
    for k in range(N + 1):
        if f[k] != 0:
            res = [r + f[k] * q for r, q in zip(res, p)]
        p = mul(p, g)
    return res

E = [Fr(0)] + [Fr(1, factorial(k)) for k in range(1, N + 1)]  # e^x - 1

# 1. Half iterate of e^x - 1 by solving a(a(x)) = E(x) coefficient by coefficient
a = [Fr(0), Fr(1)] + [Fr(0)] * (N - 1)
for k in range(2, N + 1):
    # coefficient k of a(a(x)) is 2*a_k + (terms from lower coefficients)
    a[k] = Fr(0)
    lower = compose(a, a)[k]
    a[k] = (E[k] - lower) / 2
aa = compose(a, a)
check("half iterate a(a(x)) = e^x - 1 to order 9", aa == E)
check("half iterate starts x + x^2/4 + x^3/48 + 0x^4 + x^5/3840 - 7x^6/92160",
      a[1:7] == [Fr(1), Fr(1, 4), Fr(1, 48), Fr(0), Fr(1, 3840), Fr(-7, 92160)],
      " ".join(str(c) for c in a[1:8]))
counts = [a[k] * factorial(k) for k in range(1, 8)]
check("half step of 'nonempty sets' is not a counting series: n!*coef is not a whole number at n=2 (1/2) and goes negative at n=6",
      counts[1] == Fr(1, 2) and counts[5] < 0, "n!*coef: " + ", ".join(str(c) for c in counts))

# 2. The rung flow E^t: coefficients are polynomials in t; n!*coef is integral for integer t only
def iterate_int(f, t):
    r = [Fr(0), Fr(1)] + [Fr(0)] * (N - 1)
    for _ in range(t): r = compose(f, r)
    return r
def lagrange_eval(xs, ys, t):
    s = Fr(0)
    for i, (xi, yi) in enumerate(zip(xs, ys)):
        term = Fr(yi)
        for j, xj in enumerate(xs):
            if j != i: term *= Fr(t - xj, xi - xj)
        s += term
    return s
ts = list(range(0, N + 1))
its = [iterate_int(E, t) for t in ts]
half_from_poly = [lagrange_eval(ts, [it[k] for it in its], Fr(1, 2)) for k in range(N + 1)]
check("E^t coefficients interpolated as polynomials in t give the half iterate at t=1/2", half_from_poly == a)
c2 = [lagrange_eval(ts, [it[2] for it in its], Fr(t)) * 2 for t in [Fr(1,3), Fr(1,2), Fr(2), Fr(-1)]]
check("2!*[x^2]E^t = t, so whole-number counts force t to be a whole number", c2 == [Fr(1,3), Fr(1,2), Fr(2), Fr(-1)])
integral = all(all((it[k] * factorial(k)).denominator == 1 and it[k] >= 0 for k in range(N + 1)) for it in its)
check("integer t >= 0: every n!*coef of E^t is a whole number >= 0 (counts of nested set partitions)", integral)

# 3. Finite wall: no non-decreasing f: N -> N with f(f(n)) = 2^n (only n = 0, 1 are needed)
def wall_case(f0, f1):
    # try to realise f(0)=f0, f(1)=f1 consistently with f(f(0))=1, f(f(1))=2 and monotonicity
    f = {0: f0, 1: f1}
    def setv(k, v):
        if k in f and f[k] != v: return False
        f[k] = v; return True
    if not setv(f0, 1): return False
    if not setv(f[1], 2): return False
    ks = sorted(f)
    return all(f[ks[i]] <= f[ks[i + 1]] for i in range(len(ks) - 1))
ok = not any(wall_case(x, y) for x in range(60) for y in range(60))
check("no non-decreasing f with f(f(n)) = 2^n: every choice of f(0), f(1) below 60 fails at n = 0 or 1", ok)

# 4. ... but starting at n = 2 a strictly increasing whole-number half step exists (greedy build)
def build_half(limit):
    f = {2: 3, 3: 4}
    n = 2
    changed = True
    # propagate f(f(n)) = 2^n and f(2^n) = 2^{f(n)}, filling gaps greedily
    for n in range(2, limit):
        if n not in f:
            prev = max(k for k in f if k < n)
            f[n] = f[prev] + 1
            # ensure we stay below the next known value
            nxt = [k for k in f if k > n]
            if nxt and f[n] >= f[min(nxt)]: return None
        v = f[n]
        target = 2 ** n
        if v in f and f[v] != target: return None
        f[v] = target
    return f
f = build_half(400)
ok = f is not None
if ok:
    ks = sorted(k for k in f if k < 400)
    ok = all(f[ks[i]] < f[ks[i + 1]] for i in range(len(ks) - 1)) and all(f[f[n]] == 2 ** n for n in range(2, 60))
check("strictly increasing whole-number f on n >= 2 with f(f(n)) = 2^n (checked n = 2..59)", ok,
      "f(2..12) = " + ", ".join(str(f[n]) for n in range(2, 13)) if f else "")

# 5. f(f(n)) = n + d has a strictly increasing solution exactly when d is even (then f(n) = n + d/2)
def has_solution(d, L=40):
    # search increasing f on 0..L with values bounded by L + d, depth-first
    sys_vals = {}
    def consistent(f):
        for n, v in f.items():
            if v in f and f[v] != n + d: return False
        ks = sorted(f)
        return all(f[ks[i]] < f[ks[i+1]] for i in range(len(ks)-1))
    def dfs(n, f):
        if n > 3 * d + 2: return True
        lo = f[n - 1] + 1 if n - 1 in f else 0
        for v in range(lo, n + d + 1):
            g = dict(f); g[n] = v
            if v in g and g[v] != n + d: continue
            if v <= n + d:
                g.setdefault(v, n + d)
                if g[v] != n + d: continue
            if consistent(g) and dfs(n + 1, g): return True
        return False
    return dfs(0, {})
sol = {d: has_solution(d) for d in range(1, 9)}
check("f(f(n)) = n + d: strictly increasing solution found exactly for even d (d = 1..8)",
      all(sol[d] == (d % 2 == 0) for d in sol), str(sol))

# 6. Ordered Mostowski model: subsets of the atoms with least support of size n number 2*3^n
def patterns(n):
    cnt = 0
    for bits in itertools.product([0, 1], repeat=2 * n + 1):
        # slots: I0, p1, I1, p2, ..., pn, In ; point p_i at index 2i-1
        if all(not (bits[2*i-2] == bits[2*i-1] == bits[2*i]) for i in range(1, n + 1)):
            cnt += 1
    return cnt
pc = [patterns(n) for n in range(0, 9)]
check("ordered Mostowski model: subsets of A with least support of size n = 2*3^n (brute force n<=8)",
      pc == [2 * 3 ** n for n in range(9)], str(pc))
def mobius(o, n):  # least-support-exactly-n count from 'supported by s points' count 2^{o(s)}
    return sum((-1) ** (n - s) * comb(n, s) * 2 ** o(s) for s in range(n + 1))
check("same count by inclusion-exclusion over supports: sum (-1)^(n-s) C(n,s) 2^(2s+1) = 2*3^n",
      all(mobius(lambda s: 2 * s + 1, n) == 2 * 3 ** n for n in range(30)))

# 7. orbit count of the stabiliser of s atoms on increasing pairs: C(2s+1,2) + (s+1); brute force on types
def pair_orbits(s):
    # slots 0..2s: even = open interval, odd = point. A pair u<v: slot(u) <= slot(v), same slot only if interval
    return sum(1 for i in range(2*s+1) for j in range(i, 2*s+1) if i < j or i % 2 == 0)
check("stabiliser of s atoms has C(2s+1,2) + (s+1) orbits on increasing pairs",
      all(pair_orbits(s) == comb(2*s+1, 2) + s + 1 for s in range(20)))
g2 = [mobius(lambda s: comb(2*s+1, 2) + s + 1, n) for n in range(12)]
ratios = [log2(g2[n]) / (2 * n * n) for n in range(4, 12)]
check("P of increasing pairs: least-support-n count grows like 2^(2n^2) (ratio falling toward 1) (power set turns orbit growth o(n) into about 2^o(n))",
      all(x > 0 for x in g2) and abs(ratios[-1] - 1) < 0.25, "log2(count)/(2n^2): " + ", ".join(f"{r:.3f}" for r in ratios))
# Halbeisen-Shelah ordering in the Mostowski model, read off as coefficient sequences
fin_c = [1] * 12; P_c = [2 * 3 ** n for n in range(12)]; seq11_c = [factorial(n) for n in range(12)]
check("coefficient picture matches m < fin(m) < 2^m < seq^{1-1}(m) eventually (1 < 2*3^n < n!)",
      all(fin_c[n] < P_c[n] for n in range(12)) and all(P_c[n] < seq11_c[n] for n in range(10, 12)))

# 8. ZFC obstruction, the key step checked exhaustively on a small well-ordered scale {0..6} plus a cap 7:
# if g is monotone and F = g o g is strictly inflationary (F(x) > x), then g is strictly inflationary
# and F(0) > g(0). Iterating this down the halvings gives a strictly decreasing chain of sizes.
def monotone_maps(n, top):
    def rec(i, lo, acc):
        if i == n: yield acc; return
        for v in range(lo, top + 1): yield from rec(i + 1, v, acc + [v])
    yield from rec(0, 0, [])
ok = True; tested = 0
for gl in monotone_maps(7, 7):
    g = gl + [7]                      # cap is absorbing
    if any(g[x] < x for x in range(8)) is None: pass
    F = [g[g[x]] for x in range(8)]
    if all(F[x] > x for x in range(7)):
        tested += 1
        if not (all(g[x] > x for x in range(7)) and (g[0] == 7 or F[0] > g[0])): ok = False
check("key step: monotone g with g o g strictly inflationary is itself strictly inflationary and g(x) < (g o g)(x); exhaustive on a 7-step scale",
      ok and tested > 0, f"{tested} monotone square roots tested")

# 9. Gap worlds 2^l = l^(+d): a monotone square root of P on aleph indices exists iff d is even
check("gap worlds: square root of 'jump d successors' on indices exists iff d even (same as check 5)",
      sol[2] and sol[4] and not sol[1] and not sol[3])

# 10. Rung flatness chain (finite sanity): for m >= 2, m*m >= m+m >= m+1, so flat multiplication => flat addition => flat successor
check("m*m >= 2m >= m+1 for m >= 2 (why flat multiplication forces flat addition forces flat successor)",
      all(m * m >= 2 * m >= m + 1 for m in range(2, 10000)))
check("The circle operation m+m+m*m sits above m*m, so the flat circle operation <=> flat multiplication (given m <= m*m)",
      all(m + m + m * m >= m * m >= m for m in range(1, 10000)))

print(f"\n{sum(results)} of {len(results)} checks pass")
