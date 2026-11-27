"""Two-base letter product: Hoffman-Ihara's q letter product with q a jet 1 + N, N^2 = kN - d.  (DRAFT checks)
Run: python3 two_base_checks.py   (sympy, exact)

Asked by Brennan, Nov 27: "Yes" (build the two-base letter product), "Because you could hold one base constant",
"At 1/2 perhaps?".
"""
import json, random, itertools, sympy as sp
from fractions import Fraction as Fr
random.seed(11)
res = []
def ok(name, c, note=""):
    res.append({"id": name, "pass": bool(c), "note": note}); print(("PASS " if c else "FAIL ") + name)

# ---------- jets a + bN over Q, N^2 = kN - d ----------
class Jet:
    def __init__(s, a, b, k, d): s.a, s.b, s.k, s.d = a, b, k, d
    def __add__(s, o): o = s.lift(o); return Jet(s.a + o.a, s.b + o.b, s.k, s.d)
    __radd__ = __add__
    def __neg__(s): return Jet(-s.a, -s.b, s.k, s.d)
    def __sub__(s, o): return s + (-s.lift(o))
    def __rsub__(s, o): return s.lift(o) - s
    def __mul__(s, o):
        o = s.lift(o)
        return Jet(s.a*o.a - s.d*s.b*o.b, s.a*o.b + s.b*o.a + s.k*s.b*o.b, s.k, s.d)
    __rmul__ = __mul__
    def lift(s, o): return o if isinstance(o, Jet) else Jet(o, 0, s.k, s.d)
    def norm(s): return s.a*s.a + s.k*s.a*s.b + s.d*s.b*s.b      # (a + b alpha)(a + b beta)
    def inv(s):
        nm = Fr(s.norm()); assert nm != 0
        return Jet((s.a + s.k*s.b)/nm, -s.b/nm, s.k, s.d)          # conjugate: a + b(k - N)
    def __truediv__(s, o): return s * s.lift(o).inv()
    def __pow__(s, e):
        r = Jet(1, 0, s.k, s.d)
        for _ in range(e): r = r * s
        return r
    def __eq__(s, o): o = s.lift(o); return s.a == o.a and s.b == o.b
    def shadow(s, lam): return s.a + s.b*lam                         # ring hom N -> eigenvalue lam
    def __repr__(s): return f"({s.a}) + ({s.b})N"

def qint(Nn, q):
    one = q.lift(1) if isinstance(q, Jet) else 1
    tot = one*0
    for j in range(Nn): tot = tot + q**j
    return tot
def letter(s, n, q):          # Bradley's q-MZV letter f_s(n) = q^{(s-1)n} / [n]_q^s, here with q a jet
    return q**((s-1)*n) / qint(n, q)**s
def H(word, M, q):            # truncated q-MZV: sum over M >= n1 > n2 > ... > nr >= 1
    tot = 0
    for ns in itertools.combinations(range(M, 0, -1), len(word)):
        t = 1
        for s, n in zip(word, ns): t = t * letter(s, n, q)
        tot = tot + t
    return tot

def rand_dial(real=True):
    while True:
        a = Fr(random.randint(-9, 9), random.randint(2, 7)); b = Fr(random.randint(-9, 9), random.randint(2, 7))
        if a != b and abs(1 + a) not in (0, 1) and abs(1 + b) not in (0, 1): return a, b

# T1 jet letters obey z_i <> z_j = z_(i+j) - N z_(i+j-1): Hoffman-Ihara's (1 - q) is now the jet -N
T1 = True
for _ in range(4):
    al, be = rand_dial(); k, d = al + be, al*be
    q = Jet(1, 1, k, d); Nj = Jet(0, 1, k, d)
    for n in range(1, 5):
        for i in range(1, 4):
            for j in range(1, 4):
                T1 &= letter(i, n, q)*letter(j, n, q) == letter(i+j, n, q) - Nj*letter(i+j-1, n, q)
ok("T1 jet letters q^((s-1)n)/[n]_q^s with q = 1 + N satisfy z_i<>z_j = z_(i+j) - N z_(i+j-1) exactly (4 random dials)", T1)

# T2 the letter product is commutative and associative; the triple product carries d as a plain number:
#    (z_i<>z_j)<>z_l = z_S - 2N z_(S-1) + N^2 z_(S-2) = z_S - 2N z_(S-1) + kN z_(S-2) - d z_(S-2)
k_, d_, Nn = sp.symbols("k d N")
def dia(u, v):                 # u, v: dict index -> coefficient (polys in N); extended bilinearly
    out = {}
    for i, a in u.items():
        for j, b in v.items():
            out[i+j] = sp.expand(out.get(i+j, 0) + a*b); out[i+j-1] = sp.expand(out.get(i+j-1, 0) - Nn*a*b)
    return out
def red(u): return {s: sp.expand(sp.rem(sp.Poly(c, Nn), sp.Poly(Nn**2 - k_*Nn + d_, Nn)).as_expr()) for s, c in u.items()}
T2 = True
for i, j, l in itertools.product(range(1, 4), repeat=3):
    zi, zj, zl = {i: 1}, {j: 1}, {l: 1}
    T2 &= red(dia(dia(zi, zj), zl)) == red(dia(zi, dia(zj, zl))) and dia(zi, zj) == dia(zj, zi)
trip = red(dia(dia({1: 1}, {1: 1}), {1: 1}))
T2b = sp.expand(trip[1] - (k_*Nn - d_)) == 0 and sp.expand(trip[2] + 2*Nn) == 0 and trip[3] == 1
ok("T2 the jet letter product is commutative and associative (i, j, l = 1..3)", T2)
ok("T2b triple product z1<>z1<>z1 = z3 - 2N z2 + (kN - d) z1: d enters as a number at third order", T2b,
   "matches the plane law, whose first d term is the x^2 y coefficient (Q3)")

# T3 truncated harmonic products over the jet ring (quasi-shuffle with letter product <>)
def qsh(u, v):                 # quasi-shuffle of words; coefficients are polynomials in N
    if not u: return {v: sp.Integer(1)}
    if not v: return {u: sp.Integer(1)}
    out = {}
    def add(w, c): out[w] = sp.expand(out.get(w, 0) + c)
    for w, c in qsh(u[1:], v).items(): add((u[0],) + w, c)
    for w, c in qsh(u, v[1:]).items(): add((v[0],) + w, c)
    for w, c in qsh(u[1:], v[1:]).items():
        add((u[0] + v[0],) + w, c); add((u[0] + v[0] - 1,) + w, -Nn*c)
    return out
def Heval(expr, M, q, Nj):
    tot = 0
    for w, c in expr.items():
        P = sp.Poly(c, Nn)
        for (p,), co in P.terms(): tot = tot + int(co) * Nj**p * H(w, M, q)
    return tot
T3 = True; ex = None
al, be = Fr(1, 3), Fr(-3, 5); k, d = al + be, al*be
q = Jet(1, 1, k, d); Nj = Jet(0, 1, k, d)
for M in range(1, 7):
    lhs = H((2,), M, q) * H((2,), M, q)
    rhs = 2*H((2, 2), M, q) + H((4,), M, q) - Nj*H((3,), M, q)
    T3 &= lhs == rhs
    T3 &= H((2,), M, q)*H((3,), M, q) == H((2, 3), M, q) + H((3, 2), M, q) + H((5,), M, q) - Nj*H((4,), M, q)
    T3 &= H((2,), M, q)*H((2, 1), M, q) == Heval(qsh((2,), (2, 1)), M, q, Nj)
    T3 &= H((2, 1), M, q)*H((3,), M, q) == Heval(qsh((2, 1), (3,)), M, q, Nj)
ok("T3 truncated sums multiply by the jet quasi-shuffle, e.g. H(2)H(2) = 2H(2,2) + H(4) - N H(3) (M = 1..6, depth up to 2x1)", T3)

# T4 shadows: N -> alpha and N -> beta are ring maps; they send the jet q-MZV to Bradley's at q1 = 1 + alpha and
#    q2 = 1 + beta, and the dial is (k, d) = (q1 + q2 - 2, (q1 - 1)(q2 - 1))
def Hs(word, M, qq):
    tot = 0
    for ns in itertools.combinations(range(M, 0, -1), len(word)):
        t = 1
        for s, n in zip(word, ns): t *= qq**((s-1)*n) / sum(qq**j for j in range(n))**s
        tot += t
    return tot
T4 = True
for _ in range(3):
    al, be = rand_dial(); k, d = al + be, al*be; q = Jet(1, 1, k, d)
    for w in [(2,), (3,), (2, 1), (3, 1, 2)]:
        J = H(w, 5, q)
        T4 &= J.shadow(al) == Hs(w, 5, 1 + al) and J.shadow(be) == Hs(w, 5, 1 + be)
    q1, q2 = 1 + al, 1 + be
    T4 &= k == q1 + q2 - 2 and d == (q1 - 1)*(q2 - 1)
ok("T4 the jet q-MZV is the pair (Bradley at q1 = 1 + alpha, Bradley at q2 = 1 + beta); k = q1 + q2 - 2, d = (q1 - 1)(q2 - 1)", T4,
   "for alpha != beta real, Q[N]/(N^2 - kN + d) = Q x Q, so off the border the jet carries exactly two q-MZVs")

# T5 hold one base at 1 (q2 = 1, beta = 0, d = 0): the second shadow is the ordinary truncated MZV (letters 1/n^s)
T5 = True
for kk in [Fr(1, 3), Fr(-2, 5), Fr(7, 4)]:
    q = Jet(1, 1, kk, 0)
    for w in [(2,), (2, 1), (3, 2)]:
        J = H(w, 5, q)
        plain = sum(Fr(1, 1) * __import__('math').prod([Fr(1, n**s) for s, n in zip(w, ns)]) for ns in itertools.combinations(range(5, 0, -1), len(w)))
        T5 &= J.shadow(0) == plain and J.shadow(kk) == Hs(w, 5, 1 + kk)
ok("T5 holding one base at 1 (d = 0): the jet q-MZV = (Bradley's q-MZV at q = 1 + k, the ordinary MZV)", T5)

# T6 holding one base at a constant c: the dial points form the line d = beta k - beta^2 (beta = c - 1),
#    tangent to the border d = k^2/4 at k = 2 beta; the border is the envelope of these lines
b_, kk_ = sp.symbols("b kk")
line = b_*kk_ - b_**2
T6 = sp.simplify(kk_**2/4 - line - (kk_/2 - b_)**2) == 0      # border - line = square, zero only at kk = 2b
ok("T6 holding one base at c gives the line d = (c-1)k - (c-1)^2, tangent to the border d = k^2/4 at k = 2(c-1)", T6,
   "c = 1: the line d = 0, tangent at the origin. c = 1/2: d = -k/2 - 1/4, tangent at (k, d) = (-1, 1/4)")

# T7 at c = 1/2: the second shadow is Bradley's q-MZV at q = 1/2, whose letters are 2^(n-s)/(2^n - 1)^s, and its
#    depth-1 weight-2 value is (1/4) * zeta_q(2) of Zudilin 2002 at q = 1/2, i.e. (1/4) sum sigma(m)/2^m
T7 = all(Fr(1, 2)**((s-1)*n) / sum(Fr(1, 2)**j for j in range(n))**s == Fr(2)**(n-s)/(2**n - 1)**s
         for n in range(1, 8) for s in range(1, 5))
M = 60
brad = sum(Fr(1, 2)**n / sum(Fr(1, 2)**j for j in range(n))**2 for n in range(1, M))
zud = sum(int(sp.divisor_sigma(m, 1)) * Fr(1, 2)**m for m in range(1, 2*M))
T7b = abs(float(brad - zud/4)) < 1e-15
EB = sum(Fr(1, 2**n - 1) for n in range(1, 80)); EB2 = sum(int(sp.divisor_sigma(m, 0))*Fr(1, 2)**m for m in range(1, 160))
T7c = abs(float(EB - EB2)) < 1e-20
ok("T7 at q2 = 1/2 Bradley's letters are 2^(n-s)/(2^n - 1)^s", T7)
ok("T7b at q = 1/2, Bradley's zeta_q(2) = (1/4) sum sigma(m)/2^m = (1/4) Zudilin's zeta_q(2): %.15f" % float(brad), T7b,
   "Zudilin 2002 credits Duverney with its irrationality for q = 1/p")
ok("T7c the other normalization at q = 1/2: sum 1/(2^n - 1) = sum d(m)/2^m, the Erdos-Borwein constant %.12f" % float(EB), T7c)

# T8 on the border alpha = beta (d = k^2/4) the jet ring is the dual numbers: a + bN = (a + b alpha) + b eps,
#    so the jet q-MZV is (Bradley's value, its q-derivative) at q = 1 + alpha: a literal first-order jet in q
T8 = True; qs = sp.symbols("q")
def Hsym(word, M):
    tot = 0
    for ns in itertools.combinations(range(M, 0, -1), len(word)):
        t = 1
        for s, n in zip(word, ns): t *= qs**((s-1)*n) / sum(qs**j for j in range(n))**s
        tot += t
    return tot
for al in [Fr(1, 3), Fr(-1, 4), Fr(3, 2)]:
    k, d = 2*al, al**2; q = Jet(1, 1, k, d)
    for w in [(2,), (2, 1), (3, 1)]:
        J = H(w, 4, q); f = Hsym(w, 4)
        T8 &= sp.Rational(J.a + J.b*al) == f.subs(qs, sp.Rational(1 + al)) and sp.Rational(J.b) == sp.diff(f, qs).subs(qs, sp.Rational(1 + al))
ok("T8 on the border d = k^2/4 the jet q-MZV is (zeta_q, d zeta_q/dq) at q = 1 + k/2: the value and its derivative", T8)

# T9 complex pair (d > k^2/4): N acts as a complex eigenvalue, the jet q-MZV is Bradley's at a complex q
T9 = True
k, d = 0, 1; q = Jet(1, 1, k, d)            # N^2 = -1, N -> i, q -> 1 + i
for w in [(2,), (2, 1)]:
    J = H(w, 4, q); f = Hsym(w, 4).subs(qs, 1 + sp.I)
    T9 &= sp.simplify(sp.expand_complex(sp.Rational(J.a) + sp.I*sp.Rational(J.b) - f)) == 0
ok("T9 complex pair (k, d) = (0, 1): the jet q-MZV is Bradley's at q = 1 + i (real part a, imaginary part b)", T9)

# T10 negative control: with d != 0 the jet letters do NOT obey the scalar product with -k in place of -N
al, be = Fr(1, 3), Fr(-3, 5); k, d = al + be, al*be; q = Jet(1, 1, k, d)
T10 = not (letter(1, 3, q)*letter(1, 3, q) == letter(2, 3, q) - k*letter(1, 3, q))
ok("T10 negative control: replacing the jet -N by the number -k breaks the letter product when d != 0", T10)

json.dump(res, open("two_base_checks.json", "w"), indent=1)
print(f"{sum(r['pass'] for r in res)}/{len(res)} pass")
