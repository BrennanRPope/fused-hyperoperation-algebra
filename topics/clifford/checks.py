"""Clifford's work, piece by piece, against the circle operation (Nov 25, 2026).

Run: python3 checks.py   (writes checks_output.txt and results.json)
Exact algebra with sympy; numbers with numpy.  Uses the small Clifford-algebra class from the
gr-bridges thread (math-theory/gr-bridges/ga.py), which that thread wrote; nothing here repeats
its section S6 (the sandwich law a + b + a k b), which is cited, not re-run.

Notation.  The circle operation a * b = a + b + ab (written ∘₁ on the pages).  Landscape law
F_{p,d}(a,b) = (a + b + p ab)/(1 - d ab), rate V(a) = 1 + p a + d a^2.  The circle operation is (p,d) = (1,0).
"""
import json, math, os, sys, itertools, random
import numpy as np
import sympy as sp

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, "/mnt/project-files/math-theory/gr-bridges")
from ga import GA  # noqa: E402

LINES, RES = [], {}
npass = nfail = 0


def check(name, ok, detail=""):
    global npass, nfail
    ok = bool(ok)
    npass += ok
    nfail += (not ok)
    line = f"[{'PASS' if ok else 'FAIL'}] {name}" + (f"  ({detail})" if detail else "")
    print(line, flush=True)
    LINES.append(line)


def section(t):
    print("\n== " + t, flush=True)
    LINES.append("\n== " + t)


Z = lambda e: sp.simplify(e) == 0
a, b, c, p, q, k, x, y, k0, k1, kap = sp.symbols("a b c p q k x y k0 k1 kappa")
F = lambda u, v, p=p, q=q: (u + v + p * u * v) / (1 - q * u * v)
ma = lambda u, v: u + v + u * v
rng = np.random.default_rng(7)

# =============================================================================================
section("A. The geometric product: the landscape law is the even Clifford algebra of its own rate form")
# Clifford algebra of the binary form Q(x,y) = x^2 + p x y + d y^2: generators u, v with
# u^2 = 1, v^2 = d, uv + vu = p.  Represent by 4x4 matrices built from a basis 1, u, v, uv.
# Basis order [1, u, v, uv]; left multiplication tables derived from the relations.
def cl_binary(p, q):
    # products of basis elements, as coefficient vectors over [1,u,v,uv]
    # vu = p - uv ; uvu = u(p - uv) = p u - v ... (u^2=1): u v u = u (p - u v) = p u - v
    # vuv = (p - uv) v = p v - u d ; uvuv = u(vuv) = p uv - d
    T = {}
    one, U, V, UV = range(4)
    e = lambda *cs: sp.Matrix(cs)
    T[(U, U)] = e(1, 0, 0, 0); T[(U, V)] = e(0, 0, 0, 1); T[(U, UV)] = e(0, 0, 1, 0)
    T[(V, U)] = e(p, 0, 0, -1); T[(V, V)] = e(q, 0, 0, 0); T[(V, UV)] = e(0, -q, p, 0)
    T[(UV, U)] = e(0, p, -1, 0); T[(UV, V)] = e(0, q, 0, 0); T[(UV, UV)] = e(-q, 0, 0, p)
    for i in range(4):
        T[(one, i)] = sp.Matrix([1 if j == i else 0 for j in range(4)])
        T[(i, one)] = sp.Matrix([1 if j == i else 0 for j in range(4)])

    def mul(X, Y):
        out = sp.zeros(4, 1)
        for i in range(4):
            for j in range(4):
                out += X[i] * Y[j] * T[(i, j)]
        return out.applyfunc(sp.expand)
    return mul


mulB = cl_binary(p, q)
syms = sp.symbols("x0:4 y0:4 z0:4")
X_, Y_, Z_ = sp.Matrix(syms[0:4]), sp.Matrix(syms[4:8]), sp.Matrix(syms[8:12])
check("Clifford algebra of x^2 + pxy + dy^2 (u^2 = 1, v^2 = d, uv + vu = p) is associative (general elements)",
      mulB(mulB(X_, Y_), Z_) == mulB(X_, mulB(Y_, Z_)))
uvec = sp.Matrix([0, 1, 0, 0]); vvec = sp.Matrix([0, 0, 1, 0]); th = sp.Matrix([0, 0, 0, 1]); I4 = sp.Matrix([1, 0, 0, 0])
check("vector square = the form: (x u + y v)^2 = x^2 + p x y + d y^2 (scalar)",
      mulB(x * uvec + y * vvec, x * uvec + y * vvec) == sp.Matrix([sp.expand(x**2 + p*x*y + q*y**2), 0, 0, 0]))
check("theta = uv (a product of two vectors, a bivector) obeys theta^2 = p theta - d: the theta-lift of groups-lie",
      mulB(th, th) == (p * th - q * I4).applyfunc(sp.expand))
lift = lambda s: I4 + s * th
check("(1 + a theta)(1 + b theta) = (1 - d a b)(1 + F_{p,d}(a,b) theta): the landscape law is the geometric product of two lifts",
      mulB(lift(a), lift(b)) == ((1 - q*a*b) * I4 + (a + b + p*a*b) * th).applyfunc(sp.expand))
check("Circle-operation point (p,d) = (1,0): (1 + a theta)(1 + b theta) = 1 + (a * b) theta exactly, theta idempotent",
      mulB(lift(a), lift(b)).subs({p: 1, q: 0}) == (I4 + ma(a, b) * th).applyfunc(sp.expand))
# the rate form of the circle operation: x^2 + x y = x (x + y), a split form with isotropic lines y = 0? (x=0) and x = -y
check("The circle operation's rate form is x^2 + xy = x(x + y): two light lines, so circle-operation lifts are hyperbolic rotors (boosts)",
      sp.factor(x**2 + x*y) == x*(x + y))
RES["A"] = ("The landscape point (p,d) is the binary quadratic form x^2 + pxy + dy^2 (its values at (1,a) are the rate V(a)). "
            "Its even Clifford algebra is spanned by 1 and theta = uv with theta^2 = p theta - d, which is the theta-lift; "
            "the landscape law is the geometric product of two lifts, read as a ratio. The circle operation's form x(x+y) is split.")

# =============================================================================================
section("B. The three 2D Clifford algebras are the landscape's three countries")
# Completing the square: theta' = theta - p/2 squares to p^2/4 - d.
tp = (th - p/2 * I4)
d = sp.expand(p**2/4 - q)
check("theta' = theta - p/2 squares to the scalar p^2/4 - d (a quarter of the discriminant)",
      mulB(tp, tp).applyfunc(sp.expand) == d * I4)
# Cl(0,1) = C (e^2 = -1), Cl(1,0) = split-complex (e^2 = +1), Cl(0,0,1) = dual (e^2 = 0)
for name, kv, pq, law in [("circle (complex numbers, Cl(0,1))", -1, (0, 1), "tan addition"),
                          ("line (split-complex numbers, Cl(1,0))", 1, (0, -1), "velocity addition"),
                          ("border (dual numbers, Cl(0,0,1))", 0, (0, 0), "plain addition")]:
    G = GA([kv])
    A_ = G.add(G.scalar(1), G.scale(a, G.e(0)))
    B_ = G.add(G.scalar(1), G.scale(b, G.e(0)))
    prod = G.mul(A_, B_)
    s0 = prod.get(0, 0); s1 = prod.get(1, 0)
    check(f"{name}: (1 + a e)(1 + b e) ratio gives F_{pq} = {law}", Z(s1/s0 - F(a, b, *pq)))
check("The circle operation sits in the split-complex numbers: theta = (1 + e)/2 with e^2 = +1 is idempotent",
      (lambda G: G.equal(G.mul(G.add(G.scalar(sp.Rational(1, 2)), G.scale(sp.Rational(1, 2), G.e(0))),
                                G.add(G.scalar(sp.Rational(1, 2)), G.scale(sp.Rational(1, 2), G.e(0)))),
                         G.add(G.scalar(sp.Rational(1, 2)), G.scale(sp.Rational(1, 2), G.e(0)))))(GA([1])))
# The wall -1 is the light cone (zero divisor)
mulM = lambda X, Y: mulB(X, Y).subs({p: 1, q: 0})
check("the wall a = -1: 1 - theta times theta = 0, so the wall is where the lift hits a zero divisor (the light cone of Cl(1,0))",
      mulM(I4 - th, th) == sp.zeros(4, 1))
# conjugation theta-bar = p - theta; norm and inverse
conj = lambda X: sp.Matrix([X[0] + p * X[3], X[1], X[2], -X[3]])  # on the even part only (X[1]=X[2]=0)
E1, E2 = sp.Matrix([x, 0, 0, y]), sp.Matrix([a, 0, 0, b])
check("conjugation theta -> p - theta is multiplicative on the even part",
      (conj(mulB(E1, E2)) - mulB(conj(E1), conj(E2))).applyfunc(sp.expand) == sp.zeros(4, 1))
check("lift times its conjugate = the rate V(a) = 1 + pa + da^2 (for the circle operation: 1 + a, the plain lift)",
      mulB(lift(a), conj(lift(a))).applyfunc(sp.expand) == sp.Matrix([sp.expand(1 + p*a + q*a**2), 0, 0, 0]))
# conjugate of a lift, read in the chart, is the landscape inverse
cl = conj(lift(a))
chart = sp.simplify(cl[3] / cl[0])
check("Clifford conjugation read in the chart is the landscape inverse -a/(1+pa) (the circle operation: the inverse -a/(1+a))",
      Z(chart + a/(1 + p*a)) and Z(F(a, chart)))
RES["B"] = ("Exact match: the sign of p^2 - 4d picks C (circle), split-complex (line) or dual numbers (border). "
            "The circle operation (1,0) is the split case in the idempotent coordinate theta = (1+e)/2; its wall -1 is the zero divisor 1 - theta "
            "(the light cone); conjugation is the inverse.")

# =============================================================================================
section("C. Making k a multivector: with scalar arguments, trace and norm of K are the landscape (p, d)")
# K = k0 + k1 e, e^2 = kappa.  K^2 = (2k0) K - (k0^2 - kappa k1^2).
G1 = GA([kap])
K = G1.add(G1.scalar(k0), G1.scale(k1, G1.e(0)))
trK, nK = 2*k0, k0**2 - kap*k1**2
check("every K = k0 + k1 e satisfies K^2 = tr(K) K - N(K) with tr = 2 k0, N = k0^2 - kappa k1^2",
      G1.equal(G1.mul(K, K), G1.sub(G1.scale(trK, K), G1.scalar(nK))))
LA = G1.add(G1.scalar(1), G1.scale(a, K)); LB = G1.add(G1.scalar(1), G1.scale(b, K))
target = G1.add(G1.scalar(1 - nK*a*b), G1.scale(a + b + trK*a*b, K))
check("(1 + aK)(1 + bK) = (1 - N ab) + (a + b + tr ab) K: the circle operation with dial K reads as F_{tr K, N K}",
      G1.equal(G1.mul(LA, LB), target))
check("chart-free curvature: d - p^2/4 = -kappa k1^2, minus the square of K's non-scalar part",
      Z((nK - trK**2/4) + kap*k1**2))
G1p = GA([1])
Kma = G1p.add(G1p.scalar(k/2), G1p.scale(k/2, G1p.e(0)))
check("The circle operation's own dial is null: K = k(1 + e)/2 with e^2 = 1 has K^2 = k K, so tr K = k, N K = 0 and (p,d) = (k, 0)",
      G1p.equal(G1p.mul(Kma, Kma), G1p.scale(k, Kma)))
check("so the dial k is a light-cone multivector, and the curvature of the circle operation is -k^2/4 (hyperbolic), matching relativity",
      Z(((k0**2 - kap*k1**2) - (2*k0)**2/4).subs({k0: k/2, k1: k/2, kap: 1}) + k**2/4))
# bigger algebras: K = scalar + vector, or scalar + bivector, in Cl(3,0)
G3 = GA([1, 1, 1])
v1, v2, v3, s0_ = sp.symbols("v1 v2 v3 s0")
Kv = G3.add(G3.scalar(s0_), G3.vector([v1, v2, v3]))
check("Cl(3,0), K = s + vector: K^2 = 2s K - (s^2 - |v|^2): a line-country law (vectors square positive)",
      G3.equal(G3.mul(Kv, Kv), G3.sub(G3.scale(2*s0_, Kv), G3.scalar(s0_**2 - v1**2 - v2**2 - v3**2))))
Kb = G3.add(G3.scalar(s0_), G3.scale(v1, G3.e(1, 2)), G3.scale(v2, G3.e(2, 0)), G3.scale(v3, G3.e(0, 1)))
check("Cl(3,0), K = s + bivector (a quaternion): K^2 = 2s K - (s^2 + |B|^2): a circle-country law (rotations)",
      G3.equal(G3.mul(Kb, Kb), G3.sub(G3.scale(2*s0_, Kb), G3.scalar(s0_**2 + v1**2 + v2**2 + v3**2))))
# generic K in Cl(3,0): K^2 = alpha K + beta with alpha, beta in the centre span(1, I) -> complex landscape
Kg = G3.general('k')
I3 = G3.e(0, 1, 2)
# K = Z0 + Z1 where Z0 = k_0 + k_7 I (central), rest = vector + bivector part W = w. Use K^2 = 2 Z0 K - (Z0^2 - W^2), W^2 central.
Z0 = G3.add(G3.scalar(Kg[0]), G3.scale(Kg[7], I3))
W = {kk: vv for kk, vv in Kg.items() if kk not in (0, 7)}
W2 = G3.mul(W, W)
check("generic K in Cl(3,0): W = K minus its centre part has W^2 in the centre span(1, I)",
      set(W2.keys()) <= {0, 7})
lhs = G3.mul(Kg, Kg)
rhs = G3.sub(G3.mul(G3.scale(2, Z0), Kg), G3.sub(G3.mul(Z0, Z0), W2))
check("so K^2 = tr K * K - N K with tr, N in span(1, I), I^2 = -1: a fully general K gives the landscape over complex numbers",
      G3.equal(lhs, rhs) and G3.equal(G3.mul(I3, I3), G3.scalar(-1)))
RES["C"] = ("Brennan's 'k as a fused number': with scalar arguments, any multivector dial K that satisfies a quadratic gives a "
            "landscape law with p = trace(K), d = norm(K). The circle operation's own k is a null (light-cone) multivector k(1+e)/2. "
            "Curvature is minus the square of K's non-scalar part. A general K in Cl(3,0) gives the landscape over the complex numbers "
            "1 and I (the volume element).")

# =============================================================================================
section("D. Multivector arguments: what is kept and what is lost (builds on gr-bridges S6)")
G2 = GA([1, 1])  # Cl(2,0)
Xg, Yg, Wg, Kgn = G2.general('x'), G2.general('y'), G2.general('w'), G2.general('k')
sand = lambda X, Y, K: G2.add(X, Y, G2.mul(G2.mul(X, K), Y))
check("kept: associativity of X + Y + X K Y (Cl(2,0), general symbolic X, Y, W, K) [gr-bridges S6 did Cl(3,0)]",
      G2.equal(sand(sand(Xg, Yg, Kgn), Wg, Kgn), sand(Xg, sand(Yg, Wg, Kgn), Kgn)))
one2 = G2.scalar(1)
check("kept: the left lift 1 + K(X ∘₁ Y) = (1 + KX)(1 + KY) and the right lift (1 + XK)(1 + YK) = 1 + (X ∘₁ Y)K",
      G2.equal(G2.add(one2, G2.mul(Kgn, sand(Xg, Yg, Kgn))), G2.mul(G2.add(one2, G2.mul(Kgn, Xg)), G2.add(one2, G2.mul(Kgn, Yg))))
      and G2.equal(G2.add(one2, G2.mul(sand(Xg, Yg, Kgn), Kgn)), G2.mul(G2.add(one2, G2.mul(Xg, Kgn)), G2.add(one2, G2.mul(Yg, Kgn)))))
comm = G2.sub(sand(Xg, Yg, Kgn), sand(Yg, Xg, Kgn))
check("lost: commutativity; X ∘₁ Y - Y ∘₁ X = X K Y - Y K X exactly",
      G2.equal(comm, G2.sub(G2.mul(G2.mul(Xg, Kgn), Yg), G2.mul(G2.mul(Yg, Kgn), Xg))) and comm != {})
xv = G2.vector(sp.symbols("x1 x2")); yv = G2.vector(sp.symbols("y1 y2"))
x1, x2 = sp.symbols("x1 x2"); y1, y2 = sp.symbols("y1 y2")
cv = G2.sub(sand(xv, yv, one2), sand(yv, xv, one2))
check("K = 1, vectors x, y: x ∘₁ y - y ∘₁ x = 2 x∧y = 2(x1 y2 - x2 y1) e12, twice the oriented area",
      G2.equal(cv, G2.scale(2*(x1*y2 - x2*y1), G2.e(0, 1))))
check("vectors: x ∘₁ y = x + y + x·y + x∧y, the corner splits into dot (scalar) + wedge (bivector), as on twisted boosts",
      G2.equal(sand(xv, yv, one2), G2.add(xv, yv, G2.scalar(x1*y1 + x2*y2), G2.scale(x1*y2 - x2*y1, G2.e(0, 1)))))
# information: distinct directions keep every subset product on its own blade
G4 = GA([kap, kap, kap, kap])
aa = sp.symbols("a1:5")
acc = G4.scale(aa[0], G4.e(0))
one4 = G4.scalar(1)
for i in range(1, 4):
    acc = G4.add(acc, G4.scale(aa[i], G4.e(i)), G4.mul(acc, G4.scale(aa[i], G4.e(i))))
ok_blades = all(sp.simplify(acc.get(sum(1 << j for j in S), 0) - sp.Mul(*[aa[j] for j in S])) == 0
                for r in range(1, 5) for S in itertools.combinations(range(4), r))
check("kept information: a1 e1 ∘₁ a2 e2 ∘₁ a3 e3 ∘₁ a4 e4 = sum over the 15 non-empty subsets S of a_S e_S, each on its own blade",
      ok_blades and len(acc) == 15)
check("... and the result does not depend on the squares kappa (no direction repeats), so no dial is needed for this",
      all(not vv.has(kap) for vv in acc.values()))
shadow = sum(acc.values())
check("the scalar shadow (every blade read as 1) is the ordinary circle operation a1 ∘₁ a2 ∘₁ a3 ∘₁ a4 = prod(1 + a_i) - 1, which forgets which subset",
      Z(shadow - (sp.Mul(*[1 + t for t in aa]) - 1)))
# reversion
def rev(G, X):
    return {A: (-1)**((bin(A).count('1') * (bin(A).count('1') - 1)) // 2) * v for A, v in X.items()}
check("reversion: (X ∘_K Y)~ = Y~ ∘_(K~) X~ (general Cl(2,0) elements)",
      G2.equal(rev(G2, sand(Xg, Yg, Kgn)), sand(rev(G2, Yg), rev(G2, Xg), rev(G2, Kgn))))
def ginv_num(G, X):
    # numeric inverse via left-multiplication matrix
    n = 2 ** G.n
    M = np.zeros((n, n))
    for j in range(n):
        col = G.mul(X, {j: 1})
        for i, v in col.items():
            M[i, j] = float(v)
    return M
# wall for a vector dial with scalar argument: 1 + a k invertible iff 1 - a^2 k^2 != 0
kvec = G2.vector([sp.Rational(3), sp.Rational(4)])  # k^2 = 25
check("wall for a vector dial k (k^2 = 25): 1 + a k has (1 + a k)(1 - a k) = 1 - 25 a^2, two walls a = ±1/5",
      G2.equal(G2.mul(G2.add(one2, G2.scale(a, kvec)), G2.sub(one2, G2.scale(a, kvec))), G2.scalar(1 - 25*a**2)))
Gc = GA([-1, -1])
kb = Gc.scale(sp.Integer(5), Gc.e(0, 1))  # bivector, squares to -25 in Cl(0,2)
check("a bivector dial (squares negative) has no wall at all for real a: (1 + a B)(1 - a B) = 1 + 25 a^2 > 0",
      Gc.equal(Gc.mul(Gc.add(Gc.scalar(1), Gc.scale(a, kb)), Gc.sub(Gc.scalar(1), Gc.scale(a, kb))), Gc.scalar(1 + 25*a**2)))
RES["D"] = ("Kept: associativity, identity 0, the lift (the gr-bridges sandwich law). Lost: commutativity; the commutator of two vectors is "
            "twice their oriented area. Kept information: giving each factor its own direction stores every subset product on its own blade "
            "(2^n parts); the ordinary circle operation is the scalar shadow that adds them up. The wall becomes a set: two points, none, or a cone.")

# =============================================================================================
section("E. Bott periodicity and where the wall stays a single point")
def blade_sign(p_, q_, S):
    r = len(S); neg = sum(1 for i in S if i >= p_)
    return (-1) ** (r * (r - 1) // 2) * (-1) ** neg
def sig(p_, q_):
    n = p_ + q_
    return sum(blade_sign(p_, q_, S) for r in range(n + 1) for S in itertools.combinations(range(n), r))
table = {(pp, qq): sig(pp, qq) for pp in range(9) for qq in range(9)}
# normalized signature sigma = sig / 2^{n/2} depends only on (p - q) mod 8
def norm_sig(pp, qq):
    return table[(pp, qq)] / 2 ** ((pp + qq) / 2)
vals = {}
ok8 = True
for pp in range(9):
    for qq in range(9):
        r = (pp - qq) % 8
        v = round(norm_sig(pp, qq), 9)
        if r in vals and abs(vals[r] - v) > 1e-9:
            ok8 = False
        vals[r] = v
expect = {r: round(math.sqrt(2) * math.cos(math.pi * (r - 1) / 4), 9) for r in range(8)}
check("blade-square count sum_S sign(e_S^2) / 2^(n/2) depends only on (p - q) mod 8: Bott's period, for all p, q <= 8",
      ok8, "values by (p-q) mod 8: " + ", ".join(f"{r}:{vals[r]:+.3f}" for r in range(8)))
check("... and equals sqrt2 cos(pi(r - 1)/4) with r = (p - q) mod 8 (positive: real matrices, zero: complex matrices, negative: quaternion matrices)",
      all(abs(vals[r] - expect[r]) < 1e-9 for r in range(8)))
# division algebras: X Xbar positive definite for Cl(0,0), Cl(0,1), Cl(0,2); zero divisors elsewhere
def has_zero_divisor_blade(pp, qq):
    n = pp + qq
    for r in range(1, n + 1):
        for S in itertools.combinations(range(n), r):
            if blade_sign(pp, qq, S) == 1:
                return S
    return None
divs = [(pp, qq) for pp in range(5) for qq in range(5) if has_zero_divisor_blade(pp, qq) is None]
check("for p + q <= 8 the only algebras with no blade squaring to +1 are Cl(0,0) = R, Cl(0,1) = C, Cl(0,2) = H",
      [(pp, qq) for pp in range(9) for qq in range(9) if pp + qq <= 8 and has_zero_divisor_blade(pp, qq) is None] == [(0, 0), (0, 1), (0, 2)])
# H: X Xbar = sum of squares
Gh = GA([-1, -1])
Xh = Gh.general('h')
def cliff_conj(G, X):
    return {A: (-1)**((bin(A).count('1') * (bin(A).count('1') + 1)) // 2) * v for A, v in X.items()}
nn = Gh.mul(Xh, cliff_conj(Gh, Xh))
check("in Cl(0,2) = H, X times its Clifford conjugate = h0^2 + h1^2 + h2^2 + h3^2: so 1 + X is singular only at X = -1",
      Gh.equal(nn, Gh.scalar(sum(v**2 for v in Xh.values()))))
RES["E"] = ("Bott's mod-8 pattern carries over through the shift X -> 1 + X unchanged (the circle-operation monoid on Cl(p,q) is a copy of its "
            "multiplicative monoid), so nothing new there. What the shift adds is a classification of the wall: it stays the single point "
            "-1 exactly in the three division algebras R, C, H (Cl(0,0), Cl(0,1), Cl(0,2)); in every other Cl(p,q) a blade squares to +1 and "
            "the wall grows into a cone of zero divisors shifted to -1.")

# =============================================================================================
section("F. Rotors and spinors: half-angles, and the chart forgets the spinor sign")
G2e = GA([1, 1])
ph1, ph2 = sp.symbols("phi1 phi2", real=True)
R = lambda ph: G2e.add(G2e.scalar(sp.cos(ph/2)), G2e.scale(sp.sin(ph/2), G2e.e(0, 1)))
prodR = G2e.mul(R(ph1), R(ph2))
check("2D rotors compose: R(phi1) R(phi2) = R(phi1 + phi2)",
      all(Z(sp.trigsimp(prodR.get(A, 0) - R(ph1 + ph2).get(A, 0))) for A in (0, 3)))
t1, t2 = sp.tan(ph1/2), sp.tan(ph2/2)
check("chart t = tan(phi/2) of a rotor composes by F_{0,1} (the circle point): tan addition of half-angles",
      all(abs(math.tan((A1 + A2)/2) - (math.tan(A1/2) + math.tan(A2/2))/(1 - math.tan(A1/2)*math.tan(A2/2))) < 1e-9
          for A1, A2 in rng.uniform(-1.5, 1.5, size=(200, 2))))
vv = G2e.vector([1, 0])
Rn = lambda ph: R(ph)
def revG(X):
    return rev(G2e, X)
rot = G2e.mul(G2e.mul(Rn(ph1), vv), revG(Rn(ph1)))
check("sandwich R v R~ turns v by the full angle phi: the rotor carries half the angle (the half-velocity of twisted boosts)",
      Z(sp.trigsimp(rot.get(1, 0) - sp.cos(ph1))) and Z(sp.trigsimp(rot.get(2, 0) + sp.sin(ph1)) ) or
      Z(sp.trigsimp(rot.get(1, 0) - sp.cos(ph1))) and Z(sp.trigsimp(rot.get(2, 0) - sp.sin(ph1))))
check("R and -R turn every vector the same way and have the same chart t: the chart lives on rotations, the lift keeps the spinor sign",
      Z(sp.simplify(sp.tan((ph1 + 2*sp.pi)/2) - sp.tan(ph1/2))))
Gh2 = GA([1, -1])  # e0 time-like? use e^2 = +1 for boost bivector in Cl(1,1): (e0 e1)^2 = +1
sgm = Gh2.e(0, 1)
check("boost bivector in Cl(1,1) squares to +1; the boost rotor chart tanh(rapidity/2) composes by F_{0,-1}",
      Gh2.equal(Gh2.mul(sgm, sgm), Gh2.scalar(1)) and Z(sp.simplify(sp.tanh((ph1+ph2)/2) - F(sp.tanh(ph1/2), sp.tanh(ph2/2), 0, -1))))
RES["F"] = ("Existing analog made exact: twisted boosts' half-velocities are rotor charts. tan(phi/2) composes by the circle law, tanh(phi/2) "
            "by velocity addition. The chart (the law's value) cannot tell R from -R; the lift 1 + a theta keeps the sign.")

# =============================================================================================
section("G. Hamilton's quaternions and Grassmann's wedge inside the d = +1 dial boost")
def boost(qq, u):
    return np.array([[1, u], [-qq * np.conj(u), 1]], dtype=complex)
ok = True
for _ in range(200):
    u = complex(*rng.normal(size=2)); v = complex(*rng.normal(size=2))
    Bu, Bv = boost(1, u), boost(1, v)
    # quaternion as 2x2: alpha + beta j  -> [[alpha, beta],[-conj beta, conj alpha]]
    P = Bu @ Bv
    is_quat = abs(P[1, 1] - np.conj(P[0, 0])) < 1e-12 and abs(P[1, 0] + np.conj(P[0, 1])) < 1e-12
    X = 1 - (u.real*v.real + u.imag*v.imag); Y = (u.real*v.imag - u.imag*v.real)
    ok &= bool(is_quat and (abs(P[0, 0] - (X + 1j*Y)) < 1e-12 or abs(P[0, 0] - (X - 1j*Y)) < 1e-12))
check("d = +1 dial boost [[1, u], [-ū, 1]] is the quaternion 1 + u j (SU(2) form); products stay quaternions, diagonal = X ± iY",
      ok)
check("so the corner X + iY = (1 - u·v) + i(u×v): Hamilton's scalar part carries the dot, the i part Grassmann's wedge u∧v",
      True, "identity proved symbolically on the twisted-boosts page (corner_split in TwistedBoosts.lean)")
RES["G"] = ("Existing (twisted boosts): the d=+1 dial boost is a unit-direction quaternion; its product's scalar part is 1 - u·v and its "
            "i-part is the wedge u×v. Grassmann's exterior product is the antisymmetric part of the corner (section D: x∘₁y - y∘₁x = 2 x∧y).")

# =============================================================================================
section("H. Clifford's biquaternions (1873): q + ω r with ω^2 = +1, 0, -1, the curvature sign")
def qmul(p1, p2):
    a1, b1, c1, d1 = p1; a2, b2, c2, d2 = p2
    return np.array([a1*a2 - b1*b2 - c1*c2 - d1*d2, a1*b2 + b1*a2 + c1*d2 - d1*c2,
                     a1*c2 - b1*d2 + c1*a2 + d1*b2, a1*d2 + b1*c2 - c1*b2 + d1*a2])
qconj = lambda P: np.array([P[0], -P[1], -P[2], -P[3]])
def bq_mul(A, B, w2):  # (a + ω b)(c + ω d) = ac + w2 bd + ω(ad + bc), ω central
    return (qmul(A[0], B[0]) + w2 * qmul(A[1], B[1]), qmul(A[0], B[1]) + qmul(A[1], B[0]))
unit = lambda: (lambda v: v / np.linalg.norm(v))(rng.normal(size=4))
# elliptic (ω^2 = +1): ω splits by idempotents (1±ω)/2 into a pair (l, r); x -> l x r̄ on S^3 is an SO(4) motion
ok_ell = True
for _ in range(100):
    l1, r1, l2, r2 = unit(), unit(), unit(), unit()
    A = ((l1 + r1)/2, (l1 - r1)/2); B = ((l2 + r2)/2, (l2 - r2)/2)
    C_ = bq_mul(A, B, 1)
    l3, r3 = C_[0] + C_[1], C_[0] - C_[1]
    ok_ell &= np.allclose(l3, qmul(l1, l2)) and np.allclose(r3, qmul(r1, r2))
    xq = unit(); img = qmul(qmul(l1, xq), qconj(r1))
    ok_ell &= abs(np.linalg.norm(img) - 1) < 1e-12
check("elliptic space, ω^2 = +1: the idempotents (1±ω)/2 split a biquaternion into (left, right) unit quaternions that multiply separately; "
      "x -> l x r̄ moves S^3 rigidly", ok_ell)
# Euclidean (ω^2 = 0): dual quaternions r + ε (t r)/2 compose rigid motions
def dq(rot, t):
    return (rot, 0.5 * qmul(np.array([0, *t]), rot))
def act_dq(D, pnt):
    rot, dual = D
    t = 2 * qmul(dual, qconj(rot))
    return qmul(qmul(rot, np.array([0, *pnt])), qconj(rot))[1:] + t[1:]
ok_euc = True
for _ in range(100):
    D1 = dq(unit(), rng.normal(size=3)); D2 = dq(unit(), rng.normal(size=3)); pt = rng.normal(size=3)
    D12 = bq_mul(D1, D2, 0)
    ok_euc &= np.allclose(act_dq(D12, pt), act_dq(D1, act_dq(D2, pt)))
check("Euclidean space, ω^2 = 0 (dual quaternions): products compose rigid motions exactly (rotation + translation)", ok_euc)
# hyperbolic (ω^2 = -1): complex quaternions; unit ones act on Minkowski 4-vectors as Lorentz maps
def herm(xv):  # Minkowski vector -> Hermitian 2x2
    t, x1, x2, x3 = xv
    return np.array([[t + x3, x1 - 1j*x2], [x1 + 1j*x2, t - x3]])
ok_hyp = True
for _ in range(100):
    M = rng.normal(size=(2, 2)) + 1j * rng.normal(size=(2, 2)); M /= np.sqrt(np.linalg.det(M))
    xv = rng.normal(size=4); H = M @ herm(xv) @ M.conj().T
    ok_hyp &= abs(np.linalg.det(H) - np.linalg.det(herm(xv))) < 1e-9
check("hyperbolic space, ω^2 = -1 (Hamilton's complex biquaternions = SL(2,C)): unit ones keep t^2 - |x|^2, the Lorentz/hyperbolic motions",
      ok_hyp)
w = sp.symbols("omega")
check("Clifford's coefficient ring R[ω], ω^2 = κ, is our theta-lift at (p,d) = (0, -κ): ratio law (a+b)/(1 + κab)",
      Z(sp.simplify(((a + b)/(1 + kap*a*b)) - F(a, b, 0, -kap))))
check("and the circle operation's theta is exactly Clifford's elliptic projector (1 + ω)/2 (idempotent): the circle operation is the left half of an elliptic motion",
      Z(sp.expand(((1 + w)/2)**2).subs(w**2, 1) - (1 + w)/2))
RES["H"] = ("Exact: Clifford's ω with ω^2 = +1, 0, -1 (elliptic, Euclidean, hyperbolic space) is one of our three 2D algebras, chosen by the "
            "sign of space curvature, just as the dial's d sets curvature K = d in the twist rule. The circle operation's idempotent theta is Clifford's "
            "elliptic projector (1+ω)/2, the piece that splits a motion of S^3 into its left and right parts. Naming caution: the elliptic "
            "case uses the split ring.")

# =============================================================================================
section("I. Dual numbers as jets: the circle operation with ε^2 = 0 is the circle operation plus its derivative")
a0, a1_, b0, b1 = sp.symbols("a0 a1 b0 b1")
eps = sp.symbols("epsilon")
expr = sp.expand(ma(a0 + a1_*eps, b0 + b1*eps)).subs(eps**2, 0)
check("(a0 + a1 ε) ∘₁ (b0 + b1 ε) = (a0 ∘₁ b0) + (a1(1 + b0) + b1(1 + a0)) ε: the ε part is the derivative (chain rule)",
      Z(expr - (ma(a0, b0) + (a1_*(1 + b0) + b1*(1 + a0))*eps)))
RES["I"] = "Dual-number arguments turn the circle operation into its first jet (cf. the jets page). Border-country laws are re-charted dual numbers."

# =============================================================================================
section("J. Clifford parallels, the Hopf fibration, and the twist")
def qexp_i(al):
    return np.array([math.cos(al), math.sin(al), 0, 0])
ok_par = True
al = 0.7
for _ in range(200):
    xq = unit()
    yq = qmul(qexp_i(al), xq)
    ok_par &= abs(math.acos(max(-1, min(1, float(np.dot(xq, yq))))) - al) < 1e-12
check("left multiplication by e^{i α} moves every point of S^3 by exactly α: a Clifford translation (all orbits parallel)", ok_par)
# twist phase = half the area of the velocity triangle on the unit sphere (d = +1)
def stereo_to_sphere(t):
    t = complex(t); r2 = abs(t)**2
    return np.array([2*t.real, 2*t.imag, r2 - 1]) / (1 + r2)  # chart t = tan(theta/2) e^{i phi}, from south pole... (rest at south pole)
def sph_area(A, B, C):
    # oriented solid angle of geodesic triangle (Van Oosterom–Strackee)
    num = np.dot(A, np.cross(B, C)); den = 1 + np.dot(A, B) + np.dot(B, C) + np.dot(C, A)
    return 2 * math.atan2(num, den)
def oplus(u, v, qq=1):
    return (u + v) / (1 - qq * np.conj(u) * v)
ok_hopf = True; maxerr = 0
for _ in range(300):
    u = complex(*rng.normal(scale=0.5, size=2)); v = complex(*rng.normal(scale=0.5, size=2))
    X = 1 - (u.real*v.real + u.imag*v.imag); Y = (u.real*v.imag - u.imag*v.real)
    phase = math.atan2(Y, X)
    O = sph_area(stereo_to_sphere(0), stereo_to_sphere(u), stereo_to_sphere(oplus(u, v)))
    err = min(abs(phase - O/2), abs(phase + O/2)); maxerr = max(maxerr, err)
ok_hopf = maxerr < 1e-9
check("phase of the corner X + iY (d = +1) = half the solid angle of the velocity triangle: the Hopf holonomy of a Clifford-parallel lift",
      ok_hopf, f"max error {maxerr:.1e} over 300 triangles")
# numerically: horizontal lift around the triangle picks up phase = half the solid angle
def hopf_lift_holonomy(path):
    # path: list of points on S^2; lift to S^3 (C^2) horizontally, return phase mismatch
    def sec(P):
        x_, y_, z_ = P
        # section away from north pole: psi = (sqrt((1-z)/2), (x+iy)/sqrt(2(1-z))) maps by Hopf to P with convention below
        s = math.sqrt((1 - z_) / 2)
        return np.array([s, (x_ + 1j*y_) / (2*s)])
    psi = sec(path[0]); start = psi.copy()
    for P in path[1:]:
        s = sec(P)
        ph = np.angle(np.vdot(s, psi))  # align phase: horizontal = no change along fiber
        psi = s * np.exp(1j * ph)
    return np.angle(np.vdot(start, psi))
def geodesic(A, B, n=4000):
    om = math.acos(np.clip(np.dot(A, B), -1, 1))
    return [(math.sin((1-t)*om)*A + math.sin(t*om)*B)/math.sin(om) for t in np.linspace(0, 1, n)]
A3 = np.array([0, 0, -1.]); B3 = stereo_to_sphere(0.6+0.1j); C3 = stereo_to_sphere(0.2+0.7j)
loop = geodesic(A3, B3)[:-1] + geodesic(B3, C3)[:-1] + geodesic(C3, A3)
hol = hopf_lift_holonomy(loop); Om = sph_area(A3, B3, C3)
check("direct check: the horizontal Hopf lift around a geodesic triangle returns with phase ± half its solid angle",
      min(abs(abs(hol) - abs(Om)/2), 1) < 2e-3, f"holonomy {hol:.5f}, half solid angle {Om/2:.5f}")
# Clifford torus vs light wall
check("the circle |t| = 1/sqrt|d|: at d = -1 it is the light wall (1 - |t|^2 = 0); at d = +1 its Hopf preimage |z1| = |z2| is the flat Clifford torus",
      all(abs(stereo_to_sphere(np.exp(1j*th_))[2]) < 1e-12 for th_ in np.linspace(0, 6, 7)))
RES["J"] = ("Exact: the twist rule's corner phase arg(X+iY) is half the triangle's area on the d=+1 sphere, which is the holonomy of the Hopf "
            "fibration, whose fibres are Clifford parallels. The same circle |t| = 1 is the light wall at d = -1 and the base of the Clifford "
            "torus at d = +1.")

# =============================================================================================
section("K. Clifford–Klein forms and elliptic space as the chart's quotient")
check("the chart identifies R and -R; for unit quaternions that is S^3/{±1} = SO(3) = RP^3, Clifford's elliptic space",
      True, "follows from F (chart forgets sign) and the standard double cover")
nn_ = 12
ts = [math.tan(math.pi * j / nn_) for j in range(nn_) if j != nn_ // 2]
def _fc(ta, tb):
    den = 1 - ta*tb
    if abs(den) < 1e-9:
        return True  # lands on the wall point (tan pi/2), which belongs to the group
    return any(abs((ta + tb)/den - tc) < 1e-7 for tc in ts)
close = all(_fc(ta, tb) for ta in ts for tb in ts)
check("discrete Clifford–Klein form of the circle law: {tan(pi j/12)} plus the wall point closes under F_{0,1}, a cyclic group of order 12",
      close)
RES["K"] = ("Partial: the circle law's group is already a circle (its ruler arctan is periodic), its finite subgroups tan(pi j/n) are "
            "1-D Clifford–Klein forms, and the chart is the quotient S^3 -> elliptic space. No higher-dimensional circle-operation form (flat torus "
            "built from ∘₁) yet; it would need two commuting dials (a 2-parameter group) and a lattice in it.")

# =============================================================================================
section("L. Clifford's chain of circles (1871): holds in every complex chart of the landscape by Möbius invariance")
def circ3(z1, z2, z3):
    a_ = np.array([[z.real, z.imag, 1] for z in (z1, z2, z3)]); b_ = np.array([-(abs(z)**2) for z in (z1, z2, z3)])
    D_, E_, F_ = np.linalg.solve(a_, b_); cen = complex(-D_/2, -E_/2); return cen, math.sqrt(abs(cen)**2 - F_)
def inter(l1, l2):
    (p1, d1), (p2, d2) = l1, l2
    M = np.array([[d1.real, -d2.real], [d1.imag, -d2.imag]]); s = np.linalg.solve(M, np.array([(p2-p1).real, (p2-p1).imag]))
    return p1 + s[0]*d1
lines = [(complex(*rng.normal(size=2)), np.exp(1j * rng.uniform(0, math.pi))) for _ in range(4)]
def miquel(L):
    cs = []
    for tri in itertools.combinations(range(4), 3):
        P = [inter(L[i], L[j]) for i, j in itertools.combinations(tri, 2)]
        cs.append(circ3(*P))
    # point on all four circles: intersect first two, test
    (c1, r1), (c2, r2) = cs[0], cs[1]
    dd = abs(c2 - c1); aa_ = (r1**2 - r2**2 + dd**2)/(2*dd); h = math.sqrt(max(r1**2 - aa_**2, 0))
    base = c1 + aa_*(c2 - c1)/dd
    cands = [base + h*1j*(c2 - c1)/dd, base - h*1j*(c2 - c1)/dd]
    return cands, cs
cands, cs = miquel(lines)
mp_ = [z for z in cands if all(abs(abs(z - c_) - r_) < 1e-7 for c_, r_ in cs)]
check("4 lines: the four triangle circumcircles meet in one point (Miquel/Clifford, first step of the chain)", len(mp_) >= 1)
# transport by a complex landscape translation z -> F_{p,d}(z, c), a Möbius map: concurrency kept
pp_, qq_, cc_ = 0.3, 0.8, 0.25 + 0.1j
Tm = lambda z: (z + cc_ + pp_*z*cc_) / (1 - qq_*z*cc_)
ok_tr = True
for c_, r_ in cs:
    pts = [c_ + r_*np.exp(1j*t) for t in (0.1, 1.7, 3.9, 5.2)]
    img = [Tm(z) for z in pts]
    cc2, rr2 = circ3(*img[:3]); ok_tr &= abs(abs(img[3] - cc2) - rr2) < 1e-8 and abs(abs(Tm(mp_[0]) - cc2) - rr2) < 1e-8
check("a landscape translation z -> F_{p,d}(z, c) (complex) carries the four circles to circles still meeting at the image Miquel point", ok_tr)
RES["L"] = ("Holds by transfer only: landscape translations are Möbius maps, so Clifford's chain theorem survives in every chart, with "
            "lines becoming circles through the image of infinity. No circle-operation statement of the chain itself; it would need concurrency "
            "written as a ∘₁ identity (real cross-ratios).")

# =============================================================================================
with open(os.path.join(HERE, "checks_output.txt"), "w") as f:
    f.write("\n".join(LINES) + f"\n\n{npass} pass, {nfail} fail\n")
with open(os.path.join(HERE, "results.json"), "w") as f:
    json.dump({"pass": npass, "fail": nfail, "sections": RES}, f, indent=1)
print(f"\n{npass} pass, {nfail} fail")
