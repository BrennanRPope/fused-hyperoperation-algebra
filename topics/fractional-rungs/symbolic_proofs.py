"""Exact (symbolic) proofs for the fractional and negative rungs.

Every check here is an exact algebraic identity verified by sympy (expand/simplify to 0),
so each PASS is a proof of that identity, not a numeric sample.

Notation
  E(x) = e^x - 1   (the circle operation's exponential: E(a+b) = E(a) (+) E(b))
  L(x) = log(1+x)  (its inverse, the circle operation's logarithm)
  dial:  a (+)_k b = a + b + k ab
  s-ladder: B_s(a,b) = E^s( E^-s(a) + E^-s(b) ),  E^s = s-th iterate of E
            B_0 = +, B_1 = the circle operation, B_-1 = log(e^a + e^b - 1)
"""
import sympy as sp

a, b, c, d, k, s, t, x, y, h = sp.symbols('a b c d k s t x y h')
results = []

def check(name, expr):
    ok = sp.simplify(sp.expand(expr)) == 0
    results.append((name, ok))
    print(('PASS ' if ok else 'FAIL ') + name)
    return ok

dial = lambda u, v, kk=k: u + v + kk*u*v
madd = lambda u, v: u + v + u*v

# ---------- 1. The dial ----------
check("dial commutative", dial(a, b) - dial(b, a))
check("dial associative", dial(dial(a, b), c) - dial(a, dial(b, c)))
check("dial identity 0", dial(a, 0) - a)
check("dial inverse -a/(1+ka)", sp.together(dial(a, -a/(1 + k*a))))
check("dial = the circle operation rescaled: k (a (+)_k b) = (ka)(+)(kb)", k*dial(a, b) - madd(k*a, k*b))
check("dial lift: 1 + k(a (+)_k b) = (1+ka)(1+kb)", 1 + k*dial(a, b) - (1 + k*a)*(1 + k*b))
check("mirror: -((-a)(+)_k(-b)) = a (+)_{-k} b", -dial(-a, -b) - dial(a, b, -k))
check("k=-1 is probabilistic OR: 1-(1-a)(1-b)", dial(a, b, -1) - (1 - (1 - a)*(1 - b)))
three = dial(dial(a, b), c)
check("three-way dial: layers weighted 1, k, k^2 (dimension d has weight k^(d-1))",
      three - (a + b + c + k*(a*b + b*c + c*a) + k**2*a*b*c))
four = dial(three, d)
check("four-way dial: lift = product of (1+k a_i)", 1 + k*four - (1+k*a)*(1+k*b)*(1+k*c)*(1+k*d))

# ---------- 2. Iterating a tangent-to-identity map ----------
# Formal iterate of f(x) = x + c2 x^2 + c3 x^3 + c4 x^4 + ...: coefficients are polynomials in s.
# We find them from integer iterates and then PROVE the group law P_s o P_t = P_{s+t}.
N = 7  # order kept
def compose(f, g, order=N):
    """f(g(x)) truncated, f,g as polynomials in x."""
    r = sp.expand(f.subs(x, g))
    return sum(r.coeff(x, i)*x**i for i in range(1, order + 1))

Eser = sum(x**i/sp.factorial(i) for i in range(1, N + 1))       # e^x - 1
iters = [x]
for n in range(1, N + 2):
    iters.append(compose(Eser, iters[-1]))
coeffs = []
for i in range(1, N + 1):
    pts = [(n, iters[n].coeff(x, i)) for n in range(0, N + 1)]
    coeffs.append(sp.expand(sp.interpolate(pts, s)))
Es = sum(coeffs[i - 1]*x**i for i in range(1, N + 1))
print("E^s(x) =", sp.collect(sum(sp.factor(coeffs[i-1])*x**i for i in range(1, 5)), x), "+ ...")
check("E^s at s=1 is e^x - 1 (to order 7)", Es.subs(s, 1) - Eser)
check("E^s at s=0 is x", Es.subs(s, 0) - x)
check("group law E^s o E^t = E^(s+t) (to order 7, all s,t)",
      compose(Es, Es.subs(s, t)) - Es.subs(s, s + t))
check("quadratic coefficient of E^s is s/2", coeffs[1] - s/2)
check("cubic coefficient of E^s is s/6 + s(s-1)/4", coeffs[2] - (s/6 + s*(s - 1)/4))
Eminus = Es.subs(s, -s)
check("E^-s o E^s = identity (to order 7)", compose(Eminus, Es) - x)
check("E^-1 is log(1+x) (to order 7)",
      Es.subs(s, -1) - sp.series(sp.log(1 + x), x, 0, N + 1).removeO())

# ---------- 3. The s-ladder near the identity ----------
# B_s(a,b) = E^s(E^-s a + E^-s b), expanded in a, b to total degree 4.
def trunc(expr, deg):
    p = sp.Poly(sp.expand(expr), a, b)
    return sum(cf*a**m[0]*b**m[1] for m, cf in p.terms() if sum(m) <= deg)
Es4 = sum(coeffs[i - 1]*x**i for i in range(1, 5)); Em4 = Es4.subs(s, -s)
u = Em4.subs(x, a); v = Em4.subs(x, b)
w = trunc(u + v, 4)
Bs = sum(coeffs[i - 1]*trunc(w**i, 4) for i in range(1, 5))
B4 = trunc(Bs, 4)
B4c = sp.collect(sp.expand(B4), [a, b])
print("B_s(a,b) to degree 3 =", sp.factor(trunc(Bs, 3) - a - b))
check("B_s tangent to addition, corner weight exactly s: B_s = a + b + s ab + O(3)",
      trunc(Bs, 2) - (a + b + s*a*b))
cubic = sp.factor(sp.Poly(trunc(Bs, 3) - trunc(Bs, 2), a, b).coeff_monomial(a**2*b))
print("cubic coefficient of a^2 b in B_s:", cubic)
check("cubic bend of B_s is s(s-1)/4 * (a^2 b + a b^2): zero only at s = 0 and s = 1",
      trunc(Bs, 3) - trunc(Bs, 2) - (s*(s - 1)/4)*(a**2*b + a*b**2))
check("B_1 (the circle operation) has no terms beyond ab (to degree 4)", B4.subs(s, 1) - madd(a, b))
check("B_0 is plain addition (to degree 4)", B4.subs(s, 0) - (a + b))
check("B_s commutative (to degree 4)", B4 - B4.subs({a: b, b: a}, simultaneous=True))
# associativity of the truncated law, degree 4
B4f = lambda p, q: B4.subs({a: p, b: q}, simultaneous=True)
lhs = sp.Poly(sp.expand(B4f(B4f(a, b), c)), a, b, c)
rhs = sp.Poly(sp.expand(B4f(a, B4f(b, c))), a, b, c)
diff = sum(cf*sp.prod([g**e for g, e in zip((a, b, c), m)]) for m, cf in (lhs - rhs).terms() if sum(m) <= 4)
check("B_s associative (to degree 4, all s)", diff)

# ---------- 4. The negative rung B_-1 ----------
Bm1 = lambda p, q: sp.log(sp.exp(p) + sp.exp(q) - 1)
check("B_-1(a,b) = L(E(a) + E(b))", sp.simplify(sp.exp(Bm1(a, b)) - sp.exp(sp.log(1 + (sp.exp(a)-1) + (sp.exp(b)-1)))))
check("B_-1 identity 0", sp.simplify(sp.exp(Bm1(a, 0)) - sp.exp(a)))
check("B_-1 associative: both sides = log(e^a+e^b+e^c-2)",
      sp.simplify(sp.exp(Bm1(Bm1(a, b), c)) - sp.exp(Bm1(a, Bm1(b, c)))))
check("B_-1 zoomed in is addition: h B_-1(a/h, b/h) -> a + b",
      sp.limit(h*Bm1(a/h, b/h), h, sp.oo) - (a + b))
check("B_-1 corner weight is -1 (series)", sp.series(Bm1(x, y).subs(y, x), x, 0, 3).removeO() - (2*x - x**2))

# ---------- 5. Fusing does nothing below addition ----------
# Fusing = conjugating by x -> 1 + x:  fused op(a,b) = op(1+a, 1+b) - 1.
fuse = lambda op: (lambda p, q: op(1 + p, 1 + q) - 1)
lse = lambda p, q: sp.log(sp.exp(p) + sp.exp(q))
succ = lambda p, q: q + 1
check("fused successor = successor", fuse(succ)(a, b) - succ(a, b))
check("fused log-sum-exp = log-sum-exp", sp.simplify(sp.exp(fuse(lse)(a, b)) - sp.exp(lse(a, b))))
check("fused addition = a + b + 1 (the first rung the shift touches)", fuse(lambda p, q: p + q)(a, b) - (a + b + 1))
check("fused multiplication = the circle operation", fuse(lambda p, q: p*q)(a, b) - madd(a, b))
# general reason: op(a+c, b+c) = op(a,b) + c  =>  fused op = op. Shown for max via cases in Lean.

# ---------- 6. Negative and fractional right inputs of the fused exponent ----------
fexp = lambda p, q: (1 + p)**(1 + q) - 1
check("a (*)3 (-1) = 0 (identity)", sp.simplify(fexp(a, -1)))
check("a (*)3 (-2) = circle inverse of a", sp.simplify(madd(a, fexp(a, -2))))
r = fexp(a, sp.Rational(-1, 2))
check("a (*)3 (-1/2) is the circle-operation square root: r (+) r = a", sp.simplify(madd(r, r) - a))
r3 = fexp(a, sp.Rational(-2, 3))
check("a (*)3 (-2/3) is the circle-operation cube root", sp.simplify(madd(madd(r3, r3), r3) - a))
ap, bp, cp = sp.symbols('ap bp cp', positive=True)
check("dimension law: (a (*)3 b) (*)3 c = a (*)3 (b (+) c)  [tower law; dimensions multiply]",
      sp.powsimp(sp.expand(fexp(fexp(ap, bp), cp) - fexp(ap, madd(bp, cp))), force=True))

# ---------- 7. Hirota: the KdV corner weight is a dial position ----------
k1, k2, X, T = sp.symbols('k1 k2 X T', positive=True)
A = ((k1 - k2)/(k1 + k2))**2
e1 = sp.exp(k1*X - k1**3*T); e2 = sp.exp(k2*X - k2**3*T)
f = 1 + e1 + e2 + A*e1*e2
check("KdV tau function f - 1 = e1 (+)_A e2  (dial with k = A)", (f - 1) - dial(e1, e2, A))
fx = lambda g, n=1: sp.diff(g, X, n)
bil = (f*sp.diff(f, X, T) - fx(f)*sp.diff(f, T)) + (f*fx(f, 4) - 4*fx(f)*fx(f, 3) + 3*fx(f, 2)**2)
check("Hirota bilinear KdV (DxDt + Dx^4) f.f = 0 with that corner", sp.simplify(sp.expand(bil)))
for kk_name, kk in [("A = 0 (corner empty)", 0), ("A = 1 (the full circle operation)", 1)]:
    g = 1 + e1 + e2 + kk*e1*e2
    bad = (g*sp.diff(g, X, T) - fx(g)*sp.diff(g, T)) + (g*fx(g, 4) - 4*fx(g)*fx(g, 3) + 3*fx(g, 2)**2)
    nz = sp.simplify(sp.expand(bad)) != 0
    results.append(("KdV fails with " + kk_name + " when k1 != k2", nz))
    print(('PASS ' if nz else 'FAIL ') + "KdV fails with " + kk_name + " when k1 != k2")
check("A lies in [0,1): 1 - A = 4 k1 k2/(k1+k2)^2", sp.simplify((1 - A) - 4*k1*k2/(k1 + k2)**2))

print()
print(f"{sum(ok for _, ok in results)} of {len(results)} exact checks pass")
