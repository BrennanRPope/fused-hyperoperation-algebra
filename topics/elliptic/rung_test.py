"""Verdict test: does an elliptic curve's law sit at a fractional rung of the fused tower?

Two formal group laws (written in the same coordinate) are equal iff their logarithms are equal.
  rung s:   B_s(a,b) = E^s(E^-s a + E^-s b),   log = E^-s,  E(x) = e^x - 1 (regular iterate)
  elliptic: y^2 + a1 xy + a3 y = x^3 + a2 x^2 + a4 x + a6, parameter z = -x/y,
            log = integral of the invariant differential  w(z) dz = (1 + d1 z + d2 z^2 + ...) dz,
            d_n is a polynomial in a1..a6 of weight n.  There is no a5, so d5 is forced by a1..a4.
Fit a1..a4 to the rung's log through layer 5, then see what layer 6 says (and a6 / layer 8).
"""
import sympy as sp
N = 9  # series order
x, s, z = sp.symbols('x s z')
a1, a2, a3, a4, a6 = sp.symbols('a1 a2 a3 a4 a6')

def trunc(e, var, n=N):
    return sp.series(sp.expand(e), var, 0, n).removeO()

# ---- E^m for whole m by composition, then interpolate polynomial in s (exact)
E = trunc(sp.exp(x) - 1, x)
def compose(f, g):  # f(g(x))
    return sp.expand(trunc(f.subs(x, g), x))
iters = {0: x}
cur = x
for m in range(1, N):
    cur = compose(E, cur); iters[m] = cur
Es = x
for n in range(2, N):
    pts = [(m, iters[m].coeff(x, n)) for m in range(0, N)]
    poly = sp.interpolate(pts, s)
    Es += sp.factor(poly) * x**n
Es = sp.expand(Es)
assert sp.expand(Es.coeff(x, 2) - s/2) == 0
logB = sp.expand(Es.subs(s, -s))            # log of rung s = E^-s

# ---- invariant differential of the Weierstrass curve in z
w = z**3
for _ in range(N + 2):
    w = trunc(z**3 + a1*z*w + a2*z**2*w + a3*w**2 + a4*z*w**2 + a6*w**3, z, N + 3)
X = sp.expand(z / w)           # Laurent series
Y = sp.expand(-1 / w)
# omega = dx / (2y + a1 x + a3); compute with series in z
num = sp.diff(X, z)
den = 2*Y + a1*X + a3
omega = trunc(sp.series(num/den, z, 0, N), z)
d = [sp.expand(omega.coeff(z, n)) for n in range(N)]
assert d[0] == 1
logE = sp.expand(sum(d[n] * z**(n+1) / (n+1) for n in range(N-1)))

out = []
def p(*a):
    print(*a); out.append(' '.join(str(t) for t in a))

for n in range(1, 8):
    p(f"d{n} =", sp.factor(d[n]))

# ---- fit layer by layer
rung = [sp.expand(logB.coeff(x, n+1)) for n in range(N-1)]  # rung log coeff of z^(n+1)
ell = [sp.expand(logE.coeff(z, n+1)) for n in range(N-1)]
sol = {}
for n, var in [(1, a1), (2, a2), (3, a3), (4, a4)]:
    eq = sp.expand((ell[n] - rung[n]).subs(sol))
    v = sp.solve(eq, var)[0]
    sol[var] = sp.factor(v)
    p(f"layer {n+1}: fit {var} =", sol[var])
gap6 = sp.factor(sp.expand((ell[5] - rung[5]).subs(sol)))
p("layer 6 mismatch (elliptic minus rung), no free coefficient left:", gap6)
p("  roots in s:", sp.solve(sp.numer(sp.together(gap6)), s))
eq7 = sp.expand((ell[6] - rung[6]).subs(sol))
sol[a6] = sp.factor(sp.solve(eq7, a6)[0])
p("layer 7: fit a6 =", sol[a6])
gap8 = sp.factor(sp.expand((ell[7] - rung[7]).subs(sol)))
p("layer 8 mismatch:", gap8)
p("  roots in s:", sp.solve(sp.numer(sp.together(gap8)), s))

# whole rungs and samples
for sv in [sp.Rational(1, 2), sp.Rational(1, 3), 2, -1, 3]:
    p(f"s = {sv}: layer-6 gap {sp.nsimplify(gap6.subs(s, sv))},  layer-8 gap {sp.nsimplify(gap8.subs(s, sv))}")

# discriminant of the fitted curve as a function of s (is it a real elliptic curve or singular?)
b2 = a1**2 + 4*a2; b4 = 2*a4 + a1*a3; b6 = a3**2 + 4*a6
b8 = a1**2*a6 + 4*a2*a6 - a1*a3*a4 + a2*a3**2 - a4**2
Delta = -b2**2*b8 - 8*b4**3 - 27*b6**2 + 9*b2*b4*b6
Dfit = sp.factor(Delta.subs(sol))
p("discriminant of the curve fitted through layer 7:", Dfit)
import json
json.dump({"fit": {str(k): str(v) for k, v in sol.items()}, "gap6": str(gap6), "gap8": str(gap8),
           "disc": str(Dfit), "d": [str(t) for t in d[:8]]}, open("rung_test.json", "w"), indent=1)
open("rung_test_output.txt", "w").write("\n".join(out) + "\n")
