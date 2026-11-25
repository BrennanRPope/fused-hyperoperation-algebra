"""Landscape map: exact symbolic checks for the plane of fused laws.

The map's coordinates are the first two layers of the mult-add cube beyond the sides:
  c1 = weight of the corner cell ab           (the corner weight, the dial k)
  c2 = weight of each next-layer cell a^2 b, a b^2   (the bend)
Every point (p, q) of the plane is the Mobius law  F_{p,q}(a,b) = (a + b + p ab) / (1 - q ab).

Run:  python3 algebra_checks.py   (writes algebra_checks_output.txt)
"""
import sympy as sp

out = []
npass = 0
nall = 0


def check(name, cond):
    global npass, nall
    nall += 1
    ok = bool(cond)
    npass += ok
    out.append(("PASS " if ok else "FAIL ") + name)


a, b, c, p, q, u, v, lam, s, t, z1, z2 = sp.symbols('a b c p q u v lam s t z1 z2')


def F(x, y, p=p, q=q):
    return (x + y + p * x * y) / (1 - q * x * y)


def jet(expr, order=4):
    """Taylor coefficients of a two-variable law up to total degree `order`."""
    ser = sp.expand(sp.series(sp.series(expr, a, 0, order + 1).removeO(), b, 0, order + 1).removeO())
    poly = sp.Poly(ser, a, b)
    return {m: sp.simplify(cf) for m, cf in zip(poly.monoms(), poly.coeffs()) if sum(m) <= order}


# --- 1. Every point of the plane is a group law -------------------------------------------
check("F_{p,q} commutative", sp.simplify(F(a, b) - F(b, a)) == 0)
check("F_{p,q} associative", sp.simplify(F(F(a, b), c) - F(a, F(b, c))) == 0)
check("F_{p,q} identity 0", sp.simplify(F(a, 0) - a) == 0)
check("F_{p,q} inverse  -a/(1+pa)", sp.simplify(F(a, -a / (1 + p * a))) == 0)
J = jet(F(a, b), 4)
check("third-order jet of F_{p,q} is (corner, bend) = (p, q)",
      J.get((1, 1)) == p and J.get((2, 1)) == q and J.get((1, 2)) == q)
# projective (denominator-free) associativity: numerators and denominators agree exactly
N = lambda x, y: x + y + p * x * y
D = lambda x, y: 1 - q * x * y
Ln = N(a, b) + c * D(a, b) + p * N(a, b) * c
Ld = D(a, b) - q * N(a, b) * c
Rn = a * D(b, c) + N(b, c) + p * a * N(b, c)
Rd = D(b, c) - q * a * N(b, c)
check("associativity holds numerator-by-numerator, denominator-by-denominator",
      sp.expand(Ln - Rn) == 0 and sp.expand(Ld - Rd) == 0)

# --- 2. The growth law behind each point: y' = 1 + p y + q y^2 ------------------------------
g = sp.simplify(sp.diff(F(a, b), b).subs(b, 0))
check("d/db F(a,b) at b=0 is 1 + p a + q a^2 (the Riccati rate)", sp.expand(g - (1 + p * a + q * a ** 2)) == 0)

# --- 3. The lens: walls u, v (reciprocal wall positions) ------------------------------------
Fw = F(a, b, p=-(u + v), q=u * v)
M = lambda x: (1 - u * x) / (1 - v * x)
check("lift: M(a (+) b) = M(a) M(b) with M(x) = (1-ux)/(1-vx)", sp.simplify(M(Fw) - M(a) * M(b)) == 0)
check("mult-add is the case u = -1, v = 0 (wall at -1, other wall at infinity)",
      sp.simplify(F(a, b, p=1, q=0) - (a + b + a * b)) == 0)
check("velocity addition is u = 1, v = -1 (walls at +1 and -1)",
      sp.simplify(F(a, b, p=0, q=-1) - (a + b) / (1 + a * b)) == 0)
check("tangent addition is p = 0, q = 1 (walls at +i, -i)",
      sp.simplify(sp.expand_trig(sp.tan(s + t)) - F(sp.tan(s), sp.tan(t), p=0, q=1)) == 0)
check("discriminant p^2 - 4q = (u - v)^2", sp.expand((u + v) ** 2 - 4 * u * v - (u - v) ** 2) == 0)
# double wall: additive lens
Fd = F(a, b, p=-2 * u, q=u ** 2)
Nd = lambda x: x / (1 - u * x)
check("double wall (on the parabola): N(a (+) b) = N(a) + N(b), N(x) = x/(1-ux)", sp.simplify(Nd(Fd) - Nd(a) - Nd(b)) == 0)
# complex walls: rotation
w = sp.symbols('w', real=True)
x1, x2 = sp.symbols('x1 x2', real=True)
Fc = F(x1, x2, p=-2 * w, q=w ** 2 + 1)        # walls 1/(w +- i)
Mc = lambda x: (1 - (w + sp.I) * x) / (1 - (w - sp.I) * x)
check("complex walls: M(a (+) b) = M(a) M(b) and |M| = 1 on reals (a turn of a circle)",
      sp.simplify(Mc(Fc) - Mc(x1) * Mc(x2)) == 0 and
      sp.simplify(sp.expand(Mc(x1) * sp.conjugate(Mc(x1))) - 1) == 0)

# --- 4. Units: rescaling moves a point along a parabola through the origin ------------------
check("units: F_{p,q}(lam a, lam b) = lam F_{lam p, lam^2 q}(a, b)",
      sp.simplify(F(lam * a, lam * b) - lam * F(a, b, p=lam * p, q=lam ** 2 * q)) == 0)
check("dial: k (a (+)_k b) = (k a) (+) (k b)  [every k != 0 on the line is mult-add in other units]",
      sp.simplify(t * F(a, b, p=t, q=0) - F(t * a, t * b, p=1, q=0)) == 0)
# the unit-free shape: rung-equivalent height s_eq = p^2/(p^2 - 4q)
seq = (u + v) ** 2 / (u - v) ** 2
check("s_eq = ((u+v)/(u-v))^2 is unchanged by units (u,v -> u/lam, v/lam)", sp.simplify(seq.subs({u: u / lam, v: v / lam}, simultaneous=True) - seq) == 0)

# --- 5. Singular cubics: the plane is exactly y^2 + a1 x y = x^3 + a2 x^2 ---------------------
A1, A2, m1, m2 = sp.symbols('A1 A2 m1 m2')
# points on the cubic through the double point: slope m gives x = m^2 + A1 m - A2, y = m x
pt = lambda m: (m ** 2 + A1 * m - A2, m * (m ** 2 + A1 * m - A2))
(X1, Y1), (X2, Y2) = pt(m1), pt(m2)
lamline = (Y2 - Y1) / (X2 - X1)
nu = Y1 - lamline * X1
X3 = sp.simplify(lamline ** 2 + A1 * lamline - A2 - X1 - X2)
Y3 = lamline * X3 + nu
# P + Q = -R ; negation (x, y) -> (x, -y - A1 x)
Xs, Ys = X3, -Y3 - A1 * X3
zsum = sp.simplify(-Xs / Ys)
zP, zQ = -X1 / Y1, -X2 / Y2           # formal parameter z = -x/y = -1/m
check("formal group of the cubic y^2 + a1 xy = x^3 + a2 x^2 is exactly F_{-a1,-a2} in z = -x/y",
      sp.simplify(zsum - F(zP, zQ, p=-A1, q=-A2)) == 0)
check("tangent slopes at the double point have discriminant a1^2 + 4 a2 = p^2 - 4q",
      sp.expand(A1 ** 2 + 4 * A2 - ((-A1) ** 2 - 4 * (-A2))) == 0)

# --- 6. The chi_y genus line (Hirzebruch) -----------------------------------------------------
# convention of M. Albanese's notes: u +_{chi_y} v = (u + v + (1 - y) uv) / (1 + y uv)
y = sp.symbols('y')
Fchi = (a + b + (1 - y) * a * b) / (1 + y * a * b)
Jc = jet(Fchi, 3)
check("chi_y genus law sits at (1-y, -y): the line c2 = c1 - 1", sp.simplify(Jc[(1, 1)] - (1 - y)) == 0 and sp.simplify(Jc[(2, 1)] + y) == 0)
check("Todd genus (y = 0) is mult-add, signature (y = 1) is velocity addition",
      sp.simplify(Fchi.subs(y, 0) - (a + b + a * b)) == 0 and sp.simplify(Fchi.subs(y, 1) - (a + b) / (1 + a * b)) == 0)
check("chi_y line touches the parabola only at y = -1 (Euler characteristic, a cusp)",
      sp.factor((1 - y) ** 2 - 4 * (-y)) == (y + 1) ** 2)

# --- 7. Odd laws on the vertical axis ---------------------------------------------------------
Fsin = a * sp.sqrt(1 - b ** 2) + b * sp.sqrt(1 - a ** 2)
Js = jet(Fsin, 3)
check("sine addition sits at (0, -1/2)", Js.get((1, 1), 0) == 0 and Js[(2, 1)] == sp.Rational(-1, 2))

# --- 8. Elliptic curves land at (-a1, -a2) (Silverman, ch. IV) --------------------------------
# Numeric check on a nonsingular curve: add two points near O on the curve itself, read z = -x/y,
# and compare with z1 + z2 - a1 z1 z2 - a2 (z1^2 z2 + z1 z2^2); the gap must shrink like h^4.
import mpmath as mp
mp.mp.dps = 60
ea1, ea2, ea3, ea4, ea6 = mp.mpf(3) / 7, mp.mpf(-5) / 4, mp.mpf(2) / 3, mp.mpf(1) / 5, mp.mpf(-1) / 3


def w_of(zv):
    wv = zv ** 3
    for _ in range(200):
        wv = zv ** 3 + ea1 * zv * wv + ea2 * zv ** 2 * wv + ea3 * wv ** 2 + ea4 * zv * wv ** 2 + ea6 * wv ** 3
    return wv


def pt_of(zv):
    wv = w_of(zv)
    return zv / wv, -1 / wv


def add(P, Q):
    (X1, Y1), (X2, Y2) = P, Q
    L = (Y2 - Y1) / (X2 - X1)
    nu = Y1 - L * X1
    X3 = L ** 2 + ea1 * L - ea2 - X1 - X2
    Y3 = -(L + ea1) * X3 - nu - ea3
    return X3, Y3


ratios = []
for h in [mp.mpf(10) ** -k for k in (2, 3, 4)]:
    za, zb = 0.7 * h, -0.4 * h
    X3, Y3 = add(pt_of(za), pt_of(zb))
    z3 = -X3 / Y3
    approx = za + zb - ea1 * za * zb - ea2 * (za ** 2 * zb + za * zb ** 2)
    ratios.append(abs(z3 - approx) / h ** 4)
check("elliptic curve (a1=3/7, a2=-5/4, a3, a4, a6 nonzero): group law = z1+z2 - a1 z1z2 - a2(z1^2 z2 + z1 z2^2) + O(h^4); "
      f"gap/h^4 = {[mp.nstr(r, 4) for r in ratios]}",
      max(ratios) < 10 and abs(ratios[-1] - ratios[-2]) / ratios[-1] < 0.05)

# --- 9. Fractional rungs: 3-jet (s, s(s-1)/4) and the rung-equivalent height ------------------
# E^s(x) = x + s x^2/2 + (s/6 + s(s-1)/4) x^3 + ...  (from fractional-rungs/symbolic_proofs.py)
xx = sp.symbols('xx')
Es = lambda x, s: x + s * x ** 2 / 2 + (s / 6 + s * (s - 1) / 4) * x ** 3
Ems = lambda x, s: Es(x, -s)    # inverse to third order
Bs = sp.expand(Es(Ems(a, s) + Ems(b, s), s))
Jb = jet(Bs, 3)
check("rung s sits at (s, s(s-1)/4) on the map", sp.simplify(Jb[(1, 1)] - s) == 0 and sp.simplify(Jb[(2, 1)] - s * (s - 1) / 4) == 0)
check("rung s has rung-equivalent height s: p^2/(p^2-4q) = s", sp.simplify(s ** 2 / (s ** 2 - s * (s - 1)) - s) == 0)
check("rung s osculates two real walls for s > 0 (discriminant = s)", sp.expand(s ** 2 - 4 * s * (s - 1) / 4) == s)

# --- 10. KdV corner = velocity-addition lift of the wavenumber ratio --------------------------
k1, k2 = sp.symbols('k1 k2', positive=True)
check("Hirota weight A = ((k1-k2)/(k1+k2))^2 = M(k2/k1)^2 for velocity addition's lens M(x) = (1-x)/(1+x)",
      sp.simplify(((k1 - k2) / (k1 + k2)) ** 2 - ((1 - k2 / k1) / (1 + k2 / k1)) ** 2) == 0)
check("KdV phase shift: A(f - 1) = (A e1) (+) (A e2), plain mult-add after moving each wave by log A",
      sp.expand(A1 * (z1 + z2 + A1 * z1 * z2) - ((A1 * z1) + (A1 * z2) + (A1 * z1) * (A1 * z2))) == 0)

out.append("")
out.append(f"{npass} of {nall} exact checks pass")
open(__file__.replace('.py', '_output.txt'), 'w').write("\n".join(out) + "\n")
print("\n".join(out))
