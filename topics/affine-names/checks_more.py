"""Checks for the second pass: which other affine transformations add something to the dials.

Plane numbers: eps^2 = k*eps - d, kappa = d - k^2/4, norm N(x + y eps) = x^2 + k x y + d y^2.
"""
import numpy as np
import sympy as sp

results = []


def check(name, ok, detail=""):
    results.append((name, bool(ok)))
    print(f"{'PASS' if ok else 'FAIL'}  {name}  {detail}")


k, d, x, y, t = sp.symbols('k d x y t', real=True)
kap = d - k**2 / 4

# multiplication by x + y eps as a matrix on (1, eps): eps*1 = eps, eps*eps = -d + k eps
def mat(a, b):
    return sp.Matrix([[a, -d * b], [b, a + k * b]])

E = mat(0, 1)          # multiplication by eps
E0 = E - k / 2 * sp.eye(2)

# CA1: equi-centro-affine curvature of the unit circle N = 1 equals kappa.
# Curve z(t) = exp(t E0) z0 stays on N = 1 (E0 is trace free and preserves N up to scale... check directly).
# curvature formula (SL2 acting linearly): det(z', z'') / det(z, z')^3
for kv, dv in [(0, 1), (1, 1), (1.3, 2.1), (0.5, -0.7), (2, 0.2), (-1.1, 0.9)]:
    E0n = np.array(E0.subs({k: kv, d: dv}), dtype=float)
    kapn = dv - kv**2 / 4
    z0 = np.array([1.0, 0.0])            # N(1) = 1
    from scipy.linalg import expm
    def z(tt):
        return expm(tt * E0n) @ z0
    h = 1e-4
    tt = 0.37
    z1, zp, zm = z(tt), z(tt + h), z(tt - h)
    zd = (zp - zm) / (2 * h)
    zdd = (zp - 2 * z1 + zm) / h**2
    det = lambda u, v: u[0] * v[1] - u[1] * v[0]
    curv = det(zd, zdd) / det(z1, zd)**3
    Nval = z1[0]**2 + kv * z1[0] * z1[1] + dv * z1[1]**2
    check(f"CA1 unit circle k={kv}, d={dv}", abs(curv - kapn) < 1e-5 and abs(Nval - 1) < 1e-9,
          f"curvature {curv:.6f} vs kappa {kapn:.6f}, N {Nval:.9f}")

# CA2: symbolic: for the ellipse N = 1, area = pi / sqrt(kappa), so (pi/area)^2 = kappa
# (equi-centro-affine curvature of a centred ellipse is (pi/area)^2).
M = sp.Matrix([[1, k / 2], [k / 2, d]])
area = sp.pi / sp.sqrt(M.det())
check("CA2 ellipse N=1 has (pi/area)^2 = kappa", sp.simplify((sp.pi / area)**2 - kap) == 0)

# CA3: the curvature dial's one-parameter group exp(t E0) has det 1 and is a rotation / shear / squeeze
for kv, dv, kind in [(1, 1, "rotation"), (2, 1, "shear"), (1, -1, "squeeze")]:
    E0n = np.array(E0.subs({k: kv, d: dv}), dtype=float)
    from scipy.linalg import expm
    g = expm(0.8 * E0n)
    tr = np.trace(g)
    kind_found = "rotation" if abs(tr) < 2 - 1e-9 else ("shear" if abs(abs(tr) - 2) < 1e-9 else "squeeze")
    check(f"CA3 kappa={dv - kv**2/4:+.2f} gives {kind}", abs(np.linalg.det(g) - 1) < 1e-12 and kind_found == kind,
          f"det {np.linalg.det(g):.12f}, trace {tr:.4f}")

# CA4: the 1-D moves T_a(x) = (1 + k a) x + a fix the wall -1/k, not 0; in u = 1 + kx they are u -> (1+ka) u
a = sp.symbols('a', real=True)
T = (1 + k * a) * x + a
check("CA4 T_a fixes the wall -1/k", sp.simplify(T.subs(x, -1 / k) + 1 / k) == 0)
u = 1 + k * x
check("CA4b in u = 1 + kx the move is u -> (1+ka)u", sp.simplify((1 + k * T) - (1 + k * a) * u) == 0)

# EQ1: area growth rate of 1 + t eps at t = 0 is k (trace); the trace-free part E0 keeps area
check("EQ1 d/dt det(1 + t eps) at 0 = k", sp.simplify(sp.diff((sp.eye(2) + t * E).det(), t).subs(t, 0) - k) == 0)
check("EQ1b trace(E0) = 0, so exp(t E0) is equi-affine", sp.simplify(E0.trace()) == 0)

# RF1: conjugation eps -> k - eps is an affine (oblique) reflection: det -1, squares to identity, keeps N
C = sp.Matrix([[1, k], [0, -1]])   # (x, y) -> (x + k y, -y)
Nf = lambda v: v[0]**2 + k * v[0] * v[1] + d * v[1]**2
v = sp.Matrix([x, y])
check("RF1 conjugation is a reflection keeping N",
      C.det() == -1 and sp.simplify(C * C - sp.eye(2)) == sp.zeros(2) and sp.expand(Nf(C * v) - Nf(v)) == 0)
# RF2: the half-turn about the wall is the element a = -2/k
check("RF2 a = -2/k gives ratio -1 (half-turn about the wall)", sp.simplify(1 + k * (-2 / k) + 1) == 0)

# PR1: the chord conic x^2 + kxy + dy^2 = y meets the line at infinity where x^2 + kxy + dy^2 = 0:
# discriminant k^2 - 4d = -4 kappa, so 0, 1, 2 points for kappa > 0, = 0, < 0
disc = k**2 - 4 * d
check("PR1 points at infinity counted by -4 kappa", sp.simplify(disc + 4 * kap) == 0)

# PR2: Laguerre: the angle of a plane number is (1/2i) log of the cross ratio with the two null directions.
def laguerre(kv, dv, phi):
    kapn = dv - kv**2 / 4
    E0n = np.array(E0.subs({k: kv, d: dv}), dtype=complex)
    from scipy.linalg import expm
    th = phi / np.sqrt(abs(kapn))
    w = expm(th * E0n) @ np.array([1, 0], dtype=complex)
    # slopes y/x of the four directions: 1 (slope 0), w, and the null directions 1 + k m + d m^2 = 0
    m1, m2 = np.roots([dv, kv, 1]).astype(complex)
    mw = w[1] / w[0]
    cr = ((0 - m1) * (mw - m2)) / ((0 - m2) * (mw - m1))
    return cr, kapn
cr, kapn = laguerre(1.0, 1.5, 0.6)
ang = np.log(cr) / 2j
check("PR2 Laguerre angle, kappa > 0", abs(abs(ang.real) - 0.6) < 1e-9, f"angle {ang} (expect ±0.6)")
cr, kapn = laguerre(1.0, -0.5, 0.6)
ang = np.log(cr) / 2
check("PR2b Laguerre angle, kappa < 0 (hyperbolic, log/2)", abs(abs(ang.real) - 0.6) < 1e-9, f"angle {ang} (expect ±0.6)")

# PR3: projective maps take any ellipse to any hyperbola (so projective geometry forgets kappa):
# x^2 + y^2 = z^2 and x^2 - y^2 = ... : the circle x^2 + y^2 = 1 becomes the hyperbola under (x, y, z) -> (z, y, x)
X, Y, Z = sp.symbols('X Y Z')
circle = X**2 + Y**2 - Z**2
hyper = circle.subs({X: Z, Z: X}, simultaneous=True)   # Z^2 + Y^2 - X^2  -> X^2 - Y^2 = 1 in z = 1
check("PR3 one projective swap turns the circle into a hyperbola", sp.expand(-hyper.subs(Z, 1)) == sp.expand(X**2 - Y**2 - 1))

# DL1: the discriminant reading, with the repo's sign (eps^2 = k eps - d, d > 0 round at k = 0)
Delta = k**2 - 4 * d
check("DL1 curvature = -Delta/4 = d - k^2/4", sp.simplify(-Delta / 4 - kap) == 0)
check("DL2 trace-free part squares to Delta/4", sp.simplify(E0 * E0 - Delta / 4 * sp.eye(2)) == sp.zeros(2))
check("DL3 k = 0, d = 1 (round, like i) has Delta < 0; d = -1 (relativity side) has Delta > 0",
      Delta.subs({k: 0, d: 1}) < 0 and Delta.subs({k: 0, d: -1}) > 0)

n = sum(ok for _, ok in results)
print(f"\n{n} of {len(results)} pass")
