"""DRAFT checks: do the homothety k and the curvature dial d come from one affine sphere?

Plane law F(a,b) = (a+b+kab)/(1-dab).  Inputs: k, d.  Readings out: Delta = k^2-4d, curvature c = d-k^2/4 = -Delta/4.
Shape tested: the chord surface  x^2 + kxy + dy^2 + w^2 = y  (the chord conic with one extra round direction),
rationally parametrised through the identity O=(0,0,0) by (u, u^2+v^2, v)/q, q = 1 + k u + d (u^2+v^2).
"""
import sympy as sp, random, time
from blaschke import Surface

random.seed(7)
u, v, s, t = sp.symbols('u v s t', real=True)
K_, D_ = sp.symbols('k d', real=True)
results = []


def check(name, ok, detail=""):
    results.append((name, bool(ok)))
    print(("PASS " if ok else "FAIL ") + name + ("  | " + detail if detail else ""))


def chord_surface(k, d):
    q = 1 + k * u + d * (u ** 2 + v ** 2)
    return [u / q, (u ** 2 + v ** 2) / q, v / q]


def chord_curve(k, d):
    q = 1 + k * u + d * u ** 2
    return [u / q, u ** 2 / q]


# ---------- S1: the parametrisation lies on the surface and passes through O, for symbolic k, d
x, y, w = chord_surface(K_, D_)
lhs = sp.simplify(x ** 2 + K_ * x * y + D_ * y ** 2 + w ** 2 - y)
check("S1 (u,u^2+v^2,v)/q lies on x^2+kxy+dy^2+w^2=y for all k,d, and u=v=0 is O", lhs == 0,
      "residual " + str(lhs))

# ---------- S2: the k-shear removes k exactly (it is the Cayley change of reading)
c_ = D_ - K_ ** 2 / 4
X = x + K_ * y / 2
res = sp.simplify(X ** 2 + c_ * y ** 2 + w ** 2 - y)
shear = sp.Matrix([[1, K_ / 2, 0], [0, 1, 0], [0, 0, 1]])
xc, yc = chord_curve(K_, D_)
T = sp.simplify((yc) / (xc + K_ * yc / 2))
check("S2 shear X=x+(k/2)y (det 1) carries the (k,d) surface onto X^2+(d-k^2/4)y^2+w^2=y",
      res == 0 and shear.det() == 1, "residual " + str(res))
check("S2b on the curve the shear is the Cayley map T = a/(1+(k/2)a)",
      sp.simplify(T - u / (1 + K_ * u / 2)) == 0, "T = " + str(T))

# ---------- S3/S4/S5: affine normal, affine mean curvature, cubic form on the surface
cases = [(sp.Rational(3, 5), sp.Rational(7, 10)),    # c > 0  round side
         (sp.Rational(-6, 5), sp.Rational(1, 2)),    # c > 0
         (sp.Rational(2, 1), sp.Rational(1, 1)),     # c = 0  wall (parabolic)
         (sp.Rational(1, 2), sp.Rational(1, 16)),    # c = 0
         (sp.Rational(1, 1), sp.Rational(-1, 2)),    # c < 0  two-wall side
         (sp.Rational(-3, 2), sp.Rational(1, 5))]    # c < 0
pts = [{u: 0, v: 0}, {u: sp.Rational(1, 3), v: sp.Rational(-1, 4)}, {u: sp.Rational(-1, 5), v: sp.Rational(2, 7)}]
rows = []
for (k, d) in cases:
    c = d - k ** 2 / 4
    S = Surface(chord_surface(k, d), [u, v], {u: 0, v: 0})
    vals = [S.at(p) for p in pts]
    H_pred = 2 ** sp.Rational(3, 2) * c
    xi0 = vals[0][1]
    if c != 0:
        center = sp.Matrix([-k, 2, 0]) / (4 * d - k ** 2)
        err = max(max(abs(e) for e in (xi + H_pred * (Xp - center))) for (Xp, xi, h, K2) in vals)
        what = "xi = -H (x - centre), H = 2^(3/2)(d-k^2/4)"
    else:
        err = max(max(abs(e) for e in (xi - xi0)) for (Xp, xi, h, K2) in vals)
        what = "xi constant (improper sphere), H = 0"
    lean_err = max(abs(e) for e in (xi0 - sp.sqrt(2) * sp.Matrix([-k / 2, 1, 0])))
    K2max = max(abs(K2) for (_, _, _, K2) in vals)
    rows.append((k, d, c, err, lean_err, K2max))
    check(f"S3 (k,d)=({k},{d}) c={c}: {what}", err < 1e-20, f"max err {sp.N(err, 3)}")
    check(f"S4 (k,d)=({k},{d}): affine normal at O = sqrt2*(-k/2, 1, 0), no d in it", lean_err < 1e-20,
          f"xi(O) = {[sp.N(e, 8) for e in xi0]}")
    check(f"S5 (k,d)=({k},{d}): cubic form vanishes (|K|^2_h = 0)", K2max < 1e-20, f"max |K|^2 {sp.N(K2max, 3)}")

# ---------- S6: read (k, d) back off the shape alone (normal at O + mean curvature)
recov = 0
for (k, d) in cases:
    c = d - k ** 2 / 4
    S = Surface(chord_surface(k, d), [u, v], {u: 0, v: 0})
    X1, xi1, _, _ = S.at(pts[1])
    X0, xi0, _, _ = S.at(pts[0])
    k_read = -2 * xi0[0] / xi0[1]
    # H from two points: xi1 - xi0 = -H (X1 - X0)
    H_read = -(xi1 - xi0).dot(X1 - X0) / (X1 - X0).dot(X1 - X0)
    d_read = H_read / 2 ** sp.Rational(3, 2) + k_read ** 2 / 4
    recov = max(recov, abs(k_read - k), abs(d_read - d))
check("S6 k = -2 xi_x/xi_y at O and d = H/2^(3/2) + k^2/4 recover both dials from the shape", recov < 1e-20,
      f"max error {sp.N(recov, 3)}")

# ---------- S7: same in one dimension (the chord conic is an affine circle)
curve_err = 0
for (k, d) in cases:
    c = d - k ** 2 / 4
    S = Surface(chord_curve(k, d), [u], {u: 0})
    for p in [{u: 0}, {u: sp.Rational(2, 5)}, {u: sp.Rational(-1, 3)}]:
        Xp, xi, h, _ = S.at(p)
        if c != 0:
            center = sp.Matrix([-k, 2]) / (4 * d - k ** 2)
            curve_err = max(curve_err, max(abs(e) for e in (xi + 2 ** sp.Rational(4, 3) * c * (Xp - center))))
        if p[u] == 0:
            curve_err = max(curve_err, max(abs(e) for e in (xi - 2 ** sp.Rational(1, 3) * sp.Matrix([-k / 2, 1]))))
check("S7 chord conic: xi = -kappa(x - centre), kappa = 2^(4/3)(d-k^2/4); xi(O) = 2^(1/3)(-k/2,1)", curve_err < 1e-20,
      f"max err {sp.N(curve_err, 3)}")
# general n: H = 4^((n+1)/(n+2)) c from semi-axes 1/(2 sqrt c) (n times) and 1/(2c)
n_, cc = sp.symbols('n c', positive=True)
prod = (1 / (2 * sp.sqrt(cc))) ** n_ / (2 * cc)
Hn = sp.simplify(sp.expand_power_base(prod ** (-2 / (n_ + 2)), force=True))
check("S7b any dimension n: H = 4^((n+1)/(n+2)) (d-k^2/4)",
      sp.simplify(sp.expand_power_base(Hn / (4 ** ((n_ + 1) / (n_ + 2)) * cc), force=True)) == 1, f"H = {Hn}")

# ---------- S8: a shape with a live cubic form: the Titeica surface xyz = 1
Tit = Surface([sp.exp(s), sp.exp(t), sp.exp(-s - t)], [s, t], {s: 0, t: 0})
tv = [Tit.at(p) for p in [{s: 0, t: 0}, {s: sp.Rational(1, 2), t: sp.Rational(-1, 3)}, {s: -1, t: sp.Rational(3, 4)}]]
ratios = [[xi[i] / Xp[i] for i in range(3)] for (Xp, xi, _, _) in tv]
r0 = ratios[0][0]
sphere_err = max(abs(r - r0) for row in ratios for r in row)
K2s = [K2 for (_, _, _, K2) in tv]
check("S8 xyz=1 is an affine sphere centred at 0 (xi = -H x, H constant)", sphere_err < 1e-20,
      f"H = {sp.N(-r0, 10)}")
check("S8b xyz=1 has a nonzero constant cubic form (|K|^2_h same at 3 points)",
      abs(K2s[0]) > 0.1 and max(abs(a - K2s[0]) for a in K2s) < 1e-20, f"|K|^2 = {sp.N(K2s[0], 10)}")

# ---------- S9: the matrix-orbit alternative gives the right split but the wrong types
Xg = sp.Matrix([[K_, 1], [-D_, 0]])
X0 = Xg - K_ / 2 * sp.eye(2)
a1, b1 = sp.symbols('a b')
lhs = (sp.eye(2) + a1 * Xg) * (sp.eye(2) + b1 * Xg)
rhs = (1 - D_ * a1 * b1) * sp.eye(2) + (a1 + b1 + K_ * a1 * b1) * Xg
check("S9 generator X=[[k,1],[-d,0]]: (I+aX)(I+bX) = (1-dab)I + (a+b+kab)X; tr k, det d, det X0 = d-k^2/4",
      sp.simplify(lhs - rhs) == sp.zeros(2) and Xg.trace() == K_ and Xg.det() == D_
      and sp.simplify(X0.det() - c_) == 0)
# det on traceless matrices [[p, q+r],[q-r, -p]] = -p^2 - q^2 + r^2 : signature (1 plus, 2 minus)
p_, q_, r_ = sp.symbols('p q r')
quad = sp.expand(sp.Matrix([[p_, q_ + r_], [q_ - r_, -p_]]).det())
check("S9b det on sl(2) is -p^2-q^2+r^2: round side (det>0) is a TWO-sheeted hyperboloid (hyperbolic type), "
      "two-wall side a one-sheeted one (indefinite): types do not match the plane dial",
      quad == -p_ ** 2 - q_ ** 2 + r_ ** 2, str(quad))

# ---------- S10 polish: fit the slope of H against c from float samples, then re-evaluate at 40 digits
samples = [(sp.Rational(random.randint(-20, 20), 10), sp.Rational(random.randint(-10, 20), 10)) for _ in range(8)]
samples = [(k, d) for (k, d) in samples if d - k ** 2 / 4 != 0][:6]
t0 = time.time()
hs, cs = [], []
for (k, d) in samples:
    S = Surface(chord_surface(k, d), [u, v], {u: 0, v: 0})
    X0p, xi0, _, _ = S.at(pts[0], prec=15)
    X1p, xi1, _, _ = S.at({u: sp.Rational(1, 10), v: sp.Rational(1, 20)}, prec=15)
    hs.append(float(-(xi1 - xi0).dot(X1p - X0p) / (X1p - X0p).dot(X1p - X0p)))
    cs.append(float(d - k ** 2 / 4))
slope15 = sum(a * b for a, b in zip(hs, cs)) / sum(b * b for b in cs)
t1 = time.time()
hs40 = []
for (k, d) in samples:
    S = Surface(chord_surface(k, d), [u, v], {u: 0, v: 0})
    X0p, xi0, _, _ = S.at(pts[0], prec=45)
    X1p, xi1, _, _ = S.at({u: sp.Rational(1, 10), v: sp.Rational(1, 20)}, prec=45)
    hs40.append(-(xi1 - xi0).dot(X1p - X0p) / (X1p - X0p).dot(X1p - X0p))
cs40 = [d - k ** 2 / 4 for (k, d) in samples]
slope40 = sum(a * b for a, b in zip(hs40, cs40)) / sum(b * b for b in cs40)
t2 = time.time()
target = 2 ** sp.Rational(3, 2)
e15, e40 = abs(slope15 - float(target)), abs(sp.N(slope40 - target, 50))
check("S10 polish: slope of H against d-k^2/4 is 2^(3/2)", e40 < 1e-35,
      f"15-digit fit err {e15:.1e} ({t1 - t0:.0f} s) -> 45-digit re-evaluation err {sp.N(e40, 3)} ({t2 - t1:.0f} s)")

print()
print("table k, d, c, sphere err, lean err, |K|^2")
for r in rows:
    print("  ", [str(r[0]), str(r[1]), str(r[2])] + [f"{float(e):.1e}" for e in r[3:]])
print(f"\n{sum(ok for _, ok in results)}/{len(results)} checks pass")
