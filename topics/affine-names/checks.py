"""Checks for the affine names of the four dials (DRAFT, Nov 26 2026).

Run: python3 checks.py > checks_output.txt
"""
import numpy as np
import sympy as sp
from scipy.integrate import quad

results = []


def check(name, ok, msg):
    results.append(ok)
    print(f"{name}: {'PASS' if ok else 'FAIL'}  {msg}")


a, b, x, k, d, t, th = sp.symbols('a b x k d t theta', real=True)
circ = a + b + k * a * b

# H1 the circle operation is composition of homotheties about the wall
T = lambda s: (1 + k * s) * x + s
comp = sp.expand(T(a).subs(x, T(b)) - T(circ))
fixed = sp.simplify(T(a).subs(x, -1 / k) + 1 / k)
check("H1", comp == 0 and fixed == 0,
      "T_a(x) = (1+ka)x + a gives T_a∘T_b = T_(a∘ₖb) (sympy exact), and every T_a fixes x = -1/k: "
      "the circle operation is the law of homotheties with centre at the wall -1/k and ratio 1+ka")

# H2 k = 0: homotheties become translations (centre at infinity)
check("H2", sp.simplify(T(a).subs(k, 0) - (x + a)) == 0,
      "at k = 0, T_a(x) = x + a: translations, the homotheties whose centre has gone to infinity")

# H3 rescaling the ruler: k is not an affine invariant on its own, only its sign/zero is
lam = sp.symbols('lambda', positive=True)
resc = sp.simplify((lam * a + lam * b + (k / lam) * (lam * a) * (lam * b)) / lam - circ)
check("H3", resc == 0,
      "a -> λa turns ∘ₖ into ∘_(k/λ): every k ≠ 0 is the same homothety law seen with another ruler; "
      "only k = 0 (translations) versus k ≠ 0 (a finite centre) is affine")

# C1 plane numbers: eps^2 = k eps - d. trace k, det d, trace-free part squares to -kappa
E = sp.Matrix([[0, -d], [1, k]])  # multiplication by eps on (1, eps)
E0 = E - (k / 2) * sp.eye(2)
kappa = d - k**2 / 4
ok = (sp.simplify(E**2 - (k * E - d * sp.eye(2))) == sp.zeros(2)
      and sp.simplify(E.trace() - k) == 0 and sp.simplify(E.det() - d) == 0
      and sp.simplify(E0 * E0 + kappa * sp.eye(2)) == sp.zeros(2))
check("C1", ok,
      "multiplying by ε has trace k and determinant d; splitting off the homothety (k/2)·I leaves an "
      "area-preserving part ε₀ with ε₀² = -κ, κ = d - k²/4 (sympy exact)")

# C2 the area factor of multiplying by 1 + t eps grows at rate k (trace = homothety rate)
N = sp.expand((sp.Matrix([[1, 0], [0, 1]]) + t * E).det())
check("C2", sp.simplify(N - (1 + k * t + d * t**2)) == 0 and sp.diff(N, t).subs(t, 0) == k,
      "area factor of 1 + tε is 1 + kt + dt², whose rate at t = 0 is k: k is the homothety (area-growth) part")


# C3 equi-affine curvature of the chord conic x^2 + kxy + dy^2 = y is 2^(4/3) kappa
def affine_curvature(kk, dd, s0):
    s = sp.symbols('s', real=True)
    q = 1 + kk * s + dd * s**2
    P = sp.Matrix([s / q, s**2 / q])
    P1, P2 = P.diff(s), P.diff(s, 2)
    det2 = lambda u, v: u[0] * v[1] - u[1] * v[0]
    w = det2(P1, P2)  # dσ = w^(1/3) ds
    # curvature = [x'', x'''] in affine arclength, computed by reparametrising numerically
    g = sp.lambdify(s, [P[0], P[1]], 'mpmath')
    import mpmath as mp
    mp.mp.dps = 30
    wf = sp.lambdify(s, w, 'mpmath')

    def X(sig_s):
        return mp.matrix(g(sig_s))

    # derivatives with respect to affine arclength via chain rule, numerically
    ds = lambda f, s_: mp.diff(f, s_)
    sp_ = lambda s_: mp.cbrt(wf(s_))  # dσ/ds
    D = lambda f: (lambda s_: ds(f, s_) / sp_(s_))
    fx = lambda s_: X(s_)[0]
    fy = lambda s_: X(s_)[1]
    x2, y2 = D(D(fx)), D(D(fy))
    x3, y3 = D(x2), D(y2)
    return float(x2(s0) * y3(s0) - y2(s0) * x3(s0))


ok = True
msgs = []
for kk, dd in [(0, 1), (1, 1), (0, -1), (1, -0.5), (2, 1), (-1, 0)]:
    kap = affine_curvature(sp.nsimplify(kk), sp.nsimplify(dd), 0.3)
    target = 2 ** (4 / 3) * (dd - kk**2 / 4)
    ok &= abs(kap - target) < 1e-6
    msgs.append(f"({kk},{dd}): {kap:.6f} vs {target:.6f}")
check("C3", ok, "equi-affine curvature of x²+kxy+dy² = y equals 2^(4/3)·κ, κ = d − k²/4 (sign: >0 ellipse, 0 parabola, <0 hyperbola): "
      + "; ".join(msgs))

# C4 sign convention: k = 0 relativity d = -1/c^2 is hyperbolic, circle d = 1 is round
print("Note (convention, not a check): κ > 0 round (the sphere's 1/R², eps₀² = −κ like i² = −1), κ = 0 flat, κ < 0 hyperbolic (relativity κ = −1/c²). "
      "The repo's d already has this sign when k = 0; the i² = d convention of Yaglom's tables is the mirror")

# P1 an apolar (trace-free) cubic form in the plane is pure cos 3θ, sin 3θ: the trefoil pattern
A, B, Cc, Ee = sp.symbols('A B C E', real=True)
u1, u2 = sp.cos(th), sp.sin(th)
cubic = A * u1**3 + 3 * B * u1**2 * u2 + 3 * Cc * u1 * u2**2 + Ee * u2**3
apolar = cubic.subs({Cc: -A, Ee: -B})  # C_i11 + C_i22 = 0
ok = sp.simplify(sp.expand_trig(apolar - (A * sp.cos(3 * th) + B * sp.sin(3 * th)))) == 0
check("P1", ok,
      "Blaschke's apolarity (the cubic form is trace-free against the affine metric) leaves, in 2D, only "
      "A cos 3θ + B sin 3θ (sympy exact): the cubic form is a three-lobed field, the same angular pattern as Zernike's trefoil")

# P2 an odd tensor vanishes at a centre of symmetry: centrally symmetric shapes have C = 0 at the centre
# (C(-u,-u,-u) = -C(u,u,u) and symmetry maps C to itself) -> checked on the cubic above
ok = sp.simplify(cubic.subs(th, th + sp.pi) + cubic) == 0
check("P2", ok,
      "a cubic form is odd under u -> -u, so a shape with a centre of symmetry (every ℓp ball) has cubic form 0 at its centre; "
      "this is why the simplex page found every ℓp ball and even polygon 'round at the centre' and only the triangle not")


# L1 affine isoperimetric ratio of the l^p circle along p, and its symmetry under p <-> p/(p-1)
def lp_ratio(p):
    # first quadrant: x = c^(2/p), y = s^(2/p), t in (0, pi/2)
    def xy(tt):
        c, s = np.cos(tt), np.sin(tt)
        return c ** (2 / p), s ** (2 / p)

    def derivs(tt):
        c, s = np.cos(tt), np.sin(tt)
        e = 2 / p
        x1 = -e * c ** (e - 1) * s
        y1 = e * s ** (e - 1) * c
        x2 = -e * ((e - 1) * c ** (e - 2) * (-s) * s + c ** (e - 1) * c)
        y2 = e * ((e - 1) * s ** (e - 2) * c * c + s ** (e - 1) * (-s))
        return x1, y1, x2, y2

    def integrand_aff(tt):
        x1, y1, x2, y2 = derivs(tt)
        return abs(x1 * y2 - y1 * x2) ** (1 / 3)

    def integrand_area(tt):
        x0, y0 = xy(tt)
        x1, y1, _, _ = derivs(tt)
        return 0.5 * (x0 * y1 - y0 * x1)

    om = 4 * quad(integrand_aff, 0, np.pi / 2, limit=400)[0]
    ar = 4 * quad(integrand_area, 0, np.pi / 2, limit=400)[0]
    return om**3 / (8 * np.pi**2 * ar)


ps = [1.05, 1.25, 1.5, 2.0, 3.0, 5.0, 21.0]
vals = {p: lp_ratio(p) for p in ps}
ok = abs(vals[2.0] - 1) < 1e-8 and all(vals[p] < 1 for p in ps if p != 2.0) and vals[1.05] < vals[1.5] and vals[21.0] < vals[3.0]
check("L1", ok, "Blaschke's affine isoperimetric ratio Ω³/(8π²A) of the ℓp circle: "
      + ", ".join(f"p={p}: {v:.4f}" for p, v in vals.items())
      + "; 1 only at p = 2 (ellipse), falling toward 0 at both polygon ends")
pairs = [(1.5, 3.0), (1.25, 5.0), (1.05, 21.0)]
diffs = [(p, q, lp_ratio(p), lp_ratio(q)) for p, q in pairs]
sym = all(abs(r1 - r2) < 1e-6 for _, _, r1, r2 in diffs)
print("L2 (measurement, not counted): polarity test, p versus its partner p/(p−1): "
      + "; ".join(f"{p}↔{q}: {r1:.6f} vs {r2:.6f}" for p, q, r1, r2 in diffs)
      + (" — equal: the ratio cannot tell a ball from its polar" if sym else
         " — not equal: affine perimeter does tell a ball from its polar, so p and p/(p−1) are different points of the dial (reported, not a failure)"))


# L3 the l^p circle departs from round in a four-lobed pattern (cos 4θ), not three-lobed
def lp_radius(p, thv):
    c, s_ = np.abs(np.cos(thv)), np.abs(np.sin(thv))
    return (c**p + s_**p) ** (-1 / p)


thv = np.linspace(0, 2 * np.pi, 4096, endpoint=False)
msgs = []
ok = True
for p in [1.8, 2.2, 3.0]:
    r = lp_radius(p, thv)
    co = np.abs(np.fft.rfft(r - r.mean())) / len(thv)
    low = co[1:4].max()
    ok &= co[4] > 100 * max(low, 1e-18) and co[4] > co[5:].max()
    msgs.append(f"p={p}: |c4| = {co[4]:.2e}, max |c1..c3| = {low:.1e}")
check("L3", ok, "Fourier harmonics of the ℓp circle's radius: the first non-round harmonic is 4 (a four-lobed, quartic departure), "
      "never 3: " + "; ".join(msgs) + ". The cubic form is the three-lobed departure, so the two shape dials are different directions")

print(f"\n{sum(results)} of {len(results)} pass")
