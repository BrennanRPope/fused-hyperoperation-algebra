"""Precession through the circle operation (DRAFT, Nov 25 2026).

Sections
  A  Swap = turn: the circle operation on quaternions (exact, sympy + numpy)
  B  One composition law, two ends: rotations (d = +1) and boosts (d = -1)
  C  Twist = d x area: Foucault's pendulum from the d = +1 law, Thomas from d = -1
  D  Free rigid body: Montgomery's phase, orientation built only by composing turns
  E  The torqued gyroscope (heavy symmetric top)
  F  Larmor precession and Berry's half angle
  G  Frame dragging in the spinning (rotating) BTZ black hole

Writes checks_output.txt (via stdout redirect) and results.json (for figures.py).
"""
import json, math, sys
import numpy as np
import sympy as sp
from scipy.integrate import solve_ivp

rng = np.random.default_rng(20261125)
RESULTS = {}
PASS = []


def check(name, ok, detail=""):
    PASS.append(bool(ok))
    print(("PASS " if ok else "FAIL ") + name + (f"   [{detail}]" if detail else ""))


# ---------------------------------------------------------------- quaternion tools
def qmul(a, b):
    w1, x1, y1, z1 = a; w2, x2, y2, z2 = b
    return np.array([w1*w2 - x1*x2 - y1*y2 - z1*z2,
                     w1*x2 + x1*w2 + y1*z2 - z1*y2,
                     w1*y2 - x1*z2 + y1*w2 + z1*x2,
                     w1*z2 + x1*y2 - y1*x2 + z1*w2])


def qconj(a):
    return np.array([a[0], -a[1], -a[2], -a[3]])


def qinv(a):
    return qconj(a) / np.dot(a, a)


def circ(a, b, k=1.0):
    """The circle operation a o_k b = a + b + k ab (any associative algebra)."""
    return a + b + k * qmul(a, b)


ONE = np.array([1.0, 0, 0, 0])


def pure(v):
    return np.array([0.0, *v])


def rotate(q, v):
    """Rotate 3-vector v by unit quaternion q (q v q*)."""
    return qmul(qmul(q, pure(v)), qconj(q))[1:]


def rotor_gibbs(axis_angle):
    """Rotation by |w| about w, as 1 + g with g = tan(|w|/2) w_hat (Gibbs / Rodrigues vector)."""
    th = np.linalg.norm(axis_angle)
    if th == 0:
        return np.zeros(3)
    return math.tan(th / 2) * axis_angle / th


def rotmat(axis_angle):
    th = np.linalg.norm(axis_angle)
    if th == 0:
        return np.eye(3)
    n = axis_angle / th
    K = np.array([[0, -n[2], n[1]], [n[2], 0, -n[0]], [-n[1], n[0], 0]])
    return np.eye(3) + math.sin(th) * K + (1 - math.cos(th)) * K @ K


print("== A. Swap = turn: the circle operation on quaternions")
# A1-A3 exact with sympy symbols
a = sp.symbols('a0:4'); b = sp.symbols('b0:4'); c = sp.symbols('c0:4'); k = sp.symbols('k')
def smul(x, y):
    w1, x1, y1, z1 = x; w2, x2, y2, z2 = y
    return [w1*w2 - x1*x2 - y1*y2 - z1*z2, w1*x2 + x1*w2 + y1*z2 - z1*y2,
            w1*y2 - x1*z2 + y1*w2 + z1*x2, w1*z2 + x1*y2 - y1*x2 + z1*w2]
def scirc(x, y):
    p = smul(x, y)
    return [x[i] + y[i] + k * p[i] for i in range(4)]
def slift(x):
    return [1 + k * x[0]] + [k * x[i] for i in range(1, 4)]
lhs = slift(scirc(a, b)); rhs = smul(slift(a), slift(b))
check("A1 lift survives non-commutativity: 1 + k(a o b) = (1 + ka)(1 + kb) for quaternions (sympy, exact)",
      all(sp.expand(lhs[i] - rhs[i]) == 0 for i in range(4)))
l1 = scirc(scirc(a, b), c); l2 = scirc(a, scirc(b, c))
check("A2 the circle operation stays associative on quaternions (exact)", all(sp.expand(l1[i] - l2[i]) == 0 for i in range(4)))
d1 = [scirc(a, b)[i] - scirc(b, a)[i] for i in range(4)]
d2 = [k * (smul(a, b)[i] - smul(b, a)[i]) for i in range(4)]
check("A3 swap difference: a o b - b o a = k (ab - ba) (exact)", all(sp.expand(d1[i] - d2[i]) == 0 for i in range(4)))
u = sp.symbols('u1:4'); v = sp.symbols('v1:4')
pu = [0, *u]; pv = [0, *v]
comm = [smul(pu, pv)[i] - smul(pv, pu)[i] for i in range(4)]
cr = sp.Matrix(u).cross(sp.Matrix(v))
check("A4 two directions: uv - vu = 2 (u x v), real part 0 (exact): the swap points at 90 deg to both inputs",
      sp.expand(comm[0]) == 0 and all(sp.expand(comm[i + 1] - 2 * cr[i]) == 0 for i in range(3)))

errs_turn, errs_norm, angles = [], [], []
for _ in range(2000):
    A = rng.normal(size=4) * 0.6; B = rng.normal(size=4) * 0.6; kk = rng.uniform(-2, 2)
    X = ONE + kk * A
    ab = circ(A, B, kk); ba = circ(B, A, kk)
    errs_turn.append(np.max(np.abs(ba - qmul(qmul(qinv(X), ab), X))))
    errs_norm.append(abs(np.dot(ONE + kk * ab, ONE + kk * ab) - np.dot(ONE + kk * ba, ONE + kk * ba)))
check("A5 swap = turn: b o a = X^-1 (a o b) X with X = 1 + ka (2000 random quaternions and k)",
      max(errs_turn) < 1e-12, f"max err {max(errs_turn):.1e}")
check("A6 swapping keeps the size: |1 + k(a o b)| = |1 + k(b o a)| (four-square identity)", max(errs_norm) < 1e-11,
      f"max err {max(errs_norm):.1e}")
cvec = rng.normal(size=4)
A = 0.7 * cvec; B = -1.3 * cvec
check("A7 inputs along one line commute, so no turn: one-dimensional pages have no precession",
      np.max(np.abs(circ(A, B, 0.8) - circ(B, A, 0.8))) < 1e-14)

print("== B. One composition law, two ends")
# B1: rotations. Gibbs vectors g, h; (1+g)(1+h) = (1 - g.h) + (g + h + g x h) => g' = (g + h + g x h)/(1 - g.h)
err = 0
for _ in range(2000):
    w1 = rng.normal(size=3); w2 = rng.normal(size=3)
    g = rotor_gibbs(w1); h = rotor_gibbs(w2)
    prod = qmul(ONE + pure(g), ONE + pure(h))
    gp = prod[1:] / prod[0]
    gform = (g + h + np.cross(g, h)) / (1 - g @ h)
    Rc = rotmat(w1) @ rotmat(w2)
    th = 2 * math.atan(np.linalg.norm(gp))
    Rg = rotmat(th * gp / np.linalg.norm(gp))
    err = max(err, np.max(np.abs(gp - gform)), np.max(np.abs(Rc - Rg)))
check("B1 turns compose by o_1 on quaternions: g' = (g + h + g x h)/(1 - g.h) (Rodrigues 1840), matches SO(3) "
      "for 2000 random pairs; the dot part is the plane dial d = +1, the cross part is the tilt", err < 1e-9, f"max err {err:.1e}")

# B2: boosts in the Pauli algebra: (1 + u.sigma)(1 + v.sigma) = (1 + u.v) + (u + v).sigma + i (u x v).sigma
sx = np.array([[0, 1], [1, 0]], complex); sy = np.array([[0, -1j], [1j, 0]]); sz = np.array([[1, 0], [0, -1]], complex)
SIG = [sx, sy, sz]; I2 = np.eye(2)
def vs(v):
    return v[0] * sx + v[1] * sy + v[2] * sz
def lorentz_from_sl2(M):
    """4x4 Lorentz matrix of X -> M X M^dagger on X = t + x.sigma."""
    basis = [I2] + SIG
    L = np.zeros((4, 4))
    for j in range(4):
        Y = M @ basis[j] @ M.conj().T
        for i in range(4):
            L[i, j] = 0.5 * np.real(np.trace(basis[i] @ Y))
    return L
def boost4(beta):
    bb = np.linalg.norm(beta); gam = 1 / math.sqrt(1 - bb**2); n = beta / bb
    L = np.eye(4); L[0, 0] = gam; L[0, 1:] = L[1:, 0] = gam * beta
    L[1:, 1:] = np.eye(3) + (gam - 1) * np.outer(n, n)
    return L
errs, errw = [], []
for _ in range(500):
    b1 = rng.normal(size=3); b1 *= rng.uniform(0.05, 0.95) / np.linalg.norm(b1)
    b2 = rng.normal(size=3); b2 *= rng.uniform(0.05, 0.95) / np.linalg.norm(b2)
    # half-velocity: t = tanh(rapidity / 2) along the direction
    t1 = math.tanh(math.atanh(np.linalg.norm(b1)) / 2) * b1 / np.linalg.norm(b1)
    t2 = math.tanh(math.atanh(np.linalg.norm(b2)) / 2) * b2 / np.linalg.norm(b2)
    P = (I2 + vs(t1)) @ (I2 + vs(t2))
    form = (1 + t1 @ t2) * I2 + vs(t1 + t2) + 1j * vs(np.cross(t1, t2))
    errs.append(np.max(np.abs(P - form)))
    # Wigner rotation from the corner: tan(theta/2) = |t1 x t2| / (1 + t1.t2)
    th_corner = 2 * math.atan2(np.linalg.norm(np.cross(t1, t2)), 1 + t1 @ t2)
    Lc = boost4(b1) @ boost4(b2)
    bt = Lc[1:, 0] / Lc[0, 0]
    Rw = np.linalg.inv(boost4(bt)) @ Lc
    th_lor = math.acos(np.clip((np.trace(Rw[1:, 1:]) - 1) / 2, -1, 1))
    errw.append(abs(th_corner - th_lor))
check("B2 boosts compose by o_1 in the Pauli algebra: (1 + u)(1 + v) = (1 + u.v) + (u + v) + i(u x v) "
      "(the dial flips: d = -1 puts + in front of u.v)", max(errs) < 1e-12, f"max err {max(errs):.1e}")
check("B3 the cross part of the boost corner is a turn: Thomas-Wigner angle tan(theta/2) = |u x v|/(1 + u.v) "
      "matches 4x4 Lorentz matrices (500 random pairs)", max(errw) < 1e-9, f"max err {max(errw):.1e}")
g1, g2 = 2.0, 3.0
t1 = math.tanh(math.acosh(g1) / 2) * np.array([1, 0, 0]); t2 = math.tanh(math.acosh(g2) / 2) * np.array([0, 1, 0])
thc = 2 * math.atan2(np.linalg.norm(np.cross(t1, t2)), 1 + t1 @ t2)
check("B4 perpendicular boosts: corner angle obeys the textbook cos(theta) = (g1 + g2)/(1 + g1 g2)",
      abs(math.cos(thc) - (g1 + g2) / (1 + g1 * g2)) < 1e-12)
RESULTS['B'] = {"perp_example": {"g1": g1, "g2": g2, "theta_deg": math.degrees(thc)}}

print("== C. Twist = d x area: Foucault from the d = +1 law, Thomas from d = -1")
def oplus(u, v, d):
    return (u + v) / (1 - d * np.conj(u) * v)
def twist(u, v, d):
    cu = np.conj(u) * v
    X = 1 - d * cu.real; Y = d * cu.imag
    return 2 * math.atan2(Y, X)
def chart_radius(rho, d):
    if d > 0:
        return math.tan(math.sqrt(d) * rho / 2) / math.sqrt(d)
    if d < 0:
        return math.tanh(math.sqrt(-d) * rho / 2) / math.sqrt(-d)
    return rho / 2
def fan_holonomy(points, d):
    """Sum of twists of the fan of triangles (centre, w_j, w_{j+1}), each step written as w_{j+1} = w_j (+) v_j."""
    tot = 0.0
    for j in range(len(points)):
        w = points[j]; wn = points[(j + 1) % len(points)]
        vj = oplus(-w, wn, d)
        tot += twist(w, vj, d)
    return tot
# C1 left cancellation and single-triangle area at d = +1 (sphere) and d = -1 (hyperbolic)
def sph_point(w):
    th = 2 * math.atan(abs(w)); ph = np.angle(w)
    return np.array([math.sin(th) * math.cos(ph), math.sin(th) * math.sin(ph), math.cos(th)])
def sph_area(p, q, r):
    num = p @ np.cross(q, r)
    den = 1 + p @ q + q @ r + r @ p
    return 2 * math.atan2(num, den)
def hyp_point(w):  # Poincare disk (curvature -1) -> hyperboloid
    r2 = abs(w)**2
    return np.array([(1 + r2) / (1 - r2), 2 * w.real / (1 - r2), 2 * w.imag / (1 - r2)])
def hyp_area(p, q, r):
    # signed area via angle defect from side lengths, sign from orientation in the disk
    def dist(x, y):
        return math.acosh(max(1.0, x[0] * y[0] - x[1] * y[1] - x[2] * y[2]))
    a_, b_, c_ = dist(q, r), dist(p, r), dist(p, q)
    def ang(opp, s1, s2):
        return math.acos(np.clip((math.cosh(s1) * math.cosh(s2) - math.cosh(opp)) / (math.sinh(s1) * math.sinh(s2)), -1, 1))
    A = ang(a_, b_, c_); B = ang(b_, a_, c_); C = ang(c_, a_, b_)
    return math.pi - (A + B + C)
cancel_err, tri_err_s, tri_err_h = 0, 0, 0
for _ in range(500):
    for d in (1.0, -1.0):
        rmax = 2.0 if d > 0 else 0.9
        w = (rng.uniform(-1, 1) + 1j * rng.uniform(-1, 1)) * rmax / 1.5
        x = (rng.uniform(-1, 1) + 1j * rng.uniform(-1, 1)) * rmax / 1.5
        if d < 0 and (abs(w) >= 0.95 or abs(x) >= 0.95):
            continue
        cancel_err = max(cancel_err, abs(oplus(w, oplus(-w, x, d), d) - x))
        vv = oplus(-w, x, d)
        om = twist(w, vv, d)
        if d > 0:
            ar = sph_area(sph_point(0j), sph_point(w), sph_point(x))
            tri_err_s = max(tri_err_s, abs(om - ar))
        else:
            ar = hyp_area(hyp_point(0j), hyp_point(w), hyp_point(x))
            sign = np.sign((np.conj(w) * x).imag)
            tri_err_h = max(tri_err_h, abs(om - d * sign * ar))
check("C1 the plane law cancels on the left: w (+) ((-w) (+) x) = x, at d = +1 and d = -1", cancel_err < 1e-10,
      f"max err {cancel_err:.1e}")
check("C2 one triangle: twist = +1 x (spherical area) at d = +1 (500 random triangles; Gauss-Bonnet)", tri_err_s < 1e-10,
      f"max err {tri_err_s:.1e}")
check("C3 one triangle: twist = -1 x (hyperbolic area) at d = -1 (Thomas side, as on the twisted-boosts page)",
      tri_err_h < 1e-8, f"max err {tri_err_h:.1e}")

# C4 one curve family in d: holonomy around a circle of radius rho = 2 pi (1 - cos(sqrt(d) rho))
curve = []
fam_err = 0
for d in (1.0, 0.5, 0.25, 0.0, -0.25, -0.5, -1.0):
    rho = 1.0
    t = chart_radius(rho, d)
    pts = [t * np.exp(1j * 2 * math.pi * j / 4000) for j in range(4000)]
    hol = fan_holonomy(pts, d)
    if d > 0:
        pred = 2 * math.pi * (1 - math.cos(math.sqrt(d) * rho))
    elif d < 0:
        pred = 2 * math.pi * (1 - math.cosh(math.sqrt(-d) * rho))
    else:
        pred = 0.0
    fam_err = max(fam_err, abs(hol - pred))
    curve.append({"d": d, "rho": rho, "holonomy": hol, "pred": pred})
check("C4 one law for every dial: fan of o-twists around a circle of radius rho gives 2 pi (1 - cos(sqrt(d) rho)) "
      "(d from +1 to -1; cosh for d < 0; 0 at d = 0, Galilean: no precession)", fam_err < 1e-5, f"max err {fam_err:.1e}")
RESULTS['C_family'] = curve

# C5 Foucault: latitude lambda is the circle of colatitude rho = pi/2 - lambda on the unit sphere (d = +1).
fou = []
fe = 0
for lat in [0, 10, 20, 30, 40, 48.85, 50, 60, 70, 80, 89.9]:
    rho = math.radians(90 - lat)
    t = chart_radius(rho, 1.0)
    pts = [t * np.exp(1j * 2 * math.pi * j / 4000) for j in range(4000)]
    hol = fan_holonomy(pts, 1.0)
    per_day = 2 * math.pi - hol                      # local frame turns 2 pi per day; the swing plane lags by the holonomy
    fe = max(fe, abs(per_day - 2 * math.pi * math.sin(math.radians(lat))))
    fou.append({"lat": lat, "per_day_deg": math.degrees(per_day)})
check("C5 Foucault: swing plane turns 2 pi - (cap area) = 2 pi sin(latitude) per sidereal day, from the d = +1 twists alone",
      fe < 1e-5, f"max err {fe:.1e} rad")
paris = [f for f in fou if f["lat"] == 48.85][0]["per_day_deg"] / 23.9345
check("C6 Paris (48.85 N, Foucault 1851): 11.3 deg per hour (about 31.8 h per full turn); pole 15.04 deg/h, equator 0",
      abs(paris - 11.32) < 0.01, f"{paris:.3f} deg/h")
# C7 direct simulation of a pendulum on the rotating Earth (Coriolis), exaggerated spin for speed
def foucault_sim(lat, Om=0.02, w0=1.0, T=300.0):
    Oz = Om * math.sin(math.radians(lat))
    def f(t, s):
        x, y, vx, vy = s
        return [vx, vy, -w0**2 * x + 2 * Oz * vy, -w0**2 * y - 2 * Oz * vx]
    sol = solve_ivp(f, (0, T), [1, 0, 0, 0], rtol=1e-11, atol=1e-12, dense_output=True)
    ts = np.linspace(T - 2 * math.pi / w0, T, 400)
    X = sol.sol(ts)[:2]
    C = X @ X.T
    evals, evecs = np.linalg.eigh(C)
    ax = evecs[:, -1]
    ang = math.atan2(ax[1], ax[0])
    # axis is defined mod pi; plane rotates clockwise (negative) in the north
    pred = -Oz * (T - math.pi / w0)
    diff = (ang - pred + math.pi / 2) % math.pi - math.pi / 2
    return ang, pred, diff
sim_err = max(abs(foucault_sim(lat)[2]) for lat in (30, 48.85, 70))
check("C7 a pendulum integrated on the rotating Earth (Coriolis force) turns at the same rate (3 latitudes)",
      sim_err < 2e-3, f"max miss {sim_err:.1e} rad")
RESULTS['C_foucault'] = fou
# C8 Thomas from the same fan, d = -1: per lap 2 pi (1 - gamma), magnitude 2 pi (gamma - 1)
th_err = 0; thom = []
for beta in (0.1, 0.3, 0.5, 0.7, 0.9):
    eta = math.atanh(beta); gam = 1 / math.sqrt(1 - beta**2)
    t = chart_radius(eta, -1.0)
    pts = [t * np.exp(1j * 2 * math.pi * j / 4000) for j in range(4000)]
    hol = fan_holonomy(pts, -1.0)
    th_err = max(th_err, abs(hol + 2 * math.pi * (gam - 1)))
    thom.append({"beta": beta, "per_lap_deg": -math.degrees(hol)})
check("C8 Thomas from the same fan code with the dial flipped (d = -1): 2 pi (gamma - 1) per lap, retrograde",
      th_err < 1e-4, f"max err {th_err:.1e}")
RESULTS['C_thomas'] = thom

print("== D. Free rigid body: Montgomery's phase, orientation built only by composing turns with o_1")
I = np.array([1.0, 2.0, 3.0])
def euler_rhs(t, w):
    return [(I[1] - I[2]) / I[0] * w[1] * w[2], (I[2] - I[0]) / I[1] * w[2] * w[0], (I[0] - I[1]) / I[2] * w[0] * w[1]]
def body_period(w0, axis):
    """Period of the body-frame angular momentum: time for its angle around 'axis' to advance 2 pi."""
    def ev(t, w):
        L = I * w
        o = [i for i in range(3) if i != axis]
        return math.atan2(L[o[1]], L[o[0]]) - math.atan2(I[o[1]] * w0[o[1]], I[o[0]] * w0[o[0]])
    sol = solve_ivp(euler_rhs, (0, 200), w0, rtol=1e-13, atol=1e-14, dense_output=True, max_step=0.01)
    ts = np.linspace(0, 200, 400001)
    L = (I[:, None] * sol.sol(ts))
    o = [i for i in range(3) if i != axis]
    ang = np.unwrap(np.arctan2(L[o[1]], L[o[0]]))
    ang -= ang[0]
    target = 2 * math.pi * np.sign(ang[-1])
    idx = np.argmax(np.abs(ang) >= abs(target))
    # refine by root finding on dense output
    from scipy.optimize import brentq
    def fz(t):
        Lt = I * sol.sol(t)
        a = math.atan2(Lt[o[1]], Lt[o[0]]) - math.atan2(L[o[1], 0], L[o[0], 0])
        return (a + math.pi) % (2 * math.pi) - math.pi
    T = brentq(fz, ts[idx - 3], ts[idx + 1], xtol=1e-14)
    return T, sol
def compose_orientation(sol, T, n=40000):
    """q(t+h) = q(t) (1 + g) with g from 4th-order Magnus on body angular velocity; i.e. a_{n+1} = a_n o_1 g."""
    h = T / n
    aq = np.zeros(4)                    # q = 1 + a, start at identity
    c1, c2 = 0.5 - math.sqrt(3) / 6, 0.5 + math.sqrt(3) / 6
    for j in range(n):
        w1 = sol.sol(j * h + c1 * h); w2 = sol.sol(j * h + c2 * h)
        Om = (h / 4) * (w1 + w2) + (math.sqrt(3) / 24) * h * h * np.cross(w1, w2)   # half-angle vector
        th = np.linalg.norm(Om)
        g = math.tan(th) * Om / th
        aq = circ(aq, pure(g), 1.0)
        q = ONE + aq
        q = q / np.linalg.norm(q)      # keep 1 + a a unit rotor (scale does not change the rotation)
        aq = q - ONE
    return ONE + aq
def solid_angle(sol, T, axis, n=20000):
    """Solid angle enclosed by the body-frame L curve, by the d = +1 fan of twists about the body axis."""
    ts = np.linspace(0, T, n, endpoint=False)
    L = (I[:, None] * sol.sol(ts)).T
    L = L / np.linalg.norm(L, axis=1)[:, None]
    o = [i for i in range(3) if i != axis]
    s = np.sign(L[0, axis])
    # stereographic chart about the pole s*e_axis
    pts = (L[:, o[0]] + 1j * L[:, o[1]]) / (1 + s * L[:, axis])
    return fan_holonomy(list(pts), 1.0)
mont = []
merr = 0
for w0, axis in ((np.array([0.1, 0.2, 1.0]), 2), (np.array([0.3, 0.05, 0.9]), 2), (np.array([1.0, 0.2, 0.15]), 0), (np.array([0.9, 0.3, 0.3]), 0)):
    T, sol = body_period(w0, axis)
    q = compose_orientation(sol, T)
    L0b = I * w0; E = 0.5 * np.sum(I * w0**2); Lm = np.linalg.norm(L0b)
    Lhat_space = L0b / Lm            # q(0) = 1, so space L = body L at t = 0
    # q(T) = rotation about space L by dtheta (since body L(T) = body L(0))
    vec = q[1:]; dth = 2 * math.atan2(vec @ Lhat_space, q[0])
    Lam = solid_angle(sol, T, axis)
    pred = 2 * E * T / Lm - Lam
    diff = (dth - pred + math.pi) % (2 * math.pi) - math.pi
    alt = (dth - (2 * E * T / Lm + Lam) + math.pi) % (2 * math.pi) - math.pi
    off_axis = np.linalg.norm(vec - (vec @ Lhat_space) * Lhat_space)
    merr = max(merr, abs(diff))
    mont.append({"w0": w0.tolist(), "axis": axis, "T": T, "dtheta": dth, "twoET_over_L": 2 * E * T / Lm,
                 "Lambda": Lam, "miss": diff, "miss_other_sign": alt, "off_axis": off_axis})
    print(f"   axis {axis+1}, w0 {w0}: T = {T:.4f}, dtheta = {dth:+.6f}, 2ET/L = {2*E*T/Lm:.6f}, Lambda = {Lam:+.6f}, miss {diff:.1e} (other sign {alt:.2f})")
check("D1 after one body period the composed rotation is a pure turn about the fixed angular momentum (4 runs)",
      max(m["off_axis"] for m in mont) < 1e-7, f"max off-axis {max(m['off_axis'] for m in mont):.1e}")
check("D2 Montgomery's rigid-body phase: turn = 2ET/|L| - Lambda, with Lambda the solid angle from the d = +1 fan "
      "(4 runs, both stable axes)", merr < 1e-6, f"max miss {merr:.1e} rad")
check("D3 the other sign (+Lambda) fails, so the twist term is really there",
      min(abs(m["miss_other_sign"]) for m in mont) > 0.05)
RESULTS['D_montgomery'] = mont

print("== E. The torqued gyroscope (heavy symmetric top)")
m_, g_, l_ = 1.0, 1.0, 1.0
I1, I3 = 1.5, 1.0                      # about the pivot (I1 includes m l^2)
Ib = np.array([I1, I1, I3])
def top_omega(L, q):
    Lb = rotate(qconj(q), L)
    return rotate(q, Lb / Ib)
def top_torque(q):
    ax = rotate(q, np.array([0, 0, 1.0]))
    return np.cross(l_ * ax, np.array([0, 0, -m_ * g_]))
def run_top(L0, q0, T, h):
    """Splitting integrator: half kick of L by torque, turn q by composing (1 + g) on the left, half kick."""
    L = L0.copy(); q = q0.copy()
    n = int(round(T / h))
    axes = np.zeros((n + 1, 3)); axes[0] = rotate(q, [0, 0, 1])
    Es = []; perp = 0.0
    for j in range(n):
        tau = top_torque(q)
        ax = rotate(q, [0, 0, 1])
        perp = max(perp, abs(tau @ ax), abs(tau[2]))
        L = L + 0.5 * h * tau
        w = top_omega(L, q)
        for _ in range(2):            # midpoint turn
            g = rotor_gibbs(w * h)
            qn = qmul(ONE + pure(g), q); qn /= np.linalg.norm(qn)
            qm = qmul(ONE + pure(rotor_gibbs(w * h / 2)), q); qm /= np.linalg.norm(qm)
            w = top_omega(L, qm)
        g = rotor_gibbs(w * h)
        q = qmul(ONE + pure(g), q); q /= np.linalg.norm(q)
        L = L + 0.5 * h * top_torque(q)
        axes[j + 1] = rotate(q, [0, 0, 1])
        if j % 200 == 0:
            wv = top_omega(L, q)
            Es.append(0.5 * wv @ L + m_ * g_ * l_ * axes[j + 1][2])
    return axes, np.array(Es), perp
w3 = 20.0
Wp = m_ * g_ * l_ / (I3 * w3)
# body axis 3 = axle; start horizontal along x: q0 rotates z to x (about y by +90 deg)
q0 = np.array([math.cos(math.pi / 4), 0, math.sin(math.pi / 4), 0])
L0 = I3 * w3 * np.array([1.0, 0, 0]) + I1 * Wp * np.array([0, 0, 1.0])
Tp = 2 * math.pi / Wp
axes, Es, perp = run_top(L0, q0, Tp / 4, 2e-3)
phi = np.unwrap(np.arctan2(axes[:, 1], axes[:, 0]))
rate = (phi[-1] - phi[0]) / (Tp / 4)
check("E1 the push is at 90 deg: gravity's torque is horizontal and perpendicular to the axle at every step",
      perp < 1e-12, f"max |tau.axle|, |tau_z| = {perp:.1e}")
check("E2 steady precession, axle horizontal: rate = m g l /(I3 w3) from composing turns (quarter lap)",
      abs(rate / Wp - 1) < 2e-4, f"rate {rate:.6f} vs {Wp:.6f}")
check("E3 the axle stays level in steady precession", np.max(np.abs(axes[:, 2])) < 2e-4, f"max height {np.max(np.abs(axes[:,2])):.1e}")
check("E4 energy conserved", np.ptp(Es) / abs(Es[0]) < 1e-5, f"rel drift {np.ptp(Es)/abs(Es[0]):.1e}")
# released with no precession: nutation loops
L1 = I3 * w3 * np.array([1.0, 0, 0])
axes2, Es2, _ = run_top(L1, q0, 0.8 * Tp / 4, 1e-3)
tt = np.arange(len(axes2)) * 1e-3
phi2 = np.unwrap(np.arctan2(axes2[:, 1], axes2[:, 0]))
dip = -np.min(axes2[:, 2])
pred_dip = 2 * I1 * m_ * g_ * l_ / (I3 * w3)**2    # fast-top nutation depth (sin theta = 1)
mean_rate = np.polyfit(tt, phi2, 1)[0]
check("E5 released without precession it dips and loops (nutation); depth matches 2 I1 m g l/(I3 w3)^2 (fast top)",
      abs(math.asin(dip) / pred_dip - 1) < 0.03, f"dip {math.asin(dip):.5f} rad vs {pred_dip:.5f}")
check("E6 its average precession rate is still m g l /(I3 w3) to within 3%", abs(mean_rate / Wp - 1) < 0.03,
      f"{mean_rate:.5f} vs {Wp:.5f}")
nut_f = I3 * w3 / I1
RESULTS['E'] = {"Wp": Wp, "w3": w3, "I1": I1, "I3": I3, "steady_rate": rate, "dip": float(math.asin(dip)), "pred_dip": pred_dip,
                "mean_rate": mean_rate, "nutation_freq": nut_f,
                "trace_x": axes2[::40, 0].tolist(), "trace_y": axes2[::40, 1].tolist(), "trace_z": axes2[::40, 2].tolist(),
                "steady_x": axes[::200, 0].tolist(), "steady_y": axes[::200, 1].tolist()}

print("== F. Larmor precession and Berry's half angle (spin 1/2)")
# F1: rotor steps composed with o_1; Bloch vector precesses at gamma B about B (dS/dt = gamma S x B)
gam_, B0 = 2.0, 1.5
Bhat = np.array([0.3, -0.2, 1.0]); Bhat /= np.linalg.norm(Bhat)
S0 = np.array([1.0, 0, 0])
h = 1e-3; n = 3000
aq = np.zeros(4)
gstep = pure(rotor_gibbs(-gam_ * B0 * h * Bhat))
for _ in range(n):
    aq = circ(aq, gstep, 1.0)
qL = (ONE + aq) / np.linalg.norm(ONE + aq)
S_comp = rotate(qL, S0)
S_exact = rotmat(-gam_ * B0 * n * h * Bhat) @ S0
check("F1 Larmor: spin turned by composing small turns with o_1 matches rotation at gamma B about B",
      np.max(np.abs(S_comp - S_exact)) < 1e-10, f"err {np.max(np.abs(S_comp - S_exact)):.1e}")
check("F2 precession keeps |S| and the angle to B (the turn is at 90 deg to both)",
      abs(np.linalg.norm(S_comp) - 1) < 1e-12 and abs(S_comp @ Bhat - S0 @ Bhat) < 1e-10)
# F3: Berry phase for a slowly rotating field around a cone, spinor steps (1 + i tan(b h) n.sigma) composed
def berry_run(theta0, T, b=1.0, nsteps=None):
    nsteps = nsteps or int(T / 0.02)
    hh = T / nsteps; nu = 2 * math.pi / T
    def nvec(t):
        return np.array([math.sin(theta0) * math.cos(nu * t), math.sin(theta0) * math.sin(nu * t), math.cos(theta0)])
    n0 = nvec(0.0)
    # eigenvector of n0.sigma with eigenvalue +1
    evals, evecs = np.linalg.eigh(vs(n0)); psi0 = evecs[:, 1]
    U = I2.astype(complex)
    for j in range(nsteps):
        nn = nvec((j + 0.5) * hh)
        step = I2 + 1j * math.tan(b * hh) * vs(nn)        # 1 + g form; scalar cos(b h) dropped, restored below
        U = step @ U
    U = U * math.cos(b * hh)**nsteps
    psiT = U @ psi0
    ov = np.vdot(psi0, psiT)
    total = np.angle(ov)
    dyn = b * T                                            # H = -b n.sigma, E = -b
    berry = (total - dyn + math.pi) % (2 * math.pi) - math.pi
    return berry, abs(ov)
berr = []; bres = []
for th0 in (0.4, 0.9, 1.3):
    br, mag = berry_run(th0, 3000.0)
    pred = -math.pi * (1 - math.cos(th0))
    berr.append(abs(br - pred)); bres.append({"theta0": th0, "berry": br, "pred": pred, "overlap": mag})
check("F3 Berry: spin 1/2 following a field around a cone picks up -1/2 x (solid angle), cone angles 0.4, 0.9, 1.3 rad",
      max(berr) < 3e-3, f"max miss {max(berr):.1e} (adiabatic error ~ nu/b)")
hol = fan_holonomy([chart_radius(0.9, 1.0) * np.exp(1j * 2 * math.pi * j / 4000) for j in range(4000)], 1.0)
check("F4 the Berry phase is exactly half the d = +1 twist of the field's loop: the spinor is the half-angle chart",
      abs(-0.5 * hol - bres[1]["pred"]) < 1e-6)
RESULTS['F'] = bres

print("== G. Frame dragging in the spinning BTZ black hole (2+1 gravity, l = 1, 8G = 1)")
t, r, ph = sp.symbols('t r phi', real=True)
Mb, Jb = sp.symbols('M J', positive=True)
X3 = [t, r, ph]
f = r**2 - Mb + Jb**2 / (4 * r**2)
Np = -Jb / (2 * r**2)
g3 = sp.Matrix([[-f + r**2 * Np**2, 0, r**2 * Np], [0, 1 / f, 0], [r**2 * Np, 0, r**2]])
g3 = sp.simplify(g3)
check("G1 rotating BTZ: g_tt = -(r^2 - M) and g_t phi = -J/2 exactly (static observers exist outside r^2 = M)",
      sp.simplify(g3[0, 0] + (r**2 - Mb)) == 0 and sp.simplify(g3[0, 2] + Jb / 2) == 0)
detg = sp.simplify(g3.det())
def vorticity(Om):
    xi = sp.Matrix([1, 0, Om])
    nrm = sp.simplify((xi.T * g3 * xi)[0])
    ul = (g3 * xi) / sp.sqrt(-nrm)
    w = 0
    for m in range(3):
        for nn in range(3):
            for p in range(3):
                e = sp.LeviCivita(m, nn, p)
                if e != 0:
                    w += e * ul[m] * sp.diff(ul[p], X3[nn])
    return sp.simplify(-w / (2 * sp.sqrt(-detg)))
w_static = sp.simplify(vorticity(0).subs(sp.Abs(r), r))
w_zamo = vorticity(Jb / (2 * r**2))
check("G2 a gyroscope held still beside the hole turns, relative to the static frame locked to the far boundary, at "
      "w = J / (2 (r^2 - M)) per unit proper time; a zero-angular-momentum observer's gyroscope does not",
      sp.simplify(w_static - Jb / (2 * (r**2 - Mb))) == 0 and sp.simplify(w_zamo) == 0, f"w_static = {w_static}")
# G3 independent numeric Fermi-Walker transport
Gam = [[[sp.simplify(sum(g3.inv()[a_, d_] * (sp.diff(g3[d_, b_], X3[c_]) + sp.diff(g3[d_, c_], X3[b_]) - sp.diff(g3[b_, c_], X3[d_]))
                         for d_ in range(3)) / 2) for c_ in range(3)] for b_ in range(3)] for a_ in range(3)]
Mv, Jv, r0 = 1.0, 0.8, 1.6
Gn = np.array([[[float(Gam[a_][b_][c_].subs({Mb: Mv, Jb: Jv, r: r0})) for c_ in range(3)] for b_ in range(3)] for a_ in range(3)])
gn = np.array(g3.subs({Mb: Mv, Jb: Jv, r: r0}).evalf(), dtype=float)
ut = 1 / math.sqrt(r0**2 - Mv); uvec = np.array([ut, 0, 0])
acc = np.einsum('abc,b,c->a', Gn, uvec, uvec)
def fw(tau, S):
    aS = (gn @ acc) @ S
    return -np.einsum('abc,b,c->a', Gn, uvec, S) + uvec * aS
e1 = np.array([0, math.sqrt(r0**2 - Mv + Jv**2 / (4 * r0**2)), 0])
dphi = np.array([0, 0, 1.0]); P = dphi + (uvec @ gn @ dphi) * uvec; e2 = P / math.sqrt(P @ gn @ P)
sol = solve_ivp(fw, (0, 5), e1, rtol=1e-12, atol=1e-13, dense_output=True)
taus = np.linspace(0, 5, 201)
angs = np.unwrap([math.atan2(sol.sol(s) @ gn @ e2, sol.sol(s) @ gn @ e1) for s in taus])
rate_num = np.polyfit(taus, angs, 1)[0]
w_pred = Jv / (2 * (r0**2 - Mv))
check("G3 independent check: Fermi-Walker transport integrated along a static worldline turns at the same rate",
      abs(abs(rate_num) - w_pred) < 1e-7, f"numeric {rate_num:+.9f} vs {w_pred:.9f}")
check("G4 no spin, no drag: J = 0 gives w = 0", sp.simplify(w_static.subs(Jb, 0)) == 0)
check("G5 far out the static gyroscope's rate approaches the frame-dragging angular velocity J/(2 r^2)",
      sp.limit(w_static / (Jb / (2 * r**2)), r, sp.oo) == 1)
rp, rm = sp.symbols('r_p r_m', positive=True)
TL = (rp - rm) / (2 * sp.pi); TR = (rp + rm) / (2 * sp.pi)
OmH = rm / rp
TH = (rp**2 - rm**2) / (2 * sp.pi * rp)
check("G6 horizon drag speed Omega_H l = r_-/r_+ = (T_R - T_L)/(T_R + T_L); Hawking temperature = harmonic mean of T_L, T_R",
      sp.simplify(OmH - (TR - TL) / (TR + TL)) == 0 and sp.simplify(TH - 2 / (1 / TL + 1 / TR)) == 0)
uu = OmH; rho_ = TR / TL
check("G7 its relativity-page coordinate 2u/(1 - u) (Doppler factor squared minus one) is T_R/T_L - 1: "
      "the temperature ratio is the lift", sp.simplify(2 * uu / (1 - uu) - (rho_ - 1)) == 0)
st = 0
for _ in range(200):
    r1, r2 = rng.uniform(0.1, 10, 2)
    u1, u2 = (r1 - 1) / (r1 + 1), (r2 - 1) / (r2 + 1)
    st = max(st, abs((u1 + u2) / (1 + u1 * u2) - (r1 * r2 - 1) / (r1 * r2 + 1)))
    x1, x2 = r1 - 1, r2 - 1
    st = max(st, abs((x1 + x2 + x1 * x2) - (r1 * r2 - 1)))
check("G8 stacking two drags (Einstein addition) multiplies the temperature ratios, i.e. o_1 on T_R/T_L - 1", st < 1e-12)
check("G9 the spinning bath of the probability-entanglement thread (left movers 4x hotter): horizon drags at 3/5 of light speed",
      abs((0.25 - 1) / (0.25 + 1) + 0.6) < 1e-15)
rs = np.linspace(1.02, 4, 200)
RESULTS['G'] = {"M": 1.0, "r": rs.tolist(), "J": [0.3, 0.6, 0.9],
                "w_static": {str(J): (J / (2 * (rs**2 - 1))).tolist() for J in (0.3, 0.6, 0.9)},
                "zamo": {str(J): (J / (2 * rs**2)).tolist() for J in (0.3, 0.6, 0.9)},
                "fw_check": {"r0": r0, "J": Jv, "numeric": rate_num, "pred": w_pred}}

print()
print(f"{sum(PASS)}/{len(PASS)} checks pass")
json.dump(RESULTS, open(__file__.replace('checks.py', 'results.json'), 'w'), indent=1, default=float)
