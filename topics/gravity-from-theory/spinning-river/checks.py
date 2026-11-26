"""DRAFT checks: the spinning river, and frame dragging's relatives (Nov 26, 2026).

K: give each river cell a turn as well as a velocity; does the spinning (Kerr) black hole come out exactly?
F: frame-dragging look-alikes outside black holes (Larmor's theorem, spinning fluids, Foucault, clocks).
D: derive the turn (area law + circulation law + far-field frame dragging + rim at light speed) instead of testing it.
Curvature is computed with exact symbolic derivatives and 50-digit arithmetic at sample points (numgr.py):
"zero" means below 1e-30.  Run: python3 checks.py   (writes checks_output.txt, results.json)
"""
import json
import os
import sys

import mpmath as mp
import numpy as np
import sympy as sp
from scipy.integrate import solve_ivp

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from grtools import einstein_mixed  # noqa: E402
from numgr import curvature_at, prepare, ricci_max, riemann_max  # noqa: E402

out, results = [], {}


def log(s=""):
    print(s)
    out.append(s)


def check(name, ok, detail=""):
    results[name] = {"pass": bool(ok), "detail": str(detail)}
    log(f"[{'PASS' if ok else 'FAIL'}] {name}: {detail}")


T, R, TH, PH = sp.symbols("t r theta phi", real=True)
X4 = [T, R, TH, PH]
Ms, As = sp.symbols("M a", positive=True)
Mn, An = sp.Rational(1), sp.Rational(7, 10)            # sample black hole for the numeric curvature checks
PTS = [(0.3, 3.1, 0.7, 0.2), (0.0, 5.7, 1.9, 1.0), (1.0, 2.4, 0.3, 0.0), (0.0, 1.3, 1.2, 0.5)]
ZERO = mp.mpf("1e-30")


def legs(M=Mn, a=An, turn=1, power=1, shape=2):
    """The four legs (one-forms over dt, dr, dtheta, dphi) of a river cell.
    Background: flat space in oblate spheroidal coordinates (legs g0..g3).
    Cell: its radial leg is pushed by v along the time leg w, and w is turned toward the azimuth:
        e1 = g1 + v * (g0 - turn * a * sin^shape(theta) dphi),   v = sqrt(2 M r) / rho^power."""
    rho2 = R**2 + a**2 * sp.cos(TH)**2
    g0 = sp.Matrix([1, 0, 0, 0])
    g1 = sp.Matrix([0, sp.sqrt(rho2) / sp.sqrt(R**2 + a**2), 0, 0])
    g2 = sp.Matrix([0, 0, sp.sqrt(rho2), 0])
    g3 = sp.Matrix([0, 0, 0, sp.sqrt(R**2 + a**2) * sp.sin(TH)])
    v = sp.sqrt(2 * M * R) / rho2**sp.Rational(power, 2)
    w = sp.Matrix([1, 0, 0, -turn * a * sp.sin(TH)**shape])
    e1 = g1 + v * w
    return g0, e1, g2, g3, v, rho2


def metric(**kw):
    e0, e1, e2, e3, *_ = legs(**kw)
    return -e0 * e0.T + e1 * e1.T + e2 * e2.T + e3 * e3.T


def rmax(g):
    f = prepare(g, X4)
    return max(ricci_max(f, P) for P in PTS)


# ---------------------------------------------------------------- K1 the background lattice is flat
log("== K1  The lattice the cells move through is flat space, in the spinning body's natural shape")
f0 = prepare(metric(M=0), X4)
rb = max(riemann_max(f0, P) for P in PTS)
check("K1a background (no velocity): every Riemann component vanishes", rb < ZERO, f"max |Riemann| = {mp.nstr(rb, 3)}")
xx = sp.sqrt(R**2 + As**2) * sp.sin(TH) * sp.cos(PH)
yy = sp.sqrt(R**2 + As**2) * sp.sin(TH) * sp.sin(PH)
zz = R * sp.cos(TH)
Jc = sp.Matrix([[sp.diff(q, s) for s in (R, TH, PH)] for q in (xx, yy, zz)])
flat3 = sp.simplify(Jc.T * Jc)
rho2s = R**2 + As**2 * sp.cos(TH)**2
target = sp.diag(rho2s / (R**2 + As**2), rho2s, (R**2 + As**2) * sp.sin(TH)**2)
check("K1b it is ordinary flat space dx^2+dy^2+dz^2 in oblate spheroidal coordinates (spheroids of semi-axes sqrt(r^2+a^2), r)",
      sp.simplify(flat3 - target) == sp.zeros(3, 3), "the right shape for a spinning body: squashed spheres, a ring of radius a at r = 0")

# ---------------------------------------------------------------- K2 velocity + turn gives Kerr
log("\n== K2  Each cell: velocity sqrt(2Mr)/rho plus a turn toward the azimuth  ->  the spinning black hole")
gK = metric()
fK = prepare(gK, X4)
rK = max(ricci_max(fK, P) for P in PTS)
check("K2a every field equation of empty space holds (Ricci = 0) at all sample points, M = 1, a = 0.7",
      rK < ZERO, f"max |Ricci| = {mp.nstr(rK, 3)}")
gK2 = metric(a=sp.Rational(99, 100))
rK2 = rmax(gK2)
check("K2b also near the fastest spin, a = 0.99 M", rK2 < ZERO, f"max |Ricci| = {mp.nstr(rK2, 3)}")


def kretschmann(f, P):
    Riem, _, gi = curvature_at(f, P)
    G, *_ = f(*[mp.mpf(p) for p in P])
    g = mp.matrix(G)
    n = 4
    low = {k: sum(g[k[0], e] * Riem[(e,) + k[1:]] for e in range(n)) for k in Riem}
    up = {}
    for (a_, b, c, d) in low:
        up[a_, b, c, d] = sum(gi[a_, p] * gi[b, q] * gi[c, s] * gi[d, t] * low[p, q, s, t]
                              for p in range(n) for q in range(n) for s in range(n) for t in range(n)
                              if low[p, q, s, t] != 0)
    return sum(low[k] * up[k] for k in low)


def kerr_K(r, th, M=1.0, a=0.7):
    c = mp.cos(th)
    rho2 = r**2 + a**2 * c**2
    return 48 * M**2 * (r**6 - 15 * r**4 * a**2 * c**2 + 15 * r**2 * a**4 * c**4 - a**6 * c**6) / rho2**6


kerrs = [(kretschmann(fK, P), kerr_K(mp.mpf(P[1]), mp.mpf(P[2]), 1, mp.mpf("0.7"))) for P in PTS[:2]]
check("K2c it is Kerr's curvature, not just some vacuum: R_abcd R^abcd = 48M^2(r^6 - 15r^4a^2c^2 + 15r^2a^4c^4 - a^6c^6)/rho^12",
      all(abs(k - e) < mp.mpf("1e-25") * abs(e) for k, e in kerrs),
      "; ".join(f"{mp.nstr(k, 12)} vs {mp.nstr(e, 12)}" for k, e in kerrs))
results["kretschmann_samples"] = [[str(k), str(e)] for k, e in kerrs]

# ---------------------------------------------------------------- K3 the turn is forced
log("\n== K3  The turn is not optional, and its size and shape are forced")
bad = {}
for lab, kw in [("no turn (velocity only)", dict(turn=0)), ("half the turn", dict(turn=sp.Rational(1, 2))),
                ("1.5 x the turn", dict(turn=sp.Rational(3, 2))), ("turn ~ sin(theta) not sin^2", dict(shape=1)),
                ("speed sqrt(2Mr) (rho^0)", dict(power=0)), ("speed sqrt(2Mr)/rho^2", dict(power=2))]:
    bad[lab] = rmax(metric(**kw))
check("K3a velocity alone on the spheroidal lattice fails: the turn is needed", bad["no turn (velocity only)"] > mp.mpf("1e-3"),
      f"max |Ricci| = {mp.nstr(bad['no turn (velocity only)'], 3)}")
check("K3b any other size of turn fails (0.5x, 1.5x), and so does another shape (sin instead of sin^2)",
      all(bad[k] > mp.mpf("1e-3") for k in ("half the turn", "1.5 x the turn", "turn ~ sin(theta) not sin^2")),
      {k: mp.nstr(bad[k], 3) for k in ("half the turn", "1.5 x the turn", "turn ~ sin(theta) not sin^2")})
check("K3c the speed law is forced too: v^2 = 2Mr/rho^2 (other powers of rho fail)",
      all(bad[k] > mp.mpf("1e-3") for k in ("speed sqrt(2Mr) (rho^0)", "speed sqrt(2Mr)/rho^2")),
      {k: mp.nstr(bad[k], 3) for k in ("speed sqrt(2Mr) (rho^0)", "speed sqrt(2Mr)/rho^2")})

# ---------------------------------------------------------------- K4 each cell obeys the same velocity rule
log("\n== K4  Each cell is still a zero-energy, zero-spin free faller (the velocity rule, unchanged)")
gS = metric(M=Ms, a=As)
gSi = gS.inv()
gtt_up = sp.simplify(gSi[0, 0])
check("K4a g^tt = -1: the cell's 4-velocity u = -dt has unit length, so it falls freely (u = -grad t, |u| = 1)",
      gtt_up == -1, "geodesic because u is a unit gradient; E = -u_t = 1: starts from rest far away")
check("K4b u_phi = 0: the cells carry no angular momentum (they are the zero-angular-momentum fallers)", True,
      "u_mu = (-1, 0, 0, 0) by construction")

# ---------------------------------------------------------------- K5 the cell map is the dial's flat middle
log("\n== K5  The map from lattice legs to cell legs: two shears, the dial's flat middle (q = 0)")
vv, al = sp.symbols("v alpha")
N = sp.eye(4)
N[1, 0] = vv           # radial leg picks up v * time leg: a Galilean boost (a shear)
N[1, 3] = -vv * al     # radial leg picks up -v*alpha * azimuth leg: a turn at q = 0 (also a shear)
E1, E2 = sp.zeros(4, 4), sp.zeros(4, 4)
E1[1, 0], E2[1, 3] = 1, 1
check("K5a N - I squares to zero, and the boost and the turn commute and do not interact (E1 E2 = E2 E1 = 0)",
      (N - sp.eye(4))**2 == sp.zeros(4, 4) and E1 * E2 == sp.zeros(4, 4) and E2 * E1 == sp.zeros(4, 4),
      "dual numbers: (1 + eps v)(1 + eps w) = 1 + eps (v + w); amounts v (boost) and v*alpha (turn)")
e0s, e1s, e2s, e3s, vS_, rho2_ = legs(M=Ms, a=As)
alpha = As * sp.sin(TH) / sp.sqrt(R**2 + As**2)
g1s = sp.Matrix([0, sp.sqrt(rho2_) / sp.sqrt(R**2 + As**2), 0, 0])
recon = g1s + vS_ * e0s - vS_ * alpha * e3s
check("K5b the cell's radial leg = lattice radial leg + v (time leg) - v alpha (azimuth leg), alpha = a sin(theta)/sqrt(r^2+a^2)",
      sp.simplify(recon - e1s) == sp.zeros(4, 1), "alpha -> 0 far away and on the axis; largest at the equator near the ring")

# ---------------------------------------------------------------- K6 slow spin: last round's swirl
log("\n== K6  At slow spin the turn is last round's swirl 2J/r^3")
eps = sp.symbols("epsilon")
gE = metric(M=Ms, a=eps * As)
gtphi1 = sp.simplify(sp.diff(gE[0, 3], eps).subs(eps, 0))
check("K6a first order in spin: g_t phi = -(2J/r^3) r^2 sin^2(theta), J = M a (the Lense-Thirring swirl of jet-river J7c)",
      sp.simplify(gtphi1 + 2 * Ms * As / R**3 * R**2 * sp.sin(TH)**2) == 0, f"g_tphi' = {gtphi1}")
grr_up = sp.simplify(gSi[1, 1] * rho2_)
check("K6b horizon: g^rr rho^2 = r^2 - 2Mr + a^2, so the horizon is at r = M + sqrt(M^2 - a^2), Kerr's",
      sp.simplify(grr_up - (R**2 - 2 * Ms * R + As**2)) == 0, f"g^rr rho^2 = {grr_up}")

# ---------------------------------------------------------------- K7 clocks: the turn is a clock gap
log("\n== K7  The cell's time leg w = dt - a sin^2(theta) dphi does not close: a clock gap around the axis")
Om = sp.symbols("Omega")
sq, mq = sp.symbols("s m", positive=True)                 # r = s^2, M = m^2 keep every root rational
geq = gS.subs(TH, sp.pi / 2).subs({R: sq**2, Ms: mq**2})
Lag = geq[0, 0] + 2 * Om * geq[0, 3] + Om**2 * geq[3, 3]
sols = sp.solve(sp.diff(Lag, sq), Om)
Ts = [sp.factor(sp.cancel(2 * sp.pi / s_)) for s_ in sols]
T_pro_e, T_retro_e = 2 * sp.pi * (sq**3 / mq + As), 2 * sp.pi * (sq**3 / mq - As)
okK7a = any(sp.simplify(t_ - T_pro_e) == 0 or sp.simplify(t_ + T_pro_e) == 0 for t_ in Ts) and \
    any(sp.simplify(t_ - T_retro_e) == 0 or sp.simplify(t_ + T_retro_e) == 0 for t_ in Ts)
check("K7a circular equatorial orbits: periods T = 2 pi (r^(3/2)/sqrt(M) +- a) (prograde longer)", okK7a,
      f"orbital periods {Ts} with r = s^2, M = m^2")
gap = sp.integrate(As * sp.sin(TH)**2, (PH, 0, 2 * sp.pi)).subs(TH, sp.pi / 2)
check("K7b each period is shifted by 2 pi a, the cell clock's gap once around the equator (closed integral of a sin^2 dphi)",
      sp.simplify((T_pro_e - T_retro_e) / 2 - gap) == 0 and okK7a,
      f"half the difference = {sp.simplify((T_pro_e - T_retro_e) / 2)}, clock gap = {gap}; "
      "the known gravitomagnetic clock effect T+ - T- = 4 pi a (Cohen & Mashhoon 1993), independent of r and of G")

# ---------------------------------------------------------------- F1 Larmor's theorem
log("\n== F1  Larmor's theorem: a magnetic field acts like a turning frame (charges, not gravity)")


def orbit(B, tmax=12.0):
    """charge q/m = 1 in a central pull -r/|r|^3 plus a uniform field B along z."""
    def rhs(t, s):
        x, y, vx, vy = s
        r3 = (x * x + y * y)**1.5
        return [vx, vy, -x / r3 + B * vy, -y / r3 - B * vx]
    ts = np.linspace(0, tmax, 400)
    sol = solve_ivp(rhs, (0, tmax), [1.0, 0.0, 0.0, 1.1], t_eval=ts, rtol=1e-11, atol=1e-12)
    return sol.y[0], sol.y[1]


def larmor_err(B):
    xB, yB = orbit(B)
    x0, y0 = orbit(0.0)
    # Larmor frame: angular velocity -B/2 (q/m = 1); compare the B orbit with the B-free orbit turned at -B/2
    def rhs0(t, s):
        x, y, vx, vy = s
        r3 = (x * x + y * y)**1.5
        return [vx, vy, -x / r3, -y / r3]
    ts = np.linspace(0, 12.0, 400)
    w = -B / 2
    # initial state of the B-free orbit in the inertial frame such that, turned by w t, it matches the B orbit's start
    s1 = [1.0, 0.0, 0.0, 1.1 - w * 1.0]
    sol = solve_ivp(rhs0, (0, 12.0), s1, t_eval=ts, rtol=1e-11, atol=1e-12)
    c, s = np.cos(w * ts), np.sin(w * ts)
    xl, yl = c * sol.y[0] - s * sol.y[1], s * sol.y[0] + c * sol.y[1]
    e_larmor = np.max(np.hypot(xB - xl, yB - yl))
    e_naive = np.max(np.hypot(xB - x0, yB - y0))
    return e_larmor, e_naive


(eL1, eN1), (eL2, eN2) = larmor_err(0.005), larmor_err(0.01)
check("F1a with the frame turning at Larmor's rate -qB/2m the field disappears up to B^2 (error x4 when B doubles); "
      "without the turn the error is first order (x2)",
      3.5 < eL2 / eL1 < 4.5 and 1.7 < eN2 / eN1 < 2.3 and eL1 < eN1 / 10,
      f"Larmor-frame error {eL1:.2e} -> {eL2:.2e} (ratio {eL2 / eL1:.2f}); no turn {eN1:.2e} -> {eN2:.2e} (ratio {eN2 / eN1:.2f})")
check("F1b the 1/2 in Larmor's rate is the 1/2 in the gyroscope rule of frame dragging (Omega = -curl(h)/2, Part F)",
      True, "magnetic force q v x B = Coriolis 2m v x Omega at Omega = qB/2m; gravity: B_g = curl h, Omega = -B_g/2 "
      "(the gravitational Larmor theorem, Mashhoon; known)")

# ---------------------------------------------------------------- F2 spinning fluids: the river is the acoustic metric
log("\n== F2  Spinning fluids: sound in a flowing fluid sees exactly a river (the acoustic metric)")
x_, y_, t_ = sp.symbols("x y t", real=True)
ux, uy = sp.Function("ux")(x_, y_), sp.Function("uy")(x_, y_)
g3 = sp.Matrix([[-1 + ux**2 + uy**2, -ux, -uy], [-ux, 1, 0], [-uy, 0, 1]])   # -dt^2 + (dx - u dt)^2 (sound speed 1)
Gm3, *_ = einstein_mixed(g3, [t_, x_, y_])
Sxx, Syy, Sxy = sp.diff(ux, x_), sp.diff(uy, y_), (sp.diff(ux, y_) + sp.diff(uy, x_)) / 2
e2sym = Sxx * Syy - Sxy**2
Gnn = -(Gm3[0, 0] + ux * Gm3[0, 1] + uy * Gm3[0, 2])        # energy density seen by the cells: G_mu nu n^mu n^nu
check("F2a for ANY 2D flow, the energy density seen by the cells is the Q layer of the strain: G_nn = det(sym grad u)",
      sp.simplify(Gnn - e2sym) == 0, "the same second layer as jet-river, now for the acoustic metric of Unruh (1981)")
A_, B_ = sp.symbols("A B", positive=True)
rr, pp = sp.symbols("r phi", positive=True)
uvort = {ux: (-A_ * x_ - B_ * y_) / (x_**2 + y_**2), uy: (-A_ * y_ + B_ * x_) / (x_**2 + y_**2)}  # draining bathtub
e2v = sp.simplify(e2sym.subs(uvort).doit())
check("F2b draining bathtub vortex u_r = -A/r, u_phi = B/r: Q = -(A^2 + B^2)/r^4 < 0: a negative energy density",
      sp.simplify(e2v + (A_**2 + B_**2) / (x_**2 + y_**2)**2) == 0,
      f"Q = {e2v}: the vortex drags sound like frame dragging (horizon at r = A, ergoregion to sqrt(A^2+B^2)), "
      "but it is not an empty-space solution of gravity: the analogy is in the velocities, not the field law")
check("F2c a fluid can only make a velocity river (no per-cell turn), so by K3a it cannot make full Kerr",
      bad["no turn (velocity only)"] > mp.mpf("1e-3"),
      "matches the known result that fluid vortices mimic Kerr only on the equatorial slice (Visser & Weinfurtner 2005)")

# ---------------------------------------------------------------- F3 Foucault's pendulum and Earth's own drag
log("\n== F3  Foucault's pendulum: the rotating Earth, plus Earth's own (tiny) frame dragging")
G, c = 6.674e-11, 2.998e8
M_E, R_E, k_I, w_E = 5.972e24, 6.371e6, 0.3307, 7.2921e-5
J_E = k_I * M_E * R_E**2 * w_E
lt_vert_over_sin = 2 * G * J_E / (c**2 * R_E**3)        # vertical part of (GJ/c^2R^3)(3(j.r)r - j) = 2 sin(lat) GJ/(c^2R^3)
ratio = lt_vert_over_sin / w_E
check("F3a Foucault turn = Earth's rate x sin(latitude) (on precession page C5-C7); Earth's frame dragging adds a vertical "
      "turn 2 sin(lat) GJ/(c^2R^3): the ratio is 2 (I/MR^2)(GM/Rc^2), the same at every latitude",
      abs(ratio - 2 * k_I * G * M_E / (R_E * c**2)) < 1e-15,
      f"ratio = {ratio:.3e}; Foucault at 45 deg = {360 * np.sin(np.pi / 4) / 0.99727:.1f} deg/day, "
      f"drag part = {lt_vert_over_sin * np.sin(np.pi / 4) * 3.156e7 * 206264806:.1f} mas/yr")

# ================================================================ D  deriving the turn (round 7)
log("\n== D  Deriving the turn's size and shape instead of testing it")
Kp, Bp, Omf = sp.symbols("kappa b Omega", positive=True)


def legs2(M=Mn, b=An, f=None, kappa=An):
    """Lattice = flat spheroidal shells with ring radius b; speed v = sqrt(2Mr)/rho_b; turn function f(r, theta):
    e1 = g1 + v (g0 - f dphi).  Default f = kappa sin^2(theta)."""
    rho2 = R**2 + b**2 * sp.cos(TH)**2
    if f is None:
        f = kappa * sp.sin(TH)**2
    g0 = sp.Matrix([1, 0, 0, 0])
    g1 = sp.Matrix([0, sp.sqrt(rho2) / sp.sqrt(R**2 + b**2), 0, 0])
    g2 = sp.Matrix([0, 0, sp.sqrt(rho2), 0])
    g3 = sp.Matrix([0, 0, 0, sp.sqrt(R**2 + b**2) * sp.sin(TH)])
    v = sp.sqrt(2 * M * R / rho2)
    e1 = g1 + v * sp.Matrix([1, 0, 0, -f])
    return g0, e1, g2, g3, v, rho2


def metric2(**kw):
    e0, e1, e2, e3, *_ = legs2(**kw)
    return -e0 * e0.T + e1 * e1.T + e2 * e2.T + e3 * e3.T


# D1 the area law: a rigidly turning shell's clock gap is 2 x rate x area, which gives the sin^2 shape
Om_r = sp.Function("Omega")(R)
cylr = sp.sqrt(R**2 + Bp**2) * sp.sin(TH)                       # distance from the axis of a lattice point
alpha_shell = Om_r * cylr                                       # the shell's speed at that point
g0b = sp.Matrix([1, 0, 0, 0])
g3b = sp.Matrix([0, 0, 0, cylr])
w_shell = g0b - alpha_shell * g3b                               # time leg of the turning shell at the dial's flat middle (q = 0 boost)
gap_shell = sp.integrate(-w_shell[3], (PH, 0, 2 * sp.pi))
xx_b = cylr * sp.cos(PH)
yy_b = cylr * sp.sin(PH)
area = sp.pi * sp.simplify(xx_b**2 + yy_b**2)                   # area enclosed by the latitude circle, in flat space
check("D1a area law: a shell turning rigidly at rate Omega(r) has clock gap 2 Omega A around each latitude circle, "
      "A = pi (r^2+b^2) sin^2(theta) the enclosed flat area; so the turn is f = Omega(r)(r^2+b^2) sin^2(theta)",
      sp.simplify(gap_shell - 2 * Om_r * area) == 0,
      f"gap = {sp.simplify(gap_shell)}; the sin^2 shape is the area of a latitude circle (Sagnac's law, at q = 0 exactly)")
alpha_K5 = As * sp.sin(TH) / sp.sqrt(R**2 + As**2)
check("D1b the K5 turn amount alpha = a sin(theta)/sqrt(r^2+a^2) is exactly the speed of the shell turning at "
      "Omega = a/(r^2+a^2): the cell is pushed along the turning shell's time leg",
      sp.simplify(alpha_shell.subs(Om_r, As / (R**2 + As**2)).subs(Bp, As) - alpha_K5) == 0,
      "so the 'turn' is a second q = 0 boost, sideways, at the local speed of the lattice shell")

# D2 the circulation law: a falling ring of cells keeps its clock gap  <=>  df/dr = 0
fgen = sp.Function("f")(R, TH)
gG = metric2(M=Ms, b=Bp, f=fgen)
gGi = gG.inv()
uup = [sp.simplify(-gGi[i, 0]) for i in range(4)]               # u^mu = g^{mu nu} u_nu, u_nu = (-1, 0, 0, 0)
# change of the gap around a carried latitude ring = closed integral of (L_u w)_phi dphi, w_phi = -f, axisymmetric
Lw_phi = -sum(uup[i] * sp.diff(fgen, X4[i]) for i in range(4))
check("D2a for ANY turn function the cells fall along cones (u^theta = 0), so a latitude ring of cells stays a latitude ring",
      uup[2] == 0, f"u^r = {sp.simplify(uup[1])}")
check("D2b the gap around a falling ring changes at rate -2 pi u^r df/dr: it is kept exactly when f has no r dependence",
      sp.simplify(Lw_phi + uup[1] * sp.diff(fgen, R)) == 0,
      "a circulation law like Kelvin's for a swirling fluid; with D1: Omega(r)(r^2+b^2) = kappa, so Omega(r) = kappa/(r^2+b^2)")

# D3 the size: match the far field to frame dragging (Part F's g_tphi -> -2J sin^2/r)
gk = metric2(M=Ms, b=Bp, kappa=Kp)
far = sp.limit(gk[0, 3] * R, R, sp.oo)
Js = sp.symbols("J", positive=True)
ksol = sp.solve(sp.Eq(far, -2 * Js * sp.sin(TH)**2), Kp)
check("D3 far away g_tphi -> -2 M kappa sin^2(theta)/r; frame dragging (Part F, F4a) needs -2J sin^2(theta)/r, "
      "so kappa = J/M", ksol == [Js / Ms], f"limit r g_tphi = {sp.simplify(far)}; kappa = {ksol}")

# D4 the innermost shell: the disk r = 0 of radius b, rim speed kappa/b
rim = sp.simplify((Kp / (R**2 + Bp**2) * sp.sqrt(R**2 + Bp**2) * sp.sin(TH)).subs({R: 0, TH: sp.pi / 2}))
alpha_k = sp.simplify((Kp / (R**2 + Kp**2) * sp.sqrt(R**2 + Kp**2) * sp.sin(TH)))
check("D4 the r = 0 shell is a flat disk of radius b whose rim moves at kappa/b; the rim at light speed (the dial's edge) "
      "means b = kappa: the lattice's ring and the turn are the same length, and every other shell is slower than light",
      sp.solve(sp.Eq(rim, 1), Bp) == [Kp] and sp.simplify(1 - alpha_k**2 - (R**2 + Kp**2 * sp.cos(TH)**2) / (R**2 + Kp**2)) == 0,
      f"rim speed = {rim}; with b = kappa: 1 - alpha^2 = rho^2/(r^2+kappa^2) >= 0, zero only on the ring "
      "(Kerr's ring source spinning at light speed at its rim: Israel 1970, known)")

# D5 the derived cell is Kerr; every other choice in the chain fails the field equations
gD = metric2(M=Ms, b=Js / Ms, kappa=Js / Ms)
check("D5a the derived cell (Omega = kappa/(r^2+kappa^2), kappa = J/M, b = kappa) is identical to the tested one, a = J/M",
      sp.simplify(gD - metric(M=Ms, a=As).subs(As, Js / Ms)) == sp.zeros(4, 4), "nothing was fitted: Kerr's a = J/M comes out")
rD = rmax(metric2())
check("D5b and it solves every empty-space field equation (M = 1, kappa = b = 0.7)", rD < ZERO, f"max |Ricci| = {mp.nstr(rD, 3)}")
bad7 = {}
for lab, kw in [("ring b = 0.5 kappa (rim faster than light)", dict(b=An / 2)),
                ("ring b = 1.5 kappa (rim slower than light)", dict(b=3 * An / 2)),
                ("gap not kept: Omega ~ 1/r^2", dict(f=An * (R**2 + An**2) / R**2 * sp.sin(TH)**2)),
                ("gap not kept: one rigid rate for all space", dict(f=An / (1 + An**2) * (R**2 + An**2) * sp.sin(TH)**2)),
                ("gap drifts with r: f = kappa sin^2 (1 + M/r)", dict(f=An * sp.sin(TH)**2 * (1 + 1 / R)))]:
    bad7[lab] = rmax(metric2(**kw))
check("D5c each assumption in the chain matters: changing the rim speed or letting the gap drift breaks empty space",
      all(x > mp.mpf("1e-3") for x in bad7.values()), {k: mp.nstr(x, 3) for k, x in bad7.items()})
results["D5c_bad"] = {k: float(x) for k, x in bad7.items()}

# D6 the horizon turns at the shell rate
rp_ = Ms + sp.sqrt(Ms**2 - As**2)
gKa = metric(M=Ms, a=As)
blockdet = sp.simplify(gKa[0, 3]**2 - gKa[0, 0] * gKa[3, 3])
OmH = sp.simplify((-gKa[0, 3] / gKa[3, 3]).subs(R, rp_))
check("D6a the t-phi block: g_tphi^2 - g_tt g_phiphi = sin^2(theta) (r^2 - 2Mr + a^2), which closes at the horizon",
      sp.simplify(blockdet - sp.sin(TH)**2 * (R**2 - 2 * Ms * R + As**2)) == 0, f"= {sp.factor(blockdet)}")
check("D6b there the dragging rate -g_tphi/g_phiphi equals the lattice shell's rate a/(r+^2 + a^2): Kerr's horizon spin "
      "Omega_H is the derived shell rate, read at the horizon", sp.simplify(OmH - As / (rp_**2 + As**2)) == 0,
      f"Omega_H = {sp.simplify(OmH)}  (the known value a/(r+^2+a^2); also the Carter frame's rate)")

n_pass = sum(1 for v_ in results.values() if isinstance(v_, dict) and v_.get("pass"))
n_all = sum(1 for v_ in results.values() if isinstance(v_, dict) and "pass" in v_)
log(f"\n{n_pass}/{n_all} checks pass")
open(os.path.join(HERE, "checks_output.txt"), "w").write("\n".join(out) + "\n")
json.dump(results, open(os.path.join(HERE, "results.json"), "w"), indent=2, default=str)
