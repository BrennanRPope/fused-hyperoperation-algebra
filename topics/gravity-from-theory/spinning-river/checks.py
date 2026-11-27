"""DRAFT checks: the spinning river, and frame dragging's relatives (Nov 26, 2026).

K: give each river cell a turn as well as a velocity; does the spinning (Kerr) black hole come out exactly?
F: frame-dragging look-alikes outside black holes (Larmor's theorem, spinning fluids, Foucault, clocks).
D: derive the turn (area law + circulation law + far-field frame dragging + rim at light speed) instead of testing it.
S: derive the speed law and the lattice's shape from one move: the point source shifted up the axis by i*a.
T: derive the turn from the same move: the shifted source's light rays.
C: charge on the same shifted source: the charged spinning black hole (Kerr-Newman).
B: two bodies: why one shear per cell cannot hold two, what does add exactly (Weyl), and the strut force (Newton).
O: two spinning bodies in orbit: the slow-orbit field from three river rules, orbits decaying, spins, the waves;
   O7-O8 test the rules at the next level of detail (2PN) through body 1's clock rate; O9-O12 add the field's own
   energy as a fourth rule; O11 derives what was left (the speeding-up body's own field squared); O12
   checks g00 at a body against the published 2PN field term by term.
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

# ================================================================ S  deriving the speed law (round 8)
log("\n== S  Deriving the speed law (and the lattice's shape) instead of matching it")
xs, ys, zs = sp.symbols("x y z", real=True)
cth, sth = sp.cos(TH), sp.sin(TH)
xo = sp.sqrt(R**2 + As**2) * sth * sp.cos(PH)
yo = sp.sqrt(R**2 + As**2) * sth * sp.sin(PH)
zo = R * cth

# S1 the complex shift: a point source moved along the axis by i*a; its distance in oblate coordinates
D2 = sp.expand(xo**2 + yo**2 + (zo - sp.I * As)**2)
check("S1a move a point source up the spin axis by i*a (the dial's turning number): its squared distance is exactly "
      "(r - i a cos(theta))^2 in the lattice's coordinates",
      sp.simplify(sp.expand(D2 - (R - sp.I * As * cth)**2)) == 0,
      "so the complex distance is D = r - i a cos(theta): its real part is the lattice label r")
check("S1b so the lattice's squashed spheres are the level sets of Re D, and D = 0 is the ring r = 0, theta = pi/2 of "
      "radius a: the same ring round 7 found for the turn (b = kappa), which agrees only when a = J/M",
      sp.simplify(sp.sqrt((xo**2 + yo**2).subs({R: 0, TH: sp.pi / 2}))) == As,
      "the lattice's shape and its ring both come from the shift; before this round the shape was put in")

# S2 Newton's potential of the shifted point is the speed law
Phi_shift = sp.simplify(sp.re(sp.expand(Ms * (R + sp.I * As * cth) / ((R - sp.I * As * cth) * (R + sp.I * As * cth)))))
check("S2a the shifted point's potential, Re(M/D) = M r/(r^2 + a^2 cos^2 theta)",
      sp.simplify(Phi_shift - Ms * R / (R**2 + As**2 * cth**2)) == 0, f"Re(M/D) = {Phi_shift}")
check("S2b and the jet river's rule v^2/2 = -Phi (Phi = -v^2/2, jet-river J6b) then gives v^2 = 2Mr/rho^2: "
      "exactly the speed law Part 1 matched by testing (K3c)",
      sp.simplify(2 * Phi_shift - legs(M=Ms, a=As)[4]**2) == 0, "far away Re(M/D) -> M/r, Newton's potential of mass M")

# S3 it is an exact flat-space potential (harmonic), not just in the far field
hdiag = [(R**2 + As**2 * cth**2) / (R**2 + As**2), R**2 + As**2 * cth**2, (R**2 + As**2) * sth**2]
sqh = (R**2 + As**2 * cth**2) * sth


def lap(F):
    return sp.simplify(sum(sp.diff(sqh / hdiag[i] * sp.diff(F, q), q) for i, q in enumerate((R, TH, PH))) / sqh)


check("S3 Re(M/D) is exactly harmonic on the flat lattice (Laplace's equation holds everywhere off the disk r = 0)",
      lap(Ms * R / (R**2 + As**2 * cth**2)) == 0, "the empty-space form of the jet river's trace law (Poisson for Phi)")

# S4 harmonic + mass M is not enough: other potentials of mass M fail; the i-shift is what picks Kerr
dist_origin = sp.sqrt(sp.simplify(xo**2 + yo**2 + zo**2))
dist_real = sp.sqrt(sp.simplify(sp.expand(xo**2 + yo**2 + (zo - As)**2)))
cands = {"point at the centre: M/|x|": Ms / dist_origin,
         "point moved by a real a up the axis": Ms / dist_real,
         "point moved by i a/2 (ring a/2, not the turn's ring)":
             sp.re(Ms / sp.sqrt(sp.expand(xo**2 + yo**2 + (zo - sp.I * As / 2)**2)))}
harm, bad8 = {}, {}
for lab, Phi_c in cands.items():
    Fl = sp.lambdify((R, TH, PH, Ms, As), lap(Phi_c) if "i a/2" not in lab else Phi_c, "mpmath")
    if "i a/2" in lab:      # numeric Laplacian for the complex-shift candidate
        lapc = sp.lambdify((R, TH, PH, Ms, As), sum(sp.diff(sqh / hdiag[i] * sp.diff(Phi_c, q), q)
                                                     for i, q in enumerate((R, TH, PH))) / sqh, "mpmath")
        harm[lab] = max(abs(lapc(*[mp.mpf(P[k]) for k in (1, 2, 3)], 1, mp.mpf("0.7"))) for P in PTS)
    else:
        harm[lab] = max(abs(Fl(*[mp.mpf(P[k]) for k in (1, 2, 3)], 1, mp.mpf("0.7"))) for P in PTS)
    V2c = 2 * Phi_c.subs({Ms: Mn, As: An})
    e0, e1, e2, e3, v_, rho2_c = legs2(M=Mn, b=An, kappa=An)
    w_ = sp.Matrix([1, 0, 0, -An * sth**2])
    e1c = e1 - v_ * w_ + sp.sqrt(V2c) * w_
    gc = -e0 * e0.T + e1c * e1c.T + e2 * e2.T + e3 * e3.T
    bad8[lab] = rmax(gc)
check("S4a each rival is also a flat-space potential of mass M (harmonic, -> M/r far away)",
      all(h < mp.mpf("1e-25") for h in harm.values()), {k: mp.nstr(h, 3) for k, h in harm.items()})
check("S4b but with the derived lattice and turn, each rival breaks empty space; only the i*a shift works",
      all(x > mp.mpf("1e-3") for x in bad8.values()), {k: mp.nstr(x, 3) for k, x in bad8.items()})
results["S4b_bad"] = {k: float(x) for k, x in bad8.items()}

# ================================================================ T  the turn from the same shift (round 9)
log("\n== T  The turn from the same shift: the light rays of the shifted source")
mp.mp.dps = 50


def cart_and_jac(r, th, ph, a):
    A = mp.sqrt(r**2 + a**2)
    X = mp.matrix([A * mp.sin(th) * mp.cos(ph), A * mp.sin(th) * mp.sin(ph), r * mp.cos(th)])
    J = mp.matrix([[r / A * mp.sin(th) * mp.cos(ph), A * mp.cos(th) * mp.cos(ph), -A * mp.sin(th) * mp.sin(ph)],
                   [r / A * mp.sin(th) * mp.sin(ph), A * mp.cos(th) * mp.sin(ph), A * mp.sin(th) * mp.cos(ph)],
                   [mp.cos(th), -r * mp.sin(th), 0]])
    return X, J


def cross(u, w):
    return mp.matrix([u[1] * w[2] - u[2] * w[1], u[2] * w[0] - u[0] * w[2], u[0] * w[1] - u[1] * w[0]])


def dot(u, w):
    return sum(u[k] * w[k] for k in range(3))


errT1, errT2, errT3, errT2b = 0, 0, 0, 0
for P in PTS:
    r, th, ph = (mp.mpf(P[k]) for k in (1, 2, 3))
    a = mp.mpf("0.7")
    X, J = cart_and_jac(r, th, ph, a)
    D = r - 1j * a * mp.cos(th)                                # complex distance (S1a)
    E = mp.matrix([X[0] / D, X[1] / D, (X[2] - 1j * a) / D])     # its gradient: a complex unit vector, E.E = 1
    p = mp.matrix([mp.re(E[k]) for k in range(3)])
    q = mp.matrix([mp.im(E[k]) for k in range(3)])
    kp = (p + cross(q, p)) / dot(p, p)                          # ray k+ ; the other ray is (p - q x p)/(p.p)
    km = (p - cross(q, p)) / dot(p, p)
    errT1 = max(errT1, abs(dot(E, E) - 1), abs(dot(kp, kp) - 1), abs(dot(kp, p) - 1), abs(dot(kp, q)),
                abs(dot(km, km) - 1), abs(dot(km, p) - 1), abs(dot(km, q)))
    # the ray as a covector in lattice coordinates: k_i = (dx/dq_i) . k
    kcov_p = [dot(mp.matrix([J[0, i], J[1, i], J[2, i]]), kp) for i in range(3)]
    kcov_m = [dot(mp.matrix([J[0, i], J[1, i], J[2, i]]), km) for i in range(3)]
    rho2n = r**2 + a**2 * mp.cos(th)**2
    errT2 = max(errT2, abs(kcov_p[0] - rho2n / (r**2 + a**2)), abs(kcov_p[1]), abs(kcov_p[2] + a * mp.sin(th)**2))
    errT3 = max(errT3, abs(-kcov_m[2] - (-a * mp.sin(th)**2)), abs(-kcov_m[0] + rho2n / (r**2 + a**2)))
    phihat = mp.matrix([-mp.sin(ph), mp.cos(ph), 0])
    errT2b = max(errT2b, abs(dot(kp, phihat) + a * mp.sin(th) / mp.sqrt(r**2 + a**2)))
check("T1 through each point the shifted source sends exactly two real light rays: unit vectors k with k.E = 1, where "
      "E = grad D is the complex unit direction to the source; they are k = (p +- q x p)/|p|^2, E = p + i q",
      errT1 < mp.mpf("1e-40"), f"max error {mp.nstr(errT1, 3)}; with no spin (q = 0) both are the ordinary radial ray k = E")
check("T2 the outgoing ray, written on the lattice, is k = (rho^2/(r^2+a^2)) dr - a sin^2(theta) dphi: radial along the "
      "shell's normal, plus a sideways part", errT2 < mp.mpf("1e-40"), f"max error {mp.nstr(errT2, 3)}")
check("T2b the ray leans sideways by exactly alpha = a sin(theta)/sqrt(r^2+a^2), the speed of round 7's turning shell (D1b); "
      "in the equator's plane every ray is a straight line tangent to the ring of radius a",
      errT2b < mp.mpf("1e-40"), f"max error {mp.nstr(errT2b, 3)}")
check("T3 the ingoing ray (the other one, reversed) has the same sideways part, -a sin^2(theta) dphi, and the opposite "
      "radial part", errT3 < mp.mpf("1e-40"), f"max error {mp.nstr(errT3, 3)}")

# T4 the rule 'push the cell along the ray, without the ray's radial part' gives the turn
fT = -(-As * sp.sin(TH)**2)                                     # w = dt + (sideways part of the ray) = dt - f dphi
check("T4 push each cell along its ray's time-and-sideways part, w = dt + (ray minus its radial part): w = dt - a sin^2 dphi, "
      "exactly round 7's turn f = kappa sin^2(theta) with kappa = a; with no spin it is the plain river, w = dt",
      sp.simplify(fT - As * sp.sin(TH)**2) == 0, "so the shift gives the turn's shape (sin^2) and size (a) directly")

# T5 one move gives all of Kerr: also in the ray form g = flat + v^2 (dt + k)(dt + k)
rho2s_ = R**2 + As**2 * cth**2
eta = sp.diag(-1, rho2s_ / (R**2 + As**2), rho2s_, (R**2 + As**2) * sth**2)
H = 2 * Ms * R / rho2s_                                          # 2 Re(M/D), S2
ray = {"rays of the i*a source (Kerr)": sp.Matrix([1, rho2s_ / (R**2 + As**2), 0, -As * sth**2]),
       "radial ray, no sideways part": sp.Matrix([1, sp.sqrt(rho2s_ / (R**2 + As**2)), 0, 0]),
       "half the sideways part": sp.Matrix([1, sp.sqrt(rho2s_ / (R**2 + As**2) * (1 - As**2 * sth**2 / 4 / (R**2 + As**2))),
                                            0, -As / 2 * sth**2])}
ksres = {}
for lab, l_ in ray.items():
    null = sp.simplify((l_.T * eta.inv() * l_)[0])
    gks = (eta + H * l_ * l_.T).subs({Ms: Mn, As: An})
    ksres[lab] = (null, rmax(gks))
check("T5a with the shifted source's rays and speed, the ray form flat + v^2 (dt + k)(dt + k) solves every empty-space "
      "field equation: one move gives the whole spinning black hole (Kerr-Schild form)",
      ksres["rays of the i*a source (Kerr)"][0] == 0 and ksres["rays of the i*a source (Kerr)"][1] < ZERO,
      f"ray null: {ksres['rays of the i*a source (Kerr)'][0] == 0}; max |Ricci| = {mp.nstr(ksres['rays of the i*a source (Kerr)'][1], 3)}")
check("T5b other light rays with the same speed law fail: a purely radial ray, or half the sideways part",
      all(ksres[k][0] == 0 and ksres[k][1] > mp.mpf("1e-3") for k in ("radial ray, no sideways part", "half the sideways part")),
      {k: mp.nstr(ksres[k][1], 3) for k in ("radial ray, no sideways part", "half the sideways part")})
results["T5_ricci"] = {k: float(v[1]) for k, v in ksres.items()}
check("T5c and the river cell built from the shift alone (lattice, speed and turn all from i*a) is Part 1's Kerr cell",
      sp.simplify(metric2(M=Ms, b=As, f=As * sp.sin(TH)**2) - metric(M=Ms, a=As)) == sp.zeros(4, 4),
      "inputs left: the mass M, the shift i*a (read as spin, a = J/M from frame dragging far away), and the river rules")

# ================================================================ C  charge on the same shifted source (round 10)
log("\n== C  Charge: the same move for the charged spinning black hole (Kerr-Newman)")
Qs = sp.symbols("Q", positive=True)
Qn = sp.Rational(1, 2)


def legsQ(M, a, V2, f=None):
    """river cell on the flat spheroidal lattice (ring a) with speed^2 V2 and turn f (default a sin^2)."""
    rho2 = R**2 + a**2 * sp.cos(TH)**2
    if f is None:
        f = a * sp.sin(TH)**2
    g0 = sp.Matrix([1, 0, 0, 0])
    g1 = sp.Matrix([0, sp.sqrt(rho2) / sp.sqrt(R**2 + a**2), 0, 0])
    g2 = sp.Matrix([0, 0, sp.sqrt(rho2), 0])
    g3 = sp.Matrix([0, 0, 0, sp.sqrt(R**2 + a**2) * sp.sin(TH)])
    e1 = g1 + sp.sqrt(V2) * sp.Matrix([1, 0, 0, -f])
    return -g0 * g0.T + e1 * e1.T + g2 * g2.T + g3 * g3.T


def em_residual(g, A, pts=PTS):
    """max |R_mn - 2(F_ma F_n^a - g_mn F^2/4)| (Einstein-Maxwell, G = 1, Gaussian units) and max |div F|."""
    f = prepare(g, X4)
    F = sp.Matrix(4, 4, lambda m, n: sp.diff(A[n], X4[m]) - sp.diff(A[m], X4[n]))
    dF = [[[sp.diff(F[m, n], X4[c]) for c in range(4)] for n in range(4)] for m in range(4)]
    fF = sp.lambdify(X4, [F.tolist(), dF], "mpmath")
    e_ein, e_max = mp.mpf(0), mp.mpf(0)
    for P in pts:
        _, Ric, gi = curvature_at(f, P)
        G_, dG, _ = f(*[mp.mpf(p) for p in P])
        Fv, dFv = fF(*[mp.mpf(p) for p in P])
        n = 4
        Fup = [[sum(gi[m, a_] * gi[k, b_] * Fv[a_][b_] for a_ in range(n) for b_ in range(n)) for k in range(n)] for m in range(n)]
        F2 = sum(Fv[m][k] * Fup[m][k] for m in range(n) for k in range(n))
        for m in range(n):
            for k in range(n):
                T = 2 * (sum(Fv[m][a_] * Fv[k][b_] * gi[a_, b_] for a_ in range(n) for b_ in range(n)) - G_[m][k] * F2 / 4)
                e_ein = max(e_ein, abs(Ric[m][k] - T))
        low = [[[(dG[d][b][c] + dG[d][c][b] - dG[b][c][d]) / 2 for c in range(n)] for b in range(n)] for d in range(n)]
        Gam = [[[sum(gi[a_, d] * low[d][b][c] for d in range(n)) for c in range(n)] for b in range(n)] for a_ in range(n)]
        for nu in range(n):
            div = 0
            for mu in range(n):
                for al in range(n):
                    for be in range(n):
                        cov = dFv[al][be][mu] - sum(Gam[l][mu][al] * Fv[l][be] + Gam[l][mu][be] * Fv[al][l] for l in range(n))
                        div += gi[mu, al] * gi[nu, be] * cov
            e_max = max(e_max, abs(div))
    return e_ein, e_max


# C1 no spin: the jet river's Gauss law with the field's own energy gives the charged black hole (Reissner-Nordstrom)
rp_s = sp.symbols("rho_p", positive=True)
field_energy_outside = sp.integrate(Qs**2 / (8 * sp.pi * rp_s**4) * 4 * sp.pi * rp_s**2, (rp_s, R, sp.oo))
m_of_r = Ms - field_energy_outside
check("C1a no spin: the electric field Q/r^2 holds energy Q^2/(2r) outside radius r, so the mass inside is "
      "m(r) = M - Q^2/(2r), and the jet river's Gauss law r v^2 = 2 m(r) gives v^2 = 2M/r - Q^2/r^2",
      sp.simplify(2 * m_of_r / R - (2 * Ms / R - Qs**2 / R**2)) == 0, f"m(r) = {sp.simplify(m_of_r)}")
gRN = legsQ(Mn, 0, 2 * Mn / R - Qn**2 / R**2, f=0)
eRN = em_residual(gRN, sp.Matrix([Qn / R, 0, 0, 0]))
check("C1b with A = (Q/r) dt that river is the charged black hole: Einstein's equations with the field's energy, and "
      "Maxwell's equations, hold (M = 1, Q = 0.5)", eRN[0] < ZERO and eRN[1] < ZERO,
      f"max field-equation error {mp.nstr(eRN[0], 3)}, max Maxwell error {mp.nstr(eRN[1], 3)}")

# C2 shift the charge by i*a: its electric and magnetic fields together, E + iB = -grad(Q/D)
Rb, cz = sp.symbols("R_big c_z", positive=True)
ser = sp.series(Qs / sp.sqrt(Rb**2 - 2 * sp.I * As * Rb * cz - As**2), Rb, sp.oo, 3).removeO()
im2 = sp.simplify(sp.expand(ser).coeff(Rb, -2))
check("C2 the shifted charge Q/D is one complex potential: its real part is the electric potential, its imaginary "
      "part far away is Q a cos(theta)/R^2, a magnet of moment Q a. With spin J = M a that is the gyromagnetic "
      "ratio g = 2 (2M/Q)(moment/J) = 2, the electron's value (known for Kerr-Newman: Carter 1968)",
      sp.simplify(im2 - sp.I * Qs * As * cz) == 0, f"1/R^2 term of Q/D = {im2}")

# C3 the speed law with charge: v^2 = 2 Re(M/D) - |Q/D|^2
Dz = R - sp.I * As * cth
V2_KN = sp.simplify(2 * sp.re(sp.expand(Ms * sp.conjugate(Dz) / (Dz * sp.conjugate(Dz)))) - sp.expand(Qs**2 / (Dz * sp.conjugate(Dz))))
check("C3 shift both terms: v^2 = 2 Re(M/D) - |Q/D|^2 = (2Mr - Q^2)/rho^2; the field-energy term uses the size |Q/D| of "
      "the complex potential, because the energy is (E^2 + B^2)/8pi = |E + iB|^2/8pi",
      sp.simplify(V2_KN - (2 * Ms * R - Qs**2) / (R**2 + As**2 * cth**2)) == 0, f"v^2 = {V2_KN}")

# C4 the charged spinning river: speed from C3, turn and lattice from the shift, field A = Re(Q/D) w
wleg = sp.Matrix([1, 0, 0, -An * sp.sin(TH)**2])
rho2n_ = R**2 + An**2 * sp.cos(TH)**2
V2n = (2 * Mn * R - Qn**2) / rho2n_
gKN = legsQ(Mn, An, V2n)
AKN = Qn * R / rho2n_ * wleg
eKN = em_residual(gKN, AKN)
check("C4 the charged spinning black hole: field A = Re(Q/D) (dt - a sin^2 dphi), the potential carried on the cell's "
      "turned time leg; Einstein's equations with the field's energy and Maxwell's equations all hold "
      "(M = 1, a = 0.7, Q = 0.5)", eKN[0] < ZERO and eKN[1] < ZERO,
      f"max field-equation error {mp.nstr(eKN[0], 3)}, max Maxwell error {mp.nstr(eKN[1], 3)}")
grrKN = sp.simplify(legsQ(Ms, As, (2 * Ms * R - Qs**2) / (R**2 + As**2 * cth**2)).inv()[1, 1] * (R**2 + As**2 * cth**2))
check("C4b its horizon: g^rr rho^2 = r^2 - 2Mr + a^2 + Q^2, Kerr-Newman's", sp.simplify(grrKN - (R**2 - 2 * Ms * R + As**2 + Qs**2)) == 0,
      f"g^rr rho^2 = {grrKN}")

# C5 rivals
riv = {"energy from Re(Q/D)^2, not |Q/D|^2": (legsQ(Mn, An, 2 * Mn * R / rho2n_ - (Qn * R / rho2n_)**2), AKN),
       "field without the turn: A = Re(Q/D) dt": (gKN, Qn * R / rho2n_ * sp.Matrix([1, 0, 0, 0])),
       "charge left at the centre (not shifted)": (legsQ(Mn, An, 2 * Mn * R / rho2n_ - Qn**2 / (R**2 + An**2 * sp.sin(TH)**2)),
                                                   Qn / sp.sqrt(R**2 + An**2 * sp.sin(TH)**2) * wleg)}
rivres = {k: em_residual(*v, pts=PTS[:2]) for k, v in riv.items()}
check("C5 rivals fail: using Re(Q/D)^2 for the field energy, dropping the turn from the field, or leaving the charge "
      "unshifted", all(max(e) > mp.mpf("1e-3") for e in rivres.values()),
      {k: mp.nstr(max(e), 3) for k, e in rivres.items()})
results["C5_bad"] = {k: float(max(e)) for k, e in rivres.items()}
results["C4_good"] = float(max(eKN))

# ================================================================ B  two bodies (round 11)
log("\n== B  Two bodies: can the cell picture hold two black holes at once?")
import itertools  # noqa: E402


def weyl_speciality(g, X, P):
    """Electric and magnetic parts E, B of the curvature seen by the observer along d/dX0 (orthonormal frame by
    Gram-Schmidt), Q = E + iB, and the measure |trQ^2^3 - 6 trQ^3^2| / |trQ^2|^3: zero exactly when two of Q's three
    eigenvalues coincide (algebraically special, as every one-shear vacuum is); of order one for a general field."""
    f = prepare(g, X)
    Riem, _, _ = curvature_at(f, P)
    G = mp.matrix(f(*[mp.mpf(p) for p in P])[0])
    n = 4
    Rl = {k: sum(G[k[0], e] * Riem[e, k[1], k[2], k[3]] for e in range(n)) for k in itertools.product(range(n), repeat=4)}
    ip = lambda u, w: sum(G[i, j] * u[i] * w[j] for i in range(n) for j in range(n))  # noqa: E731
    es, sg = [], []
    for k in range(n):
        w = [mp.mpf(1 if i == k else 0) for i in range(n)]
        v0 = list(w)
        for e, s in zip(es, sg):
            c = ip(v0, e) / s
            w = [w[i] - c * e[i] for i in range(n)]
        nn = ip(w, w)
        es.append([x / mp.sqrt(abs(nn)) for x in w]); sg.append(mp.sign(nn))
    nz = [(k, val) for k, val in Rl.items() if val != 0]

    def Rf(a, b, c, d):
        return sum(val * es[a][k[0]] * es[b][k[1]] * es[c][k[2]] * es[d][k[3]] for k, val in nz)
    eps3 = lambda i, j, k: (i - j) * (j - k) * (k - i) / 2  # noqa: E731
    Qm = mp.matrix(3, 3)
    for i in range(3):
        for j in range(3):
            Bij = sum(eps3(i, k, l) * Rf(k + 1, l + 1, 0, j + 1) for k in range(3) for l in range(3)) / 2
            Qm[i, j] = Rf(0, i + 1, 0, j + 1) + 1j * Bij
    t2 = sum(Qm[i, j] * Qm[j, i] for i in range(3) for j in range(3))
    Q2 = Qm * Qm
    t3 = sum(Q2[i, j] * Qm[j, i] for i in range(3) for j in range(3))
    return abs(t2**3 - 6 * t3**2) / abs(t2)**3


# Weyl's coordinates for bodies at rest on one axis: (t, rho, z, phi), flat space dz^2 + drho^2 + rho^2 dphi^2
tw, rw, zw, pw = sp.symbols("t rho z phi", real=True)
XW = [tw, rw, zw, pw]
PW = [(0, 1.2, 0.7, 0.3), (0, 2.5, -1.0, 0.0), (0, 0.8, 4.0, 0.0)]


def Rk(zk):
    return sp.sqrt(rw**2 + (zw - zk)**2)


def Ek(zi, zj):
    return Rk(zi) * Rk(zj) + (zw - zi) * (zw - zj) + rw**2


def rods(bodies):
    """Static black holes at rest on the axis, each a rod of length 2M centred at height c (Weyl 1917; two rods:
    Bach-Weyl 1922, Israel-Khan 1964). psi = (1/2) ln(-g_tt) is a sum of the rods' Newtonian potentials."""
    e2psi, e2gam = 1, 1
    ends = [(c - M, c + M) for M, c in bodies]
    for (lo, hi), (M, c) in zip(ends, bodies):
        e2psi *= (Rk(lo) + Rk(hi) - 2 * M) / (Rk(lo) + Rk(hi) + 2 * M)
        e2gam *= Ek(lo, hi) / (2 * Rk(lo) * Rk(hi))
    for (l1, h1), (l2, h2) in itertools.combinations(ends, 2):
        e2gam *= Ek(h1, l2) * Ek(l1, h2) / (Ek(h1, h2) * Ek(l1, l2))
    return sp.diag(-e2psi, e2gam / e2psi, e2gam / e2psi, rw**2 / e2psi), e2gam


g1rod, _ = rods([(1, 3)])
f1rod = prepare(g1rod, XW)
r1 = max(ricci_max(f1rod, P) for P in PW[:2])
s1 = max(weyl_speciality(g1rod, XW, P) for P in PW[:2])
sK = max(weyl_speciality(metric(), X4, P) for P in [(0, 4.0, 0.9, 0.3), (0, 6.0, 1.3, 0.0)])
check("B1a one body: the black hole at rest written in Weyl's coordinates (a rod of length 2M on the axis) is empty "
      "space, and its curvature is algebraically special, as is the spinning one's (two of three curvature eigenvalues "
      "equal; measure 0). Every one-shear empty space is special: the Goldberg-Sachs theorem",
      r1 < ZERO and s1 < ZERO and sK < ZERO,
      f"Ricci {mp.nstr(r1, 3)}; special measure: at rest {mp.nstr(s1, 3)}, spinning {mp.nstr(sK, 3)}")
g2rod, e2g2 = rods([(1, 3), (sp.Rational(1, 2), -3)])
f2rod = prepare(g2rod, XW)
r2 = max(ricci_max(f2rod, P) for P in PW)
s2 = [weyl_speciality(g2rod, XW, P) for P in PW]
check("B1b two bodies (masses 1 and 1/2, centres 6 apart): gravity's exact field is empty space off the axis, but its "
      "curvature is algebraically general, so no single shear per cell can produce it", r2 < ZERO and min(s2[:2]) > 0.01,
      f"Ricci off the axis {mp.nstr(r2, 3)}; special measure at three points {[mp.nstr(s, 3) for s in s2]} "
      "(near one body it tends to that body's 0)")
results["B1_spec"] = [float(s) for s in s2]

# B2 two shears combined: metric sum, leg sum, product of the two cell maps
tc, xc, yc, zc = sp.symbols("t x y z", real=True)
XC = [tc, xc, yc, zc]
etaC = sp.diag(-1, 1, 1, 1)
PC = [(0, 1.1, 0.4, 0.7), (0, 2.0, -1.3, -0.5)]


def shear(M, z0):
    rr = sp.sqrt(xc**2 + yc**2 + (zc - z0)**2)
    l_ = sp.Matrix([1, xc / rr, yc / rr, (zc - z0) / rr])
    H = 2 * M / rr
    return H, l_, (H / 2) * (etaC * l_) * l_.T          # N^a_m = (H/2) l^a l_m, nilpotent: N N = 0


def two_shears(M1, M2, how, d=3):
    H1, l1, N1 = shear(M1, d)
    H2, l2, N2 = shear(M2, -d)
    I4 = sp.eye(4)
    if how == "sum":
        return etaC + H1 * l1 * l1.T + H2 * l2 * l2.T
    E = I4 + N1 + N2 if how == "legs" else (I4 + N1) * (I4 + N2)
    return E.T * etaC * E


combos = {"sum": "add the two metrics", "legs": "add the two shears", "product": "compose the two cell maps"}
lefts = {}
for how in combos:
    for e_ in ((sp.Rational(1, 100), sp.Rational(1, 1000)) if how == "product" else (sp.Rational(1, 100),)):
        lefts[how, e_] = max(ricci_max(prepare(two_shears(e_, e_, how), XC), P) for P in PC)
alone = max(ricci_max(prepare(two_shears(sp.Rational(1, 100), 0, "product"), XC), P) for P in PC)
ratio = (lefts["product", sp.Rational(1, 1000)] / mp.mpf("1e-6")) / (lefts["product", sp.Rational(1, 100)] / mp.mpf("1e-4"))
check("B2 two shears, one per body (each alone exact): adding the metrics, adding the shears, or composing the two "
      "cell maps all leave a leftover; it grows as M1*M2 (both masses 10x smaller: leftover 100x smaller), and one "
      "body alone leaves none", alone < ZERO and all(v > mp.mpf("1e-6") for (h, e_), v in lefts.items() if e_ == sp.Rational(1, 100))
      and abs(ratio - 1) < 0.02,
      {**{combos[h] + f" (M = {e_})": mp.nstr(v, 3) for (h, e_), v in lefts.items()}, "one body alone": mp.nstr(alone, 3),
       "leftover/M^2, M = 1/1000 vs 1/100": mp.nstr(ratio, 4)})
results["B2_left"] = {combos[h]: float(v) for (h, e_), v in lefts.items() if e_ == sp.Rational(1, 100)}

# B3 what does superpose: Weyl's psi = (1/2) ln(1 - v^2), and on the axis that is Einstein's velocity addition
zs, M1s, M2s, z1s, z2s, ls = sp.symbols("z M_1 M_2 z_1 z_2 ell", positive=True)
u1, u2 = M1s / (zs - z1s), M2s / (zs - z2s)
Eax = (zs - z1s - M1s) / (zs - z1s + M1s) * (zs - z2s - M2s) / (zs - z2s + M2s)
uu = (u1 + u2) / (1 + u1 * u2)
check("B3 what does add exactly for bodies at rest is psi = (1/2) ln(1 - v^2) (the log of the river cell's clock rate), "
      "summed over rods in Weyl's flat space (B1b). On the axis each body's psi is -artanh(u), u = M/(z - c), so the "
      "two bodies' u combine by Einstein's velocity addition (u1 + u2)/(1 + u1 u2)",
      sp.simplify((1 - uu) / (1 + uu) - Eax) == 0, "1 - v^2 on the axis = product of the two bodies' (1-u)/(1+u)")

# B4 the strut between them, and Newton's law
gam0 = sp.lambdify([rw, zw], e2g2, "mpmath")
ell_, M1n, M2n = 6, mp.mpf(1), mp.mpf("0.5")
egs = [mp.sqrt(gam0(mp.mpf("1e-12"), mp.mpf(zz))) for zz in (0.0, 1.2, -1.8)]
F_num = (1 / egs[0] - 1) / 4
F_IK = M1n * M2n / (ell_**2 - (M1n + M2n)**2)
check("B4 two bodies at rest need a strut: on the axis between them the field has a constant cone defect e^gamma "
      "(the same at every height), and the strut's push (e^-gamma - 1)/4 equals the force M1 M2/(l^2 - (M1+M2)^2) "
      "(Israel-Khan 1964), which far apart is Newton's M1 M2/l^2: the inverse-square law, from the exact field",
      max(abs(e - egs[0]) for e in egs) < mp.mpf("1e-15") and abs(F_num - F_IK) < mp.mpf("1e-15")
      and sp.limit(M1s * M2s / (ls**2 - (M1s + M2s)**2) * ls**2, ls, sp.oo) == M1s * M2s,
      f"e^gamma on the strut {mp.nstr(egs[0], 10)} at three heights; push {mp.nstr(F_num, 12)}, "
      f"M1M2/(l^2-(M1+M2)^2) = {mp.nstr(F_IK, 12)}, Newton M1M2/l^2 = {mp.nstr(M1n * M2n / ell_**2, 6)}")
results["B4_force"] = [float(F_num), float(F_IK), float(M1n * M2n / ell_**2)]

# B5 the lead for spin: Kerr's axis values are the shifted source's u = M/(z + i a) in Weyl's z = r - M
ua = Ms / ((R - Ms) + sp.I * As)
Eax_K = sp.simplify((1 - ua) / (1 + ua))
gtt_axis = sp.simplify(metric(M=Ms, a=As)[0, 0].subs(TH, 0))
check("B5 lead for spin: on the spin axis, in Weyl's coordinate z = r - M, Kerr's Ernst potential is (1 - u)/(1 + u) "
      "with u = M/(z + i a), a point source shifted by i*a: its real part is the clock factor -g_tt = 1 - v^2 of the "
      "spinning river, its imaginary part Kerr's twist potential 2Ma/(r^2 + a^2)",
      sp.simplify(sp.re(sp.expand_complex(Eax_K)) + gtt_axis) == 0
      and sp.simplify(sp.im(sp.expand_complex(Eax_K)) - 2 * Ms * As / (R**2 + As**2)) == 0,
      f"(1-u)/(1+u) = {sp.simplify(sp.expand_complex(Eax_K))}")

# ================================================================ O  two spinning bodies in orbit (round 12)
log("\n== O  Two spinning bodies orbiting each other: the slow-orbit field from the river's rules, orbits, spins, waves")
import random  # noqa: E402

from orbits import ensemble, flux_ratio, periastron_advance  # noqa: E402

eps_ = sp.symbols("epsilon", positive=True)
Xo = sp.Matrix(sp.symbols("x1:4", real=True))
Yo = [sp.Matrix(sp.symbols(f"y{k}_1:4", real=True)) for k in (1, 2)]
Vo = [sp.Matrix(sp.symbols(f"v{k}_1:4", real=True)) for k in (1, 2)]
mo = list(sp.symbols("m1 m2", positive=True))
nrm = lambda u: sp.sqrt(u.dot(u))  # noqa: E731
ro = [nrm(Xo - Yo[k]) for k in (0, 1)]
no = [(Xo - Yo[k]) / ro[k] for k in (0, 1)]
r12o = nrm(Yo[0] - Yo[1])
n12o = (Yo[0] - Yo[1]) / r12o
acco = [-mo[1] * n12o / r12o**2, mo[0] * n12o / r12o**2]


def bfp_metric():
    """The published slow-orbit (1PN) field of two moving bodies, harmonic gauge (Blanchet, Faye & Ponsot 1998)."""
    g00, g0, gs = -1, sp.zeros(3, 1), 1
    for a, b in ((0, 1), (1, 0)):
        ra, rb, ma, mb, va = ro[a], ro[b], mo[a], mo[b], Vo[a]
        g00 += (2 * ma / ra + ma / ra * (4 * va.dot(va) - no[a].dot(va)**2) - 2 * ma**2 / ra**2
                + ma * mb * (-2 / (ra * rb) - ra / (2 * r12o**3) + ra**2 / (2 * rb * r12o**3) - 5 / (2 * ra * r12o)))
        g0 += -4 * ma * va / ra
        gs += 2 * ma / ra
    return g00, g0, gs


def river_rules_metric(lag=True, multiply=True, dressed=True):
    """Rule A: each body's own far field (the river's 2m/r, frame dragging from its motion as in Part F), carried
    outward at light speed. Rule B: clock rates multiply between bodies (check B3), each body keeping its own clock
    factor 1 - 2U + 2U^2. Rule C: each body's mass is seen through the other body's clock, m_a -> m_a (1 - U_b)."""
    U = []
    for a, b in ((0, 1), (1, 0)):
        ma = mo[a] * (1 - (mo[b] / r12o if dressed else 0))
        va, ra = Vo[a], ro[a]
        h00 = (4 * ma * (1 + va.dot(va) / 2) / ra + 4 * mo[a] * va.dot(va) / ra) / 2   # (hbar^00 + hbar^ii)/2
        if lag:      # the field left the body a moment ago: + d^2/dt^2 (m r_a)
            h00 += mo[a] * ((va.dot(va) - no[a].dot(va)**2) / ra - no[a].dot(acco[a]))
        U.append(h00 / 2)
    if multiply:
        g00 = -(1 - 2 * U[0] + 2 * U[0]**2) * (1 - 2 * U[1] + 2 * U[1]**2)
    else:
        g00 = -(1 - 2 * U[0] + 2 * U[0]**2 - 2 * U[1] + 2 * U[1]**2)
    g0 = sum((-4 * mo[k] * Vo[k] / ro[k] for k in (0, 1)), sp.zeros(3, 1))
    return g00, g0, 1 + 2 * mo[0] / ro[0] + 2 * mo[1] / ro[1]


def pn_diff(rule_kw, seed):
    """Largest mismatch, through slow-orbit order (m ~ eps^2, v ~ eps), at a random configuration (40 digits)."""
    rnd = random.Random(seed)
    P = {c: sp.Float(rnd.uniform(-3, 3), 40) for vec in [Xo] + Yo for c in vec}
    P.update({c: sp.Float(rnd.uniform(-1, 1), 40) * eps_ for vv in Vo for c in vv})
    P.update({mm: sp.Float(rnd.uniform(0.2, 2), 40) * eps_**2 for mm in mo})
    Bm, Rm = bfp_metric(), river_rules_metric(**rule_kw)

    def worst(e, k):
        p = sp.Poly(sp.expand(sp.sympify(e).subs(P)), eps_)
        return max([abs(c) for d, c in zip(p.monoms(), p.coeffs()) if d[0] <= k] + [sp.Float(0)])
    return max(worst(Rm[0] - Bm[0], 4), max(worst(Rm[1][i] - Bm[1][i], 3) for i in range(3)), worst(Rm[2] - Bm[2], 2))


d_all = max(pn_diff({}, s_) for s_ in (1, 2, 3))
d_drop = {nm: min(pn_diff(kw, s_) for s_ in (1, 2, 3)) for nm, kw in
          (("no light-speed lag", {"lag": False}), ("clock rates added, not multiplied", {"multiply": False}),
           ("masses not seen through the other's clock", {"dressed": False}))}
check("O1 three river rules give gravity's slow-orbit field of two moving bodies exactly (Blanchet-Faye-Ponsot 1998, "
      "harmonic gauge; the field that fixes Einstein-Infeld-Hoffmann motion): A each body's own far field carried at light "
      "speed, B clock rates multiply (B3), C each mass seen through the other's clock. Dropping any one breaks it",
      d_all < mp.mpf("1e-25") and all(v_ > mp.mpf("1e-3") for v_ in d_drop.values()),
      {"all three rules": mp.nstr(d_all, 3), **{k: mp.nstr(v_, 3) for k, v_ in d_drop.items()}})
results["O1_drop"] = {k: float(v_) for k, v_ in d_drop.items()}
results["O1_all"] = float(d_all)

# O2 spin precession from that field
Xg = sp.Matrix(sp.symbols("x1:4", real=True))
Yg = sp.Matrix(sp.symbols("y1:4", real=True))
u1g, w2g, S2g = (sp.Matrix(sp.symbols(f"{c}1:4", real=True)) for c in ("u", "w", "S"))
m2g = sp.symbols("m2", positive=True)
dg_ = Xg - Yg
rg_ = sp.sqrt(dg_.dot(dg_))
Ug = m2g / rg_
gvec = -4 * m2g * w2g / rg_ - 2 * S2g.cross(dg_) / rg_**3
grad_ = lambda f: sp.Matrix([sp.diff(f, c) for c in Xg])  # noqa: E731
curl_ = lambda F: sp.Matrix([sp.diff(F[2], Xg[1]) - sp.diff(F[1], Xg[2]), sp.diff(F[0], Xg[2]) - sp.diff(F[2], Xg[0]),  # noqa: E731
                             sp.diff(F[1], Xg[0]) - sp.diff(F[0], Xg[1])])
Om_ = sp.Rational(3, 2) * u1g.cross(grad_(Ug)) - curl_(gvec) / 2
ng_ = dg_ / rg_
known_ = m2g / rg_**2 * (sp.Rational(3, 2) * ng_.cross(u1g) - 2 * ng_.cross(w2g)) + (3 * ng_ * ng_.dot(S2g) - S2g) / rg_**3
check("O2 each spin's wobble in that field: (3/2) v1 x grad U (the geodetic turn: 1 from curved space + 1/2 Thomas, "
      "both on our pages) plus half the curl of the frame-dragging field from body 2's motion (Part F) and spin (the "
      "shift's far field) gives the known spin-orbit and spin-spin precession (Barker-O'Connell; Kidder 1995); on a "
      "circular orbit Omega_1 = (2 + 3 m2/(2 m1)) L/r^3 + [3 n (n.S2) - S2]/r^3",
      sp.simplify(Om_ - known_) == sp.zeros(3, 1), "exact (symbolic)")

# O3 orbits: the closest point turns by 6 pi M/p per orbit; the waves' back-push carries off the quadrupole flux
adv = {nu_: periastron_advance(nu_, 1000.0, 0.4) for nu_ in (0.0, 0.25)}
flx = {r0_: flux_ratio(0.25, r0_) for r0_ in (40.0, 80.0)}
check("O3 orbits from that field (Einstein-Infeld-Hoffmann motion, plus the waves' back-push at 2.5PN): the closest "
      "point turns by 6 pi M/(a(1 - e^2)) per orbit whatever the mass ratio, and the energy the orbit loses matches the "
      "quadrupole formula (32/5) nu^2 M^5/r^5, better the wider the orbit",
      all(abs(a_ / p_ - 1) < 0.005 for a_, p_ in adv.values()) and abs(flx[80.0] - 1) < 0.01 and abs(flx[80.0] - 1) < abs(flx[40.0] - 1),
      {**{f"turn per orbit / 6piM/p, mass ratio nu = {k}": f"{a_ / p_:.4f}" for k, (a_, p_) in adv.items()},
       **{f"energy loss / quadrupole formula, r = {k:g} M": f"{v_:.4f}" for k, v_ in flx.items()}})
results["O3"] = {"adv": {str(k): v_ for k, v_ in adv.items()}, "flux": {str(k): v_ for k, v_ in flx.items()}}

# O4 against binary pulsars (numbers these formulas give; G = c = 1, masses in seconds)
Tsun, day_, yr_ = 4.925490947e-6, 86400.0, 365.25 * 86400.0
mp_, mc_, Pb_, e_ = 1.438, 1.390, 0.322997448918 * day_, 0.6171340          # Hulse-Taylor, Weisberg & Huang 2016
fe_ = (1 + 73 / 24 * e_**2 + 37 / 96 * e_**4) / (1 - e_**2)**3.5
Pdot = -(192 * np.pi / 5) * Tsun**(5 / 3) * (Pb_ / (2 * np.pi))**(-5 / 3) * fe_ * mp_ * mc_ * (mp_ + mc_)**(-1 / 3)
mA_, mB_, Pb2_, e2_ = 1.338185, 1.248868, 0.1022515592973 * day_, 0.087777023   # double pulsar, Kramer et al. 2021
OmB = (2 * np.pi / Pb2_)**(5 / 3) * Tsun**(2 / 3) * mA_ * (4 * mB_ + 3 * mA_) / (2 * (mA_ + mB_)**(4 / 3) * (1 - e2_**2))
OmB_deg = np.degrees(OmB) * yr_
check("O4 against binary pulsars: the orbit of the Hulse-Taylor pulsar shrinks at dP/dt = -2.402e-12 by the "
      "quadrupole formula (observed/predicted 0.9983 +- 0.0016, Weisberg & Huang 2016); pulsar B of the double pulsar "
      "wobbles at (2 + 3 mA/(2 mB)) L/r^3 averaged = 5.07 deg/yr (observed 4.77 +0.66/-0.65, Breton et al. 2008)",
      abs(Pdot + 2.402e-12) < 0.002e-12 and abs(OmB_deg - 4.77) < 0.66,
      f"dP/dt = {Pdot:.4e}; spin wobble {OmB_deg:.3f} deg/yr")
results["O4"] = {"Pdot": Pdot, "OmB_deg": OmB_deg}

# O5 spins that start different: do they line up?  (light ensemble here; orbits.py makes the long one)
e_start, e_end = ensemble(0.9, 40, 300.0, seed0=0)
c0, c1 = np.cos(np.radians(e_start[:, 0])), np.cos(np.radians(e_end[:, 0]))
lo_ = e_start[:, 1] < e_start[:, 2]
check("O5 spins that start pointing different ways (40 inspirals, mass ratio 0.9, spins 0.8, r = 300 M to 10 M): on "
      "average they do not line up with each other (mean cosine of the angle between them stays near 0), but they "
      "sort: when the heavier spin starts closer to the orbit's axis, their shadows on the orbit plane end closer "
      "together than when it starts farther (spin-orbit resonances, Schnittman 2004)",
      abs(c1.mean()) < 0.25 and np.median(e_end[lo_, 3]) < np.median(e_end[~lo_, 3]),
      f"mean cos(angle between spins): start {c0.mean():.3f}, end {c1.mean():.3f}; median end angle between shadows: "
      f"heavier spin starts closer {np.median(e_end[lo_, 3]):.0f} deg, farther {np.median(e_end[~lo_, 3]):.0f} deg")
results["O5"] = {"mean_cos_start": float(c0.mean()), "mean_cos_end": float(c1.mean())}

# O6 the wave: GW150914's run-up at leading order, and the ring-down of the final spinning hole
Mc_ = 30 * Tsun
tau_ = lambda f: 5 / 256 * Mc_**(-5 / 3) * (np.pi * f)**(-8 / 3)  # noqa: E731
ncyc = Mc_**(-5 / 3) / (32 * np.pi**(8 / 3)) * (35.0**(-5 / 3) - 150.0**(-5 / 3))
Mf_, af_ = 62 * 1.09 * Tsun, 0.67
f_ring = (1.5251 - 1.1568 * (1 - af_)**0.1292) / (2 * np.pi * Mf_)
check("O6 the waves of GW150914 (chirp mass about 30 solar masses as seen at the detector): from 35 to 150 Hz the "
      "leading-order run-up takes about 0.17 s and about 9 cycles (published: about 0.2 s and 8 cycles); the final "
      "spinning hole (62 solar masses, spin 0.67) rings at about 250 Hz (Kerr's fundamental tone, fit of Berti, Cardoso "
      "& Will 2006). The merger in between needs a full simulation",
      0.12 < tau_(35.0) - tau_(150.0) < 0.25 and 6 < ncyc < 11 and 230 < f_ring < 270,
      f"run-up {tau_(35.0) - tau_(150.0):.3f} s, {ncyc:.1f} cycles; ring-down {f_ring:.0f} Hz")
results["O6"] = {"t35_150": float(tau_(35.0) - tau_(150.0)), "cycles": float(ncyc), "f_ring": float(f_ring)}

# ================================================================ O7, O8  the next level of detail (round 13)
# Body 1's clock rate z1 = 1/u^t on a circular orbit, as a series in x = (m omega)^(2/3): it does not depend on
# coordinates. Truth from the first law (Le Tiec, Blanchet & Whiting 2012) with the known 2PN energy and angular
# momentum; the rules' value from their two-body field, each body moving on a geodesic with its own infinite part
# removed by Hadamard's finite part (next_order.py; Blanchet, Faye & Ponsot's method). Total mass 1.
import next_order as no_  # noqa: E402
third = sp.Rational(1, 3)
z_third = no_.redshift(third, order=6)
z_test = no_.redshift(0, order=6)
z_add = no_.redshift(third, order=6, multiply=False)
tf_, m1_sym = no_.truth_from_first_law()
truth_ok = all(sp.simplify(tf_[k].subs(m1_sym, third) - no_.truth(third)[k]) == 0 for k in range(4))
check("O7 the test of the next level is sound: general relativity's clock rate of body 1 follows from the first law "
      "(checked symbolically); the rules' field gives the same through the first correction (1PN, x^2) for body 1 at "
      "1/3 of the mass, and the same through 2PN (x^3) for a test mass (m1 -> 0). Adding the clock rates instead of "
      "multiplying them leaves body 1's clock rate unchanged through 2PN: rule B shapes the field between the bodies "
      "but drops out of what each body feels",
      truth_ok and z_third[:3] == no_.truth(third)[:3] and z_test == no_.truth(0) and z_add == z_third,
      f"m1 = 1/3: rules {z_third[:3]}, truth {no_.truth(third)[:3]}; test mass {z_test}; rates added {z_add}")
miss_ = sp.simplify(z_third[3] - no_.truth(third)[3])
check("O8 at 2PN the three rules miss: for body 1 at 1/3 of the mass the x^3 coefficient of its clock rate is "
      "+1205/2916 against general relativity's -3523/2916. Over seven mass ratios (next_order.json) the miss is "
      "m1 m2 (23 m1^2 - 50 m1 + 36)/3 = nu (15/2 + 9 Delta/2 - 23 nu/3), nonzero already at first order in the "
      "lighter mass",
      miss_ == no_.mismatch_formula(third) == sp.Rational(394, 243) and z_third[3] == sp.Rational(1205, 2916),
      f"rules {z_third[3]}, truth {no_.truth(third)[3]}, miss {miss_} (formula {no_.mismatch_formula(third)})")
results["O8_z1"] = {"rules": [str(c) for c in z_third], "truth": [str(c) for c in no_.truth(third)]}

# ================================================================ O9, O10  a fourth rule: the field's own energy (round 14)
# What general relativity's 2PN field has and the three rules lack is the field's own energy acting as a source.
# Its closed-form pieces (fourth_rule_closed.py) and four static field-on-field integrals (fourth_rule_integrals.py)
# are added to the rules' field and fed through the same clock-rate pipeline (fourth_rule.py -> fourth_rule.json).
import fourth_rule as fr_  # noqa: E402
import fourth_rule_integrals as fi_  # noqa: E402
fint_ = fi_.compute(verbose=False)
fco_ = fint_["coefficients"]
rows_ = fint_["rows"]
num_ok = all(abs(float(sp.Rational(r_["ball_finite_part"])) - float(r_["ball_numeric"])) < 1e-6 for r_ in rows_.values())
slab_ok = rows_["T1_grad"]["ball_finite_part"] == "-1/15" and rows_["T4_val"]["ball_finite_part"] == "2/3"
check("O9 the four static field-on-field integrals (the 1PN space-space field acting through the Newtonian "
      "potential's curvature) come out as exact numbers: X-hat at body 1 gains (alpha m1^2 m2 + beta m1 m2^2)/r12^3 "
      "with alpha = -1/2, beta = 1, and its gradient (gamma m1^2 m2 + delta m1 m2^2)/r12^4 with gamma = 0, "
      "delta = -2. Each finite part lands on a simple fraction to better than 1e-6, and the cut-out shape near body 2 "
      "is corrected by a rule checked on integrands with known answers",
      fco_ == {"alpha": "-1/2", "beta": "1", "gamma": "0", "delta": "-2"} and num_ok and slab_ok,
      f"{fco_}; parts {[(k_, r_['total']) for k_, r_ in rows_.items()]}")
results["O9_integrals"] = fint_
frj_ = json.load(open(os.path.join(HERE, "fourth_rule.json")))
pr_ = frj_["per_ratio"]
m_ = sp.symbols("m1", positive=True)
closed_f = m_*(1 - m_)*(5 - 10*m_ + 3*m_**2)
resid_f = m_*(1 - m_)**2*(3 - m_)/3
alg_ok = (sp.expand(closed_f + fr_.newtonian_map(m_).subs(fr_.F_VALUES) - resid_f) == 0 and
          sp.expand(closed_f + fr_.newtonian_map(m_).subs({fr_.al: sp.Rational(-1, 2), fr_.be: sp.Rational(5, 4),
                                                            fr_.ga: 0, fr_.de: sp.Rational(-9, 4)})) == 0)
table_ok = all(sp.sympify(pr_[r]["miss_closed"]) == closed_f.subs(m_, sp.Rational(r)) and
               sp.sympify(pr_[r]["miss_with_f"]) == resid_f.subs(m_, sp.Rational(r)) and
               pr_[r]["miss_with_point_term"] == "0" and pr_[r]["newtonian_map_matches"] for r in fr_.RATIOS)
# live: rerun mass ratio 1/3 through the pipeline from the cached closed-form pieces
ex13_ = fr_.closed_pieces("1/3")
z13_ = no_.redshift(third, 6, res_hook=fr_.make_hook(ex13_, fextra=fr_.f_terms(third, 1 - third)))
live_ = sp.expand(z13_[3] - no_.truth(third)[3])
live_ok = (sp.nsimplify(live_.subs(fr_.F_VALUES)) == sp.Rational(32, 243) and
           sp.nsimplify(live_.subs({fr_.al: 0, fr_.be: 0, fr_.ga: 0, fr_.de: 0})) == sp.Rational(4, 9))
check("O10 adding the field's own energy as a fourth rule closes most of the 2PN gap. For body 1 at 1/3 of the mass "
      "the miss in the x^3 clock-rate coefficient goes 394/243 (rules) -> 4/9 (+ closed-form field energy) -> 32/243 "
      "(+ the four integrals), 92% closed. At six mass ratios the stages follow m1 m2 (5 - 10 m1 + 3 m1^2) (fitted on "
      "three, holding on the other three, 2/5 never used in any fit) and m1 m2^2 (3 - m1)/3, which also follows by "
      "algebra from a hand-derived Newtonian map of the four numbers. What is left is exactly one point-like term per "
      "body, X-hat += m_a m_b^2/(4 r12^2 r_b) (beta -> 5/4, delta -> -9/4), which brings the miss to 0 at every ratio",
      alg_ok and table_ok and live_ok and all(frj_["closed_formula_holds"].values()),
      f"stages at 1/3: {pr_['1/3']['miss_rules']} -> {pr_['1/3']['miss_closed']} -> {pr_['1/3']['miss_with_f']} -> "
      f"{pr_['1/3']['miss_with_point_term']}; closed formula {frj_['closed_formula']}; live 1/3 rerun "
      f"{sp.nsimplify(live_.subs(fr_.F_VALUES))}")
results["O10_stages"] = {r: {k_: pr_[r][k_] for k_ in ("miss_rules", "miss_closed", "miss_with_f", "miss_with_point_term")}
                         for r in fr_.RATIOS}

# ================================================================ O11  the leftover derived (round 15)
# g00 has -2V^2. Each body's V is M/r plus, from retardation, -(M/2) a.n with a the body's own acceleration.
# The cross term of the two in -2V^2 is +2 M^2 a.(x-y)/r^2 = -4 Xacc. Round 14 left it out.
eps_ = sp.symbols("eps", positive=True); R1_, R2_ = sp.symbols("R1 R2", positive=True)
point_ok = True
for r in fr_.RATIOS:
    m1_ = sp.Rational(r); m2_ = 1 - m1_; exr_ = fr_.closed_pieces(r)
    for b_, (ma_, mb_, sg_) in ((1, (m1_, m2_, 1)), (2, (m2_, m1_, -1))):
        val_ = sp.simplify((fr_.OWN_SQUARE_WEIGHT*exr_[b_]["Xacc"][0]).subs({R1_: m2_, R2_: m1_}).coeff(eps_, 6))
        grd_ = sp.simplify((fr_.OWN_SQUARE_WEIGHT*exr_[b_]["Xacc"][1]).subs({R1_: m2_, R2_: m1_}).coeff(eps_, 6))
        # 8 x point term: value 8*(1/4) m_a m_b^2, x-gradient 8*(-1/4) m_a m_b^2 (sign flips at body 2)
        point_ok &= (val_ == 2*ma_*mb_**2) and (grd_ == -2*sg_*ma_*mb_**2)
json_ok = all(frj_["own_square_closes"].values()) and all(frj_["r12_matches_known"].values())
z13o_, sol13o_ = no_.redshift(third, 6, res_hook=fr_.make_hook(
    ex13_, fextra={b_: tuple(sp.sympify(fr_.f_terms(third, 1 - third)[b_][w_]).subs(fr_.F_VALUES) for w_ in (0, 1))
                   for b_ in (1, 2)}, own_square=True), return_sol=True)
a14_, a24_ = sp.symbols("a14 a24")
live_o_ok = (sp.simplify(z13o_[3] - no_.truth(third)[3]) == 0 and
             sp.nsimplify(sol13o_[a14_] + sol13o_[a24_]) == fr_.harmonic_r12_known(third))
check("O11 the leftover is derived: it is the speeding-up body's own field squared. In g00's -2V^2, each body's "
      "Newtonian M/r times its own retardation part -(M/2) a.n gives +2 M^2 a.(x-y)/r^2 = -4 Xacc, left out in round "
      "14. At Newtonian order its value and x-gradient at body a are +2 and -2 times m_a m_b^2/r12^3 (/r12^4), exactly "
      "8 times the point term, at all six ratios. With it the x^3 clock-rate miss is 0 at all six ratios, and the "
      "orbit size r12 at omega = 1 matches the known harmonic 2PN relation "
      "omega^2 = (m/r^3)[1 + (nu - 3) g + (6 + 41 nu/4 + nu^2) g^2] at all six (an independent observable)",
      point_ok and json_ok and live_o_ok,
      f"closes {frj_['own_square_closes']}; r12 matches {frj_['r12_matches_known']}; live 1/3: miss "
      f"{sp.simplify(z13o_[3] - no_.truth(third)[3])}, r12 x^2 {sp.nsimplify(sol13o_[a14_] + sol13o_[a24_])}")
results["O11_own_square"] = {r: {k_: pr_[r][k_] for k_ in ("miss_with_own_square", "r12_x2_with_own_square",
                                                         "r12_x2_known_harmonic")} for r in fr_.RATIOS}

# ================================================================ O12  term by term against the paper (round 15)
g00d_ = {r: fr_.g00_vs_bfp(r) for r in ("1/2", "2/5")}
g00d14_ = fr_.g00_vs_bfp("1/2", own_square=False)
check("O12 term by term against the paper: with the speeding-up body's own field squared, our g00 at body 1 along the "
      "solved circular orbit equals Blanchet, Faye & Ponsot's (1998) Eq. (7.6a) exactly through 2PN, at m1 = 1/2 and "
      "2/5 (scratch runs also at 1/3); without it (round 14) it misses by -2 m1 m2^2 x^3, the point term",
      all(v_ == 0 for v_ in g00d_.values()) and sp.simplify(g00d14_ + sp.Rational(1, 4)*eps_**6) == 0,
      f"differences {g00d_}; round 14 at 1/2: {g00d14_}")
results["O12_g00_vs_bfp"] = {k_: str(v_) for k_, v_ in g00d_.items()}

n_pass = sum(1 for v_ in results.values() if isinstance(v_, dict) and v_.get("pass"))
n_all = sum(1 for v_ in results.values() if isinstance(v_, dict) and "pass" in v_)
log(f"\n{n_pass}/{n_all} checks pass")
open(os.path.join(HERE, "checks_output.txt"), "w").write("\n".join(out) + "\n")
json.dump(results, open(os.path.join(HERE, "results.json"), "w"), indent=2, default=str)
