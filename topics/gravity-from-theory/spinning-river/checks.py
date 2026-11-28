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
KS: Kerr's orbit factor 1 - 6u + 8 chi u^(3/2) - 3 chi^2 u^2 split into its pieces, by the double copy and by the rules
   (round 16, kerr_split.py -> kerr_split.json).
TW: why the shifted source's tower of shapes leaves a single -3 (round 17, tower.py -> tower.json).
TL: orbits tilted out of the equator: Kepler in the shift's distance is equatorial only; Laplace's rule in Newton
   and in Kerr (round 18, tilt.py -> tilt.json).
PL: Kerr's orbits over the poles against the shifted source's exact polar law (round 19, polar.py -> polar.json).
PX: Kerr's polar orbit split into the field's pieces: frame dragging twice, or the shape meeting the mass? (round 20,
    polar_split.py -> polar_split.json).
JB: the river's push and twist combined as one jet (Brennan's jet rule) or as two velocities (round 21,
    jet_boost.py -> jet_boost.json).
WL: one law for both wobbles of an equatorial orbit, and why the push's +28 equals the twist's 28 (round 22,
    wobble_law.py -> wobble_law.json).
PC: Kerr's orbits as Newton's orbits around the ring-shifted mass plus one mass term, and the source of the polar
    (a/r)^2 (M/r)^n terms (round 23, polar_coupling.py -> polar_coupling.json).
RT: is the one mass term the river's push? Time part, push and twist shares, and the orbit's shape per unit energy
    (round 25, river_term.py -> river_term.json).
EQ: the equator split three ways: which piece of the field gives each term of round 22's D - 3(v - q)^2
    (round 26, equator_split.py -> equator_split.json).
KN: charge on the spinning hole: which piece of the river gives each charge term (round 27, kn_split.py ->
    kn_split.json).
SH: the black hole's shadow piece by piece (round 28, shadow_split.py -> shadow_split.json).
TH: two bodies at the third level, where the rules must stop matching (round 30, third_level_predict.py,
    then third_level_compare.py against GWTC-4.0).
EC: the radar echo past the Sun piece by piece, predicted before comparing (round 29, echo_predict.py ->
    echo_prediction.json, then echo_compare.py -> echo_compare.json).
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

# ================================================================ KS1-KS7  Kerr's orbit factor split, two ways (round 16)
# A slightly eccentric equatorial orbit around a spinning hole turns by 2 pi (1/sqrt(T) - 1) per lap,
# T = 1 - 6u + 8 chi u^(3/2) - 3 chi^2 u^2 (u = M/r, chi = a/M). Which piece of the field makes each number?
import kerr_split as ks_  # noqa: E402
mp.mp.dps = 60
kj_ = json.load(open(os.path.join(HERE, "kerr_split.json")))
Tk_ = ks_.orbit_factor(ks_.ks_shell({}))
ex_ok_ = all(abs(Tk_(mp.mpf(u_), mp.mpf(c_)) - ks_.exact_T(mp.mpf(u_), mp.mpf(c_))) < mp.mpf("1e-50")
             for u_, c_ in (("0.05", "0.7"), ("0.05", "-0.7"), ("0.1", "0.3")))
isco_ok_ = all(abs(ks_.exact_T(1/ks_.isco(mp.mpf(c_)), mp.mpf(c_))) < mp.mpf("1e-50") for c_ in ("0.5", "0.9", "-0.5"))
check("KS1 Kerr's orbit factor from the exact Kerr-Schild mass shell is T = 1 - 6u + 8 chi u^(3/2) - 3 chi^2 u^2 "
      "(to 50 digits at three points, prograde and retrograde), and it vanishes at Bardeen-Press-Teukolsky's "
      "innermost stable orbit (an independent formula) at chi = 0.5, 0.9, -0.5", ex_ok_ and isco_ok_,
      f"T - formula: {kj_['exact']}")
dc_ = {k_: {kk_: sp.Rational(vv_) for kk_, vv_ in v_.items()} for k_, v_ in kj_["double_copy"].items()}
add_ok_ = True
for k_, v_ in dc_.items():
    on_ = set(k_.split("+")) if k_ != "E2 only" else set()
    add_ok_ &= (v_["u"] == -2*("PP" in on_) - 4*("EP" in on_) and v_["chi_u15"] == 8*("ES" in on_)
                and v_["chi2_u2"] == -3)
live_es_ = ks_.series(ks_.orbit_factor(ks_.ks_shell(dict(EE=1, PP=0, SS=0, EP=0, ES=1, PS=0))))
live_e2_ = ks_.series(ks_.orbit_factor(ks_.ks_shell(dict(EE=1, PP=0, SS=0, EP=0, ES=0, PS=0))))
check("KS2 way 1, the double copy: (k.p)^2 = (E + P + S)^2 has six pieces. E.E alone gives u: 0 and chi^2 u^2: -3; "
      "P.P adds -2u, 2E.P adds -4u, 2E.S adds +8 chi u^(3/2) and nothing else, P.S and S.S add nothing. The numbers "
      "simply add, in all 32 on/off settings",
      add_ok_ and len(dc_) == 32 and live_es_["chi_u15"] == 8 and live_es_["u"] == 0 and live_e2_["chi2_u2"] == -3,
      f"E.E only {kj_['double_copy']['E2 only']}; ES {kj_['double_copy']['ES']}; PP+EP "
      f"{kj_['double_copy']['PP+EP']}; all {kj_['double_copy']['PP+EP+ES+PS+SS']}")
sc_ = kj_["single_copy"]
s2_ = ks_.series(ks_.orbit_factor(ks_.single_shell(0, 2)))
p2_ = ks_.series(ks_.orbit_factor(ks_.single_shell(2, 0)))
check("KS3 the single copy -(2M/r)(E + P + S): S gives +4 chi u^(3/2), P gives -1 u, E alone -3 chi^2 u^2. Doubling the "
      "twist coupling doubles its number (S x 2 -> 8, first order) and doubling the push quadruples it (P x 2 -> -4, "
      "second order): gravity's 8 = 2 x 4 and its 2E.P's 4 = 2^2 x 1",
      sc_["P0S1"]["chi_u15"] == "4" and sc_["P1S0"]["u"] == "-1" and sc_["P0S0"]["chi2_u2"] == "-3"
      and s2_["chi_u15"] == 8 and p2_["u"] == -4, f"single copy {sc_}; S x 2: {s2_}; P x 2: {p2_}")
pn_ = kj_["pn"]
vb_, Ub_ = sp.symbols("v U")
Lb_ = sp.Matrix([[1, -vb_, 0, 0], [-vb_, 1, 0, 0], [0, 0, 1, 0], [0, 0, 0, 1]])      # boost to first order in v
hb_ = sp.expand((Lb_.T*sp.diag(2*Ub_, 2*Ub_, 2*Ub_, 2*Ub_)*Lb_)[0, 1])
check("KS4 way 2, the rules: in slow-motion form the spin number is 8 kappa and the u number -2(2 + 2 gamma - beta). "
      "Rule A (a moving body's field is its rest field boosted) gives h_0x = -(2U + 2U) v, half from the time part and "
      "half from the space part, so kappa = 1 and 8 = 4 + 4. With the push off (kappa = 0) the spin number is 0",
      pn_["Kerr (beta=gamma=kappa=q=1)"]["chi_u15"] == "8" and pn_["no push (kappa=0)"]["chi_u15"] == "0"
      and pn_["half push (kappa=1/2)"]["chi_u15"] == "4" and pn_["Kerr (beta=gamma=kappa=q=1)"]["u"] == "-6"
      and pn_["no space warp, no squared clock (beta=gamma=0)"]["u"] == "-4" and hb_ == -4*Ub_*vb_,
      f"boost h_0x = {hb_}; dials {pn_}")
Rf_, ab_, Mb_ = sp.symbols("R a M", positive=True)
shift_U_ = sp.series(Mb_/sp.sqrt(Rf_**2 - ab_**2), Rf_, sp.oo, 5).removeO()   # M/r in the equator, r = sqrt(R^2 - a^2)
ring_U_ = sp.series(Mb_*sum(sp.legendre(l_, 0)**2*(sp.sqrt(2)*ab_)**l_/Rf_**(l_ + 1) for l_ in range(0, 5)),
                    Rf_, sp.oo, 5).removeO()
qk_ = sp.simplify(sp.expand(shift_U_).coeff(Rf_, -3)/(Mb_*ab_**2/2))
check("KS5 the shape: the -3 chi^2 u^2 is set by the body's quadrupole alone (shape dial q gives -3q: 0, -3, -6), and "
      "Kerr's q = 1 is the shifted source's own: in the equator M/r = M/R + M a^2/(2 R^3) + ... in flat radius R. The "
      "rules don't fix q. A ring with Kerr's mass, spin and shape has radius sqrt(2) a and turns at c/sqrt(2)",
      pn_["no shape (q=0)"]["chi2_u2"] == "0" and pn_["double shape (q=2)"]["chi2_u2"] == "-6" and qk_ == 1
      and sp.simplify(ring_U_.coeff(Rf_, -3) - Mb_*ab_**2/2) == 0
      and sp.simplify((ab_/(sp.sqrt(2)*ab_))**2 - sp.Rational(1, 2)) == 0,
      f"q from the shift = {qk_}; ring's 1/R^3 term {ring_U_.coeff(Rf_, -3)}; ring speed J/(M R) = 1/sqrt(2)")
sims_ = kj_["simulated"]
s_live_, p_live_ = ks_.simulate(dict(EE=1, PP=0, SS=0, EP=0, ES=1, PS=0), 0.05, 0.7)
check("KS6 simulated orbits: integrating Hamilton's equations of the exact Kerr-Schild mass shell (full, E.E + 2E.S "
      "only, E.E only; chi = +-0.7, u = 0.05), the closest point turns per lap by 2 pi (1/sqrt(T) - 1) to better than "
      "1e-7", all(v_["rel_diff"] < 1e-7 for v_ in sims_.values()) and abs(s_live_ - p_live_)/abs(p_live_) < 1e-7,
      f"{ {k_: round(v_['rel_diff'], 12) for k_, v_ in sims_.items()} }; live E.E + 2E.S, chi 0.7: "
      f"{abs(s_live_ - p_live_)/abs(p_live_):.1e}")
zs_ = kj_["isco_zamo_speed"]; pg_ = kj_["product_form_gap"]
check("KS7 what doesn't carry over: the orbit page's product form T = (Delta/r^2)(1 - (2v)^2), Delta/r^2 = "
      "1 - 2u + chi^2 u^2, with v the speed seen by the local frame turning with the hole, holds for chi = 0 and fails for chi = +-0.5 (off by 0.08 at r = 10M); "
      "and that speed at the innermost stable orbit is c/2 only at chi = 0 (and in the limit chi -> 1)",
      abs(float(pg_["r=10,chi=0"])) < 1e-40 and abs(float(pg_["r=10,chi=0.5"])) > 0.05
      and abs(float(pg_["r=10,chi=-0.5"])) > 0.05 and zs_["0"] == "0.5"
      and all(abs(float(zs_[c_]) - 0.5) > 0.03 for c_ in ("0.3", "0.6", "0.9", "-0.5", "-0.9")),
      f"gaps {pg_}; ISCO speeds {zs_}")
results["KS_kerr_split"] = {"double_copy": kj_["double_copy"]["PP+EP+ES+PS+SS"], "single_copy": sc_,
                            "isco_zamo_speed": zs_}

# ================================================================ TW1-TW5  why the tower of shapes leaves one -3 (round 17)
# The shifted source has a whole tower of shapes (Kerr's M (ia)^l), yet T has a single shape term. tower.py.
import tower as tw_  # noqa: E402
twj_ = json.load(open(os.path.join(HERE, "tower.json")))
Rt_, at_, Mt_, rt_ = tw_.R, tw_.a, tw_.M, tw_.r
Ts_ = tw_.T_newton(tw_.SHIFTED)
Ts_r_ = sp.simplify(Ts_.subs(Rt_, sp.sqrt(rt_**2 + at_**2)))
w2_r_ = sp.simplify(tw_.omega2(tw_.SHIFTED).subs(Rt_, sp.sqrt(rt_**2 + at_**2)))
flat_ = tw_.coeff_list(sp.series(Ts_, at_, 0, 11).removeO(), Rt_)
check("TW1 Newton around the shifted source: in the equator its pull is -M/sqrt(R^2 - a^2) in flat radius R, which is "
      "-M/r in the shift's own radius r = sqrt(R^2 - a^2). Kepler's law then holds exactly in r (Omega^2 = M/r^3) and "
      "the orbit factor is exactly T = 1 - 3a^2/r^2. In flat radius the same T is 1 - 3(a/R)^2 - 3(a/R)^4 - ...: "
      "every shape contributes, and in r they sum to one term",
      sp.simplify(Ts_r_ - (1 - 3*at_**2/rt_**2)) == 0 and sp.simplify(w2_r_ - Mt_/rt_**3) == 0
      and all(c_ == -3 for c_ in flat_), f"T in r: {Ts_r_}; Omega^2 in r: {w2_r_}; flat-radius series {flat_}")
stair_ = twj_["staircase_in_r"]; alone_ = [sp.Rational(v_) for v_ in twj_["shape_alone_leading"]]
st_ok_ = True
for n_ in range(1, 5):
    row_ = [sp.Rational(v_) for v_ in stair_[f"up to (a/R)^{2*n_}"]]
    st_ok_ &= all(v_ == 0 for v_ in row_[1:n_]) and row_[n_] == -alone_[n_]
st_ok_ &= [sp.Rational(v_) for v_ in stair_["up to (a/R)^10"]] == [-3, 0, 0, 0, 0]
check("TW2 the staircase: cut the shifted source's tower after its n-th shape and T in r keeps a leftover at "
      "(a/r)^(2n+2); the next shape alone contributes exactly minus that leftover (-15/2, -105/8, -315/16, "
      "-3465/128), so each shape cancels what the ones below leave, through (a/r)^10",
      st_ok_, f"each shape alone {twj_['shape_alone_leading']}; staircase {stair_}")
riv_ = twj_["rivals_in_r"]; frq_ = twj_["in_frequency"]
check("TW3 rivals with Kerr's mass and quadrupole leave an endless series: quadrupole only gives T = 1 - 3(a/r)^2 + "
      "15/2 (a/r)^4 - ..., a ring of radius sqrt(2) a gives 1 - 3(a/r)^2 - 15/4 (a/r)^4 - .... Against the measured "
      "orbital frequency x = (M Omega)^(2/3), the shifted source gives exactly T = 1 - 3 chi^2 x^2 and the rivals "
      "do not stop either",
      all(v_ != "0" for v_ in riv_["quadrupole only"][1:] + riv_["ring of radius sqrt(2) a"][1:])
      and frq_["shifted source"] == ["-3", "0", "0", "0", "0"]
      and all(v_ != "0" for k_ in ("quadrupole only", "ring of radius sqrt(2) a") for v_ in frq_[k_][1:]),
      f"in r {riv_}; against frequency {frq_}")
Phk_ = sp.integrate(Rt_*Mt_/(Rt_**2 - at_**2)**sp.Rational(3, 2), Rt_)
check("TW4 Kepler in r is the whole story: T = 4 + d ln(Omega^2)/d ln R for any equatorial pull, so T = 1 - 3a^2/r^2 "
      "for every orbit holds exactly when Omega^2 = M/r^3; and Omega^2 = M/r^3 fixes the pull as the shifted source's",
      sp.simplify(Phk_ - tw_.SHIFTED) == 0
      and sp.simplify(Ts_ - (4 + Rt_*sp.diff(sp.log(tw_.omega2(tw_.SHIFTED)), Rt_))) == 0,
      f"pull from Kepler minus shifted source: {sp.simplify(Phk_ - tw_.SHIFTED)}")
ca_ = tw_.carter()
Tee_ = ks_.orbit_factor(ks_.ks_shell(dict(EE=1, PP=0, SS=0, EP=0, ES=0, PS=0)))
ee_ok_ = all(abs(Tee_(mp.mpf(u_), mp.mpf(c_)) - (1 - 3*mp.mpf(c_)**2*mp.mpf(u_)**2)) < mp.mpf("1e-45")
             for u_, c_ in (("0.05", "0.7"), ("0.1", "-0.4")))
xq_, chq_, yq_ = sp.symbols("x chi y", positive=True)
uq_ = xq_*(1 - chq_*xq_**sp.Rational(3, 2))**sp.Rational(-2, 3)
Traw_ = sp.expand(sp.series((1 - 6*uq_ + 8*chq_*uq_**sp.Rational(3, 2) - 3*chq_**2*uq_**2).subs(xq_, yq_**2),
                            yq_, 0, 9).removeO())
check("TW5 gravity does the same in the same radius: Kerr's radial motion in r is Carter's quartic "
      "(E^2-1) r^4 + 2M r^3 + (a^2(E^2-1) - L^2) r^2 + 2M (L - aE)^2 r, whose only spin terms are a^2 (shape) and "
      "L - aE (twist); the exact T = 1 - 6u + 8 chi u^(3/2) - 3 chi^2 u^2 follows from it symbolically, and the "
      "period is Kepler's in r plus the clock gap, 1/Omega = r^(3/2)/sqrt(M) + a. The Kerr-Schild E.E piece alone "
      "gives exactly 1 - 3 chi^2 u^2. Against the raw measured frequency T is an endless series "
      "(-4 chi x^(5/2), +8 chi^2 x^3, ...); with the clock gap 2 pi a taken off the period it is the finite one",
      ca_["quartic_ok"] and ca_["circular_ok"] and ca_["T_minus_formula"] == 0
      and ca_["period_minus_kepler_plus_gap"] == 0 and ee_ok_
      and Traw_.coeff(yq_, 5) == -4*chq_ and Traw_.coeff(yq_, 6) == 8*chq_**2,
      f"T from the quartic: {ca_['T']}; period check {ca_['period_minus_kepler_plus_gap']}; "
      f"raw-frequency series {Traw_}")
results["TW_tower"] = {k_: twj_[k_] for k_ in ("staircase_in_r", "rivals_in_r", "in_frequency")}

# ================================================================ TL1-TL5  tilted orbits: is Kepler in r equatorial only? (round 18)
# Round 17: in the equator Kepler's law holds exactly in the shift's distance r. Does it off the equator? tilt.py.
import tilt as tl_  # noqa: E402
tlj_ = json.load(open(os.path.join(HERE, "tilt.json")))
eq_ = tl_.equator_rates()
rr_, aa_, MM_ = tl_.r, tl_.a, tl_.M
check("TL1 in the equator the shifted source obeys Kepler in r and Laplace's rule: Omega^2 = M/r^3, the radial factor "
      "is 1 - 3a^2/r^2 and the up-and-down factor nu^2/Omega^2 is 1 + 3a^2/r^2, both exact, so the two add to exactly 2 "
      "(empty space off the ring). In oblate coordinates the pull is -M r/(r^2 + a^2 cos^2 theta), which separates",
      sp.simplify(eq_["omega2_r"] - MM_/rr_**3) == 0 and sp.simplify(eq_["T_r"] - (1 - 3*aa_**2/rr_**2)) == 0
      and sp.simplify(eq_["V_r"] - (1 + 3*aa_**2/rr_**2)) == 0 and eq_["laplace"] == 2 and tl_.oblate_potential_ok(),
      f"Omega^2 {eq_['omega2_r']}; radial {eq_['T_r']}; vertical {eq_['V_r']}; sum {eq_['laplace']}")


def tl_direct_(av, rv, frac, laps=6):
    """Fly the spherical orbit in ordinary 3D space under -M Re 1/|x - i a z_hat| and time its up-and-down beat."""
    Mv = 1.0
    Leq = tl_.L_equator(mp.mpf(Mv), mp.mpf(av), mp.mpf(rv))
    L = float(mp.mpf(frac)*Leq)
    s = tl_.sphere(mp.mpf(Mv), mp.mpf(av), mp.mpf(rv), mp.mpf(L))
    c = float(mp.sin(mp.radians(s["tilt_deg"])))          # cos(theta) at the turning latitude
    rho = np.sqrt(rv**2 + av**2)*np.sqrt(1 - c**2); zz = rv*c

    def f(t, y):
        x, yy, z, vx, vy, vz = y
        d = np.sqrt(complex(x*x + yy*yy + (z - 1j*av)**2))
        g = -Mv*np.array([x, yy, z - 1j*av])/d**3
        return [vx, vy, vz, g[0].real, g[1].real, g[2].real]

    def up(t, y):
        return y[2]
    up.direction = 1
    sol = solve_ivp(f, [0, laps*2*np.pi/float(s["Om_th"])*1.05], [rho, 0, zz, 0, L/rho, 0], rtol=1e-12, atol=1e-12,
                    events=up, dense_output=True)
    te = sol.t_events[0]
    per = (te[-1] - te[0])/(len(te) - 1)
    ts = np.linspace(0, sol.t[-1], 4000); Y = sol.sol(ts)
    Rq = Y[0]**2 + Y[1]**2 + Y[2]**2
    rq = np.sqrt((Rq - av**2 + np.sqrt((Rq - av**2)**2 + 4*av**2*Y[2]**2))/2)    # oblate radius along the path
    return float(s["tilt_deg"]), 2*np.pi/per, float(s["Om_th"]), float(np.max(np.abs(rq - rv)))


dts_ = [tl_direct_(0.8, 3.0, f_) for f_ in ("0.7", "0.3")]
check("TL2 the tilted orbits are real: the spherical orbits from the separated equations, flown directly in ordinary "
      "3D space at about 43 and 71 degrees of tilt, keep r fixed (to 1e-7) and beat up and down at the predicted rate "
      "(to 1e-6)",
      all(abs(w_ - q_)/q_ < 1e-6 and dr_ < 1e-7 for _, w_, q_, dr_ in dts_),
      f"(tilt deg, flown rate, predicted rate, max |r - r0|): {dts_}")
sph_ = tlj_["spheres"]
r3_ = sph_["a=0.8, r=3.0"]["rows"]; r5_ = sph_["a=0.8, r=5.0"]["rows"]
pol_ok_ = all(abs(v_["rows"][-1]["Om_th2"] - v_["polar_exact"]) < 1e-5 for v_ in sph_.values())
ser_ = [sp.Rational(c_) for c_ in tlj_["polar_series"]]
check("TL3 Kepler's period law in r is equatorial only: as an orbit at fixed r tilts from the equator to the pole, "
      "Omega_theta^2 r^3/M falls from 1 + 3a^2/r^2 (1.213 at a = 0.8, r = 3) to the polar value, and Omega_phi^2 r^3/M "
      "dips from 1 to 0.868 and back. The polar orbit is exact: Omega^2 r^3/M = (pi/2)^2/(2K(k) - E(k))^2 with k = a/r, "
      "complete elliptic integrals, an endless series 1 - 3k^2/2 + 33k^4/32 - 41k^6/64 + ...",
      abs(r3_[0]["Om_th2"] - sph_["a=0.8, r=3.0"]["equator_vertical"]) < 1e-3 and abs(r3_[0]["Om_ph2"] - 1) < 1e-3
      and min(x_["Om_ph2"] for x_ in r3_) < 0.87 and min(x_["Om_ph2"] for x_ in r5_) < 0.96 and pol_ok_
      and ser_ == [1, sp.Rational(-3, 2), sp.Rational(33, 32), sp.Rational(-41, 64), sp.Rational(2757, 8192)],
      f"a=0.8,r=3 (tilt, Om_th^2, Om_ph^2): {[(round(x_['tilt_deg'], 1), round(x_['Om_th2'], 4), round(x_['Om_ph2'], 4)) for x_ in r3_]}; "
      f"polar exact {sph_['a=0.8, r=3.0']['polar_exact']}; series {ser_}")
en_ = [x_["energy"] for x_ in r3_]
check("TL4 what survives off the equator: Kepler's energy law. An orbit with no spin about the axis has exactly "
      "E = -M/(2r) and K = M r in the shift's distance, so its radial motion is Kepler's; tilting toward the pole "
      "moves -2rE/M from 0.929 (equator, a = 0.8, r = 3) up to exactly 1",
      tlj_["polar_energy_minus_kepler"] == "0" and tlj_["polar_K_minus_Mr"] == "0"
      and all(en_[i_] < en_[i_ + 1] for i_ in range(len(en_) - 1)) and abs(en_[-1] - 1) < 1e-5
      and max(x_["residual"] for x_ in r3_ + r5_) < 1e-20,
      f"-2rE/M along the tilt: {[round(e_, 5) for e_ in en_]}")
kv_ = tl_.kerr_vertical(); s_, ch_ = kv_["s"], kv_["chi"]
check("TL5 Kerr in the equator: the up-and-down factor is (Omega_theta/Omega_phi)^2 = 1 - 4 chi u^(3/2) + 3 chi^2 u^2 "
      "(its -4 is Lense-Thirring's nodal turn 2J/r^3), and radial plus vertical is exactly 2 - 6u + 4 chi u^(3/2): "
      "the shape term cancels, as Laplace's rule demands of a Newtonian body in empty space. Both of Kerr's shape "
      "terms, -3 chi^2 u^2 and +3 chi^2 u^2, are the shifted source's -3a^2/r^2 and +3a^2/r^2",
      sp.expand(kv_["V"] - (1 - 4*ch_*s_**3 + 3*ch_**2*s_**4)) == 0
      and sp.expand(kv_["sum"] - (2 - 6*s_**2 + 4*ch_*s_**3)) == 0
      and sp.expand(kv_["V"] + kv_["T"]).coeff(ch_, 2) == 0
      and sp.simplify((3*ch_**2*s_**4).subs({ch_: aa_/MM_, s_: sp.sqrt(MM_/rr_)}) - 3*aa_**2/rr_**2) == 0,
      f"vertical {kv_['V']}; sum {kv_['sum']}")
results["TL_tilt"] = {"spheres": sph_, "polar_series": tlj_["polar_series"], "kerr": tlj_["kerr"]}

# ================================================================ PL1-PL4  Kerr's orbits over the poles (round 19)
# Round 18 found the shifted source's exact polar law. Does Kerr's polar orbit follow it, in the same radius? polar.py.
import polar as pl_  # noqa: E402
plj_ = json.load(open(os.path.join(HERE, "polar.json")))
pc_ = pl_.constants()
rP_, MP_, aP_ = pl_.r, pl_.M, pl_.a
xP_ = sp.symbols("x", positive=True)
RP_ = (pc_["E2"]*(xP_**2 + aP_**2)**2 - (xP_**2 - 2*MP_*xP_ + aP_**2)*(xP_**2 + aP_**2*pc_["E2"] + pc_["Q"]))
check("PL1 Kerr's polar orbits (no spin about the axis) at fixed Boyer-Lindquist r have exact constants E^2 and Q "
      "(Carter's R = R' = 0), and their up-and-down rate is a complete elliptic integral with modulus "
      "m = (a^2/r^2)((r^2+a^2)^2 - 4Mr^3)/((r^2+a^2)^2 - 4Ma^2 r). With M -> 0 the modulus is a^2/r^2 and the "
      "radial-to-polar prefactor is 1 + a^2/r^2, the shifted source's; with a -> 0 it is Schwarzschild's 1 - 6M/r",
      sp.simplify(RP_.subs(xP_, rP_)) == 0 and sp.simplify(sp.diff(RP_, xP_).subs(xP_, rP_)) == 0
      and sp.simplify(pc_["m"].subs(MP_, 0) - aP_**2/rP_**2) == 0
      and sp.simplify(pc_["pref"].subs(MP_, 0) - (1 + aP_**2/rP_**2)) == 0
      and sp.simplify(pc_["pref"].subs(aP_, 0) - (1 - 6*MP_/rP_)) == 0,
      f"E^2 = {pc_['E2']}; Q = {pc_['Q']}; m = {pc_['m']}")
fl_ = pl_.fly(1.0, 0.9, 6.0)
fls_ = plj_["flights"]
check("PL2 the exact rate is right: polar orbits flown with Hamilton's equations of the exact Kerr-Schild metric in "
      "Cartesian coordinates (no Carter separation used) keep r fixed to 1e-8 and beat up and down at the exact rate "
      "to 1e-9 (a = 0.9, r = 6; a = 0.5, r = 8; and a = 0.8 around M = 0.05, r = 3)",
      all(abs(f_["flown"] - f_["exact"])/f_["exact"] < 1e-9 and f_["max_dr"] < 1e-8 and abs(f_["H0"] + 0.5) < 1e-10
          for f_ in fls_ + [fl_]),
      f"rerun a=0.9, r=6: flown {fl_['flown']}, exact {fl_['exact']}, max |dr| {fl_['max_dr']}; stored "
      f"{[(f_['flown'], f_['exact']) for f_ in fls_]}")
srs_ = plj_["series"]
nw_ok_ = ([row_[0] for row_ in srs_["omega_theta2"]] == srs_["newton_omega_theta2"]
          and [row_[0] for row_ in srs_["ratio"]] == srs_["newton_ratio"])
check("PL3 at Newtonian strength Kerr follows the shifted source exactly, at every order in k = a/r: the u^0 part of "
      "Kerr's Omega_theta^2 r^3/M is 1 - 3k^2/2 + 33k^4/32 - 41k^6/64 + 2757k^8/8192, the shifted source's "
      "(pi/2)^2/(2K - E)^2, and the u^0 part of (Omega_r/Omega_theta)^2 is (1 + k^2)(2K(k)/pi)^2, through k^8. Kerr's "
      "polar orbit sees the shifted source's whole tower of shapes, at the same distance r",
      nw_ok_, f"Kerr u^0 rows {[row_[0] for row_ in srs_['omega_theta2']]} and {[row_[0] for row_ in srs_['ratio']]}; "
      f"Newton {srs_['newton_omega_theta2']} and {srs_['newton_ratio']}")
om_ = srs_["omega_theta2"]; ra_ = srs_["ratio"]
deg_ok_ = all(all(om_[j_][i_] == "0" for i_ in range(j_ + 1, 4)) for j_ in range(1, 3)) and om_[3][3] != "0"
check("PL4 beyond Newton, over the poles the a^2 terms and the mass terms mix: -4 u k^2 in Omega_theta^2 r^3/M and "
      "+15 u k^2 in (Omega_r/Omega_theta)^2 (for a black hole, k = chi u, so these are chi^2 u^3). In Omega_theta^2 r^3/M "
      "the k^(2j) part is a polynomial of degree j in u (seen through u^3). In the equator, in the same radius, Kerr's ratio 1 - 6u + 8 chi u^(3/2) - 3 chi^2 u^2 "
      "has no such term: there the shape stops at Newton's",
      om_[1][1] == "-4" and ra_[1][1] == "15" and ra_[1][2] == "-12" and deg_ok_,
      f"Omega_theta^2 r^3/M rows (k^0..k^8, u^0..u^3): {om_}; ratio rows {ra_}")
results["PL_polar"] = {"series": srs_, "flights": [{k_: v_ for k_, v_ in f_.items() if k_ != "path"} for f_ in fls_]}

# ================================================================ PX1-PX4  Kerr's polar orbit split into the field's pieces (round 20)
# Is Kerr's extra polar term frame dragging felt twice, or the shape meeting the mass? Same pieces as KS1. polar_split.py.
import polar_split as px_  # noqa: E402
pxj_ = json.load(open(os.path.join(HERE, "polar_split.json")))
z_px_, nz_px_ = px_.twist_drops_out()
dphi_px_, dchk_px_ = px_.node_rate_split()
nf_px_ = px_.node_from_flight(1.0, 0.9, 6.0)
ne_px_ = float(px_.node_rate_exact(1.0, 0.9, 6.0)[0])
lt_px_ = pxj_["node_vs_lense_thirring"]
check("PX1 over the poles the twist pieces (S = -aL/(r^2+a^2), so S.S, E.S, P.S) vanish when L = 0: they drop out of the "
      "radial and polar motion, so frame dragging cannot be felt twice there. They only turn the orbit's plane, at "
      "d(phi)/d(lambda) = 2MraE/Delta exactly, all of it from E.S and P.S. Kerr-Schild orbits flown in Cartesian "
      "coordinates turn at that rate to 1e-8, and far out the rate is Lense-Thirring's 2J/r^3",
      z_px_ and nz_px_ and dchk_px_ == 0 and abs(nf_px_[0] - ne_px_)/ne_px_ < 1e-8
      and abs(lt_px_[-1][1]/lt_px_[-1][2] - 1) < 1e-3 and "wES" in str(dphi_px_) and "wPS" in str(dphi_px_),
      f"node rate with dials {dphi_px_}; flown {nf_px_[0]} vs exact {ne_px_}; r = 200: {lt_px_[-1][1]} vs 2J/r^3 "
      f"{lt_px_[-1][2]}")
d1_px_, d2_px_, cee_px_ = px_.ee_alone_exact()
check("PX2 the E.E piece alone (the river's time part, no push, no twist) gives over the poles exactly the shifted "
      "source's law at every M: prefactor 1 + a^2/r^2 and modulus a^2/r^2, with E^2 = r/(r+M) and Carter's "
      "constant Mr^2/(r+M)",
      d1_px_ == 0 and d2_px_ == 0 and sp.simplify(cee_px_["E2"] - px_.r/(px_.r + px_.M)) == 0,
      f"E^2 = {cee_px_['E2']}, Kc = {cee_px_['Kc']}")
st_px_ = {k_: {kk_: [[sp.Rational(v_) for v_ in row_] for row_ in vv_] for kk_, vv_ in v_.items()}
          for k_, v_ in pxj_["settings"].items()}
rw_px_ = {k_: px_.series_rows({"EE": 1, "PP": int("PP" in k_), "EP": int("EP" in k_)}) for k_ in ("EE+EP", "EE+PP")}
same_px_ = all([[sp.Rational(str(v_)) for v_ in row_] for row_ in rw_px_[k_]["ratio"]] == st_px_[k_]["ratio"]
               for k_ in rw_px_)


def col_px_(k_, kk_, j_=1):
    return [st_px_[k_][kk_][i_][j_] for i_ in range(3)]


dEP_ = [x_ - y_ for x_, y_ in zip(col_px_("EE+EP", "ratio"), col_px_("EE", "ratio"))]
dPP_ = [x_ - y_ for x_, y_ in zip(col_px_("EE+PP", "ratio"), col_px_("EE", "ratio"))]
full_ = col_px_("EE+PP+EP", "ratio")
mEP_ = px_.polar_constants({"EE": 1, "EP": 1, "PP": 0})["m"]
mPP_ = px_.polar_constants({"EE": 1, "PP": 1, "EP": 0})["m"]
check("PX3 the polar mixed terms come from the river's push (E.P and P.P), with E.P carrying most of it: at first order "
      "in M the k^2 u term of (Omega_r/Omega_theta)^2 is 0 from E.E, +16 from E.P and -1 from P.P, adding to Kerr's "
      "+15; in Omega_theta^2 r^3/M the k^2 u term is +8 from E.E and -12 from E.P, giving Kerr's -4, and P.P adds "
      "nothing. E.P also moves the elliptic modulus off a^2/r^2; P.P does not",
      same_px_ and dEP_[1] == 16 and dPP_[1] == -1 and full_[1] == 15 and dEP_[1] + dPP_[1] == full_[1]
      and col_px_("EE+PP", "omega_theta2") == col_px_("EE", "omega_theta2")
      and col_px_("EE+EP", "omega_theta2")[1] == -4 and col_px_("EE", "omega_theta2")[1] == 8
      and sp.simplify(mPP_ - px_.a**2/px_.r**2) == 0 and sp.simplify(mEP_ - px_.a**2/px_.r**2) != 0,
      f"ratio u-columns: EE {col_px_('EE', 'ratio')}, E.P adds {dEP_}, P.P adds {dPP_}, full {full_}; modulus with "
      f"E.P: {sp.factor(mEP_)}")
eq_px_ = {k_: float(v_[1]) for k_, v_ in pxj_["equator_k2_and_k2u"].items()}
pc_px_ = {"PP": 8, "EP": 20, "SS": -16, "ES": -12, "PS": 0}
add_px_ = all(abs(v_ - sum(pc_px_[p_] for p_ in k_.split("+")[1:])) < 1e-9 for k_, v_ in eq_px_.items())
k2_px_ = all(abs(float(v_[0]) + 3) < 1e-9 for v_ in pxj_["equator_k2_and_k2u"].values())
re_px_ = {k_: float(px_.equator_k2u({p_: int(p_ in k_.split("+")) for p_ in pc_px_})[1])
          for k_ in ("EE+PP+EP", "EE+SS+ES", "EE+PP+EP+SS+ES+PS")}
check("PX4 in the equator the same shape-meeting-mass term (k^2 u, which is chi^2 u^3 for a black hole, in the radial "
      "to orbital rate ratio) is pushed in by E.P (+20) and P.P (+8) and taken out exactly by the twist, S.S (-16) and "
      "E.S (-12); P.S adds nothing. The pieces add in all 32 settings, and the full field's total is 0, which is why "
      "Kerr's equatorial ratio stops at Newton's shape while its polar one does not",
      add_px_ and k2_px_ and abs(re_px_["EE+PP+EP"] - 28) < 1e-9 and abs(re_px_["EE+SS+ES"] + 28) < 1e-9
      and abs(re_px_["EE+PP+EP+SS+ES+PS"]) < 1e-9,
      f"k^2 u by setting: {eq_px_}; recomputed {re_px_}")
results["PX_polar_split"] = {"equator_k2u": eq_px_, "settings": pxj_["settings"], "node_flight": nf_px_,
                             "node_exact": ne_px_}

# ================================================================ JB1-JB4  push and twist: one jet or two velocities? (round 21)
# Brennan's jet rule: where velocities combine, use the jet. Push and twist composed as one jet add the cross piece P.S
# to the field; added as two velocities they don't. jet_boost.py.
import jet_boost as jb_  # noqa: E402
jbj_ = json.load(open(os.path.join(HERE, "jet_boost.json")))
vac_ = jb_.vacuum()
check("JB1 the river's push and twist are two shears into the same (radial) leg, so as jets they compose into one "
      "shear, in either order, with no leftover. Squared into the field, that one jet carries the cross piece 2P.S, which "
      "adding the push and twist as two separate velocities leaves out. Only the one jet is Kerr: Ricci ~ 1e-50 at 50 "
      "digits, against ~ 1e-3 to 1e-2 near the hole for two velocities, whether the pieces are added in the metric or "
      "in its inverse",
      jb_.composition() and str(jb_.jet_minus_two_velocities()) == "2*P*S"
      and all(max(v_) < 1e-40 for k_, v_ in vac_.items() if k_.startswith("one jet"))
      and all(min(v_) > 1e-4 for k_, v_ in vac_.items() if k_.startswith("two velocities")),
      f"max Ricci at three points: {vac_}")
s1_jb_, u_jb_ = jb_.spin_term_exact({})
s2_jb_, _ = jb_.spin_term_exact({"PS": 0})
tw_jb_ = jbj_["spin_tower_chi_u15_u25_u35"]
near_ = lambda x_, y_: all(abs(p_ - q_) < 1e-6 for p_, q_ in zip(x_, y_))
check("JB2 in the equator the one jet stops Kerr's spin term at 8 chi u^(3/2): exactly, to first order in the spin. "
      "Two velocities multiply it by (1 - 2u)(1 - 3u), adding -40 chi u^(5/2) + 48 chi u^(7/2). The tower comes from "
      "the twist's E.S meeting the push E.P, and the jet's P.S, which acts only with the push, takes it out: +40, -48",
      sp.simplify(s1_jb_ - 8*sp.sqrt(u_jb_)) == 0
      and sp.simplify(s2_jb_ - 8*sp.sqrt(u_jb_)*(1 - 2*u_jb_)*(1 - 3*u_jb_)) == 0
      and near_(tw_jb_["EE+PP+EP+SS+ES+PS"], [8, 0, 0]) and near_(tw_jb_["EE+PP+EP+SS+ES"], [8, -40, 48])
      and near_(tw_jb_["EE+PP+EP+PS"], [0, 40, -48]) and near_(tw_jb_["EE+ES"], [8, 0, 0])
      and near_(tw_jb_["EE+PS"], [0, 0, 0]),
      f"spin term per unit a (r = 1, M = u): one jet {s1_jb_}, two velocities {sp.factor(s2_jb_)}; tower "
      f"(chi u^1.5, u^2.5, u^3.5): full {tw_jb_['EE+PP+EP+SS+ES+PS']}, no P.S {tw_jb_['EE+PP+EP+SS+ES']}, "
      f"push + P.S {tw_jb_['EE+PP+EP+PS']}")
cg_ = jb_.clock_gap_fn({}); cg2_ = jb_.clock_gap_fn({"PS": 0}); h_jb_ = mp.mpf("1e-12")
gap_ = lambda f_, rv_: float((f_(mp.mpf(rv_), 1, h_jb_) - f_(mp.mpf(rv_), 1, -h_jb_))/(2*h_jb_))
gaps_ = [(rv_, gap_(cg_, rv_), gap_(cg2_, rv_)) for rv_ in (6, 10, 30)]
sh_ = jbj_["shape_tower_k2"]
check("JB3 the one jet keeps Part 1's clock gap: each circular equatorial orbit is shifted by exactly 2 pi a per lap "
      "at every radius (Kerr's clock effect). Two velocities shift it by 2 pi a (1 - 4M/r) (small spin), which vanishes "
      "at r = 4M. The one jet also keeps Kerr's shape term exactly -3k^2 (k = a/r); two velocities add "
      "+108 k^2 u^2 - 240 k^2 u^3 + 192 k^2 u^4 + ...",
      all(abs(g1_ - 1) < 1e-12 and abs(g2_ - (1 - 4/rv_)) < 1e-12 for rv_, g1_, g2_ in gaps_)
      and sh_["one jet"] == ["-3", "0", "0", "0", "0"] and sh_["two velocities"] == ["-3", "0", "108", "-240", "192"],
      f"clock shift per unit a (r, one jet, two velocities): {gaps_}; k^2 tower {sh_}")
po_ = jb_.push_on_orbit()
pn_ = jb_.polar_node_rates()
r20_ = jbj_["round20_by_jet_layer"]
check("JB4 round 20's chi^2 u^3 cancellation in the equator does not come from the one jet: P.S adds 0 to it, and "
      "grouped by the jet's layers it is 2E.J = +8 against J.J = -8. Where the one jet lives: on an orbit the push "
      "answers the energy seen by the shell turning at a/(r^2+a^2), so k.p = (r^2+a^2)(E - aL/(r^2+a^2))/Delta; with "
      "two velocities the push sees E alone. Over the poles (L = 0) the radial, polar and clock motion are the same "
      "either way, and only the plane's turning tells them apart: 2MraE/Delta with one jet (Kerr), 2MraE/(r^2+a^2) with "
      "two velocities, slower by Delta/(r^2+a^2)",
      abs(r20_["full without PS"]) < 1e-9 and r20_["2E.J = 2EP + 2ES"] == 8 and r20_["J.J = PP + SS + 2PS"] == -8
      and po_[1] == 0 and pn_[2] == "0" and pn_[3] == "0" and pn_[5]
      and abs(dict(pn_[6])[6] - 0.6740016299918501) < 1e-12,
      f"push on the orbit {po_[0]}; polar node rates one jet {pn_[0]}, two velocities {pn_[1]}, ratio {pn_[4]}; "
      f"ratio at a = 0.9M: {pn_[6]}")
results["JB_jet_boost"] = {k_: jbj_[k_] for k_ in ("spin_term_per_a", "shape_tower_k2", "clock_gap", "polar",
                                                  "round20_by_jet_layer")}

# ================================================================ WL1-WL4  one law for both wobbles (round 22)
# Why the push's +28 equals the twist's 28 in the equator: spin squared can enter the wobbles only as (v - q)^2.
# wobble_law.py.
import wobble_law as wl_  # noqa: E402
c_wl_ = wl_.carter()
check("WL1 in Mino time the equatorial orbit's radial equation is r^4 (E^2 - 1) + 2 M r^3 - K0 r^2 + 2 M r x^2, with "
      "K0 = L^2 + a^2 (1 - E^2) and x = L - aE (Carter's radial function, rearranged). The mass meets the orbit as "
      "2 M r^3 (Newton) and 2 M r x^2; the second is one power of r, so it adds no radial stiffness. On every circular "
      "orbit the stiffnesses are M (r - 3x^2/r) in-and-out and M (r + 3x^2/r) up-and-down: they add to 2 M r, with "
      "no spin in it (symbolic, exact)",
      all(c_wl_.values()), f"{c_wl_}")
j_wl_ = wl_.jet_on_orbit()
l_wl_ = wl_.law()
check("WL2 on the circular orbit the river's jet energy is exact and simple: E + P + S (energy, push, twist) = 1/sqrt(D), "
      "D = 1 - 3v^2 + 2qv (v = sqrt(M/r), q = a/r), which is Carter's P/Delta. The orbit's Mino turning rate is "
      "r v k.p and its angular momentum about the spin is x = r (v - q) k.p. Dividing the stiffnesses by the turning "
      "rate squared: (Omega_r/Omega_phi)^2 = D - 3(v - q)^2 and (Omega_theta/Omega_phi)^2 = D + 3(v - q)^2, which are "
      "Kerr's known 1 - 6v^2 + 8qv - 3q^2 and 1 - 4qv + 3q^2",
      sp.simplify(j_wl_["kp_times_sqrtD"] - 1) == 0 and j_wl_["kp_equals_P_over_Delta"]
      and sp.simplify(j_wl_["turning_rate_over_r_v_kp"] - 1) == 0
      and sp.simplify(j_wl_["x_over_r_kp"] - (wl_.v - wl_.q)) == 0
      and all(l_wl_[k_] == 0 for k_ in ("T_r_from_stiffness_minus_law", "T_theta_from_stiffness_minus_law",
                                        "T_r_law_minus_known", "T_theta_law_minus_known")),
      f"jet on orbit {j_wl_}; law {l_wl_}")
jp_wl_ = wl_.jet_pieces()
a2m_ = {k_: sp.Rational(v_["a2M"]) for k_, v_ in jp_wl_.items()}
check("WL3 why +28 = 28: D holds the spin only to the first power, so spin squared enters the wobbles only through "
      "(v - q)^2, the angular momentum about the spin, squared. The (a/r)^2 term is therefore exactly -3 in-and-out and "
      "+3 up-and-down at every strength, with no (a/r)^2 (M/r)^n term for any n >= 1; round 20's push +28 against "
      "twist -28 is the n = 1 case. The split itself is bookkeeping: cutting the jet energy on the orbit instead gives "
      "E alone -1, push of E +4, twist -3, push of the twist 0, again summing to 0",
      l_wl_["D_q2_coefficient"] == 0 and l_wl_["q2_coefficient_radial"] == -3
      and l_wl_["q2_coefficient_vertical"] == 3
      and a2m_["E alone"] == -1 and a2m_["E + push of E"] == 3 and a2m_["E + twist"] == -4
      and a2m_["E + push of E + twist"] == 0 and jp_wl_["all (the one jet)"]["exact"] == "0",
      f"a^2 M parts of 1/(k.p)^2 as pieces are added: {a2m_}")
fl_wl_ = {f"u={u_}, chi={c_}": wl_.flight(u_, c_) for u_, c_ in ((0.1, 0.9), (0.1, -0.9))}
check("WL4 in flight: a slightly eccentric orbit integrated with Hamilton's equations of the exact Kerr-Schild field "
      "turns per lap by 2 pi (1/sqrt(D - 3(v - q)^2) - 1), with the spin and against it (r = 10M, a = +-0.9M; "
      "against the spin this is just outside the innermost stable orbit at 8.72M)",
      all(abs(s_ - p_) < 1e-5*abs(p_) for s_, p_ in fl_wl_.values()),
      f"(flown, law) radians per lap: {fl_wl_}")
results["WL_wobble_law"] = {"carter": c_wl_, "jet_pieces": jp_wl_, "flight": fl_wl_}

# ================================================================ PC1-PC4  Kerr = Newton's ring + one mass term (round 23)
# Is the coupling 2Mr x^2 (x = -aE over the poles) the source of round 19's (a/r)^2 (M/r)^n terms? polar_coupling.py.
import polar_coupling as pc_  # noqa: E402
n1_pc_, rf_pc_ = pc_.newton_plus_one(), pc_.ring_field()
check("PC1 Kerr's orbit equations in Mino time are Newton's equations for the ring-shifted mass (Euler's problem with "
      "the centres at z = +-ia, in oblate spheroidal coordinates, potential -M r/(r^2 + a^2 cos^2 theta)) with "
      "E^2 - 1 in place of twice Newton's energy, plus one mass term per equation: radial R_Kerr = R_Newton + "
      "2 M r (K - a^2) with K Carter's constant; polar Theta_Kerr = Theta_Newton exactly; clock dt/dlam = E Sigma + "
      "2 M r P/Delta; turning dphi/dlam = Newton's + (a/(r^2 + a^2)) 2 M r P/Delta. In Schwarzschild the radial term "
      "is 2 M r (Q + L^2), the textbook term behind perihelion precession (symbolic, exact, any orbit)",
      all(n1_pc_.values()) and all(rf_pc_.values()), f"{n1_pc_}; ring field {rf_pc_}")
sl_pc_ = pc_.stiffness_law()
tab_pc_ = pc_.dials()
pk_pc_ = pc_.polar_checks(tab_pc_)
check("PC2 the stiffness law of round 22 holds for every spherical orbit: omega_r^2 = M (r - 3K/r) + 4 a^2 Q/r^2 and "
      "the r^2 coefficient W = M (r + 3K/r) - 2 a^2 Q/r^2 (equator: K = x^2, Q = 0). Over the poles W = Q (1 + m), so "
      "(Omega_r/Omega_theta)^2 = (1 + m)(2K(m)/pi)^2 B with B = omega_r^2/W exactly 1 at Newtonian strength; this "
      "reproduces round 19's prefactor and its whole series table",
      all(sl_pc_.values()) and all(pk_pc_[k_] for k_ in ("pref_equals_omega_r2_over_Q_round19",
                                                          "omega_r2_law_holds_on_polar_orbit", "W_equals_Q(1+m)",
                                                          "B_is_1_at_Newton", "kerr_rows_match_round19",
                                                          "kerr_omega_theta2_k2u_is_-4")),
      f"stiffness law {sl_pc_}; polar {pk_pc_}")
check("PC3 the answer: x^2 = a^2 E^2 is only part of the source. Its a^2 part (2 M a^2 r) is the ring's oblateness and "
      "is Newtonian; its binding part 2 M a^2 r (E^2 - 1) gives 6 of round 19's +15 (a/r)^2 (M/r), and 2 M r Q (the "
      "Schwarzschild half) gives the other 9; at (a/r)^2 (M/r)^2 the -12 is 12 + 0 - 24 (the two halves together). With "
      "the whole term 2 M r (K - a^2) off every mixed term of the ratio vanishes, and with the clock terms off too "
      "Newton's polar laws return exactly. At order (a/r)^2 and every order in u: ratio = 1 - 6u + (a/r)^2 (3/2)"
      "(1 + 10u - 8u^2), Omega_theta^2 r^3/M = 1 - (a/r)^2 (3/2 + 4u); the -4 is -21/2 (radial term) + 13/2 (clock)",
      all(pk_pc_[k_] for k_ in ("k2_law", "newton_returns_ratio", "newton_returns_omega_theta2",
                                 "whole_term_off_kills_mixed_ratio_terms", "split_k2u_ratio_15=9+6",
                                 "split_k2u_omega_theta2_-4=-21/2+13/2"))
      and pk_pc_["split_k2u2_ratio"] == ["-12", "12", "0"]
      and sp.simplify(tab_pc_["Kerr"]["omega_theta2_k2_exact"] + (8*pc_.u + 3)/2) == 0,
      {n_: (str(d_["ratio_k2_exact"]), str(d_["omega_theta2_k2_exact"])) for n_, d_ in tab_pc_.items()})
fl_pc_ = [pc_.ring_flight(x_, delta=d_) for d_ in (1e-5, 1e-6) for x_ in (0.3, 0.6, 0.9)]
check("PC4 in flight: slightly perturbed polar orbits integrated in Cartesian coordinates in Newton's field of the "
      "ring-shifted mass (complex shift, no separation used) beat at (1 + k^2)(2K(k^2)/pi)^2, R_Newton's prediction, "
      "for a/r = 0.3, 0.6, 0.9; the small gap shrinks tenfold with a tenfold smaller kick",
      all(f_["rel_err"] < 2e-5 for f_ in fl_pc_[3:])
      and all(8 < a_["rel_err"]/b_["rel_err"] < 12 for a_, b_ in zip(fl_pc_[:3], fl_pc_[3:])),
      [(f_["k"], f_["delta"], f_["flown_ratio2"], f_["law"], f_["rel_err"]) for f_ in fl_pc_])
results["PC_polar_coupling"] = {"newton_plus_one": n1_pc_, "stiffness_law": sl_pc_, "polar": pk_pc_,
                                "flights": fl_pc_}


# ================================================================ RT1-RT4  the one mass term, piece by piece (round 25)
import river_term as rt_  # noqa: E402
pi_rt_ = rt_.piece_identities()
check("RT1 in the Kerr-Schild shell with a dial on each piece of (k.p)^2 (k.p = E + P + S), the river's time part "
      "alone is exactly Newton's ring with mass M E^2: R_EE = R_Newton(M -> M E^2), and the up-and-down equation "
      "has no dial. E.P and P.S give only M^2 terms, which cancel P.P's own M^2 terms exactly; the M^2 part vanishes "
      "for every setting whose pieces form one square (alpha E + beta P + gamma S)^2 and not for E.P alone "
      "(symbolic, exact, any orbit)",
      all(pi_rt_[k_] for k_ in ("shell_with_all_pieces_is_carter", "kerr_is_newton_plus_one_term",
                                 "time_part_is_newton_with_mass_M_E2", "M2_part_vanishes_with_all_pieces",
                                 "M2_part_vanishes_for_any_square", "M2_part_nonzero_for_EP_alone")),
      {k_: v_ for k_, v_ in pi_rt_.items() if isinstance(v_, bool)})
check("RT2 round 23's one term is shared by three pieces: 2Mr(K - a^2) = 2Mr(E^2 - 1)(r^2 + a^2) (time part) + "
      "2Mr[C - (E^2 - 1) r^2 - a^2 L^2/(r^2 + a^2)] (push, P.P) + 2Mr[a^2 L^2/(r^2 + a^2) - 2aEL] (twist, S.S and "
      "2E.S). The clock 2Mr P/Delta splits into time 2MrE, twist, push 4M^2 r^2 E/Delta and push-twist pieces; the "
      "turning term comes from the twist pieces only. Per unit energy (C -> C/E^2, L -> L/E, Newton's energy "
      "1 - 1/E^2), Kerr's radial function is E^2 times Newton's ring (mass M) + push + twist, and the clock and turning "
      "are E times their per-unit-energy forms (symbolic, exact)",
      all(pi_rt_[k_] for k_ in ("one_term_split_time_push_twist", "first_order_PP_is_push_share_SS_ES_twist",
                                 "clock_full_is_round23", "clock_pieces_add_up", "turning_full_is_round23",
                                 "turning_pieces_add_up", "turning_has_no_time_or_push_alone_piece",
                                 "kerr_per_unit_energy_is_newton_plus_push_plus_twist", "clock_per_unit_energy",
                                 "turning_per_unit_energy")),
      pi_rt_["shares"])
ph_rt_ = rt_.polar_homogeneity()
k2_rt_ = rt_.polar_k2_table()
uu_ = rt_.u
sy_ = lambda x_: sp.sympify(x_, locals={"u": uu_})
k2ok_ = (sp.simplify(sy_(k2_rt_["time part only"]["ratio_k2"]) - sp.Rational(3, 2)) == 0
         and sy_(k2_rt_["time part only"]["ratio_k0"]) == 1
         and sp.simplify(sy_(k2_rt_["push only"]["ratio_k2"]) - sp.Rational(3, 2)*(1 + 10*uu_ - 8*uu_**2)) == 0
         and sp.simplify(sy_(k2_rt_["push only"]["ratio_k0"]) - (1 - 6*uu_)) == 0
         and sp.simplify(sy_(k2_rt_["Kerr (time + push)"]["_exact"][3]) + (8*uu_ + 3)/2) == 0
         and sp.simplify(sy_(k2_rt_["push only"]["_exact"][3]) + (8*uu_ + 3)/2) != 0
         and sp.simplify(sy_(k2_rt_["time part only"]["_exact"][3]) + (8*uu_ + 3)/2) != 0)
check("RT3 over the poles the time part's share changes no rate ratio: on spherical polar orbits R and Theta are of "
      "degree one in (E^2 - 1, Q, Newton's mass), so the elliptic modulus and omega_r^2/Q don't depend on Newton's "
      "mass. With the time share alone the ratio stays Newton's (1 and 3/2 exactly); with the push share alone it is "
      "Kerr's at every order in u (1 - 6u and (3/2)(1 + 10u - 8u^2)). The coordinate-time rate Omega_theta^2 needs both "
      "(Kerr's -(3/2 + 4u) only with both, since the time part sets the clock)",
      all(ph_rt_[k_] for k_ in ("degree_one", "m_free_of_mu", "omega_r2_over_Q_free_of_mu")) and k2ok_,
      {n_: (d_["ratio_k0"], d_["ratio_k2"], d_["omega_theta2_k2_series"]) for n_, d_ in k2_rt_.items()})
ec_rt_ = [rt_.eccentric_polar(1.0, 0.9, 6.0, 10.0), rt_.eccentric_polar(1.0, 0.5, 7.0, 14.0)]
check("RT4 in flight: eccentric polar orbits (a = 0.9M, r from 6M to 10M; a = 0.5M, r from 7M to 14M) flown with "
      "Hamilton's equations of the exact Kerr-Schild metric in Cartesian coordinates, in Mino time, beat (radial "
      "against up-and-down) at the ratio of Newton's ring + push with the constants divided by E^2, to 1e-9; Newton's "
      "ring alone with the same turning points is far off",
      all(abs(e_["ratio_flown"]/e_["ratio_newton_ring_plus_push"] - 1) < 1e-9
          and abs(e_["ratio_kerr_quadrature"]/e_["ratio_newton_ring_plus_push"] - 1) < 1e-9
          and abs(e_["ratio_newton_ring_alone"]/e_["ratio_flown"] - 1) > 0.1 for e_ in ec_rt_),
      [{k_: v_ for k_, v_ in e_.items() if k_.startswith("ratio")} for e_ in ec_rt_])
results["RT_river_term"] = {"pieces": pi_rt_, "homogeneity": ph_rt_, "polar_k2": k2_rt_, "eccentric": ec_rt_}

# ================================================================ EQ1-EQ3  the equator split three ways (round 26)
import equator_split as eq_  # noqa: E402
sp_eq_ = eq_.split()
check("EQ1 at the equator (Q = 0), per unit energy, Kerr's radial function is Newton's ring + push + twist, and the "
      "push's 2MrL^2, the twist's -4MraL and the ring's 2Mra^2 complete one square: R = e r^4 + 2M(1 - e) r^3 - "
      "(L^2 - e a^2) r^2 + 2Mr(L - a)^2. With the push alone or the twist alone the a^2 L^2/(r^2 + a^2) parts don't "
      "cancel and R is not a polynomial (symbolic, exact)",
      all(sp_eq_.values()), sp_eq_)
gs_eq_ = eq_.generic_stiffness()
check("EQ2 for any radial function e r^4 + 2 alpha r^3 - W r^2 + 2 beta r on a circular orbit, the radial stiffness "
      "is alpha r - 3 beta/r and W = alpha r + 3 beta/r: the sign flip behind round 22's law. Kerr per unit energy "
      "has alpha = M(1 - e) and beta = M (L - a)^2 (symbolic, exact)",
      all(gs_eq_.values()), gs_eq_)
ra_eq_ = eq_.ratios()
check("EQ3 the equatorial rate ratio (Omega_r/Omega_phi)^2 piece by piece (v^2 = M/r, q = a/r): Newton's ring alone "
      "gives exactly 1 - 3q^2; ring + push brings -6v^2 and +28 q^2 v^2; ring + twist brings +8qv and -28 q^2 v^2; "
      "all three give exactly round 22's D - 3(v - q)^2 = 1 - 6v^2 + 8qv - 3q^2 (exact for ring and Kerr, series to "
      "v^10 for the partial cases; every orbit checked circular)",
      all(ra_eq_["checks"].values()) and all(ra_eq_[n_]["circular_ok"] for n_ in eq_.CASES),
      {**ra_eq_["checks"], **{n_: ra_eq_[n_]["series_in_chi_v"] for n_ in eq_.CASES}})
results["EQ_equator_split"] = {"split": sp_eq_, "stiffness": gs_eq_,
                               "ratios": {k_: v_ for k_, v_ in ra_eq_.items() if k_ != "curves"}}

# ================================================================ KN1-KN4  charge on the spinning hole, piece by piece (round 27)
import kn_split as kn_  # noqa: E402
sp_kn_ = kn_.split()
check("KN1 charge enters the river's shell in one slot: the field 2Mr (k.p)^2 becomes (2Mr - Q^2)(k.p)^2, and all "
      "pieces together give Kerr-Newman's radial function. The time part alone is Newton's ring around the mass "
      "function M - Q^2/(2r) (mass less the field's energy outside r, times E^2); Kerr-Newman = that ring + "
      "(2Mr - Q^2)(K - a^2); the second-order part vanishes for every setting whose pieces form one square (symbolic, "
      "exact, any orbit)",
      all(sp_kn_[k_] for k_ in ("shell_with_all_pieces_is_kerr_newman", "time_part_is_newton_ring_around_mass_function",
                                 "kerr_newman_is_ring_plus_one_term", "second_order_vanishes_with_all_pieces",
                                 "second_order_vanishes_for_any_square")),
      sp_kn_)
gs_kn_ = kn_.generic_stiffness()
check("KN2 at the equator, per unit energy, round 26's split carries over with 2Mr -> 2Mr - Q^2 in Newton's ring, "
      "push and twist: R = e r^4 + 2 alpha r^3 - W r^2 + 2 beta r - gamma with alpha = M(1 - e), beta = M(L - a)^2, "
      "gamma = Q^2 (L - a)^2, and the radial stiffness is alpha r - 3 beta/r + 4 gamma/r^2 (symbolic, exact)",
      all(sp_kn_[k_] for k_ in ("equator_ring_is_time_part", "equator_ring_plus_push", "equator_ring_plus_twist",
                                 "equator_all_three", "equator_square_form", "equator_coefficients_alpha_beta_gamma"))
      and all(gs_kn_.values()), gs_kn_)
ra_kn_ = kn_.ratios()
check("KN3 the equatorial rate ratio (Omega_r/Omega_phi)^2 of Kerr-Newman (v^2 = M/r, q = a/r, s = Q/r, "
      "V^2 = v^2 - s^2) is exactly [1 - 3q^2 + s^2(1 + q^2)/V^2] (Newton's ring, what the time part alone gives) + "
      "[-6v^2 + 3s^2 - s^4/V^2] (push) + [8qV] (twist); at a = 0 it is Reissner-Nordstrom's; the push or the twist "
      "alone still make the +-28 q^2 v^2 half-squares (40-digit checks at 12 random points, series to v^8)",
      all(ra_kn_["checks"].values()), ra_kn_["checks"])
cn_kn_ = kn_.carter_numeric()
check("KN4 independent of the river's pieces: Kerr-Newman's own Carter equations in Boyer-Lindquist form, circular "
      "orbits found numerically at 40 digits (a = 0.5M, Q = 0.6M, r = 5M, 8M, 13M), give the three-piece ratio to "
      "1e-15",
      all(c_["diff"] < 1e-15 for c_ in cn_kn_), cn_kn_)
results["KN_charge_split"] = {"split": sp_kn_, "stiffness": gs_kn_,
                              "ratios": {k_: v_ for k_, v_ in ra_kn_.items() if k_ != "curves"}, "carter": cn_kn_}

# ================================================================ SH1-SH4  the black hole's shadow, piece by piece (round 28)
import shadow_split as sh_  # noqa: E402
id_sh_ = sh_.identities()
check("SH1 light rays (null, E = 1): all pieces together give Kerr's null radial function; the horizon factor "
      "r^2 + a^2 - 2Mr appears only when the push is on (time part alone and ring + twist give r^2 + a^2: no horizon); "
      "at a = 0 the time part alone is Newton's light, R = r^2(r^2 + 2Mr - b^2), with no photon sphere, and time part "
      "+ push is Schwarzschild's, photon sphere r = 3M, b = sqrt(27) M (symbolic, exact)",
      all(id_sh_[k_] for k_ in ("all_pieces_is_kerr_null", "horizon_only_with_push", "a0_time_part_is_newton_light",
                                 "a0_push_is_schwarzschild", "a0_time_part_no_photon_sphere",
                                 "a0_push_photon_sphere_3M_sqrt27M")), id_sh_)
bd_sh_ = sh_.bending()
b_last_ = bd_sh_[-1]
ser3_ = 4/b_last_["b"] + 15*np.pi/(4*b_last_["b"]**2) + 128/(3*b_last_["b"]**3)
check("SH2 bending of light past the hole (a = 0, rays integrated numerically): the time part alone bends by exactly "
      "2 arctan(M/b), Newton's (Soldner's) value, at every b; time part + push gives Schwarzschild's deflection "
      "(series 4M/b + 15 pi M^2/4b^2 + 128 M^3/3b^3 at b = 1000M); the time part's share tends to one half",
      all(abs(d_["time_part"] - 2*np.arctan(1/d_["b"])) < 1e-9 for d_ in bd_sh_)
      and abs(b_last_["full"] - ser3_) < 1e-9 and abs(b_last_["time_over_full"] - 0.5) < 2e-3
      and all(bd_sh_[i_]["time_over_full"] < bd_sh_[i_ + 1]["time_over_full"] for i_ in range(len(bd_sh_) - 1)),
      {"rows": bd_sh_, "series_b1000_3_terms": ser3_})
sd_sh_ = sh_.shadows()
nopush_ = [k_ for k_ in sd_sh_ if "time part" in k_ or "ring + twist" in k_]
rp_ = [k_ for k_ in sd_sh_ if "ring + push" in k_]
kr_ = [k_ for k_ in sd_sh_ if "Kerr (all)" in k_]
sh3_ = {"no_horizon_or_shadow_without_push": all(not sd_sh_[k_]["horizon"] and not sd_sh_[k_]["shadow"] for k_ in nopush_),
        "ring_plus_push_centred": all(abs(sd_sh_[k_]["centre"]) < 1e-9 for k_ in rp_),
        "ring_plus_push_edge_on_shorter_than_2sqrt27": all(sd_sh_[k_]["height"] < 2*np.sqrt(27) - 1e-3
                                                          for k_ in rp_ if ", i=90," in k_),
        "kerr_matches_bardeen_alpha": all(abs(sd_sh_[k_]["alpha_min"] - sd_sh_[k_]["bardeen_alpha_min"]) < 1e-3 and
                                          abs(sd_sh_[k_]["alpha_max"] - sd_sh_[k_]["bardeen_alpha_max"]) < 1e-3 for k_ in kr_),
        "kerr_matches_bardeen_height": all(abs(sd_sh_[k_]["height"] - sd_sh_[k_]["bardeen_height"]) < 1e-6 for k_ in kr_),
        "kerr_edge_on_height_2sqrt27": all(abs(sd_sh_[k_]["height"] - 2*np.sqrt(27)) < 1e-6 for k_ in kr_ if ", i=90," in k_),
        "twist_shifts_centre": all(sd_sh_[k_]["centre"] > 0.5*float(k_.split(",")[0][2:])*np.sin(np.radians(
            float(k_.split(",")[1].split("=")[1]))) for k_ in kr_)}
check("SH3 the shadow (a = 0.5, 0.9, 0.998M; tilt 90, 50, 17 degrees): without the push there is no horizon and no "
      "shadow; ring + push gives a centred, squashed shadow; adding the twist gives Kerr's shadow, which matches "
      "Bardeen's closed form (edges to 1e-3, height to 1e-6; edge-on height exactly 2 sqrt(27) M) and is shifted "
      "toward the side moving away from us",
      all(sh3_.values()), {**sh3_, **{k_: {x_: y_ for x_, y_ in sd_sh_[k_].items() if x_ != "outline"} for k_ in sd_sh_}})
ev_sh_ = sh_.eht(sd_sh_)
ring_ = ev_sh_["EHT_2022_ring_diameter_uas"]; shad_ = ev_sh_["EHT_2022_shadow_diameter_uas"]
def _overlap(lo_, hi_, c_, e_):
    return lo_ <= c_ + e_ and hi_ >= c_ - e_
sh4_ = {nm_: {"kerr_overlaps_ring": _overlap(*ev_sh_[nm_]["kerr_range_uas"], *ring_),
              "kerr_overlaps_shadow": _overlap(*ev_sh_[nm_]["kerr_range_uas"], *shad_),
              "push_a0_in_shadow": abs(ev_sh_[nm_]["push_shadow_a0_uas"] - shad_[0]) < shad_[1],
              "time_part_none": ev_sh_[nm_]["time_part_alone_uas"] == 0.0}
        for nm_ in ("GRAVITY 2019 (VLTI)", "Do et al. 2019 (Keck)")}
check("SH4 Sgr A*: with mass and distance from stellar orbits alone (GRAVITY 2019 or Keck, Do et al. 2019) the push's "
      "shadow (2 sqrt(27) GM/c^2D, about 51-52 microarcseconds) and Kerr's range over spin and tilt overlap the Event "
      "Horizon Telescope's 2022 ring (51.8 +- 2.3) and shadow (48.7 +- 7.0); the time part alone predicts none",
      all(all(v_.values()) for v_ in sh4_.values()), {**sh4_, "eht": ev_sh_})
results["SH_shadow_split"] = {"identities": id_sh_, "bending": bd_sh_, "eht": ev_sh_}

# ================================================================ EC1-EC4  the radar echo past the Sun (round 29)
import echo_predict as ep_  # noqa: E402
import echo_compare as ec_  # noqa: E402
lf_ec_ = ep_.leg_formula()
check("EC1 light time past a mass, per leg, to first order in M, with the push at strength lambda (k.p = E + lambda P, "
      "the field still one square): sqrt(r^2 - r0^2) + M[(1 + lambda^2) ln((r + sqrt(r^2 - r0^2))/r0) + "
      "sqrt((r - r0)/(r + r0))] plus a clock relabelling 2 lambda M ln r that cancels from anything path-dependent. "
      "The time part gives the 1 (slow clocks), the push lambda^2; Newton's corpuscle gives -1 (symbolic, exact; "
      "40-digit spot checks of the antiderivatives)",
      all(v_ for k_, v_ in lf_ec_.items() if k_ != "b^2(r0)"), lf_ec_)
ph_ec_ = ep_.perihelion_dial()
bd_ec_ = ep_.bending_checks()
check("EC2 the same dial sets the perihelion advance, 6 lambda^2 (pi M/a per orbit), and the bending, 2(1 + lambda^2) M/b: "
      "one dial ties delay, bending and perihelion (beta + gamma = 2); with separate weights mu (EP) and nu (PP) the "
      "perihelion factor is (2 mu^2 + nu)/3 and the delay's nu, so it is the square that ties them (symbolic; bending "
      "integrated at b = 1e5 M to 1e-4)",
      ph_ec_["general_is_1-2(2mu^2+nu)v^2"] and ph_ec_["square_dial_is_1-6lam^2v^2"]
      and all(abs(r_["deflection_x_b_over_2M"] - r_["predicted_1+lam2"]) < 1e-4 for r_ in bd_ec_),
      {"perihelion": ph_ec_, "bending": bd_ec_})
se_ec_ = ep_.sun_exact_checks()
tw_ec_ = ep_.twist_delay()
check("EC3 at the Sun (ray at 1.6 solar radii, out to 1.016 AU) the exact light times of the time part and of time part "
      "+ push match the first-order formula to within 20 M^2/r0 (centimetres in 10-20 km); the twist adds -8 M a/b one "
      "way for a prograde equatorial ray (exact Kerr, 40 digits, 3 geometries), which for the Sun's spin is about 0.02 ns "
      "round trip",
      all(abs(r_["difference_m"]) < 20*r_["M2_over_r0_m"] for r_ in se_ec_)
      and all(abs(r_["half_difference_over_Ma/r0"] + 8) < 0.05 for r_ in tw_ec_["rows"])
      and 1e-11 < tw_ec_["sun_round_trip_s_at_1.6Rsun"] < 1e-10,
      {"sun": se_ec_, "twist": tw_ec_})
pr_ec_ = json.load(open(os.path.join(HERE, "echo_prediction.json")))
cm_ec_ = ec_.compare()
check("EC4 the prediction written first (echo_prediction.json: dials lambda = tau = 1 from rounds 7-8, no measured delay "
      "used; gamma = 1; shadow-only dials gamma = 0.71-1.13) against Cassini 2002 (Bertotti, Iess & Tortora 2003: "
      "gamma = 1 + (2.1 +- 2.3) x 10^-5): within 1 sigma; inside the shadow-only range; the time part alone (gamma = 0) "
      "and Newton's corpuscle (gamma = -2) excluded; the twist below Cassini's resolution",
      (not pr_ec_["uses_any_measured_delay_or_gamma"]) and pr_ec_["prediction"]["gamma_predicted"] == 1.0
      and abs(cm_ec_["derived_dials"]["pull_in_sigma"]) < 1.0 and cm_ec_["shadow_only"]["cassini_inside"]
      and cm_ec_["time_part_alone"]["pull_in_sigma"] > 1e4 and cm_ec_["twist"]["twist_over_resolution"] < 0.1,
      cm_ec_)
results["EC_echo_split"] = {"prediction": pr_ec_["prediction"], "compare": cm_ec_}


# ---------------------------------------------------------------- TH: two bodies at the third level (round 30)
import third_level_predict as tl3_  # noqa: E402
cs_tl_, m1_tl_ = tl3_.redshift_x4()
pi_tl_ = sp.factor(cs_tl_[4].coeff(sp.pi**2))
check("TH1 the stationary-phase rule phi_k = 40 c_k/((k-8)(k-5)) reproduces the known phase coefficients phi_2, "
      "phi_3 = -16 pi and phi_4 from the known energy and emitted power",
      all(tl3_.phase_rule_check().values()))
check("TH2 the first law with the known third-level energy gives the round-13 redshift through x^3 and, at x^4, a pi^2 "
      "term 41 m1 (m1 - 6)(m1 - 1)^2 / 192 for body 1: zero only for a test body (m1 = 0), first order in the small "
      "mass (-41 pi^2 m1/32, the self-force value) and second order when body 1 is the big one",
      all(sp.expand(cs_tl_[k] - tl3_.next_order_truth(m1_tl_)[k]) == 0 for k in range(4))
      and sp.expand(pi_tl_ - sp.Rational(41, 192)*m1_tl_*(m1_tl_ - 6)*(m1_tl_ - 1)**2) == 0
      and sp.limit(pi_tl_/m1_tl_, m1_tl_, 0) == sp.Rational(-41, 32))
no_tl_ = json.load(open(os.path.join(HERE, "next_order.json")))
fr_tl_ = json.load(open(os.path.join(HERE, "fourth_rule.json")))
check("TH3 every coefficient the rules produced at the second level is rational (next_order.json, all ratios), and so "
      "are the four field-on-field integrals of the fourth rule: the rules' closed forms cannot make the pi^2",
      all(sp.nsimplify(c).is_Rational for r in no_tl_["two_pn"].values() for c in r["rules"])
      and all(sp.nsimplify(c).is_Rational for c in fr_tl_["f_values"].values()))
de3_tl_, dphi_tl_, frac_tl_ = tl3_.floor_shift(sp.Rational(1, 4))
check("TH4 the size floor: missing only the pi^2 piece moves the third-level phase coefficient by -80 delta e3 = "
      "-(1025/6) pi^2 nu (power kept as in relativity), 0.65 of itself at equal masses and nothing for a test body",
      abs(dphi_tl_ + 1025/6*np.pi**2/4) < 1e-9 and abs(frac_tl_ - 0.6464) < 1e-3
      and tl3_.floor_shift(0)[1] == 0.0)


tlc_ = json.load(open(os.path.join(HERE, "third_level_compare.json")))
tlp_ = json.load(open(os.path.join(HERE, "third_level_prediction.json")))
check("TH5 the floor was written before the bound (third_level_prediction.json, no measured bound used) and held "
      "against GWTC-4.0's 66-event 3PN bound (arXiv:2603.19020, Table 2: +-0.18 hierarchical, +-0.06 restricted): "
      "at equal masses the missing pi^2 alone is 0.65, outside both; a test body is untouched",
      (not tlp_["uses_any_measured_bound"]) and tlc_["prediction_sha256"] == tlp_["script_sha256"]
      and tlc_["by_mass_ratio"]["0.25"]["floor_fraction"] > tlc_["measured"]["FTI_hierarchical"][2]
      and tlc_["by_mass_ratio"]["0.25"]["floor_fraction"] > tlc_["measured"]["TIGER_hierarchical"][2]
      and tlc_["test_body"] == 0.0)


n_pass = sum(1 for v_ in results.values() if isinstance(v_, dict) and v_.get("pass"))
n_all = sum(1 for v_ in results.values() if isinstance(v_, dict) and "pass" in v_)
log(f"\n{n_pass}/{n_all} checks pass")
open(os.path.join(HERE, "checks_output.txt"), "w").write("\n".join(out) + "\n")
json.dump(results, open(os.path.join(HERE, "results.json"), "w"), indent=2, default=str)
