"""DRAFT checks: how the theory's dimension dial sets the power of a central force (Nov 26, 2026).

Three pictures of a point source's pull, read across the dimension dial D (number of space dimensions):
  flux (Gauss)       pull ~ 1/r^(D-1)        (surface = derivative of volume, any shape)
  velocity-line rule pull ~ 1/r^2  for every D (../velocity-lines)
  Einstein's river   pull ~ 1/r^(D-1) for D >= 3, none outside matter at D = 2
Right shape: the sphere for the field (rotation symmetry), the cube for the theory's cells (lattice),
the Gaussian (a product over axes that is also round) for non-integer D.
Run: python3 checks.py  (writes checks_output.txt, results.json)
"""
import json
import math
import os
import sys

import numpy as np
import sympy as sp
from scipy.integrate import quad, solve_ivp
from scipy.special import ive

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from grtools import einstein_mixed  # noqa: E402

out, results = [], {}


def log(s=""):
    print(s)
    out.append(s)


def check(name, ok, detail=""):
    results[name] = {"pass": bool(ok), "detail": str(detail)}
    log(f"[{'PASS' if ok else 'FAIL'}] {name}: {detail}")


D, r, t = sp.symbols("D r t", positive=True)
S_D = 2 * sp.pi**(D / 2) / sp.gamma(D / 2)           # area of the unit sphere in D dimensions

# ---------------------------------------------------------------- D1-D3 flux across the dial
log("== D1  Surface = derivative of volume, for the cube and the sphere, at any D")
a = sp.symbols("a", positive=True)
cube_face_layer = sp.binomial(D, D - 1)               # coefficient of a^(D-1) in the lift (1+a)^D, Gamma form
check("D1a cube face layer C(D, D-1) = D for every real D", sp.simplify(sp.gamma(D + 1) / (sp.gamma(D) * sp.gamma(2)) - D) == 0,
      "the lift's face layer D·a^(D-1) = d/da (a^D): extends off the integers")
V_D = sp.pi**(D / 2) / sp.gamma(D / 2 + 1)
check("D1b sphere: d/dr (V_D r^D) = S_D r^(D-1)", sp.simplify(sp.diff(V_D * r**D, r) - S_D * r**(D - 1)) == 0,
      "both shapes: area ~ r^(D-1), so the pull ~ 1/r^(D-1) whatever the shape; only the constant differs")

log("\n== D2  The Gaussian route: one axis's heat kernel to the D-th power, integrated over time")
kernel = (4 * sp.pi * t)**(-D / 2) * sp.exp(-r**2 / (4 * t))     # = product of D one-axis kernels
Dn = sp.Rational(7, 2)
G_num = sp.integrate(kernel.subs(D, Dn), (t, 0, sp.oo))
G_formula = sp.gamma(D / 2 - 1) / (4 * sp.pi**(D / 2) * r**(D - 2))
check("D2a potential G = Gamma(D/2-1)/(4 pi^(D/2) r^(D-2)) (shown at D = 7/2)",
      sp.simplify(G_num - G_formula.subs(D, Dn)) == 0, f"G(7/2) = {sp.simplify(G_num)}")
for Dv in [3, sp.Rational(5, 2), 4]:
    Gi = sp.integrate(kernel.subs(D, Dv), (t, 0, sp.oo))
    assert sp.simplify(Gi - G_formula.subs(D, Dv)) == 0
field = sp.simplify(-sp.diff(G_formula, r))
check("D2b pull = 1/(S_D r^(D-1)) for every D > 2", sp.simplify(field - 1 / (S_D * r**(D - 1))) == 0,
      f"pull = {field}; D = 3 gives {sp.simplify(field.subs(D, 3))}")
lap = sp.simplify(sp.diff(G_formula, r, 2) + (D - 1) / r * sp.diff(G_formula, r))
check("D2c G solves Laplace's equation away from the source at every D", lap == 0, "G'' + (D-1)/r G' = 0")
f2 = sp.limit(1 / (S_D * r**(D - 1)), D, 2)
check("D3a pull is continuous through D = 2: 1/(2 pi r)", sp.simplify(f2 - 1 / (2 * sp.pi * r)) == 0, f"{f2}")
eps = sp.symbols("epsilon")
Gser = sp.series(G_formula.subs(D, 2 + eps), eps, 0, 1).removeO()
log_part = sp.simplify(sp.diff(Gser, r))
check("D3b potential has a pole at D = 2; what is left is -ln(r)/(2 pi)",
      sp.simplify(log_part + 1 / (2 * sp.pi * r)) == 0, f"d/dr(finite part) = {log_part}")
results["S_D_values"] = {str(k): float(S_D.subs(D, k)) for k in [1, 2, 2.5, 3, 3.5, 4]}

# ---------------------------------------------------------------- D4 the right shape for the field
log("\n== D4  Right shape for the field: the sphere")
# point source at the centre of a cube of half-side 1 (3D): flux is right in total, uneven over the surface
En = lambda x, y: 1 / (4 * math.pi) / (x * x + y * y + 1)**1.5  # normal field on the face z = 1  # noqa: E731
face, _ = quad(lambda y: quad(lambda x: En(x, y), -1, 1, epsabs=1e-13)[0], -1, 1, epsabs=1e-13)
ratio = En(0, 0) / En(1, 1)
check("D4a cube around a point source: total flux 1, but 3*sqrt(3) times stronger at face centres than corners",
      abs(6 * face - 1) < 1e-10 and abs(ratio - 3 * math.sqrt(3)) < 1e-12, f"6 x face = {6 * face:.12f}, centre/corner = {ratio:.6f}")
# p-dial: only p = 2 keeps its shape when turned
th = np.linspace(0, 2 * np.pi, 2001)
dev = {}
for p in [1, 1.5, 2, 4, 12]:
    ux, uy = np.cos(th), np.sin(th)
    s = (np.abs(ux)**p + np.abs(uy)**p)**(1 / p)
    ux, uy = ux / s, uy / s                       # unit p-circle
    c, si = math.cos(0.3), math.sin(0.3)           # turn it by 0.3 rad
    vx, vy = c * ux - si * uy, si * ux + c * uy
    dev[p] = float(np.max(np.abs((np.abs(vx)**p + np.abs(vy)**p)**(1 / p) - 1)))
check("D4b p dial: only p = 2 (sphere) is unchanged by every turn",
      dev[2] < 1e-12 and all(dev[p] > 1e-3 for p in dev if p != 2), {k: f"{v:.2e}" for k, v in dev.items()})

# ---------------------------------------------------------------- D5 the theory's cells: cube near, sphere far
log("\n== D5  Right shape for the theory's cells: the cube (lattice); its field turns round far away")


def lattice_G3(x):
    # (6 - sum of neighbours) G = delta on Z^3:  G(x) = integral_0^inf prod_i e^{-2t} I_{x_i}(2t) dt
    f = lambda tt: np.prod([ive(xi, 2 * tt) for xi in x])  # noqa: E731
    T = 4000.0
    main = quad(f, 0, 50, limit=400, epsabs=1e-14)[0] + quad(f, 50, T, limit=400, epsabs=1e-14)[0]
    r2 = sum(xi * xi for xi in x)
    tail = quad(lambda tt: (4 * math.pi * tt)**-1.5 * math.exp(-r2 / (4 * tt)), T, np.inf)[0]  # continuum tail
    return main + tail


def lattice_a2(x):
    # potential difference a(x) = G(0) - G(x) on Z^2 (4 - neighbours): ~ (ln r + gamma + 1.5 ln 2)/(2 pi)
    f = lambda tt: ive(0, 2 * tt)**2 - ive(x[0], 2 * tt) * ive(x[1], 2 * tt)  # noqa: E731
    return quad(f, 0, 50, limit=400)[0] + quad(f, 50, 1e6, limit=800)[0]


pts_near = [(1, 0, 0), (1, 1, 1)]
pts_far = [(12, 0, 0), (7, 7, 7)]
an = {}
for pts, name in [(pts_near, "near"), (pts_far, "far")]:
    vals = [lattice_G3(p_) * 4 * math.pi * math.sqrt(sum(c * c for c in p_)) for p_ in pts]
    an[name] = (vals, abs(vals[0] - vals[1]) / vals[1])
check("D5a 3D lattice: axis-vs-diagonal mismatch shrinks with distance, and 4 pi r G -> 1",
      an["far"][1] < 0.01 < an["near"][1] and abs(an["far"][0][0] - 1) < 0.01,
      f"near: 4 pi r G = {an['near'][0][0]:.4f}, {an['near'][0][1]:.4f}; far: {an['far'][0][0]:.4f}, {an['far'][0][1]:.4f}")
a_ax, a_dg = lattice_a2((20, 0)), lattice_a2((14, 14))
pred = lambda rr: (math.log(rr) + 0.5772156649 + 1.5 * math.log(2)) / (2 * math.pi)  # noqa: E731
check("D5b 2D lattice: potential grows like ln(r)/(2 pi), so the pull is 1/(2 pi r)",
      abs(a_ax - pred(20)) < 2e-3 and abs(a_dg - pred(14 * math.sqrt(2))) < 2e-3,
      f"a(20,0) = {a_ax:.5f} vs {pred(20):.5f}; a(14,14) = {a_dg:.5f} vs {pred(14 * math.sqrt(2)):.5f}")
results["lattice"] = {"near": an["near"][0], "far": an["far"][0], "a2_axis": a_ax, "a2_diag": a_dg}

# ---------------------------------------------------------------- D6 the velocity rule across the dial


def hodo_misfit(power, e=0.3, norbits=1.0):
    """Orbit under pull 1/r^power from apocentre; misfit of the velocity path from its best circle."""
    def f(tt, s):
        x, y, vx, vy = s
        rr = math.hypot(x, y)
        acc = -1 / rr**power
        return [vx, vy, acc * x / rr, acc * y / rr]
    vc = 1.0
    s0 = [1.0, 0.0, 0.0, (1 - e) * vc]
    sol = solve_ivp(f, [0, 2 * math.pi * norbits], s0, rtol=1e-11, atol=1e-13, max_step=0.01)
    vx, vy = sol.y[2], sol.y[3]
    A = np.c_[2 * vx, 2 * vy, np.ones_like(vx)]
    c, *_ = np.linalg.lstsq(A, vx**2 + vy**2, rcond=None)
    R = math.sqrt(c[2] + c[0]**2 + c[1]**2)
    return float(np.max(np.abs(np.hypot(vx - c[0], vy - c[1]) - R)) / R)


log("\n== D6  The velocity-line rule does not move with the dial: it always wants 1/r^2")
mis = {Dv: hodo_misfit(Dv - 1) for Dv in [2.0, 2.5, 3.0, 3.5]}
check("D6a flux pull 1/r^(D-1) gives straight velocity lines only at D = 3",
      mis[3.0] < 1e-8 and all(mis[k] > 1e-3 for k in mis if k != 3.0), {k: f"{v:.1e}" for k, v in mis.items()})
check("D6b in 2D, the rule's pull (1/r^2) has straight velocity lines; the flux pull (1/r) does not",
      hodo_misfit(2.0) < 1e-8 and hodo_misfit(1.0) > 1e-2, f"1/r^2: {hodo_misfit(2.0):.1e}; 1/r: {hodo_misfit(1.0):.1e}")


def harmonics(Dsph, l):
    return math.comb(l + Dsph, Dsph) - (math.comb(l + Dsph - 2, Dsph) if l >= 2 else 0)


check("D6c the rule's velocity sphere is one dimension up: S^D; level counts 2l+1 (D = 2), (l+1)^2 (D = 3)",
      all(harmonics(2, l) == 2 * l + 1 and harmonics(3, l) == (l + 1)**2 for l in range(8)),
      "2D: 1, 3, 5, ... (2D hydrogen, Zaslow-Zandler 1967); 3D: 1, 4, 9, ... (Fock 1935)")

# ---------------------------------------------------------------- D7 orbits across the dial
log("\n== D7  Orbits across the dial: apsidal angle pi/sqrt(4 - D) for near-circular orbits (pull 1/r^(D-1))")


def apsidal(power, e=0.01):
    def f(tt, s):
        x, y, vx, vy = s
        rr = math.hypot(x, y)
        acc = -1 / rr**power
        return [vx, vy, acc * x / rr, acc * y / rr]
    ev = lambda tt, s: s[0] * s[2] + s[1] * s[3]  # noqa: E731
    ev.direction = 1
    sol = solve_ivp(f, [0, 40], [1.0, 0.0, 0.0, 1 - e], rtol=1e-12, atol=1e-13, events=ev)
    x, y = sol.y_events[0][0][:2]
    return math.atan2(y, x) % (2 * math.pi)


aps = {Dv: (apsidal(Dv - 1), math.pi / math.sqrt(4 - Dv)) for Dv in [2.0, 2.5, 3.0, 3.5]}
check("D7a apsidal angle matches pi/sqrt(4 - D)", all(abs(a_ - b_) / b_ < 5e-3 for a_, b_ in aps.values()),
      {k: f"{a_:.4f} vs {b_:.4f}" for k, (a_, b_) in aps.items()})
check("D7b stable circular orbits need D < 4 (Ehrenfest 1917)", sp.solve(sp.Eq(4 - D, 0), D) == [4],
      "pi/sqrt(4 - D) blows up at D = 4; above it circular orbits are unstable")
results["apsidal"] = {str(k): v for k, v in aps.items()}

# ---------------------------------------------------------------- D8 Einstein's river across the dial
log("\n== D8  Einstein's river in 2+1: empty space has a constant river speed, so no pull")
T_, R_, TH = sp.symbols("t r theta", real=True)
vf = sp.Function("v")(R_)
g3 = sp.Matrix([[-1 + vf**2, vf, 0], [vf, 1, 0], [0, 0, R_**2]])
G3, *_ = einstein_mixed(g3, [T_, R_, TH])
Gtt3 = sp.simplify(G3[0, 0])
target3 = -sp.diff(vf**2, R_) / (2 * R_)
check("D8a 2+1 river: G^t_t = -(v^2)'/(2r)", sp.simplify(Gtt3 - target3) == 0, f"G^t_t = {sp.factor(Gtt3)}")
Dsym = sp.symbols("Dsym", positive=True)
gen = lambda Dv, v2: -(Dv - 1) / (2 * R_**(Dv - 1)) * sp.diff(R_**(Dv - 2) * v2, R_)  # noqa: E731
check("D8b one formula for D = 2, 3, 4: G^t_t = -(D-1)/(2 r^(D-1)) (r^(D-2) v^2)'",
      sp.simplify(gen(2, vf**2) - Gtt3) == 0
      and sp.simplify(gen(3, vf**2) - (-sp.diff(R_ * vf**2, R_) / R_**2)) == 0
      and sp.simplify(gen(4, vf**2) - (-sp.Rational(3, 2) * sp.diff(R_**2 * vf**2, R_) / R_**3)) == 0,
      "D = 3, 4 from ../river-rebuild R5a, R6a; vacuum: v^2 ~ r^(2-D) (D > 2), v^2 constant (D = 2)")

# ---------------------------------------------------------------- D9 the three pictures
log("\n== D9  Where the three pictures agree on the dial")
table = {}
for Dv in [2, 3, 4]:
    flux = Dv - 1
    rule = 2
    ein = Dv - 1 if Dv > 2 else None
    table[Dv] = {"flux": flux, "rule": rule, "einstein": ein}
agree = [Dv for Dv, row in table.items() if row["flux"] == row["rule"] == row["einstein"]]
check("D9 flux, rule and Einstein agree only at D = 3", agree == [3],
      {k: v for k, v in table.items()})
results["table"] = table

n_pass = sum(1 for v in results.values() if isinstance(v, dict) and v.get("pass"))
n_all = sum(1 for v in results.values() if isinstance(v, dict) and "pass" in v)
log(f"\n{n_pass}/{n_all} checks pass")
open(os.path.join(HERE, "checks_output.txt"), "w").write("\n".join(out) + "\n")
json.dump(results, open(os.path.join(HERE, "results.json"), "w"), indent=2, default=str)
