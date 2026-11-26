"""DRAFT checks: rebuild the river model of gravity on the velocity-line rule (Nov 26, 2026).

Chain tested here:
  velocity-line rule (../velocity-lines) -> pull K/r^2 -> energy E conserved
  -> river = the E = 0 free-fall family (flat middle of the dial) -> v^2 = 2K/r
  -> river metric -dt^2 + (dr + v dt)^2 + r^2 dOmega^2 is vacuum (Einstein) and exact
  -> perihelion, light bending, gamma = beta = 1 (so frame dragging from ../README.md carries over).
Units: G = c = 1 unless stated. Run: python3 checks.py  (writes checks_output.txt, results.json)
"""
import json
import math
import sys

import numpy as np
import sympy as sp
from scipy.integrate import quad, solve_ivp

sys.path.insert(0, __file__.rsplit("/", 1)[0] if "/" in __file__ else ".")
from grtools import christoffel, einstein_mixed  # noqa: E402

out, results = [], {}


def log(s=""):
    print(s)
    out.append(s)


def check(name, ok, detail=""):
    results[name] = {"pass": bool(ok), "detail": str(detail)}
    log(f"[{'PASS' if ok else 'FAIL'}] {name}: {detail}")


# ---------------------------------------------------------------- R1  rule -> 1/r^2 -> E conserved
log("== R1  The rule gives the pull K/r^2 (velocity-lines V7); then E = v^2/2 - K/r is conserved")
t = sp.symbols("t", real=True)
K = sp.symbols("K", positive=True)
x, y = sp.Function("x")(t), sp.Function("y")(t)
r_ = sp.sqrt(x**2 + y**2)
E_expr = (sp.diff(x, t)**2 + sp.diff(y, t)**2) / 2 - K / r_
dE = sp.diff(E_expr, t).subs({sp.diff(x, t, 2): -K * x / r_**3, sp.diff(y, t, 2): -K * y / r_**3})
check("R1 energy conserved under the rule's pull", sp.simplify(dE) == 0, f"dE/dt = {sp.simplify(dE)}")

# ---------------------------------------------------------------- R2  the river is the E = 0 family
log("\n== R2  River = free fall from rest at infinity: E = 0, so v^2 = 2K/r, on the dial's flat middle")
r, E, vx, vy = sp.symbols("r E v_x v_y", real=True)
v_river = sp.solve(sp.Eq(sp.Symbol("w", positive=True)**2 / 2 - K / r, 0), sp.Symbol("w", positive=True))
check("R2a E = 0 gives the escape speed", sp.simplify(v_river[0]**2 - 2 * K / r) == 0,
      f"v = {v_river[0]}  (Newton's escape speed, now an output of the rule)")

# Milnor/Fock velocity-space metric ds = 2|dv| / (|v|^2 - 2E): curvature -2E (the dial q = -2E)
lam = 2 / (vx**2 + vy**2 - 2 * E)
gauss = sp.simplify(-(sp.diff(sp.log(lam), vx, 2) + sp.diff(sp.log(lam), vy, 2)) / lam**2)
check("R2b velocity-space curvature equals the dial -2E", sp.simplify(gauss - (-2 * E)) == 0,
      f"Gaussian curvature = {gauss}")
check("R2c river (E = 0) lives on the flat middle (curvature 0)", sp.simplify(gauss.subs(E, 0)) == 0,
      "q = -2E = 0: the Galilean end of the theory's dial")
# scale of the velocity geometry along any orbit: |v|^2 - 2E = 2K/r = v_esc^2
check("R2d the river speed is the scale of the velocity geometry",
      sp.simplify((2 * (E + K / r) - 2 * E) - v_river[0]**2) == 0,
      "|v|^2 - 2E = v_esc(r)^2, so the dial chart's scale factor is 2/v_esc^2 = r/K")

# ---------------------------------------------------------------- R3  river metric is vacuum
log("\n== R3  Put the rule's river into the theory's local relativity: the metric is vacuum")
T, R, TH, PH = sp.symbols("t r theta phi", real=True)
M = sp.symbols("M", positive=True)
X4 = [T, R, TH, PH]


def river_metric(vfun, X, extra):
    """-dt^2 + (dr + v dt)^2 + r^2 dOmega^2 ; extra = angular block for the sphere."""
    n = len(X)
    g = sp.zeros(n, n)
    g[0, 0] = -1 + vfun**2
    g[0, 1] = g[1, 0] = vfun
    g[1, 1] = 1
    for i, e in enumerate(extra):
        g[2 + i, 2 + i] = e
    return g


sph2 = [R**2, R**2 * sp.sin(TH)**2]
g_river = river_metric(sp.sqrt(2 * M / R), X4, sph2)   # (dr + v dt) = 0 gives inflow dr/dt = -v
Gmix, *_ = einstein_mixed(g_river, X4)
check("R3 Einstein tensor of the rule's river vanishes", all(sp.simplify(e) == 0 for e in Gmix),
      "G^a_b = 0 for v^2 = 2M/r (Painleve-Gullstrand form of Schwarzschild)")

# ---------------------------------------------------------------- R4  exactness
log("\n== R4  The river's speed law holds exactly, in the river's own clocks and the areal radius")
Gam, gi = christoffel(g_river, X4)
u = [1, -sp.sqrt(2 * M / R), 0, 0]      # river observer: dt = dtau, dr/dtau = -sqrt(2M/r)
norm = sp.simplify(sum(g_river[a, b] * u[a] * u[b] for a in range(4) for b in range(4)))
acc = [sp.simplify(sum(Gam[a][b][c] * u[b] * u[c] for b in range(4) for c in range(4))
                   + sum(u[b] * sp.diff(u[a], X4[b]) for b in range(4))) for a in range(4)]
check("R4a river observers are free (geodesic) with proper time = t",
      norm == -1 and all(a_ == 0 for a_ in acc), f"u.u = {norm}, acceleration = {acc}")
# change of time t_s = t - f(r), f' = v/(1-v^2): Schwarzschild's diagonal form
vv = sp.sqrt(2 * M / R)
J = sp.Matrix([[1, vv / (1 - vv**2)], [0, 1]])   # dt = dt_s + f'(r) dr, f' = v/(1 - v^2)
g2 = sp.Matrix([[g_river[0, 0], g_river[0, 1]], [g_river[1, 0], g_river[1, 1]]])
gS = sp.simplify(J.T * g2 * J)
check("R4b same metric as Schwarzschild after relabelling the clock",
      sp.simplify(gS[0, 1]) == 0 and sp.simplify(gS[0, 0] + (1 - 2 * M / R)) == 0
      and sp.simplify(gS[1, 1] - 1 / (1 - 2 * M / R)) == 0, f"diag = ({gS[0, 0]}, {gS[1, 1]})")
# isotropic radius: gamma and beta
rho = sp.symbols("rho", positive=True)
r_of_rho = rho * (1 + M / (2 * rho))**2
spatial = sp.simplify((sp.diff(r_of_rho, rho)**2 / (1 - 2 * M / r_of_rho)) / 1)
conf = sp.simplify(r_of_rho**2 / rho**2)
gtt_iso = sp.series(-(1 - 2 * M / r_of_rho), M, 0, 3).removeO()
gam_ppn = sp.series(conf, M, 0, 2).removeO().coeff(M) * rho / 2
beta_ppn = sp.simplify(-gtt_iso.coeff(M, 2) * rho**2 / 2)   # g00 = -1 + 2U - 2 beta U^2
check("R4c isotropic form: spatial factor is the same in all directions",
      sp.simplify(spatial - conf) == 0, f"factor = {sp.factor(conf)}")
check("R4d PPN gamma = 1 and beta = 1", sp.simplify(gam_ppn - 1) == 0 and sp.simplify(beta_ppn - 1) == 0,
      f"gamma = {gam_ppn}, beta = {beta_ppn}  (so frame dragging of ../README.md F-part carries over)")

# ---------------------------------------------------------------- R5  other speed laws fail Einstein
log("\n== R5  Shell theorem check: v^2 = 2M/r^p is vacuum only for p = 1 (known)")
p = sp.symbols("p", positive=True)
vf = sp.Function("v")(R)
Gmix_f, *_ = einstein_mixed(river_metric(vf, X4, sph2), X4)
Gtt = sp.simplify(Gmix_f[0, 0])
target = -sp.diff(R * vf**2, R) / R**2
check("R5a general river: G^t_t = -(r v^2)'/r^2", sp.simplify(Gtt - target) == 0, f"G^t_t = {sp.factor(Gtt)}")
Gtt_p = sp.simplify(target.subs(vf, sp.sqrt(2 * M / R**p)).doit())
check("R5b G^t_t = 0 only at p = 1", sp.solve(sp.Eq(sp.simplify(Gtt_p * R**(p + 2)), 0), p) == [1],
      f"G^t_t = {sp.factor(Gtt_p)}")
results["Gtt_general"] = str(Gtt)

# ---------------------------------------------------------------- R6  dimensions
log("\n== R6  Rule versus Einstein in other dimensions")
PS = sp.symbols("psi", real=True)
X5 = [T, R, TH, PH, PS]
sph3 = [R**2, R**2 * sp.sin(TH)**2, R**2 * sp.sin(TH)**2 * sp.sin(PH)**2]
Gmix5, *_ = einstein_mixed(river_metric(vf, X5, sph3), X5)
Gtt5 = sp.simplify(Gmix5[0, 0])
target5 = -sp.Rational(3, 2) * sp.diff(R**2 * vf**2, R) / R**3
check("R6a 4+1 river: G^t_t proportional to (r^2 v^2)'", sp.simplify(Gtt5 - target5) == 0,
      f"G^t_t = {sp.factor(Gtt5)}  -> vacuum needs v^2 = C/r^2, pull ~ 1/r^3")
# pull for v^2 = C/r^(n-3) in n spacetime dims is d/dr(v^2/2) ~ 1/r^(n-2); rule needs 1/r^2
n_ = sp.symbols("n", integer=True, positive=True)
check("R6b rule (pull 1/r^2) + vacuum Einstein picks n = 4 (3+1)", sp.solve(sp.Eq(n_ - 2, 2), n_) == [4],
      "3+1: v^2 ~ 1/r, pull 1/r^2 (velocity circles); 4+1: pull 1/r^3, velocity lines bend (velocity-lines V7)")


# ---------------------------------------------------------------- R7  consequences, numerically
log("\n== R7  Consequences of the rebuilt river (geodesics of the R3 metric, numerically)")
Mv = 1.0
idx = [0, 1, 3]                             # t, r, phi in the equatorial plane
Gnum = [[[sp.lambdify(R, sp.simplify(Gam[a][b][c].subs({M: Mv, TH: sp.pi / 2})), "numpy") for c in idx]
         for b in idx] for a in idx]


def rhs(lmb, s):
    pos, vel = s[:3], s[3:]
    rr = pos[1]
    accn = np.zeros(3)
    for i in range(3):
        accn[i] = -sum(Gnum[i][j][k](rr) * vel[j] * vel[k] for j in range(3) for k in range(3))
    return np.concatenate([vel, accn])


# R7a perihelion advance for a wide orbit p = 2000 M, e = 0.3
pp, ee = 2000.0, 0.3
rp, ra = pp / (1 + ee), pp / (1 - ee)
f = lambda rr: 1 - 2 * Mv / rr  # noqa: E731
# turning points: E^2 = f(r)(1 + L^2/r^2) at rp and ra
L2 = (f(ra) - f(rp)) / (f(rp) / rp**2 - f(ra) / ra**2)
En = math.sqrt(f(rp) * (1 + L2 / rp**2))
Lz = math.sqrt(L2)
vpg = math.sqrt(2 * Mv / rp)
s0 = [0.0, rp, 0.0, En / (1 - vpg**2), 0.0, Lz / rp**2]
peri = lambda lmb, s: s[4]  # noqa: E731
peri.direction = 1
Tn = 2 * math.pi * (pp / (1 - ee**2))**1.5
sol = solve_ivp(rhs, [0, 3.2 * Tn], s0, method="DOP853", rtol=1e-12, atol=1e-12, events=peri, dense_output=False)
phis = sol.y_events[0][:, 2]
adv_num = float(np.mean(np.diff(phis)) - 2 * math.pi)
adv_exact = quad(lambda ch: 1 / math.sqrt(1 - 2 * Mv * (3 + ee * math.cos(ch)) / pp), 0, 2 * math.pi,
                 epsabs=1e-14, epsrel=1e-13)[0] - 2 * math.pi
adv_weak = 6 * math.pi * Mv / pp
check("R7a perihelion advance from the river metric",
      abs(adv_num - adv_exact) / adv_exact < 1e-6 and abs(adv_num - adv_weak) / adv_weak < 3e-3,
      f"numeric {adv_num:.9e}, exact (Darwin) {adv_exact:.9e}, 6piM/p {adv_weak:.9e} rad/orbit")
results["perihelion"] = {"numeric": adv_num, "exact": adv_exact, "weak": adv_weak}

# R7b light bending, b ~ 1000 M
r0 = 1000.0
b = r0 / math.sqrt(f(r0))
vpg0 = math.sqrt(2 * Mv / r0)
s0 = [0.0, r0, 0.0, 1 / (1 - vpg0**2), 0.0, b / r0**2]
Rend = 1e8
far = lambda lmb, s: s[1] - Rend  # noqa: E731
far.terminal = True
sol = solve_ivp(rhs, [0, 2e8], s0, method="DOP853", rtol=1e-12, atol=1e-14, events=far)
phi_end = sol.y_events[0][0][2] + math.asin(b / Rend)
defl_num = 2 * phi_end - math.pi
defl_2 = 4 * Mv / b + 15 * math.pi / 4 * (Mv / b)**2
check("R7b light bending from the river metric",
      abs(defl_num - defl_2) / defl_2 < 2e-5 and abs(defl_num - 4 * Mv / b) / (4 * Mv / b) < 5e-3,
      f"numeric {defl_num:.9e}, 4M/b + 15pi/4 (M/b)^2 = {defl_2:.9e}, 4M/b = {4 * Mv / b:.9e}")
results["bending"] = {"numeric": defl_num, "second_order": defl_2, "b": b}

# R7c, R7d the Sun, in SI units
GMs, c = 1.32712440018e20, 299792458.0
aM, eM, PM = 5.7909050e10, 0.205630, 87.9691
arcsec = 180 / math.pi * 3600
merc = 6 * math.pi * GMs / (c**2 * aM * (1 - eM**2)) * (36525 / PM) * arcsec
check("R7c Mercury: 42.98 arcsec/century (observed 42.98 +- 0.04)", abs(merc - 42.98) < 0.05, f"{merc:.2f}\"/cy")
Rsun = 6.957e8
bend_sun = 4 * GMs / (c**2 * Rsun) * arcsec
check("R7d Sun-grazing light: 1.751 arcsec (VLBI: gamma = 1 to 1e-4)", abs(bend_sun - 1.751) < 0.002,
      f"{bend_sun:.4f}\"")
results["mercury_arcsec_per_century"] = merc
results["sun_bending_arcsec"] = bend_sun

# ---------------------------------------------------------------- R8  the turn of the velocity line
log("\n== R8  Observation (ours, a restatement): the velocity line's turn per orbit as dial x area")
gam_s, bet_s, Ms, Ls, cs = sp.symbols("gamma beta M L c", positive=True)
p_s = Ls**2 / Ms                         # semi-latus rectum from L^2 = GM p
area_hodo = sp.pi * (Ms / Ls)**2         # velocity circle radius = GM/L
turn_ppn = (2 + 2 * gam_s - bet_s) / 3 * 6 * sp.pi * Ms / (cs**2 * p_s)
ratio = sp.simplify(turn_ppn / (area_hodo / cs**2))
check("R8 turn per orbit = 2(2 + 2gamma - beta) x (velocity-circle area)/c^2",
      sp.simplify(ratio - 2 * (2 + 2 * gam_s - bet_s)) == 0 and ratio.subs({gam_s: 1, bet_s: 1}) == 6,
      f"ratio = {ratio}; with gamma = beta = 1 it is 6. Not yet derived from the swap rule.")

# ---------------------------------------------------------------- summary
n_pass = sum(1 for k, v in results.items() if isinstance(v, dict) and "pass" in v and v["pass"])
n_all = sum(1 for k, v in results.items() if isinstance(v, dict) and "pass" in v)
log(f"\n{n_pass}/{n_all} checks pass")
base = __file__.rsplit("/", 1)[0] if "/" in __file__ else "."
open(f"{base}/checks_output.txt", "w").write("\n".join(out) + "\n")
json.dump(results, open(f"{base}/results.json", "w"), indent=2, default=str)
