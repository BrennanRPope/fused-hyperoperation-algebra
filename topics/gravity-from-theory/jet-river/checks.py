"""DRAFT checks: the velocity rule for the whole river, as a jet field (Nov 26, 2026).

The river's state at each cell is a jet: velocity v, its gradient S (first derivative), and the gradient of S
(second derivative). The jet-fields page's box of one time step, det(I + h S) = 1 + hP + h^2 Q + h^3 R, gives
the layers. Question: does a law on these layers force exact Schwarzschild without borrowing Einstein's
equation or the "proper time and areal radius" bridge?
Run: python3 checks.py   (writes checks_output.txt, results.json)
"""
import json
import os
import sys

import sympy as sp

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from grtools import christoffel, einstein_mixed, riemann  # noqa: E402
from lintools import lin_einstein  # noqa: E402

out, results = [], {}


def log(s=""):
    print(s)
    out.append(s)


def check(name, ok, detail=""):
    results[name] = {"pass": bool(ok), "detail": str(detail)}
    log(f"[{'PASS' if ok else 'FAIL'}] {name}: {detail}")


T, R, TH, PH = sp.symbols("t r theta phi", real=True)
X4 = [T, R, TH, PH]
h, M, G_, rho, H = sp.symbols("h M G rho H", positive=True)
v = sp.Function("v")(R)


def river(vv, extra_tphi=0):
    """-dt^2 + (dr + v dt)^2 + r^2 dOmega^2 ; flat slices, unit lapse. (dr + v dt = 0: inflow dr/dt = -v.)"""
    g = sp.zeros(4, 4)
    g[0, 0] = -1 + vv**2
    g[0, 1] = g[1, 0] = vv
    g[1, 1] = 1
    g[2, 2] = R**2
    g[3, 3] = R**2 * sp.sin(TH)**2
    return g


# ---------------------------------------------------------------- J1 the jet's box for the river
log("== J1  The river cell's one-step box: det(I + h S) = 1 + hP + h^2 Q + h^3 R")
a_, b_ = sp.diff(v, R), v / R            # inflow -v r_hat: strain eigen-rates -v' (radial), -v/r (sideways, twice)
S = -sp.diag(a_, b_, b_)
box = sp.expand(sp.det(sp.eye(3) + h * S))
P, Q, Rr = [sp.simplify(box.coeff(h, k)) for k in (1, 2, 3)]
check("J1a layers of the river's strain box", sp.simplify(Q - sp.diff(R * v**2, R) / R**2) == 0,
      f"P = {P},  Q = {sp.factor(Q)} = (r v^2)'/r^2,  R = {Rr}")

# ---------------------------------------------------------------- J2 Einstein's tensor is the Q layer
log("\n== J2  Every Einstein component of the river is the Q layer or its derivative")
Gm, *_ = einstein_mixed(river(v), X4)
Gtt, Grr, Gthth = [sp.simplify(Gm[i, i]) for i in (0, 1, 2)]
check("J2a G^t_t = -Q", sp.simplify(Gtt + Q) == 0, f"G^t_t = {sp.factor(Gtt)}")
check("J2b G^r_r = -Q", sp.simplify(Grr + Q) == 0, f"G^r_r = {sp.factor(Grr)}")
check("J2c G^theta_theta = G^phi_phi = -(r^2 Q)'/(2r)",
      sp.simplify(Gthth + sp.diff(R**2 * Q, R) / (2 * R)) == 0 and sp.simplify(Gm[3, 3] - Gthth) == 0,
      f"G^th_th = {sp.factor(Gthth)}")
offd = [sp.simplify(Gm[i, j]) for i in range(4) for j in range(4) if i != j]
check("J2d all off-diagonal components vanish (static river)", all(e == 0 for e in offd),
      "so vacuum <=> Q = 0 and (r^2 Q)' = 0 <=> Q = 0")
results["einstein_components"] = {"Gtt": str(Gtt), "Grr": str(Grr), "Gthth": str(Gthth)}

# ---------------------------------------------------------------- J3 Q is a divergence: exact Gauss law
log("\n== J3  Q is a null Lagrangian (jet-fields), so it integrates to an exact Gauss law")
x, y, z = sp.symbols("x y z", real=True)
u = [sp.Function(f"u{i}")(x, y, z) for i in range(3)]
Xs = [x, y, z]
J = sp.Matrix(3, 3, lambda i, j: sp.diff(u[i], Xs[j]))
e2J = J[0, 0] * J[1, 1] - J[0, 1] * J[1, 0] + J[0, 0] * J[2, 2] - J[0, 2] * J[2, 0] + J[1, 1] * J[2, 2] - J[1, 2] * J[2, 1]
flux = [sp.Rational(1, 2) * (u[i] * sum(J[j, j] for j in range(3)) - sum(u[j] * J[i, j] for j in range(3))) for i in range(3)]
div = sum(sp.diff(flux[i], Xs[i]) for i in range(3))
check("J3a Q(grad u) = div F for any 3D field, F = (u (div u) - (u.grad)u)/2", sp.simplify(sp.expand(div - e2J)) == 0,
      "the second layer of the box is a pure divergence")
rr = sp.symbols("rr", positive=True)
Fr = sp.simplify(sp.Rational(1, 2) * (-v * P - v * sp.diff(v, R)))      # radial F for inflow u = -v(r) r_hat
check("J3b over a ball: integral of Q = 4 pi r^2 F_r = 4 pi r v^2",
      sp.simplify(4 * sp.pi * R**2 * Fr - 4 * sp.pi * R * v**2) == 0
      and sp.simplify(sp.integrate((R**2 * Q).subs(v, R**sp.Rational(3, 4)).doit().subs(R, rr), (rr, 0, R))
                      - R * R**sp.Rational(3, 2)) == 0,
      "so Q = 8 pi G rho integrates to r v^2 = 2 G m(r): the river speed law with the enclosed mass, exactly")

# ---------------------------------------------------------------- J4 every cell at zero energy
log("\n== J4  The velocity rule for every cell at once: each cell is a zero-energy free faller")
m = sp.Function("m")(R)
cell_energy = v**2 / 2 - G_ * m / R
vsol = sp.sqrt(2 * G_ * m / R)
Qcell = sp.simplify(Q.subs(v, vsol).doit())
check("J4a all cells at E = 0 (v^2 = 2Gm(r)/r) <=> Q = 8 pi G rho, with m' = 4 pi r^2 rho",
      sp.simplify(Qcell - 2 * G_ * sp.diff(m, R) / R**2) == 0, f"Q = {Qcell} = 8 pi G rho")
vac = sp.dsolve(sp.Eq(sp.diff(R * sp.Function("w")(R), R), 0))
check("J4b empty space: Q = 0 gives r v^2 = constant = 2GM (Schwarzschild's river)", "C1" in str(vac.rhs),
      f"r v^2 = {vac.rhs}")
check("J4c uniform ball: v = H r, Q = 3H^2 = 8 pi G rho (Friedmann's equation)",
      sp.simplify(Q.subs(v, H * R).doit() - 3 * H**2) == 0, "inside a uniform dust ball the river is a contracting Friedmann universe")

# ---------------------------------------------------------------- J5 Birkhoff: m cannot change in empty space
log("\n== J5  In empty space the enclosed mass cannot change (Birkhoff)")
mt = sp.Function("m")(T)
vt = sp.sqrt(2 * mt / R)
gt = river(vt)
Gt, *_ = einstein_mixed(gt, X4)
nonzero = {f"G[{i},{j}]": sp.factor(sp.simplify(Gt[i, j])) for i in range(4) for j in range(4) if sp.simplify(Gt[i, j]) != 0}
only_mdot = all(sp.simplify(e.subs(sp.Derivative(mt, T), 0)) == 0 for e in nonzero.values())
check("J5 a river v^2 = 2m(t)/r breaks the field equations unless dm/dt = 0", len(nonzero) > 0 and only_mdot,
      nonzero)

# ---------------------------------------------------------------- J6 second jet: tides are curvature
log("\n== J6  The second derivative of the jet is the tide (curvature)")
g = river(v)
Gam, gi = christoffel(g, X4)
Riem = riemann(Gam, X4)
uvec = [1, -v, 0, 0]                    # river cell 4-velocity (proper time = t)
er = [0, 1, 0, 0]                       # radial unit vector (orthogonal to u, unit length)
eth = [0, 0, 1 / R, 0]


def tidal(e):
    # E(e,e) = R_{a b c d} e^a u^b e^c u^d  with R_{abcd} = g_{a s} R^s_{bcd}
    tot = 0
    for A in range(4):
        for B in range(4):
            for C in range(4):
                for D in range(4):
                    if e[A] == 0 or uvec[B] == 0 or e[C] == 0 or uvec[D] == 0:
                        continue
                    low = sum(g[A, s] * Riem[s][B][C][D] for s in range(4))
                    tot += low * e[A] * uvec[B] * e[C] * uvec[D]
    return sp.simplify(tot)


Err, Ethth = tidal(er), tidal(eth)
vp, vpp = sp.diff(v, R), sp.diff(v, R, 2)
jet_rr = -(v * vpp + vp**2)                 # -(v.grad S + S^2) radial: D S/Dt = -v d/dr for inflow cells
jet_thth = -(v * sp.diff(v / R, R) + (v / R)**2)
check("J6a Riemann tide R(e,u,e,u) = -(DS/Dt + S^2) following the cell, radial and sideways, any v(r)",
      sp.simplify(Err - jet_rr) == 0 and sp.simplify(Ethth - jet_thth) == 0,
      f"R_(r u r u) = {sp.factor(Err)};  R_(th u th u) = {sp.factor(Ethth)}")
trace = sp.simplify(Err + 2 * Ethth)
lap = sp.simplify(sp.diff(R**2 * sp.diff(v**2 / 2, R), R) / R**2)
check("J6b trace of the tide R_uu = Laplacian of Phi with Phi = -v^2/2 (Newton's Poisson equation, exact)",
      sp.simplify(trace + lap) == 0, f"R_uu = {sp.factor(trace)}")
vS = sp.sqrt(2 * M / R)
check("J6c Schwarzschild river: tides -2M/r^3 along r, +M/r^3 sideways; trace 0",
      sp.simplify(Err.subs(v, vS).doit() + 2 * M / R**3) == 0 and sp.simplify(Ethth.subs(v, vS).doit() - M / R**3) == 0,
      "geodesic deviation a = -R(.,u,.,u) xi: stretched along r by 2M/r^3, squeezed sideways by M/r^3")

# ---------------------------------------------------------------- J7 the vector law and frame dragging
log("\n== J7  The Q layer's gradient gives the vector law; a swirl gives frame dragging; full spin bends the lattice")
s = sp.symbols("s11 s12 s13 s22 s23 s33")
Sg = sp.Matrix([[s[0], s[1], s[2]], [s[1], s[3], s[4]], [s[2], s[4], s[5]]])
e2S = sp.Rational(1, 2) * (Sg.trace()**2 - (Sg * Sg).trace())
grad = sp.Matrix(3, 3, lambda i, j: sp.diff(e2S, s[[0, 1, 2, 3, 4, 5][[(0, 0), (0, 1), (0, 2), (1, 1), (1, 2), (2, 2)].index(tuple(sorted((i, j))))]]) * (1 if i == j else sp.Rational(1, 2)))
check("J7a dQ/dS = (tr S) I - S", sp.simplify(grad - (Sg.trace() * sp.eye(3) - Sg)) == sp.zeros(3, 3),
      "momentum constraint (vacuum): div((tr S) I - S) = 0, the divergence of Q's gradient")
Jsym, a_ = sp.symbols("J a", positive=True)
vS = sp.sqrt(2 * M / R)
w = sp.Function("w")(R)


def swirl(vv, ww, twist=0):
    """first-order pieces added to the river: cells also circle at angular speed w (dphi -> dphi - w dt);
    optional twist g_r phi (the cell's own rotation in Doran's Kerr river, first order in spin)."""
    g1 = sp.zeros(4, 4)
    g1[0, 3] = g1[3, 0] = -ww * R**2 * sp.sin(TH)**2
    g1[1, 3] = g1[3, 1] = twist
    return g1


nz = lambda Gm: {f"G[{i},{j}]": Gm[i, j] for i in range(4) for j in range(4) if Gm[i, j] != 0}
GW = lin_einstein(river(v), swirl(v, w), X4)
check("J7b vector law for a swirl w(r) on any river: G^t_phi = -r (r w'' + 4 w') sin^2/2, so (r^4 w')' = 0, w = 2J/r^3",
      sp.simplify(GW[0, 3] + R * (R * sp.diff(w, R, 2) + 4 * sp.diff(w, R)) * sp.sin(TH)**2 / 2) == 0,
      "the Lense-Thirring frame-dragging rate, from the Q layer's vector law (momentum constraint)")
fA = nz(lin_einstein(river(vS), swirl(vS, 2 * Jsym / R**3), X4))
check("J7c Schwarzschild river + swirl w = 2J/r^3 satisfies every field equation at first order in spin",
      len(fA) == 0, "frame dragging is a velocity field: inflow v = sqrt(2M/r) plus swirl 2J/r^3 (slow-spin Kerr)")
fB = nz(lin_einstein(river(vS), swirl(vS, 2 * Jsym / R**3, -(Jsym / M) * vS * sp.sin(TH)**2), X4))
check("J7d Doran's twist term g_r phi = -a v sin^2 also passes: at first order it is only a relabeling phi -> phi + h(r)",
      len(fB) == 0, "r^2 sin^2 h' = -a v; the cell's own rotation shows up only at second order in spin")

# second order in spin: the Kerr river's lattice is no longer flat
e_ = sp.symbols("epsilon", positive=True)
A_ = e_ * a_
s2 = sp.sin(TH)**2
rho2 = R**2 + A_**2 * sp.cos(TH)**2
gD = sp.zeros(3, 3)                                     # Doran's Kerr river at t = const
gD[0, 0] = rho2 / (R**2 + A_**2)
gD[0, 2] = gD[2, 0] = -sp.sqrt(2 * M * R) * A_ * s2 / sp.sqrt(R**2 + A_**2)
gD[1, 1] = rho2
gD[2, 2] = (R**2 + A_**2) * s2 + 2 * M * R * A_**2 * s2**2 / rho2
rp = sp.symbols("r", positive=True)
gD = gD.subs(R, rp)
X3 = [rp, TH, PH]
tr2 = lambda x: sp.series(x, e_, 0, 3).removeO()
gD = gD.applyfunc(tr2)
gDi = gD.inv().applyfunc(lambda x: tr2(sp.simplify(x)))
C3 = [[[tr2(sp.expand(sum(gDi[i, d] * (sp.diff(gD[d, j], X3[k]) + sp.diff(gD[d, k], X3[j]) - sp.diff(gD[j, k], X3[d]))
                         for d in range(3)) / 2)) for k in range(3)] for j in range(3)] for i in range(3)]
Ric3 = sp.Matrix(3, 3, lambda b, d: tr2(sp.expand(sum(sp.diff(C3[i][d][b], X3[i]) - sp.diff(C3[i][i][b], X3[d]) for i in range(3))
                                                  + sum(C3[i][i][k] * C3[k][d][b] - C3[i][d][k] * C3[k][i][b]
                                                        for i in range(3) for k in range(3)))))
Rs3 = sp.simplify(tr2(sp.expand(sum(gDi[i, j] * Ric3[i, j] for i in range(3) for j in range(3)))))
first_flat = all(sp.simplify(x.coeff(e_, 1)) == 0 for x in Ric3)
check("J7e Doran's Kerr river: lattice flat at first order in spin, curved at second: 3R = 2 M a^2 (3 sin^2 - 2)/r^5",
      first_flat and sp.simplify(Rs3 - 2 * M * a_**2 * e_**2 * (3 * s2 - 2) / rp**5) == 0,
      f"3R = {Rs3}; known theorem (Garat & Price 2000): no conformally flat slicing of Kerr exists at all")
results["J7e_scalar_curvature"] = str(Rs3)

# ---------------------------------------------------------------- summary
n_pass = sum(1 for v_ in results.values() if isinstance(v_, dict) and v_.get("pass"))
n_all = sum(1 for v_ in results.values() if isinstance(v_, dict) and "pass" in v_)
log(f"\n{n_pass}/{n_all} checks pass")
open(os.path.join(HERE, "checks_output.txt"), "w").write("\n".join(out) + "\n")
json.dump(results, open(os.path.join(HERE, "results.json"), "w"), indent=2, default=str)
