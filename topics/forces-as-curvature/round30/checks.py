"""DRAFT round 30: a lump's own tidal pull, built from the law (round 23's gap; round 25's tide = 3/2 (qbar - q_loc)).
Question: rounds 23 and 25 found the law accounts for the local, direction-blind part of curvature (q set by the matter
here) but not the tide, and concluded 'the law needs Einstein's equation itself'. Here: split Einstein's constraint on
a slice at rest into (i) a local trace rule, which is the law's own (q_loc = q + 8 pi rho/3), and (ii) the contracted
Bianchi identity, which every geometry obeys, and check whether those two, with the dial allowed to differ by direction,
build the tide; and count what is left over. Exact parts in sympy; numbers with numpy/scipy. Units G = c = 1, Lambda = 3q."""
import json, io, time, math
import numpy as np
import sympy as sp
from scipy.integrate import solve_ivp, quad

T0 = time.time(); _log = io.StringIO(); passed = failed = 0; results = {}
def say(s): print(s); _log.write(s + "\n")
def check(name, ok, detail):
    global passed, failed
    passed += bool(ok); failed += (not ok)
    say(f"{'PASS' if ok else 'FAIL'} {name}  | {detail}")

def curvature(g, X):
    """Christoffel, Riemann R^a_bcd, Ricci, scalar for a metric matrix g in coordinates X."""
    n = len(X); gi = g.inv()
    G = [[[sp.simplify(sum(gi[a, d] * (sp.diff(g[d, b], X[c]) + sp.diff(g[d, c], X[b]) - sp.diff(g[b, c], X[d]))
                           for d in range(n)) / 2) for c in range(n)] for b in range(n)] for a in range(n)]
    def Riem(a, b, c, d):
        return sp.diff(G[a][b][d], X[c]) - sp.diff(G[a][b][c], X[d]) + sum(G[a][c][e] * G[e][b][d] - G[a][d][e] * G[e][b][c] for e in range(n))
    Ric = sp.Matrix(n, n, lambda b, d: sp.simplify(sum(Riem(a, b, a, d) for a in range(n))))
    Rs = sp.simplify(sum(gi[i, j] * Ric[i, j] for i in range(n) for j in range(n)))
    return G, Riem, Ric, Rs, gi

# ---------------------------------------------------------------- L1
say("== L1. A round lump's slice at rest is the law in every plane, with a dial that depends on direction")
r, th, ph = sp.symbols('r theta phi', positive=True)
Qb = sp.Function('qbar')(r)
X = [r, th, ph]
g = sp.diag(1 / (1 - Qb * r**2), r**2, r**2 * sp.sin(th)**2)
Gm, Riem, Ric, Rs, gi = curvature(g, X)
def sec(i, j):  # sectional curvature of the coordinate plane (i, j), lowered R_ijij / (g_ii g_jj)
    low = sum(g[i, a] * Riem(a, j, i, j) for a in range(3))
    return sp.simplify(low / (g[i, i] * g[j, j]))
K_across = sec(1, 2); K_along = sec(0, 1); K_along2 = sec(0, 2)
q_loc_expr = sp.simplify(Rs / 6)
tide = sp.simplify(K_across - K_along)
ok1 = (sp.simplify(K_across - Qb) == 0 and sp.simplify(K_along - (Qb + r * Qb.diff(r) / 2)) == 0
       and sp.simplify(K_along2 - K_along) == 0
       and sp.simplify(q_loc_expr - (2 * K_along + K_across) / 3) == 0
       and sp.simplify(tide - sp.Rational(3, 2) * (Qb - q_loc_expr)) == 0)
check("L1 (exact) on the slice at rest of any round lump, dr^2/(1 - qbar(r) r^2) + r^2 dOmega^2, every plane bends like a piece of the law's 3-sphere, but with its own dial: planes across the radius at qbar, planes along it at qbar + r qbar'/2. The mean of the three bendings is the local dial q_loc = R/6 (the law's direction-blind part, rounds 23 and 25); the split across - along = -r qbar'/2 = 3/2 (qbar - q_loc) is the tide. In three dimensions the curvature is exactly this set of plane bendings (Riemann is fixed by Ricci), so a direction-dependent dial is all the slice's curvature there is",
      ok1, f"across {K_across}; along {K_along}; scalar/6 {q_loc_expr}; tide {tide}")
results["L1"] = dict(across=str(K_across), along=str(K_along), q_loc=str(q_loc_expr), tide=str(tide))

# ---------------------------------------------------------------- L2
say("== L2. Fitting (the Bianchi identity) plus the law's local rule builds the tide")
# mixed Einstein tensor G^i_j = Ric^i_j - R/2 delta; divergence nabla_i G^i_r
Rmix = gi * Ric
Gmix = sp.simplify(Rmix - Rs / 2 * sp.eye(3))
div_r = sp.simplify(sum(sp.diff(Gmix[i, 0], X[i]) for i in range(3))
                    + sum(Gm[i][i][k] * Gmix[k, 0] for i in range(3) for k in range(3))
                    - sum(Gm[k][i][0] * Gmix[i, k] for i in range(3) for k in range(3)))
Grr = sp.simplify(Gmix[0, 0]); Gthth = sp.simplify(Gmix[1, 1])
bianchi_form = sp.simplify(Grr.diff(r) + 2 / r * (Grr - Gthth))
# with G^r_r = -K_across and G^th_th = -K_along, Bianchi reads qbar' = -(2/r)(K_across - K_along); with the law's
# local rule (mean of bendings = q_loc) this is qbar' = (3/r)(q_loc - qbar): qbar is the volume average of q_loc.
ql = sp.Function('q_loc')(r)
ode = sp.Eq(Qb.diff(r), 3 / r * (ql - Qb))
sol_check = sp.simplify(sp.diff((Qb * r**3), r) - 3 * r**2 * ql)          # (r^3 qbar)' = 3 r^2 q_loc  <=>  ode
ok2a = div_r == 0 and sp.simplify(Grr + Qb) == 0 and sp.simplify(Gthth + K_along) == 0 and sp.simplify(bianchi_form) == 0
# numbers: a smooth lump on the law's background; build qbar from the ODE alone, compare with Einstein's m(r)
q0, rho0 = 0.01, 0.05
rho = lambda x: rho0 * math.exp(-x * x)
qloc = lambda x: q0 + 8 * math.pi * rho(x) / 3
sol = solve_ivp(lambda x, y: [3 / x * (qloc(x) - y[0])], (1e-6, 6.0), [qloc(1e-6)], rtol=1e-11, atol=1e-14, dense_output=True)
rr = np.linspace(0.05, 6.0, 120)
m_of = lambda x: quad(lambda s: 4 * math.pi * s * s * rho(s), 0, x, epsabs=1e-14)[0]
qbar_E = np.array([q0 + 2 * m_of(x) / x**3 for x in rr])
qbar_ode = sol.sol(rr)[0]
err_q = float(np.max(np.abs(qbar_ode - qbar_E)))
M = m_of(30.0)
tide_ode = 1.5 * (qbar_ode - np.array([qloc(x) for x in rr]))
out = rr > 4.5
err_out = float(np.max(np.abs(tide_ode[out] - 3 * M / rr[out]**3) / (3 * M / rr[out]**3)))
same_ode = sp.simplify(sol_check - r**3 * (Qb.diff(r) - ode.rhs)) == 0       # (r^3 qbar)' - 3 r^2 q_loc = r^3 (qbar' - rhs)
ok2 = ok2a and same_ode and err_q < 1e-8 and err_out < 1e-6
check("L2 (exact + numbers) the contracted Bianchi identity, which every geometry obeys whatever the matter, holds identically here and says qbar' = -(2/r)(across - along): wherever the across-dial changes with radius, the along-dial must split from it. Add the law's local rule (mean bending = q_loc = q + 8 pi rho/3) and it becomes (r^3 qbar)' = 3 r^2 q_loc, so qbar is q_loc averaged over the ball inside, and the tide follows. For a smooth lump, building qbar from the law's rule and fitting alone matches Einstein's mass function, and the tide outside is 3M/r^3",
      ok2, f"nabla_i G^i_r = {div_r}; G^r_r = -qbar, G^th_th = -along; qbar from fitting vs q + 2m(r)/r^3: max error {err_q:.1e}; tide outside vs 3M/r^3: {err_out:.1e} (M = {M:.4f})")
results["L2"] = dict(div=str(div_r), err_qbar=err_q, err_tide_out=err_out, M=M, q0=q0, rho0=rho0,
                     r=rr.tolist(), qbar=qbar_ode.tolist(), qloc=[qloc(x) for x in rr], tide=tide_ode.tolist())

# ---------------------------------------------------------------- L3
say("== L3. Without round symmetry: what the law's rule and fitting leave free is a gravitational wave")
x, y, z, k, A_, B_, eps = sp.symbols('x y z k A B epsilon', real=True)
h = sp.Matrix([[A_ * sp.cos(k * z), B_ * sp.cos(k * z), 0], [B_ * sp.cos(k * z), -A_ * sp.cos(k * z), 0], [0, 0, 0]])
Xc = [x, y, z]
# linearised curvature of delta + eps h: R_ij = (d_k d_i h_kj + d_k d_j h_ki - d^2 h_ij - d_i d_j h)/2, R = d_i d_j h_ij - d^2 h
trh = h.trace()
Rlin = sp.Matrix(3, 3, lambda i, j: sp.simplify((sum(sp.diff(h[kk, j], Xc[kk], Xc[i]) + sp.diff(h[kk, i], Xc[kk], Xc[j]) for kk in range(3))
                                                  - sum(sp.diff(h[i, j], c, 2) for c in Xc) - sp.diff(trh, Xc[i], Xc[j])) / 2))
Rscal_lin = sp.simplify(sum(sp.diff(h[i, j], Xc[i], Xc[j]) for i in range(3) for j in range(3)) - sum(sp.diff(trh, c, 2) for c in Xc))
count = 6 - 3 - 1
ok3 = Rscal_lin == 0 and Rlin != sp.zeros(3, 3) and count == 2 and sp.simplify(Rlin.trace()) == 0
check("L3 (exact, linearised; textbook York 1973) away from round symmetry a slice's shape has 6 functions; 3 are only the choice of coordinates and the law's local rule fixes 1, leaving 2 per point. A transverse ripple in those 2 (both polarisations, A and B) keeps every local dial at its law value (scalar curvature unchanged) yet bends planes differently by direction: it is a tide carried by a gravitational wave (round 26). So matter's own tide is built by the law's rule plus fitting, but the 2 wave functions are free data the law does not fix",
      ok3, f"delta R = {Rscal_lin}; delta Ric = {[str(Rlin[0,0]), str(Rlin[0,1]), str(Rlin[2,2])]} (traceless); free functions per point = 6 - 3 - 1 = {count}")
results["L3"] = dict(dR=str(Rscal_lin), dRic=str(Rlin), free=count)

# ---------------------------------------------------------------- L4
say("== L4. Two lumps at rest: no round symmetry, tide built with no free wave data (Brill-Lindquist)")
m1, m2, d = 1.0, 0.6, 3.0
psi_s = 1 + m1 / (2 * sp.sqrt(x**2 + y**2 + (z - d / 2)**2)) + m2 / (2 * sp.sqrt(x**2 + y**2 + (z + d / 2)**2))
# Ricci of g = psi^4 delta in 3D (exact formula for a conformally flat metric)
dpsi = [sp.diff(psi_s, c) for c in Xc]
lap = sum(sp.diff(psi_s, c, 2) for c in Xc)
gradsq = sum(v**2 for v in dpsi)
RicBL = sp.Matrix(3, 3, lambda i, j: -2 / psi_s * sp.diff(psi_s, Xc[i], Xc[j]) + 6 / psi_s**2 * dpsi[i] * dpsi[j]
                  - (2 / psi_s**2 * gradsq + 2 / psi_s * lap) * (1 if i == j else 0))
RsBL = sum(RicBL[i, i] for i in range(3)) / psi_s**4
fR = sp.lambdify((x, y, z), RsBL, "numpy"); fRic = sp.lambdify((x, y, z), RicBL, "numpy"); fpsi = sp.lambdify((x, y, z), psi_s, "numpy")
rng = np.random.default_rng(7)
pts = rng.uniform(-6, 6, (400, 3)); pts = pts[(np.linalg.norm(pts - [0, 0, d / 2], axis=1) > 0.8) & (np.linalg.norm(pts + [0, 0, d / 2], axis=1) > 0.8)]
vac = float(np.max(np.abs([fR(*p) for p in pts])))
def eigs(p):
    Rm = np.array(fRic(*p), dtype=float) / fpsi(*p)**4        # mixed Ricci^i_j (orthonormal up to the conformal factor)
    return np.sort(np.linalg.eigvalsh(Rm))
# far field: eigenvalues -> (-2, 1, 1) (m1 + m2)/R^3 with areal radius R = rho psi^2
Rfar = 200.0; p = np.array([Rfar * 0.6, 0.0, Rfar * 0.8]); ev = eigs(p)
Rareal = Rfar * fpsi(*p)**2; Mt = m1 + m2
far_err = float(np.max(np.abs(ev - np.array([-2, 1, 1]) * Mt / Rareal**3)) / (Mt / Rareal**3))
# grid of principal stretch directions in the plane y = 0 for the picture
gx = np.linspace(-5, 5, 29); gz = np.linspace(-5, 5, 29); field = []
for xx in gx:
    for zz in gz:
        if min(math.hypot(xx, zz - d / 2), math.hypot(xx, zz + d / 2)) < 0.7: continue
        Rm = np.array(fRic(xx, 1e-9, zz), dtype=float) / fpsi(xx, 1e-9, zz)**4
        w, V = np.linalg.eigh(Rm); i = int(np.argmin(w))          # most negative Ricci = strongest stretch direction
        field.append([xx, zz, float(V[0, i]), float(V[2, i]), float(w[i]), float(w.max())])
ok4 = vac < 1e-9 and far_err < 0.02
check("L4 (exact + numbers; textbook Brill-Lindquist 1963) two lumps of mass 1 and 0.6 at rest, with no round symmetry: the slice psi^4 (flat) with psi = 1 + m1/2r1 + m2/2r2 keeps the local dial at the vacuum value everywhere between them (scalar curvature 0), while the planes bend differently by direction around both: the pair's tide, built with the wave data set to zero. Far away it becomes one lump's tide of the total mass, Ricci eigenvalues (-2, 1, 1) M/R^3",
      ok4, f"max |R| at {len(pts)} points outside the lumps: {vac:.1e}; far field at R = {Rareal:.0f}: eigenvalues {ev} vs (-2,1,1) x {Mt/Rareal**3:.3e}, relative error {far_err:.1e}")
results["L4"] = dict(m1=m1, m2=m2, d=d, vacuum_maxR=vac, far_err=far_err, field=field)

say(f"\n{passed} passed, {failed} failed  ({time.time() - T0:.0f} s)")
results["summary"] = dict(passed=passed, failed=failed)
json.dump(results, open("results.json", "w"), indent=1, default=str)
open("checks_output.txt", "w").write(_log.getvalue())
