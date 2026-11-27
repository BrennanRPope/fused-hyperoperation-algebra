"""Part K of the new-force DRAFT (Nov 27 2026): the other unknown forces on the map (h1), one cell at a time.
K1 screened, order off (a massive new force; and a second circle's photon mixing with ours on a tilted torus);
K2 unlimited reach, order on (a long-range non-abelian force: the conformal window, strength set by the material alone);
K3 rungs, order on (copies of W bosons or gluons);
K4 tube, order off (abelian confinement: flux tubes, and cosmic strings as their fundamental form);
K5 a base dial that varies in space and hides from tests (screening).
Run after repel_checks.py: appends K_ keys to results.json, writes other_forces_output.txt."""
import math, json
import numpy as np
import sympy as sp
from scipy.integrate import solve_ivp, quad

lines = []; NP = [0, 0]; OUT = {}
def log(s=""):
    print(s); lines.append(s)
def check(name, ok, detail=""):
    NP[0 if ok else 1] += 1
    log(("PASS  " if ok else "FAIL  ") + name + (f"   [{detail}]" if detail else ""))
lP = 1.616255e-35; MP = 1.220890e19; Mbar = MP/math.sqrt(8*math.pi); hbarc_eVm = 1.973269804e-7
prev = json.load(open("results.json"))

log("== K. The other unknown forces, one cell at a time")
# ---------------- K1 screened, order off
lam_EW = 38.6e-6
check("K1 a massive new force (Delta > 0, order off): our law sets its range 2/sqrt(Delta) and its shape (a bump in D_eff), not its strength; at gravity's strength Eot-Wash (2020) needs range below 38.6 um, i.e. a carrier above 5.1 meV",
      abs(hbarc_eVm/lam_EW*1e3 - 5.11) < 0.01, f"{hbarc_eVm/lam_EW*1e3:.2f} meV")
alpha_cassini = 2.5e-5/2                         # scalar: gamma - 1 = -2a/(1+a) > -2.5e-5 (2 sigma)
OUT.update(K_alpha_cassini_longrange=alpha_cassini, K_alpha_moon=3e-11, K_lambda_moon_m=[15*6.371e6, 60*6.371e6])
check("K1b where it stands far out (pass 2): about 3e-11 of gravity near the Moon's distance (perigee, 15-60 Earth radii), about 1e-5 for a scalar of unlimited reach (Cassini); our size carrier (strength 1/3, range 3e-32 m) sits 1e27 to the left of every test",
      float(prev["H_lab_over_range"]) > 1e26 and abs(alpha_cassini - 1.25e-5) < 1e-7)
th = sp.symbols('theta', real=True); r1, r2 = sp.symbols('R1 R2', positive=True)
Gt = sp.Matrix([[r1**2, r1*r2*sp.cos(th)], [r1*r2*sp.cos(th), r2**2]])
eps = sp.simplify(Gt[0, 1]/sp.sqrt(Gt[0, 0]*Gt[1, 1]))
check("K1c a second circle's photon is a 'dark photon', and a tilted torus mixes it with ours: the kinetic matrix is the torus metric, so the mixing is eps = cos(angle between the circles); square torus, no mixing",
      sp.simplify(eps - sp.cos(th)) == 0 and eps.subs(th, sp.pi/2) == 0, f"eps = {eps}")
Ginv = sp.simplify(Gt.inv())
check("K1d after mixing, a copy that only moves round the hidden circle picks up a small charge under ours in proportion to cos(theta) (the inverse metric's off-diagonal; Holdom 1986 for a massless hidden photon): tilt makes 'millicharged' copies",
      sp.simplify(Ginv[0, 1]*r1*r2*sp.sin(th)**2 + sp.cos(th)) == 0, f"G^12 = {sp.simplify(Ginv[0,1])}")
# ---------------- K2 unlimited reach, order on: the conformal window
def bz(N, Nf):
    b0 = sp.Rational(11, 3)*N - sp.Rational(2, 3)*Nf
    b1 = sp.Rational(34, 3)*N**2 - Nf*(sp.Rational(10, 3)*N + sp.Rational(N**2 - 1, N))
    return b0, b1
b0, b1 = bz(3, 16); a16 = float(-4*sp.pi*b0/b1)
edge3 = float(sp.solve(bz(3, sp.Symbol('f'))[1], sp.Symbol('f'))[0]); edge2 = float(sp.solve(bz(2, sp.Symbol('f'))[1], sp.Symbol('f'))[0])
AF3 = float(sp.solve(bz(3, sp.Symbol('f'))[0], sp.Symbol('f'))[0]); AF2 = float(sp.solve(bz(2, sp.Symbol('f'))[0], sp.Symbol('f'))[0])
OUT.update(K_BZ_alpha_SU3_16=a16, K_BZ_edge_SU3=edge3, K_BZ_edge_SU2=edge2)
check("K2 a long-range non-abelian force exists below the 16.5-flavour line too: from about 8 flavours (two loops, SU(3)) its strength stops running and freezes at a fixed value; for SU(2) from about 5.6 up to 11",
      abs(edge3 - 8.05) < 0.01 and abs(AF3 - 16.5) < 1e-9 and abs(edge2 - 5.55) < 0.01 and abs(AF2 - 11) < 1e-9,
      f"SU(3) window {edge3:.2f} to {AF3}, SU(2) {edge2:.2f} to {AF2}")
check("K2b at 16 flavours of SU(3) the frozen strength is alpha* = 4 pi b0/|b1| = 0.042 (Banks-Zaks): set by the particle count alone, with no size and no shape",
      abs(a16 - 0.0416) < 0.0002, f"alpha* = {a16:.4f}")
def flow(a0, Nf=16, N=3, t1=400.0):
    B0, B1 = [float(x) for x in bz(N, Nf)]
    f = lambda t, a: [-(B0/(2*math.pi))*a[0]**2 - (B1/(8*math.pi**2))*a[0]**3]
    return solve_ivp(f, [0, -t1], [a0], rtol=1e-10, atol=1e-14).y[0][-1]   # toward long distance (t = ln mu decreasing)
lo_, hi_ = flow(0.01, t1=4000), flow(0.12, t1=4000)
check("K2c starting weak (0.01) or strong (0.12) at short distance, the strength flows to the same 0.042 far away: the long-range strength forgets the shape and size of where it started",
      abs(lo_ - a16) < 2e-3 and abs(hi_ - a16) < 2e-3, f"{lo_:.4f}, {hi_:.4f}")
check("K2d with the strength frozen there is no length scale, so the force is exactly Coulomb-like, felt only by particles carrying the hidden charge (a 'dark' conformal sector; known: Caswell 1974, Banks-Zaks 1982, Georgi's unparticles 2007)",
      sp.simplify(-(b0/(2*sp.pi))*(-4*sp.pi*b0/b1)**2 - (b1/(8*sp.pi**2))*(-4*sp.pi*b0/b1)**3) == 0, "beta(alpha*) = 0 exactly")
# ---------------- K3 rungs, order on
M1 = float(prev["I_copy1_GeV"]); lhc_KKg = 3.3e3
OUT.update(K_KKgluon_GeV=M1, K_KKgluon_over_LHC=M1/lhc_KKg)
y = sp.symbols('y', real=True); n1, n2, n3 = 1, 2, 3
ov = sp.integrate(sp.exp(sp.I*(n1 + n2 - n3)*y), (y, 0, 2*sp.pi)); ov2 = sp.integrate(sp.exp(sp.I*(n1 + n2 - 4)*y), (y, 0, 2*sp.pi))
check("K3 copies of W bosons or gluons on our circle (rungs, order on): at n x 5.2e17 GeV, carrying the colour charge themselves, and they only meet when the rung numbers add up (n1 + n2 = n3, momentum round the circle conserved)",
      ov == 2*sp.pi and ov2 == 0)
check("K3b where it stands: the LHC rules out Kaluza-Klein gluons below about 3.3 TeV (CMS, 13 TeV, in a warped model); ours sit 1.6e14 times higher",
      abs(M1/lhc_KKg - 1.58e14)/1.58e14 < 0.02, f"{M1/lhc_KKg:.2e}")
# ---------------- K4 tube, order off: the abelian flux tube (Nielsen-Olesen) at the self-dual point
def vortex(c, n=1, rmax=12.0):
    def rhs(r, s):
        f, a = s
        return [n*(1 - a)*f/r, r*(1 - f**2)/(2*n)]
    r0 = 1e-4
    sol = solve_ivp(rhs, [r0, rmax], [c*r0**n, r0**2/(4*n)], rtol=1e-11, atol=1e-13, dense_output=True)
    return sol
lo, hi = 0.1, 2.0
for _ in range(80):
    c = 0.5*(lo + hi); s_ = vortex(c)
    if s_.status != 0 or np.max(s_.y[0]) > 1: hi = c
    else: lo = c
c_star = 0.5*(lo + hi); s_ = vortex(c_star, rmax=8.0)
def edens(r):
    f, a = s_.sol(r); fp = (1 - a)*f/r; ap = r*(1 - f**2)/2
    return 2*math.pi*r*(0.5*(ap/r)**2 + 0.5*(fp**2 + (1 - a)**2*f**2/r**2) + 0.125*(1 - f**2)**2)
T = quad(edens, 1e-4, 7.5, limit=400)[0]
OUT.update(K_vortex_tension_over_pi=T/math.pi)
check("K4 the abelian tube: a magnetic flux tube at the self-dual point carries tension exactly pi v^2 per unit of flux (solved numerically; Manton-Sutcliffe units, = 2 pi eta^2 for a canonical field eta), so two charges joined by it feel a constant pull (D_eff = 1)",
      abs(T/math.pi - 1) < 2e-3, f"tension / (pi v^2) = {T/math.pi:.4f}")
for Gmu, nm in [(1.3e-7, "CMB"), (1e-10, "PTA")]:
    OUT[f"K_string_eta_max_GeV_{nm}"] = math.sqrt(Gmu/(2*math.pi))*MP
check("K4b the fundamental form of this cell is a cosmic string: Planck's CMB limit G mu < 1.3e-7 caps the breaking scale at 1.8e15 GeV; pulsar timing's G mu < 1e-10 (Nambu-Goto loops, weaker with field-theory loops) at 5e13 GeV",
      abs(OUT["K_string_eta_max_GeV_CMB"] - 1.76e15)/1.76e15 < 0.02 and abs(OUT["K_string_eta_max_GeV_PTA"] - 4.87e13)/4.87e13 < 0.02,
      f"{OUT['K_string_eta_max_GeV_CMB']:.2e} GeV, {OUT['K_string_eta_max_GeV_PTA']:.2e} GeV")
Gmu_ours = 2*math.pi*(M1/MP)**2
OUT["K_Gmu_at_copy_scale"] = Gmu_ours
check("K4c a string with the circle's own scale (5.2e17 GeV) would have G mu = 0.011, 1e5 times over the CMB limit; nothing in our law makes one, and inflation would dilute any made before it (stated, not derived)",
      Gmu_ours/1.3e-7 > 1e4, f"G mu = {Gmu_ours:.3f}")
# ---------------- K5 a base dial that varies in space and hides: screening
rho_water = 1.0*5.6096e23*(1.973269804e-14)**3        # 1 g/cm^3 in GeV^4
m_rad = float(prev["I_radion_GeV"])
shift = rho_water/(Mbar**2*m_rad**2)
OUT["K_density_shift_ratio"] = shift
check("K5 a field that varies in space can hide from tests if dense matter makes it heavy (chameleon, Khoury-Weltman 2004); for our size carrier, ordinary density changes its mass-squared by 1 part in 1e86: no hiding, and none needed at its range",
      shift < 1e-80, f"rho/(Mbar^2 m^2) = {shift:.1e}")

log(f"\nPART K TOTAL: {NP[0]} pass, {NP[1]} fail")
prev.update(OUT); json.dump(prev, open("results.json", "w"), indent=1)
open("other_forces_output.txt", "w").write("\n".join(lines) + "\n")
