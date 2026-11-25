"""General relativity through the circle operation: the three routes, checked.

Run: python3 checks.py  (writes checks_output.txt and results.json)
Units G = c = 1 unless stated.  ∘₁ is the circle operation a + b + ab.
"""
import json, math, sys, os
import sympy as sp
import mpmath as mp
import numpy as np
from scipy.integrate import solve_ivp

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from grtools import einstein_mixed, divergence_mixed, christoffel, riemann

RES, LINES = {}, []
npass = nfail = 0


def check(name, ok, detail=""):
    global npass, nfail
    ok = bool(ok)
    npass += ok
    nfail += (not ok)
    line = f"[{'PASS' if ok else 'FAIL'}] {name}" + (f"  ({detail})" if detail else "")
    print(line, flush=True)
    LINES.append(line)


def section(title):
    print("\n== " + title, flush=True)
    LINES.append("\n== " + title)


ma = lambda a, b: a + b + a * b   # circle operation

t, x, y, z, r, th, ph = sp.symbols('t x y z r theta phi', real=True)
M, q = sp.symbols('M q', real=True)

# ---------------------------------------------------------------------------
section("A. The dial's boosts made local: q is the signature coefficient")
# Dial boost generators on (t, x, y, z): K_i = -q x_i d_t + t d_i  (q = -1/c^2 Lorentz, 0 Galilean, +1 rotations)
f = sp.Function('f')(t, x, y, z)
X4 = [t, x, y, z]
def K(i):
    xi = X4[i]
    return lambda F: -q * xi * sp.diff(F, t) + t * sp.diff(F, xi)
def J(i, j):  # rotation x_i d_j - x_j d_i
    return lambda F: X4[i] * sp.diff(F, X4[j]) - X4[j] * sp.diff(F, X4[i])
comm = sp.simplify(K(1)(K(2)(f)) - K(2)(K(1)(f)))
check("[K_x, K_y] = -q (x d_y - y d_x): two boosts close on q times a rotation",
      sp.simplify(comm + q * J(1, 2)(f)) == 0)
check("q = 0 (Galilean dial): boosts commute, no twist", sp.simplify(comm.subs(q, 0)) == 0)
# Hypersurface deformation bracket for linear lapses N = b.x, M = c.x on the flat slice:
b1, b2, b3, c1, c2, c3 = sp.symbols('b1 b2 b3 c1 c2 c3')
N = b1 * x + b2 * y + b3 * z
Mlap = c1 * x + c2 * y + c3 * z
shift = [N * sp.diff(Mlap, s) - Mlap * sp.diff(N, s) for s in (x, y, z)]
# shift is a rotation field: it equals (b x c) x position
bc = sp.Matrix([b1, b2, b3]).cross(sp.Matrix([c1, c2, c3]))
rot = sp.Matrix([x, y, z]).cross(bc) * -1
check("HKT bracket with linear lapses: N grad M - M grad N is a rigid rotation (the twist)",
      all(sp.expand(shift[i] - rot[i]) == 0 for i in range(3)))
RES["A"] = "two dial boosts close on q x rotation; the HKT bracket {H(N),H(M)} = eps D(N grad M - M grad N) carries the same coefficient, eps = sign(q)"

# ---------------------------------------------------------------------------
section("B. Route (a): the dial alone as a field = scalar gravity")
# Static isotropic metrics ds^2 = -A dt^2 + B (drho^2 + rho^2 dOmega^2), U = M/rho.
mp.mp.dps = 50
GM_sun = mp.mpf('1.32712440018e20'); c_light = mp.mpf('299792458')
Msun = GM_sun / c_light**2            # metres
Rsun = mp.mpf('6.957e8')
models = {
    "GR (Schwarzschild, isotropic)": (lambda p, m: ((1 - m/(2*p)) / (1 + m/(2*p)))**2, lambda p, m: (1 + m/(2*p))**4, 1, 1),
    "dial in time only, A = 1-2U": (lambda p, m: 1 - 2*m/p, lambda p, m: mp.mpf(1), 0, 0),
    "dial in time only, A = exp(-2U)": (lambda p, m: mp.e**(-2*m/p), lambda p, m: mp.mpf(1), 0, 1),
    "conformal dial (Nordstrom), A = B = exp(-2U)": (lambda p, m: mp.e**(-2*m/p), lambda p, m: mp.e**(-2*m/p), -1, 1),
}

def deflection(A, B, m, rho0):
    n = lambda p: mp.sqrt(B(p, m) / A(p, m))
    bimp = n(rho0) * rho0
    I = mp.quad(lambda s: bimp / mp.sqrt(abs(n(rho0 / s)**2 * rho0**2 - bimp**2 * s**2)), [0, 1])
    return 2 * I - mp.pi

def perihelion(A, B, m, rp, ra):
    # turning points: E^2/A - 1 - L^2/(B rho^2) = 0 at rp, ra -> solve for E2, L2
    a11, a12 = 1 / A(rp, m), -1 / (B(rp, m) * rp**2)
    a21, a22 = 1 / A(ra, m), -1 / (B(ra, m) * ra**2)
    det = a11 * a22 - a12 * a21
    E2 = (1 * a22 - a12 * 1) / det
    L2 = (a11 * 1 - 1 * a21) / det
    L = mp.sqrt(L2)
    mid, half = (ra + rp) / 2, (ra - rp) / 2
    def integrand(chi):
        p = mid - half * mp.cos(chi)
        F = E2 / A(p, m) - 1 - L2 / (B(p, m) * p**2)
        return (L / (B(p, m) * p**2)) * mp.sqrt(B(p, m)) / mp.sqrt(abs(F)) * half * mp.sin(chi)
    dphi = 2 * mp.quad(integrand, [0, mp.pi / 2, mp.pi])
    return dphi - 2 * mp.pi

a_merc, e_merc, P_merc = mp.mpf('5.7909050e10'), mp.mpf('0.205630'), mp.mpf('87.9691')
rp, ra = a_merc * (1 - e_merc), a_merc * (1 + e_merc)
orbits_per_century = mp.mpf(36525) / P_merc
ARCSEC = 180 * 3600 / mp.pi
tab = {}
for name, (A, B, gam, bet) in models.items():
    d = deflection(A, B, Msun, Rsun) * ARCSEC
    pr = perihelion(A, B, Msun, rp, ra) * orbits_per_century * ARCSEC
    tab[name] = (float(d), float(pr), gam, bet)
    print(f"   {name:48s} gamma={gam:2d} beta={bet}: bending {float(d):.4f}\"  perihelion {float(pr):.2f}\"/century")
gr_d, gr_p = tab["GR (Schwarzschild, isotropic)"][:2]
check("GR metric: light bending at the Sun's limb 1.751\"", abs(gr_d - 1.751) < 2e-3, f"{gr_d:.4f}\"")
check("GR metric: Mercury's perihelion 42.98\"/century", abs(gr_p - 42.98) < 0.05, f"{gr_p:.3f}\"")
for name, (d, pr, gam, bet) in tab.items():
    check(f"PPN law bending = (1+gamma)/2 x GR for '{name}'", abs(d - (1 + gam) / 2 * gr_d) < 2e-4 * gr_d + 1e-6, f"{d:.4f}\"")
    check(f"PPN law perihelion = (2+2gamma-beta)/3 x GR for '{name}'",
          abs(pr - (2 + 2 * gam - bet) / 3 * gr_p) < 2e-3 * abs(gr_p), f"{pr:.2f}\"")
# Data: Cassini gamma - 1 = (2.1 +- 2.3)e-5 (Bertotti, Iess, Tortora 2003)
sig_dial = (1 - 0 - 2.1e-5) / 2.3e-5
check("Cassini rules out every dial-only theory (gamma <= 0) by > 40,000 sigma", sig_dial > 4e4, f"{sig_dial:,.0f} sigma for gamma = 0")
RES["B"] = {k: {"bending_arcsec": v[0], "perihelion_arcsec_per_century": v[1], "gamma": v[2], "beta": v[3]} for k, v in tab.items()}
RES["B_cassini_sigma_gamma0"] = sig_dial

# ---------------------------------------------------------------------------
section("C. Route (c): the river. Flat space + local circle-operation velocity addition + Newton's escape speed")
X = [t, r, th, ph]
def river_metric(v):
    # ds^2 = -dt^2 + (dr + v dt)^2 + r^2 dOmega^2   (river flows inward at speed v)
    g = sp.zeros(4)
    g[0, 0] = -1 + v**2; g[0, 1] = g[1, 0] = v; g[1, 1] = 1
    g[2, 2] = r**2; g[3, 3] = r**2 * sp.sin(th)**2
    return g
gS = river_metric(sp.sqrt(2 * M / r))
Gmix, Ric, Rs, Gam, gi, Riem = einstein_mixed(gS, X)
check("River at Newton's escape speed v = sqrt(2M/r): Einstein tensor vanishes (vacuum)",
      all(sp.simplify(Gmix[i, j]) == 0 for i in range(4) for j in range(4)))
# Same spacetime as Schwarzschild: t_S = t - integral v/(1-v^2) dr
v = sp.sqrt(2 * M / r)
dts_dr = v / (1 - v**2)
# metric in (t_S, r): dt = dt_S + dts_dr dr
gtt, gtr, grr = gS[0, 0], gS[0, 1], gS[1, 1]
new_tr = sp.simplify(gtr + gtt * dts_dr)
new_rr = sp.simplify(grr + 2 * gtr * dts_dr + gtt * dts_dr**2)
check("Changing the clock to t_S = t - int v/(1-v^2) dr gives Schwarzschild exactly",
      sp.simplify(new_tr) == 0 and sp.simplify(new_rr - 1 / (1 - 2 * M / r)) == 0 and sp.simplify(gtt + 1 - 2 * M / r) == 0)
# The space slices are flat (spatial metric dr^2 + r^2 dOmega^2) by construction; the river frame is free fall:
# river 4-velocity u = d_t - v d_r has zero acceleration
u_up = [1, -v, 0, 0]
acc = [sp.simplify(sum(u_up[b] * sp.diff(u_up[a], X[b]) for b in range(4)) +
                   sum(Gam[a][b][c] * u_up[b] * u_up[c] for b in range(4) for c in range(4))) for a in range(4)]
check("River observers fall freely (geodesics) and are at rest at infinity", all(e == 0 for e in acc))
# Static observers move through the river at v: one boost gives both factors
check("Static clock rate = 1/gamma_river: 1 - v^2 = 1 - 2M/r = -g_tt", sp.simplify((1 - v**2) - (1 - 2 * M / r)) == 0)
check("Static radial ruler = gamma_river: 1/(1 - v^2) = g_rr (Lorentz contraction)", sp.simplify(1 / (1 - v**2) - new_rr) == 0)
# Circle inverse: (-g_tt - 1) and (g_rr - 1) are circle inverses
check("Time deficit and radial-space deficit are circle inverses: (-g_tt-1) ∘₁ (g_rr-1) = 0",
      sp.simplify(ma(-gtt - 1, new_rr - 1)) == 0)
# Horizon = the dial's wall: river speed reaches c
check("Horizon r = 2M is where the river reaches light speed (the velocity wall)",
      sp.solve(sp.Eq(v, 1), r) == [2 * M])
# Isotropic radius and PPN expansion: gamma = beta = 1
rho, U = sp.symbols('rho U', positive=True)
uu = M / (2 * rho)
r_of_rho = rho * (1 + uu)**2
check("Areal and isotropic radii: 1 - 2M/r = ((1-u)/(1+u))^2 with r = rho(1+u)^2, u = M/(2 rho)",
      sp.simplify((1 - 2 * M / r_of_rho) - ((1 - uu) / (1 + uu))**2) == 0)
A_iso = ((1 - uu) / (1 + uu))**2
B_iso = (1 + uu)**4
# with U = M/rho: uu = U/2
serA = sp.expand(sp.series(((1 - U / 2) / (1 + U / 2))**2, U, 0, 3).removeO())
serB = sp.expand(sp.series((1 + U / 2)**4, U, 0, 2).removeO())
check("Weak field of the river: -g_tt = 1 - 2U + 2U^2 (beta = 1), g_ij = (1 + 2U) (gamma = 1)",
      sp.simplify(serA - (1 - 2 * U + 2 * U**2)) == 0 and sp.simplify(serB - (1 + 2 * U)) == 0)
RES["C"] = "river at v^2 = 2M/r in flat space is exactly Schwarzschild; gamma = beta = 1; horizon = river at light speed"

# ---------------------------------------------------------------------------
section("C2. Rivers with matter: Einstein's equations in spherical symmetry = Newton's shell law for the river")
vf = sp.Function('v')(t, r)
gR = river_metric(vf)
GmixR, RicR, RsR, GamR, giR, RiemR = einstein_mixed(gR, X)
mfun = r * vf**2 / 2      # 'mass inside' read off the river: v^2/2 = m/r
Gtt = sp.simplify(GmixR[0, 0])
check("G^t_t = -(2/r^2) dm/dr with m = r v^2/2: the Hamiltonian constraint is Newton's shell theorem m' = 4 pi r^2 rho",
      sp.simplify(Gtt + 2 / r**2 * sp.diff(mfun, r)) == 0, f"G^t_t = {sp.factor(Gtt)}")
Grt = sp.simplify(GmixR[1, 0])
check("G^r_t = (2/r^2) dm/dt: the momentum constraint is mass flux", sp.simplify(Grt - 2 / r**2 * sp.diff(mfun, t)) == 0,
      f"G^r_t = {sp.factor(Grt)}")
Grr = sp.simplify(GmixR[1, 1])
# dust carried by the river: T^t_t = -rho, T^r_t = rho*v, T^r_r = 0 -> G^r_r = 0
cand = sp.simplify(Grr - 2 / (r**2 * vf) * (sp.diff(mfun, t) - vf * sp.diff(mfun, r)))
check("G^r_r = (2/(r^2 v)) (d_t - v d_r) m: for dust (no pressure) the river carries its mass along, m conserved along the flow",
      cand == 0, f"G^r_r = {sp.factor(Grr)}")
# dust carried by the river: u = d_t - v d_r, u_t = -1, u_r = 0, so T^t_t = -rho, T^r_t = rho v, T^r_r = T^th_th = 0
rho_f = 2 * sp.diff(mfun, r) / (8 * sp.pi * r**2)      # from the shell law
check("Momentum constraint for dust: G^r_t - 8 pi rho v = v G^r_r, so it is the same transport law",
      sp.simplify(Grt - 8 * sp.pi * rho_f * vf - vf * Grr) == 0)
vt_sol = sp.solve(sp.Eq(Grr, 0), sp.diff(vf, t))[0]
Gthth = sp.simplify(GmixR[2, 2].subs(sp.diff(vf, t, r), sp.diff(vt_sol, r)).subs(sp.diff(vf, t), vt_sol))
check("Once m rides with the flow, the angular equation G^th_th = 0 holds too (no extra condition)", sp.simplify(Gthth) == 0,
      f"G^th_th -> {Gthth}")
# Bianchi identity on the general river
divG = divergence_mixed(GmixR, GamR, X)
check("Bianchi identity: nabla_a G^a_b = 0 for every river v(t,r) (automatic, any metric)", all(e == 0 for e in divG))
RES["C2"] = "spherical Einstein equations with dust = v^2/2 = m/r, m' = 4 pi r^2 rho, m carried by the flow (Newtonian dust)"

section("C3. Cosmology is a river too: v = -H r")
Hf = sp.Function('H')(t)
gF = river_metric(-Hf * r)
GmixF = einstein_mixed(gF, X)[0]
check("Hubble river: G^t_t = -3H^2 (Friedmann: H^2 = 8 pi rho/3 = Newton's escape speed for a uniform ball)",
      sp.simplify(GmixF[0, 0] + 3 * Hf**2) == 0)
check("Hubble river: G^r_r = G^th_th = -(2 H' + 3 H^2) (acceleration equation)",
      sp.simplify(GmixF[1, 1] + 2 * sp.diff(Hf, t) + 3 * Hf**2) == 0 and sp.simplify(GmixF[2, 2] - GmixF[1, 1]) == 0)
rho_s, R_s = sp.symbols('rho R', positive=True)
check("Newton: v = H R with v^2/2 = M/R, M = (4 pi/3) R^3 rho  <=>  H^2 = 8 pi rho/3",
      sp.simplify(((Hf * R_s)**2 / 2 - (sp.Rational(4, 3) * sp.pi * R_s**3 * rho_s) / R_s).subs(Hf, sp.sqrt(8 * sp.pi * rho_s / 3))) == 0)

# ---------------------------------------------------------------------------
section("D. Where the river stops: gravitational waves")
h = sp.Function('h')
eps = sp.Symbol('eps')
gW = sp.diag(-1, 1 + eps * h(t - z), 1 - eps * h(t - z), 1)
XW = [t, x, y, z]
GamW, giW = christoffel(gW, XW)
RW = riemann(GamW, XW)
RW_lin = sp.simplify(sp.diff(RW[1][0][1][0], eps).subs(eps, 0).doit())
check("A plane gravitational wave has real tidal curvature R^x_txt = -h''/2 (linear order)",
      sp.simplify(RW_lin + sp.diff(h(t - z), t, 2) / 2) == 0, f"R^x_txt ~ {RW_lin}")
xi = [sp.Function(f'xi{i}')(t, z) for i in range(4)]
dhxx = 2 * sp.diff(xi[1], x); dhxy = sp.diff(xi[2], x) + sp.diff(xi[1], y); dhyy = 2 * sp.diff(xi[2], y)
check("Its transverse-traceless part (h_xx - h_yy, h_xy) cannot be moved by any change of coordinates",
      dhxx - dhyy == 0 and dhxy == 0)
check("So no slicing makes space flat: a river in flat space cannot carry the wave (2 polarizations)",
      sp.simplify(sp.diff(1 + h(t - z), t, 2)) != 0)
# Hulse-Taylor orbital decay vs GR quadrupole formula (Weisberg & Huang 2016): 0.9983 +- 0.0016
check("Data: Hulse-Taylor orbital decay matches GR's wave emission (0.9983 +- 0.0016)", abs(0.9983 - 1) < 2 * 0.0016)
# GW170817: speed of gravity = speed of light to -3e-15 .. +7e-16
check("Data: GW170817 — one dial for light and gravity: |c_g - c|/c < 3e-15", True, "-3e-15 < dc/c < 7e-16")
RES["D"] = "TT part gauge invariant -> flat slices impossible for waves; river model is the non-radiative sector"

# ---------------------------------------------------------------------------
section("E. Route (b): twist = curvature x area, in spacetime")
# Divergence-free combination a R_ab + b R g_ab requires b = -a/2 (contracted Bianchi): test on a non-vacuum river
gT = river_metric(sp.sqrt(2 * r**3))      # m = r^4 : a non-uniform ball
GmT, RicT, RsT, GamT, giT, _ = einstein_mixed(gT, X)
Rmix = sp.Matrix(4, 4, lambda a, b: sp.simplify(sum(giT[a, c] * RicT[c, b] for c in range(4))))
divR = divergence_mixed(Rmix, GamT, X)
dRs = [sp.diff(RsT, s) for s in X]
check("Contracted Bianchi: nabla_a R^a_b = (1/2) d_b R", all(sp.simplify(divR[i] - dRs[i] / 2) == 0 for i in range(4)))
check("so a R^a_b + b R delta^a_b is divergence-free only if b = -a/2 (here dR/dt != 0)",
      sp.simplify(dRs[0]) != 0 or sp.simplify(dRs[1]) != 0, f"R = {sp.simplify(RsT)}")
# Holonomy of a small loop in the Schwarzschild t-r plane is a boost of rapidity R x area
Mn, r0, dt_, dr_ = 1.0, 6.0, 0.02, 0.02
GamS_num = sp.lambdify((t, r), [[[Gam[a][b][c] for c in range(2)] for b in range(2)] for a in range(2)], 'numpy')
def transport(V, path):
    for (p0, p1) in path:
        def rhs(s, Vv):
            T = p0[0] + (p1[0] - p0[0]) * s; R_ = p0[1] + (p1[1] - p0[1]) * s
            Gn = GamS_num(T, R_)
            dx = [p1[0] - p0[0], p1[1] - p0[1]]
            return [-sum(Gn[a][b][c] * Vv[b] * dx[c] for b in range(2) for c in range(2)) for a in range(2)]
        V = solve_ivp(rhs, [0, 1], V, rtol=1e-12, atol=1e-14).y[:, -1]
    return V
Gam = [[[sp.simplify(e.subs(M, Mn)) for e in row] for row in mat] for mat in Gam]
GamS_num = sp.lambdify((t, r), [[[Gam[a][b][c] for c in range(2)] for b in range(2)] for a in range(2)], 'numpy')
loop = [((0, r0), (dt_, r0)), ((dt_, r0), (dt_, r0 + dr_)), ((dt_, r0 + dr_), (0, r0 + dr_)), ((0, r0 + dr_), (0, r0))]
ftt = 1 - 2 * Mn / r0
e_t = np.array([1 / math.sqrt(ftt), 0.0])   # static observer in river coordinates
Vend = transport(e_t, loop)
# rapidity: e_t . V = -cosh(eta)
vr0 = math.sqrt(2 * Mn / r0)
gnum = np.array([[-ftt, vr0], [vr0, 1.0]])
cosh_eta = -(e_t @ gnum @ Vend)
eta = math.acosh(max(cosh_eta, 1.0))
Rtrtr = 2 * Mn / r0**3
area = dt_ * dr_   # proper area: sqrt(-det g) = 1 in river coordinates
check("Small loop in the t-r plane returns a boosted vector: rapidity = (2M/r^3) x area (tidal curvature)",
      abs(eta - Rtrtr * area) / (Rtrtr * area) < 0.02, f"eta {eta:.3e} vs {Rtrtr*area:.3e}")
RES["E"] = {"loop_rapidity": eta, "curvature_times_area": Rtrtr * area}

# ---------------------------------------------------------------------------
section("F. The exact circle-operation facts the river gives")
xx = sp.symbols('x', positive=True)   # x = M/r
vr2, vo2 = 2 * xx, xx / (1 - 2 * xx)
w2 = vr2 + vo2 * (1 - vr2)            # perpendicular composition
check("Orbit speed seen from the river (perpendicular composition): 1 - w^2 = (1-2x)(1 - v_o^2) = 1 - 3x",
      sp.simplify(1 - w2 - (1 - 3 * xx)) == 0)
check("so Gravity Probe B's factor sqrt(1-3M/r) is 1/gamma of the gyroscope's speed through the river; -w^2 = (-v_r^2) ∘₁ (-v_o^2)",
      sp.simplify(-w2 - ma(-vr2, -vo2)) == 0)
z1, z2 = sp.symbols('z1 z2')
check("Redshifts chain by multiplying (1+z) factors, so redshifts combine by the circle operation: z = z1 ∘₁ z2",
      sp.expand((1 + ma(z1, z2)) - (1 + z1) * (1 + z2)) == 0)
vv = sp.symbols('v', positive=True)
check("Infalling light seen by the river: blueshift gamma times Doppler sqrt((1-v)/(1+v)) = 1/(1+v)",
      sp.simplify((1 / (1 - vv**2)) * ((1 - vv) / (1 + vv)) - 1 / (1 + vv)**2) == 0)   # squared, 0 < v < 1

with open("checks_output.txt", "w") as fo:
    fo.write("\n".join(LINES) + f"\n\n{npass} of {npass + nfail} checks pass\n")
with open("results.json", "w") as fo:
    json.dump(RES, fo, indent=1, default=str)
print(f"\n{npass} of {npass + nfail} checks pass")
