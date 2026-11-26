"""Does one dial set the reach of every force? Checks for the DRAFT page.

Run: python3 -I checks.py  ->  checks_output.txt, results.json
Each check prints PASS/FAIL. Part P is the polishing step: key numbers again at
30-digit precision, reporting how far they moved.

Parts (fundamental to derived):
 A  the law's own linear form: the generator v = 1 + kx + dx^2 is a Riccati flow;
    linearised it is w'' = (Delta/4) w, the Jacobi equation with curvature K = -Delta/4
 B  carrier mass: m^2 = Delta/4 = minus the curvature reading; Yukawa = flux x screening
 C  circle class: a closed direction gives a massless zero mode plus massive rungs n sqrt(-Delta)
 D  confinement: flux that keeps one open direction (a tube) gives a constant pull;
    no single (k, d) makes a tube from a point source
 E  the dynamical tube (London / Abrikosov): sideways screening with the same mass m
    sets the tube's width (2/m) and its tension (~ m^2 ln)
 F  the forces placed: reach = hbar c / (m c^2) from measured masses
 P  polishing
"""
import json, math
import numpy as np
from scipy import integrate, special
import sympy as sp
import mpmath as mp

OUT, RES = [], {}
npass = nfail = 0


def check(name, ok, detail=""):
    global npass, nfail
    npass += bool(ok)
    nfail += (not ok)
    line = f"{'PASS' if ok else 'FAIL'}  {name}" + (f"   [{detail}]" if detail else "")
    print(line)
    OUT.append(line)


def say(s):
    print(s)
    OUT.append(s)


# ------------------------------------------------------------------ A
say("== A. The law's own linear form")
L, k, d, K, m, r, s, q = sp.symbols('L k d K m r s q', real=True)
u = sp.Function('u')
w = sp.Function('w')
# A1: x = -u'/(d u) with u'' - k u' + d u = 0 gives x' = 1 + k x + d x^2
x = -sp.diff(u(L), L) / (d * u(L))
xp = sp.diff(x, L).subs(sp.diff(u(L), L, 2), k * sp.diff(u(L), L) - d * u(L))
check("A1 Riccati: x = -u'/(d u), u'' - k u' + d u = 0  =>  x' = 1 + k x + d x^2",
      sp.simplify(xp - (1 + k * x + d * x ** 2)) == 0)
# A2: u = e^{kL/2} w removes the drift: w'' = (Delta/4) w
uu = sp.exp(k * L / 2) * w(L)
expr = sp.diff(uu, L, 2) - k * sp.diff(uu, L) + d * uu
w2 = sp.solve(sp.Eq(expr, 0), sp.diff(w(L), L, 2))[0]
check("A2 drift removed: w'' = (k^2 - 4d)/4 * w = (Delta/4) w",
      sp.simplify(w2 - (k ** 2 - 4 * d) / 4 * w(L)) == 0)


# A3-A5: the three classes, numerically, from the flow x' = v(x)
def flow(kk, dd, x0, Lmax, n=4001):
    sol = integrate.solve_ivp(lambda t, y: [1 + kk * y[0] + dd * y[0] ** 2], [0, Lmax], [x0],
                              t_eval=np.linspace(0, Lmax, n), rtol=1e-12, atol=1e-14)
    return sol.t, sol.y[0]


# line class: k=3, d=1, Delta=5; walls (-3 +- sqrt5)/2; flow runs up to the upper wall (attracting)
kk, dd = 3.0, 1.0
Dl = kk * kk - 4 * dd
wall_hi = (-kk + math.sqrt(Dl)) / (2 * dd)
wall_lo = (-kk - math.sqrt(Dl)) / (2 * dd)
t, xs = flow(kk, dd, -1.0, 12.0)
# x' = v < 0 between the walls, so x falls toward wall_lo
dist = np.abs(xs - wall_lo)
sel = (t > 4) & (t < 10)
rate = -np.polyfit(t[sel], np.log(dist[sel]), 1)[0]
check("A3 line class (k=3, d=1): flow reaches its wall exponentially, rate sqrt(Delta)",
      abs(rate - math.sqrt(Dl)) < 1e-4, f"rate {rate:.6f} vs sqrt5 = {math.sqrt(Dl):.6f}")
RES['A3_rate'] = rate
# border: k=2, d=1, Delta=0; double wall at -1; approach like 1/L (a power, no length scale)
t, xs = flow(2.0, 1.0, -3.0, 400.0, 8001)
dist = np.abs(xs + 1)
sel = (t > 50) & (t < 400)
slope = np.polyfit(np.log(t[sel]), np.log(dist[sel]), 1)[0]
check("A4 border (k=2, d=1, Delta=0): approach is a power law, slope -1 (no reach scale)",
      abs(slope + 1) < 0.01, f"log-log slope {slope:.5f}")
# circle class: k=1, d=1, Delta=-3: ruler closes, period 2 pi / sqrt(3)
P = integrate.quad(lambda z: 1 / (1 + z + z * z), -np.inf, np.inf)[0]
check("A5 circle class (k=1, d=1): ruler closes with length 2 pi / sqrt(-Delta)",
      abs(P - 2 * math.pi / math.sqrt(3)) < 1e-9, f"{P:.12f}")
# A6 Jacobi / Raychaudhuri: spread R of field lines: R'' = -K R. Expansion theta = R'/R obeys
# theta' = -K - theta^2. With x = -theta/K this is x' = 1 + K x^2: our law with k = 0, d = K,
# whose curvature reading -Delta/4 = -(0 - 4K)/4 = K.
th = sp.Function('th')
xx = -th(L) / K
lhs = sp.diff(xx, L).subs(sp.diff(th(L), L), -K - th(L) ** 2)
check("A6 Raychaudhuri theta' = -K - theta^2 is the law x' = 1 + K x^2 (k=0, d=K)",
      sp.simplify(lhs - (1 + K * xx ** 2)) == 0)
Dk = 0 ** 2 - 4 * K
check("A6b its curvature reading -Delta/4 equals the Jacobi curvature K", sp.simplify(-Dk / 4 - K) == 0)
# A7 the three spreads from a point (R(0)=0, R'(0)=1)
Rpos = sp.sin(sp.sqrt(K) * r) / sp.sqrt(K)
Rneg = sp.sinh(sp.sqrt(-K) * r) / sp.sqrt(-K)
ok = all(sp.simplify(sp.diff(R, r, 2) + K * R) == 0 for R in (Rpos, Rneg))
check("A7 point-source spreads: sin (K>0), r (K=0), sinh (K<0) solve R'' = -K R", ok)

# ------------------------------------------------------------------ B
say("== B. Carrier mass is a reading of the same dials")
# B1: Yukawa in 3D: phi = e^{-mr}/(4 pi r) solves (nabla^2 - m^2) phi = 0 for r > 0
phi = sp.exp(-m * r) / (4 * sp.pi * r)
lap = sp.diff(r ** 2 * sp.diff(phi, r), r) / r ** 2
check("B1 e^{-mr}/(4 pi r) solves nabla^2 phi = m^2 phi (r>0); reduced u = r phi obeys u'' = m^2 u",
      sp.simplify(lap - m ** 2 * phi) == 0)
# B2: Fourier: G(r) = (1/(2 pi^2 r)) int_0^inf q sin(qr)/(q^2+m^2) dq = e^{-mr}/(4 pi r)
mp.mp.dps = 20
errs = []
for mm, rr in [(1.0, 0.5), (1.0, 2.0), (0.3, 3.0), (2.5, 0.7)]:
    I = mp.quadosc(lambda qq: qq * mp.sin(qq * rr) / (qq ** 2 + mm ** 2), [0, mp.inf], omega=rr)
    G = I / (2 * mp.pi ** 2 * rr)
    errs.append(abs(G - mp.e ** (-mm * rr) / (4 * mp.pi * rr)) / (mp.e ** (-mm * rr) / (4 * mp.pi * rr)))
check("B2 Fourier of 1/(q^2 + m^2) is e^{-mr}/(4 pi r) (4 settings)", max(errs) < 1e-12,
      f"max rel err {float(max(errs)):.1e}")


# B3: D dimensions: G_D = (2pi)^{-D/2} (m/r)^{D/2-1} K_{D/2-1}(mr); m -> 0 gives the dimension-dial result
def GD(D, mm, rr):
    nu = D / 2 - 1
    return (2 * math.pi) ** (-D / 2) * (mm / rr) ** nu * special.kv(nu, mm * rr)


def Gflux(D, rr):
    return math.gamma(D / 2 - 1) / (4 * math.pi ** (D / 2) * rr ** (D - 2))


ok3 = abs(GD(3, 0.8, 1.7) - math.exp(-0.8 * 1.7) / (4 * math.pi * 1.7)) < 1e-14
lim = [abs(GD(D, 1e-12, 1.3) / Gflux(D, 1.3) - 1) for D in (2.5, 3, 3.5, 4, 5)]
check("B3 D-dim massive Green's function: exact e^{-mr}/(4 pi r) at D=3; m->0 gives Gamma(D/2-1)/(4 pi^{D/2} r^{D-2})",
      ok3 and max(lim) < 1e-5, f"max rel dev at m=1e-12: {max(lim):.1e} (approach is slow, like (mr)^(D-2), at D=2.5)")
# B3b: far away every D decays at the same rate m: d ln G / dr -> -m
rates = []
for D in (2.5, 3.0, 4.0, 5.0):
    r1, r2 = 40.0, 41.0
    rates.append(-(math.log(GD(D, 1.0, r2)) - math.log(GD(D, 1.0, r1))) - (D - 1) / 2 * math.log(r2 / r1))
check("B3b far decay rate is m in every dimension (after the r^{-(D-1)/2} power)",
      max(abs(x - 1) for x in rates) < 5e-3, ", ".join(f"{x:.4f}" for x in rates))
# B4 3D: log slope adds: d ln G/dr = -1/r (flux, from D) - m (screening, from Delta)
check("B4 3D: d ln G / dr = -1/r - m exactly: the power (dimension dial) and the screening add",
      sp.simplify(sp.diff(sp.log(phi), r) - (-1 / r - m)) == 0)
# B5 effective dimension reading: D_eff = 1 - d ln F / d ln r, F = -G'
F = -sp.diff(phi, r)
Deff = sp.simplify(1 - r * sp.diff(sp.log(F), r))
check("B5 Yukawa: D_eff(r) = 3 + m^2 r^2/(1 + m r): 3 near, grows without bound far",
      sp.simplify(Deff - (3 + m ** 2 * r ** 2 / (1 + m * r))) == 0, str(Deff))
check("B5b massless: D_eff = 3 at every r", sp.simplify(Deff.subs(m, 0) - 3) == 0)
# B6 reach from base dials: in ruler units m = sqrt(Delta)/2, reach 2/sqrt(Delta); massless on d = k^2/4
pts = [(3.0, 1.0), (2.0, 0.5), (4.0, -1.0), (1.0, 0.1)]
okb6 = True
detail = []
h = 0.01
Lg = np.arange(-30, 30 + h / 2, h)
n_ = len(Lg)
for kk, dd in pts:
    Dl = kk * kk - 4 * dd
    # Green's function of -w'' + (Delta/4) w on a grid (finite differences), source at L = 0
    main = np.full(n_, 2 / h ** 2 + Dl / 4)
    off = np.full(n_ - 1, -1 / h ** 2)
    from scipy.linalg import solve_banded
    ab = np.zeros((3, n_)); ab[0, 1:] = off; ab[1] = main; ab[2, :-1] = off
    rhs = np.zeros(n_); rhs[n_ // 2] = 1 / h
    G = solve_banded((1, 1), ab, rhs)
    i1, i2 = n_ // 2 + int(4 / h), n_ // 2 + int(6 / h)
    rr_ = -(math.log(G[i2]) - math.log(G[i1])) / 2
    mex = 2 / h * math.asinh(h * math.sqrt(Dl) / 4)  # exact grid mass
    okb6 &= abs(rr_ - mex) < 1e-6 and abs(rr_ - math.sqrt(Dl) / 2) < 1e-4
    detail.append(f"({kk},{dd}): {rr_:.5f} vs {math.sqrt(Dl)/2:.5f}")
check("B6 Green's function of -w'' + (Delta/4) w falls at m = sqrt(Delta)/2 (reach 2/sqrt Delta), on a grid", okb6,
      "; ".join(detail))
# B7 geometric screening: hyperbolic space of curvature -kappa (kappa = s^2): massless potential
phiH = s * (sp.cosh(s * r) / sp.sinh(s * r) - 1) / (4 * sp.pi)
lapH = sp.diff(sp.sinh(s * r) ** 2 * sp.diff(phiH, r), r) / sp.sinh(s * r) ** 2
check("B7 H^3 (curvature -s^2): s(coth(sr) - 1)/(4 pi) is harmonic (r>0)", sp.simplify(lapH) == 0)
near = sp.limit(phiH * 4 * sp.pi * r, r, 0)
check("B7b near the source it is Coulomb 1/(4 pi r)", sp.simplify(near - 1) == 0)
ss = 0.7
fH = lambda rr: ss * 2 / math.expm1(2 * ss * rr) / (4 * math.pi)  # coth x - 1 = 2/(e^{2x} - 1)
rateH = -(math.log(fH(30.0)) - math.log(fH(29.0)))
check("B7c far away it falls like e^{-2 s r}: negative curvature screens with mass (D-1) sqrt(kappa)",
      abs(rateH - 2 * ss) < 1e-8, f"{rateH:.10f} vs {2*ss}")
# B7d reduced form: phi = u / sinh(sr) with u'' = s^2 u: curved space acts as mass s in reduced form
uH = sp.exp(-s * r)
phi2 = uH / sp.sinh(s * r)
lap2 = sp.diff(sp.sinh(s * r) ** 2 * sp.diff(phi2, r), r) / sp.sinh(s * r) ** 2
check("B7d e^{-sr}/sinh(sr) is also harmonic on H^3 (same function up to a factor 2 e^{...})",
      sp.simplify(lap2) == 0)

# ------------------------------------------------------------------ C
say("== C. A closed direction: massless zero mode plus massive rungs")
# 4D massless field on R^3 x circle of length Pc. Image sum vs rung (Kaluza-Klein) sum.
Pc = 2 * math.pi / math.sqrt(3)   # the circle-class ruler length for k=1, d=1
mu = 2 * math.pi / Pc             # = sqrt(-Delta)
errs = []
for rr in (0.05, 0.3, 1.0, 3.0):
    img = sum(1 / (4 * math.pi ** 2 * (rr * rr + (j * Pc) ** 2)) for j in range(-20000, 20001))
    kk_ = sum(math.exp(-abs(n) * mu * rr) for n in range(-400, 401)) / (4 * math.pi * rr * Pc)
    closed = 1 / math.tanh(math.pi * rr / Pc) / (4 * math.pi * rr * Pc)
    errs.append(max(abs(img - closed) / closed, abs(kk_ - closed) / closed))
check("C1 images on the circle = sum over rungs n of e^{-|n| sqrt(-Delta) r}/(4 pi r P) = coth(pi r/P)/(4 pi r P)",
      max(errs) < 1e-4, f"max rel err {max(errs):.1e}")
check("C1b rung masses are n sqrt(-Delta): mu = 2 pi / P", abs(mu - math.sqrt(3)) < 1e-12, f"mu = {mu:.12f}")
rr = 6.0
ratio = (1 / math.tanh(math.pi * rr / Pc)) - 1
check("C2 far away only rung 0 survives (1/(4 pi r P): 3D Coulomb, infinite reach); rung 1 is e^{-sqrt3 r}",
      abs(ratio - 2 * math.exp(-mu * rr)) / (2 * math.exp(-mu * rr)) < 1e-4, f"coth-1 = {ratio:.3e}")

# ------------------------------------------------------------------ D
say("== D. Confinement: flux that keeps one open direction")
# D1 dimension-dial pull 1/(S_D r^{D-1}) at D = 1: S_1 = 2, constant pull 1/2, linear potential
SD = lambda D: 2 * math.pi ** (D / 2) / math.gamma(D / 2)
check("D1 dimension dial at D=1: pull 1/(S_1 r^0) = 1/2 per unit charge, constant: potential grows linearly",
      abs(SD(1) - 2) < 1e-14)
# D2 a charge in a square tube of side a with walls that keep flux in (Neumann). Pull on the axis.
a = 1.0


def pull_modes(z, nmax=60):
    tot = 1 / (2 * a * a)
    for mi in range(0, nmax, 2):
        for ni in range(0, nmax, 2):
            if mi == 0 and ni == 0:
                continue
            em = 1 if mi == 0 else 2
            en = 1 if ni == 0 else 2
            kmn = math.pi * math.hypot(mi, ni) / a
            tot += em * en / (a * a) * math.exp(-kmn * z) / 2
    return tot


def pull_images(z, N=400):
    i = np.arange(-N, N + 1)
    I, J = np.meshgrid(i, i)
    rho2 = (a * I) ** 2 + (a * J) ** 2
    terms = z / (4 * math.pi * (z * z + rho2) ** 1.5)
    # tail beyond the box treated as a smooth sheet: int_{rho>R} 2 pi rho z/(4 pi (z^2+rho^2)^1.5) drho / a^2
    R = (N + 0.5) * a
    tail = z / (2 * a * a * math.sqrt(z * z + R * R))
    return terms.sum() + tail


errs = [abs(pull_modes(z) - pull_images(z)) / pull_modes(z) for z in (0.3, 0.5, 0.8)]
check("D2 tube pull: wall-mode sum and image sum agree where both converge", max(errs) < 1e-3,
      f"max rel diff {max(errs):.1e}")
zs = np.logspace(-2, 0.7, 60)
Fz = np.array([pull_images(z) if z < 0.4 else pull_modes(z) for z in zs])
Deff_tube = 1 - np.gradient(np.log(Fz), np.log(zs))
RES['tube_z'] = zs.tolist()
RES['tube_F'] = Fz.tolist()
check("D2b tube: D_eff runs from 3 near the charge to 1 far away",
      abs(Deff_tube[0] - 3) < 0.02 and abs(Deff_tube[-1] - 1) < 1e-3,
      f"D_eff near {Deff_tube[0]:.4f}, far {Deff_tube[-1]:.6f}")
check("D2c tube far pull = flux / (2 x area): constant, so the potential is linear (string tension)",
      abs(pull_modes(5.0) - 1 / (2 * a * a)) < 1e-12, f"{pull_modes(5.0):.15f}")
# D2d the transverse wall modes are themselves screened: the lowest falls like e^{-2 pi z/a}
corr = lambda z: pull_modes(z) - 1 / (2 * a * a) if z < 1 else sum(
    (1 if mi == 0 else 2) * (1 if ni == 0 else 2) / (a * a) * math.exp(-math.pi * math.hypot(mi, ni) / a * z) / 2
    for mi in range(0, 60, 2) for ni in range(0, 60, 2) if (mi, ni) != (0, 0))
rate_t = -(math.log(corr(7.0)) - math.log(corr(6.0)))
check("D2d tube: the sideways rungs are massive (rate 2 pi / a for a centred charge): only rung 0 reaches far",
      abs(rate_t - 2 * math.pi / a) < 1e-5, f"{rate_t:.10f}")
# D3 no constant K turns a point source into a tube: R' = cos(sqrt K r), 1, cosh(sqrt(-K) r)
#    never stays at 0; on the sphere R' = 0 only at the equator, then the area shrinks.
Rp_pos = sp.diff(Rpos, r)
eq = sp.solve(sp.Eq(Rp_pos.subs(K, 1), 0), r)
check("D3 no constant curvature keeps a point source's area fixed: R' vanishes only at isolated points (sphere's equator)",
      len(eq) >= 1 and sp.simplify(sp.diff(Rp_pos.subs(K, 1), r).subs(r, sp.pi / 2)) != 0,
      f"sphere K=1: R'=0 at r = {eq}")

# ------------------------------------------------------------------ E
say("== E. A tube made by screening sideways (London / Abrikosov)")
Phi = 1.0
for mm in (0.5, 1.0, 2.0):
    B = lambda rho: Phi * mm ** 2 / (2 * math.pi) * special.k0(mm * rho)
    flux = integrate.quad(lambda rho: B(rho) * 2 * math.pi * rho, 0, np.inf)[0]
    rms2 = integrate.quad(lambda rho: rho ** 2 * B(rho) * 2 * math.pi * rho, 0, np.inf)[0] / flux
    if mm == 1.0:
        RES['E_flux'] = flux
    check(f"E1 m={mm}: B = Phi m^2 K0(m rho)/(2 pi) carries flux Phi and has rms radius 2/m",
          abs(flux - Phi) < 1e-9 and abs(math.sqrt(rms2) - 2 / mm) < 1e-8,
          f"flux {flux:.12f}, rms {math.sqrt(rms2):.10f} vs {2/mm}")
rho = sp.symbols('rho', positive=True)
Bs = sp.besselk(0, m * rho)
check("E2 K0(m rho) solves (nabla^2 - m^2) B = 0 off the core (2D radial)",
      sp.simplify(sp.diff(rho * sp.diff(Bs, rho), rho) / rho - m ** 2 * Bs) == 0)
tail = [-(math.log(special.k0(x + 1) * math.sqrt(x + 1)) - math.log(special.k0(x) * math.sqrt(x))) for x in (30.0,)]
check("E3 sideways tail e^{-m rho}/sqrt(rho): the same mass m as a Yukawa field", abs(tail[0] - 1) < 2e-3,
      f"rate {tail[0]:.5f}")
# E4 London energy per length with core xi: closed form (Phi^2 m^2/4pi)(m xi)K0(m xi)K1(m xi)
okE = True
dets = []
for mm, xi in [(1.0, 0.01), (0.5, 0.02), (2.0, 0.005)]:
    B = lambda rr_: Phi * mm ** 2 / (2 * math.pi) * special.k0(mm * rr_)
    dB = lambda rr_: -Phi * mm ** 3 / (2 * math.pi) * special.k1(mm * rr_)
    E = 0.5 * integrate.quad(lambda rr_: (B(rr_) ** 2 + dB(rr_) ** 2 / mm ** 2) * 2 * math.pi * rr_, xi, np.inf,
                             limit=400)[0]
    Ec = Phi ** 2 * mm ** 2 / (4 * math.pi) * (mm * xi) * special.k0(mm * xi) * special.k1(mm * xi)
    okE &= abs(E - Ec) / Ec < 1e-7
    dets.append(f"{E:.8f}/{Ec:.8f}")
check("E4 tube tension sigma = (Phi^2 m^2 / 4 pi) (m xi) K0(m xi) K1(m xi) ~ (Phi^2 m^2/4 pi) ln(1/(m xi))", okE,
      "; ".join(dets))
# E5 doubling m (same core) multiplies sigma by about 4 (times a slowly falling log)
s1 = (1.0) ** 2 * (0.01) * special.k0(0.01) * special.k1(0.01)
s2 = (2.0) ** 2 * (0.02) * special.k0(0.02) * special.k1(0.02)
check("E5 one mass sets both: width 2/m and tension ~ m^2 (x log): doubling m gives ~4x tension",
      3.0 < s2 / s1 < 4.0, f"ratio {s2/s1:.4f}")

# ------------------------------------------------------------------ F
say("== F. The forces placed (reach = hbar c / mc^2)")
hbarc_MeV_fm = 197.3269804
masses = {  # MeV; sources in README
    'W': 80369.2, 'Z': 91188.0, 'pi+-': 139.57039, 'pi0': 134.9768,
}
reach = {kname: hbarc_MeV_fm / v for kname, v in masses.items()}
for kname, v in reach.items():
    say(f"   reach of {kname}: {v:.6g} fm")
RES['reach_fm'] = reach
hbarc_eV_m = 197.3269804e-9
for kname, bound in (('photon', 1e-18), ('graviton', 2.42e-23)):
    v = hbarc_eV_m / bound
    say(f"   {kname}: mass below {bound} eV  ->  reach above {v:.3g} m")
    RES[f'{kname}_min_reach_m'] = v
sqrt_sigma = 440.0  # MeV, typical lattice value
say(f"   quark string tension sigma ~ (440 MeV)^2 = {sqrt_sigma**2/hbarc_MeV_fm:.0f} MeV/fm "
    f"= {sqrt_sigma**2/hbarc_MeV_fm*1.602176634e-13/1e-15/1e3:.3g} kN")
check("F1 W reach is under 1/500 of the pion's", 500 < reach['pi+-'] / reach['W'] < 700,
      f"{reach['pi+-']/reach['W']:.1f}")

# ------------------------------------------------------------------ P
say("== P. Polishing (30 digits)")
mp.mp.dps = 30
mm, rr = mp.mpf(1), mp.mpf(2)
I = mp.quadosc(lambda qq: qq * mp.sin(qq * rr) / (qq ** 2 + mm ** 2), [0, mp.inf], omega=rr)
G = I / (2 * mp.pi ** 2 * rr)
ex = mp.e ** (-mm * rr) / (4 * mp.pi * rr)
check("P1 Fourier Yukawa at 30 digits", abs(G - ex) / ex < mp.mpf('1e-25'), mp.nstr(abs(G - ex) / ex, 3))
# P2 decay rate of line-class flow, exact: x - wall_lo = C e^{-sqrt(Delta) L}(1 + ...); log-derivative limit
Dl = mp.mpf(5)
sq = mp.sqrt(Dl)
wl, wh = (-3 - sq) / 2, (-3 + sq) / 2
# exact solution of x' = (x - wl)(x - wh): (x - wl)/(x - wh) = A e^{(wl - wh) L}
R0 = (mp.mpf(-1) - wh) / (mp.mpf(-1) - wl)   # (x - wh)/(x - wl) = R0 e^{(wh - wl) L}
xL = lambda LL: (wh - wl * R0 * mp.e ** ((wh - wl) * LL)) / (1 - R0 * mp.e ** ((wh - wl) * LL))
gap = lambda LL: (wh - wl) / (1 - R0 * mp.e ** ((wh - wl) * LL))   # = x(L) - wl, exactly
rate30 = -(mp.log(abs(gap(40))) - mp.log(abs(gap(39))))
check("P2 line-class wall rate = sqrt(Delta) to 30 digits (moved from fit)", abs(rate30 - sq) < mp.mpf('1e-15'),
      f"{mp.nstr(rate30, 20)}; fit was {RES['A3_rate']:.6f}, moved {float(abs(rate30 - RES['A3_rate'])):.1e}")
# P3 London tension closed form vs quadrature at 30 digits
mm, xi = mp.mpf(1), mp.mpf('0.01')
Bf = lambda t: mm ** 2 / (2 * mp.pi) * mp.besselk(0, mm * t)
dBf = lambda t: -mm ** 3 / (2 * mp.pi) * mp.besselk(1, mm * t)
E = mp.quad(lambda t: (Bf(t) ** 2 + dBf(t) ** 2 / mm ** 2) * mp.pi * t, [xi, 1, 10, mp.inf])
Ec = mm ** 2 / (4 * mp.pi) * (mm * xi) * mp.besselk(0, mm * xi) * mp.besselk(1, mm * xi)
check("P3 tube tension closed form at 30 digits", abs(E - Ec) / Ec < mp.mpf('1e-25'), mp.nstr(abs(E - Ec) / Ec, 3))
# P4 tube flux exactly Phi
fl = mp.quad(lambda t: mm ** 2 / (2 * mp.pi) * mp.besselk(0, mm * t) * 2 * mp.pi * t, [0, 1, mp.inf])
check("P4 tube flux = Phi at 30 digits", abs(fl - 1) < mp.mpf('1e-25'), mp.nstr(abs(fl - 1), 3))

say(f"\nTOTAL: {npass} pass, {nfail} fail")
RES['pass'], RES['fail'] = npass, nfail
open('checks_output.txt', 'w').write("\n".join(OUT) + "\n")
json.dump(RES, open('results.json', 'w'), indent=1, default=float)
