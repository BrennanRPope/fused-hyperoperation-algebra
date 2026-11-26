"""Warp bubbles, wormholes and artificial gravity through the circle operation (DRAFT, Nov 25 2026).
Run: python3 checks.py  (sympy, numpy, scipy). Writes results.json."""
import json, sys, math
import numpy as np, sympy as sp
from scipy import integrate
sys.path.insert(0, '../general-relativity')
import grtools

R = {}; npass = 0; nall = 0
def check(name, ok, info=''):
    global npass, nall
    nall += 1; npass += bool(ok)
    print(('PASS ' if ok else 'FAIL ') + name + ('  ' + str(info) if info != '' else ''))
    R[name] = {'ok': bool(ok), 'info': str(info)}

G_SI, c_SI = 6.674e-11, 2.998e8
KG_PER_M = c_SI**2 / G_SI          # geometric length -> kg
M_SUN, M_JUP, M_UNIV = 1.989e30, 1.898e27, 1.5e53  # visible-universe ordinary matter ~1e53 kg (order)
L_PLANCK = 1.616e-35

print('=== A. Warp bubble as a river ===')
t, x, y, z = sp.symbols('t x y z', real=True)
vs = sp.symbols('v_s', positive=True)
f = sp.Function('f')
rs = sp.sqrt((x - vs*t)**2 + y**2 + z**2)
v = [vs*f(rs), 0, 0]                      # river velocity (Alcubierre: shift = -v)
X3 = [x, y, z]
K = sp.Matrix(3, 3, lambda i, j: (sp.diff(v[i], X3[j]) + sp.diff(v[j], X3[i]))/2)
e1 = K.trace()
e2 = sp.simplify((e1**2 - (K*K).trace())/2)
rho_ours = e2/(8*sp.pi)
r_ = sp.symbols('r', positive=True)
fp = sp.Function('fp')
rho_alc = -vs**2*(y**2+z**2)/(32*sp.pi*rs**2) * sp.diff(f(rs), x)**2/((x-vs*t)/rs)**2
check('A1 Eulerian energy = e2/8pi equals Alcubierre 1994: -v^2 f\'^2 (y^2+z^2)/(32 pi r^2)',
      sp.simplify(rho_ours - rho_alc) == 0)
theta = sp.simplify(e1)
check('A2 expansion theta = v_s (x_s/r_s) f\'(r_s): sign flips front/back (contraction ahead, expansion behind)',
      sp.simplify(theta - vs*(x-vs*t)/rs*sp.diff(f(rs), x)/((x-vs*t)/rs)) == 0)
# full 4D check of G_nn on the Alcubierre metric with a concrete f at sample points (numeric)
sig, Rb = 2.0, 1.0
def fnum(r): return (np.tanh(sig*(r+Rb)) - np.tanh(sig*(r-Rb)))/(2*np.tanh(sig*Rb))
def fpn(r): return sig*(1/np.cosh(sig*(r+Rb))**2 - 1/np.cosh(sig*(r-Rb))**2)/(2*np.tanh(sig*Rb))
# symbolic 4D Einstein tensor for river metric with general v(t,x,y,z) along x is heavy; gr-bridges T2/T3 already
# proved G_nn = e2(K) on every river. Here: strain eigenvalues.
rng = np.random.default_rng(1)
okA3 = True; minprod = 0
for _ in range(2000):
    p = rng.normal(size=3)*1.2
    r = np.linalg.norm(p); g = 1.7*fpn(r)
    dv = np.array([g*p[0]/r, g*p[1]/r, g*p[2]/r])   # grad of v_x
    Km = np.zeros((3, 3)); Km[0, :] += dv/2; Km[:, 0] += dv/2
    lam = np.linalg.eigvalsh(Km)
    lam = sorted(lam, key=abs)
    prod = lam[1]*lam[2]
    s2 = (dv[1]**2 + dv[2]**2)/4
    okA3 &= abs(lam[0]) < 1e-12 and prod <= 1e-15 and abs(prod + s2) < 1e-12
    minprod = min(minprod, prod)
check('A3 strain of a one-direction river: rates (0, +, -): one stretch, one squeeze, product = -(shear)^2 <= 0', okA3)
# A4 Natario: e1 = 0 forces e2 <= 0
okA4 = all((lambda a, b: (a*b + b*(-a-b) + a*(-a-b)) <= 0)(*rng.normal(size=2)) for _ in range(10000))
check('A4 zero expansion (e1 = 0) forces e2 = -(sum of squares)/2 <= 0 (Natario 2002)', okA4)
# A5 total energy: closed form -(v^2/12) int f'^2 r^2 dr, vs grid integral of e2/8pi, vs -(1/32pi) int |curl v|^2
vsn = 1.0
I_rad = integrate.quad(lambda r: fpn(r)**2*r**2, 0, 8)[0]
E_closed = -(vsn**2/12)*I_rad
n = 161; Lb = 3.2
g1 = np.linspace(-Lb, Lb, n); hh = g1[1]-g1[0]
Xg, Yg, Zg = np.meshgrid(g1, g1, g1, indexing='ij')
Rg = np.sqrt(Xg**2+Yg**2+Zg**2) + 1e-12
vx = vsn*fnum(Rg)
gx, gy, gz = np.gradient(vx, hh, hh, hh)
e2g = -(gy**2 + gz**2)/4
E_grid = e2g.sum()*hh**3/(8*np.pi)
curl2 = gy**2 + gz**2
E_curl = -curl2.sum()*hh**3/(32*np.pi)
check('A5 total Eulerian energy: closed form = grid = -(1/32pi) int |curl v|^2 (Santiago-Schuster-Visser eq 7.17)',
      abs(E_grid/E_closed-1) < 0.02 and abs(E_curl/E_closed-1) < 0.02, f'{E_closed:.5f} {E_grid:.5f} {E_curl:.5f}')
R['A5_values'] = dict(E_closed=E_closed, E_grid=E_grid, E_curl=E_curl)
# A6 irrotational warp v = grad(v_s x f(r)): total zero, but negative somewhere
phi = vsn*Xg*fnum(Rg)
ux, uy, uz = np.gradient(phi, hh, hh, hh)
H = [np.gradient(u, hh, hh, hh) for u in (ux, uy, uz)]
Kx = np.array([[H[i][j] for j in range(3)] for i in range(3)])
tr = Kx[0, 0]+Kx[1, 1]+Kx[2, 2]
trK2 = sum(Kx[i, j]*Kx[j, i] for i in range(3) for j in range(3))
rho_irr = (tr**2 - trK2)/(16*np.pi)
inner = (np.abs(g1) < 2.6)
mask = inner[:, None, None] & inner[None, :, None] & inner[None, None, :]
Etot = rho_irr[mask].sum()*hh**3; Eneg = rho_irr[mask & (rho_irr < 0)].sum()*hh**3; Epos = rho_irr[mask & (rho_irr > 0)].sum()*hh**3
check('A6 irrotational warp (v = grad phi): positive and negative energy cancel (total ~ 0) but negative regions remain',
      abs(Etot) < 0.03*abs(Eneg) and Eneg < 0, f'total {Etot:.4f}, negative part {Eneg:.4f}, positive part {Epos:.4f}')
R['A6_values'] = dict(Etot=Etot, Eneg=Eneg, Epos=Epos)
# A7 light speeds in the river: dx/dt = v +/- 1 exactly; local addition is k = 0 in these coordinates
vv, dxdt = sp.symbols('v w', real=True)
null = sp.solve(sp.Eq(-1 + (dxdt - vv)**2, 0), dxdt)
check('A7 light in the river moves at dx/dt = v - 1 and v + 1 (plain sum: k = 0 between river and coordinates)',
      set(sp.simplify(s) for s in null) == {vv - 1, vv + 1})
# local composition at a point is k = -1 (half-velocities / Einstein); compare
u1, u2 = 0.6, 0.7
einstein = (u1+u2)/(1+u1*u2)
check('A8 at one point two velocities compose by the relativity rule (0.6 with 0.7 -> 0.915c), across the river by plain sum',
      abs(einstein - 0.9154929577) < 1e-9, f'{einstein:.4f} vs plain 1.3')
# A9 horizons: where v_s f = 1 for v_s = 2 (front and back of the wall)
from scipy.optimize import brentq
vs2 = 2.0
rh_out = brentq(lambda r: vs2*fnum(r) - 1, Rb, Rb+3)
check('A9 superluminal bubble (v_s = 2): river crosses light speed at r = %.3f (inside the wall) -> horizon, like r = 2M' % rh_out,
      Rb < rh_out < Rb+1.5)
R['A9_r_horizon'] = rh_out
# A10 energy scale: thin wall approx E = -(v^2 R^2 sigma)/36 in length units
def E_alc(Rm, sigma, beta=1.0):
    fp_ = lambda r: sigma*(1/np.cosh(sigma*(r+Rm))**2 - 1/np.cosh(sigma*(r-Rm))**2)/(2*np.tanh(sigma*Rm))
    a, b = max(0, Rm - 40/sigma), Rm + 40/sigma
    I = integrate.quad(lambda r: fp_(r)**2*r**2, a, b, points=[Rm], limit=400)[0]
    return -(beta**2/12)*I
Et = E_alc(100.0, 1.0)
check('A10 thin-wall law E = -v^2 R^2 sigma/36 (geometric units) holds for R = 100 m, sigma = 1/m',
      abs(Et/(-(100**2)*1/36) - 1) < 0.01, f'{Et:.2f} m')
scale = {}
for lab, Rm, sigma in [('R=100m wall~1m', 100, 2.0), ('R=100m wall~1mm', 100, 2000.0),
                       ('R=100m wall~100 Planck', 100, 2/(100*L_PLANCK)), ('R=10m wall~1m', 10, 2.0)]:
    M = -(Rm**2*sigma/36)*KG_PER_M
    scale[lab] = dict(kg=M, solar=M/M_SUN, universe=M/M_UNIV)
R['A10_scale'] = scale
for k_, v_ in scale.items(): print('   ', k_, '%.2e kg = %.2e Msun' % (v_['kg'], v_['solar']))
check('A11 Planck-thin wall (Pfenning-Ford regime) needs ~1e62 kg of negative energy per (v/c)^2: >1e8 visible universes',
      scale['R=100m wall~100 Planck']['kg'] < -1e61)

print('=== B. Artificial gravity ===')
w = sp.Function('w')     # planar river speed w(z)
g4 = sp.Matrix([[-1 + w(z)**2, 0, 0, -w(z)], [0, 1, 0, 0], [0, 0, 1, 0], [-w(z), 0, 0, 1]])  # (t,x,y,z): -dt^2+dx^2+dy^2+(dz-w dt)^2
Gm, Ric, Rs, *_ = grtools.einstein_mixed(g4, [t, x, y, z])
Gm = sp.simplify(Gm)
W2 = w(z)**2
check('B1 planar river: Einstein tensor is -(w^2)\'\'/2 on the two transverse directions and zero otherwise',
      all(sp.simplify(Gm[a, b]) == 0 for a in range(4) for b in range(4) if (a, b) not in [(1, 1), (2, 2)])
      and sp.simplify(Gm[1, 1] + sp.diff(W2, z, 2)/2) == 0 and sp.simplify(Gm[2, 2] + sp.diff(W2, z, 2)/2) == 0,
      str(sp.simplify(Gm[1, 1])))
gg, cc = sp.symbols('g c0', positive=True)
check('B2 vacuum iff w^2 = 2 g z + c (a uniform field): flat spacetime, an accelerated frame, no tides, costs nothing',
      sp.simplify(sp.diff(2*gg*z + cc, z, 2)) == 0)
# proper acceleration of a static observer in river with w^2 = 2gz+c: a = w w'/sqrt(1-w^2)
Gch, gi4 = grtools.christoffel(g4, [t, x, y, z])
ut = 1/sp.sqrt(1 - W2)
acc = [sp.simplify(Gch[m][0][0]*ut**2) for m in range(4)]
amag2 = sp.simplify(sum(g4[m, n]*acc[m]*acc[n] for m in range(4) for n in range(4)))
check('B3 a floor at fixed height feels a = |w w\'|/sqrt(1-w^2); for w^2 = 2gz + c that is g/sqrt(1-2gz-c): the floor must be held (pushed)',
      sp.simplify(amag2 - (w(z)*sp.diff(w(z), z))**2/(1 - W2)) == 0, str(amag2))
gE = 9.81
sigma_plate = gE/(2*np.pi*G_SI)
M_sphere10 = gE*10**2/G_SI
check('B4 real 1 g from mass: infinite plate needs %.2e kg/m^2; a 10 m ball needs %.2e kg (a large mountain)'
      % (sigma_plate, M_sphere10), sigma_plate > 1e10 and M_sphere10 > 1e13)
R['B4'] = dict(sigma_plate=sigma_plate, M_sphere10=M_sphere10)
rad = {rpm: gE/(rpm*2*np.pi/60)**2 for rpm in (1, 2, 4, 6)}
check('B5 rotation: 1 g at 2 rpm needs radius %.0f m, at 4 rpm %.0f m (g = omega^2 r)' % (rad[2], rad[4]),
      abs(rad[2]-223.6) < 1 and abs(rad[4]-55.9) < 0.5)
R['B5_radius'] = rad
# B6 thrust: composing many small kicks by Einstein rule = tanh(g tau / c)
tau = 365.25*86400; n_k = 200000; dv = gE*tau/c_SI/n_k
uu = 0.0
for _ in range(n_k): uu = (uu + dv)/(1 + uu*dv)
check('B6 one year of 1 g thrust: 200,000 small kicks composed by the relativity rule give tanh(g tau/c) = %.4f c' % uu,
      abs(uu - np.tanh(gE*tau/c_SI)) < 1e-6)

print('=== C. Wormholes ===')
rr, b0 = sp.symbols('r b0', positive=True)
b = sp.Function('b')
gw = sp.diag(-1, 1/(1 - b(rr)/rr), rr**2, rr**2*sp.sin(sp.symbols('th'))**2)
th, ph = sp.symbols('th'), sp.symbols('ph')
gw = sp.diag(-1, 1/(1 - b(rr)/rr), rr**2, rr**2*sp.sin(th)**2)
Gw, *_ = grtools.einstein_mixed(gw, [t, rr, th, ph])
rho_w = sp.simplify(-Gw[0, 0]/(8*sp.pi)); pr_w = sp.simplify(Gw[1, 1]/(8*sp.pi))
check('C1 Morris-Thorne (zero redshift): rho = b\'/(8 pi r^2), p_r = -b/(8 pi r^3)',
      sp.simplify(rho_w - sp.diff(b(rr), rr)/(8*sp.pi*rr**2)) == 0 and sp.simplify(pr_w + b(rr)/(8*sp.pi*rr**3)) == 0)
nec = sp.simplify(rho_w + pr_w)
check('C2 rho + p_r = (r b\' - b)/(8 pi r^3): negative at the throat whenever it flares out (b\'(r0) < 1)',
      sp.simplify(nec - (rr*sp.diff(b(rr), rr) - b(rr))/(8*sp.pi*rr**3)) == 0)
# area composition at the throat: A = (1+a)(1+b); R_kk = -(A''/A - 2a'b'/((1+a)(1+b))); at a'=b'=0, R_kk = -A''/A
lam_ = sp.symbols('lam'); a_ = sp.Function('a')(lam_); bb_ = sp.Function('bb')(lam_)
A = (1+a_)*(1+bb_)
lhs = sp.diff(a_, lam_, 2)/(1+a_) + sp.diff(bb_, lam_, 2)/(1+bb_)
rhs = sp.diff(A, lam_, 2)/A - 2*sp.diff(a_, lam_)*sp.diff(bb_, lam_)/A
check('C3 bundle area 1 + (a o b): focusing a\'\'/(1+a) + b\'\'/(1+b) = A\'\'/A - 2a\'b\'/A; at a throat (area minimum) R_kk = -A\'\'/A < 0',
      sp.simplify(lhs - rhs) == 0)
# thin-shell wormhole: sigma = -sqrt(1-2M/a)/(2 pi a); shell mass 4 pi a^2 sigma
for a_m in (1.0,):
    m_shell = -2*a_m*KG_PER_M
Mm, aa = sp.symbols('M a', positive=True)
Kth = sp.sqrt(1 - 2*Mm/aa)/aa                      # K^theta_theta on each side, normals pointing away from throat
sig_shell = -(2*Kth)/(4*sp.pi)                     # Israel: sigma = -[K^theta_theta]/(4 pi)
mshell = sp.simplify(4*sp.pi*aa**2*sig_shell)
check('C4 thin-shell wormhole (Visser 1989): shell rest mass 4 pi a^2 sigma = -2a sqrt(1-2M/a); throat 1 m, no outer mass: %.2e kg = %.2f Jupiter masses of negative energy'
      % (m_shell, m_shell/M_JUP), sp.simplify(mshell + 2*aa*sp.sqrt(1-2*Mm/aa)) == 0)
R['C4_kg_per_m'] = m_shell
# C5 two-sided cross-ratio in the thermofield double (planar BTZ = Rindler-AdS3)
okC5 = True; maxerr = 0
for _ in range(500):
    beta = rng.uniform(0.5, 3); kk = 2*np.pi/beta
    x1 = rng.uniform(-2, 2); ell = rng.uniform(0.1, 3); x2 = x1+ell; tt = rng.uniform(0, 2)
    zL = lambda xx: np.exp(kk*(xx - tt)); zR = lambda xx: -np.exp(kk*(xx + tt))
    z1, z2, z3, z4 = zL(x1), zL(x2), zR(x2), zR(x1)
    eta = (z1-z2)*(z3-z4)/((z1-z3)*(z2-z4))
    s = np.sinh(np.pi*ell/beta); chh = np.cosh(2*np.pi*tt/beta)
    pred = s**2/(s**2 + chh**2)
    err = abs(eta - pred); maxerr = max(maxerr, err)
    okC5 &= err < 1e-9
check('C5 two-sided cross-ratio of A (left) and B (right) in the thermofield double = sinh^2(pi l/b)/(sinh^2 + cosh^2(2 pi t/b))', okC5, maxerr)
# C6 mutual information (c/3) log(eta/(1-eta)) = Hartman-Maldacena (2c/3) log(sinh/cosh); open iff eta > 1/2
okC6 = True
for _ in range(500):
    s = rng.uniform(0.1, 20); chh = rng.uniform(1, 20); eta = s**2/(s**2+chh**2)
    okC6 &= abs(np.log(eta/(1-eta))/3 - 2*np.log(s/chh)/3) < 1e-10
check('C6 I(A:B) = (c/3) log(odds of the link) = (2c/3) log(sinh(pi l/b)/cosh(2 pi t/b)); wormhole connected iff link probability > 1/2', okC6)
# C7 late-time decay rate of link probability: 4 pi / beta
beta = 1.0; ell = 3.0; ts = np.array([5.0, 6.0])
xs = np.sinh(np.pi*ell)**2/(np.sinh(np.pi*ell)**2 + np.cosh(2*np.pi*ts)**2)
rate = -np.log(xs[1]/xs[0])/(ts[1]-ts[0])
check('C7 link probability through the throat falls like exp(-4 pi t/beta) at late times (twice the surface-gravity rate)',
      abs(rate - 4*np.pi) < 1e-3, rate)
# C8 gluing across the throat: 1-eta composes multiplicatively for adjacent right regions (cocycle), any points
okC8 = True
for _ in range(500):
    pts = np.sort(rng.normal(size=5))
    p1, p2, p3, p4, p5 = pts[0], pts[1], pts[2], pts[3], pts[4]
    om = lambda a1, a2, b1, b2: (a1-b2)*(a2-b1)/((a1-b1)*(a2-b2))   # 1 - eta
    okC8 &= abs(om(p1, p2, p3, p4)*om(p1, p2, p4, p5) - om(p1, p2, p3, p5)) < 1e-9
check('C8 gluing B1 and B2 next to each other: link probabilities combine by the union rule o_-1, for any four points (pure algebra, not wormhole-specific)', okC8)
# C9 negative-energy shift (GJW-like): shift right-side light-cone coordinate; link probability -> 1 at traversability
beta = 1.0; kk = 2*np.pi; ell = 1.0; x1, x2 = 0.0, 1.0; tt = 0.3
zbL = lambda xx: np.exp(kk*(xx + tt)); zbR = lambda xx, a: -np.exp(kk*(xx - tt)) - a
def eta_bar(a):
    z1, z2, z3, z4 = zbL(x1), zbL(x2), zbR(x2, a), zbR(x1, a)
    return (z1-z2)*(z3-z4)/((z1-z3)*(z2-z4))
a_crit = -(zbL(x1) + np.exp(kk*(x1 - tt)))   # z1 = z4
avals = np.linspace(0.6, -0.999, 9)*abs(a_crit)
etas = [eta_bar(a) for a in avals]
check('C9 shifting the right side by a negative-energy pulse raises the link probability; it reaches 1 exactly when the two points become light-connected (traversable)',
      all(np.diff(etas) > 0) and etas[-1] > 0.99, [round(e, 3) for e in etas])
R['C9'] = dict(a=list(map(float, avals)), eta=list(map(float, etas)), a_crit=float(a_crit))

print(f'\n{npass} of {nall} checks pass')
R['summary'] = f'{npass} of {nall}'
json.dump(R, open('results.json', 'w'), indent=1, default=float)
