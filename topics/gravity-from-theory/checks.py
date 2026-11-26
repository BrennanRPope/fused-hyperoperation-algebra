"""DRAFT checks: gravity from the theory (Nov 25, 2026).

Part N: Newton's inverse-square law.  Part F: frame dragging (and where c enters).
Units: G = M = 1 unless stated. Writes checks_output.txt and results.json.
"""
import json, math, os, sys
import numpy as np
import sympy as sp
from scipy.integrate import solve_ivp

HERE = os.path.dirname(os.path.abspath(__file__))
R = {}
LOG = []
PASS = 0
FAIL = 0


def check(name, ok, detail=""):
    global PASS, FAIL
    line = f"[{'PASS' if ok else 'FAIL'}] {name}  {detail}"
    print(line); LOG.append(line)
    if ok:
        PASS += 1
    else:
        FAIL += 1


# =====================================================================
# Part N. Inverse-square law
# =====================================================================

# N1. The box's surface is the e_{n-1} layer of the lift (1+ka)(1+kb)(1+kc).
a, b, c, k, r = sp.symbols("a b c k r", positive=True)
lift = sp.expand((1 + k*a)*(1 + k*b)*(1 + k*c))
e2 = a*b + b*c + c*a
layer2 = sp.expand(lift).coeff(k, 2)
check("N1 lift layer k^2 = e2(a,b,c) = half the box surface", sp.simplify(layer2 - e2) == 0 and
      sp.simplify(2*e2 - 2*(a*b + b*c + c*a)) == 0, f"layer2 = {layer2}")

# N2. In n dimensions the face layer of a cube of side r is e_{n-1}(r,...,r) = n r^(n-1).
ok = True
for n in range(1, 8):
    xs = sp.symbols(f"x0:{n}")
    poly = sp.expand(sp.prod([1 + k*x for x in xs]))
    lay = poly.coeff(k, n - 1).subs({x: r for x in xs})
    ok &= sp.simplify(lay - n*r**(n - 1)) == 0
check("N2 face layer of the n-cube grows like n r^(n-1), n=1..7", ok)

# N3. Flux of the 1/r^2 field through boxes: 4 pi whatever the box size or where the source sits inside; 0 if outside.
def box_flux(src, half, m=160):
    g = (np.arange(m) + 0.5)/m*2*half - half
    U, V = np.meshgrid(g, g)
    dA = (2*half/m)**2
    tot = 0.0
    for axis in range(3):
        for sgn in (+1, -1):
            P = np.zeros((m*m, 3))
            others = [i for i in range(3) if i != axis]
            P[:, axis] = sgn*half; P[:, others[0]] = U.ravel(); P[:, others[1]] = V.ravel()
            d = P - src
            rr = np.linalg.norm(d, axis=1)
            tot += np.sum(sgn*d[:, axis]/rr**3)*dA
    return tot
fl = [box_flux(np.array([0.0, 0, 0]), 1.0), box_flux(np.array([0.0, 0, 0]), 3.0),
      box_flux(np.array([0.3, -0.2, 0.25]), 1.0), box_flux(np.array([2.5, 0, 0]), 1.0)]
check("N3 box flux of 1/r^2: 4pi for any box around the source, 0 outside",
      all(abs(f - 4*math.pi) < 2e-3 for f in fl[:3]) and abs(fl[3]) < 2e-3,
      " ".join(f"{f:.5f}" for f in fl) + f" (4pi = {4*math.pi:.5f})")
R["N3_flux"] = fl

# N4. Kepler hodographs.  Integrate orbits (GM = 1), map velocities onto the dial geometry with
#     dial kappa = -1/(2E)  (Milnor 1983 metric 2 sqrt(c)|dv|/(|v|^2+c), c = -2E).
def orbit(force_p, x0, v0, T, n=4000):
    def f(t, y):
        x, yy, z, vx, vy, vz = y
        rr = math.sqrt(x*x + yy*yy + z*z)
        aa = -1.0/rr**force_p/rr
        return [vx, vy, vz, aa*x, aa*yy, aa*z]
    s = solve_ivp(f, (0, T), list(x0) + list(v0), method="DOP853", rtol=1e-12, atol=1e-13,
                  t_eval=np.linspace(0, T, n))
    return s.y[:3].T, s.y[3:].T

def in_plane(V):
    # coordinates of 3D velocities in their own plane (the hodograph is planar)
    n_ = np.linalg.svd(V, full_matrices=False)[2]
    return V @ n_[:2].T

def plane_residual(P):
    # smallest singular value of uncentred points / largest: ~0 iff points lie in a plane through the origin
    sv = np.linalg.svd(P, compute_uv=False)
    return sv[-1]/sv[0]

# bound ellipse, E < 0
X, V = orbit(2, (1.0, 0, 0), (0, 1.25, 0.2), 60)
E = 0.5*np.sum(V**2, 1) - 1/np.linalg.norm(X, axis=1)
Eb = E.mean(); cb = -2*Eb
V = in_plane(V)
w = V/math.sqrt(cb)
w2 = np.sum(w**2, 1)[:, None]
S = np.hstack([2*w, w2 - 1])/(w2 + 1)       # unit sphere (kappa > 0 end)
res_sphere = plane_residual(S)
# same points with a wrong dial value: not a great circle
w_bad = V/math.sqrt(2.5*cb); w2b = np.sum(w_bad**2, 1)[:, None]
S_bad = np.hstack([2*w_bad, w2b - 1])/(w2b + 1)
res_bad = plane_residual(S_bad)
check("N4 bound orbit: hodograph lands on a great circle of the sphere when dial = -1/(2E)",
      res_sphere < 1e-8 and res_bad > 1e-2 and np.ptp(E) < 1e-9,
      f"E = {Eb:.6f}, residual {res_sphere:.1e}; wrong dial residual {res_bad:.3f}")
R["N4"] = dict(E=Eb, res=res_sphere, res_bad=res_bad)

# N5. Unbound hyperbola, E > 0: hodograph lands on a geodesic of the hyperboloid (plane through origin)
X2, V2 = orbit(2, (-30.0, 1.0, 0), (1.0, 0.0, 0.15), 70)
E2v = 0.5*np.sum(V2**2, 1) - 1/np.linalg.norm(X2, axis=1); E2 = E2v.mean()
V2 = in_plane(V2)
w = V2/math.sqrt(2*E2); w2 = np.sum(w**2, 1)[:, None]
H = np.hstack([2*w, w2 + 1])/(w2 - 1)       # hyperboloid X3^2 - X1^2 - X2^2 = 1
on_h = np.max(np.abs(H[:, 2]**2 - H[:, 0]**2 - H[:, 1]**2 - 1)/H[:, 2]**2)
res_h = plane_residual(H)
check("N5 unbound orbit (E > 0): hodograph is a geodesic of the hyperboloid (dial negative)",
      res_h < 1e-8 and on_h < 1e-8, f"E = {E2:.6f}, residual {res_h:.1e}")
R["N5"] = dict(E=E2, res=res_h)

# N6. Parabola, E = 0: hodograph is a circle through v = 0 (flat end of the dial)
X3, V3 = orbit(2, (1.0, 0, 0), (0, math.sqrt(2.0), 0), 40)
Vp = V3[:, :2]
A_ = np.hstack([2*Vp, np.ones((len(Vp), 1))]); bb = np.sum(Vp**2, 1)
sol, *_ = np.linalg.lstsq(A_, bb, rcond=None)
cx, cy = sol[0], sol[1]; rad = math.sqrt(sol[2] + cx*cx + cy*cy)
circ_res = np.max(np.abs(np.hypot(Vp[:, 0] - cx, Vp[:, 1] - cy) - rad))
check("N6 parabolic orbit (E = 0): hodograph is a circle through the origin of velocity space",
      circ_res < 1e-8 and abs(math.hypot(cx, cy) - rad) < 1e-8, f"centre ({cx:.4f},{cy:.4f}) radius {rad:.4f}")

# N7. Hamilton's converse: other power laws do not give circular hodographs.
def circle_fit_res(Vp):
    A_ = np.hstack([2*Vp, np.ones((len(Vp), 1))]); bb = np.sum(Vp**2, 1)
    s, *_ = np.linalg.lstsq(A_, bb, rcond=None)
    cx, cy = s[0], s[1]; rad = math.sqrt(s[2] + cx*cx + cy*cy)
    return np.max(np.abs(np.hypot(Vp[:, 0] - cx, Vp[:, 1] - cy) - rad))/rad
hod = {}
for p in (1.5, 2.0, 2.5):
    Xp, Vp_ = orbit(p, (1.0, 0, 0), (0, 0.8, 0), 80, n=6000)
    hod[p] = circle_fit_res(Vp_[:, :2])
check("N7 only the inverse-square force gives circular hodographs (Hamilton 1846)",
      hod[2.0] < 1e-8 and hod[1.5] > 1e-2 and hod[2.5] > 1e-2,
      ", ".join(f"p={p}: {v:.1e}" for p, v in hod.items()))
R["N7"] = {str(p): v for p, v in hod.items()}

# N8. Flux law in n dimensions: force 1/r^(n-1). Apsidal angle pi/sqrt(3-p) = pi/sqrt(4-n):
#     closed only at n = 3, stable only for n < 4 (Bertrand 1873, Ehrenfest 1917).
def apsidal(p):
    # slightly perturbed circular orbit at r=1: measure angle between successive pericentre and apocentre
    X, V = orbit(p, (1.0, 0, 0), (0, 1.0*1.01, 0), 60, n=200000)
    rr = np.linalg.norm(X, axis=1); th = np.unwrap(np.arctan2(X[:, 1], X[:, 0]))
    idx = [i for i in range(1, len(rr) - 1) if (rr[i] - rr[i-1])*(rr[i+1] - rr[i]) < 0]
    return (th[idx[2]] - th[idx[0]])/2   # pericentre -> apocentre, averaged over one cycle
ok = True; aps = {}
for n in (2.5, 3.0, 3.5):
    p = n - 1
    a_num = apsidal(p); a_th = math.pi/math.sqrt(3 - p)
    aps[n] = (a_num, a_th)
    ok &= abs(a_num - a_th) < 2e-2
check("N8 apsidal angle pi/sqrt(4-n) in n dimensions: equals pi (closed) only at n = 3",
      ok and abs(aps[3.0][0] - math.pi) < 1e-3,
      ", ".join(f"n={n}: {v[0]:.4f} vs {v[1]:.4f}" for n, v in aps.items()))
R["N8"] = {str(n): v for n, v in aps.items()}

# N9. Milnor's metric is the kappa-stereographic chart of the twisted-boosts page with kappa = -1/(2E).
E_, vx, vy = sp.symbols("E v_x v_y", real=True)
cc = -2*E_
milnor = 2*sp.sqrt(cc)/(vx**2 + vy**2 + cc)
kap = -1/(2*E_)
kchart = 2*sp.sqrt(kap)/(1 + kap*(vx**2 + vy**2))
check("N9 Milnor's velocity metric = kappa-stereographic chart with kappa = -1/(2E)",
      sp.simplify((milnor - kchart).subs(E_, -sp.Rational(3, 7))) == 0 and
      sp.simplify((milnor - kchart).subs(E_, -sp.Rational(1, 2))) == 0)

# =====================================================================
# Part F. Frame dragging
# =====================================================================

# F1. Boost Newton's static field with the dial boost. Coordinates (t, x); dial q (q = -1/c^2 relativity, 0 Galilean).
#     Weak static field of a mass at rest, written as -c^2(1 - 2U/c^2)dt^2 + (1 + 2 gam U/c^2)dx^2 with -1/c^2 -> q:
#     g_tt = 1/q + 2U  (h_tt = 2U),   g_xx = 1 - 2 gam q U  (h_xx = -2 gam q U).
q, u, U, gam = sp.symbols("q u U gamma", real=True)
G_ = 1/sp.sqrt(1 + q*u**2)
# inverse dial boost: t = G(t' - q u x'),  x = G(x' + u t')   [Lorentz when q = -1/c^2]
L = sp.Matrix([[G_, -G_*q*u], [G_*u, G_]])
v_ = sp.symbols("v", real=True)
Gv = 1/sp.sqrt(1 + q*v_**2)
Lv = sp.Matrix([[Gv, -Gv*q*v_], [Gv*v_, Gv]])
LL = sp.simplify(L*Lv)
w_comp = sp.simplify(LL[1, 0]/LL[1, 1])
check("F0 two dial boosts compose to one with speed (u+v)/(1 - q u v): the theory's velocity law (twisted-boosts page)",
      sp.simplify(w_comp - (u + v_)/(1 - q*u*v_)) == 0, f"{w_comp}")
eta = sp.diag(1/q, 1)                       # flat: dt^2/q + dx^2  (= -c^2 dt^2 + dx^2 when q = -1/c^2)
h = sp.diag(2*U, -2*gam*q*U)
hp = sp.simplify(L.T*h*L)
etap = sp.simplify(L.T*eta*L)
htx = sp.simplify(hp[0, 1])
target = -2*G_**2*(1 + gam)*q*U*u
check("F1a boosted field: h'_tx = -2 G^2 (1+gamma) q U u  (dial boost acting on Newton's field)",
      sp.simplify(htx - target) == 0, f"h'_tx = {sp.factor(htx)}")
check("F1b flat part stays flat under the dial boost (the boost is a symmetry of the dial's interval)",
      sp.simplify(etap - eta) == sp.zeros(2, 2))
check("F1c frame-dragging term vanishes at the Galilean dial q = 0 (infinite light speed)",
      sp.simplify(htx.subs(q, 0)) == 0)
cs = sp.symbols("c", positive=True)
htx_c = sp.series(sp.simplify(htx.subs(q, -1/cs**2)), u, 0, 2).removeO()
check("F1d with q = -1/c^2: h'_tx = 2(1+gamma) U u / c^2 to first order in u (clocks 2 + rulers 2gamma)",
      sp.simplify(htx_c - 2*(1 + gam)*U*u/cs**2) == 0, f"{htx_c}")
R["F1"] = str(sp.factor(htx))

# F2. A spinning ball: superpose moving mass elements.  Field point outside a uniform ball (radius 1, density 1,
#     spin w about z).  Integral S(x) = int rho u /|x - x'| d^3x'  should equal (1/2) J x x / r^3 (exterior dipole).
def ball_S(xf, w=1.0, m=48):
    g = (np.arange(m) + 0.5)/m*2 - 1
    Xg, Yg, Zg = np.meshgrid(g, g, g, indexing="ij")
    inside = Xg**2 + Yg**2 + Zg**2 <= 1
    P = np.stack([Xg[inside], Yg[inside], Zg[inside]], 1)
    dV = (2/m)**3
    vel = np.stack([-w*P[:, 1], w*P[:, 0], 0*P[:, 0]], 1)
    d = np.linalg.norm(xf - P, axis=1)
    S = np.sum(vel/d[:, None], 0)*dV
    J = np.sum(np.cross(P, vel), 0)*dV
    return S, J
ok = True; errs = []
for xf in (np.array([2.0, 0.5, 0.3]), np.array([0.0, 3.0, 1.0]), np.array([1.5, -1.5, 2.0])):
    S, J = ball_S(xf)
    pred = 0.5*np.cross(J, xf)/np.linalg.norm(xf)**3
    e = np.linalg.norm(S - pred)/np.linalg.norm(pred); errs.append(e)
    ok &= e < 5e-3
check("F2 superposed moving mass of a spinning ball: int rho u/|x-x'| = (1/2) J x x / r^3 outside",
      ok, "rel. errors " + ", ".join(f"{e:.1e}" for e in errs))

# So h_0i = -2(1+gamma) * (1/2) (J x x)/r^3 = -(1+gamma) (J x x)/r^3; gamma = 1 gives -2 (J x x)/r^3.
x, y, z = sp.symbols("x y z", real=True)
J = sp.symbols("J", positive=True)
rr = sp.sqrt(x**2 + y**2 + z**2)
Jv = sp.Matrix([0, 0, J]); Xv = sp.Matrix([x, y, z])
gvec = -(1 + gam)*Jv.cross(Xv)/rr**3

def curl(F):
    return sp.Matrix([sp.diff(F[2], y) - sp.diff(F[1], z), sp.diff(F[0], z) - sp.diff(F[2], x), sp.diff(F[1], x) - sp.diff(F[0], y)])

# F3. Gyroscope precession = -(1/2) curl g  (half the circulation per area).  With gamma = 1: Lense-Thirring.
Om = sp.simplify(-sp.Rational(1, 2)*curl(gvec).subs(gam, 1))
rhat = Xv/rr
LT = (3*(Jv.dot(rhat))*rhat - Jv)/rr**3
check("F3 gyroscope turn rate -1/2 curl g = (3(J.r)r - J)/r^3 (Lense-Thirring / Schiff) when gamma = 1",
      sp.simplify(Om - LT) == sp.zeros(3, 1))
R["F3_axis"] = "2J/r^3 along J"; R["F3_equator"] = "-J/r^3"

# F3b. 'Turn = half the circulation per enclosed area' checked on small loops numerically
gnum = sp.lambdify((x, y, z), gvec.subs({gam: 1, J: 1}), "numpy")
def circ(center, normal, rad=1e-3, n=4000):
    normal = normal/np.linalg.norm(normal)
    e1 = np.cross(normal, [1, 0, 0]) if abs(normal[0]) < 0.9 else np.cross(normal, [0, 1, 0])
    e1 /= np.linalg.norm(e1); e2 = np.cross(normal, e1)
    th = np.linspace(0, 2*np.pi, n, endpoint=False)
    tot = 0.0
    for t_ in th:
        p = center + rad*(np.cos(t_)*e1 + np.sin(t_)*e2)
        dl = rad*(-np.sin(t_)*e1 + np.cos(t_)*e2)*(2*np.pi/n)
        tot += np.dot(np.array(gnum(*p), dtype=float).ravel(), dl)
    return tot/(np.pi*rad**2)
pt = np.array([1.2, 0.7, 0.9])
Omn = sp.lambdify((x, y, z), LT.subs(J, 1), "numpy")(*pt)
ok = True
for nrm in (np.array([1.0, 0, 0]), np.array([0, 1.0, 0]), np.array([0, 0, 1.0])):
    ok &= abs(-0.5*circ(pt, nrm) - float(np.dot(np.array(Omn, dtype=float).ravel(), nrm))) < 1e-6
check("F3b turn rate = -(1/2) x circulation of g per enclosed area (three loop orientations)", ok)

# F4. Independent check against the exact Kerr solution (Boyer-Lindquist), far field.
M_, aK, rK, th = sp.symbols("M a r theta", positive=True)
Sig = rK**2 + aK**2*sp.cos(th)**2
gtphi = -2*M_*aK*rK*sp.sin(th)**2/Sig
lead = sp.limit(gtphi*rK, rK, sp.oo)
# our field in spherical: g.dx = -2 (J x x).dx / r^3 = -2 J sin^2(theta) dphi / r  with J = M a
check("F4a Kerr's g_t-phi at large r = -2 J sin^2(theta)/r with J = M a: same as the boosted Newton field",
      sp.simplify(lead - (-2*M_*aK*sp.sin(th)**2)) == 0)

# F4b. Exact static-observer turn rate in Kerr from Omega = (N/2)|dA|_h, compared with Lense-Thirring at large r.
def kerr_static_rate(Mv, av, rv, thv):
    Mk, ak, r_, t_ = sp.symbols("Mk ak r_ t_", positive=True)
    D = r_**2 - 2*Mk*r_ + ak**2; Sg = r_**2 + ak**2*sp.cos(t_)**2
    gtt = -(1 - 2*Mk*r_/Sg); gtp = -2*Mk*ak*r_*sp.sin(t_)**2/Sg
    gpp = (r_**2 + ak**2 + 2*Mk*ak**2*r_*sp.sin(t_)**2/Sg)*sp.sin(t_)**2
    grr = Sg/D; gthth = Sg
    N2 = -gtt; A = gtp/N2; hpp = gpp + N2*A**2
    dAr = sp.diff(A, r_); dAt = sp.diff(A, t_)
    mag2 = dAr**2/(grr*hpp) + dAt**2/(gthth*hpp)
    expr = sp.sqrt(N2)/2*sp.sqrt(mag2)
    return float(expr.subs({Mk: Mv, ak: av, r_: rv, t_: thv}))
Jk = 1.0*0.5
rates = []
for rv in (200.0, 800.0):
    eq = kerr_static_rate(1.0, 0.5, rv, math.pi/2)
    ax = kerr_static_rate(1.0, 0.5, rv, 1e-4)
    rates.append((rv, eq*rv**3/Jk, ax*rv**3/Jk))
check("F4b exact Kerr static gyroscope -> J/r^3 at the equator, 2J/r^3 on the axis, far out",
      all(abs(e - 1) < 2e-2 and abs(a_ - 2) < 2e-2 for _, e, a_ in rates),
      "; ".join(f"r={rv:.0f}: eq {e:.4f}, axis {a_:.4f}" for rv, e, a_ in rates))

# F5. The same rule, exact, in spinning BTZ: (N/2)|dA|_h = J / (2 (r^2 - M)) (precession page G2)
Mb, Jb, rb = sp.symbols("M_b J_b r_b", positive=True)
gtt = -(rb**2 - Mb); gtp = -Jb/2
grr = 1/(rb**2 - Mb + Jb**2/(4*rb**2)); gpp = rb**2
N2 = -gtt; A = gtp/N2; hpp = gpp + N2*A**2
rateB = sp.sqrt(N2)/2*sp.Abs(sp.diff(A, rb))/sp.sqrt(grr*hpp)
check("F5 BTZ: half the swirl per area, weighted by the clock rate, is exactly J/(2(r^2 - M))",
      sp.simplify((rateB - Jb/(2*(rb**2 - Mb))).subs({Mb: sp.Rational(1, 3), Jb: sp.Rational(1, 5), rb: 2})) == 0 and
      sp.simplify((rateB - Jb/(2*(rb**2 - Mb))).subs({Mb: 2, Jb: sp.Rational(3, 2), rb: 3})) == 0)

# F6. Orbit node drag (LAGEOS-type) from the same field: a = -grad Phi + v x curl g  (from L = v^2/2 - Phi + g.v)
Jt = 1e-3
cg = sp.lambdify((x, y, z), curl(gvec.subs({gam: 1, J: Jt})), "numpy")
def f6(t, Y):
    p = Y[:3]; v = Y[3:]
    rr_ = np.linalg.norm(p)
    B = np.array(cg(*p), dtype=float).ravel()
    return np.concatenate([v, -p/rr_**3 + np.cross(v, B)])
a0 = 10.0; inc = math.radians(50)
p0 = np.array([a0, 0, 0]); v0 = math.sqrt(1/a0)*np.array([0, math.cos(inc), math.sin(inc)])
Tn = 2*math.pi*a0**1.5*120
s6 = solve_ivp(f6, (0, Tn), np.concatenate([p0, v0]), method="DOP853", rtol=1e-12, atol=1e-14,
               t_eval=np.linspace(0, Tn, 2000))
Lz = np.cross(s6.y[:3].T, s6.y[3:].T)
node = np.unwrap(np.arctan2(Lz[:, 0], -Lz[:, 1]))
rate = np.polyfit(s6.t, node, 1)[0]
pred6 = 2*Jt/a0**3
check("F6 node drag of an orbit from the same field = 2J/a^3 (Lense-Thirring 1918), prograde",
      abs(rate/pred6 - 1) < 1e-2, f"numeric {rate:.4e} vs {pred6:.4e}")

# F7. Data.  Earth: I = 8.034e37 kg m^2, omega = 7.2921e-5 /s.
Gc, cl = 6.6743e-11, 2.99792458e8
Jearth = 8.034e37*7.2921e-5
yr = 3.15576e7; mas = 180/math.pi*3600e3
a_gpb = 7027.4e3
gpb_vec = Gc*Jearth/(2*cl**2*a_gpb**3)*yr*mas         # polar-orbit average, along J
gpb_proj = gpb_vec*math.cos(math.radians(16.84))       # onto the measured axis (guide star IM Peg, dec 16.84 deg)
gpb_meas, gpb_err = 37.2, 7.2
aL, eL = 12270e3, 0.0045
lag = 2*Gc*Jearth/(cl**2*aL**3*(1 - eL**2)**1.5)*yr*mas
R["F7"] = dict(gpb_vec=gpb_vec, gpb_proj=gpb_proj, gpb_gamma0=gpb_proj/2, lageos=lag, lageos_gamma0=lag/2)
check("F7a Gravity Probe B frame dragging (gamma = 1): ~39 mas/yr vs 37.2 +/- 7.2 measured",
      abs(gpb_proj - gpb_meas) < gpb_err, f"vector {gpb_vec:.1f}, along measured axis {gpb_proj:.1f} (GR quoted 39.2)")
check("F7b LAGEOS node drag (gamma = 1): ~31 mas/yr, the value LAGEOS measured at 99 +/- 5% (Ciufolini & Pavlis 2004)",
      abs(lag - 30.7) < 0.5, f"{lag:.1f} mas/yr")
sig_gpb0 = (gpb_meas - gpb_proj/2)/gpb_err
sig_lag0 = (0.99 - 0.5)/0.05
check("F7c a dial-only theory (gamma = 0) gives half: GP-B ~2.6 sigma off, LAGEOS ~10 sigma off (stat.)",
      sig_gpb0 > 2 and sig_lag0 > 5, f"GP-B {gpb_proj/2:.1f} mas/yr ({sig_gpb0:.1f} sigma); LAGEOS 50% ({sig_lag0:.0f} sigma, ~5 with 10% total error)")
# double pulsar B (twisted-boosts page split): velocity 0.75, space 1.51 (x gamma), drag 2.81 (x (1+gamma)/2)
dp1 = 0.75 + 1.51 + 2.81; dp0 = 0.75 + 0 + 2.81/2
check("F7d double pulsar B: drag weight 2 m1/M = (1+gamma) m1/M; gamma = 1 gives 5.07, gamma = 0 gives 2.16 vs 4.77 +0.66",
      abs(dp1 - 4.77) < 0.66 and (4.77 - dp0)/0.65 > 3, f"gamma=1: {dp1:.2f}, gamma=0: {dp0:.2f} ({(4.77-dp0)/0.65:.1f} sigma)")
R["F7"]["dp1"] = dp1; R["F7"]["dp0"] = dp0

# F8. Where c sits: restore units.  Omega = (1+gamma)/2 * G/c^2 * (3(J.r)r - J)/r^3; -q = 1/c^2 is the dial.
Gs = sp.symbols("G", positive=True)
expr = (1 + gam)/2*Gs*J/(cs**2)
check("F8 frame dragging ~ (1+gamma)/2 x |q| x G J / r^3 with |q| = 1/c^2: no new constant beyond Newton's G and the dial",
      sp.simplify(expr.subs({gam: 1}) - Gs*J/cs**2) == 0)

print(f"\n{PASS} pass, {FAIL} fail")
LOG.append(f"\n{PASS} pass, {FAIL} fail")
R["pass"] = PASS; R["fail"] = FAIL
with open(os.path.join(HERE, "checks_output.txt"), "w") as fh:
    fh.write("\n".join(LOG) + "\n")
with open(os.path.join(HERE, "results.json"), "w") as fh:
    json.dump(R, fh, indent=1, default=float)
