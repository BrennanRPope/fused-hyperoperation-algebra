"""Quantum action: numerical checks (DRAFT, Nov 26 2026).

Run:  python3 checks.py   (numpy, scipy, sympy, mpmath)  -> checks_output.txt, results.json
Sections
  A  action is area: the swap rule in phase space, spin spheres, contraction to hbar
  B  action is length: one dial for Fermat, Newton and Einstein; Kepler; fish-eye; Snell
  C  action shapes: cube (spins), simplex (oscillators, Born weights), lattice points = states,
     Bohr-Sommerfeld and the Maslov half
  E  gravity: Regge deficit = swap-rule turn, 2D action is topological,
     Ponzano-Regge (6j symbol ~ cos(Regge action + pi/4))
"""
import json, math, os, sys
from fractions import Fraction
import numpy as np
from scipy.linalg import expm, eigh_tridiagonal
from scipy.integrate import solve_ivp
from scipy.optimize import minimize_scalar, brentq
from scipy import stats
import sympy as sp
import mpmath as mp

HERE = os.path.dirname(os.path.abspath(__file__))
out_lines, results = [], {}
npass = nfail = 0

def check(name, ok, detail=""):
    global npass, nfail
    tag = "PASS" if ok else "FAIL"
    npass += ok; nfail += (not ok)
    line = f"[{tag}] {name}: {detail}"
    print(line); out_lines.append(line)

def wrap(a):
    return (a + np.pi) % (2*np.pi) - np.pi

rng = np.random.default_rng(20261126)

# ---------------------------------------------------------------- A. area
# A1: displacements in phase space, applied to a real wavefunction on a grid.
hbar = 0.5
N, L = 4096, 40.0
x = np.linspace(-L/2, L/2, N, endpoint=False); dx = x[1]-x[0]
p = 2*np.pi*hbar*np.fft.fftfreq(N, d=dx)
psi0 = np.exp(-x**2/(2*hbar)) ; psi0 = psi0/np.sqrt(np.sum(abs(psi0)**2)*dx)

def displace(psi, x0, p0):
    # D(x0,p0) = exp(i(p0 X - x0 P)/hbar) = e^{-i x0 p0/2hbar} e^{i p0 X/hbar} e^{-i x0 P/hbar}
    psi = np.fft.ifft(np.exp(-1j*x0*p/hbar)*np.fft.fft(psi))
    return np.exp(-1j*x0*p0/(2*hbar))*np.exp(1j*p0*x/hbar)*psi

def loop_phase(verts):
    psi = psi0.astype(complex)
    for a, b in zip(verts, verts[1:]+verts[:1]):
        psi = displace(psi, b[0]-a[0], b[1]-a[1])
    ov = np.sum(np.conj(psi0)*psi)*dx
    return np.angle(ov), abs(ov)

def shoelace(verts):
    v = np.array(verts); return 0.5*np.sum(v[:,0]*np.roll(v[:,1],-1)-np.roll(v[:,0],-1)*v[:,1])

A1 = []
polys = [[(0,0),(1,0),(0,1)], [(0,0),(1.2,0),(1.2,0.7),(0,0.7)], [(0,0),(0.5,-0.3),(1,0.4),(0.2,0.9),(-0.4,0.3)],
         [(0,0),(0,1),(1,0)], [(0.3,0.2),(1.5,0.2),(0.9,1.4)]]
for P in polys:
    ph, mag = loop_phase([tuple(map(float,v)) for v in P])
    A = shoelace(P)
    A1.append((A, ph, A/hbar, mag))
err = max(abs(wrap(ph - A/hbar)) for A, ph, _, _ in A1)
check("A1 a closed loop of phase-space moves turns the quantum state by (enclosed area)/hbar",
      err < 1e-6 and min(m for *_, m in A1) > 1-1e-9,
      f"5 polygons incl. reversed and pentagon; worst phase error {err:.1e} rad; |overlap|=1")
results["A1"] = [dict(area=A, phase=float(ph), area_over_hbar=r) for A, ph, r, _ in A1]

# A1b: swapping two moves costs exactly their parallelogram area / hbar
a, b = (0.7, 0.2), (-0.3, 0.9)
psiab = displace(displace(psi0.astype(complex), *a), *b)
psiba = displace(displace(psi0.astype(complex), *b), *a)
sw = np.angle(np.sum(np.conj(psiba)*psiab)*dx)
area_ab = a[0]*b[1]-a[1]*b[0]
check("A1b swap rule in phase space: D(a)D(b) = e^{i(a x b)/hbar} D(b)D(a)",
      abs(wrap(sw - area_ab/hbar)) < 1e-8, f"swap phase {sw:.6f}, area/hbar {area_ab/hbar:.6f}")

# A2: spin-j coherent states, Berry phase of a loop = j * solid angle
def spin_mats(j):
    m = np.arange(j, -j-1, -1)
    Jz = np.diag(m)
    jp = np.sqrt(j*(j+1) - m[1:]*(m[1:]+1))
    Jp = np.diag(jp, 1)
    Jx = (Jp + Jp.T)/2; Jy = (Jp - Jp.T)/(2j)
    return Jx, Jy, Jz.astype(complex)

def coherent(j, n, mats):
    Jx, Jy, Jz = mats
    th = math.acos(max(-1,min(1,n[2]))); ph = math.atan2(n[1], n[0])
    top = np.zeros(int(round(2*j+1)), complex); top[0] = 1
    return expm(-1j*ph*Jz) @ expm(-1j*th*Jy) @ top

def bargmann(j, pts, mats):
    st = [coherent(j, q, mats) for q in pts]
    prod = 1+0j
    for u, v in zip(st, st[1:]+st[:1]):
        prod *= np.vdot(u, v)
    return -np.angle(prod)

A2 = []
for j in [0.5, 1, 1.5, 3, 7.5]:
    mats = spin_mats(j)
    for th0 in [0.3, 0.9, 1.4]:
        M = 400
        pts = [np.array([math.sin(th0)*math.cos(t), math.sin(th0)*math.sin(t), math.cos(th0)]) for t in np.linspace(0, 2*np.pi, M, endpoint=False)]
        g = bargmann(j, pts, mats)
        Om = 2*np.pi*(1-math.cos(th0))
        A2.append((j, th0, float(g), float(wrap(-j*Om))))
errA2 = max(abs(wrap(g - e)) for *_, g, e in A2)
check("A2 spin j: a coherent state carried round a loop picks up j x (solid angle) (Berry 1984)",
      errA2 < 2e-3, f"15 loops, j = 1/2..15/2; worst error {errA2:.1e} rad (400-step loops)")
results["A2"] = A2

# A3: contraction. Same small square in flat phase space vs on the spin sphere of radius^2 = j*hbar.
def sph_from_tangent(u, v, R):
    r = math.hypot(u, v); c = r/R
    if r == 0: return np.array([0, 0, 1.0])
    return np.array([math.sin(c)*u/r, math.sin(c)*v/r, math.cos(c)])

def slerp_path(P, Q, n):
    om = math.acos(max(-1, min(1, P @ Q)))
    return [ (math.sin((1-t)*om)*P + math.sin(t*om)*Q)/math.sin(om) for t in np.linspace(0, 1, n, endpoint=False)]

def sph_poly_area(V):
    # Girard: sum of interior angles - (n-2) pi, for a small convex geodesic polygon on the unit sphere
    n = len(V); s = 0
    for i in range(n):
        A, B, C = V[i-1], V[i], V[(i+1) % n]
        t1 = A - (A@B)*B; t2 = C - (C@B)*B
        s += math.acos(max(-1, min(1, t1@t2/np.linalg.norm(t1)/np.linalg.norm(t2))))
    return s - (n-2)*math.pi

hb = 1.0; side = 1.2; Aflat = side**2
sq = [(-side/2, -side/2), (side/2, -side/2), (side/2, side/2), (-side/2, side/2)]
A3 = []
for j in [2, 5, 10, 25, 60]:
    R = math.sqrt(j*hb)
    V = [sph_from_tangent(u, v, R) for u, v in sq]
    path = sum((slerp_path(V[i], V[(i+1) % 4], 60) for i in range(4)), [])
    mats = spin_mats(j)
    phase = bargmann(j, path, mats)            # quantum phase (radians)
    Omega = sph_poly_area(V)                    # holonomy turn = curvature x area on unit sphere
    A3.append(dict(j=j, d=1/R**2, phase=float(phase), turn=float(Omega), flat_area_over_hbar=Aflat/hb,
                   ratio_phase_turn=float(-phase/Omega)))
ok_ratio = all(abs(r["ratio_phase_turn"] - r["j"]) < 0.02*r["j"] for r in A3)
dev = [abs(-r["phase"] - Aflat/hb) for r in A3]
check("A3a quantum phase = j x holonomy turn on every spin sphere (phase per turn is the spin)",
      ok_ratio, ", ".join(f"j={r['j']}: {r['ratio_phase_turn']:.3f}" for r in A3))
slope = np.polyfit(np.log([r["j"] for r in A3]), np.log(dev), 1)[0]
check("A3b as the sphere flattens (d = 1/(j hbar) -> 0) the turn vanishes but the phase -> flat area/hbar, gap ~ 1/j",
      dev[-1] < 0.01 and abs(slope + 1) < 0.15,
      f"phase gap from A/hbar: {', '.join(f'{d:.4f}' for d in dev)}; log-log slope {slope:.2f}")
results["A3"] = A3

# A4: [J_x/sqrt(j), J_y/sqrt(j)] -> i near the pole  (Inonu-Wigner / Holstein-Primakoff)
j = 40; Jx, Jy, Jz = spin_mats(j)
X, P_ = Jx/math.sqrt(j), Jy/math.sqrt(j)
C = X@P_ - P_@X
vals = [ (C[m, m]/1j).real for m in range(4)]
check("A4 near the top of a big spin, [X, P] = i(1 - m/j): the Heisenberg rule is the flattened swap rule",
      all(abs(v - (1 - m/j)) < 1e-12 for m, v in enumerate(vals)), f"j=40, first diagonal entries {np.round(vals, 4).tolist()}")

# ---------------------------------------------------------------- B. length
def S_dial(t, xx, m, d):
    dt, dxx = np.diff(t), np.diff(xx)
    if d == 0:
        return np.sum(0.5*m*dxx**2/dt)
    return (m/d)*np.sum(np.sqrt(dt**2 + d*dxx**2))

T = 2.0; tt = np.linspace(0, T, 20001); x_line = 0.3*tt
B1 = {}
for d in [0.5, 0.0, -0.5]:
    base = S_dial(tt, x_line, 1.0, d)
    eps_vals = [-0.1, -0.05, 0.05, 0.1]
    dS = [S_dial(tt, x_line + e*np.sin(np.pi*tt/T), 1.0, d) - base for e in eps_vals]
    if d == 0:
        geo = None
    else:
        geo = [np.sum(np.sqrt(np.diff(tt)**2 + d*np.diff(x_line + e*np.sin(np.pi*tt/T))**2)) - np.sum(np.sqrt(np.diff(tt)**2 + d*np.diff(x_line)**2)) for e in eps_vals]
    B1[d] = dict(dS=dS, dlength=geo)
okB1 = all(min(v["dS"]) > 0 for v in B1.values()) and min(B1[0.5]["dlength"]) > 0 and max(B1[-0.5]["dlength"]) < 0
check("B1a one action S_d = (m/d) int sqrt(1 + d v^2) dt is least on the straight path at d = +, 0, -",
      okB1, "every wiggle raises S; the length it measures goes UP for d>0 (Fermat) and DOWN for d<0 (proper time: twin paradox)")
results["B1"] = {str(k): v for k, v in B1.items()}

# B1b: d -> 0 limit is Newton
xw = 0.3*tt + 0.2*np.sin(2*np.pi*tt/T)
lim = [S_dial(tt, xw, 1.0, d) - 1.0*T/d for d in [-1e-2, -1e-3, -1e-4]]
newt = S_dial(tt, xw, 1.0, 0)
check("B1b S_d - mT/d -> Newton's int (1/2) m v^2 dt as the dial d -> 0",
      abs(lim[-1] - newt) < 1e-4 and abs(lim[0]-newt) > abs(lim[-1]-newt), f"{[round(v,6) for v in lim]} -> {newt:.6f}")

# B1c: the symmetry of S_d composes velocities by our plane law (u+v)/(1-d u v)
def boost(t, xx, u, d):
    g = 1/math.sqrt(1 + d*u*u)
    return g*(t + d*u*xx), g*(xx - u*t)
okc = True; worst = 0
for d in [0.7, -0.6]:
    for u in [0.3, -0.5]:
        t2, x2 = boost(tt, xw, u, d)
        s1 = np.sum(np.sqrt(np.diff(tt)**2 + d*np.diff(xw)**2)); s2 = np.sum(np.sqrt(np.diff(t2)**2 + d*np.diff(x2)**2))
        worst = max(worst, abs(s1 - s2))
uu, vv, dd = sp.symbols("u v d")
# velocity of the line x = v t seen after boost by u is (v-u)/(1+d u v); so boosts compose by w = (u+v)/(1 - d u v)
t_ = sp.symbols("t")
tb = (t_ + dd*uu*vv*t_); xb = (vv*t_ - uu*t_)
vel = sp.simplify(xb/tb)
comp_ok = sp.simplify(vel - (vv-uu)/(1+dd*uu*vv)) == 0
check("B1c the action's symmetry (rotation for d>0, Lorentz boost for d<0) adds velocities by (u+v)/(1-duv), our plane law at k=0",
      worst < 1e-9 and comp_ok, f"length unchanged to {worst:.1e}; velocity law checked symbolically")

# B1d: k is a slide. (x+y+kxy)/(1-dxy) is conjugate to (x+y)/(1-d'xy) by x -> x/(1+ax)
xs, ys, kk, a_ = sp.symbols("x y k a")
law = lambda X, Y, k, d: (X+Y+k*X*Y)/(1-d*X*Y)
f = lambda X: X/(1+a_*X)
lhs = f(law(xs, ys, kk, dd))
dp = sp.symbols("dp")
rhs = law(f(xs), f(ys), 0, dp)
num = sp.expand(sp.numer(sp.together(lhs - rhs)))
sol = sp.solve(sp.Poly(num, xs, ys).coeffs(), [a_, dp], dict=True)
check("B1d the k part of the full law is only a change of velocity ruler: x -> x/(1 + (k/2)x) gives d' = d - k^2/4",
      any(sp.simplify(s[a_] - kk/2) == 0 and sp.simplify(s[dp] - (dd - kk**2/4)) == 0 for s in sol), f"solutions {sol}")
results["B1d"] = str(sol)

# B2: Kepler at fixed energy: every orbit has the same action (all ellipses, one action)
def kepler_action(e, E=-0.5, k=1.0, m=1.0):
    a = k/(-2*E); rp = a*(1-e)
    vp = math.sqrt(2*(E + k/rp)/m)
    def f(t, y):
        r3 = (y[0]**2+y[1]**2)**1.5
        return [y[2], y[3], -k*y[0]/r3/m, -k*y[1]/r3/m]
    Pd = 2*np.pi*math.sqrt(m*a**3/k)
    sol = solve_ivp(f, [0, Pd], [rp, 0, 0, vp], rtol=1e-11, atol=1e-12, dense_output=True)
    ts = np.linspace(0, Pd, 20001); Y = sol.sol(ts)
    integrand = m*(Y[2]**2 + Y[3]**2)
    S = np.trapezoid(integrand, ts)
    Lz = m*rp*vp
    return S, 2*np.pi*Lz, Y
B2 = []
for e in [0.0, 0.3, 0.6, 0.9]:
    S, Jth, _ = kepler_action(e)
    B2.append(dict(e=e, action=S, angular_part=Jth, radial_part=S - Jth))
check("B2 Kepler: at one energy every ellipse, round or thin, has the same action 2*pi*k*sqrt(m/-2E) (hydrogen's degeneracy)",
      all(abs(r["action"] - 2*np.pi) < 1e-5 for r in B2),
      "; ".join(f"e={r['e']}: {r['action']:.6f} (angular {r['angular_part']:.3f} + radial {r['radial_part']:.3f})" for r in B2))
results["B2"] = B2

# B3: Maxwell's fish-eye n = 2/(1+r^2): rays are circles, and lift to great circles of the sphere
def fisheye_ray(x0, y0, ang, s_max=6.0):
    n = lambda X, Y: 2/(1+X*X+Y*Y)
    def f(s, y):
        X, Y, px, py = y
        nn = n(X, Y)
        gx = -4*X/(1+X*X+Y*Y)**2; gy = -4*Y/(1+X*X+Y*Y)**2
        return [px/nn**2*nn, py/nn**2*nn, gx/nn*nn, gy/nn*nn]
    # Hamiltonian form: dr/ds = p/n, dp/ds = grad n, |p| = n
    def g(s, y):
        X, Y, px, py = y; nn = n(X, Y)
        gx = -4*X/(1+X*X+Y*Y)**2; gy = -4*Y/(1+X*X+Y*Y)**2
        return [px/nn, py/nn, gx, gy]
    nn0 = n(x0, y0)
    sol = solve_ivp(g, [0, s_max], [x0, y0, nn0*math.cos(ang), nn0*math.sin(ang)], rtol=1e-11, atol=1e-12, max_step=0.01)
    return sol.y[0], sol.y[1]
worst_circ = worst_plane = 0
rays = []
for (x0, y0, ang) in [(0.5, 0, 1.2), (-0.2, 0.7, -0.4), (1.5, -0.5, 2.5), (0, 0, 0.7)]:
    X, Y = fisheye_ray(x0, y0, ang); rays.append((X.tolist()[::20], Y.tolist()[::20]))
    A_ = np.c_[2*X, 2*Y, np.ones_like(X)]; bb = X**2 + Y**2
    cfit, *_ = np.linalg.lstsq(A_, bb, rcond=None)
    if abs(x0) + abs(y0) > 0:
        resid = np.sqrt((X-cfit[0])**2 + (Y-cfit[1])**2) - math.sqrt(cfit[2] + cfit[0]**2 + cfit[1]**2)
        worst_circ = max(worst_circ, np.max(abs(resid)))
    rr = X**2 + Y**2
    S3 = np.c_[2*X/(1+rr), 2*Y/(1+rr), (1-rr)/(1+rr)]
    sv = np.linalg.svd(S3, compute_uv=False)
    worst_plane = max(worst_plane, sv[2]/sv[0])
check("B3 Maxwell's fish-eye lens: light rays are circles, and lifted to the sphere they are great circles",
      worst_circ < 1e-6 and worst_plane < 1e-6, f"circle fit error {worst_circ:.1e}; third singular value ratio {worst_plane:.1e} (plane through centre)")
results["B3_rays"] = rays

# B4: Snell from least optical path
n1, n2 = 1.0, 1.5; Apt, Bpt = (0, 1), (1, -1)
opl = lambda xi: n1*math.hypot(xi-Apt[0], Apt[1]) + n2*math.hypot(Bpt[0]-xi, Bpt[1])
xi = minimize_scalar(opl, bounds=(0, 1), method="bounded", options={"xatol": 1e-12}).x
s1 = (xi-Apt[0])/math.hypot(xi-Apt[0], Apt[1]); s2 = (Bpt[0]-xi)/math.hypot(Bpt[0]-xi, Bpt[1])
check("B4 Fermat: the least optical path through a surface obeys Snell's law n1 sin a1 = n2 sin a2",
      abs(n1*s1 - n2*s2) < 1e-6, f"{n1*s1:.10f} vs {n2*s2:.10f}")

# ---------------------------------------------------------------- C. shapes
# C1: Archimedes. Points spread evenly on a product of spheres land evenly in the cube of their heights.
Np, Nsp = 200000, 3
v = rng.normal(size=(Np, Nsp, 3)); v /= np.linalg.norm(v, axis=2, keepdims=True)
z = v[:, :, 2]
ks = [stats.kstest(z[:, i], "uniform", args=(-1, 2)).pvalue for i in range(Nsp)]
hist, _ = np.histogramdd(z, bins=(4, 4, 4), range=[(-1, 1)]*3)
chi = stats.chisquare(hist.ravel()).pvalue
check("C1a three spins' actions (heights J_z) fill the cube evenly (Archimedes' hat-box, Duistermaat-Heckman)",
      min(ks) > 1e-3 and chi > 1e-3, f"KS p-values {np.round(ks, 3).tolist()}, 64-cell chi-square p {chi:.3f}")
counts = [(int(2*j+1))**3 for j in [0.5, 1, 1.5, 2]]
check("C1b quantum states of 3 spins = whole-step lattice points of the cube [-j, j]^3: (2j+1)^3",
      counts == [8, 27, 64, 125], f"{counts}")

# C2: random pure quantum states: Born weights |psi_i|^2 fill the simplex evenly
Nl = 4
psi = rng.normal(size=(Np, Nl)) + 1j*rng.normal(size=(Np, Nl)); psi /= np.linalg.norm(psi, axis=1, keepdims=True)
w = abs(psi)**2
kspv = stats.kstest(w[:, 0], stats.beta(1, Nl-1).cdf).pvalue
pair = stats.kstest(w[:, 0] + w[:, 1], stats.beta(2, Nl-2).cdf).pvalue
check("C2a Born weights of random 4-level pure states fill the simplex evenly (flat Dirichlet)",
      kspv > 1e-3 and pair > 1e-3, f"KS p (one weight ~ Beta(1,3)) {kspv:.3f}, (two weights ~ Beta(2,2)) {pair:.3f}")
from math import comb
okC2 = True; tab = []
for Nosc in [2, 3, 4]:
    for nq in [0, 1, 5, 12]:
        cnt = sum(1 for mm in np.ndindex(*(nq+1,)*Nosc) if sum(mm) == nq)
        okC2 &= cnt == comb(nq+Nosc-1, Nosc-1); tab.append((Nosc, nq, cnt))
check("C2b states of N oscillators sharing n quanta = lattice points on the simplex face: C(n+N-1, N-1)",
      okC2, f"{tab[:6]} ...")
# phase-space volume / (2 pi hbar)^N vs count (Weyl law with the Maslov shift)
Nosc = 3; E = 40
vol_states = (E)**Nosc/math.factorial(Nosc)         # {sum H_i <= E}, hbar=omega=1: volume/(2pi)^N = E^N/N!
cnt_le = comb(int(E - Nosc/2) + Nosc, Nosc)          # states with sum(n_i + 1/2) <= E
check("C2c the action simplex's volume counts the states (Weyl), once each oscillator's half quantum is paid",
      abs(cnt_le/vol_states - 1) < 0.02, f"states {cnt_le} vs volume {vol_states:.1f}")

# C3: Bohr-Sommerfeld with the Maslov half: quartic well V = x^4, m = 1/2 (H = p^2 + x^4), hbar = 1
Ng, Lg = 6000, 12.0
xg = np.linspace(-Lg/2, Lg/2, Ng); h = xg[1]-xg[0]
diag = 2/h**2 + xg**4; off = -np.ones(Ng-1)/h**2
Eex = eigh_tridiagonal(diag, off, select="i", select_range=(0, 11))[0]
def loop_action(E):  # closed-loop area  = 2 * int p dx, p = sqrt(E - x^4)
    xt = E**0.25
    return 2*float(mp.quad(lambda X: mp.sqrt(max(E - X**4, 0)), [-xt, xt]))
def bs_level(n, shift):
    return brentq(lambda E: loop_action(E) - 2*np.pi*(n + shift), 1e-6, 1e4)
E_half = np.array([bs_level(n, 0.5) for n in range(12)]); E_zero = np.array([bs_level(n, 0.0) if n > 0 else np.nan for n in range(12)])
rel_half = abs(E_half - Eex)/Eex; rel_zero = abs(E_zero - Eex)/Eex
check("C3a quartic well: (loop area) = 2 pi hbar (n + 1/2) gives the true levels, error falling every step (18% at n=0, as usual for WKB)",
      all(np.diff(rel_half) < 0) and rel_half[11] < 1e-3, f"relative errors n=0..4: {np.round(rel_half[:5], 4).tolist()}, n=11: {rel_half[11]:.1e}")
check("C3b without the half (the loop's own turning) the levels are badly off",
      np.all(rel_zero[1:12] > 10*rel_half[1:12]), f"errors without 1/2, n=1..5: {np.round(rel_zero[1:6], 3).tolist()}")
results["C3"] = dict(exact=Eex.tolist(), bs_half=E_half.tolist(), bs_zero=E_zero.tolist())
# the half = 2 quarter turns: the loop's tangent is vertical twice per lap (Maslov index 2)
th = np.linspace(0, 2*np.pi, 4001)[:-1]
Ec = 5.0
rr_ = np.array([brentq(lambda r: r*r*math.sin(t)**2 + r**4*math.cos(t)**4 - Ec, 1e-9, 10) for t in th])
Xl, Pl = rr_*np.cos(th), rr_*np.sin(th)
dX = np.roll(Xl, -1) - Xl
sg = np.sign(dX); vert = int(np.sum(sg != np.roll(sg, -1)))
check("C3c the energy loop has exactly 2 vertical tangents per lap (Maslov index 2 = two quarter turns = the half quantum)",
      vert == 2, f"sign changes of dx along the loop: {vert}")

# ---------------------------------------------------------------- E. gravity
# E1: Regge. On a geodesic polyhedron approximating a sphere, the deficit angle at a vertex = curvature x its area.
def icosphere(nsub):
    t = (1+5**0.5)/2
    V = [(-1,t,0),(1,t,0),(-1,-t,0),(1,-t,0),(0,-1,t),(0,1,t),(0,-1,-t),(0,1,-t),(t,0,-1),(t,0,1),(-t,0,-1),(-t,0,1)]
    V = [np.array(v)/np.linalg.norm(v) for v in V]
    F = [(0,11,5),(0,5,1),(0,1,7),(0,7,10),(0,10,11),(1,5,9),(5,11,4),(11,10,2),(10,7,6),(7,1,8),(3,9,4),(3,4,2),(3,2,6),(3,6,8),(3,8,9),(4,9,5),(2,4,11),(6,2,10),(8,6,7),(9,8,1)]
    for _ in range(nsub):
        cache = {}; NF = []
        def mid(i, j):
            key = (min(i,j), max(i,j))
            if key not in cache:
                m = V[i]+V[j]; V.append(m/np.linalg.norm(m)); cache[key] = len(V)-1
            return cache[key]
        for a, b, c in F:
            ab, bc, ca = mid(a,b), mid(b,c), mid(c,a)
            NF += [(a,ab,ca),(b,bc,ab),(c,ca,bc),(ab,bc,ca)]
        F = NF
    return np.array(V), F
E1 = []
for R in [1.0, 2.0]:
    for nsub in [1, 2, 3]:
        V, F = icosphere(nsub); V = V*R
        ang = np.zeros(len(V)); area = np.zeros(len(V))
        for tri in F:
            P = V[list(tri)]
            Ar = 0.5*np.linalg.norm(np.cross(P[1]-P[0], P[2]-P[0]))
            for i in range(3):
                u = P[(i+1)%3]-P[i]; w_ = P[(i+2)%3]-P[i]
                ang[tri[i]] += math.acos(u@w_/np.linalg.norm(u)/np.linalg.norm(w_)); area[tri[i]] += Ar/3
        deficit = 2*np.pi - ang
        E1.append(dict(R=R, nsub=nsub, total=float(deficit.sum()), ratio=float(np.median(deficit/area)*R**2)))
check("E1a Regge: the deficit angles of any triangulated sphere add to 4 pi (Gauss-Bonnet, exact)",
      all(abs(r["total"] - 4*np.pi) < 1e-9 for r in E1), f"totals {[round(r['total'],12) for r in E1]}")
check("E1b each deficit (the turn a vector makes going round the vertex) = (1/R^2) x the vertex's area: the swap rule's d x area",
      all(abs(r["ratio"] - 1) < 0.02 for r in E1 if r["nsub"] >= 2), f"median deficit/area x R^2: {[round(r['ratio'],4) for r in E1]}")
results["E1"] = E1

# E1c: 2D gravity's action does not move: random edge flips change every deficit but not their sum
V, F = icosphere(1); F = [tuple(f) for f in F]
def deg_deficits(F, nV):
    deg = np.zeros(nV)
    for f_ in F:
        for i in f_: deg[i] += 1
    return 2*np.pi - (np.pi/3)*deg       # equilateral triangles (dynamical triangulations)
d0 = deg_deficits(F, len(V)); sums = [d0.sum()]; changed = 0
for it in range(300):
    fi = rng.integers(len(F)); a, b, c = F[fi]
    e = (a, b)
    others = [k for k, g in enumerate(F) if k != fi and a in g and b in g]
    if not others: continue
    g = F[others[0]]; dpt = [q for q in g if q not in (a, b)][0]
    deg = np.zeros(len(V))
    for f_ in F:
        for i in f_: deg[i] += 1
    if deg[a] <= 3 or deg[b] <= 3: continue
    if any(c in f_ and dpt in f_ for f_ in F): continue
    F[fi] = (a, dpt, c); F[others[0]] = (b, c, dpt)
    dd_ = deg_deficits(F, len(V)); sums.append(dd_.sum()); changed += np.any(dd_ != d0)
check("E1c 2D: after hundreds of random flips of the triangulated universe, the Regge action (sum of deficits) stays 4 pi",
      max(abs(np.array(sums) - 4*np.pi)) < 1e-9 and changed > 50, f"{len(sums)-1} flips, local curvature changed {changed} times, total never moved")

# E2: Ponzano-Regge. Wigner 6j ~ cos(sum (j+1/2) theta + pi/4)/sqrt(12 pi V)
def fact(n): return math.factorial(n)
def delta(a, b, c):
    return Fraction(fact(a+b-c)*fact(a-b+c)*fact(-a+b+c), fact(a+b+c+1))
def sixj(j1, j2, j3, j4, j5, j6):
    # integer spins only here (all sums integral); Racah formula, exact
    tri = [(j1,j2,j3),(j1,j5,j6),(j4,j2,j6),(j4,j5,j3)]
    pref = mp.mpf(1)
    for t_ in tri:
        dl = delta(*t_); pref *= mp.sqrt(mp.mpf(dl.numerator)/dl.denominator)
    a = [sum(t_) for t_ in tri]
    b = [j1+j2+j4+j5, j2+j3+j5+j6, j3+j1+j6+j4]
    s = 0
    for t in range(max(a), min(b)+1):
        den = 1
        for ai in a: den *= fact(t-ai)
        for bi in b: den *= fact(bi-t)
        s += (-1)**t * Fraction(fact(t+1), den)
    return pref*mp.mpf(s.numerator)/s.denominator

def tetra(l):
    # l = (BC, CA, AB, AD, BD, CD) as j1..j6 (+1/2); returns volume and interior dihedral at each edge
    BC, CA, AB, AD, BD, CD = l
    A = np.zeros(3); B = np.array([AB, 0, 0])
    cx = (AB**2 + CA**2 - BC**2)/(2*AB); Cp = np.array([cx, math.sqrt(CA**2 - cx**2), 0])
    dxx = (AB**2 + AD**2 - BD**2)/(2*AB)
    dy = (AD**2 - CD**2 + (Cp[0])**2 + Cp[1]**2 - 2*Cp[0]*dxx)/(2*Cp[1])
    dz2 = AD**2 - dxx**2 - dy**2
    if dz2 <= 0: return None
    D = np.array([dxx, dy, math.sqrt(dz2)])
    P = {"A": A, "B": B, "C": Cp, "D": D}
    vol = abs(np.dot(B - A, np.cross(Cp - A, D - A)))/6
    def dihedral(e, others):
        p, q = P[e[0]], P[e[1]]; ax = (q-p)/np.linalg.norm(q-p)
        r1, r2 = P[others[0]] - p, P[others[1]] - p
        r1 -= (r1@ax)*ax; r2 -= (r2@ax)*ax
        return math.acos(r1@r2/np.linalg.norm(r1)/np.linalg.norm(r2))
    edges = [("BC","AD"), ("CA","BD"), ("AB","CD"), ("AD","BC"), ("BD","CA"), ("CD","AB")]
    th = [dihedral(e, o) for e, o in edges]
    return vol, th

def ponzano_regge(js):
    l = [jj + 0.5 for jj in js]
    tv = tetra(l)
    if tv is None: return None
    vol, th = tv
    S = sum(li*(math.pi - ti) for li, ti in zip(l, th))   # Regge action: length x exterior angle
    return math.cos(S + math.pi/4)/math.sqrt(12*math.pi*vol), math.cos(S)/math.sqrt(12*math.pi*vol), S

mp.mp.dps = 60
base = (6, 7, 8, 7, 6, 8)
E2 = []
for lam in [1, 2, 3, 4, 6, 8]:
    js = tuple(lam*q for q in base)
    ex = float(sixj(*js)); pr, pr0, S = ponzano_regge(js)
    E2.append(dict(lam=lam, js=js, exact=ex, pr=pr, pr_no_maslov=pr0, regge_action=S))
rel = [abs(r["exact"] - r["pr"])/abs(r["pr"]) for r in E2 if abs(r["pr"]) > 0.3*max(abs(q["pr"]) for q in E2)]
check("E2a Ponzano-Regge: 6j symbols of big spins = cos(Regge action + pi/4)/sqrt(12 pi V): 3D quantum gravity's sum gives Einstein's action back",
      max(rel[-2:]) < 0.05, "; ".join(f"x{r['lam']}: exact {r['exact']:.3e}, formula {r['pr']:.3e}" for r in E2))
# scan one spin to see the oscillation and the need for the quarter turn
scan = []
for j6 in range(14, 60):
    js = (30, 31, 32, 31, 30, j6)
    if (30+31+j6) % 1: continue
    pr = ponzano_regge(js)
    if pr is None: continue
    scan.append((j6, float(sixj(*js)), pr[0], pr[1]))
sc = np.array(scan)
mid = (sc[:, 0] > 18) & (sc[:, 0] < 56)
corr_m = np.corrcoef(sc[mid, 1], sc[mid, 2])[0, 1]; corr_0 = np.corrcoef(sc[mid, 1], sc[mid, 3])[0, 1]
check("E2b the quarter turn pi/4 (a Maslov phase) is needed: with it the formula tracks the exact values, without it it doesn't",
      corr_m > 0.98 and corr_0 < 0.8, f"correlation with pi/4: {corr_m:.4f}; without: {corr_0:.4f}  (j = (30,31,32,31,30,j6), {mid.sum()} values)")
results["E2"] = E2; results["E2_scan"] = scan


# E3: Turaev-Viro / Mizoguchi-Tada. Turn on q = e^{i pi/r}: the quantum 6j symbol follows the Regge action of a
# SPHERICAL tetrahedron plus a volume term: 3D quantum gravity with a cosmological constant Lambda ~ (2 pi / r)^2.
def qn(n, r): return mp.sin(n*mp.pi/r)/mp.sin(mp.pi/r)
def qf(n, r):
    pr_ = mp.mpf(1)
    for i in range(1, n+1): pr_ *= qn(i, r)
    return pr_
def q6j(j1, j2, j3, j4, j5, j6, r):
    tri = [(j1,j2,j3),(j1,j5,j6),(j4,j2,j6),(j4,j5,j3)]
    pref = mp.mpf(1)
    for (a1, b1, c1) in tri: pref *= mp.sqrt(qf(a1+b1-c1, r)*qf(a1-b1+c1, r)*qf(-a1+b1+c1, r)/qf(a1+b1+c1+1, r))
    a = [sum(t_) for t_ in tri]; b = [j1+j2+j4+j5, j2+j3+j5+j6, j3+j1+j6+j4]
    s_ = mp.mpf(0)
    for t in range(max(a), min(b)+1):
        den = mp.mpf(1)
        for ai in a: den *= qf(t-ai, r)
        for bi in b: den *= qf(bi-t, r)
        s_ += (-1)**t*qf(t+1, r)/den
    return pref*s_
def sph_tet(Ls):
    BC, CA, AB, AD, BD, CD = Ls
    dist = {(0,1): AB, (0,2): CA, (0,3): AD, (1,2): BC, (1,3): BD, (2,3): CD}
    G = np.eye(4)
    for (i, k), val in dist.items(): G[i, k] = G[k, i] = math.cos(val)
    detG = np.linalg.det(G)
    if detG <= 0: return None
    Ci = np.linalg.inv(G)
    pairs = [(1,2),(0,2),(0,1),(0,3),(1,3),(2,3)]
    th = []
    for (i, k) in pairs:
        rest = [q for q in range(4) if q not in (i, k)]
        th.append(math.acos(-Ci[rest[0], rest[1]]/math.sqrt(Ci[rest[0], rest[0]]*Ci[rest[1], rest[1]])))
    return th, detG
def sph_volume(Ls, n=300):
    # Schlafli: dV = (1/2) sum L_i d(theta_i), integrated along t*L from tiny to 1
    ts = np.linspace(0.05, 1, n); vol = 0; prev = None
    for t in ts:
        st = sph_tet([t*l for l in Ls])
        if st is None: prev = None; continue
        if prev is not None:
            vol += 0.5*sum(l*(t + tp)/2*(a - b) for l, a, b in zip(Ls, st[0], prev))
        prev, tp = st[0], t
    return vol
E3 = {}
for r in [160, 400]:
    rows = []
    for j6 in range(14, 60):
        js = (30, 31, 32, 31, 30, j6)
        if max(js) > (r-2)/2: continue
        Ls = [2*np.pi*(jj + 0.5)/r for jj in js]
        st = sph_tet(Ls)
        if st is None: continue
        th, detG = st
        vol = sph_volume(Ls)
        regge_sph = sum(Li*(np.pi - t) for Li, t in zip(Ls, th))
        flat = ponzano_regge(js)
        rows.append(dict(j6=j6, exact=float(q6j(*js, r)),
                         with_volume=math.cos(r/(2*np.pi)*(regge_sph + 2*vol) + np.pi/4)*detG**-0.25,
                         no_volume=math.cos(r/(2*np.pi)*regge_sph + np.pi/4)*detG**-0.25,
                         flat_q1=(flat[0] if flat else float("nan")), volume=vol))
    ex = np.array([q["exact"] for q in rows])
    cw = np.corrcoef(ex, [q["with_volume"] for q in rows])[0, 1]
    cn = np.corrcoef(ex, [q["no_volume"] for q in rows])[0, 1]
    fl = np.array([q["flat_q1"] for q in rows]); okf = ~np.isnan(fl)
    cf = np.corrcoef(ex[okf], fl[okf])[0, 1]
    sign = float(np.mean(np.sign(ex) == np.sign([q["with_volume"] for q in rows])))
    E3[r] = dict(rows=rows, corr_with_volume=cw, corr_no_volume=cn, corr_flat=cf, sign_agree=sign)
check("E3 turning on q = e^{i pi/r} curves the tetrahedra: quantum 6j = cos(spherical Regge action + volume term + pi/4) (Turaev-Viro; Mizoguchi-Tada)",
      all(v["corr_with_volume"] > 0.99 and v["sign_agree"] == 1.0 and v["corr_no_volume"] < 0.5 and v["corr_flat"] < 0.8 for v in E3.values()),
      "; ".join(f"r={r}: with volume {v['corr_with_volume']:.4f} (signs {v['sign_agree']:.0%}), no volume {v['corr_no_volume']:.3f}, flat q=1 formula {v['corr_flat']:.3f}" for r, v in E3.items()))
results["E3"] = {str(r): v for r, v in E3.items()}

# ---------------------------------------------------------------- summary
summary = f"{npass} passed, {nfail} failed"
print(summary); out_lines.append(summary)
with open(os.path.join(HERE, "checks_output.txt"), "w") as fh: fh.write("\n".join(out_lines) + "\n")
with open(os.path.join(HERE, "results.json"), "w") as fh: json.dump(results, fh, indent=1, default=float)
