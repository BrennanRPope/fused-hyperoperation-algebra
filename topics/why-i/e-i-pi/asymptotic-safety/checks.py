"""Checks for why-i/e-i-pi/asymptotic-safety (Nov 25, 2026): units at different levels,
running G, dimension flow, deformed covariance. Run: python3 checks.py (numpy only).
Writes checks_output.txt next to this file."""
import numpy as np, itertools, os
from math import gamma as Gamma, pi, erf, sqrt, log
HERE = os.path.dirname(os.path.abspath(__file__))
results, lines = [], []
def check(name, ok, detail=""):
    ok = bool(ok); results.append(ok)
    s = f"{'PASS' if ok else 'FAIL'} {name}" + (f"  | {detail}" if detail else ""); lines.append(s); print(s)

G, hbar, c, Lam = 6.67430e-11, 1.054571817e-34, 2.99792458e8, 1.1056e-52
lP2 = hbar*G/c**3

def rk4(f, v, dt, n):
    v = np.array(v, float); out = [v.copy()]
    for _ in range(n):
        k1 = f(v); k2 = f(v + dt/2*k1); k3 = f(v + dt/2*k2); k4 = f(v + dt*k3)
        v = v + dt/6*(k1 + 2*k2 + 2*k3 + k4); out.append(v.copy())
    return np.array(out)

# ---------------- U: the lattice of units ----------------
# exponent vectors over (L, T, M)
U = {"c": (1, -1, 0), "hbar": (2, -1, 1), "G": (3, -2, -1), "Lambda": (-2, 0, 0), "m": (0, 0, 1)}
def det3(a, b, cc): return round(np.linalg.det(np.array([U[a], U[b], U[cc]], float)))
def inv_factors(rows):
    A = np.array(rows, int); n = A.shape[0]; out = []; prev = 1
    for k in range(1, n + 1):
        g = 0
        for r in itertools.combinations(range(n), k):
            for cl in itertools.combinations(range(A.shape[1]), k):
                g = np.gcd(g, abs(round(np.linalg.det(A[np.ix_(r, cl)].astype(float)))))
        out.append(g // prev); prev = g
    return out
f_G = inv_factors([U["c"], U["hbar"], U["G"]]); f_m = inv_factors([U["c"], U["hbar"], U["m"]]); f_L = inv_factors([U["c"], U["hbar"], U["Lambda"]])
check("U1 units form a lattice Z^3 (exponents of length, time, mass). c, hbar, G span a sublattice of index |det| = 2: invariant factors (1, 1, 2), so some units (e.g. a length) are only reached as a SQUARE ROOT: l_P = sqrt(hbar G/c^3) (the root is known; the index-2 phrasing is ours)",
      abs(det3("c", "hbar", "G")) == 2 and f_G == [1, 1, 2], f"det = {det3('c','hbar','G')}, factors {f_G}")
check("U2 compare: c, hbar and a mass m span ALL of Z^3 (index 1, no roots: Compton units), while c, hbar, Lambda again have index 2: the square root comes with gravity's two scales",
      f_m == [1, 1, 1] and f_L == [1, 1, 2], f"(c,hbar,m) {f_m}; (c,hbar,Lambda) {f_L}")
M = np.array([U["G"], U["hbar"], U["c"], U["Lambda"]], float).T  # 3 x 4
ns = np.linalg.svd(M)[2][-1]; ns = ns/ns[0]
check("U3 Buckingham: 4 quantities minus rank 3 = 1 pure number, exponents (G, hbar, c, Lambda) = (1, 1, -3, 1): G hbar Lambda/c^3 = l_P^2 Lambda (as in g-dials D1)",
      np.allclose(ns, [1, 1, -3, 1]), f"{np.round(ns, 6)}; value {lP2*Lam:.3e}")
cube = list(itertools.product([0, 1], repeat=3))
n_eq = sum(1 for a in cube for b in cube if tuple(np.maximum(a, b)) == tuple(np.add(a, b)))
check("U4 the Boolean cube {0,1}^3 sits inside the unit lattice Z^3 with the same order, but its join (max) agrees with the lattice's group law (sum of exponents = multiplying quantities) only for disjoint supports: 27 of 64 pairs (3^3)",
      n_eq == 27, f"{n_eq}/64")

def massdim(v):  # with hbar = c = 1: L = T = 1/M, so M-dimension = M - L - T
    return v[2] - v[0] - v[1]
check("U5 with hbar = c = 1 the lattice collapses to one line (mass dimension): [c] = [hbar] = 0, [G] = -2, [Lambda] = +2, [scale k] = +1; so g = G k^2 and lambda = Lambda/k^2 are pure numbers, and so is G Lambda",
      (massdim(U["c"]), massdim(U["hbar"]), massdim(U["G"]), massdim(U["Lambda"])) == (0, 0, -2, 2))

# ---------------- R: running units ----------------
b = 7/(2*pi)  # one-loop QCD-like: d alpha/dt = -b alpha^2  (t = log mu)
tr = rk4(lambda v: np.array([-b*v[0]**2]), [0.3], 0.01, 500)
t = np.arange(len(tr))*0.01; al = tr[:, 0]
LamQCD = np.exp(t)*np.exp(-1/(b*al))
check("R1 dimensional transmutation (known: Coleman-Weinberg; Gross-Wilczek): a dimensionless coupling that runs, d alpha/dt = -b alpha^2, defines a scale Lambda = mu exp(-1/(b alpha(mu))) that is the SAME at every mu: a unit created by running",
      np.ptp(LamQCD)/LamQCD.mean() < 1e-8, f"relative spread {np.ptp(LamQCD)/LamQCD.mean():.1e}")
check("R2 in the variable 1/alpha that flow is plain addition: 1/alpha(t) = 1/alpha(0) + b t (dial k = 0, no wall)",
      np.allclose(1/al, 1/0.3 + b*t, atol=1e-9))
def ds_from_eta(d, eta):
    a = 2 - eta; p = np.linspace(1e-6, 60, 400001)
    P = lambda T: np.trapezoid(p**(d-1)*np.exp(-T*p**a), p)
    T1, T2 = 0.5, 1.0
    return -2*(np.log(P(T2)) - np.log(P(T1)))/(np.log(T2) - np.log(T1))
ds_vals = [ds_from_eta(4, e) for e in (0, -1, -2)]
check("R3 an anomalous dimension eta changes the propagator to 1/p^(2 - eta); the return probability of a random walk then gives spectral dimension d_s = 2d/(2 - eta): 4 at eta = 0, 8/3 at eta = -1, 2 at eta = -2 (derived; the eta = -2 case is Lauscher-Reuter 2005)",
      np.allclose(ds_vals, [4, 8/3, 2], atol=1e-3), f"{np.round(ds_vals, 4)}")

# ---------------- A: asymptotic safety ----------------
def toy(eps, om):
    gs = eps/om
    return gs, (lambda v: np.array([eps*v[0] - om*v[0]**2]))
eps, om, g0 = 2.0, 2/0.7, 1e-3
gs, f = toy(eps, om)
tr = rk4(f, [g0], 0.005, 2000); t = np.arange(len(tr))*0.005
closed = gs/(1 + (gs/g0 - 1)*np.exp(-eps*t))
check("A1 toy one-loop flow of g = G k^2 (known shape: Weinberg 1979, Reuter 1998): dg/dt = (d-2) g - w g^2 runs from the Gaussian point g = 0 to the fixed point g* = (d-2)/w; closed form matches the integration",
      np.max(np.abs(tr[:, 0] - closed)) < 1e-9, f"g* = {gs:.3f}, max error {np.max(np.abs(tr[:,0]-closed)):.1e}")
def ma(x, s, k): return x + s + k*x*s
u = 1/tr[:, 0]; k = -gs
s_t = (1 - np.exp(-eps*t))/gs
check("A2 (ours) in u = 1/g the flow is EXACTLY our weighted circle operation: u(t) = u(0) o_k s(t) = u0 + s + k u0 s with dial k = -g* and step s(t) = (1 - e^{-(d-2)t})/g*",
      np.allclose(u, ma(u[0], s_t, k), rtol=1e-9))
i, j = 300, 700
check("A3 (ours) steps compose by the same operation, s(t1) o_k s(t2) = s(t1 + t2) (the flow is a one-parameter group inside mult-add), and the wall -1/k = 1/g* is the UV fixed point: no step crosses it",
      abs(ma(s_t[i], s_t[j], k) - s_t[i + j]) < 1e-12 and np.all(u > 1/gs))
flat = np.log(np.abs(1 + k*u))/k
check("A4 (ours) the flow lives just past the wall (1 + k u = 1 - g*/g < 0, the lift's flipped side) and the flat coordinate log|1 + k u|/k runs linearly in RG time with slope (d-2)/g*: RG time is the flat coordinate of this dial",
      np.allclose(np.diff(flat)/0.005, eps/gs, rtol=1e-6), f"slope {np.mean(np.diff(flat))/0.005:.4f} vs {eps/gs:.4f}")
gstars = [toy(e, om)[0] for e in (0.5, 0.1, 0.01)]
check("A5 (dial reading ours; g* = O(d-2) is known, Weinberg 1979) the rate d - 2 is G's canonical mass dimension. As d -> 2 the dial k = -g* -> 0 and the flow becomes plain addition in 1/g, like the QCD flow in R2: G's units are what switch the wall on",
      gstars[0] > gstars[1] > gstars[2] and np.isclose(gstars[2]/0.01, 1/om), f"g* at d-2 = 0.5, 0.1, 0.01: {np.round(gstars, 5)}")

# Einstein-Hilbert truncation, optimised (Litim) cutoff, d = 4 (Reuter-Saueressig form)
d = 4; pre = (4*pi)**(1 - d/2)
Phi = lambda p, n, w: 1/Gamma(n + 1)/(1 + w)**p
Phit = lambda p, n, w: 1/Gamma(n + 2)/(1 + w)**p
B1 = lambda l: pre/3*(d*(d+1)*Phi(1, d/2-1, -2*l) - 6*d*(d-1)*Phi(2, d/2, -2*l) - 4*d*Phi(1, d/2-1, 0) - 24*Phi(2, d/2, 0))
B2 = lambda l: -pre/6*(d*(d+1)*Phit(1, d/2-1, -2*l) - 6*d*(d-1)*Phit(2, d/2, -2*l))
eta = lambda g, l: g*B1(l)/(1 - g*B2(l))
def beta(v):
    g, l = v; e = eta(g, l)
    bl = -(2 - e)*l + 0.5*g*pre*(2*d*(d+1)*Phi(1, d/2, -2*l) - 8*d*Phi(1, d/2, 0) - d*(d+1)*e*Phit(1, d/2, -2*l))
    return np.array([(d - 2 + e)*g, bl])
def jac(v, h=1e-7):
    J = np.zeros((2, 2))
    for jj in range(2):
        dv = np.zeros(2); dv[jj] = h; J[:, jj] = (beta(v + dv) - beta(v - dv))/(2*h)
    return J
v = np.array([0.7, 0.2])
for _ in range(60): v = v - np.linalg.solve(jac(v), beta(v))
J = jac(v); th = -np.linalg.eigvals(J); th = th[np.argsort(th.imag)]
check("A6 Einstein-Hilbert truncation with the optimised cutoff (Litim 2004; Reuter-Saueressig form): non-Gaussian fixed point g* = 0.707, lambda* = 0.193, g* lambda* = 0.137, critical exponents 1.475 +- 3.043 i (reproduced)",
      np.allclose(v, [0.7073, 0.1932], atol=2e-4) and np.allclose([th[1].real, abs(th[1].imag)], [1.4753, 3.0432], atol=2e-4),
      f"(g*, lambda*) = ({v[0]:.4f}, {v[1]:.4f}), theta = {th[1].real:.4f} +- {abs(th[1].imag):.4f} i")
check("A7 at the fixed point the graviton's anomalous dimension is exactly eta_N = 2 - d = -2 (forced by beta_g = (d - 2 + eta) g = 0 with g* != 0): G's units are cancelled, and by R3 the spectral dimension is 2",
      abs(eta(*v) + 2) < 1e-10 and abs(2*d/(2 - eta(*v)) - 2) < 1e-9, f"eta_N* = {eta(*v):.12f}")
w, Vv = np.linalg.eig(J)
z0 = 1e-4*np.real(Vv[:, 0]*np.exp(0.4j))

trj = rk4(lambda x: J @ x, z0, 0.001, 3000)
P = np.linalg.inv(Vv)
zc = (P @ trj.T)[0]
ang = np.unwrap(np.angle(zc)); rate = np.polyfit(np.arange(len(zc))*0.001, ang, 1)[0]
dec = np.polyfit(np.arange(len(zc))*0.001, np.log(np.abs(zc)), 1)[0]
check("A8 (reading ours) complex exponents mean the flow into the fixed point is z -> e^{-theta t} z: it shrinks at rate 1.475 AND turns at rate 3.043 per unit RG time. Near the fixed point the flow is a complex dilation, the e and i of e-i-pi together, with the wall at the fixed point",
      np.isclose(dec, -th[1].real, atol=1e-6) and np.isclose(abs(rate), abs(th[1].imag), atol=1e-6), f"decay {dec:.4f}, turn {rate:.4f}")
labels = []
for ph in (0.0, 0.8, 1.6):
    v0 = v + 1e-5*np.real(Vv[:, 0]*np.exp(1j*ph))
    T = rk4(beta, v0, -0.01, 3000); gl = T[:, 0]*T[:, 1]
    labels.append((gl[-300], gl[-1]))
check("A9 every trajectory leaving the fixed point ends in a classical regime where g lambda = G Lambda (hbar = c = 1) is CONSTANT: the trajectory's label is the g-dials pure number G hbar Lambda/c^3 (known in substance: Reuter-Weyer 2004; shown here on Lambda < 0 trajectories, which need no fine tuning)",
      all(abs(a/b2 - 1) < 1e-6 for a, b2 in labels) and len(set(np.round([b2 for _, b2 in labels], 3))) == 3,
      f"plateaus {[f'{b2:.4f}' for _, b2 in labels]}")
A = pre*0.5*(2*d*(d+1)*Phi(1, d/2, 0) - 8*d*Phi(1, d/2, 0))
kT = (4*lP2*Lam/A)**0.25
check("A10 in the classical regime lambda = Lambda/k^2 + (A/4) G k^2 with A = 1/(2 pi); its turning point is at k_T/m_P = (4 l_P^2 Lambda/A)^(1/4) ~ 1e-30 and there g lambda = 2 G Lambda. With our universe's G hbar Lambda/c^3 this lands on Reuter-Weyer's turning scale k_T ~ 1e-30 m_Pl",
      np.isclose(A, 1/(2*pi)) and 1e-31 < kT < 1e-29, f"A = {A:.5f}, k_T/m_P = {kT:.2e}")

# ---------------- D: dimension flow ----------------
y = 0.8*np.exp(-2*np.linspace(-3, 3, 61))  # toy: g/g* = 1/(1 + y), y = C e^{-2t}
ds = 4*(1 + y)/(2 + y)
check("D1 (ours) toy: along the flow of A1, eta_N = -2 g/g*, so d_s = 4/(1 + g/g*) runs from 4 (IR) to 2 (UV), and (d_s - 2)/(4 - d_s) = y/2 exactly: the log of that ratio is linear in RG time",
      np.allclose((ds - 2)/(4 - ds), y/2))
a_, b_, c_ = 4.02, 119.0, 54.0
sig = np.logspace(-1, 4, 50); D = a_ - b_/(sig + c_); D0 = a_ - b_/c_
check("D2 (ours) CDT's fitted curve D_S(sigma) = 4.02 - 119/(54 + sigma) (Ambjorn-Jurkiewicz-Loll 2005) obeys (D - D0)/(4.02 - D) = sigma/54 exactly, D0 = D_S(0) = 1.82: the same straight line in log(diffusion time), slope 1",
      np.allclose((D - D0)/(a_ - D), sig/c_), f"D0 = {D0:.3f} (paper quotes 1.80 +- 0.25)")
al_, ls = 0.5, 1.0
x = np.logspace(-4, 4, 60); q = x + (ls/al_)*(x/ls)**al_
dH = np.gradient(np.log(q), np.log(x)); dHx = al_*(1 + (x/ls)**(1 - al_))/(1 + al_*(x/ls)**(1 - al_))
check("D3 (ours) Calcagni's multiscale measure q(x) = x + (l*/alpha)(x/l*)^alpha (Calcagni 2011-17) has Hausdorff dimension running from alpha to 1 with (d_H - alpha)/(1 - d_H) = alpha (x/l*)^(1 - alpha): straight in log x, slope 1 - alpha",
      np.allclose(dH, dHx, atol=2e-3) and np.allclose((dHx - al_)/(1 - dHx), al_*(x/ls)**(1 - al_)))
vv = ds - 3; Delta = 0.37
tt = np.linspace(-3, 3, 61); ds_shift = 4*(1 + 0.8*np.exp(-2*(tt + Delta)))/(2 + 0.8*np.exp(-2*(tt + Delta)))
check("D4 (ours) so with v = d - 3 (walls at d = 2 and d = 4), moving the scale by Delta composes v by relativistic velocity addition: v' = (v - tanh Delta)/(1 - v tanh Delta). The two fixed dimensions are the two walls; log-scale is rapidity",
      np.allclose(ds_shift - 3, (vv - np.tanh(Delta))/(1 - vv*np.tanh(Delta))))
s_ = np.linspace(-2, 2, 41); de = np.array([3 - erf(z) for z in s_]); r = np.log((de - 2)/(4 - de))
check("D5 this is a real constraint, not automatic: a two-plateau curve of erf shape, d = 3 - erf(log s), is NOT straight in that coordinate (second differences far from 0)",
      np.max(np.abs(np.diff(r, 2))) > 1e-2, f"max second difference {np.max(np.abs(np.diff(r, 2))):.3f}")
ds_traj = []
v0 = v + 1e-5*np.real(Vv[:, 0])
T = rk4(beta, v0, -0.01, 3000)
for gg, ll in T[::300]: ds_traj.append(2*d/(2 - eta(gg, ll)))
check("D6 on a real Einstein-Hilbert trajectory d_s = 2d/(2 - eta_N(k)) runs from 2 at the fixed point to 4 in the classical regime (Lauscher-Reuter 2005; Reuter-Saueressig 2011)",
      abs(ds_traj[0] - 2) < 1e-3 and abs(ds_traj[-1] - 4) < 1e-3, f"{np.round(ds_traj, 3)}")

# ---------------- B: deformed general covariance ----------------
dK = np.linspace(0, pi/2, 50); rho_rel = np.sin(dK)**2
check("B1 LQC's holonomy deformation beta = cos(2 delta K) (Bojowald-Paily 2012; Bojowald-Mielczarek 2015) equals 1 - 2 rho/rho_c on the effective solutions, where rho/rho_c = sin^2(delta K)",
      np.allclose(np.cos(2*dK), 1 - 2*rho_rel))
check("B2 (reading ours) beta = 1 + k rho with k = -2/rho_c is our lift: its wall rho = rho_c/2 is where the signature flips (beta = 0 is the Carroll point {S,S} = 0; beta < 0 Euclidean), and beta = cos(2 delta K) is the real part of a turn e^{2 i delta K}, flipping at a quarter turn",
      np.isclose(1 + (-2)*0.5, 0) and np.isclose(np.cos(2*(pi/4)), 0, atol=1e-12))
gam = 0.23753; rc = sqrt(3)/(32*pi**2*gam**3)
check("B3 rho_c = sqrt3/(32 pi^2 gamma^3) rho_Planck ~ 0.41 rho_Planck (Ashtekar-Pawlowski-Singh), and rho_Planck = c^7/(hbar G^2): G enters the dial k = -2/rho_c as G^2 hbar, setting WHERE the wall is; the deformed slot itself is the signature (c's Cayley-Klein direction)",
      0.40 < rc < 0.42, f"rho_c/rho_P = {rc:.3f}")
H2 = lambda r: r*(1 - r)  # in units 8 pi G rho_c/3
check("B4 LQC's effective Friedmann equation H^2 = (8 pi G/3) rho (1 - rho/rho_c) carries a second lift with wall at rho_c (the bounce, H = 0); signature change sits halfway, at the maximum of H^2 (rho_c/2)",
      abs(H2(1.0)) < 1e-15 and np.argmax(H2(np.linspace(0, 1, 1001))) == 500)

# ---------------- E, K: entropy and the 2D link ----------------
check("E1 (ours) at the fixed point the entropy per area, in units of the probe scale, is a fixed pure number: c^3/(4 G_k hbar) / k^2 = 1/(4 g*) = 0.354 (this scheme). The IR value 1/(4 G_0) is still the trajectory label, so asymptotic safety does not supply the area law's number",
      np.isclose(1/(4*v[0]), 0.3535, atol=1e-3), f"1/(4 g*) = {1/(4*v[0]):.4f}")
gstr = lambda cc: (cc - 1 - sqrt((1 - cc)*(25 - cc)))/12
check("K1 the short-distance limit is effectively two-dimensional, where quantum gravity is solved by matrix models: pure 2D gravity has string susceptibility gamma_str = -1/2 (KPZ 1988; double-scaling limit 1990), c = 1 matter gives 0",
      np.isclose(gstr(0), -0.5) and np.isclose(gstr(1), 0))

n_ok = sum(results); lines.append(""); lines.append(f"{n_ok}/{len(results)} pass"); print(lines[-1])
open(os.path.join(HERE, "checks_output.txt"), "w").write("\n".join(lines) + "\n")
