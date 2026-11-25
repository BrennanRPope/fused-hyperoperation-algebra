"""Path integrals in a stronger sense: what the three missing pieces need.  topics/why-i/path-integral-strong (Nov 25, 2026).
Each check tests an identity or a numerical fact that the write-up leans on. Published theorems are cited, not re-proved.
Run: python3 checks.py   (numpy, sympy, scipy)
"""
import numpy as np, sympy as sp
from scipy.linalg import expm, eigh
rng = np.random.default_rng(1125)
results = []
def check(name, ok, detail=""):
    results.append((name, bool(ok))); print(("PASS " if ok else "FAIL ") + name + (("  | " + detail) if detail else ""))

# ===================== GAP 1: THE DYNAMICS =====================
# K1 Bargmann: boost and shift commute up to the phase e^{i m v a}: the order cost of boost-then-shift is the mass
x, a, v, m, p = sp.symbols('x a v m p', real=True)
psi = sp.Function('psi')
shift = lambda f: f.subs(x, x + a)
boost = lambda f: sp.exp(sp.I*m*v*x)*f
lhs = shift(boost(psi(x))); rhs = boost(shift(psi(x)))
check("K1 shift-after-boost = e^(i m v a) x boost-after-shift: the boost/shift twist is set by the mass (Bargmann 1954)",
      sp.simplify(lhs/rhs - sp.exp(sp.I*m*v*a)) == 0, str(sp.simplify(lhs/rhs)))
# K2: translation invariance => H = f(P); boost covariance [K, H] = i P/m*m ... with K = m x: [m x, f(P)] = i m f'(P); demand = i P
f = sp.Function('f')
sol = sp.dsolve(sp.Eq(m*sp.diff(f(p), p), p), f(p))
check("K2 shift symmetry + boost rule [K, H] = iP (with [K, P] = i m) force H = P^2/(2m) + constant (Levy-Leblond)",
      sp.simplify(sol.rhs - p**2/(2*m)).free_symbols <= {sp.Symbol('C1')}, str(sol))
# K3: check [m x, P^2/2m] = i P on a test function (P = -i d/dx, hbar = 1)
g = sp.exp(-x**2)*sp.sin(3*x)
P = lambda h: -sp.I*sp.diff(h, x)
comm = m*x*P(P(g))/(2*m) - P(P(m*x*g))/(2*m)
check("K3 operator check: [m x, P^2/2m] = i P", sp.simplify(comm - sp.I*P(g)) == 0)

# Q: the quadratic world. Poisson brackets of monomials: {q^a p^b, q^c p^d} = (a d - b c) q^(a+c-1) p^(b+d-1)
q_, p_ = sp.symbols('q p')
def pb(F, G): return sp.expand(sp.diff(F, q_)*sp.diff(G, p_) - sp.diff(F, p_)*sp.diff(G, q_))
ok = True
for A_, B_, C_, D_ in [(2, 0, 0, 2), (1, 1, 2, 0), (3, 1, 0, 2), (1, 2, 2, 1), (0, 3, 3, 0)]:
    ok &= sp.simplify(pb(q_**A_*p_**B_, q_**C_*p_**D_) - (A_*D_ - B_*C_)*q_**(A_+C_-1)*p_**(B_+D_-1)) == 0
check("Q1 the Poisson bracket of two monomials has coefficient a*d - b*c: the order dial's wedge, on exponents", ok)
quad = [sp.Integer(1), q_, p_, q_**2, p_**2, q_*p_]
deg = lambda e: sp.Poly(e, q_, p_).total_degree() if e != 0 else -1
check("Q2 brackets of polynomials of degree <= 2 stay of degree <= 2 (a finite-dimensional algebra: Heisenberg + sp(2))",
      all(deg(pb(F, G)) <= 2 for F in quad for G in quad))
check("Q3 cubics do not close: {q^3, p^3} = 9 q^2 p^2 has degree 4 (why Groenewold-van Hove stops at quadratics)",
      pb(q_**3, p_**3) == 9*q_**2*p_**2)
check("Q4 the order dial's affine pair sits inside: {q, q p} = q, the affine rule [X, Y] = Y (Klauder's Q and D)",
      pb(q_, q_*p_) == q_)
# Q5: the oscillator path integral as a finite composition: kick-drift steps in Sp(2) -> exact flow
w, T = 1.3, 2.0
exact = np.array([[np.cos(w*T), np.sin(w*T)/w], [-w*np.sin(w*T), np.cos(w*T)]])
errs = []
for N in [10, 100, 1000]:
    dt = T/N
    drift = np.array([[1, dt], [0, 1.]]); kick = np.array([[1, 0], [-w*w*dt, 1.]])
    half = np.array([[1, 0], [-w*w*dt/2, 1.]])
    step = half @ drift @ half
    errs.append(np.linalg.norm(np.linalg.matrix_power(step, N) - exact))
check("Q5 for a quadratic Hamiltonian the time-sliced chain is a product of 2x2 symplectic matrices converging to the exact flow (error ~ 1/N^2)",
      errs[0]/errs[1] > 50 and errs[1]/errs[2] > 50, f"{[f'{e:.1e}' for e in errs]}")
check("Q6 each step keeps oriented phase-space area (det = 1): the chain never loses the twist",
      np.isclose(np.linalg.det(step), 1))

# ===================== GAP 2: WHAT FIXES hbar =====================
# B1 Bohr-Sommerfeld with the Maslov half: area(E_n) = 2 pi hbar (n + 1/2) gives the exact oscillator levels
hbar, mass = 1.0, 1.0
L, Ng = 24.0, 1200
xs = np.linspace(-L/2, L/2, Ng); dx = xs[1] - xs[0]
lap = (np.diag(np.ones(Ng-1), 1) - 2*np.eye(Ng) + np.diag(np.ones(Ng-1), -1))/dx**2
H = -hbar**2/(2*mass)*lap + np.diag(0.5*mass*w**2*xs**2)
E = eigh(H, eigvals_only=True)[:8]
bs = hbar*w*(np.arange(8) + 0.5)
check("B1 Bohr-Sommerfeld: orbit area 2*pi*E/w set to 2*pi*hbar*(n + 1/2) reproduces the computed oscillator levels",
      np.allclose(E, bs, rtol=2e-3), f"computed {np.round(E[:4], 4)}, area rule {bs[:4]}")
# B2 anharmonic: area rule is only approximate (an honest limit)
H4 = -0.5*lap + np.diag(xs**4)
E4 = eigh(H4, eigvals_only=True)[:6]
from scipy.integrate import quad as nquad
from scipy.optimize import brentq
def area4(En):
    xt = En**0.25
    return 2*nquad(lambda X: np.sqrt(max(0, 2*(En - X**4))), -xt, xt)[0]
bs4 = [brentq(lambda En: area4(En) - 2*np.pi*(n + .5), 1e-3, 50) for n in range(6)]
rel = np.abs(np.array(bs4) - E4)/E4
check("B2 for x^4 the same area rule is off by 18% at n = 0 and improves fast with n (it is the leading semiclassical term)",
      rel[0] < 0.2 and rel[5] < rel[0] and rel[5] < 0.01, f"relative errors {np.round(rel, 4)}")
# B3 Kostant-Souriau on the sphere: area 4*pi*s must be 2*pi*hbar*n, so s = n*hbar/2; with the shifted radius
# s = hbar*(j + 1/2) the 2j + 1 states fill the sphere one area-cell 2*pi*hbar each
ok = all(np.isclose(4*np.pi*hbar*(j + .5)/(2*np.pi*hbar), 2*j + 1) for j in [0, .5, 1, 1.5, 2])
check("B3 spin: a sphere of area 4*pi*hbar*(j + 1/2) holds exactly 2j + 1 cells of area 2*pi*hbar (Kostant-Souriau integrality)", ok)
# B4 a closed chain in the Heisenberg group composes to its enclosed area; requiring the identity forces area in 2*pi*hbar*Z
def heis_loop(vs_):
    c = 0.; tot = np.zeros(2)
    for v_ in vs_:
        c += 0.5*(tot[0]*v_[1] - tot[1]*v_[0]); tot += v_
    return tot, c
sq = np.array([[2., 0], [0, np.pi], [-2, 0], [0, -np.pi]])
tot, c = heis_loop(sq)
check("B4 a closed loop of steps returns to the start shifted only by its area (here 2*pi); with hbar = 1 its phase e^(i*area) = 1 exactly",
      np.allclose(tot, 0) and np.isclose(c, 2*np.pi) and np.isclose(np.exp(1j*c/hbar), 1))
check("B5 the same loop with hbar = 0.7 gives phase != 1: the rule fixes area in units of hbar, not hbar itself",
      not np.isclose(np.exp(1j*c/0.7), 1), f"phase {np.exp(1j*c/0.7):.3f}")

# ===================== GAP 3: THE MEASURE =====================
# M1 Feynman-Kac (real dial): Brownian paths weighted by exp(-int V) reproduce e^{-tH} f
Vf = lambda X: 0.5*X**2
t_, x0 = 1.0, 0.4
fgrid = np.exp(-(xs - 1.0)**2)
Hfk = -0.5*lap + np.diag(Vf(xs))
u_exact = np.interp(x0, xs, expm(-t_*Hfk) @ fgrid)
Npaths, Nst = 40000, 200
dt = t_/Nst
B = x0 + np.cumsum(rng.normal(scale=np.sqrt(dt), size=(Npaths, Nst)), axis=1)
weights = np.exp(-dt*(0.5*Vf(x0) + Vf(B[:, :-1]).sum(axis=1) + 0.5*Vf(B[:, -1])))
u_mc = np.mean(weights*np.exp(-(B[:, -1] - 1.0)**2))
se = np.std(weights*np.exp(-(B[:, -1] - 1.0)**2))/np.sqrt(Npaths)
check("M1 Feynman-Kac: averaging Brownian paths with weight exp(-int V dt) matches e^(-tH) f (rigorous, Kac 1949)",
      abs(u_mc - u_exact) < 4*se + 2e-3, f"paths {u_mc:.4f} +/- {se:.4f}, exact {u_exact:.4f}")
# M2 the affine coherent-state label space with Klauder's metric is our half-plane with the Poincare metric
qq, pp = sp.symbols('q p', positive=True)
A_, B_ = 1/qq, pp                                      # group coords: x -> A x + B  (U[p,q] composition (q,p)(q',p') = (q q', p + p'/q))
J = sp.Matrix([[sp.diff(B_, pp), sp.diff(B_, qq)], [sp.diff(A_, pp), sp.diff(A_, qq)]])
poinc = sp.diag(1/A_**2, 1/A_**2)
pulled = sp.simplify(J.T*poinc*J)
check("M2 Klauder's affine phase-space metric q^2 dp^2 + dq^2/q^2 is exactly the Poincare length on the order dial's half-plane (a = 1/q, b = p)",
      sp.simplify(pulled - sp.diag(qq**2, 1/qq**2)) == sp.zeros(2), str(pulled))
# composition rule check
def Uqp(q1, p1, q2, p2): return (q1*q2, p1 + p2/q1)
ok = True
for _ in range(20):
    q1, q2 = rng.uniform(.2, 3, 2); p1, p2 = rng.normal(size=2)
    qn, pn = Uqp(q1, p1, q2, p2)
    ok &= np.isclose(1/qn, (1/q1)*(1/q2)) and np.isclose(pn, (1/q1)*p2 + p1)
check("M3 Klauder's labels compose like the order dial's maps: (a, b) = (1/q, p) follows x -> a x + b", ok)
# M4 Brownian motion on the half-plane: homogeneous (same law from every start, after moving by the group)
def hyp_bm(z0, T_, n_, N_):
    xx = np.full(N_, z0.real); yy = np.full(N_, z0.imag); dt_ = T_/n_
    for _ in range(n_):
        dW = rng.normal(scale=np.sqrt(dt_), size=(2, N_))
        xx = xx + yy*dW[0]
        yy = yy*np.exp(dW[1] - 0.5*dt_)                  # exact log-step for dy = y dW
    return xx + 1j*yy
def hdist(z1, z2): return np.arccosh(1 + np.abs(z1 - z2)**2/(2*z1.imag*z2.imag))
Z1 = hyp_bm(1j, 1.0, 400, 20000); Z2 = hyp_bm(3 + 0.25j, 1.0, 400, 20000)
d1 = hdist(Z1, 1j); d2 = hdist(Z2, 3 + 0.25j)
qs = [.25, .5, .75, .9]
q1_, q2_ = np.quantile(d1, qs), np.quantile(d2, qs)
check("M4 hyperbolic Brownian motion looks the same from every start: distance-from-start quantiles agree (start i vs start 3 + 0.25i)",
      np.allclose(q1_, q2_, rtol=0.04), f"{np.round(q1_, 3)} vs {np.round(q2_, 3)}")
check("M5 its stretch coordinate drifts like E[log lam_t] = -t/2 (generator (1/2) lam^2 Laplacian)",
      abs(np.mean(np.log(Z1.imag)) + 0.5) < 0.02, f"{np.mean(np.log(Z1.imag)):.4f}")

# M6 curvature of Klauder's metric with the scale put back: (beta*hbar)^-1 q^2 dp^2 + beta*hbar q^-2 dq^2
beta, hb = sp.symbols('beta hbar', positive=True)
E_, G_ = (beta*hb)**-1*qq**2, beta*hb/qq**2            # metric E dp^2 + G dq^2, coordinates (p, q)
sE, sG = sp.sqrt(E_), sp.sqrt(G_)
Kgauss = sp.simplify(-1/(sE*sG)*(sp.diff(sp.diff(sE, qq)/sG, qq) + sp.diff(sp.diff(sG, pp)/sE, pp)))
check("M6 with its scale put back, Klauder's metric has Gaussian curvature -1/(beta*hbar) (scalar curvature -2/(beta*hbar), the figure Klauder quotes)",
      sp.simplify(Kgauss + 1/(beta*hb)) == 0, f"K = {Kgauss}")
# B6 its area element is dp dq (the symplectic area), so on a closed surface of genus g with that curvature,
# Gauss-Bonnet gives area 2*pi*(2g-2)*beta*hbar and the number of hbar-cells is (2g-2)*beta: hbar cancels
area_el = sp.simplify(sp.sqrt(E_*G_))
g_ = sp.symbols('g', positive=True, integer=True)
Area = 2*sp.pi*(2*g_ - 2)/(1/(beta*hb))
cells = sp.simplify(Area/(2*sp.pi*hb))
check("B6 [ours, arithmetic] Klauder's area element is dp dq, so a closed genus-g surface of that curvature holds (2g-2)*beta cells of size 2*pi*hbar: integrality fixes beta, not hbar",
      area_el == 1 and sp.simplify(cells - (2*g_ - 2)*beta) == 0, f"area element {area_el}, cells {cells}")

# ===================== SHAPES ELSEWHERE =====================
# S1 Duistermaat-Heckman on the sphere: int e^{i t z} dA equals its stationary-phase answer (two poles) exactly
from scipy.integrate import quad as q1d
ok = True
for tt in [0.3, 1.7, 5.0, 12.0]:
    exact_val = 2*np.pi*q1d(lambda z: np.cos(tt*z), -1, 1)[0] + 2j*np.pi*q1d(lambda z: np.sin(tt*z), -1, 1)[0]
    stat = 2*np.pi*(np.exp(1j*tt) - np.exp(-1j*tt))/(1j*tt)     # north pole + south pole terms
    ok &= np.isclose(exact_val, stat)
check("S1 Duistermaat-Heckman: on the sphere, int exp(i t z) dA equals the sum of its two stationary-phase terms exactly, for every t", ok)
# S2 Feynman checkerboard vs Kac telegraph: the same two-direction chain, turning amplitude imaginary vs real
pmom, mm, tt = 1.1, 0.8, 1.5
sz = np.diag([1., -1.]); sx = np.array([[0, 1.], [1, 0]])
def chain(turn, N):
    eps = tt/N
    move = np.diag([np.exp(-1j*pmom*eps), np.exp(1j*pmom*eps)])   # right movers / left movers in momentum space
    flip = np.eye(2) + turn*eps*sx
    return np.linalg.matrix_power(move @ flip, N)
errD = abs(chain(1j*mm, 4000) - expm(tt*(-1j*pmom*sz + 1j*mm*sx))).max()
errT = abs(chain(mm, 4000) - expm(tt*(-1j*pmom*sz + mm*sx))).max()
check("S2 checkerboard/telegraph: one two-way chain; turn weight i*m*eps gives the 1+1 Dirac propagator, real turn rate gives Kac's telegraph kernel",
      errD < 1e-3 and errT < 1e-3, f"errors {errD:.1e}, {errT:.1e}")
# S3 Suzuki-Trotter: a quantum spin's partition function = a classical Ising ring one dimension up
hT, gT, betaT, NT = 0.7, 0.4, 2.0, 12
eps = betaT/NT
Kc = -0.5*np.log(np.tanh(eps*hT)); Ac = np.sqrt(np.sinh(eps*hT)*np.cosh(eps*hT))
import itertools
Zc = 0.
for s in itertools.product([1, -1], repeat=NT):
    s = np.array(s)
    Zc += Ac**NT*np.exp(Kc*np.sum(s*np.roll(s, -1)) + eps*gT*np.sum(s))
Zq_trot = np.trace(np.linalg.matrix_power(expm(eps*hT*sx) @ expm(eps*gT*sz), NT))
Zq = np.trace(expm(betaT*(hT*sx + gT*sz)))
check("S3 Suzuki-Trotter: the time-sliced quantum trace equals a classical Ising ring sum exactly, and tends to Tr e^(beta H)",
      np.isclose(Zc, Zq_trot.real) and abs(Zq_trot - Zq)/Zq < 0.01, f"ring {Zc:.5f}, sliced {Zq_trot.real:.5f}, exact {Zq:.5f}")
# S4 lattice gauge: the product of link phases around a big loop = product of plaquette phases (discrete Stokes)
Lg = 5
th_x = rng.uniform(-np.pi, np.pi, (Lg + 1, Lg + 1)); th_y = rng.uniform(-np.pi, np.pi, (Lg + 1, Lg + 1))
loop = 0.
for i in range(Lg): loop += th_x[i, 0]
for j in range(Lg): loop += th_y[Lg, j]
for i in range(Lg): loop -= th_x[i, Lg]
for j in range(Lg): loop -= th_y[0, j]
plaq = sum(th_x[i, j] + th_y[i + 1, j] - th_x[i, j + 1] - th_y[i, j] for i in range(Lg) for j in range(Lg))
check("S4 Wilson loop: the phase around a loop equals the product of the plaquette phases inside (discrete Stokes)", np.isclose(np.exp(1j*loop), np.exp(1j*plaq)))
# S5 a closed geodesic on a hyperbolic surface: element [[2,1],[1,1]] moves points on its axis by 2*arccosh(tr/2)
Mg = np.array([[2., 1], [1, 1]])
ev = np.linalg.eigvals(Mg)
fp = [(Mg[0, 0] - Mg[1, 1] + s*np.sqrt((Mg[0, 0] + Mg[1, 1])**2 - 4))/(2*Mg[1, 0]) for s in (1, -1)]
cen, rad = (fp[0] + fp[1])/2, abs(fp[0] - fp[1])/2
z0 = cen + rad*np.exp(1j*1.1)
gz = (Mg[0, 0]*z0 + Mg[0, 1])/(Mg[1, 0]*z0 + Mg[1, 1])
check("S5 Selberg/Gutzwiller shape: a hyperbolic element slides points along its axis by length 2*arccosh(|tr|/2) (here tr = 3)",
      np.isclose(hdist(z0, gz), 2*np.arccosh(3/2)) and np.isclose(abs(gz - cen), rad), f"{hdist(z0, gz):.5f} vs {2*np.arccosh(1.5):.5f}")
# S6 Lefschetz thimble: the oscillating Fresnel integral becomes a decaying Gaussian on the turned line
from scipy.special import fresnel
def fres_R(R):                                   # int_{-R}^{R} exp(i x^2) dx, via Fresnel integrals
    S_, C_ = fresnel(R*np.sqrt(2/np.pi)); return 2*np.sqrt(np.pi/2)*(C_ + 1j*S_)
fres = fres_R(400.)                              # the real-line integral only converges slowly (error ~ 1/R)
rot = np.exp(1j*np.pi/4)*q1d(lambda s: np.exp(-s**2), -np.inf, np.inf)[0]
check("S6 thimble: int exp(i x^2) dx along the real line = int along the line turned by 45 degrees (a Gaussian) = sqrt(pi) e^(i pi/4)",
      abs(fres - rot) < 5e-3 and np.isclose(rot, np.sqrt(np.pi)*np.exp(1j*np.pi/4)), f"{fres:.4f} vs {rot:.4f}")
# S7 Magnus/BCH: log(e^A e^B) = A + B + [A,B]/2 + ([A,[A,B]] - [B,[A,B]])/12 + O(4)
Am = 0.05*rng.normal(size=(3, 3)); Bm = 0.05*rng.normal(size=(3, 3))
from scipy.linalg import logm
cm = lambda X, Y: X @ Y - Y @ X
bch3 = Am + Bm + cm(Am, Bm)/2 + (cm(Am, cm(Am, Bm)) - cm(Bm, cm(Am, Bm)))/12
err2 = abs(logm(expm(Am) @ expm(Bm)) - (Am + Bm + cm(Am, Bm)/2)).max()
err3 = abs(logm(expm(Am) @ expm(Bm)) - bch3).max()
check("S7 Magnus/BCH: the log of an ordered product is the sum plus nested twists; the third-order term cuts the error by a lot",
      err3 < err2/10, f"error to 2nd order {err2:.1e}, to 3rd {err3:.1e}")
print(f"\n{sum(o for _, o in results)}/{len(results)} pass")
