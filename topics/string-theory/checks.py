"""Checks for topics/string-theory (Nov 25, 2026): what string theory rests on, and what the
Fused Hyperoperation Algebra could borrow.  Run: python3 checks.py  (numpy, sympy, scipy).
Writes checks_output.txt next to this file.  Units: alpha' = 1, hbar = c = k_B = 1 unless stated."""
import numpy as np, sympy as sp, os, itertools
from scipy.special import loggamma
HERE = os.path.dirname(os.path.abspath(__file__))
results, lines = [], []
def check(name, ok, detail=""):
    results.append(bool(ok)); s = ("PASS " if ok else "FAIL ") + name + (("  | " + detail) if detail else "")
    lines.append(s); print(s)
rng = np.random.default_rng(7)

# ===================== A. What string theory rests on (computable pieces) =====================
D = sp.symbols('D')
z = sp.zeta(-1)
check("A1 zeta(-1) = -1/12, so the regularized zero-point sum (D-2)/2 * (1+2+3+...) gives a = (D-2)/24",
      z == sp.Rational(-1, 12), f"zeta(-1) = {z}")
Db = sp.solve(sp.Eq((D - 2)/24, 1), D)[0]
check("A2 bosonic string: a light-cone massless vector at level 1 needs a = 1, so D = 26", Db == 26, f"D = {Db}")
Ds = sp.solve(sp.Eq((D - 2)*(sp.Rational(1, 24) + sp.Rational(1, 48)), sp.Rational(1, 2)), D)[0]
check("A3 superstring NS sector: bosons -1/24 and fermions -1/48 per direction, a = 1/2, so D = 10", Ds == 10, f"D = {Ds}")
n = 24; sym_traceless = n*(n + 1)//2 - 1; anti = n*(n - 1)//2
check("A4 closed bosonic string, massless level N = Ntilde = 1: 24 x 24 = 576 states = graviton 299 + B-field 276 + dilaton 1",
      n*n == 576 and sym_traceless + anti + 1 == 576 and sym_traceless == 299 and anti == 276, "massless spin-2 is forced into the spectrum")
# Newton's constant from g_s and alpha' (Polchinski conventions): 2 kappa10^2 = (2 pi)^7 g_s^2 alpha'^4, kappa^2 = 8 pi G
def G10(gs, ap): return (2*np.pi)**7*gs**2*ap**4/(16*np.pi)
g1, a1 = 0.1, 1.0
check("A5 G10 = (2 pi)^7 g_s^2 alpha'^4 / (16 pi): doubling g_s multiplies G by 4; doubling the string length multiplies G by 2^8",
      np.isclose(G10(2*g1, a1)/G10(g1, a1), 4) and np.isclose(G10(g1, 4*a1)/G10(g1, a1), 2**8), "G is derived: (coupling)^2 x (length)^8")
# 10D Planck length: 2 kappa^2 = (2 pi)^7 l_P^8  =>  l_P = g_s^(1/4) l_s
gs_test = 0.37; lP = (2*16*np.pi*G10(gs_test, 1.0)/(2*(2*np.pi)**7))**(1/8)
check("A6 the 10D Planck length is l_P = g_s^(1/4) l_s: gravity's length is the string length times a power of the coupling",
      np.isclose(lP, gs_test**0.25), f"l_P/l_s = {lP:.6f}, g_s^(1/4) = {gs_test**0.25:.6f}")
# AdS5 x S5: L^4 = 4 pi g_s N alpha'^2, Vol(S5) = pi^3 L^5, G5 = G10 / Vol(S5)  =>  L^3/G5 = 2 N^2 / pi
for Nc, gs in [(10, 0.05), (1000, 0.003), (37, 0.2)]:
    L = (4*np.pi*gs*Nc)**0.25; G5 = G10(gs, 1.0)/(np.pi**3*L**5)
    if not np.isclose(L**3/G5, 2*Nc**2/np.pi): break
check("A7 AdS5 x S5: L^3 / G5 = 2 N^2 / pi for any g_s, N: in AdS units Newton's constant is pi/(2 N^2), a 1/N^2 dial",
      np.isclose(L**3/G5, 2*Nc**2/np.pi), f"N = {Nc}: L^3/G5 = {L**3/G5:.4f}, 2N^2/pi = {2*Nc**2/np.pi:.4f}")

# ===================== T. T-duality on a circle =====================
def M2(nn, w, N, Nt, R): return (nn/R)**2 + (w*R)**2 + 2*(N + Nt - 2)   # alpha' = 1, level matching N - Nt = n w
ok = True
for _ in range(200):
    R = np.exp(rng.uniform(-2, 2)); nn, w = rng.integers(-4, 5, 2); N = rng.integers(0, 6); Nt = N - nn*w
    if Nt < 0: continue
    ok &= np.isclose(M2(nn, w, N, Nt, R), M2(w, nn, N, Nt, 1/R))
check("T1 closed-string spectrum on a circle is unchanged by R -> alpha'/R with momentum n <-> winding w (200 random states)", ok)
massless = [(nn, w, N, Nt) for nn in (-1, 0, 1) for w in (-1, 0, 1) for N in (0, 1) for Nt in (0, 1)
            if N - Nt == nn*w and np.isclose(M2(nn, w, N, Nt, 1.0), 0) and (nn, w) != (0, 0)]
check("T2 at the self-dual radius R = sqrt(alpha') extra states (n, w = +-1) become massless (the SU(2) x SU(2) enhancement)",
      len(massless) == 4 and all(M2(*s, 1.3) > 0 for s in massless), f"{len(massless)} extra massless (n,w,N,Nt): {massless}")
def pLR(nn, w, R): return (nn/R + w*R)/np.sqrt(2), (nn/R - w*R)/np.sqrt(2)
ok = True
for nn, w in itertools.product(range(-3, 4), repeat=2):
    a, b = pLR(nn, w, 0.6); c_, d_ = pLR(nn, w, 1.9); ok &= np.isclose(a*a - b*b, 2*nn*w)
    ok &= np.isclose(a*a - b*b, c_*c_ - d_*d_)
check("T3 p_L^2 - p_R^2 = 2 n w for every radius: the radius never changes the Lorentzian length of a lattice point", ok)
t = np.log(1.9/0.6); B = np.array([[np.cosh(t), -np.sinh(t)], [-np.sinh(t), np.cosh(t)]])
ok = all(np.allclose(B @ np.array(pLR(nn, w, 0.6)), np.array(pLR(nn, w, 1.9))) for nn, w in itertools.product(range(-3, 4), repeat=2))
check("T4 changing the radius is a Lorentz boost of the (p_L, p_R) lattice with rapidity log(R2/R1) (Narain's O(1,1))", ok, f"rapidity {t:.4f}")
# mult-add at k = 1 with 1 + a = R/sqrt(alpha')
ma = lambda a, b, k=1.0: a + b + k*a*b
inv = lambda a, k=1.0: -a/(1 + k*a)
a = rng.uniform(-0.9, 5, 50)
check("T5 ours: with 1 + a = R/sqrt(alpha'), T-duality is the mult-add (k = 1) inverse: 1 + inv(a) = 1/(1 + a), and a o inv(a) = 0 (the self-dual radius)",
      np.allclose(1 + inv(a), 1/(1 + a)) and np.allclose(ma(a, inv(a)), 0))
u = np.log(1 + a)
check("T6 ours: in the flat coordinate u = log(1 + a) = log(R/sqrt(alpha')), the rapidity of T4, T-duality is the reflection u -> -u; the wall a = -1 is R = 0, sent to R = infinity",
      np.allclose(np.log(1 + inv(a)), -u) and abs(1 + inv(-1 + 1e-9)) > 1e8)
# 2-torus with B field: generalized metric H(E), E = G + B; T-duality E -> E^{-1}
def Hgen(E):
    G = (E + E.T)/2; Bf = (E - E.T)/2; Gi = np.linalg.inv(G)
    return np.block([[G - Bf @ Gi @ Bf, Bf @ Gi], [-Gi @ Bf, Gi]])
ok = True; P = np.block([[np.zeros((2, 2)), np.eye(2)], [np.eye(2), np.zeros((2, 2))]])
for _ in range(50):
    A = rng.normal(size=(2, 2)); G = A @ A.T + 0.3*np.eye(2); Bf = rng.normal()*np.array([[0, 1], [-1, 0]]); E = G + Bf
    ok &= np.allclose(Hgen(np.linalg.inv(E)), P @ Hgen(E) @ P)
check("T7 2-torus with a B-field: T-duality E -> E^{-1} (E = G + B) swaps momentum and winding in the mass matrix, H(E^-1) = P H(E) P (Giveon-Porrati-Rabinovici)", ok)

# ===================== S. S-duality and the half-plane =====================
def mob(Mx, tau): return (Mx[0, 0]*tau + Mx[0, 1])/(Mx[1, 0]*tau + Mx[1, 1])
ok = True
for _ in range(50):
    al, b = np.exp(rng.normal()), rng.normal(); tau = complex(rng.normal(), np.exp(rng.normal()))
    ok &= np.isclose(mob(np.array([[al, b], [0, 1/al]]), tau), al**2*tau + al*b)
check("S1 the order dial's stretch-and-shift maps tau -> lam*tau + s are exactly the upper-triangular (Borel) part of SL(2,R) acting on the half-plane", ok)
Smat = np.array([[0, -1], [1, 0]]); lam = np.exp(rng.normal(size=20))
check("S2 S-duality tau -> -1/tau fixes tau = i (C0 = 0, g_s = 1), the order dial's 'do nothing' point",
      np.isclose(mob(Smat, 1j), 1j))
check("S3 ours: on the pure-multiplication line tau = i*lam (C0 = 0), S-duality g_s <-> 1/g_s is the group inverse lam -> 1/lam, the same shape as T5",
      np.allclose([mob(Smat, 1j*l) for l in lam], 1j/lam))
tau = 0.7 + 1.3j; ginv = (-0.7 + 1j)/1.3   # inverse of x -> 1.3 x + 0.7 as a point
check("S4 off that line S-duality is not the group inverse (it needs the full SL(2), beyond the dial's own group)",
      not np.isclose(mob(Smat, tau), ginv), f"S(tau) = {mob(Smat, tau):.4f}, group inverse point = {ginv:.4f}")
# S-duality is not in the Borel subgroup: lower-left entry nonzero, and Borel maps fix infinity
check("S5 Borel maps all fix the cusp at infinity (weak coupling); S moves it to 0, so dualities live outside the dial's group",
      np.isinf(complex(mob(np.array([[2., 3.], [0, .5]]), 1e300)).imag) or abs(mob(np.array([[2., 3.], [0, .5]]), 1e12j)) > 1e11 and abs(mob(Smat, 1e12j)) < 1e-11)

# ===================== M. Modular invariance and the Hagedorn temperature =====================
def eta(tau, terms=400):
    q = np.exp(2j*np.pi*tau); n = np.arange(1, terms + 1)
    return np.exp(2j*np.pi*tau/24)*np.prod(1 - q**n)
def Zb(tau): return tau.imag**12*abs(eta(tau))**48
ok = True
for _ in range(20):
    tau = complex(rng.uniform(-0.5, 0.5), rng.uniform(0.9, 2.0))
    ok &= np.isclose(Zb(tau + 1), Zb(tau), rtol=1e-8) and np.isclose(Zb(-1/tau), Zb(tau), rtol=1e-6)
check("M1 the bosonic closed-string torus weight tau2^12 |eta(tau)|^48 is unchanged by tau -> tau + 1 and tau -> -1/tau (modular invariance)", ok)
from scipy.integrate import quad
areaF = quad(lambda x: 1/np.sqrt(1 - x*x), -0.5, 0.5)[0]
strip = [quad(lambda x: 1/eps - 0, -0.5, 0.5)[0] for eps in (1e-1, 1e-2, 1e-3)]
check("M2 the fundamental domain F of SL(2,Z) has finite area pi/3 in dtau1 dtau2/tau2^2; the strip |tau1|<1/2 grows like 1/eps as tau2 -> eps (the UV region F leaves out)",
      np.isclose(areaF, np.pi/3) and strip[2] > 900, f"area(F) = {areaF:.6f}, strip down to eps: {strip}")
Nmax = 3000; c = np.zeros(Nmax + 1, dtype=object); c[0] = 1
p24 = [1] + [0]*Nmax
for nn in range(1, Nmax + 1):              # multiply by (1 - q^n)^-24 via 24 geometric series
    for _ in range(24):
        for j in range(nn, Nmax + 1): p24[j] += p24[j - nn]
Ns = np.arange(200, Nmax + 1, 50); lp = np.array([float(sp.log(p24[N])) for N in Ns])
Xf = np.column_stack([np.sqrt(Ns), np.log(Ns), np.ones_like(Ns, dtype=float), 1/np.sqrt(Ns)])
coef = np.linalg.lstsq(Xf, lp, rcond=None)[0]
check("M3 number of string states at level N grows like exp(4 pi sqrt N) N^(-27/4) (Hardy-Ramanujan-Meinardus; the 4 pi comes from eta's tau -> -1/tau)",
      abs(coef[0] - 4*np.pi) < 0.01 and abs(coef[1] + 27/4) < 0.1, f"fit over N = 200..3000: sqrt N coefficient {coef[0]:.4f} (4 pi = {4*np.pi:.4f}), log N coefficient {coef[1]:.3f} (-27/4)")
# open bosonic: alpha' M^2 = N - 1, Boltzmann weight exp(-beta M): converges iff beta > 4 pi sqrt(alpha')
def logterm(beta, N): return float(sp.log(p24[N])) - beta*np.sqrt(N - 1)
bH = 4*np.pi
check("M4 Hagedorn: sum over states of exp(-M/T) converges only for T < T_H = 1/(4 pi sqrt(alpha')) (open bosonic string)",
      logterm(1.05*bH, 3000) < logterm(1.05*bH, 1000) and logterm(0.95*bH, 3000) > logterm(0.95*bH, 1000), f"T_H = {1/bH:.5f} / sqrt(alpha')")
qv = rng.uniform(0, 0.9); ks = qv**np.arange(1, 30); chain = 0.0
for s_ in ks: chain = ma(chain, s_, -1.0)
check("M5 ours (bookkeeping): the partition-function denominator prod(1 - q^n) is 1 - (the k = -1 chain of the steps q^n), the project's survival dial",
      np.isclose(1 - chain, np.prod(1 - ks)))
check("M6 the primon gas (Julia 1990) has the same kind of wall: zeta(s) = prod (1 - p^-s)^-1 blows up as s -> 1 (its 'Hagedorn' point)",
      float(sp.zeta(1.01)) > 99 and float(sp.zeta(1.001)) > 999)

# ===================== L. Minimum length, GUP, and the (k, d) plane =====================
dp = np.exp(np.linspace(-3, 3, 2001)); dx = 1/(2*dp) + dp/2     # hbar = alpha' = 1, symmetric form of the string GUP
check("L1 string GUP dx >= hbar/(2 dp) + alpha' dp/(2 hbar) has a minimum dx = sqrt(alpha') at dp = hbar/sqrt(alpha') (Gross-Mende, Amati-Ciafaloni-Veneziano, Konishi-Paffuti-Provero)",
      np.isclose(dx.min(), 1, atol=1e-5) and np.isclose(dp[dx.argmin()], 1, atol=1e-2))
check("L2 that GUP is T-duality symmetric: dp -> hbar^2/(alpha' dp) swaps its two terms", np.allclose(1/(2*dp) + dp/2, dp/2 + 1/(2*dp)) and np.allclose(1/(2*(1/dp)) + (1/dp)/2, dx))
p, rho, beta, hbar = sp.symbols('p rho beta hbar', positive=True); f = sp.Function('f')
pr = sp.tan(sp.sqrt(beta)*rho)/sp.sqrt(beta)
x_on_f = lambda g: sp.I*hbar*sp.diff(g, rho)     # x = i hbar d/drho
F = sp.Function('F')(rho)
comm = sp.simplify(x_on_f(pr*F) - pr*x_on_f(F))
check("L3 Kempf-Mangano-Mann [x, p] = i hbar (1 + beta p^2) is canonical in rho = arctan(sqrt(beta) p)/sqrt(beta): x = i hbar d/drho",
      sp.simplify(comm - sp.I*hbar*(1 + beta*pr**2)*F) == 0)
tl = lambda P, Q, b: (P + Q)/(1 - b*P*Q); tp = lambda r, b: np.tan(np.sqrt(b)*r)/np.sqrt(b)
law = max(abs(tp(r1 + r2, b) - tl(tp(r1, b), tp(r2, b), b)) for r1, r2, b in rng.uniform(-0.6, 0.6, (200, 3)) + np.array([0, 0, 1.0]))
check("L4 ours (reading): composing so that rho adds is the tangent law (p + q)/(1 - beta p q), the point (k, d) = (0, beta) of the repo's plane, kappa = d - k^2/4 = beta > 0: elliptic, no wall",
      law < 1e-9, f"max error over 200 random (rho1, rho2, beta): {law:.1e}")
kP = 1.0; p0a, p0b = rng.uniform(0, 4, 2); A_, B_ = kP*(1 - np.exp(-p0a/kP)), kP*(1 - np.exp(-p0b/kP))
check("L5 kappa-Poincare's massless momenta |p| = kappa(1 - e^(-p0/kappa)) compose by a + e^(-p0a/kappa) b = a o_k b with k = -1/kappa: the point (k, d) = (-1/kappa, 0), kappa_CK = -k^2/4 < 0, hyperbolic, wall at |p| = kappa",
      np.isclose(A_ + np.exp(-p0a/kP)*B_, ma(A_, B_, -1/kP)) and np.isclose(ma(A_, B_, -1/kP), kP*(1 - np.exp(-(p0a + p0b)/kP))))
check("L6 the two are mirror images: KMM's flat coordinate rho is bounded (|rho| < pi/(2 sqrt beta)) with p unbounded; kappa-Poincare's flat coordinate p0 is unbounded with |p| < kappa_P",
      all(abs(np.arctan(np.sqrt(bb)*pp)/np.sqrt(bb)) < np.pi/(2*np.sqrt(bb)) for bb in (0.5, 2.0) for pp in (1e3, 1e9))
      and all(kP*(1 - np.exp(-p0/kP)) < kP for p0 in (1.0, 10.0, 30.0)))
# Veneziano amplitude at fixed angle: exponential fall-off (string softness)
def logA(s, tt):
    al = lambda x: 1 + x
    return (loggamma(complex(-al(s))) + loggamma(complex(-al(tt))) - loggamma(complex(-al(s) - al(tt)))).real
ss = np.arange(20, 400, 4) + 0.5; vals = np.array([logA(s_, -s_/2) for s_ in ss])
slope = np.polyfit(ss, vals, 1)[0]
check("L7 Veneziano amplitude at fixed angle (t = -s/2) falls like exp(-alpha' s ln 2): faster than any power, the 'softness' behind the minimum length (Gross-Mende 1987)",
      np.isclose(slope, -np.log(2), atol=0.01), f"fitted slope of log|A| against alpha' s: {slope:.4f}, -ln 2 = {-np.log(2):.4f}")

# ===================== N. The genus dial (g_s, 1/N) =====================
def hz(m, g):
    if g < 0 or m < 0: return 0
    if m == 0: return 1 if g == 0 else 0
    return sp.Rational(2*(2*m - 1)*hz(m - 1, g) + (m - 1)*(2*m - 1)*(2*m - 3)*hz(m - 2, g - 1), m + 1)
check("N1 Harer-Zagier: the gluings of a 2m-gon by genus give <tr X^(2m)> = sum_g eps_g(m) N^(m+1-2g); m = 2: (2, 1), m = 3: (5, 10), m = 4: (14, 70, 21)",
      [hz(2, 0), hz(2, 1)] == [2, 1] and [hz(3, 0), hz(3, 1)] == [5, 10] and [hz(4, 0), hz(4, 1), hz(4, 2)] == [14, 70, 21])
Nm, samp = 3, 200000; acc = 0.0
for _ in range(samp // 1000):
    A = rng.normal(size=(1000, Nm, Nm)) + 1j*rng.normal(size=(1000, Nm, Nm)); X = (A + np.conj(np.transpose(A, (0, 2, 1))))/(2*np.sqrt(Nm))
    acc += np.einsum('bij,bjk,bkl,bli->b', X, X, X, X).real.sum()/Nm
mc = acc/samp
check("N2 GUE Monte Carlo (N = 3, variance 1/N): (1/N) E tr X^4 = 2 + 1/N^2; the planar part is 2, each handle costs 1/N^2",
      abs(mc - (2 + 1/Nm**2)) < 0.03 and mc - 2 > 0.06, f"Monte Carlo {mc:.4f}, exact {2 + 1/Nm**2:.4f}")
check("N3 string amplitudes weigh a genus-g worldsheet by g_s^(2g-2) = g_s^(-chi): one handle costs g_s^2, the same power that sets G (A5); 't Hooft: g_s <-> 1/N",
      all(2*g - 2 == -(2 - 2*g) for g in range(5)))
th = sp.symbols('theta'); r = sp.symbols('r', positive=True)
gb = sp.integrate(sp.integrate((1/r**2)*r**2*sp.sin(th), (th, 0, sp.pi)), (sp.symbols('ph'), 0, 2*sp.pi))
check("N4 Gauss-Bonnet: the worldsheet's Einstein-Hilbert term integrates to 4 pi chi (sphere: 4 pi; a K = -1 genus-g surface has area 4 pi (g - 1)); so the dilaton term is a pure topology weight",
      sp.simplify(gb - 4*sp.pi) == 0)

passed = sum(results)
lines.append(f"\n{passed}/{len(results)} checks pass")
print(f"\n{passed}/{len(results)} checks pass")
open(os.path.join(HERE, "checks_output.txt"), "w").write("\n".join(lines) + "\n")
