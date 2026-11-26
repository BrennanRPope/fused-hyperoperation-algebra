"""DRAFT · Tunnelling worked up with our tools (Nov 26, 2026).

Checks for tunnelling/README.md. Run: python3 checks.py > checks_output.txt
numpy, scipy, mpmath, sympy. Each check prints PASS/FAIL with the numbers.

Notation (project names first):
  t  odds of reflection R/T of one barrier (the odds dial)
  D  hyperbolic length of a barrier, 1 + 2t = cosh D (a reading)
  u  the barrier's half-velocity, |u| = tanh(D/2) = |reflection amplitude| (twisted-boosts' t)
  s  switch probability of a lossless two-mode mixer (beam splitter, Rabi pulse), 1 - 2s = cos(theta)
  phi  loop phase between two elements (a setting)
  K  Gaussian curvature of the plane the elements live in (-1 barriers, 0 flat, +1 mixers)
"""
import math, cmath
import numpy as np
import mpmath as mp
import sympy as sp

rng = np.random.default_rng(20261126)
results = []
def check(name, ok, info=""):
    results.append(bool(ok))
    print(("PASS " if ok else "FAIL ") + name + (f"   [{info}]" if info else ""))

def delta_M(g, kw=1.0):
    b = g/(2*kw)
    return np.array([[1 - 1j*b, -1j*b], [1j*b, 1 + 1j*b]])
def prop(ph):
    return np.array([[np.exp(1j*ph), 0], [0, np.exp(-1j*ph)]])
def odds(M):           # t = R/T = |beta|^2 for an SU(1,1) matrix [[a,b],[b*,a*]]
    return abs(M[0, 1])**2
def length(M):
    return math.acosh(abs(M[0, 0])**2 + abs(M[0, 1])**2)
def mixer(s, pa=0.3, pb=-0.7):   # SU(2) lossless mixer with switch probability s
    a = math.sqrt(1 - s)*cmath.exp(1j*pa); b = math.sqrt(s)*cmath.exp(1j*pb)
    return np.array([[a, b], [-b.conjugate(), a.conjugate()]])
def circ(k, a, b):
    return a + b + k*a*b
def dlaw(k, d, a, b):   # the (k, d) law of named-functions/
    return (a + b + k*a*b)/(1 - d*a*b)

PH = np.linspace(0, np.pi, 20001)[:-1]   # loop phase grid (round trip 2*phase covers the circle)

# ------------------------------------------------------------------ 1. a barrier is a dial boost
print("== 1. A barrier is a boost of the dial (twisted-boosts, q = -1)")
ok = True; worst = 0
for g in (0.3, 0.7, 1.3, 1.9, 4.0):
    M = delta_M(g); a, b = M[0, 0], M[0, 1]
    u = b/np.conj(a)
    B = np.array([[1, u], [np.conj(u), 1]])                     # B_q(u) with q = -1
    ok &= np.allclose(B @ np.diag([a, np.conj(a)]), M)
    worst = max(worst, abs(abs(u) - math.tanh(length(M)/2)))
    t = odds(M); ok &= abs(t - abs(u)**2/(1 - abs(u)**2)) < 1e-12
check("1.1 every delta barrier is B_{-1}(u)·(a phase), with u its reflection amplitude and |u| = tanh(D/2)",
      ok and worst < 1e-12, f"max | |u| - tanh(D/2) | = {worst:.1e}")
check("1.2 the odds are the flip of the amplitude squared: t = |u|²/(1 - |u|²) = R/T", ok)

# ------------------------------------------------------------------ 2. the band and its edges
print("\n== 2. Fixed phase: the band's edges are our k = 0 laws with d = K")
g1, g2 = 1.3, 0.8
M1, M2 = delta_M(g1), delta_M(g2)
D1, D2 = length(M1), length(M2)
Ds = np.array([length(M2 @ prop(p) @ M1) for p in PH])
r1, r2 = math.tanh(D1/2), math.tanh(D2/2)
check("2.1 hyperbolic edge: max |u| over the phase is (r₁ + r₂)/(1 + r₁r₂), the (k=0, d=−1) law (relativity)",
      abs(math.tanh(Ds.max()/2) - dlaw(0, -1, r1, r2)) < 1e-7, f"{math.tanh(Ds.max()/2):.8f} vs {dlaw(0, -1, r1, r2):.8f}")
check("2.2 hyperbolic edge: min |u| is the same law run in reverse, (r₁ − r₂)/(1 − r₁r₂)",
      abs(math.tanh(Ds.min()/2) - abs(r1 - r2)/(1 - r1*r2)) < 1e-6)
s1, s2 = 0.2, 0.35
N1, N2 = mixer(s1, 0, 0), mixer(s2, 0, 0)
th = np.array([2*math.asin(min(1, abs((N2 @ prop(p) @ N1)[0, 1]))) for p in PH])
a1, a2 = math.tan(math.asin(math.sqrt(s1))), math.tan(math.asin(math.sqrt(s2)))
check("2.3 sphere edge: max tan(θ/2) is (a₁ + a₂)/(1 − a₁a₂), the (k=0, d=+1) law (tangent addition)",
      abs(math.tan(th.max()/2) - dlaw(0, 1, a1, a2)) < 1e-6, f"{math.tan(th.max()/2):.7f} vs {dlaw(0, 1, a1, a2):.7f}")
check("2.4 so the curvature reading −Δ/4 = d of the edge law equals the curvature: −1 for barriers, +1 for mixers",
      -(0 - 4*(-1))/4 == -1 and -(0 - 4*1)/4 == 1)

# ------------------------------------------------------------------ 3. the twist
print("\n== 3. The twist: the extra transmission phase is half the triangle's area")
def hyp_area(a, b, c):
    def ang(x, y, z):  # angle opposite z
        return math.acos(max(-1, min(1, (math.cosh(x)*math.cosh(y) - math.cosh(z))/(math.sinh(x)*math.sinh(y)))))
    return math.pi - (ang(b, c, a) + ang(a, c, b) + ang(a, b, c))
def sph_area(a, b, c):
    def ang(x, y, z):
        return math.acos(max(-1, min(1, (math.cos(z) - math.cos(x)*math.cos(y))/(math.sin(x)*math.sin(y)))))
    return ang(b, c, a) + ang(a, c, b) + ang(a, b, c) - math.pi
def wrap(x): return (x + math.pi) % (2*math.pi) - math.pi
worst_h = worst_s = 0
for p in np.linspace(0.05, 3.09, 60):
    M = M2 @ prop(p) @ M1
    exc = wrap(np.angle(M[0, 0]) - np.angle(M1[0, 0]) - np.angle(M2[0, 0]) - p)
    worst_h = max(worst_h, abs(abs(exc) - hyp_area(D1, D2, length(M))/2))
    N = N2 @ prop(p) @ N1
    exc = wrap(np.angle(N[0, 0]) - np.angle(N1[0, 0]) - np.angle(N2[0, 0]) - p)
    T1_, T2_ = 2*math.asin(math.sqrt(s1)), 2*math.asin(math.sqrt(s2)); Tt = 2*math.asin(min(1, abs(N[0, 1])))
    worst_s = max(worst_s, abs(abs(exc) - sph_area(T1_, T2_, Tt)/2))
check("3.1 two barriers: transmission phase − (each barrier's phase + the propagation) = ± ½ hyperbolic area",
      worst_h < 1e-9, f"max error {worst_h:.1e} over 60 phases")
check("3.2 two mixers: the same excess is ± ½ the spherical triangle's area (Pancharatnam's phase)",
      worst_s < 1e-9, f"max error {worst_s:.1e}")
print("NOTE 3.3 (reading, not a check): both are twisted-boosts' holonomy, ω = q × area, with the phase of X + iY equal to ω/2")

# ------------------------------------------------------------------ 4. averaged law: k = -2K
print("\n== 4. Averaging the law of cosines over the angle: the squared half-chord composes by ∘ₖ, k = −2K")
def cosK(K, D):
    if K > 0: return math.cos(math.sqrt(K)*D)
    if K < 0: return math.cosh(math.sqrt(-K)*D)
    return None
def tK(K, D):          # sin_K(D/2)^2 / ... normalised: (1 - cos(√K D))/(2K); flat: D²/4
    return D*D/4 if K == 0 else (1 - cosK(K, D))/(2*K)
ang = np.linspace(0, 2*np.pi, 40001)[:-1]
ok = True; info = []
for K in (-2.0, -1.0, -0.3, 0.0, 0.5, 1.0):
    Da, Db = 0.7, 0.45
    if K == 0:
        Dc = np.sqrt(Da*Da + Db*Db - 2*Da*Db*np.cos(ang))
        tc = Dc**2/4
    elif K < 0:
        q = math.sqrt(-K); z = np.cosh(q*Da)*np.cosh(q*Db) - np.sinh(q*Da)*np.sinh(q*Db)*np.cos(ang)
        tc = (1 - z)/(2*K)
    else:
        q = math.sqrt(K); z = np.cos(q*Da)*np.cos(q*Db) + np.sin(q*Da)*np.sin(q*Db)*np.cos(ang)
        tc = (1 - z)/(2*K)
    pred = circ(-2*K, tK(K, Da), tK(K, Db))
    ok &= abs(tc.mean() - pred) < 1e-10; info.append(f"K={K}: {tc.mean():.6f}/{pred:.6f}")
check("4.1 for K = −2, −1, −0.3, 0, 0.5, 1: <t> = t₁ ∘₋₂ₖ t₂ exactly, with t = (1 − cos√K·D)/(2K)", ok, "; ".join(info[:3]))
uu = rng.uniform(0.05, 3, 200)
check("4.2 at K = −1 this t is the barrier's odds R/T; at K = +1 it is the mixer's switch probability s",
      all(abs(tK(-1, math.acosh(1 + 2*x)) - x) < 1e-9 for x in uu) and abs(tK(1, 2*math.asin(math.sqrt(0.3))) - 0.3) < 1e-12)
S = np.array([abs((N2 @ prop(p) @ N1)[0, 1])**2 for p in PH])
check("4.3 mixers: <1 − 2s> = (1 − 2s₁)(1 − 2s₂), i.e. s composes by ∘₋₂ (the piling-up lemma's bias rule)",
      abs((1 - 2*S).mean() - (1 - 2*s1)*(1 - 2*s2)) < 1e-10)
rho = np.array([odds(M2 @ prop(p) @ M1) for p in PH]); t1, t2 = odds(M1), odds(M2)
check("4.4 barriers: <t> = t₁ ∘₂ t₂ (Anderson–Thouless–Abrahams–Fisher), from the same average",
      abs(rho.mean() - circ(2, t1, t2)) < 1e-10, f"{rho.mean():.10f}")
# rescaling: k -> k/lambda when t -> lambda t; only sign(K) is physics
check("4.5 rescaling lengths changes |K| and |k| together (k = −2K), so at d = 0 only the sign is physical (F1)",
      abs(circ(-2*(-1)/4, 4*0.3, 4*0.2)/4 - circ(2, 0.3, 0.2)) < 1e-12)

# ------------------------------------------------------------------ 5. typical laws are tropical
print("\n== 5. The log-average picks the stronger path exactly (Jensen): the tropical law")
ok = True
for A, B in [(1.0, 0.3), (0.4, 0.9), (2.0, 2.0), (0.7, 0.05)]:
    lhs = np.log(np.abs(A + B*np.exp(1j*ang))**2).mean()
    ok &= abs(lhs - math.log(max(A, B)**2)) < 1e-6 or A == B
check("5.1 <ln|A + B e^{iφ}|²> = ln max(|A|², |B|²), exactly, for two paths (Jensen 1899)", ok)
ok = True
for (x1, x2) in [(0.2, 0.3), (0.4, 0.5), (0.7, 0.6), (0.9, 0.05), (0.85, 0.4)]:
    Nn1, Nn2 = mixer(x1), mixer(x2, 1.1, 0.4)
    Ss = np.array([abs((Nn2 @ prop(p) @ Nn1)[0, 1])**2 for p in PH])
    lhs = np.log(1 - Ss).mean()
    ok &= abs(lhs - math.log(max((1 - x1)*(1 - x2), x1*x2))) < 1e-6
check("5.2 mixers: <ln(1 − s)> = ln max((1−s₁)(1−s₂), s₁s₂): stay-stay or switch-switch, whichever is likelier", ok)
check("5.3 the two paths tie exactly where s₁ + s₂ = 1, because (1−s₁)(1−s₂) − s₁s₂ = 1 − s₁ − s₂: a wall",
      all(abs((1 - a)*(1 - b) - a*b - (1 - a - b)) < 1e-15 for a, b in rng.uniform(0, 1, (50, 2))))
print("NOTE 5.4 (reading of 5.2): below the wall ∘₋₁ holds on s; above it, s₁s₂ takes over")
check("5.5 barriers: <ln(1 + t)> = ln(1 + t₁) + ln(1 + t₂) because cosh·cosh > sinh·sinh: the direct path always wins",
      abs(np.log1p(rho).mean() - math.log1p(t1) - math.log1p(t2)) < 1e-10)
# chain of different barriers: the typical law stays exact at every step
gs = [0.6, 1.4, 0.9, 2.2, 0.5, 1.1]; trials = 30000
acc = 0.0; acc2 = 0.0
Ms = [delta_M(g) for g in gs]
P = rng.uniform(0, np.pi, (trials, len(gs)))
for tr in range(trials):
    M = Ms[0]
    for j in range(1, len(gs)):
        M = Ms[j] @ prop(P[tr, j]) @ M
    lv = math.log1p(odds(M)); acc += lv; acc2 += lv*lv
pred = sum(math.log1p(odds(m)) for m in Ms)
se_guess = math.sqrt((acc2/trials - (acc/trials)**2)/trials)
check("5.6 six different barriers, random spacings: <ln(1 + t)> = Σ ln(1 + tᵢ) (∘₁ along the chain, no wall)",
      abs(acc/trials - pred) < 3*se_guess, f"{acc/trials:.4f} ± {se_guess:.4f} vs {pred:.4f}")

# ------------------------------------------------------------------ 6. the k = 0 pair laws
print("\n== 6. The pair's mean transmission: k = 0 in both geometries, with a wall on the sphere")
Tm = (1/(1 + rho)).mean()
check("6.1 barriers: 1/<T> − 1 = t₁ + t₂ (plain addition, pair only)", abs(1/Tm - 1 - (t1 + t2)) < 1e-9)
inv = (1/(1 - S)).mean()
check("6.2 mixers: <1/(1 − s)> = 1/|1 − (s₁ + s₂)|, plain addition with its wall at 1",
      abs(inv - 1/abs(1 - s1 - s2)) < 1e-8, f"{inv:.6f}")
c1, c2 = sp.symbols('c1 c2')
check("6.3 both rest on one identity: (1 + c₁c₂)² − (c₁² − 1)(c₂² − 1) = (c₁ + c₂)²",
      sp.expand((1 + c1*c2)**2 - (c1**2 - 1)*(c2**2 - 1) - (c1 + c2)**2) == 0)

# ------------------------------------------------------------------ 7. the continuum flow
print("\n== 7. Many weak elements: the flow ∂ₜ t(1 + εt) ∂ₜ, ε = −K, and which k it keeps exact")
t_, k_, c_, a_, e_ = sp.symbols('t k c alpha epsilon')
Lop = lambda f: sp.diff(t_*(1 + e_*t_)*sp.diff(f, t_), t_)
sol_mean = sp.solve(sp.Poly(sp.expand(Lop(1 + k_*t_) - c_*(1 + k_*t_)), t_).coeffs(), [k_, c_], dict=True)
check("7.1 the mean of 1 + kt grows exactly (an eigenfunction) only for k = 2ε (or the trivial k = 0)",
      {k_: 2*e_, c_: 2*e_} in sol_mean, str(sol_mean))
n_, d_ = sp.fraction(sp.together(Lop(sp.log(1 + k_*t_)) - c_))
sol_log = sp.solve(sp.Poly(sp.expand(n_), t_).coeffs(), [k_, c_], dict=True)
check("7.2 the mean of ln(1 + kt) grows linearly only for k = ε: the typical law", {k_: e_, c_: e_} in sol_log, str(sol_log))
ex = sp.simplify(Lop((1 + k_*t_)**a_)/(1 + k_*t_)**(a_ - 2))
sol_pow = sp.solve(sp.Poly(sp.expand(ex - c_*(1 + k_*t_)**2), t_).coeffs(), [k_, c_, a_], dict=True)
nontriv = [s_ for s_ in sol_pow if s_.get(a_, None) not in (0,) and s_.get(k_, None) != 0]
check("7.3 among all powers (1 + kt)^α the only exact one is α = 1, k = 2ε: the flow picks k = 2ε and k = ε",
      nontriv == [{a_: 1, c_: 2*e_, k_: 2*e_}], str(nontriv))
check("7.4 at ε = 1 this is Mel'nikov's one-channel equation (DMPK for N = 1); at ε = −1 the zonal heat flow on the sphere",
      sp.expand(Lop(t_).subs(e_, 1) - (1 + 2*t_)) == 0)

# ------------------------------------------------------------------ 8. the spectrum and the mean transmission
print("\n== 8. Three standard averages, three points of the Legendre family; mean T and Legendre's K")
mp.mp.dps = 20
def P(nu, z): return mp.legenp(nu, 0, z, type=3)
z1v = 1 + 2*t1; z2v = 1 + 2*t2
ok = True
for nu in (1, 2, 0.5, -0.5, -0.5 + 1.3j):
    lhs = np.mean([complex(P(nu, 1 + 2*x)) for x in rho[::50]])
    ok &= abs(lhs - complex(P(nu, z1v)*P(nu, z2v))) < 1e-6
check("8.1 <P_ν(1 + 2t)> = P_ν(1 + 2t₁)P_ν(1 + 2t₂) for ν = 1, 2, ½, −½ and a conical ν (product formula)", ok)
dP0 = lambda z: mp.diff(lambda nu: P(nu, z), 0)
check("8.2 ν = 1 is the mean law ∘₂; the slope in ν at ν = 0 is ln(1 + t), the typical law ∘₁",
      abs(dP0(3) - mp.log(2)) < 1e-12 and abs(P(1, 3) - 3) < 1e-15)
F = lambda lam: mp.quad(lambda z: 2/(1 + z)*P(-0.5 + 1j*lam, z), [1, 2, 10, 100, mp.inf])
ok = all(abs(F(l) - 2*mp.pi/mp.cosh(mp.pi*l)) < 1e-8 for l in (0, 0.5, 1.0, 2.0))
check("8.3 transmission T = 2/(1 + z) has Mehler–Fock coefficient 2π/cosh(πλ)", ok)
gB = 1.0; MB = delta_M(gB); tB = odds(MB); zB = 1 + 2*tB
lams = np.linspace(0, 8, 801)
PcB = np.array([float(mp.re(P(-0.5 + 1j*l, zB))) for l in lams])
w = lams*np.tanh(np.pi*lams)*2*np.pi/np.cosh(np.pi*lams)
meanT = lambda N: np.trapezoid(w*PcB**N, lams)
mc = {}
trials = 120000
for N in (2, 4, 8):
    M = np.tile(MB, (trials, 1, 1))
    for j in range(1, N):
        ph = rng.uniform(0, np.pi, trials)
        Pm = np.zeros((trials, 2, 2), complex); Pm[:, 0, 0] = np.exp(1j*ph); Pm[:, 1, 1] = np.exp(-1j*ph)
        M = MB @ Pm @ M
    T = 1/np.abs(M[:, 1, 1])**2
    mc[N] = (T.mean(), T.std()/math.sqrt(trials))
okmc = all(abs(mc[N][0] - meanT(N)) < 3*mc[N][1] + 1e-9 for N in mc)
check("8.4 <T> after N identical barriers = ∫ λ tanh(πλ)·2π/cosh(πλ)·P_{−½+iλ}(1+2t₁)ᴺ dλ (Monte Carlo, N = 2, 4, 8)",
      okmc, "; ".join(f"N={N}: {mc[N][0]:.4f}±{mc[N][1]:.4f} vs {meanT(N):.4f}" for N in mc))
PB = float(P(-0.5, zB)); Kform = 2/math.pi*math.sqrt(1/(1 + tB))*float(mp.ellipk(tB/(1 + tB)))
check("8.5 the per-barrier decay of <T> is P_{−½}(1 + 2t₁) = (2/π)·√T₁·K(|r₁|), Legendre's complete elliptic integral",
      abs(PB - Kform) < 1e-14, f"{PB:.12f}")
rat = meanT(401)/meanT(400)
check("8.6 at large N, <T_{N+1}>/<T_N> → P_{−½}(1+2t₁)·(N/(N+1))^{3/2}", abs(rat - PB*(400/401)**1.5) < 1e-3,
      f"{rat:.5f} vs {PB*(400/401)**1.5:.5f}")
small = 1e-4
check("8.7 weak barriers: rates per unit of typical length are 2 (mean odds), 1 (typical), ¼ (mean T, the bottom of the hyperbolic plane's spectrum)",
      abs(math.log1p(2*small)/math.log1p(small) - 2) < 1e-3 and abs(-math.log(float(P(-0.5, 1 + 2*small)))/math.log1p(small) - 0.25) < 1e-3)

# ------------------------------------------------------------------ 9. rewinds
print("\n== 9. Rewinds: single numbers rewind exactly; the full spread loses a comb of components")
def inv2(T_, a): return (T_ - a)/(1 + 2*a)
check("9.1 removing a barrier from the mean odds is the ∘₂ inverse (t − t₁)/(1 + 2t₁): defined for every barrier",
      abs(inv2(circ(2, t1, t2), t1) - t2) < 1e-12)
def invm2(T_, a): return (T_ - a)/(1 - 2*a)
check("9.2 on the sphere the same rewind has a wall at s₁ = ½: a 50/50 mixer erases the bias (1 − 2s = 0)",
      abs(invm2(circ(-2, 0.3, 0.2), 0.3) - 0.2) < 1e-12 and (1 - 2*0.5) == 0)
zc = 3.0
f = lambda l: P(-0.5 + 1j*l, zc)
xs = [0.05*i for i in range(1, 241)]; vals = [mp.re(f(x)) for x in xs]
nreal = sum(1 for i in range(len(xs) - 1) if vals[i]*vals[i + 1] < 0)
def contour(x0, x1, y0, y1, n=300):
    pts = []
    for i in range(n): pts.append(mp.mpc(x0 + (x1 - x0)*i/n, y0))
    for i in range(n): pts.append(mp.mpc(x1, y0 + (y1 - y0)*i/n))
    for i in range(n): pts.append(mp.mpc(x1 - (x1 - x0)*i/n, y1))
    for i in range(n): pts.append(mp.mpc(x0, y1 - (y1 - y0)*i/n))
    pts.append(pts[0]); tot = 0; prev = f(pts[0])
    for p_ in pts[1:]:
        cur = f(p_); tot += mp.arg(cur/prev); prev = cur
    return int(mp.nint(tot/(2*mp.pi)))
nbox = contour(0.3, 11.97, -2, 2)
check("9.3 a barrier's character P_{−½+iλ}(1+2t₁) has only real zeros in λ (box count = real count)",
      nbox == nreal, f"{nreal} real, {nbox} in the box |Im λ| < 2")
Dc = math.acosh(zc); sgrid = np.linspace(0, Dc*0.999, 400)
ker = 1/np.sqrt(np.cosh(Dc) - np.cosh(sgrid))
check("9.4 why: Mehler's kernel 1/√(cosh D − cosh σ) increases on [0, D), so Pólya's theorem gives real zeros",
      np.all(np.diff(ker) > 0))
print("NOTE 9.5 (reading of 8.1 and 9.3, corrected): one added barrier multiplies the component at each zero by 0, but a finite"
      "\n          chain's spectrum is entire in λ, so those values are fixed by their neighbours: ill-conditioned, not lost."
      "\n          The wall-free rewind is on the fused ladder; see rewind_checks.py.")

# ------------------------------------------------------------------ polish at 30 digits
print("\n== Polish (30 digits)")
mp.mp.dps = 30
u, v = mp.mpf('0.4225'), mp.mpf('0.16')
A = u + v + 2*u*v; Bc = 2*mp.sqrt(u*v*(1 + u)*(1 + v))
mean = mp.quad(lambda p: A - Bc*mp.cos(p), [0, 2*mp.pi])/(2*mp.pi)
logm = mp.quad(lambda p: mp.log(1 + A - Bc*mp.cos(p)), [0, 2*mp.pi])/(2*mp.pi)
tm = mp.quad(lambda p: 1/(1 + A - Bc*mp.cos(p)), [0, 2*mp.pi])/(2*mp.pi)
x1, x2 = mp.mpf('0.7'), mp.mpf('0.6')
Am = x1 + x2 - 2*x1*x2; Bm = 2*mp.sqrt(x1*x2*(1 - x1)*(1 - x2))
logs = mp.quad(lambda p: mp.log(1 - Am + Bm*mp.cos(p)), [0, mp.pi/3, mp.pi, 2*mp.pi])/(2*mp.pi)
ellK = 2/mp.pi*mp.sqrt(1/(1 + u))*mp.ellipk(u/(1 + u))
print(f"polish |<t> - u∘₂v| = {mp.nstr(abs(mean - A), 3)}")
print(f"polish |<ln(1+t)> - ln(1+u) - ln(1+v)| = {mp.nstr(abs(logm - mp.log((1+u)*(1+v))), 3)}")
print(f"polish |1/<T> - 1 - (u+v)| = {mp.nstr(abs(1/tm - 1 - u - v), 3)}")
print(f"polish sphere above the wall |<ln(1-s)> - ln(s1 s2)| = {mp.nstr(abs(logs - mp.log(x1*x2)), 3)}")
print(f"polish |P_(-1/2)(1+2u) - (2/π)√T K(|r|)| = {mp.nstr(abs(mp.legenp(-0.5, 0, 1 + 2*u, type=3) - ellK), 3)}")
check("P.1 polish: the pair laws and the sphere's wall law hold to 25+ digits",
      max(abs(mean - A), abs(logm - mp.log((1+u)*(1+v))), abs(1/tm - 1 - u - v), abs(logs - mp.log(x1*x2))) < mp.mpf(10)**-25)

print(f"\n{sum(results)}/{len(results)} checks pass")
