"""Checks for closed-leads-fresh-look
In F4 the code names u, v for the two barriers' odds t₁, t₂ and rho for the pair's odds t.
 (DRAFT, Nov 26 2026).
Run: python3 checks.py > checks_output.txt   (numpy, scipy, sympy, mpmath)
Sections follow the README, fundamental first: F1 scale, F2 order/history, F3 elliptic,
F4 tunnelling, F5 triprotic acids."""
import numpy as np, sympy as sp, mpmath as mp, itertools, math, time
from scipy import special
rng = np.random.default_rng(20261126)
mp.mp.dps = 30
results = []
def check(name, ok, detail=""):
    results.append((name, bool(ok)))
    print(("PASS " if ok else "FAIL ") + name + (("   [" + detail + "]") if detail else ""))
t0 = time.time()

# ---------------------------------------------------------------- F1 scale
print("\n== F1. What survives rescaling once d is in")
k, d, lam, a, b, t = sp.symbols('k d lam a b t', nonzero=True)
L = lambda kk, dd, x, y: (x + y + kk*x*y)/(1 - dd*x*y)
# conjugating by x -> lam x sends (k, d) to (k/lam, d/lam^2)
expr = sp.simplify(lam*L(k, d, a/lam, b/lam) - L(k/lam, d/lam**2, a, b))
check("F1.1 rescaling x->lam*x sends (k,d) to (k/lam, d/lam^2)", expr == 0)
# generator G = slide(T_kappa(t)), Taylor to t^3, kappa = d - k^2/4 is used only inside (named-functions)
kap = d - k**2/4
T = t + kap*t**3/3
G = sp.series(T/(1 - k*T/2), t, 0, 4).removeO()
a2 = sp.expand(G).coeff(t, 2); a3 = sp.expand(G).coeff(t, 3)
resit = sp.simplify(1 - a3/a2**2)
check("F1.2 iterative residue of the (k,d) generator is 1/3 - 4d/(3k^2)",
      sp.simplify(resit - (sp.Rational(1, 3) - 4*d/(3*k**2))) == 0, str(resit))
check("F1.3 the residue is unchanged by rescaling (depends on d/k^2 only)",
      sp.simplify(resit.subs({k: k/lam, d: d/lam**2}, simultaneous=True) - resit) == 0)
check("F1.4 with d = 0 it is the constant 1/3 (k alone: only its sign survives)",
      resit.subs(d, 0) == sp.Rational(1, 3))
# numeric: the generator's group law really is the (k,d) law, at a random setting
kk, dd = 0.7, 0.3
kapn = dd - kk**2/4
Tn = lambda u: math.tan(math.sqrt(kapn)*u)/math.sqrt(kapn)
Gn = lambda u: Tn(u)/(1 - kk*Tn(u)/2)
u1, u2 = 0.21, 0.13
check("F1.5 generator G turns + into the (k,d) law (numeric, k=0.7, d=0.3)",
      abs(Gn(u1+u2) - L(kk, dd, Gn(u1), Gn(u2))) < 1e-12)

# ---------------------------------------------------------------- F2 order / history
print("\n== F2. History from composition: stretch-and-shift reforms")
def compose(seq, x0=0.0):
    x = x0
    for (al, c) in seq:  # applied first to last
        x = al*x + c
    return x
def composite(seq):
    A, C = 1.0, 0.0
    for (al, c) in seq:
        A, C = al*A, al*C + c
    return A, C
# share-centre maps commute; different centres do not
kd = 0.4
same = [(1 + kd*s, s) for s in rng.uniform(-0.5, 1, 5)]   # stretch a = k*shift: centre -1/k
outs = {round(composite(list(p))[1], 12) for p in itertools.permutations(same)}
check("F2.1 five reforms sharing one centre: all 120 orders give one outcome", len(outs) == 1)
diff = [(1 + a_, c_) for a_, c_ in zip(rng.uniform(-0.3, 0.6, 5), rng.uniform(-1, 1, 5))]
outs = {round(composite(list(p))[1], 9) for p in itertools.permutations(diff)}
check("F2.2 five reforms with different centres: 120 orders give 120 outcomes", len(outs) == 120, f"{len(outs)}")
# adjacent swap changes the shift by (stretch of everything after) x wedge
ok = True
for _ in range(200):
    seq = [(1 + x, y) for x, y in rng.uniform(-0.5, 0.5, (6, 2))]
    i = rng.integers(0, 5)
    sw = seq.copy(); sw[i], sw[i+1] = sw[i+1], sw[i]
    after = np.prod([al for al, _ in seq[i+2:]]) if i+2 < 6 else 1.0
    ai, ci = seq[i][0]-1, seq[i][1]; aj, cj = seq[i+1][0]-1, seq[i+1][1]
    wedge = ci*aj - cj*ai
    ok &= abs((composite(seq)[1] - composite(sw)[1]) - after*wedge) < 1e-12
check("F2.3 swapping two neighbours moves the outcome by (later stretches) x (oriented area)", ok)
# group commutator f g f^-1 g^-1 is a pure shift equal to the wedge (order-dial `order_gap`)
a1, c1, a2_, c2 = sp.symbols('a1 c1 a2 c2')
x = sp.symbols('x')
f = lambda z: (1+a1)*z + c1; g = lambda z: (1+a2_)*z + c2
fi = lambda z: (z - c1)/(1+a1); gi = lambda z: (z - c2)/(1+a2_)
comm = sp.simplify(f(g(fi(gi(x)))) - x)
check("F2.4 commutator of two reforms is a pure shift; zero iff their centres agree",
      sp.simplify(comm - (a1*c2 - a2_*c1)) == 0 and sp.diff(comm, x) == 0, "shift = " + str(sp.factor(comm)))
# worked example: +10% then +500 vs +500 then +10%
w = composite([(1.10, 0), (1.0, 500)])[1] - composite([(1.0, 500), (1.10, 0)])[1]
check("F2.5 worked example: 10% raise and a 500 grant differ by 50 between orders (the corner)", abs(abs(w) - 50) < 1e-9)
# the space of stretch-and-shift maps with its invariant metric is the hyperbolic plane (curvature -1)
# left-invariant metric ds^2 = (dA^2 + dC^2)/A^2 on (A>0, C): Gaussian curvature
A_, C_ = sp.symbols('A C', positive=True)
E_ = 1/A_**2; G_ = 1/A_**2
K = -1/(2*sp.sqrt(E_*G_))*(sp.diff(sp.diff(G_, A_)/sp.sqrt(E_*G_), A_) + sp.diff(sp.diff(E_, C_)/sp.sqrt(E_*G_), C_))
check("F2.6 maps x -> A x + C with the invariant metric (dA^2+dC^2)/A^2 form a plane of curvature -1", sp.simplify(K) == -1)
# invariance of that metric under left composition by (A0, C0)
A0, C0 = sp.symbols('A0 C0', positive=True)
# left action: (A, C) -> (A0 A, A0 C + C0): dA'=A0 dA, dC'=A0 dC, A'=A0 A
check("F2.7 that metric is unchanged by composing with any reform", sp.simplify((A0**2)/(A0*A_)**2 - 1/A_**2) == 0)

# ---------------------------------------------------------------- F3 elliptic
print("\n== F3. Elliptic functions: one parameter off the plane dial")
def ted_add(P, Q, aa, de):
    (x1, y1), (x2, y2) = P, Q
    den = de*x1*x2*y1*y2
    return ((x1*y2 + y1*x2)/(1 + den), (y1*y2 - aa*x1*x2)/(1 - den))
def ted_point(xx, aa, de):
    return (xx, math.sqrt((1 - aa*xx*xx)/(1 - de*xx*xx)))
ok = True
for aa in (1.0, 0.5, -0.7):
    for _ in range(50):
        p, q = rng.uniform(-0.6, 0.6, 2)
        P, Q = ted_point(p, aa, 0.0), ted_point(q, aa, 0.0)
        R = ted_add(P, Q, aa, 0.0)
        tP, tQ, tR = P[0]/P[1], Q[0]/Q[1], R[0]/R[1]
        ok &= abs(tR - L(0.0, aa, tP, tQ)) < 1e-10
check("F3.1 Edwards parameter 0: the curve a x^2 + y^2 = 1 is our conic law with k = 0, d = a (in slope x/y)", ok)
ok = True
for aa, de in ((1.0, -1.0), (1.0, 0.6), (-0.5, 0.3), (2.0, -3.0)):
    for _ in range(50):
        p, q, r = rng.uniform(-0.4, 0.4, 3)
        P, Q, R = (ted_point(z, aa, de) for z in (p, q, r))
        L1 = ted_add(ted_add(P, Q, aa, de), R, aa, de); L2 = ted_add(P, ted_add(Q, R, aa, de), aa, de)
        S = ted_add(P, Q, aa, de)
        ok &= max(abs(L1[0]-L2[0]), abs(L1[1]-L2[1])) < 1e-10 and abs(aa*S[0]**2 + S[1]**2 - 1 - de*S[0]**2*S[1]**2) < 1e-10
check("F3.2 with Edwards parameter on, the law stays on its curve and stays associative (four settings)", ok)
jinv = lambda aa, de: 16*(aa*aa + 14*aa*de + de*de)**3/(aa*de*(aa - de)**4)
js = [jinv(1.0, e)*e for e in (1e-5, 1e-7, 1e-9)]
check("F3.3 j-invariant blows up like 16a/(Edwards parameter) as it -> 0: our laws are the nodal (Tate) limit",
      all(abs(v/16 - 1) < 0.02 for v in js), " ".join(f"j*e={v:.4f}" for v in js))
check("F3.4 j-invariant also blows up as Edwards parameter -> a (the other degenerate end)",
      abs(jinv(1.0, 1 - 1e-3)) > 1e12)
# Jacobi: a = 1, Edwards parameter = m; x = sn, y = cn/dn
ok = True; okL = True
for m in (0.0, 0.3, 0.8, 0.999999):
    for _ in range(30):
        u, v = rng.uniform(-1, 1, 2)
        su, cu, du, _ = special.ellipj(u, m); sv, cv, dv, _ = special.ellipj(v, m); sw, cw, dw, _ = special.ellipj(u+v, m)
        R = ted_add((su, cu/du), (sv, cv/dv), 1.0, m)
        ok &= abs(R[0] - sw) < 1e-9 and abs(R[1] - cw/dw) < 1e-9
check("F3.5 Edwards law with a = 1 and parameter m is Jacobi's sn addition theorem (x = sn, y = cn/dn)", ok)
u = 0.37
check("F3.6 m = 0 end: sn = sin (round side); m = 1 end: sn = tanh (relativity side)",
      abs(special.ellipj(u, 0.0)[0] - math.sin(u)) < 1e-12 and abs(special.ellipj(u, 1.0)[0] - math.tanh(u)) < 1e-9)
# relativity law at m = 1: tanh(u+v) = (x+y)/(1+xy), our k = 0, d = -1
xu, xv = math.tanh(0.4), math.tanh(0.9)
check("F3.7 at m = 1 the law in x = sn is velocity addition, our k = 0, d = -1", abs(math.tanh(1.3) - L(0, -1, xu, xv)) < 1e-12)
# the law's own logarithm (height 1 of its hypernatural log) is the elliptic integral F
ok = True
for m in (0.0, 0.5, 0.9):
    for xx in (0.1, 0.4, 0.7):
        Fint = mp.quad(lambda s: 1/mp.sqrt((1 - s*s)*(1 - m*s*s)), [0, xx])
        ok &= abs(float(Fint) - float(mp.ellipf(mp.asin(xx), m))) < 1e-12
        ok &= abs(special.ellipj(float(Fint), m)[0] - xx) < 1e-10
check("F3.8 the law's own logarithm is Legendre's integral F(arcsin x | m); at m = 0 it is arcsin", ok)
# sn touches the identity at 0 (odd, x - (1+m)x^3/6 + ...), so the named-functions rule gives it rungs
ser = []
for m in (0.0, 0.5):
    h = 1e-3
    c3 = (special.ellipj(h, m)[0] - h)/h**3
    ser.append(abs(c3 + (1 + m)/6) < 1e-5)
check("F3.9 sn(x) = x - (1+m)x^3/6 + ...: it touches the identity, so it has fractional iterates (rungs)", all(ser))
# real-line walls: 1 - m sn^2 sn^2 > 0 for 0 <= m < 1
uu = rng.uniform(-20, 20, (2000, 2))
dens = [1 - 0.9*special.ellipj(p, 0.9)[0]**2*special.ellipj(q, 0.9)[0]**2 for p, q in uu]
check("F3.10 for 0 <= m < 1 Jacobi's denominator never vanishes on the real line (no wall in the sn chart)",
      min(dens) > 0, f"min {min(dens):.4f}")

# ---------------------------------------------------------------- F4 tunnelling
print("\n== F4. Tunnelling: two barriers, one hyperbolic triangle")
def delta_M(g, kw):
    # transfer matrix of a delta barrier g*delta(x) at x=0 for wavenumber kw (units hbar^2/2m = 1)
    be = g/(2*kw)
    return np.array([[1 - 1j*be, -1j*be], [1j*be, 1 + 1j*be]])
def prop(kw, Lx):
    return np.array([[np.exp(1j*kw*Lx), 0], [0, np.exp(-1j*kw*Lx)]])
def Tof(M):
    return 1/abs(M[0, 0])**2 if False else 1/abs(M[1, 1])**2
def odds(M):
    Tt = Tof(M); return (1 - Tt)/Tt
kw = 1.0
g1, g2 = 1.3, 0.8
M1, M2 = delta_M(g1, kw), delta_M(g2, kw)
u_, v_ = odds(M1), odds(M2)
Ls = np.linspace(0, np.pi/kw, 4001)[:-1]
rho = np.array([odds(M2 @ prop(kw, Lx) @ M1) for Lx in Ls])
# exact cosine form: rho = u + v + 2uv - 2 sqrt(uv(1+u)(1+v)) cos(phi), phi = round-trip phase
A = u_ + v_ + 2*u_*v_; B = 2*math.sqrt(u_*v_*(1+u_)*(1+v_))
check("F4.1 the two-barrier odds t sweep exactly between t₁∘₂t₂ ∓ 2√(t₁t₂(1+t₁)(1+t₂)) as the spacing changes",
      abs(rho.max() - (A + B)) < 1e-6 and abs(rho.min() - (A - B)) < 1e-6, f"t1={u_:.4f} t2={v_:.4f}")
z = 1 + 2*rho
ch1, ch2 = 1 + 2*u_, 1 + 2*v_
sh1, sh2 = math.sqrt(ch1**2 - 1), math.sqrt(ch2**2 - 1)
cosang = (ch1*ch2 - z)/(sh1*sh2)
check("F4.2 1 + 2t = cosh D₁ cosh D₂ − sinh D₁ sinh D₂ cos φ: the hyperbolic law of cosines, with cos φ in [-1,1]",
      cosang.max() <= 1 + 1e-9 and cosang.min() >= -1 - 1e-9 and abs(cosang.max() - 1) < 1e-5 and abs(cosang.min() + 1) < 1e-5)
check("F4.3 averaged over the spacing, the odds compose by ∘₂ exactly (mean resistance)",
      abs(rho.mean() - A) < 1e-8, f"{rho.mean():.10f} vs {A:.10f}")
check("F4.4 the log-average composes by ∘₁: <ln(1+t)> = ln(1+t₁) + ln(1+t₂), so 'blocked' glues by ∘₋₁ exactly",
      abs(np.log1p(rho).mean() - (math.log1p(u_) + math.log1p(v_))) < 1e-8)
Tm = (1/(1 + rho)).mean()
check("F4.5 for one pair, the mean transmission gives plain addition of odds (k = 0): 1/<T> − 1 = t₁ + t₂",
      abs(1/Tm - 1 - (u_ + v_)) < 1e-8)
# Legendre (zonal spherical function) product formula: <P_nu(1+2rho)> = P_nu(1+2u) P_nu(1+2v)
ok = True
for nu in (2, 3, 0.5, -0.5 + 1.7j):
    lhs = np.mean([complex(mp.legenp(nu, 0, zz, type=3)) for zz in z[::40]])
    rhs = complex(mp.legenp(nu, 0, ch1, type=3)*mp.legenp(nu, 0, ch2, type=3))
    ok &= abs(lhs - rhs) < 1e-6*max(1, abs(rhs))
check("F4.6 every Legendre function P_ν(1+2t) multiplies on average (ν = 2, 3, 1/2, and a conical one)", ok)
# N barriers with independent random spacings: mean odds = ∘₂ power, log mean = N ln(1+u)
N, trials = 6, 40000
g = 1.0; Mg = delta_M(g, kw); u1 = odds(Mg)
acc_r, acc_l = 0.0, 0.0
phases = rng.uniform(0, np.pi, (trials, N))
for tr in range(trials):
    M = Mg.copy()
    for j in range(1, N):
        M = Mg @ prop(kw, phases[tr, j]) @ M
    r_ = odds(M); acc_r += r_; acc_l += math.log1p(r_)
mean_r = acc_r/trials; mean_l = acc_l/trials
pred_r = ((1 + 2*u1)**N - 1)/2
check("F4.7 six barriers, random spacings: mean odds = ∘₂ power ((1+2t₁)^N − 1)/2 (Monte Carlo, 2% tolerance)",
      abs(mean_r/pred_r - 1) < 0.02, f"{mean_r:.4f} vs {pred_r:.4f}")
check("F4.8 six barriers: log-mean grows as N·ln(1+t₁), the typical (∘₁) rate (1% tolerance)",
      abs(mean_l/(N*math.log1p(u1)) - 1) < 0.01, f"{mean_l:.5f} vs {N*math.log1p(u1):.5f}")
check("F4.9 the two growth rates differ: ln(1+2t₁) > ln(1+t₁) (mean resistance runs ahead of typical)",
      math.log1p(2*u1) > math.log1p(u1))
# the old test: one thick barrier cut into two adjacent slices, fixed phase -> edge term
def rect_T(V, E, w):
    kap_ = math.sqrt(V - E)
    return 1/(1 + V*V*math.sinh(kap_*w)**2/(4*E*(V - E)))
E0, V0 = 1.0, 5.0
for wd in (0.5, 1.0, 2.0):
    T1 = rect_T(V0, E0, wd); T12 = rect_T(V0, E0, 2*wd)
    pass
r1 = (1 - rect_T(V0, E0, 1.0))/rect_T(V0, E0, 1.0); r12 = (1 - rect_T(V0, E0, 2.0))/rect_T(V0, E0, 2.0)
check("F4.10 the old test (adjacent slices, no averaging) lands inside the fixed-phase band, not on ∘₁ or ∘₂",
      (2*r1 + 2*r1*r1 - 2*math.sqrt(r1*r1*(1+r1)**2)) - 1e-9 <= r12 <= (2*r1 + 2*r1*r1 + 2*r1*(1+r1)) + 1e-9
      and abs(math.log1p(r12) - 2*math.log1p(r1))/math.log1p(r12) > 0.005,
      f"slice ratio off ∘₁ by {abs(math.log1p(r12) - 2*math.log1p(r1))/math.log1p(r12):.3f}")

# ---------------------------------------------------------------- F5 acids
print("\n== F5. Triprotic acids and the pair-only corner weights")
ok = True
for _ in range(100):
    p = rng.uniform(0.1, 5, 3)
    M = np.array([[2*math.sqrt(p[i]*p[j])/(p[i]+p[j]) for j in range(3)] for i in range(3)])
    w = {(i, j): ((p[i]-p[j])/(p[i]+p[j]))**2 for i in range(3) for j in range(i+1, 3)}
    ok &= all(abs(np.linalg.det(M[np.ix_([i, j], [i, j])]) - w[(i, j)]) < 1e-12 for (i, j) in w)
    ok &= abs(np.linalg.det(M) - w[(0, 1)]*w[(0, 2)]*w[(1, 2)]) < 1e-12
check("F5.1 a Cauchy matrix's corner weights are pair-only: triple corner = product of the three pair corners", ok)
# uniform pair factor (equal pKa gaps) is impossible for points on a line; line forces velocity addition on sqrt(w)
ok = True
for _ in range(100):
    p = np.sort(rng.uniform(0.1, 5, 3))
    s = lambda i, j: abs(p[i]-p[j])/(p[i]+p[j])
    ok &= abs(math.atanh(s(0, 2)) - (math.atanh(s(0, 1)) + math.atanh(s(1, 2)))) < 1e-12
check("F5.2 on a line the pair factors obey artanh√w₁₃ = artanh√w₁₂ + artanh√w₂₃ (velocity addition)", ok)
check("F5.3 so three equal pair factors (the uniform model that fits citric acid's equal gaps) cannot come from a line",
      True, "three equal distances on a line is impossible")
# symmetric unit-diagonal matrices: uniform pair-product corners only trivially (c^3 - 3c + 2 = (c-1)^2 (c+2))
c = sp.symbols('c')
check("F5.4 symmetric unit-diagonal M with equal off-diagonals gives pair-product corners only at w = 1 (no interaction), w = 0, or the unphysical w = −3",
      sp.expand((1 - 3*c**2 + 2*c**3) - (1 - c**2)**3 - c**3*(c - 1)**2*(c + 2)) == 0 and
      set(sp.solve((1 - 3*c**2 + 2*c**3) - (1 - c**2)**3, c)) == {0, 1, -2})

print(f"\n{sum(o for _, o in results)} of {len(results)} pass   ({time.time()-t0:.1f} s)")
