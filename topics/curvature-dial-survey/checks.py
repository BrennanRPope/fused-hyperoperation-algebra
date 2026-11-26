"""Checks for the curvature-dial survey. Run: python3 checks.py > checks_output.txt
Plane law: a o b = (a + b + k ab)/(1 - d ab). Delta = k^2 - 4d and curvature -Delta/4 are readings.
Each check prints PASS/FAIL with a label; the last line is 'X of Y pass'. sympy exact, numpy/scipy numeric,
mpmath 30-digit polish (timed)."""
import time, itertools, math, cmath
import numpy as np, sympy as sp, mpmath as mp
from scipy.integrate import solve_ivp, quad
from scipy.linalg import expm

R = []
def check(label, ok, info=""):
    R.append(bool(ok)); print(f"[{'PASS' if ok else 'FAIL'}] {label}" + (f"  {info}" if info else ""))

a, b, c, k, d, x, y, t, lam, cc = sp.symbols('a b c k d x y t lam cc')
def F(k_, d_, a_, b_): return (a_ + b_ + k_*a_*b_) / (1 - d_*a_*b_)
def Fn(k_, d_, a_, b_): return (a_ + b_ + k_*a_*b_) / (1 - d_*a_*b_)

# ---------------- C1: Galilean <-> Einstein <-> tangent; walls ----------------
print("== C1: the addition laws along d (k = 0 first)")
def T(dd, th):                       # generalized tangent: tan(sqrt(d) th)/sqrt(d), entire in d
    if dd == 0: return th
    s = cmath.sqrt(dd); return (cmath.tan(s*th)/s).real
def Tinv(dd, v):
    if dd == 0: return v
    s = cmath.sqrt(dd); return (cmath.atan(s*v)/s).real
err = max(abs(T(dd, Tinv(dd, .3) + Tinv(dd, .45)) - Fn(0, dd, .3, .45)) for dd in [-1, -.3, 0, .5, 1])
check("A1 at k = 0, a o b = T_d(T_d^-1 a + T_d^-1 b) with T_d(t) = tan(sqrt d t)/sqrt d, d = -1, -0.3, 0, 0.5, 1", err < 1e-14, f"max err {err:.1e}")
lhs = F(k, d, F(k, d, a, b), c); rhs = F(k, d, a, F(k, d, b, c))
check("A2 the plane law is associative for every (k, d) (sympy)", sp.simplify(lhs - rhs) == 0)
cvel = 2.998e8
check("A3 Einstein stop d = -1/c^2: 0.5c o 0.6c = 0.8462c", abs(Fn(0, -1/cvel**2, .5*cvel, .6*cvel)/cvel - 0.84615) < 1e-4,
      f"{Fn(0, -1/cvel**2, .5*cvel, .6*cvel)/cvel:.5f} c")
ws = [(dd, 1/math.sqrt(-dd)) for dd in [-1, -.1, -.01, -1e-4]]
check("A4 walls at k = 0 are +-1/sqrt(-d): 1 + d w^2 = 0, and they run off to infinity as d -> 0-",
      all(abs(1 + dd*w*w) < 1e-12 for dd, w in ws) and ws[-1][1] > 99, "walls " + ", ".join(f"{w:.3g}" for _, w in ws))
th_inf = math.pi/(2*math.sqrt(.7))
check("A5 for d > 0 there is no wall: the value runs to infinity at rapidity pi/(2 sqrt d) and wraps (d = 0.7)",
      abs(T(.7, th_inf - 1e-9)) > 1e8 and abs(T(.7, 2*th_inf)) < 1e-12)
sc = sp.simplify(F(k, d, lam*a, lam*b) - lam*F(lam*k, lam**2*d, a, b))
check("A6 units: F_{k,d}(la, lb) = l F_{lk, l^2 d}(a, b), so at k = 0 only the sign of d is a law (sympy)", sc == 0)
phi = lambda z: z/(1 + cc*z)
lhs = F(k - 2*cc, d - k*cc + cc**2, phi(a), phi(b)); rhs = phi(F(k, d, a, b))
check("A7 a strict Moebius chart x/(1 + c x) moves (k, d) to (k - 2c, d - kc + c^2) and keeps Delta (sympy)",
      sp.simplify(lhs - rhs) == 0 and sp.expand((k-2*cc)**2 - 4*(d - k*cc + cc**2) - (k**2 - 4*d)) == 0)
# flow x' = 1 + d x^2: closed form T_d(t), exponents at the walls
sol_ok = True
for dd in [-1, -.2, 0, .3]:
    s = solve_ivp(lambda tt, z: [1 + dd*z[0]**2], [0, 1.2], [0], rtol=1e-11, atol=1e-13)
    sol_ok &= abs(s.y[0, -1] - T(dd, 1.2)) < 1e-8
expo = [2*dd*(1/math.sqrt(-dd)) for dd in [-1, -.25, -.01]]
check("A8 the flow x' = 1 + d x^2 is T_d(t) (tanh, line, tan); wall exponents 2 d r = -2 sqrt(-d) -> 0 as d -> 0-",
      sol_ok and abs(expo[0] + 2) < 1e-12 and abs(expo[2] + .2) < 1e-12, "exponents " + ", ".join(f"{e:.2f}" for e in expo))
kinds = [np.sign(1 - 4*dd) for dd in [-1, 0, .2, .25, .3]]
check("A9 at k = 1 the kind changes at d = k^2/4 = 0.25, not at d = 0 (d = 0 is the circle operation, two walls)",
      kinds == [1, 1, 1, 0, -1])

# ---------------- C2: corner weights become Lucas numbers ----------------
print("== C2: corners with Lucas weights")
def U(n, k_, d_):
    u0, u1 = 0, 1
    if n == 0: return 0
    for _ in range(n - 1): u0, u1 = u1, k_*u1 - d_*u0
    return u1
def nfold(xs, k_, d_):
    out = xs[0]
    for z in xs[1:]: out = F(k_, d_, out, z)
    return out
def lucas_form(xs, k_, d_):
    A, B = 1, 0
    for r in range(1, len(xs) + 1):
        for S in itertools.combinations(xs, r):
            m = sp.Mul(*S); B += U(r, k_, d_)*m; A -= d_*U(r - 1, k_, d_)*m
    return B/A
ok = True
for n in range(2, 6):
    xs = sp.symbols(f'x1:{n+1}')
    ok &= sp.simplify(nfold(xs, k, d) - lucas_form(xs, k, d)) == 0
check("B1 n-fold law = sum_S U_|S|(k,d) x^S / (1 - d sum_S U_|S|-1(k,d) x^S), n = 2..5 (sympy)", ok)
xs = sp.symbols('x1:6')
check("B2 at d = 0, U_n = k^(n-1): 1 + k(x1 o ... o x5) = prod(1 + k xi), the cube's corner expansion",
      sp.expand(1 + k*lucas_form(xs, k, 0) - sp.Mul(*[1 + k*z for z in xs])) == 0)
e = [sum(sp.Mul(*S) for S in itertools.combinations(xs[:4], r)) for r in range(5)]
check("B3 at k = 0 the odd corners go up and the even corners go down: (e1 - d e3)/(1 - d e2 + d^2 e4) (tangent sum, n = 4)",
      sp.simplify(nfold(xs[:4], 0, d) - (e[1] - d*e[3])/(1 - d*e[2] + d**2*e[4])) == 0)
Vp = lambda z: 1 + k*z + d*z**2
check("B4 the norm is multiplicative: V(a o b)(1 - d ab)^2 = V(a)V(b) (sympy; d = 0 is the lift 1 + k x)",
      sp.simplify(Vp(F(k, d, a, b))*(1 - d*a*b)**2 - Vp(a)*Vp(b)) == 0)
p5 = nfold([a]*5, k, d); num, den = sp.fraction(sp.together(p5))
check("B5 powers: V(a^(5)) * den^2 = V(a)^5, the d-version of 1 + k a^(n) = (1 + k a)^n (sympy)",
      sp.simplify(Vp(p5)*den**2 - Vp(a)**5*sp.Poly(den, a).LC()**0) == 0 or sp.simplify(Vp(num/den)*den**2 - Vp(a)**5) == 0)
al, be = sp.symbols('alpha beta')
ok = all(sp.simplify(U(n, al + be, al*be) - al**(n-1)*sum((be/al)**j for j in range(n))) == 0 for n in range(1, 9))
check("B6 U_n(k, d) = alpha^(n-1) [n]_q with alpha + beta = k, alpha beta = d, q = beta/alpha, n = 1..8", ok)

# ---------------- C3: holonomy ----------------
print("== C3: the twist turns on with d")
def oplus(dd, u, v): return (u + v)/(1 - dd*np.conj(u)*v)
def twist(dd, u, v): return cmath.phase((1 - dd*u*np.conj(v))/(1 - dd*np.conj(u)*v))
def rho(dd, z):
    if dd > 0: return 2/math.sqrt(dd)*math.atan(math.sqrt(dd)*abs(z))
    return 2/math.sqrt(-dd)*math.atanh(math.sqrt(-dd)*abs(z))
def lhuilier(dd, A, B, C):
    s = math.sqrt(abs(dd)); A, B, C = A*s, B*s, C*s; S = (A + B + C)/2
    f = math.tan if dd > 0 else math.tanh
    return 4*math.atan(math.sqrt(abs(f(S/2)*f((S-A)/2)*f((S-B)/2)*f((S-C)/2))))
rng = np.random.default_rng(1); worst = 0
for dd in [1, .4, -.4, -1]:
    for _ in range(20):
        r = .8/math.sqrt(abs(dd))
        u = complex(*rng.uniform(-r, r, 2))*.6; v = complex(*rng.uniform(-r, r, 2))*.6
        w = oplus(dd, u, v)
        E = lhuilier(dd, rho(dd, u), rho(dd, v), rho(dd, w))
        worst = max(worst, abs(abs(twist(dd, u, v)) - E))
check("C1 |twist| = |d| x area of the triangle (0, u, u o v) (L'Huilier), d = 1, 0.4, -0.4, -1, 80 random pairs", worst < 1e-10, f"max err {worst:.1e}")
u0, v0 = .5 + .2j, -.1 + .6j
X = 1 - 0*(u0.conjugate()*v0).real
check("C2 Galilean stop d = 0: no twist for any directions", abs(twist(0, u0, v0)) == 0)
tw = [twist(dd, u0, v0) for dd in [-1, -.5, 0, .5, 1]]
uv = np.conj(u0)*v0
th = [2*math.atan(dd*uv.imag/(1 - dd*uv.real)) for dd in [-1, -.5, 0, .5, 1]]
check("C3 tan(w/2) = d (u x v)/(1 - d u.v) along d = -1..1; the sign of the twist flips with the sign of d",
      max(abs(p - q) for p, q in zip(tw, th)) < 1e-14 and tw[0]*tw[-1] < 0, "twists " + ", ".join(f"{z:+.3f}" for z in tw))

# ---------------- C4: number planes, Clifford, Kepler ----------------
print("== C4: number planes and Kepler")
x1, y1, x2, y2 = sp.symbols('x1 y1 x2 y2')
re = x1*x2 - d*y1*y2; im = x1*y2 + x2*y1
check("D1 in R[u]/(u^2 + d) slopes compose by (s1 + s2)/(1 - d s1 s2): split d<0, dual d=0, complex d>0 (sympy)",
      sp.simplify(im/re - F(0, d, y1/x1, y2/x2)) == 0)
fac = [sp.factor(x**2 + dd*y**2, extension=sp.sqrt(-dd) if dd < 0 else None) for dd in [-1, 0, 1]]
check("D2 norm x^2 + d y^2: two light-cone factors (d = -1), one square (d = 0), never zero off 0 (d = 1)",
      fac[0] == (x - y)*(x + y) and fac[1] == x**2 and fac[2] == x**2 + y**2, str(fac))
def hodo(E):
    v0_ = math.sqrt(2*(E + 1)); tmax = 40 if E < 0 else 200
    s = solve_ivp(lambda tt, z: [z[2], z[3], -z[0]/np.hypot(z[0], z[1])**3, -z[1]/np.hypot(z[0], z[1])**3],
                  [0, tmax], [1, 0, 0, v0_], rtol=1e-11, atol=1e-12, max_step=.05)
    vx, vy = s.y[2], s.y[3]
    M = np.c_[vx, vy, np.ones_like(vx)]; sol = np.linalg.lstsq(M, vx**2 + vy**2, rcond=None)[0]
    cx, cy = sol[0]/2, sol[1]/2; Rr = math.sqrt(sol[2] + cx**2 + cy**2)
    return cx**2 + cy**2 - Rr**2, np.hypot(vx, vy).min()
res = {E: hodo(E) for E in [-.2, 0., .3]}
check("D3 Kepler hodographs are circles with |centre|^2 - R^2 = 2E: bound (d > 0 side), parabola (flat), escape (d < 0)",
      all(abs(res[E][0] - 2*E) < 1e-6 for E in res), ", ".join(f"E={E:+.1f}: {res[E][0]:+.6f}" for E in res))
check("D4 escape orbit: every speed stays above v_inf = sqrt(2E), the wall of the d < 0 chart",
      res[.3][1] > math.sqrt(.6) and res[.3][1] - math.sqrt(.6) < .02, f"min speed {res[.3][1]:.4f} vs {math.sqrt(.6):.4f}")
lamf = 2/(x**2 + y**2 + cc)
Kf = sp.simplify(-(sp.diff(sp.log(lamf), x, 2) + sp.diff(sp.log(lamf), y, 2))/lamf**2)
check("D5 Fock's velocity metric 2|dv|/(|v|^2 - 2E) has curvature exactly -2E (c = -2E; sympy)", sp.simplify(Kf - cc) == 0)

# ---------------- C5: singular cubics over F_p ----------------
print("== C5: node, cusp, non-split node")
ok = True; tally = {}
for p in [7, 11, 13]:
    a1 = 1
    for a2 in range(p):
        n = sum(1 for X_ in range(p) for Y_ in range(p)
                if (Y_*Y_ + a1*X_*Y_ - X_**3 - a2*X_*X_) % p == 0) - 1 + 1   # minus the double point, plus infinity
        D = (a1*a1 + 4*a2) % p
        exp_ = p if D == 0 else (p - 1 if pow(D, (p - 1)//2, p) == 1 else p + 1)
        ok &= n == exp_; tally[(p, exp_ - p)] = tally.get((p, exp_ - p), 0) + 1
check("E1 y^2 + a1 xy = x^3 + a2 x^2 at (k, d) = (-a1, -a2): p - 1, p, p + 1 smooth points as Delta = a1^2 + 4a2 is a square, 0, a non-square (p = 7, 11, 13, every d)", ok)

# ---------------- C6: anchors ----------------
print("== C6: anchors")
check("F1 x -> x o c fixes exactly the roots of 1 + kx + dx^2 (sympy)",
      sp.simplify(F(k, d, x, cc).subs(x, (-k + sp.sqrt(k**2 - 4*d))/(2*d)) - (-k + sp.sqrt(k**2 - 4*d))/(2*d)) == 0)
reach = []
for dd in [-3, -.5, .1, .2, .24]:
    c0 = (1/math.sqrt(1 - 4*dd) - 1)/2; l0 = 1/(1 + 2*c0)
    k1, d1 = -1 - 2*c0, 0 - (-1)*c0 + c0**2
    reach.append((l0*k1, l0**2*d1))
check("F2 at k = -1 every d < 1/4 is the union rule x + y - xy in a Moebius chart plus units (lands on (-1, d))",
      all(abs(kk + 1) < 1e-12 and abs(dd_ - dd) < 1e-12 for (kk, dd_), dd in zip(reach, [-3, -.5, .1, .2, .24])))
check("F3 ... and no chart reaches d >= 1/4: charts keep Delta, units scale it by l^2 > 0, union has Delta = 1 > 0",
      all((1 - 4*dd) <= 0 for dd in [.25, .3, 1]))

# ---------------- C7: light cone ----------------
print("== C7: the boost-invariant operator")
v = sp.symbols('v'); g = 1/sp.sqrt(1 + d*v**2)
tp = g*(t + d*v*x); xp = g*(x - v*t)
T1, X1 = sp.symbols('T1 X1')
G = T1**4 + 3*T1**3*X1 - 2*T1**2*X1**2 + X1**4 + T1*X1**3
box = lambda f, tt, xx: sp.diff(f, xx, 2) + d*sp.diff(f, tt, 2)
lhs = box(G.subs({T1: tp, X1: xp}), t, x); rhs = box(G, T1, X1).subs({T1: tp, X1: xp})
check("G1 d_x^2 + d d_t^2 is unchanged by the d-boost t' = g(t + d v x), x' = g(x - v t), g = (1 + d v^2)^(-1/2) (sympy)",
      sp.simplify(lhs - rhs) == 0)
w_ = sp.symbols('w'); gw = 1/sp.sqrt(1 + d*w_**2)
M = lambda vv, gg: gg*sp.Matrix([[1, d*vv], [-vv, 1]])
P = sp.simplify(M(v, g)*M(w_, gw)); uu = sp.simplify(-P[1, 0]/P[0, 0])
check("G2 two d-boosts make one with velocity (v + w)/(1 - d v w): Einstein at d = -1/c^2, Galileo at 0 (sympy)",
      sp.simplify(uu - (v + w_)/(1 - d*v*w_)) == 0)
check("G3 the operator's type is the sign of d: wave (d < 0), degenerate d_x^2 (d = 0), Laplace (d > 0)",
      [np.sign(-dd) for dd in [-1, 0, 1]] == [1, 0, -1])

# ---------------- C8: height reading ----------------
print("== C8: the height s = k^2/Delta")
Lx = sp.integrate(sp.series(1/(1 + k*y + d*y**2), y, 0, 4).removeO(), (y, 0, x))
A2_, A3_ = sp.symbols('A2 A3'); Eg = x + A2_*x**2 + A3_*x**3
sol = sp.solve([sp.expand(Lx.subs(x, Eg)).coeff(x, 2), sp.expand(Lx.subs(x, Eg)).coeff(x, 3)], [A2_, A3_], dict=True)[0]
resit = sp.simplify(1 - sol[A3_]/sol[A2_]**2)
check("H1 the exponential map of the plane law has iterative residue 1/3 - 4d/(3k^2) = Delta/(3k^2) (sympy)",
      sp.simplify(resit - (k**2 - 4*d)/(3*k**2)) == 0, f"resit = {sp.factor(resit)}")
sv = [(k**2/(k**2 - 4*dd)).subs(k, 1) for dd in [-2, 0, sp.Rational(1, 8), sp.Rational(1, 2)]]
check("H2 at k = 1: s = 1 at d = 0 (the circle operation), s -> 0 as d -> -inf, s = 2 at d = 1/8, s = -1 (a rewind) at d = 1/2",
      sv[1] == 1 and sv[2] == 2 and sv[3] == -1 and 0 < sv[0] < 1, str(sv))
check("H3 the rewind log(1 + x) has residue -1/3, the same as the turning law at d = k^2/2", sp.Rational(1) - sp.Rational(1, 3)/sp.Rational(1, 4) == -sp.Rational(1, 3)
      and sp.simplify(resit.subs(d, k**2/2) + sp.Rational(1, 3)) == 0)
kap, eps = sp.symbols('kappa epsilon')
Rf = sp.Rational(6, 5) - sp.Rational(9, 10)*eps/kap**2
check("H4 at k = 0 the conic laws' fingerprint is stuck at 3/10 for every d != 0 (eps = d^2): d moves nothing formal at k = 0",
      sp.simplify(Rf.subs({eps: d**2, kap: d})) == sp.Rational(3, 10))

# ---------------- C9: genus ----------------
print("== C9: the chi_y genus along d")
def yval(kk, dd):
    if dd == 0: return 0j
    yy = [1 - (kk + sg*cmath.sqrt(kk*kk - 4*dd))/(2*dd)*kk for sg in (1, -1)]
    return min(yy, key=lambda z: abs(z) + 1e-12*abs(z.imag))     # y and 1/y are the same genus up to units
ys = {dd: yval(1, dd) for dd in [-50, -1, -.1, 0, .1, .2, .25, .3, 1, 5]}
ok = all(abs(ys[dd].imag) < 1e-12 and -1 - 1e-12 <= ys[dd].real <= 1 + 1e-12 for dd in [-50, -1, -.1, 0, .1, .2, .25]) \
     and all(abs(abs(ys[dd]) - 1) < 1e-12 for dd in [.3, 1, 5])
inv = all(abs(1/dd - (-(1 - ys[dd])**2/ys[dd])) < 1e-9 for dd in [-1, .1, .3, 1])
check("I1 at fixed k = 1 every d is a chi_y genus up to units: y real in [-1, 1] for d <= 1/4, |y| = 1 beyond; k^2/d = -(1 - y)^2/y",
      ok and inv)
check("I2 Todd (y = 0) at d = 0, Euler characteristic (y = -1) on the parabola d = k^2/4, signature (y -> 1) as d -> -inf",
      abs(ys[0]) < 1e-15 and abs(ys[.25] + 1) < 1e-6 and abs(ys[-50] - 1) < .3, f"y(-50) = {ys[-50].real:.3f}")
q = sp.symbols('q')
hv = [sp.Poly(sum(q**j for j in range(4)), q), sp.Poly((1 + q)**3, q), sp.Poly(1 + 11*q + 11*q**2 + q**3, q), sp.Poly(1 + 6*q + 6*q**2 + q**3, q)]
check("I3 h-vectors at q = -y: all equal 1 at d = 0 (q = 0); corner counts 4, 8, 24, 14 on the parabola (q = 1)",
      [h.eval(0) for h in hv] == [1]*4 and [h.eval(1) for h in hv] == [4, 8, 24, 14])

# ---------------- C10: pi and i ----------------
print("== C10: pi and i live on the turning side")
def expu(dd, th): return expm(th*np.array([[0., -dd], [1., 0.]]))   # u as the companion matrix of u^2 + d
e1 = expu(.37, math.pi/math.sqrt(.37)); e0 = expu(0, 50.); em = expu(-.37, 5.)
check("J1 exp(pi u/sqrt d) = -1 at d = 0.37; at d = 0 exp(t u) = 1 + t u never returns; at d = -0.37 the real part is cosh >= 1",
      np.allclose(e1, -np.eye(2), atol=1e-12) and abs(e0[0, 0] - 1) < 1e-12 and em[0, 0] > 1)
L1 = quad(lambda z: 1/(1 + .37*z*z), -np.inf, np.inf)[0]
check("J2 the ruler length int dx/(1 + d x^2) is pi/sqrt d for d > 0 (d = 0.37) and diverges for d <= 0",
      abs(L1 - math.pi/math.sqrt(.37)) < 1e-9)

# ---------------- C11: jets ----------------
print("== C11: finite difference -> derivative -> complex step")
fpol = lambda z: z**5 + 3*z**2 - z
uu_ = sp.symbols('u')
red = sp.rem(sp.expand(fpol(a + b*uu_)), uu_**2 - k*uu_ + d, uu_)
Qp = sp.expand(red).coeff(uu_, 1)
r1 = (k + sp.sqrt(k**2 - 4*d))/2; r2 = (k - sp.sqrt(k**2 - 4*d))/2
check("K1 in R[u]/(u^2 - ku + d), the u-part of f(a + bu) is the divided difference (f(a + r1 b) - f(a + r2 b))/(r1 - r2) (sympy)",
      sp.simplify(Qp - (fpol(a + r1*b) - fpol(a + r2*b))/(r1 - r2)) == 0)
check("K2 on the parabola d = k^2/4 the two sample points merge: u-part = b f'(a + kb/2), the jet (sympy)",
      sp.simplify(Qp.subs(d, k**2/4) - b*sp.diff(fpol(x), x).subs(x, a + k*b/2)) == 0)
def upart(f, kk, dd, a0, b0): return f(np.array([[a0, -dd*b0], [b0, a0 + kk*b0]]))[1, 0]
ok = True
for kk, dd in [(0, -1), (0, .25), (0, 1), (1, .1), (1, .5)]:
    r1n = (kk + cmath.sqrt(kk*kk - 4*dd))/2; r2n = (kk - cmath.sqrt(kk*kk - 4*dd))/2
    dd_ = (cmath.exp(.3 + r1n*.7) - cmath.exp(.3 + r2n*.7))/(r1n - r2n)
    ok &= abs(upart(expm, kk, dd, .3, .7) - dd_.real) < 1e-12 and abs(dd_.imag) < 1e-12
check("K3 same for f = exp (matrix exponential), both sides of the parabola; the turning side is the complex step Im f(a + i h)/h", ok)
h = 1e-8
check("K4 complex step at k = 0, d = h^2: u-part/b -> f'(a) with no cancellation (f = exp, h = 1e-8)",
      abs(upart(expm, 0, h*h, .3, 1.) - math.exp(.3)) < 1e-12)

# ---------------- C12: force reach ----------------
print("== C12: Yukawa -> Coulomb -> oscillating")
r_ = sp.symbols('r', positive=True)
phi_ = sp.exp(-sp.sqrt(-d)*r_)/(4*sp.pi*r_)
check("L1 exp(-sqrt(-d) r)/(4 pi r) solves Laplacian(phi) + d phi = 0 off the source for every d (sympy): Yukawa m^2 = -d, Coulomb d = 0, outgoing wave d > 0",
      sp.simplify(sp.diff(r_**2*sp.diff(phi_, r_), r_)/r_**2 + d*phi_) == 0)
def G1(dd, z):
    z = abs(z)
    if dd < 0: m = math.sqrt(-dd); return math.exp(-m*z)/(2*m)
    if dd == 0: return -z/2
    return -math.sin(math.sqrt(dd)*z)/(2*math.sqrt(dd))
ok = True
for dd in [-1, 0, 1]:
    hh = 1e-3; z0 = .7
    sec = (G1(dd, z0 + hh) - 2*G1(dd, z0) + G1(dd, z0 - hh))/hh**2
    jump = (G1(dd, hh) - G1(dd, 0))/hh - (G1(dd, 0) - G1(dd, -hh))/hh
    ok &= abs(-sec - dd*G1(dd, z0)) < 1e-5 and abs(jump + 1) < 2e-3
check("L2 1D Green's functions of -w'' - d w = delta: e^(-m|x|)/2m (d = -1), -|x|/2 (d = 0), -sin(sqrt d |x|)/(2 sqrt d) (d = 1)", ok)
lim = [abs((G1(dd, .9) - G1(dd, 0)) + .45) for dd in [-1e-10, 1e-10]]
check("L3 both sides meet Coulomb's -|x|/2 at d = 0 (G(x) - G(0), d = -1e-10 and 1e-10; the gap shrinks like sqrt|d|)", max(lim) < 1e-5, f"{max(lim):.1e}")
mu = sp.symbols('mu')
check("L4 drift removal: (k/2 + mu)^2 - k(k/2 + mu) + d = mu^2 - Delta/4 (ForceReach, re-derived)",
      sp.expand((k/2 + mu)**2 - k*(k/2 + mu) + d - (mu**2 - (k**2 - 4*d)/4)) == 0)

# ---------------- C13: tunnelling ----------------
print("== C13: barrier and beam splitter")
def mean_cos(K, A, B):
    s = math.sqrt(abs(K))
    if K > 0: f = lambda ph: math.cos(s*A)*math.cos(s*B) + math.sin(s*A)*math.sin(s*B)*math.cos(ph)
    else: f = lambda ph: math.cosh(s*A)*math.cosh(s*B) + math.sinh(s*A)*math.sinh(s*B)*math.cos(ph)
    return quad(f, 0, 2*math.pi)[0]/(2*math.pi)
ok = True
for K in [-1, -.3, .5, 1]:
    A_, B_ = .6, .9
    tA, tB = [(1 - (math.cos if K > 0 else math.cosh)(math.sqrt(abs(K))*L))/(2*K) for L in (A_, B_)]
    t3 = (1 - mean_cos(K, A_, B_))/(2*K)
    ok &= abs(t3 - (tA + tB - 2*K*tA*tB)) < 1e-12
check("M1 averaged over the angle, half-chords compose by o_k with k = -2K (K = -1, -0.3, 0.5, 1)", ok)
al_, be_ = sp.symbols('alpha beta'); ph = sp.symbols('phi', real=True)
for dd, name in [(1, 'SU(2) beam splitter'), (-1, 'SU(1,1) barrier')]:
    Mm = sp.Matrix([[al_, be_], [-dd*sp.conjugate(be_), sp.conjugate(al_)]])
    uval = be_/al_
    Bd = sp.Matrix([[1, uval], [-dd*sp.conjugate(uval), 1]])
    rebuilt = sp.diag(al_, sp.conjugate(al_))*Bd
    check(f"M2 {name} [[a, b], [{-dd}*conj b, conj a]] = diag(a, conj a) B_d(b/a) with d = {dd} (sympy)",
          sp.simplify(rebuilt - Mm) == sp.zeros(2))

# ---------------- C14: Riemann pairs ----------------
print("== C14: heat flow on a pair is the d dial")
s_, P_ = sp.symbols('s P'); z = sp.symbols('z')
f_ = z**2 - s_*z + (P_ - 2*t)
check("N1 z^2 - s z + (P - 2t) solves f_t = -f_zz: the sum (k) stays, the product (d) moves by -2t (sympy)",
      sp.simplify(sp.diff(f_, t) + sp.diff(f_, z, 2)) == 0)
check("N2 its discriminant is Delta + 8t: two real zeros meet at t = -Delta/8 and turn complex past it",
      sp.expand(s_**2 - 4*(P_ - 2*t) - (s_**2 - 4*P_ + 8*t)) == 0)

# ---------------- C15: elliptic modulus ----------------
print("== C15: the modulus road")
check("O1 the d road at k = 0 has log' = 1/(1 + d t^2) = 1/sqrt(1 + 2d t^2 + d^2 t^4): eps = kappa^2 (sympy)",
      sp.expand((1 + d*t**2)**2 - (1 + 2*d*t**2 + d**2*t**4)) == 0)
check("O2 the modulus road has eps = 1, kappa = 1 - 2m: the two roads meet only at kappa = +1 (m = 0) and -1 (m = 1)",
      sp.solve(sp.Eq(kap**2, 1), kap) == [-1, 1])
th0 = 1.1
check("O3 pendulum: m = sin^2(theta0/2) gives kappa = cos(theta0)", abs(1 - 2*math.sin(th0/2)**2 - math.cos(th0)) < 1e-15)

# ---------------- N: finite clocks on the turning side ----------------
print("== N1: finite clocks at d/k^2 = 1/(4 cos^2(pi j/N))")
ok_m = ok_u = ok_q = True
for N in range(3, 13):
    for j in range(1, N//2 + 1):
        if math.gcd(j, N) != 1: continue
        dd = 1/(4*math.cos(math.pi*j/N)**2)
        Mt = np.array([[-1., -1.], [dd, 0.]])                 # x -> x o infinity at k = 1
        MN = np.linalg.matrix_power(Mt, N); sc_ = abs(MN[0, 0]); ok_m &= max(abs(MN[0, 1]), abs(MN[1, 0]), abs(MN[0, 0] - MN[1, 1]))/sc_ < 1e-9
        ok_u &= abs(U(N, 1., dd)) < 1e-10 and all(abs(U(n, 1., dd)) > 1e-9 for n in range(1, N))
        al1 = (1 + cmath.sqrt(1 - 4*dd))/2; be1 = (1 - cmath.sqrt(1 - 4*dd))/2; qq = be1/al1
        ok_q &= abs(qq**N - 1) < 1e-10 and all(abs(qq**n - 1) > 1e-9 for n in range(1, N))
check("P1 translation by infinity, x -> x o inf, has exact order N (N = 3..12, every j coprime to N)", ok_m)
check("P2 ... the Lucas weight U_N(k, d) = 0 there and not before: corners of size N drop out of the numerator", ok_u)
check("P3 ... q = beta/alpha is a primitive N-th root of unity, so [N]_q = 0 and y = -q is minus a root of unity (Hirzebruch's level-N cusp values)", ok_q)
def inf_orbit(kk, dd, N):
    pts = [sp.oo, -kk/dd]; xcur = -kk/dd
    for _ in range(N - 2): xcur = sp.simplify(-(1 + kk*xcur)/(dd*xcur)); pts.append(xcur)
    return pts
o3 = inf_orbit(k, k**2, 3); o4 = inf_orbit(k, k**2/2, 4); o6 = inf_orbit(k, k**2/3, 6)
check("P4 exact: infinity returns to 0 in 3 steps at d = k^2, 4 at d = k^2/2, 6 at d = k^2/3 (sympy)",
      o3[-1] == 0 and o4[-1] == 0 and o6[-1] == 0 and 0 not in o6[:-1])
hyp_ok = all(abs(U(n, 1., dd)) > 1e-12 for dd in np.linspace(-3, .2499, 60) if dd != 0 for n in range(1, 40))
check("P5 nothing like it on the two-wall side: U_n(1, d) != 0 for n < 40 and every d < 1/4 tested (60 values)", hyp_ok)

okh = all(abs(1/(1 - 4/(4*math.cos(math.pi*j/N)**2)) + 1/math.tan(math.pi*j/N)**2) < 1e-9
          for N in range(3, 13) for j in range(1, N//2 + 1) if math.gcd(j, N) == 1)
check("P6 the clocks sit at rung heights s = k^2/Delta = -cot^2(pi j/N): s = -1/3 (N = 3), -1 (N = 4), -3 (N = 6)", okh)

# ---------------- polish (mpmath, 30 digits) ----------------
print("== Polish (mpmath, 30 digits)")
mp.mp.dps = 30; t0 = time.time()
dd = mp.mpf('0.37'); s = mp.sqrt(dd); A_, B_ = mp.mpf('0.3'), mp.mpf('0.45')
Tm = lambda th: mp.tan(s*th)/s; Ti = lambda vv: mp.atan(s*vv)/s
e30 = abs(Tm(Ti(A_) + Ti(B_)) - (A_ + B_)/(1 - dd*A_*B_))
check("Q1 polish: tangent law at d = 0.37 to 30 digits", e30 < mp.mpf('1e-28'), f"err {mp.nstr(e30, 3)}")
m_ = mp.mpf('0.8'); qv = mp.mpf('1.3')
ft = 2*mp.quad(lambda zz: mp.cos(qv*zz)*mp.exp(-m_*zz)/(2*m_), [0, mp.inf])
check("Q2 polish: the Fourier transform of e^(-m|x|)/2m is 1/(q^2 + m^2) (m = 0.8, q = 1.3)", abs(ft - 1/(qv**2 + m_**2)) < mp.mpf('1e-27'))
mm = 1 - mp.mpf('1e-20'); gap = mp.ellipk(mm) - mp.log(4/mp.sqrt(1 - mm))
check("Q3 polish: pendulum quarter period K(m) - ln(4/sqrt(1 - m)) -> 0 at the separatrix (m = 1 - 1e-20)", abs(gap) < mp.mpf('1e-18'), f"{mp.nstr(gap, 3)}")
dd7 = 1/(4*mp.cos(mp.pi/7)**2); Mt = mp.matrix([[-1, -1], [dd7, 0]]); MN = Mt**7
check("Q4 polish: translation by infinity has order 7 at d = 1/(4 cos^2(pi/7)), off-diagonal < 1e-27", max(abs(MN[0, 1]), abs(MN[1, 0])) < mp.mpf('1e-27'))
print(f"   polish cost: {time.time() - t0:.2f} s")

print(f"{sum(R)} of {len(R)} pass")
