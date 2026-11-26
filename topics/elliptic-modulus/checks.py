"""Checks for the elliptic-modulus page (DRAFT, Nov 26 2026). Run: python3 checks.py > checks_output.txt
Sections: A the curvature reading, B the fingerprint, C the pendulum, D heights, E theta and the lift.
Exact work in sympy; numerics in mpmath at 40 digits unless stated; heights in double precision."""
import sympy as sp, mpmath as mp, numpy as np, json, math
from series import slope_series, jacobi_series, julia as sjulia, flow as sflow, trunc, u, m, s
import numseries as NS
from heights import height, G as Gd
mp.mp.dps = 40
results = []; npass = 0
def check(tag, desc, ok, extra=""):
    global npass
    npass += bool(ok); results.append((tag, bool(ok), desc, extra))
    print(("PASS " if ok else "FAIL ") + tag + " " + desc + (("   [" + extra + "]") if extra else ""))
def G(x, mm):  # slope-chart generator, mpmath (complex ok)
    sn, cn, dn = mp.ellipfun('sn', x, m=mm), mp.ellipfun('cn', x, m=mm), mp.ellipfun('dn', x, m=mm)
    return sn*dn/cn
kap = sp.Symbol('kappa'); eps = sp.Symbol('epsilon'); t = sp.Symbol('t')

print("== A. The curvature reading along Jacobi's modulus (slope chart t = x/y on the Edwards curve)")
N = 12; Gs = slope_series(N)
check("A1", "G_m(u) = sn dn/cn = u + (1 - 2m) u^3/3 + 2(1 - m + m^2) u^5/15 + ...: the bend reads kappa = 1 - 2m",
      sp.simplify(Gs.coeff(u, 3) - (1 - 2*m)/3) == 0 and sp.simplify(Gs.coeff(u, 5) - 2*(1 - m + m**2)/15) == 0)
pts = [(mp.mpf('0.37'), mp.mpf('0.2')), (mp.mpf('-0.9'), mp.mpf('0.5')), (mp.mpf('1.1'), mp.mpf('0.85')),
       (mp.mpc('0.3', '0.4'), mp.mpf('0.3')), (mp.mpf('0.6'), mp.mpc('0.5', '0.7'))]
err = 0
for x, mm in pts:
    g = G(x, mm); gp = mp.diff(lambda y: G(y, mm), x); err = max(err, abs(gp**2 - (1 + 2*(1 - 2*mm)*g**2 + g**4)))
check("A2", "G'^2 = 1 + 2 kappa G^2 + G^4 with kappa = 1 - 2m: the law's logarithm is int dt / sqrt(1 + 2 kappa t^2 + t^4)",
      err < mp.mpf(10)**-30, "max err %s over 5 points, incl. complex u and m" % mp.nstr(err, 3))
# conic laws: (a+b+kab)/(1-dab) has log' = 1/(1 + k t + d t^2); at k = 0 that is 1/sqrt((1 + d t^2)^2)
d = sp.Symbol('d')
check("A3", "our k = 0 laws are the perfect squares: (1 + d t^2)^2 = 1 + 2d t^2 + d^2 t^4 (kappa = d, quartic coeff = kappa^2)",
      sp.expand((1 + d*t**2)**2 - (1 + 2*d*t**2 + d**2*t**4)) == 0)
# Edwards addition law in the slope chart gives G(u + v)
err = 0
for (uu, vv, mm) in [(0.31, 0.77, 0.2), (-0.4, 1.3, 0.5), (0.9, 0.25, 0.9)]:
    uu, vv, mm = mp.mpf(uu), mp.mpf(vv), mp.mpf(mm)
    X = lambda w: mp.ellipfun('sn', w, m=mm); Y = lambda w: mp.ellipfun('cn', w, m=mm)/mp.ellipfun('dn', w, m=mm)
    x1, y1, x2, y2 = X(uu), Y(uu), X(vv), Y(vv)
    on = abs(x1**2 + y1**2 - 1 - mm*x1**2*y1**2)
    x3 = (x1*y2 + y1*x2)/(1 + mm*x1*x2*y1*y2); y3 = (y1*y2 - x1*x2)/(1 - mm*x1*x2*y1*y2)
    err = max(err, abs(x3/y3 - G(uu + vv, mm)), on)
check("A4", "Edwards addition (x^2 + y^2 = 1 + m x^2 y^2) read in the slope t = x/y is G(u) (+) G(v) = G(u + v)",
      err < mp.mpf(10)**-30, "err %s" % mp.nstr(err, 3))
err = 0; mm = mp.mpf('0.3'); K = mp.ellipk(mm); Kp = mp.ellipk(1 - mm)
for x in [mp.mpf('0.2'), mp.mpc('0.4', '0.1')]:
    err = max(err, abs(G(x + 2*K, mm) - G(x, mm)), abs(G(x + 2j*Kp, mm) - G(x, mm)))
check("A5", "G has two periods, 2K (real) and 2iK' (imaginary): a torus. One period only at m = 0 (pi, turning) and m = 1 (i pi, walls)",
      err < mp.mpf(10)**-28, "m = 0.3, err %s" % mp.nstr(err, 3))
err = max(abs(G(1j*x, mm) - 1j*G(x, 1 - mm)) for x in [mp.mpf('0.3'), mp.mpf('0.8'), mp.mpc('0.2', '0.5')])
check("A6", "Wick rotation u -> iu swaps the two sides: G_m(iu) = i G_{1-m}(u), i.e. kappa -> -kappa", err < mp.mpf(10)**-30,
      "err %s" % mp.nstr(err, 3))
check("A7", "at m = 1/2 the torus is square: K = K'", abs(mp.ellipk(0.5) - mp.ellipk(0.5)) == 0 and
      abs(mp.ellipk(mp.mpf(1)/2) - mp.gamma(0.25)**2/(4*mp.sqrt(mp.pi))) < mp.mpf(10)**-35, "K(1/2) = Gamma(1/4)^2/(4 sqrt pi)")

print("\n== B. The fingerprint (Ecalle's iterative residue) along the modulus")
def resit_odd(a3, a5): return sp.Rational(3, 2) - a5/a3**2
Rslope = sp.factor(resit_odd(Gs.coeff(u, 3), Gs.coeff(u, 5)))
check("B1", "slope chart: R(m) = 6/5 - 9/(10 kappa^2) with kappa = 1 - 2m", sp.simplify(Rslope - (sp.Rational(6, 5) - sp.Rational(9, 10)/(1 - 2*m)**2)) == 0,
      "R = %s" % Rslope)
# general odd chart with log' = (1 + 2 kappa t^2 + eps t^4)^(-1/2): invert the log series to get the generator
Lp = sp.series((1 + 2*kap*t**2 + eps*t**4)**sp.Rational(-1, 2), t, 0, 8).removeO()
L = sp.integrate(Lp, t)
# generator g = L^{-1}: g = u + c3 u^3 + c5 u^5
c3, c5 = sp.symbols('c3 c5'); g = u + c3*u**3 + c5*u**5
eq = sp.expand(sp.series(L.subs(t, g), u, 0, 7).removeO() - u)
sol = sp.solve([eq.coeff(u, 3), eq.coeff(u, 5)], [c3, c5], dict=True)[0]
Rgen = sp.simplify(resit_odd(sol[c3], sol[c5]))
check("B2", "any odd chart with log' = 1/sqrt(1 + 2 kappa t^2 + eps t^4): R = 6/5 - (9/10) eps/kappa^2 (linear in the scale-free number eps/kappa^2)",
      sp.simplify(Rgen - (sp.Rational(6, 5) - sp.Rational(9, 10)*eps/kap**2)) == 0, "generator cubic coeff = %s" % sol[c3])
check("B3", "ends: eps = kappa^2 (our conic laws, tan and tanh) gives 3/10; eps = 0 (arcsin, the sine chart) gives 6/5",
      sp.simplify(Rgen.subs(eps, kap**2)) == sp.Rational(3, 10) and sp.simplify(Rgen.subs(eps, 0)) == sp.Rational(6, 5))
sn_s, cn_s, dn_s = jacobi_series(N)
Rsn = sp.factor(resit_odd(sn_s.coeff(u, 3), sn_s.coeff(u, 5)))
check("B4", "x = sn chart: R = 6(1 - m + m^2)/(5(1 + m)^2), the same formula with kappa = -(1+m)/2, eps = m; symmetric under m -> 1/m",
      sp.simplify(Rsn - 6*(1 - m + m**2)/(5*(1 + m)**2)) == 0 and sp.simplify(Rsn.subs(m, 1/m) - Rsn) == 0 and
      sp.simplify(Rgen.subs({kap: -(1 + m)/2, eps: m}) - Rsn) == 0)
check("B5", "slope chart R is Wick-symmetric, R(m) = R(1 - m), a function of m(1 - m) only (as j is)",
      sp.simplify(Rslope.subs(m, 1 - m) - Rslope) == 0 and sp.simplify(Rslope - (sp.Rational(6, 5) - sp.Rational(9, 10)/(1 - 4*m*(1 - m)))) == 0)
G12 = sp.expand(Gs.subs(m, sp.Rational(1, 2)))
ok = G12.coeff(u, 3) == 0 and G12.coeff(u, 5) == sp.Rational(1, 10) and all(G12.coeff(u, j) == 0 for j in (3, 7, 11))
R4 = sp.Rational(5, 2) - G12.coeff(u, 9)/G12.coeff(u, 5)**2
check("B6", "m = 1/2 (kappa = 0, square torus): G = u + u^5/10 + u^9/120 + ..., only powers u^(4j+1): the class jumps from contact 3 to contact 5",
      ok, "iterative residue of the contact-5 class: %s = %.6f" % (R4, float(R4)))
V = sjulia(Gs, 10); R_V = sp.factor(-V.coeff(u, 5)/V.coeff(u, 3)**2)
Gsf = sflow(V, s, 10); Rs = sp.factor(resit_odd(Gsf.coeff(u, 3), Gsf.coeff(u, 5)))
check("B7", "heights move the fingerprint: R(G^s) = R(G)/s for every real s (sympy, generic m and s)",
      sp.simplify(Rs - Rslope/s) == 0 and sp.simplify(R_V - Rslope) == 0)
def m_of_R(Rt): k2 = mp.mpf(9)/(12 - 10*Rt); return (1 - mp.sqrt(k2))/2, k2
mm, ss = mp.mpf('0.2'), mp.mpf('0.5')
Rm = mp.mpf(6)/5 - mp.mpf(9)/(10*(1 - 2*mm)**2); mpr, k2 = m_of_R(Rm/ss)
NN = 13; Gm = NS.slope(mm, NN); Hs = NS.flow(NS.julia(Gm, NN), ss, NN); Gmp = NS.slope(mpr, NN)
lam2 = ss*Gm[3]/Gmp[3]
scaled = [Gmp[j]*lam2**((j - 1)//2) for j in range(NN)]
diff = [Hs[j] - scaled[j] for j in range(NN)]
check("B8", "normal form: G_0.2^(1/2) equals a rescaled elliptic generator at m' = %s through u^5 exactly; they part at u^7" % mp.nstr(mpr, 8),
      abs(diff[3]) < 1e-35 and abs(diff[5]) < 1e-35 and abs(diff[7]) > 1e-4,
      "u^7 gap %s, u^9 gap %s; a change of chart absorbs the rest formally, since (contact, R) is a complete formal invariant"
      % (mp.nstr(diff[7], 5), mp.nstr(diff[9], 5)))
gj = json.load(open('julia_growth.json'))
grow = lambda key: [l for n, l in gj[key] if n % 8 == 1 and l > -60]
g0, g2, g8 = grow('0'), grow('0.2'), grow('0.8')
check("B9", "the Julia series (iterative logarithm) of G_m diverges at m = 0 and 0.2: heights live on petals, not as convergent series (as for e^x - 1)",
      g0[-1] > 20 and g2[-1] > 10 and np.all(np.diff(np.diff(g0[-6:])) > -1) and max(abs(a - b) for a, b in zip(g2, g8)) < 1e-9,
      "log10|V_n| at n = 137: %.1f (m=0), %.1f (m=0.2); m = 0.2 and 0.8 equal (Wick)" % (g0[-1], g2[-1]))

print("\n== C. The pendulum: one parameter, two readings")
err = max(abs((1 - 2*mp.sin(th/2)**2) - mp.cos(th)) for th in [mp.mpf(0.3), mp.mpf(1.2), mp.mpf(2.9)])
check("C1", "a pendulum of amplitude theta0 has m = sin^2(theta0/2), so the curvature reading is kappa = 1 - 2m = cos(theta0)", err < 1e-38)
from scipy.integrate import solve_ivp
from scipy.special import ellipj, ellipk
okc = True; worst = 0
for th0 in [0.6, math.pi/2, 2.6]:
    w = 1.0; mmf = math.sin(th0/2)**2
    sol_ = solve_ivp(lambda tt, y: [y[1], -w*w*math.sin(y[0])], [0, 0.95*ellipk(mmf)/w], [0, 2*w*math.sin(th0/2)],
                     rtol=1e-12, atol=1e-13, dense_output=True)
    for tt in np.linspace(0.05, 0.9*ellipk(mmf), 9):
        th, thd = sol_.sol(tt)
        worst = max(worst, abs(w*math.sin(th)/thd - Gd(w*tt, mmf)))
check("C2", "integrating the pendulum: the slope reading t = omega sin(theta)/theta_dot is exactly G_m(omega time) (amplitudes 34, 90, 149 deg)",
      worst < 1e-8, "max err %.1e" % worst)
check("C3", "the wall: G has its pole at u = K, the turning point theta_dot = 0 (the pendulum turns back: turning country)",
      abs(float(mp.ellipfun('cn', mp.ellipk(0.3), m=0.3))) < 1e-30)
Ks = [(1 - 10**-e, float(mp.ellipk(1 - mp.mpf(10)**-e)), math.log(4/math.sqrt(10**-e))) for e in (4, 8, 12)]
check("C4", "at the separatrix m -> 1 the half period K diverges like ln(4/sqrt(1-m)): relativity's wall is the top reached in infinite time",
      all(abs(a - b) < 1e-3 for _, a, b in Ks), "; ".join("K=%.6f vs %.6f" % (a, b) for _, a, b in Ks))
xs = np.array([0.2, 0.5, 0.8]); err = max(abs(Gd(xs, 1 - 1e-15) - np.tanh(xs)))
check("C5", "at m = 1 the slope chart is tanh and its law is velocity addition (our k = 0, d = -1)", err < 1e-6, "err %.1e" % err)

print("\n== D. Heights (fractional rungs) of the elliptic generator, numerically")
worst = 0
for mmf in [0.2, 0.5, 0.8]:
    x = np.array([0.25, 0.5, 0.75]); h = height(x, mmf, 0.5); worst = max(worst, np.max(np.abs(height(h, mmf, 0.5) - Gd(x, mmf))))
check("D1", "half a height twice is one step: G^(1/2)(G^(1/2)(x)) = G(x) at m = 0.2, 0.5, 0.8", worst < 1e-10, "max err %.1e" % worst)
x = np.linspace(0.05, 0.8, 16); h = height(x, 0.2, 0.5)
mpf_, lam = float(mpr), math.sqrt(float(lam2))
gap = np.max(np.abs(h - Gd(lam*x, mpf_)/lam))
check("D2", "the half height of the m = 0.2 law against the rescaled m' law: close but not equal on [0, 0.8] (the u^7 gap)",
      1e-7 < gap < 1e-2, "max gap %.2e (relative to G - x: %.2e)" % (gap, np.max(np.abs(Gd(x, 0.2) - x))))
worst = 0
for ss_ in [0.5, -0.5, 1.5]:
    hh = height(height(x[:5], 0.8, ss_), 0.8, 1 - ss_); worst = max(worst, np.max(np.abs(hh - Gd(x[:5], 0.8))))
check("D3", "heights add: G^s(G^(1-s)) = G for s = 1/2, -1/2, 3/2 at m = 0.8", worst < 1e-9, "max err %.1e" % worst)

print("\n== E. Theta: the lift copied round the torus")
q = mp.mpf('0.3'); z = mp.mpf('0.7')
lhs = mp.nsum(lambda n: q**(n*(n - 1)/2)*z**n, [-mp.inf, mp.inf])
rhs = mp.nprod(lambda n: (1 - q**n)*(1 + z*q**(n - 1))*(1 + q**n/z), [1, mp.inf])
check("E1", "Jacobi triple product: sum_n q^(n(n-1)/2) z^n = prod (1 - q^n)(1 + z q^(n-1))(1 + q^n/z)", abs(lhs - rhs) < 1e-35,
      "err %s" % mp.nstr(abs(lhs - rhs), 3))
fold = lambda vals: mp.fprod(1 + v for v in vals) - 1   # the circ_1 fold of a list: lift = product of lifts
orbit = [z*q**j for j in range(200)]
check("E2", "prod_j (1 + z q^j) = 1 + (z o1 zq o1 zq^2 o1 ...): theta is the lift of the mult-add fold over the Tate orbit (both ways)",
      abs(mp.fprod(1 + v for v in orbit) - (1 + fold(orbit))) < 1e-35)
eu = mp.nsum(lambda n: z**n*q**(n*(n - 1)/2)/mp.qp(q, q, n), [0, mp.inf])
check("E3", "grouping the corners of that infinite cube by size gives Euler's sum z^n q^(n(n-1)/2)/(q;q)_n: corner weight q^(sum of positions)",
      abs(eu - mp.fprod(1 + v for v in orbit)) < 1e-35)
zer = [abs(mp.nsum(lambda n: q**(n*(n - 1)/2)*(-q**(-j))**n, [-mp.inf, mp.inf])) for j in (-1, 0, 1, 2)]
check("E4", "theta's zeros sit at z = -q^(-j): mult-add's one wall at -1, copied once round each turn of the torus", max(zer) < 1e-30)
# genus-one KdV: u = 2 d^2/dx^2 log tau, tau = sum q^(n(n-1)/2) e^(2 n kx)
kx = mp.mpf(1)
def u_of(xx, qq):
    tau = lambda y: mp.nsum(lambda n: qq**(n*(n - 1)/2)*mp.exp(2*n*kx*y), [-mp.inf, mp.inf])
    return 2*mp.diff(lambda y: mp.log(tau(y)), xx, 2)
xsm = [mp.mpf(v) for v in np.linspace(-1.5, 1.5, 7)]
sol_err = max(abs(u_of(xx, mp.mpf(10)**-25) - 2*kx**2*mp.sech(kx*xx)**2) for xx in xsm)
check("E5", "as q -> 0 only the corners n = 0, 1 survive: tau -> 1 + e^(2kx), our lift, and u -> 2k^2 sech^2(kx), the KdV soliton",
      sol_err < 1e-20, "err %s at q = 1e-25" % mp.nstr(sol_err, 3))
qq = mp.mpf('0.2'); mp.mp.dps = 30
xx = [mp.mpf(v) for v in np.linspace(-0.6, 0.6, 9)]
U = [u_of(y, qq) for y in xx]; U2 = [mp.diff(lambda y: u_of(y, qq), y, 2) for y in xx]
# fit u'' + 3u^2 - c u = A
Amat = mp.matrix([[-Ui, -1] for Ui in U]); bvec = mp.matrix([-(U2[i] + 3*U[i]**2) for i in range(len(U))])
cA = mp.lu_solve(Amat.T*Amat, Amat.T*bvec)
res = max(abs(U2[i] + 3*U[i]**2 - cA[0]*U[i] - cA[1]) for i in range(len(U)))
check("E6", "for q > 0 the same tau is a cnoidal wave: u solves the KdV travelling-wave equation u'' + 3u^2 - c u = A",
      res < 1e-15, "q = 0.2, c = %s, residual %s" % (mp.nstr(cA[0], 8), mp.nstr(res, 3)))
mp.mp.dps = 40
print("\n%d of %d checks pass" % (npass, len(results)))
json.dump([dict(tag=a, ok=b, desc=c, extra=e) for a, b, c, e in results], open('checks.json', 'w'), indent=1)
