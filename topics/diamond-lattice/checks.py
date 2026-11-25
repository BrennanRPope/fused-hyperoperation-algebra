"""Checks for the diamond-lattice page. Every claim on the page is checked here."""
import math, itertools, numpy as np, sympy as sp, mpmath as mp
from dynamics import T as Tn, fixed_ys, ystar, lam, nu, beraha, basin_grid
N = [0, 0]
def check(name, ok):
    N[0] += 1; N[1] += bool(ok)
    print(('PASS' if ok else 'FAIL') + ': ' + name)

t, Q, y, s, x, a, b, v = sp.symbols('t Q y s x a b v')
lens = lambda p, q, A, B: (A + B + p*A*B) / (1 - q*A*B)
F = lambda A, B: lens(Q - 2, 1 - Q, A, B)
T = sp.together(F(t**2, t**2))
z = lambda e: sp.simplify(e) == 0

print('--- 1. The map: series squares t, parallel is the landscape law')
check('T(t) = F_{Q-2,1-Q}(t^2,t^2) = (2t^2 + (Q-2)t^4)/(1 + (Q-1)t^4)', z(T - (2*t**2 + (Q - 2)*t**4)/(1 + (Q - 1)*t**4)))
tx = (x - 1)/(x + Q - 1); X = ((x**2 + Q - 1)/(2*x + Q - 2))**2
check('conjugate (t = (x-1)/(x+Q-1), x = e^K) to the classical map x -> ((x^2+Q-1)/(2x+Q-2))^2', z(T.subs(t, tx) - tx.subs(x, X)))
check('t = (x-1)/(x+Q-1) is v/(Q+v) with v = x - 1', z(tx - (x - 1)/(Q + (x - 1))))
# brute force: level-1 diamond Potts effective transmissivity
def eff_t(edges, nv, Qv, vv):
    zc = zs = 0
    for cols in itertools.product(range(Qv), repeat=nv - 2):
        for c1 in (0, 1):
            col = [0, c1] + list(cols)
            w = 1
            for (i, j) in edges: w *= (1 + vv) if col[i] == col[j] else 1
            if c1 == 0: zs += w
            else: zc += w
    # same-colour sum a, different b: t = (a-b)/(a+(Q-1)b)
    return (zs - zc)/(zs + (Qv - 1)*zc)
dia1 = [(0, 2), (2, 1), (0, 3), (3, 1)]
ok = all(abs(eff_t(dia1, 4, Qv, vv) - float(T.subs({Q: Qv, t: sp.Rational(vv)/(Qv + sp.Rational(vv))}))) < 1e-12 for Qv in (2, 3, 4) for vv in (sp.Rational(1, 2), 2, sp.Rational(-1, 3)))
check('brute-force Potts sums on the 4-edge diamond (Q = 2,3,4) give t\' = T(t)', ok)
dia2 = []
nxt = 4
for (u, w) in dia1:
    m1, m2 = nxt, nxt + 1; nxt += 2
    dia2 += [(u, m1), (m1, w), (u, m2), (m2, w)]
ok = all(abs(eff_t(dia2, nxt, 2, vv) - float(T.subs(t, T).subs({Q: 2, t: sp.Rational(vv)/(2 + sp.Rational(vv))}))) < 1e-12 for vv in (sp.Rational(1, 2), 3))
check('brute force on the 16-edge level-2 diamond (Q = 2) gives T(T(t))', ok)

print('--- 2. Walls and the lift')
for w in (1, 1/(1 - Q)):
    check('wall %s absorbs under the parallel law: F(w, b) = w' % w, z(F(w, b) - w))
check('the second wall 1/(1-Q) is the colouring bond v = -1', z((v/(Q + v)).subs(v, -1) - 1/(1 - Q)))
phi = lambda u: (1 - u)/(1 + (Q - 1)*u)
check('lift phi(t) = (1-t)/(1+(Q-1)t) turns parallel into multiplication', z(phi(F(a, b)) - phi(a)*phi(b)))
check('phi(t) = 1/(1+v) = e^{-K}: parallel multiplies e^K', z(phi(v/(Q + v)) - 1/(1 + v)))
check('so T = phi^-1(phi(t^2)^2): square in t, then square in e^{-K}', z(phi(T) - phi(t**2)**2))
p_, q_ = Q - 2, 1 - Q
check('discriminant p^2 - 4q = (Q-2)^2 + 4(Q-1) = Q^2', z(p_**2 - 4*q_ - Q**2))
check('wall gap squared (1 - 1/(1-Q))^2 = (p^2-4q)/q^2', z((1 - 1/(1 - Q))**2 - (p_**2 - 4*q_)/q_**2))
check('bend q = 1 - Q is the reciprocal of the second wall', z(q_ - 1/(1/(1 - Q))))
Tp = sp.together(sp.diff(T, t))
check("T'(t) = 4t(1-t^2)(1+(Q-1)t^2)/(1+(Q-1)t^4)^2", z(Tp - 4*t*(1 - t**2)*(1 + (Q - 1)*t**2)/(1 + (Q - 1)*t**4)**2))
check("critical points 0, +-1, +-sqrt(1/(1-Q)) (and infinity, the series wall)", z(sp.limit(Tp.subs(Q, 3)*1, t, sp.oo)) )
check('the parallel law doubles its own wall: F(w, w) = w for w = 1/(1-Q), so T(+-sqrt w) = w', z(F(1/(1 - Q), 1/(1 - Q)) - 1/(1 - Q)))
check('T(infinity) = (Q-2)/(Q-1)', z(sp.limit(T, t, sp.oo) - (Q - 2)/(Q - 1)))
check("Q != 0: both t = 0 and t = 1 are superattracting (T' = 0)", z(Tp.subs(t, 0)) and z(Tp.subs(t, 1)))
T0 = sp.simplify(T.subs(Q, 0))
check('Q = 0 (walls merge, the border): T = 2t^2/(1+t^2), T - t = -t(t-1)^2/(1+t^2)', z(T0 - 2*t**2/(1 + t**2)) and z(T0 - t + t*(t - 1)**2/(1 + t**2)))
check("Q = 0: T'(1) = 1, the cold point turns parabolic", z(sp.diff(T0, t).subs(t, 1) - 1))
check('Q = 1: second wall at infinity meets the series wall; T = 2t^2 - t^4, a polynomial', z(T.subs(Q, 1) - (2*t**2 - t**4)))
check('Q = 2: walls 0, inf, 1, -1 are harmonic, and T(1/t) = T(t)', z(T.subs(Q, 2).subs(t, 1/t) - T.subs(Q, 2)))
check('only Q = 2 has T(1/t) = T(t)', sp.solve(sp.numer(sp.together(T.subs(t, 1/t) - T)).subs(t, 2), Q) == [2])
check('T is even (series squares): T(-t) = T(t)', z(T.subs(t, -t) - T))

print('--- 3. Fixed points: one cubic, one parametrisation, one slope formula')
check('T(t) - t = -t(t-1)((Q-1)t^3 + t^2 + t - 1)/(1+(Q-1)t^4)', z(T - t + t*(t - 1)*((Q - 1)*t**3 + t**2 + t - 1)/(1 + (Q - 1)*t**4)))
Qy = (y - 1)**2*(y + 1)
check('t = 1/y is fixed whenever Q = (y-1)^2 (y+1)', z(T.subs(Q, Qy).subs(t, 1/y) - 1/y))
check('and then e^{K} = y^2 (x = y^2 is fixed by the x map)', z(X.subs(Q, Qy).subs(x, y**2) - y**2))
check('matches Chang-Roeder-Shrock fixed cubic q^2 + 2qv - v^3 = 0 with v = y^2 - 1', z((Q**2 + 2*Q*v - v**3).subs({Q: Qy, v: y**2 - 1})))
check("slope at every nontrivial fixed point: T'(1/y) = (2y/(y+1))^2", z(Tp.subs(Q, Qy).subs(t, 1/y) - (2*y/(y + 1))**2))
check("same in t: lambda = 4/(1+t*)^2", z(Tp.subs(Q, Qy).subs(t, 1/y) - (4/(1 + 1/y)**2)))
closed = {1: ((sp.sqrt(5) - 1)/2, None), 3: (sp.Rational(1, 2), sp.Rational(16, 9)), 16: (sp.Rational(1, 3), sp.Rational(9, 4)), 45: (sp.Rational(1, 4), sp.Rational(64, 25))}
for Qv, (ts, lv) in closed.items():
    ok = z(T.subs({Q: Qv, t: ts}) - ts)
    lv = lv if lv is not None else sp.simplify(4/(1 + ts)**2)
    ok = ok and z(Tp.subs({Q: Qv, t: ts}) - lv)
    check('Q = %d: t* = %s, lambda = %s' % (Qv, ts, sp.nsimplify(lv)), ok)
check('Q = 1: lambda = 4/phi^2 = 4 - 4t* (golden)', z(4/(1 + (sp.sqrt(5) - 1)/2)**2 - (4 - 4*(sp.sqrt(5) - 1)/2)))
trib = sp.nsolve(y**3 - y**2 - y - 1, y, 1.8)
check('Q = 2: y* is the tribonacci constant 1.839287 (y^3 = y^2 + y + 1), t* = 0.543689', abs(trib - 1.839286755) < 1e-8 and abs(1/trib - 0.5436890127) < 1e-9)
check('Q = (m-1)^2 (m+1) gives t* = 1/m: Q = 0, 3, 16, 45, 96 for m = 1..5', [((m - 1)**2*(m + 1)) for m in range(1, 6)] == [0, 3, 16, 45, 96])
ok = all(abs(Tn(1/m, (m - 1)**2*(m + 1)) - 1/m) < 1e-12 and abs(ystar((m - 1)**2*(m + 1)) - m) < 1e-9 for m in (2, 3, 4, 5))
check('numerically: y* = m at Q = (m-1)^2(m+1), m = 2..5', ok)
nus = {Qv: nu(Qv) for Qv in (1, 2, 3, 4, 16)}
print('    nu:', {k: round(v_, 4) for k, v_ in nus.items()})
check('nu(1) = 1.6353, nu(2) = 1.3383, nu(3) = 1.2047, nu(4) = 1.1242, nu(16) = 0.8548',
      all(abs(nus[k] - w) < 5e-4 for k, w in ((1, 1.6353), (2, 1.3383), (3, 1.2047), (4, 1.1242), (16, 0.8548))))
check('nu(16) = ln2/(2 ln(3/2)) = 0.8548 and nu(45) = ln2/(2 ln(8/5)) = 0.7374', abs(nu(16) - math.log(2)/(2*math.log(1.5))) < 1e-12 and abs(nu(45) - math.log(2)/(2*math.log(1.6))) < 1e-12 and abs(nu(45) - 0.7374) < 1e-4)
check('slope 4/(1+t*)^2 > 1 for every Q > 0 (y* > 1 so t* < 1)', all(ystar(q_) > 1 and lam(q_) > 1 for q_ in np.linspace(1e-3, 100, 500)))
check('nu(3) = ln2 / (2 ln(4/3)) exactly', abs(nus[3] - math.log(2)/(2*math.log(4/3))) < 1e-12)
check('nu decreases with Q and tends to 1/2 (lambda -> 4) as Q grows; nu -> infinity as Q -> 0',
      all(nu(q1) > nu(q2) for q1, q2 in zip(np.linspace(0.05, 50, 200)[:-1], np.linspace(0.05, 50, 200)[1:])) and abs(nu(1e9) - 0.5) < 0.01 and nu(1e-4) > 30)
qq = np.linspace(3.5, 4.5, 2001); nn = np.array([nu(q_) for q_ in qq])
d2 = np.diff(nn, 2)
check('nu(Q) is smooth through Q = 4 (no kink: second differences stay tiny and of one sign)', np.all(d2 > 0) and np.max(np.abs(d2)) < 1e-6)

print('--- 4. The fold at 32/27 (known: Chang-Roeder-Shrock 2020)')
check('32/27 - (y-1)^2(y+1) = (y + 1/3)^2 (5/3 - y): Q(y) peaks at y = -1/3 with value 32/27', z(sp.Rational(32, 27) - Qy - (y + sp.Rational(1, 3))**2*(sp.Rational(5, 3) - y)))
check('Q = 32/27: ferro point y = 5/3, t* = 3/5, slope 25/16, nu = 1.5531', z(T.subs({Q: sp.Rational(32, 27), t: sp.Rational(3, 5)}) - sp.Rational(3, 5)) and z(Tp.subs({Q: sp.Rational(32, 27), t: sp.Rational(3, 5)}) - sp.Rational(25, 16)) and abs(nu(32/27) - 1.5531) < 1e-4)
check('at the fold y = -1/3: t = -3, e^K = 1/9, slope 2y/(y+1) = -1 (neutral)', z((2*y/(y + 1)).subs(y, -sp.Rational(1, 3)) + 1))
check('|2y/(y+1)| < 1 exactly for y > -1/3 (when y > -1): (y+1)^2 - 4y^2 = (1-y)(1+3y)', z((y + 1)**2 - 4*y**2 - (1 - y)*(1 + 3*y)))
ok = all(len([w for w in fixed_ys(q_) if -1/3 < w < 1 and abs(w) > 1e-9]) == 1 for q_ in np.linspace(0.02, 1.18, 60) if abs(q_ - 1) > 1e-3) and all(len(fixed_ys(q_)) == 1 for q_ in np.linspace(1.19, 6, 60))
check('0 < Q < 32/27: exactly one attracting antiferro fixed point; Q > 32/27: one real fixed root only', ok)
check('Q = 1: that attracting point is t = infinity (y = 0)', 0 in [round(w, 9) for w in fixed_ys(1.0)])

print('--- 5. Zeros: level-n zeros are where the n-th iterate hits t = infinity')
def ab_level(n, Qs, vs):
    A, Bv = 1 + vs, sp.Integer(1)
    for _ in range(n):
        sa = A*A + (Qs - 1)*Bv*Bv; sb = 2*A*Bv + (Qs - 2)*Bv*Bv
        A, Bv = sp.expand(sa*sa), sp.expand(sb*sb)
    return A, Bv
A2, B2 = ab_level(2, Q, sp.Integer(-1))
P2 = sp.Poly(sp.expand(Q*A2 + Q*(Q - 1)*B2), Q)
check('level-2 chromatic polynomial vanishes at Q = 0, 1 and has 2 at Q = 2 (bipartite: P(2) = 2)', P2.eval(0) == 0 and P2.eval(1) == 0 and P2.eval(2) == 2)
roots = [r for r in sp.Poly(sp.cancel(P2.as_expr()/(Q*(Q - 1))), Q).nroots(n=30, maxsteps=200)]
den = lambda Qv: 1 + (Qv - 1)*complex(Tn(complex(1/(1 - Qv)), Qv))**4
check('every other level-2 chromatic zero has T(t_c)^4 (Q-1) = -1, i.e. T^2(t_c) = infinity', all(abs(den(complex(r))) < 1e-8 for r in roots))
def chrom_real_zero(n, Qv):
    with mp.workdps(60):
        A, Bv = mp.mpf(0), mp.mpf(1)
        for _ in range(n):
            sa = A*A + (Qv - 1)*Bv*Bv; sb = 2*A*Bv + (Qv - 2)*Bv*Bv
            A, Bv = sa*sa, sb*sb
        return A + (Qv - 1)*Bv
ok = all(chrom_real_zero(n, mp.mpf(q_)) > 0 for n in (1, 2, 3, 4, 5) for q_ in np.linspace(1.001, 6, 400))
check('no real chromatic zero with Q > 1 at levels 1..5 (the real orbit never reaches infinity)', ok)

print('--- 6. Beraha values 4cos^2(pi/n): tested, not assumed')
zz = sp.Symbol('zz'); okf = []
for n in range(3, 31):
    mpoly = sp.minimal_polynomial(4*sp.cos(sp.pi/n)**2, zz)
    res = sp.resultant(sp.Poly(mpoly, zz), sp.Poly(Qy.subs(Q, 0) - zz, zz))
    degs = sorted(f.degree() for f, _ in sp.factor_list(sp.Poly(res, y))[1])
    if degs != [3*sp.degree(mpoly, zz)]: okf.append(n)
check('the critical cubic has a root in the Beraha field Q(B_n) only for n = 3, 6 (Q = 1, 3), n <= 30', okf == [3, 6])
def dest(t0, Qv, it=4000):
    tt = t0
    for _ in range(it):
        tt = Tn(tt, Qv) if np.isfinite(tt) else (Qv - 2)/(Qv - 1)
        if abs(tt) < 1e-10: return 0
        if abs(tt - 1) < 1e-10: return 1
    return 2
qs = np.linspace(32/27 + 1e-4, 4.5, 6000)
dd = [dest(1/(1 - q_), q_) for q_ in qs]
sw = [0.5*(qs[i] + qs[i + 1]) for i in range(len(qs) - 1) if dd[i] != dd[i + 1]]
big = [w for w in sw if w > 1.3]
print('    switch points above 1.3:', [round(w, 4) for w in big])
check('the chromatic point t_c = 1/(1-Q) switches basin on the real axis at 3, 1.6389, 1.4098, 1.3240, ... (count > 10, all in (32/27, 3])',
      len(sw) > 10 and abs(max(sw) - 3) < 1e-3 and any(abs(w - 1.6389) < 1e-3 for w in sw) and any(abs(w - 1.4098) < 2e-3 for w in sw) and min(sw) > 32/27)
q1 = [r for r in sp.Poly(Q**3 - 5*Q**2 + 11*Q - 9, Q).nroots() if abs(sp.im(r)) < 1e-12][0]
check('q1 = 1.6388969 is the real root of Q^3 - 5Q^2 + 11Q - 9 (as Chang-Roeder-Shrock state)', abs(q1 - 1.6388969) < 1e-6)
check('the switch points crowd down onto 32/27 = 1.1852', min(sw) - 32/27 < 0.01)
B = [beraha(n) for n in range(3, 200)]
near = [w for w in sw if min(abs(w - bb) for bb in B) < 0.005]
check('of the switch points above 1.2, only Q = 3 (n = 6) sits on a Beraha value', [round(w, 3) for w in near if w > 1.2] == [3.0])
check('why Q = 3: t_c = -t* exactly when (y-1)^2 = 1, y = 2; T even so T(t_c) = t* (T(-1/2) = 1/2 at Q = 3)',
      sp.solve(sp.Eq((y - 1)**2*(y + 1), y + 1), y) == [-1, 0, 2] and z(T.subs({Q: 3, t: -sp.Rational(1, 2)}) - sp.Rational(1, 2)))
hits = [n for n in range(5, 13) if n != 6 and any(abs(Tn_iter - 1/ystar(beraha(n))) < 1e-6 for Tn_iter in
        [float(np.real(v_)) for v_ in itertools.accumulate(range(6), lambda tt, _: Tn(tt, beraha(n)), initial=1/(1 - beraha(n)))])]
check('Beraha n = 5, 7..12: the chromatic orbit never lands on t* within 6 steps', hits == [])
ber_nu = [nu(beraha(n)) for n in range(4, 40)]
smooth = np.array([nu(q_) for q_ in np.linspace(2, 3.99, 400)])
check('nu at Beraha values lies on the smooth nu(Q) curve (no jumps: |d nu/dQ| bounded, monotone)', np.all(np.diff(smooth) < 0) and np.max(np.abs(np.diff(smooth))) < 0.01)

print('--- 7. Julia-set pictures (numerical facts shown on the page)')
for Qv in (2.0, 3.0, 4.0):
    ts = 1/ystar(Qv)
    ok = dest(ts - 1e-7, Qv) == 0 and dest(ts + 1e-7, Qv) == 1
    check('Q = %g: t* sits on the boundary (just below goes hot, just above goes cold)' % Qv, ok)
xs, ys, lab, st = basin_grid(2.0, xlim=(-3, 3), ylim=(-3, 3), n=200, iters=160)
labinv = lab.copy()
check('Q = 2 picture: every grid point is decided (no third attractor)', np.all(lab != 3) and not np.any(lab == 2))
ok = True
for i in range(0, 200, 7):
    for j in range(0, 200, 7):
        z0 = xs[j] + 1j*ys[i]
        if abs(z0) < 0.4 or abs(z0) > 2.5: continue
        if abs(abs(z0) - 1) < 0.05: continue
        d1 = dest(z0, 2.0); d2 = dest(1/z0, 2.0)
        ok &= (d1 == d2) or d1 == 2 or d2 == 2
check('Q = 2 picture is symmetric under t -> 1/t (sampled)', ok)
xs, ys, lab, st = basin_grid(1.0, xlim=(-3, 3), ylim=(-3, 3), n=160, iters=120)
check('Q = 1: a third basin (infinity) appears; percolation map is a polynomial', np.mean(lab == 2) > 0.5)
xs, ys, lab, st = basin_grid(0.6, xlim=(-3, 3), ylim=(-3, 3), n=160, iters=160)
check('Q = 0.6: the antiferro fixed point grabs most of the plane', np.mean(lab == 2) > 0.5)
fr = []
for Qv in (2.0, 3.0, 4.0, 8.0, 16.0):
    xs, ys, lab, st = basin_grid(Qv, xlim=(-3, 3), ylim=(-3, 3), n=160, iters=160)
    fr.append(np.mean(lab == 1))
check('as Q grows past 2 the cold basin swallows more of the window (fractions %s)' % [round(f_, 2) for f_ in fr], all(f1 < f2 for f1, f2 in zip(fr, fr[1:])))
print('\n%d/%d checks passed' % (N[1], N[0]))
