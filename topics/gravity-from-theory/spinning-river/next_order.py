"""DRAFT. Round 13: do the three rules still hold at the next level of detail (2PN)?

The test quantity is the redshift of body 1 on a circular orbit, z1 = 1/u^t, written in
x = (m omega)^(2/3). It does not depend on coordinates, so the rules' answer can be held
against general relativity's answer number for number.

Truth: the first law of binary mechanics (Le Tiec, Blanchet and Whiting 2012),
z1 = dM/dm1 - omega dJ/dm1 at fixed omega, applied to the known 2PN energy and angular
momentum of a circular binary.

Rules: each body's own field (linear part with its light-speed lag kept to the order needed,
plus the nonlinear part of a moving Schwarzschild body), clock factors multiply (rule B),
the frame-drag and space parts add, masses dressed by the local field (rule C, full form
mu = m / (sqrt(-g) sqrt(-g_uv v^u v^v))). Each body moves on a geodesic of the field with its
own infinite self-part removed by Hadamard's finite part (the method of Blanchet, Faye and
Ponsot 1998).

Units G = c = 1, total mass m = 1, omega = eps, so x = eps^2.

    python3 next_order.py          # writes next_order.json (about 5 minutes)
    from next_order import redshift; redshift(Rational(1,3), order=6)
"""
import functools, json, sys, time
import sympy as sp

e = sp.symbols('eps', positive=True)


def tr(x, o):
    """Drop every term above eps^o (after expanding, every term is a monomial in eps)."""
    x = sp.expand(x)
    return sp.Add(*[u for u in sp.Add.make_args(x) if u.as_powers_dict().get(e, 0) <= o])


# ---------------------------------------------------------------- truth (first law)
def truth(m1):
    """z1 coefficients of x^0..x^3 from the first law with the known 2PN E and J (m = 1)."""
    m1 = sp.nsimplify(m1); m2 = 1 - m1
    return [sp.Integer(1),
            -m2*(2*m1 + 3*m2)/2,
            -m2*(36*m1**3 + 100*m1**2*m2 + 96*m1*m2**2 + 27*m2**3)/24,
            -m2*(54*m1**5 + 224*m1**4*m2 + 356*m1**3*m2**2 + 284*m1**2*m2**3
                 + 124*m1*m2**4 + 27*m2**5)/16]


def truth_from_first_law():
    """Derive truth() symbolically from E(omega) and J(omega); returns the four coefficients."""
    m1, m2, w, X = sp.symbols('m1 m2 omega X', positive=True)
    m = m1 + m2; nu = m1*m2/m**2; mu = m1*m2/m
    x = (m*w)**sp.Rational(2, 3)
    E = -(mu*x/2)*(1 + (-sp.Rational(3, 4) - nu/12)*x
                   + (-sp.Rational(27, 8) + sp.Rational(19, 8)*nu - nu**2/24)*x**2)
    J = (mu*m/sp.sqrt(x))*(1 + (sp.Rational(3, 2) + nu/6)*x
                           + (sp.Rational(27, 8) - sp.Rational(19, 8)*nu + nu**2/24)*x**2)
    z1 = 1 + sp.diff(E, m1) - w*sp.diff(J, m1)
    z1X = sp.simplify(z1.subs(w, X**sp.Rational(3, 2)/m))
    s = sp.series(z1X, X, 0, 4).removeO()
    return [sp.factor(sp.simplify(s.coeff(X, k).subs(m2, 1 - m1))) for k in range(4)], m1


# ---------------------------------------------------------------- moving Schwarzschild, nonlinear part
@functools.lru_cache(None)
def boosted_nonlinear():
    """Harmonic Schwarzschild boosted to velocity v, expanded in eps; returns its m^2 and
    higher parts per metric component as functions of symbols m, r, rx.., vx.. ."""
    m, r = sp.symbols('m r', positive=True)
    rx, ry, rz, vx, vy, vz = sp.symbols('rx ry rz vx vy vz', real=True)
    rv = sp.Matrix([rx, ry, rz]); v = sp.Matrix([vx, vy, vz])*e
    v2 = sp.expand((v.T*v)[0]); rdv = sp.expand((rv.T*v)[0])
    g2 = 1 + v2 + v2**2 + v2**3
    g = 1 + v2/2 + sp.Rational(3, 8)*v2**2 + sp.Rational(5, 16)*v2**3
    k = sp.Rational(1, 2) + sp.Rational(3, 8)*v2 + sp.Rational(5, 16)*v2**2
    Xv = rv + k*rdv*v
    q = tr(g2*rdv**2/r**2, 6)

    def Rinv(n, order):
        s = sum(sp.binomial(-sp.Rational(n, 2), j)*q**j for j in range(4))
        return tr(s, order)/r**n
    M = m*e**2; O = 6
    G00 = -(1 - 2*M*Rinv(1, O) + 2*M**2*Rinv(2, O) - 2*M**3*Rinv(3, O))

    def Gij(i, j):
        d = 1 if i == j else 0
        return (1 + 2*M*Rinv(1, O) + M**2*Rinv(2, O))*d + M**2*Rinv(4, O)*Xv[i]*Xv[j]
    L = sp.zeros(4, 4); L[0, 0] = g
    for i in range(3):
        L[0, i+1] = -g*v[i]; L[i+1, 0] = -g*v[i]
        for j in range(3):
            L[i+1, j+1] = (1 if i == j else 0) + k*v[i]*v[j]
    Gr = sp.zeros(4, 4); Gr[0, 0] = tr(G00, O)
    for i in range(3):
        for j in range(3):
            Gr[i+1, j+1] = tr(Gij(i, j), O)
    need = {(0, 0): 6, (0, 1): 5, (0, 2): 5, (0, 3): 5}
    for i in range(3):
        for j in range(i, 3):
            need[(i+1, j+1)] = 4
    out = {}
    for (a, b), o in need.items():
        expr = 0
        for al in range(4):
            for be in range(4):
                if L[al, a] != 0 and L[be, b] != 0 and Gr[al, be] != 0:
                    expr += tr(L[al, a]*L[be, b], o)*Gr[al, be]
        P = sp.Poly(tr(expr, o), m)
        out[(a, b)] = sum(c*m**d[0] for d, c in zip(P.monoms(), P.coeffs()) if d[0] >= 2)
    return out


# ---------------------------------------------------------------- the rules' field and its regularized value at each body
def regularized_field(m1, order=6, kappa=1, multiply=True):
    """Build the rules' two-body metric for a circular orbit and take its Hadamard finite part
    (and that of its x-gradient, and of sqrt(-g)) at each body."""
    X, Y, Z, t, s = sp.symbols('X Y Z t s', real=True)
    nx, ny, nz = sp.symbols('nx ny nz', real=True)
    R1, R2 = sp.symbols('R1 R2', positive=True)
    r1, r2 = sp.symbols('r1 r2', positive=True)
    p1, p2, q1, q2 = sp.symbols('p1 p2 q1 q2')
    m1 = sp.nsimplify(m1); m2 = 1 - m1
    o00, o0i, oij = order, order - 1, order - 2
    sgn = {1: 1, 2: -1}; Rs = {1: R1, 2: R2}; rs = {1: r1, 2: r2}
    mus = {1: e**2*(m1 + e**2*p1 + e**4*q1), 2: e**2*(m2 + e**2*p2 + e**4*q2)}
    xs = sp.Matrix([X, Y, Z])

    def traj(a):
        return sgn[a]*Rs[a]*sp.Matrix([sp.cos(e*t), sp.sin(e*t), 0])

    def at0(expr, a):
        ex = sp.expand(expr.subs(t, 0))
        B = sp.expand((X - sgn[a]*Rs[a])**2 + Y**2 + Z**2)
        ex = ex.replace(lambda u: u.is_Pow and sp.expand(u.base) == B,
                        lambda u: rs[a]**(2*u.exp))
        return sp.expand(ex)

    def pot(a):
        y = traj(a); d = xs - y; ra = sp.sqrt((d.T*d)[0])
        vt = y.diff(t)
        v2 = sp.expand(sp.simplify((vt.T*vt)[0]))
        mu = mus[a]; sig = mu*(1 + v2)
        V = sig*(1/rs[a] + at0(sp.diff(ra, t, 2), a)/2 + at0(sp.diff(ra**3, t, 4), a)/24)
        Vi = [mu*(at0(vt[i]/ra, a) + at0(sp.diff(vt[i]*ra, t, 2), a)/2) for i in range(3)]
        v0 = [vt[i].subs(t, 0) for i in range(3)]
        W = [[mu*(v0[i]*v0[j] - (1 if i == j else 0)*v2)/rs[a] for j in range(3)] for i in range(3)]
        return tr(V, o00), [tr(u, o0i) for u in Vi], [[tr(u, oij) for u in row] for row in W], v0, v2

    NL = boosted_nonlinear()
    mS, rS = sp.symbols('m r', positive=True)
    rx, ry, rz, vx, vy, vz = sp.symbols('rx ry rz vx vy vz', real=True)

    def nl(a, v0, v2):
        meff = mus[a]*(1 - v2/2 - v2**2/8)          # mu sqrt(1 - v^2)
        sub = {mS: meff/e**2, rS: rs[a], rx: X - sgn[a]*Rs[a], ry: Y, rz: Z,
               vx: v0[0]/e, vy: v0[1]/e, vz: v0[2]/e}
        out = {}
        for kk, vv in NL.items():
            o = o00 if kk == (0, 0) else (o0i if kk[0] == 0 else oij)
            out[kk] = tr(vv.subs(sub), o)
        return out

    own = {}
    for a in (1, 2):
        V, Vi, W, v0, v2 = pot(a)
        own[a] = dict(V=V, Vi=Vi, W=W, N=nl(a, v0, v2))
    g = sp.zeros(4, 4)
    c1 = 1 - 2*own[1]['V'] - own[1]['N'][(0, 0)]
    c2 = 1 - 2*own[2]['V'] - own[2]['N'][(0, 0)]
    # rule B: clock factors multiply (multiply=False adds them instead, to test whether B matters)
    g[0, 0] = -tr(c1*c2, o00) if multiply else -tr(c1 + c2 - 1, o00)
    for i in range(3):
        g[0, i+1] = g[i+1, 0] = tr(sum(-4*own[a]['Vi'][i] + own[a]['N'][(0, i+1)] for a in (1, 2)), o0i)
        for j in range(i, 3):
            d = 1 if i == j else 0
            val = d + sum(2*own[a]['V']*d + 4*own[a]['W'][i][j] + own[a]['N'][(i+1, j+1)]
                          for a in (1, 2)) + kappa*4*own[1]['V']*own[2]['V']*d
            g[i+1, j+1] = g[j+1, i+1] = tr(val, oij)

    Dsym = sp.symbols('Dsym', positive=True)

    @functools.lru_cache(None)
    def taylor(qq, K, sig):
        u = sp.symbols('u')
        ser = sp.series((1 + u)**sp.Rational(qq, 2), u, 0, K + 1).removeO()
        out = sp.expand(ser.subs(u, 2*sig*s*nx/Dsym + s**2/Dsym**2))
        out = sp.Add(*[w for w in sp.Add.make_args(out) if sp.degree(w, s) <= K])
        return sp.expand(out*Dsym**qq)

    @functools.lru_cache(None)
    def avg(a, b, c):
        if a % 2 or b % 2 or c % 2:
            return 0
        df = lambda k: sp.factorial2(k - 1) if k > 0 else 1
        return sp.Rational(df(a)*df(b)*df(c), sp.factorial2(a + b + c + 1))

    def angavg(poly):
        P = sp.Poly(sp.expand(poly), nx, ny, nz)
        return sp.expand(sum(c*avg(*mon) for mon, c in zip(P.monoms(), P.coeffs())))

    def finite(expr, b):
        o = 3 - b
        xb = sgn[b]*Rs[b]; sig = 1 if b == 1 else -1
        groups = {}
        for term in sp.Add.make_args(sp.expand(expr)):
            pw = term.as_powers_dict()
            eb = pw.get(rs[b], 0); eo = pw.get(rs[o], 0)
            groups.setdefault((eb, eo), []).append(term/(rs[b]**eb*rs[o]**eo))
        total = 0
        for (eb, eo), lst in groups.items():
            C = sp.expand(sp.expand(sp.Add(*lst).subs({X: xb + s*nx, Y: s*ny, Z: s*nz}))*s**eb)
            if C == 0:
                continue
            pows = [w.as_powers_dict().get(s, 0) for w in sp.Add.make_args(C)]
            prod = sp.expand(C*taylor(int(eo), max(0, -min(pows)), sig))
            c0 = sp.Add(*[w for w in sp.Add.make_args(prod) if w.as_powers_dict().get(s, 0) == 0])
            total += angavg(c0)
        return sp.expand(sp.sympify(total).subs(Dsym, R1 + R2))

    def dX(expr):
        return sp.expand(sp.diff(expr, X) + sp.diff(expr, r1)*(X - R1)/r1 + sp.diff(expr, r2)*(X + R2)/r2)

    res = {}
    for b in (1, 2):
        G = sp.zeros(4, 4); dG = sp.zeros(4, 4)
        for i in range(4):
            for j in range(i, 4):
                G[i, j] = G[j, i] = finite(g[i, j], b)
                dG[i, j] = dG[j, i] = finite(dX(g[i, j]), b)
        h = [[tr(g[i+1, j+1], 4) for j in range(3)] for i in range(3)]
        T = lambda a, c: tr(sp.expand(a*c), 4)
        detg = (T(h[0][0], T(h[1][1], h[2][2]) - T(h[1][2], h[2][1]))
                - T(h[0][1], T(h[1][0], h[2][2]) - T(h[1][2], h[2][0]))
                + T(h[0][2], T(h[1][0], h[2][1]) - T(h[1][1], h[2][0])))
        mg = T(-tr(g[0, 0], 4), detg)
        sqmg = finite(tr(1 + (mg - 1)/2 - T(mg - 1, mg - 1)/8, 4), b)
        res[b] = (G, dG, sqmg)
    return res


# ---------------------------------------------------------------- circular orbit and the redshift
def redshift(m1, order=6, kappa=1, multiply=True, verbose=False, res_hook=None):
    """z1 coefficients of x^0..x^(order/2) from the rules, for body mass m1 (total mass 1)."""
    m1 = sp.nsimplify(m1); m2 = 1 - m1
    res = regularized_field(m1, order, kappa, multiply)
    if res_hook is not None:      # add extra regularized field pieces (round 14: the field-energy rule)
        res = res_hook(res)
    R1, R2 = sp.symbols('R1 R2', positive=True)
    p1, p2, q1, q2 = sp.symbols('p1 p2 q1 q2')
    a12, a14, a22, a24 = sp.symbols('a12 a14 a22 a24')
    d1s, d2s = sp.symbols('d1s d2s')

    def RS(x):
        x = sp.sympify(x); out = 0
        for i in range(3):
            for j in range(3 - i):
                der = sp.diff(x, R1, i, R2, j) if (i or j) else x
                out += der.subs({R1: m2, R2: m1})*d1s**i*d2s**j/(sp.factorial(i)*sp.factorial(j))
        return sp.expand(out.subs({d1s: e**2*a12 + e**4*a14, d2s: e**2*a22 + e**4*a24}))

    def series_fn(x, o, coef):
        u = tr(x - 1, o); out = 1; term = 1
        for k in range(1, o//2 + 1):
            term = tr(term*u, o); out += coef(k)*term
        return tr(out, o)
    psqrt = lambda x, o: series_fn(x, o, lambda k: sp.binomial(sp.Rational(1, 2), k))
    pinv = lambda x, o: series_fn(x, o, lambda k: (-1)**k)

    vel = {1: [1, 0, e*R1, 0], 2: [1, 0, -e*R2, 0]}
    eqs = {}; Ns = {}; Mu = {}
    for b in (1, 2):
        G, dG, sqmg = res[b]
        G = G.applyfunc(lambda u: tr(RS(u), order)); dG = dG.applyfunc(lambda u: tr(RS(u), order))
        sqmg = tr(RS(sqmg), order - 2)
        vb = [RS(x) for x in vel[b]]
        N = psqrt(-tr(sum(G[i, j]*vb[i]*vb[j] for i in range(4) for j in range(4)), order), order)
        Py = tr(sum(G[2, j]*vb[j] for j in range(4)), order)
        Fx = tr(sum(dG[i, j]*vb[i]*vb[j] for i in range(4) for j in range(4)), order)/2
        eqs[b] = tr(-e*Py - Fx, order)       # radial geodesic equation on the circle
        Ns[b] = N
        Mu[b] = tr((m1 if b == 1 else m2)*pinv(tr(sqmg*N, order - 2), order - 2), order - 2)
    M1 = Mu[1] - (m1 + e**2*p1 + e**4*q1); M2 = Mu[2] - (m2 + e**2*p2 + e**4*q2)
    E1, E2 = eqs[1], eqs[2]
    newt = [sp.simplify(E1.coeff(e, 2)), sp.simplify(E2.coeff(e, 2))]
    assert newt == [0, 0], newt
    sol = {}
    steps = [[E1.coeff(e, 4), E2.coeff(e, 4), M1.coeff(e, 2), M2.coeff(e, 2)], [a12, a22, p1, p2]],
    if order >= 6:
        steps = steps + ([[E1.coeff(e, 6), E2.coeff(e, 6), M1.coeff(e, 4), M2.coeff(e, 4)], [a14, a24, q1, q2]],)
    for eq, unk in steps:
        eq = [sp.simplify(x.subs(sol)) for x in eq]
        sol.update(sp.solve(eq, unk, dict=True)[0])
    for kk in (q1, q2, a14, a24, p1, p2, a12, a22):
        sol.setdefault(kk, 0)
    z1 = sp.expand(tr(Ns[1].subs(sol), order).subs(sol))
    cs = [sp.nsimplify(sp.simplify(z1.coeff(e, 2*k))) for k in range(order//2 + 1)]
    if verbose:
        print('m1 =', m1, 'rules:', cs, 'truth:', truth(m1)[:len(cs)], flush=True)
    return cs


def mismatch_formula(m1):
    """Fitted rules - truth in the x^3 coefficient (m = 1): m1 m2 (23 m1^2 - 50 m1 + 36)/3."""
    m1 = sp.nsimplify(m1)
    return m1*(1 - m1)*(23*m1**2 - 50*m1 + 36)/3


if __name__ == '__main__':
    t0 = time.time()
    fit_masses = [sp.Integer(0), sp.Rational(1, 5), sp.Rational(1, 4), sp.Rational(1, 3),
                  sp.Rational(1, 2), sp.Rational(2, 3)]
    out = {'one_pn': {}, 'two_pn': {}, 'kappa0': {}}
    out['one_pn']['1/3'] = [str(c) for c in redshift(sp.Rational(1, 3), order=4, verbose=True)]
    rows = {}
    for mm in fit_masses + [sp.Rational(2, 5)]:
        cs = redshift(mm, order=6, verbose=True)
        rows[mm] = cs
        out['two_pn'][str(mm)] = {'rules': [str(c) for c in cs],
                                  'truth': [str(c) for c in truth(mm)],
                                  'difference_x3': str(cs[3] - truth(mm)[3])}
    out['kappa0']['1/3'] = [str(c) for c in redshift(sp.Rational(1, 3), order=6, kappa=0, verbose=True)]
    # rule B tested on the motion: add the clock rates instead of multiplying them
    out['added_rates'] = {'1PN_1/3': [str(c) for c in redshift(sp.Rational(1, 3), order=4, multiply=False, verbose=True)],
                          '2PN_1/3': [str(c) for c in redshift(sp.Rational(1, 3), order=6, multiply=False, verbose=True)]}
    # fit a quartic in m1 to the x^3 difference on six masses, then test it on m1 = 2/5
    a = sp.symbols('a0:5'); mm = sp.symbols('mm')
    poly = sum(a[k]*mm**k for k in range(5))
    fit = sp.solve([poly.subs(mm, k) - (rows[k][3] - truth(k)[3]) for k in fit_masses], a, dict=True)[0]
    fitted = sp.factor(poly.subs(fit))
    out['fit'] = str(fitted)
    out['fit_matches_formula'] = bool(sp.simplify(fitted - mismatch_formula(mm)) == 0)
    k = sp.Rational(2, 5)
    out['held_out_2_5'] = {'measured': str(rows[k][3] - truth(k)[3]),
                           'formula': str(mismatch_formula(k))}
    out['seconds'] = round(time.time() - t0)
    json.dump(out, open('next_order.json', 'w'), indent=1)
    print(json.dumps(out, indent=1))
