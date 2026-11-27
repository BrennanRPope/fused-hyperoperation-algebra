"""Second-order (2PN) periastron advance for two bodies, from Gupta's Einstein-Maxwell-dilaton Lagrangian
(arXiv:2609.07413, Eqs. 4.2, 5.9, 7.1; transcribed in sources/gupta-2609.07413.md on branch sources-gupta).
Route: order-reduce the accelerations with the Newtonian equations of motion, Legendre-transform to H,
centre-of-mass frame = total momentum zero, solve H = E for p_r^2 as a polynomial in u = 1/r, radial action
i_r to second order around Kepler (finite parts by d/dA), turn = 2 pi (-d i_r/dJ - 1).
Reported as turn/pi = T k^2/J^2 + (T2 + T2e w) k^4/J^4 with w = -E J^2/(mu k^2), as for a small visitor."""
import os, sympy as sp
HERE = os.path.dirname(os.path.abspath(__file__))
s = sp.symbols('s', positive=True)            # PN bookkeeping: 1/c
G, m1, m2, q1, q2, a, r = sp.symbols('G m1 m2 q1 q2 a r', positive=True)
al1, al2, be1, be2, bb1, bb2, d1, d2, e1, e2, z12, gb, G12 = sp.symbols('al1 al2 be1 be2 bb1 bb2 d1 d2 e1 e2 z12 gb G12')
names = 'v1v1 v2v2 v1v2 nv1 nv2 na1 na2 a1a2 v1a2 a1v2 v1a1 v2a2'.split()
K = {n: sp.Symbol(n) for n in names}
L1s, L2s = sp.symbols('L1 L2')
def halves():
    ns = dict(G=G, m1=m1, m2=m2, q1=q1, q2=q2, a=a, r=r, al1=al1, al2=al2, be1=be1, bb1=bb1, bb2=bb2, d2=d2, e1=e1,
              z12=z12, gb=gb, G12=G12, L1=L1s, L2=L2s, **K)
    out = {}
    for line in open(os.path.join(HERE, 'gupta_halves.txt')):
        if '=' in line and not line.startswith('#'):
            k, v = line.split('=', 1); out[k.strip()] = sp.sympify(v.strip(), locals=ns)
    return out['ell_1PN'], out['ell_2PN']
def image(ex):
    ex = ex.subs({L1s: 7 + 10*gb + sp.Rational(7, 2)*gb**2 - 2*d2 + 2*bb1, L2s: 8 + 10*gb + sp.Rational(7, 2)*gb**2 - 2*d2 + 4*bb1})
    sw = {m1: m2, m2: m1, q1: q2, q2: q1, al1: al2, al2: al1, be1: be2, be2: be1, bb1: bb2, bb2: bb1, d1: d2, d2: d1, e1: e2, e2: e1,
          K['v1v1']: K['v2v2'], K['v2v2']: K['v1v1'], K['nv1']: -K['nv2'], K['nv2']: -K['nv1'], K['na1']: -K['na2'], K['na2']: -K['na1'],
          K['v1a2']: K['a1v2'], K['a1v2']: K['v1a2'], K['v1a1']: K['v2a2'], K['v2a2']: K['v1a1']}
    return ex, ex.subs(sw, simultaneous=True)
def full_L():
    l1, l2 = halves(); A1, B1 = image(l1); A2, B2 = image(l2)
    return A1 + B1, A2 + B2
def params(case, **kw):
    """Substitutions for the body parameters. case 'GR' or 'EMD-extremal-opposite' (coupling a, alpha = a, beta = beta' = 0,
    q^2 = G m^2 (1 + a^2), opposite charges)."""
    if case == 'GR':
        return {q1: 0, q2: 0, al1: 0, al2: 0, be1: 0, be2: 0, bb1: 0, bb2: 0, d1: 0, d2: 0, e1: 0, e2: 0, z12: 0, gb: 0, G12: G, a: 0}
    av = sp.sympify(kw['a']); f = 1 + av**2
    return {a: av, al1: av, al2: av, be1: 0, be2: 0, bb1: 0, bb2: 0, e1: 0, e2: 0, z12: 0, d1: av**2/f**2, d2: av**2/f**2,
            gb: -2*av**2/f, G12: G*f, q1: sp.sqrt(G*f)*m1, q2: -sp.sqrt(G*f)*m2}
def hamiltonian(L1, L2, sub):
    """H(r, pr, J) in the centre-of-mass frame (p1 = p, p2 = -p), 2D, n = x axis."""
    L1 = L1.subs(sub); L2 = L2.subs(sub)
    k = sp.nsimplify(sp.expand((G12*m1*m2 - q1*q2).subs(sub)))
    v1x, v1y, v2x, v2y = sp.symbols('v1x v1y v2x v2y')
    a1x, a2x = -k/(m1*r**2), k/(m2*r**2)            # Newtonian accelerations along n
    comp = {K['v1v1']: v1x**2 + v1y**2, K['v2v2']: v2x**2 + v2y**2, K['v1v2']: v1x*v2x + v1y*v2y, K['nv1']: v1x, K['nv2']: v2x,
            K['na1']: a1x, K['na2']: a2x, K['a1a2']: a1x*a2x, K['v1a2']: v1x*a2x, K['a1v2']: a1x*v2x, K['v1a1']: v1x*a1x, K['v2a2']: v2x*a2x}
    L1c = sp.expand(L1.subs(comp)); L2c = sp.expand(L2.subs(comp))
    px, py = sp.symbols('px py')
    at = {v1x: px/m1, v1y: py/m1, v2x: -px/m2, v2y: -py/m2}
    cross = sum((sp.diff(L1c, v)**2) for v in (v1x, v1y))/(2*m1) + sum((sp.diff(L1c, v)**2) for v in (v2x, v2y))/(2*m2)
    H0 = (px**2 + py**2)*(1/m1 + 1/m2)/2 - k/r
    H1 = -L1c.subs(at); H2 = (-L2c + cross).subs(at)
    ms = {m1: sub.get(m1, m1), m2: sub.get(m2, m2)}
    return H0.subs(ms), sp.expand(H1.subs(ms)), sp.expand(H2.subs(ms)), k.subs(ms)
def turn_coeffs(H0, H1, H2, k, mu):
    pr, u, E = sp.symbols('pr u E'); J = sp.Symbol('J', positive=True)
    X = sp.Symbol('X')                                 # X = pr^2
    def toX(H):
        H = sp.expand(H.subs({sp.Symbol('py'): J*u, r: 1/u}))
        px = sp.Symbol('px'); H = sp.expand(H)
        assert all(m % 2 == 0 for m in sp.Poly(H, px).monoms()[0]) if H.has(px) else True
        return sp.expand(H.subs(px, sp.sqrt(X)))
    h1, h2 = toX(H1), toX(H2)
    X0 = 2*mu*(E + k*u) - J**2*u**2
    X1 = -2*mu*h1.subs(X, X0)
    X2 = -2*mu*(h2.subs(X, X0) + sp.diff(h1, X).subs(X, X0)*X1)
    P = sp.expand(X0 + s**2*X1 + s**4*X2)
    c = [P.coeff(u, j) for j in range(sp.degree(P, u) + 1)]
    A, B, C = sp.symbols('A B C')
    pp, q2s = B/(-C), (B/(-C))**2 + A/(-C)
    def I(kk):  # integral_0^pi (p + q cos)^k
        return sp.pi*sum(sp.binomial(kk, j)*pp**(kk - j)*q2s**(j//2)*sp.factorial2(j - 1)/sp.factorial2(j) for j in range(0, kk + 1, 2))
    d = {j: c[j] for j in range(3, len(c))}
    d1 = {j: sp.expand(cj).coeff(s, 2) for j, cj in d.items()}
    first = sum(cj*I(j - 2) for j, cj in d.items())/(2*sp.pi*sp.sqrt(-C))
    sq = sum(d1[i]*d1[j]*I(i + j - 2) for i in d1 for j in d1)
    second = s**4*sp.diff(sq/sp.sqrt(-C), A)/(4*sp.pi)
    ir = B/sp.sqrt(-A) - sp.sqrt(-C) + first + second
    ir = ir.subs({A: c[0], B: c[1]/2, C: c[2]})
    Kad = -sp.diff(ir, J)
    ser = sp.series(Kad, s, 0, 5).removeO()
    t1 = sp.cancel(sp.expand(2*ser.coeff(s, 2)*J**2/k**2))
    t2 = sp.expand(2*ser.coeff(s, 4)*J**4/k**4)
    # E enters as w = -E J^2/(mu k^2)
    w = sp.Symbol('w'); t2w = sp.expand(sp.cancel(sp.expand(t2.subs(E, -w*mu*k**2/J**2))))
    return sp.nsimplify(t1), sp.nsimplify(t2w.coeff(w, 0)), sp.nsimplify(t2w.coeff(w, 1)), sp.expand(t2w - t2w.coeff(w, 0) - w*t2w.coeff(w, 1))
def run(case, M1=1, M2=1, **kw):
    sub = params(case, **kw); sub = {kk: (v.subs({G: 1, m1: M1, m2: M2}) if hasattr(v, 'subs') else v) for kk, v in sub.items()}
    sub.update({G: 1, m1: M1, m2: M2})
    L1, L2 = full_L(); H0, H1, H2, k = hamiltonian(L1, L2, sub)
    T, T2, T2e, rest = turn_coeffs(H0, H1, H2, k, sp.Rational(M1*M2)/(M1 + M2) if isinstance(M1, int) else M1*M2/(M1 + M2))
    nu = sp.Symbol('nu')
    return T, T2, T2e, rest
if __name__ == "__main__":
    import time; t = time.time()
    print("GR equal:", run('GR'), time.time() - t)
