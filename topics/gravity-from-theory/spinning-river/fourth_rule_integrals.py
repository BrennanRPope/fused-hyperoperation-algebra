"""Round 14 (DRAFT): the static 'field energy acting on field energy' pieces of GR's 2PN potential X-hat.

GR's 2PN time-time potential X-hat has a source term W-hat_ij d_ij V: the 1PN space-space field W-hat (itself made
of field energy) acting through the curvature of the Newtonian potential V. With both pieces static, four products
appear at third order in the masses (units r12 = 1, body 1 at y1, body 2 at y2):

  T1 = W_self1 : dd V2   (m1^2 m2)      T2 = W_self2 : dd V1   (m1 m2^2)
  T3 = W_cross : dd V1   (m1^2 m2)      T4 = W_cross : dd V2   (m1 m2^2)

with W_self = m^2 (n n - delta)/(4 r^2) and W_cross = -m1 m2 (A + A^T), A_ij = d/dy1_i d/dy2_j ln(r1 + r2 + r12).
X-hat at body 1 is the Poisson integral  X(y1) = -(1/4 pi) Pf Int F / r1,  and its x-gradient uses c1/r1^2.
Result:  X-hat(y1) = (alpha m1^2 m2 + beta m1 m2^2)/r12^3,  dX/dx(y1) = (gamma m1^2 m2 + delta m1 m2^2)/r12^4.

Method:
  * bipolar coordinates, d^3x = 2 pi r1 r2 dr1 dr2 / r12: the inner r2 integral exactly (sympy), the outer r1
    integral numerically after subtracting its singular series at each end, whose finite part is added exactly;
  * this orders the integral around body 1, so the cut-out near body 2 is a thin slab, not a ball. For a
    leading term A(u)/rho^3 near body 2 (u = n2.n12), ball minus slab is  a2/3 + b/4  in X units, where
    A = c0 + b u + c2 u^2 and a2 = 2 c2/3 (checked on test integrands with known answers);
  * d_ij(1/r) = Pf(3 n_i n_j - delta_ij)/r^3 - (4 pi/3) delta_ij delta^3: the point part adds
    (1/3) m_b W_kk(y_b) times the kernel at y_b, with W_kk(y_b) taken as its direction average (Hadamard).

Answer: alpha = -1/2, beta = 1, gamma = 0, delta = -2.  About 1-2 minutes.
Cross-checks (python3 fourth_rule_integrals.py --cross, a few more minutes): T4 re-ordered around body 2, and T3's
value by direct 3D integration with a small ball cut out at body 1.
"""
import sympy as sp, mpmath as mp, json, os, sys
HERE = os.path.dirname(os.path.abspath(__file__))
r1, r2, t, eps, u = sp.symbols('r1 r2 t eps u', positive=True)
D = 1
S = r1 + r2 + D
c1 = (r2**2 - r1**2 - D**2)/(2*r1*D)       # n1.n12
c2 = (r2**2 + D**2 - r1**2)/(2*r2*D)       # n2.n12
c12 = (r1**2 + r2**2 - D**2)/(2*r1*r2)     # n1.n2
Akk = -2/(D*S) + (c2 - c12 + 1 - c1)/S**2
n1An1 = -(1 - c1**2)/(D*S) + (c12 + c1)*(c1 - 1)/S**2
n2An2 = -(1 - c2**2)/(D*S) + (1 + c2)*(c2 - c12)/S**2
# F for each product, after contracting with the trace-free part 3 n n - delta of d_ij(1/r_b)
P = {'T1': (3*c12**2 - 1)/(4*r1**2*r2**3),
     'T2': (3*c12**2 - 1)/(4*r2**2*r1**3),
     'T3': -(6*n1An1 - 2*Akk)/r1**3,
     'T4': -(6*n2An2 - 2*Akk)/r2**3}
MASS = {'T1': 'm1^2 m2', 'T2': 'm1 m2^2', 'T3': 'm1^2 m2', 'T4': 'm1 m2^2'}
KERN = {'val': sp.Integer(1), 'grad': c1/r1}          # times 1/r1 from the bipolar measure
mp.mp.dps = 60


def _inner(expr, var, lo, hi):
    e = sp.apart(sp.cancel(sp.together(expr)), var)
    Fi = sp.integrate(e, var)
    return sp.simplify(Fi.subs(var, hi) - Fi.subs(var, lo))


def _neg_power(term):
    q = sp.powsimp(term.subs(sp.log(t), 1))
    return (q*t**8).is_polynomial(t) and sp.Poly(q*t**8, t).monoms()[0][0] - 8 < 0


def _pf_piece(f, var, a, b, end, side):
    """Hadamard finite part of Int_a^b f d(var), singular at var = end (approached from side +1 or -1)."""
    ser = None
    for n in (1, 2, 3):
        try:
            ser = sp.expand(sp.series(f.subs(var, end + side*t), t, 0, n).removeO()); break
        except NotImplementedError:
            pass
    sing = sp.Add(*[term for term in sp.Add.make_args(ser) if _neg_power(term)])
    g = f - sing.subs(t, side*(var - end))
    num = mp.quad(sp.lambdify(var, g, 'mpmath'), mp.linspace(a, b, 9), method='gauss-legendre', maxdegree=10)
    Fi = sp.integrate(sing, (t, eps, abs(b - a)))
    fin = sp.expand(sp.series(Fi, eps, 0, 1).removeO()).as_independent(eps)[0]
    return num + mp.mpf(sp.N(fin.subs(sp.log(eps), 0), 40))


def _bipolar(expr, outer):
    """-(1/2) Pf Int Int expr r2 dr1 dr2, outer variable r1 (cut-outs: ball at y1, slab at y2) or r2 (the reverse)."""
    o, i = (r1, r2) if outer == 'r1' else (r2, r1)
    low = _inner(expr*r2, i, 1 - o, o + 1); high = _inner(expr*r2, i, o - 1, o + 1)
    tot = (_pf_piece(low, o, 0, sp.Rational(1, 2), 0, 1) + _pf_piece(low, o, sp.Rational(1, 2), 1, 1, -1)
           + _pf_piece(high, o, 1, 2, 1, 1)
           + mp.quad(sp.lambdify(o, high, 'mpmath'), [2, 10, 100, mp.inf], method='gauss-legendre'))
    return -tot/2


def slab_to_ball(expr3d):
    """Ball-minus-slab correction at body 2 for a 3D integrand expr3d (already including 1/r1 or c1/r1^2)."""
    f = expr3d.subs(r1, sp.sqrt(r2**2 + 1 - 2*r2*u))
    lead = sp.expand(sp.series(f*r2**3, r2, 0, 1).removeO().subs(r2, 0))
    p = sp.Poly(lead, u)
    if p.degree() > 2: raise ValueError('leading angular term beyond l = 2: correction rule not derived')
    cc = [p.coeff_monomial(u**k) for k in range(3)]
    return sp.Rational(2, 3)*cc[2]/3 + cc[1]/4, lead


def point_parts():
    """(1/3) m_b W_kk(y_b) x kernel(y_b) from the -(4 pi/3) delta part of d_ij(1/r_b). Only y2 contributes at y1."""
    # W_self1_kk at y2 = m1^2 (1 - 3)/(4 r12^2);  W_cross_kk = -2 m1 m2 A_kk, direction-averaged at y2
    Wself1_kk = sp.Rational(-1, 2)
    lim = sp.limit(Akk.subs(r1, sp.sqrt(r2**2 + 1 - 2*r2*u)), r2, 0)
    Wc_kk = sp.integrate(-2*lim, (u, -1, 1))/2            # direction average
    k_val, k_grad = 1, -1                                   # kernels 1/r1 and c1/r1^2 at y2 (c1 = -1 there)
    third = sp.Rational(1, 3)
    return {('T1', 'val'): third*Wself1_kk*k_val, ('T1', 'grad'): third*Wself1_kk*k_grad,
            ('T4', 'val'): third*Wc_kk*k_val, ('T4', 'grad'): third*Wc_kk*k_grad,
            ('T2', 'val'): 0, ('T2', 'grad'): 0, ('T3', 'val'): 0, ('T3', 'grad'): 0}, Wc_kk


def compute(verbose=True):
    pts, Wc_kk = point_parts()
    rows = {}
    for name in ('T1', 'T2', 'T3', 'T4'):
        for kn in ('val', 'grad'):
            shell = _bipolar(P[name]*KERN[kn], 'r1')
            corr, lead = slab_to_ball(P[name]*KERN[kn]/r1) if name in ('T1', 'T4') else (sp.Integer(0), 0)
            ball = shell + mp.mpf(sp.N(corr, 40))
            ballq = sp.nsimplify(float(ball), tolerance=1e-6, rational=True)
            tot = ballq + pts[(name, kn)]
            rows[f'{name}_{kn}'] = {'mass': MASS[name], 'shell': mp.nstr(shell, 12), 'slab_to_ball': str(corr),
                                    'leading_near_body2': str(lead), 'ball_finite_part': str(ballq),
                                    'ball_numeric': mp.nstr(ball, 12), 'point_part': str(pts[(name, kn)]),
                                    'total': str(tot)}
            if verbose: print(name, kn, rows[f'{name}_{kn}'], flush=True)
    g = lambda k: sp.Rational(rows[k]['total'])
    coeff = {'alpha': g('T1_val') + g('T3_val'), 'beta': g('T2_val') + g('T4_val'),
             'gamma': g('T1_grad') + g('T3_grad'), 'delta': g('T2_grad') + g('T4_grad')}
    return {'rows': rows, 'Wcross_kk_at_body2': str(Wc_kk), 'coefficients': {k: str(v) for k, v in coeff.items()}}


def cross_checks():
    out = {}
    # T4 ordered around body 2: ball at y2 directly, no slab correction needed there
    for kn in ('val', 'grad'):
        v = _bipolar(P['T4']*KERN[kn], 'r2')
        out[f'T4_{kn}_ordered_at_body2'] = mp.nstr(v, 12)
    # T3 by direct spherical integration around body 1 with a ball of radius 1e-3 cut out
    r, c = sp.symbols('r c', real=True)
    sub = {r1: r, r2: sp.sqrt(1 + 2*r*c + r**2)}
    mp.mp.dps = 30
    for kn, kk in (('val', 1/r),):        # the gradient's 1/r^3 piece needs the finite part, not a hard cut
        f = sp.lambdify((r, c), (P['T3']*kk).subs(sub)*r**2, 'mpmath')
        I = 2*mp.pi*mp.quad(lambda rr: mp.quad(lambda cc: f(rr, cc), [-1, 0, 1]),
                            [mp.mpf('1e-3'), 0.01, 0.1, 0.5, 0.9, 0.99, 1, 1.01, 1.1, 2, 10, mp.inf])
        out[f'T3_{kn}_direct_ball'] = mp.nstr(-I/(4*mp.pi), 10)
    mp.mp.dps = 60
    return out


if __name__ == '__main__':
    res = compute()
    if '--cross' in sys.argv:
        res['cross_checks'] = cross_checks(); print(res['cross_checks'])
    print('coefficients:', res['coefficients'])
    json.dump(res, open(os.path.join(HERE, 'fourth_rule_integrals.json'), 'w'), indent=1)
