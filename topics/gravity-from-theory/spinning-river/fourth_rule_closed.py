"""Round 14 (DRAFT): the closed-form part of what general relativity's 2PN field has and the three rules lack.

D = GR - rules, for two bodies on a circle (total mass 1, units of next_order.py):
  g_ij += 4 Wc ;  g_0i += -8 Rcross ;  g_00 += 8 Xcross + 16 V1k V2k + 8 Xacc_self (+ the static field-on-field
  pieces, done by integration in fourth_rule_integrals.py).
Every closed-form piece is a derivative, with respect to a body's position, of Fock's function g = ln(r1 + r2 + r12),
whose Laplacian is 1/(r1 r2). finite() takes Hadamard's finite part at a body, as next_order.py does.
"""
import sympy as sp, functools
from next_order import e
X, Y, Z, s = sp.symbols('X Y Z s', real=True)
nx, ny, nz = sp.symbols('nx ny nz', real=True)
R1, R2 = sp.symbols('R1 R2', positive=True)
r1s, r2s, Ds = sp.symbols('r1 r2 Dsym', positive=True)
P1 = sp.symbols('p1x p1y p1z', real=True); P2 = sp.symbols('p2x p2y p2z', real=True)
xv = [X, Y, Z]
_r1 = sp.sqrt(sum((xv[i] - P1[i])**2 for i in range(3)))
_r2 = sp.sqrt(sum((xv[i] - P2[i])**2 for i in range(3)))
_r12 = sp.sqrt(sum((P1[i] - P2[i])**2 for i in range(3)))
g = sp.log(_r1 + _r2 + _r12)
CONF = {P1[0]: R1, P1[1]: 0, P1[2]: 0, P2[0]: -R2, P2[1]: 0, P2[2]: 0}
_B1 = sp.expand((X - R1)**2 + Y**2 + Z**2); _B2 = sp.expand((X + R2)**2 + Y**2 + Z**2); _B12 = sp.expand((R1 + R2)**2)


def conf(expr):
    ex = expr.subs(CONF)
    def rep(u):
        b = sp.expand(u.base)
        if b == _B1: return r1s**(2*u.exp)
        if b == _B2: return r2s**(2*u.exp)
        if b == _B12: return Ds**(2*u.exp)
        return u
    ex = ex.replace(lambda u: u.is_Pow and not u.exp.is_Integer, rep)
    ex = ex.replace(lambda u: u.is_Pow and u.exp.is_Integer and u.exp < 0 and sp.expand(u.base) in (_B1, _B2, _B12), rep)
    return ex


@functools.lru_cache(None)
def _dg(idx):
    ex = g
    for who, i in idx:
        ex = sp.diff(ex, (P1 if who == '1' else P2)[i])
    return ex


def D(*idx):
    """derivative of g with respect to body positions, e.g. D(('1',0),('2',1)) = d^2 g / dy1_x dy2_y."""
    return conf(_dg(tuple(idx)))


def build(m1, m2):
    """The closed-form pieces of D for the circular orbit (positions y1 = (R1,0,0), y2 = (-R2,0,0) at t = 0)."""
    M1, M2 = e**2*m1, e**2*m2
    v1 = [0, e*R1, 0]; v2 = [0, -e*R2, 0]; a1 = [-e**2*R1, 0, 0]; a2 = [e**2*R2, 0, 0]
    rng = range(3); dl = lambda i, j: 1 if i == j else 0
    Wc = [[-M1*M2*(D(('1', i), ('2', j)) + D(('1', j), ('2', i))) for j in rng] for i in rng]
    Rc = []
    for i in rng:
        comp = M1*M2*(v1[i] - v2[i])/(Ds*r1s) + M1*M2*(v2[i] - v1[i])/(Ds*r2s)
        nl1 = -2*M1*M2*sum(v2[k]*D(('1', k), ('2', i)) + v1[k]*D(('1', i), ('2', k)) for k in rng)
        nl2 = sp.Rational(3, 2)*M1*M2*sum(v1[k]*D(('1', k), ('2', i)) + v2[k]*D(('2', k), ('1', i)) for k in rng)
        Rc.append(comp + nl1 + nl2)
    vsq1 = sum(q**2 for q in v1); vsq2 = sum(q**2 for q in v2)
    P = dict(Wc=Wc, Rc=Rc)
    P['Xa'] = M1*M2*vsq1/(Ds*r1s) + M1*M2*vsq2/(Ds*r2s)                                   # -4 pi V sigma_ii
    P['Xb'] = -2*M1*M2*sum(v1[i]*v2[k]*D(('2', i), ('2', k)) + v2[i]*v1[k]*D(('1', i), ('1', k)) for i in rng for k in rng)
    P['Xc'] = M1*M2*(sum(a2[k]*D(('2', k)) for k in rng) + sum(v2[k]*v2[l]*D(('2', k), ('2', l)) for k in rng for l in rng)
                     + sum(a1[k]*D(('1', k)) for k in rng) + sum(v1[k]*v1[l]*D(('1', k), ('1', l)) for k in rng for l in rng))
    P['Xd'] = 3*M1*M2*sum(v1[k]*v2[l]*D(('1', k), ('2', l)) for k in rng for l in rng)       # (3/2)(dV/dt)^2
    P['Xe'] = -4*M1*M2*sum(v1[j]*v2[i]*D(('1', i), ('2', j)) for i in rng for j in rng)      # -2 dVj dVi
    P['Xf1c'] = M1*M2*sum((v1[i]*v1[j] - dl(i, j)*vsq1)*D(('2', i), ('2', j)) for i in rng for j in rng)
    P['Xf2c'] = M1*M2*sum((v2[i]*v2[j] - dl(i, j)*vsq2)*D(('1', i), ('1', j)) for i in rng for j in rng)
    P['Xacc'] = (-(M1**2/2)*sum(a1[k]*(xv[k] - [R1, 0, 0][k]) for k in rng)/r1s**2
                 - (M2**2/2)*sum(a2[k]*(xv[k] - [-R2, 0, 0][k]) for k in rng)/r2s**2)          # speeding-up body
    P['VV'] = M1*M2*sum(v1[k]*v2[k] for k in rng)/(r1s*r2s)
    return P


@functools.lru_cache(None)
def _avg(a, b, c):
    if a % 2 or b % 2 or c % 2: return 0
    df = lambda k: sp.factorial2(k - 1) if k > 0 else 1
    return sp.Rational(df(a)*df(b)*df(c), sp.factorial2(a + b + c + 1))


def _angavg(poly):
    Pl = sp.Poly(sp.expand(poly), nx, ny, nz)
    return sp.expand(sum(c*_avg(*mon) for mon, c in zip(Pl.monoms(), Pl.coeffs())))


def dX(expr):
    return sp.diff(expr, X) + sp.diff(expr, r1s)*(X - R1)/r1s + sp.diff(expr, r2s)*(X + R2)/r2s


def finite(expr, b):
    """Hadamard finite part at body b: the direction average of the r^0 term."""
    if b == 1: sub = {X: R1 + s*nx, Y: s*ny, Z: s*nz}; rb, ro, sig = r1s, r2s, 1
    else: sub = {X: -R2 + s*nx, Y: s*ny, Z: s*nz}; rb, ro, sig = r2s, r1s, -1
    ex = sp.sympify(expr).subs(sub)
    ex = ex.subs(R1, Ds - R2) if b == 1 else ex.subs(R2, Ds - R1)
    ex = ex.subs({rb: s, ro: Ds*sp.sqrt(1 + 2*sig*s*nx/Ds + s**2/Ds**2)})
    total = 0
    for term in sp.Add.make_args(sp.expand(ex)):
        c0 = sp.expand(sp.series(term, s, 0, 1).removeO()).coeff(s, 0)
        if c0 != 0: total += _angavg(c0)
    return sp.expand(sp.expand(total).subs(Ds, R1 + R2))


def at_bodies(m1):
    """(value, x-gradient) of every closed-form piece at each body. About 6 minutes."""
    m1 = sp.nsimplify(m1); P = build(m1, 1 - m1); out = {}
    for b in (1, 2):
        res = {}
        both = lambda name, ex: res.__setitem__(name, (finite(ex, b), finite(dX(ex), b)))
        for i in range(3):
            for j in range(i, 3): both(('Wc', i, j), P['Wc'][i][j])
            both(('Rc', i), P['Rc'][i])
        for k in ('Xa', 'Xb', 'Xc', 'Xd', 'Xe', 'Xf1c', 'Xf2c', 'Xacc', 'VV'): both(k, P[k])
        out[b] = res
    return out
