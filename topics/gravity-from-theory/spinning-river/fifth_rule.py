"""DRAFT. Round 31: a fifth rule, gravity's field energy acting on its own field energy. Does it make the pi^2?

Part 24 showed the four rules can only give rational numbers at the third level, while relativity's third-level
energy has -(205 pi^2/96) nu x^4. Here the question is where a pi^2 can come from at all, and what the fifth rule
needs to make relativity's.

1. Census. Every two-body field piece reduces to overlaps (1/4pi) int d^3x r1^a r2^b = R(a,b) r12^(a+b+3)
   (Riesz's formula, continued in a and b as dimensional and Hadamard regularization do when no logarithm appears).
   R is rational whenever a or b is odd, and the only finite even-even value that is not zero is
   R(-2,-2) = pi^2/4: the overlap of body 1's field squared with body 2's field squared. A field squared is
   field energy's own size (U^2 ~ (grad U)^2 r^2), so pi^2 first appears when field energy acts on field energy,
   with both masses squared (G^3 m1^2 m2^2), which is exactly the slot relativity's pi^2 sits in (its energy term is
   proportional to mu nu = m1^2 m2^2/m^3, and through the first law it is first order in a small body's mass).
2. Static. With the bodies at rest, field energy acting on field energy, int U^2 (grad U)^2, has a zero pi^2 part:
   its m1^2 m2^2 overlaps are R(-4,-2), R(-2,-4) (both 0).
3. Moving. With the bodies moving, field energy's time part, int U^2 (dU/dt)^2, does make pi^2:
   (pi^2/4)(G^3 m1^2 m2^2/r^3)[ alpha v1^2 + beta v2^2 + gamma v1.v2 + (n.v) terms ],
   alpha = beta = -1/8, gamma = 1. Relativity's pi^2 can depend on the bodies' velocities only through
   v1 - v2 (moving the whole pair must not change it, up to a total time derivative: 2 alpha + gamma = 0 = 2 beta + gamma),
   which on circles puts pi^2 in the nu^1 term alone. The time part alone fails that (2 alpha + gamma = 3/4): it
   would also put pi^2 into a nu^2 term that relativity doesn't have. The light-speed lag of field energy (rule A applied
   to field energy) is what has to restore it.
4. Size. With the boost-safe shape, a circle sees only one number K: L_pi = K (pi^2/4) G^3 m1^2 m2^2 v12^2 / r^3.
   Relativity's -(205 pi^2/96) nu needs K = K_GR (computed below).

    python3 fifth_rule.py       # writes fifth_rule.json (seconds)
"""
import json, os
import sympy as sp, mpmath as mp

HERE = os.path.dirname(os.path.abspath(__file__))
a_, b_, t_ = sp.symbols('a b t')


def riesz(A, B):
    """(1/4pi) int d^3x r1^A r2^B = riesz(A,B) r12^(A+B+3); oo marks a pole (a logarithm in regularization)."""
    e = (sp.pi**sp.Rational(3, 2)*sp.gamma((a_ + 3)/2)*sp.gamma((b_ + 3)/2)*sp.gamma(-(a_ + b_ + 3)/2)
         / (sp.gamma(-a_/2)*sp.gamma(-b_/2)*sp.gamma((a_ + b_ + 6)/2))/(4*sp.pi))
    return sp.simplify(sp.limit(e.subs({a_: A + t_, b_: B + 2*t_}), t_, 0))


def census(lo=-6):
    out = {}
    for A in range(lo, 0):
        for B in range(lo, 0):
            v = riesz(A, B)
            kind = 'pole' if v.has(sp.oo, -sp.oo, sp.zoo) else ('pi^2' if v.has(sp.pi) else 'rational')
            out[f'{A},{B}'] = (str(v), kind)
    return out


def numeric_overlap():
    """int d^3x /(r1^2 r2^2) at r12 = 2 in prolate spheroidal coordinates; should be pi^3/2."""
    mp.mp.dps = 20
    val = mp.quad(lambda s: mp.quad(lambda t: 2*mp.pi/(s*s - t*t), [-1, 1]), [1, mp.inf])
    return float(val), float(mp.pi**3/2)


# positions and velocities as symbols
y1 = sp.Matrix(sp.symbols('y1x y1y y1z', real=True)); y2 = sp.Matrix(sp.symbols('y2x y2y y2z', real=True))
v1 = sp.Matrix(sp.symbols('v1x v1y v1z', real=True)); v2 = sp.Matrix(sp.symbols('v2x v2y v2z', real=True))
r12 = sp.sqrt(((y1 - y2).T*(y1 - y2))[0])
F = sp.pi**2/4/r12                       # (1/4pi) int r1^-2 r2^-2


def vd(v, y, f):
    return sum(v[i]*sp.diff(f, y[i]) for i in range(3))


def identity_check():
    """(v.n)^2/r^4 = [ (v.grad_y)^2 r^-2 + 2 v^2 r^-4 ] / 8, with r = |x - y|, n = (x - y)/r."""
    x = sp.Matrix(sp.symbols('x y z', real=True)); r = sp.sqrt(((x - y2).T*(x - y2))[0])
    n = (x - y2)/r
    lhs = (v2.T*n)[0]**2/r**4
    rhs = (vd(v2, y2, vd(v2, y2, 1/r**2)) + 2*(v2.T*v2)[0]/r**4)/8
    return sp.simplify(lhs - rhs) == 0


def moving_time_part():
    """pi^2 part of (1/4pi) int U^2 (dU/dt)^2, m1^2 m2^2 terms, U = m1/r1 + m2/r2, dU_a/dt = v_a . grad_{y_a} U_a.
    U1^2 (dU2/dt)^2 -> (v2.grad2)^2 F / 8 (the 2 v^2 R(-2,-4) part is zero); U2^2 (dU1/dt)^2 likewise;
    4 U1 U2 dU1/dt dU2/dt = (v1.grad1)(v2.grad2)(r1^-2 r2^-2) -> (v1.grad1)(v2.grad2) F."""
    T = vd(v2, y2, vd(v2, y2, F))/8 + vd(v1, y1, vd(v1, y1, F))/8 + vd(v1, y1, vd(v2, y2, F))
    # circular orbit at an instant: separation along x, velocities along y
    R_ = sp.symbols('R', positive=True); w1, w2 = sp.symbols('w1 w2', real=True)
    circ = {y1[0]: R_, y1[1]: 0, y1[2]: 0, y2[0]: 0, y2[1]: 0, y2[2]: 0,
            v1[0]: 0, v1[1]: w1, v1[2]: 0, v2[0]: 0, v2[1]: w2, v2[2]: 0}
    Tc = sp.simplify(T.subs(circ)*R_**3/(sp.pi**2/4))
    P = sp.Poly(sp.expand(Tc), w1, w2)
    alpha, beta, gamma = P.coeff_monomial(w1**2), P.coeff_monomial(w2**2), P.coeff_monomial(w1*w2)
    return alpha, beta, gamma


def static_part():
    """m1^2 m2^2 part of (1/4pi) int U^2 (grad U)^2: (grad U1)^2 U2^2 + (grad U2)^2 U1^2 + grad(U1^2).grad(U2^2)
    -> R(-4,-2) + R(-2,-4) - 2 R(-2,-4) (the last by parts: -int U1^2 lap(U2^2), lap(r^-2) = 2 r^-4)."""
    return riesz(-4, -2) + riesz(-2, -4) - 2*riesz(-2, -4)


def energy_shift_per_K():
    """Circular-orbit energy at fixed omega from L = mu v^2/2 + mu m/r + K (pi^2/4) m1^2 m2^2 v^2/r^3 (G = c = 1),
    to first order in K; returned as the coefficient c in E = -(mu x/2)(1 + ... + c K nu x^3)."""
    m1, m2, w, r, K, eps = sp.symbols('m1 m2 omega r K epsilon', positive=True)
    m = m1 + m2; mu = m1*m2/m; lam = eps*K*sp.pi**2/4*m1**2*m2**2
    # radial equation at rdot = 0 with v = r omega
    eq = mu*r*w**2 - mu*m/r**2 - lam*w**2/r**2
    r0 = (m/w**2)**sp.Rational(1, 3)
    r1 = sp.symbols('r1')
    sol = sp.solve(sp.series(eq.subs(r, r0 + eps*r1), eps, 0, 2).removeO().coeff(eps, 1), r1)[0]
    rr = r0 + eps*sol
    v2 = rr**2*w**2
    E = mu*v2/2 - mu*m/rr + lam*v2/rr**3          # v dL/dv - L for this L
    dE = sp.simplify(sp.series(E, eps, 0, 2).removeO().coeff(eps, 1))
    x = (m*w)**sp.Rational(2, 3); nu = m1*m2/m**2
    c = sp.simplify(dE/(-(mu*x/2)*K*nu*x**3))
    return sp.nsimplify(sp.simplify(c))


def circle_bracket(al, be, ga):
    """alpha v1^2 + beta v2^2 + gamma v1.v2 on a circle, in units of v12^2 (v1 = m2 v, v2 = -m1 v, m = 1)."""
    nu = sp.symbols('nu'); m1 = sp.symbols('m1')
    expr = sp.expand(al*(1 - m1)**2 + be*m1**2 - ga*m1*(1 - m1))
    return sp.simplify(expr.subs(m1**2, m1 - nu).subs(m1**2, m1 - nu))


def main():
    cen = census()
    pis = [k for k, (v, kind) in cen.items() if kind == 'pi^2']
    al, be, ga = moving_time_part()
    cK = energy_shift_per_K()
    target = -sp.Rational(205, 96)*sp.pi**2       # relativity's pi^2 coefficient of nu x^3 in E/(-mu x/2)
    K_GR = sp.simplify(target/cK)
    num, exact = numeric_overlap()
    out = {
        'census_pi2_entries': pis,
        'R(-2,-2)': str(riesz(-2, -2)),
        'numeric_overlap_r12_2': [num, exact],
        'identity_ok': identity_check(),
        'static_pi2_part': str(static_part()),
        'moving_time_part_alpha_beta_gamma': [str(al), str(be), str(ga)],
        'boost_safe_time_part': bool(2*al + ga == 0 and 2*be + ga == 0),
        'time_part_on_circles': str(circle_bracket(al, be, ga)),
        'energy_coefficient_per_K': str(cK),
        'K_GR': str(K_GR),
    }
    json.dump(out, open(os.path.join(HERE, 'fifth_rule.json'), 'w'), indent=1)
    print(json.dumps(out, indent=1))


if __name__ == '__main__':
    main()
