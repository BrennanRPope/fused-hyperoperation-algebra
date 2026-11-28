"""DRAFT. Round 32: give field energy the light-speed lag. Does the pi^2 come out at relativity's size, K = 41/16?

Answer: no. With the lag, field energy's own channel makes no pi^2 at all, so 41/16 has to come from elsewhere.

1. Identity. For any field U, with dU.dU = (grad U)^2 - (dU/dt)^2 and box = -d^2/dt^2 + lap,
       U^2 dU.dU = (1/3) d_mu(U^3 d^mu U) - (1/3) U^3 box U.
   Integrated over space, the first term is an exact total time derivative of a function of the bodies' positions and
   velocities, which changes no orbit and no energy. So field energy acting on field energy is, up to that, -(1/3) int U^3 box U.
2. With the lag. The lagged potential U + dU, dU = m d^2 r/dt^2 / 2 (the first term of the retarded m/r), satisfies
   box(U + dU) = the point sources, to this order. Then -(1/3) int U^3 box U only samples fields at the bodies:
   powers of r12 with rational coefficients. A pi^2 needs the overlap R(-2,-2) = pi^2/4 of round 31, which a
   point sample never makes. So the lagged m1^2 m2^2 field energy has no pi^2: K_time = 0.
3. The bookkeeping, done explicitly with F = (pi^2/4)/r12 (the pi^2 part of (1/4pi) int U1^2 U2^2 / (m1^2 m2^2)),
   D_a = v_a . grad_a and A = a1 . grad1 + a2 . grad2 (the accelerations):
     time part alone, -int U^2 (dU/dt)^2        -> -(1/8)(D1^2 + D2^2) F - D1 D2 F         (round 31)
     the lag's change to the static part         -> -(3/8)(D1^2 + D2^2) F - (1/2) A F
     sum                                         -> -(1/2)[(D1 + D2)^2 + A] F = -(1/2) d^2F/dt^2
   The sum is exactly minus half the second time derivative of F: a total time derivative. Round 31's statement that
   the lag "restores the boost-safe shape" was wrong: the lag removes the time part's pi^2 entirely.
4. So what's left. In a single field with two-derivative field energy, f(U)(dU)^2, a relabeling of the field makes it
   free, and the identity above makes its lagged self-energy local. Relativity's 41/16 has to come from field energy
   between different parts of the field: the time part meeting the push (the spatial part, relativity's change to the
   path) and the twist (frame dragging). Round 33 (pi2_channels.py) narrows this further.

    python3 lagged_field_energy.py     # writes lagged_field_energy.json (seconds)
"""
import json, os
import sympy as sp
from fifth_rule import riesz, y1, y2, v1, v2, F, vd, moving_time_part

HERE = os.path.dirname(os.path.abspath(__file__))
a1 = sp.Matrix(sp.symbols('a1x a1y a1z', real=True)); a2 = sp.Matrix(sp.symbols('a2x a2y a2z', real=True))
D1 = lambda f: vd(v1, y1, f)
D2 = lambda f: vd(v2, y2, f)
A = lambda f: vd(a1, y1, f) + vd(a2, y2, f)


def identity_check():
    """U^2 dU.dU - [(1/3) d_mu(U^3 d^mu U) - (1/3) U^3 box U] == 0 for a generic U(t, x, y, z)."""
    t, x, y, z = sp.symbols('t x y z', real=True)
    U = sp.Function('U')(t, x, y, z)
    sig = {t: -1, x: 1, y: 1, z: 1}
    dd = sum(sig[c]*sp.diff(U, c)**2 for c in (t, x, y, z))
    div = sum(sig[c]*sp.diff(U**3*sp.diff(U, c), c) for c in (t, x, y, z))
    box = sum(sig[c]*sp.diff(U, c, 2) for c in (t, x, y, z))
    return sp.simplify(U**2*dd - (div/3 - U**3*box/3)) == 0


def time_derivative_of(G):
    """d/dt of G(y1, y2, v1, v2) along the motion (dy/dt = v, dv/dt = a)."""
    out = 0
    for i in range(3):
        out += v1[i]*sp.diff(G, y1[i]) + v2[i]*sp.diff(G, y2[i]) + a1[i]*sp.diff(G, v1[i]) + a2[i]*sp.diff(G, v2[i])
    return out


def pieces():
    """pi^2 parts (units of m1^2 m2^2, after the 1/4pi), as functions of the positions, velocities, accelerations."""
    # time part alone: -int U^2 (dU/dt)^2, the negative of round 31's moving part
    T_time = -(D2(D2(F))/8 + D1(D1(F))/8 + D1(D2(F)))
    # J = pi^2 of int (U1^2 U2 d^2U2/dt^2 + U1 U2^2 d^2U1/dt^2):
    #   U d^2U/dt^2 = (1/2) d^2(U^2)/dt^2 - (dU/dt)^2, d^2U_a/dt^2 = D_a^2 U_a + a_a . grad_a U_a,
    #   (dU/dt)^2 = [D^2 r^-2 + 2 v^2 r^-4]/8 and R(-2,-4) = 0, so U1^2 (dU2/dt)^2 -> D2^2 F / 8
    J = (D2(D2(F))/2 - D2(D2(F))/8 + vd(a2, y2, F)/2) + (D1(D1(F))/2 - D1(D1(F))/8 + vd(a1, y1, F)/2)
    # the lag's change to the static part int U^2 (grad U)^2 = -(1/3) int U^3 lap U:
    #   -int U^2 dU lap U (point samples, no pi^2) - (1/3) int U^3 lap dU, lap dU = d^2U/dt^2, m1^2 m2^2 part -> -J
    lag = -J
    return T_time, J, lag


def lap_lag_check():
    """lap_x of dU = (m/2) d^2|x - y|/dt^2 equals d^2(m/|x - y|)/dt^2 (velocity and acceleration terms); box(U + dU) = source."""
    x = sp.Matrix(sp.symbols('X Y Z', real=True)); r = sp.sqrt(((x - y2).T*(x - y2))[0])
    d2 = lambda f: vd(v2, y2, vd(v2, y2, f)) + vd(a2, y2, f)
    dU = d2(r)/2
    lap = sum(sp.diff(dU, x[i], 2) for i in range(3))
    return sp.simplify(lap - d2(1/r)) == 0


def operator_weights(expr):
    """Coefficients of D1^2 F, D2^2 F, D1 D2 F, A F in expr, read off on generic numbers."""
    basis = [D1(D1(F)), D2(D2(F)), D1(D2(F)), A(F)]
    c = sp.symbols('c0:4')
    import random
    rng = random.Random(32)
    syms = list(y1) + list(y2) + list(v1) + list(v2) + list(a1) + list(a2)
    pts = [{s: sp.Rational(rng.randint(-30, 30), rng.randint(5, 13)) for s in syms} for k in range(6)]
    M = sp.Matrix([[sp.N(b.subs(p)/sp.pi**2, 40) for b in basis] for p in pts])
    rhs = sp.Matrix([sp.N(expr.subs(p)/sp.pi**2, 40) for p in pts])
    sol = M[:4, :].LUsolve(rhs[:4, :])
    resid = max(abs(x) for x in (M*sol - rhs))           # the last two points check the fit
    assert resid < 1e-25, resid
    return [sp.nsimplify(round(float(x), 12), rational=True) for x in sol]


def main():
    T_time, J, lag = pieces()
    G = D1(F) + D2(F)                                   # dF/dt
    d2F = time_derivative_of(G)
    total = T_time + lag
    out = {
        'identity_ok': identity_check(),
        'lap_lag_ok': lap_lag_check(),
        'round31_time_part_alpha_beta_gamma': [str(s) for s in moving_time_part()],
        'unlagged_equals_J_minus_half_d2F': sp.simplify(T_time - (J - d2F/2)) == 0,
        'lagged_total_equals_minus_half_d2F': sp.simplify(total + d2F/2) == 0,
        'weights_D1D1_D2D2_D1D2_A': {
            'time_part': [str(w) for w in operator_weights(T_time)],
            'lag': [str(w) for w in operator_weights(lag)],
            'sum': [str(w) for w in operator_weights(total)],
            'minus_half_d2F_dt2': [str(w) for w in operator_weights(-d2F/2)],
        },
        'K_time_with_lag': 0,
        'K_GR_needed': '41/16',
        'R(-2,-2)': str(riesz(-2, -2)),
    }
    json.dump(out, open(os.path.join(HERE, 'lagged_field_energy.json'), 'w'), indent=1)
    print(json.dumps(out, indent=1))


if __name__ == '__main__':
    main()
