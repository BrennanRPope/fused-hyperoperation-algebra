"""DRAFT. Round 33: where the push and the twist meet the time part, and what relativity's pi^2 actually is.

The river's pieces, in the split of Einstein's field that matches them (Kaluza-Klein, as in Kol and Smolkin):
  time part  -> phi    (g00 = -e^{2 phi}: Newton's pull, the clocks)
  push       -> sigma  (the spatial metric gamma_ij: relativity's change to the path)
  twist      -> A_i    (g0i: frame dragging)

1. What relativity's pi^2 is, pinned by two invariants. Write the pi^2 part of the two-body Lagrangian (G = c = 1,
   per unit reduced mass, m = 1) as  nu pi^2 [ a v^2/r^3 + b rdot^2/r^3 + s/r^4 ]. Coordinates move it around, but two
   things don't: the energy of circular orbits at a given frequency, and the periastron advance of nearly circular ones.
   Relativity's are -(205/96) pi^2 nu x^3 in E/(-mu x/2) (round 30) and +(123/32) pi^2 nu x^3 in the advance k
   (Damour, Jaranowski, Schafer 2000). ADM's pi^2 terms (DJS: p^2/64, -3 (n.p)^2/64, -21/(32 r^4)) give both, and so do
   two simpler forms:
       harmonic-like:  (41/64) pi^2 G^3 m1^2 m2^2 (v12^2 - 3 rdot^2)/r^3 = -(41/16) (D1 + D2)^2 F
       static:         (41/64) pi^2 G^4 m1^2 m2^2 m / r^4             = +(41/16) (a . grad) F
   with F = (pi^2/4) G^3 m1^2 m2^2 / r12 (round 31's overlap of the two bodies' fields squared), D_a = v_a . grad_a and
   a the accelerations. The two differ by (41/16) d^2F/dt^2, a total time derivative. So relativity's pi^2 is one
   thing: F felt through the bodies' acceleration, with weight 41/16. Round 31's boost-safe K v12^2 (K = 41/16) gets
   circles right but not the periastron (369/64 instead of 123/32): it is not the same physics.
2. The time part meeting the push, with the bodies at rest, makes no pi^2. Static, phi enters Einstein's action only
   as gamma^{ij} d_i phi d_j phi, so every static meeting has an even number of time-part legs. A third-level static
   term has five masses (odd), so it always contains a direct body-to-body Newtonian link: it is a second-level piece
   (rational, round 14) times m/r12. No pi^2.
3. The time part meeting itself while moving makes no pi^2 once it lags (round 32). So relativity's pi^2 has to come
   from the moving meetings: the twist (frame dragging) with the time part, and the push with the time part while
   the bodies move. Those are relativity's G^3 m1^2 m2^2 v^2 terms; this script does not compute them.

    python3 pi2_channels.py      # writes pi2_channels.json (about a minute)
"""
import json, os
import sympy as sp

HERE = os.path.dirname(os.path.abspath(__file__))
w, eps, J, x = sp.symbols('omega epsilon J x', positive=True)
r, r0, R = sp.symbols('r r0 R', positive=True)

FORMS = {   # (a, b, s) in nu pi^2 [a v^2/r^3 + b rdot^2/r^3 + s/r^4], Lagrangian sign
    'ADM (DJS, from -H)': (-sp.Rational(1, 64), sp.Rational(3, 64), sp.Rational(21, 32)),
    'harmonic-like': (sp.Rational(41, 64), -sp.Rational(123, 64), 0),
    'static only': (0, 0, sp.Rational(41, 64)),
    'round 31 boost-safe K': (sp.Rational(41, 64), 0, 0),
}
GR_E, GR_K = -sp.Rational(205, 96), sp.Rational(123, 32)


def circle_energy(a, b, s):
    """pi^2 nu x^3 coefficient of E/(-mu x/2) at fixed frequency, first order in the pi^2 terms."""
    Lc = R**2*w**2/2 + 1/R + eps*(a*R**2*w**2/R**3 + s/R**4)
    rr0 = w**sp.Rational(-2, 3); r1 = sp.symbols('r1')
    sol = sp.solve(sp.series(sp.diff(Lc, R).subs(R, rr0 + eps*r1), eps, 0, 2).removeO().coeff(eps, 1), r1)[0]
    E = w*sp.diff(Lc, w) - Lc
    dE = sp.series(E.subs(R, rr0 + eps*sol), eps, 0, 2).removeO().coeff(eps, 1)
    X = w**sp.Rational(2, 3)
    return sp.nsimplify(sp.simplify(dE/(-X/2*X**3)))


def periastron(a, b, s):
    """pi^2 nu x^3 coefficient of k = omega_phi/omega_r - 1 for nearly circular orbits (Newtonian k = 0 exactly)."""
    mr = 1 + 2*eps*(a + b)/r**3                       # mass for rdot
    mt = 1 + 2*eps*a/r**3                             # mass for r dphi/dt
    V = J**2/(2*r**2*mt) - 1/r - eps*s/r**4           # effective potential at fixed angular momentum J
    Js = [j for j in sp.solve(sp.series(sp.diff(V, r).subs(r, r0), eps, 0, 2).removeO(), J)
          if j.subs({eps: 0, r0: 4}).is_positive][0]
    Js = sp.series(Js, eps, 0, 2).removeO()
    wr2 = (sp.diff(V, r, 2)/mr).subs(r, r0).subs(J, Js)
    wphi = (J/(r**2*mt)).subs(r, r0).subs(J, Js)
    k1 = sp.series(wphi/sp.sqrt(wr2), eps, 0, 2).removeO().coeff(eps, 1)
    return sp.nsimplify(sp.simplify(k1.subs(r0, 1/x)/x**3))


def static_parity():
    """Static bulk vertices of Einstein's action in the time-part/push split: phi only through
    sqrt(gamma) gamma^{ij} d_i phi d_j phi, so each vertex has 0 or 2 phi legs; sigma (push) has no static coupling
    to a point mass. A connected bulk piece's outside phi legs = sum of even numbers minus 2 per inner phi line: even.
    Returns, for 2..5 masses, whether a static diagram can avoid a direct body-to-body link (needs an even count)."""
    phi_degrees = {'sigma phi phi': 2, 'sigma sigma phi phi': 2, 'sigma^n': 0}
    assert all(d % 2 == 0 for d in phi_degrees.values())
    return {n: (n % 2 == 0) for n in range(2, 6)}


def gr_weights():
    """relativity's pi^2 in round 32's basis (D1^2 F, D2^2 F, D1 D2 F, (a.grad) F), both forms."""
    import lagged_field_energy as lg
    F, D1, D2, A = lg.F, lg.D1, lg.D2, lg.A
    harm = -sp.Rational(41, 16)*(D1(D1(F)) + D2(D2(F)) + 2*D1(D2(F)))
    stat = sp.Rational(41, 16)*A(F)
    d2F = lg.time_derivative_of(D1(F) + D2(F))
    same_up_to_td = sp.simplify(harm - stat + sp.Rational(41, 16)*d2F) == 0
    # check the harmonic-like form is (41/64) pi^2 (v12^2 - 3 rdot^2)/r^3 with v1 = v12, v2 = 0 at a point
    y1, y2, v1, v2 = lg.y1, lg.y2, lg.v1, lg.v2
    pt = {y1[0]: 3, y1[1]: 1, y1[2]: 2, y2[0]: 0, y2[1]: 0, y2[2]: 0, v1[0]: 1, v1[1]: -2, v1[2]: 5,
          v2[0]: 0, v2[1]: 0, v2[2]: 0}
    rr = sp.sqrt(14); rdot = sp.Rational(3 - 2 + 10)/rr; v2_ = 30
    harm_ok = sp.simplify(harm.subs(pt) - sp.Rational(41, 64)*sp.pi**2*(v2_ - 3*rdot**2)/rr**3) == 0
    return {'harmonic_like': [str(w_) for w_ in lg.operator_weights(harm)],
            'static': [str(w_) for w_ in lg.operator_weights(stat)],
            'differ_by_41_16_d2F_dt2': same_up_to_td, 'harmonic_form_matches_v2_minus_3rdot2': harm_ok}


def main():
    inv = {k: {'circle_energy': str(circle_energy(*v)), 'periastron': str(periastron(*v))} for k, v in FORMS.items()}
    out = {
        'relativity': {'circle_energy': str(GR_E), 'periastron': str(GR_K)},
        'forms': {k: [str(c) for c in v] for k, v in FORMS.items()},
        'invariants': inv,
        'gr_forms_all_match': all(inv[k]['circle_energy'] == str(GR_E) and inv[k]['periastron'] == str(GR_K)
                                  for k in ('ADM (DJS, from -H)', 'harmonic-like', 'static only')),
        'round31_misses_periastron': inv['round 31 boost-safe K']['periastron'] != str(GR_K),
        'gr_weights_D1D1_D2D2_D1D2_A': gr_weights(),
        'static_can_avoid_direct_link_by_mass_count': static_parity(),
        'time_part_alone_lagged_K': 0,
    }
    json.dump(out, open(os.path.join(HERE, 'pi2_channels.json'), 'w'), indent=1)
    print(json.dumps(out, indent=1))


if __name__ == '__main__':
    main()
