"""DRAFT. Round 30: where must the river's two-body rules stop matching general relativity?
Written and run BEFORE any measured gravitational-wave bound was looked up in this round.

1. The structural claim. Every piece of the rules' two-body field (Parts 9-10: the light-speed
   potentials, multiplied clocks, dressed masses, and the field-energy pieces built from Fock's
   ln(r1+r2+r12) and four Hadamard integrals) is a closed form whose finite part at a body is a
   rational number at a rational mass ratio. So the rules give rational coefficients at every level.
   General relativity's third-level (3PN) circular-orbit energy has the term -(205 pi^2/96) nu x^4
   (Blanchet, Living Reviews in Relativity 17, 2 (2014), Eq. 233), and through the first law of
   binary mechanics (Le Tiec, Blanchet and Whiting 2012) the redshift z1 inherits a pi^2 term at
   x^4. A rational number never equals rational + c pi^2 with c != 0, so the rules, as they stand,
   miss at the third level at every mass ratio with nu > 0, and only the test-mass limit (nu = 0,
   one exact Schwarzschild shell) is safe. The pi^2 comes from gravity's field energy acting on
   itself (cubic and higher field-on-field terms), which none of the four rules contains.

2. A size floor. The pi^2 piece is the part the rules are sure to lack. If their rational part
   happened to agree with relativity, the miss in the energy would be exactly
   delta e3 = +205 pi^2 nu / 96. Keeping the emitted power as in relativity (the rules' power is
   untested), the stationary-phase rule phi_k = 40 c_k / ((k-8)(k-5)) turns that into a shift of
   the 3PN phase coefficient, Delta phi_6 = -80 delta e3, and the fractional shift that
   gravitational-wave tests report, d_phi6 = Delta phi_6 / phi_6(non-log part).

    python3 third_level_predict.py      # writes third_level_prediction.json (a few seconds)
"""
import hashlib, json, os, datetime
import sympy as sp

HERE = os.path.dirname(os.path.abspath(__file__))
nu, x, v = sp.symbols('nu x v', positive=True)
pi2 = sp.pi**2

# known energy of a circular binary, E = -(mu x / 2)(1 + e1 x + e2 x^2 + e3 x^3)   (Blanchet 2014, Eq. 233)
e1 = -sp.Rational(3, 4) - nu/12
e2 = -sp.Rational(27, 8) + sp.Rational(19, 8)*nu - nu**2/24
e3 = (-sp.Rational(675, 64) + (sp.Rational(34445, 576) - sp.Rational(205, 96)*pi2)*nu
      - sp.Rational(155, 96)*nu**2 - sp.Rational(35, 5184)*nu**3)
# emitted power through 2PN (Blanchet 2014, Eq. 314), used only to check the phase rule on phi_2..phi_4
f2 = -sp.Rational(1247, 336) - sp.Rational(35, 12)*nu
f3 = 4*sp.pi
f4 = -sp.Rational(44711, 9072) + sp.Rational(9271, 504)*nu + sp.Rational(65, 18)*nu**2
# the 3PN phase coefficient of the standard stationary-phase waveform, non-log part
# (Buonanno, Iyer, Ochsner, Pan and Sathyaprakash, PRD 80, 084043 (2009), Eq. 3.3)
phi6_nonlog = (sp.Rational(11583231236531, 4694215680) - sp.Rational(640, 3)*pi2
               - sp.Rational(6848, 21)*sp.EulerGamma
               + (-sp.Rational(15737765635, 3048192) + sp.Rational(2255, 12)*pi2)*nu
               + sp.Rational(76055, 1728)*nu**2 - sp.Rational(127825, 1296)*nu**3)


def redshift_x4():
    """z1 through x^4 from the first law, z1 = 1 + dE/dm1 - omega dJ/dm1 at fixed omega, with J from dE = omega dJ."""
    m1, m2, w, X = sp.symbols('m1 m2 omega X', positive=True)
    m = m1 + m2; n = m1*m2/m**2; mu = m1*m2/m
    xx = (m*w)**sp.Rational(2, 3)
    sub = lambda c: c.subs(nu, n)
    Ex = -(mu*X/2)*(1 + sub(e1)*X + sub(e2)*X**2 + sub(e3)*X**3)
    # circular orbits: dE/dX = omega dJ/dX with omega = X^(3/2)/m  ->  J = m * int X^(-3/2) dE/dX dX
    J = m*sp.integrate(sp.expand(sp.diff(Ex, X)*X**sp.Rational(-3, 2)), X)
    E = Ex.subs(X, xx); J = J.subs(X, xx)
    z1 = 1 + sp.diff(E, m1) - w*sp.diff(J, m1)
    z1X = z1.subs(w, X**sp.Rational(3, 2)/m).subs(m2, 1 - m1)
    s = sp.series(sp.simplify(z1X), X, 0, 5).removeO()
    return [sp.expand(sp.simplify(s.coeff(X, k))) for k in range(5)], m1


def phase_rule_check():
    """phi_k = 40 c_k/((k-8)(k-5)), c_k the v^k coefficient of E'(v)/F(v) normalised to 1; checked on phi_2..phi_4."""
    E = v**2*(1 + e1*v**2 + e2*v**4 + e3*v**6)
    F = 1 + f2*v**2 + f3*v**3 + f4*v**4
    r = sp.series(sp.diff(E, v)/(2*v)/F, v, 0, 5).removeO()
    got = {k: sp.simplify(40*r.coeff(v, k)/((k - 8)*(k - 5))) for k in (2, 3, 4)}
    known = {2: sp.Rational(3715, 756) + sp.Rational(55, 9)*nu, 3: -16*sp.pi,
             4: sp.Rational(15293365, 508032) + sp.Rational(27145, 504)*nu + sp.Rational(3085, 72)*nu**2}
    return {k: bool(sp.simplify(got[k] - known[k]) == 0) for k in got}


def floor_shift(n):
    """Fractional shift of phi_6 if the rules miss only the pi^2 piece of e3 (power kept as in relativity)."""
    de3 = sp.Rational(205, 96)*pi2*n
    # E'(v) = 2v(1 + 2 e1 v^2 + 3 e2 v^4 + 4 e3 v^6): a change de3 moves c_6 by 4 de3, and phi_6 = 40 c_6/((-2)(1))
    dphi6 = -20*4*de3
    return float(de3), float(dphi6), float(dphi6/phi6_nonlog.subs(nu, n))


def main():
    cs, m1 = redshift_x4()
    pi_part_z1 = sp.factor(cs[4].coeff(sp.pi**2))
    rational_check = [sp.factor(cs[k] - next_order_truth(m1)[k]) for k in range(4)]
    grid = [0, 0.01, 0.02, 0.05, 0.08, 0.1, 0.15, 0.2, 0.2222, 0.24, 0.25]
    rows = []
    for n in grid:
        de3, dphi6, frac = floor_shift(sp.nsimplify(n))
        q = float(((1 - 2*n) - (1 - 4*n)**0.5)/(2*n)) if n > 0 else 0.0
        rows.append({'nu': n, 'mass_ratio_q': round(q, 4), 'delta_e3': de3, 'delta_phi6': dphi6,
                     'fractional_d_phi6': frac, 'phi6_nonlog': float(phi6_nonlog.subs(nu, sp.nsimplify(n)))})
    out = {
        'written_utc': datetime.datetime.now(datetime.timezone.utc).strftime('%Y-%m-%dT%H:%M:%SZ'),
        'uses_any_measured_bound': False,
        'claim': ('The rules give rational third-level coefficients; relativity has a pi^2 term proportional to nu; '
                  'so the rules miss at the third level for every nu > 0 and match only at nu = 0.'),
        'z1_x4_pi2_coefficient_body1': str(pi_part_z1),
        'z1_lower_orders_match_round13_truth': all(c == 0 for c in rational_check),
        'phase_rule_checks_phi2_phi4': phase_rule_check(),
        'floor_prediction': rows,
        'prediction_summary': ('If the rules miss only the pi^2 piece, the 3PN phase coefficient shifts by '
                               f"{rows[-1]['fractional_d_phi6']:.3f} of itself at equal masses and by "
                               f"{rows[4]['fractional_d_phi6']:.3f} at nu = 0.08 (about 1:9); zero at nu = 0."),
    }
    out['script_sha256'] = hashlib.sha256(open(__file__, 'rb').read()).hexdigest()
    json.dump(out, open(os.path.join(HERE, 'third_level_prediction.json'), 'w'), indent=1)
    print(json.dumps(out, indent=1))


def next_order_truth(m1):
    m2 = 1 - m1
    return [sp.Integer(1), -m2*(2*m1 + 3*m2)/2,
            -m2*(36*m1**3 + 100*m1**2*m2 + 96*m1*m2**2 + 27*m2**3)/24,
            -m2*(54*m1**5 + 224*m1**4*m2 + 356*m1**3*m2**2 + 284*m1**2*m2**3 + 124*m1*m2**4 + 27*m2**5)/16]


if __name__ == '__main__':
    main()
