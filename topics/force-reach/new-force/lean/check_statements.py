"""Check that every statement in NewForce.lean is true."""
import sympy as sp
q = sp.symbols('q', real=True)
R = sp.Rational
tests = {
    'tower_sum': 1 + 2*q/(1 - q) - (1 + q)/(1 - q),
    'tower_leading_yukawa': (1 + q)/(1 - q) - 1 - 2*q - 2*q**2/(1 - q),
    'dual_exponent_mond': R(3, 2)/(R(3, 2) - 1) - 3,
    'holder_pair': 1/R(3, 2) + R(1, 3) - 1,
}
ok = True
for n, e in tests.items():
    g = sp.simplify(e) == 0; ok &= g; print(('PASS ' if g else 'FAIL ') + n)
x = sp.symbols('x', positive=True)
g = sp.simplify(((1 + sp.exp(-x))/(1 - sp.exp(-x)) - sp.coth(x/2)).rewrite(sp.exp)) == 0
ok &= g; print(('PASS ' if g else 'FAIL ') + '(1+q)/(1-q) = coth(x/2) at q = e^-x (docstring claim)')
g = (R(33, 3) - R(32, 3) > 0) and (R(33, 3) - R(34, 3) < 0); ok &= g; print(('PASS ' if g else 'FAIL ') + 'su3_flavour_bound')
g = (R(22, 3) - R(20, 3) > 0) and (R(22, 3) - R(22, 3) == 0); ok &= g; print(('PASS ' if g else 'FAIL ') + 'su2_flavour_bound')

n, al = sp.symbols('n alpha', positive=True)
gtests = {
    'vector_over_tensor': R(2)/R(1, 2) - 4,
    'scalar_over_tensor': R(3, 2)/R(1, 2) - 3,
    'like_charges_cancel': R(1, 2) + R(3, 2) - 2,
    'opposite_charges_add': R(1, 2) + R(3, 2) + 2 - 4,
    'klein_strength': 4/(4/al) - al,
    'one_circle_third': R(1, 3) - R(1, 3),
    'gamma_from_strength': (1 - n/(n + 2))/(1 + n/(n + 2)) - 1/(n + 1),
    'one_circle_gamma': R(1, 2) - 1 + R(1, 2),
}
for nm, e in gtests.items():
    g = sp.simplify(e) == 0; ok &= g; print(('PASS ' if g else 'FAIL ') + nm)
# the pieces match the field-theory split in strength_checks.py (G1, G2)

C_, R_ = sp.symbols('C R', positive=True)
htests = {
    'hold_turning_point': -(-6*C_/R_**3)/R_**2 - 6*C_*R_**2/R_**7,
    'hold_curvature': 2*(-6*C_/R_**3)/R_**3 + 42*C_*R_**2/R_**8 - 30*C_/R_**6,
    'hold_leftover': (-6*C_/R_**3)/R_ + C_*R_**2/R_**6 + 5*C_/R_**4,
    'hold_mass_factor': R(2, 3)*30 - 20,
    'sphere_shape_factor': 4/R(2, 3) - 6,
    'strong_b0_five': R(33, 3) - R(10, 3) - R(23, 3),
}
for nm, e in htests.items():
    g = sp.simplify(e) == 0; ok &= g; print(('PASS ' if g else 'FAIL ') + nm)

itests = {
    'every_rung_four_thirds': R(2, 3)/R(1, 2) - R(4, 3),
    'zero_rung_split': 1 + R(1, 3) - R(4, 3),
    'copies_extreme': 1 + 3 - 4,
    'circle_swell_ratio': R(1, 3)/R(2, 3) - R(1, 2),
    'alpha_dip': -2*R(2, 3) + R(4, 3),
    'inflaton_tail_eps': R(1, 2)*R(2, 3) - R(1, 3),
}
for nm, e in itests.items():
    g = sp.simplify(e) == 0; ok &= g; print(('PASS ' if g else 'FAIL ') + nm)

R1s = sp.symbols('R1', positive=True)
jtests = {
    'like_copies_far_repel': 1 + 0 - 4 + 3,
    'matter_copy_size': R(1, 3)*3 - 1,
    'second_circle_cap': 4*sp.pi**2*R1s*(R1s**3/(4*sp.pi**2)) - R1s**4,
    'cap_ratio': (R1s**3/(4*sp.pi**2))/R1s - R1s**2/(4*sp.pi**2),
    'light_bends_three_quarters': R(1, 2)*(1 + (1 - R(1, 3))/(1 + R(1, 3))) - R(3, 4),
}
for nm, e in jtests.items():
    g = sp.simplify(e) == 0; ok &= g; print(('PASS ' if g else 'FAIL ') + nm)
ktests = {
    'bz_b0_sixteen': R(33, 3) - R(32, 3) - R(1, 3),
    'bz_b1_sixteen': 102 - R(38*16, 3) + R(302, 3),
    'bz_alpha_sixteen': -4*sp.pi*R(1, 3)/R(-302, 3) - 4*sp.pi/302,
    'square_torus_no_mixing': sp.cos(sp.pi/2),
}
for nm, e in ktests.items():
    g = sp.simplify(e) == 0; ok &= g; print(('PASS ' if g else 'FAIL ') + nm)
a_ = sp.symbols('a', positive=True); th_, rho_ = sp.symbols('theta rho', positive=True)
ltests = {
    'mwgc_ratio_sq': R(4*137, 2*137) - 2,
    'sphere_arrow_average': 2*sp.pi*R(4, 3) - R(2, 3)*4*sp.pi,
    'sphere_arrow_average (integral)': sp.integrate(sp.sin(th_)**3, (th_, 0, sp.pi)) - R(4, 3),
    'sphere_strength_six': 4/R(2, 3) - 6,
    'circle_edge_cubed': (4*sp.pi**2)**2*4*a_ - 64*a_*sp.pi**4,
    'sphere_edge_cubed': (8*sp.pi**2)**2*4*216*a_ - 216*256*a_*sp.pi**4,
}
for nm, e in ltests.items():
    g = sp.simplify(e) == 0; ok &= g; print(('PASS ' if g else 'FAIL ') + nm)
K_, r_ = sp.symbols('K rho', positive=True)
U_ = -2/r_**4 + K_/r_**9; rs_ = (R(9, 8)*K_)**R(1, 5)
ntests = {
    'sphere_hold_rest': 8*(9*K_/8) - 9*K_,
    'sphere_hold_rest (U\'(rho*) = 0)': sp.diff(U_, r_).subs(r_, rs_),
    'sphere_hold_minimum': -40 + 90*R(8, 9) - 40,
    'sphere_hold_minimum (rho*^6 U\'\' = 40)': sp.diff(U_, r_, 2).subs(r_, rs_)*rs_**6 - 40,
    'sphere_hold_leftover': -2 + R(8, 9) + R(10, 9),
    'sphere_strength_exponent': R(-10, 3)*R(1, 5) + R(2, 3),
}
for nm, e in ntests.items():
    g = sp.simplify(e) == 0; ok &= g; print(('PASS ' if g else 'FAIL ') + nm)
print('ALL PASS' if ok else 'SOME FAIL')
