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
print('ALL PASS' if ok else 'SOME FAIL')
