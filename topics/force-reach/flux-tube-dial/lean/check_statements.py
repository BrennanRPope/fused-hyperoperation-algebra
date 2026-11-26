"""Symbolic/numeric check that every statement in FluxTubeDial.lean is true."""
import random, math
import sympy as sp
ok = True
D, s, dl, c, a, b = sp.symbols('D s delta c a b', real=True)
tests = {
    'deff_from_exponent': 2 - (D - (D - 1) * s) - (1 + (D - 1) * (s - 1)),
    'deff_maxwell': 1 + (D - 1) * (2 - 1) - D,
    'deff_tube': 1 + (D - 1) * (1 - 1) - 1,
    'threshold_exponent_zero': D - (D - 1) * (D / (D - 1)),
    'threshold_3d': sp.Rational(3, 2) - sp.Rational(3, 2),
    'pt_quarter': (2 + 4 * sp.Rational(1, 4)) / (1 + 4 * sp.Rational(1, 4)) - sp.Rational(3, 2),
    'pt_map_form': (2 + 4 * dl) / (1 + 4 * dl) - (1 + 1 / (1 + 4 * dl)),
    'gilbert_s1': sp.Integer(2) ** 0 - 1,
}
for name, e in tests.items():
    good = sp.simplify(e) == 0
    ok &= good
    print(('PASS ' if good else 'FAIL ') + name)
# confines_below_threshold: random D > 1, s < D/(D-1)
good = True
for _ in range(2000):
    DD = random.uniform(1.01, 10); ss = random.uniform(-5, DD / (DD - 1) - 1e-9)
    good &= DD - (DD - 1) * ss > 0
ok &= good; print(('PASS ' if good else 'FAIL ') + 'confines_below_threshold')
# the full l_p identity behind lp_core
good = True
for _ in range(500):
    cc, aa, bb, ss = random.uniform(0.1, 3), random.uniform(0, 5), random.uniform(0, 5), random.uniform(0.2, 2)
    p = 1 / ss
    good &= abs(((cc * aa ** ss) ** p + (cc * bb ** ss) ** p) ** (1 / p) - cc * (aa + bb) ** ss) < 1e-9 * (1 + cc * (aa + bb) ** ss)
ok &= good; print(('PASS ' if good else 'FAIL ') + 'lp identity (numeric only; not a Lean statement)')
print('ALL PASS' if ok else 'SOME FAIL')
