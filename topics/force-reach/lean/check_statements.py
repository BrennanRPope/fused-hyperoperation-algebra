"""Numeric/symbolic check that every statement in ForceReach.lean is true."""
import random, math
import sympy as sp
ok = True
k, d, m, r, mu, K, x, r1, r2 = sp.symbols('k d m r mu K x r1 r2', real=True)
tests = {
    'disc_is_spread_sq': (r1 - r2) ** 2 - ((r1 + r2) ** 2 - 4 * r1 * r2),
    'drift_removal': (k / 2 + mu) ** 2 - k * (k / 2 + mu) + d - (mu ** 2 - (k ** 2 - 4 * d) / 4),
    'raychaudhuri_is_law': -K * (1 + K * x ** 2) - (-K - (-K * x) ** 2),
    'curvature_reading_k0': -(0 - 4 * K) / 4 - K,
    'yukawa_deff': 1 - r * (-m - 2 / r + m / (1 + m * r)) - (3 + m ** 2 * r ** 2 / (1 + m * r)),
}
for name, e in tests.items():
    good = sp.simplify(e) == 0
    ok &= good
    print(('PASS ' if good else 'FAIL ') + name)
# double_root_iff_border: random rationals
good = True
for _ in range(200):
    kk = random.uniform(-5, 5); dd = kk * kk / 4 if random.random() < 0.5 else random.uniform(-5, 5)
    lhs = abs((kk / 2) ** 2 - kk * kk / 2 + dd) < 1e-12
    rhs = abs(kk * kk - 4 * dd) < 1e-12
    good &= (lhs == rhs)
ok &= good; print(('PASS ' if good else 'FAIL ') + 'double_root_iff_border')
D = -3.7
good = abs(2 * math.pi / (2 * math.pi / math.sqrt(-D)) - math.sqrt(-D)) < 1e-14
ok &= good; print(('PASS ' if good else 'FAIL ') + 'rung_mass')
s = sp.symbols('s', real=True)
good = sp.simplify(sp.diff(-m * sp.exp(-m * s), s) - m ** 2 * sp.exp(-m * s)) == 0
ok &= good; print(('PASS ' if good else 'FAIL ') + 'yukawa_reduced')
print('ALL PASS' if ok else 'SOME FAIL')
