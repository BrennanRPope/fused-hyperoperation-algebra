"""Data for index.html: flow polynomials P_n(t) (n<=10, float coefficients) and the
exact sign grid of P_n(t) for n = 1..16, t = 0..3 in steps of 1/40. Prints JSON."""
import json, math
import sympy as sp
from fractions import Fraction as F
t = sp.symbols('t'); Nb = 16; fact = math.factorial
Eb = [F(0)] + [F(1, fact(i)) for i in range(1, Nb + 1)]
def mulb(a, b):
    c = [F(0)] * (Nb + 1)
    for i, ai in enumerate(a):
        if ai:
            for j, bj in enumerate(b):
                if i + j <= Nb: c[i + j] += ai * bj
    return c
Mb = sp.zeros(Nb, Nb); pw = [F(1)] + [F(0)] * Nb
for j in range(1, Nb + 1):
    pw = mulb(pw, Eb)
    for i in range(1, Nb + 1): Mb[i - 1, j - 1] = sp.Rational(pw[i].numerator, pw[i].denominator)
Nm = Mb - sp.eye(Nb); L = sp.zeros(Nb, Nb); P = sp.eye(Nb)
for m in range(1, Nb + 1): P = P * Nm; L += ((-1)**(m + 1)) * P / m
v = sp.zeros(Nb, 1); v[0] = 1; cur = v; res = v.copy()
for m in range(1, Nb + 1): cur = L * cur / m; res += cur * t**m
Pb = [sp.Poly(sp.expand(res[i] * fact(i + 1)), t) for i in range(Nb)]
coef = [[float(c) for c in p.all_coeffs()[::-1]] for p in Pb[:10]]
grid = []
for p in Pb:
    row = ''
    for j in range(0, 121):
        val = p.eval(sp.Rational(j, 40))
        row += '-' if val < 0 else ('0' if val == 0 else '+')
    grid.append(row)
json.dump({'coef': coef, 'grid': grid}, open('figures.json', 'w'))
print(json.dumps({'coef': coef})[:300]); print('\n'.join(grid))
