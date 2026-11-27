"""Jensen polynomials of the Riemann Xi function and the circle degree of their real Newton maps.

Xi(2t) = 4 * integral_0^inf Phi(u) cos(2 t u) du, where Xi(t) = xi(1/2 + i t), with the kernel
Phi(u) = 2 sum_n (2 pi^2 n^4 e^{9u} - 3 pi n^2 e^{5u}) exp(-pi n^2 e^{4u}) (Polya 1926).
With gamma(m) = m! * 4/(2m)! * integral_0^inf u^{2m} Phi(u) du, the series sum gamma(m) x^m / m! is
Xi(2 sqrt(-x)) (x = -t^2), so gamma(0) = xi(1/2). Rescaling x does not move roots off the real line. Polya: RH holds iff every Jensen polynomial
J^{d,n}(X) = sum_j C(d, j) gamma(n + j) X^j has only real roots (Griffin, Ono, Rolen and Zagier 2019
prove it for each d once n is large). Writes jensen.json.
"""
import json
from math import comb
import mpmath as mp
import numpy as np
from dialtop import newton_parts, rp1_degree_alg, distinct_real_roots

mp.mp.dps = 60


def Phi(u):
    s = mp.mpf(0)
    for n in range(1, 8):  # terms beyond n = 7 are below 1e-60 here
        s += (2 * mp.pi ** 2 * n ** 4 * mp.e ** (9 * u) - 3 * mp.pi * n ** 2 * mp.e ** (5 * u)) * mp.e ** (-mp.pi * n ** 2 * mp.e ** (4 * u))
    return 2 * s


def gammas(M):
    out = []
    for m in range(M + 1):
        mom = mp.quad(lambda u: u ** (2 * m) * Phi(u), [0, 0.25, 0.5, 1, 2])
        out.append(mp.factorial(m) * 4 / mp.factorial(2 * m) * mom)
    return out


def xi_half():
    s = mp.mpf(1) / 2
    return s * (s - 1) / 2 * mp.pi ** (-s / 2) * mp.gamma(s / 2) * mp.zeta(s)


if __name__ == '__main__':
    D, NMAX = 8, 30
    g = gammas(D + NMAX)
    rows = []
    for d in range(1, D + 1):
        for n in range(0, NMAX + 1):
            # J(X) = sum_j C(d,j) g(n+j) X^j; rescale X so coefficients are of similar size
            c = [comb(d, j) * g[n + j] for j in range(d + 1)]
            ratio = c[0] / c[d] if c[d] != 0 else 1
            scale = mp.power(abs(ratio), mp.mpf(1) / d)
            cs = [float(c[j] * scale ** j / c[0]) for j in range(d + 1)]
            coeffs = cs[::-1]                       # highest power first
            r = len(distinct_real_roots(coeffs, tol=1e-6))
            num, den, _ = newton_parts(coeffs)
            deg = rp1_degree_alg(num, den) if d >= 2 else 1
            rows.append({'d': d, 'n': n, 'real_roots': r, 'newton_degree': deg})
    json.dump({'gamma': [mp.nstr(x, 30) for x in g], 'xi_half': mp.nstr(xi_half(), 30), 'rows': rows},
              open('jensen.json', 'w'), indent=1)
    print('xi(1/2)', mp.nstr(xi_half(), 20), 'gamma(0)', mp.nstr(g[0], 20))
    bad = [r for r in rows if r['real_roots'] != r['d'] or (r['d'] >= 2 and r['newton_degree'] != 2 - r['d'])]
    print(len(rows), 'Jensen polynomials;', len(bad), 'not hyperbolic or degree != 2 - d')
