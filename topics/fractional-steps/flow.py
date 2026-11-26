"""Fractional steps of the critical Poisson(1) branching process.

Offspring pgf f(s) = e^(s-1); in x = s - 1 this is the rung step g(x) = e^x - 1.
Abel function on the unit disk:  alpha(s) = lim [phi(g^n(s-1)) - n],
generator (only possible one, see README):  u(s) = 1/alpha'(s),
fractional step:  F_t(s) = alpha^{-1}(alpha(s) + t)  (push n steps, shift, pull back n steps).
"""
import mpmath as mp
from fatou_series import c as C_FR, a as A_FR

mp.mp.dps = 40
C = mp.mpf(C_FR.numerator) / C_FR.denominator
A = [mp.mpf(v.numerator) / v.denominator for v in A_FR[:24]]
NSTEP = 1500


def phi(x):
    return -2 / x + C * mp.log(-x) + sum(A[k] * x ** (k + 1) for k in range(len(A)))


def dphi(x):
    return 2 / x ** 2 + C / x + sum((k + 1) * A[k] * x ** k for k in range(len(A)))


def push(s, n=NSTEP):
    """x_n = g^n(s-1) and sum_{j<n} x_j (log of the chain-rule product)."""
    x = mp.mpc(s) - 1
    S = mp.mpf(0)
    for _ in range(n):
        S += x
        x = mp.expm1(x)
    return x, S


def alpha(s):
    x, _ = push(s)
    return phi(x) - NSTEP


def dalpha(s):
    x, S = push(s)
    return dphi(x) * mp.exp(S)


def u(s):
    return 1 / dalpha(s)


def inv_phi(z, guess):
    y = guess
    for _ in range(60):
        dy = (phi(y) - z) / dphi(y)
        y -= dy
        if abs(dy) < mp.mpf(10) ** (-mp.mp.dps + 5):
            break
    return y


def F(t, s, n=NSTEP):
    x, _ = push(s, n)
    y = inv_phi(phi(x) + t, x)
    for _ in range(n):
        y = mp.log1p(y)
    return 1 + y


def taylor(fun, M=128, rho=0.9):
    """Taylor coefficients at s = 0 from M samples on |s| = rho."""
    pts = [rho * mp.expjpi(2 * mp.mpf(j) / M) for j in range(M)]
    vals = [fun(p) for p in pts]
    coeffs = []
    for k in range(M // 2):
        ck = mp.fsum(vals[j] * mp.expjpi(-2 * mp.mpf(j * k) / M) for j in range(M)) / M
        coeffs.append(ck / rho ** k)
    return coeffs
