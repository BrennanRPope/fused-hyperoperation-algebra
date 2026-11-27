"""Vectorised Hurwitz zeta (Euler-Maclaurin) and the L-functions used in round 7."""
import numpy as np
from scipy.special import loggamma, bernoulli

B2K = [float(b) for b in bernoulli(40)[2::2]]   # B_2, B_4, ...
import math
FACT = [float(math.factorial(2*k)) for k in range(1, 21)]

def hurwitz(s, a, K=14):
    """zeta(s, a) for complex array s (|Im s| large ok), real a in (0,1]."""
    s = np.asarray(s, dtype=complex)
    tmax = np.max(np.abs(s.imag)) if s.size else 0
    N = int(tmax/(2*np.pi)*1.6) + 30
    n = np.arange(N) + a
    logn = np.log(n)
    out = np.zeros(s.shape, dtype=complex)
    flat = s.ravel(); res = np.zeros(flat.shape, dtype=complex)
    CH = max(1, int(4e6 // N))
    for i in range(0, flat.size, CH):
        ss = flat[i:i+CH]
        res[i:i+CH] = np.exp(-np.outer(ss, logn)).sum(axis=1)
    X = N + a
    lx = np.log(X)
    res += np.exp((1-flat)*lx)/(flat-1) + 0.5*np.exp(-flat*lx)
    poch = flat.copy()            # (s)_1
    term_pow = np.exp((-flat-1)*lx)
    for k in range(1, K+1):
        res += B2K[k-1]/FACT[k-1] * poch * term_pow
        poch = poch*(flat+2*k-1)*(flat+2*k)
        term_pow = term_pow/(X*X)
    return res.reshape(s.shape)

def dirichlet_series_periodic(s, coeffs):
    """sum a_n n^-s with a_n periodic of period q, coeffs = [a_1..a_q]."""
    q = len(coeffs)
    s = np.asarray(s, dtype=complex)
    tot = np.zeros(s.shape, dtype=complex)
    for r, c in enumerate(coeffs, start=1):
        if c != 0:
            tot += c*hurwitz(s, r/q)
    return tot*np.exp(-s*np.log(q))

def kron(D, n):
    from sympy.ntheory import factorint
    from sympy.functions.combinatorial.numbers import jacobi_symbol as js
    val = 1
    for p, e in factorint(n).items():
        if p == 2:
            if D % 2 == 0: return 0
            v = 1 if D % 8 in (1, 7) else -1
        else:
            v = int(js(D % p, p)) if D % p else 0
        val *= v**e
    return int(val)

def char_coeffs(D):
    q = abs(D)
    return [kron(D, n) for n in range(1, q+1)]

zeta = lambda s: hurwitz(s, 1.0)
