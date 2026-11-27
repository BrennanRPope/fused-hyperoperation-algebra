"""Torus, Klein-bottle and odd-half spectral zetas (square torus R^2/Z^2, glide g(x,y)=(x+1/2,-y)).

Dropping the common factor (2 pi)^(-2s) (round 8, klein.py):
  torus      Z_T(s) = 4 zeta(s) L(s, chi_-4)        = sum r2(n) n^-s   (Epstein zeta of the square lattice)
  Klein      Z_K(s) = Z_T(s)/2 - eta(2s)            eta(2s) = (1 - 2^(1-2s)) zeta(2s) = sum (-1)^(k+1) k^-2s
  odd half   Z_O(s) = Z_T(s)/2 + eta(2s)            (tones that flip sign round the twist)
  Z_K + Z_O = Z_T.
Completed torus zeta: pi^-s Gamma(s) Z_T(s) is symmetric under s -> 1-s, so
  theta_T(t) = Im log Gamma(1/2 + i t) - t log pi,   Z_T(1/2+it) e^{i theta_T(t)} is real.
Riemann: theta_R(t) = Im log Gamma(1/4 + i t/2) - (t/2) log pi.
"""
import numpy as np
from scipy.special import loggamma
from fastl import hurwitz, dirichlet_series_periodic

def zeta(s):
    return hurwitz(s, 1.0)

def L4(s):
    return dirichlet_series_periodic(s, [1, 0, -1, 0])

def eta2(s):
    s = np.asarray(s, dtype=complex)
    return (1 - np.exp((1 - 2*s)*np.log(2)))*hurwitz(2*s, 1.0)

def parts(s):
    """(Z_T(s), eta(2s)) on an array."""
    s = np.asarray(s, dtype=complex)
    return 4*zeta(s)*L4(s), eta2(s)

def Z(s, kind):
    zt, e = parts(s)
    if kind == 'torus': return zt
    if kind == 'klein': return zt/2 - e
    if kind == 'odd':   return zt/2 + e
    raise ValueError(kind)

def theta_T(t):
    t = np.asarray(t, dtype=float)
    return loggamma(0.5 + 1j*t).imag - t*np.log(np.pi)

def theta_R(t):
    t = np.asarray(t, dtype=float)
    return loggamma(0.25 + 0.5j*t).imag - 0.5*t*np.log(np.pi)

def gram_points(theta, t0, t1, n_from=None):
    """All Gram points g with theta(g) = n pi in [t0, t1] (theta increasing there), by bisection."""
    lo, hi = theta(np.array([t0, t1]))
    n0 = int(np.ceil(lo/np.pi)); n1 = int(np.floor(hi/np.pi))
    ns = np.arange(n0, n1 + 1)
    a = np.full(ns.shape, float(t0)); b = np.full(ns.shape, float(t1))
    for _ in range(60):
        m = 0.5*(a + b)
        up = theta(m) > ns*np.pi
        b = np.where(up, m, b); a = np.where(up, a, m)
    return ns, 0.5*(a + b)
