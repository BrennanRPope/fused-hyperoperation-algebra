"""Katz-Pavlovic / Obukhov dyadic model in the Cheskidov normalization.
a_j' = -nu*lam^(2*gam*j)*a_j + lam^(j-1)*a_{j-1}^2 - lam^j*a_j*a_{j+1},  j=0..N
lam = 2^theta, theta = 1 + n/2 for an n-dimensional NS-like model (fully intermittent
Bernstein scaling). For true dissipation |xi|^(2 alpha): gam = alpha/theta.
"""
import numpy as np
from scipy.integrate import solve_ivp

def rhs_factory(N, lam, nu, gam):
    lj = lam ** np.arange(N + 1)
    diss = nu * lam ** (2 * gam * np.arange(N + 1))
    def f(t, a):
        ap = np.zeros_like(a)
        ap[1:] += lj[:-1] * a[:-1] ** 2
        ap[:-1] -= lj[:-1] * a[:-1] * a[1:]
        ap -= diss * a
        return ap
    def jac(t, a):
        J = np.zeros((N + 1, N + 1))
        idx = np.arange(N + 1)
        J[idx, idx] = -diss
        J[idx[:-1], idx[:-1]] -= lj[:-1] * a[1:]
        J[idx[:-1], idx[1:]] -= lj[:-1] * a[:-1]
        J[idx[1:], idx[:-1]] += 2 * lj[:-1] * a[:-1]
        return J
    return f, jac

def run(N, lam, nu=0.0, gam=0.5, T=5.0, a0=None, method='DOP853', rtol=1e-10, atol=1e-30, max_pts=200000):
    if a0 is None:
        a0 = np.zeros(N + 1); a0[0] = 1.0
    f, jac = rhs_factory(N, lam, nu, gam)
    kw = dict(rtol=rtol, atol=atol, dense_output=False)
    if method in ('Radau', 'BDF', 'LSODA'):
        kw['jac'] = jac
    sol = solve_ivp(f, (0, T), a0, method=method, **kw)
    return sol
