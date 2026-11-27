"""Hardy Z(t) for Riemann zeta at large t by the Riemann-Siegel formula (main sum + C0, C1 terms),
used for Gram statistics of zeta at heights whose zero spacing matches the torus window."""
import numpy as np
from zetas import theta_R

def Z_rs(t):
    t = np.atleast_1d(np.asarray(t, dtype=float))
    a = np.sqrt(t/(2*np.pi)); N = np.floor(a).astype(int); p = a - N
    th = theta_R(t)
    out = np.zeros_like(t)
    Nmax = N.max()
    for n in range(1, Nmax + 1):
        m = N >= n
        out[m] += np.cos(th[m] - t[m]*np.log(n))/np.sqrt(n)
    out *= 2
    c0 = np.cos(2*np.pi*(p*p - p - 1/16))/np.cos(2*np.pi*p)
    # C1 = -Psi'''(p)/(96 pi^2), Psi(p) = c0; numerical third derivative
    h = 1e-3
    def psi(q): return np.cos(2*np.pi*(q*q - q - 1/16))/np.cos(2*np.pi*q)
    d3 = (psi(p + 2*h) - 2*psi(p + h) + 2*psi(p - h) - psi(p - 2*h))/(2*h**3)
    c1 = -d3/(96*np.pi**2)
    out += (-1.0)**(N - 1)*(t/(2*np.pi))**(-0.25)*(c0 + c1/a)
    return out
