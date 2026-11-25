"""Published flows of the dimensionless Newton coupling g = G k^(d-2) (and lambda = Lambda/k^2),
shared by checks.py and figures.py. numpy only. Sources in README.md."""
import numpy as np
from math import gamma as Gamma, pi

def rk4(f, v, dt, n):
    v = np.array(v, float); out = [v.copy()]
    for _ in range(n):
        k1 = f(v); k2 = f(v + dt/2*k1); k3 = f(v + dt/2*k2); k4 = f(v + dt*k3)
        v = v + dt/6*(k1 + 2*k2 + 2*k3 + k4); out.append(v.copy())
    return np.array(out)

# ---- Einstein-Hilbert truncation, optimised (Litim) cutoff, Reuter-Saueressig form (as in ../checks.py) ----
d = 4; pre = (4*pi)**(1 - d/2)
Phi  = lambda p, n, w: 1/Gamma(n + 1)/(1 + w)**p
Phit = lambda p, n, w: 1/Gamma(n + 2)/(1 + w)**p
B1 = lambda l: pre/3*(d*(d+1)*Phi(1, d/2-1, -2*l) - 6*d*(d-1)*Phi(2, d/2, -2*l) - 4*d*Phi(1, d/2-1, 0) - 24*Phi(2, d/2, 0))
B2 = lambda l: -pre/6*(d*(d+1)*Phit(1, d/2-1, -2*l) - 6*d*(d-1)*Phit(2, d/2, -2*l))

def eh_beta(X=0.0, Y=0.0):
    """EH (g, lambda) flow. Optional matter: X = N_S + 2N_D - 4N_V adds X/(6 pi) to B1,
    Y = N_S - 4N_D + 2N_V adds Y g/(4 pi) to beta_lambda (coefficients of Dona-Eichhorn-Percacci 2013,
    eqs. 35-36; the hybrid with the Litim graviton part is illustrative)."""
    def beta(v):
        g, l = v; b1 = B1(l) + X/(6*pi); e = g*b1/(1 - g*B2(l))
        bl = -(2 - e)*l + 0.5*g*pre*(2*d*(d+1)*Phi(1, d/2, -2*l) - 8*d*Phi(1, d/2, 0) - d*(d+1)*e*Phit(1, d/2, -2*l)) + Y*g/(4*pi)
        return np.array([(d - 2 + e)*g, bl])
    return beta

def jac(beta, v, h=1e-7):
    J = np.zeros((2, 2))
    for j in range(2):
        dv = np.zeros(2); dv[j] = h; J[:, j] = (beta(v + dv) - beta(v - dv))/(2*h)
    return J

def newton(beta, v, n=80):
    v = np.array(v, float)
    for _ in range(n): v = v - np.linalg.solve(jac(beta, v), beta(v))
    return v

# ---- one-dimensional flows of g ----
def one_loop(g, gs=1.0, eps=2.0): return eps*g*(1 - g/gs)                       # Weinberg 1979 shape; Bonanno-Reuter with omega' = omega
def br_full(g, w, wp, b2): return 2*g*(1 - wp*g)/(1 - b2*g)                      # Bonanno-Reuter 2000, beta = 2g(1 - omega' g)/(1 - B2 g)
B1_0, B2_0 = B1(0.0), B2(0.0)                                                    # Litim, lambda = 0: -11/(3 pi), 1/(12 pi)
def litim0(g): return (2 + g*B1_0/(1 - g*B2_0))*g
def ak_two_loop(G, eps, c=0.0):                                                  # Aida-Kitazawa 1997, eq. (4.18)
    return eps*G - (25 - c)/6*(G/(4*pi))*G - 5*(25 - c)/12*(G/(4*pi))**2*G

def root(f, a, b, n=200):
    fa = f(a)
    for _ in range(n):
        m = 0.5*(a + b); fm = f(m)
        if np.sign(fm) == np.sign(fa): a, fa = m, fm
        else: b = m
    return 0.5*(a + b)
def deriv(f, x, h=1e-6): return (f(x + h) - f(x - h))/(2*h)
