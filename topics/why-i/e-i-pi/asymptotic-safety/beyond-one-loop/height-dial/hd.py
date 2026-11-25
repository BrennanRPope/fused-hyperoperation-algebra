"""Shared pieces for the height-dial test. numpy only.
Imports the fractional-rung machinery (regular iteration of E(x) = e^x - 1) and the published flows."""
import numpy as np, os, sys
from math import pi, e, log
HERE = os.path.dirname(os.path.abspath(__file__)); UP = os.path.dirname(HERE)
ROOT = os.path.abspath(os.path.join(HERE, "../../../../.."))          # math-theory/
sys.path[:0] = [UP, os.path.join(ROOT, "fractional-rungs")]
from regular_iteration import Es                                     # Es(s, x): E iterated s times (ours, fractional-rungs)
from flows import eh_beta, newton, jac, ak_two_loop, root, deriv, B1_0, B2_0   # published flows (see ../README.md)

def e_s(s):  return 1 + float(Es(s, np.array([1.0]))[0])             # e at height s (ours, e-every-height): e_0 = 2, e_1 = e
def wall(s): return float(Es(s, np.array([-60.0]))[0])               # finite end of the height-s law, E^s(-inf)
def V(s, a, h=1e-5):
    """Generator of B_s(a,b) = E^s(E^-s a + E^-s b): V_s(a) = (E^s)'(E^-s(a)). V_1(a) = 1 + a."""
    x = float(Es(-s, np.array([a]))[0])
    return float((Es(s, np.array([x + h]))[0] - Es(s, np.array([x - h]))[0])/(2*h))
def wall_rate(s, x1=-20.0, x2=-19.0):
    """Exponent at the wall: E^s(x) - w_s ~ C e^(mu x) as x -> -inf; returns mu."""
    w = wall(s); a = Es(s, np.array([x1, x2])); return float(np.log((a[1] - w)/(a[0] - w))/(x2 - x1))
def s_for(L, lo=-1.5, hi=3.0):
    """Height whose e_s has ln e_s = L (e_s is increasing in s)."""
    for _ in range(70):
        m = (lo + hi)/2
        if log(e_s(m)) < L: lo = m
        else: hi = m
    return (lo + hi)/2

# ---- litim cutoff, lambda = 0, in d dimensions (Reuter-Saueressig form, same as ../flows.py at d = 4) ----
from math import gamma as Gamma
def litim0_ratio(d):
    pre = (4*pi)**(1 - d/2)
    Phi  = lambda p, n, w: 1/Gamma(n + 1)/(1 + w)**p
    Phit = lambda p, n, w: 1/Gamma(n + 2)/(1 + w)**p
    b1 = pre/3*(d*(d+1)*Phi(1, d/2-1, 0) - 6*d*(d-1)*Phi(2, d/2, 0) - 4*d*Phi(1, d/2-1, 0) - 24*Phi(2, d/2, 0))
    b2 = -pre/6*(d*(d+1)*Phit(1, d/2-1, 0) - 6*d*(d-1)*Phit(2, d/2, 0))
    gs = -(d - 2)/(b1 - (d - 2)*b2)                    # eta(g*) = -(d-2)
    return 1/(1 - gs*b2)                               # theta/(d-2)
def ak_ratio(eps, c=0.0):
    f = lambda G: ak_two_loop(G, eps, c); Gs = root(lambda G: f(G)/G, 1e-9, 400.0)
    return -deriv(f, Gs)/eps
