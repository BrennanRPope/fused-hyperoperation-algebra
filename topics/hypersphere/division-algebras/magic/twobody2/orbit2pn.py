"""Integrate the 2PN centre-of-mass Hamiltonian (from pn2.py) for equal masses and measure the turn per lap."""
import sys, os; sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np, sympy as sp
from scipy.integrate import solve_ivp
from scipy.optimize import brentq
from pn2 import full_L, params, hamiltonian, G, m1, m2, r
def build(case, **kw):
    sub = params(case, **kw); sub = {k: (v.subs({G: 1, m1: 1, m2: 1}) if hasattr(v, 'subs') else v) for k, v in sub.items()}; sub.update({G: 1, m1: 1, m2: 1})
    L1, L2 = full_L(); H0, H1, H2, k = hamiltonian(L1, L2, sub)
    pr, J = sp.symbols('pr J', real=True)
    H = (H0 + H1 + H2).subs({sp.Symbol('px'): pr, sp.Symbol('py'): J/r})
    f = sp.lambdify((r, pr, J), H, 'numpy'); dr = sp.lambdify((r, pr, J), sp.diff(H, pr), 'numpy')
    dphi = sp.lambdify((r, pr, J), sp.diff(H, J), 'numpy'); dpr = sp.lambdify((r, pr, J), -sp.diff(H, r), 'numpy')
    return f, dr, dphi, dpr, float(k)
def measure(model, Jk, ecc=0.3, laps=8):
    f, dr, dphi, dpr, k = model; mu = 0.5; J = Jk*k   # so that k/J = 1/Jk, as K/L for a small visitor
    w = (1 - ecc**2)/2; E = -w*mu*k**2/J**2
    rp0 = J**2/(mu*k*(1 + ecc)); rp = brentq(lambda x: f(x, 0.0, J) - E, 0.6*rp0, 1.3*rp0)
    rhs = lambda t, y: [dr(y[0], y[1], J), dpr(y[0], y[1], J), dphi(y[0], y[1], J)]
    ev = lambda t, y: y[1]; ev.direction = 1
    Tk = 2*np.pi*(J**2/(mu*k*(1 - ecc**2)))**1.5*np.sqrt(mu/k)
    sol = solve_ivp(rhs, [0, laps*1.2*Tk], [rp, 0.0, 0.0], rtol=1e-12, atol=1e-12, events=ev, dense_output=True, max_step=Tk/400)
    te = [t for t in sol.t_events[0] if t > 0.1*Tk]; ph = [sol.sol(t)[2] for t in te]
    turns = np.diff([0.0] + ph) - 2*np.pi
    return float(np.mean(turns)), J, E, sol
