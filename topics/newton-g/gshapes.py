"""Newton's G lead, reopened (Nov 25, 2026). Shared numbers for checks.py and figures.py. numpy only."""
import numpy as np
from math import pi, log, exp

# ---------- Litim cutoff, lambda = 0 (from ../why-i/e-i-pi/asymptotic-safety/beyond-one-loop/flows.py) ----------
B1, B2 = -11/(3*pi), 1/(12*pi)
gstar = 2/(2*B2 - B1)                     # = 12 pi / 23
theta_litim = 23/11
def litim0(g): return (2 + g*B1/(1 - g*B2))*g

# ---------- one-loop shape and Sakharov / Bonanno-Reuter closed form ----------
def one_loop(g, gs=1.0, eps=2.0): return eps*g*(1 - g/gs)
def G_sakharov(k, G0, w): return 1/(1/G0 + w*k**2)      # 1/G = 1/G0 + w k^2 (additive induced stiffness)

# ---------- 2+eps gravity, Aida-Kitazawa (1997) eq. 4.18 ----------
def ak_beta(G, eps, c=0.0):
    return eps*G - (25 - c)/6*(G/(4*pi))*G - 5*(25 - c)/12*(G/(4*pi))**2*G
def rho_gravity(c=0.0): return 15/(25 - c)

# ---------- critical branching: F(x) = f(1+x) - 1 for offspring pgf f ----------
LAWS = {
    "geometric (linear fractional)": (lambda x: x/(1 - x), 2.0, 6.0),     # f = 1/(2-s): sigma^2 = 2, f'''(1) = 6
    "Poisson(1)  [E(x) = e^x - 1]":  (lambda x: np.expm1(x), 1.0, 1.0),   # f = e^(s-1)
    "binary split (0 or 2)":          (lambda x: x + x**2/2, 1.0, 0.0),    # f = (1+s^2)/2
}
def rho_branching(sig2, f3): return 1 - 2*f3/(3*sig2**2)

def bisect(fun, a, b, n=200):
    fa = fun(a)
    for _ in range(n):
        m = 0.5*(a + b); fm = fun(m)
        if np.sign(fm) == np.sign(fa): a, fa = m, fm
        else: b = m
    return 0.5*(a + b)

def theta_map(F, eps):
    """Near-critical family x -> F(m x), m = e^eps. Gaussian multiplier m (rate eps per step).
    Returns -log of the multiplier at the other fixed point (its exponent)."""
    m = exp(eps)
    h = lambda x: F(m*x) - x
    xs = bisect(h, -0.999/m if F is LAWS["geometric (linear fractional)"][0] else -50.0, -1e-14)
    dF = (F(m*xs + 1e-9) - F(m*xs - 1e-9))/2e-9
    return -log(m*dF)

def theta_flow(beta, eps, a, b):
    gs = bisect(lambda g: beta(g, eps), a, b)
    d = (beta(gs + 1e-9, eps) - beta(gs - 1e-9, eps))/2e-9
    return -d

# ---------- shapes: induced share R = g/g* along the flow ----------
def R_logistic(t):  return 1/(1 + np.exp(-2*t))
def R_richards(t, p): return (1/(1 + np.exp(-2*p*t)))**(1/p)
def litim_t_of_R(R, th=theta_litim):        # exact: R^(1/2) (1-R)^(-1/th) = e^t
    return 0.5*np.log(R) - np.log(1 - R)/th
def R_litim(t, th=theta_litim):
    out = np.empty_like(np.asarray(t, float))
    for i, ti in enumerate(np.atleast_1d(t)):
        out.flat[i] = bisect(lambda R: litim_t_of_R(R, th) - ti, 1e-300, 1 - 1e-16)
    return out
def align(Rfun, t):          # shift so R = 1/2 at t = 0
    t0 = bisect(lambda s: Rfun(np.array([s]))[0] - 0.5, -20, 20)
    return Rfun(t + t0)
