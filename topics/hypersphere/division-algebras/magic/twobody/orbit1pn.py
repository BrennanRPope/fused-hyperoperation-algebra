"""Integrate the 1PN two-body equations (pn1.py Lagrangian, relative motion) and measure the turn per lap."""
import numpy as np, sympy as sp
from scipy.integrate import solve_ivp
import pn1
rs, rds, phs, phds = sp.symbols('r_ rd_ ph_ phd_')
def stu_params(kA, sA, kB, sB, TA, TB):
    """Strengths (G = total mass = 1 units) for two extremal STU bodies, and the twice-sourced terms C from the probe turns."""
    mA, mB = sum(kA), sum(kB); wA = [sp.S(k)/mA for k in kA]; wB = [sp.S(k)/mB for k in kB]
    s = 4*sum(a*b for a, b in zip(wA, wB)) - 1
    v = 4*sum(a*b*x*y for a, b, x, y in zip(wA, wB, sA, sB))
    C = lambda T: 7 + 6*s - s**2 + v**2 - 8*v - T*(1 + s - v)**2
    return dict(s=s, v=v, CA=C(TA), CB=C(TB), XA=sp.S(mA)/(mA + mB))
def measure(par, Jscale=40.0, ecc=0.3, laps=6):
    nu = par['XA']*(1 - par['XA'])
    L1 = pn1.L1_rel(nu, par['XA']).subs({pn1.s: par['s'], pn1.v: par['v'], pn1.CA: par['CA'], pn1.CB: par['CB']})
    L1 = L1.subs({pn1.r: rs, pn1.rd: rds, pn1.ph: phds})
    K = float(1 + par['s'] - par['v'])
    Lag = (rds**2 + rs**2*phds**2)/2 + K/rs + L1
    # Euler-Lagrange: d/dt dL/dq' - dL/dq = 0, linear in accelerations
    rdd, phdd = sp.symbols('rdd phdd')
    def ddt(expr):
        return sp.diff(expr, rs)*rds + sp.diff(expr, rds)*rdd + sp.diff(expr, phds)*phdd
    eqs = [ddt(sp.diff(Lag, rds)) - sp.diff(Lag, rs), ddt(sp.diff(Lag, phds))]
    sol = sp.solve(eqs, [rdd, phdd], dict=True)[0]
    f = sp.lambdify((rs, rds, phds), [sol[rdd], sol[phdd]], 'numpy')
    J = Jscale*K; p = J**2/K; r0 = p/(1 + ecc)
    y0 = [r0, 0.0, 0.0, J/r0**2]
    def rhs(t, y):
        a, b = f(y[0], y[1], y[3]); return [y[1], a, y[3], b]
    Tk = 2*np.pi*(p/(1 - ecc**2))**1.5/np.sqrt(K)
    ev = lambda t, y: y[1]; ev.direction = 1   # rdot crosses 0 upward: periastron
    out = solve_ivp(rhs, [0, (laps + 0.5)*Tk], y0, rtol=1e-12, atol=1e-12, events=ev, dense_output=True)
    keep = out.t_events[0] > 0.1*Tk; ph_peri = np.concatenate([[0.0], out.y_events[0][keep, 2]])
    adv = np.mean(np.diff(ph_peri)) - 2*np.pi
    Jact = float(sp.diff(Lag, phds).subs({rs: r0, rds: 0, phds: y0[3]}))
    return adv/(np.pi*K**2/Jact**2), out
