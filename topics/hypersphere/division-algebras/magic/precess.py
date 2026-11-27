"""Turn per lap of a probe orbiting a static source, first order, from the exact background.
W(u) = (C^2/B)[(E - q A_t)^2/A - m(phi)^2 - L^2 u^2/C],  u = 1/r,  dphi = L du/sqrt(W).
Series: F = W/L^2 = c0 + c1 u - c2 u^2 + c3 u^3;  turn = 2pi(1/sqrt(c2)-1) + (3pi/2) c3 c1 / c2^(5/2).
Reported as T = turn / (pi K^2/L^2), where K is the Newtonian pull (c1 = 2K/L^2), so GR gives 6 and Sommerfeld 1."""
import sympy as sp, mpmath as mp
from fields import ghs, jnw, r, th
u,L,E=sp.symbols('u L E',positive=True)
def W_of(g,phi,At,q,b):
    A=-g[0,0]; B=g[1,1]; C=g[2,2]/r**2; m=sp.exp(b*phi)
    W=(C**2/B)*((E-q*At)**2/A-m**2-L**2/(r**2*C))
    return sp.simplify(W.subs(r,1/u))
def T_series(W):
    F=sp.series((W/L**2).subs(E,1),u,0,4).removeO()
    c1=sp.limit(F.coeff(u,1)*L**2,L,sp.oo)/2   # = K
    c2=-F.coeff(u,2); c3=sp.limit(F.coeff(u,3),L,sp.oo)
    d=sp.limit((c2-1)*L**2,L,sp.oo)          # c2 = 1 + d/L^2
    K=sp.simplify(c1)
    turn_coef=sp.simplify(-d + 3*c3*K)        # turn = pi*(-d + 3 c3 K)/L^2
    return sp.nsimplify(sp.simplify(turn_coef/K**2)), K
def T_numeric(W,K,Lv=300,e=0.3):
    Wf=sp.lambdify((u,E,L),W,'mpmath'); mp.mp.dps=40
    Kf=float(K); pp=mp.mpf(Lv)**2/Kf; u1=(1+e)/pp; u2=(1-e)/pp
    E0=mp.findroot(lambda Ev: Wf(u1,Ev,Lv), 1-Kf/(2*pp)); E0=mp.re(E0)
    u2=mp.re(mp.findroot(lambda uu: Wf(uu,E0,Lv), u2))
    lo,hi=min(u1,u2),max(u1,u2); mid=(lo+hi)/2; half=(hi-lo)/2
    def f(t):
        w=mp.re(Wf(mid-half*mp.cos(t),E0,Lv))
        return Lv*half*mp.sin(t)/mp.sqrt(w) if w>0 else mp.mpf(1)
    turn=2*mp.quad(f,[0,mp.pi])-2*mp.pi
    return float(turn/(mp.pi*Kf**2/Lv**2))
