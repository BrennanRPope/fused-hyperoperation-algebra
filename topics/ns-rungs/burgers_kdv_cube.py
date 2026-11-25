"""Rung 2 (mult-add) warm-up: the two-wave cube  tau = 1 + e1 + e2 + A*e1*e2.
Every term is w*exp(a x + b t). With softmax weights p_j over the terms:
  (log tau)_x = sum p_j a_j ,  (log tau)_xx = sum p_j a_j^2 - (sum p_j a_j)^2   (exact)
Burgers (Cole-Hopf, u = -2 nu (log phi)_x): phi must solve phi_t = nu phi_xx. A term solves it
iff b = nu a^2; the corner term has a=a1+a2, b=b1+b2, so b - nu a^2 = -2 nu a1 a2 != 0.
-> corner weight A must be 0 -> shocks MERGE (a plateau is lost).
KdV (Hirota, u = 2 (log tau)_xx): corner weight A = ((k1-k2)/(k1+k2))^2 != 0 ->
solitons pass through each other and both survive, shifted."""
import numpy as np
def soft(terms):
    T=np.array(terms); M=T.max(0); W=np.exp(T-M); return W/W.sum(0)
x=np.linspace(-30,90,6001); dx=x[1]-x[0]
nu=0.05
def burgers(t,A=0.0):
    # plateaus u=2 | 1 | 0, shocks start at x=-4 and x=+4
    a=[0.0,-1/(2*nu),-2/(2*nu)]; b=[nu*ai*ai for ai in a]
    th=[0*x, a[1]*(x-4)+b[1]*t, a[2]*(x+4)+b[2]*t + a[1]*(-8)]
    terms=list(th); slopes=list(a); rates=list(b)
    if A>0: terms.append(np.log(A)+th[1]+th[2]); slopes.append(a[1]+a[2]); rates.append(b[1]+b[2])
    p=soft(terms); S=np.array(slopes)[:,None]; R=np.array(rates)[:,None]
    u=-2*nu*(p*S).sum(0)
    heat_rel=(p*(R-nu*S**2)).sum(0)          # (phi_t - nu phi_xx)/phi, exact
    return u, heat_rel
def shocks(u):
    g=np.abs(np.gradient(u,x)); i=np.where((g[1:-1]>g[:-2])&(g[1:-1]>=g[2:])&(g[1:-1]>1))[0]+1
    return [round(float(x[j]),2) for j in i]
out=[]
pr=lambda s: (print(s), out.append(s))
pr("BURGERS (nu=0.05). Shock 2|1 moves at 1.5, shock 1|0 at 0.5; they meet at t=8.")
for A in (0.0,1.0):
    _,h=burgers(3.0,A); pr(f"  corner weight A={A}: max |phi_t - nu phi_xx|/phi = {np.max(np.abs(h)):.3e}")
for t in (0,4,8,16,30):
    u,_=burgers(t); pr(f"  t={t:2d}: shocks at {shocks(u)}  ({len(shocks(u))} shock(s))")
# KdV  u_t + 6 u u_x + u_xxx = 0
k1,k2=1.0,0.6
def kdv(t,A='hirota'):
    if A=='hirota': A=((k1-k2)/(k1+k2))**2
    e1=k1*(x+10)-k1**3*t; e2=k2*x-k2**3*t
    terms=[0*x,e1,e2]; sl=[0,k1,k2]
    if A>0: terms.append(np.log(A)+e1+e2); sl.append(k1+k2)
    p=soft(terms); S=np.array(sl,float)[:,None]
    return 2*((p*S*S).sum(0)-((p*S).sum(0))**2)
def kdv_res(A='hirota',t=8.0,h=1e-4):
    u=kdv(t,A); ut=(kdv(t+h,A)-kdv(t-h,A))/(2*h)
    ux=np.gradient(u,x); uxxx=np.gradient(np.gradient(ux,x),x)
    return np.max(np.abs(ut+6*u*ux+uxxx)[20:-20])/np.max(np.abs(ut))
A12=((k1-k2)/(k1+k2))**2
pr(f"\nKdV, k1={k1}, k2={k2}: Hirota corner weight A12=((k1-k2)/(k1+k2))^2 = {A12:.4f}")
pr(f"  relative PDE residual mid-collision (t=8): corner A12: {kdv_res():.1e}   corner 1 (plain mult-add): {kdv_res(1.0):.1e}   corner 0: {kdv_res(0.0):.1e}")
def peaks(u):
    i=np.where((u[1:-1]>u[:-2])&(u[1:-1]>=u[2:])&(u[1:-1]>0.05))[0]+1
    return [(round(float(x[j]),2),round(float(u[j]),3)) for j in i]
for t in (0,8,25,50):
    pr(f"  t={t:2d}: peaks (x, height) {peaks(kdv(t))}")
pr(f"  both heights k^2/2 = {k1**2/2}, {k2**2/2} survive. Phase shifts: fast +{-np.log(A12)/k1:.3f}, slow {np.log(A12)/k2:.3f}")
open('burgers_kdv_output.txt','w').write("\n".join(out)+"\n")
np.savez('cube_waves.npz',x=x)
