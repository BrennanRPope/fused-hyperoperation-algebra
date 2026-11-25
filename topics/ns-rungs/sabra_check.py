"""Control run: classic one-chain Sabra model (a=1, b, c=-1-b), lambda=2, same integrator,
to see whether the arrival-time test detects finite-time blow-up where it is known to occur."""
import numpy as np, sys
from scipy.integrate import solve_ivp
def run(b,N=30,tmax=6,theta=1e-2,lam=2.0):
    k=lam**np.arange(N); a=1.0; c=-1-b
    c0,c1,c2=a*lam,b*lam,-c*lam
    def rhs(t,y):
        u=y[:N]+1j*y[N:]; du=np.zeros(N,complex)
        X,Y,Z=u[:-2],u[1:-1],u[2:]; kk=k[:-2]
        du[:-2]+=1j*kk*c0*np.conj(Y)*Z; du[1:-1]+=1j*kk*c1*np.conj(X)*Z; du[2:]+=1j*kk*c2*X*Y
        return np.concatenate([du.real,du.imag])
    r=np.random.default_rng(0); u=np.zeros(N,complex); u[:3]=np.array([1,.7,.5])*np.exp(2j*np.pi*r.random(3)); u[3:]=1e-12*np.exp(2j*np.pi*r.random(N-3))
    tgt=N-3
    def ev(t,y): return abs(y[tgt]+1j*y[N+tgt])**2-theta*k[tgt]**(-2/3)
    ev.terminal=True; ev.direction=1
    s=solve_ivp(rhs,(0,tmax),np.concatenate([u.real,u.imag]),method='DOP853',rtol=1e-9,atol=1e-14,events=ev,max_step=0.02)
    En=s.y[:N]**2+s.y[N:]**2
    tn=[next((float(s.t[i]) for i in range(len(s.t)) if En[n,i]>=theta*k[n]**(-2/3)),None) for n in range(N)]
    return s.t[-1],tn
for b in [float(v) for v in sys.argv[1:]]:
    te,tn=run(b); print(f"b={b}: t_end={te:.4f}  t_n=",[None if v is None else round(v,4) for v in tn])
