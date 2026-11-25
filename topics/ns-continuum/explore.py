import numpy as np, time
from scipy.integrate import solve_ivp
from dyadic import rhs_factory
n=3.0; lam=2**(1+n/2); N=28
f,_=rhs_factory(N,lam,0.0,0.5)
a0=np.zeros(N+1); a0[0]=1.0
jt=N-4
ev=lambda t,a: a[jt]-1e-2*lam**(-jt/3); ev.terminal=True
t0=time.time()
sol=solve_ivp(f,(0,10),a0,method='DOP853',rtol=1e-11,atol=1e-60,events=ev)
print('time',time.time()-t0,'npts',len(sol.t),'tend',sol.t[-1])
t=sol.t; A=sol.y
for frac in [0.5,0.8,0.9,0.95,0.98,0.99,0.995,0.999,1.0]:
    k=np.searchsorted(t,frac*t[-1]); k=min(k,len(t)-1)
    prof=np.log(np.abs(A[:,k])+1e-300)/np.log(lam)
    print(f"t={t[k]:.10f}", ' '.join(f"{p:6.2f}" for p in prof[:N-2]))
