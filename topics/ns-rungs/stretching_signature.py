"""Does blow-up in the shell models look like a tetration-type flow or a self-feeding pole?
Test on the arrival time t_n of each shell (k_n = lam^n):
  pole / self-similar blow-up: t* - t_n = C r^n  (straight line on a log scale, constant ratio r < 1)
  rung-3 flow (Y' = Y log Y):  steps shrink like 1/n   (ratio -> 1, no finite t*)
  rung-4 flow (tetration):     steps shrink like 1/(n log n) (ratio -> 1, no finite t*)
We fit both a geometric model and a tetration-type model to the late shells and compare errors."""
import numpy as np, json, math
from multiprocessing import Pool
from helical_shell import run
from sabra_check import run as sabra_run
def slog(x):
    n=0
    while x>1: x=math.log(x); n+=1
    return n+(x-1)
def fits(tn, lam, n0, n1):
    n=np.arange(len(tn)); T=np.array([np.nan if v is None else v for v in tn],float)
    T=np.fmax.accumulate(np.where(np.isnan(T),-np.inf,T))   # arrival front
    m=(n>=n0)&(n<=n1)&np.isfinite(T); n,T=n[m],T[m]
    # geometric: t_n = t* - C r^n  -> scan t*, linear fit of log(t*-t_n)
    best=None
    for ts in T[-1]+np.logspace(-7,0,400):
        y=np.log(ts-T); A=np.polyfit(n,y,1); pred=ts-np.exp(np.polyval(A,n)); e=np.sqrt(np.mean((pred-T)**2))
        if best is None or e<best[0]: best=(e,ts,math.exp(A[0]))
    eg,ts,r=best
    # tetration-type: t_n = a + b*slog(n * 2 ln lam + c), scan c
    bt=None
    for c in np.linspace(-0.9*n[0]*2*math.log(lam),20,300):
        x=np.array([slog(v*2*math.log(lam)+c) for v in n]); A=np.polyfit(x,T,1)
        if A[0]<=0: continue
        e=np.sqrt(np.mean((np.polyval(A,x)-T)**2))
        if bt is None or e<bt: bt=e
    # rung-3 type: t_n = a + b*log(n + c)
    b3=None
    for c in np.linspace(-0.9*n[0],50,300):
        x=np.log(n+c); A=np.polyfit(x,T,1)
        if A[0]<=0: continue
        e=np.sqrt(np.mean((np.polyval(A,x)-T)**2))
        if b3 is None or e<b3: b3=e
    span=T[-1]-T[0]
    return dict(geo_err=eg/span, tet_err=bt/span, r3_err=b3/span, ratio=r, tstar=ts, n=n.tolist(), T=T.tolist())
def job(a):
    kind,eps,seed=a
    if kind=='sabra':
        te,tn=sabra_run(-0.5,N=30); lam=2.0; n0,n1=8,26
    else:
        r=run(eps,N=35,tmax=40,seed=seed); tn=r['tn']; lam=(1+5**.5)/2; n0,n1=12,31
    f=fits(tn,lam,n0,n1); f.update(kind=kind,eps=eps,seed=seed); return f
if __name__=="__main__":
    jobs=[('sabra',None,0)]+[('helical',e,s) for e in (0.5,4.0,float('inf')) for s in range(4)]
    with Pool(4) as p: res=p.map(job,jobs)
    json.dump(res,open('stretching_signature.json','w'),indent=1,default=str)
    for f in res:
        print(f"{f['kind']:8s} eps={f['eps']} seed={f['seed']}: step ratio r={f['ratio']:.3f}  relative error  geometric {f['geo_err']:.4f} | tetration-type {f['tet_err']:.4f} | rung-3 type {f['r3_err']:.4f}")
