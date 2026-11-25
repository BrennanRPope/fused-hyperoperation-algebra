"""Scan the cross-handed dial eps (and N, seeds) in the helical shell model.
Blow-up time t*: the time energy reaches the Kolmogorov level at shell N-3 (None if not by tmax)."""
import numpy as np, json, sys, itertools
from multiprocessing import Pool
from helical_shell import run
def job(a):
    eps,N,seed,tmax=a
    r=run(eps,N=N,tmax=tmax,seed=seed)
    return dict(eps=eps,N=N,seed=seed,tmax=tmax,tstar=(r['tend'] if r['status']==1 else None),
                dE=abs(r['E1']-r['E0']),dH=abs(r['H1']-r['H0']),
                last_shell=max([n for n,v in enumerate(r['tn']) if v is not None]))
if __name__=="__main__":
    which=sys.argv[1]
    if which=="eps":
        jobs=[(e,30,s,40.0) for e in (0.0,0.05,0.1,0.25,0.5,1.0,2.0,4.0,float('inf')) for s in range(6)]
    else:
        jobs=[(e,N,s,40.0) for e in (0.5,4.0) for N in (20,25,30,35) for s in range(4)]
    with Pool(8) as p: res=p.map(job,jobs)
    json.dump(res,open(f"scan_{which}.json","w"),indent=1,default=str)
    for r in res: print(r)
