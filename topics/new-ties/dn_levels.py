"""Level-n estimates D_n for the three no-shared-wall examples (distinct composite maps of n letters, derivative at
an attracting wall). Slow O(1/n) convergence; used only as a cross-check of the proved brackets."""
import json, numpy as np
from fractions import Fraction as F
from scipy.optimize import brentq
import planelaws as P
out={}
for q in (2,3,12):
    maps=[P.map_with(F(1),F(q),s)[0] for s in (F(1,5),F(3,20))]+[P.map_with(F(-1),F(-q),F(1,4))[0]]
    mf=[np.array([[float(x) for x in r] for r in M]) for M in maps]
    def key(A): return tuple(np.round((A/A[1,1]).ravel(),10))
    lv={key(np.eye(2)):np.eye(2)}; rows=[]
    for n in range(1,15):
        new={}
        for A in lv.values():
            for g in mf:
                B=g@A; new.setdefault(key(B),B)
        lv=new
        g=np.array([abs(np.linalg.det(A))/(A[1,0]+A[1,1])**2 for A in lv.values()])
        rows.append((n,len(lv),brentq(lambda D: np.sum(g**D)-1,1e-6,3)))
    out[str(q)]=rows; print(q,[(n,round(d,4)) for n,c,d in rows])
json.dump(out,open("dn_levels.json","w"),indent=1)
