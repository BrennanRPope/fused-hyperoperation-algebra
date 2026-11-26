"""DRAFT: neighbour corners (order-sensitive) vs full circle operation of all runs (order-free) vs free pair table."""
import json,numpy as np
from scipy.optimize import least_squares
R=json.load(open('words.json'))
T=np.array([r['Tstar'] for r in R]); RU=[[min(k,5) for k in r['runs']] for r in R]
def pred(par,i,kind):
    p=par[:5]; kap=par[5] if len(par)>5 else 0; rr=RU[i]; v=np.array([p[k-1] for k in rr])
    if kind=='sum': return v.sum()
    if kind=='nbr': return v.sum()+kap*np.sum(v*np.roll(v,-1)) if len(v)>1 else v.sum()
    if kind=='alln': # all pairs, order-free, second order only
        return v.sum()+kap*(v.sum()**2-(v**2).sum())/2
    if kind=='nbr_shuffled':
        return v.sum()+kap*np.sum(v*np.roll(v,-1))
def fit(idx,kind,rr=None):
    x0=np.array([1.4,4.6,7.4,10.0,12.4]+([0.0] if kind!='sum' else []))
    f=lambda par: np.array([(pred(par,i,kind)-T[i])/T[i] for i in idx])
    return least_squares(f,x0).x
n=len(T); folds=np.array_split(np.random.default_rng(1).permutation(n),20)
out={}
for kind in ['sum','nbr','alln']:
    e=np.zeros(n)
    for f in folds:
        par=fit(np.setdiff1d(np.arange(n),f),kind); e[f]=[(pred(par,i,kind)-T[i])/T[i] for i in f]
    out[kind]=dict(median=float(np.median(abs(e))*100),p90=float(np.percentile(abs(e),90)*100),max=float(abs(e).max()*100),err=e.tolist(),par=fit(np.arange(n),kind).tolist())
    print(kind,{k:v for k,v in out[kind].items() if k not in('err',)})
# control: shuffle run order within each word, refit nbr: if adjacency matters, shuffled should do worse
rng=np.random.default_rng(7); orig=[list(r) for r in RU]
meds=[]
for rep in range(5):
    for i in range(n): RU[i]=list(rng.permutation(orig[i]))
    e=np.zeros(n)
    for f in folds:
        par=fit(np.setdiff1d(np.arange(n),f),'nbr'); e[f]=[(pred(par,i,'nbr')-T[i])/T[i] for i in f]
    meds.append(float(np.median(abs(e))*100))
RU[:]=orig
out['nbr_shuffled_medians']=meds; print('shuffled-order neighbour corners medians',meds)
json.dump(out,open('corner2.json','w'))
