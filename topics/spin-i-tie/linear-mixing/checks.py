# DRAFT. Round 8: does linear mixing pin the information dial to k = 0?
# Test: "an irrelevant coin flip doesn't change what's left to learn."
# leftover(B|A) := S(AB) - S(A)  (no chain-rule form assumed).
# Global coin: A' = (A, C), C independent of everything.  Local coin: split one outcome a into two
# halves with the same state of B (the mixture of B is unchanged; that is what linear mixing says).
# Tsallis S_q composes independent pieces by the law at k = 1 - q; Renyi H_a adds (k = 0) but is not Shannon.
import json, numpy as np
from scipy.linalg import eigvalsh
rng=np.random.default_rng(8)
def ts(p,q):
    p=np.asarray(p,float).ravel(); p=p[p>1e-300]
    return float(-(p*np.log(p)).sum()) if abs(q-1)<1e-12 else float((1-(p**q).sum())/(q-1))
def re(p,a):
    p=np.asarray(p,float).ravel(); p=p[p>1e-300]
    return float(-(p*np.log(p)).sum()) if abs(a-1)<1e-12 else float(np.log((p**a).sum())/(1-a))
def joint(nA=3,nB=4): P=rng.random((nA,nB)); return P/P.sum()
def left(P,f,q): return f(P,q)-f(P.sum(1),q)
def gcoin(P,c): return np.einsum('ab,c->acb',P,c).reshape(-1,P.shape[1])   # A' = (A, C)
def lcoin(P,t=0.5,i=0): return np.vstack([P[:i],t*P[i:i+1],(1-t)*P[i:i+1],P[i+1:]])  # split outcome i
res={}; Ps=[joint() for _ in range(40)]; coins=[np.array([x,1-x]) for x in (0.5,0.3,0.1)]
# G1: global coin multiplies the Tsallis leftover by the law's lift 1 + k*S(C), exactly
r=[]
for q in [0.3,0.6,1.0,1.5,2.0,3.0]:
    for P in Ps:
        for c in coins: r.append(abs(left(gcoin(P,c),ts,q)-(1+(1-q)*ts(c,q))*left(P,ts,q)))
res['G1']=max(r)<1e-12
# G2-G4: scan the dial; worst change of the leftover under each coin
grid=np.linspace(0.2,3.0,57); g_ts=[];l_ts=[];g_re=[];l_re=[]
for q in grid:
    g_ts.append(max(abs(left(gcoin(P,c),ts,q)-left(P,ts,q)) for P in Ps for c in coins))
    l_ts.append(max(abs(left(lcoin(P,t),ts,q)-left(P,ts,q)) for P in Ps for t in (0.5,0.2)))
    g_re.append(max(abs(left(gcoin(P,c),re,q)-left(P,re,q)) for P in Ps for c in coins))
    l_re.append(max(abs(left(lcoin(P,t),re,q)-left(P,re,q)) for P in Ps for t in (0.5,0.2)))
i1=int(np.argmin(abs(grid-1))); away=[j for j in range(len(grid)) if abs(grid[j]-1)>0.04]
res['G2']=g_ts[i1]<1e-12 and all(g_ts[j]>1e-4 for j in away)     # global coin pins the law's dial to k = 0
res['G3']=max(g_re)<1e-12                                          # Renyi passes the global coin at every order
res['G4']=l_re[i1]<1e-12 and all(l_re[j]>1e-4 for j in away) and l_ts[i1]<1e-12 and all(l_ts[j]>1e-4 for j in away)  # local coin leaves only Shannon
# G5: concavity does not pick the dial: Tsallis is concave at every q > 0 (random midpoints)
ok=True
for q in grid:
    for _ in range(200):
        p=rng.dirichlet(np.ones(5)); s=rng.dirichlet(np.ones(5))
        if ts((p+s)/2,q) < (ts(p,q)+ts(s,q))/2-1e-12: ok=False
res['G5']=ok
# G6: quantum. cq-state rho_AB = sum_a p_a |a><a| (x) rho_a. Von Neumann leftover is unchanged by both coins;
# quantum Tsallis leftover is multiplied by 1 + k S(C) under the global coin.
def rdm(d=3):
    X=rng.normal(size=(d,d))+1j*rng.normal(size=(d,d)); r=X@X.conj().T; return r/np.trace(r).real
def cq(p,rhos): return np.concatenate([pa*eigvalsh(r) for pa,r in zip(p,rhos)])  # spectrum of rho_AB
r6=[]
for q in [0.5,1.0,2.0]:
    for _ in range(20):
        p=rng.dirichlet(np.ones(3)); rhos=[rdm() for _ in range(3)]; c=coins[1]
        L=ts(cq(p,rhos),q)-ts(p,q)
        pg=np.outer(p,c).ravel(); Lg=ts(cq(pg,[r for r in rhos for _ in c]),q)-ts(pg,q)
        pl=np.concatenate([[p[0]/2,p[0]/2],p[1:]]); Ll=ts(cq(pl,[rhos[0]]+rhos),q)-ts(pl,q)
        if q==1.0: r6+= [abs(Lg-L),abs(Ll-L)]
        else: r6.append(abs(Lg-(1+(1-q)*ts(c,q))*L))
res['G6']=max(r6)<1e-10
for k,v in res.items(): print(k,"PASS" if v else "FAIL")
print(f"{sum(res.values())}/{len(res)} passed")
json.dump({'checks':{k:bool(v) for k,v in res.items()},'q_grid':grid.tolist(),'global_tsallis':g_ts,'local_tsallis':l_ts,'global_renyi':g_re,'local_renyi':l_re},open('results.json','w'),indent=1)
