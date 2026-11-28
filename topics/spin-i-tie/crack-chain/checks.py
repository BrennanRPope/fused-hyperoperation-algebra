# DRAFT. Round 7: can the law a o_k b = a + b + k*a*b give Shannon's chain rule instead of assuming it?
# Entropy at dial k: Tsallis S_q with q = 1 - k (natural-log units), S_q(p) = (1 - sum p^q)/(q - 1); q = 1 is Shannon.
# Three ways to write "joint = first + what's left" at dial k:
#  (a) additive, conditionals weighted by p(a)^q         (escort weights)
#  (b) the law,  conditionals weighted by p(a) honestly  (law-composed, honest average)
#  (c) additive, conditionals weighted by p(a) honestly  (Shannon's form)
import json, numpy as np
rng=np.random.default_rng(7)
def S(p,q):
    p=np.asarray(p,float); p=p[p>0]
    return float(-(p*np.log(p)).sum()) if abs(q-1)<1e-12 else float((1-(p**q).sum())/(q-1))
def law(a,b,k): return a+b+k*a*b
def joint(nA=3,nB=4): P=rng.random((nA,nB)); return P/P.sum()
def parts(P,q):
    pA=P.sum(1); cond=[S(P[i]/pA[i],q) for i in range(len(pA))]
    return pA,np.array(cond),S(P.ravel(),q),S(pA,q)
qs=[0.3,0.6,1.0,1.5,2.0,3.0]; res={}; resid={}
# C1: (a) holds exactly at every dial (Furuichi 2006's chain rule)
r=[]
for q in qs:
    for _ in range(50):
        pA,c,SAB,SA=parts(joint(),q); r.append(abs(SAB-(SA+(pA**q*c).sum())))
res['C1']=max(r)<1e-12
# C2: the escort weights add up to the law's lift: sum p^q = 1 + k*S_q, k = 1 - q
r=[]
for q in qs:
    for _ in range(50):
        p=rng.random(5); p/=p.sum(); r.append(abs((p**q).sum()-(1+(1-q)*S(p,q))))
res['C2']=max(r)<1e-12
# C3: for independent A, B, (a) reduces to the law: S(AB) = S(A) o_k S(B)
r=[]
for q in qs:
    for _ in range(50):
        a=rng.random(3); a/=a.sum(); b=rng.random(4); b/=b.sum()
        r.append(abs(S(np.outer(a,b).ravel(),q)-law(S(a,q),S(b,q),1-q)))
res['C3']=max(r)<1e-12
# C4 and C5: (b) and (c) fail for dependent A, B at every dial except k = 0
grid=np.linspace(0.2,3.0,57); Ps=[joint() for _ in range(40)]
rb=[];rc=[]
for q in grid:
    eb=ec=0
    for P in Ps:
        pA,c,SAB,SA=parts(P,q); avg=(pA*c).sum()
        eb=max(eb,abs(SAB-law(SA,avg,1-q))); ec=max(ec,abs(SAB-(SA+avg)))
    rb.append(eb); rc.append(ec)
rb=np.array(rb); rc=np.array(rc); i1=np.argmin(abs(grid-1))
res['C4']=rb[i1]<1e-12 and all(rb[j]>1e-4 for j in range(len(grid)) if abs(grid[j]-1)>0.04)
res['C5']=rc[i1]<1e-12 and all(rc[j]>1e-4 for j in range(len(grid)) if abs(grid[j]-1)>0.04)
# C6: (a) chains consistently over three variables at every dial (grouping is associative)
r=[]
for q in qs:
    for _ in range(30):
        P=rng.random((2,3,3)); P/=P.sum(); pA=P.sum((1,2)); pAB=P.sum(2)
        t=S(pA,q)
        for i in range(2):
            t+=pA[i]**q*S(pAB[i]/pA[i],q)
            for j in range(3): t+=pAB[i,j]**q*S(P[i,j]/pAB[i,j],q)
        r.append(abs(S(P.ravel(),q)-t))
res['C6']=max(r)<1e-12
for k,v in res.items(): print(k,"PASS" if v else "FAIL")
print(f"{sum(res.values())}/{len(res)} passed")
json.dump({'checks':{k:bool(v) for k,v in res.items()},'q_grid':grid.tolist(),'resid_law_honest':rb.tolist(),'resid_add_honest':rc.tolist()},open('results.json','w'),indent=1)
