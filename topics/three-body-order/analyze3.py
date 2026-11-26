import json,numpy as np
from words import rho
R=json.load(open('words.json'))
T=np.array([r['Tstar'] for r in R])
K=np.zeros((len(R),5))
for i,r in enumerate(R):
    for k in r['runs']: K[i,min(k,5)-1]+=1
dbl=np.array([sum(1 for j in range(len(r['word'])) if r['word'][j]==r['word'][(j+1)%len(r['word'])]) for r in R],float)
# adjacent run-length pairs
pairs=[(a,b) for a in range(2,5) for b in range(a,5)]
P=np.zeros((len(R),len(pairs)))
for i,r in enumerate(R):
    rr=r['runs']
    for j in range(len(rr)):
        a,b=sorted((min(rr[j],4),min(rr[(j+1)%len(rr)],4)))
        if (a,b) in pairs: P[i,pairs.index((a,b))]+=1
def loo(X,y):
    n=len(y); e=np.zeros(n)
    for i in range(n):
        m=np.ones(n,bool); m[i]=False; c=np.linalg.lstsq(X[m],y[m],rcond=None)[0]; e[i]=(X[i]@c-y[i])/y[i]
    return np.median(abs(e))*100, abs(e).max()*100
lr8=np.array([np.log(rho(r['word'],np.exp(0.8j*np.pi))) for r in R])
lr1=np.array([r['ell']/2 for r in R])
print('runs           ',loo(K,T))
print('runs+doubles   ',loo(np.column_stack([K,dbl]),T))
print('runs+pairs     ',loo(np.column_stack([K,P]),T))
print('runs+logrho(.8)',loo(np.column_stack([K,lr8]),T))
print('runs+logrho(-1)',loo(np.column_stack([K,lr1]),T))
print('runs+doubles+logrho(.8)',loo(np.column_stack([K,dbl,lr8]),T))
print('runs+pairs+logrho(-1)',loo(np.column_stack([K,P,lr1]),T))
c=np.linalg.lstsq(K,T,rcond=None)[0]; tres=T-K@c
for nm,f in [('doubles',dbl),('logrho.8',lr8),('logrho-1',lr1)]:
    cc=np.linalg.lstsq(K,f,rcond=None)[0]; print(nm,'partial corr',np.corrcoef(f-K@cc,tres)[0,1])
print('doubles total',dbl.sum(),'words with doubles',(dbl>0).sum())
