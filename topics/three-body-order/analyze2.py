import json,numpy as np
from words import rho
R=json.load(open('words.json'))
T=np.array([r['Tstar'] for r in R]); L=np.array([r['L'] for r in R]); ell=np.array([r['ell'] for r in R])
K=np.zeros((len(R),5))
for i,r in enumerate(R):
    for k in r['runs']: K[i,min(k,5)-1]+=1
def fit(X,y,m):
    c,*_=np.linalg.lstsq(X[m],y[m],rcond=None); return c
def looerr(X,y):
    n=len(y); e=np.zeros(n)
    for i in range(n):
        m=np.ones(n,bool); m[i]=False; c=fit(X,y,m); e[i]=(X[i]@c-y[i])/y[i]
    return e
ths=np.linspace(0.5,1.0,51)*np.pi
LR={}
for th in ths:
    LR[th]=np.array([np.log(rho(r['word'],np.exp(1j*th))) for r in R])
base=looerr(K,T)
print('runs only: median %.4f%% p90 %.4f%% max %.2f%%'%(np.median(abs(base))*100,np.percentile(abs(base),90)*100,abs(base).max()*100))
worst=np.argsort(-abs(base))[:6]
for i in worst: print('  worst',R[i]['cls'],R[i]['word'][:40],R[i]['runs'][:12],'%.2f%%'%(base[i]*100))
best=None
for th in ths:
    e=looerr(np.column_stack([K,LR[th]]),T)
    med=np.median(abs(e))*100
    if best is None or med<best[1]: best=(th,med,e)
    print('theta/pi %.2f med %.4f%% max %.2f%%'%(th/np.pi,med,abs(e).max()*100)) if abs((th/np.pi*100)%5)<1e-6 else None
print('best theta/pi',best[0]/np.pi,best[1])
# nested: choose theta on training only
n=len(T); en=np.zeros(n); chosen=[]
idx=np.arange(n); rng=np.random.default_rng(0); folds=np.array_split(rng.permutation(n),10)
for f in folds:
    tr=np.setdiff1d(idx,f)
    bestth=None
    for th in ths[::2]:
        X=np.column_stack([K,LR[th]])
        # inner LOO approx by 5-fold
        inner=np.array_split(rng.permutation(tr),5); errs=[]
        for g in inner:
            t2=np.setdiff1d(tr,g); m=np.zeros(n,bool); m[t2]=True; c=fit(X,T,m); errs+=list(abs((X[g]@c-T[g])/T[g]))
        s=np.median(errs)
        if bestth is None or s<bestth[1]: bestth=(th,s)
    chosen.append(bestth[0]/np.pi)
    X=np.column_stack([K,LR[bestth[0]]]); m=np.zeros(n,bool); m[tr]=True; c=fit(X,T,m); en[f]=(X[f]@c-T[f])/T[f]
print('nested 10-fold: chosen theta/pi',np.round(chosen,2),'median %.4f%% max %.2f%%'%(np.median(abs(en))*100,abs(en).max()*100))
# runs-only 10-fold same folds
eb=np.zeros(n)
for f in folds:
    tr=np.setdiff1d(idx,f); m=np.zeros(n,bool); m[tr]=True; c=fit(K,T,m); eb[f]=(K[f]@c-T[f])/T[f]
print('runs-only 10-fold median %.4f%% max %.2f%%'%(np.median(abs(eb))*100,abs(eb).max()*100))
# partial: does arrangement matter? ell residual after run counts vs period residual
for th in [np.pi, best[0]]:
    lr=LR[th]; c=np.linalg.lstsq(K,lr,rcond=None)[0]; lres=lr-K@c
    c2=np.linalg.lstsq(K,T,rcond=None)[0]; tres=T-K@c2
    print('theta/pi %.2f partial corr(log rho | runs, T | runs) = %.3f ; sd of log-rho residual %.3f'%(th/np.pi,np.corrcoef(lres,tres)[0,1],lres.std()))
json.dump(dict(best=best[0]/np.pi),open('best_theta.json','w'))
