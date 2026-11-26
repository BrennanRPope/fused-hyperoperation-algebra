"""DRAFT control: does WHICH run sits next to which matter, or only how many of each?"""
import json,numpy as np
R=json.load(open('words.json')); T=np.array([r['Tstar'] for r in R])
def feats(RU):
    K=np.zeros((len(RU),5)); pairs=[(a,b) for a in range(1,6) for b in range(a,6)]; P=np.zeros((len(RU),len(pairs)))
    for i,rr in enumerate(RU):
        for k in rr: K[i,min(k,5)-1]+=1
        for j in range(len(rr)):
            a,b=sorted((min(rr[j],5),min(rr[(j+1)%len(rr)],5))); P[i,pairs.index((a,b))]+=1
    return K,P
def loo(X,y):
    n=len(y); e=np.zeros(n); H=X@np.linalg.pinv(X); r=y-H@y; e=r/(1-np.diag(H))/y
    return float(np.median(abs(e))*100), float(abs(e).max()*100)
RU=[r['runs'] for r in R]; K,P=feats(RU)
print('counts only',loo(K,T)); print('counts + neighbour pairs (real order)',loo(np.column_stack([K,P]),T))
rng=np.random.default_rng(3); sh=[]
for rep in range(20):
    K2,P2=feats([list(rng.permutation(rr)) for rr in RU]); sh.append(loo(np.column_stack([K,P2]),T)[0])
print('counts + neighbour pairs (runs shuffled within word) median: mean %.4f min %.4f max %.4f'%(np.mean(sh),min(sh),max(sh)))
json.dump(dict(counts=loo(K,T),real=loo(np.column_stack([K,P]),T),shuffled=sh),open('pairs_ctrl.json','w'))
