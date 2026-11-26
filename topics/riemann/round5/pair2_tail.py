import math, numpy as np, itertools
def qint(n,q): return sum(q**i for i in range(n))
def qfact(n,q):
    o=1.0
    for i in range(1,n+1): o*=qint(i,q)
    return o
def qmulti(es,q): return qfact(sum(es),q)/math.prod(qfact(e,q) for e in es)
N=800_000
spf=np.zeros(N+1,np.int64)
sv=np.ones(N+1,bool);sv[:2]=False
for i in range(2,int(N**.5)+1):
    if sv[i]: sv[i*i::i]=False
P=np.nonzero(sv)[0]
for p in P[::-1]: spf[p::p]=p
def evec(n):
    es=[]
    while n>1:
        p=spf[n];e=0
        while n%p==0: n//=p;e+=1
        es.append(e)
    return es
form={'0.3':1.6786848175642866,'0.6':1.7239302392989948,'0.9':1.791972311037411}
for q in [0.3,0.6,0.9]:
    W=np.zeros(N+1); 
    for n in range(2,N+1): W[n]=qmulti(evec(n),q)
    n=np.arange(N+1,dtype=float); n[0]=1
    terms=W*n**-2.0
    for M in [200_000,400_000,800_000]:
        d=1+terms[2:M+1].sum(); half=terms[M//2+1:M+1].sum()  # tail ~ last-octave sum (density ~ constant x log-ish)
        print(q,M,d,d+half,form[str(q)], (d+half)/form[str(q)]-1)
