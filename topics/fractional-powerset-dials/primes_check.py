import math
N=10**7
# smallest prime factor sieve
spf=list(range(N+1))
for i in range(2,int(N**0.5)+1):
    if spf[i]==i:
        for j in range(i*i,N+1,i):
            if spf[j]==j: spf[j]=i
primes=[p for p in range(2,N+1) if spf[p]==p]
print("primes",len(primes))
# zeta identity
for s,target in [(2,15/math.pi**2),(3,1.2020569031595942/1.0173430619844491)]:
    pr=1.0
    for p in primes: pr*=1+p**-s
    print("s",s,pr,target)
# counts: squarefree n<=X with weights
import array
def stats(Xs):
    # iterate n, factor via spf
    acc={k:0.0 for k in ['sf','corner','centre','unionhalf']}
    out=[]
    xi=0
    for n in range(1,N+1):
        m=n; sf=True; ok1=True; w_c=1.0; w_u=1.0; last=0
        while m>1:
            p=spf[m]; m//=p
            if p==last: sf=False; break
            last=p
            if p%4!=1: ok1=False
            w_c*=0.5
            w_u*=p*(math.sqrt(1+1/p)-1)
        if sf:
            acc['sf']+=1; acc['centre']+=w_c; acc['unionhalf']+=w_u
            if ok1: acc['corner']+=1
        if n==Xs[xi]:
            out.append((n,dict(acc))); xi+=1
            if xi==len(Xs): break
    return out
Xs=[10**k for k in range(3,8)]
res=stats(Xs)
for X,a in res:
    L=math.log(X)
    print(X, {k:round(v/(X/math.sqrt(L)),4) for k,v in a.items() if k!='sf'}, 'sf/X',round(a['sf']/X,5))
# local exponent between consecutive decades
for (X1,a1),(X2,a2) in zip(res,res[1:]):
    for k in ['corner','centre','unionhalf']:
        e=1+ (math.log(a2[k]/X2)-math.log(a1[k]/X1))/(math.log(math.log(X2))-math.log(math.log(X1)))
        print(X1,X2,k,"alpha est",round(e,3))
