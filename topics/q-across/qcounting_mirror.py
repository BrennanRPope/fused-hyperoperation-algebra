"""Python copy of QCounting.lean's finite definitions, to check the stated values."""
from math import comb
from knots_mirror import padd, trim
def pzeros(k,p): return [0]*k+p
def gauss(n,k):
    if k==0: return [1]
    if n==0: return []
    return padd(gauss(n-1,k-1),pzeros(k,gauss(n-1,k)))
def peval(x,p): return sum(c*x**i for i,c in enumerate(p))
assert trim(gauss(4,2))==[1,1,2,1,1] and peval(2,gauss(4,2))==35 and peval(0,gauss(4,2))==1 and peval(-1,gauss(4,2))==2
for n in range(9):
    for k in range(n+1):
        assert peval(1,gauss(n,k))==comb(n,k)
        exp=0 if (n%2==0 and k%2==1) else comb(n//2,k//2)
        assert peval(-1,gauss(n,k))==exp,(n,k)
def padd2(P,Q):
    L=max(len(P),len(Q)); P=P+[[]]*(L-len(P)); Q=Q+[[]]*(L-len(Q)); return [padd(a,b) for a,b in zip(P,Q)]
def mulLin(i,P): return padd2(P,[[]]+[pzeros(i,p) for p in P])
for n in range(7):
    P=[[1]]
    for i in range(n): P=mulLin(i,P)
    R=[pzeros(k*(k-1)//2,gauss(n,k)) for k in range(n+1)]
    assert [trim(p) for p in P]==[trim(r) for r in R]
print('QCounting.lean values all match')
