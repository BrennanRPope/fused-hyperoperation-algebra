"""Exact checks (plain Python integers) for every list the Lean files state."""
import json
N=30
def mul(p,q,n=None):
    r=[0]*(len(p)+len(q)-1)
    for i,a in enumerate(p):
        for j,b in enumerate(q): r[i+j]+=a*b
    return r[:n] if n else r
def trim(p):
    while p and p[-1]==0: p=p[:-1]
    return p
out={}
# Euler product
e=[1]+[0]*(N-1)
for k in range(1,N):
    f=[1]+[0]*(k-1)+[-1]; e=mul(e,f,N)
pent=[0]*N
for k in range(-10,11):
    g=k*(3*k-1)//2
    if 0<=g<N: pent[g]+=(-1)**k
assert e==pent; out['euler30']=e
# partitions to 50
M=50
p=[1]+[0]*M
for k in range(1,M+1):
    for i in range(k,M+1): p[i]+=p[i-k]
out['p']=p
assert all(p[5*n+4]%5==0 for n in range(10)); assert all(p[7*n+5]%7==0 for n in range(7)); assert all(p[11*n+6]%11==0 for n in range(5))
# check series product
pp=mul(e,p[:N],N); assert pp==[1]+[0]*(N-1)
# theta
th=[0]*N
for k in range(-6,7):
    if k*k<N: th[k*k]+=1
out['theta']=th
th4=mul(mul(th,th,N),mul(th,th,N),N)
r4=[1]+[8*sum(d for d in range(1,n+1) if n%d==0 and d%4) for n in range(1,N)]
assert th4==r4; out['theta4']=th4
# triple product
tp=[1]+[0]*(N-1)
for n in range(1,N):
    if 2*n<N: tp=mul(tp,[1]+[0]*(2*n-1)+[-1],N)
    if 2*n-1<N:
        f=[1]+[0]*(2*n-2)+[1]; tp=mul(mul(tp,f,N),f,N)
assert tp==th, (tp,th)
# gauss binomials
from functools import lru_cache
@lru_cache(None)
def gauss(n,k):
    if k==0: return (1,)
    if n==0: return ()
    a=list(gauss(n-1,k-1)); b=[0]*k+list(gauss(n-1,k))
    L=max(len(a),len(b)); a+= [0]*(L-len(a)); b+=[0]*(L-len(b))
    return tuple(trim([x+y for x,y in zip(a,b)]))
out['gauss']={f"{n},{k}":list(gauss(n,k)) for n in range(7) for k in range(n+1)}
print('gauss 4 2',gauss(4,2), 'at q=2', sum(c*2**i for i,c in enumerate(gauss(4,2))))
print('euler',e); print('p',p); print('theta4',th4)
json.dump(out,open('series_results.json','w'))
print('all series checks passed')
