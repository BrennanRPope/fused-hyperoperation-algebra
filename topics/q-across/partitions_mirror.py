"""Line-for-line Python copy of Partitions.lean's definitions, to check the stated lists."""
import json
from knots_mirror import padd, pmul, trim
def pzeros(k,p): return [0]*k+p
def one(N): return [1]+[0]*(N-1)
def tmul(N,p,q): return pmul(p,q)[:N]
def oneMinus(N,k,p): return padd(p,[-c for c in pzeros(k,p)])[:N]
def onePlus(N,k,p): return padd(p,pzeros(k,p))[:N]
def divAux(a,buf):
    out=[]
    for x in a:
        b=x+(buf[0] if buf else 0); out.append(b); buf=buf[1:]+[b]
    return out
def divOneMinus(k,p): return divAux(p,[0]*k)
def r1(N): return range(1,1+max(N-1,0))
def eulerProd(N):
    acc=one(N)
    for k in r1(N): acc=oneMinus(N,k,acc)
    return acc
def neg(k): return (-1)**k
def pentCoeff(m):
    acc=0
    for k in range(m+2):
        acc+= (neg(k) if (k*max(3*k-1,0))//2==m else 0) + (neg(k) if k>=1 and k*(3*k+1)//2==m else 0)
    return acc
def tma(N,x,y): return padd(padd(x,y),pmul(x,y))[:N]
def eulerMA(N):
    acc=[]
    for k in r1(N): acc=tma(N,acc,pzeros(k,[-1]))
    return acc
def partGen(N):
    acc=one(N)
    for k in r1(N): acc=divOneMinus(k,acc)
    return acc
def theta(N): return [1 if i==0 else (2 if any(j*j==i for j in range(i+1)) else 0) for i in range(N)]
def tripleProd(N):
    acc=one(N)
    for k in r1(N): acc=oneMinus(N,2*k,onePlus(N,2*k-1,onePlus(N,2*k-1,acc)))
    return acc
def sigmaStar(n): return sum(d for d in range(1,1+n) if n%d==0 and d%4!=0)
def r4(n): return 1 if n==0 else 8*sigmaStar(n)
R=json.load(open('series_results.json'))
assert eulerProd(30)==[pentCoeff(m) for m in range(30)]==R['euler30']
assert trim(padd([1],eulerMA(30)))==trim(eulerProd(30))
assert partGen(51)==R['p']
assert tmul(30,eulerProd(30),partGen(30))==one(30)
assert tripleProd(30)==theta(30)
t2=tmul(30,theta(30),theta(30)); assert tmul(30,t2,t2)==[r4(n) for n in range(30)]
assert [r4(p) for p in [2,3,5,7,11,13,17,19,23,29]]==[24,32,48,64,96,112,144,160,192,240]
print('Partitions.lean lists all match')
