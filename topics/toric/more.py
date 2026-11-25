import sympy as sp, math, numpy as np, itertools, json
from sympy import primerange, divisor_sigma, totient, n_order
out={}
# 1. elliptic curves: wobble a_p within 2 sqrt p
def ec_count(p,A,B):
    sq={}
    for y in range(p): sq[y*y%p]=sq.get(y*y%p,0)+1
    return 1+sum(sq.get((x**3+A*x+B)%p,0) for x in range(p))
ec=[]
for p in primerange(5,110):
    n1=ec_count(p,-1,0); n2=ec_count(p,1,1)
    ec.append((p,p+1-n1,p+1-n2,round(2*math.sqrt(p),2)))
out['ec']=ec
print('a_p for y^2=x^3-x and y^2=x^3+x+1 (p, a1, a2, bound):'); print(ec)
# 2. Circle n-th power vs ordinary power mod n (AKS identity)
a=sp.symbols('a')
aks=[]
for n in range(2,200):
    ok=all(math.comb(n,k)%n==0 for k in range(1,n))
    aks.append(ok)
assert [n for n,ok in zip(range(2,200),aks) if ok]==list(primerange(2,200))
print('AKS check: n | every middle Pascal entry exactly for primes up to 199: ok')
out['pascal_mod']={n:[math.comb(n,k)%n for k in range(n+1)] for n in [4,5,6,7,8,9]}
print(out['pascal_mod'])
# 3. q -> 1 limits
q=sp.symbols('q')
def qint(n): return sum(q**i for i in range(n))
def qfact(n): return sp.prod([qint(i) for i in range(1,n+1)])
def qbin(n,k): return sp.cancel(qfact(n)/(qfact(k)*qfact(n-k)))
lim=[]
for n,k in [(3,1),(4,2),(5,2),(6,3)]:
    b=sp.expand(qbin(n,k)); lim.append((n,k,str(b),b.subs(q,1)))
for n in range(1,6):
    GL=sp.prod([q**n-q**i for i in range(n)])
    lim.append(('GL',n,str(sp.factor(GL)),sp.limit(GL/(q-1)**n,q,1)))
print(lim)
# Circle-operation coordinates: P^n count = ((1+a)^(n+1)-1)/a
for n in range(1,5):
    e=sp.expand(((1+a)**(n+1)-1)/a); print('P^%d'%n,e, ' at a=0:',e.subs(a,0))
# 4. sublattices of Z^2, Z^3 by Hermite normal form, vs zeta products
def count_hnf(n,m):
    # number of upper-triangular HNF matrices with det m: diag d1..dn, entries above diag in [0,d_j)
    tot=0
    for ds in itertools.product(*[sp.divisors(m)]*n):
        if math.prod(ds)!=m: continue
        c=1
        for j,d in enumerate(ds): c*=d**j   # column j has j entries above, each in [0,d_j)
        tot+=c
    return tot
def coeff(n,m):  # Dirichlet coefficient of zeta(s)...zeta(s-n+1) at m
    # multiplicative convolution of n^0, n^1, ..., n^(n-1)
    f={1:1}
    for i in range(n):
        g={}
        for x,v in f.items():
            for k in range(1,m//x+1):
                if x*k<=m: g[x*k]=g.get(x*k,0)+v*k**i
        f=g
    return f.get(m,0)
lat={}
for n in [1,2,3]:
    row=[(m,count_hnf(n,m),coeff(n,m)) for m in range(1,31)]
    assert all(x==y for _,x,y in row); lat[n]=[r[1] for r in row]
    print('Z^%d sublattices index 1..30:'%n,lat[n])
print('prime index in Z^2: p+1, Z^3: p^2+p+1 ->',[(p,lat[2][p-1],lat[3][p-1]) for p in [2,3,5,7,11,13]])
out['lat']=lat
# 5. Frobenius on roots of unity; Dirichlet counts
for nn in [7,12]:
    for p in [2,3,5,7,11,13,29,31]:
        if nn%p==0: print(nn,p,'ramified'); continue
        o=n_order(p,nn); print('n=%d p=%d p mod n=%d order %d splits into %d'%(nn,p,p%nn,o,totient(nn)//o))
ps=list(primerange(2,10**6))
for nn in [12,7,3]:
    c={}
    for p in ps: c[p%nn]=c.get(p%nn,0)+1
    print('primes < 1e6 mod',nn,sorted(c.items()))
    out['dir%d'%nn]=sorted(c.items())
# 6. zeros for toric zeta over Z
z=np.load('/mnt/project-files/math-theory/primes/zeta_zeros_1000.npy')
out['zeros']=[round(float(x),3) for x in z[:14]]
print(out['zeros'])
# check: local factor of P1xP1 over Z at p equals product of shifted zeta local factors
s=sp.symbols('s')
for p in [2,3,5]:
    Np=lambda m:(p**m+1)**2
    # log Z_p(T) = sum N_m T^m/m with T=p^-s ; compare to -sum h_i log(1-p^i T), h=(1,2,1)
    T=sp.Rational(1,7)
    lhs=sum(sp.Rational(Np(m),m)*T**m for m in range(1,60))
    rhs=-sum(h*sp.log(1-p**i*T) for i,h in enumerate([1,2,1])) if p*p*T<1 else None
    if rhs is not None: print('P1xP1 local zeta p=%d T=1/7:'%p, float(lhs), float(rhs))
json.dump(out,open('more.json','w'),default=str)
