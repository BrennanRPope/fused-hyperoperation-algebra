import math
def li2(z,terms=200000):
    return sum(z**k/k**2 for k in range(1,terms))
# naturals additive dial: log GF of prod(1+c q^n) ~ L(c)/(1-q), L(c) = -Li2(-c)
full=-li2(-1.0); print("L(1)",full, math.pi**2/12)
target=full/2
lo,hi=0,1
for _ in range(60):
    m=(lo+hi)/2
    if -li2(-m,5000)<target: lo=m
    else: hi=m
print("participation c giving half of log-mass:",lo, " L(1/2)=",-li2(-0.5,5000)," vs half",target)
# union-centre dials on naturals: c_n solves (1+c q^n)^2=1+q^n  -> log GF exactly half
# corner: odd parts, L = pi^2/24
print("odd-parts L", math.pi**2/24)
# rationals: sum over coprime a,b of (ab)^-2 = zeta(2)^2/zeta(4)=2.5
s=0.0; M=3000
from math import gcd
for a in range(1,M):
    for b in range(1,M//a+1 if a<M else 1):
        pass
# faster: prod over primes (1+p^-s)/(1-p^-s)
N=10**6; sieve=bytearray([1])*(N+1); sieve[0]=sieve[1]=0
for i in range(2,int(N**.5)+1):
    if sieve[i]: sieve[i*i::i]=bytearray(len(sieve[i*i::i]))
pr=[p for p in range(N+1) if sieve[p]]
P=1.0
for p in pr: P*=(1+p**-2)/(1-p**-2)
print("Q>0 Euler product at s=2",P," zeta(2)^2/zeta(4)=",(math.pi**2/6)**2/(math.pi**4/90))
# Gaussian prime ideals: split p=1 mod 4 -> two ideals norm p; inert p=3 mod4 -> one ideal norm p^2; 2 ramified -> norm 2
P=1.0
for p in pr:
    if p==2: P*=1+2**-2
    elif p%4==1: P*=(1+p**-2)**2
    else: P*=1+p**-4
# zeta_K(s)=zeta(s)L(s,chi4); L(2,chi4)=Catalan G
G=0.915965594177219
zK2=math.pi**2/6*G
L4=sum((-1)**k/(2*k+1)**4 for k in range(200000))
zK4=math.pi**4/90*L4
print("Q(i) prime ideals product at s=2",P," zeta_K(2)/zeta_K(4)=",zK2/zK4)
