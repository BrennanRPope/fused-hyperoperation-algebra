"""Molecular biology through mult-add: checks. Run: python3 -I checks.py"""
import numpy as np, itertools, math
rng=np.random.default_rng(1)
ok=0; tot=0
def check(name,cond,detail=""):
    global ok,tot; tot+=1; ok+=bool(cond); print(("PASS " if cond else "FAIL ")+name+(" | "+detail if detail else ""))
wadd=lambda k,a,b: a+b+k*a*b
ma=lambda a,b: a+b+a*b
inv=lambda u: -u/(1+u)

# ---- 1. Helix-coil (Zimm-Bragg) is the dial gas with k = 1/sigma, z = sigma*s
def zb_ring(s,sig,N):
    M=np.array([[s,1.0],[sig*s,1.0]])   # rows: previous (h,c); cols: next (h,c)
    return np.trace(np.linalg.matrix_power(M,N))
def dial_ring(z,k,N):
    T=np.array([[1.0,z],[1.0,k*z]])
    return np.trace(np.linalg.matrix_power(T,N))
err=0
for _ in range(200):
    s=rng.uniform(0.2,3); sig=10**rng.uniform(-4,0); N=int(rng.integers(2,40))
    a=zb_ring(s,sig,N); b=dial_ring(sig*s,1/sig,N); err=max(err,abs(a-b)/abs(a))
check("Zimm-Bragg ring = dial-gas ring with k=1/sigma, z=sigma*s (200 random cases)",err<1e-10,f"max rel err {err:.1e}")
# brute force: ZB weights by enumeration on an open chain vs matrix
def zb_brute(s,sig,N):
    Z=0
    for st in itertools.product([0,1],repeat=N):  # 1 = helix
        w=1.0
        for i,x in enumerate(st):
            if x==1:
                w*=s*(sig if (i==0 or st[i-1]==0) else 1.0)
        Z+=w
    return Z
def zb_open(s,sig,N):
    M=np.array([[s,1.0],[sig*s,1.0]]); v=np.array([0.0,1.0])  # start as if preceded by coil
    for _ in range(N): v=v@M
    return v.sum()
e=max(abs(zb_brute(1.3,0.01,n)-zb_open(1.3,0.01,n))/zb_open(1.3,0.01,n) for n in range(1,13))
check("Zimm-Bragg transfer matrix matches brute-force sum over all chains (N<=12)",e<1e-12,f"{e:.1e}")
# odd rings vanish on the wall z=-1/k  <=> s = -1
v=[abs(zb_ring(-1.0,sig,N)) for sig in (1e-4,1e-2,0.3) for N in (3,5,7,9,11)]
check("every odd Zimm-Bragg ring vanishes at s = -1 (the wall), any sigma",max(v)<1e-9,f"max {max(v):.1e}")
vv=[abs(zb_ring(-1.0,0.01,N)) for N in (2,4,6)]
check("even rings do not vanish there",min(vv)>0.1,str(np.round(vv,3)))
# sharpness: helix fraction slope at s=1 for long chain = 1/(4 sqrt(sigma))
for sig in (1e-2,1e-3,1e-4):
    lam=lambda s: ((1+s)+math.sqrt((1-s)**2+4*sig*s))/2
    h=1e-6; th=lambda s: s*(math.log(lam(s+h))-math.log(lam(s-h)))/(2*h)
    slope=(th(1+1e-4)-th(1-1e-4))/2e-4
    check(f"helix-fraction slope at s=1 is 1/(4 sqrt sigma), sigma={sig}",abs(slope-1/(4*math.sqrt(sig)))/slope<1e-3,f"{slope:.2f} vs {1/(4*math.sqrt(sig)):.2f}")
print(" measured-scale dial values: sigma 1e-2..1e-4 -> k = 1/sigma = 100..10000")

# ---- 2. Two-site binding: Z = 1 + a + b + w ab = 1 + (a (+)_w b)
for w in (0,0.5,1,4,25):
    a,b=rng.uniform(0,3,2)
    Z=1+a+b+w*a*b
    check(f"two-site binding polynomial = 1 + wAdd(w,a,b), w={w}",abs(Z-(1+wadd(w,a,b)))<1e-12)
check("w=1 (independent sites): Z factorises (1+a)(1+b)",abs(1+wadd(1,.7,1.9)-(1.7*2.9))<1e-12)
check("w=0 (exclusive sites): Z = 1+a+b, plain addition",abs(1+wadd(0,.7,1.9)-3.6)<1e-12)
# Hill coefficient at half saturation for identical sites
def hill_half(w):
    Z=lambda x:1+2*x+w*x*x; Y=lambda x:(2*x+2*w*x*x)/(2*Z(x))
    from scipy.optimize import brentq
    x50=brentq(lambda x:Y(x)-.5,1e-9,1e9); h=1e-6
    f=lambda lx: math.log(Y(math.exp(lx))/(1-Y(math.exp(lx))))
    return (f(math.log(x50)+h)-f(math.log(x50)-h))/(2*h), x50
for w in (0.25,1,4,100):
    nh,x50=hill_half(w)
    check(f"two-site Hill coefficient = 2 sqrt(w)/(1+sqrt(w)), w={w}",abs(nh-2*math.sqrt(w)/(1+math.sqrt(w)))<1e-6,f"nH={nh:.4f}, x50={x50:.4f}=1/sqrt(w)")
# wall: weight -1/w on one site makes Z forget the other
check("a site at the wall -1/w makes Z independent of the other site",abs((1+wadd(4,0.3,-0.25))-(1+wadd(4,7.0,-0.25)))<1e-12)

# ---- 3. MWC hemoglobin = sum of two mult-add powers (dial 1 for R, dial c for T)
L,c,KR,n=9000,0.014,2.8,4   # LibreTexts teaching fit
def Y(p):
    a=p/KR
    return (a*(1+a)**(n-1)+L*c*a*(1+c*a)**(n-1))/((1+a)**n+L*(1+c*a)**n)
mpow=lambda a,n:(1+a)**n-1
def wpow(k,a,n):
    x=0.0
    for _ in range(n): x=wadd(k,x,a)
    return x
a=1.7
check("(1+c a)^4 = 1 + c*wpow(c,a,4): T state is the weighted power at dial c",abs((1+c*a)**4-(1+c*wpow(c,a,4)))<1e-12)
check("(1+a)^4 = 1 + mpow(a,4): R state is the plain mult-add power",abs((1+a)**4-(1+mpow(a,4)))<1e-12)
from scipy.optimize import brentq
p50=brentq(lambda p:Y(p)-.5,0.1,500)
h=1e-5; lg=lambda lp: math.log(Y(math.exp(lp))/(1-Y(math.exp(lp))))
nH=(lg(math.log(p50)+h)-lg(math.log(p50)-h))/(2*h)
check("MWC (L=9000,c=0.014,K_R=2.8 torr) gives P50 near human 26-27 torr",20<p50<32,f"P50={p50:.1f} torr")
check("MWC gives Hill coefficient near the fitted 2.8",2.5<nH<3.2,f"nH={nH:.2f}")
# limits: L->0 pure R, L->inf pure T: no cooperativity
for LL in (1e-9,1e12):
    L0=L; L=LL
    p=brentq(lambda p:Y(p)-.5,1e-3,1e6); nh=(lg(math.log(p)+h)-lg(math.log(p)-h))/(2*h)
    check(f"one dial alone (L={LL:g}) has Hill coefficient 1: cooperativity needs two dials mixed",abs(nh-1)<1e-3,f"nH={nh:.4f}")
    L=L0

# ---- 4. Fidelity: catches combine by the 'or' law (k=-1); Eigen threshold = mult-add inverse
from fractions import Fraction as F
eps=[F(1,10**5),F(1,10**2),F(1,10**3)]   # base selection, proofreading, mismatch repair (textbook orders)
catch=[1-e for e in eps]; tot_c=catch[0]
for x in catch[1:]: tot_c=wadd(-1,tot_c,x)
check("three fidelity stages: total catch = c1 (+)_-1 c2 (+)_-1 c3, error 1e-10",(1-tot_c)==F(1,10**10),f"error {float(1-tot_c):.2e} (exact rational)")
# Eigen: master grows iff sigma*(1-mu)^L > 1  <=>  s (*) l > 0, l = mpow(-mu,L)
for _ in range(500):
    mu=10**rng.uniform(-6,-2); Lg=int(10**rng.uniform(2,5)); sig=1+10**rng.uniform(-3,1)
    l=mpow(-mu,Lg); s=sig-1
    lhs=sig*(1-mu)**Lg-1; rhs=ma(s,l)
    assert abs(lhs-rhs)<1e-9*max(1,abs(lhs))
check("Eigen: master's excess growth = s (*) l with l = mpow(-mu, L) (500 cases)",True)
l=mpow(-1e-4,10000); sth=inv(l)
check("threshold advantage is the mult-add inverse of the copy loss",abs(ma(sth,l))<1e-12,f"l={l:.4f}, s*={sth:.4f}")
for name,Ug in (("DNA microbes (Drake ~0.003/genome)",0.003),("lytic RNA viruses (~1/genome)",1.0)):
    l=math.exp(-Ug)-1
    print(f"   {name}: copy loss l={l:.4f}, needs advantage s > {inv(l):.4f} (sigma > {1+inv(l):.3f})")
check("Drake: DNA microbes sit far from the wall (need sigma>1.003)",abs((1+inv(math.exp(-0.003)-1))-math.exp(0.003))<1e-12)
check("Drake: RNA viruses sit near it (need sigma>e=2.718)",abs(1+inv(math.exp(-1)-1)-math.e)<1e-12)
print(f"\n{ok}/{tot} checks passed")
