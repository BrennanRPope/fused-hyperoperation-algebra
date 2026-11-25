import numpy as np
from scipy.optimize import minimize_scalar
rng=np.random.default_rng(7); res=[]
def check(name,ok,info=""):
    res.append(ok); print(("PASS " if ok else "FAIL ")+name+(" | "+info if info else ""))
wadd=lambda k,a,b:a+b+k*a*b
lens=lambda k,q,a,b:(a+b+k*a*b)/(1-q*a*b)

# 1. Haldane: Poisson crossovers, recombination = odd number of crossovers
def simulate_r(lengths,m=0,N=400000):
    """Chi-square model: chiasma points Poisson rate 2(m+1) per Morgan on the 4-strand bundle;
    every (m+1)-th point is a crossover (random start); each crossover involves the gamete with prob 1/2.
    Returns recombination fractions for consecutive intervals and their unions."""
    L=sum(lengths); edges=np.cumsum([0]+list(lengths))
    counts=np.zeros((N,len(lengths)),dtype=np.int64)
    rate=2*(m+1)
    for i in range(N):
        n=rng.poisson(rate*L); pts=np.sort(rng.uniform(0,L,n))
        start=rng.integers(0,m+1); co=pts[start::m+1]
        co=co[rng.random(len(co))<0.5]   # crossover seen by this gamete
        counts[i]=np.histogram(co,bins=edges)[0]
    return counts
lengths=[0.15,0.25]
C=simulate_r(lengths,m=0,N=200000)
odd=lambda x:(x%2==1)
r1=odd(C[:,0]).mean(); r2=odd(C[:,1]).mean(); r12=odd(C[:,0]+C[:,1]).mean()
check("no interference (Poisson crossovers): r₁₂ = r₁ + r₂ − 2r₁r₂, the dial at k = −2",abs(r12-wadd(-2,r1,r2))<0.004,
      f"simulated {r12:.4f}, dial {wadd(-2,r1,r2):.4f}")
x=np.linspace(0.01,3,50); r=(1-np.exp(-2*x))/2
check("Haldane's map function r = (1 − e^(−2x))/2 makes 1 − 2r multiply: u = −2r adds by plain ∘₁",
      np.allclose(-2*r[:-1]*0+(-2*((1-np.exp(-2*(x[:-1]+x[1:])))/2)),wadd(1,-2*r[:-1],-2*r[1:])))
check("the wall of dial −2 is r = 1/2: adding an unlinked interval leaves the pair unlinked",all(abs(wadd(-2,a,0.5)-0.5)<1e-15 for a in (0,0.1,0.33,0.49)))

# 2. Interference: coincidence c gives dial k = -2c
for m in (1,4):
    C=simulate_r(lengths,m=m,N=200000)
    r1=odd(C[:,0]).mean(); r2=odd(C[:,1]).mean(); r12=odd(C[:,0]+C[:,1]).mean()
    dbl=(odd(C[:,0])&odd(C[:,1])).mean(); c=dbl/(r1*r2)
    check(f"with interference (chi-square model m = {m}): r₁₂ = r₁ ∘_k r₂ with k = −2c, by bookkeeping (c the coincidence)",abs(r12-wadd(-2*c,r1,r2))<1e-12,
          f"c = {c:.3f}, k = {-2*c:.3f}")
# c depends on interval size: no single dial
def coincidence(m,l,N=150000):
    C=simulate_r([l,l],m=m,N=N); a=odd(C[:,0]); b=odd(C[:,1]); return (a&b).mean()/(a.mean()*b.mean())
cs=[coincidence(4,l) for l in (0.1,0.25,0.5)]
check("under interference the coincidence, and so the dial, changes with interval length (m = 4)",max(cs)-min(cs)>0.3,
      "c at 10, 25, 50 cM: "+", ".join(f"{c:.2f}" for c in cs))

# 3. Lens laws for map functions: 'unlinked absorbs' forces 2k + q = -4
ok=True
for K in np.linspace(-0.5,1.5,9):
    k,q=-2*K,-4+4*K
    ok&=all(abs(lens(k,q,a,0.5)-0.5)<1e-12 for a in (0.05,0.2,0.4))
check("a lens law (a + b + k·ab)/(1 − q·ab) keeps r = 1/2 absorbing exactly when 2k + q = −4",ok)
k,q=-1.0,2.0
bad=lens(k,q,0.3,0.5)
check("off that line r = 1/2 is not absorbing (example k = −1, q = 2)",abs(bad-0.5)>0.01,f"gives {bad:.3f}")
check("Haldane is the point (k, q) = (−2, 0) on the line",2*-2+0==-4)
R1,R2=rng.uniform(0,0.98,2); r1,r2=R1/2,R2/2
check("Kosambi's rule r₁₂ = (r₁ + r₂)/(1 + 4r₁r₂) is the point (0, −4)",abs(lens(0,-4,r1,r2)-(r1+r2)/(1+4*r1*r2))<1e-15)
# associativity of the lens line
a,b,c_=rng.uniform(0,0.45,3); errs=[]
for K in (0,0.3,0.7,1):
    k,q=-2*K,-4+4*K; errs.append(abs(lens(k,q,lens(k,q,a,b),c_)-lens(k,q,a,lens(k,q,b,c_))))
check("every lens law on that line is associative, so it defines a map function",max(errs)<1e-14,f"{max(errs):.1e}")
w=lambda R:2*R/(1-R)
R12=(R1+R2)/(1+R1*R2)
check("Kosambi is plain ∘₁ after the change w = 2R/(1 − R), R = 2r",abs(w(R12)-wadd(1,w(R1),w(R2)))<1e-10*(1+w(R12)))
# small-interval coincidence on the line is K = -k/2
K=0.4;k,q=-2*K,-4+4*K; e=1e-5; r12=lens(k,q,e,e); cc=(2*e-r12)/(2*e*e)
check("on the line the coincidence of two tiny intervals is K = −k/2 (Haldane 1, Kosambi 0)",abs(cc-K)<1e-3,f"{cc:.4f}")
# does the chi-square model follow a lens law? fit K
def mapfun_r(m,l,N=150000):
    return None
pairs=[(0.1,0.1),(0.1,0.3),(0.2,0.2),(0.3,0.5),(0.5,0.5),(0.05,0.6),(0.15,0.15),(0.4,0.4)]
obs=[]
for l1,l2 in pairs:
    C=simulate_r([l1,l2],m=4,N=300000)
    obs.append((odd(C[:,0]).mean(),odd(C[:,1]).mean(),odd(C[:,0]+C[:,1]).mean()))
obs=np.array(obs)
f=lambda K: np.sum((lens(-2*K,-4+4*K,obs[:,0],obs[:,1])-obs[:,2])**2)
best=minimize_scalar(f,bounds=(-4,1.5),method="bounded")
miss=lambda K: np.max(np.abs(lens(-2*K,-4+4*K,obs[:,0],obs[:,1])-obs[:,2]))
noise=np.sqrt(0.25/300000)
print(f"   chi-square m=4 vs lens line: Kosambi (K=0) worst miss {miss(0):.4f}; best K={best.x:.2f} worst miss {miss(best.x):.4f}; noise per r about {noise:.4f}")
check("the chi-square interference model (m = 4) is not exactly a lens law: Kosambi, its small-interval match, misses by many noise widths",miss(0)>5*noise,
      f"Kosambi miss {miss(0):.4f}, noise {noise:.4f}")
GENE_FIT=dict(K=float(best.x),miss_best=float(miss(best.x)),miss_kos=float(miss(0)))
# 4. Gene–environment: factor V Leiden and oral contraceptives (Vandenbroucke 1994 counts)
cells={'both':(25,2),'fvl':(10,4),'oc':(84,63),'none':(36,100)}
def ors(c):
    ref=c['none'][0]/c['none'][1]
    return {k:(v[0]/v[1])/ref for k,v in c.items()}
o=ors(cells); e1,e2,e12=o['oc']-1,o['fvl']-1,o['both']-1
kd=(e12-e1-e2)/(e1*e2)
check("the excess risks combine on the dial: OR₁₁ − 1 = (OR₁₀ − 1) ∘_k (OR₀₁ − 1)",abs(wadd(kd,e1,e2)-e12)<1e-12,
      f"ORs {o['oc']:.1f}, {o['fvl']:.1f}, {o['both']:.1f} give k = {kd:.2f}")
check("k = 0 is additive risk (RERI = 0), k = 1 multiplicative; RERI = k·ERR₁·ERR₂",abs((e12-e1-e2)-kd*e1*e2)<1e-12)
# bootstrap uncertainty of k
ks=[]
for _ in range(20000):
    c={}
    for kk,(a,b) in cells.items():
        c[kk]=(rng.binomial(155,a/155),rng.binomial(169,b/169))
    if min(min(v) for v in c.values())==0: continue
    oo=ors(c); a1,a2,a12=oo['oc']-1,oo['fvl']-1,oo['both']-1
    if a1*a2!=0: ks.append((a12-a1-a2)/(a1*a2))
lo,hi=np.percentile(ks,[2.5,97.5])
check("but the joint cell has only 2 controls: the dial's 95% range spans both additive and multiplicative",lo<0 and hi>1,
      f"k = {kd:.2f}, bootstrap 95% {lo:.1f} to {hi:.1f}")
# 5. affected children
check("at least one affected child of two carriers among n: 1 − (3/4)ⁿ = wpow(−1, 1/4, n)",
      all(abs((1-0.75**n)-(lambda n:(lambda f:f(f,n))(lambda f,n:0 if n==0 else wadd(-1,f(f,n-1),0.25)))(n))<1e-15 for n in range(1,8)))
print(f"\n{sum(res)}/{len(res)} checks passed")
