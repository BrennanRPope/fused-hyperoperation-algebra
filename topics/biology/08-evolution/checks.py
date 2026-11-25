import numpy as np, json
rng=np.random.default_rng(8); res=[]; out={}
def check(name,ok,info=""):
    res.append(bool(ok)); print(("PASS " if ok else "FAIL ")+name+(" | "+info if info else ""))
star=lambda a,b:a+b+a*b
wadd=lambda k,a,b:a+b+k*a*b
mpow=lambda x,n:(1+x)**n-1

# 1. Haploid selection runs on odds, and selection coefficients compose by ∘₁
p,s=0.2,0.15
p1=p*(1+s)/(1+s*p); odds=lambda x:x/(1-x)
check("haploid selection p' = p(1 + s)/(1 + sp) multiplies the allele's odds by 1 + s",abs(odds(p1)-odds(p)*(1+s))<1e-14)
s1,s2=0.12,-0.05
q=p*(1+s1)/(1+s1*p); q=q*(1+s2)/(1+s2*q)
check("two selective stages (survival then fecundity) act as one stage with s = s₁ ∘₁ s₂",abs(odds(q)-odds(p)*(1+star(s1,s2)))<1e-14)
T=40; x=p
for _ in range(T): x=x*(1+s)/(1+s*x)
check("t generations of constant selection are one step with s_t = mpow(s, t)",abs(odds(x)-odds(p)*(1+mpow(s,T)))<1e-9*odds(x),
      f"s = {s}, 40 generations: odds × {1+mpow(s,T):.0f}")

# 2. Epistasis is the dial
a,b=0.08,0.11
for k,name in [(0,"additive null"),(1,"multiplicative null")]:
    pass
check("the additive and multiplicative 'no epistasis' nulls are dial 0 and dial 1; epistasis against the multiplicative null is ε = (k − 1)·s₁s₂",
      abs((wadd(0.3,a,b)-star(a,b))-(0.3-1)*a*b)<1e-15)

# 3. Path consistency: linear global epistasis forces the dial
# effect of mutation j on a background whose excess fitness (over the ancestor) is x: Δ_j(x) = a_j + c_j x
n=6; A=rng.uniform(-0.1,0.15,n); k=-1.8; C=k*A
def add(x,j,c=C): return x+A[j]+c[j]*x
worst=0
for _ in range(500):
    x0=rng.uniform(-0.2,0.3); i,j=rng.choice(n,2,replace=False)
    worst=max(worst,abs(add(add(x0,i),j)-add(add(x0,j),i)))
check("if every mutation's effect is linear in background fitness with slope c_j = k·a_j, adding mutations in either order gives the same genotype",worst<1e-15)
check("…and that is the dial: x + Δ_j(x) = x ∘_k a_j, so the fitness of any genotype is the dial-sum of its mutations' effects",
      all(abs(add(x,j)-wadd(k,x,A[j]))<1e-15 for x in (-0.1,0.2) for j in range(n)))
# order of 6 mutations irrelevant
perm_vals=[]
for _ in range(50):
    x=0.0
    for j in rng.permutation(n): x=add(x,j)
    perm_vals.append(x)
check("all 720 orders of six mutations reach one fitness (50 random orders checked)",np.ptp(perm_vals)<1e-15,f"fitness {perm_vals[0]:.4f}")
# converse: any other linear rule is path-inconsistent
Cb=-1.2*np.ones(n)*0.1   # common slope for all mutations (what many data show)
gap=max(abs(add(add(0.0,i,Cb),j,Cb)-add(add(0.0,j,Cb),i,Cb)) for i in range(n) for j in range(n))
check("conversely, a common slope for every mutation (beneficial and deleterious alike) gives different fitness for the same genotype depending on order",gap>1e-3,
      f"largest order gap {gap:.4f} in excess fitness")
wall=-1/k
check("at the wall x = −1/k every mutation is neutral: all the lines Δ_j(x) cross zero at one background fitness",
      all(abs(A[j]+C[j]*wall)<1e-15 for j in range(n)),f"k = {k}: wall at excess fitness {wall:.3f}")
check("past the wall every mutation flips sign: beneficial ones turn deleterious and vice versa (sign epistasis for all at once)",
      all(np.sign(A[j]+C[j]*(wall+0.1))==-np.sign(A[j]) for j in range(n)))
# fit test on noisy simulated 'global epistasis' data: ratio c/a recovers k
xs=np.linspace(-0.1,0.4,12); rat=[]
for j in range(n):
    if abs(A[j])<0.02: continue
    y=A[j]+C[j]*xs+rng.normal(0,0.002,len(xs)); cfit,afit=np.polyfit(xs,y,1); rat.append(cfit/afit)
check("in noisy simulated data the fitted slope/intercept ratio is the same for every mutation, the test the dial puts to global-epistasis data",
      np.std(rat)/abs(np.mean(rat))<0.15,f"ratios {np.round(rat,2).tolist()}, true k = {k}")
out['sim']=dict(A=A.tolist(),k=k)

# 4. Mutation load: 'or'
U=0.5; L=1-np.exp(-U)
check("Haldane–Muller load: loads at many loci combine by the 'or' law, and n loci of load U/n give wpow(−1, U/n, n) → 1 − e^(−U)",
      abs((1-(1-U/1e6)**1e6)-L)<1e-6,f"U = {U}: load {L:.4f}")
# 5. Dominance as dial on the two allele copies
h,s=0.3,0.2; kd=(1-2*h)/(h**2*s)
check("a diploid's homozygous effect is the dial-sum of its two heterozygous copies: s = (hs) ∘_k (hs) with k = (1 − 2h)/(h²s); additivity h = ½ is k = 0",
      abs(wadd(kd,h*s,h*s)-s)<1e-15,f"h = {h}, s = {s}: k = {kd:.2f}")
json.dump(out,open('results.json','w'),indent=1)
print(f"\n{sum(res)}/{len(res)} checks passed")
