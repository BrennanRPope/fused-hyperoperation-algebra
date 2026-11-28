# DRAFT. Round 12: how far from 0 can the law's k sit, given measured three-slit limits?
# Sorkin parameter (Sinha et al. 2010): kappa = I3 / delta, delta = |I_AB| + |I_BC| + |I_CA| (pairwise terms).
# Measured caps: Sinha et al., Science 329, 418 (2010): |kappa| < ~1e-2.
#                Kauten et al., New J. Phys. 19, 033017 (2017): more than 4 orders of magnitude below pairwise, ~1e-4.
# The law as a measure (Sugeno lambda-measure): I2 = k*a*b, I3 = k^2*a*b*c (round 11).
# Reading A (the law alone runs the odds): kappa = |k| abc / (ab + bc + ca)  -> equal slits: kappa = |k| w / 3.
# Reading B (the law rides on top of ordinary quantum interference): delta comes from quantum terms
#   |I_ij| = 2 sqrt(w_i w_j)|cos|, at most 6w for equal slits, so kappa >= k^2 w^3 / (6w) = (k w)^2 / 6.
import json, itertools, math, numpy as np
rng=np.random.default_rng(12)
def I(mu,parts):
    n=len(parts); return sum((-1)**(n-r)*mu(frozenset().union(*S)) for r in range(1,n+1) for S in itertools.combinations(parts,r))
def law_mu(w,k): return lambda A: (np.prod([1+k*w[i] for i in A])-1)/k
P=[frozenset([i]) for i in range(3)]; res={}
# B1: reading A formula for kappa, checked from the measure itself
r=[]
for _ in range(100):
    w=rng.random(3)*0.4; k=rng.uniform(-2,2); mu=law_mu(w,k)
    d=sum(abs(I(mu,[P[i],P[j]])) for i,j in ((0,1),(1,2),(0,2)))
    r.append(abs(abs(I(mu,P))/d-abs(k)*np.prod(w)/(w[0]*w[1]+w[1]*w[2]+w[0]*w[2])))
res['B1']=max(r)<1e-12
# B2: equal slits give kappa = |k| w / 3
w=np.full(3,1/3); r=[abs(abs(I(law_mu(w,k),P))/sum(abs(I(law_mu(w,k),[P[i],P[j]])) for i,j in ((0,1),(1,2),(0,2)))-abs(k)/9) for k in (-1,-0.1,0.01,0.5)]
res['B2']=max(r)<1e-9   # small k loses digits in (prod - 1)/k, so a looser tolerance
# B3: reading B: the largest quantum pairwise sum over all phases is 6w for equal slits
best=0
for _ in range(20000):
    ph=rng.uniform(0,2*np.pi,3); psi=np.sqrt(1/3)*np.exp(1j*ph); mu=lambda A: abs(sum(psi[i] for i in A))**2
    best=max(best,sum(abs(I(mu,[P[i],P[j]])) for i,j in ((0,1),(1,2),(0,2))))
res['B3']=best<=6/3+1e-12 and best>6/3-1e-3
# B4: the caps (equal slits, w = 1/3, i.e. each slit carries a third of the flux)
caps={}
for name,kap in (('Sinha 2010',1e-2),('Kauten 2017',1e-4)):
    caps[name]={'kappa':kap,'A_kw':3*kap,'A_k':9*kap,'B_kw':math.sqrt(6*kap),'B_k':3*math.sqrt(6*kap)}
res['B4']=abs(caps['Kauten 2017']['A_k']-9e-4)<1e-15 and abs(caps['Kauten 2017']['B_k']-3*math.sqrt(6e-4))<1e-15
for k,v in res.items(): print(k,"PASS" if v else "FAIL")
print(f"{sum(res.values())}/{len(res)} passed")
for name,c in caps.items():
    print(f"{name}: kappa < {c['kappa']:.0e}.  law alone: |k w| < {c['A_kw']:.1e}, |k| < {c['A_k']:.1e} at w = 1/3.  "
          f"law on top of quantum: |k w| < {c['B_kw']:.3f}, |k| < {c['B_k']:.3f} at w = 1/3")
json.dump({'checks':{k:bool(v) for k,v in res.items()},'caps':caps},open('results.json','w'),indent=1)
