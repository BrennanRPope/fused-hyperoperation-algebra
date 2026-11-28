# DRAFT. Round 11: is the interference rule the law in disguise, and what does averaging do to it?
# The law as a measure on disjoint pieces: mu(A u B) = mu(A) + mu(B) + k mu(A) mu(B)  (Sugeno's lambda-measure).
# Quantum measure on paths: mu(A) = |sum of amplitudes in A|^2  (Sorkin's quantum measure).
# Interference orders (Sorkin): I1 = mu(A); I2(A,B) = mu(AB) - mu(A) - mu(B);
# I3(A,B,C) = mu(ABC) - mu(AB) - mu(BC) - mu(AC) + mu(A) + mu(B) + mu(C); I4 likewise.
import json, itertools, math, numpy as np
rng=np.random.default_rng(11)
def I(mu,parts):  # Sorkin interference term of order len(parts), by inclusion-exclusion
    n=len(parts); return sum((-1)**(n-r)*mu(frozenset().union(*S)) for r in range(1,n+1) for S in itertools.combinations(parts,r))
def law_mu(w,k):   # law measure from single-path weights w: lift 1 + k mu multiplies over disjoint pieces
    return lambda A: (np.prod([1+k*w[i] for i in A])-1)/k if abs(k)>1e-15 else sum(w[i] for i in A)
def q_mu(psi): return lambda A: abs(sum(psi[i] for i in A))**2
res={}
# Q1: two paths: |sqrt(p) e^{ia} + sqrt(q) e^{ib}|^2 = p + q + 2 sqrt(pq) cos(phi): the law of cosines
r=[]
for _ in range(200):
    p,q=rng.random(2); a,b=rng.uniform(0,2*np.pi,2)
    r.append(abs(abs(math.sqrt(p)*np.exp(1j*a)+math.sqrt(q)*np.exp(1j*b))**2-(p+q+2*math.sqrt(p*q)*math.cos(a-b))))
res['Q1']=max(r)<1e-12
# Q2: not the law at any fixed k: the k it would need, 2 cos(phi)/sqrt(pq), changes with p and q at a fixed phase
phi=math.pi/3; ks=[2*math.cos(phi)/math.sqrt(p*q) for p,q in rng.random((200,2))]
res['Q2']=np.std(ks)>1.0
# Q3: the law at k has interference at EVERY order: I_n = k^(n-1) * product of the pieces
r=[]
for k in (-1.0,-0.5,0.7,2.0):
    for _ in range(50):
        w=rng.random(4)*0.4; mu=law_mu(w,k); P=[frozenset([i]) for i in range(4)]
        for n in (2,3,4): r.append(abs(I(mu,P[:n])-k**(n-1)*np.prod(w[:n])))
res['Q3']=max(r)<1e-12
# Q4: quantum stops at order 2: I2 is generally nonzero, I3 = I4 = 0 exactly (Sorkin)
i2=[];i34=[]
for _ in range(200):
    psi=rng.normal(size=4)+1j*rng.normal(size=4); mu=q_mu(psi); P=[frozenset([i]) for i in range(4)]
    i2.append(abs(I(mu,P[:2]))); i34+= [abs(I(mu,P[:3])),abs(I(mu,P[:4]))]
res['Q4']=np.median(i2)>0.1 and max(i34)<1e-10
# Q5: averaging over phases (aggregation): coherence g scales I2 by g; full averaging (g = 0) leaves the k = 0 law
# Exact average over a uniform random phase on path 2 of 2 is p + q; partial dephasing multiplies the cross term by g.
p,q=0.3,0.5; gs=np.linspace(0,1,21); I2g=[]; I2mc=[]
for g in gs:
    I2g.append(2*math.sqrt(p*q)*g)   # expected: cross term times the coherence g
    if g<=0: th=rng.uniform(0,2*np.pi,200000)
    else: th=rng.normal(0,math.sqrt(-2*math.log(g)) if g<1 else 0.0,200000)   # wrapped Gaussian noise, E[e^{i th}] = g
    I2mc.append(float(np.mean(np.abs(math.sqrt(p)+math.sqrt(q)*np.exp(1j*th))**2))-(p+q))
th=rng.uniform(0,2*np.pi,400000); mc=np.mean(np.abs(math.sqrt(p)+math.sqrt(q)*np.exp(1j*th))**2)
res['Q5']=abs(mc-(p+q))<2e-3 and max(abs(np.array(I2mc)-np.array(I2g)))<5e-3
# Q6: so averaging quantum measures can never give the law at k != 0: any average keeps I3 = 0,
# but the law at k != 0 has I3 = k^2 * a*b*c != 0. Check: random averages of quantum measures have I3 = 0.
r=[]
for _ in range(50):
    psis=[rng.normal(size=3)+1j*rng.normal(size=3) for _ in range(20)]; wts=rng.dirichlet(np.ones(20))
    mu=lambda A: sum(wt*abs(sum(ps[i] for i in A))**2 for wt,ps in zip(wts,psis)); P=[frozenset([i]) for i in range(3)]
    r.append(abs(I(mu,P)))
w=np.array([0.2,0.3,0.25]); lawI3=abs(I(law_mu(w,0.5),[frozenset([i]) for i in range(3)]))
res['Q6']=max(r)<1e-10 and lawI3>1e-3
for k,v in res.items(): print(k,"PASS" if v else "FAIL")
print(f"{sum(res.values())}/{len(res)} passed")
print(f"Monte Carlo phase average {mc:.4f} vs p+q {p+q:.4f}; law k=0.5 gives I3 = {lawI3:.5f}")
# data for the picture: |I_n| for one example of each kind
w=np.array([0.3,0.25,0.2,0.15]); P=[frozenset([i]) for i in range(4)]
psi=np.sqrt(w)*np.exp(1j*np.array([0,1.0,2.2,4.0]))
bars={'classical (k = 0)':[abs(I(law_mu(w,0.0),P[:n])) for n in (1,2,3,4)],
      'quantum':[abs(I(q_mu(psi),P[:n])) for n in (1,2,3,4)],
      'the law at k = 2':[abs(I(law_mu(w,2.0),P[:n])) for n in (1,2,3,4)]}
json.dump({'checks':{k:bool(v) for k,v in res.items()},'bars':bars,'g':gs.tolist(),'I2g':I2g,'I2mc':I2mc},open('results.json','w'),indent=1)
