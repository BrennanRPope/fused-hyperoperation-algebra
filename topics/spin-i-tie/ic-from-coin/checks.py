# DRAFT. Round 9: does information causality (IC) follow from the coin rule plus no-signalling?
# Coin rule (round 8): an irrelevant coin flip doesn't change what's left to learn -> Shannon on records.
# Test: PR boxes obey no-signalling, records obey Shannon, so if IC still fails, those two aren't enough.
import json, itertools, math, numpy as np
def h(p): return 0.0 if p<=0 or p>=1 else -(p*math.log2(p)+(1-p)*math.log2(1-p))
def box(E):  # P(A,B|x,y): A xor B = x*y with probability (1+E)/2, uniform marginals
    return lambda A,B,x,y: ((1+E)/2 if (A^B)==(x&y) else (1-E)/2)/2
res={}
# N1: PR box (E = 1) is no-signalling: each side's marginal ignores the other side's input
P=box(1.0); ok=True
for x,y,A in itertools.product((0,1),repeat=3):
    m=[sum(P(A,B,x,yy) for B in (0,1)) for yy in (0,1)]; ok&=abs(m[0]-m[1])<1e-15
for x,y,B in itertools.product((0,1),repeat=3):
    m=[sum(P(A,B,xx,y) for A in (0,1)) for xx in (0,1)]; ok&=abs(m[0]-m[1])<1e-15
res['N1']=ok
# N2: the depth-1 IC game, computed exactly. Alice holds a0, a1, inputs x = a0 xor a1, sends m = a0 xor A (1 bit).
# Bob wants bit b, inputs y = b, guesses m xor B. IC says I(a0:guess|b=0) + I(a1:guess|b=1) <= 1 bit.
def ic_sum(E):
    P=box(E); tot=0
    for b in (0,1):
        J=np.zeros((2,2))
        for a0,a1,A,B in itertools.product((0,1),repeat=4):
            J[(a0,a1)[b], a0^A^B]+=0.25*P(A,B,a0^a1,b)
        pa=J.sum(1); pg=J.sum(0)
        tot+=sum(J[i,j]*math.log2(J[i,j]/(pa[i]*pg[j])) for i in (0,1) for j in (0,1) if J[i,j]>0)
    return tot
Es=np.linspace(0,1,101); S=[ic_sum(E) for E in Es]
res['N2']=abs(ic_sum(1.0)-2)<1e-12 and abs(ic_sum(1/math.sqrt(2))-2*(1-h((1+1/math.sqrt(2))/2)))<1e-12 and ic_sum(1/math.sqrt(2))<1
# N3: the records here are classical bits, so Shannon's coin rule holds for them (coin on everything and on one outcome)
rng=np.random.default_rng(9)
def H(p): p=np.asarray(p).ravel(); p=p[p>0]; return float(-(p*np.log2(p)).sum())
ok=True
for _ in range(50):
    J=rng.random((3,4)); J/=J.sum(); L=H(J)-H(J.sum(1)); c=np.array([.3,.7])
    G=np.einsum('ab,c->acb',J,c).reshape(-1,4); Lo=np.vstack([.4*J[:1],.6*J[:1],J[1:]])
    ok&=abs(H(G)-H(G.sum(1))-L)<1e-12 and abs(H(Lo)-H(Lo.sum(1))-L)<1e-12
res['N3']=ok
# So: no-signalling (N1) + coin rule on records (N3) still allow 2 bits through a 1-bit message (N2). Not enough.
# N4: the missing piece (Barnum et al. 2009): the theory's own states need ONE entropy for "how mixed" and
# "how unpredictable". A gbit (the local system of a PR box; state space = a square) fails this.
def meas_ent(x,y): return min(h((1+x)/2),h((1+y)/2))     # its two fine-grained measurements
def mix_ent(x,y):  # least Shannon entropy over mixtures of the 4 corners with mean (x,y)
    best=9
    for a in np.linspace(0,1,4001):
        b=(1+x)/2-a; c=(1+y)/2-a; d=1-a-b-c
        w=np.array([a,b,c,d])
        if (w>=-1e-12).all(): best=min(best,H(np.clip(w,0,1)))
    return best
res['N4']=abs(meas_ent(0.5,0)-h(0.75))<1e-12 and abs(mix_ent(0.5,0)-1.5)<1e-3
# N5: a qubit (state space = a ball) passes: both entropies equal von Neumann's for random states
ok=True
for _ in range(50):
    r=rng.random(); th=np.linspace(0,math.pi,721)
    me=min(h((1+r*math.cos(t))/2) for t in th); mx=h((1+r)/2)   # mixing: eigen-decomposition is optimal
    ok&=abs(me-mx)<1e-6
res['N5']=ok
for k,v in res.items(): print(k,"PASS" if v else "FAIL")
print(f"{sum(res.values())}/{len(res)} passed")
print(f"IC sum: classical E=0.5 -> {ic_sum(0.5):.4f}, quantum E=0.7071 -> {ic_sum(1/math.sqrt(2)):.4f}, PR E=1 -> {ic_sum(1.0):.4f}")
print(f"gbit at (0.5,0): measurement entropy {meas_ent(0.5,0):.4f} bits, mixing entropy {mix_ent(0.5,0):.4f} bits")
json.dump({'checks':{k:bool(v) for k,v in res.items()},'E':Es.tolist(),'ic_sum':S},open('results.json','w'),indent=1)
