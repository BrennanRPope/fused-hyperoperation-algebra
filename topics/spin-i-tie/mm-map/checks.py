# DRAFT. Round 13: the law against Masanes and Mueller's five requirements
# (Ll. Masanes, M. P. Mueller, "A derivation of quantum theory from physical requirements", NJP 13, 063001 (2011), arXiv 1004.1483).
# R1 Finiteness, R2 Local tomography, R3 Equivalence of subspaces, R4 Symmetry (any pure state reversibly to any other),
# R5 All measurements allowed. Their starting frame already assumes convex mixing with affine probabilities.
# Counting fact used (Hardy 2001, quant-ph/0101012): a system needs K fiducial numbers; classical K = d, complex quantum K = d^2,
# real quantum K = d(d+1)/2, quaternionic K = d(2d-1). Local tomography <=> K_AB = K_A * K_B.
import json, itertools
law=lambda a,b,k: a+b+k*a*b
fams={'classical':lambda d:d,'complex quantum':lambda d:d*d,'real quantum':lambda d:d*(d+1)//2,'quaternionic':lambda d:d*(2*d-1)}
res={}; table={}
# M1: with N = K - 1 (free numbers once the total is fixed), local tomography says N_AB = N_A o N_B at k = +1
ok=True
for name in ('classical','complex quantum'):
    K=fams[name]
    for da,db in itertools.product(range(2,6),repeat=2):
        ok&= K(da*db)-1 == law(K(da)-1,K(db)-1,1)
res['M1']=ok
# M2: real and quaternionic counts break the k = +1 law (real two-bit systems need more numbers than the law gives,
#     the extra ones only joint measurements see; the quaternionic count falls short, and that composite is not well defined anyway)
for name,K in fams.items():
    table[name]={'N_A':K(2)-1,'N_AB actual':K(4)-1,'N_AB by law k=1':law(K(2)-1,K(2)-1,1)}
res['M2']=table['real quantum']['N_AB actual']>table['real quantum']['N_AB by law k=1'] and \
          table['quaternionic']['N_AB actual']<table['quaternionic']['N_AB by law k=1']
# M3: the k = +1 law holds for every K = d^r, so it does not pick r = 2 over r = 1
res['M3']=all((da*db)**r-1==law(da**r-1,db**r-1,1) for r in range(1,5) for da in range(2,5) for db in range(2,5))
# M4: the mixing frame. Tsallis leftover after an independent fair coin is (1 + k S_coin) S_B, so it equals S_B only at k = 0
import math
def Sq(p,q): return -sum(x*math.log(x) for x in p if x>0) if abs(q-1)<1e-12 else (1-sum(x**q for x in p))/(q-1)
ok=True
for q in (0.5,1.0,2.0,3.0):
    k=1-q; c=[.5,.5]; B=[.2,.3,.5]; AB=[x*y for x in c for y in B]
    left=Sq(AB,q)-Sq(c,q); ok&= abs(left-(1+k*Sq(c,q))*Sq(B,q))<1e-12; ok&= (abs(left-Sq(B,q))<1e-12)==(q==1.0)
res['M4']=ok
for k,v in res.items(): print(k,"PASS" if v else "FAIL")
print(f"{sum(res.values())}/{len(res)} passed")
for n,t in table.items(): print(f"{n}: one bit N = {t['N_A']}, two bits N = {t['N_AB actual']}, k=1 law says {t['N_AB by law k=1']}")
json.dump({'passed':res,'two_bit_counts':table},open('results.json','w'),indent=1)
