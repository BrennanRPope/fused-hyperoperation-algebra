# DRAFT. Round 6 check: information causality (IC) with a Renyi entropy, per tree depth n.
# Bit-guessing with bias E at each level: Bob's guess of each of the 2^n bits has bias E^n.
# One simple Renyi version (ours, one of several definitions): per bit I_a(e) = 1 - H_a((1+e)/2),
# H_a the binary Renyi entropy in bits (a=1 is Shannon). IC at depth n: 2^n * I_a(E^n) <= 1.
# Threshold E_a(n) = largest E obeying it. Quantum boxes reach E = 1/sqrt2 = 0.7071.
import json, math
from scipy.optimize import brentq
def H(a,p):
    q=1-p
    if abs(a-1)<1e-12:
        return -sum(x*math.log2(x) for x in (p,q) if x>0)
    if a==math.inf: return -math.log2(max(p,q))
    return math.log2(p**a+q**a)/(1-a)
def I(a,e): return 1-H(a,(1+e)/2)
def thr(a,n):
    f=lambda E: 2**n*I(a,E**n)-1
    return brentq(f,1e-9,1-1e-12)
alphas=[0.5,1,2,3,math.inf]; ns=list(range(1,13))
tab={str(a):[thr(a,n) for n in ns] for a in alphas}
res={}
# R1: Renyi-2 at n=1 equals sqrt(sqrt2-1)=0.6436 (matches Oughton-Timpson's 0.644)
r1=abs(tab['2'][0]-math.sqrt(math.sqrt(2)-1))<1e-9; res['R1']=r1
# R2: Shannon thresholds stay above 1/sqrt2 at every depth and approach it
s=tab['1']; r2=all(x>1/math.sqrt(2) for x in s) and s[-1]-1/math.sqrt(2)<0.02; res['R2']=r2
# R3: every alpha>1 dips below 1/sqrt2 at some depth (would forbid quantum boxes)
r3=all(min(tab[str(a)])<1/math.sqrt(2) for a in [2,3,math.inf]); res['R3']=r3
# R4: alpha=0.5 stays above 1/sqrt2 at every depth computed (consistent with QM, like Shannon)
r4=all(x>1/math.sqrt(2) for x in tab['0.5']); res['R4']=r4
# R5: the edge order is 2 ln2 = 1.386: small-bias I_a ~ a E^2/(2 ln2), so deep levels sit above
# 1/sqrt2 only when a < 2 ln2. Check 1.38 stays above at depths 1..39 and 1.40 dips below.
ok=lambda a: min(thr(a,n) for n in range(1,40))>1/math.sqrt(2)
res['R5']=ok(1.38) and not ok(1.40)
for a in alphas:
    print(f"alpha={a}: "+" ".join(f"{x:.4f}" for x in tab[str(a)]))
for k,v in res.items(): print(k,"PASS" if v else "FAIL")
print(f"{sum(res.values())}/{len(res)} passed")
json.dump({'n':ns,'thresholds':tab,'checks':res},open('results.json','w'),indent=1)
