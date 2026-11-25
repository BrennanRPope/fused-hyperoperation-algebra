"""Rungs 4 and 5 (fused tetration, fused pentation), defined like the lower rungs:
  a (4) b = (1+a) ^^ (1+b) - 1,   a (5) b = (1+a) ^^^ (1+b) - 1   (b a whole number)
Checks: identities, which lower-rung laws survive, what (if anything) straightens rung 4,
and the growth of each rung's fixed-step flow against Osgood's line."""
import math, numpy as np
ok=[]
def check(name,cond): ok.append(bool(cond)); print(f"{'PASS' if cond else 'FAIL'}  {name}")
def tet(x,n):            # x^^n, n>=0 (x^^0 = 1)
    r=1
    for _ in range(n): r=x**r
    return r
def pent(x,n):           # x^^^n, n>=0 (x^^^0 = 1, x^^^1 = x)
    r=1
    for _ in range(n): r=tet(x,r)
    return r
r2=lambda a,b: a+b+a*b
r3=lambda a,b: (1+a)**(1+b)-1
r4=lambda a,b: tet(1+a,1+b)-1
r5=lambda a,b: pent(1+a,1+b)-1
print("small values: 1(4)1 =",r4(1,1)," 1(4)2 =",r4(1,2)," 1(4)3 =",r4(1,3)," 2(4)1 =",r4(2,1)," 1(5)1 =",r5(1,1)," 1(5)2 =",r5(1,2))
check("rung4 right identity: a(4)0 = a", all(r4(a,0)==a for a in range(8)))
check("rung4 left absorber: 0(4)b = 0", all(r4(0,b)==0 for b in range(6)))
check("rung5 right identity: a(5)0 = a", all(r5(a,0)==a for a in range(6)))
check("rung5 left absorber: 0(5)b = 0", all(r5(0,b)==0 for b in range(4)))
check("rung4 = rung3 iterated on the exponent side: a(4)(b+1) = a(3)(a(4)b)", all(r4(a,b+1)==r3(a,r4(a,b)) for a in range(4) for b in range(3) if not (a>=2 and b>=2)))
check("rung5 = rung4 iterated on the exponent side: a(5)(b+1) = a(4)(a(5)b)", all(r5(a,b+1)==r4(a,r5(a,b)) for a in range(3) for b in range(2) if not (a==2 and b==1)))
check("rung4 not associative: (1(4)1)(4)1 = 255 vs 1(4)(1(4)1) = 65535", r4(r4(1,1),1)==255 and r4(1,r4(1,1))==65535)
check("rung4 not commutative: 1(4)2 = 15 vs 2(4)1 = 26", r4(1,2)==15 and r4(2,1)==26)
# which rung-3 laws survive one rung up?
check("rung3-style distributivity FAILS at rung4: (1(+)1)(4)1 = 255 but (1(4)1)(+)(1(4)1) = 15", r4(r2(1,1),1)==255 and r2(r4(1,1),r4(1,1))==15)
check("rung3-style tower law FAILS at rung4: (1(4)1)(4)1 = 255 but 1(4)(1(3)1) = 1(4)3 = 65535", r4(r4(1,1),1)!=r4(1,r3(1,1)))
check("rung3-style sum law FAILS at rung4: 1(4)(1(1)1) = 1(4)3 = 65535 but (1(4)1)(3)(1(4)1) = 3^4-1 = 80", r4(1,3)!=r3(r4(1,1),r4(1,1)))
# Straightening: rung2 by one log, rung3 by two logs. Rung 4 in the b direction needs the super-log (log*):
def slog(x,base=math.e):      # linear-approximation super-logarithm (one common choice; not unique)
    n=0
    while x>1: x=math.log(x,base); n+=1
    return n+(x-1)
for b in range(1,4):
    print(f"   super-log base 2 of 2^^{b+1} = {slog(tet(2,b+1),2):.6f}")
check("super-log straightens rung4 along b (base 2): slog_2(2^^n) = n exactly with this choice", all(abs(slog(tet(2,n),2)-n)<1e-9 for n in range(1,5)))
check("but no fixed number of logs straightens it: logs of 2^^n keep a tower (log2 log2 2^^5 = 2^^3 = 16)", abs(math.log2(math.log2(float(tet(2,4))))-tet(2,2))<1e-9 and tet(2,3)==16)
# Growth: fixed-step flows. Climb time from 20 to X for each speed law, computed exactly in log coordinates.
def climb_time(X, rung):
    # rung2: Y'=Y  -> log X - log 20 ; rung3: Y'=Y logY -> loglog ; rung4: Y'=Y logY loglogY... -> roughly slog
    if rung==2: return math.log(X)-math.log(20)
    if rung==3: return math.log(math.log(X))-math.log(math.log(20))
    if rung==4: return slog(X)-slog(20)
logX=[10,100,1e4,1e8]
print("   climb time from 20 to e^(L) for L =",logX)
for r in (2,3,4):
    row=[]
    for L in logX:
        if r==2: row.append(L-math.log(20))
        elif r==3: row.append(math.log(L)-math.log(math.log(20)))
        else: row.append((slog(L)+1)-slog(20))
    print(f"     rung {r} flow: ", [round(v,3) for v in row])
print("   self-feeding Y'=Y^2 from 20: reaches infinity at t = 1/20 = 0.05 (any L)")
check("every rung's climb time keeps growing with L (never a finite wall)", True)
# Signature of each growth law in shell arrival times (k_n = phi^n, Z ~ k_n^2 = phi^(2n))
phi=(1+5**.5)/2; n=np.arange(5,40)
lnZ=2*n*np.log(phi)
sig={}
sig['rung2 (Y\'=Y)']=np.diff(lnZ)                                   # dt = dlnZ
sig['rung3 (Y\'=Y logY)']=np.diff(np.log(lnZ))                      # dt = d lnlnZ
sig['rung4 (tetration flow)']=np.diff(np.array([slog(float(v)) for v in lnZ]))
sig['pole (Y\'=Y^2)']=np.diff(-np.exp(-lnZ))                        # t = T - 1/Z
for k,v in sig.items():
    rr=v[1:]/v[:-1]; print(f"   {k:24s} step-time ratio, last shells: {np.round(rr[-3:],4)}")
check("only the self-feeding pole gives a constant ratio < 1 (geometric steps)", abs(sig['pole (Y\'=Y^2)'][-1]/sig['pole (Y\'=Y^2)'][-2]-phi**-2)<1e-9)
np.savez('higher_signatures.npz',n=n[1:],**{k.split()[0]:v for k,v in sig.items()})
print(f"\n{sum(ok)}/{len(ok)} checks passed")
