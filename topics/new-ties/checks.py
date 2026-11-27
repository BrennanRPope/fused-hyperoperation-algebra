"""Checks for the new-ties draft. Run: python3 checks.py"""
import json, numpy as np, sympy as sp
from fractions import Fraction as F
import union_ifs as U
out=[]; ok=lambda name,cond,info="": out.append((name,bool(cond),info))

# U1 union map: x o_k a with k=-1/M is the union of two independent slabs in a box M; fixed point M
x,a,M=sp.symbols('x a M',positive=True); k=-1/M
f=a+x+k*a*x
ok("U1 union size: a+x-ax/M = M*(1-(1-a/M)(1-x/M))", sp.simplify(f-M*(1-(1-a/M)*(1-x/M)))==0)
ok("U2 centre of every union map in box M is M (full box)", sp.simplify(f.subs(x,M)-M)==0)
# U3 commute iff same box (slopes != 1)
a1,a2,M1,M2=sp.symbols('a1 a2 M1 M2',positive=True)
g=lambda y,aa,MM: aa+y-aa*y/MM
gap=sp.simplify(g(g(x,a2,M2),a1,M1)-g(g(x,a1,M1),a2,M2))
ok("U3 order gap = a1*a2*(1/M2-1/M1)... zero iff M1=M2", sp.simplify(gap - a1*a2*(M1-M2)/(M1*M2))==0 or sp.simplify(gap + a1*a2*(M1-M2)/(M1*M2))==0, str(sp.factor(gap)))
# U4-U7 exact counts of distinct composites match Cartier-Foata trace counts (exact rational arithmetic)
res=json.load(open("union_ifs.json"))
for i,(name,r) in enumerate(res.items()):
    ok(f"U{4+i} distinct composites = trace counts, {name}", r["match"], f"last {r['counts'][-1]} of {r['words']} words")
# U8 Moran equation is the special case of one map per box
for name,r in res.items():
    if name.startswith("C"): ok("U8 one map per box: trace root = Moran root", abs(r["D_trace"]-r["D_moran"])<1e-9)
# U9 upper bound: D_n (pressure over distinct level-n maps) decreases and stays above D*
for name,r in res.items():
    Dn=np.array(r["Dn"]); ok(f"U9 D_n decreasing to D* from above, {name}", np.all(np.diff(Dn)<=1e-12) and Dn[-1]>=r["D_trace"]-1e-9, f"D_11={Dn[-1]:.4f}, D*={r['D_trace']:.4f}")
# U10 extrapolate D_n ~ D* + c/n
for name,r in res.items():
    Dn=np.array(r["Dn"]); n=np.arange(1,len(Dn)+1); A=np.vstack([np.ones(5),1/n[-5:]]).T
    c=np.linalg.lstsq(A,Dn[-5:],rcond=None)[0]
    ok(f"U10 1/n extrapolation of D_n lands on D* (within 0.005), {name}", abs(c[0]-r["D_trace"])<0.005, f"extrap {c[0]:.4f} vs D* {r['D_trace']:.4f}, Moran {r['D_moran']:.4f}")
# U11 box counting is nearer D* than naive Moran where they differ
for name,r in res.items():
    if name.startswith("C"): continue
    eps=np.array(r["eps"]); N=np.array(r["N"]); sl=np.polyfit(np.log(1/eps[6:12]),np.log(N[6:12]),1)[0]
    ok(f"U11 box count nearer D* than Moran, {name}", abs(sl-r["D_trace"])<abs(sl-r["D_moran"]), f"box {sl:.4f}, D* {r['D_trace']:.4f}, Moran {r['D_moran']:.4f}")

# C1-C4 Carnot parabola: discriminant 4 means every heat-type engine closes its gap at Carnot's rate e^{-beta}
b,g_=sp.symbols('beta gamma',positive=True)
etaC=1-sp.exp(-b); eta=etaC/(2-g_*etaC)
kk=2*g_-2; dd=-g_*(2-g_)
ok("C1 discriminant k^2-4d = 4 for every gamma (curvature reading -Delta/4 = -1)", sp.simplify(kk**2-4*dd-4)==0)
rhs=sp.Rational(1,2)*(1+kk*eta+dd*eta**2)
ok("C2 Riccati: d eta/d beta = (1/2)(1+k eta+d eta^2)", sp.simplify(sp.diff(eta,b)-rhs)==0)
wall=1/(2-g_)
ratio=sp.simplify((wall-eta)*sp.exp(b))
lim=sp.limit(ratio,b,sp.oo)
ok("C3 gap to the wall ~ C*e^{-beta}: (wall-eta)e^beta -> 2/(2-gamma)^2", sp.simplify(lim-2/(2-g_)**2)==0, str(lim))
ca=1-sp.exp(-b/2)
ok("C4 Curzon-Ahlborn closes at half rate: (1-eta_CA)e^{beta/2} = 1", sp.simplify((1-ca)*sp.exp(b/2)-1)==0)
# C5 general fact: rate at a wall of (1/2)P(eta) is (1/2)sqrt(Delta)
e_=sp.symbols('e'); K,D=sp.symbols('K D',real=True)
P=1+K*e_+D*e_**2; r1,r2=sp.solve(P,e_)
ok("C5 linear rate at a wall = (1/2)|P'(wall)| = (1/2)sqrt(Delta)", sp.simplify((sp.Rational(1,2)*sp.diff(P,e_).subs(e_,r1))**2-(K**2-4*D)/4)==0)

passed=sum(c for _,c,_ in out)
with open("checks_output.txt","w") as fh:
    for n,c,i in out:
        line=f"{'PASS' if c else 'FAIL'}  {n}  {i}"; print(line); fh.write(line+"\n")
    fh.write(f"{passed}/{len(out)} pass\n")
print(f"{passed}/{len(out)} pass")
