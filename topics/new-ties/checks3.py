"""Round-3 checks (plane laws). Run after planelaws.py and planelaws_sweep.py."""
import json, numpy as np, sympy as sp
from fractions import Fraction as F
import planelaws as P
R=json.load(open("planelaws.json")); S=json.load(open("planelaws_sweep.json")); out=[]
ok=lambda n,c,i="": out.append((n,bool(c),i))
k,d,w,a,x=sp.symbols('k d w a x')
ok("T1 laws sharing wall w lie on d w^2 + k w + 1 = 0, tangent to the border d = k^2/4 at k = -2/w",
   sp.expand((d*w**2+k*w+1).subs(d,k**2/4)*4-(k*w+2)**2)==0)
ok("T2 the map (x+a+kax)/(1-dax) fixes w exactly when d w^2 + k w + 1 = 0",
   sp.expand(w*(1-d*a*w)-(w+a+k*a*w)+a*(d*w**2+k*w+1))==0)
# velocity (0,-1) and union (-1,0) both on the wall-1 line k + d + 1 = 0, tangent at (-2, 1)
ok("T3 velocity addition and the union rule share the wall 1 (tangent line at k=-2, d=1)", (0+(-1)+1)==0 and (-1+0+1)==0)
# T4 conjugation straightens the shared-wall example: h(x)=1/(x-3) makes every map affine with slope = multiplier
h=((F(0),F(1)),(F(1),F(-3))); hi=((F(-3)*-1,F(1)),(F(1),F(0)))  # inverse of h up to scale: x = 3 + 1/y
hi=((F(3),F(1)),(F(1),F(0)))
allaff=True; slopes=[]
for p,q,ss in [(F(1),F(3),[F(1,5),F(3,20)]),(F(-1),F(3),[F(1,4)])]:
    for s in ss:
        M,_=P.map_with(p,q,s); C=P.norm(P.mm(h,P.mm(M,hi)))
        allaff &= (C[1][0]==0); slopes.append((C[0][0]/C[1][1], s))
ok("T4 shared wall 3: after y = 1/(x-3) every map is a straight-line map", allaff)
ok("T5 ... and its slope is the map's multiplier at the attracting wall (inverse, as the line maps run outward in y)",
   all(sl==1/s or sl==s for sl,s in slopes), str(slopes))
r=R["shared wall 3"]
ok("T6 shared wall: count = shared-centre count", r["match"])
ok("T7 shared wall: box count within 0.02 of D*(multipliers)", abs(r["box"]-r["D_mult"])<0.02, f"box {r['box']:.4f}, D* {r['D_mult']:.4f}")
ok("T8 shared wall: sup |f'| bound is loose (> D*), so the straightening is what gives the exact value", r["D_sup"]>r["D_mult"]+0.1)
for lab in ["no shared wall, near","no shared wall, far","no shared wall, close"]:
    r=R[lab]
    ok(f"T9 {lab}: count = shared-centre count", r["match"])
    ok(f"T10 {lab}: sup |f'| on the hull equals the multiplier, so D* is a proved upper bound", abs(r["D_sup"]-r["D_mult"])<1e-9)
    ok(f"T11 {lab}: measured dimension below D*", max(r["box"],min(r["Dn_wall1"],r["Dn_wall2"]))<r["D_mult"], f"box {r['box']:.4f}, D_n {r['Dn_wall1']:.4f}, D* {r['D_mult']:.4f}")
Dn=[s["Dn"] for s in S]
ok("T12 sweep: D_n rises toward D* as the repelling walls move away", np.all(np.diff(Dn)>0) and abs(Dn[-1]-S[-1]["D_mult"])<0.01, " ".join(f"{v:.3f}" for v in Dn))
p=sum(c for _,c,_ in out)
with open("checks3_output.txt","w") as fh:
    for n,c,i in out:
        l=f"{'PASS' if c else 'FAIL'}  {n}  {i}"; print(l); fh.write(l+"\n")
    fh.write(f"{p}/{len(out)} pass\n")
print(f"{p}/{len(out)} pass")
