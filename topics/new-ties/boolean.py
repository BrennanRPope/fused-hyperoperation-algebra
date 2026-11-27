"""Round 5 (DRAFT): is there a family?  Union, intersect, symmetric-difference ("XOR") and disjoint-union fractals.
(README section 9.)  Reading: X is a random set of density x in a box of size 1, A an independent random set of
density a.  Any Boolean combination f(X, A) has density affine in x.  A fractal arises from repeating such steps."""
import itertools, json, sympy as sp, numpy as np
from fractions import Fraction as F
import union_ifs as U
x,a,k,M=sp.symbols('x a k M')

def classify():
    """all 16 Boolean functions f(X,A): density map x -> P(f), its slope, and its fixed point (centre)"""
    rows=[]
    for bits in itertools.product([0,1],repeat=4):          # f(0,0), f(0,1), f(1,0), f(1,1)
        f=dict(zip([(0,0),(0,1),(1,0),(1,1)],bits))
        dens=sum(f[(X,Y)]*(x if X else 1-x)*(a if Y else 1-a) for X in (0,1) for Y in (0,1))
        dens=sp.expand(dens); slope=sp.simplify(sp.diff(dens,x))
        if slope==0: rows.append((bits,str(dens),"constant (forgets x)",None)); continue
        cen=sp.simplify(sp.solve(sp.Eq(dens,x),x)[0]) if slope!=1 else None
        fixed=(cen is not None) and sp.diff(cen,a)==0
        rows.append((bits,str(dens),str(slope),str(cen) if cen is not None else "none (slope 1)",fixed))
    return rows

NAMES={(0,0,0,1):"X AND A (intersect)",(0,1,1,1):"X OR A (union)",(0,1,1,0):"X XOR A (symmetric difference)",
       (0,0,1,0):"X AND NOT A (difference)",(1,0,1,1):"X OR NOT A",(1,0,0,1):"X XNOR A",
       (0,0,1,1):"X (do nothing)",(1,1,0,0):"NOT X (complement)",(0,1,0,0):"NOT X AND A",(1,1,0,1):"NOT X OR A",
       (1,1,1,0):"NOT (X AND A)",(1,0,0,0):"NOT (X OR A)"}

if __name__=="__main__":
    out={}
    print("Boolean steps with a centre that does not depend on a:")
    rows=classify(); out["classification"]=[]
    for r in rows:
        if len(r)==5:
            bits,dens,slope,cen,fixed=r; name=NAMES.get(bits,str(bits))
            out["classification"].append(dict(name=name,density=dens,slope=slope,centre=cen,fixed_centre=bool(fixed)))
            print(f"  {name:32s} density {dens:22s} slope {slope:10s} centre {cen:16s} {'FIXED' if fixed else ''}")
    # fused law: x o_k a = a + x + k a x : centre -1/k, slope 1 + k a
    law=a+x+k*a*x
    out["fused_centre"]=str(sp.solve(sp.Eq(law,x),x)); out["fused_slope"]=str(sp.diff(law,x))
    # union in a box of size M is k=-1/M; XOR in box 1 equals union in box 1/2
    union_box=a+x-a*x/M
    out["xor_is_union_in_half_box"]=bool(sp.simplify(union_box.subs(M,sp.Rational(1,2))-(a+x-2*a*x))==0)
    # De Morgan: AND with a = complement o OR with (1-a) o complement
    out["and_is_conjugate_of_or"]=bool(sp.simplify(1-((1-a)+(1-x)-(1-a)*(1-x))-a*x)==0)
    # every fused law is multiplication in the coordinate u = 1 + k x
    out["fused_is_multiplication"]=bool(sp.simplify((1+k*law)-(1+k*a)*(1+k*x))==0)
    # k = 0 is plain adding: a translation, no fixed point, no attractor
    out["k0_translation"]=str(sp.diff(law.subs(k,0),x))
    print("XOR in box 1 = union in box 1/2:",out["xor_is_union_in_half_box"],
          "| AND = De Morgan conjugate of OR:",out["and_is_conjugate_of_or"],
          "| every fused law is multiplication in u = 1 + kx:",out["fused_is_multiplication"])
    # fair coins: AND 1/2 and OR 1/2 give binary expansion; the attractor is all of [0,1]
    fair=[(F(1),[F(1,2)]),(F(0),[F(1,2)])]
    out["fair_coins_Dstar"]=U.trace_root(fair)
    # a three-centre Boolean fractal: OR with 4/5 and 17/20 (centre 1), AND with 3/13 (centre 0),
    # XOR with 9/22 and 15/26 (centre 1/2, slopes 2/11 and -2/13: the second one flips).
    # Simple unit-fraction slopes collide by accident (e.g. OR 1/5 with XOR 1/6, 1/7 at length 4: r s1 s2 - r s2 + s1 - s2 = 0)
    boxes=[(F(1),[F(1,5),F(3,20)]),(F(0),[F(3,13)]),(F(1,2),[F(2,11),F(-2,13)])]
    maps=U.maps_from(boxes); n=9
    counts,levels=U.distinct_composites(maps,n)
    tc=U.trace_counts(boxes,n)
    absb=[(c,[abs(s) for s in sl]) for c,sl in boxes]
    Ds=U.trace_root(absb); Dm=U.moran_root(absb)
    pts=U.attractor_points(maps,levels[-1]); eps,Ns=U.boxdim(pts,pts.min(),pts.max(),range(4,17))
    sl=np.polyfit(np.log(1/np.array(eps))[4:],np.log(np.array(Ns))[4:],1)[0]
    out["boolean_fractal"]=dict(boxes=[(str(c),[str(s) for s in s_]) for c,s_ in boxes],counts=counts,trace_counts=tc,
        match=counts==tc,D_star=Ds,D_moran=Dm,box=float(sl),words=len(maps)**n)
    print(f"three-centre Boolean fractal: counts {counts[-3:]} vs trace {tc[-3:]} match {counts==tc}; "
          f"{len(maps)**n} words; D* {Ds:.4f}, Moran {Dm:.4f}, box {sl:.4f}")
    # the same without XOR is round 1's example A moved into one box (centres 0 and 1)
    json.dump(out,open("boolean.json","w"),indent=1)
