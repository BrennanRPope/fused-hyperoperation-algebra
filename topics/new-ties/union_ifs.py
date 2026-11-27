"""Union fractal: maps x -> x o_k a = a + (1+k a) x, k = -1/M (slab of size a joined in a box of size M).
Slope s = 1 - a/M, fixed point (centre) = M, the full box. Maps in one box commute; maps in different boxes don't.
Predicted dimension D*: 1 + sum_boxes [prod_{i in box}(1 - s_i^D) - 1] = 0  (Cartier-Foata trace count, weighted).
Naive Moran (ignores commuting): sum s_i^D = 1.
"""
from fractions import Fraction as F
import numpy as np, json, sys
from scipy.optimize import brentq

def maps_from(boxes):
    # boxes: list of (M, [slopes]) with exact Fractions
    out=[]
    for bi,(M,sl) in enumerate(boxes):
        for s in sl:
            out.append((bi, s, M*(1-s)))   # x -> s x + c, c = a = M(1-s)
    return out

def trace_root(boxes):
    f=lambda D: 1+sum(np.prod([1-float(s)**D for s in sl])-1 for M,sl in boxes)
    return brentq(f,1e-9,5)
def moran_root(boxes):
    f=lambda D: sum(float(s)**D for M,sl in boxes for s in sl)-1
    return brentq(f,1e-9,5)

def distinct_composites(maps, n):
    # exact: map = (slope, intercept) as Fractions. compose w then m: x -> m(w(x))
    level={(F(1),F(0))}
    counts=[]; levels=[]
    for _ in range(n):
        new=set()
        for (s,c) in level:
            for (_,ms,mc) in maps:
                new.add((ms*s, ms*c+mc))
        level=new; counts.append(len(level)); levels.append(level)
    return counts, levels

def trace_counts(boxes,n):
    # coefficients of 1/mu(t), mu = 1 + sum[(1-t)^m - 1]
    import sympy as sp
    t=sp.symbols('t'); mu=1+sum((1-t)**len(sl)-1 for M,sl in boxes)
    ser=sp.series(1/mu,t,0,n+1).removeO()
    return [int(ser.coeff(t,i)) for i in range(1,n+1)]

def pressure_D(level):
    sl=np.array([float(s) for s,c in level]); n=None
    return brentq(lambda D: np.sum(sl**D)-1,1e-9,5)

def attractor_points(maps,level):
    # images of every box centre under all composites of the last level
    cents=sorted(set(float(m[2]/(1-m[1])) for m in maps))
    pts=np.array([float(s)*x+float(c) for (s,c) in level for x in cents])
    return pts

def boxdim(pts, lo, hi, js):
    Ns=[]; eps=[]
    for j in js:
        e=(hi-lo)/2**j
        Ns.append(len(np.unique(np.floor((pts-lo)/e)))); eps.append(e)
    return np.array(eps),np.array(Ns)

EXAMPLES={
 "A: two in box 1, one in box 3": [(F(1),[F(1,5),F(3,20)]),(F(3),[F(1,4)])],
 "B: two boxes, two maps each": [(F(1),[F(1,5),F(1,7)]),(F(4),[F(1,6),F(1,8)])],
 "C: three boxes, one map each": [(F(1),[F(1,5)]),(F(2),[F(1,6)]),(F(4),[F(1,7)])],
 "D: three in box 1, one in box 5": [(F(1),[F(1,5),F(1,7),F(1,9)]),(F(5),[F(1,4)])],
}
if __name__=="__main__":
    res={}
    for name,boxes in EXAMPLES.items():
        maps=maps_from(boxes); n=int(sys.argv[1]) if len(sys.argv)>1 else 11
        counts,levels=distinct_composites(maps,n)
        tc=trace_counts(boxes,n)
        Dn=[pressure_D(L) for L in levels]
        Ds=trace_root(boxes); Dm=moran_root(boxes)
        pts=attractor_points(maps,levels[-1])
        lo,hi=pts.min(),pts.max()
        eps,Ns=boxdim(pts,lo,hi,range(4,18))
        res[name]=dict(counts=counts,trace_counts=tc,match=counts==tc,Dn=Dn,D_trace=Ds,D_moran=Dm,
                       eps=eps.tolist(),N=Ns.tolist(),words=len(maps)**n)
        print(name, "counts match trace monoid:",counts==tc, counts[-3:], "words:",len(maps)**n)
        print("   D_trace=%.5f  D_moran=%.5f  pressure D_n: %s"%(Ds,Dm," ".join("%.4f"%d for d in Dn[-4:])))
        sl=np.polyfit(np.log(1/eps[6:12]),np.log(Ns[6:12]),1)[0]
        print("   box-count slope (mid scales): %.4f"%sl)
    json.dump(res,open("union_ifs.json","w"),indent=1)
