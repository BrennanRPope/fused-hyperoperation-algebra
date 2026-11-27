"""Round 2: lower bound construction, generalizations, and a learning-rule scenario. Saves round2.json."""
from fractions import Fraction as F
import numpy as np, json, itertools, math
from scipy.optimize import brentq
import union_ifs as U
R={}

# ---------- R1 lower bound: first-return generators S = {u v : u in box-1 monoid \ 0, v in (other boxes' free product) \ 0}
def box_monoid_slopes(sl, L):
    # all nonzero multisets of the box's maps of total length <= L -> slope products (distinct monoid elements)
    out=[]
    m=len(sl)
    for exps in itertools.product(range(L+1),repeat=m):
        t=sum(exps)
        if 1<=t<=L: out.append(math.prod(float(s)**e for s,e in zip(sl,exps)))
    return out
def other_slopes(boxes, L):
    # distinct elements (traces) of the free product of the other boxes, length <= L, via exact composites
    maps=U.maps_from(boxes)
    level={(F(1),F(0))}; allm=set()
    for _ in range(L):
        level={(ms*s, ms*c+mc) for (s,c) in level for (_,ms,mc) in maps}; allm|=level
    return [float(s) for s,c in allm]
def sdim_S(boxes, L):
    A=box_monoid_slopes(boxes[0][1],L); B=other_slopes(boxes[1:],L)
    A=np.array(A); B=np.array(B)
    return brentq(lambda D: np.sum(A**D)*np.sum(B**D)-1,1e-6,5), len(A)*len(B)
def free_check(boxes, Lgen, depth):
    # generators of S_L as exact maps; check that distinct generator words of given depth give distinct maps
    b1=U.maps_from([boxes[0]]); bo=U.maps_from(boxes[1:])
    def elems(maps,L,free_product):
        level={(F(1),F(0))}; allm=set()
        for _ in range(L):
            level={(ms*s, ms*c+mc) for (s,c) in level for (_,ms,mc) in maps}; allm|=level
        return allm
    A=elems(b1,Lgen,False); B=elems(bo,Lgen,True)
    gens=[(sa*sb, sa*cb+ca) for (sa,ca) in A for (sb,cb) in B]   # u after v: x -> u(v(x))
    gens=list(set(gens))
    words=1; level={(F(1),F(0)):1}
    seen=[(F(1),F(0))]
    lvl=[(F(1),F(0))]
    for _ in range(depth):
        lvl=[(gs*s, gs*c+gc) for (s,c) in lvl for (gs,gc) in gens]
    return len(gens), len(lvl), len(set(lvl))
R["R1"]={}
for name in ["A: two in box 1, one in box 3","B: two boxes, two maps each","D: three in box 1, one in box 5"]:
    boxes=U.EXAMPLES[name]; Ds=U.trace_root(boxes)
    seq=[]
    for L in [1,2,3,4,6,8,10]:
        d,ng=sdim_S(boxes,L); seq.append((L,d,ng))
    g,nw,nd=free_check(boxes,2,3)
    R["R1"][name]=dict(D_trace=Ds,sdim_S=seq,free_check=dict(generators=g,words=nw,distinct=nd))
    print(name,"D*=%.5f"%Ds," s-dim(S_L):"," ".join("L%d:%.4f"%(L,d) for L,d,_ in seq)," free:",g,nw,nd)

# ---------- R2 any slopes, including a reflection (past the wall: slope < 0)
boxes_ref=[(F(1),[F(1,5),F(-1,6)]),(F(3),[F(1,4)])]
maps=U.maps_from(boxes_ref)
counts,levels=U.distinct_composites(maps,10); tc=U.trace_counts(boxes_ref,10)
f=lambda D: 1+sum(np.prod([1-abs(float(s))**D for s in sl])-1 for M,sl in boxes_ref)
Ds=brentq(f,1e-9,5); Dn=[U.pressure_D({(abs(s),c) for s,c in L}) for L in levels]
R["R2"]=dict(counts=counts,trace=tc,match=counts==tc,D_trace=Ds,Dn=Dn)
print("R2 reflection: match",counts==tc,"D*=%.4f D_10=%.4f"%(Ds,Dn[-1]))

# ---------- R3 the plane: complex similarities z -> c + lam (z - c); commute iff same centre
class G:  # Gaussian rationals as pairs of Fractions
    pass
def cm(a,b): return (a[0]*b[0]-a[1]*b[1], a[0]*b[1]+a[1]*b[0])
def ca(a,b): return (a[0]+b[0],a[1]+b[1])
def cs(a,b): return (a[0]-b[0],a[1]-b[1])
Z0=(F(0),F(0)); ONE=(F(1),F(0))
plane=[ # (centre, lam)
  ((F(0),F(0)), (F(1,4),F(1,4))),   # |lam| = 0.354, turns 45 deg
  ((F(0),F(0)), (F(3,10),F(0))),    # same centre, no turn
  ((F(1),F(0)), (F(0),F(1,3))),     # centre 1, turns 90 deg
  ((F(0),F(1)), (F(1,3),F(0))),     # centre i
]
pmaps=[(lam, cm(cs(ONE,lam),c)) for c,lam in plane]   # z -> lam z + (1-lam)c
lvl={(ONE,Z0)}; pc=[]
for n in range(9):
    lvl={(cm(l,s), ca(cm(l,t),b)) for (s,t) in lvl for (l,b) in pmaps}; pc.append(len(lvl))
plane_boxes=[(0,[0.25*2**0.5,0.3]),(1,[1/3]),(2,[1/3])]
ptc=U.trace_counts([(F(0),[F(1),F(2)]),(F(1),[F(3)]),(F(2),[F(5)])],9)  # only the pattern matters
fD=lambda D: 1+sum(np.prod([1-s**D for s in sl])-1 for _,sl in plane_boxes)
Dp=brentq(fD,1e-9,5); Dmor=brentq(lambda D: sum(s**D for _,sl in plane_boxes for s in sl)-1,1e-9,5)
# box counting in 2D: images of the three centres under all level-9 maps
cents=[(F(0),F(0)),(F(1),F(0)),(F(0),F(1))]
P=np.array([[float(x[0]),float(x[1])] for (s,t) in lvl for c in cents for x in [ca(cm(s,c),t)]])
lo=P.min(0); span=(P.max(0)-lo).max()
eps=[];Ns=[]
for j in range(3,12):
    e=span/2**j; eps.append(e); Ns.append(len({tuple(r) for r in np.floor((P-lo)/e).astype(int)}))
bs=np.polyfit(np.log(1/np.array(eps[3:8])),np.log(np.array(Ns[3:8])),1)[0]
R["R3"]=dict(counts=pc,trace=ptc,match=pc==ptc,D_trace=Dp,D_moran=Dmor,box=bs,points=P[::7].tolist())
print("R3 plane: match",pc==ptc,"D*=%.4f Moran=%.4f box=%.4f"%(Dp,Dmor,bs))

# ---------- R4 plane laws (k,d): Möbius maps x -> (x + a + k a x)/(1 - d a x). Same law commute (same walls).
def mob(k,d,a):  # matrix [[1+k a, a],[-d a, 1]]
    return ((1+k*a, a),(-d*a, F(1)))
def mm(A,B): return tuple(tuple(sum(A[i][t]*B[t][j] for t in range(2)) for j in range(2)) for i in range(2))
def norm(A):
    for v in (A[1][1],A[1][0],A[0][1],A[0][0]):
        if v!=0: return tuple(tuple(x/v for x in row) for row in A)
LAWSETS={"union x2 + law k=0,d=-4 (walls +-1/2)":[mob(F(-1),F(0),F(1,2)),mob(F(-1),F(0),F(2,3)),mob(F(0),F(-4),F(1,5))],
          "union x2 + velocity addition":[mob(F(-1),F(0),F(3,5)),mob(F(-1),F(0),F(5,7)),mob(F(0),F(-1),F(5,11))],
          "collision example: union x2 + velocity, slabs 1/2, 2/3, 1/2":[mob(F(-1),F(0),F(1,2)),mob(F(-1),F(0),F(2,3)),mob(F(0),F(-1),F(1,2))]}
mtc=U.trace_counts([(F(1),[F(1),F(2)]),(F(2),[F(3)])],10)
R["R4"]={}
for lab,laws in LAWSETS.items():
    lv={norm(((F(1),F(0)),(F(0),F(1))))}; mc=[]
    for n in range(10):
        lv={norm(mm(g,A)) for A in lv for g in laws}; mc.append(len(lv))
    first=next((i+1 for i,(a,b) in enumerate(zip(mc,mtc)) if a!=b),None)
    R["R4"][lab]=dict(counts=mc,trace=mtc,match=mc==mtc,first_collision_length=first)
    print("R4",lab,"match:",mc==mtc,"first collision at length",first)

# ---------- R5 scenario: Rescorla-Wagner / moving-average learning, V <- V + rate*(target - V)
# trial types: reward with salience 0.8, reward with salience 0.85 (target 1), no reward with rate 0.75 (target 0)
trials=[("reward, bright cue",F(4,5),F(1)),("reward, dim cue",F(17,20),F(1)),("no reward",F(3,4),F(0))]
rmaps=[(None,1-r,r*t) for _,r,t in trials]
n=11
dist={(F(1),F(0)):1}; Hs=[]; cnt=[]
for i in range(n):
    nd={}
    for (s,c),w in dist.items():
        for (_,ms,mc_) in rmaps:
            key=(ms*s, ms*c+mc_); nd[key]=nd.get(key,0)+w
    dist=nd; tot=3**(i+1); p=np.array(list(dist.values()))/tot
    Hs.append(float(-(p*np.log2(p)).sum())); cnt.append(len(dist))
rboxes=[(F(1),[F(1,5),F(3,20)]),(F(0),[F(1,4)])]
Dr=U.trace_root(rboxes); growth=1/brentq(lambda t: 1+(1-t)**2-1+(1-t)-1,1e-9,0.5)
lam=-np.mean([math.log(float(1-r)) for _,r,_ in trials])
# entropy rate of the random walk (uniform trial order) from increments of H_n
hrw_bits=Hs[-1]-Hs[-2]
R["R5"]=dict(counts=cnt,H_bits=Hs,growth=growth,bits_per_trial_max=math.log2(growth),
             bits_if_all_differ=math.log2(3),hrw_bits_increment=hrw_bits,D_set=Dr,
             D_measure_est=hrw_bits*math.log(2)/lam, lyapunov=lam)
print("R5 learning: distinct states",cnt[-3:],"growth/trial %.4f (%.3f bits), all-differ %.3f bits, H_n increment %.3f bits"%(growth,math.log2(growth),math.log2(3),hrw_bits))
print("   set dimension %.4f, measure dimension est %.4f"%(Dr,hrw_bits*math.log(2)/lam))
json.dump(R,open("round2.json","w"),indent=1,default=float)
