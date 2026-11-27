"""Round 3: union fractal for the plane laws (k, d). Maps x -> (x + a + k a x)/(1 - d a x), matrix [[1+ka, a],[-da, 1]].
Walls of law (k,d): roots of d x^2 + k x + 1 = 0  <=>  k = -(p+q)/(pq), d = 1/(pq).
Laws sharing a wall w lie on the straight line d w^2 + k w + 1 = 0 in the (k,d) plane."""
from fractions import Fraction as F
import numpy as np, json, math
from scipy.optimize import brentq
import union_ifs as U

def law(p,q): return (-(p+q)/(p*q), 1/(p*q))
def mat(k,d,a): return ((1+k*a, a),(-d*a, F(1)))
def mm(A,B): return tuple(tuple(sum(A[i][t]*B[t][j] for t in range(2)) for j in range(2)) for i in range(2))
def norm(A):
    for v in (A[1][1],A[1][0],A[0][1],A[0][0]):
        if v!=0: return tuple(tuple(x/v for x in row) for row in A)
def app(A,x): return (A[0][0]*x+A[0][1])/(A[1][0]*x+A[1][1])
def deriv(A,x):
    a,b=A[0]; c,d=A[1]; return abs(float(a*d-b*c))/float(c*x+d)**2
def map_with(p,q,s):
    """the map of law (p,q) whose multiplier at the attracting wall p is s (s rational): find a exactly.
    For Möbius fixing p,q the multiplier at p is lambda; in coordinate y=(x-p)/(x-q): y -> lambda y.
    Build matrix = H^-1 diag(lambda,1) H and check it has the law's form."""
    H=((F(1),-p),(F(1),-q)); Hi=((-q,p),(F(-1),F(1)))  # inverse up to scale
    M=norm(mm(Hi,mm(((s,F(0)),(F(0),F(1))),H)))
    k,d=law(p,q); a=M[0][1]
    assert M==norm(mat(k,d,a)), "not in the law's family"
    return M,a
def run(classes,n=10,name=""):
    # classes: list of (attracting wall p, repelling wall q, [multipliers s])
    maps=[]; cls=[]
    for ci,(p,q,ss) in enumerate(classes):
        for s in ss:
            M,a=map_with(p,q,s); maps.append(M); cls.append(ci)
    ps=[float(p) for p,q,ss in classes]; lo,hi=min(ps),max(ps)
    lv={norm(((F(1),F(0)),(F(0),F(1))))}; counts=[]
    for i in range(n):
        lv={norm(mm(g,A)) for A in lv for g in maps}; counts.append(len(lv))
    tc=U.trace_counts([(F(i),[F(j) for j in range(len(ss))]) for i,(p,q,ss) in enumerate(classes)],n)
    # D* with multipliers; upper bound with sup |f'| on the hull [lo,hi]
    mult=[[float(s) for s in ss] for p,q,ss in classes]
    xs=np.linspace(lo,hi,2001)
    sup=[];inf=[]
    for M in maps:
        dv=[deriv(M,F(x).limit_denominator(10**6)) for x in xs[::50]]+[deriv(M,F(lo)),deriv(M,F(hi))]
        sup.append(max(dv)); inf.append(min(dv))
    def Dstar(groups): return brentq(lambda D: 1+sum(np.prod([1-s**D for s in g])-1 for g in groups),1e-9,5)
    def group(vals):
        out=[[] for _ in classes]
        for v,c in zip(vals,cls): out[c].append(v)
        return out
    D_mult=Dstar(mult); D_sup=Dstar(group(sup)); D_inf=Dstar(group(inf))
    # Bowen-type estimates from the distinct level-n maps: derivative at an attractor point (midpoint of walls' images)
    x0=F(lo+hi)/2 if False else F(ps[0]).limit_denominator()  # a point of the attractor: an attracting wall
    def Dn_at(level,pt):
        g=np.array([deriv(A,pt) for A in level]); return brentq(lambda D: np.sum(g**D)-1,1e-9,5)
    Dn=Dn_at(lv,F(classes[0][0])); Dn2=Dn_at(lv,F(classes[-1][0]))
    # box counting: images of the attracting walls
    pts=np.array([float(app(A,p)) for A in lv for p,q,ss in classes])
    a0,b0=pts.min(),pts.max(); eps=[];N=[]
    for j in range(4,18):
        e=(b0-a0)/2**j; eps.append(e); N.append(len(np.unique(np.floor((pts-a0)/e))))
    box=np.polyfit(np.log(1/np.array(eps[6:12])),np.log(np.array(N[6:12])),1)[0]
    r=dict(counts=counts,trace=tc,match=counts==tc,D_mult=D_mult,D_sup=D_sup,D_inf=D_inf,Dn_wall1=Dn,Dn_wall2=Dn2,box=box,
           laws=[[float(x) for x in law(p,q)] for p,q,ss in classes],points=pts[::5].tolist())
    print(f"{name}: count match {counts==tc}; D*(multipliers)={D_mult:.4f}, upper D*(sup|f'|)={D_sup:.4f}, D*(inf)={D_inf:.4f}, D_n={Dn:.4f}/{Dn2:.4f}, box={box:.4f}")
    return r

R={}
s1,s2,s3=F(1,5),F(3,20),F(1,4)
R["shared wall 3"]=run([(F(1),F(3),[s1,s2]),(F(-1),F(3),[s3])],name="shared repelling wall 3")
R["shared wall inf (union fractal A)"]=None
R["no shared wall, near"]=run([(F(1),F(3),[s1,s2]),(F(-1),F(-3),[s3])],name="no shared wall (3 and -3)")
R["no shared wall, far"]=run([(F(1),F(12),[s1,s2]),(F(-1),F(-12),[s3])],name="no shared wall (12 and -12)")
R["no shared wall, close"]=run([(F(1),F(2),[s1,s2]),(F(-1),F(-2),[s3])],name="no shared wall (2 and -2)")
R.pop("shared wall inf (union fractal A)")
json.dump(R,open("planelaws.json","w"),indent=1)
