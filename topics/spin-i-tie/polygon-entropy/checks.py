# DRAFT. Round 10: does "one entropy for mixing and measuring" (monoentropy) round a polygon into a disk?
# Regular n-gon state space (inscribed circle radius 1; n = 4 is the gbit square [-1,1]^2 of round 9).
# Mixing entropy: least Shannon entropy over ways to write the state as a mix of corners (pure states).
# Measurement entropy: least Shannon entropy over fine-grained measurements (effects on the edge rays).
# Both are least values of a concave function over a polytope, so a vertex wins: at most 3 nonzero terms.
# We enumerate every set of up to 3 corners / edge effects, so both entropies are exact.
import json, itertools, math, numpy as np
def verts(n):
    R=1/math.cos(math.pi/n); t=math.pi/n+2*math.pi*np.arange(n)/n
    return np.c_[R*np.cos(t),R*np.sin(t)]
def facets(n):  # affine functional a.x + c, zero on edge i, max 1 over the polygon (rows: ax, ay, c)
    V=verts(n); F=[]
    for i in range(n):
        p,q=V[i],V[(i+1)%n]; nrm=np.array([q[1]-p[1],p[0]-q[0]]); c=-nrm@p
        if c<0: nrm,c=-nrm,-c
        vals=V@nrm+c; F.append(np.r_[nrm,c]/vals.max())
    return np.array(F)
def Hrows(P):
    P=np.clip(P,0,1); with_=np.where(P>1e-15,P*np.log2(np.where(P>1e-15,P,1)),0); return -with_.sum(1)
def meas_ent(n,X):  # X: states (m,2)
    F=facets(n); Xh=np.c_[X,np.ones(len(X))]; best=np.full(len(X),np.inf); u=np.array([0,0,1.0])
    for r in (1,2,3):
        for S in itertools.combinations(range(n),r):
            A=F[list(S)].T
            c,res,rk,_=np.linalg.lstsq(A,u,rcond=None)
            if np.abs(A@c-u).max()>1e-9 or (c<-1e-12).any(): continue
            P=(Xh@A)*c; best=np.minimum(best,Hrows(P))
    return best
def mix_ent(n,X):
    V=verts(n); best=np.full(len(X),np.inf)
    for r in (1,2,3):
        for S in itertools.combinations(range(n),r):
            M=np.vstack([V[list(S)].T,np.ones(r)])
            W=np.linalg.lstsq(M,np.c_[X,np.ones(len(X))].T,rcond=None)[0].T
            ok=(np.abs(W@M.T-np.c_[X,np.ones(len(X))]).max(1)<1e-9)&(W>-1e-12).all(1)
            h=Hrows(W); best=np.where(ok,np.minimum(best,h),best)
    return best
def inside(n,X):
    F=facets(n); return ((np.c_[X,np.ones(len(X))]@F.T)>=-1e-12).all(1)
g=np.linspace(-1.5,1.5,121); GX=np.array([(x,y) for y in g for x in g])
ns=list(range(3,31)); out={}; res={}
for n in ns:
    X=GX[inside(n,GX)]; gap=mix_ent(n,X)-meas_ent(n,X)
    V=verts(n); mid=(V+np.roll(V,-1,0))/2; gm=mix_ent(n,mid)-meas_ent(n,mid)
    pt=np.array([[0.5,0.2]]); gp=float(abs(mix_ent(n,pt)-meas_ent(n,pt))[0])
    out[n]=dict(max=float(max(abs(gap).max(),abs(gm).max())),mean=float(abs(gap).mean()),frac=float((abs(gap)>0.1).mean()),mid=float(gm.min()),pt=gp)
    print(f"n={n:2d}  largest |gap| {out[n]['max']:.4f}  edge-midpoint gap {out[n]['mid']:.4f}  average |gap| {out[n]['mean']:.4f}  share |gap|>0.1 bit {out[n]['frac']:.3f}  |gap| at (0.5,0.2) {out[n]['pt']:.4f}")
res['P1']=out[3]['max']<1e-9                                   # triangle (a classical trit) is monoentropic
X=np.array([[0.5,0.0]]); res['P2']=abs(meas_ent(4,X)[0]-0.811278)<1e-5 and abs(mix_ent(4,X)[0]-1.5)<1e-9  # matches round 9
# P3: the worst disagreement never shrinks: even n give exactly 1 bit at edge midpoints; odd n climb toward 1 bit
res['P3']=all(abs(out[n]['mid']-1)<1e-9 for n in ns if n%2==0) and all(np.diff([out[n]['mid'] for n in ns if n%2==1 and n>=5])>0) and out[29]['mid']>0.94
# P4: on average the disagreement fades as sides are added (both parities), so the bulk rounds off
res['P4']=out[30]['mean']<0.2*out[4]['mean'] and out[29]['mean']<0.75*out[5]['mean']
# P5: same trend further out, on a coarser grid: 31 -> 61 sides (odd) and 30 -> 62 sides (even)
gc=np.linspace(-1.2,1.2,41); GC=np.array([(x,y) for y in gc for x in gc]); far={}
for n in (30,31,61,62):
    X=GC[inside(n,GC)]; far[n]=float(np.abs(mix_ent(n,X)-meas_ent(n,X)).mean())
    print(f"coarse grid n={n}: average |gap| {far[n]:.4f}")
res['P5']=far[62]<far[30] and far[61]<far[31]
for k,v in res.items(): print(k,"PASS" if v else "FAIL")
print(f"{sum(res.values())}/{len(res)} passed")
json.dump({'checks':{k:bool(v) for k,v in res.items()},'per_n':out,'coarse_far':far},open('results.json','w'),indent=1)
