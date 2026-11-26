"""DRAFT polish step: put each line-up exactly on the equator (bisection on time) before testing the swap rule."""
import sys,json,numpy as np
sys.path.insert(0,'../three-body'); import tbpy
from physics import solid_angle, middle, line_angle
d={o['cls']:o for o in json.load(open('../three-body-word-formula/data/li-liao-695.json'))}
def exact_cross(X,t,dt):
    a,b=0.0,dt; fa=tbpy.shape(X[None])[0,2]
    for _ in range(60):
        m=(a+b)/2; Xm=tbpy.flow(X,m) if m>0 else X; fm=tbpy.shape(Xm[None])[0,2]
        if np.sign(fm)==np.sign(fa): a,fa=m,fm
        else: b=m
    return tbpy.flow(X,(a+b)/2)
out={}
for c in ['I.A1','I.A2','I.B1','II.C2','I.A9','II.C4']:
    o=d[c]; n=40000; ts,X=tbpy.traj(tbpy.suvakov(o['v1'],o['v2']),o['T'],n); S=tbpy.shape(X); dt=ts[1]-ts[0]
    idx=np.where(np.sign(S[:-1,2])*np.sign(S[1:,2])<0)[0]
    Xc={i:exact_cross(X[i],ts[i],dt) for i in idx}
    gaps=[];rows=[]
    for k,i in enumerate(idx):
        mi=middle(Xc[i])
        for j in idx[k+1:]:
            if middle(Xc[j])==mi: break
        else: continue
        arc=np.vstack([tbpy.shape(Xc[i][None]),S[i+1:j+1],tbpy.shape(Xc[j][None])]); arc[0,2]=arc[-1,2]=0
        arc/=np.linalg.norm(arc,axis=1)[:,None]
        a0=np.arctan2(arc[-1,1],arc[-1,0]); a1=np.arctan2(arc[0,1],arc[0,0]); da=(a1-a0+np.pi)%(2*np.pi)-np.pi
        g=a0+da*np.linspace(0,1,400); eq=np.stack([np.cos(g),np.sin(g),0*g],1)
        Om=solid_angle(np.vstack([arc,eq[1:-1]])); turn=(line_angle(Xc[j])-line_angle(Xc[i]))%np.pi
        gm=(turn+Om/2)%np.pi; gm=min(gm,np.pi-gm); gaps.append(gm); rows.append([float(Om),float(turn)])
    out[c]=dict(worst=float(max(gaps)),median=float(np.median(gaps)),rows=rows)
    print(c,'pairs',len(gaps),'after polish: worst |turn + Omega/2| mod pi = %.1e, median %.1e'%(max(gaps),np.median(gaps)))
json.dump(out,open('swap_polished.json','w'))
