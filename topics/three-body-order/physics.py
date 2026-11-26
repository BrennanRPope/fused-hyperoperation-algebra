"""DRAFT: physical order readings of real orbits.
1) swap rule on the shape sphere: the line-up's turn between two line-ups with the same middle body = half the enclosed solid angle (mod pi)
2) pairwise windings per period (the whole-number part of the braid, which includes the overall twist the word leaves out)"""
import sys, json, numpy as np
sys.path.insert(0,'../three-body'); import tbpy
d=json.load(open('../three-body-word-formula/data/li-liao-695.json')); byc={o['cls']:o for o in d}
def sample(p1,p2,T,n):
    X0=tbpy.suvakov(p1,p2); ts,X=tbpy.traj(X0,T,n); return ts,X
def solid_angle(P):
    p=np.array([0,0,1.0]); tot=0
    for a,b in zip(P,np.roll(P,-1,axis=0)):
        tot+=2*np.arctan2(p@np.cross(a,b),1+p@a+p@b+a@b)
    return tot
def middle(Xrow):
    r=Xrow[:6].reshape(3,2); u=r[1]-r[0]; 
    # project on the line, find which body lies between the other two
    ax=(r[np.argmax([np.linalg.norm(r[i]-r[j]) for i,j in [(0,1),(0,2),(1,2)] ]==0)] if False else None)
    s=[np.dot(r[i]-r.mean(0), (r[0]-r[2] if np.linalg.norm(r[0]-r[2])>1e-9 else u)) for i in range(3)]
    return int(np.argsort(s)[1])
def line_angle(Xrow):
    z1,z2=tbpy.jacobi(Xrow[None,:]); z=z1[0] if abs(z1[0])>abs(z2[0]) else z2[0]
    return np.angle(z)%np.pi
def swap_rule(p1,p2,T,n=40000):
    ts,X=sample(p1,p2,T,n); S=tbpy.shape(X)
    idx=np.where(np.sign(S[:-1,2])*np.sign(S[1:,2])<0)[0]
    rows=[]
    for k,i in enumerate(idx):
        # next line-up with the same middle body
        mi=middle(X[i]); 
        for j in idx[k+1:]:
            if middle(X[j])==mi: break
        else: continue
        arc=S[i:j+1].copy(); arc[0,2]=0; arc[-1,2]=0; arc/=np.linalg.norm(arc,axis=1)[:,None]
        a0=np.arctan2(arc[-1,1],arc[-1,0]); a1=np.arctan2(arc[0,1],arc[0,0])
        # equator path back from end to start, short way (both points lie on the same collision-free arc)
        da=(a1-a0+np.pi)%(2*np.pi)-np.pi; g=a0+da*np.linspace(0,1,200)
        eq=np.stack([np.cos(g),np.sin(g),0*g],1)
        Om=solid_angle(np.vstack([arc,eq[1:-1]]))
        turn=(line_angle(X[j])-line_angle(X[i]))%np.pi
        gap=min(((turn-s*Om/2)%np.pi, np.pi-(turn-s*Om/2)%np.pi) for s in (1,-1))
        gp=min((turn-Om/2)%np.pi, np.pi-(turn-Om/2)%np.pi); gm=min((turn+Om/2)%np.pi, np.pi-(turn+Om/2)%np.pi)
        rows.append(dict(t0=float(ts[i]),t1=float(ts[j]),Omega=float(Om),turn=float(turn),gap_plus=float(gp),gap_minus=float(gm)))
    return rows
def windings(p1,p2,T,n=200000):
    ts,X=sample(p1,p2,T,n); r=X[:,:6].reshape(-1,3,2); out={}; worst=0
    for i,j in [(0,1),(0,2),(1,2)]:
        v=r[:,j]-r[:,i]; a=np.unwrap(np.arctan2(v[:,1],v[:,0])); worst=max(worst,np.abs(np.diff(a)).max())
        out[f'{i+1}{j+1}']=(a[-1]-a[0])/(2*np.pi)
    return out,worst
if __name__=='__main__':
    res={'swap':{},'wind':{},'twins':{}}
    for c in ['I.A1','I.A2','I.B1','II.C2','I.A9','II.C4']:
        o=byc[c]; rows=swap_rule(o['v1'],o['v2'],o['T'])
        gp=max(r['gap_plus'] for r in rows); gm=max(r['gap_minus'] for r in rows)
        res['swap'][c]=dict(n=len(rows),worst_gap_plus=gp,worst_gap_minus=gm,rows=rows[:40])
        print(c,'line-up pairs',len(rows),'worst |turn - Omega/2| mod pi = %.2e   |turn + Omega/2| = %.2e'%(gp,gm))
    for o in d[:120]:
        w,worst=windings(o['v1'],o['v2'],o['T'],n=max(20000,4000*len(o['word'])))
        res['wind'][o['cls']]=dict(w=w,worst_step=worst,word=o['word'])
    tw=json.load(open('../three-body-word-formula/twins.json'))
    for c,e in tw.items():
        a,_=windings(*e['cat'],n=400000); b,_=windings(*e['twin'],n=400000)
        res['twins'][c]=dict(cat=a,twin=b); print('twin',c,{k:round(v,3) for k,v in a.items()},{k:round(v,3) for k,v in b.items()})
    json.dump(res,open('physics.json','w'))
