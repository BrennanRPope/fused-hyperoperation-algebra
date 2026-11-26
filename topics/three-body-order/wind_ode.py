"""DRAFT: pairwise windings by integrating d(angle)/dt = (r x v)/|r|^2 alongside the motion (no sampling gaps)."""
import numpy as np, json, sys
from scipy.integrate import solve_ivp
sys.path.insert(0,'../three-body'); import tbpy
PAIRS=[(0,1),(0,2),(1,2)]
def rhs(t,Y):
    r=Y[:6].reshape(3,2); v=Y[6:12].reshape(3,2); a=np.zeros((3,2)); th=np.zeros(3)
    for i,j in PAIRS:
        d=r[j]-r[i]; q=d@d; f=d/q**1.5; a[i]+=f; a[j]-=f
    for k,(i,j) in enumerate(PAIRS):
        d=r[j]-r[i]; u=v[j]-v[i]; th[k]=(d[0]*u[1]-d[1]*u[0])/(d@d)
    return np.concatenate([v.ravel(),a.ravel(),th])
def windings(p1,p2,T,rtol=1e-12):
    Y0=np.concatenate([tbpy.suvakov(p1,p2),np.zeros(3)])
    s=solve_ivp(rhs,(0,T),Y0,method='DOP853',rtol=rtol,atol=1e-13)
    w=s.y[12:,-1]/(2*np.pi); close=np.abs(s.y[:12,-1]-Y0[:12]).max()
    return dict(zip(['12','13','23'],w.tolist())),float(close)
if __name__=='__main__':
    tw=json.load(open('../three-body-word-formula/twins.json')); out={}
    for c,e in tw.items():
        a,ca=windings(*e['cat']); b,cb=windings(*e['twin'])
        out[c]=dict(cat=a,twin=b,close_cat=ca,close_twin=cb)
        print(c,{k:round(v,4) for k,v in a.items()},{k:round(v,4) for k,v in b.items()},'closure %.1e %.1e'%(ca,cb),flush=True)
    json.dump(out,open('twin_windings.json','w'),indent=1)
