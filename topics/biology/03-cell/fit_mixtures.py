"""Fit the interaction dial to three real binary-mixture toxicity data sets (drc package)."""
import numpy as np, pandas as pd
from scipy.optimize import least_squares, brentq
from scipy.stats import f as fdist

def odds_single(D,e,b): return np.where(D>0,(np.maximum(D,1e-300)/e)**b,0.0)

def pred_odds_dial(p,df):
    d,b1,e1,b2,e2,k=p
    D=df.dose.values; pct=df.pct.values.astype(float)
    D1=D*pct/100; D2=D*(100-pct)/100
    o1=odds_single(D1,e1,b1); o2=odds_single(D2,e2,b2)
    o=o1+o2+k*o1*o2
    return d/(1+np.maximum(o,-0.999))

def ct_odds(D1,D2,e1,b1,e2,b2,a):
    # Chou-Talalay: v1 + v2 + a v1 v2 = 1, v_i = D_i/(e_i o^(1/b_i)); solve for effect odds o
    if D1<=0 and D2<=0: return 0.0
    if D2<=0: return (D1/e1)**b1
    if D1<=0: return (D2/e2)**b2
    def g(lo):
        v1=D1/(e1*np.exp(lo/b1)); v2=D2/(e2*np.exp(lo/b2)); return v1+v2+a*v1*v2-1
    # breakpoints: where each drug alone reaches the effect (v=1) or its wall (v=-1/a)
    pts=[b1*np.log(D1/e1),b2*np.log(D2/e2)]
    if a!=0 and -1/a>0: pts+=[b1*np.log(D1/(e1*(-1/a))),b2*np.log(D2/(e2*(-1/a)))]
    grid=np.unique(np.concatenate([np.linspace(60,-60,2401),np.clip(pts,-60,60)]))[::-1]
    gv=g(grid)
    for i in range(len(grid)-1):
        if gv[i]<0<=gv[i+1]:
            if gv[i+1]==0: return np.exp(grid[i+1])
            return np.exp(brentq(g,grid[i+1],grid[i],xtol=1e-14))
    return np.exp(-60.0)  # combination never reaches additivity line (antagonism past the wall)

def pred_ct(p,df):
    d,b1,e1,b2,e2,a=p
    out=[]
    for D,pct in zip(df.dose.values,df.pct.values.astype(float)):
        o=ct_odds(D*pct/100,D*(100-pct)/100,e1,b1,e2,b2,a); out.append(d/(1+o))
    return np.array(out)

def fit(df,model,fix=None):
    df=df.copy(); df.loc[df.pct==999,'pct']=50; df.loc[df.dose==0,'dose']=0.0
    pr=pred_odds_dial if model=='odds' else pred_ct
    y=df.rgr.values
    s1=df[df.pct==100]; s2=df[df.pct==0]
    e1=np.median(s1.dose); e2=np.median(s2.dose)
    p0=[y[df.dose==0].mean(),2,e1,2,e2,0.0]
    lo=[0,0.2,e1/1e3,0.2,e2/1e3,-0.9 if model=='odds' else -3.0]; hi=[10,20,e1*1e3,20,e2*1e3,50]
    if fix is not None:
        f=lambda q: pr(list(q)+[fix],df)-y
        r=least_squares(f,p0[:5],bounds=(lo[:5],hi[:5]),x_scale='jac')
        return list(r.x)+[fix], np.sum(r.fun**2), len(y)
    f=lambda q: pr(q,df)-y
    r=least_squares(f,p0,bounds=(lo,hi),x_scale='jac')
    return list(r.x), np.sum(r.fun**2), len(y)

def profile_ci(df,model,p,rss,n,grid):
    out=[]
    for k in grid:
        _,r,_=fit(df,model,fix=k)
        F=(r-rss)/(rss/(n-6)); out.append((k,F))
    crit=fdist.ppf(0.95,1,n-6)
    inside=[k for k,F in out if F<=crit]
    return (min(inside),max(inside)) if inside else (np.nan,np.nan), out

if __name__=="__main__":
    for name in ['glymet','acidiq','mecter']:
        df=pd.read_csv(f'data/{name}.csv')
        for model in ['odds','ct']:
            p,rss,n=fit(df,model)
            print(name,model,'k=%.3f'%p[5],'slopes %.2f %.2f'%(p[1],p[3]),'rss %.4g n %d'%(rss,n))
