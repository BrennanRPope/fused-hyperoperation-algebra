import sys; sys.path.insert(0,'.')
import numpy as np, pandas as pd
from scipy.optimize import least_squares
from fit_mixtures import ct_odds
from scipy.stats import f as fdist
def pred(p,df):
    d,c,b1,e1,b2,e2,a=p
    return np.array([c+(d-c)/(1+ct_odds(D*pc/100,D*(100-pc)/100,e1,b1,e2,b2,a)) for D,pc in zip(df.dose,df.pct.astype(float))])
for nm in ['glymet','acidiq','mecter']:
    df=pd.read_csv(f'data/{nm}.csv'); df.loc[df.pct==999,'pct']=50; y=df.rgr.values
    e1=np.median(df[df.pct==100].dose); e2=np.median(df[df.pct==0].dose)
    def fitp(fix=None):
        p0=[y[df.dose==0].mean(),0.01,2,e1,2,e2,0.0]; lo=[0,-0.5,0.2,e1/1e3,0.2,e2/1e3,-3]; hi=[10,1,20,e1*1e3,20,e2*1e3,50]
        if fix is None:
            r=least_squares(lambda q:pred(q,df)-y,p0,bounds=(lo,hi),x_scale='jac'); return r.x,np.sum(r.fun**2)
        r=least_squares(lambda q:pred(list(q)+[fix],df)-y,p0[:6],bounds=(lo[:6],hi[:6]),x_scale='jac'); return list(r.x)+[fix],np.sum(r.fun**2)
    p,rss=fitp(); n=len(y)
    out=[]
    for a in (-1,-0.5,0):
        _,r=fitp(a); F=(r-rss)/(rss/(n-7)); out.append((a,1-fdist.cdf(F,1,n-7)))
    print(nm,'floor c=%.3f'%p[1],'alpha=%.2f'%p[6],'p(HSA,-0.5,Loewe)=',[f'{q:.2g}' for _,q in out])
