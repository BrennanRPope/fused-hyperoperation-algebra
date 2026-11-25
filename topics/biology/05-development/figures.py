import sys, json, numpy as np, pandas as pd
sys.path.insert(0,"../tools"); from style import *
from scipy.optimize import least_squares
setup()
exec('\n'.join(open('checks.py').read().split('# 2. Real growth data')[1].split('out={}')[0].split('\n')[1:]))  # richards, fit, sets
fits=json.load(open('fits.json'))
fig,ax=plt.subplots(1,4,figsize=(14,3.9))
t=np.linspace(-4,4,300)
for nu,c,lab in ((-1/3,S3,"von Bertalanffy, power 3"),(1e-9,S4,"Gompertz, ν → 0"),(1,S1,"logistic, power −1 (inverse)"),(2.5,S2,"ν = 2.5, power −0.4")):
    ax[0].plot(t,richards([1,1,0,nu],t),color=c,label=lab)
ax[0].set_title("Richards family: 1 + mpow(x, −1/ν)"); ax[0].set_xlabel("time (scaled)"); ax[0].set_ylabel("size / asymptote"); ax[0].legend(fontsize=7.8,loc="upper left")
for a,(nm,(tc,yc,lab)) in zip(ax[1:],sets.items()):
    df=pd.read_csv(f'data/{nm}.csv'); tt=df[tc].values.astype(float); yy=df[yc].values.astype(float)
    a.plot(tt,yy,"o",ms=3.5,color=INK,alpha=0.35)
    T=np.linspace(tt.min(),tt.max(),200); f=fits[nm]
    a.plot(T,richards(f['p'],T),color=S2,label=f"best ν = {f['nu']:.2f}")
    pL,_=fit(tt,yy,fixnu=1); a.plot(T,richards(pL,T),"--",color=S1,label="logistic")
    pG,_=fit(tt,yy,fixnu=1e-9); a.plot(T,richards(pG,T),":",color=S4,label="Gompertz")
    a.set_title(lab,fontsize=10.5); a.set_xlabel(tc); a.legend(fontsize=8,loc="upper left")
save(fig,"figures/f1_growth_family.png")
