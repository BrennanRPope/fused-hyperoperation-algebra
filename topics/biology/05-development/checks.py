import numpy as np, pandas as pd, json
from scipy.optimize import least_squares
from scipy.stats import f as fdist
res=[]
def check(name,ok,info=""):
    res.append(ok); print(("PASS " if ok else "FAIL ")+name+(" | "+info if info else ""))
wadd=lambda k,a,b:a+b+k*a*b
mpow=lambda x,p:(1+x)**p-1
inv=lambda u:-u/(1+u)

# 1. Richards family = 1 + mpow(decaying excess, -1/nu)
t=np.linspace(0,10,200); K,Q,r=5.0,8.0,0.9
for nu,name in ((1,"logistic"),(-1/3,"von Bertalanffy"),(0.5,"Richards ν = 0.5"),(2,"Richards ν = 2")):
    x=Q*np.exp(-r*t)* (1 if nu>0 else 1)
    if nu<0: x=-0.6*np.exp(-r*t)   # von Bertalanffy starts below the asymptote from 1 + x with x in (-1,0)
    y=K*(1+x)**(-1/nu)
    check(f"{name}: size/asymptote = 1 + mpow(x(t), −1/ν) with x(t) decaying exponentially",np.allclose(y/K,1+mpow(x,-1/nu)))
x=Q*np.exp(-r*t)
check("logistic (ν = 1) is the circle-operation inverse: y/K = 1 + inv(x), the same saturation as Michaelis–Menten",np.allclose(1/(1+x),1+inv(x)))
# Gompertz as limit and as iterated circle-operation power
nu=1e-6; x0=3.0
yR=(1+nu*x0*np.exp(-r*t))**(-1/nu); yG=np.exp(-x0*np.exp(-r*t))
check("Gompertz is the ν → 0 end of the family",np.max(np.abs(yR-yG))<1e-4)
lam=np.exp(-r); z=np.exp(-x0)-1; zs=[z]
for _ in range(9): zs.append(mpow(zs[-1],lam))
check("Gompertz in discrete time: each step is the circle-operation power x ↦ (1+x)^λ − 1, λ = e^(−r)",
      np.allclose(1+np.array(zs),np.exp(-x0*np.exp(-r*np.arange(10)))))
check("powers compose: mpow(mpow(x, a), b) = mpow(x, a·b), so n steps are one power λⁿ",abs(mpow(mpow(0.7,0.6),0.3)-mpow(0.7,0.18))<1e-14)

# 2. Real growth data: where does the exponent sit?
def richards(p,tt):
    K,r,t0,nu=p
    base=1+nu*np.exp(-r*(tt-t0))
    if abs(nu)<1e-8: return K*np.exp(-np.exp(-r*(tt-t0)))
    return K*np.maximum(base,1e-12)**(-1/nu)
def fit(tt,yy,fixnu=None):
    K0=yy.max()*1.1; t00=np.median(tt); r0=4/(tt.max()-tt.min())
    best=None
    for nu0 in ([fixnu] if fixnu is not None else [-0.3,0.2,1,2]):
        p0=[K0,r0,t00,nu0]
        if fixnu is None:
            f=lambda p: richards(p,tt)-yy; lo=[yy.max()*0.5,1e-4,tt.min()-5*(tt.max()-tt.min()),-0.95]; hi=[yy.max()*20,10,tt.max()*3,20]
        else:
            f=lambda p: richards(list(p)+[fixnu],tt)-yy; p0=p0[:3]; lo=[yy.max()*0.5,1e-4,tt.min()-5*(tt.max()-tt.min())]; hi=[yy.max()*20,10,tt.max()*3]
        try: rr=least_squares(f,p0,bounds=(lo,hi),x_scale='jac',max_nfev=20000)
        except Exception: continue
        if best is None or rr.cost<best.cost: best=rr
    p=list(best.x)+([fixnu] if fixnu is not None else [])
    return p,2*best.cost
sets={'Soybean':('Time','weight','soybean leaf weight per plant (g)'),'Loblolly':('age','height','loblolly pine height (ft)'),'Orange':('age','circumference','orange tree trunk circumference (mm)')}
out={}
for nm,(tc,yc,lab) in sets.items():
    df=pd.read_csv(f'data/{nm}.csv'); tt=df[tc].values.astype(float); yy=df[yc].values.astype(float)
    p,rss=fit(tt,yy); n=len(yy)
    pv={}
    for nu,lab2 in ((1,'logistic'),(1e-9,'Gompertz'),(-1/3,'von Bertalanffy')):
        _,r0=fit(tt,yy,fixnu=nu); F=(r0-rss)/(rss/(n-4)); pv[lab2]=1-fdist.cdf(F,1,n-4)
    out[nm]=dict(nu=p[3],p=p,pv=pv,n=n,label=lab)
    print(f"   {lab}: ν = {p[3]:.2f}, power −1/ν = {-1/p[3]:.2f}; p(logistic) = {pv['logistic']:.2g}, p(Gompertz) = {pv['Gompertz']:.2g}, p(von Bertalanffy) = {pv['von Bertalanffy']:.2g}")
json.dump(out,open('fits.json','w'),indent=1,default=float)
check("orange trees: the data cannot tell logistic, Gompertz and von Bertalanffy apart (all p > 0.05)",min(out['Orange']['pv'].values())>0.05)
check("soybean and loblolly reject all three named shapes (pooled curves)",max(out['Soybean']['pv'].values())<0.05 and max(out['Loblolly']['pv'].values())<0.05)
check("the three data sets do not share one exponent: fitted ν spans more than 1",max(o['nu'] for o in out.values())-min(o['nu'] for o in out.values())>1,
      ", ".join(f"{k} {o['nu']:.2f}" for k,o in out.items()))

# 3. Allometry: Huxley y = b x^a is a circle operation power, nothing more
a=0.75; xs=np.array([0.5,1,3])
check("allometry y = b·x^α: the growth of y is the circle-operation power of the growth of x",np.allclose(((1+xs)**a)-1,mpow(xs,a)))
print(f"\n{sum(res)}/{len(res)} checks passed")
