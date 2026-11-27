import json, numpy as np
from scipy.optimize import curve_fit
def load(*fs):
    d={}
    for f in fs: d.update({int(k):v for k,v in json.load(open(f)).items()})
    n=np.array(sorted(d)); return n, np.array([d[k] for k in n])
out={}
for name,fs in [("symmetrized",("threshold.json","threshold_hi.json")),("doubled",("doubled_threshold.json",))]:
    n,b=load(*fs); m=n>=10; n,b=n[m],b[m]
    r={}
    f=lambda n,L,C,a: L-C*n**(-a)
    p,_=curve_fit(f,n,b,p0=[2,3,1],maxfev=20000); r["power"]=dict(L=p[0],C=p[1],a=p[2],rms=float(np.sqrt(np.mean((f(n,*p)-b)**2))))
    g=lambda n,L,C,q: L-C*q**n
    p,_=curve_fit(g,n,b,p0=[1.9,3,0.9],maxfev=20000); r["geometric"]=dict(L=p[0],C=p[1],q=p[2],rms=float(np.sqrt(np.mean((g(n,*p)-b)**2))))
    h=lambda n,A,B: A+B*np.log(n)
    p,_=curve_fit(h,n,b); r["log_unbounded"]=dict(A=p[0],B=p[1],rms=float(np.sqrt(np.mean((h(n,*p)-b)**2))))
    # power with L fixed at 2
    f2=lambda n,C,a: 2-C*n**(-a)
    p,_=curve_fit(f2,n,b,p0=[3,1]); r["power_L2"]=dict(C=p[0],a=p[1],rms=float(np.sqrt(np.mean((f2(n,*p)-b)**2))))
    # leave-last-out stability of power L
    p,_=curve_fit(f,n[:-2],b[:-2],p0=[2,3,1],maxfev=20000); r["power_L_dropping_last2"]=p[0]
    r["ndata"]=[int(n[0]),int(n[-1])]
    out[name]=r
print(json.dumps(out,indent=1,default=float)); json.dump(out,open("fitlim.json","w"),indent=1,default=float)
