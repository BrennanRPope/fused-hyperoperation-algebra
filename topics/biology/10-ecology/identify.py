import numpy as np, json
from scipy.optimize import least_squares
from checks import feed_exact, pars
rng=np.random.default_rng(11)
out={}
for reps in (5,20):
  Ns=np.repeat([0.5,1,2,4,8,16],reps); Ps=np.tile(np.repeat([1,2,4,8,12],1),6*reps//5)
  for sd in (0.1,0.25):
    ks=[]
    for b in range(200):
        y=feed_exact(Ns,Ps,k=0.5,**pars)*np.exp(rng.normal(0,sd,len(Ns)))
        def model(p):
            al,h,c=np.exp(p[:3]); x=al*Ns*h; z=c*Ps; return al*Ns/(1+x+z+p[3]*x*z)
        f=least_squares(lambda p:np.log(model(p))-np.log(y),[0,0,0,0.3],bounds=([-5,-5,-5,-2],[5,5,5,10]))
        ks.append(f.x[3])
    lo,hi=np.percentile(ks,[2.5,97.5]); out[f"{len(Ns)} trials, {int(sd*100)}% noise"]=[lo,hi]
    print(len(Ns),sd,round(lo,2),round(hi,2))
json.dump(out,open('identify.json','w'),indent=1)
