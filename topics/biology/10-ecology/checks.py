import numpy as np, json
from scipy.optimize import least_squares
rng=np.random.default_rng(10); res=[]; out={}
def check(name,ok,info=""):
    res.append(bool(ok)); print(("PASS " if ok else "FAIL ")+name+(" | "+info if info else ""))
wadd=lambda k,a,b:a+b+k*a*b

# 1. Functional response: predator time budget as a Markov chain
def feed_exact(N,P,al,h,ga,tau,k):
    x=al*N*h; y=ga*P*tau
    return al*N/(1+wadd(k,x,y))
def feed_markov(N,P,al,h,ga,tau,k):
    # states: S searching, H handling, IS interfering (from search), IH interfering while handling
    Q=np.zeros((4,4))
    Q[0,1]=al*N; Q[1,0]=1/h; Q[0,2]=ga*P; Q[2,0]=1/tau; Q[1,3]=k*ga*P; Q[3,1]=1/tau
    Q-=np.diag(Q.sum(1)); A=np.vstack([Q.T,np.ones(4)]); b=np.r_[np.zeros(4),1]
    pi=np.linalg.lstsq(A,b,rcond=None)[0]; return al*N*pi[0]
pars=dict(al=0.8,h=0.3,ga=0.5,tau=0.6)
worst=max(abs(feed_markov(N,P,k=k,**pars)-feed_exact(N,P,k=k,**pars))/feed_exact(N,P,k=k,**pars)
          for N in (0.5,3,20) for P in (1,4,10) for k in (0,0.4,1,2.5))
check("a predator that searches, handles, and is drawn into interference from either state feeds at αN/(1 + x ∘ₖ y), x = αNh, y = γPτ",worst<1e-12,
      f"steady state of the 4-state chain vs formula: largest relative gap {worst:.1e}")
check("k = 0 (no interference while handling) is Beddington–DeAngelis, 1 + x + y",
      abs(feed_exact(3,4,k=0,**pars)-pars['al']*3/(1+pars['al']*3*pars['h']+pars['ga']*4*pars['tau']))<1e-14)
check("k = 1 (handlers interfere as often as searchers) is Crowley–Martin, (1 + x)(1 + y)",
      abs(feed_exact(3,4,k=1,**pars)-pars['al']*3/((1+pars['al']*3*pars['h'])*(1+pars['ga']*4*pars['tau'])))<1e-14)
# Gillespie simulation of one predator: realised feeding rate matches
def gillespie(N,P,al,h,ga,tau,k,T=2e5):
    s=0; t=0; eaten=0
    rates={0:[(1,al*N),(2,ga*P)],1:[(0,1/h),(3,k*ga*P)],2:[(0,1/tau)],3:[(1,1/tau)]}
    while t<T:
        rs=rates[s]; tot=sum(r for _,r in rs); t+=rng.exponential(1/tot)
        u=rng.uniform(0,tot); acc=0
        for nxt,r in rs:
            acc+=r
            if u<=acc:
                if s==0 and nxt==1: eaten+=1
                s=nxt; break
    return eaten/T
g=gillespie(5,6,k=0.5,**pars); f=feed_exact(5,6,k=0.5,**pars)
check("a stochastic simulation of one predator (200,000 time units) eats at the formula's rate",abs(g-f)/f<0.02,f"simulated {g:.4f}, formula {f:.4f}")
# Fitting: with k = 0.5 truth, BD and CM are biased; the dial recovers k
Ns=np.repeat([0.5,1,2,4,8,16],5); Ps=np.tile([1,2,4,8,12],6); ktrue=0.5
y=feed_exact(Ns,Ps,k=ktrue,**pars)*np.exp(rng.normal(0,0.04,len(Ns)))
def model(p,k=None):
    al,h,c=np.exp(p[:3]); kk=p[3] if k is None else k
    x=al*Ns*h; z=c*Ps; return al*Ns/(1+x+z+kk*x*z)
fd=least_squares(lambda p:np.log(model(p))-np.log(y),[0,0,0,0.3])
fb=least_squares(lambda p:np.log(model(p,0))-np.log(y),[0,0,0])
fc=least_squares(lambda p:np.log(model(p,1))-np.log(y),[0,0,0])
n=len(y); aic=lambda r,kp:n*np.log(np.mean(r**2))+2*kp
A=dict(dial=aic(fd.fun,4),BD=aic(fb.fun,3),CM=aic(fc.fun,3))
check("in a simulated feeding experiment crossing prey and predator density (true k = 0.5, 4% noise) the dial fit recovers k, and both classic forms fit worse",
      abs(fd.x[3]-ktrue)<0.25 and A['dial']<min(A['BD'],A['CM']),
      f"fitted k = {fd.x[3]:.2f}; AIC dial {A['dial']:.1f}, BD {A['BD']:.1f}, CM {A['CM']:.1f}")
out['fr']=dict(kfit=fd.x[3],aic=A)

# 2. Diversity: Tsallis entropy is the q-logarithm of the Hill number, and combines by the dial k = 1 − q
lnq=lambda x,q:(x**(1-q)-1)/(1-q)
def tsallis(p,q): p=np.asarray(p); return (1-np.sum(p**q))/(q-1)
def hill(p,q): p=np.asarray(p); return np.sum(p**q)**(1/(1-q))
pA=rng.dirichlet(np.ones(6)); pB=rng.dirichlet(np.ones(4)); pAB=np.outer(pA,pB).ravel()
q=2.0
check("Tsallis (HCDT) entropy is the q-logarithm of Hill's effective number: S_q = ln_q(D_q)",abs(tsallis(pA,q)-lnq(hill(pA,q),q))<1e-14)
check("two independent communities: Hill numbers multiply and Tsallis entropies combine by the dial with k = 1 − q",
      abs(hill(pAB,q)-hill(pA,q)*hill(pB,q))<1e-12 and abs(tsallis(pAB,q)-wadd(1-q,tsallis(pA,q),tsallis(pB,q)))<1e-14)
vals=[tsallis(np.ones(m)/m,q) for m in (10,100,10000)]
check("for q > 1 the wall −1/k = 1/(q − 1) is the ceiling: Tsallis entropy approaches it as richness grows and never passes it",
      all(v<1/(q-1) for v in vals) and 1/(q-1)-vals[-1]<1e-3,f"q = 2: {', '.join(f'{v:.4f}' for v in vals)} → wall {1/(q-1):.1f}")
check("q = 1 is dial 0: Shannon entropies of independent communities add",
      abs(-np.sum(pAB*np.log(pAB))-(-np.sum(pA*np.log(pA))-np.sum(pB*np.log(pB))))<1e-12)

# 3. Multiple stressors: the two nulls are dial 0 and dial −1 on proportional losses
l1,l2=0.3,0.2
check("multiple-stressor nulls: 'simple addition' of losses is dial 0, 'multiplicative' (independent action) is the 'or' law, dial −1",
      abs((1-(1-l1)*(1-l2))-wadd(-1,l1,l2))<1e-15)
# 4. Island colonisation
check("an island receiving colonists from n sources, each succeeding with chance c, is colonised with wpow(−1, c, n)",abs((1-0.95**8)-(lambda f:f(f,8))(lambda f,m:0 if m==0 else wadd(-1,f(f,m-1),0.05)))<1e-15)
json.dump(out,open('results.json','w'),indent=1,default=float)
print(f"\n{sum(res)}/{len(res)} checks passed")
