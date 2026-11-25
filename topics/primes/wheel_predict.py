import numpy as np, matplotlib, math
matplotlib.use('Agg')
import matplotlib.pyplot as plt
N=10**8
s=np.ones(N+200,bool); s[:2]=False
for i in range(2,int((N+200)**.5)+1):
    if s[i]: s[i*i::i]=False
P=np.nonzero(s)[0]
P=P[(P>10**6)&(P<N)]
nxt=np.nonzero(s)[0]; idx=np.searchsorted(nxt,P); nextp=nxt[idx+1]
out=[]
M=1
for p in [2,3,5,7,11,13,17,19,23]:
    M*=p
    if M<6: continue
    # survivors mod M as bool
    if M>2*10**8: break
    w=np.ones(M,bool)
    for q in [2,3,5,7,11,13,17,19,23]:
        if M%q==0: w[::q]=False
    surv=np.nonzero(w)[0]
    r=P%M; j=np.searchsorted(surv,r,side='right')
    ext=np.concatenate([surv,[surv[0]+M]])
    nextsurv=P-r+ext[j]
    hit=(nextsurv==nextp).mean()
    out.append((p,M,len(surv),hit)); print(p,M,len(surv),round(hit,4))
# naive: next odd
hit_odd=(P+2==nextp).mean(); print('next odd',hit_odd)
np.save('wheel_predict.npy',np.array(out))
fig,ax=plt.subplots(figsize=(10,5.2))
labs=['next odd\nnumber']+[f'wheel\nafter {o[0]}' for o in out[1:]]
vals=[hit_odd]+[o[3] for o in out[1:]]
ax.bar(range(len(vals)),vals,color='#7b5ea7')
import math
pred=[]
for o in out[1:]:
    from sympy import primerange
    prod=np.prod([1-1/q for q in primerange(2,o[0]+1)])
    pred.append(np.mean(1/np.log(P[::50]))/prod)
ax.plot(range(1,len(vals)),pred,'o--',color='#e76f51',ms=6,label='prediction: (1/log x) ÷ fraction of spots the wheel leaves open')
ax.legend(frameon=False,fontsize=9,loc='upper left')
for i,v in enumerate(vals): ax.text(i,v-0.025,f'{100*v:.1f}%',ha='center',fontsize=10,color='w',weight='bold')
ax.set_xticks(range(len(vals))); ax.set_xticklabels(labs,fontsize=9)
ax.set_ylabel('how often the guess is the next prime')
ax.set_ylim(0,max(vals)*1.25)
ax.set_title('Guessing where the next prime is born: "the next open spot on the wheel".\nTested on all 5.7 million primes between 1 million and 100 million.\nEach new prime on the wheel adds a few points. Certainty needs every prime up to √x: the full sieve.',fontsize=11)
for sp in ['top','right']: ax.spines[sp].set_visible(False)
fig.tight_layout(); fig.savefig('figures/p10_wheel_guess.png',dpi=150)
