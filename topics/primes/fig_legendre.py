import numpy as np, matplotlib, itertools, math
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from sympy import primerange, primepi
plt.rcParams.update({'font.size':11})
x=1000; ps=list(primerange(2,32)); a=len(ps)
layer=[0]*(a+1); cnt=[0]*(a+1)
for mask in range(1<<a):
    d=1; k=0
    for i in range(a):
        if mask>>i&1: d*=ps[i]; k+=1
    layer[k]+=(-1)**k*(x//d); cnt[k]+=1
S=np.cumsum(layer); true=int(primepi(x))-a+1
print('layers',layer,'cum',S.tolist(),'true',true, 'corners',sum(cnt))
fig,axs=plt.subplots(1,2,figsize=(15,5.8))
ax=axs[0]; ks=np.arange(a+1)
ax.axhline(true,color='k',lw=1,ls='--'); ax.text(a+0.1,true+25,f'exact: {true} numbers ≤ 1000\nuntouched by primes up to 31\n(= the primes from 37 to 997, plus 1)',fontsize=9,ha='right')
ax.plot(ks,S,'-',color='#999',lw=1)
for k in ks:
    ax.plot(k,S[k],'o',ms=10,color='#2a9d8f' if k%2==0 else '#e76f51',mec='k')
ax.set_xticks(ks); ax.set_xlabel('layers of the 11-dimensional cube added so far (primes 2…31, 2048 corners)')
ax.set_ylabel('running count')
for k in range(4): ax.text(k+0.15,S[k],f'{S[k]}',fontsize=9,va='center')
ax.set_title('Legendre\'s count is the cube summed layer by layer.\nEven layers overshoot, odd layers undershoot, and it settles exactly.\n(Layers 5 to 11 add nothing: those corners are products bigger than 1000.)',fontsize=11)
for s in ['top','right']: ax.spines[s].set_visible(False)
ax=axs[1]
xs=np.unique(np.logspace(3,8,40).astype(int))
from sympy import sieve
sieve.extend(10**8)
ratio=[]
for X in xs:
    pr=np.array(list(primerange(2,int(math.isqrt(X))+1)))
    smooth=X*np.prod(1-1/pr)
    exact=int(primepi(X))-len(pr)+1
    ratio.append(smooth/exact)
ax.semilogx(xs,ratio,'-o',ms=4,color='#7b5ea7')
c=2*math.exp(-0.5772156649)
ax.axhline(c,color='#d62728',ls='--'); ax.text(xs[0],c+0.004,f'2e^(−γ) ≈ {c:.4f}',color='#d62728',fontsize=10)
ax.axhline(1,color='k',lw=0.8)
ax.set_xlabel('x'); ax.set_ylabel('smooth cube estimate ÷ true count')
ax.set_title('Replace each corner\'s count ⌊x/d⌋ with the smooth x/d and the cube\ndrifts above the truth, creeping toward a 12% overcount (Mertens 2e^(−γ)).\nThe rounding at each corner is not noise: it carries real information.',fontsize=11)
for s in ['top','right']: ax.spines[s].set_visible(False)
fig.tight_layout(); fig.savefig('figures/p6_legendre.png',dpi=150); plt.close(fig)
print(ratio[-3:])
