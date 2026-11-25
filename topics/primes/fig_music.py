import numpy as np, matplotlib, math
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from sympy import primerange, factorint
plt.rcParams.update({'font.size':11})
G=np.load('zeta_zeros_1000.npy')
def psi_true(x):
    # chebyshev psi
    tot=0.0
    for p in primerange(2,int(x)+1):
        pk=p
        while pk<=x: tot+=math.log(p); pk*=p
    return tot
def psi_N(x,N):
    rho=0.5+1j*G[:N]
    s=np.sum(2*np.real(np.exp(rho[None,:]*np.log(x)[:,None])/rho[None,:]),axis=1)
    return x - s - math.log(2*math.pi) - 0.5*np.log(1-x**-2.0)
fig,axs=plt.subplots(2,1,figsize=(15,10))
ax=axs[0]
x=np.linspace(1.5,50,20000); lx=np.log(x)
f=np.zeros_like(x)
for g in G: f-=np.cos(g*lx)
ax.plot(x,f,color='#4b2e83',lw=0.7)
pp=[n for n in range(2,51) if len(factorint(n))==1]
for n in pp:
    isp=factorint(n)[list(factorint(n))[0]]==1
    ax.axvline(n,color='goldenrod' if isp else '#aaa',lw=0.8,ls='-' if isp else ':',zorder=0)
    ax.text(n,f.max()*1.05,str(n),ha='center',fontsize=8,color='k' if isp else '#888')
ax.set_xlim(1.5,50); ax.set_ylim(f.min()*1.1,f.max()*1.15)
ax.set_xlabel('x'); ax.set_yticks([])
ax.set_title('The music of the primes: add up 1,000 pure waves, one per zeta zero  (−Σ cos(γ·log x)).\nNothing about primes is put in by hand, yet the sum spikes at every prime (gold) and prime power (dotted grey).',fontsize=11.5)
for s in ['top','right','left']: ax.spines[s].set_visible(False)
ax=axs[1]
x=np.linspace(2.0,60,6000)
tr=np.array([psi_true(v) for v in x])
ax.plot(x,tr,color='k',lw=2.2,label='the true prime staircase ψ(x): jumps by log p at each prime power')
for N,c in [(10,'#e9c46a'),(100,'#e76f51'),(1000,'#2a9d8f')]:
    ax.plot(x,psi_N(x,N),color=c,lw=1.1,label=f'x minus {N} zero-waves')
ax.plot(x,x,'--',color='#999',lw=1,label='smooth trend x (no zeros)')
ax.legend(frameon=False,fontsize=9.5,loc='upper left')
ax.set_xlabel('x'); ax.set_ylabel('ψ(x)')
ax.set_title('The explicit formula: start with the straight line x and subtract one wave per zero. With 1,000 zeros\nthe line has bent itself into the staircase, step for step. Each step is a prime (or prime power) being born.',fontsize=11.5)
for s in ['top','right']: ax.spines[s].set_visible(False)
fig.tight_layout(); fig.savefig('figures/p7_music.png',dpi=150); plt.close(fig)
xx=np.array([29.5,30.5,31.5]); print(psi_N(xx,1000),[psi_true(v) for v in xx])
