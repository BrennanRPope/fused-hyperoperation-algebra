import numpy as np, matplotlib, math
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.integrate import quad
plt.rcParams.update({'font.size':11})
N=10**8
s=np.ones(N+1,bool); s[:2]=False
for i in range(2,int(N**.5)+1):
    if s[i]: s[i*i::i]=False
P=np.nonzero(s)[0]; print('pi',len(P))
small=P[P<5000]
def S(H):
    k=len(H); val=1.0
    for p in small:
        nu=len(set(h%p for h in H))
        val*=(1-nu/p)/(1-1/p)**k
    return val
pats={'twins\n(n, n+2)':(0,2),'cousins\n(n, n+4)':(0,4),'"sexy"\n(n, n+6)':(0,6),'(n, n+2,\nn+6)':(0,2,6),'(n, n+4,\nn+6)':(0,4,6),'quadruplets\n(n,n+2,n+6,n+8)':(0,2,6,8)}
rows=[]
for name,H in pats.items():
    m=max(H); ok=np.ones(len(s)-m,bool)
    for h in H: ok&=s[h:len(s)-m+h]
    actual=int(ok.sum())
    k=len(H); integ=quad(lambda t: 1/math.log(t)**k,2,N,limit=500)[0]
    rows.append((name,H,actual,S(H)*integ,integ))
    print(H,actual,round(S(H),4),round(S(H)*integ),round(integ))
fig,axs=plt.subplots(1,3,figsize=(18,6.2),gridspec_kw={'width_ratios':[1.25,1,1]})
ax=axs[0]; xs=np.arange(len(rows)); w=0.27
ax.bar(xs-w,[r[2] for r in rows],w,color='#264653',label='actual count up to 10⁸')
ax.bar(xs,[r[3] for r in rows],w,color='#e9c46a',label='prediction: product over primes × smooth density')
ax.bar(xs+w,[r[4] for r in rows],w,color='#ccc',label='naive coin-flip model (no wheel correction)')
ax.set_yscale('log'); ax.set_xticks(xs); ax.set_xticklabels([r[0] for r in rows],fontsize=8.5)
for i,r in enumerate(rows): ax.text(i,r[2]*1.5,f'{r[2]/r[3]:.3f}',ha='center',fontsize=9,color='#264653')
ax.legend(fontsize=8.5,frameon=False,loc='upper right'); ax.set_ylim(500,2e7)
ax.set_title('Prime patterns: the wheel-corrected prediction lands within a fraction of a percent.\nNumbers above bars: actual ÷ prediction.',fontsize=10.5)
for sp in ['top','right']: ax.spines[sp].set_visible(False)
ax=axs[1]
pp=small[:15]
fac=[(1-2/p)/(1-1/p)**2 if p>2 else 2.0 for p in pp]
cum=np.cumprod(fac)
ax.bar(range(len(pp)),fac,color=['#e76f51' if f>1 else '#2a9d8f' for f in fac])
ax.plot(range(len(pp)),cum,'-o',color='k',ms=4,label='running product')
ax.axhline(1,color='#999',lw=0.8)
ax.set_xticks(range(len(pp))); ax.set_xticklabels([str(p) for p in pp],fontsize=9); ax.set_xlabel('prime p')
ax.set_title('Twin primes: each prime\'s wheel casts a vote.\nBar = (fraction of spots mod p where n and n+2 both survive) ÷ (naive guess).\nThe running product settles at 2 × 0.6602 = 1.3203.',fontsize=10.5)
ax.legend(frameon=False,fontsize=9)
for sp in ['top','right']: ax.spines[sp].set_visible(False)
ax=axs[2]
g=np.diff(P); i=np.arange(1,len(P))
rec=[];best=0
for k in range(len(g)):
    if g[k]>best: best=g[k]; rec.append((P[k],g[k]))
rp=np.array(rec)
ax.loglog(rp[:,0],rp[:,1],'o',color='#7b5ea7',ms=6,label='record gaps (each bigger than all before it)')
xx=np.logspace(1,8,100); ax.loglog(xx,np.log(xx)**2,'--',color='#d62728',label='Cramér: (log p)²')
ax.loglog(xx,np.log(xx),':',color='#555',label='average gap: log p')
ax.set_xlabel('p'); ax.set_ylabel('gap after p'); ax.legend(frameon=False,fontsize=9)
ax.set_title('Gaps: on average about log p, and the record gaps\nhug Cramér\'s (log p)² line from his random model.',fontsize=10.5)
for sp in ['top','right']: ax.spines[sp].set_visible(False)
fig.tight_layout(); fig.savefig('figures/p8_patterns_gaps.png',dpi=150); plt.close(fig)
print(rec[-5:])
