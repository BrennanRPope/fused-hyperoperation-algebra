import numpy as np, matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Arc
plt.rcParams.update({'font.size':11,'font.family':'DejaVu Sans'})
def primes_upto(n):
    s=np.ones(n+1,bool); s[:2]=False
    for i in range(2,int(n**.5)+1):
        if s[i]: s[i*i::i]=False
    return np.nonzero(s)[0]
cols={2:'#1f77b4',3:'#ff7f0e',5:'#2ca02c',7:'#d62728'}

# ---------- F1: hops ----------
N=60
fig,ax=plt.subplots(figsize=(14,5.2))
ax.axhline(0,color='#333',lw=1)
hit=set()
for p,c in cols.items():
    h=p*0.32
    for m in range(p, N-p+1, p):
        ax.add_patch(Arc((m+p/2,0),p,2*h*1.0,theta1=0,theta2=180,color=c,lw=1.6,alpha=0.85))
    for m in range(2*p,N+1,p): hit.add(m)
P=set(primes_upto(N))
for n in range(2,N+1):
    if n in hit:
        ax.plot(n,0,'o',ms=5,color='#bbbbbb',zorder=3)
    else:
        ax.plot(n,0,'*',ms=15,color='gold',mec='k',mew=0.8,zorder=4)
    ax.text(n,-0.55,str(n),ha='center',va='top',fontsize=8,color=('k' if n not in hit else '#999'))
for p,c in cols.items():
    yy={2:-1.2,3:-1.75,5:-2.3,7:-2.85}[p]
    ax.annotate(f'{p} is born and starts hopping',xy=(p,-0.75),xytext=(p+0.6,yy),fontsize=8.5,color=c,ha='left',va='center',
                arrowprops=dict(arrowstyle='-',color=c,lw=0.8))
for p,c in cols.items():
    ax.plot([p*p,p*p],[0,p*0.32+0.4],ls=':',color=c,lw=1.2) if p*p<=N else None
ax.annotate('11: the first spot no hop lands on.\nEvery earlier wave missed it, so it is born.',xy=(11,0.1),xytext=(18,3.1),fontsize=10,
            arrowprops=dict(arrowstyle='->',color='k'))
ax.annotate('49 = 7×7: the first number only 7 can catch\n(dotted lines mark p×p for each prime)',xy=(49,2.6),xytext=(40,4.0),fontsize=9,color=cols[7],
            arrowprops=dict(arrowstyle='->',color=cols[7]))
ax.set_xlim(0.5,N+0.8); ax.set_ylim(-3.4,4.8); ax.axis('off')
ax.set_title('The birth of a prime: each prime hops along the line with stride p, landing on its multiples.\nGold stars are where no hop has landed. Each one is a new prime.',fontsize=12.5)
handles=[plt.Line2D([],[],color=c,lw=2,label=f'prime {p}, stride {p}') for p,c in cols.items()]
ax.legend(handles=handles,loc='upper left',frameon=False,ncol=4,fontsize=9)
fig.tight_layout(); fig.savefig('figures/p1_hops.png',dpi=160); plt.close(fig)

# ---------- F2: product of waves ----------
x=np.linspace(1,120.5,40000)
W=np.ones_like(x)
fig,axs=plt.subplots(5,1,figsize=(14,8),sharex=True,gridspec_kw={'height_ratios':[0.8,0.8,0.8,0.8,3]})
for a,(p,c) in zip(axs,cols.items()):
    w=np.abs(np.sin(np.pi*x/p)); W*=w
    a.plot(x,w,color=c,lw=1); a.set_yticks([]); a.set_ylabel(f'|sin(πx/{p})|',rotation=0,ha='right',va='center',fontsize=9)
    a.set_ylim(0,1.15)
    for s in ['top','right','left']: a.spines[s].set_visible(False)
a=axs[-1]
a.fill_between(x,W,color='#7b5ea7',alpha=0.35,lw=0); a.plot(x,W,color='#4b2e83',lw=0.9)
pr=[n for n in primes_upto(121) if n>7]
ints=np.arange(2,121)
Wi=np.prod([np.abs(np.sin(np.pi*ints/p)) for p in cols],axis=0)
for n,v in zip(ints,Wi):
    if v>1e-9:
        a.plot([n,n],[0,v],color='goldenrod',lw=2); a.plot(n,v,'*',color='gold',mec='k',mew=0.6,ms=11)
        a.text(n,v+0.045,str(n),ha='center',fontsize=7.5)
a.set_ylim(0,Wi.max()*1.18); a.set_yticks([]); a.set_ylabel('product\nof all four',rotation=0,ha='right',va='center',fontsize=9)
for s in ['top','right','left']: a.spines[s].set_visible(False)
a.set_xlabel('x')
fig.suptitle('Multiply the four waves together. The product is zero on every multiple of 2, 3, 5 or 7 (the composites).\nThe whole numbers where it stays above zero are exactly the primes from 11 to 113: the "vacuum" the products never reach.',fontsize=12)
fig.tight_layout(rect=[0,0,1,0.93]); fig.savefig('figures/p2_wave_product.png',dpi=150); plt.close(fig)
print('nonzero ints', [int(n) for n,v in zip(ints,Wi) if v>1e-9])
