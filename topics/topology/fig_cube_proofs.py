"""Pictures for the AKS identity, Fermat and Lucas on the mult-add cube."""
import numpy as np, matplotlib, itertools
from math import comb
matplotlib.use('Agg')
import matplotlib.pyplot as plt
B='#1d2bd8'; R='#cf3519'; K='#15131a'; M='#b9b7c2'
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':12,'savefig.dpi':160,'savefig.bbox':'tight','savefig.facecolor':'white'})
cols3=['#ffffff',B,R]
def necklace(ax,cx,cy,s,r=0.32,colors=None):
    n=len(s); th=np.pi/2-2*np.pi*np.arange(n)/n
    ax.add_patch(plt.Circle((cx,cy),r,fill=False,color=M,lw=1))
    for t,v in zip(th,s):
        c=(colors or ['#ffffff',B])[v]
        ax.add_patch(plt.Circle((cx+r*np.cos(t),cy+r*np.sin(t)),0.085,facecolor=c,edgecolor=K,lw=0.9,zorder=3))
def orbits(n,k=None,a=2):
    seen=set(); out=[]
    for s in itertools.product(range(a),repeat=n):
        if k is not None and sum(s)!=k: continue
        if s in seen: continue
        orb=[]; t=s
        for _ in range(n):
            if t not in orb: orb.append(t)
            t=t[1:]+t[:1]
        seen.update(orb); out.append(orb)
    return out

# Figure 8: AKS, prime 5 vs composite 4
fig,axs=plt.subplots(1,2,figsize=(13,2.9))
for ax,n,title in [(axs[0],5,'5-cube, layer 2: 10 corners = 2 rotation groups of 5'),(axs[1],4,'4-cube, layer 2: 6 corners = a group of 4 and a group of 2')]:
    orbs=orbits(n,2)
    for row,orb in enumerate(orbs):
        for col,s in enumerate(orb):
            necklace(ax,col*0.95,-row*0.95,s)
        ax.text(len(orb)*0.95-0.35,-row*0.95,f'{len(orb)} corners',va='center',fontsize=11,color=R if len(orb)<n else B)
    ax.set_xlim(-0.5,6.0); ax.set_ylim(-len(orbs)*0.95+0.35,0.5); ax.set_aspect('equal'); ax.axis('off')
    ax.set_title(title,fontsize=12)
plt.savefig('figures/t8_aks_necklaces.png'); plt.close()

# Figure 9: Fermat with p = 3, a = 3 colours
fig,ax=plt.subplots(figsize=(12,4.4))
orbs=orbits(3,None,3)
fixed=[o for o in orbs if len(o)==1]; rest=[o for o in orbs if len(o)==3]
for i,o in enumerate(rest):
    cx=(i%4)*2.9; cy=-(i//4)*1.0
    for j,s in enumerate(o): necklace(ax,cx+j*0.85,cy,s,r=0.3,colors=cols3)
for i,o in enumerate(fixed):
    necklace(ax,12.2,-i*0.75+0.2,o[0],r=0.3,colors=cols3)
ax.text(4.6,0.75,'24 corners in 8 groups of 3',ha='center',color=B,fontsize=12)
ax.text(12.2,0.75,'3 diagonal corners',ha='center',color=R,fontsize=12)
ax.set_xlim(-0.6,13.0); ax.set_ylim(-1.9,1.1); ax.set_aspect('equal'); ax.axis('off')
ax.set_title('Fermat on the 3 × 3 × 3 grid cube: 3³ − 3 = 24 is a multiple of 3',fontsize=12)
plt.savefig('figures/t9_fermat_necklaces.png'); plt.close()

# Figure 10: Lucas, n = 7 = 2*3 + 1, p = 3
n,p=7,3; a,b=divmod(n,p)
ks=np.arange(n+1)
tot=[comb(n,k) for k in ks]
fix=[comb(a,k//p)*comb(b,k%p) for k in ks]
fig,ax=plt.subplots(figsize=(10,4.4))
ax.bar(ks,fix,color=R,width=0.65,label='corners whose blocks are all-on or all-off')
ax.bar(ks,[t-f for t,f in zip(tot,fix)],bottom=fix,color=B,alpha=0.75,width=0.65,label='the rest, moving in groups of 3')
for k,t,f in zip(ks,tot,fix):
    ax.text(k,t+0.8,f'{t} ≡ {f}',ha='center',fontsize=11)
ax.set_xticks(ks); ax.set_xlabel('layer k of the 7-cube')
ax.set_ylim(0,42); ax.legend(frameon=False,loc='upper right')
for s_ in ['top','right']: ax.spines[s_].set_visible(False)
ax.set_title('Lucas: axes in blocks {1,2,3} {4,5,6} {7}; mod 3 only the red corners count: C(2, k div 3)·C(1, k mod 3)',fontsize=11)
plt.savefig('figures/t10_lucas_blocks.png'); plt.close()
print('ok')
