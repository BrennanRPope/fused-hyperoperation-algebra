import numpy as np, matplotlib, itertools
matplotlib.use('Agg')
import matplotlib.pyplot as plt
plt.rcParams.update({'font.size':11})
fig=plt.figure(figsize=(17,6.7))
# panel 1: cube 2,3,5 in oblique projection
ax=fig.add_subplot(1,3,1)
proj=lambda a,b,c: (a*1.0+c*0.45, b*1.0+c*0.38)
pr=[2,3,5]
corners=list(itertools.product([0,1],repeat=3))
for u,v in itertools.combinations(corners,2):
    if sum(abs(np.array(u)-np.array(v)))==1:
        x1,y1=proj(*u); x2,y2=proj(*v); ax.plot([x1,x2],[y1,y2],color='#888',lw=1.3,zorder=1)
for c in corners:
    n=int(np.prod([p**e for p,e in zip(pr,c)])); k=sum(c)
    x,y=proj(*c); col='#2a9d8f' if k%2==0 else '#e76f51'
    ax.plot(x,y,'o',ms=34,color=col,mec='k',zorder=2)
    ax.text(x,y+0.02,f'{n}',ha='center',va='center',fontsize=12,weight='bold',color='w',zorder=3)
    ax.text(x+0.13,y-0.15,'+' if k%2==0 else '−',fontsize=15,weight='bold',color=col,zorder=3)
ax.annotate('×2',xy=(0.5,-0.1),ha='center',fontsize=11,color='#555')
ax.annotate('×3',xy=(-0.14,0.5),ha='center',fontsize=11,color='#555')
ax.annotate('×5',xy=(0.12,0.27),fontsize=11,color='#555')
ax.set_xlim(-0.35,1.75); ax.set_ylim(-0.35,1.65); ax.set_aspect('equal'); ax.axis('off')
ax.set_title('Three primes combined by the circle operation: the cube.\nEach corner is a combination; layers 1, 3, 3, 1.\nGreen = even layer (+), red = odd layer (−): that sign is μ(n).',fontsize=10.5)
# panel 2: 4-cube as layers (Hasse diagram)
ax=fig.add_subplot(1,3,2)
pr4=[2,3,5,7]; subs=list(itertools.product([0,1],repeat=4)); pos={}
for k in range(5):
    lay=[s for s in subs if sum(s)==k]
    for i,s in enumerate(lay):
        pos[s]=(k, (i-(len(lay)-1)/2)*0.9)
for u,v in itertools.combinations(subs,2):
    d=np.array(v)-np.array(u)
    if abs(d).sum()==1:
        ax.plot([pos[u][0],pos[v][0]],[pos[u][1],pos[v][1]],color='#ccc',lw=0.8,zorder=1)
for s,(x,y) in pos.items():
    n=int(np.prod([p**e for p,e in zip(pr4,s)])); k=sum(s)
    col='#2a9d8f' if k%2==0 else '#e76f51'
    ax.plot(x,y,'o',ms=26,color=col,mec='k',zorder=2); ax.text(x,y,str(n),ha='center',va='center',fontsize=8.5,color='w',weight='bold',zorder=3)
for k,cnt in enumerate([1,4,6,4,1]):
    ax.text(k,-3.0,f'layer {k}\n{cnt} corner{"s" if cnt>1 else ""}\nsign {"+" if k%2==0 else "−"}',ha='center',fontsize=9)
ax.set_xlim(-0.5,4.5); ax.set_ylim(-3.4,3.0); ax.axis('off')
ax.set_title('Four primes: the 4-cube, laid out layer by layer.\n16 corners = every squarefree number built from 2, 3, 5, 7.',fontsize=10.5)
# panel 3: lattice
ax=fig.add_subplot(1,3,3)
K=5
for a in range(K):
    for b in range(K):
        n=2**a*3**b; inside=(a<=1 and b<=1)
        if inside:
            col='#2a9d8f' if (a+b)%2==0 else '#e76f51'
            ax.plot(a,b,'s',ms=30,color=col,mec='k',zorder=2); ax.text(a,b,str(n),ha='center',va='center',color='w',weight='bold',fontsize=10,zorder=3)
        else:
            ax.plot(a,b,'s',ms=30,color='#eee',mec='#bbb',zorder=2); ax.text(a,b,str(n),ha='center',va='center',fontsize=8 if n<1000 else 6.5,color='#666',zorder=3)
ax.add_patch(plt.Rectangle((-0.45,-0.45),1.9,1.9,fill=False,lw=2.5,color='k'))
ax.set_xticks(range(K)); ax.set_yticks(range(K))
ax.set_xlabel('power of 2'); ax.set_ylabel('power of 3')
ax.set_xlim(-0.7,K-0.3); ax.set_ylim(-0.7,K-0.3); ax.set_aspect('equal')
for s in ['top','right']: ax.spines[s].set_visible(False)
ax.set_title('Euler\'s ζ(s) sums 1/nˢ over the whole grid of prime powers.\n1/ζ(s) lives only on the unit square (cube) with ± signs.\nOff the cube μ(n) = 0.',fontsize=10.5)
fig.tight_layout(); fig.savefig('figures/p5_cube.png',dpi=150); plt.close(fig)
