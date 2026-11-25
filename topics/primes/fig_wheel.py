import numpy as np, matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from math import gcd
plt.rcParams.update({'font.size':11})
fig,axs=plt.subplots(1,2,figsize=(15,7.6))
# left: mod 30 wheel
ax=axs[0]; M=30
S=[r for r in range(M) if gcd(r,M)==1]
th=lambda r: np.pi/2-2*np.pi*r/M
ax.add_patch(plt.Circle((0,0),1,fill=False,color='#999',lw=1))
for r in range(M):
    x,y=np.cos(th(r)),np.sin(th(r))
    if r in S:
        ax.plot([0,x],[0,y],color='goldenrod',lw=2,alpha=0.7)
        ax.plot(x,y,'o',ms=15,color='gold',mec='k',zorder=3)
        ax.text(1.13*x,1.13*y,str(r),ha='center',va='center',fontsize=12,weight='bold')
    else:
        c='#1f77b4' if r%2==0 else ('#ff7f0e' if r%3==0 else '#2ca02c')
        ax.plot(x,y,'o',ms=7,color=c,alpha=0.6)
        ax.text(1.13*x,1.13*y,str(r),ha='center',va='center',fontsize=8,color='#999')
for a,b in zip(S,S[1:]+[S[0]+M]):
    m=(a+b)/2; x,y=0.8*np.cos(th(m)),0.8*np.sin(th(m))
    ax.text(x,y,f'+{b-a}',ha='center',va='center',fontsize=9,color='#7b5ea7')
ax.set_title('After 2, 3 and 5 are born, the wheel of 30 has 8 open spokes.\nGaps between them (purple) repeat forever: 6,4,2,4,2,4,6,2',fontsize=11.5)
ax.set_aspect('equal'); ax.axis('off'); ax.set_xlim(-1.3,1.3); ax.set_ylim(-1.3,1.35)
ax.text(0,-1.32,'blue: hit by 2   orange: hit by 3   green: hit by 5',ha='center',fontsize=9,color='#555')
# right: spiral of 7 laps
ax=axs[1]; M2=210
kills=[(7*s)%M2 for s in S]
for n in range(M2):
    lap=n//M; r=0.35+0.11*lap+0.11*(n%M)/M
    t=th(n%M); x,y=r*np.cos(t),r*np.sin(t)
    if gcd(n,M)!=1:
        ax.plot(x,y,'.',ms=3,color='#ccc')
    elif n%7==0:
        ax.plot(x,y,'X',ms=10,color='#d62728',mec='k',mew=0.5,zorder=3)
    else:
        ax.plot(x,y,'o',ms=7,color='gold',mec='k',mew=0.6,zorder=2)
tt=np.linspace(0,7,2000); rr=0.35+0.11*tt
ax.plot(rr*np.cos(np.pi/2-2*np.pi*tt),rr*np.sin(np.pi/2-2*np.pi*tt),color='#ddd',lw=0.8,zorder=0)
n=11; r=0.35+0.11*(n%M)/M; x,y=r*np.cos(th(n)),r*np.sin(th(n))
ax.annotate('11, the first open spot after 1,\nis the next prime to be born',xy=(x,y),xytext=(0.55,1.25),fontsize=10,arrowprops=dict(arrowstyle='->'))
ax.text(0,0,'7 laps\n= 210',ha='center',va='center',fontsize=10,color='#777')
ax.set_title('Adding 7: wind the 30-wheel round 7 times (210 spots).\n7 lands on exactly 8 of the 56 open spots (red X), one per spoke. 48 stay open.',fontsize=11.5)
ax.set_aspect('equal'); ax.axis('off'); lim=1.25; ax.set_xlim(-lim,lim); ax.set_ylim(-lim,lim+0.1)
fig.tight_layout(); fig.savefig('figures/p3_wheel.png',dpi=150); plt.close(fig)
print(S, sorted(kills))

# gap populations across stages
from sympy import primerange
