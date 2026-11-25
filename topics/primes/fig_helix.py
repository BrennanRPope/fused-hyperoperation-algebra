import numpy as np, matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
plt.rcParams.update({'font.size':10})
G=np.load('zeta_zeros_1000.npy')
cols={2:'#1f77b4',3:'#ff7f0e',5:'#2ca02c',7:'#d62728'}
fig=plt.figure(figsize=(16,7))
ax=fig.add_subplot(1,2,1,projection='3d')
X=40; x=np.linspace(0,X,6000)
for p,c in cols.items():
    r=0.4+0.18*list(cols).index(p)
    y=r*np.sin(2*np.pi*x/p); z=-r*np.cos(2*np.pi*x/p)
    ax.plot(x,y,z,color=c,lw=1.3,label=f'prime {p}: one turn every {p}')
    m=np.arange(p,X+1,p); ax.scatter(m,0*m,-r+0*m,color=c,s=14,depthshade=False)
from sympy import isprime
for n in range(2,X+1):
    if all(n%p for p in cols) :
        ax.scatter([n],[0],[-1.3],marker='*',s=120,color='gold',edgecolor='k',depthshade=False)
        ax.text(n,0,-1.55,str(n),fontsize=8,ha='center')
ax.plot([0,X],[0,0],[-1.3,-1.3],color='#999',lw=0.8)
ax.set_title('The sieve as helices: each prime coils around the number line once every p steps,\ntouching bottom (dots) on its multiples. Gold stars: no coil touches bottom there.',fontsize=10.5)
ax.view_init(elev=14,azim=-62); ax.set_yticks([]); ax.set_zticks([]); ax.set_xlabel('x')
ax.legend(loc='upper left',fontsize=8,frameon=False)
ax=fig.add_subplot(1,2,2,projection='3d')
t=np.linspace(0,np.log(30),6000)
for g,c in zip(G[:3],['#264653','#e76f51','#2a9d8f']):
    rho=0.5+1j*g; w=np.exp(rho*t)/abs(rho)
    ax.plot(t,w.real,w.imag,color=c,lw=1.2,label=f'zero γ = {g:.3f}: x^ρ/ρ, ρ = ½ + {g:.2f}i')
ax.set_xlabel('log x'); ax.set_yticks([]); ax.set_zticks([])
ax.view_init(elev=14,azim=-62)
ax.legend(loc='upper left',fontsize=8,frameon=False)
ax.set_title('The zeros as helices: each zeta zero is a cone-shaped coil along log x,\nwidening like √x and turning at rate γ. Their shadows, added up, are the waves in the music picture.',fontsize=10.5)
fig.tight_layout(); fig.savefig('figures/p9_helices.png',dpi=150); plt.close(fig)
