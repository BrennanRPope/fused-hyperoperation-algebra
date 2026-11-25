import numpy as np, matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from math import comb
B='#1d2bd8'; R='#cf3519'; G='#a86a00'; K='#15131a'; M='#8a8794'
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':12,'axes.edgecolor':K,'axes.labelcolor':K,'xtick.color':K,'ytick.color':K,'savefig.dpi':160,'savefig.bbox':'tight','savefig.facecolor':'white'})
out='figures/'

# T1 helix = universal cover of the circle
fig=plt.figure(figsize=(7,6)); ax=fig.add_subplot(111,projection='3d')
s=np.linspace(-0.2,3.2*2*np.pi,1200)
ax.plot(np.cos(s),np.sin(s),s/(2*np.pi),color=B,lw=2.4,label='helix: the log, unrolled')
c=np.linspace(0,2*np.pi,300); ax.plot(np.cos(c),np.sin(c),-0.9*np.ones_like(c),color=K,lw=2.4,label='Circle-operation circle |1+a| = 1')
for k in range(4):
    th=0.9; ax.scatter([np.cos(th)],[np.sin(th)],[th/(2*np.pi)+k],color=R,s=40,zorder=5)
    ax.plot([np.cos(th)]*2,[np.sin(th)]*2,[th/(2*np.pi)+k,-0.9],color=R,lw=0.8,ls=':')
ax.scatter([np.cos(0.9)],[np.sin(0.9)],[-0.9],color=R,s=60)
ax.text(-0.2,-1.0,-1.35,'one point on the circle',color=R,fontsize=11,ha='center')
ax.text(-0.9,0.3,3.6,'above it, one point per whole turn\n(that ladder is \u03c0\u2081 = \u2124)',color=R,fontsize=11,ha='center')
ax.set_axis_off(); ax.view_init(18,-58); ax.legend(loc='lower right',frameon=False)
plt.savefig(out+'t1_helix_cover.png'); plt.close()

# T2 torus cells = Pascal = Betti
fig,axs=plt.subplots(1,2,figsize=(11,4.6),gridspec_kw={'width_ratios':[1.1,1]})
ax=axs[0]
rows=[[1,1],[1,2,1],[1,3,3,1],[1,4,6,4,1]]
for d,row in enumerate(rows,1):
    for k,v in enumerate(row):
        xx=k-len(row)/2+0.5
        ax.add_patch(plt.Circle((xx,-d),0.36,color=[B,B,G,R,K][k] if False else B,alpha=0.12+0.17*k))
        ax.text(xx,-d,str(v),ha='center',va='center',fontsize=14,color=K)
    ax.text(-3.3,-d,f'T{"¹²³⁴"[d-1]}',ha='left',va='center',fontsize=14,color=B)
for k,lab in enumerate(['points','loops','surfaces','solids','4-solids']):
    pass
ax.text(0,-0.25,'Holes of the n-torus, by dimension\n(computed: 1,1 | 1,2,1 | 1,3,3,1 | 1,4,6,4,1)',ha='center',fontsize=11,color=K)
ax.set_xlim(-3.4,2.8); ax.set_ylim(-4.6,0.3); ax.set_aspect('equal'); ax.axis('off')
ax=axs[1]
u=np.linspace(0,2*np.pi,80); v=np.linspace(0,2*np.pi,40); U,V=np.meshgrid(u,v)
Rr,r=2,0.8; X=(Rr+r*np.cos(V))*np.cos(U); Y=(Rr+r*np.cos(V))*np.sin(U); Z=r*np.sin(V)
ax.remove(); ax=fig.add_subplot(122,projection='3d')
ax.plot_surface(X,Y,Z,color='#c9ccf5',alpha=0.35,linewidth=0)
ax.plot((Rr+r)*np.cos(u),(Rr+r)*np.sin(u),0*u,color=B,lw=2.6,label='loop a (grade 1)')
ax.plot(Rr+r*np.cos(u),0*u,r*np.sin(u),color=R,lw=2.6,label='loop b (grade 1)')
ax.scatter([Rr+r],[0],[0],color=K,s=40,label='the point (grade 0)')
ax.set_axis_off(); ax.view_init(32,-60); ax.set_box_aspect((1,1,0.4))
ax.legend(loc='lower left',frameon=False,fontsize=10)
ax.set_title('T² = circle × circle: 1 point, 2 loops, 1 surface (ab)',fontsize=11)
plt.savefig(out+'t2_torus_pascal.png'); plt.close()

# T3 Chern number by plaquette products
def qwz(kx,ky,m):
    d=np.array([np.sin(kx),np.sin(ky),m+np.cos(kx)+np.cos(ky)])
    H=np.array([[d[2],d[0]-1j*d[1]],[d[0]+1j*d[1],-d[2]]]); return np.linalg.eigh(H)[1][:,0]
def chern(m,N=36):
    ks=np.linspace(0,2*np.pi,N,endpoint=False); u=np.array([[qwz(a,b,m) for b in ks] for a in ks]); F=np.zeros((N,N))
    for i in range(N):
        for j in range(N):
            u1,u2,u3,u4=u[i,j],u[(i+1)%N,j],u[(i+1)%N,(j+1)%N],u[i,(j+1)%N]
            F[i,j]=np.angle(np.vdot(u1,u2)*np.vdot(u2,u3)*np.vdot(u3,u4)*np.vdot(u4,u1))
    return F.sum()/(2*np.pi),F
fig,axs=plt.subplots(1,2,figsize=(11,4.4))
_,F=chern(1.0,48)
im=axs[0].imshow(F.T,origin='lower',extent=[0,2*np.pi,0,2*np.pi],cmap='RdBu_r',vmin=-abs(F).max(),vmax=abs(F).max())
axs[0].set_xlabel('kx'); axs[0].set_ylabel('ky'); axs[0].set_title('Twist in each little square (m = 1)\nadds up to exactly 1 turn',fontsize=11)
ms=np.linspace(-3.5,3.5,57); Cs=[chern(m,24)[0] for m in ms]
axs[1].step(ms,Cs,where='mid',color=B,lw=2.4); axs[1].axhline(0,color=M,lw=0.8)
axs[1].set_xlabel('mass m'); axs[1].set_ylabel('Chern number'); axs[1].set_yticks([-1,0,1])
axs[1].set_title('Whole turns only: the number jumps, never slides',fontsize=11)
for sp_ in ['top','right']: axs[1].spines[sp_].set_visible(False)
plt.savefig(out+'t3_chern.png'); plt.close()

# T4 Hopf invariant one: 2-adic valuation of 3^n - 1 vs n
n=np.arange(1,41); v=[]
for k in n:
    m=3**int(k)-1; c=0
    while m%2==0: m//=2; c+=1
    v.append(c)
fig,ax=plt.subplots(figsize=(9,4.2))
ax.plot(n,n,color=M,lw=1.2,label='n (what is needed)')
ok=[k for k,vv in zip(n,v) if vv>=k]
ax.scatter(n,v,color=B,s=28,label='how many times 2 divides 3ⁿ − 1',zorder=3)
ax.scatter(ok,[v[k-1] for k in ok],color=R,s=90,zorder=4,label='enough: n = 1, 2, 4')
for k in ok: ax.annotate({1:'complex numbers\nS¹',2:'quaternions\nS³',4:'octonions\nS⁷'}[k],(k,v[k-1]),xytext=({1:11,2:13,4:15}[k],{1:7.6,2:10.6,4:13.6}[k]),color=R,fontsize=10,arrowprops=dict(arrowstyle='-',color=R,lw=0.8))
ax.set_xlabel('n'); ax.set_ylim(0,20); ax.legend(frameon=False,loc='upper right')
for sp_ in ['top','right']: ax.spines[sp_].set_visible(False)
ax.set_title('Circle powers ψ² and ψ³ must commute, and that forces 2ⁿ | 3ⁿ − 1',fontsize=11)
plt.savefig(out+'t4_hopf_invariant.png'); plt.close()

# T5 trefoil = (2,3) circle powers on the torus
fig=plt.figure(figsize=(6.4,5.2)); ax=fig.add_subplot(111,projection='3d')
ax.plot_surface(X,Y,Z,color='#c9ccf5',alpha=0.22,linewidth=0)
s=np.linspace(0,2*np.pi,800); p,q=2,3
ph=p*s; th=q*s
ax.plot((Rr+r*np.cos(th))*np.cos(ph),(Rr+r*np.cos(th))*np.sin(ph),r*np.sin(th),color=R,lw=3)
ax.set_axis_off(); ax.view_init(55,-60); ax.set_box_aspect((1,1,0.4))
ax.set_title('Trefoil: z ↦ (z², z³) on the torus',fontsize=12)
plt.savefig(out+'t5_trefoil.png'); plt.close()

# T6 Riemann-Roch Pascal
N=7; KK=7
tab=np.array([[comb(nn+k,nn) for k in range(KK+1)] for nn in range(1,N+1)])
fig,ax=plt.subplots(figsize=(8.6,4.6))
ax.imshow(np.log(tab),cmap='Blues',aspect='auto')
for i in range(N):
    for j in range(KK+1):
        ax.text(j,i,str(tab[i,j]),ha='center',va='center',fontsize=10,color='white' if np.log(tab[i,j])>np.log(tab).max()*0.6 else K)
ax.set_xticks(range(KK+1)); ax.set_yticks(range(N)); ax.set_yticklabels([f'CP{"¹²³⁴⁵⁶⁷"[i]}' for i in range(N)])
ax.set_xlabel('twist k of the line bundle O(k)')
ax.set_title('χ(O(k)) from the circle-operation log (Todd class) = Pascal C(n+k, n)',fontsize=11)
plt.savefig(out+'t6_riemann_roch.png'); plt.close()
print('done')
