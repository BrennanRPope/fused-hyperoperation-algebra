import numpy as np, matplotlib, itertools
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from math import comb, factorial, e, log
B='#1d2bd8'; R='#cf3519'; G='#a86a00'; K='#15131a'; M='#8a8794'
LAY=[K,B,G,R,'#6a2bb0','#0b7a5a','#555']
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':12,'axes.edgecolor':K,'axes.labelcolor':K,'xtick.color':K,'ytick.color':K,'savefig.dpi':160,'savefig.bbox':'tight','savefig.facecolor':'white'})
out='figures/'

# E1: the 3-cube, every switch weighs 1/3; layer totals vs 1/k!
fig=plt.figure(figsize=(12,5.2))
ax=fig.add_subplot(121,projection='3d')
corners=list(itertools.product([0,1],repeat=3))
for a,b in itertools.combinations(corners,2):
    if sum(x!=y for x,y in zip(a,b))==1: ax.plot(*zip(a,b),color=M,lw=1.2)
for c in corners:
    k=sum(c); ax.scatter(*c,s=260,color=LAY[k],zorder=5)
    lab={0:'1',1:'1/3',2:'1/9',3:'1/27'}[k]
    ax.text(c[0]+0.07,c[1]+0.07,c[2]+0.09,lab,color=LAY[k],fontsize=12,weight='bold')
ax.set_axis_off(); ax.view_init(20,-62)
ax.set_title('Mult-add 1/3 three times: the 3-cube\nevery corner weighs (1/3)^(switches on)',fontsize=12)
ax2=fig.add_subplot(122)
n=3; ks=np.arange(4)
lay=[comb(n,k)/n**k for k in ks]; inv=[1/factorial(k) for k in ks]
ax2.bar(ks-0.18,lay,0.36,color=[LAY[k] for k in ks],label='layer k total: C(3,k)/3^k')
ax2.bar(ks+0.18,inv,0.36,color='none',edgecolor=K,hatch='//',label='1/k!  (the layer of e)')
for k in ks:
    ax2.text(k-0.18,lay[k]+0.03,['1','1','1/3','1/27'][k],ha='center',fontsize=11,color=LAY[k])
    ax2.text(k+0.18,inv[k]+0.03,['1','1','1/2','1/6'][k],ha='center',fontsize=11,color=K)
ax2.set_xticks(ks); ax2.set_xticklabels([f'k={k}' for k in ks]); ax2.set_ylim(0,1.25)
ax2.set_title(f'Sum of corners = (4/3)³ = {(4/3)**3:.4f}.  Sum of 1/k! for k ≤ 3 = {sum(inv):.4f}',fontsize=11)
ax2.legend(frameon=False,fontsize=10); ax2.spines[['top','right']].set_visible(False)
plt.savefig(out+'e1_cube_weights.png'); plt.close()

# E2: layers converge to 1/k!
fig,ax=plt.subplots(figsize=(10,5))
ks=np.arange(0,7); ns=[2,4,8,32,1000]; cols=['#c6c9f2','#9aa0ec','#6c74e3','#3c47dc',B]
w=0.14
for i,(n,c) in enumerate(zip(ns,cols)):
    vals=[comb(n,k)/n**k if k<=n else 0 for k in ks]
    ax.bar(ks+(i-2)*w,vals,w,color=c,label=f'n = {n}')
ax.plot(ks,[1/factorial(k) for k in ks],'o',color=R,ms=9,label='1/k!',zorder=5)
ax.set_xticks(ks); ax.set_xticklabels([f'layer {k}' for k in ks]); ax.set_ylabel('total weight of the layer')
ax.set_title('Each layer of the n-cube settles on 1/k! as n grows',fontsize=12)
ax.legend(frameon=False,ncol=3); ax.spines[['top','right']].set_visible(False)
plt.savefig(out+'e2_layers_converge.png'); plt.close()

# E3: collisions, k = 2 (grid) and k = 3 (cube)
fig=plt.figure(figsize=(13,4.6))
for i,n in enumerate([4,8,16]):
    ax=fig.add_subplot(1,4,i+1)
    for a in range(n):
        for b in range(n):
            ax.add_patch(plt.Rectangle((a,b),0.92,0.92,color=R if a==b else B,alpha=0.85 if a==b else 0.25))
    ax.set_xlim(0,n); ax.set_ylim(0,n); ax.set_aspect('equal'); ax.axis('off')
    ax.set_title(f'n = {n}: {n} of {n*n} picks collide\n(share {1/n:.3f})',fontsize=11)
ax=fig.add_subplot(1,4,4,projection='3d')
n=6
pts=np.array(list(itertools.product(range(n),repeat=3)))
coll=np.array([len(set(p))<3 for p in pts])
ax.scatter(*pts[~coll].T,color=B,alpha=0.25,s=14)
ax.scatter(*pts[coll].T,color=R,s=18)
ax.set_axis_off(); ax.view_init(22,-55)
ax.set_title(f'k = 3, n = {n}: {coll.sum()} of {n**3} collide\n(at most C(3,2)/n = {3/n:.2f})',fontsize=11)
fig.suptitle('Ordered picks of k axes out of n. Red picks reuse an axis; their share shrinks like C(k,2)/n',fontsize=12,y=1.02)
plt.savefig(out+'e3_collisions.png'); plt.close()

# E4: staircase climbs, trapped by s_N (from below, in the limit) and u_N (from above)
fig,axs=plt.subplots(1,2,figsize=(14,5.2))
sN=np.cumsum([1/factorial(k) for k in range(12)])
uN=[None]+[sN[k]+1/(factorial(k)*k) for k in range(1,12)]
ax=axs[0]; N=np.arange(1,41); st=(1+1/N)**N
ax.step(N,st,where='post',color=B,lw=2.4,label='staircase (1+1/n)ⁿ')
for k in (2,3):
    ax.axhline(sN[k],color=G,lw=1); ax.text(40.5,sN[k],f's{k}',color=G,va='center',fontsize=10)
ax.axhline(uN[2],color=R,lw=1,ls='--'); ax.text(40.5,uN[2],'u2',color=R,va='center',fontsize=10)
ax.axhline(e,color=K,lw=1.6); ax.text(1.5,e+0.012,'e = 2.71828…',color=K)
ax.set_ylim(1.95,2.85); ax.set_xlim(1,42.5); ax.set_xlabel('n (number of mult-add steps)')
ax.set_title('The staircase only climbs',fontsize=12)
ax.legend(loc='lower right',frameon=False); ax.spines[['top','right']].set_visible(False)
ax=axs[1]; N=np.unique(np.logspace(0,5,400).astype(int)); st=(1+1/N)**N
ax.step(N,st,where='post',color=B,lw=2.4)
for k in (4,5):
    ax.axhline(sN[k],color=G,lw=1); ax.text(1.3e5,sN[k],f's{k}',color=G,va='center',fontsize=10)
for k in (4,):
    ax.axhline(uN[k],color=R,lw=1,ls='--'); ax.text(1.3e5,uN[k],f'u{k}',color=R,va='center',fontsize=10)
ax.axhline(e,color=K,lw=1.6)
ax.set_xscale('log'); ax.set_ylim(2.70,2.721); ax.set_xlim(10,2.5e5)
ax.set_xlabel('n (log scale)')
ax.set_title('Zoomed in: under every u N, past every s N in the end',fontsize=12)
ax.spines[['top','right']].set_visible(False)
plt.savefig(out+'e4_staircase.png'); plt.close()

# E5: why e is the natural base
fig,axs=plt.subplots(1,3,figsize=(15,4.8))
ax=axs[0]; x=np.linspace(-1.2,1.3,300)
for b,c,lab in [(2,M,'2ˣ'),(e,B,'eˣ'),(3,G,'3ˣ')]:
    ax.plot(x,b**x,color=c,lw=2.2,label=lab); ax.plot(x,1+log(b)*x,color=c,lw=1,ls='--')
ax.scatter([0],[1],color=K,zorder=5)
ax.set_title('Only eˣ leaves the unit corner at slope 1\n(dashed: tangent at x = 0; slopes ln2, 1, ln3)',fontsize=11)
ax.set_ylim(0,3.6); ax.legend(frameon=False); ax.spines[['top','right']].set_visible(False)
ax=axs[1]; a=np.linspace(-0.6,1.6,300)
ax.plot(a,np.log1p(a),color=B,lw=2.4,label='L(a) = ln(1 + a)')
ax.plot(a,a,color=M,lw=1,ls='--',label='slope 1 line')
for n,c in [(1,R),(2,G),(5,K)]:
    ax.scatter([1/n],[np.log1p(1/n)],color=c,zorder=5)
    ax.text(1/n-0.6 if n==1 else 1/n+0.03,np.log1p(1/n)-0.14,f'n={n}: n·L(1/n) = {n*np.log1p(1/n):.3f}',color=c,fontsize=10)
ax.set_title('L turns ⊕ into +:  L(a ⊕ b) = L(a) + L(b)\nn steps of 1/n weigh n·L(1/n) → 1',fontsize=11)
ax.legend(frameon=False,loc='upper left'); ax.spines[['top','right']].set_visible(False)
ax=axs[2]
for n,c in [(8,G),(16,B),(64,K)]:
    z=[1+0j]
    for _ in range(n): z.append(z[-1]*(1+1j*np.pi/n))
    ax.plot(np.real(z),np.imag(z),'-o',color=c,ms=3,lw=1.4,label=f'n = {n}')
th=np.linspace(0,np.pi,200); ax.plot(np.cos(th),np.sin(th),color=R,lw=1,ls=':')
z=[1+0j]
for _ in range(6): z.append(z[-1]*(1+1/6))
ax.plot(np.real(z),np.zeros(7),'-o',color=B,ms=4,lw=1.4)
ax.annotate('real steps 1/n: out to e',(np.real(z[-1]),0),(1.35,-0.5),arrowprops=dict(arrowstyle='->',color=B),color=B,fontsize=10)
ax.text(-1.9,1.7,'imaginary steps iπ/n: round to −1',color=K,fontsize=10)
ax.set_aspect('equal'); ax.set_xlim(-2.0,2.9); ax.set_ylim(-0.8,1.95)
ax.set_title('Same mult-add, two directions:\n(1+1/n)ⁿ → e  and  (1+iπ/n)ⁿ → −1',fontsize=11)
ax.legend(frameon=False,loc='lower left',fontsize=9); ax.spines[['top','right']].set_visible(False)
plt.savefig(out+'e5_natural_base.png'); plt.close()
print('ok')
