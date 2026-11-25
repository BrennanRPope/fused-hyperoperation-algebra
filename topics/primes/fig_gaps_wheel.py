import numpy as np, matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
stages=[3,5,7,11,13,17,19]; P=[2,3,5,7,11,13,17,19]
res={}
M=1; surv=np.array([1])  # residues mod M coprime
M=2; surv=np.array([1])
for p in P[1:]:
    # lift to M*p
    allr=(surv[None,:]+M*np.arange(p)[:,None]).ravel()
    allr=np.sort(allr[allr%p!=0]); M*=p; surv=allr
    g=np.diff(np.concatenate([surv,[surv[0]+M]]))
    vals,cnt=np.unique(g,return_counts=True)
    res[p]=(M,len(surv),dict(zip(vals.tolist(),cnt.tolist())))
    print(p,M,len(surv),dict(list(zip(vals.tolist(),cnt.tolist()))[:8]), 'max gap',g.max())
fig,axs=plt.subplots(1,2,figsize=(15,6))
ax=axs[0]
gs=[2,4,6,8,10,12,14,18]
cmap=plt.get_cmap('tab10')
xs=list(range(len(res)))
for i,g in enumerate(gs):
    y=[res[p][2].get(g,0)/res[p][1] for p in res]
    ax.plot(xs,y,'-o',color=cmap(i%10),lw=2 if g in(2,6) else 1.2,label=f'gap {g}')
ax.set_xticks(xs); ax.set_xticklabels([f'after {p}\n(mod {res[p][0]:,})' for p in res],fontsize=8.5)
ax.set_ylabel('share of all gaps in the wheel')
ax.legend(ncol=3,fontsize=9,frameon=False)
ax.set_title('Each new prime reshapes the wheel\'s gap pattern.\nGaps 2 and 4 stay exactly tied at every stage; 6 overtakes them once 11 is born.',fontsize=11)
for s in ['top','right']: ax.spines[s].set_visible(False)
ax=axs[1]
p=19; d=res[p][2]; ks=sorted(d)
ax.bar(ks,[d[k] for k in ks],width=1.6,color=['#d62728' if k==6 else '#7b5ea7' for k in ks])
ax.set_yscale('log'); ax.set_xlabel('gap size'); ax.set_ylabel('how many times it appears (log scale)')
ax.set_title(f'The full gap pattern of the wheel after 19 (mod {res[19][0]:,}, {res[19][1]:,} open spokes).\nOnly even gaps appear. 6 (red) is now the most common. The longest gap is {max(ks)}.',fontsize=11)
for s in ['top','right']: ax.spines[s].set_visible(False)
fig.tight_layout(); fig.savefig('figures/p4_wheel_gaps.png',dpi=150); plt.close(fig)
