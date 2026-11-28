# DRAFT figure for round 15
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
fig,ax=plt.subplots(figsize=(13,4.2)); ax.axis('off'); ax.set_xlim(0,13); ax.set_ylim(0,4.2)
steps=[('law_lift\nlaw_comm, law_assoc','the law and its lift'),('tsallis_indep','Tsallis composes\nby the law, k = 1 - q'),
       ('coin_forces_k_zero','coin rule\nforces k = 0'),('tsallis_breaks_coin_rule','every q ≠ 1 fails'),('Faddeev 1956 (cited)','honest weights\ngive Shannon')]
for i,(name,what) in enumerate(steps):
    x=0.2+i*2.6; col='#ef6c00' if i<4 else '#9e9e9e'
    ax.add_patch(plt.Rectangle((x,1.3),2.1,1.7,fc='white',ec=col,lw=2.5))
    ax.text(x+1.05,2.55,name,ha='center',va='center',fontsize=(9 if len(name)<20 else 7.6),family='monospace',color=col)
    ax.text(x+1.05,1.8,what,ha='center',va='center',fontsize=10)
    if i<4: ax.annotate('',(x+2.55,2.15),(x+2.1,2.15),arrowprops=dict(arrowstyle='->',lw=2))
ax.text(0.2,0.6,'orange: written in CoinRule.lean and checked symbolically (5 of 5), not yet compiled in Lean     grey: cited, not proved here',fontsize=10)
ax.text(0.2,3.6,'DRAFT  round 15: the coin-rule chain as a Lean file',fontsize=12,color='#b71c1c')
fig.savefig('lp1_chain.png',dpi=130,bbox_inches='tight')
