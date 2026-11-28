# DRAFT figure for round 6.
import json, math, matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
d=json.load(open('results.json')); n=d['n']
fig,ax=plt.subplots(figsize=(7.5,4.6))
lab={'0.5':'Rényi 0.5','1':'Shannon','2':'Rényi 2','3':'Rényi 3','inf':'Rényi ∞'}
for k,v in d['thresholds'].items(): ax.plot(n,v,'o-',ms=4,label=lab[k])
ax.axhline(1/math.sqrt(2),color='k',ls='--',lw=1); ax.text(1.1,0.713,'quantum boxes 0.707',ha='left',fontsize=9)
ax.annotate('Rényi 2 fails at the first level: 0.644\n(Oughton & Timpson 2024)',xy=(1,0.6436),xytext=(3.2,0.56),arrowprops=dict(arrowstyle='->'),fontsize=9)
ax.set_xlabel('tree depth n'); ax.set_ylabel('largest bias E that information causality allows')
ax.set_title('DRAFT: information causality with a Rényi entropy, level by level',fontsize=11)
ax.set_ylim(0.38,0.95); ax.legend(fontsize=8,loc='upper right'); fig.tight_layout(); fig.savefig('rb1_depth.png',dpi=130)
