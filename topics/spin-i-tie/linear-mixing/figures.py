# DRAFT figure for round 8.
import json, numpy as np, matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch
d=json.load(open('results.json')); q=np.array(d['q_grid']); f=lambda v: np.maximum(np.array(v),1e-16)
fig,ax=plt.subplots(1,2,figsize=(12,4.8),gridspec_kw={'width_ratios':[1,1.15]})
a=ax[0]
a.semilogy(q,f(d['global_tsallis']),'o-',ms=3,label='Tsallis (the law at k = 1 − q), coin on everything')
a.semilogy(q,f(d['local_tsallis']),'s-',ms=3,label='Tsallis, coin on one outcome')
a.semilogy(q,f(d['local_renyi']),'^-',ms=3,label='Rényi, coin on one outcome')
a.semilogy(q,f(d['global_renyi']),'k-',lw=2,label='Rényi, coin on everything: never changes')
a.axvline(1,color='gray',ls=':'); a.text(1.05,1e-10,'Shannon\n(k = 0)',fontsize=9)
a.set_ylim(1e-17,1e7); a.set_xlabel('order q   (dial k = 1 − q)')
a.set_ylabel('change in what is left to learn (log)')
a.set_title('Does an irrelevant coin flip change what is left to learn?',fontsize=10); a.legend(fontsize=7.5,loc='upper right')
b=ax[1]; b.axis('off'); b.set_xlim(0,10); b.set_ylim(0,10)
boxes=[(5,9.0,'Mixtures add linearly\n(quantum mechanics already assumes this)','#cfe3f7'),
       (5,7.0,'An irrelevant coin flip changes nothing\n(assumed: our reading of linear mixing)','#cfe3f7'),
       (5,5.0,'Information dial k = 0, and the entropy is Shannon\n(this round, 6 of 6: coin on everything gives k = 0; coin on one outcome gives Shannon)','#d7f0d2'),
       (5,3.0,'Shannon chain rule\n(rounds 5 and 7)','#fde3c8'),
       (5,1.0,'Information causality gives 2√2\n(round 3; "no more than the bits sent" is still assumed)','#fde3c8')]
for x,y,t,c in boxes:
    b.add_patch(FancyBboxPatch((x-4.4,y-0.7),8.8,1.4,boxstyle='round,pad=0.05',fc=c,ec='gray'))
    b.text(x,y,t,ha='center',va='center',fontsize=8.5)
for y in (8.3,6.3,4.3,2.3): b.annotate('',xy=(5,y-0.55),xytext=(5,y),arrowprops=dict(arrowstyle='->',lw=1.5))
b.set_title('Where the dial gets pinned (blue: assumed, green: this round, orange: earlier)',fontsize=9.5)
fig.suptitle('DRAFT: linear mixing pins the information dial to k = 0',fontsize=11)
fig.tight_layout(); fig.savefig('lm1_coin.png',dpi=130)
