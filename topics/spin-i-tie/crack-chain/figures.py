# DRAFT figures for round 7.
import json, numpy as np, matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
d=json.load(open('results.json')); q=np.array(d['q_grid'])
fig,ax=plt.subplots(1,2,figsize=(11,4.4))
a=ax[0]
f=lambda v: np.maximum(np.array(v),1e-16)
a.semilogy(q,f(d['resid_add_honest']),'o-',ms=3,label='plain sum, honest average (Shannon form)')
a.semilogy(q,f(d['resid_law_honest']),'s-',ms=3,label='the law, honest average')
a.semilogy(q,np.full_like(q,1e-16),'k-',lw=2,label='plain sum, weights p$^q$: exact (rounding floor)')
a.axvline(1,color='gray',ls=':'); a.text(1.05,1e-9,'k = 0\n(Shannon)',fontsize=9)
a.set_xlabel('q   (dial k = 1 − q)'); a.set_ylabel('worst chain-rule error (log), 40 random joints')
a.set_title('Which chain rules hold at each dial',fontsize=10); a.legend(fontsize=7.5,loc='upper right'); a.set_ylim(1e-17,1e3)
b=ax[1]
rng=np.random.default_rng(3); p=rng.random(5); p/=p.sum()
qq=np.linspace(0.2,3,100)
S=lambda q: -(p*np.log(p)).sum() if abs(q-1)<1e-9 else (1-(p**q).sum())/(q-1)
b.plot(qq,[(p**x).sum() for x in qq],lw=2,label='total weight  Σ p$^q$')
b.plot(qq,[1+(1-x)*S(x) for x in qq],'--',lw=2,label="the law's lift  1 + k·S")
b.axhline(1,color='gray',ls=':'); b.axvline(1,color='gray',ls=':')
b.annotate('weights are honest\nprobabilities only here',xy=(1,1),xytext=(1.6,1.8),arrowprops=dict(arrowstyle='->'),fontsize=9)
b.set_xlabel('q   (dial k = 1 − q)'); b.set_ylabel('sum of conditional weights')
b.set_title("The weights add up to the law's lift",fontsize=10); b.legend(fontsize=8)
fig.suptitle('DRAFT: the chain rule holds at every dial; only k = 0 averages honestly',fontsize=11)
fig.tight_layout(); fig.savefig('cc1_weights.png',dpi=130)
