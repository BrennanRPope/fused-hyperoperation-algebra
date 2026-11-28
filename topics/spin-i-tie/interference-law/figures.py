# DRAFT figure for round 11.
import json, numpy as np, matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
d=json.load(open('results.json'))
fig,ax=plt.subplots(1,3,figsize=(15,4.6))
a=ax[0]; p,q=0.3,0.5; ph=np.linspace(0,2*np.pi,300)
a.plot(ph,p+q+2*np.sqrt(p*q)*np.cos(ph),lw=2,label='quantum: p + q + 2√(pq)·cos φ')
for k,c in [(0,'k'),(-1,'C2'),(2,'C3')]:
    a.axhline(p+q+k*p*q,color=c,ls='--',lw=1.2,label=f'the law at k = {k}: flat in φ')
a.set_xlabel('phase difference φ'); a.set_ylabel('odds for "either path"'); a.set_xticks([0,np.pi,2*np.pi]); a.set_xticklabels(['0','π','2π'])
a.set_title('Two paths: the dial is an angle, not a k',fontsize=10); a.legend(fontsize=7.5,loc='upper center')
b=ax[1]; names=list(d['bars']); x=np.arange(4); wdt=0.27
for i,nm in enumerate(names):
    b.bar(x+(i-1)*wdt,np.maximum(d['bars'][nm],1e-6),wdt,label=nm)
b.set_yscale('log'); b.set_ylim(1e-4,1); b.set_xticks(x); b.set_xticklabels(['1 path','2 paths','3 paths','4 paths'])
b.set_ylabel('size of the interference term (log; no bar = exactly 0)'); b.legend(fontsize=8)
b.set_title('Quantum stops at 2 paths; the law at k ≠ 0 never stops',fontsize=10)
c=ax[2]; g=np.array(d['g'])
c.plot(g,d['I2g'],'k-',lw=1.5,label='expected: 2√(pq) × coherence'); c.plot(g,d['I2mc'],'o',ms=4,label='simulated phase noise')
c.annotate('fully averaged:\nthe k = 0 law (plain adding)',xy=(0,0),xytext=(0.25,0.12),arrowprops=dict(arrowstyle='->'),fontsize=9)
c.set_xlabel('coherence left after averaging over phases'); c.set_ylabel('two-path interference term')
c.set_title('Aggregation wipes out interference and lands on k = 0',fontsize=10); c.legend(fontsize=8,loc='upper left')
fig.suptitle('DRAFT: interference is not the law in disguise; averaging over phases turns it into the k = 0 law',fontsize=11)
fig.tight_layout(); fig.savefig('il1_orders.png',dpi=125)
