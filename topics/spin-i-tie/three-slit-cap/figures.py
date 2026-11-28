# DRAFT figure for round 12.
import json, numpy as np, matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
d=json.load(open('results.json'))['caps']; k=np.logspace(-5,0.5,300); w=1/3
fig,a=plt.subplots(figsize=(8.2,5))
a.loglog(k,k*w/3,lw=2,label='the law alone runs the odds:  κ = |k|·w/3')
a.loglog(k,(k*w)**2/6,lw=2,label='the law rides on top of quantum:  κ = (k·w)²/6')
for nm,c,ls in (('Sinha 2010',d['Sinha 2010'],'--'),('Kauten 2017',d['Kauten 2017'],'-.')):
    a.axhline(c['kappa'],color='k',ls=ls,lw=1); a.text(1.3e-5,c['kappa']*1.3,f"{nm}: measured κ below {c['kappa']:.0e}",fontsize=9)
ka=d['Kauten 2017']['A_k']; kb=d['Kauten 2017']['B_k']
a.plot([ka,kb],[1e-4,1e-4],'ro',ms=7)
a.annotate(f'|k| < {ka:.0e}',xy=(ka,1e-4),xytext=(ka*0.12,3e-6),arrowprops=dict(arrowstyle='->'),fontsize=9)
a.annotate(f'|k| < {kb:.2f}',xy=(kb,1e-4),xytext=(kb*1.5,1e-6),arrowprops=dict(arrowstyle='->'),fontsize=9)
a.axvspan(kb,k[-1],color='red',alpha=0.08); a.text(0.35,3e-9,'ruled out\neither way',fontsize=9,color='darkred')
a.set_xlabel('|k|, the law\'s dial (each slit carries a third of the flux, w = 1/3)'); a.set_ylabel('three-path term ÷ two-path terms  (κ)')
a.set_ylim(1e-10,1); a.set_title('DRAFT: three-slit experiments cap how far the law\'s k can sit from 0',fontsize=11)
a.legend(fontsize=8.5,loc='lower right'); fig.tight_layout(); fig.savefig('ts1_cap.png',dpi=130)
