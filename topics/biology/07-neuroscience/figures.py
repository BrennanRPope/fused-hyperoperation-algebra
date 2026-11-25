import sys, json, numpy as np
sys.path.insert(0,'../tools')
from style import setup, save, S1, S2, S3, S4, INK, MUTE
import matplotlib.pyplot as plt
from dendrite import build, respond, lens_constants
setup(); r=json.load(open('results.json'))
G=build(); EE,EI=65.0,-5.0; iE,iI=40,10
k,q=r['hao']['k'],r['hao']['q']
fig,ax=plt.subplots(1,2,figsize=(11,4.3))
# left: shunting term against EPSP*IPSP
a0=ax[0]
cols=[S1,S2,S3,S4]
for c,gI in zip(cols,[0.2,0.6,1.2,2.0]):
    b=respond(G,[(iI,gI,EI)]); xs=[];ys=[]
    for gE in np.linspace(0.005,0.35,25):
        a=respond(G,[(iE,gE,EE)]); ab=respond(G,[(iE,gE,EE),(iI,gI,EI)])
        xs.append(a*b); ys.append(ab-a-b)
    a0.plot(xs,ys,'o',ms=4,color=c,label=f"IPSP {b:.1f} mV")
xx=np.linspace(-15,0,50); a0.plot(xx,k*xx,color=INK,lw=1.2,ls='--',label=f"dial: slope k = {k:.3f}/mV")
a0.set_xlabel("EPSP × IPSP (mV²)"); a0.set_ylabel("together − (EPSP + IPSP)  (mV)")
a0.set_title("The shunting term is the dial's k·ab"); a0.legend(fontsize=9,loc='upper left')
# right: k against location of excitation
a1=ax[1]
ks=np.array(r['location']['ks']); a1.plot(np.arange(1,61),ks,color=S1,label="inhibition on trunk (node 10)")
a1.axvline(iI,color=S1,lw=1,ls=':')
ktr=r['oblique']['ktr']; kin=r['oblique']['kin']
a1.plot(np.arange(1,61),ktr,color=S2,label="inhibition on oblique branch: E on trunk")
a1.plot(np.arange(61,86),kin,color=S2,ls='--',label="inhibition on oblique branch: E on branch")
a1.axhline(r['location']['kdist'],color=MUTE,lw=1,ls='--')
a1.text(62,r['location']['kdist']+0.015,"Z_II/(|E_I|·Z_sI)",color=MUTE,fontsize=9)
a1.set_xlabel("where the excitation sits (node; 1–60 trunk, 61–85 oblique)"); a1.set_ylabel("k (per mV)")
a1.set_title("k is set by geometry, not by synapse strength"); a1.legend(fontsize=8.5,loc='upper left')
save(fig,'figures/f1_dendrite_dial.png')
