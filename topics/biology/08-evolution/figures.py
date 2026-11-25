import sys, json, numpy as np
sys.path.insert(0,'../tools')
from style import setup, save, S1, S2, INK, MUTE
import matplotlib.pyplot as plt
setup(); r=json.load(open('results.json')); A=np.array(r['sim']['A']); k=r['sim']['k']; wall=-1/k
xs=np.linspace(-0.15,0.75,100)
fig,ax=plt.subplots(1,2,figsize=(11,4.2),sharey=True)
for a in A:
    c=S1 if a>0 else S2
    ax[0].plot(xs,a+k*a*xs,color=c)
    ax[1].plot(xs,a-0.12*xs,color=c)
for a0 in ax: a0.axhline(0,color=MUTE,lw=1); a0.set_xlabel("background fitness (excess over ancestor)")
ax[0].axvline(wall,color=INK,ls='--',lw=1); ax[0].text(wall+0.01,0.12,"wall −1/k:\nevery mutation\nneutral",fontsize=9,color=INK)
ax[0].set_ylabel("fitness effect of the mutation"); ax[0].set_title("Path-consistent: slope ∝ effect (the dial)")
ax[1].set_title("Common slope: order changes the fitness")
ax[0].plot([],[],color=S1,label="beneficial"); ax[0].plot([],[],color=S2,label="deleterious"); ax[0].legend(fontsize=9,loc='lower left')
save(fig,'figures/f1_global_epistasis.png')
