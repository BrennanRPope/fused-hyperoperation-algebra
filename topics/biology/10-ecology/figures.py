import sys, json, numpy as np
sys.path.insert(0,'../tools')
from style import setup, save, S1, S2, S3, S4, INK, MUTE
import matplotlib.pyplot as plt
setup()
al,h,ga,tau=0.8,0.3,0.5,0.6
fr=lambda N,P,k:al*N/(1+al*N*h+ga*P*tau+k*al*N*h*ga*P*tau)
fig,ax=plt.subplots(1,2,figsize=(11,4.2))
P=np.linspace(0,12,100)
for c,k,lab in zip([S1,S3,S2,S4],[0,0.5,1,2.5],["k = 0 Beddington–DeAngelis","k = 0.5","k = 1 Crowley–Martin","k = 2.5 (handlers attract interference)"]):
    ax[0].plot(P,fr(16,P,k),color=c,label=lab)
ax[0].set_xlabel("predator density P"); ax[0].set_ylabel("feeding rate per predator")
ax[0].set_title("Prey plentiful (N = 16): the dial sets how fast\ninterference bites"); ax[0].legend(fontsize=8.5)
d=json.load(open('identify.json')); labs=list(d.keys())
for i,l in enumerate(labs):
    lo,hi=d[l]; ax[1].plot([lo,hi],[i,i],color=S1,lw=6,solid_capstyle='butt')
ax[1].axvline(0,color=MUTE,ls='--',lw=1); ax[1].axvline(1,color=MUTE,ls='--',lw=1); ax[1].axvline(0.5,color=INK,lw=1)
ax[1].set_ylim(-0.8,len(labs)-0.5); ax[1].text(0.03,-0.65,"BD",color=MUTE,fontsize=9); ax[1].text(1.03,-0.65,"CM",color=MUTE,fontsize=9)
ax[1].set_yticks(range(len(labs))); ax[1].set_yticklabels(labs,fontsize=9)
ax[1].set_xlim(-0.2,3.3); ax[1].set_xlabel("fitted k, 95% of 200 simulated experiments (true k = 0.5)")
ax[1].set_title("How big an experiment tells BD from CM"); ax[1].grid(axis='y',visible=False)
save(fig,'figures/f1_functional_response.png')
