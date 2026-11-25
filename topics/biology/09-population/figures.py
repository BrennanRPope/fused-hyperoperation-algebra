import sys, json, numpy as np
sys.path.insert(0,'../tools')
from style import setup, save, S1, S2, S3, S4, INK, MUTE
import matplotlib.pyplot as plt
setup(); r=json.load(open('results.json'))
fig,ax=plt.subplots(1,2,figsize=(11,4.2))
N=2000; R0=r['rf']['R0']
for c,p in zip([S1,S2,S3],[0.0008,0.0012,0.002]):
    S=N-1.0; C=1.0; cum=[1.0]; t=0
    while C>1e-6 and t<60:
        newC=S*(1-(1-p)**C); S-=newC; C=newC; cum.append(cum[-1]+C); t+=1
    ax[0].plot(np.array(cum)/N,color=c,marker='o',ms=3,label=f"R₀ = {p*N:.1f}")
ax[0].set_xlabel("generation of infection"); ax[0].set_ylabel("fraction ever infected")
ax[0].set_title("Reed–Frost: each step is wpow(−1, p, I)"); ax[0].legend(fontsize=9)
Ns=np.linspace(0,250,300); lam=3.0; cc=0.02
for col,b in zip([S1,S2,S3,S4],[1,2,5,20]):
    ax[1].plot(Ns,lam*Ns/(1+cc*Ns/b)**b,color=col,label=f"b = {b}" + (" (Beverton–Holt)" if b==1 else ""))
ax[1].plot(Ns,lam*Ns*np.exp(-cc*Ns),color=INK,ls='--',lw=1.2,label="Ricker (b → ∞)")
ax[1].set_xlabel("population this year N"); ax[1].set_ylabel("next year N'")
ax[1].set_title("Hassell: denominator 1 + mpow(aN, b)"); ax[1].legend(fontsize=9)
save(fig,'figures/f1_epidemic_density.png')
