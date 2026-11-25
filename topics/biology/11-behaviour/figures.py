import sys, json, numpy as np
sys.path.insert(0,'../tools')
from style import setup, save, S1, S2, S3, S4, INK, MUTE
import matplotlib.pyplot as plt
setup(); r=json.load(open('results.json'))
fig,ax=plt.subplots(1,2,figsize=(11,4.2))
D=np.linspace(0,30,200); k=0.2
ax[0].plot(D,1/(1+k*D),color=S1,label="κ → 0: hyperbolic 1/(1 + kD)")
for col,kap in zip([S3,S2],[0.3,1.0]):
    c=k*kap   # c/κ = k: the same short-delay slope as the hyperbola
    u=(np.exp(c*D)-1)/kap
    ax[0].plot(D,1/(1+u),color=col,label=f"κ = {kap}" + (": exponential" if kap==1 else ""))
ax[0].set_xlabel("delay D"); ax[0].set_ylabel("discount factor d(D)")
ax[0].set_title("Discounts that compose over split delays by dial κ"); ax[0].legend(fontsize=9)
Vs=r['rw']; ax[1].plot(range(len(Vs)),Vs,'o-',color=S1,ms=4,label="Rescorla–Wagner, αβ = 0.15")
n=np.arange(31); ax[1].plot(n,1-(1-0.15)**n,color=INK,ls='--',lw=1.2,label="λ·wpow(−1, αβ, n)")
ax[1].set_xlabel("trial"); ax[1].set_ylabel("associative strength V/λ")
ax[1].set_title("The learning curve is the 'or' power"); ax[1].legend(fontsize=9)
save(fig,'figures/f1_discount_learning.png')
