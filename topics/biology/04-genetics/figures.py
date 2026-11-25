import sys, numpy as np
sys.path.insert(0,"../tools"); from style import *
setup()
lens=lambda k,q,a,b:(a+b+k*a*b)/(1-q*a*b)
fig,ax=plt.subplots(1,2,figsize=(11,4.3))
r=np.linspace(0,0.5,300)
for K,c,lab in ((1,S1,"Haldane, no interference (k = −2)"),(0.5,S3,"halfway (K = 0.5)"),(0,S2,"Kosambi (K = 0)")):
    ax[0].plot(r,lens(-2*K,-4+4*K,r,0.2),color=c,label=lab)
ax[0].plot(r,np.minimum(r+0.2,0.5),":",color=MUTE,label="map lengths simply add, capped")
ax[0].axhline(0.5,color=INK,lw=1); ax[0].text(0.01,0.47,"r = ½: unlinked, the wall",fontsize=9.5,color=INK)
ax[0].set_xlabel("recombination fraction r₁ of the first interval"); ax[0].set_ylabel("r over both intervals (r₂ = 0.2)")
ax[0].set_title("Every curve ends at the wall r = ½"); ax[0].legend(fontsize=8.5,loc="lower right")
K=np.linspace(-0.5,1.5,50); ax[1].plot(-2*K,-4+4*K,color=S4)
for (k,q),lab in (((-2,0),"Haldane"),((0,-4),"Kosambi")):
    ax[1].plot(k,q,"o",ms=9,color=S4,mec="white",mew=2); ax[1].annotate(lab,(k,q),xytext=(10,4),textcoords="offset points",fontsize=10,color=INK)
ax[1].plot(-1,2,"x",ms=9,color=MUTE,mew=2); ax[1].annotate("off the line: r = ½ is not\nabsorbing (½ ↦ 0.93)",(-1,2),xytext=(10,-6),textcoords="offset points",fontsize=9,color=MUTE)
ax[1].set_xlabel("dial k"); ax[1].set_ylabel("denominator weight q"); ax[1].set_title("Map-function lens laws lie on 2k + q = −4")
save(fig,"figures/f1_linkage_wall.png")
# f2 G x E
fig,ax=plt.subplots(figsize=(7.2,4.2))
e1,e2=2.7,5.9; ks=np.linspace(-0.3,2.2,100)
ax.plot(ks,1+e1+e2+ks*e1*e2,color=S1)
for k,t in ((0,"additive risks"),(1,"multiplicative risks")):
    ax.plot(k,1+e1+e2+k*e1*e2,"o",ms=8,color=S1,mec="white",mew=2); ax.annotate(t,(k,1+e1+e2+k*e1*e2),xytext=(8,-14),textcoords="offset points",fontsize=9.5,color=INK)
ax.axhline(34.7,color=S2,lw=1.5); ax.text(-0.25,36,"observed: both factors, odds ratio 34.7",fontsize=9.5,color=S2)
ax.axvspan(-0.1,2.2,color=S2,alpha=0.08); ax.text(1.55,4,"bootstrap 95%\nrange runs to 13",fontsize=9,color=MUTE)
ax.set_xlabel("dial k"); ax.set_ylabel("predicted odds ratio with both factors")
ax.set_title("Factor V Leiden with the pill: k ≈ 1.6, poorly pinned")
save(fig,"figures/f2_gene_environment.png")
