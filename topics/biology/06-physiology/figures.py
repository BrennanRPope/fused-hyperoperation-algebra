import sys, numpy as np
sys.path.insert(0,"../tools"); from style import *
setup()
fig,ax=plt.subplots(1,2,figsize=(11,4.3))
for g,c in ((0.1,S3),(0.25,S1),(0.6,S2)):
    F=np.linspace(0,1,200); a=g; b=g; v=b*(1-F)/(F+a)
    ax[0].plot(F,v,color=c,label=f"a/F₀ = {g}")
    fo=a*(np.sqrt(1+1/a)-1); ax[0].plot(fo,b*(1-fo)/(fo+a),"o",ms=7,color=c,mec="white",mew=2)
ax[0].set_xlabel("force F / F₀"); ax[0].set_ylabel("shortening velocity v / v_max")
ax[0].set_title("Hill's muscle law: the level set X ∘₁ Y = X₀"); ax[0].legend(fontsize=9)
ax[0].text(0.45,0.75,"dots: peak power, where X = Y",fontsize=9,color=MUTE)
b1=np.logspace(-1,2,200); b2=10
for om,c,lab in ((0,S1,"same site, ω = 0: ratios add"),(1,S2,"separate sites, ω = 1: ratios multiply"),(3,S4,"sites that help each other, ω = 3")):
    ax[1].loglog(b1+1,1+b1+b2+om*b1*b2,color=c,label=lab)
ax[1].set_xlabel("dose ratio of antagonist 1 alone"); ax[1].set_ylabel("dose ratio of the pair (antagonist 2 alone: 11)")
ax[1].set_title("Two antagonists: the dial reads the mechanism"); ax[1].legend(fontsize=9,loc="upper left")
save(fig,"figures/f1_muscle_antagonists.png")
