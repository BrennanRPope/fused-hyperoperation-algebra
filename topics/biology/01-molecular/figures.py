import sys, math, numpy as np
sys.path.insert(0,"../tools"); from style import *
setup()
from scipy.optimize import brentq
# f1 helix-coil sharpening
fig,ax=plt.subplots(1,2,figsize=(10.5,4.2),gridspec_kw={"width_ratios":[3,2]})
s=np.linspace(0.6,1.6,500)
for sig,col in ((1e-1,S3),(1e-2,S2),(1e-4,S1)):
    lam=lambda s: ((1+s)+np.sqrt((1-s)**2+4*sig*s))/2
    h=1e-6; th=s*(np.log(lam(s+h))-np.log(lam(s-h)))/(2*h)
    ax[0].plot(s,th,color=col,label=f"σ = {sig:g}  →  dial k = 1/σ = {1/sig:g}")
ax[0].set_xlabel("s, helix propagation weight"); ax[0].set_ylabel("helix fraction (long chain)")
ax[0].set_title("Helix–coil: the dial k = 1/σ sets how sharp the switch is"); ax[0].legend(loc="upper left",fontsize=9.5)
ks=np.logspace(-0.5,5.2,200)
ax[1].semilogx(ks,1/(4*np.sqrt(1/ks)),color=S1)
ax[1].axvspan(1e2,1e4,color=S1,alpha=0.12); ax[1].text(1e3,ax[1].get_ylim()[1]*0.85 if False else 50,"measured\npeptides\nk ≈ 10²–10⁴",ha="center",color=INK,fontsize=9.5)
ax[1].set_xlabel("dial k = 1/σ"); ax[1].set_ylabel("steepness at s = 1  (= √k / 4)"); ax[1].set_title("Steepness grows as √k")
save(fig,"figures/f1_helix_coil_dial.png")
# f2 binding
fig,ax=plt.subplots(1,2,figsize=(10.5,4.2))
w=np.logspace(-2,3,300); ax[0].semilogx(w,2*np.sqrt(w)/(1+np.sqrt(w)),color=S1)
for wv,lab in ((0.01,"near-exclusive"),(1,"independent: plain ✳"),(100,"strongly cooperative")):
    ax[0].plot(wv,2*np.sqrt(wv)/(1+np.sqrt(wv)),"o",ms=8,color=S2,mec="white",mew=2)
    ax[0].annotate(lab,(wv,2*np.sqrt(wv)/(1+np.sqrt(wv))),xytext=(8,-4) if wv>0.05 else (8,4),textcoords="offset points",fontsize=9.5,color=INK)
ax[0].set_xlabel("cooperativity ω  (the dial in Z = 1 + a ✳_ω b)"); ax[0].set_ylabel("Hill coefficient at half saturation")
ax[0].set_title("Two sites: the dial sets cooperativity")
L,c,KR,n=9000,0.014,2.8,4
p=np.logspace(-0.3,2.6,400); a=p/KR
Y=(a*(1+a)**3+L*c*a*(1+c*a)**3)/((1+a)**4+L*(1+c*a)**4)
ax[1].semilogx(p,a/(1+a),"--",color=S3,label="R state alone (dial 1)")
ax[1].semilogx(p,c*a/(1+c*a),"--",color=S2,label="T state alone (dial c = 0.014)")
ax[1].semilogx(p,Y,color=S1,label="MWC mix: (1+a)⁴ + L(1+ca)⁴")
ax[1].axvline(26.8,color=MUTE,lw=1); ax[1].text(24,0.62,"human P50 ≈ 27 torr",fontsize=9,color=MUTE,ha="right")
ax[1].set_xlabel("oxygen pressure (torr)"); ax[1].set_ylabel("fraction of sites filled"); ax[1].set_title("Hemoglobin: two mult-add powers mixed")
ax[1].legend(fontsize=9,loc="upper left")
save(fig,"figures/f2_binding_dial.png")
# f3 Eigen threshold
fig,ax=plt.subplots(figsize=(7.5,4.3))
U=np.logspace(-3.5,0.7,300); l=np.exp(-U)-1; sstar=-l/(1+l)
ax.loglog(U,sstar,color=S1,label="threshold s* = mult-add inverse of the copy loss")
ax.fill_between(U,sstar,1e-4,color=S2,alpha=0.12)
ax.text(0.5,2e-3,"master sequence lost\n(error catastrophe)",color=INK,fontsize=10)
for u,lab in ((0.003,"DNA microbes\n(Drake: ~0.003)"),(0.1,"retroviruses\n(~0.1)"),(1.0,"lytic RNA viruses\n(~1)")):
    ax.plot(u,math.exp(u)-1,"o",ms=9,color=S2,mec="white",mew=2); ax.annotate(lab,(u,math.exp(u)-1),xytext=(-10,12),textcoords="offset points",fontsize=9,ha="right")
ax.set_xlabel("mutations per genome per copy, U"); ax.set_ylabel("fitness advantage needed, s")
ax.set_title("Eigen's threshold is the mult-add inverse"); ax.set_ylim(1e-4,30); ax.legend(loc="upper left",fontsize=9.5)
save(fig,"figures/f3_eigen_inverse.png")
print("ok")
