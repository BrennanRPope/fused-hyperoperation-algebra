import sys, numpy as np
sys.path.insert(0,"../tools"); from style import *
setup()
cols=(S1,S2,S3)
# f1 double-reciprocal lines on the dial
fig,ax=plt.subplots(1,3,figsize=(12,4.2),sharey=False)
for a,(k,title) in zip(ax,((0,"Ping-pong: dial k = 0"),(1/6,"Malate dehydrogenase: k = 1/6"),(1,"Random, independent: k = 1"))):
    x=np.linspace(-1.6 if k==0 else (-1/k-0.6 if k>0.5 else -7.5),3,300)
    for uB,c in zip((0.25,1,4),cols):
        a.plot(x,1+x+uB+k*x*uB,color=c,label=f"u_B = {uB:g}")
    if k>0:
        a.plot(-1/k,1-1/k,"o",ms=9,color=INK,mec="white",mew=2)
        a.annotate("wall u_A = −1/k",(-1/k,1-1/k),xytext=(10,-14),textcoords="offset points",fontsize=9.5,color=INK)
    else:
        a.text(-1.5,7.2,"parallel: the wall\nis at infinity",fontsize=9.5,color=INK)
    a.axhline(0,color=MUTE,lw=0.8); a.axvline(0,color=MUTE,lw=0.8)
    a.set_xlabel("u_A = K_a / A   (reciprocal of first substrate)"); a.set_title(title)
ax[0].set_ylabel("V / v   (reciprocal rate)"); ax[0].legend(fontsize=9,loc="upper left",bbox_to_anchor=(0,0.82))
save(fig,"figures/f1_lines_meet_at_wall.png")
# f2 two-nutrient growth: dial 1, dial 0, and Liebig
fig,ax=plt.subplots(1,2,figsize=(10.5,4.2))
u=np.logspace(-2,2,300)
for k,c,lab in ((1,S1,"multiplicative Monod (k = 1)"),(0,S2,"additive resistance (k = 0)")):
    ax[0].loglog(u,1/(1+u+u+k*u*u),color=c,label=lab)
ax[0].loglog(u,1/(1+u),"--",color=S3,label="Liebig minimum (not on the dial)")
ax[0].set_xlabel("u = K / S, both nutrients equally scarce"); ax[0].set_ylabel("growth rate / maximum")
ax[0].set_title("Two scarce nutrients: three laws"); ax[0].legend(fontsize=9,loc="lower left")
ks=np.linspace(0,1.4,200); g=np.logspace(-2,2,400)
miss=[np.max(np.abs(1/(1+2*g+k*g*g)-1/(1+g))) for k in ks]
ax[1].plot(ks,miss,color=S4); ib=int(np.argmin(miss))
ax[1].plot(ks[ib],miss[ib],"o",ms=8,color=S4,mec="white",mew=2)
ax[1].annotate(f"best k ≈ {ks[ib]:.2f}, still {miss[ib]*100:.0f}% off",(ks[ib],miss[ib]),xytext=(25,30),textcoords="offset points",fontsize=9.5,color=INK)
ax[1].set_xlabel("dial k   (k < 0 blows up: growth would hit a pole)"); ax[1].set_ylabel("worst gap to Liebig's minimum"); ax[1].set_title("No dial setting reproduces Liebig")
save(fig,"figures/f2_nutrient_dial.png")
