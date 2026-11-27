import json, numpy as np, matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch
from fractions import Fraction as F
import union_ifs as U
plt.rcParams.update({"font.size":11})
GREEN="#2e9e5b"; BLUE="#2f6db5"; ORANGE="#e08a1e"; GREY="#999999"; RED="#c0392b"
res=json.load(open("union_ifs.json"))

# ---- Fig 1: tie map
fig,ax=plt.subplots(figsize=(13,8)); ax.axis("off")
N={ # name: (x,y,colour,label)
 "setops":(0.12,0.80,BLUE,"Union at the core\nk = -1/(box size)\n(set operations)"),
 "order":(0.12,0.50,ORANGE,"Maps commute exactly\nwhen they share a centre\n(order dial)"),
 "trace":(0.40,0.66,BLUE,"Entropy from shared centres\n1 + sum[(1-t)^m - 1]\n(chaos, Cartier-Foata)"),
 "union":(0.70,0.66,GREEN,"BUILT 1: union fractal\ndimension from the\nshared-centre rule"),
 "zoom":(0.92,0.86,ORANGE,"Zoom = dimension D\n(dial audit: r's own number)"),
 "moran":(0.92,0.46,ORANGE,"Moran equation\nsum s^D = 1 (1946)"),
 "carnot":(0.40,0.24,BLUE,"Low-dissipation engines on\nthe parabola d = k^2/4 - 1\n(Carnot cascade)"),
 "curv":(0.12,0.18,BLUE,"Curvature -Delta/4\n(derived-dial catalogue)"),
 "speed":(0.70,0.24,GREEN,"BUILT 2: discriminant 4\n= every engine closes its\ngap at Carnot's rate"),
 "fuzzy":(0.40,0.95,BLUE,"Fuzzy monogamy\n(holographic cone)"),
 "pair":(0.70,0.95,BLUE,"Pair dial and triple overlap\n= I3 (set operations)"),
 "qsw":(0.40,0.43,BLUE,"q switch: no partial chaos\nfrom a smooth dial"),
}
E=[("setops","union",GREEN),("order","trace",GREY),("trace","union",GREEN),("order","union",GREEN),("union","zoom",GREEN),
   ("union","moran",GREEN),("carnot","speed",GREEN),("curv","speed",GREEN),("qsw","trace",GREY),("fuzzy","pair",BLUE)]
for a,b,c in E:
    (x1,y1,_,_),(x2,y2,_,_)=N[a],N[b]
    ax.annotate("",xy=(x2,y2),xytext=(x1,y1),arrowprops=dict(arrowstyle="-|>",color=c,lw=2.5 if c==GREEN else 1.5,shrinkA=38,shrinkB=38))
for k,(x,y,c,l) in N.items():
    ax.text(x,y,l,ha="center",va="center",fontsize=10,color="white" if c!=GREY else "black",
            bbox=dict(boxstyle="round,pad=0.5",fc=c,ec="none",alpha=0.95))
ax.text(0.55,0.83,"reading only: holographic states always have\na negative triple overlap (no build)",fontsize=9,color=BLUE,ha="center")
ax.set_xlim(0,1.04); ax.set_ylim(0.05,1.05)
ax.set_title("Figure 1. What ties the new pieces together. Green = built and tested here; blue = our pieces; orange = known",fontsize=12)
plt.savefig("figures/fig1-tie-map.png",dpi=130,bbox_inches="tight"); plt.close()

# ---- Fig 2: the union fractal, example A
boxes=U.EXAMPLES["A: two in box 1, one in box 3"]; maps=U.maps_from(boxes)
counts,levels=U.distinct_composites(maps,9)
fig,axs=plt.subplots(2,1,figsize=(12,6.5),gridspec_kw={"height_ratios":[1,1.3]})
ax=axs[0]; ax.set_xlim(0.9,3.1); ax.set_ylim(-0.6,2.2); ax.axis("off")
ax.plot([1,3],[1.7,1.7],"k",lw=6); ax.text(2,1.95,"the whole picture lives between the two full-box points",ha="center")
cols=[GREEN,"#7ccf9a",ORANGE]
for i,(bi,s,c) in enumerate(maps):
    lo,hi=float(s)*1+float(c),float(s)*3+float(c)
    y=1.0-0.45*i; ax.plot([lo,hi],[y,y],color=cols[i],lw=6)
    ax.annotate("",xy=((lo+hi)/2,y+0.08),xytext=(2,1.62),arrowprops=dict(arrowstyle="->",color=cols[i]))
    lab=["box 1, slab 4/5 (shrink 1/5)","box 1, slab 17/20 (shrink 3/20)","box 3, slab 9/4 (shrink 1/4)"][i]
    ax.text(hi+0.04 if i<2 else lo-0.04,y,lab,va="center",ha="left" if i<2 else "right",color=cols[i])
ax.text(1,-0.45,"centre 1 = box 1 full",ha="center",color=GREEN); ax.text(3,-0.45,"centre 3 = box 3 full",ha="center",color=ORANGE)
ax.set_title("Figure 2. The union fractal. Each map joins a slab to x inside a box; it pulls everything toward 'box full'.")
ax=axs[1]
pts=U.attractor_points(maps,levels[-1])
ax.eventplot(np.sort(pts),lineoffsets=0,colors="k",linewidths=0.3)
ax.set_yticks([]); ax.set_xlabel("x (size of the union)")
ax.set_title("The attractor: one green pair commutes (same box), so many words land on the same point")
plt.tight_layout(); plt.savefig("figures/fig2-union-fractal.png",dpi=130,bbox_inches="tight"); plt.close()

# ---- Fig 3: counts
fig,ax=plt.subplots(figsize=(9,5.5))
for name,c in zip(res,[GREEN,BLUE,GREY,ORANGE]):
    r=res[name]; n=np.arange(1,len(r["counts"])+1)
    ax.semilogy(n,r["counts"],"o",color=c,label=name.split(":")[0]+": distinct maps")
    ax.semilogy(n,r["trace_counts"],"-",color=c,alpha=0.6)
    ax.semilogy(n,[len(U.maps_from(U.EXAMPLES[name]))**i for i in n],":",color=c,alpha=0.6)
ax.set_xlabel("length of the word (number of unions)"); ax.set_ylabel("how many")
ax.set_title("Figure 3. Dots: distinct maps, counted exactly. Lines: the shared-centre formula.\nDotted: all words. Dots sit on the lines in every case.")
ax.legend(fontsize=9); plt.savefig("figures/fig3-counts.png",dpi=130,bbox_inches="tight"); plt.close()

# ---- Fig 4: dimension
fig,axs=plt.subplots(1,2,figsize=(13,5))
ax=axs[0]
for name,c in zip(res,[GREEN,BLUE,GREY,ORANGE]):
    r=res[name]; n=np.arange(1,len(r["Dn"])+1)
    ax.plot(1/n,r["Dn"],"o-",color=c,label=name.split(":")[0])
    ax.axhline(r["D_trace"],color=c,ls="--",lw=1); ax.plot(0,r["D_trace"],"*",color=c,ms=14)
    ax.plot(0,r["D_moran"],"x",color=c,ms=10,mew=2)
ax.set_xlabel("1 / word length"); ax.set_ylabel("dimension")
ax.set_title("Upper bounds D_n run straight to the star (our formula),\nnot to the x (Moran's equation, ignores commuting)")
ax.legend(fontsize=9)
ax=axs[1]; names=list(res); xs=np.arange(len(names)); w=0.27
box=[np.polyfit(np.log(1/np.array(res[n]["eps"][6:12])),np.log(np.array(res[n]["N"][6:12])),1)[0] for n in names]
ax.bar(xs-w,[res[n]["D_moran"] for n in names],w,color=GREY,label="Moran (ignores commuting)")
ax.bar(xs,[res[n]["D_trace"] for n in names],w,color=GREEN,label="our formula D*")
ax.bar(xs+w,box,w,color=BLUE,label="box counting of the picture")
ax.set_xticks(xs); ax.set_xticklabels([n.split(":")[0] for n in names]); ax.set_ylim(0.5,0.85)
ax.set_ylabel("dimension"); ax.legend(fontsize=9); ax.set_title("Measured dimension follows our formula")
plt.suptitle("Figure 4. The dimension of the union fractal",y=1.02); plt.savefig("figures/fig4-dimension.png",dpi=130,bbox_inches="tight"); plt.close()

# ---- Fig 5: the commuting discount as one slab varies
from scipy.optimize import brentq
s2=np.linspace(0.01,0.4,80)
Dt=[U.trace_root([(F(1),[F(1,5),F(str(round(v,4)))]),(F(3),[F(1,4)])]) for v in s2]
Dm=[U.moran_root([(F(1),[F(1,5),F(str(round(v,4)))]),(F(3),[F(1,4)])]) for v in s2]
Dsep=[U.moran_root([(F(1),[F(1,5)]),(F(2),[F(str(round(v,4)))]),(F(3),[F(1,4)])]) for v in s2]
fig,ax=plt.subplots(figsize=(9,5))
ax.plot(s2,Dm,color=GREY,lw=2,label="second slab in its own box (= Moran)")
ax.plot(s2,Dt,color=GREEN,lw=2.5,label="second slab in the SAME box as the first (our formula)")
ax.fill_between(s2,Dt,Dm,color=GREEN,alpha=0.15)
ax.text(0.25,0.66,"commuting discount:\ndimension lost because\norder stops mattering",color=GREEN)
ax.axhline(U.moran_root([(F(1),[F(1,5)]),(F(3),[F(1,4)])]),color=ORANGE,ls="--",label="without the second slab")
ax.set_xlabel("shrink factor of the second slab's map"); ax.set_ylabel("dimension")
ax.set_title("Figure 5. Same slabs, same shrink factors: putting two in one box lowers the dimension"); ax.legend(fontsize=9)
plt.savefig("figures/fig5-commuting-discount.png",dpi=130,bbox_inches="tight"); plt.close()

# ---- Fig 6: Carnot speed
fig,axs=plt.subplots(1,2,figsize=(13,5))
ax=axs[0]; kk=np.linspace(-2.6,0.6,200)
ax.plot(kk,kk**2/4,color="k",lw=2,label="border, curvature 0")
ax.plot(kk,kk**2/4-1,color=GREEN,lw=1,ls="--",label="curvature -1 leaf (discriminant 4)")
g=np.linspace(0,1,50); ax.plot(2*g-2,-g*(2-g),color=GREEN,lw=4,label="heat-type engines (Schmiedl-Seifert)")
ax.plot([-2],[0],"o",color=ORANGE,ms=10); ax.text(-2.5,0.15,"parity\n(lower bound)",color=ORANGE)
ax.plot([0],[-1],"o",color=BLUE,ms=10); ax.text(0.05,-0.95,"velocity addition\n(upper bound)",color=BLUE)
ax.plot([-1],[0],"s",color=RED,ms=9); ax.text(-0.95,0.12,"Curzon-Ahlborn: curvature -1/4",color=RED)
ax.set_xlabel("k"); ax.set_ylabel("d"); ax.legend(fontsize=9,loc="lower left")
ax.set_title("Every engine on the parabola has curvature exactly -1")
ax=axs[1]; b=np.linspace(0,6,200)
for gg,c in zip([0,0.33,0.66,1],["#9bd8b0","#5fc286","#2e9e5b","#1b6b3c"]):
    etaC=1-np.exp(-b); eta=etaC/(2-gg*etaC); wall=1/(2-gg)
    ax.semilogy(b,(wall-eta)*(2-gg)**2/2,color=c,lw=2,label=f"engine gamma={gg}")
ax.semilogy(b,np.exp(-b),"k--",label="Carnot's rate e^-beta (waste fraction)")
ax.semilogy(b,np.exp(-b/2),color=RED,ls=":",lw=2,label="Curzon-Ahlborn, half rate")
ax.set_xlabel("beta = ln(T_hot / T_cold)"); ax.set_ylabel("gap to the engine's limit (scaled)"); ax.legend(fontsize=9)
ax.set_title("Discriminant 4 means: every engine closes in at Carnot's rate")
plt.suptitle("Figure 6. Build 2: the meaning of the Carnot parabola, read with the curvature dial",y=1.02)
plt.savefig("figures/fig6-carnot-speed.png",dpi=130,bbox_inches="tight"); plt.close()
print("ok")
