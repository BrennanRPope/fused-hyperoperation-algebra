import json, numpy as np, matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from fractions import Fraction as F
import union_ifs as U
plt.rcParams.update({"font.size":11})
GREEN="#2e9e5b"; BLUE="#2f6db5"; ORANGE="#e08a1e"; GREY="#999999"; RED="#c0392b"; PURPLE="#8e44ad"
# Fig 15: where each Boolean step pulls a set's density, against the other set's density a
a=np.linspace(0.02,0.98,300)
fig,ax=plt.subplots(figsize=(10,6.2))
fixed=[("X AND A, X AND NOT A  (intersect, difference): pull toward empty",0*a,BLUE),
       ("X XOR A, X XNOR A  (symmetric difference): pull toward half",0.5+0*a,PURPLE),
       ("X OR A, X OR NOT A  (union): pull toward full",1+0*a,GREEN)]
for lab,y,c in fixed: ax.plot(a,y,color=c,lw=3.5,label=lab)
moving=[("NOT X AND A",a/(1+a)),("NOT (X AND A)",1/(1+a)),("NOT (X OR A)",(1-a)/(2-a)),("NOT X OR A",1/(2-a))]
for i,(lab,y) in enumerate(moving): ax.plot(a,y,color=GREY,lw=1.5,ls="--",label=lab+": centre moves with a" if i==0 else lab)
ax.set_xlabel("density a of the independent set you combine with"); ax.set_ylabel("density the step pulls toward (its centre)")
ax.set_ylim(-0.08,1.12)
ax.set_title("Figure 15. Every Boolean step that keeps x pulls a set's density toward a point.\n"
             "Only three points stay put whatever a is: empty (AND), half (XOR), full (OR). Those steps commute and make union fractals.")
ax.legend(fontsize=8.5,loc="center right"); plt.savefig("figures/fig15-boolean-centres.png",dpi=130,bbox_inches="tight"); plt.close()
# Fig 16: two Boolean fractals
cases=[("OR with 4/5 and 17/20, AND with 1/4 (= round 1's example A, moved)",[(F(1),[F(1,5),F(3,20)]),(F(0),[F(1,4)])],GREEN,10),
       ("OR with 4/5 and 17/20, AND with 3/13, XOR with 9/22 and 15/26",[(F(1),[F(1,5),F(3,20)]),(F(0),[F(3,13)]),(F(1,2),[F(2,11),F(-2,13)])],PURPLE,7)]
fig,axs=plt.subplots(len(cases),1,figsize=(11,4.6))
for ax,(lab,boxes,c,n) in zip(axs,cases):
    maps=U.maps_from(boxes); counts,levels=U.distinct_composites(maps,n); pts=U.attractor_points(maps,levels[-1])
    ax.eventplot(np.sort(pts),colors=c,linewidths=0.25,linelengths=0.8)
    absb=[(cc,[abs(s) for s in sl]) for cc,sl in boxes]
    ax.set_title(f"{lab}\nshared-centre formula D* = {U.trace_root(absb):.3f};  after {n} steps: {counts[-1]:,} different densities from {len(maps)**n:,} step sequences",fontsize=9.5)
    for cx,name in [(0,"empty"),(0.5,"half"),(1,"full")]:
        ax.axvline(cx,color=GREY,lw=0.8,ls=":"); ax.text(cx,1.25,name,ha="center",fontsize=8,color=GREY)
    ax.set_yticks([]); ax.set_xlim(-0.03,1.03); ax.set_ylim(0.4,1.5)
axs[-1].set_xlabel("density of the set after many random Boolean steps")
plt.suptitle("Figure 16. Boolean fractals: the densities a set can reach by repeated AND / OR / XOR with independent sets",y=1.02)
plt.tight_layout(); plt.savefig("figures/fig16-boolean-fractals.png",dpi=130,bbox_inches="tight"); plt.close(); print("ok")
