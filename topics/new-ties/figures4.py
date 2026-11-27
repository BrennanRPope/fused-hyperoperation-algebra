import json, numpy as np, matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from offline import bounds, zc, rho
plt.rcParams.update({"font.size":11})
GREEN="#2e9e5b"; BLUE="#2f6db5"; ORANGE="#e08a1e"; GREY="#999999"; RED="#c0392b"
O=json.load(open("offline.json")); DN=json.load(open("dn_levels.json")); S=json.load(open("planelaws_sweep.json"))
# Fig 13: the proved bracket against the distance of the repelling walls
fig,ax=plt.subplots(figsize=(10,6))
q=[r[0] for r in O["sweep"]]; lo=[r[1] for r in O["sweep"]]; up=[r[2] for r in O["sweep"]]; ds=O["sweep"][0][3]
ax.fill_between(q,lo,up,color=BLUE,alpha=0.35,label="proved bracket for the true dimension (width ≤ 0.002)")
ax.semilogx(q,lo,color=BLUE,lw=1); ax.semilogx(q,up,color=BLUE,lw=1)
ax.axhline(ds,color=GREEN,lw=2.5,ls="--",label=f"shared-centre formula D* = {ds:.3f} (the line case, and the ceiling)")
ax.semilogx([s["q"] for s in S],[s["box"] for s in S],"s",color=ORANGE,ms=6,label="box counting (finite resolution)")
for qq in ("2","3","12"):
    d=dict((int(n),v) for n,c,v in DN[qq])
    ax.plot([float(qq)]*2,[d[7],d[14]],"v-",color=GREY,ms=5)
ax.plot([],[],"v-",color=GREY,label="level-n estimates D₇ → D₁₄ (still falling)")
ax.set_xlabel("distance of the repelling walls (first law repels from +q, second from −q)"); ax.set_ylabel("dimension")
ax.set_title("Figure 13. Off the tangent line the true dimension is pinned down: a thin proved band,\nalways below the formula, closing on it as the repelling walls move away")
ax.legend(fontsize=9,loc="lower right"); plt.savefig("figures/fig13-proved-band.png",dpi=130,bbox_inches="tight"); plt.close()
# Fig 14: why it is smaller: where each law's maps receive their inputs, and the factor they pay there
fig,axs=plt.subplots(2,1,figsize=(11,7),gridspec_kw=dict(height_ratios=[1,1.4]))
cl=[(1.0,3.0,[0.2,0.15]),(-1.0,-3.0,[0.25])]
lo3,up3,_,info=bounds(cl,nb=3); inp=info["inputs"]
ax=axs[0]
ax.plot([-3.4,3.4],[0,0],color="k",lw=1)
for x,lab,c in [(1,"attracting wall\nof law 1",GREEN),(-1,"attracting wall\nof law 2",BLUE),(3,"repelling wall\nof law 1",GREEN),(-3,"repelling wall\nof law 2",BLUE)]:
    ax.plot([x],[0],"o" if abs(x)==1 else "x",color=c,ms=10,mew=2); ax.text(x,0.35,lab,ha="center",fontsize=9,color=c)
ax.plot(inp[0],[-0.35,-0.35],color=GREEN,lw=7,solid_capstyle="butt"); ax.text(np.mean(inp[0]),-0.75,"law 1 blocks receive points only here",ha="center",fontsize=9,color=GREEN)
ax.plot(inp[1],[-1.1,-1.1],color=BLUE,lw=7,solid_capstyle="butt"); ax.text(np.mean(inp[1]),-1.5,"law 2 blocks receive points only here",ha="center",fontsize=9,color=BLUE)
ax.set_ylim(-1.8,1.0); ax.set_xlim(-3.4,3.4); ax.axis("off")
ax.set_title("Figure 14. Why bent maps give a smaller fractal (walls ±3).\nTop: after a switch, each law gets its input away from its own attracting wall")
ax=axs[1]; y=np.linspace(-1,1,400)
for lam,ls in [(0.2,"-"),(0.04,"--")]:
    ax.plot(y,rho(zc(y,1.0,3.0),lam),color=GREEN,ls=ls,label=f"law 1, block multiplier {lam}")
    ax.plot(y,rho(zc(y,-1.0,-3.0),lam),color=BLUE,ls=ls,label=f"law 2, block multiplier {lam}")
for i,c in [(0,GREEN),(1,BLUE)]: ax.axvspan(*inp[i],color=c,alpha=0.08)
ax.axhline(1,color=GREY,lw=1); ax.set_ylim(0,1.05)
ax.set_xlabel("where the block is applied"); ax.set_ylabel("stretch it achieves ÷ its multiplier")
ax.set_title("Bottom: a block stretches by its full multiplier only at its own wall; where it actually works it does less",fontsize=10)
ax.legend(fontsize=8,loc="lower center",ncol=2)
plt.tight_layout(); plt.savefig("figures/fig14-why-smaller.png",dpi=130,bbox_inches="tight"); plt.close(); print("ok")
