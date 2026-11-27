import json, numpy as np, matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
plt.rcParams.update({"font.size":11})
GREEN="#2e9e5b"; BLUE="#2f6db5"; ORANGE="#e08a1e"; GREY="#999999"; RED="#c0392b"
R=json.load(open("planelaws.json")); S=json.load(open("planelaws_sweep.json"))
# Fig 11: tangent lines
fig,ax=plt.subplots(figsize=(10,7))
k=np.linspace(-4.2,4.2,400); ax.plot(k,k**2/4,"k",lw=2.5,label="border: two walls merge (d = k²/4)")
for wv,c in [(1,GREEN),(-1,BLUE),(2,ORANGE),(0.5,RED),(-2,"#8e44ad")]:
    ax.plot(k,(-1-k*wv)/wv**2,color=c,lw=1.6,label=f"laws sharing the wall {wv}")
    kt=-2/wv; ax.plot(kt,kt**2/4,"o",color=c,ms=7)
ax.plot(k,0*k,color=GREY,lw=1.6,label="laws sharing the wall at infinity (d = 0): straight-line maps")
for (x,y,lab) in [(0,0,"plain addition"),(-1,0,"union rule"),(0,-1,"velocity addition"),(-2,0,"parity"),(-2,1,"double wall 1")]:
    ax.plot(x,y,"ks",ms=6); ax.text(x+0.08,y+0.12,lab,fontsize=9)
ax.set_xlim(-4.2,4.2); ax.set_ylim(-3,4.5); ax.set_xlabel("k"); ax.set_ylabel("d")
ax.set_title("Figure 11. Laws that share a wall lie on one tangent line of the border.\nAlong a tangent line the union fractal straightens out and our formula is exact.")
ax.legend(fontsize=8,loc="upper center"); plt.savefig("figures/fig11-tangent-lines.png",dpi=130,bbox_inches="tight"); plt.close()
# Fig 12
fig,axs=plt.subplots(1,2,figsize=(13,5))
ax=axs[0]
for lab,y,c in [("shared wall 3",1.0,GREEN),("no shared wall, near",0.0,ORANGE)]:
    P=np.array(R[lab]["points"]); ax.eventplot(np.sort(P),lineoffsets=y,linelengths=0.6,colors=c,linewidths=0.3)
    ax.text(-1.05,y+0.42,f"{lab}: measured {R[lab]['box']:.3f}, formula {R[lab]['D_mult']:.3f}",color=c,fontsize=10)
ax.set_yticks([]); ax.set_ylim(-0.6,1.7); ax.set_xlabel("x")
ax.set_title("Same multipliers, same walls 1 and -1.\nOnly where the repelling walls sit differs")
ax=axs[1]; q=[s["q"] for s in S]
ax.semilogx(q,[s["D_mult"] for s in S],color=GREEN,lw=2.5,label="our formula = proved upper bound")
ax.semilogx(q,[s["Dn"] for s in S],"o-",color=BLUE,label="dimension from the distinct maps (D_n)")
ax.semilogx(q,[s["box"] for s in S],"s",color=ORANGE,label="box counting")
ax.set_xlabel("distance of the repelling walls (±q)"); ax.set_ylabel("dimension"); ax.legend(fontsize=9)
ax.set_title("Without a shared wall the formula is an upper bound;\nthe gap closes as the repelling walls move away")
plt.suptitle("Figure 12. Plane laws: exact on a tangent line, an upper bound off it",y=1.02)
plt.savefig("figures/fig12-plane-laws.png",dpi=130,bbox_inches="tight"); plt.close(); print("ok")
