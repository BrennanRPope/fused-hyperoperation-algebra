import json, numpy as np, matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
plt.rcParams.update({"font.size":11})
GREEN="#2e9e5b"; BLUE="#2f6db5"; ORANGE="#e08a1e"; GREY="#999999"; RED="#c0392b"
R=json.load(open("round2.json")); R1=json.load(open("union_ifs.json"))

# Fig 7 sandwich
fig,ax=plt.subplots(figsize=(10,5.5))
for (name,r),c in zip(R["R1"].items(),[GREEN,BLUE,ORANGE]):
    L=[x[0] for x in r["sdim_S"]]; d=[x[1] for x in r["sdim_S"]]
    ax.plot(L,d,"o-",color=c,label=name.split(":")[0]+": lower bounds (free sub-fractal)")
    Dn=R1[name]["Dn"]; ax.plot(range(1,len(Dn)+1),Dn,"s--",color=c,alpha=0.7,label=name.split(":")[0]+": upper bounds D_n")
    ax.axhline(r["D_trace"],color=c,lw=0.8)
ax.set_ylim(0.45,0.82); ax.set_xlabel("length L (generators up to length L, or words of length n)"); ax.set_ylabel("dimension")
ax.set_title("Figure 7. The squeeze: lower bounds climb and upper bounds fall to the same line, our formula D*")
ax.legend(fontsize=8,ncol=2); plt.savefig("figures/fig7-squeeze.png",dpi=130,bbox_inches="tight"); plt.close()

# Fig 8 plane
r=R["R3"]; P=np.array(r["points"])
fig,axs=plt.subplots(1,2,figsize=(13,5.5),gridspec_kw={"width_ratios":[1.3,1]})
ax=axs[0]; ax.scatter(P[:,0],P[:,1],s=0.3,c="k"); ax.set_aspect("equal")
for (x,y),lab,c in [((0,0),"centre 0: two maps\n(one turns 45 degrees)",GREEN),((1,0),"centre 1: turns 90 degrees",BLUE),((0,1),"centre i",ORANGE)]:
    ax.plot(x,y,"o",color=c,ms=10); ax.text(x+0.03,y+0.03,lab,color=c,fontsize=9)
ax.set_title("The union fractal in the plane: four turning, shrinking maps")
ax=axs[1]; ax.bar([0,1,2],[r["D_moran"],r["D_trace"],r["box"]],color=[GREY,GREEN,BLUE])
ax.set_xticks([0,1,2]); ax.set_xticklabels(["Moran","our formula","box counting"]); ax.set_ylim(1,1.3)
for i,v in enumerate([r["D_moran"],r["D_trace"],r["box"]]): ax.text(i,v+0.005,f"{v:.3f}",ha="center")
ax.set_ylabel("dimension"); ax.set_title("Same formula in two dimensions")
plt.suptitle("Figure 8. Generalization to the plane (complex similarities)",y=1.0)
plt.savefig("figures/fig8-plane.png",dpi=130,bbox_inches="tight"); plt.close()

# Fig 9 where it holds
rows=["straight-line maps, shrinking\n(the union fractal)","straight-line maps with flips\n(past the wall)","plane: turning + shrinking","plane laws (k, d): Möbius maps","maps that collide by accident"]
cols=["count of distinct\nmaps","upper bound\non dimension","exact\ndimension"]
val=[["yes","proved","proved*"],["yes","proved","proved*"],["yes","proved","numerical"],["yes","needs\nderivatives","open"],["no","proved","lower"]]
colr={"yes":GREEN,"proved":GREEN,"proved*":GREEN,"numerical":BLUE,"needs\nderivatives":ORANGE,"open":GREY,"no":RED,"lower":RED}
fig,ax=plt.subplots(figsize=(10,5)); ax.axis("off")
for j,c in enumerate(cols): ax.text(0.42+0.2*j,0.95,c,ha="center",va="center",fontweight="bold")
for i,rw in enumerate(rows):
    y=0.8-0.17*i; ax.text(0.0,y,rw,va="center")
    for j in range(3):
        v=val[i][j]; ax.text(0.42+0.2*j,y,v,ha="center",va="center",color="white",fontsize=10,
            bbox=dict(boxstyle="round,pad=0.4",fc=colr[v],ec="none"))
ax.text(0,0.0,"* for rational or algebraic numbers, using Hochman's 2014 theorem, and assuming no accidental collisions",fontsize=9)
ax.set_xlim(0,1); ax.set_ylim(-0.05,1.02)
ax.set_title("Figure 9. Where the formula holds")
plt.savefig("figures/fig9-where-it-holds.png",dpi=130,bbox_inches="tight"); plt.close()

# Fig 10 learning
r=R["R5"]; n=np.arange(1,len(r["counts"])+1)
fig,axs=plt.subplots(1,2,figsize=(13,5))
ax=axs[0]; ax.semilogy(n,3.0**n,":",color=GREY,label="trial histories (3^n)")
ax.semilogy(n,r["counts"],"o-",color=GREEN,label="different learned strengths (exact)")
ax.set_xlabel("number of trials"); ax.set_ylabel("how many"); ax.legend()
ax.set_title("Two reward cues toward the same target commute,\nso the learner forgets their order")
ax=axs[1]; vals=[r["bits_if_all_differ"],r["bits_per_trial_max"],r["hrw_bits_increment"]]
ax.bar([0,1,2],vals,color=[GREY,GREEN,BLUE])
for i,v in enumerate(vals): ax.text(i,v+0.02,f"{v:.2f}",ha="center")
ax.set_xticks([0,1,2]); ax.set_xticklabels(["if every order\nmattered","most it can keep\n(our count)","kept on average\n(random trials)"])
ax.set_ylabel("bits of history per trial"); ax.set_title("How much history the learner's state carries")
plt.suptitle("Figure 10. Scenario: Rescorla–Wagner learning (also any moving average with several rates)",y=1.02)
plt.savefig("figures/fig10-learning.png",dpi=130,bbox_inches="tight"); plt.close()
print("ok")
