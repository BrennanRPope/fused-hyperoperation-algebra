import numpy as np, json, math, matplotlib
matplotlib.use("Agg"); import matplotlib.pyplot as plt
plt.rcParams.update({"font.size":11,"axes.spines.top":False,"axes.spines.right":False,"figure.dpi":130})
C=["#2a6f97","#e07a5f","#3d405b","#81b29a","#9c6ade","#c1121f"]
def slog(x):
    n=0
    while x>1: x=math.log(x); n+=1
    return n+(x-1)
fig,ax=plt.subplots(1,3,figsize=(15,4.4))
# (a) height reached vs time for fixed-step flows, measured in "number of exponentials" = slog(Y)
t=np.linspace(0,6,601); Y0=20.0
h2=[slog(Y0*math.exp(tt)) for tt in t]
h3=[slog(math.exp(math.log(Y0)*math.exp(tt))) if math.log(Y0)*math.exp(tt)<700 else slog(math.log(Y0)*math.exp(tt))+1 for tt in t]
h4=[slog(Y0)+tt for tt in t]
h5=[]  # pentation flow: slog grows like tetration of t; draw slog(slog-height) linear -> height = exp-tower in t, cap the plot
for tt in t:
    v=slog(Y0)
    for _ in range(int(tt)): v=min(v*3,60)
    h5.append(min(60, v*(3**(tt-int(tt))) if v<60 else 60))
ax[0].plot(t,h2,color=C[0],lw=2,label="rung 2 (exponential)")
ax[0].plot(t,h3,color=C[1],lw=2,label="rung 3 (double exponential)")
ax[0].plot(t,h4,color=C[4],lw=2,label="rung 4 (tetration): one more exp per step")
T=1/Y0; ax[0].axvline(T,color="k",ls="--",lw=2); ax[0].text(T+0.08,8.4,"Y′ = Y²: infinite at t = 0.05",fontsize=9)
ax[0].axhline(4,color="#bbb",lw=1,ls=":"); ax[0].text(3.2,4.15,"height 4 ≈ 10^1,650,000",fontsize=8.5,color="#555")
ax[0].set_ylim(0,9); ax[0].set_xlabel("time"); ax[0].set_ylabel("height = number of stacked exponentials (log*)")
ax[0].set_title("Higher rungs climb absurdly fast, but every height is finite",fontsize=10.5); ax[0].legend(frameon=False,fontsize=8.5,loc="lower right")
# (b) step-time ratio signatures
d=np.load("higher_signatures.npz"); n=d["n"]
for key,lab,c in [("rung2","rung 2: steps stay equal",C[0]),("rung3","rung 3: steps shrink like 1/n",C[1]),("rung4","rung 4: like 1/(n log n)",C[4]),("pole","self-feeding pole: steps shrink by a fixed fraction",C[2])]:
    v=d[key]; ax[1].semilogy(n,v/v[0],color=c,lw=2,label=lab)
ax[1].set_xlabel("shell n"); ax[1].set_ylabel("time to reach the next shell (relative)")
ax[1].set_title("Signature of each growth law in a cascade",fontsize=10.5); ax[1].legend(frameon=False,fontsize=8.5,loc="lower left"); ax[1].set_ylim(1e-12,3)
# (c) data: arrival fronts, time left vs shell
res=json.load(open("stretching_signature.json"))
for f in res:
    T=np.array(f["T"]); nn=np.array(f["n"]); ts=f["tstar"]
    if f["kind"]=="sabra":
        ax[2].semilogy(nn,ts-T,"o-",color=C[5],ms=4,lw=2,label=f"Sabra control: fixed fraction r = {f['ratio']:.2f}",zorder=4)
    else:
        col={0.5:C[4],4.0:C[1]}.get(f["eps"],C[2]); ax[2].semilogy(nn,ts-T,"-",color=col,lw=1,alpha=.7)
for e,c in [(0.5,C[4]),(4.0,C[1]),(float('inf'),C[2])]:
    ax[2].plot([],[],color=c,label=f"helical toy, ε = {'cross-only' if e==float('inf') else e} (4 seeds)")
ax[2].set_xlabel("shell n"); ax[2].set_ylabel("time left before blow-up, t* − tₙ")
ax[2].set_title("Data: straight lines on this scale = self-feeding pole",fontsize=10.5); ax[2].legend(frameon=False,fontsize=8.5,loc="lower left")
plt.tight_layout(); plt.savefig("figures/n8_higher_rungs.png"); plt.close(); print("ok")
