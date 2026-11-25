import numpy as np, json, collections, matplotlib
matplotlib.use("Agg"); import matplotlib.pyplot as plt
plt.rcParams.update({"font.size":11,"axes.spines.top":False,"axes.spines.right":False,"figure.dpi":130})
C=["#2a6f97","#e07a5f","#3d405b","#81b29a","#9c6ade","#c1121f"]
theta=1e-2
def arrivals(f):
    d=np.load(f); En=d['En']; k=d['k']; t=d['t']; out=[]
    for n in range(len(k)):
        h=np.where(En[n]>=theta*k[n]**(-2/3))[0]; out.append(t[h[0]] if len(h) else np.nan)
    return np.array(out), t[-1]
fig,ax=plt.subplots(1,3,figsize=(14,4.3),gridspec_kw={"width_ratios":[1.1,1.1,0.9]})
for f,lab,c in [("shell_eps0.00.npz","ε = 0 (same-handed only)",C[3]),("shell_eps0.25.npz","ε = 0.25",C[0]),
                ("shell_eps0.50.npz","ε = 0.5",C[4]),("shell_eps4.00.npz","ε = 4",C[1]),("shell_epsinf.npz","cross-handed only",C[2])]:
    tn,te=arrivals(f); n=np.arange(len(tn))
    # running max of the first arrival to show the front
    front=np.fmax.accumulate(np.where(np.isnan(tn),np.nan,tn)); last=np.where(~np.isnan(tn))[0].max(); front[last+1:]=np.nan
    ax[0].plot(n,tn,"o",ms=3.5,color=c,alpha=.55); ax[0].plot(n,front,"-",color=c,lw=2,label=lab)
ax[0].set_yscale("log"); ax[0].set_ylim(1e-2,60); ax[0].set_xlabel("shell n   (size of eddy shrinks as 1.618⁻ⁿ)"); ax[0].set_ylabel("time energy first reaches shell n")
ax[0].axhline(12,color="#999",ls=":",lw=1); ax[0].text(0.3,13.5,"ε = 0 run stops here, front stalled at shell 9",fontsize=8.5,color="#555")
ax[0].legend(frameon=False,fontsize=8.5,loc="lower right"); ax[0].set_title("A wall flattening out = all scales reached in finite time",fontsize=10.5)
r=json.load(open("scan_eps.json"))+json.load(open("scan_eps005_long.json"))
d=collections.defaultdict(list)
for x in r:
    if x['tstar'] is not None and str(x['eps'])!='inf': d[float(x['eps'])].append(float(x['tstar']))
for e,v in d.items():
    col=C[5] if e==1.0 else C[0]
    ax[1].scatter([e]*len(v),v,color=col,s=28,alpha=.8,zorder=3)
es=np.array(sorted(e for e in d if e<=0.5)); med=np.array([np.median(d[e]) for e in es])
A=np.polyfit(np.log(es),np.log(med),1)
xx=np.logspace(np.log10(0.04),np.log10(0.6),50); ax[1].plot(xx,np.exp(np.polyval(A,np.log(xx))),color=C[0],lw=1.5,ls="--",label=f"fit t* ∝ ε^{A[0]:.2f}  (ε ≤ 0.5)")
ax[1].scatter([],[],color=C[5],label="ε = 1: hidden extra invariant\n(only 2 of 6 seeds by t = 40)")
ax[1].set_xscale("log"); ax[1].set_yscale("log"); ax[1].set_xlabel("cross-handed dial ε"); ax[1].set_ylabel("blow-up time t*  (6 seeds each)")
ax[1].set_title("Smoothness breaks at every ε tried; later as ε → 0",fontsize=10.5); ax[1].legend(frameon=False,fontsize=8.5,loc="upper right",bbox_to_anchor=(1.0,0.72))
ax[1].text(0.042,0.13,"ε = 0: never\n(conserved Σ k|u⁺|² caps it)",fontsize=8.5,color=C[3])
rN=json.load(open("scan_N.json")); dN=collections.defaultdict(list)
for x in rN: dN[(x['eps'],x['N'])].append(float(x['tstar']))
for e,c in [(0.5,C[4]),(4.0,C[1])]:
    Ns=sorted({k[1] for k in dN if k[0]==e}); med=[np.median(dN[(e,N)]) for N in Ns]
    for N in Ns: ax[2].scatter([N]*len(dN[(e,N)]),dN[(e,N)],color=c,s=20,alpha=.6)
    ax[2].plot(Ns,med,"-",color=c,lw=2,label=f"ε = {e}: median")
ax[2].set_xlabel("number of shells N (smallest eddy 1.618⁻ᴺ)"); ax[2].set_ylabel("t*"); ax[2].set_ylim(0,1.6)
ax[2].set_title("Adding 15 smaller shells costs no extra time",fontsize=10.5); ax[2].legend(frameon=False,fontsize=9)
plt.tight_layout(); plt.savefig("figures/n4_shell_dial.png"); plt.close()
print("fit exponent",A[0], {e:round(float(np.median(v)),3) for e,v in sorted(d.items())})
# spectra heat maps
fig,ax=plt.subplots(1,2,figsize=(12,4),sharey=True)
for a,f,ttl in [(ax[0],"shell_eps0.00.npz","ε = 0: same-handed only (front creeps and fades)"),(ax[1],"shell_eps0.50.npz","ε = 0.5: cross-handed on (energy runs to the smallest shell)")]:
    d=np.load(f); t=d['t']; En=d['En']; k=d['k']
    comp=np.log10(En*k[:,None]**(2/3)+1e-30)
    tt=np.linspace(0,t[-1],600); M=np.array([np.interp(tt,t,row) for row in comp])
    m=a.pcolormesh(tt,np.arange(len(k)),M,vmin=-8,vmax=1,cmap="magma",shading="auto"); a.set_title(ttl,fontsize=10.5); a.set_xlabel("time")
ax[0].set_ylabel("shell n"); fig.colorbar(m,ax=ax,label="log₁₀ (shell energy ÷ Kolmogorov level)",pad=0.01)
plt.savefig("figures/n6_spectra.png",bbox_inches="tight"); plt.close()
