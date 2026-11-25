import numpy as np, matplotlib
matplotlib.use("Agg"); import matplotlib.pyplot as plt
plt.rcParams.update({"font.size":11,"axes.spines.top":False,"axes.spines.right":False,"figure.dpi":130})
C=["#2a6f97","#e07a5f","#3d405b","#81b29a","#9c6ade"]
fig,ax=plt.subplots(1,2,figsize=(12,4.4))
# left: model ODEs in rung-3 coordinates w = log log Y
t=np.linspace(0,4,600); Y0=np.e**2
w0=np.log(np.log(Y0))
ax[0].plot(t,w0+t,color=C[1],lw=2.5,label="Y' = Y log Y  (rung 3): straight line, slope 1")
ax[0].plot(t,np.log(np.log(Y0)+t),color=C[0],lw=2,label="Y' = Y  (rung 2): bends down")
for p,c,ls in [(1.2,C[4],"-."),(2,"k","--")]:
    T=Y0**(1-p)/(p-1); tt=np.linspace(0,T*0.99999999,500); Y=(Y0**(1-p)-(p-1)*tt)**(1/(1-p))
    ax[0].plot(tt,np.log(np.log(Y)),color=c,ls=ls,lw=2,label=f"Y' = Y^{p}: turns vertical at t = {T:.3g}")
ax[0].set_ylim(0.3,5.2); ax[0].set_xlim(0,4); ax[0].set_xlabel("time"); ax[0].set_ylabel("w = log log Y")
ax[0].set_title("Rung-3 coordinates: the safe edge is a straight line",fontsize=10.5); ax[0].legend(frameon=False,fontsize=8.5,loc="upper left",bbox_to_anchor=(0.06,1.0))
# right: helical shell model enstrophy in the same coordinates
for f,lab,c in [("shell_eps0.00.npz","ε = 0",C[3]),("shell_eps0.25.npz","ε = 0.25",C[0]),("shell_eps0.50.npz","ε = 0.5",C[4]),("shell_epsinf.npz","cross-handed only",C[2])]:
    d=np.load(f); Z=d['Z']; tt=d['t']; w=np.log(np.log(np.e*Z/Z[0]))
    ax[1].plot(tt,w,color=c,lw=1.6,label=lab)
ax[1].set_xscale("symlog",linthresh=1); ax[1].set_xlabel("time (linear to 1, then log)"); ax[1].set_ylabel("log log (e·Z/Z₀),  Z = enstrophy")
ax[1].set_title("Shell model in rung-3 coordinates: cross-handed runs turn vertical",fontsize=10.5); ax[1].legend(frameon=False,fontsize=9)
plt.tight_layout(); plt.savefig("figures/n7_rung3_lens.png"); plt.close(); print("ok")
