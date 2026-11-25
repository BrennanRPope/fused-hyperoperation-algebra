import numpy as np, matplotlib, json, os
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D  # noqa
plt.rcParams.update({"font.size":11,"axes.spines.top":False,"axes.spines.right":False,"figure.dpi":130})
C=["#2a6f97","#e07a5f","#3d405b","#81b29a","#f2cc8f","#9c6ade"]
os.makedirs("figures",exist_ok=True)

# ---- f1: the ladder of fixed-step flows and Osgood's line
fig,ax=plt.subplots(1,2,figsize=(11,4.3))
t=np.linspace(0,6,600)
Y0=np.e**np.e*1.0001
curves=[("rung 1: Y' = 1  (line)", np.log(np.log(np.log(Y0+t)))),
        ("rung 2: Y' = Y  (exponential)", np.log(np.log(np.log(Y0)+t))),
        ("rung 3: Y' = Y log Y  (double exp)", np.log(np.log(Y0)*np.exp(t)) .__array__()*0 + np.log(np.log(Y0))+t),
        ]
# rung 4: Y' = Y logY loglogY  -> logloglog Y grows linearly
r4=np.log(np.log(np.log(Y0)))+t
lab_y="log log log Y"
ax[0].plot(t,curves[0][1],color=C[3],lw=2,label=curves[0][0])
ax[0].plot(t,curves[1][1],color=C[0],lw=2,label=curves[1][0])
ax[0].plot(t,np.log(curves[2][1]),color=C[1],lw=2,label=curves[2][0])
ax[0].plot(t,r4,color=C[5],lw=2,label="rung 4: Y' = Y logY loglogY  (tower)")
# self-feeding Y' = Y^2 from same start: Y = Y0/(1-Y0 t)
T=1/Y0; tt=np.linspace(0,T*0.999999,400); Ys=Y0/(1-Y0*tt)
ax[0].plot(tt,np.log(np.log(np.log(Ys))),color="k",lw=2.5,ls="--",label="self-feeding Y' = Y²: ends at t = 1/Y₀")
ax[0].set_xlabel("time"); ax[0].set_ylabel(lab_y); ax[0].set_ylim(-0.1,6); ax[0].legend(fontsize=8.5,frameon=False,loc="upper left")
ax[0].set_title("Every rung's flow lives forever; only self-feeding ends",fontsize=11)
X=np.logspace(2,300,300); lX=np.log(X)
I3=np.log(lX)-np.log(np.log(20)); I4=np.log(np.log(lX))-np.log(np.log(np.log(20)))
I2=1/20-1/X; I15=2/np.sqrt(20)-2/np.sqrt(X)
ax[1].plot(lX,I3,color=C[1],lw=2,label="∫ dY/(Y log Y)  → ∞ (slowly)")
ax[1].plot(lX,I4,color=C[5],lw=2,label="∫ dY/(Y logY loglogY)  → ∞ (slower)")
ax[1].plot(lX,I15,color="k",lw=2,ls=":",label="∫ dY/Y^1.5  → 0.447")
ax[1].plot(lX,I2*10,color="k",lw=2,ls="--",label="∫ dY/Y²  → 0.05 (×10)")
ax[1].set_xscale("log"); ax[1].set_xlabel("log X  (how far up you go)"); ax[1].set_ylabel("time to climb from 20 to X")
ax[1].set_title("Osgood's line: infinite climb time = no blow-up",fontsize=11); ax[1].legend(fontsize=8.5,frameon=False)
plt.tight_layout(); plt.savefig("figures/n1_ladder.png"); plt.close()

# ---- f2: the two-wave cube: Burgers (corner 0) vs KdV (corner A12)
import importlib.util, sys
spec=importlib.util.spec_from_file_location("bk","burgers_kdv_cube.py")
import io,contextlib
src=open("burgers_kdv_cube.py").read()
ns={}
with contextlib.redirect_stdout(io.StringIO()): exec(src,ns)
x=ns["x"]
fig=plt.figure(figsize=(11.5,7.2))
def cube(ax,corner,title,col):
    pts={"1":(0,0),"e₁":(1,0),"e₂":(0,1),corner:(1,1)}
    for a,b in [("1","e₁"),("1","e₂"),("e₁",corner),("e₂",corner)]:
        ax.plot(*zip(pts[a],pts[b]),color="#999",lw=1.5,ls="--" if corner.startswith("0") and corner in (a,b) else "-")
    for k,(px,py) in pts.items():
        empty=k.startswith("0")
        ax.scatter(px,py,s=900,color="white" if empty else col,edgecolor=col,lw=2,zorder=3)
        ax.text(px,py,k,ha="center",va="center",fontsize=10 if len(k)<6 else 8,color=col if empty else "white",zorder=4,weight="bold")
    ax.set_xlim(-.45,1.45); ax.set_ylim(-.45,1.45); ax.axis("off"); ax.set_title(title,fontsize=10.5)
a1=fig.add_axes([0.02,0.55,0.18,0.38]); cube(a1,"0·e₁e₂","Burgers (Cole–Hopf)\ntop corner must be 0",C[1])
a2=fig.add_axes([0.02,0.07,0.18,0.38]); cube(a2,"A·e₁e₂","KdV (Hirota)\ncorner A = ((k₁−k₂)/(k₁+k₂))²",C[0])
ts=np.linspace(0,20,201); U=np.array([ns["burgers"](tt)[0] for tt in ts])
b=fig.add_axes([0.27,0.55,0.4,0.38]); m=b.pcolormesh(x,ts,U,cmap="magma",shading="auto"); b.set_xlim(-10,25)
b.set_ylabel("time"); b.set_title("Burgers u(x,t): two shocks meet at t = 8 and fuse into one",fontsize=10.5); fig.colorbar(m,ax=b,pad=0.01)
b2=fig.add_axes([0.73,0.55,0.25,0.38])
for tt,cc in [(0,C[2]),(6,C[0]),(16,C[1])]:
    b2.plot(x,ns["burgers"](tt)[0],color=cc,lw=2,label=f"t={tt}")
b2.set_xlim(-10,25); b2.legend(frameon=False,fontsize=9); b2.set_title("plateaus 2 | 1 | 0  →  2 | 0",fontsize=10.5)
ts=np.linspace(0,50,251); K=np.array([ns["kdv"](tt) for tt in ts])
k=fig.add_axes([0.27,0.07,0.4,0.38]); m=k.pcolormesh(x,ts,K,cmap="viridis",shading="auto"); k.set_xlim(-15,50)
k.set_ylabel("time"); k.set_xlabel("x"); k.set_title("KdV u(x,t): solitons cross, both survive, both shifted",fontsize=10.5); fig.colorbar(m,ax=k,pad=0.01)
k2=fig.add_axes([0.73,0.07,0.25,0.38])
for tt,cc in [(0,C[2]),(8,C[0]),(50,C[1])]:
    k2.plot(x,ns["kdv"](tt),color=cc,lw=2,label=f"t={tt}")
k2.set_xlim(-15,50); k2.legend(frameon=False,fontsize=9); k2.set_xlabel("x"); k2.set_title("heights 0.5 and 0.18 kept",fontsize=10.5)
plt.savefig("figures/n2_cube_waves.png"); plt.close()

# ---- f3: triad-sign cube
fig=plt.figure(figsize=(6.4,5.6)); ax=fig.add_subplot(111,projection="3d")
phi=(1+5**.5)/2
import itertools
for s in itertools.product((1,-1),repeat=3):
    p=[(1-v)/2 for v in s]; layer=sum(1 for v in s if v<0)
    G=s[0]+s[1]*phi+s[2]*phi**2
    homo=layer in (0,3)
    col=C[3] if homo else (C[1])
    if abs(G)<1e-9: col="#bbbbbb"
    ax.scatter(*p,s=420,color=col,edgecolor="k",depthshade=False)
    lab="".join("+" if v>0 else "−" for v in s)
    ax.text(p[0],p[1],p[2]+0.13,lab,ha="center",fontsize=11,weight="bold")
for s in itertools.product((0,1),repeat=3):
    for d in range(3):
        if s[d]==0:
            q=list(s); q[d]=1; ax.plot(*zip(s,q),color="#888",lw=1)
ax.plot([0,1],[0,1],[0,1],color=C[3],lw=3,alpha=.6)
ax.set_axis_off(); ax.view_init(18,-58)
ax.set_title("The triad-sign cube: helicity signs of the 3 waves in a triad\n"
             "green = same-handed poles (layers 0 and 3), orange = cross-handed middle (layers 1, 2)\n"
             "grey = coupling vanishes on golden-ratio shells (flat triangle)",fontsize=9.5)
plt.tight_layout(); plt.savefig("figures/n3_sign_cube.png"); plt.close()

# ---- f5: restricted Euler 2D vs 3D
d=np.load("restricted_euler.npz")
fig,ax=plt.subplots(1,1,figsize=(6,4))
ax.semilogy(d["t"],d["n"],color=C[1],lw=2,label="3D: |A| (self-feeding, Y' ~ Y²)")
ax.axvline(float(d["T"]),color="k",ls=":",lw=1); ax.text(float(d["T"]),2,"  blow-up",fontsize=9)
ax.semilogy([0,float(d["T"])],[d["n"][0]]*2,color=C[3],lw=2,label="2D: |A| stays put (right side is 0)")
ax.set_xlabel("time"); ax.set_ylabel("size of velocity gradient"); ax.legend(frameon=False,fontsize=9)
ax.set_title("Restricted Euler: the square cancels, the cube runs away",fontsize=10.5)
plt.tight_layout(); plt.savefig("figures/n5_restricted_euler.png"); plt.close()
print("figures 1,2,3,5 done")
