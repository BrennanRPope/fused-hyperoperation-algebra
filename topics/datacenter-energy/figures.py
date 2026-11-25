import math, numpy as np, matplotlib
matplotlib.use("Agg"); import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
D = "/mnt/project-files/math-theory/datacenter-energy/figures/"
import os; os.makedirs(D, exist_ok=True)
plt.rcParams.update({"font.size": 11, "axes.spines.top": False, "axes.spines.right": False})
BL, OR, GR, RD, GY = "#3b6fb6", "#e08a2e", "#4a9a5b", "#c0504d", "#888"

# d1: the PUE square. IT = unit square, electrical strip p, cooling strip c, corner p*c
p, c = 0.0684, 0.4602
fig, ax = plt.subplots(figsize=(6.4, 5.2))
ax.add_patch(Rectangle((0, 0), 1, 1, fc="#dfe8f5", ec=BL, lw=1.5)); ax.text(.5, .5, "IT power\n1", ha="center", va="center", fontsize=13)
ax.add_patch(Rectangle((1, 0), p, 1, fc=OR, ec="k", alpha=.8)); ax.text(1+p/2, .5, f"power chain losses p = {p:.3f}", rotation=90, ha="center", va="center", fontsize=9)
ax.add_patch(Rectangle((0, 1), 1, c, fc="#9fd0aa", ec="k")); ax.text(.5, 1+c/2, f"cooling the IT heat  c = {c:.3f}", ha="center", va="center")
ax.add_patch(Rectangle((1, 1), p, c, fc=RD, ec="k")); ax.annotate(f"corner p·c = {p*c:.3f}\ncooling the power chain's own heat\n(5.6% of all overhead)", xy=(1+p/2, 1+c/2), xytext=(1.18, 1.25), arrowprops=dict(arrowstyle="->"), fontsize=9)
ax.set_xlim(-.05, 1.75); ax.set_ylim(-.05, 1.55); ax.set_aspect("equal"); ax.axis("off")
ax.set_title("PUE 1.56 as a mult-add square: (1+p)(1+c) = 1 + p + c + p·c")
fig.tight_layout(); fig.savefig(D+"d1_pue_square.png", dpi=150); plt.close(fig)

# d2: neper shares (order-free) of the full chip-to-meter overhead
names = ["Distribution", "UPS", "Server PSU", "Board VRM", "Cooling & other"]
eff = [0.975, 0.96, 0.95, 0.90, 1/(1+c)]
L = np.array([-math.log(e) for e in eff])
fig, ax = plt.subplots(figsize=(7.2, 3.6))
left = 0
cols = [OR, OR, "#b07acb", "#b07acb", GR]
for n, l, col in zip(names, L, cols):
    ax.barh(0, l, left=left, color=col, ec="white"); ax.text(left+l/2, 0, f"{n}\n{100*l/L.sum():.0f}%", ha="center", va="center", fontsize=8.5, rotation=0 if l > .09 else 90)
    left += l
ax.axvline(-math.log(.975*.96), color="k", ls=":"); ax.text(-math.log(.975*.96), .48, "facility side", ha="right", fontsize=8)
ax.axvline(-math.log(.975*.96*.95*.90), color="k", ls="--"); ax.text(-math.log(.975*.96*.95*.90)+.005, .48, "IT meter: PUE can't see right of here", fontsize=8)
ax.set_yticks([]); ax.set_xlabel(f"nepers = ln(input/output); total {L.sum():.3f} = ln {math.exp(L.sum()):.3f} watts drawn per chip watt")
ax.set_title("Where a chip watt's overhead goes, in log units that add in any order")
ax.set_ylim(-.6, .7); fig.tight_layout(); fig.savefig(D+"d2_neper_shares.png", dpi=150); plt.close(fig)

# d3: unit-count rule with the Titanium-minimum fit (fixed + proportional + square)
a, b, cq = 0.01762, -0.06468, 0.14581
X = np.linspace(0.05, 2.6, 400)
fig, ax = plt.subplots(figsize=(7.2, 4.4))
for k, col in zip([1, 2, 3], [BL, OR, GR]):
    loss = k*a + b*X + cq*X**2/k
    ax.plot(X, 100*loss/(X+loss), color=col, label=f"{k} unit{'s' if k>1 else ''} sharing")
xs = math.sqrt(a/cq)
for k in [1, 2]:
    xk = xs*math.sqrt(k*(k+1)); ax.axvline(xk, color=GY, ls=":"); ax.text(xk+.02, 7.6 if k == 1 else 3.0, f"switch at x*·√{k*(k+1)} = {xk:.2f}", fontsize=8.5)
ax.set_ylim(2.5, 8); ax.set_xlabel("total load, in units of one PSU's rating"); ax.set_ylabel("% of input lost")
ax.set_title("Paralleled units: add one when load passes x*·√(k(k+1))  (x* = √(a/c))")
ax.legend(frameon=False); fig.tight_layout(); fig.savefig(D+"d3_unit_rule.png", dpi=150); plt.close(fig)

# d4: fan cube and spare-unit savings
fig, axs = plt.subplots(1, 2, figsize=(9.6, 4))
d = np.linspace(-.4, .2, 200)
axs[0].plot(100*d, 100*((1+d)**3-1), color="k", lw=2, label="(1+d)³ − 1, all layers")
axs[0].plot(100*d, 100*3*d, color=BL, ls="--", label="3d (edges only)")
axs[0].plot(100*d, 100*(3*d+3*d*d), color=OR, ls=":", label="3d + 3d² (edges + faces)")
axs[0].axhline(0, color=GY, lw=.6); axs[0].axvline(-20, color=GY, ls=":"); axs[0].text(-39, 10, "−20% speed → −49% power", fontsize=9)
axs[0].set_xlabel("fan speed change d (%)"); axs[0].set_ylabel("fan power change (%)"); axs[0].legend(frameon=False, fontsize=8.5)
axs[0].set_title("The mult-add cube: layers 3d, 3d², d³")
N = np.arange(1, 13)
for m, lab, col in [(3, "cube law", "k"), (2.75, "turbulent duct, Δp∝Q^1.75", RD), (2, "laminar, Δp∝Q", BL)]:
    axs[1].plot(N, 100*(1-(N/(N+1))**(m-1)), "o-", color=col, ms=4, label=lab)
axs[1].set_xlabel("N units needed"); axs[1].set_ylabel("fan energy saved (%)"); axs[1].legend(frameon=False, fontsize=8.5)
axs[1].set_title("Run the spare too: N+1 at N/(N+1) speed")
fig.tight_layout(); fig.savefig(D+"d4_fans.png", dpi=150); plt.close(fig)
print("ok")
