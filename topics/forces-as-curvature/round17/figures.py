"""DRAFT picture for round 17. Run after checks.py. Blue = the law, orange = loop read-back, grey = reference."""
import json, math
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

OR, BL, GR = "#e07b00", "#1f6fc2", "#777777"
R = json.load(open("results.json"))
plt.rcParams.update({"font.size": 10.5, "figure.dpi": 130})
fig = plt.figure(figsize=(21, 5.6))
ax = [fig.add_subplot(1, 4, 1), fig.add_subplot(1, 4, 2), fig.add_subplot(1, 4, 3, projection="3d"), fig.add_subplot(1, 4, 4)]

a = ax[0]
for L, mk in zip(range(3), ("s", "^", "o")):
    lv = R["spectrum"][str(L)]["levels"]
    a.plot([L] * 3, lv, mk, color=OR, ms=9, mfc="none" if L < 2 else OR)
    for v, mult in zip(lv, (6, 16, 30)):
        a.annotate(f"{v:.3f}" + (f"  ({mult} modes)" if L == 2 else ""), (L, v), textcoords="offset points", xytext=(9, -3), fontsize=8.5, color=OR)
for n2 in (4, 9, 16):
    a.axhline(n2, color=GR, ls="--", lw=0.9)
    a.text(3.55, n2 + 0.25, f"(n + 1)² = {n2}", color=GR, fontsize=8.5, ha="right")
a.set_xticks(range(3)); a.set_xticklabels([f"{L}\n({600 * 8**L:,})" for L in range(3)])
a.set_xlim(-0.4, 3.6); a.set_ylim(2.5, 19.5)
a.set_xlabel("refinement level of the law's 3-sphere (pieces)"); a.set_ylabel("frequency² × radius²")
a.set_title("A. Light on the law's closed space:\nlevels of 6, 16, 30 modes at (n + 1)²")

b = ax[1]
rows = [r for r in R["hands"] if r["level"] == 2]
for i, r in enumerate(rows):
    w = math.sqrt(r["cluster"])
    b.barh(i + 0.18, r["plus"], height=0.34, color=BL, label="right-handed (curl > 0)" if i == 0 else None)
    b.barh(i - 0.18, -r["minus"], height=0.34, color=OR, label="left-handed (curl < 0)" if i == 0 else None)
    b.text(r["plus"] + 0.4, i + 0.18, f"{r['plus']} at {r['mean_plus']:+.3f}", va="center", fontsize=8.5)
    b.text(-r["minus"] - 0.4, i - 0.18, f"{r['minus']} at {r['mean_minus']:+.3f}", va="center", ha="right", fontsize=8.5)
b.axvline(0, color="k", lw=0.8)
b.set_yticks(range(3)); b.set_yticklabels([f"level n + 1 = {int(math.sqrt(r['cluster']))}" for r in rows])
b.set_xlim(-26, 26); b.set_xlabel("number of modes (left ← → right)")
b.legend(fontsize=8.5, loc="lower right")
b.set_title("B. Each level splits evenly into two hands,\ncurl = ±(n + 1) (level 2)")

c = ax[2]
def fibre(theta, phi, N=300):          # the Hopf circle over the point (theta, phi) of the velocity sphere, stereographic to 3D
    s = np.linspace(0, 2 * np.pi, N)
    z1 = np.cos(theta / 2) * np.exp(1j * (s + phi)); z2 = np.sin(theta / 2) * np.exp(1j * s)
    x = np.stack([z1.real, z1.imag, z2.real, z2.imag]); return x[:3] / (1 - x[3])
for th, col in ((0.55, BL), (1.05, OR)):
    for ph in np.linspace(0, 2 * np.pi, 11)[:-1]:
        X = fibre(th, ph); c.plot(X[0], X[1], X[2], color=col, lw=1.0)
c.view_init(elev=28, azim=35)
c.set_box_aspect((1, 1, 1)); c.set_xlim(-1.6, 1.6); c.set_ylim(-1.6, 1.6); c.set_zlim(-1.6, 1.6); c.set_axis_off()
h2 = R["hopf_field"][2]
c.text2D(0.02, 0.02, f"on the law's pieces (level 2): overlap with the lowest level {h2['hopf']['overlap']:.4f},\n"
         f"frequency² {h2['hopf']['freq2']:.3f}, curl {h2['hopf']['curl']:+.3f}; mirror image curl {h2['mirror']['curl']:+.3f}",
         transform=c.transAxes, fontsize=8.5, color=BL)
c.text2D(0.02, 0.92, "the law's twist turns along these circles: each lies over one\nvelocity (blue: one ring of velocities, orange: another); every two link", transform=c.transAxes, fontsize=8.5, color=GR)
c.set_title("C. The law's twist on the 3-sphere is the Hopf field:\na lowest photon of one hand")

d = ax[3]
caps = R["hopf_law"]["caps"]
ar = np.array([r["area"] for r in caps]) / 2
d.plot([0, 6.5], [0, 6.5], color=GR, ls="--", lw=1, label="half the enclosed area")
d.plot(ar, [abs(r["law"]) for r in caps], "o", color=BL, ms=11, mfc="none", label="the law's charge-½ phase (round 10)")
d.plot(ar, [abs(r["hopf"]) for r in caps], "x", color=OR, ms=9, mew=2, label="the Hopf connection's phase")
d.set_xlabel("half the area a loop of velocities encloses"); d.set_ylabel("phase picked up round the loop")
d.legend(fontsize=8.5, loc="upper left")
ins = d.inset_axes([0.58, 0.1, 0.38, 0.32])
fe = [r["flux_err"] for r in R["hopf_field"]]
ins.plot(range(3), fe, "o-", color=OR); ins.set_xticks(range(3)); ins.set_ylim(0, 0.25)
ins.set_title("per triangle: flux vs\nhalf image area, mismatch", fontsize=7.5); ins.tick_params(labelsize=7.5); ins.set_xlabel("level", fontsize=7.5)
d.set_title("D. Going round loops: the law's phase is\nthe Hopf phase, half the area, loop by loop")

fig.tight_layout()
fig.savefig("figures/f24_light_as_loops.png", bbox_inches="tight")
