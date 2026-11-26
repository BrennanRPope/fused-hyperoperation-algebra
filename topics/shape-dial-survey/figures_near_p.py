"""Figure 24: the Lee-Yang near zone boxed for every roundness 1.1 <= p <= 5 (run after near_zone_p.py)."""
import json, math
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from matplotlib.collections import PatchCollection
D = json.load(open("near_zone_p.json")); FAR = json.load(open("far_bounds_p.json"))["F2_strips"]
BLUE, ORANGE, GREY, RED, GREEN, PURPLE = "#2b6cb0", "#dd6b20", "#888888", "#c53030", "#2f855a", "#6b46c1"
plt.rcParams.update({"font.size": 11, "axes.spines.top": False, "axes.spines.right": False})
fig = plt.figure(figsize=(15.5, 11.5))
gs = fig.add_gridspec(2, 2, wspace=0.22, hspace=0.34)

# A: every certified box in (p, b)
ax = fig.add_subplot(gs[0, :])
boxes = D["boxes"]
zs = [b["zeros"] for b in boxes]; zmin, zmax = min(zs), max(zs)
cmap = plt.get_cmap("viridis")
pats, cols = [], []
for b in boxes:
    pats.append(Rectangle((b["p"][0], b["b"][0]), b["p"][1] - b["p"][0], b["b"][1] - b["b"][0]))
    cols.append(cmap((b["zeros"] - zmin)/max(1, zmax - zmin)) if b["ok"] else RED)
ax.add_collection(PatchCollection(pats, facecolor=cols, edgecolor="white", linewidth=0.15))
pp = np.linspace(1.1, 5, 50)
ax.plot(pp, pp - 1, color="k", lw=1.6)
ax.fill_between(pp, pp - 1, 4.2, color=GREY, alpha=0.12)
ax.text(1.6, 2.6, "a = b + 1 > p: past the rule's edge,\nnot needed (the rule says the zeros\nturn complex up here)", fontsize=10, color=GREY)
ax.text(3.45, 2.6, "b = p − 1  (a = p)", rotation=17, fontsize=10)
split = [b for b in boxes if b["depth"] > 0]
ax.set_xlim(1.1, 5); ax.set_ylim(0, 4.05)
ax.set_xlabel("roundness p"); ax.set_ylabel("step b (the climb goes b → a = b + 1)")
sm = plt.cm.ScalarMappable(cmap=cmap, norm=plt.Normalize(zmin, zmax)); sm.set_array([])
cb = fig.colorbar(sm, ax=ax, pad=0.01, fraction=0.03); cb.set_label("zeros isolated in the box's s-range")
nbad = sum(not b["ok"] for b in boxes)
ax.set_title(f"Every one of {len(boxes):,} boxes passes ({len(split)} came from splitting a box in four"
             + (f"; {nbad} fail, red" if nbad else "") + "): zeros isolated, and the next step's sign alternates")
for x in (2, 3, 4):
    ax.axvline(x, color="white", lw=0.8, ls=":")

# B: one hard box, p near 5, top step
ax = fig.add_subplot(gs[1, 0])
xg, wg = np.polynomial.legendre.leggauss(600); xg = (xg + 1)/2; wg = wg/2
def G(p, c, s):
    return np.array([np.sum(wg*(1 - xg**p)**c*np.cos(si*xg)) for si in s])
p, b = 4.925, 3.625
s = np.linspace(0, 40, 1600)
gb, ga = G(p, b, s), G(p, b + 1, s)
sc = np.sign(s)*np.abs(s)**1.0
ax.plot(s, gb*(1 + s)**1.5, color=BLUE, lw=2, label=f"G at step b = {b}")
ax.plot(s, ga*(1 + s)**1.5, color=ORANGE, lw=2, label=f"G at a = b + 1 = {b + 1}")
zi = np.where(np.diff(np.sign(gb)) != 0)[0]
for i in zi:
    z = s[i]; sg = np.sign(ga[i])
    ax.plot(z, 0, "o", color=BLUE, ms=6)
    ax.annotate("+" if sg > 0 else "−", (z, ga[i]*(1 + z)**1.5), textcoords="offset points",
                xytext=(0, 8 if sg > 0 else -16), ha="center", color=ORANGE, fontsize=13, fontweight="bold")
ax.axhline(0, color="k", lw=0.8)
ax.set_xlim(0, 40); ax.set_xlabel("s"); ax.set_ylabel("G × (1 + s)^1.5  (scaled so the wiggles stay visible)")
ax.legend(fontsize=9.5, loc="upper right")
ax.set_title(f"p = {p}: at each zero of the blue curve the orange one\nhas a sign, and the signs alternate + − + − …")

# C: near and far zones meet at every p
ax = fig.add_subplot(gs[1, 1])
for k, st in enumerate(FAR):
    if st["p"][0] >= 5: break
    ax.fill_between(st["p"], [st["T"]]*2, 95, color=GREEN, alpha=0.25, lw=0, label="far zone, proved by formulas in p (s ≥ T*)" if k == 0 else None)
    ax.plot(st["p"], [st["T"]]*2, color=GREEN, lw=1.6)
TN = {(1.1, 2): 42, (2, 3): 31, (3, 4): 50, (4, 5): 76}
for k, ((lo, hi), t) in enumerate(TN.items()):
    ax.fill_between([lo, hi], 0, t, color=BLUE, alpha=0.22, lw=0, label="near zone, boxed here (0 ≤ s ≤ end of grid)" if k == 0 else None)
    ax.plot([lo, hi], [t, t], color=BLUE, lw=1.6)
ax.set_xlim(1.1, 5); ax.set_ylim(0, 95)
ax.set_xlabel("roundness p"); ax.set_ylabel("s")
ax.text(2.05, 91, "where the two overlap, each box finds\na zero of G past T*, so the alternation\ncarries from one zone into the other", fontsize=9.5, color=GREY, va="top")
ax.legend(fontsize=9.5, loc="lower right")
ax.set_title("The two zones overlap at every roundness from 1.1 to 5")

fig.suptitle("Figure 24. Lee–Yang near zone boxed for every roundness 1.1 ≤ p ≤ 5 (DRAFT, computer-assisted)", fontsize=14, y=0.995)
fig.savefig("figures/fig24-lee-yang-near-p.png", dpi=110, bbox_inches="tight"); print("saved")
