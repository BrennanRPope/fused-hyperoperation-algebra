"""Figure 25: the Lee-Yang near zone boxed for every roundness 5 <= p <= 8, shown with the 1.1-5 run (run after
near_zone_p.py with NZ_RANGE=5-8)."""
import json
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from matplotlib.collections import PatchCollection
from scipy.special import roots_jacobi
D1 = json.load(open("near_zone_p.json")); D2 = json.load(open("near_zone_p58.json"))
FAR = json.load(open("far_bounds_p.json"))["F2_strips"]
BLUE, ORANGE, GREY, RED, GREEN = "#2b6cb0", "#dd6b20", "#888888", "#c53030", "#2f855a"
plt.rcParams.update({"font.size": 11, "axes.spines.top": False, "axes.spines.right": False})
fig = plt.figure(figsize=(15.5, 11.5))
gs = fig.add_gridspec(2, 2, wspace=0.22, hspace=0.34)

# A: every certified box, both runs
ax = fig.add_subplot(gs[0, :])
boxes = D1["boxes"] + D2["boxes"]
zs = [b["zeros"] for b in boxes]; zmin, zmax = min(zs), max(zs)
cmap = plt.get_cmap("viridis")
pats, cols = [], []
for b in boxes:
    pats.append(Rectangle((b["p"][0], b["b"][0]), b["p"][1] - b["p"][0], b["b"][1] - b["b"][0]))
    cols.append(cmap((b["zeros"] - zmin)/max(1, zmax - zmin)) if b["ok"] else RED)
ax.add_collection(PatchCollection(pats, facecolor=cols, edgecolor="none"))
pp = np.linspace(1.1, 8, 50)
ax.plot(pp, pp - 1, color="k", lw=1.6)
ax.fill_between(pp, pp - 1, 7.4, color=GREY, alpha=0.12)
for pe in (4, 6, 8):
    ax.plot([pe, pe], [pe - 1, 7.4], color=RED, lw=3, solid_capstyle="butt")
ax.text(1.4, 5.2, "a = b + 1 > p: the rule says zeros turn complex here,\nso no climb is needed. Except at even p (red lines),\nwhere the rule says real zeros for every a: open\nbetween whole numbers", fontsize=10, color=GREY)
ax.text(6.2, 6.0, "b = p − 1  (a = p)", rotation=27, fontsize=10)
ax.axvline(5, color="k", lw=0.8, ls="--")
ax.text(5.03, 0.15, "new run →", fontsize=10)
ax.text(4.97, 0.15, "← earlier run", fontsize=10, ha="right")
ax.set_xlim(1.1, 8); ax.set_ylim(0, 7.15)
ax.set_xlabel("roundness p"); ax.set_ylabel("step b (the climb goes b → a = b + 1)")
sm = plt.cm.ScalarMappable(cmap=cmap, norm=plt.Normalize(zmin, zmax)); sm.set_array([])
cb = fig.colorbar(sm, ax=ax, pad=0.01, fraction=0.03); cb.set_label("zeros isolated in the box's s-range")
n2 = len(D2["boxes"]); nbad = sum(not b["ok"] for b in boxes)
ax.set_title(f"All {len(boxes):,} boxes pass, {n2:,} of them new for 5 ≤ p ≤ 8"
             + (f" ({nbad} fail, red)" if nbad else "") + ": zeros isolated, next step's sign alternates")

# B: one hard box near p = 8, top step
ax = fig.add_subplot(gs[1, 0])
def G(p, c, s, N=500):
    x, w = roots_jacobi(N, c, 0.0); x = (x + 1)/2; w = w/2**(c + 1)   # weight (1 - x)^c on [0, 1]
    f = ((1 - x**p)/(1 - x))**c
    return (w*f*np.cos(np.outer(s, x))).sum(1)
p, b = 7.95, 6.925
s = np.linspace(0.05, 110, 3000)
gb, ga = G(p, b, s), G(p, b + 1, s)
from scipy.ndimage import maximum_filter1d
def env(g):  # local size of the wiggles, over a window of about one period
    return np.maximum(maximum_filter1d(np.abs(g), size=int(2*np.pi/(s[1] - s[0]))), 1e-300)
scb, sca = gb/env(gb), ga/env(ga)
ax.plot(s, scb, color=BLUE, lw=1.8, label=f"G at step b = {b}")
ax.plot(s, sca, color=ORANGE, lw=1.8, label=f"G at a = b + 1 = {b + 1}")
zi = np.where(np.diff(np.sign(gb)) != 0)[0]
for i in zi:
    z = s[i]; sg = np.sign(ga[i])
    ax.plot(z, 0, "o", color=BLUE, ms=4)
    ax.annotate("+" if sg > 0 else "−", (z, sca[i]), textcoords="offset points",
                xytext=(0, 6 if sg > 0 else -14), ha="center", color=ORANGE, fontsize=11, fontweight="bold")
ax.axhline(0, color="k", lw=0.8)
ax.set_xlim(0, 110); ax.set_xlabel("s"); ax.set_ylabel("G divided by its local wiggle size")
ax.set_ylim(-1.35, 1.6)
ax.legend(fontsize=9.5, loc="upper right")
ax.set_title(f"p = {p}, top step: {len(zi)} zeros by s = 110, and at each one\nthe orange curve's sign alternates + − + − …")

# C: near and far zones meet at every p
ax = fig.add_subplot(gs[1, 1])
for k, st in enumerate(FAR):
    ax.fill_between(st["p"], [st["T"]]*2, 230, color=GREEN, alpha=0.25, lw=0, label="far zone, proved by formulas in p (s ≥ T*)" if k == 0 else None)
    ax.plot(st["p"], [st["T"]]*2, color=GREEN, lw=1.4)
TN = {(1.1, 2): 42, (2, 3): 31, (3, 4): 50, (4, 5): 76, (5, 6): 112, (6, 7): 152, (7, 8): 205}
for k, ((lo, hi), t) in enumerate(TN.items()):
    ax.fill_between([lo, hi], 0, t, color=BLUE, alpha=0.22, lw=0, label="near zone, boxed (0 ≤ s ≤ end of grid)" if k == 0 else None)
    ax.plot([lo, hi], [t, t], color=BLUE, lw=1.6)
ax.set_xlim(1.1, 8); ax.set_ylim(0, 230)
ax.set_xlabel("roundness p"); ax.set_ylabel("s")
ax.legend(fontsize=9.5, loc="upper left")
ax.set_title("The two zones overlap at every roundness from 1.1 to 8")

fig.suptitle("Figure 25. Lee–Yang near zone boxed out to roundness 8 (DRAFT, computer-assisted)", fontsize=14, y=0.995)
fig.savefig("figures/fig25-lee-yang-near-p8.png", dpi=110, bbox_inches="tight"); print("saved")
