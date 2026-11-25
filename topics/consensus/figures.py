"""Figures for the consensus page. Run: python3 figures.py  (writes figures/*.png)"""
import math
import os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle, Polygon
from scipy.stats import binom

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "figures")
os.makedirs(OUT, exist_ok=True)
SKY, INK, ORG, GREY, PALE = "#4fb3e8", "#10223a", "#ff7a1a", "#8a97a8", "#dcecf7"
plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 12, "axes.edgecolor": INK,
                     "axes.labelcolor": INK, "xtick.color": INK, "ytick.color": INK,
                     "axes.spines.top": False, "axes.spines.right": False})


def save(fig, name):
    fig.savefig(os.path.join(OUT, name), dpi=150, bbox_inches="tight", facecolor="white")
    plt.close(fig)


# f1  BFT: the quorum region in the square (quorum Q, faulty f)
fig, ax = plt.subplots(figsize=(6, 5))
Q = np.linspace(0, 1, 400)
ax.fill_between(Q, 0, np.minimum(2 * Q - 1, 1 - Q).clip(0), color=SKY, alpha=.55, lw=0)
ax.plot(Q, 2 * Q - 1, color=INK, lw=1.6)
ax.plot(Q, 1 - Q, color=INK, lw=1.6, ls="--")
ax.plot([2 / 3], [1 / 3], "o", color=ORG, ms=10, zorder=5)
ax.annotate("best corner: Q = 2/3, f = 1/3", (2 / 3, 1 / 3), (0.03, 0.38), color=INK,
            arrowprops=dict(arrowstyle="-", color=ORG))
ax.text(0.40, 0.93, "solid: safety, overlap 2Q − 1 > f", color=INK, fontsize=11)
ax.text(0.03, 0.12, "dashed: liveness, Q ≤ 1 − f", color=INK, fontsize=11)
ax.text(0.62, 0.07, "safe and live", color=INK, fontsize=11, weight="bold")
ax.set_xlim(0, 1); ax.set_ylim(0, 1)
ax.set_xlabel("quorum size Q (share of validators)"); ax.set_ylabel("faulty share f")
save(fig, "f1_bft_region.png")

# f2  committee: layer masses of a 100-cube, failure shaded
c, f = 100, 0.25
j = np.arange(c + 1)
pm = binom.pmf(j, c, f)
fig, ax = plt.subplots(figsize=(6.4, 4.2))
cols = [ORG if k >= math.ceil(c / 3) else SKY for k in j]
ax.bar(j, pm, color=cols, width=1.0, lw=0)
ax.axvline(c / 3, color=INK, lw=1.2, ls="--")
ax.set_xlim(5, 50)
ax.set_xlabel("layer j = faulty members in a committee of 100")
ax.set_ylabel("mass of layer j")
ax.text(34.5, pm.max() * .75, f"committee fails:\n{binom.sf(33, c, f):.1%} of the cube", color=ORG, fontsize=11)
ax.text(6, pm.max() * .95, "population 25% faulty", color=INK, fontsize=11)
save(fig, "f2_committee_layers.png")

# f3  coin age: attacker share of the next 20 blocks (from checks.py C5 output)
rows = {"plain stake": [0.101, 0.195, 0.295], "coin age s·(1+t)": [0.389, 0.532, 0.642],
        "fused (1+s)^(1+t) − 1": [0.415, 0.682, 0.860]}
x = np.arange(3)
fig, ax = plt.subplots(figsize=(6.4, 4.2))
for i, (name, vals) in enumerate(rows.items()):
    ax.bar(x + (i - 1) * 0.27, vals, 0.25, label=name, color=[GREY, SKY, ORG][i])
ax.axhline(0.5, color=INK, lw=1.2, ls="--")
ax.text(2.42, 0.52, "half", color=INK, fontsize=11)
ax.set_xticks(x, ["10% stake", "20% stake", "30% stake"])
ax.set_ylabel("attacker's share of next 20 blocks")
ax.set_ylim(0, 1)
ax.legend(frameon=False, fontsize=10, loc="upper left")
save(fig, "f3_coin_age.png")

# f4  Sybil corner: (1+s1)(1+s2) = 1 + s1 + s2 + s1 s2
fig, ax = plt.subplots(figsize=(5.2, 5.2))
s1, s2 = 0.9, 0.6
ax.add_patch(Rectangle((0, 0), 1, 1, color=PALE))
ax.add_patch(Rectangle((1, 0), s1, 1, color=SKY, alpha=.75))
ax.add_patch(Rectangle((0, 1), 1, s2, color=SKY, alpha=.75))
ax.add_patch(Rectangle((1, 1), s1, s2, color=ORG))
for (xx, yy, t) in [(.5, .5, "1"), (1 + s1 / 2, .5, "s₁"), (.5, 1 + s2 / 2, "s₂"), (1 + s1 / 2, 1 + s2 / 2, "s₁·s₂\nSybil gain")]:
    ax.text(xx, yy, t, ha="center", va="center", color=INK, fontsize=14 if len(t) < 4 else 12, weight="bold")
ax.set_xlim(-0.05, 2.0); ax.set_ylim(-0.05, 1.75); ax.set_aspect("equal"); ax.axis("off")
save(fig, "f4_sybil_corner.png")

# f5  energy, log scale
fig, ax = plt.subplots(figsize=(6.4, 3.4))
names = ["Bitcoin (work), 2025", "Ethereum (work), 2022", "Ethereum (stake), 2022"]
vals = [138.0, 23.0, 0.0026]
ax.barh(names[::-1], vals[::-1], color=[SKY, GREY, ORG][::-1])
ax.set_xscale("log"); ax.set_xlim(1e-3, 1e3)
for i, v in enumerate(vals[::-1]):
    ax.text(v * 1.3, i, f"{v:g} TWh/yr", va="center", color=INK, fontsize=11)
ax.set_xlabel("electricity per year (TWh, log scale)")
save(fig, "f5_energy.png")

# f6  churn: share of an old validator set that has left (k = -1 power), vs 1/3
staked, cap = 36.0e6, 256.0
cc = cap / staked
days = np.linspace(0, 400, 400)
t = days * 1440 / 6.4
fig, ax = plt.subplots(figsize=(6.4, 3.8))
ax.plot(days, np.minimum(t * cc, 1), color=SKY, lw=2.2, label="fixed exit cap (linear)")
ax.plot(days, 1 - (1 - cc) ** t, color=ORG, lw=2.2, label="proportional churn: 1 − (1 − c)^t")
ax.axhline(1 / 3, color=INK, ls="--", lw=1.2)
ax.text(5, 0.345, "1/3 of an old set gone: checkpoints needed", color=INK, fontsize=10.5)
ax.set_xlabel("days"); ax.set_ylabel("share of old validator set that has left")
ax.set_ylim(0, 0.6); ax.legend(frameon=False, fontsize=10, loc="lower right")
save(fig, "f6_churn.png")
print("wrote", sorted(os.listdir(OUT)))
