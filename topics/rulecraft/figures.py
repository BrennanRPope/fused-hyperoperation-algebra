"""Figures for the Rulecraft x Fused Hyperoperation Algebra page. Run: python3 figures.py"""
import os
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from math import comb
import numpy as np

OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "figures")
BLUE, ORANGE, AQUA, GRAY, INK, INK2 = "#257de0", "#f47230", "#00b894", "#c9c9c9", "#191919", "#666666"
plt.rcParams.update({"font.size": 11, "axes.edgecolor": "#bbbbbb", "axes.labelcolor": INK2,
                     "xtick.color": INK2, "ytick.color": INK2, "axes.spines.top": False,
                     "axes.spines.right": False})

def box(ax, x, y, w, h, fc, label=None, tc=INK, fs=10):
    ax.add_patch(Rectangle((x, y), w, h, facecolor=fc, edgecolor="white", linewidth=2))
    if label: ax.text(x + w/2, y + h/2, label, ha="center", va="center", color=tc, fontsize=fs)

# r1: the boundary square, river numbers
fig, axs = plt.subplots(1, 2, figsize=(11, 4.6), gridspec_kw={"width_ratios": [1, 1.15]})
ax = axs[0]; XS, XT = 0.4, 0.96
box(ax, 0, 0, 1, 1, "#eeeeee", "1\n(no change)")
box(ax, 1, 0, XS, 1, BLUE, f"inside\n{XS}", "white", fs=9)
box(ax, 0, 1, 1, XT, AQUA, f"outside\nX_T = {XT}", "white")
box(ax, 1, 1, XS, XT, ORANGE, f"cross\n{XS*XT:.3f}", "white", fs=9)
ax.set_xlim(-0.05, 1.5); ax.set_ylim(-0.05, 2.05); ax.set_aspect("equal"); ax.axis("off")
ax.set_title("The boundary square: (1+X_S)(1+X_T)\nthe orange corner is the spillover", fontsize=11, color=INK)
ax = axs[1]
parts = [("Upper itself", 4.0, BLUE), ("Middle town", 9.6 - 4.0, AQUA), ("Lower town", 17.44 - 4.0 - 5.6, ORANGE)]
# Upper plant's value split by beneficiary (everyone built)
up = {"Upper": 4.0, "Middle": 0.4*1.4*10, "Lower": 0.4*1.4*1.4*10}
left = 0
for (name, v), c in zip(up.items(), (BLUE, AQUA, ORANGE)):
    ax.barh([0], [v], left=left, color=c, edgecolor="white", linewidth=2, height=0.5)
    ax.text(left + v/2, 0, f"{name}\n{v:.2f}", ha="center", va="center", color="white", fontsize=10)
    left += v
ax.axvline(6, color=INK2, linestyle="--", linewidth=1.5)
ax.text(6, 0.42, "plant cost 6", ha="center", color=INK2, fontsize=10)
ax.set_xlim(0, 18.5); ax.set_ylim(-0.6, 0.7); ax.set_yticks([])
ax.set_xlabel("value of the Upper town's plant (everyone built), by who gains")
ax.set_title("Upper sees 4 of 17.44: 77% spills downstream", fontsize=11, color=INK)
fig.tight_layout(); fig.savefig(f"{OUT}/r1_boundary.png", dpi=150); plt.close(fig)

# r2: panel grids
fig, axs = plt.subplots(1, 2, figsize=(11, 5.2))
for ax, (title, mem) in zip(axs, [("Expert + resident: area 100", [(1, 9), (9, 1)]),
                                   ("Two experts: area 36", [(1, 9), (1, 9)])]):
    ex = [m[0] for m in mem]; ins = [m[1] for m in mem]
    x0 = 0
    for a, ea in enumerate(ex):
        y0 = 0
        for b, ib in enumerate(ins):
            own = a == b
            box(ax, x0, y0, ea, ib, BLUE if own else ORANGE,
                f"{ea}x{ib}={ea*ib}" if min(ea, ib) >= 2 else None, "white")
            y0 += ib
        x0 += ea
    ax.set_xlim(-0.3, 10.3); ax.set_ylim(-0.3, 18.3); ax.set_aspect("equal")
    ax.set_xlabel("experience (summed)"); ax.set_ylabel("intellect (summed)")
    ax.set_title(title, fontsize=11, color=INK)
axs[1].text(2.8, 9, ha="left", s= "blue: each member's own corner\norange: cross corners, one member's\nexperience joined to another's intellect",
            fontsize=9.5, color=INK2)
fig.tight_layout(); fig.savefig(f"{OUT}/r2_panels.png", dpi=150); plt.close(fig)

# r3: event horizon
fig, axs = plt.subplots(1, 2, figsize=(11, 4.4))
ax = axs[0]
loss = np.linspace(0, 0.97, 400); back = loss/(1-loss)
ax.plot(loss*100, back*100, color=BLUE, linewidth=2)
ax.plot(loss*100, loss*100, color=GRAY, linewidth=2, linestyle="--")
for L in (0.2, 0.5, 0.9):
    ax.plot([L*100], [L/(1-L)*100], "o", color=BLUE, markersize=8, markeredgecolor="white")
    ax.annotate(f"lose {L:.0%}: need +{L/(1-L):.0%}", (L*100, L/(1-L)*100), xytext={0.2: (8, 330), 0.5: (30, 650), 0.9: (45, 1250)}[L], textcoords="data",
                fontsize=9.5, color=INK2)
ax.axvline(100, color=ORANGE, linewidth=2)
ax.text(88, 2600, "event horizon\n(-1 absorbs)", ha="right", color=INK, fontsize=10)
ax.set_xlim(0, 101); ax.set_ylim(0, 3200)
ax.set_xlabel("share lost (%)"); ax.set_ylabel("gain needed to get back (%)")
ax.set_title("The exact way back, -a/(1+a), blows up at -1", fontsize=11, color=INK)
ax.text(62, 160, "naive: same % back", color=INK2, fontsize=9.5)
ax = axs[1]
for k, (a, x) in enumerate([(0.0, 0), (-0.5, 1.7), (-0.9, 3.0), (-1.0, 3.9)]):
    side = 1 + a
    if side > 0: box(ax, x, 0, side, 1, BLUE)
    else: ax.plot([x, x], [0, 1], color=ORANGE, linewidth=3)
    ax.text(x + 0.25, -0.18, f"a = {a:g}\nside {side:g}", ha="center", va="top", fontsize=9.5, color=INK2)
ax.set_xlim(-0.1, 4.4); ax.set_ylim(-0.7, 1.2); ax.set_aspect("equal"); ax.axis("off")
ax.set_title("The side 1+a shrinks to zero: the box collapses\nand no multiplication can rebuild it", fontsize=11, color=INK)
fig.tight_layout(); fig.savefig(f"{OUT}/r3_horizon.png", dpi=150); plt.close(fig)

# r4: debt trap by cube layer
x = 0.02/1.03; years = np.arange(1, 31)
L1 = np.array([comb(n, 1)*x for n in years])*100
L2 = np.array([comb(n, 2)*x**2 for n in years])*100
L3 = np.array([sum(comb(n, k)*x**k for k in range(3, n+1)) for n in years])*100
fig, ax = plt.subplots(figsize=(9, 4.4))
ax.bar(years, L1, color=BLUE, edgecolor="white", linewidth=1, label="layer 1: simple (n·x)")
ax.bar(years, L2, bottom=L1, color=AQUA, edgecolor="white", linewidth=1, label="layer 2: pairs of years")
ax.bar(years, L3, bottom=L1+L2, color=ORANGE, edgecolor="white", linewidth=1, label="layer 3 and up")
ax.text(30, L1[-1]+L2[-1]+L3[-1]+1.5, "+78.1", ha="center", fontsize=10, color=INK)
ax.set_xlabel("years with no primary surplus (r = 5%, g = 3%)")
ax.set_ylabel("rise in debt, points of GDP (start 100)")
ax.set_title("A debt trap, split into cube layers: 58 points simple, 20 from compounding", fontsize=11, color=INK)
ax.legend(frameon=False, loc="upper left")
fig.tight_layout(); fig.savefig(f"{OUT}/r4_debt_layers.png", dpi=150); plt.close(fig)

# r5: tension held in balance
fig, ax = plt.subplots(figsize=(10, 3.8))
s = 5; x0 = 0
for d in (0, 2, 4):
    w, h = s + d, s - d
    box(ax, x0, 0, w, h, BLUE, f"{w} x {h} = {w*h}", "white")
    if d:
        box(ax, x0 + w - d, h, d, d, ORANGE, f"{d*d}", "white")
    ax.text(x0 + w/2, -0.8, "balanced" if d == 0 else f"lopsided by {d}: loses {d}² = {d*d}",
            ha="center", fontsize=10, color=INK2)
    x0 += w + 1.5
ax.set_xlim(-0.3, x0); ax.set_ylim(-1.6, 6); ax.set_aspect("equal"); ax.axis("off")
ax.set_title("Two principles sharing 10 units: the corner is largest when held in balance, (s+d)(s-d) = s² - d²",
             fontsize=11, color=INK)
fig.tight_layout(); fig.savefig(f"{OUT}/r5_balance.png", dpi=150); plt.close(fig)
print("ok")
