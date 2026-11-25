"""Figures for the information-theory page. Run: python3 figures.py"""
import math, os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle, Polygon

OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "figures")
os.makedirs(OUT, exist_ok=True)
BLUE, ORANGE, AQUA, YELLOW, VIOLET, RED = "#257de0", "#f47230", "#00b894", "#f2a600", "#6a4ce2", "#ec4647"
INK, INK2, GRID = "#191919", "#666666", "#d9d9d9"
plt.rcParams.update({"font.size": 11, "text.color": INK, "axes.labelcolor": INK,
                     "xtick.color": INK2, "ytick.color": INK2, "axes.edgecolor": GRID,
                     "font.family": "DejaVu Sans"})

# ---------- f1: four ways to count the corner ----------
def square(ax, title, sub, cells, p=0.38, q=0.30, outer=False):
    ax.set_aspect("equal"); ax.axis("off")
    if outer:  # the circle-operation square: side 1+a by 1+b, unit cell at origin
        a, b = 0.55, 0.42
        rects = {"1": (0,0,1,1), "a": (1,0,a,1), "b": (0,1,1,b), "ab": (1,1,a,b)}
        lab = {"1": "1", "a": "a", "b": "b", "ab": "ab"}
        lim = (1+a, 1+b)
    else:
        rects = {"1": (p,q,1-p,1-q), "a": (0,q,p,1-q), "b": (p,0,1-p,q), "ab": (0,0,p,q)}
        lab = {"1": "both fine\n(1−p)(1−q)", "a": "p(1−q)", "b": "(1−p)q", "ab": "pq"}
        lim = (1, 1)
    for key, (x, y, w_, h_) in rects.items():
        style = cells.get(key)
        fc = "white"; hatch = None
        if style == "+": fc = BLUE
        elif style == "++": fc = VIOLET
        elif style == "-": fc = "white"
        ax.add_patch(Rectangle((x, y), w_, h_, fc=fc, ec="white" if fc != "white" else GRID, lw=2,
                               alpha=0.85 if fc != "white" else 1))
        txt = lab[key]
        if style == "++": txt += "\n(twice)"
        ax.text(x+w_/2, y+h_/2, txt, ha="center", va="center", fontsize=9,
                color="white" if fc != "white" else INK2)
    ax.set_xlim(-0.05, lim[0]+0.05); ax.set_ylim(-0.05, lim[1]+0.05)
    ax.set_title(title, fontsize=12, color=INK, pad=6, loc="left")
    ax.text(0, -0.13, sub, transform=ax.transAxes, fontsize=9.5, color=INK2, va="top")

fig, axs = plt.subplots(1, 4, figsize=(13, 4.4))
square(axs[0], "k = 1  circle operation", "(1+a)(1+b) − 1 = a + b + ab\nvolume grows: counting, odds, Campbell",
       {"a": "+", "b": "+", "ab": "+"}, outer=True)
square(axs[1], "k = 0  plain addition", "p + q: the corner counted twice\nShannon, Rényi, code lengths",
       {"a": "+", "b": "+", "ab": "++"})
square(axs[2], "k = −1  either one", "p + q − pq: everything but 'both fine'\nfailures, erasures, Gini, Kraft slack",
       {"a": "+", "b": "+", "ab": "+"})
square(axs[3], "k = −2  exactly one", "p + q − 2pq: the two off-diagonal cells\nbit flips, parity, check nodes",
       {"a": "+", "b": "+"})
fig.suptitle("One square, four ways to count its corner: the weight k in a ∘ₖ b = a + b + k·ab",
             x=0.01, ha="left", fontsize=13)
plt.tight_layout(rect=(0, 0.10, 1, 0.93)); fig.savefig(f"{OUT}/f1_corner_square.png", dpi=150); plt.close(fig)

# ---------- f2: the circle-operation cube, union vs parity ----------
def iso(x, y, z):
    return x * 1.0 + y * 0.55, z * 1.0 + y * 0.45
def cube(ax, title, color_of, sign_of, sub):
    ax.set_aspect("equal"); ax.axis("off")
    corners = [(x, y, z) for x in (0, 1) for y in (0, 1) for z in (0, 1)]
    for c in corners:
        for d in range(3):
            if c[d] == 0:
                c2 = list(c); c2[d] = 1
                (x1, y1), (x2, y2) = iso(*c), iso(*c2)
                ax.plot([x1, x2], [y1, y2], color=GRID, lw=2, zorder=1)
    for c in corners:
        g = sum(c); X, Y = iso(*c)
        col = color_of(g)
        ax.scatter([X], [Y], s=520, color=col if col else "white", edgecolor=GRID if not col else "white",
                   lw=2, zorder=3)
        name = "".join(s for s, v in zip("pqr", c) if v) or "1"
        ax.text(X, Y, name, ha="center", va="center", fontsize=9, zorder=4,
                color="white" if col else INK2)
        s = sign_of(g)
        if s: ax.text(X - 0.2, Y + 0.13, s, fontsize=12, color=INK, zorder=4, ha="center")
    ax.set_title(title, fontsize=12, loc="left")
    ax.text(0, -0.04, sub, transform=ax.transAxes, fontsize=9.5, color=INK2, va="top")
    ax.set_xlim(-0.4, 1.9); ax.set_ylim(-0.6, 1.9)
fig, axs = plt.subplots(1, 2, figsize=(10, 4.6))
cube(axs[0], "Either part fails (k = −1)", lambda g: [None, BLUE, ORANGE, VIOLET][g],
     lambda g: ["", "+", "−", "+"][g],
     "Every corner but the empty one, signs alternating by layer:\np + q + r − pq − pr − qr + pqr  (inclusion–exclusion)")
cube(axs[1], "Odd number of flips (k = −2)", lambda g: [None, BLUE, None, VIOLET][g],
     lambda g: ["", "", "", ""][g],
     "Only the odd layers (1·3·3·1 → the 3 and the 1):\n1 − 2·P(odd) = (1 − 2p)³  (piling-up lemma)")
plt.tight_layout(); fig.savefig(f"{OUT}/f2_cube_union_parity.png", dpi=150); plt.close(fig)

# ---------- f3: the weight ladder ----------
fig, ax = plt.subplots(figsize=(11, 4.4))
ax.axhline(0, color=GRID, lw=2)
marks = [
    (1, "k = 1", "Tsallis q = 0 (count − 1)\nCampbell cost M − 1\nvariable node: odds − 1"),
    (0, "k = 0", "Shannon, Rényi\nexpected code length\nLLRs, hazards, bits"),
    (-1, "k = −1", "independent failures\nerasure channel\nGini–Simpson (Tsallis q = 2)\nKraft slack"),
    (-2, "k = −2", "binary symmetric channel\nparity, piling-up\ncheck node"),
]
for x, lab, txt in marks:
    ax.scatter([x], [0], s=120, color=BLUE, zorder=3, edgecolor="white", lw=2)
    if x == -1:
        ax.text(x, 0.08, txt, ha="center", va="bottom", fontsize=9.2, color=INK2)
        ax.text(x + 0.1, -0.1, lab, ha="left", va="top", fontsize=11, color=INK)
    else:
        ax.text(x, 0.08, lab, ha="center", fontsize=11, color=INK)
        ax.text(x, -0.14, txt, ha="center", va="top", fontsize=9.2, color=INK2)
ax.scatter([-1.5, -4/3], [0, 0], s=60, color=AQUA, zorder=3, edgecolor="white", lw=2)
ax.text(-1.5, 0.08, "−3/2: 3 letters", ha="center", va="bottom", fontsize=9, color=INK2)
ax.text(-4/3, -0.1, "−4/3: 4 letters\n(depolarizing qubit)", ha="center", va="top", fontsize=9, color=INK2)
ax.annotate("", xy=(-2.15, 0.85), xytext=(1.15, 0.85),
            arrowprops=dict(arrowstyle="->", color=INK2, lw=1.2))
ax.text(-0.5, 0.9, "more cancelling in the corner", ha="center", fontsize=10, color=INK2)
ax.text(-1.85, 0.55, "channels: k = −(1 + c), c = chance a second error\nlands back on the right letter = 1/(letters − 1)",
        ha="center", fontsize=9, color=INK2)
ax.text(0.6, 0.5, "entropy: k = 1 − q\n(Tsallis index q)", ha="center", fontsize=9, color=INK2)
ax.set_xlim(-2.5, 1.5); ax.set_ylim(-0.75, 1.05); ax.axis("off")
ax.set_title("The weight ladder: where each information law sits", loc="left", fontsize=13)
plt.tight_layout(); fig.savefig(f"{OUT}/f3_weight_ladder.png", dpi=150); plt.close(fig)

# ---------- f4: decoding, two circle operations on one segment ----------
fig, ax = plt.subplots(figsize=(7.5, 5))
r1 = np.linspace(-0.999, 0.999, 400); r2 = 0.6
ax.plot(r1, r1 * r2, color=BLUE, lw=2)
ax.plot(r1, (r1 + r2) / (1 + r1 * r2), color=ORANGE, lw=2)
ax.plot(r1, r1, color=GRID, lw=1, ls="--")
ax.text(0.55, 0.22, "check node: ρ₁·ρ₂\n(scales toward the middle)", color=INK, fontsize=10)
ax.text(-0.95, 0.55, "variable node: (ρ₁+ρ₂)/(1+ρ₁ρ₂)\n(velocity addition, pushes out)", color=INK, fontsize=10)
ax.set_xlabel("incoming belief ρ₁ = 1 − 2p = tanh(LLR/2)")
ax.set_ylabel("outgoing belief ρ")
ax.set_title("Belief propagation, other input fixed at ρ₂ = 0.6", loc="left", fontsize=12)
ax.grid(color=GRID, lw=0.6); ax.set_xlim(-1, 1); ax.set_ylim(-1, 1)
for s in ("top", "right"): ax.spines[s].set_visible(False)
plt.tight_layout(); fig.savefig(f"{OUT}/f4_decoding.png", dpi=150); plt.close(fig)

# ---------- f5: data processing along a cascade ----------
def h2(p): return 0.0 if p <= 0 or p >= 1 else -p*math.log2(p) - (1-p)*math.log2(1-p)
fig, ax = plt.subplots(figsize=(7.5, 4.8))
n = np.arange(1, 21)
series = [(0.02, "BSC, p = 0.02", BLUE), (0.05, "BSC, p = 0.05", ORANGE), (0.10, "BSC, p = 0.10", AQUA)]
for p, lab, col in series:
    C = [1 - h2((1 - (1 - 2*p)**m) / 2) for m in n]
    ax.plot(n, C, color=col, lw=2, marker="o", ms=4, label=lab)
ax.plot(n, 0.9**n, color=VIOLET, lw=2, marker="o", ms=4, label="BEC, erasure 0.10")
ax.set_xlabel("links in the cascade"); ax.set_ylabel("capacity, bits per use")
ax.set_title("A cascade only loses: capacity of n identical links", loc="left", fontsize=12)
ax.legend(frameon=False, fontsize=9.5); ax.grid(color=GRID, lw=0.6); ax.set_ylim(0, 1)
for s in ("top", "right"): ax.spines[s].set_visible(False)
plt.tight_layout(); fig.savefig(f"{OUT}/f5_cascade_capacity.png", dpi=150); plt.close(fig)
print("figures written to", OUT)
