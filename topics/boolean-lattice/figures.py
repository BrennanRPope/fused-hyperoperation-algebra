"""Figures for the Boolean lattice page. Run: python3 figures.py (writes fig*.png here)."""
import os
from itertools import combinations, permutations
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
L = "abcd"
RANKC = ["#9e9e9e", "#3b7dd8", "#2ca58d", "#e0a020", "#d0453b"]

def label(S):
    return "".join(L[i] for i in S) if S else "1"

def subsets(n):
    return [S for m in range(n + 1) for S in combinations(range(n), m)]

def pos(S, n):
    if n == 2:
        v = {0: np.array([-1.0, 1.0]), 1: np.array([1.0, 1.0])}
    elif n == 3:
        v = {0: np.array([-1.3, 1.0]), 1: np.array([0.0, 1.0]), 2: np.array([1.3, 1.0])}
    else:
        v = {0: np.array([-2.4, 1.0]), 1: np.array([-0.8, 1.0]), 2: np.array([0.8, 1.0]), 3: np.array([2.4, 1.7])}
    p = np.array([0.0, 0.0])
    for i in S:
        p = p + v[i]
    return p

def hasse(ax, n, weight=None, title=""):
    els = subsets(n)
    for S in els:
        for T in els:
            if len(T) == len(S) + 1 and set(S) <= set(T):
                a, b = pos(S, n), pos(T, n)
                ax.plot([a[0], b[0]], [a[1], b[1]], color="#bbbbbb", lw=1.0, zorder=1)
    for S in els:
        x, y = pos(S, n)
        r = len(S)
        txt = label(S)
        if weight is not None and S:
            w = weight(r)
            txt = f"{txt}\n×{w}"
            alpha = 1.0 if w != "0" else 0.25
        else:
            alpha = 1.0
        ax.scatter([x], [y], s=900 if n == 4 else 1100, color=RANKC[r], alpha=alpha, zorder=2, edgecolor="white")
        ax.text(x, y, txt, ha="center", va="center", fontsize=8 if n == 4 else 10, color="white", zorder=3, fontweight="bold")
    ax.set_title(title, fontsize=11)
    ax.set_aspect("equal"); ax.axis("off")
    ax.margins(0.12)

# Figure 1: the shape at n = 2, 3, 4
fig, axs = plt.subplots(1, 3, figsize=(16, 6.4), gridspec_kw={"width_ratios": [1, 1.3, 1.9]})
hasse(axs[0], 2, title="n = 2: a square\na ∘ₖ b = a + b + k·ab")
hasse(axs[1], 3, title="n = 3: a cube\n7 non-empty corners")
hasse(axs[2], 4, title="n = 4: the tesseract (Boolean lattice B₄)\n15 non-empty corners")
handles = [plt.Line2D([], [], marker="o", ls="", ms=11, color=RANKC[r], label=f"rank {r}: {['the 1 (empty corner)', 'single steps', 'pairs', 'triples', 'all four'][r]}") for r in range(5)]
fig.legend(handles=handles, loc="lower center", ncol=5, frameon=False, fontsize=10)
fig.suptitle("The shape: every term of a chain a₁ ∘ₖ a₂ ∘ₖ … ∘ₖ aₙ sits on one corner of the n-cube (each line adds one step)", fontsize=13)
fig.subplots_adjust(bottom=0.12, top=0.84)
fig.savefig(os.path.join(HERE, "fig1-shape-n2-n3-n4.png"), dpi=150)
plt.close(fig)

# Figure 2: the dial k weights the ranks (n = 3)
def kw(k):
    def f(r):
        w = k ** (r - 1)
        if k == 0:
            return "1" if r == 1 else "0"
        if isinstance(w, float) and w != int(w):
            return f"{w:g}"
        return f"{int(w):+d}" if k < 0 else f"{int(w)}"
    return f
fig, axs = plt.subplots(1, 4, figsize=(17, 5.6))
for ax, k, t in zip(axs, [1, 0.5, 0, -1], ["k = 1: every corner once\n(ring theory's circle operation)",
                                          "k = 1/2: deeper corners fade\n(weight k^(rank − 1))",
                                          "k = 0: only single steps survive\n(plain addition)",
                                          "k = −1: signs alternate by rank\n(inclusion–exclusion)"]):
    hasse(ax, 3, weight=kw(k), title=t)
fig.suptitle("The dial k does not change the shape. It puts a weight on each layer: a corner of rank m counts k^(m − 1).", fontsize=13)
fig.subplots_adjust(top=0.80)
fig.savefig(os.path.join(HERE, "fig2-dial-k-weights.png"), dpi=150)
plt.close(fig)

# Figure 3: the plane's two dials weight the ranks by Lucas numbers
def lucas(k, d, m):
    u0, u1 = 0, 1
    if m == 0:
        return 0
    for _ in range(m - 1):
        u0, u1 = u1, k * u1 - d * u0
    return u1
cases = [((1, 0), "(1, 0)  circle operation"), ((0, -1), "(0, −1)  velocity addition"),
         ((0, 1), "(0, 1)  tangent of a sum"), ((1, -1), "(1, −1)  Fibonacci"), ((2, 1), "(2, 1)  U_m = m")]
fig, ax = plt.subplots(figsize=(13, 5.6))
ranks = range(1, 7)
for row, ((k, d), name) in enumerate(cases):
    for m in ranks:
        num = lucas(k, d, m)
        den = -d * lucas(k, d, m - 1) if m >= 2 else 0
        for col_off, val, lab in [(0, num, "top"), (0.5, den, "bottom")]:
            x = (m - 1) * 1.2 + col_off
            c = "#3b7dd8" if val > 0 else ("#d0453b" if val < 0 else "#eeeeee")
            ax.add_patch(plt.Rectangle((x, -row - 0.45), 0.48, 0.9, color=c, alpha=0.85 if val else 1))
            ax.text(x + 0.24, -row, f"{val:+d}" if val else "0", ha="center", va="center", fontsize=10,
                    color="white" if val else "#888888", fontweight="bold")
    ax.text(-0.25, -row, f"(k, d) = {name}", ha="right", va="center", fontsize=11)
for m in ranks:
    ax.text((m - 1) * 1.2 + 0.49, 0.85, f"rank {m}", ha="center", fontsize=11)
    ax.text((m - 1) * 1.2 + 0.24, 0.55, "top", ha="center", fontsize=8, color="#555555")
    ax.text((m - 1) * 1.2 + 0.74, 0.55, "bottom", ha="center", fontsize=8, color="#555555")
ax.set_xlim(-5.6, 7.3); ax.set_ylim(-len(cases) + 0.3, 1.2); ax.axis("off")
ax.set_title("Two dials, same cube: a chain of the plane's law (x + y + k·xy)/(1 − d·xy) is a ratio of corner sums.\n"
             "A rank-m corner gets weight U_m(k, d) on top and −d·U_(m−1)(k, d) on the bottom (Lucas numbers; blue +, red −).",
             fontsize=12)
fig.tight_layout()
fig.savefig(os.path.join(HERE, "fig3-two-dials-lucas.png"), dpi=150)
plt.close(fig)

# Figure 4: chains through the cube are the orderings (permutohedron)
fig, axs = plt.subplots(1, 2, figsize=(15, 6.4), gridspec_kw={"width_ratios": [1, 1.1]})
ax = axs[0]
hasse(ax, 3, title="A path from the bottom to the top of the cube\nadds the steps in one order: here b, then a, then c")
path = [(), (1,), (0, 1), (0, 1, 2)]
for S, T in zip(path, path[1:]):
    a, b = pos(S, 3), pos(T, 3)
    ax.annotate("", xy=b, xytext=a, arrowprops=dict(arrowstyle="-|>", color="#222222", lw=2.6, shrinkA=17, shrinkB=17), zorder=4)
ax = axs[1]
orders = ["abc", "bac", "bca", "cba", "cab", "acb"]
ang = np.pi / 2 + np.arange(6) * np.pi / 3
P = np.c_[np.cos(ang), np.sin(ang)] * 2.2
for i in range(6):
    a, b = P[i], P[(i + 1) % 6]
    ax.plot([a[0], b[0]], [a[1], b[1]], color="#888888", lw=1.6)
    o1, o2 = orders[i], orders[(i + 1) % 6]
    sw = [ch for ch in o1 if o1.index(ch) != o2.index(ch)]
    mid = (a + b) / 2 * 1.2
    ax.text(mid[0], mid[1], f"swap {sw[0]},{sw[1]}", ha="center", va="center", fontsize=9, color="#555555")
for i, o in enumerate(orders):
    ax.scatter([P[i, 0]], [P[i, 1]], s=1400, color="#3b7dd8", edgecolor="white", zorder=3)
    ax.text(P[i, 0], P[i, 1], o, ha="center", va="center", color="white", fontsize=11, fontweight="bold", zorder=4)
ax.text(0, 0.25, "permutohedron\n(3! = 6 orderings)", ha="center", va="center", fontsize=11)
ax.text(0, -0.55, "with numbers, every vertex\ngives the same total;\nwith matrices, an edge changes it\nby k times a commutator", ha="center", va="center", fontsize=9, color="#444444")
ax.set_aspect("equal"); ax.axis("off"); ax.margins(0.18)
ax.set_title("The 6 bottom-to-top paths are the 6 orderings,\nthe corners of a hexagon (an edge swaps two neighbours)", fontsize=11)
fig.suptitle("When order matters the terms still sit on the cube's corners. The orderings are the cube's paths, the permutohedron.", fontsize=13)
fig.subplots_adjust(top=0.82)
fig.savefig(os.path.join(HERE, "fig4-paths-permutohedron.png"), dpi=150)
plt.close(fig)
print("wrote fig1..fig4")
