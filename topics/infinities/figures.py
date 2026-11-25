"""Figures for the infinities page.  python3 figures.py  -> figures/*.png"""
import math, sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch

sys.path.insert(0, "../fractional-rungs")
from regular_iteration import Es

BLUE = "#1d2bd8"; RED = "#cf3519"; INK = "#15131a"; MUTE = "#8a8794"; WASH = "#eef0ff"
plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 12, "axes.edgecolor": INK,
                     "axes.labelcolor": INK, "xtick.color": INK, "ytick.color": INK})

def save(fig, name):
    fig.savefig(f"figures/{name}", dpi=160, bbox_inches="tight", facecolor="white")
    plt.close(fig)

# ------------------------------------------------------------- 1. same count, different shape
def draw_tree(ax, children, depth, color_fn, title):
    # children(node) -> list of child nodes ; node = tuple of bits
    level = [()]
    pos = {(): 0.5}
    for d in range(depth):
        nxt = []
        for node in level:
            nxt += children(node, d)
        n = len(nxt)
        for i, c in enumerate(nxt):
            pos[c] = (i + 0.5)/n
        for c in nxt:
            p = c[:-1]
            ax.plot([pos[p], pos[c]], [-d, -(d+1)], color=color_fn(c), lw=1.3, solid_capstyle="round")
        level = nxt
    for node in level:
        ax.plot(pos[node], -depth, "o", ms=2.5, color=color_fn(node))
    ax.set_title(title, fontsize=13, color=INK, pad=8)
    ax.set_xlim(-0.02, 1.02); ax.set_ylim(-depth - 0.6, 0.4); ax.axis("off")

D = 8
def comb_children(node, d):
    if all(b == 0 for b in node):
        return [node + (0,), node + (1,)]
    return [node + (0,)]
split = {0, 1, 3, 7}          # uniform tree: every branch splits together, rarely
def uni_children(node, d):
    return [node + (0,), node + (1,)] if d in split else [node + (0,)]

fig, axs = plt.subplots(1, 2, figsize=(11, 5.2))
draw_tree(axs[0], comb_children, D, lambda c: BLUE if all(b == 0 for b in c) else MUTE,
          "Comb: one new dead-end per level\n1 + 1 + 1 + …  →  countable")
draw_tree(axs[1], uni_children, D, lambda c: RED,
          "Every branch keeps splitting\n(1+1)(1+1)(1+1)… − 1  →  continuum")
comb_counts = [d + 1 for d in range(D + 1)]
uni_counts = [2**sum(1 for k in split if k < d) for d in range(D + 1)]
for ax, cnt in [(axs[0], comb_counts), (axs[1], uni_counts)]:
    for d, c in enumerate(cnt):
        ax.text(1.06, -d, str(c), fontsize=10, color=INK, va="center", transform=ax.transData, clip_on=False)
    ax.text(1.06, 0.55, "pieces", fontsize=9, color=MUTE, ha="left", clip_on=False)
fig.text(0.5, 0.01, "Both trees have about as many pieces at every level. Only the shape decides the size of the limit.",
         ha="center", fontsize=11, color=MUTE)
save(fig, "f1_same_count.png")

# ------------------------------------------------------------- 2. three scales of size
fig, axs = plt.subplots(3, 1, figsize=(11, 6.4), gridspec_kw={"hspace": 0.95})
# a) cardinality ladder
ax = axs[0]
ax.axhline(0, color=INK, lw=1)
pts = [(0, r"$\aleph_0$"), (1.6, r"$\aleph_1$"), (2.6, r"$\aleph_2$"), (3.4, r"$\aleph_3$"), (4.0, ""), (4.45, ""), (4.75, ""), (6.2, r"$\mathfrak{c} = 2^{\aleph_0}$")]
for x, lab in pts:
    ax.plot(x, 0, "o", ms=9 if lab else 5, color=BLUE if x < 6 else RED)
    if lab: ax.text(x, 0.35, lab, ha="center", fontsize=15)
ax.text(5.3, -0.55, "how many steps in between: the axioms don't decide", ha="center", fontsize=10, color=MUTE)
ax.text(0.8, -0.55, "nothing in between", ha="center", fontsize=10, color=MUTE)
ax.set_xlim(-0.5, 7.2); ax.set_ylim(-0.9, 0.8); ax.axis("off")
ax.set_title("Cardinality: separate steps, never a smooth line", loc="left", fontsize=13)
# b) dimension
ax = axs[1]
xs = np.linspace(0, 1, 400)
ax.imshow(np.vstack([xs, xs]), extent=(0, 1, -0.12, 0.12), aspect="auto", cmap=matplotlib.colors.LinearSegmentedColormap.from_list("b", ["#dfe2ff", BLUE]))
for x, lab, dy in [(0, "countable sets\n(and some uncountable)", 0.32), (math.log(2)/math.log(3), "Cantor set\nlog 2 / log 3", 0.32),
                   (0.722, "density-1/5 corners\nH(1/5)", -0.62), (1, "full line", 0.32)]:
    ax.plot([x, x], [-0.16, 0.16], color=INK, lw=1.2)
    ax.text(x, dy, lab, ha="center", fontsize=9.5, va="bottom" if dy > 0 else "top")
ax.set_xlim(-0.08, 1.08); ax.set_ylim(-0.9, 0.9); ax.axis("off")
ax.set_title("Fractal dimension: smooth from 0 to 1, but only among sets that are already uncountable", loc="left", fontsize=13)
# c) growth of the finite stages
ax = axs[2]
ax.axhline(0, color=INK, lw=1)
marks = [(0, "n\nrung 0", BLUE), (1.2, "n²", MUTE), (2.0, "n¹⁰⁰", MUTE), (3.2, "s = ¼", RED), (4.1, "s = ½", RED),
         (5.0, "s = ¾", RED), (6.2, "2ⁿ\nrung 1", BLUE)]
for x, lab, c in marks:
    ax.plot(x, 0, "o", ms=8, color=c)
    ax.text(x, 0.3, lab, ha="center", va="bottom", fontsize=10.5)
ax.annotate("", xy=(5.6, -0.45), xytext=(2.6, -0.45), arrowprops=dict(arrowstyle="<->", color=RED))
ax.text(4.1, -0.62, "fractional rungs: a continuum of growth rates, all faster than every polynomial, all slower than 2ⁿ",
        ha="center", va="top", fontsize=9.5, color=RED)
ax.set_xlim(-0.5, 7.2); ax.set_ylim(-1.0, 0.9); ax.axis("off")
ax.set_title("Growth of the finite stages: the fused rungs fill the gap smoothly", loc="left", fontsize=13)
save(fig, "f2_three_scales.png")

# ------------------------------------------------------------- 3. zoomed mult-add sweeps dimension
def Dt(t):
    return (math.log1p(t) - (t*math.log(t) if t > 0 else 0)/(1+t))/math.log(2)
fig, axs = plt.subplots(1, 2, figsize=(11, 4.4), gridspec_kw={"width_ratios": [1.15, 1]})
ax = axs[0]
ts = np.linspace(0, 1, 400)
ax.plot(ts, [Dt(t) for t in ts], color=BLUE, lw=2.6)
for t in [0.1, 0.5]:
    ax.plot(t, Dt(t), "o", color=RED); ax.text(t + 0.03, Dt(t) - 0.07, f"t = {t}: dim {Dt(t):.2f}", fontsize=10)
ax.set_xlabel("zoom t in  a + b + t·ab     (0 = addition, 1 = mult-add)")
ax.set_ylabel("dimension of the corners it lands on")
ax.set_xlim(0, 1); ax.set_ylim(0, 1.05)
ax.spines[["top", "right"]].set_visible(False)
ax.set_title("D(t) = [ln(1+t) − t·ln t / (1+t)] / ln 2", fontsize=12)
ax = axs[1]
n = 60
for t, c in [(0.1, "#9aa3ff"), (0.4, "#5b67f0"), (1.0, BLUE)]:
    w = np.array([math.comb(n, k)*t**k for k in range(n+1)]); w = w/w.sum()
    ax.plot(np.arange(n+1)/n, w, color=c, lw=2, label=f"t = {t}")
    ax.axvline(t/(1+t), color=c, ls=":", lw=1)
ax.set_xlabel("share of switches on (grade / n), n = 60")
ax.set_ylabel("weight of that layer")
ax.legend(frameon=False)
ax.spines[["top", "right"]].set_visible(False)
ax.set_title("Weighted cube layers C(n,k)·tᵏ peak at t/(1+t)", fontsize=12)
save(fig, "f3_zoom_dimension.png")

# ------------------------------------------------------------- 4. fractional rung counts
def logN(s, n):
    c = float(Es(-s, np.array([1.0]))[0])
    return np.array([float(Es(s - 1, np.array([v*c]))[0]) for v in n])
ns = np.unique(np.round(np.logspace(0, 6, 60)).astype(int)).astype(float)
fig, ax = plt.subplots(figsize=(10, 5))
cols = {0.0: BLUE, 0.25: "#e48a76", 0.5: "#d85e44", 0.75: RED, 1.0: BLUE}
for s in [0.0, 0.25, 0.5, 0.75, 1.0]:
    y = logN(s, ns)/math.log(10)
    ax.plot(ns, y, color=cols[s], lw=2.6 if s in (0, 1) else 2, ls="-" if s in (0, 1) else "--")
    ax.text(ns[-1]*1.25, y[-1], {0.0: "rung 0: n", 1.0: "rung 1: 2ⁿ − 1"}.get(s, f"s = {s}"), va="center", fontsize=10.5,
            color=cols[s])
m = ns >= 2
for k, lab, xl in [(2, "n²", 30), (10, "n¹⁰", 30)]:
    ax.plot(ns[m], k*np.log10(ns[m]), color=MUTE, lw=1, ls=":")
    pass
ax.text(1.3, 3e4, "dotted grey: n² (lower) and n¹⁰ (upper)", color=MUTE, fontsize=10)
ax.set_xscale("log"); ax.set_yscale("log")
ax.set_xlabel("n (how many ones are fused together)")
ax.set_ylabel("digits in the count  (log₁₀ of 1 + count)")
ax.set_xlim(1, 1e6*6)
ax.spines[["top", "right"]].set_visible(False)
ax.set_title("Fusing n ones at a fractional rung: each rung eventually outgrows every polynomial, and always loses to 2ⁿ", fontsize=12)
save(fig, "f4_rung_counts.png")

# ------------------------------------------------------------- 5. Cichon's diagram
fig, ax = plt.subplots(figsize=(10, 4.6))
nodes = {
    r"$\aleph_1$": (-0.6, 0), "add(N)": (0.6, 0), "add(M)": (2.2, 0), "cov(M)": (3.8, 0), "non(N)": (5.4, 0),
    "cov(N)": (0.6, 2.4), "non(M)": (2.2, 2.4), "cof(M)": (3.8, 2.4), "cof(N)": (5.4, 2.4), r"$\mathfrak{c}$": (6.6, 2.4),
    r"$\mathfrak{b}$": (2.2, 1.2), r"$\mathfrak{d}$": (3.8, 1.2),
}
edges = [(r"$\aleph_1$", "add(N)"), ("add(N)", "add(M)"), ("add(M)", "cov(M)"), ("cov(M)", "non(N)"),
         ("add(N)", "cov(N)"), ("cov(N)", "non(M)"), ("non(M)", "cof(M)"), ("cof(M)", "cof(N)"), ("cof(N)", r"$\mathfrak{c}$"),
         ("non(N)", "cof(N)"), ("add(M)", r"$\mathfrak{b}$"), (r"$\mathfrak{b}$", "non(M)"), ("cov(M)", r"$\mathfrak{d}$"), (r"$\mathfrak{d}$", "cof(M)"), (r"$\mathfrak{b}$", r"$\mathfrak{d}$")]
for a, b in edges:
    (x1, y1), (x2, y2) = nodes[a], nodes[b]
    if y1 == y2:   # horizontal: leave the box edges
        p1, p2 = (x1 + 0.5, y1), (x2 - 0.5, y2)
    else:
        p1, p2 = (x1, y1 + 0.24), (x2, y2 - 0.24)
    ax.annotate("", xy=p2, xytext=p1, arrowprops=dict(arrowstyle="-|>", color=INK, lw=1.1, shrinkA=0, shrinkB=0))
for k, (x, y) in nodes.items():
    c = RED if k in (r"$\mathfrak{b}$", r"$\mathfrak{d}$") else (BLUE if k in (r"$\aleph_1$", r"$\mathfrak{c}$") else INK)
    ax.add_patch(FancyBboxPatch((x - 0.48, y - 0.22), 0.96, 0.44, boxstyle="round,pad=0.02,rounding_size=0.12",
                                fc="white", ec=c, lw=1.4))
    ax.text(x, y, k, ha="center", va="center", fontsize=16 if k.startswith("$") or k.startswith("r") or "$" in k else 12, color=c)
ax.set_xlim(-1.3, 7.3); ax.set_ylim(-0.6, 3.0); ax.axis("off")
ax.set_title("Cichoń's diagram: ten in-between sizes, arrows point to larger (all provable in ZFC)", fontsize=12)
save(fig, "f5_cichon.png")
print("figures written")
