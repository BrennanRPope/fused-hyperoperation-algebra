"""Figures for topics/hyperoperation-ladder. Run: python3 figures.py (writes figures/*.png)."""
import math, os, sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "fractional-rungs"))
OUT = os.path.join(HERE, "figures")
os.makedirs(OUT, exist_ok=True)

INK, INK2, GRID, AXIS = "#191919", "#666666", "#e7e7e7", "#c2c2c2"
BLUE, ORANGE = "#257de0", "#f47230"
GOOD, BAD = "#0aac23", "#dc3c3d"
SEQ = ["#9dc8f7", "#84b9f3", "#6aabf1", "#519ced", "#338ced", "#257de0", "#1b64be", "#0f4a96", "#0e3d7e"]
plt.rcParams.update({"font.size": 11, "text.color": INK, "axes.labelcolor": INK2, "xtick.color": INK2,
                     "ytick.color": INK2, "axes.edgecolor": AXIS, "axes.titlecolor": INK,
                     "axes.spines.top": False, "axes.spines.right": False, "figure.dpi": 150})

def save(fig, name):
    fig.savefig(os.path.join(OUT, name), bbox_inches="tight", facecolor="white")
    plt.close(fig)

# ---------------------------------------------------------------- F1: two ladders, one shared step
fig, ax = plt.subplots(figsize=(9, 5.6))
ax.set_xlim(0, 10); ax.set_ylim(-0.6, 6.4); ax.axis("off")
whole = [("rung 0", "successor  b + 1"), ("rung 1", "a + b + 1"), ("rung 2", "a ∘ b = a + b + ab"),
         ("rung 3", "(1+a)^(1+b) − 1"), ("rung 4", "fused tetration"), ("rung 5", "fused pentation")]
height = [(-1, "log(eᵃ + eᵇ − 1)"), (0, "a + b"), (0.5, "B½  (no closed form)"), (1, "a ∘ b"),
          (1.5, "B³ᐟ²"), (2, "(1+a)(1+b)e^(L(a)L(b)) − 1")]
for i, (r, f) in enumerate(whole):
    y = i
    ax.add_patch(FancyBboxPatch((0.3, y - 0.32), 3.9, 0.64, boxstyle="round,pad=0.02,rounding_size=0.08",
                                fc="#eaf3fd", ec=BLUE, lw=2 if i == 2 else 1))
    ax.text(0.5, y, r, va="center", fontsize=10, color=INK2)
    ax.text(1.45, y, f, va="center", fontsize=11)
for s, f in height:
    y = 1 + 2 * s * (1 if s >= 0 else 1)   # height 0 at y=1, height 1 at y=3, so the ∘ boxes line up
    y = {-1: -0.0, 0: 1, 0.5: 2, 1: 2, 1.5: 4, 2: 5}[s] if False else y
    ys = {-1: 0, 0: 1, 0.5: 1.5, 1: 2, 1.5: 3.5, 2: 5}[s]
    ax.add_patch(FancyBboxPatch((5.8, ys - 0.22), 3.9, 0.44, boxstyle="round,pad=0.02,rounding_size=0.08",
                                fc="#fdeee6", ec=ORANGE, lw=2 if s == 1 else 1))
    ax.text(6.0, ys, f"height {s:g}", va="center", fontsize=10, color=INK2)
    ax.text(7.35, ys, f, va="center", fontsize=10.5)
ax.annotate("", xy=(5.75, 2), xytext=(4.25, 2), arrowprops=dict(arrowstyle="<->", color=INK, lw=1.6))
ax.text(5.0, 2.22, "same operation", ha="center", fontsize=10)
ax.text(2.25, 6.0, "Whole-rung ladder\n(Goodstein's hyperoperations, moved by x ↦ x+1)", ha="center", fontsize=11, color=BLUE)
ax.text(7.75, 6.0, "Height ladder\n(conjugate + by iterates of eˣ − 1)", ha="center", fontsize=11, color=ORANGE)
ax.text(5.0, -0.55, "They share exactly one step. Above and below ∘ the two ladders are different operations.",
        ha="center", fontsize=10, color=INK2)
save(fig, "f1_two_ladders.png")

# ---------------------------------------------------------------- F2: the law board
laws = ["commutative", "associative", "identity 0\n(two-sided)", "0 right\nidentity", "0 left\nabsorber",
        "left\ndistributive\nover the\nrung below", "right\ndistributive\nover the\nrung below"]
rows = ["rung 0  b+1", "rung 1  a+b+1", "rung 2  ∘", "rung 3", "rung 4", "rung 5", "every rung ≥ 6"]
H, F, NA = 1, 0, -1
table = [
    [F, F, F, F, F, NA, NA],
    [H, H, F, F, F, H, H],
    [H, H, H, H, F, H, H],
    [F, F, F, H, H, F, H],
    [F, F, F, H, H, F, F],
    [F, F, F, H, H, F, F],
    [F, F, F, H, H, F, F],
]
fig, ax = plt.subplots(figsize=(10.5, 5.2))
for i, row in enumerate(table):
    for j, v in enumerate(row):
        c = {H: GOOD, F: BAD, NA: "#efefef"}[v]
        ax.add_patch(plt.Rectangle((j + 0.04, -i - 0.92), 0.92, 0.88, fc=c, ec="white", lw=2))
        ax.text(j + 0.5, -i - 0.48, {H: "holds", F: "fails", NA: "—"}[v], ha="center", va="center",
                color="white" if v != NA else INK2, fontsize=10, weight="bold")
for j, l in enumerate(laws):
    ax.text(j + 0.5, 0.15, l, ha="center", va="bottom", fontsize=9)
for i, r in enumerate(rows):
    ax.text(-0.1, -i - 0.48, r, ha="right", va="center", fontsize=10.5)
ax.set_xlim(-2.2, len(laws)); ax.set_ylim(-len(rows) - 0.3, 1.3); ax.axis("off")
ax.set_title("Every cell is a Lean theorem (HyperLadder.lean), for every rung n, standard axioms only",
             fontsize=11, loc="left", x=-0.17)
ax.text(-2.2, -len(rows) - 0.25, "The 'every rung ≥ 6' row is one theorem per column, not a sample. Distributivity is "
        "rung n+1 over rung n; rung 0 has no rung below.", fontsize=9, color=INK2)
save(fig, "f2_law_board.png")

# ---------------------------------------------------------------- F3: symmetry maps
def lad(n, a, b):
    if n == 0: return b + 1
    if n == 1: return a + b + 1
    if n == 2: return a + b + a * b
    return (a + 1) ** (b + 1) - 1
N = 7
fig, axs = plt.subplots(1, 4, figsize=(12, 3.4))
for n, ax in enumerate(axs):
    M = np.array([[np.sign(lad(n, a, b) - lad(n, b, a)) for a in range(N)] for b in range(N)])
    cols = {-1: "#f6c1aa", 0: "#efefef", 1: "#9dc8f7"}
    for a in range(N):
        for b in range(N):
            ax.add_patch(plt.Rectangle((a, b), 0.94, 0.94, fc=cols[M[b, a]], ec="none"))
    ax.set_xlim(0, N); ax.set_ylim(0, N); ax.set_aspect("equal")
    ax.set_xticks(np.arange(N) + 0.47, range(N)); ax.set_yticks(np.arange(N) + 0.47, range(N))
    ax.tick_params(length=0); [s.set_visible(False) for s in ax.spines.values()]
    ax.set_xlabel("a"); ax.set_ylabel("b" if n == 0 else "")
    ax.set_title(f"rung {n}: " + ("symmetric" if n in (1, 2) else "not symmetric"), fontsize=11,
                 color=GOOD if n in (1, 2) else INK)
fig.text(0.5, -0.04, "Blue: a⋆b > b⋆a.   Orange: a⋆b < b⋆a.   Gray: equal.   Rungs 4 and up repeat rung 3's pattern "
         "(proved: commutative only at rungs 1 and 2).", ha="center", fontsize=9.5, color=INK2)
save(fig, "f3_symmetry.png")

# ---------------------------------------------------------------- F4: the dial collapses at k = 0
fig, axs = plt.subplots(1, 3, figsize=(12, 3.6))
k = np.logspace(-3, 0, 200)
axs[0].loglog(k, 2 + 1 / k, color=BLUE, lw=2)
axs[0].set_xlabel("dial k"); axs[0].set_title("rung 1 at a = b = 1: 2 + 1/k\nblows up as k → 0", fontsize=10.5)
b = np.linspace(0, 3, 100)
for kk, c in zip([1, 0.5, 0.2, 0], [SEQ[7], SEQ[5], SEQ[3], INK2]):
    axs[1].plot(b, 1 + b + kk * b, color=c, lw=2, ls="--" if kk == 0 else "-")
    axs[1].text(3.05, 1 + 3 + kk * 3, f"k = {kk:g}" + (": +" if kk == 0 else ""), va="center", fontsize=9)
axs[1].set_xlim(0, 3.6); axs[1].set_xlabel("b"); axs[1].set_title("rung 2: 1 ∘ₖ b = 1 + b + kb\nflattens to + at k = 0", fontsize=10.5)
m = np.arange(1, 7)
for kk, c in zip([1, 0.5, 0.2], [SEQ[7], SEQ[5], SEQ[3]]):
    axs[2].plot(m, ((1 + kk) ** m - 1) / kk, "o-", color=c, lw=2, ms=6, label=f"k = {kk:g}")
axs[2].plot(m, m, "o--", color=INK2, lw=2, ms=6, label="k → 0: m·a")
axs[2].set_xlabel("m = 1 + kb"); axs[2].set_title("rung 3 at a = 1: ((1+k)^m − 1)/k\nflattens to m·a at k = 0", fontsize=10.5)
axs[2].legend(frameon=False, fontsize=9)
for ax in axs: ax.grid(color=GRID, lw=0.6)
fig.text(0.5, -0.1, "Every k ≠ 0 is the k = 1 ladder in rescaled coordinates (proved). k = 0 is a different, degenerate object, "
         "not the ordinary + × ^ ladder.", ha="center", fontsize=9.5, color=INK2)
save(fig, "f4_dial.png")

# ---------------------------------------------------------------- F5: heights fill in between, but never hit rung 3
import regular_iteration as RI
fig, ax = plt.subplots(figsize=(8.5, 5))
A = np.linspace(0.0, 2.0, 41)
hs = [-1, -0.5, 0, 0.5, 1, 1.5, 2]
for s, c in zip(hs, [SEQ[0], SEQ[1], SEQ[2], SEQ[4], SEQ[5], SEQ[6], SEQ[8]]):
    ys = [RI.B(s, a, 1.0) for a in A]
    ax.plot(A, ys, color=c, lw=2)
    ax.text(2.03, ys[-1] + {-1: -0.12, 0: 0.18, -0.5: 0.0}.get(s, 0), f"height {s:g}", va="center", fontsize=9)
ax.plot(A, (1 + A) ** 2 - 1, color=ORANGE, lw=2, ls="--")
ax.text(2.03, 8.0, "rung 3: (1+a)² − 1", va="center", color=INK, fontsize=9)
ax.set_xlabel("a"); ax.set_ylabel("a ⋆ 1"); ax.grid(color=GRID, lw=0.6)
ax.set_title("Height ladder B_s(a, 1): one smooth family, commutative at every height.\n"
             "Height 2 is not rung 3; the families meet only at height 1 = rung 2.", fontsize=10.5, loc="left")
ax.set_xlim(0, 2.6)
save(fig, "f5_heights.png")
print("wrote", sorted(os.listdir(OUT)))
