"""Figures for the names and notation review (writes fig*.png next to this file)."""
import matplotlib, os
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch
HERE = os.path.dirname(os.path.abspath(__file__))
plt.rcParams.update({"font.size": 11, "axes.spines.top": False, "axes.spines.right": False})

def box(ax, x, y, w, h, text, fc):
    ax.add_patch(FancyBboxPatch((x - w/2, y - h/2), w, h, boxstyle="round,pad=0.02", fc=fc, ec="#34495e", lw=1.2))
    ax.text(x, y, text, ha="center", va="center", fontsize=10.5)
def arrow(ax, a, b, text, color, off=(0, 0), ha="center"):
    ax.add_patch(FancyArrowPatch(a, b, arrowstyle="-|>", mutation_scale=16, color=color, lw=1.8))
    ax.text((a[0]+b[0])/2 + off[0], (a[1]+b[1])/2 + off[1], text, ha=ha, va="center", fontsize=10, color=color, weight="bold")

# ---------- Fig 1: the sign flip
fig, axes = plt.subplots(1, 2, figsize=(13.5, 5.6), gridspec_kw={"width_ratios": [1.15, 1]})
ax = axes[0]; ax.axis("off"); ax.set_xlim(-1.2, 10); ax.set_ylim(0, 10)
ax.set_title("The same law, read on two sides, has opposite k", weight="bold")
W, P = "#d6eaf8", "#fdebd0"
box(ax, 2.0, 7.6, 3.0, 1.3, "weights x, y\nexample: 1 and 1", W)
box(ax, 8.0, 7.6, 3.0, 1.3, "x ✳ₖ y = x + y + k·xy\nexample k = 1: 3", W)
box(ax, 2.0, 2.4, 3.0, 1.3, "probabilities μ, ν\nexample: 1/2 and 1/2", P)
box(ax, 8.0, 2.4, 3.0, 1.3, "μ ✳₋ₖ ν = μ + ν − k·μν\nexample: 3/4", P)
arrow(ax, (3.6, 7.6), (6.4, 7.6), "compose with +k", "#1f618d", off=(0, 0.45))
arrow(ax, (3.6, 2.4), (6.4, 2.4), "compose with −k", "#ba4a00", off=(0, -0.45))
arrow(ax, (2.0, 6.85), (2.0, 3.15), "μ = x / (1 + kx)", "#555555", off=(-0.2, 0), ha="right")
arrow(ax, (8.0, 6.85), (8.0, 3.15), "same map", "#555555", off=(0.2, 0), ha="left")
ax.text(5.0, 5.0, "the square commutes\n(weight_prob_flip, proved in Lean)", ha="center", va="center", fontsize=10, style="italic")
ax.text(5.0, 0.55, "Check: one fair element, weight 1, has probability 1/2. Two joined: weights 1 ✳ 1 = 3,\n"
        "probability 3/4 = 1/2 + 1/2 − 1/4, the \"either one\" rule. Also 1 − kμ = 1/(1 + kx).", ha="center", fontsize=9)

ax = axes[1]; ax.axis("off"); ax.set_xlim(0, 10); ax.set_ylim(0, 10)
ax.set_title("Which side each page is on today", weight="bold")
rows = [("Counting / power set", "weights", "+1 = at most once", W),
        ("Repetition dial (item 1)", "weights", "+1/m = at most m", W),
        ("Zoom dial (item 2)", "weights", "k = 1 (pieces r^d)", W),
        ("Dial survey, item 4", "factor (1 − kz)^(−1/k)", "−1 subsets, +1 multisets", P),
        ("Primon row", "factor (1 − k p^−s)^(−1/k)", "−1 = fermions", P),
        ("Determinantal (DPP) row", "eigenvalues (probabilities)", "−1 = fermions", P)]
ax.text(0.1, 9.3, "page", weight="bold"); ax.text(4.0, 9.3, "quantity composed", weight="bold"); ax.text(7.3, 9.3, "what k means", weight="bold")
for i, (a, b, c, fc) in enumerate(rows):
    y = 8.3 - 1.25*i
    ax.add_patch(FancyBboxPatch((0, y - 0.45), 10, 0.9, boxstyle="round,pad=0.02", fc=fc, ec="none"))
    ax.text(0.1, y, a, va="center", fontsize=9.5); ax.text(4.0, y, b, va="center", fontsize=9.5); ax.text(7.3, y, c, va="center", fontsize=9.5)
ax.text(5.0, 0.35, "Blue rows read k on weights, orange rows on probabilities or the opposite-sign factor.\n"
        "Both are right. The freeze should say which side each page uses.", ha="center", fontsize=9, style="italic")
fig.suptitle("Figure 1. Why k looks like it has two signs: weights compose by ✳ₖ, probabilities by ✳₋ₖ",
             fontsize=13, weight="bold", y=1.0)
fig.tight_layout(); fig.savefig(os.path.join(HERE, "fig1-sign-flip.png"), dpi=130, bbox_inches="tight"); plt.close(fig)

# ---------- Fig 2: symbols and letters
fig, axes = plt.subplots(1, 2, figsize=(13.5, 5.4), gridspec_kw={"width_ratios": [1, 1.1]})
ax = axes[0]
syms = ["✳", "⊞", "⊕", "⊛"]; files = [256, 129, 143, 25]; uses = [1974, 1625, 1116, 308]
cols = ["#27ae60", "#c0392b", "#c0392b", "#7f8c8d"]
notes = {"✳": "keep: free in mathlib, already the Lean infix",
         "⊞": "clash: mathlib uses ⊞ for biproducts",
         "⊕": "clash: Lean's Sum type; also the plane law",
         "⊛": "local helper in two files"}
bars = ax.barh(syms[::-1], uses[::-1], color=cols[::-1])
for b, u, f, s_ in zip(bars, uses[::-1], files[::-1], syms[::-1]):
    ax.text(b.get_width() + 30, b.get_y() + b.get_height()/2, f"{u} uses in {f} files\n{notes[s_]}", va="center", fontsize=9)
ax.set_xlim(0, 3700); ax.set_xlabel("times the symbol appears in the public repo (all text files)")
ax.set_title("Symbols for mult-add in use now", weight="bold")

ax = axes[1]; ax.axis("off"); ax.set_xlim(0, 10); ax.set_ylim(0, 10)
ax.set_title("The same dial under other names", weight="bold")
tab = [("this project", "k", "a + b + k·ab", "#d5f5e3"),
       ("Coleman–McGuinness 1991", "c", "(x + y + cxy)/(1 − dxy)", "#eaf2f8"),
       ("Tsallis 1988", "1 − q", "q-logarithm, q = 1 − k", "#eaf2f8"),
       ("Sugeno 1974", "λ", "λ-measure: A + B + λAB", "#eaf2f8"),
       ("Box–Cox 1964", "λ", "(x^λ − 1)/λ transform", "#eaf2f8"),
       ("project's own q today", "k + 1", "clashes with Tsallis q", "#fadbd8"),
       ("plane's second dial today", "q", "rename to d (as in C–M)", "#fadbd8")]
ax.text(0.1, 9.3, "source", weight="bold"); ax.text(4.3, 9.3, "letter", weight="bold"); ax.text(5.9, 9.3, "where it appears", weight="bold")
for i, (a, b, c, fc) in enumerate(tab):
    y = 8.3 - 1.12*i
    ax.add_patch(FancyBboxPatch((0, y - 0.42), 10, 0.84, boxstyle="round,pad=0.02", fc=fc, ec="none"))
    ax.text(0.1, y, a, va="center", fontsize=9.5); ax.text(4.3, y, b, va="center", fontsize=10.5, weight="bold"); ax.text(5.9, y, c, va="center", fontsize=9.5)
ax.text(5.0, 0.2, "Green: keep. Red: two letters that mean something else elsewhere.", ha="center", fontsize=9, style="italic")
fig.suptitle("Figure 2. One symbol and one letter for the freeze: ✳ for mult-add, k for the dial, d for the plane",
             fontsize=13, weight="bold", y=1.0)
fig.tight_layout(); fig.savefig(os.path.join(HERE, "fig2-symbols-and-letters.png"), dpi=130, bbox_inches="tight"); plt.close(fig)
print("ok")
