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

# ---------- Fig 3: collision map
items = [("⊞", ["mult-add (Base)", "rung-1 shifted add (Rungs)", "MARing product"], ["mathlib: biproduct"]),
         ("⊕", ["mult-add (README, Rungs)", "the (k, q) plane law"], ["Lean: sum of types", "direct sum, XOR"]),
         ("⊛", ["rung 3 (Rungs)", "mult-add (FunctionalAnalysis)"], ["mathlib: Day convolution", "convolution"]),
         ("✳", ["mult-add (most pages)"], ["none found"]),
         ("q", ["k + 1", "plane's second dial (−1/c²)"], ["Tsallis q = 1 − k", "q-analog base"]),
         ("k", ["dial on weights (+1 once)", "dial on probabilities (−1 once)"], ["Boltzmann k, a field"]),
         ("d", ["proposed plane dial", "zoom-dial dimension"], ["Coleman–McGuinness d (a match)"]),
         ("mult-add", ["a + b + ab"], ["computing: a·b + c"]),
         ("Fused", ["conjugated by x ↦ 1 + x"], ["fused multiply-add (FMA)"])]
fig, ax = plt.subplots(figsize=(13.5, 9.6)); ax.axis("off"); ax.set_xlim(0, 13.5); ax.set_ylim(-0.6, len(items) + 0.6)
ax.text(2.6, len(items) + 0.25, "meanings in this project", ha="center", weight="bold", fontsize=11)
ax.text(10.9, len(items) + 0.25, "meanings elsewhere", ha="center", weight="bold", fontsize=11)
for i, (name, ours, theirs) in enumerate(items):
    y = len(items) - 1 - i + 0.2
    inside = len(ours) > 1
    clash = theirs not in (["none found"], ["Coleman–McGuinness d (a match)"])
    fc = "#fadbd8" if (clash and inside) else ("#fdebd0" if clash or inside else "#d5f5e3")
    box(ax, 6.75, y, 1.6, 0.62, name, fc)
    for j, o in enumerate(ours):
        yy = y + (j - (len(ours) - 1)/2) * 0.3
        ax.text(5.6, yy, o, ha="right", va="center", fontsize=9.3, color="#1f618d")
        ax.plot([5.65, 5.95], [yy, y], color="#1f618d", lw=0.8)
    for j, t in enumerate(theirs):
        yy = y + (j - (len(theirs) - 1)/2) * 0.3
        ax.text(7.9, yy, t, ha="left", va="center", fontsize=9.3, color="#ba4a00" if clash else "#1e8449")
        ax.plot([7.55, 7.85], [y, yy], color="#ba4a00" if clash else "#1e8449", lw=0.8)
ax.text(6.75, -0.45, "Red: clashes both inside the project and outside.  Orange: one kind of clash.  Green: no clash found.",
        ha="center", fontsize=9.5, style="italic")
fig.suptitle("Figure 3. Naming collisions: each name, what it means here (left) and elsewhere (right), all present since v0.1.0 except d",
             fontsize=12.5, weight="bold", y=0.98)
fig.savefig(os.path.join(HERE, "fig3-collision-map.png"), dpi=130, bbox_inches="tight"); plt.close(fig)

# ---------- Fig 4: the circle operation on the dial, and the Lean glyph
fig, axes = plt.subplots(1, 2, figsize=(14, 5.6), gridspec_kw={"width_ratios": [1.25, 1]})
ax = axes[0]; ax.set_xlim(-2.3, 2.3); ax.set_ylim(-1.6, 2.2); ax.axis("off")
ax.set_title("One family, a ∘ₖ b = a + b + k·ab: both ring-theory conventions are points on it", weight="bold", fontsize=11)
ax.annotate("", xy=(2.2, 0), xytext=(-2.2, 0), arrowprops=dict(arrowstyle="-|>", color="#34495e", lw=1.5))
ax.text(2.2, -0.25, "k", fontsize=12, ha="right")
pts = [(-1, "∘₋₁: a + b − ab", "minus convention\nColeman–Easdown;\nMcConnell–Stokes", "#ba4a00"),
       (0, "∘₀: a + b", "plain addition", "#7f8c8d"),
       (1, "∘₁: a + b + ab", "plus convention\nDivinsky 1955;\nmathlib PreQuasiregular", "#1f618d")]
for x, top, bot, c in pts:
    ax.plot([x], [0], "o", ms=13, color=c); ax.text(x, -0.42, f"k = {x:+d}" if x else "k = 0", ha="center", fontsize=10)
    ax.text(x, 0.45, top, ha="center", fontsize=11, weight="bold", color=c)
    ax.text(x, 1.0, bot, ha="center", va="bottom", fontsize=9, color=c)
ax.annotate("", xy=(0.92, -0.75), xytext=(-0.92, -0.75), arrowprops=dict(arrowstyle="<->", color="#555", lw=1.2, connectionstyle="arc3,rad=0.25"))
ax.text(0, -1.35, "mirror a ↦ −a turns one into the other:  a ∘₋ₖ b = −((−a) ∘ₖ (−b))   (circle_mirror, proved in Lean)",
        ha="center", fontsize=9.5, style="italic")
ax = axes[1]; ax.axis("off"); ax.set_xlim(0, 10); ax.set_ylim(0, 10)
ax.set_title("Which glyph for Lean? (pages write ∘ₖ)", weight="bold", fontsize=11)
rows = [("∘", "Lean's function composition", "#fadbd8"), ("⊚", "mathlib: composition (TypeVec)", "#fadbd8"),
        ("○", "mathlib: composition of relations", "#fadbd8"), ("⊙", "mathlib: Hadamard product", "#fadbd8"),
        ("✳", "not on mathlib's character list; no shortcut", "#fdebd0"),
        ("◦", "free and allowed, but looks like ∘ (see the row above)", "#fdebd0"),
        ("◯", "free in mathlib, allowed, typed \\bigcirc, clearly not ∘", "#d5f5e3")]
for i, (g, note, fc) in enumerate(rows):
    y = 9.0 - 1.22*i
    ax.add_patch(FancyBboxPatch((0, y - 0.5), 10, 1.0, boxstyle="round,pad=0.02", fc=fc, ec="none"))
    ax.text(0.7, y, g, ha="center", va="center", fontsize=17); ax.text(1.6, y, note, va="center", fontsize=10)
ax.text(5, 0.25, "Recommended: a ◯ b for k = 1 and a ◯[k] b in Lean", ha="center", fontsize=10, weight="bold", color="#1e8449")
fig.suptitle("Figure 4. The circle operation, generalized: ∘ₖ on the page, ◯[k] in Lean", fontsize=13, weight="bold", y=1.0)
fig.tight_layout(); fig.savefig(os.path.join(HERE, "fig4-circle-operation.png"), dpi=130, bbox_inches="tight"); plt.close(fig)

# ---------- Fig 5: is the circle operation a "path" construction?
fig, axes = plt.subplots(1, 2, figsize=(14.5, 6.2), gridspec_kw={"width_ratios": [1, 1.15]})
ax = axes[0]; ax.axis("off"); ax.set_xlim(0, 10); ax.set_ylim(0, 10)
ax.set_title("A chain of 4 steps, composed with ∘₁, keeps all 15 ordered sub-chains", weight="bold", fontsize=10.5)
steps = ["a", "b", "c", "d"]
for i, s_ in enumerate(steps):
    box(ax, 1.6 + 2.2*i, 9.0, 1.2, 0.7, s_, "#d6eaf8")
    if i < 3: arrow(ax, (2.25 + 2.2*i, 9.0), (2.95 + 2.2*i, 9.0), "", "#1f618d")
import itertools
lv = {1: [], 2: [], 3: [], 4: []}
for r in range(1, 5):
    for c in itertools.combinations(steps, r): lv[r].append("".join(c))
ys = {1: 7.5, 2: 6.0, 3: 4.5, 4: 3.0}
for r in range(1, 5):
    ax.text(0.05, ys[r], f"length {r}:", fontsize=9, va="center", color="#555")
    ax.text(1.6, ys[r], "   ".join(lv[r]), fontsize=10.5, va="center", family="DejaVu Sans Mono")
ax.text(5, 1.9, "a ∘₁ b ∘₁ c ∘₁ d = the sum of all of them, in order\n(chain4_landscape, proved in Lean, non-commutative ring)", ha="center", fontsize=9.5)
ax.text(5, 0.6, "With dial k, a sub-chain of length m is weighted k^(m−1): k = 0 keeps only the single steps.",
        ha="center", fontsize=9, style="italic")

ax = axes[1]; ax.axis("off"); ax.set_xlim(0, 10); ax.set_ylim(0, 10)
ax.set_title("Where it sits: inside the existing \"keep the whole path\" constructions", weight="bold", fontsize=10.5)
box(ax, 5, 8.6, 6.4, 1.1, "Path algebra / free category (quivers)\nand Chen–Lyons path signatures: the whole path, nothing lost", "#d5f5e3")
box(ax, 5, 5.6, 6.4, 1.1, "Our chain composite (1+a)(1+b)(1+c)(1+d) − 1\n= strict iterated sums (discrete signature); k = dilation", "#d6eaf8")
box(ax, 5, 2.6, 6.4, 1.1, "The same composite over numbers:\none number, order forgotten (1 ∘₁ 2 = 2 ∘₁ 1 = 5)", "#fdebd0")
arrow(ax, (5, 8.0), (5, 6.2), "a special case of", "#1e8449", off=(0.2, 0), ha="left")
arrow(ax, (5, 5.0), (5, 3.2), "evaluate in numbers\n(loses order)", "#ba4a00", off=(0.2, 0), ha="left")
box(ax, 1.3, 0.8, 2.4, 0.9, "Homotopy paths:\nunrelated", "#eeeeee")
box(ax, 8.7, 0.8, 2.6, 0.9, "Path integrals:\nchecked separately", "#eeeeee")
fig.suptitle("Figure 5. \"Path theory\": ours is a case of the path constructions, not the other way round",
             fontsize=13, weight="bold", y=1.0)
fig.tight_layout(); fig.savefig(os.path.join(HERE, "fig5-path-names.png"), dpi=130, bbox_inches="tight"); plt.close(fig)
