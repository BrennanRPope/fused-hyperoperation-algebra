"""Pictures for semiring-audit/README.md (DRAFT). Run: python3 figures.py"""
import itertools
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch

OUT = "figures/"
BLUE, GREEN, RED, ORANGE, GREY, PURPLE = "#2e6fac", "#2e8b57", "#c0392b", "#d98c1f", "#777777", "#7d4fa3"
plt.rcParams.update({"font.size": 11, "axes.spines.top": False, "axes.spines.right": False})
rng = np.random.default_rng(3)


def box(ax, x, y, w, h, text, color, fs=11, tc="white"):
    ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.02,rounding_size=0.04",
                                fc=color, ec="none"))
    ax.text(x + w/2, y + h/2, text, ha="center", va="center", color=tc, fontsize=fs)


def arrow(ax, p, q, text="", color=GREY, fs=10, off=(0, 0.04), style="-|>"):
    ax.add_patch(FancyArrowPatch(p, q, arrowstyle=style, mutation_scale=16, color=color, lw=1.8))
    if text:
        ax.text((p[0]+q[0])/2 + off[0], (p[1]+q[1])/2 + off[1], text, ha="center", va="bottom",
                fontsize=fs, color=color)


# ------------------------------------------------------------------ fig1: one semiring, three rulers
fig, ax = plt.subplots(figsize=(12, 5.6)); ax.set_xlim(0, 12); ax.set_ylim(0, 5.6); ax.axis("off")
ax.set_title("One semiring read on three rulers: the pair our Lean proves distributive again and again", fontsize=13)
cols = [(0.4, "Lift ruler  y", "or:   y₁ + y₂\nand:  y₁ × y₂\n\nnothing = 0\none = 1", BLUE,
         "ordinary arithmetic\n(counting, probability)"),
        (4.4, "Fused ruler  a = y − 1", "or:   a + b + 1   (rung 1, sadd)\nand:  a + b + ab   (∘₁, mult-add)\n\nnothing = −1  (the wall)\none = 0", PURPLE,
         "our FusedPeano, AbstractAlgebra,\nPvsNP, LinearAlgebra, Rungs"),
        (8.4, "Height ruler  a = E(s) = eˢ − 1", "or:   log(eˢ¹ + eˢ²)\nand:  s₁ + s₂\n\nnothing = −∞\none = 0", GREEN,
         "the log semiring; at its cold end\nthe tropical (min, +) semiring")]
for x, title, body, colr, foot in cols:
    box(ax, x, 3.95, 3.2, 0.7, title, colr, fs=12)
    box(ax, x, 1.55, 3.2, 2.25, body, "#f3f3f3", fs=11, tc="black")
    ax.text(x + 1.6, 1.15, foot, ha="center", va="top", fontsize=10, color=colr)
arrow(ax, (3.65, 2.7), (4.35, 2.7), "x ↦ x − 1", color=GREY)
arrow(ax, (7.65, 2.7), (8.35, 2.7), "a = eˢ − 1", color=GREY)
ax.text(6, 5.15, "The same distributive law,  and(x, or(y, z)) = or(and(x, y), and(x, z)),  holds on every ruler,"
        " because each map carries both operations.", ha="center", fontsize=10.5, color="black")
fig.savefig(OUT + "fig1-one-semiring-three-rulers.png", dpi=140, bbox_inches="tight"); plt.close(fig)

# ------------------------------------------------------------------ fig2: map of hits
rows = [
    ("FusedPeano.madd_sadd, instCommSemiring", "(ℕ, a+b+1, a+b+ab) ≅ positive integers (+, ×)", "transport of structure", "known"),
    ("AbstractAlgebra.multAdd_fadd_distrib, shift_ring_hom", "ring (R, a+b+1, ∘₁), zero −1, one 0", "textbook exercise; Equiv.commRing", "known"),
    ("PvsNP.multAdd_sAdd_distrib, evalS_lift", "same law; circuits over the shifted semiring", "transport of semiring circuits", "known"),
    ("LinearAlgebra.sadd, fl_sadd, fl_scale", "semimodule maps over the shifted ring", "affine maps about the corner (−1,−1)", "known"),
    ("q-c D1/D6: ⊕ₖ = a+b+1/k with ∘ₖ", "generated semiring, generator g(x) = 1 + kx", "Pap's g-semiring / pseudo-analysis", "known; dial map ours"),
    ("sadd read on the height ruler", "log semiring (log-sum-exp, +)", "Maslov; Litvinov 2005", "link ours"),
    ("operations-catalog rescale; exploration soft min cut", "tropical limit of the log semiring", "Maslov dequantization", "known"),
    ("PvsNP.minplus_distrib, minplus_assoc", "(min, +) on ℕ", "tropical semiring; Mathlib Tropical", "known"),
    ("solitons-rungs check 3: log τ = lse of corners", "tropical curve = creases of the max", "tropical KdV / KP solitons", "known"),
    ("PvsNP.witness_fold, np/sharpP/parityP_point", "one fold, three semiring sums: OR, +, XOR", "Boolean, ℕ, GF(2): NP, #P, ⊕P", "known; fold reading ours"),
    ("OrderDial.ksum_distrib_gap, ksum_or", "AND distributes over the k-sum only on 0/1 or k=0", "Kleisli of relations vs bags", "known"),
    ("LogicSets distributivity defect", "Boolean lattice and ring on 0, 1", "Stone 1936", "known; defect formula ours"),
    ("PowersetDials.conv_W, cell_merge", "union convolution of product weights", "zeta transform; subset convolution", "known"),
    ("MeasureProbability.lam_append (and GameTheory)", "Sugeno λ-measure = g-decomposable measure", "Sugeno 1974; Weber 1984; Pap", "known"),
    ("InformationTheory.outer_sum (K5)", "Σᵢⱼ xᵢyⱼ = (Σx)(Σy)", "generalized distributive law", "known"),
    ("Rungs.fexp_multAdd_right, ftet_no_multAdd_dist", "rung 3 right-distributes, rung 4 fails", "not a semiring (one-sided)", "known"),
]
fig, ax = plt.subplots(figsize=(16.5, 8.2)); ax.axis("off")
ax.set_title("Where our files already use semiring structure, and its standard name", fontsize=13, loc="left")
hdr = ["Our file and theorem", "What it really is", "Standard name", "Known or ours"]
xs = [0.0, 0.33, 0.62, 0.86]
for x, h in zip(xs, hdr):
    ax.text(x, 0.97, h, fontsize=11.5, weight="bold", transform=ax.transAxes)
for i, r in enumerate(rows):
    y = 0.92 - i*0.057
    if i % 2 == 0:
        ax.add_patch(plt.Rectangle((-0.005, y - 0.02), 1.0, 0.052, transform=ax.transAxes, fc="#f4f6f8", ec="none"))
    for j, (x, t) in enumerate(zip(xs, r)):
        colr = "black"
        if j == 3:
            colr = GREEN if t == "known" else ORANGE
        ax.text(x, y, t, fontsize=10, transform=ax.transAxes, color=colr, va="bottom")
fig.savefig(OUT + "fig2-map-of-hits.png", dpi=130, bbox_inches="tight"); plt.close(fig)

# ------------------------------------------------------------------ fig3: distributivity matrix
E, Ei = np.expm1, np.log1p
lift = lambda op: (lambda x, y: E(op(Ei(x), Ei(y))))
ops = {
    "log-sum-exp": lambda x, y: np.logaddexp(x, y),
    "+": lambda x, y: x + y,
    "×": lambda x, y: x * y,
    "a+b+1": lambda x, y: x + y + 1,
    "∘₁": lambda x, y: x + y + x*y,
}
ops["+⁽²⁾"] = lift(ops["∘₁"])
ops["E₊(a+b+1)"] = lift(ops["a+b+1"])
names = list(ops)
M = np.zeros((len(names), len(names)))
for i, mul in enumerate(names):
    for j, add in enumerate(names):
        A, B, C = rng.uniform(0.02, 0.4, size=(3, 4000))
        l = ops[mul](A, ops[add](B, C)); r = ops[add](ops[mul](A, B), ops[mul](A, C))
        m = np.isfinite(l) & np.isfinite(r)
        err = np.max(np.abs(l[m] - r[m]) / (1 + np.abs(r[m]))) if m.any() else 1
        M[i, j] = np.log10(max(err, 1e-17))
fig, ax = plt.subplots(figsize=(8.6, 7))
im = ax.imshow(M, cmap="RdYlGn_r", vmin=-16, vmax=0)
ax.set_xticks(range(len(names)), names, rotation=35, ha="right"); ax.set_yticks(range(len(names)), names)
ax.set_xlabel("the 'or' (summed across paths)"); ax.set_ylabel("the 'and' (applied along a path)")
for i in range(len(names)):
    for j in range(len(names)):
        ax.text(j, i, "✓" if M[i, j] < -9 else "", ha="center", va="center", fontsize=14, color="white")
ax.set_title("Does the row distribute over the column?\nGreen ✓ = yes (error at rounding level); red = no", fontsize=12)
cb = fig.colorbar(im, ax=ax, shrink=0.8); cb.set_label("log₁₀ of the largest distributivity error")
fig.savefig(OUT + "fig3-who-distributes-over-whom.png", dpi=140, bbox_inches="tight"); plt.close(fig)

# ------------------------------------------------------------------ fig4: two ladders
fig, ax = plt.subplots(figsize=(12, 5.2)); ax.set_xlim(0, 12); ax.set_ylim(0, 5.2); ax.axis("off")
ax.set_title("Bennett's ladder climbs by exp and every neighbouring pair is a semiring.\n"
             "Our height ladder climbs by E = eˣ − 1: neighbours are not, but each rung has a partner that is.", fontsize=12)
ben = ["log(eᵃ+eᵇ)", "a + b", "a × b", "a^(log b)"]
for i, t in enumerate(ben):
    box(ax, 0.6, 0.5 + i*1.05, 2.4, 0.7, t, BLUE)
    if i:
        arrow(ax, (1.8, 0.5 + (i-1)*1.05 + 0.72), (1.8, 0.5 + i*1.05 - 0.02), "", color=GREEN)
        ax.text(3.1, 0.5 + i*1.05 - 0.2, "distributes ✓", color=GREEN, fontsize=10)
ax.text(1.8, 4.85, "Bennett 1915 (exp)", ha="center", fontsize=11, color=BLUE)
ours = ["a + b   (height 0)", "a ∘₁ b   (height 1)", "a +⁽²⁾ b   (height 2)"]
part = ["log(eᵃ+eᵇ)", "a + b + 1", "E₊(a + b + 1)"]
for i, (t, p) in enumerate(zip(ours, part)):
    y = 0.9 + i*1.3
    box(ax, 5.6, y, 2.6, 0.75, t, PURPLE)
    box(ax, 9.2, y, 2.4, 0.75, p, GREEN)
    arrow(ax, (8.25, y + 0.37), (9.15, y + 0.37), "", color=GREEN)
    ax.text(8.7, y + 0.55, "✓", color=GREEN, fontsize=13, ha="center")
    if i:
        arrow(ax, (6.9, 0.9 + (i-1)*1.3 + 0.77), (6.9, y - 0.02), "", color=RED)
        ax.text(5.5, y - 0.33, "✗ does not distribute", color=RED, fontsize=10, ha="right")
ax.text(6.9, 4.85, "our addition at height s", ha="center", fontsize=11, color=PURPLE)
ax.text(10.4, 4.85, "its semiring partner", ha="center", fontsize=11, color=GREEN)
ax.text(10.4, 0.45, "partner at height s = (Eˢ)₊(log-sum-exp)", ha="center", fontsize=9.5, color=GREEN)
fig.savefig(OUT + "fig4-two-ladders.png", dpi=140, bbox_inches="tight"); plt.close(fig)

# ------------------------------------------------------------------ fig5: witness fold
def foldW(kk, xs):
    r = 0
    for x in reversed(xs): r = x + r + kk*x*r
    return r
fig, axs = plt.subplots(1, 2, figsize=(13, 4.6))
ns = np.arange(0, 7)
for kk, lab, colr, off in [(-1, "k = −1: OR  → Boolean semiring (NP: is there one?)", BLUE, -0.25),
                           (0, "k = 0: count → (ℕ, +, ×) (#P: how many?)", GREEN, 0),
                           (-2, "k = −2: XOR → GF(2) (⊕P: odd or even?)", ORANGE, 0.25)]:
    axs[0].bar(ns + off, [foldW(kk, [1]*n) for n in ns], width=0.25, color=colr, label=lab)
axs[0].set_xlabel("number of yes-witnesses folded"); axs[0].set_ylabel("value of the dial fold")
axs[0].set_title("One fold (PvsNP.witness_fold), three different semiring sums")
axs[0].legend(fontsize=9, loc="upper left"); axs[0].set_ylim(0, 8)
cs = np.linspace(-1, 2, 301); ks = np.linspace(-2.5, 1, 301)
Cg, Kg = np.meshgrid(cs, ks)
G = Kg * Cg * (Cg - 1)          # gap with a = b = 1
im = axs[1].imshow(np.sign(G)*np.log10(1 + np.abs(G)), origin="lower", aspect="auto",
                   extent=[cs[0], cs[-1], ks[0], ks[-1]], cmap="coolwarm", vmin=-0.8, vmax=0.8)
for cc in (0, 1): axs[1].axvline(cc, color="black", lw=1.5)
axs[1].axhline(0, color="black", lw=1.5)
for kk, t in [(-1, "OR"), (-2, "XOR")]:
    axs[1].axhline(kk, color="black", lw=0.7, ls=":"); axs[1].text(1.55, kk + 0.05, t, fontsize=10)
axs[1].set_xlabel("c (the value we AND with)"); axs[1].set_ylabel("dial k")
axs[1].set_title("Distributivity gap k·a·b·c·(c−1) at a = b = 1 (OrderDial)\nzero only on the black lines: c = 0, c = 1, or k = 0")
fig.colorbar(im, ax=axs[1], shrink=0.85, label="signed size of the gap")
fig.tight_layout(); fig.savefig(OUT + "fig5-witness-fold.png", dpi=140, bbox_inches="tight"); plt.close(fig)

# ------------------------------------------------------------------ fig6: Maslov / tropical dial
fig, axs = plt.subplots(1, 2, figsize=(13, 4.6))
x = np.linspace(-3, 3, 400)
for h, colr in [(1, PURPLE), (0.5, BLUE), (0.2, GREEN), (0.05, ORANGE)]:
    axs[0].plot(x, h*np.logaddexp(x/h, 0), color=colr, lw=2, label=f"temperature h = {h}")
axs[0].plot(x, np.maximum(x, 0), "k--", lw=1.5, label="max(a, 0): the tropical end")
axs[0].set_xlabel("a"); axs[0].set_ylabel("h · log(e^(a/h) + e^(0/h))")
axs[0].set_title("Cooling log-sum-exp gives max (Maslov dequantization)\n(operations-catalog 'rescale'; quantum-to-classical D6)")
axs[0].legend(fontsize=9)
cuts = rng.uniform(0, 3, size=(6,))
hs = np.logspace(-2.5, 0.7, 200)
soft = [-h*np.log(np.exp(-cuts/h).sum()) for h in hs]
axs[1].semilogx(hs, soft, color=BLUE, lw=2, label="softened min cut  −h·log Σ e^(−cut/h)")
axs[1].axhline(cuts.min(), color="k", ls="--", label="ordinary min cut (holographic RT surface)")
for cc in cuts: axs[1].axhline(cc, color=GREY, lw=0.6, alpha=0.6)
axs[1].set_xlabel("softness h (log scale)"); axs[1].set_ylabel("value")
axs[1].set_title("The exploration page's soft min cut is the same dial\n(grey lines: the six candidate cuts)")
axs[1].legend(fontsize=9, loc="lower left")
fig.tight_layout(); fig.savefig(OUT + "fig6-tropical-dial.png", dpi=140, bbox_inches="tight"); plt.close(fig)
print("figures written")
