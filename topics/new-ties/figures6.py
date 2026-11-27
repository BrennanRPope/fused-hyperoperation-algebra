import json, numpy as np, matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from fractions import Fraction as F
import nocollide as N
all_words = N.all_words
plt.rcParams.update({"font.size": 11})
GREEN = "#2e9e5b"; BLUE = "#2f6db5"; ORANGE = "#e08a1e"; GREY = "#999999"; RED = "#c0392b"
# Fig 17: 2-adic levels where words ending in law 2 land, started at the repelling wall -Q
fig, axs = plt.subplots(1, 3, figsize=(13, 4.4), sharey=True)
for ax, Q, col in zip(axs, (2, 12, 3), (GREEN, GREEN, RED)):
    cnt = {}
    for n in (1, 2, 3, 4):
        for w in all_words(n):
            lv = -N.v(N.z1(N.word_apply(w, F(-Q), Q), Q), 2)       # level = log2 of the 2-adic size of z1
            cnt[lv] = cnt.get(lv, 0) + 1
    xs = sorted(cnt); tot = sum(cnt.values())
    ax.bar(xs, [cnt[x] / tot for x in xs], color=col, width=0.8)
    ax.set_xlim(-14, 12); ax.set_xlabel("2-adic level of z₁ (log₂ of its 2-adic size)")
    if len(xs) == 1:
        x0 = xs[0]
        ax.annotate("", xy=(x0 + 3, 0.45), xytext=(x0 + 0.6, 0.45), arrowprops=dict(arrowstyle="->", lw=2, color=BLUE))
        ax.text(x0 + 3.3, 0.40, "a ratio map C moves\nevery point right by ≥ 2\nlevels, off this bar", color=BLUE, fontsize=9)
        ax.set_title(f"walls ±{Q}: every word lands on ONE level\n(proved for all words, Theorem P5)", fontsize=10.5)
    else:
        ax.set_title(f"walls ±{Q}: words spread over {len(xs)} levels\n(so this 2-adic proof can't work)", fontsize=10.5)
axs[0].set_ylabel("share of words")
fig.suptitle("Figure 17. Why different words give different maps at walls ±2 and ±12. Every word that ends in the second law\n"
             f"(all {tot:,} with up to 4 runs, runs up to 2 letters / 4 c's) is evaluated at −Q and measured 2-adically.", fontsize=11.5)
plt.tight_layout(); plt.savefig("figures/fig17-two-adic-levels.png", dpi=130, bbox_inches="tight"); plt.close()
# Fig 18: what is now proved about the dimension, with no unchecked assumption of ours
fig, ax = plt.subplots(figsize=(10.5, 4.6))
rows = [("walls ±2", 0.3810, 0.3818, None, None), ("walls ±3", 0.4259, 0.4418, 0.4405, 0.4418), ("walls ±12", 0.5581, 0.5599, None, None)]
for i, (lab, lo, hi, clo, chi) in enumerate(rows):
    y = 2 - i
    ax.plot([lo, hi], [y, y], color=GREEN, lw=9, solid_capstyle="butt", label="proved now, resting only on the cited theorem" if i == 0 else None)
    if clo:
        ax.plot([clo, chi], [y - 0.22, y - 0.22], color=ORANGE, lw=5, solid_capstyle="butt",
                label="tighter range if no collisions at ±3 (checked to length 14, not proved)")
        ax.text(lo, y + 0.18, f"{lo:.3f}", ha="center", fontsize=9); ax.text(hi, y + 0.18, f"{hi:.4f}", ha="center", fontsize=9)
    else:
        ax.text((lo + hi) / 2, y + 0.18, f"{lo:.4f} to {hi:.4f}", ha="center", fontsize=9)
ax.axvline(0.6104, color=GREY, ls="--", lw=1.5); ax.text(0.6085, 2.45, "line-case formula D* = 0.610", ha="right", color=GREY, fontsize=9)
ax.set_yticks([2, 1, 0]); ax.set_yticklabels([r[0] for r in rows]); ax.set_ylim(-0.6, 2.7); ax.set_xlim(0.35, 0.63)
ax.set_xlabel("true dimension of the union fractal")
ax.set_title("Figure 18. The original examples' dimensions, now proved with no assumption of ours.\n"
             "±2 and ±12: collisions ruled out for every word. ±3: a free part of the system gives the lower end.", fontsize=11)
ax.legend(fontsize=9, loc="lower right"); plt.savefig("figures/fig18-unconditional-dimensions.png", dpi=130, bbox_inches="tight"); plt.close()
