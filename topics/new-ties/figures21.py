"""Round 21 pictures (DRAFT): fig68 where Brennan's 1, 2, 4, 8 comes from (the letters' powers of 2);
fig69 vanishing sums of 2-5-19 numbers: thousands in general, only 19 + 1 = 20 with word-shaped sides;
fig70 the all-lengths reduction (Theorem 21.1) and where it stops.  Usage: python3 figures21.py (needs explore21.json)."""
import json, math
from fractions import Fraction as F
import numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch
import explore21 as E

J = json.load(open("explore21.json"))

# ---------------------------------------------------------------- fig68
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(13, 4.9))
letters = ["N  (1/19)", "K  (1/10)", "m  (1/20)", "l  (19/200)"]
b1 = [0, 1, 2, 3]; b3 = [0, 3, 4, 7]
x = np.arange(4)
r1 = ax1.bar(x - 0.2, b1, 0.38, color="#1f77b4", label="the break of 13.4 (kappa = 1/10, mu = 1/20)")
r3 = ax1.bar(x + 0.2, b3, 0.38, color="#bbbbbb", label="the 3-adic break (kappa = 3/40, mu = 3/80; section 14)")
for i in range(4):
    ax1.text(i - 0.2, b1[i] + 0.15, "2^%d = %d" % (b1[i], 2 ** b1[i]), ha="center", fontsize=9, color="#1f77b4", weight="bold")
    ax1.text(i + 0.2, b3[i] + 0.15, "%d" % 2 ** b3[i], ha="center", fontsize=8.5, color="#555555")
ax1.set_xticks(x); ax1.set_xticklabels(["N", "K", "m", "l"], fontsize=11)
ax1.set_ylabel("power of 2 in the letter's stretch 1/s")
ax1.set_ylim(0, 8.5); ax1.legend(fontsize=8, loc="upper left")
ax1.set_title("A. Brennan's 1, 2, 4, 8 is really there: the stretches 19, 10, 20, 200\n"
              "have 2-parts 1, 2, 4, 8. The other known break gives 1, 8, 16, 128", fontsize=9.5)
A_, B_ = np.meshgrid(range(1, 7), range(0, 7))
ax2.set_xlim(0.5, 6.5); ax2.set_ylim(-0.5, 6.5)
for a in range(1, 7):
    for b in range(0, 7):
        hit = (a, b) == (1, 2)
        ax2.add_patch(plt.Rectangle((a - 0.45, b - 0.45), 0.9, 0.9, color="#1f77b4" if hit else "#f0f0f0"))
        ax2.text(a, b, "%d,%d,%d,%d" % (1, 2 ** a, 2 ** b, 2 ** (a + b)), ha="center", va="center", fontsize=6.5,
                 color="white" if hit else "#555555")
ax2.set_xlabel("a = power of 2 in 1/kappa"); ax2.set_ylabel("b = power of 2 in 1/mu")
ax2.set_title("B. On the collision surface the four 2-parts are always 1, 2^a, 2^b, 2^(a+b)\n"
              "(proved: nu = kappa/(2 - kappa), lambda = mu (2 - kappa)). Doubling needs a = 1, b = 2", fontsize=9.5)
fig.suptitle("Fig 68 (DRAFT). The 1, 2, 4, 8 comes from the example's chosen numbers, belongs to the four letters, and never to a length",
             fontsize=11)
fig.tight_layout(rect=(0, 0, 1, 0.94)); fig.savefig("figures/fig68-one-two-four-eight.png", dpi=130); plt.close(fig)

# ---------------------------------------------------------------- fig69
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(13, 4.9), gridspec_kw=dict(width_ratios=[1, 1.35]))
vals = [J["B"]["primitive"], 1]
ax1.bar([0, 1], vals, color=["#bbbbbb", "#2ca02c"], width=0.6)
ax1.set_yscale("log"); ax1.set_ylim(0.5, 1e5)
ax1.set_xticks([0, 1]); ax1.set_xticklabels(["each side 10-lacunary\n(up to 10^9, 3 terms a side)",
                                            "each side word-shaped\n(up to 10^15, 4 terms a side)"], fontsize=9)
ax1.text(0, vals[0] * 1.4, "{:,}".format(vals[0]), ha="center", fontsize=11)
ax1.text(1, 1.4, "1:  19 + 1 = 20", ha="center", fontsize=11, color="#2ca02c", weight="bold")
ax1.set_ylabel("primitive vanishing sums (up to scaling)")
ax1.set_title("A. Sums of products of 2, 5 and 19 that balance:\nplenty in general, exactly one with word-shaped sides", fontsize=9.5)
ex = [((19, 1), (20,), True), ((76, 5), (80, 1), False), ((361, 20), (380, 1), False), ((380, 19, 1), (400,), False)]
for k, (A, B, ok) in enumerate(ex):
    y = 3 - k
    for t in A:
        ax2.plot(math.log10(t), y + 0.12, "o", color="#1f77b4", ms=8)
        ax2.text(math.log10(t), y + 0.27, str(t), ha="center", fontsize=8, color="#1f77b4")
    for t in B:
        ax2.plot(math.log10(t), y - 0.12, "s", color="#ff7f0e", ms=8)
        ax2.text(math.log10(t), y - 0.42, str(t), ha="center", fontsize=8, color="#ff7f0e")
    lab = " + ".join(map(str, A)) + " = " + " + ".join(map(str, B))
    why = "both sides word-shaped" if ok else ("left side has more 2s further down" if k == 2 else
                                               "a neighbour ratio is no product of word ratios")
    ax2.text(3.0, y, lab + "\n" + why, fontsize=8.5, va="center", color="#2ca02c" if ok else "#555555")
ax2.set_xlim(-0.3, 5.4); ax2.set_ylim(-0.7, 3.7); ax2.set_yticks([])
ax2.set_xlabel("log10 of the term (blue: one side, orange: the other)")
ax2.set_title("B. A word's terms step down by word ratios 19^L / (r 200^L 20^M), r = 19 or 10:\n"
              "powers of 2 and 5 can only grow toward the large end. Most balances break that", fontsize=9.5)
fig.suptitle("Fig 69 (DRAFT). The only balance a word can use is 1 + 19 = 20, the identity behind N N = N m", fontsize=11)
fig.tight_layout(rect=(0, 0, 1, 0.94)); fig.savefig("figures/fig69-only-nineteen-plus-one.png", dpi=130); plt.close(fig)

# ---------------------------------------------------------------- fig70
fig, ax = plt.subplots(figsize=(13, 5.2)); ax.axis("off"); ax.set_xlim(0, 13); ax.set_ylim(0, 5.2)
def box(x, y, w, h, txt, col, fs=8.6):
    ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.08", fc=col, ec="#444444", lw=1))
    ax.text(x + w / 2, y + h / 2, txt, ha="center", va="center", fontsize=fs, wrap=True)
def arrow(x1, y1, x2, y2, txt=""):
    ax.annotate("", xy=(x2, y2), xytext=(x1, y1), arrowprops=dict(arrowstyle="->", lw=1.4, color="#444444"))
    if txt: ax.text((x1 + x2) / 2, (y1 + y2) / 2 + 0.15, txt, ha="center", fontsize=8, color="#444444")
G, O, R = "#d9f2d9", "#ffe9cc", "#f8d7d7"
box(0.2, 3.3, 2.6, 1.5, "Two canonical words\nwith the same fingerprint:\nS' psi(W) = S psi(W')\n(a balance of 2-5-19 terms)", G)
box(3.5, 3.3, 2.9, 1.5, "Take the smallest balancing\npiece that holds the largest\nterm. Its sides are\nword-shaped (Q5)", G)
box(7.1, 3.3, 2.8, 1.5, "Conjecture 21.2: that piece is\none term = one term, or\nc 19 + c = c 20.\nSearched: only these", O)
box(10.5, 3.3, 2.3, 1.5, "c 19 + c = c 20 needs a\nword starting N N:\nnot canonical. So the\nleading terms agree\n(Lemma 21.3)", G)
box(7.1, 0.5, 2.8, 1.6, "Cut the first letter: the\nrest, u and u', have the same\nfingerprint and one N/K\nletter fewer: induction", G)
box(3.5, 0.5, 2.9, 1.6, "Same canonical form of u, u'\nand the same t_u force\nu = u' (Lemma 21.4;\nonly lN = mK, Q6)", G)
box(0.2, 0.5, 2.6, 1.6, "So the two words are\nequal: Conjecture 21.2\nimplies Conjecture 15.3\nat every length (21.5)", G)
arrow(2.85, 4.05, 3.45, 4.05); arrow(6.45, 4.05, 7.05, 4.05); arrow(9.95, 4.05, 10.45, 4.05)
arrow(11.65, 3.25, 9.95, 1.9); arrow(7.05, 1.3, 6.45, 1.3); arrow(3.45, 1.3, 2.85, 1.3)
ax.text(8.5, 2.7, "the one open step (orange)", ha="center", fontsize=8.5, color="#b35900")
ax.text(6.5, 0.05, "Green: proved by hand in round 21 (with the earlier Proposition 15.2 and the canonical form of 19.3). "
        "Orange: open for all sizes; checked to 4 terms a side and 10^15.", ha="center", fontsize=8.5)
ax.set_title("Fig 70 (DRAFT). One 2-5-19 statement, Conjecture 21.2, would give the fingerprint conjecture for every length",
             fontsize=11)
fig.savefig("figures/fig70-all-lengths-reduction.png", dpi=130, bbox_inches="tight"); plt.close(fig)
print("ok")
