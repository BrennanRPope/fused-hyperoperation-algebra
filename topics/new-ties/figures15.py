"""Round 15 pictures (DRAFT): fig47 every path off the line sits on one straight line set by two numbers psi(P), psi(Q);
fig48 the fingerprint searches (dense to 14 letters, sparse with long runs); fig49 what is proved and what is only searched."""
import json
import numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch
from fractions import Fraction as F
from kappa15 import kappa_int
K = json.load(open("kappa15.json")); S = json.load(open("sparse15.json"))
LET = {"l": (F(19, 200), F(0)), "m": (F(1, 20), F(0)), "N": (F(1, 19), F(1, 19)), "K": (F(1, 10), F(1, 10))}
def wmap(w):
    s, t = F(1), F(0)
    for x in w: a, b = LET[x]; s, t = s * a, s * b + t
    return s, t
def psi(w): s, t = wmap(w); return t / s
def unit(q): q = F(q); return kappa_int(q.numerator) == 1 and kappa_int(q.denominator) == 1

# fig47
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(13, 5.2), gridspec_kw=dict(width_ratios=[1.6, 1]))
X = np.geomspace(0.25, 4, 400)
ax1.plot(X, X - 1, color="#555", lw=2, label="every pair, rescaled: y = x - 1")
ax1.axhline(0, color="#1f77b4", lw=2.5, label="the line (b = 0): only the point x = 1")
pairs = [("NmN", "NN"), ("NK", "NlN"), ("NlN", "NmmN"), ("NN", "KN"), ("NmmN", "NlK"), ("NKN", "NNK"), ("NmN", "KmN")]
cols = plt.cm.tab10(np.arange(len(pairs)))
for i, ((P, Q), c) in enumerate(zip(pairs, cols)):
    pP, pQ = psi(P), psi(Q); a0 = wmap(P)[0] / wmap(Q)[0]; x0 = float(a0 * pP / pQ)
    lo, hi = float(pQ / (pP + F(1, 9)) * pP / pQ), float((pQ + F(1, 9)) / pP * pP / pQ)
    ax1.plot([x0], [x0 - 1], "o", color=c, ms=7, mec="k")
    ar = pQ / pP; ok = unit(ar)
    ax1.plot([], [], "o", color=c, label="%s vs %s: psi %s, %s; return stretch %s is %s" % (P, Q, pP, pQ, ar, "a stretch" if ok else "NOT a stretch"))
ax1.set_xscale("log"); ax1.set_xlim(0.25, 4); ax1.set_ylim(-0.8, 3)
from matplotlib.ticker import FixedLocator, NullLocator
ax1.xaxis.set_major_locator(FixedLocator([0.25, 0.5, 1, 2, 4])); ax1.xaxis.set_minor_locator(NullLocator()); ax1.set_xticklabels(["1/4", "1/2", "1", "2", "4"])
ax1.set_xlabel("x = a psi(P) / psi(Q)   (log scale; a = s_P / s_Q)"); ax1.set_ylabel("y = b / psi(Q)")
ax1.set_title("A. The state of a path off the line, P against Q, is b = psi(P) a - psi(Q):\nrescaled, every pair is one line (dot = the pair). The real cut keeps a pair alive\nonly within about 1/(9 psi(Q)) of x = 1, so real sizes cannot rule out a return")
ax1.legend(fontsize=6.6, loc="upper left")
ax2.axvspan(-2.5, 2.5, color="#2ca02c", alpha=0.35)
ax2.axvspan(12.5, 20, color="#d62728", alpha=0.3)
for x, y, lab in [(-2, 0.86, "200/361 (NN vs Km)"), (-1, 0.78, "10/19 (N vs K), 20/19 (NN vs Nm)"), (0, 0.70, "1/2 and 2 (KN vs Nm)"), (1, 0.62, "19/10, 19/20"), (2, 0.54, "361/200")]:
    ax2.plot([x], [y - 0.04], "v", color="k"); ax2.text(x + 0.4, y - 0.045, lab, ha="left", va="center", fontsize=7.5)
ax2.text(0, 0.2, "every size ratio the front\nrelations can make:\n19-adic valuation -2 to 2", ha="center", fontsize=9)
ax2.text(16, 0.2, "what a bad pair's second\ncollision needs:\nvaluation a >= 13", ha="center", fontsize=9)
ax2.set_xlim(-4, 20); ax2.set_ylim(0, 1); ax2.set_yticks([])
ax2.set_xlabel("19-adic valuation of the size ratio rho = s_h / s_g")
ax2.set_title("B. If fingerprints only match through the front relations\n(Conjecture 15.3), the gap can never be crossed (Cor. 15.4)")
fig.suptitle("Fig 47 (DRAFT). Lemma 15.1 (proved by hand): real sizes give no bounded invariant, only one line per pair; the question is arithmetic", fontsize=10.5)
fig.tight_layout(); fig.savefig("figures/fig47-one-line-per-pair.png", dpi=130); plt.close(fig)

# fig48
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(13, 4.9))
L = np.arange(1, K["n"] + 1)
ax1.semilogy(L, np.cumsum(K["words_per_length"]), "o-", color="#1f77b4", label="reduced words searched, up to this length")
ax1.semilogy(L, K["words_ending_in_N_or_K_by_length"], "s-", color="#9467bd", ms=4, label="of which end in N or K")
ax1.semilogy(L, K["distinct_canonical_forms_by_length"], "-", color="#2ca02c", lw=6, alpha=0.45, label="different canonical forms")
ax1.semilogy(L, K["distinct_fingerprints_by_length"], "k--", lw=1.4, label="different fingerprints (the same curve: no clash)")
ax1.set_xlabel("letters per word"); ax1.set_ylabel("count (log scale)")
ax1.set_title("A. Dense search, every word to 14 letters (exact integers):\n{:,} words, {:,} fingerprints, 0 clashes".format(K["words"], K["distinct_fingerprints_by_length"][-1]))
ax1.legend(fontsize=8, loc="upper left")
ax1.text(K["n"], 6, "psi is a power of 2, 5, 19 only for\nN, K, NN, KN (Thm 15.5 proves it\nfor <= 2 N/K letters)", ha="right", fontsize=8)
pl = S["plan"][1:]
xs = np.arange(len(pl)); ax2.bar(xs, [p["words"] for p in pl], color="#ff7f0e")
for x, p in zip(xs, pl):
    ax2.text(x, p["words"] * 1.15, "%s words\nl, m runs up to %d\n0 clashes" % (format(p["words"], ","), p["longest_run"]), ha="center", fontsize=8)
ax2.set_yscale("log"); ax2.set_ylim(1e4, 2e7)
ax2.set_xticks(xs); ax2.set_xticklabels(["%d N/K letters" % p["NK_letters"] for p in pl])
ax2.set_ylabel("words searched (log scale)")
ax2.set_title("B. Sparse search, the shape of a bad pair's words (few N/K,\nlong l and m runs, up to %d letters): 0 clashes" % (1 + 2 * 300 + 1))
fig.suptitle("Fig 48 (DRAFT). Conjecture 15.3 searched two ways: equal fingerprint only ever comes from the front relations", fontsize=10.5)
fig.tight_layout(); fig.savefig("figures/fig48-fingerprint-search.png", dpi=130); plt.close(fig)

# fig49
fig, ax = plt.subplots(figsize=(12.5, 5.6)); ax.axis("off"); ax.set_xlim(0, 12.5); ax.set_ylim(0, 6)
def box(x, y, w, h, text, kind):
    fc = {"proved": "#cfe8cf", "searched": "#ffe7b3", "goal": "#d6e4f5"}[kind]
    ec = {"proved": "#2ca02c", "searched": "#e08a00", "goal": "#1f77b4"}[kind]
    ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.08", fc=fc, ec=ec, lw=2, ls="-" if kind != "searched" else "--"))
    ax.text(x + w / 2, y + h / 2, text, ha="center", va="center", fontsize=8.6)
def arrow(x1, y1, x2, y2, t=""):
    ax.annotate("", xy=(x2, y2), xytext=(x1, y1), arrowprops=dict(arrowstyle="-|>", lw=1.6, color="#444"))
    if t: ax.text((x1 + x2) / 2 + (0.45 if x1 == x2 else 0), (y1 + y2) / 2 + (0 if x1 == x2 else 0.15), t, fontsize=7.5, ha="center", color="#444")
box(0.2, 4.2, 3.2, 1.4, "Lemma 15.1 (by hand)\nstate = line b = psi(P) a - psi(Q);\nreturn needs equal fingerprints", "proved")
box(0.2, 2.2, 3.2, 1.4, "Prop. 15.2 (by hand)\nfront relations and the\ncanonical form keep the fingerprint", "proved")
box(0.2, 0.2, 3.2, 1.4, "Thm 15.5 (by hand)\nwith <= 2 N/K letters, psi is a unit\nonly for N, K, NN, KN", "proved")
box(4.4, 2.2, 3.4, 1.6, "Conjecture 15.3 (searched)\nequal fingerprint => same canonical form\n29.2 million words to 14 letters, plus\n3.2 million sparse words: 0 clashes", "searched")
box(4.4, 4.4, 3.4, 1.2, "Conjecture 14.4 (Round 14)\nreturns only at 9 stretches", "searched")
box(8.8, 2.2, 3.4, 1.6, "Cor. 15.4 (by hand, from 15.3)\nno second collision for any bad pair:\nfront ratios have 19-adic size <= 2 < 13", "proved")
box(8.8, 4.4, 3.4, 1.2, "Goal: the original break\nis exactly s_F, for every a", "goal")
box(4.4, 0.2, 3.4, 1.4, "Why no proof yet: an equation in powers\nof 2, 5, 19 with unboundedly many terms\n(finiteness theorems need a fixed number)", "searched")
arrow(3.4, 4.9, 4.4, 3.4, "reduces to"); arrow(3.4, 2.9, 4.4, 2.9); arrow(3.4, 0.9, 4.4, 2.4, "first case")
arrow(6.1, 3.8, 6.1, 4.4, "implies"); arrow(7.8, 3.0, 8.8, 3.0, "implies"); arrow(10.5, 3.8, 10.5, 4.4, "gives")
ax.text(0.2, 5.85, "green, solid: proved   orange, dashed: only searched or open   blue: the goal", fontsize=8.5)
ax.set_title("Fig 49 (DRAFT). Round 15: what is proved and what is only searched", fontsize=11)
fig.tight_layout(); fig.savefig("figures/fig49-proved-vs-searched.png", dpi=130); plt.close(fig)
print("wrote fig47, fig48, fig49")
