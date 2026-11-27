"""Round 14 pictures (DRAFT): fig44 the exhaustive psi-collision search, fig45 the neighbour-map plane,
fig46 excursions for all stretches at once, and the 19-adic gap a bad pair would have to cross."""
import json
import numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from fractions import Fraction as F
PS = json.load(open("psi14.json")); NB = json.load(open("nb14.json")); U = json.load(open("unif14.json"))

# fig44
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 4.8))
n = PS["n"]; lv = np.arange(1, n + 1); sizes = np.array(PS["level_sizes"])
first = {q: min(int(k) for k in r["by_max_length"]) for q, r in PS["ratios"].items()}
classes = [sum(1 for q in first if first[q] <= L) for L in lv]
ax1.semilogy(lv, np.cumsum(sizes), "o-", color="#1f77b4", label="words searched, up to this length")
ax1b = ax1.twinx(); ax1b.step(lv, classes, where="mid", color="#d62728", lw=2, label="different collision ratios found")
ax1b.set_ylim(0, 8); ax1b.set_ylabel("different ratios s_h/s_g", color="#d62728")
ax1.set_xlabel("letters per word"); ax1.set_ylabel("words (log scale)", color="#1f77b4")
ax1.set_title("A. 29 million words, still only 5 ratios")
h1, l1 = ax1.get_legend_handles_labels(); h2, l2 = ax1b.get_legend_handles_labels(); ax1.legend(h1 + h2, l1 + l2, fontsize=8, loc="upper left")
items = sorted(PS["ratios"].items(), key=lambda t: t[1]["v19"])
xs = np.arange(len(items)); cnts = [r["count"] for q, r in items]
ax2.bar(xs, cnts, color=["#2ca02c" if r["v19"] >= -2 else "#999" for q, r in items])
for x, (q, r) in zip(xs, items):
    ax2.text(x, r["count"] + 0.4, "19-adic val. %d\ne.g. %s ~ %s" % (r["v19"], r["example"][0], r["example"][1]), ha="center", fontsize=7.5)
ax2.set_xticks(xs); ax2.set_xticklabels(["s_h/s_g = " + q for q, r in items], fontsize=8)
ax2.set_ylim(0, max(cnts) * 1.35); ax2.set_ylabel("collision pairs found, up to 14 letters")
ax2.set_title("B. All are the known front relations (valuations -2 to 1);\na bad pair would need valuation a >= 13")
fig.suptitle("Fig 44 (DRAFT). Every pair of reduced words of up to 14 letters with equal psi, found exhaustively", fontsize=11)
fig.tight_layout(); fig.savefig("figures/fig44-psi-collisions.png", dpi=130); plt.close(fig)

# fig45
S = NB["states_ab"]; a = np.array([s["a"] for s in S]); b = np.array([s["b"] for s in S]); dep = np.array([s["depth"] for s in S])
lin = b == 0
fig, ax = plt.subplots(figsize=(10, 5.2))
ax.axvspan(9 / 19, 19 / 9, color="#ffd27f", alpha=0.35, label="windows where leaving the line is possible (9/19 to 19/9)")
ax.scatter(a[lin], b[lin], s=6, color="#1f77b4", label="on the line (b = 0): the stretches lam^i / mu^j  (%d)" % lin.sum())
sc = ax.scatter(a[~lin], b[~lin], s=6, c=dep[~lin], cmap="viridis", label="off the line (%d), coloured by length so far" % (~lin).sum())
ax.set_xscale("log"); ax.set_xlabel("a, the size ratio of the two prefixes (log scale)"); ax.set_ylabel("b, the offset")
ax.set_title("Fig 45 (DRAFT). Reading two words from the left, every pair up to %d letters in all (exact)" % NB["depth"])
ax.text(0.02, 0.05, "No path returns to the identity (a = 1, b = 0) except by equal letters.\n"
        "From every stretch on the line, every trip off it ends within %d letters and never comes back." % NB["per_stretch_longest_trip"],
        transform=ax.transAxes, fontsize=8.5, bbox=dict(facecolor="white", alpha=0.8, edgecolor="#999"))
plt.colorbar(sc, ax=ax, label="combined length")
ax.legend(fontsize=8, loc="center left", bbox_to_anchor=(0.0, 0.42)); fig.tight_layout(); fig.savefig("figures/fig45-neighbour-plane.png", dpi=130); plt.close(fig)

# fig46
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 4.8), gridspec_kw=dict(width_ratios=[1.25, 1]))
L = sorted(int(k) for k in U["nodes_by_length"]); nodes = [U["nodes_by_length"][str(k)] for k in L]
ra = [U["returns_by_length"].get(str(k), 0) for k in L]; ru = [U["unit_returns_by_length"].get(str(k), 0) for k in L]
fu = {q: r["first_length"] for q, r in U["unit_alphas"].items()}
du = [sum(1 for q in fu if fu[q] <= k) for k in L]
ax1.semilogy(L, nodes, "o-", color="#1f77b4", label="branches alive (real sizes only)")
ax1.semilogy(L, [max(x, 0.8) for x in ra], "s-", color="#ff7f0e", label="returns to the line at a real stretch")
ax1.semilogy(L, [max(x, 0.8) for x in ru], "^-", color="#2ca02c", label="returns whose stretch has only 2, 5, 19 in it")
ax1.semilogy(L, du, "D-", color="#d62728", lw=2, label="different stretches among those (flat at 9 from length 4)")
ax1.set_xlabel("letters in the excursion (both words)"); ax1.set_ylabel("count (log scale)")
ax1.set_title("A. Real sizes alone: ~3.5x more branches per letter.\nPrimes: only 9 stretches ever come back")
ax1.legend(fontsize=7.5, loc="upper left")
vals = sorted(U["unit_alphas"], key=lambda q: float(F(q)))
v19 = [U["unit_alphas"][q]["v19"] for q in vals]
ax2.axvspan(-30, -12.5, color="#d62728", alpha=0.12); ax2.axvspan(12.5, 30, color="#d62728", alpha=0.12)
ax2.text(-21, 1.5, "a bad pair's first trip\nstarts at stretch 1/rho,\nvaluation -a <= -13", ha="center", fontsize=8, color="#d62728")
ax2.text(21, 1.5, "(mirror image)", ha="center", fontsize=8, color="#d62728")
from collections import Counter
cnt = Counter(v19); ax2.bar(list(cnt), list(cnt.values()), color="#2ca02c", width=0.8)
for k in cnt: ax2.text(k, cnt[k] + 0.05, ",".join(q for q in vals if U["unit_alphas"][q]["v19"] == k), ha="center", fontsize=6.5, rotation=90, va="bottom")
ax2.set_xlim(-30, 30); ax2.set_ylim(0, 4.2); ax2.set_xlabel("19-adic valuation of the returning stretch"); ax2.set_ylabel("how many of the 9")
ax2.set_title("B. The gap a second collision would have to cross")
fig.suptitle("Fig 46 (DRAFT). Trips off the line for every stretch at once, up to %d letters" % U["maxlen"], fontsize=11)
fig.tight_layout(); fig.savefig("figures/fig46-all-stretches.png", dpi=130); plt.close(fig)
print("ok")
