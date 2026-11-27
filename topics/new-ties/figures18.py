"""Round 18 pictures (DRAFT): fig56 the independent recheck against the original bounds; fig57 the 19-adic size of T (proved
at most 2 for three letters, F2), and how it grows with more letters; fig58 the bounds for Theorem 18.1, round by round;
fig59 where Conjecture 15.3 stands now. Usage: python3 figures18.py (needs pairs18.json)."""
import json, itertools
import numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from pairs17 import v, runs, P_of, T3
from pairs18 import prefixes
P = json.load(open("pairs18.json")); R = P["rounds"]; E = P["exact_check"]

# fig56: my bounds against the independent referee's, stage by stage
fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(15, 5.2), gridspec_kw=dict(width_ratios=[1, 1, 1.3]))
names = ["16.1\n(Case II)", "16.2", "17.1"]
mine_start = [5.064e12, 1.532e17, 1.146e62]; ref_start = [4.939e12, 1.458e17, 1.1226e62]
x = np.arange(3); w = 0.36
ax1.bar(x - w / 2, mine_start, w, color="#1f77b4", label="this thread (Rounds 16, 17, log form corrected)")
ax1.bar(x + w / 2, ref_start, w, color="#ff7f0e", label="independent referee (Round 18)")
for xi, a, b in zip(x, mine_start, ref_start):
    ax1.text(xi - w / 2, a * 5, "%.3g" % a, ha="center", fontsize=7.5, rotation=90, va="bottom")
    ax1.text(xi + w / 2, b * 5, "%.3g" % b, ha="center", fontsize=7.5, rotation=90, va="bottom")
ax1.set_yscale("log"); ax1.set_ylim(1, 1e80); ax1.set_xticks(x); ax1.set_xticklabels(names, fontsize=8.5)
ax1.set_ylabel("Matveev's starting bound (log scale)"); ax1.legend(fontsize=7.5, loc="upper left")
ax1.set_title("A. Starting bounds from Matveev's theorem\n(they agree to within 5 percent)", fontsize=9.5)
lab = ["16.1 Case II:\nexponent E", "16.1 Case I:\ngap H", "16.2: long\nword's gap", "17.1: box\nfor exact check"]
mine_end = [10, 7, 39, 25]; ref_end = [10, 0, 9, 28]
x = np.arange(4)
ax2.bar(x - w / 2, mine_end, w, color="#1f77b4"); ax2.bar(x + w / 2, ref_end, w, color="#ff7f0e")
for xi, a, b in zip(x, mine_end, ref_end):
    ax2.text(xi - w / 2, a + 0.6, str(a), ha="center", fontsize=8)
    ax2.text(xi + w / 2, b + 0.6, str(b) if b else "none\nneeded", ha="center", fontsize=8 if b else 7)
ax2.set_xticks(x); ax2.set_xticklabels(lab, fontsize=8); ax2.set_ylim(0, 46)
ax2.set_ylabel("final bound after the lattice rounds")
ax2.set_title("B. Final bounds before the exact check\n(different routes, same conclusion: no exception)", fontsize=9.5)
mine = [416, 32, 25, 25]; ref = [506, 39, 30, 28, 28]
ax3.plot(range(1, len(mine) + 1), mine, "o-", color="#1f77b4", label="this thread (Round 17)")
ax3.plot(range(1, len(ref) + 1), ref, "s--", color="#ff7f0e", label="referee (Round 18)")
for i, val in enumerate(mine): ax3.text(i + 1 - 0.15, val * 0.8, str(val), fontsize=8, color="#1f77b4", ha="center", va="top")
for i, val in enumerate(ref): ax3.text(i + 1 + 0.15, val * 1.12, str(val), fontsize=8, color="#ff7f0e", ha="center")
ax3.set_yscale("log"); ax3.set_ylim(15, 900); ax3.set_xticks(range(1, 6))
ax3.set_xticklabels(["round %d" % i for i in range(1, 6)], fontsize=8)
ax3.set_ylabel("bound on every gap after the round (log scale)"); ax3.legend(fontsize=8)
ax3.set_title("C. Theorem 17.1 round by round, both starting near 1.1 x 10^62.\nThe referee used a 5-adic step "
              "(90 in round 1, against 105 for the 2-adic one)", fontsize=9.5)
fig.suptitle("Fig 56 (DRAFT). The independent recheck of Theorems 16.1, 16.2 and 17.1: an independent Claude session redid every step\n"
             "from the statements alone, without reading this thread's code. All three theorems hold; no conclusion changes.", fontsize=10)
fig.tight_layout(); fig.savefig("figures/fig56-independent-recheck.png", dpi=130); plt.close(fig)

# fig57: the 19-adic size of T
H = 24; grid = np.zeros((H + 1, H + 1), dtype=int)
for p in prefixes(H):
    for s in runs(H):
        h1 = p[1] + p[2]; h2 = s[1] + s[2]
        grid[h1, h2] = max(grid[h1, h2], v(T3(p, s), 19))
def gaps(Hmax):
    return [w for w in runs(Hmax)]
def T_general(ws):
    """T for N u1 x1 ... uk xk: T_1 = 1, T_j = 19^(L_1 + ... + L_(j-1)) + r_j P_j T_(j-1)"""
    T = 1; Ls = 0
    for (r, L, M) in ws:
        Ls += L; T = 19 ** Ls + r * P_of(L, M) * T
    return T
counts = {}
for k, Hm in ((1, 40), (2, 10), (3, 4), (4, 2), (5, 1)):
    mx = 0; n = 0
    for ws in itertools.product(gaps(Hm), repeat=k):
        if ws[0] == (19, 0, 0): continue          # an opening N N has one N/K letter fewer
        n += 1; mx = max(mx, v(T_general(ws), 19))
    counts[k + 1] = (mx, n, Hm)
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(13, 5.3), gridspec_kw=dict(width_ratios=[1.1, 1]))
im = ax1.imshow(grid.T, origin="lower", cmap="viridis", vmin=0, vmax=3)
cb = fig.colorbar(im, ax=ax1, ticks=[0, 1, 2, 3]); cb.set_label("largest v19(T3) over the words with these gap lengths")
ax1.set_xlabel("first gap H1"); ax1.set_ylabel("second gap H2")
ax1.set_title("A. Three letters: the 19-adic size never passes 2\n(F2, proved by hand from the leading terms mod 19)", fontsize=9.5)
ks = sorted(counts); vals = [counts[k][0] for k in ks]
cols = ["#2ca02c" if k <= 3 else "#ff7f0e" for k in ks]
ax2.bar(ks, vals, 0.6, color=cols)
for k, val in zip(ks, vals):
    ax2.text(k, val + 0.1, "%d\n(%s words,\ngaps <= %d)" % (val, format(counts[k][1], ","), counts[k][2]), ha="center", fontsize=8)
ax2.set_xticks(ks); ax2.set_xlabel("N/K letters in the word"); ax2.set_ylabel("largest 19-adic size of T found")
ax2.set_ylim(0, max(vals) + 2.2)
ax2.set_title("B. With more letters the largest 19-adic size found grows\n(green: proved, F1 and F2; orange: found by search, small gaps only)", fontsize=9.5)
fig.suptitle("Fig 57 (DRAFT). How many factors 19 the numerator T can have. Bounded for three letters, which is what lets the size\n"
             "comparison in Theorem 18.1 work without knowing the 19-adic shift; for four or more letters it is not yet bounded by any proof", fontsize=10)
fig.tight_layout(); fig.savefig("figures/fig57-nineteen-adic-size.png", dpi=130); plt.close(fig)

# fig58: the funnel for Theorem 18.1
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(13.5, 5.4), gridspec_kw=dict(width_ratios=[1, 1.35]))
XS = float(P["a_priori"]["X_star"])
seq = [XS] + [float(r["new_X"]) for r in R]
ax1.bar(range(len(seq)), seq, 0.6, color=["#7f7f7f"] + ["#1f77b4"] * len(R))
for i, val in enumerate(seq): ax1.text(i, val * 4, "%.1e" % val if val > 1e6 else "%d" % val, ha="center", fontsize=8.5)
ax1.bar([len(seq)], [1], 0.6, color="#2ca02c"); ax1.text(len(seq), 3, "exact\ncheck:\n0 clashes", ha="center", fontsize=8)
ax1.set_yscale("log"); ax1.set_ylim(0.5, 1e35); ax1.set_xticks(range(len(seq) + 1))
ax1.set_xticklabels(["a priori\n(Matveev,\n3 and 4 logs)"] + ["round %d" % (i + 1) for i in range(len(R))] + ["H1 + H2\n<= %d" % P["final_bound"]], fontsize=8)
ax1.set_ylabel("largest gap length X still possible (log scale)")
ax1.set_title("A. The bound on every gap, round by round", fontsize=10)
keys = [("m11", "first gap of one word (shared 3-log form)", "#1f77b4"),
        ("m21", "second gap, or the other's first (one 4-log form per prefix)", "#ff7f0e"),
        ("m2", "a second gap (2-adic and 5-adic D-step)", "#2ca02c"),
        ("other_word_sum", "the other word's two gaps together, from sizes", "#9467bd")]
w = 0.2; xs = np.arange(len(R))
for k, (key, lab, col) in enumerate(keys):
    vals = [r[key] for r in R]
    ax2.bar(xs + (k - 1.5) * w, vals, w, color=col, label=lab)
    for xi, val in zip(xs + (k - 1.5) * w, vals): ax2.text(xi, val + 3, str(val), ha="center", fontsize=7.5)
ax2.set_xticks(xs)
ax2.set_xticklabels(["round %d\n(start %s)" % (i + 1, "%.1e" % float(r["X"]) if r["X"] > 1e6 else r["X"]) for i, r in enumerate(R)], fontsize=8)
ax2.set_ylim(0, 340); ax2.set_ylabel("bound after the round (the new X is the largest bar)"); ax2.legend(fontsize=7.8, loc="upper right")
ax2.set_title("B. Inside each round: one of the two words is short, which caps the other", fontsize=10)
fig.suptitle("Fig 58 (DRAFT). Theorem 18.1, three N/K letters against three: from 2 x 10^29 to %d in %d rounds (the last gains nothing),\n"
             "then all %s reduced three-letter words with H1 + H2 <= %d compared exactly: every fingerprint different"
             % (P["final_bound"], len(R), format(E["three_letter_words"], ","), P["final_bound"]), fontsize=10)
fig.tight_layout(); fig.savefig("figures/fig58-three-vs-three-funnel.png", dpi=130); plt.close(fig)

# fig59: the grid
fig, ax = plt.subplots(figsize=(8.6, 6.2))
ks = ["1", "2", "3", "4+"]; G, O = "#cfe8cf", "#ffe7b3"
status = {(0, 0): ("trivial", G), (0, 1): ("proved\nThm 15.5", G), (1, 1): ("proved\nThm 16.2", G), (0, 2): ("proved\nThm 16.1", G),
          (1, 2): ("proved\nThm 17.1", G), (2, 2): ("proved\nThm 18.1\n(new)", "#9fd89f"), (0, 3): ("searched", O),
          (1, 3): ("searched", O), (2, 3): ("searched", O), (3, 3): ("searched", O)}
for (i, j), (t, c) in status.items():
    for (a, b) in ((i, j), (j, i)):
        ax.add_patch(plt.Rectangle((b, 3 - a), 1, 1, fc=c, ec="white", lw=3))
        ax.text(b + 0.5, 3 - a + 0.5, t, ha="center", va="center", fontsize=9)
ax.set_xlim(0, 4); ax.set_ylim(0, 4); ax.set_aspect("equal")
ax.set_xticks(np.arange(4) + 0.5); ax.set_xticklabels(ks); ax.set_yticks(np.arange(4) + 0.5); ax.set_yticklabels(ks[::-1])
ax.set_xlabel("N/K letters in the second word"); ax.set_ylabel("N/K letters in the first word")
ax.set_title("Fig 59 (DRAFT). Conjecture 15.3, pair type by pair type, after round 18\n(green: proved; orange: searched only). "
             "Every pair with at most three N/K letters each is now proved", fontsize=10)
fig.tight_layout(); fig.savefig("figures/fig59-where-the-conjecture-stands.png", dpi=130); plt.close(fig)
print("wrote fig56 to fig59;", "letter counts:", counts)
