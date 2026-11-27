"""Round 17 pictures (DRAFT): fig53 the bounds for Theorem 17.1, round by round; fig54 the 2-adic step (why a long second gap
would need K 19^m to be 1 to many binary places, and how close it actually comes); fig55 where Conjecture 15.3 stands now."""
import json, math
import numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from pairs17 import v, jmax, strip, T2, T3, name2, name3
P = json.load(open("pairs17.json")); R = P["rounds"]

# fig53: the funnel
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(13.5, 5.4), gridspec_kw=dict(width_ratios=[1, 1.35]))
XS = float(P["a_priori"]["X_star"]); E = P["exact_check"]
seq = [XS] + [float(r["new_X"]) for r in R]
lab1 = ["a priori\n(Matveev,\n3, 4, 5 logs)"] + ["round %d" % (i + 1) for i in range(len(R))]
ax1.bar(range(len(seq)), seq, 0.6, color=["#7f7f7f"] + ["#1f77b4"] * len(R))
for i, val in enumerate(seq): ax1.text(i, val * 4, "%.1e" % val if val > 1e6 else "%d" % val, ha="center", fontsize=8.5)
ax1.bar([len(seq)], [1], 0.6, color="#2ca02c")
ax1.text(len(seq), 3, "exact\ncheck:\n0 clashes", ha="center", fontsize=8)
ax1.set_yscale("log"); ax1.set_ylim(0.5, 1e68); ax1.set_xticks(range(len(seq) + 1)); ax1.set_xticklabels(lab1 + ["below\n%d" % P["final_bound"]], fontsize=8)
ax1.set_ylabel("largest gap length X still possible (log scale)")
ax1.set_title("A. The bound on every gap, round by round", fontsize=10)
keys = [("V_short_hV", "V short: V's gap H' (shared size Lambda3, 2-adic D-step)", "#1f77b4"),
        ("H1_bound_V_short", "V short: W's first gap H1 (one lattice per fingerprint)", "#ff7f0e"),
        ("H2_bound_V_short", "V short: W's second gap H2 (exact 2-adic logarithm)", "#2ca02c"),
        ("W_short_Hprime_bound", "W short: V's gap H' from sizes alone", "#9467bd")]
w = 0.2; xs = np.arange(len(R))
for k, (key, lab, col) in enumerate(keys):
    vals = [r[key] for r in R]
    ax2.bar(xs + (k - 1.5) * w, vals, w, color=col, label=lab)
    for xi, val in zip(xs + (k - 1.5) * w, vals): ax2.text(xi, val + 4, str(val), ha="center", fontsize=7.5)
ax2.set_xticks(xs); ax2.set_xticklabels(["round %d\n(start %s)" % (i + 1, "%.1e" % float(r["X"]) if r["X"] > 1e6 else r["X"]) for i, r in enumerate(R)], fontsize=8)
ax2.set_ylim(0, 470); ax2.set_ylabel("bound after the round (the new X is the largest bar)")
ax2.legend(fontsize=7.8, loc="upper right")
ax2.set_title("B. Inside each round: either V is short or W is short", fontsize=10)
fig.suptitle("Fig 53 (DRAFT). Theorem 17.1, three N/K letters against two: from 10^62 to %d in %d rounds (the last one gains nothing),\n"
             "then %s three-letter words against %s two-letter words compared exactly: no shared fingerprint"
             % (P["final_bound"], len(R), format(E["three_letter_words"], ","), format(E["two_letter_words"], ",")), fontsize=10)
fig.tight_layout(); fig.savefig("figures/fig53-three-vs-two-funnel.png", dpi=130); plt.close(fig)

# fig54: the 2-adic step for one fingerprint
Vw = (19, 3, 2); K = strip(T2(Vw)); M = 400
ms = np.arange(-M, M + 1)
vals = [v(abs(K * 19 ** int(m) - 1) if m >= 0 else abs(K - 19 ** int(-m)), 2) for m in ms]
fig, ax = plt.subplots(figsize=(11, 4.8))
ax.vlines(ms, 0, vals, color="#2ca02c", lw=0.8)
j = jmax(K, M)
ax.axhline(j, color="k", ls="--", lw=1)
ax.text(-M + 5, j + 0.4, "largest possible over |m| <= %d: %d binary places (exact 2-adic logarithm), so a matching\n"
        "three-letter word has second gap H2 <= %d" % (M, j, j // 2), fontsize=8.5)
m0 = 1 - Vw[1]
ax.annotate("m = g' - L' = %d: the two-letter word's own shift.\nHere K 19^m - 1 = r'P'/19^L', exactly v2(r'P') = %d places.\n"
            "(The D-step rules this m out; the bound\nincludes it anyway.)" % (m0, j), xy=(m0, j), xytext=(-150, 21.3), fontsize=8,
            arrowprops=dict(arrowstyle="->", lw=0.8))
for H2 in (6, 8, 10):
    ax.axhline(2 * H2, color="#d62728", lw=1, alpha=0.6)
    ax.text(M - 5, 2 * H2 + 0.2, "a second gap of %d letters needs %d places" % (H2, 2 * H2), ha="right", fontsize=8, color="#d62728")
ax.set_xlim(-M, M); ax.set_ylim(0, 25.5)
ax.set_xlabel("m = (19-adic size of T3) - L1 - L2, the unknown shift")
ax.set_ylabel("binary places where K 19^m and 1 agree")
ax.set_title("Fig 54 (DRAFT). The 2-adic step in Theorem 17.1, for the fingerprint K = %d of %s.\nA three-letter word with this "
             "fingerprint and second gap H2 needs K 19^m = 1 to 2 H2 binary places (red). Away from one spike, it never gets close (green)" % (K, name2(Vw)), fontsize=9.5)
fig.tight_layout(); fig.savefig("figures/fig54-two-adic-step.png", dpi=130); plt.close(fig)

# fig55: the grid
fig, ax = plt.subplots(figsize=(8.6, 6.2))
ks = ["1", "2", "3", "4+"]
G, O = "#cfe8cf", "#ffe7b3"
status = {(0, 0): ("trivial", G), (0, 1): ("proved\nThm 15.5", G), (1, 1): ("proved\nThm 16.2", G), (0, 2): ("proved\nThm 16.1", G),
          (1, 2): ("proved\nThm 17.1\n(new)", "#9fd89f"), (2, 2): ("searched\n(0 clashes)", O), (0, 3): ("searched", O), (1, 3): ("searched", O),
          (2, 3): ("searched", O), (3, 3): ("searched", O)}
for (i, j), (t, c) in status.items():
    for (a, b) in ((i, j), (j, i)):
        ax.add_patch(plt.Rectangle((b, 3 - a), 1, 1, fc=c, ec="white", lw=3))
        ax.text(b + 0.5, 3 - a + 0.5, t, ha="center", va="center", fontsize=9)
ax.set_xlim(0, 4); ax.set_ylim(0, 4); ax.set_aspect("equal")
ax.set_xticks(np.arange(4) + 0.5); ax.set_xticklabels(ks); ax.set_yticks(np.arange(4) + 0.5); ax.set_yticklabels(ks[::-1])
ax.set_xlabel("N/K letters in the second word"); ax.set_ylabel("N/K letters in the first word")
ax.set_title("Fig 55 (DRAFT). Conjecture 15.3, pair type by pair type, after round 17\n(green: proved; orange: searched only)", fontsize=10)
fig.tight_layout(); fig.savefig("figures/fig55-where-the-conjecture-stands.png", dpi=130); plt.close(fig)
print("wrote fig53, fig54, fig55")
