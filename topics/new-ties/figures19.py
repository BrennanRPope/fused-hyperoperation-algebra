"""Round 19 pictures (DRAFT): fig60 the 19-adic size of T for four letters (at most 3, F3) and how it escapes at five;
fig61 the bounds for Theorems 19.1 to 19.3, round by round; fig62 the case tree for four letters against four;
fig63 where Conjecture 15.3 stands now. Usage: python3 figures19.py (needs pairs19.json)."""
import json
import numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from pairs17 import v
from pairs19 import box_words, hsum, T_of
P = json.load(open("pairs19.json")); E = P["exact_check"]

# fig60
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(13, 4.8))
hist = {}
for w in box_words(3, (14, 14, 14)):
    if hsum(w) <= 14:
        g = v(T_of(w), 19); hist[g] = hist.get(g, 0) + 1
ks = sorted(hist)
ax1.bar(ks, [hist[k] for k in ks], color="#2ca02c")
for k in ks: ax1.text(k, hist[k] * 1.3, f"{hist[k]:,}", ha="center", fontsize=9)
ax1.axvline(3.5, color="red", ls="--"); ax1.text(3.55, max(hist.values()) / 30, "F3: never more\nthan 3", color="red", fontsize=9)
ax1.set_yscale("log"); ax1.set_xticks(range(0, 5)); ax1.set_xlim(-0.6, 4.6); ax1.set_ylim(1, max(hist.values()) * 8)
ax1.set_xlabel("number of factors 19 in T (its 19-adic size g)"); ax1.set_ylabel("four-letter words (log scale)")
ax1.set_title("A. Four letters: g <= 3 for every word (proved by hand, F3);\nshown: all %s words with gap sum <= 14" %
              f"{sum(hist.values()):,}", fontsize=9.5)
Ms = list(range(0, 6001)); gs = []
for M in Ms: gs.append(v(T_of(((10, 0, 0), (10, 1, 0), (10, 0, 0), (19, 0, M))), 19))
ax2.plot(np.array(Ms) + 1, gs, ".", color="#bbbbbb", ms=2, label="g for each M")
ax2.step(np.array(Ms) + 1, np.maximum.accumulate(gs), where="post", color="#ff7f0e", lw=2, label="largest g so far")
ax2.legend(fontsize=8, loc="upper left")
for M, g in ((0, 3), (38, 4), (5092, 5)):
    ax2.plot(M + 1, g, "ko"); ax2.annotate("m^%d: g = %d" % (M, g), (M + 1, g), textcoords="offset points", xytext=(6, -14), fontsize=9)
ax2.set_xscale("log"); ax2.set_ylim(-0.3, 5.8); ax2.set_yticks(range(6))
ax2.set_xlabel("1 + M, the number of m letters before the last N (log scale)"); ax2.set_ylabel("g = factors 19 in T")
ax2.set_title("B. Five letters: N K l K K m^M N. The leading units cancel mod 19,\nso g = 2 + v19(M - mu) keeps growing: "
              "no fixed bound", fontsize=9.5)
fig.suptitle("Fig 60 (DRAFT). The 19-adic size of T: bounded for four letters, unbounded for five", fontsize=11)
fig.tight_layout(); fig.savefig("figures/fig60-nineteen-adic-four-and-five.png", dpi=130); plt.close(fig)

# fig61
fig, axs = plt.subplots(1, 3, figsize=(15, 4.8))
titles = {"42": "A. Four letters against two (Thm 19.1)", "43": "B. Four against three (Thm 19.2)",
          "44": "C. Four against four (Thm 19.3)"}
for ax, key in zip(axs, ("42", "43", "44")):
    R = P["pairs"][key]["rounds"]; ys = [float(P["pairs"][key]["a_priori"]["X_star"])] + [float(r["new_X"]) for r in R]
    xs = list(range(len(ys)))
    ax.plot(xs, ys, "o-", color="#1f77b4")
    for x, y in zip(xs, ys): ax.text(x, y * 3, ("%.1e" % y) if y > 1e4 else str(int(y)), ha="center", fontsize=8)
    ax.set_yscale("log"); ax.set_ylim(10, max(ys) * 1e6)
    ax.set_xticks(xs); ax.set_xticklabels(["Matveev\nstart"] + ["round %d" % i for i in range(1, len(R) + 1)], fontsize=7.5)
    ax.set_ylabel("bound on every gap (log scale)"); ax.set_title(titles[key], fontsize=10)
axs[2].text(0.98, 0.95, "then every short-box word against all\n%s four-letter words with gap sum <= %d:\nno clash" %
            (f"{E['four_letter_words']:,}", E["four_letter_streamed_up_to"]), transform=axs[2].transAxes, ha="right", va="top", fontsize=8.5,
            bbox=dict(fc="#eef6ee", ec="#2ca02c"))
fig.suptitle("Fig 61 (DRAFT). From Matveev's bound to an exact check: each round of lattice and 2-adic steps cuts the bound "
             "to about its logarithm, until the size step stops it", fontsize=10.5)
fig.tight_layout(); fig.savefig("figures/fig61-four-letter-funnels.png", dpi=130); plt.close(fig)

# fig62: the case tree for four against four (final round's numbers)
R = P["pairs"]["44"]["rounds"][-1]; a, b, c, d, f1, f2 = (R[k] for k in ("a", "b", "c", "d", "f1", "f2"))
fig, ax = plt.subplots(figsize=(14, 7.2)); ax.axis("off"); ax.set_xlim(0, 14); ax.set_ylim(0, 7.4)
def box(x, y, t, c="#e8f0fa", w=2.5):
    ax.text(x, y, t, ha="center", va="center", fontsize=8.6, bbox=dict(boxstyle="round,pad=0.35", fc=c, ec="#555"))
def arrow(x1, y1, x2, y2, t=""):
    ax.annotate("", (x2, y2 + 0.28), (x1, y1 - 0.28), arrowprops=dict(arrowstyle="->", color="#555"))
    if t: ax.text((x1 + x2) / 2 + 0.08, (y1 + y2) / 2, t, fontsize=7.8, color="#333")
G = "#d8efd8"
box(7, 7.0, "W = N x1 x2 x3,  W' = N y1 y2 y3, same fingerprint\nP(1,1), one lattice: min(x1, y1) <= a = %d; call W the one with x1 <= a" % a, w=6)
box(7, 6.0, "S(3,3), the D-step: min(x3, y3) <= c = %d" % c)
arrow(7, 7.0, 7, 6.0)
box(3.5, 5.0, "A: x3 <= c\nP(2,1), lattice per x1: min(x2, y1) <= b = %d" % b); box(10.5, 5.0, "B: y3 <= c\nP(2,1): min(x2, y1) <= b = %d" % b)
arrow(7, 6.0, 3.5, 5.0); arrow(7, 6.0, 10.5, 5.0)
box(1.6, 3.8, "x2 <= b:\nW short (a, b, c)", G); box(4.8, 3.8, "y1 <= b: S(2,3), 2-adic\nwith x3 known: min(x2, y3) <= d = %d" % d)
arrow(3.5, 5.0, 1.6, 3.8); arrow(3.5, 5.0, 4.8, 3.8)
box(2.9, 2.6, "x2 <= d:\nW short (a, d, c)", G); box(6.0, 2.6, "y3 <= d: S(2,2), x3 and y3 known:\nmin(x2, y2) <= f1 = %d" % f1)
arrow(4.8, 3.8, 2.9, 2.6); arrow(4.8, 3.8, 6.0, 2.6)
box(4.4, 1.4, "W short (a, f1, c)\nor W' short (b, f1, d)", G); arrow(6.0, 2.6, 4.4, 1.4)
box(9.0, 3.8, "B1: x2 <= b: S(3,3) with y3 known:\nmin(x3, y2) <= d"); box(12.4, 3.8, "B2: y1 <= b: S(3,3):\nmin(x3, y2) <= d")
arrow(10.5, 5.0, 9.0, 3.8); arrow(10.5, 5.0, 12.4, 3.8)
box(8.1, 2.6, "x3 <= d: W short (a, b, d)\ny2 <= d: S(3,1), y2 y3 known:\nmin(x3, y1) <= f2 = %d" % f2)
box(12.4, 2.6, "y2 <= d: W' short (b, d, c)\nx3 <= d: S(2,2):\nmin(x2, y2) <= f1")
arrow(9.0, 3.8, 8.1, 2.6); arrow(12.4, 3.8, 12.4, 2.6)
box(8.1, 1.4, "W short (a, b, f2)\nor W' short (f2, d, c)", G); box(12.4, 1.4, "W short (a, f1, d)\nor W' short (b, f1, c)", G)
arrow(8.1, 2.6, 8.1, 1.4); arrow(12.4, 2.6, 12.4, 1.4)
ax.text(7, 0.35, "Green: one word is short, so its gap sum is at most %d. The size step then caps the other word's gap sum at %d: "
        "10^3 20^(sum) <= T(other) <= 548,720 T(short)." % (R["short_sum"], R["other_bound"]), ha="center", fontsize=9)
ax.set_title("Fig 62 (DRAFT). The case tree for four N/K letters against four, with the final round's bounds. P = lattice step "
             "on a known front part, S = 2-adic step on a known back part", fontsize=10)
fig.tight_layout(); fig.savefig("figures/fig62-four-vs-four-tree.png", dpi=130); plt.close(fig)

# fig63: the grid
fig, ax = plt.subplots(figsize=(9.4, 7.2))
ks = ["1", "2", "3", "4", "5+"]; G, GN, O = "#cfe8cf", "#9fd89f", "#ffe7b3"
status = {(0, 0): ("trivial", G), (0, 1): ("Thm 15.5", G), (1, 1): ("Thm 16.2", G), (0, 2): ("Thm 16.1", G),
          (1, 2): ("Thm 17.1", G), (2, 2): ("Thm 18.1", G), (0, 3): ("Thm 19.1\n(new)", GN), (1, 3): ("Thm 19.1\n(new)", GN),
          (2, 3): ("Thm 19.2\n(new)", GN), (3, 3): ("Thm 19.3\n(new)", GN)}
for j in range(5): status[(j, 4)] = ("searched;\nneeds a\n19-adic step", O)
for (i, j), (t, c) in status.items():
    for (a_, b_) in ((i, j), (j, i)):
        ax.add_patch(plt.Rectangle((b_, 4 - a_), 1, 1, fc=c, ec="white", lw=3))
        ax.text(b_ + 0.5, 4 - a_ + 0.5, t, ha="center", va="center", fontsize=8.6)
ax.set_xlim(0, 5); ax.set_ylim(0, 5); ax.set_aspect("equal")
ax.set_xticks(np.arange(5) + 0.5); ax.set_xticklabels(ks); ax.set_yticks(np.arange(5) + 0.5); ax.set_yticklabels(ks[::-1])
ax.set_xlabel("N/K letters in the second word"); ax.set_ylabel("N/K letters in the first word")
ax.set_title("Fig 63 (DRAFT). Conjecture 15.3, pair type by pair type, after round 19\n(green: proved, computer-assisted; "
             "orange: searched only; five letters need a bound the size step cannot give)", fontsize=9.5)
fig.tight_layout(); fig.savefig("figures/fig63-where-the-conjecture-stands.png", dpi=130); plt.close(fig)
print("figures written")
