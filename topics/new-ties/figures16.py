"""Round 16 pictures (DRAFT): fig50 why no 3-letter word has fingerprint 1 (how close powers of 19 come to 2^x 5^y, against
how close they would have to come); fig51 the bounds, from Matveev to LLL to the exact check; fig52 where Conjecture 15.3
stands, pair type by pair type."""
import json, math
import numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from decimal import Decimal, getcontext
getcontext().prec = 80
J = json.load(open("baker16.json")); P = json.load(open("pairs16.json"))
l2, l5, l19 = (Decimal(p).ln() for p in (2, 5, 19))

# fig50
Es = np.arange(4, 61); closest = []; need = []
for E in Es:
    best = None
    for y in range(0, int(E * 1.83) + 3):
        t = Decimal(int(E)) * l19 - y * l5
        if t < 0: break
        x = int((t / l2).to_integral_value())
        lam = abs(t - x * l2)
        if lam > 0 and (best is None or lam < best): best = lam
    closest.append(-float(best.log10())); need.append(float(E) * math.log10(19) - math.log10(2 * (1 + 380 * E * E)))
fig, ax = plt.subplots(figsize=(10, 5.2))
ax.plot(Es, closest, "o-", ms=4, color="#1f77b4", label="how close 19^E comes to the nearest 2^x 5^y: -log10 |E log 19 - x log 2 - y log 5| (illustration)")
ax.plot(Es, need, "-", lw=2.5, color="#d62728", label="how close it would have to come for a 3-letter word with psi = 19^E")
ax.axvline(J["case_II_final_bound"], color="k", ls="--", lw=1)
ax.text(J["case_II_final_bound"] + 0.5, 45, "past E = %d the red line is out of reach\nfor every E, not just these (LLL + Matveev);\nbelow it, every case is checked exactly" % J["case_II_final_bound"], fontsize=9)
ax.set_xlabel("E, the power of 19"); ax.set_ylabel("closeness, in decimal digits")
ax.set_ylim(0, 60); ax.legend(fontsize=8.5, loc="upper left")
ax.set_title("Fig 50 (DRAFT). Theorem 16.1: a 3-letter word with fingerprint 1 would need powers of 19 and of 2 and 5 to agree\nto about E digits; they never agree to more than a few (the blue points), so none exists")
fig.tight_layout(); fig.savefig("figures/fig50-why-no-three-letter-unit.png", dpi=130); plt.close(fig)

# fig51
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(13, 4.9))
labels = ["Matveev", "LLL 1", "LLL 2", "LLL 3", "exact\ncheck"]
c2 = [J["matveev_bounds"]["case_II_E"]] + [r["new_bound"] for r in J["lll_rounds"] if r["case"] == "II"][:3]
c1 = [J["matveev_bounds"]["case_I_H"]] + [r["new_bound"] for r in J["lll_rounds"] if r["case"] == "I"][:3]
while len(c1) < 4: c1.append(c1[-1])
x = np.arange(4)
ax1.bar(x - 0.2, c2, 0.4, color="#1f77b4", label="Case II: the power E (psi = 19^E)")
ax1.bar(x + 0.2, c1, 0.4, color="#ff7f0e", label="Case I: the run length H (word ends in ...N with no gap)")
for xi, a, b in zip(x, c2, c1):
    ax1.text(xi - 0.2, a * 1.5, format(a, ",") if a < 1e6 else "%.1e" % a, ha="center", fontsize=8)
    ax1.text(xi + 0.2, b * 1.5, format(b, ",") if b < 1e6 else "%.1e" % b, ha="center", fontsize=8)
ax1.bar([4], [1], 0.6, color="#2ca02c"); ax1.text(4, 2, "%d candidates,\n0 solutions" % J["exact_check"]["candidates_checked"], ha="center", fontsize=8)
ax1.set_yscale("log"); ax1.set_ylim(0.5, 1e15); ax1.set_xticks(range(5)); ax1.set_xticklabels(labels)
ax1.set_ylabel("largest value still possible (log scale)"); ax1.legend(fontsize=8)
ax1.set_title("A. Theorem 16.1 (3 letters vs fingerprint 1):\nfrom 5 trillion to 10 in three lattice reductions")
fb = [p["final_H_bound"] for p in P["worst_small_words"]]
ax2.bar([0], [P["matveev"]["H_big_bound"]], 0.5, color="#9467bd")
ax2.bar([1], [P["largest_final_bound"]], 0.5, color="#9467bd")
ax2.bar([2], [1], 0.5, color="#2ca02c")
ax2.text(0, P["matveev"]["H_big_bound"] * 1.6, "%.1e" % P["matveev"]["H_big_bound"], ha="center", fontsize=8)
ax2.text(1, P["largest_final_bound"] * 1.6, "<= %d for every one\nof the %d small words" % (P["largest_final_bound"], P["small_words"]), ha="center", fontsize=8)
ax2.text(2, 2, "%d words,\n0 shared\nfingerprints" % P["exact_check"]["words"], ha="center", fontsize=8)
ax2.set_yscale("log"); ax2.set_ylim(0.5, 1e20); ax2.set_xticks(range(3))
ax2.set_xticklabels(["Matveev (4 logs)", "LLL, one lattice\nper small word", "exact check"])
ax2.set_ylabel("largest run length H still possible (log scale)")
ax2.set_title("B. Theorem 16.2 (all pairs of 2-letter words): the 2-adic step\nmakes one word small (H <= %d), then the same funnel" % P["matveev"]["H_small_bound"])
fig.suptitle("Fig 51 (DRAFT). The bounds: linear forms in logarithms (Matveev 2000), lattice reduction (de Weger's method), then exact arithmetic", fontsize=10.5)
fig.tight_layout(); fig.savefig("figures/fig51-bounds-funnel.png", dpi=130); plt.close(fig)

# fig52
fig, ax = plt.subplots(figsize=(8.6, 6.2))
ks = ["1", "2", "3", "4+"]
status = {(0, 0): ("trivial", "#cfe8cf"), (0, 1): ("proved\nThm 15.5", "#cfe8cf"), (1, 1): ("proved\nThm 16.2", "#cfe8cf"),
          (0, 2): ("proved\nThm 16.1", "#cfe8cf"), (1, 2): ("searched\n(0 clashes)", "#ffe7b3"), (2, 2): ("searched\n(0 clashes)", "#ffe7b3"),
          (0, 3): ("searched", "#ffe7b3"), (1, 3): ("searched", "#ffe7b3"), (2, 3): ("searched", "#ffe7b3"), (3, 3): ("searched", "#ffe7b3")}
for (i, j), (t, c) in status.items():
    for (a, b) in ((i, j), (j, i)):
        ax.add_patch(plt.Rectangle((b, 3 - a), 1, 1, fc=c, ec="white", lw=3))
        ax.text(b + 0.5, 3 - a + 0.5, t, ha="center", va="center", fontsize=9)
ax.set_xlim(0, 4); ax.set_ylim(0, 4); ax.set_aspect("equal")
ax.set_xticks(np.arange(4) + 0.5); ax.set_xticklabels(ks); ax.set_yticks(np.arange(4) + 0.5); ax.set_yticklabels(ks[::-1])
ax.set_xlabel("N/K letters in the second word"); ax.set_ylabel("N/K letters in the first word")
ax.set_title("Fig 52 (DRAFT). Conjecture 15.3, pair type by pair type\n(green: proved; orange: searched only: every word to 14 letters,\nand 2-4 N/K letters with long l, m runs)", fontsize=10)
fig.tight_layout(); fig.savefig("figures/fig52-where-the-conjecture-stands.png", dpi=130); plt.close(fig)
print("wrote fig50, fig51, fig52")
