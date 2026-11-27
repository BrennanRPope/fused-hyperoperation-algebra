"""Round 20 pictures (DRAFT): fig64 our bounds against the independent referee's for Theorems 18.1 and 19.1 to 19.3;
fig65 the 19-adic step for the five-letter family (exact records against log_19 M, and the published bound);
fig66 the 12 cancelling patterns and how far the exact maximum sits above log_19 B; fig67 the road to five letters.
Usage: python3 figures20.py (needs adic20.json)."""
import json, math
import numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
import adic20 as A

J = json.load(open("adic20.json"))
LOG19 = math.log(19)

# fig64: ours against the referee's
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(13, 4.8))
names = ["18.1\n3 vs 3", "19.1\n4 vs 2", "19.2\n4 vs 3", "19.3\n4 vs 4"]
ours_start = [None, 2.142e29, 2.081e62, 4.823e112]
ref_start = [6.438e45, 6.438e45, 2.775e62, 1.550e79]
x = np.arange(4)
ax1.bar(x[1:] - 0.18, [math.log10(s) for s in ours_start[1:]], 0.36, color="#1f77b4", label="ours (round 18: no single start)")
ax1.bar(x + 0.18, [math.log10(s) for s in ref_start], 0.36, color="#ff7f0e", label="referee")
for i in range(4):
    if ours_start[i]: ax1.text(i - 0.18, math.log10(ours_start[i]) + 2, "%.0e" % ours_start[i], ha="center", fontsize=8)
    ax1.text(i + 0.18, math.log10(ref_start[i]) + 2, "%.0e" % ref_start[i], ha="center", fontsize=8)
ax1.set_xticks(x); ax1.set_xticklabels(names); ax1.set_ylabel("log10 of the starting bound (Matveev)")
ax1.set_ylim(0, 125); ax1.legend(fontsize=8, loc="upper left")
ax1.set_title("A. Starting bounds from Matveev's theorem. Different splits,\nso different starts; both are finite, which is all a start needs",
              fontsize=9.5)
rounds_ours = {"19.1": [369, 61, 50, 50], "19.2": [269, 68, 59, 55, 55], "19.3": [371, 74, 66, 61, 59, 59]}
rounds_ref = {"18.1": [422, 42, 33], "19.1": [558, 50, 38], "19.2": [765, 59, 47, 45], "19.3": [263, 61, 51, 49]}
cols = {"18.1": "#9467bd", "19.1": "#2ca02c", "19.2": "#d62728", "19.3": "#1f77b4"}
for k, r in rounds_ours.items():
    ax2.plot(range(1, len(r) + 1), r, "-o", color=cols[k], label="ours " + k)
for k, r in rounds_ref.items():
    ax2.plot(range(1, len(r) + 1), r, "--s", color=cols[k], mfc="white", label="referee " + k)
ax2.set_yscale("log"); ax2.set_xlabel("lattice round (after the first huge one)"); ax2.set_ylabel("bound (log scale)")
ax2.legend(fontsize=7, ncol=2)
ax2.set_title("B. Round by round. Ours bound every gap, the referee's the gap sum;\neach side's exact check covers its own last bound: no clash",
              fontsize=9.5)
fig.suptitle("Fig 64 (DRAFT). An independent session redid Theorems 18.1 and 19.1 to 19.3 from the statements alone: all hold",
             fontsize=11)
fig.tight_layout(); fig.savefig("figures/fig64-independent-recheck-18-19.png", dpi=130); plt.close(fig)

# fig65: the family
ties, bad = A.casework(4, 5)
fam = next(p for p in bad if p["rs"] == (10, 10, 10, 19))
A3, u = A.u_of(A.pattern_front(fam, ()))
mu = A.mu_digits(u, A.K_for(10 ** 300))
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(13, 4.8))
ks = range(1, 230)
recM = [mu % 19 ** k for k in ks]          # first M with v19(M - mu) >= k
pts = [(M, 2 + k) for M, k in zip(recM, ks) if M > 0]
# keep only genuine records (a later k can share the same M)
best = {}
for M, g in pts: best[M] = max(best.get(M, 0), g)
Ms = sorted(best); gs = [best[M] for M in Ms]
lx = [math.log10(M) for M in Ms]
ax1.plot(lx, gs, "o", color="#ff7f0e", ms=3.5, label="exact records: the first M where g reaches each value")
xx = np.linspace(0, 300, 200)
ax1.plot(xx, 2 + xx / math.log10(19), "k--", lw=1, label="2 + log_19 M (a random 19-adic mu would do this)")
ax1.set_xlabel("log10 M (the number of m letters before the last N)"); ax1.set_ylabel("g = factors 19 in T")
ax1.legend(fontsize=8, loc="upper left"); ax1.set_xlim(0, 305)
ax1.set_title("A. N K l K K m^M N: g grows, but only like log M.\nExact, from mu's 19-adic digits (checked by brute force to M = 10^5)",
              fontsize=9.5)
Bexp = [2, 13, 50, 100, 200, 300]
bl = [A.bl_bound_g(A3, u, 10 ** e) for e in Bexp]
ex = [2 + A.max_v_upto(mu, A.K_for(10 ** e), 10 ** e)[0] for e in Bexp]
ax2.semilogy(Bexp, bl, "-o", color="#7f7f7f", label="published two-log bound (Bugeaud-Laurent, as restated)")
ax2.semilogy(Bexp, [18 * b for b in bl], ":", color="#7f7f7f", label="x 18 (conservative reading of its constant)")
ax2.semilogy(Bexp, ex, "-o", color="#2ca02c", label="exact maximum over all M <= B (the digit step)")
for e, b, g in zip(Bexp, bl, ex):
    ax2.annotate("%d" % g, (e, g), textcoords="offset points", xytext=(4, -12), fontsize=8, color="#2ca02c")
ax2.set_xlabel("log10 B, the largest M allowed"); ax2.set_ylabel("bound on g (log scale)"); ax2.legend(fontsize=8)
ax2.set_title("B. The published bound holds for every M but is about 10^4 times too big;\nthe digit step gives the exact answer for M up to any B",
              fontsize=9.5)
fig.suptitle("Fig 65 (DRAFT). The 19-adic step on the family that broke the factor-19 bound", fontsize=11)
fig.tight_layout(); fig.savefig("figures/fig65-nineteen-adic-step-family.png", dpi=130); plt.close(fig)

# fig66: the 12 patterns and the excess over log_19 B
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(13, 4.8), gridspec_kw=dict(width_ratios=[1.05, 1]))
ax1.axis("off")
rows = [(p["word"].split(" ")[0].split("m^M")[0], p["A3"], p["free_front_gaps"]) for p in J["patterns"]["list"]]
ax1.text(0.0, 1.0, "The 12 five-letter patterns whose leading units cancel mod 19", fontsize=10, weight="bold", va="top")
ax1.text(0.0, 0.93, "front (N-items may also carry m's), then m^M N;   g = A3 + 1 + v19(M - mu(front))", fontsize=8.5, va="top")
for i, (w, a3, nf) in enumerate(rows):
    y = 0.85 - i * 0.068
    ax1.text(0.02, y, w + "  m^M N", family="monospace", fontsize=9, va="top",
             color="#d62728" if w == "NKlKK" else "black")
    ax1.text(0.62, y, "A3 = %d" % a3, fontsize=9, va="top")
    ax1.text(0.78, y, "%d free front gap%s" % (nf, "" if nf == 1 else "s"), fontsize=9, va="top")
ax1.text(0.02, 0.85 - 12 * 0.068, "red: the family of round 19.  All 12: tied terms 0, 1, 4; the last N bare.", fontsize=8.5,
         va="top", color="#555555")
H = J["reduction"]["by_B"]
labs = ["10^2", "10^100", "10^200", "10^300"]
keys = ["2", "100", "200", "300"]
cats = ["-1", "0", "1", "2", "3", "4"]
bottom = np.zeros(4)
colors = ["#c7e9c0", "#74c476", "#fdae6b", "#e6550d", "#a63603", "#7f2704"]
for c, col in zip(cats, colors):
    vals = np.array([H[k]["hist"].get(c, 0) for k in keys], float)
    ax2.bar(labs, vals, bottom=bottom, color=col, label="k_max - ceil(log_19 B) = %s" % c)
    bottom += vals
for i, k in enumerate(keys):
    ax2.text(i, 13900, "max g %d" % H[k]["max_g"], ha="center", fontsize=9)
ax2.set_ylim(0, 15000); ax2.set_ylabel("fronts (13,517: all with free front gaps <= 30)")
ax2.set_xlabel("B, the largest M allowed"); ax2.legend(fontsize=7.5, loc="upper left", bbox_to_anchor=(1.01, 1.0))
ax2.set_title("B. For every front, the exact maximum of v19(M - mu) over M <= B\nis within 4 of log_19 B: g only grows like log M",
              fontsize=9.5)
fig.suptitle("Fig 66 (DRAFT). Every unbounded case at five letters has one shape, and the digit step tames it", fontsize=11)
fig.tight_layout(); fig.savefig("figures/fig66-twelve-patterns.png", dpi=130); plt.close(fig)

# fig67: the road to five letters
fig, ax = plt.subplots(figsize=(15, 4.8)); ax.axis("off"); ax.set_xlim(0, 13); ax.set_ylim(0, 4.6)
boxes = [
    (0.2, "Casework:\n12 patterns,\none 19-adic\nform each", "done (round 20)", "#2ca02c"),
    (2.35, "Published\ntwo-log bound:\ng <= C x height(front)\nx log(M)^2", "done (round 20,\nconstants as restated)", "#2ca02c"),
    (4.5, "Matveev a priori\nbound for 5 vs 1..5\n(g only as a\ncoefficient)", "not done", "#7f7f7f"),
    (6.65, "Lattice rounds\non the fronts,\n2-adic steps\non the backs", "not done\n(same kind as round 19)", "#7f7f7f"),
    (8.8, "Digit step\nper front:\ng <= A3 + 1\n+ log_19 B + small", "done for the family\nand fronts <= 30", "#2ca02c"),
    (10.95, "Size step and\nexact check\n(~10^10 words,\n100-300 x round 19)", "not done\n(the big part)", "#7f7f7f"),
]
for x0, txt, st, col in boxes:
    ax.add_patch(plt.Rectangle((x0, 1.9), 1.9, 1.6, fc="white", ec=col, lw=2.2))
    ax.text(x0 + 0.95, 2.85, txt, ha="center", va="center", fontsize=8.8)
    ax.text(x0 + 0.95, 1.55, st, ha="center", va="top", fontsize=8.2, color=col)
for x0, _, _, _ in boxes[:-1]:
    ax.annotate("", xy=(x0 + 2.12, 2.7), xytext=(x0 + 1.93, 2.7), arrowprops=dict(arrowstyle="->", lw=1.5))
ax.text(6.5, 4.25, "Fig 67 (DRAFT). What a proof for five letters needs, and what this round settled", ha="center", fontsize=11)
ax.text(6.5, 3.85, "Green: done.  Grey: not done.  The new 19-adic layer is the two green boxes on either side; the grey ones are the "
        "round-19 method again, bigger.", ha="center", fontsize=9, color="#444444")
ax.text(6.5, 0.45, "Judgement, not proved: no new idea is needed for five letters, only a larger run. All lengths at once would still need "
        "something uniform.", ha="center", fontsize=9, style="italic")
fig.savefig("figures/fig67-road-to-five-letters.png", dpi=130, bbox_inches="tight"); plt.close(fig)
print("wrote fig64 to fig67")
