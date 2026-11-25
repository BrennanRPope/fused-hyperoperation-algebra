"""Figures for the quantum-to-classical DRAFT. Reads results.json from checks.py."""
import json, math, cmath
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch, Circle

R = json.load(open("results.json"))
B, C, D, E, F = R["B"], R["C"], R["D"], R["E"], R["F"]
plt.rcParams.update({"font.size": 11, "axes.titlesize": 12, "figure.dpi": 130})
BLUE, ORANGE, GREEN, RED, PURPLE, GRAY = "#2b6cb0", "#dd6b20", "#2f855a", "#c53030", "#6b46c1", "#718096"


def save(fig, name):
    fig.tight_layout()
    fig.savefig("figures/" + name, bbox_inches="tight")
    plt.close(fig)


# fig1: along one path, composing by o_k is turning on a circle
fig, ax = plt.subplots(1, 2, figsize=(11, 4.6))
a = ax[0]
a.add_patch(Circle((0, 0), 1, fill=False, color=GRAY, lw=1))
th = np.cumsum([0, 0.6, 0.9, 0.5, 1.1])
cols = [BLUE, ORANGE, GREEN, PURPLE]
for i in range(4):
    t0, t1 = th[i], th[i + 1]
    arc = np.linspace(t0, t1, 30)
    a.plot(1.08 * np.cos(arc), 1.08 * np.sin(arc), color=cols[i], lw=4)
    a.text(1.28 * math.cos((t0 + t1) / 2), 1.28 * math.sin((t0 + t1) / 2), f"segment {i+1}\nS{i+1}/ħ",
           ha="center", va="center", fontsize=9, color=cols[i])
a.annotate("", xy=(math.cos(th[-1]), math.sin(th[-1])), xytext=(0, 0), arrowprops=dict(arrowstyle="->", lw=2, color=RED))
a.annotate("", xy=(1, 0), xytext=(0, 0), arrowprops=dict(arrowstyle="->", lw=1.2, color=GRAY))
a.text(0.35 * math.cos(th[-1]) - 0.15, 0.35 * math.sin(th[-1]) + 0.1, "whole path:\nturn = (S1+S2+S3+S4)/ħ", color=RED, fontsize=9)
a.set_xlim(-1.7, 1.7); a.set_ylim(-1.5, 1.7); a.set_aspect("equal"); a.axis("off")
a.set_title("One path: chaining segments by ∘ₖ at k = i/ħ\nis just adding their turns (the lift 1 + k·a = e^{iS/ħ})")
a = ax[1]
x = np.linspace(-2.9, 2.9, 400)
a.plot(x, 2 * np.tan(x / 2), color=BLUE, lw=2)
a.set_ylim(-8, 8); a.axhline(0, color=GRAY, lw=0.6); a.axvline(0, color=GRAY, lw=0.6)
a.set_xlabel("turn θ = S/ħ"); a.set_ylabel("chart value t = 2·tan(θ/2)")
a.set_title("Same thing in a real chart: t's compose by the tan law\n(t1 + t2)/(1 − t1·t2/4), the circle class (Cayley)")
a.text(-2.7, 5.5, "real k never gives this class:\nthe tan law is ∘ₖ at imaginary k", fontsize=9, color=RED)
save(fig, "fig1-one-path-is-just-turning.png")

# fig2: paths cancel except near Newton's
fig, ax = plt.subplots(1, 2, figsize=(12, 4.8))
a = ax[0]
for hb, col in [("0.04", ORANGE), ("0.01", BLUE)]:
    d = B["integrand"][hb]
    a.plot(d["x"], d["re"], color=col, lw=1, label=f"ħ = {hb}")
a.axvline(B["xmid_newton"], color=RED, ls="--", lw=1.5, label="Newton's path, halfway")
a.set_xlabel("where a path is halfway through the trip (x at t = T/2)")
a.set_ylabel("contribution of all paths through that point\n(real part, scaled)")
a.set_title("Every midpoint contributes about as strongly,\nbut away from Newton's point the wiggles cancel")
a.legend(fontsize=9, loc="lower left")
a = ax[1]
d = B["spiral"]["0.01"]
re, im, xs = np.array(d["re"]), np.array(d["im"]), np.array(d["x"])
sc = a.scatter(re, im, c=xs, cmap="coolwarm", s=2)
j = np.argmin(abs(xs - B["xmid_newton"]))
a.plot(re[j], im[j], "o", color="k", ms=6)
a.annotate("Newton's midpoint:\nthe running sum moves\nstraight here", (re[j], im[j]), xytext=(re[j] - 0.33, im[j] + 0.15),
           arrowprops=dict(arrowstyle="->"), fontsize=9)
a.plot(re[0], im[0], "s", color=GRAY); a.text(re[0] - 0.05, im[0] + 0.06, "start (x = −1.2)", fontsize=8)
a.plot(re[-1], im[-1], "s", color=GRAY); a.text(re[-1] - 0.1, im[-1] + 0.06, "end (x = 1.5)\n= the whole amplitude", fontsize=8)
a.set_xlabel("running sum, real part"); a.set_ylabel("running sum, imaginary part")
a.set_title("Running sum over midpoints (ħ = 0.01): a Cornu spiral.\nTight curls at both ends = cancelling paths")
a.set_aspect("equal")
cb = fig.colorbar(sc, ax=a, fraction=0.04); cb.set_label("midpoint x")
save(fig, "fig2-paths-cancel-except-newtons.png")

# fig3: Ehrenfest
fig, ax = plt.subplots(1, 2, figsize=(11.5, 4.4))
a = ax[0]
d0 = B["ehrenfest"]["0.04"]
a.plot(d0["t"], d0["xn"], color="k", lw=2.5, label="Newton's law")
for hb, col in [("0.04", ORANGE), ("0.01", BLUE)]:
    d = B["ehrenfest"][hb]
    a.plot(d["t"], d["xq"], color=col, lw=1.2, ls="--", label=f"quantum wave packet centre, ħ = {hb}")
a.set_xlabel("time"); a.set_ylabel("position")
a.set_title("A wave packet in the quartic well rides Newton's path")
a.legend(fontsize=9)
a = ax[1]
hbs = [0.04, 0.02, 0.01, 0.005]
a.loglog(hbs, B["ehrenfest_dev"], "o-", color=BLUE)
a.set_xlabel("ħ"); a.set_ylabel("largest gap from Newton (t ≤ 2)")
a.set_title(f"Gap shrinks in step with ħ (slope {B['ehrenfest_slope']:.2f})")
save(fig, "fig3-wave-packet-follows-newton.png")

# fig4: the dial plane
fig, ax = plt.subplots(1, 2, figsize=(12.5, 5.6), gridspec_kw={"width_ratios": [1.5, 1]})
a = ax[0]
a.set_xlim(-5.2, 5.2); a.set_ylim(-5.2, 5.2); a.set_aspect("equal"); a.axis("off")
a.add_patch(Circle((0, 0), 4.6, fill=False, ls="--", color=RED, lw=1.5))
a.annotate("", xy=(5, 0), xytext=(-5, 0), arrowprops=dict(arrowstyle="->", color=GRAY))
a.annotate("", xy=(0, 5), xytext=(0, -5), arrowprops=dict(arrowstyle="->", color=GRAY))
a.text(4.9, -0.45, "real k", fontsize=9, color=GRAY); a.text(0.15, 4.3, "imaginary k", fontsize=9, color=GRAY)
a.plot(0, 0, "o", color="k", ms=8)
a.text(0.2, 0.25, "k = 0: every path weighs 1\n→ counting paths (ħ → ∞)", fontsize=9)
a.plot([-4.4, -0.3], [0, 0], color=GREEN, lw=5, alpha=0.6)
a.text(-4.3, -0.75, "k = −1/(temperature): Boltzmann weights\nstatistics, heat kernel (Wick-rotated)", fontsize=9, color=GREEN)
a.plot([0, 0], [0.3, 4.4], color=BLUE, lw=5, alpha=0.6)
a.text(0.2, 2.5, "k = i/ħ: quantum\nphases e^{iS/ħ}", fontsize=9, color=BLUE)
a.annotate("", xy=(-1.9, 1.9), xytext=(-0.4, 0.5), arrowprops=dict(arrowstyle="<->", color=PURPLE, connectionstyle="arc3,rad=-0.3"))
a.text(-2.6, 1.0, "Wick rotation\n(turn the dial 90°)", fontsize=9, color=PURPLE)
a.text(-5.1, 0.25, "k → −∞:\nmin-plus\nleast action", fontsize=9, color=RED)
a.text(4.0, 0.25, "k → +∞:\nmax-plus\n(Viterbi)", fontsize=9, color=RED)
a.text(-1.5, 4.75, "k → i∞: stationary phase → Newton", fontsize=9, color=RED)
a.text(-3.4, -4.1, "the dashed ring at infinity = the classical world", fontsize=10, color=RED)
a.set_title("The path dial: weights are 1 + k·a = e^{kS}.\nEvery finite k is the same arithmetic re-charted; new rules live only at 0 and ∞")
a = ax[1]
a.set_xlim(0, 10); a.set_ylim(0, 10); a.axis("off")
a.annotate("", xy=(9.5, 7), xytext=(0.5, 7), arrowprops=dict(arrowstyle="->", lw=1.5))
a.plot(9.0, 7, "o", color="k"); a.text(7.6, 7.5, "k = 0\nGalilean\n(c → ∞)", fontsize=9)
a.text(0.5, 7.5, "k = −1/c²\nrelativity", fontsize=9)
a.text(0.5, 8.9, "Velocity dial (what's composed: speeds)", fontsize=10, weight="bold")
a.annotate("", xy=(9.5, 3), xytext=(0.5, 3), arrowprops=dict(arrowstyle="->", lw=1.5))
a.plot(0.8, 3, "o", color="k"); a.text(0.4, 3.5, "k = 0\ncounting\n(ħ → ∞)", fontsize=9)
a.text(7.2, 3.5, "|k| = 1/ħ → ∞\nclassical", fontsize=9, color=RED)
a.text(0.5, 4.9, "Path dial (what's composed: weights)", fontsize=10, weight="bold")
a.text(0.3, 0.7, "Classical sits at k = 0 on the first dial but at\n|k| = ∞ on the second: two different dials,\nnot one line from quantum to relativity.", fontsize=9.5, color=RED)
save(fig, "fig4-the-path-dial.png")

# fig5: one graph many aggregates + sweep
fig, ax = plt.subplots(1, 2, figsize=(12.5, 4.8))
a = ax[0]
a.axis("off"); a.set_xlim(-0.5, 4.5); a.set_ylim(-1.6, 2.4)
pos = {"s": (0, 0.4), "a": (2, 1.6), "b": (2, -0.8), "t": (4, 0.4)}
edges = [("s", "a", "cost 1.0\nwidth 5"), ("a", "t", "1.0+... 2.0 / 2"), ("s", "b", "0.5 / 3"), ("b", "t", "2.0 / 4"),
         ("s", "t", "3.4 / 1"), ("a", "b", "0.3 / 6")]
lab = {("s", "a"): "1.0 / 5", ("a", "t"): "2.0 / 2", ("s", "b"): "0.5 / 3", ("b", "t"): "2.0 / 4", ("s", "t"): "3.4 / 1", ("a", "b"): "0.3 / 6"}
for (u, v), l in lab.items():
    (x0, y0), (x1, y1) = pos[u], pos[v]
    a.add_patch(FancyArrowPatch((x0, y0), (x1, y1), arrowstyle="-|>", mutation_scale=14, shrinkA=14, shrinkB=14, color=GRAY))
    off = {("s", "t"): (0.6, -0.18), ("a", "b"): (0.3, 0.25)}.get((u, v), (0, 0.12))
    a.text((x0 + x1) / 2 + off[0], (y0 + y1) / 2 + off[1], l, fontsize=8, ha="center", color=BLUE)
for n_, (x0, y0) in pos.items():
    a.add_patch(Circle((x0, y0), 0.22, color="white", ec="k", zorder=3)); a.text(x0, y0, n_, ha="center", va="center", zorder=4)
rows = [("count paths (+, ×)", "4"), ("shortest (min, +)", "2.5"), ("widest (max, min)", "4"), ("reachable? (OR, AND)", "yes")]
for i, (r_, v_) in enumerate(rows):
    a.text(-0.4, -1.0 - 0.0 * i, "", fontsize=9)
a.text(-0.4, -1.45, "  •  ".join(f"{r_}: {v_}" for r_, v_ in rows[:2]) + "\n  •  ".join([""] + [f"{r_}: {v_}" for r_, v_ in rows[2:]]),
       fontsize=8.5)
a.set_title("One route map (edge labels: cost / width), four aggregates.\nAll come from one provenance polynomial over the 4 paths")
a = ax[1]
k = -np.array(D["sweep"]["k"])
Sp = np.array(D["sweep"]["Spath"])
names_p = ["s→a→t", "s→b→t", "s→t", "s→a→b→t"]
W = np.exp(-np.outer(k, Sp - Sp.min()))
share = W / W.sum(axis=1, keepdims=True)
for i, (nm, col) in enumerate(zip(names_p, [ORANGE, RED, GRAY, PURPLE])):
    a.semilogx(k, share[:, i], color=col, lw=2, label=f"{nm} (cost {Sp[i]:.1f})")
a.axhline(0.25, color="k", ls=":", lw=0.8)
a.text(0.011, 0.27, "small |k|: each path 1/4 (counting)", fontsize=8.5)
a.text(30, 0.9, "large |k|: only the\nshortest path (min-plus)", fontsize=8.5, ha="center")
a.set_xlabel("|k|  (real negative dial = 1/temperature)")
a.set_ylabel("share of the path sum carried by each path")
a.set_title("Turning one dial on the same map:\nsmall |k| counts every path, large |k| keeps only the best")
a.legend(fontsize=8.5, loc="center left")
save(fig, "fig5-one-map-many-aggregates.png")

# fig6: which candidates are genuine semirings
fig, ax = plt.subplots(figsize=(11, 4.8))
names = ["(+, ×) counting", "(min, +) least action", "(max, ×) Viterbi", "(max, min) widest", "(OR, AND) reachable",
         "(⊕ₖ, ∘ₖ) our dial, any k", "(∘₋₁, ×) union as 'or'", "soft-max with soft-min, t=4", "series/parallel bonds"]
gaps = [0, 0, 0, 0, 0, 0, float(D["union_fail"]["-1"]), float(D["bottleneck_soft_err"][1]), 0.0455]
cols = [GREEN] * 6 + [RED] * 3
ax.barh(names[::-1], [max(g, 1e-3) if g > 0 else 0 for g in gaps[::-1]], color=cols[::-1])
for i, g in enumerate(gaps[::-1]):
    ax.text(max(g, 0) + 0.003, i, "genuine semiring (distributes)" if g == 0 else f"distributive law misses by {g:.3f}", va="center", fontsize=9)
ax.set_xlim(0, 0.09)
ax.set_xlabel("gap in a·(b ⊕ c) = a·b ⊕ a·c on a sample")
ax.set_title("Which path aggregations are genuine semirings?\nGreen: yes at every value. Red: only at their end points (Boolean, or the hard limit)")
save(fig, "fig6-which-are-semirings.png")

# fig7: height flow as critical branching
fig, ax = plt.subplots(1, 2, figsize=(12, 4.6))
a = ax[0]
a.axis("off")
def grow(seed):
    r_ = np.random.default_rng(seed)
    gens_ = [[None]]
    for g_ in range(5):
        kids = []
        for pi, _ in enumerate(gens_[-1]):
            kids += [pi] * r_.poisson(1.0)
        gens_.append(kids)
    return gens_
seed = next(sd for sd in range(500) if len(grow(sd)[5]) >= 2 and max(len(g) for g in grow(sd)) <= 7 and len(grow(sd)[1]) >= 2)
gens_ = grow(seed)
xpos = [[0.5]]
for gi in range(1, 6):
    n_ = len(gens_[gi])
    xs_ = [0.08 + 0.84 * (j + 0.5) / max(n_, 1) for j in range(n_)]
    xpos.append(xs_)
    for j, par in enumerate(gens_[gi]):
        a.plot([xpos[gi - 1][par], xs_[j]], [1.0 - 0.18 * (gi - 1), 1.0 - 0.18 * gi], color=BLUE)
        a.plot(xs_[j], 1.0 - 0.18 * gi, "o", color=BLUE, ms=4)
    a.text(1.02, 1.0 - 0.18 * gi, f"generation {gi}: {n_}", fontsize=8, color=GRAY)
a.plot(0.5, 1.0, "o", color="k")
a.set_xlim(0, 1.25); a.set_ylim(-0.0, 1.05)
a.set_title("Each member has Poisson(1) children (average exactly 1).\nE(x) = e^x − 1 is this family's generating map, shifted by 1")
a = ax[1]
n = np.arange(1, len(E["surv"]) + 1)
a.loglog(n, E["surv"], color=BLUE, lw=2, label="exact: −Eⁿ(−1)")
a.loglog(np.arange(1, len(E["mc"]) + 1), E["mc"], "o", color=ORANGE, ms=4, label="40,000 simulated families")
a.loglog(n, 2 / n, "--", color=RED, label="Kolmogorov: 2/n")
a.set_xlabel("generations n (= height steps)"); a.set_ylabel("chance the family line is still alive")
a.set_title("The height flow's 'fades like 2/height'\nis a critical family tree dying out")
a.legend(fontsize=9)
save(fig, "fig7-height-flow-is-a-family-tree.png")

# fig8: decoherence
fig, ax = plt.subplots(1, 2, figsize=(12, 4.5))
a = ax[0]
xs = np.array(F["x"])
for m, col in zip(F["M"], [BLUE, GREEN, ORANGE, PURPLE, RED]):
    a.plot(xs, F["fringes"][str(m)], color=col, lw=1.4, label=f"{m} environment records")
a.set_xlabel("position on the screen"); a.set_ylabel("intensity")
a.set_title("Two-path fringes fade as more of the environment\nkeeps a faint record of which path (each record: angle 0.15)")
a.legend(fontsize=8, loc="upper right")
a = ax[1]
Ms = np.arange(0, 200)
V = np.cos(0.15) ** Ms
a.plot(Ms, V, color=BLUE, lw=2, label="fringe visibility V = ∏ cos θᵢ")
a.plot(Ms, 1 - V ** 2, color=RED, lw=2, label="which-path record 1 − V² = p₁ ∘₋₁ p₂ ∘₋₁ …")
a.set_xlabel("number of environment particles M"); a.set_ylabel("value")
a.set_title("The record composes by the union rule (k = −1):\nmany weak records add up to certainty")
a.legend(fontsize=9)
save(fig, "fig8-decoherence-is-the-union-rule.png")

# fig9: scorecard
fig, ax = plt.subplots(figsize=(11.5, 5.4))
ax.axis("off")
rows = [
    ("Composing phases along a path by ∘ₖ", "just phase addition, re-charted (lift)", "known", GRAY),
    ("Paths cancel → Newton's path (quartic well)", "yes, checked; that's the sum, not ∘ₖ", "known (Feynman, stationary phase)", GRAY),
    ("k → 0 as the classical limit", "no: classical is |k| → ∞; k = 0 counts paths", "reading ours", ORANGE),
    ("One dial line quantum → classical → relativity", "no: velocity dial and path dial are different", "ours (negative)", ORANGE),
    ("∘ₖ with partner ⊕ = a + b + 1/k is a semiring", "yes, every k; the wall −1/k is 'no path'", "ours as a statement (simple)", GREEN),
    ("Dial plane: counting at 0, statistics, quantum, tropical at ∞", "yes, a clean map", "known pieces (Maslov, Wick); map ours", GREEN),
    ("Union rule ∘₋₁ as the 'or' of path sums", "only on yes/no values; else fails", "known (reliability is hard)", GRAY),
    ("Height flow = critical family tree dying out", "yes, exactly (2/n)", "known math; link to our flow ours", GREEN),
    ("Height flow = zooming-out (renormalization)", "only generically (marginal coupling)", "see running-coupling handoff", ORANGE),
    ("Decoherence record composes by the union rule", "yes, exactly", "trivially known; ∘₋₁ reading ours", GREEN),
]
for xh, lb in [(0.0, "Question"), (0.42, "Answer"), (0.76, "Known or ours")]:
    ax.text(xh, 1.0, lb, weight="bold", transform=ax.transAxes)
for i, (q, ans, src, col) in enumerate(rows):
    y = 0.92 - i * 0.092
    ax.add_patch(plt.Rectangle((-0.01, y - 0.03), 1.02, 0.08, color=col, alpha=0.15, transform=ax.transAxes))
    ax.text(0.0, y, q, fontsize=9.5, transform=ax.transAxes)
    ax.text(0.42, y, ans, fontsize=9.5, transform=ax.transAxes)
    ax.text(0.76, y, src, fontsize=9, transform=ax.transAxes)
ax.set_title("Scorecard: what our operations add to the quantum-to-classical story (green = adds a clean reading)", pad=28)
save(fig, "fig9-scorecard.png")
print("figures written")
