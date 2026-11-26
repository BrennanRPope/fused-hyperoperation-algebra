"""Pictures for the third pass (DRAFT): the shadow proof of the gluing law, and the fuzzy-cut proof.
Reads results3.json (run checks3.py first)."""
import json, math, itertools
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, Rectangle
from scipy.special import logsumexp
from figures import BLUE, ORANGE, GREEN, RED, PURPLE, GRAY
R = json.load(open("results3.json")); FIG = "figures/"
TEAL = "#2a9d8f"

# ---------------------------------------------------------------- fig 13: the shadow proof
a1, a2, b1, m, b2 = 0.0, 4.0, 4.8, 6.0, 9.0
f = lambda p: (a1 - p) / (p - a2)   # negative for p > a2; |f| is the height of p's foot                 # orientation-preserving Moebius map: a1 -> 0, a2 -> infinity
finv = lambda z: (a1 + a2 * z) / (z + 1)
def semicircle(u, w, n=300):
    th = np.linspace(0, np.pi, n); c, r = (u + w) / 2, (w - u) / 2; return c + r * np.cos(th), r * np.sin(th)
def perp(p, n=200):                               # geodesic from boundary point p meeting gamma_A at right angles
    rr = abs(f(p)); th = np.linspace(np.pi, np.pi / 2, n); z = finv(rr * np.exp(1j * th)); return z.real, z.imag
fig, axs = plt.subplots(1, 2, figsize=(16, 7.2))
ax = axs[0]
ax.axhline(0, color="black", lw=1)
x, y = semicircle(a1, a2); ax.plot(x, y, color=BLUE, lw=3, label="γ_A, the minimal surface of A")
for (u, w, col, lab) in [(b1, m, ORANGE, "B₁"), (m, b2, GREEN, "B₂")]:
    ax.plot([u, w], [0, 0], color=col, lw=7, solid_capstyle="butt"); ax.text((u + w) / 2, -0.35, lab, color=col, ha="center", fontsize=14, weight="bold")
ax.plot([a1, a2], [0, 0], color=BLUE, lw=7, solid_capstyle="butt"); ax.text(2, -0.35, "A", color=BLUE, ha="center", fontsize=14, weight="bold")
feet = {}
for p, lab in [(b1, "b₁"), (m, "m"), (b2, "b₂")]:
    px, py = perp(p); ax.plot(px, py, "--", color=GRAY, lw=1.4); feet[lab] = (px[-1], py[-1])
    ax.plot(px[-1], py[-1], "o", color="black", ms=5); ax.text(p, 0.08, lab, ha="center", fontsize=11)
def arc_between(F1, F2, col, lab, off):
    th = np.linspace(np.arctan2(F1[1], F1[0] - 2), np.arctan2(F2[1], F2[0] - 2), 100)
    ax.plot(2 + 2 * np.cos(th), 2 * np.sin(th), color=col, lw=7, alpha=0.9, solid_capstyle="butt")
    tm = th.mean(); ax.text(2 + 2.45 * np.cos(tm) + off, 2.45 * np.sin(tm), lab, color=col, fontsize=14, weight="bold", ha="center")
arc_between(feet["b₁"], feet["m"], ORANGE, "t₁", 0.05); arc_between(feet["m"], feet["b₂"], GREEN, "t₂", 0.05)
ax.set_xlim(-0.6, 9.6); ax.set_ylim(-3.2, 5.6); ax.set_aspect("equal"); ax.axis("off")
ax.set_title("(a) Drop B's endpoints onto A's minimal surface (blue) at right angles.\nB₁'s 'shadow' has length t₁, B₂'s has t₂; they meet at m's foot,\nso the shadow of B₁B₂ has length t₁ + t₂", fontsize=11)
ax.text(-0.4, -2.9, "Lemma (shadow):  1 − x(A,B) = e^(−t(A,B))\nGluing:  shadows add, so 1 − x multiplies:\n"
        "1 − x₁₂ = (1 − x₁)(1 − x₂)  ⇔  x₁₂ = x₁ ∘₋₁ x₂", fontsize=11.5, bbox=dict(boxstyle="round", fc="#f4f4f4", ec="#cccccc"))
ax = axs[1]
ax.axhline(0, color="black", lw=1); ax.plot([0, 0], [0, 7.1], color=BLUE, lw=3)
ax.text(0.2, 7.3, "A's minimal surface becomes the vertical axis\n(Möbius map f: a₁ → 0, a₂ → ∞)", color=BLUE, fontsize=11)
hs = {}
for p, lab, col in [(b1, "b₁", "black"), (m, "m", "black"), (b2, "b₂", "black")]:
    rr = abs(f(p)); th = np.linspace(np.pi, np.pi / 2, 100); ax.plot(rr * np.cos(th), rr * np.sin(th), "--", color=GRAY, lw=1.4)
    ax.plot(0, rr, "o", color="black", ms=5); ax.text(-rr, -0.45, f"f({lab})", ha="center", fontsize=11); hs[lab] = rr
for (u, w, col, lab) in [(hs["b₁"], hs["m"], ORANGE, "t₁ = log h(b₁)/h(m)"), (hs["m"], hs["b₂"], GREEN, "t₂ = log h(m)/h(b₂)")]:
    ax.plot([0, 0], [u, w], color=col, lw=7, alpha=0.85)
    ax.text(0.15, (u + w) / 2, lab, color=col, fontsize=12, va="center", weight="bold")
ax.plot([-hs["b₁"], -hs["m"]], [0, 0], color=ORANGE, lw=7, solid_capstyle="butt"); ax.plot([-hs["m"], -hs["b₂"]], [0, 0], color=GREEN, lw=7, solid_capstyle="butt")
ax.text(-7.6, -3.6, "On the vertical axis, hyperbolic length is log(height ratio).\nThe perpendicular from boundary point f(p) < 0 lands at height h(p) = |f(p)|.\n"
        "So t(A,B) = log h(b₁)/h(b₂), and the algebra\n1 − x(A,B) = h(b₂)/h(b₁)  (Ptolemy)  finishes the lemma.\n"
        "Adding lengths = multiplying height ratios (telescoping).", fontsize=11, bbox=dict(boxstyle="round", fc="#f4f4f4", ec="#cccccc"))
ax.set_xlim(-7.8, 4.5); ax.set_ylim(-3.9, 8.2); ax.set_aspect("equal"); ax.axis("off")
ax.set_title("(b) The same picture after the Möbius map.\nA shadow is a segment of the vertical axis", fontsize=11)
fig.suptitle("Fig 13. Ours: a picture proof of the gluing law. Cross-ratios glue by ∘₋₁ because shadows on A's minimal surface add",
             fontsize=13, y=1.0)
fig.savefig(FIG + "fig13_shadow_proof.png", bbox_inches="tight", dpi=130); plt.close(fig)

# ---------------------------------------------------------------- fig 14: one distance, three rungs
d = np.linspace(0.05, 4.0, 400)
x = 1 / np.cosh(d / 2) ** 2; X = 1 / np.sinh(d / 2) ** 2; t = -2 * np.log(np.tanh(d / 2)); MI = np.maximum(0, np.log(X))
dstar = 2 * np.arcsinh(1)
fig, axs = plt.subplots(1, 2, figsize=(14, 5.4))
ax = axs[0]
ax.plot(d, x, color=BLUE, lw=2.5, label="x = sech²(δ/2)   (rung ∘₋₁: cross-ratio)")
ax.plot(d, np.minimum(X, 4), color=ORANGE, lw=2.5, label="X = csch²(δ/2)   (rung ∘₁: odds)")
ax.plot(d, np.minimum(t, 4), color=GREEN, lw=2.5, label="t = −2 log tanh(δ/2)   (rung +: shadow length)")
ax.plot(d, np.minimum(MI, 4), color=RED, lw=2.5, ls="--", label="I(A:B) in units of c/3 = max(0, log X)")
ax.axvline(dstar, color=GRAY, ls=":"); ax.text(dstar + 0.05, 1.6, "δ* = 2 arcsinh 1 ≈ 1.76:\nbeyond this, A and B\nshare no entanglement\n(x = ½, X = 1, t = log 2)", fontsize=9.5)
ax.set_ylim(0, 4); ax.set_xlabel("δ = hyperbolic distance between the minimal surfaces of A and B"); ax.legend(fontsize=9.5, loc="upper right", bbox_to_anchor=(1, 0.72))
ax.set_title("(a) Every number on the ladder is a function of one distance δ", fontsize=11)
ax = axs[1]; ax.axis("off"); ax.set_xlim(0, 10); ax.set_ylim(0, 6)
rows = [(5.0, "x = sech²(δ/2)", "glue: x₁ ∘₋₁ x₂  (complements multiply)", BLUE),
        (3.0, "X = x/(1−x) = csch²(δ/2)", "glue: X₁ ∘₁ X₂  (the height-flow law)", ORANGE),
        (1.0, "t = −log(1−x)", "glue: t₁ + t₂  (shadows add)", GREEN)]
for yy, lab, law, col in rows:
    ax.add_patch(FancyBboxPatch((0.2, yy - 0.5), 4.2, 1.0, boxstyle="round,pad=0.05", fc=col, alpha=0.15, ec=col))
    ax.text(2.3, yy, lab, ha="center", va="center", fontsize=11.5); ax.text(4.7, yy, law, va="center", fontsize=11, color=col)
for y0, y1, lab in [(4.45, 3.55, "flip φ₋₁"), (2.45, 1.55, "X = eᵗ − 1")]:
    ax.annotate("", xy=(2.3, y1), xytext=(2.3, y0), arrowprops=dict(arrowstyle="<->", lw=1.6)); ax.text(2.45, (y0 + y1) / 2, lab, va="center", fontsize=10)
ax.text(0.2, 5.75, "Proved at the bottom rung (shadows add); the other two follow by the flip", fontsize=11, weight="bold")
ax.set_title("(b) Why one proof is enough", fontsize=11)
fig.suptitle("Fig 14. Ours (reading standard hyperbolic geometry): the cross-ratio, its odds and the shadow length are one distance in three coordinates", fontsize=12.5, y=1.0)
fig.savefig(FIG + "fig14_three_rungs.png", bbox_inches="tight", dpi=130); plt.close(fig)

# ---------------------------------------------------------------- fig 15: fuzzy cuts never break monogamy
rng = np.random.default_rng(7)
def soft_I3(W, nb, h):
    F = {}
    for T in itertools.product((0, 1), repeat=3):
        vals = []
        for b in itertools.product((0, 1), repeat=nb):
            lab = np.array(list(T) + [0] + list(b)); vals.append(-(W * (lab[:, None] != lab[None, :])).sum() / 2 / h)
        F[T] = -h * logsumexp(vals)
    S = {T: F[T] - F[(0, 0, 0)] for T in F}
    return sum((-1) ** (sum(T) + 1) * S[T] for T in S if sum(T))
hs = np.logspace(-2, 1.3, 40)
fig, axs = plt.subplots(1, 2, figsize=(14, 5.3))
ax = axs[0]
for i in range(12):
    nb = int(rng.integers(1, 5)); N = 4 + nb; W = rng.exponential(1, (N, N)) * (rng.random((N, N)) < 0.7); W = np.triu(W, 1); W = W + W.T
    ax.plot(hs, [soft_I3(W, nb, h) for h in hs], lw=1.6, alpha=0.85)
ax.axhline(0, color="black", lw=1.2); ax.set_xscale("log")
ax.set_xlabel("softness h   (h → 0: sharp min cut;  h = 1/log D: tensor network, bond dimension D)", fontsize=9.5)
ax.set_ylabel("I₃ of the fuzzy entropies")
ax.set_title("(a) 12 random graphs: I₃ stays ≤ 0 at every softness (now proved)", fontsize=11)
ax = axs[1]
pq = R["K"]["potts_star_best_I3"]; f5 = R["K"]["five_party_worst"]
labels = ["Ising\n(2 states)"] + [f"Potts\n{q} states" for q in ("3", "8", "20")]
vals = [pq["2"], pq["3"], pq["8"], pq["20"]]
cols = [TEAL if v < 1e-9 else RED for v in vals]
ax.bar(range(4), [max(v, 1e-16) for v in vals], color=cols)
ax.set_yscale("log"); ax.set_ylim(1e-16, 1e-1); ax.set_xticks(range(4)); ax.set_xticklabels(labels)
ax.axhline(1e-9, color=GRAY, ls=":"); ax.text(3.45, 2e-9, "numerical zero", fontsize=9, ha="right")
ax.set_ylabel("largest I₃ found (one hub, optimised weights)")
ax.set_title("(b) The two-state symmetry matters. One hub: 8 or 20 states break\nmonogamy (red). Two hubs: 3 states already break it (check K12b)", fontsize=10.5)
fig.suptitle("Fig 15. Ours: softened (fuzzy) min cuts keep the monogamy rule at every softness, for every graph with positive weights", fontsize=12.5, y=1.0)
fig.savefig(FIG + "fig15_fuzzy_cuts.png", bbox_inches="tight", dpi=130); plt.close(fig)

# ---------------------------------------------------------------- fig 16: the four-copy proof as a picture
fig = plt.figure(figsize=(16, 7.4)); ax = fig.add_axes([0, 0, 1, 0.92]); ax.axis("off"); ax.set_xlim(0, 16); ax.set_ylim(0, 7.4)
def pattern(x0, y0, s, title, col):
    for i, sg in enumerate(s):
        ax.add_patch(Rectangle((x0 + 0.42 * i, y0), 0.38, 0.38, fc=(ORANGE if sg > 0 else "#264653"), ec="white"))
    ax.text(x0 + 0.84, y0 + 0.55, title, ha="center", fontsize=9, color=col)
even = [(1, 1, 1, 1), (1, 1, -1, -1), (1, -1, 1, -1), (1, -1, -1, 1)]
odd = [(1, 1, 1, -1), (1, 1, -1, 1), (1, -1, 1, 1), (-1, 1, 1, 1)]
ax.text(0.3, 6.9, "Step 1. Monogamy ⇔ one inequality between 8 pinnings of the four terminals A, B, C, O", fontsize=12, weight="bold")
for k, s in enumerate(even): pattern(0.4 + 1.9 * k, 5.6, s, "even", BLUE)
for k, s in enumerate(odd): pattern(0.4 + 1.9 * k, 4.4, s, "odd", RED)
ax.text(0.4, 5.25, "Z(even 1)·Z(even 2)·Z(even 3)·Z(even 4)   ≤   Z(odd 1)·Z(odd 2)·Z(odd 3)·Z(odd 4)", fontsize=11)
ax.text(0.4, 3.95, "(orange = +, dark = −; the order is A, B, C, O)", fontsize=9, color=GRAY)
ax.text(8.4, 6.9, "Step 2. Run four copies at once: a product of 4 Z's\nis one Z of a 4-copy magnet", fontsize=12, weight="bold")
ax.text(8.4, 5.3, "Stack the 4 even pinnings: terminal j gets the j-th column\nof the 4×4 Hadamard matrix. Stack the odd ones: the 'odd frame'\n"
        "v₁ = (+,+,+,−), v₂ = (+,+,−,+), v₃ = (+,−,+,+), v₄ = (−,+,+,+).\nNegating copy 4 everywhere (a symmetry) turns the Hadamard frame\n"
        "into (v₁, v₂, v₃, −v₄). So the claim is:  Z(v₁,v₂,v₃,+v₄) ≥ Z(v₁,v₂,v₃,−v₄).", fontsize=10.5)
ax.text(0.3, 3.3, "Step 3. Rotate each 4-copy spin into the frame:  αᵢ = vᵢ·σ / 2", fontsize=12, weight="bold")
ax.text(0.4, 1.9, "The energy σ·σ' = Σ αᵢαᵢ' is unchanged (orthogonal frame).\nTerminal j now pushes only coordinate αⱼ (it sits at 2eⱼ).\n"
        "The 16 spin values become 8 'even' sign vectors (α₁α₂α₃α₄ = +1)\nand the 8 points ±2eᵢ. Every moment of that set is ≥ 0.", fontsize=10.5)
ax.text(8.4, 3.3, "Step 4. Griffiths' first inequality finishes it", fontsize=12, weight="bold")
ax.text(8.4, 1.5, "Z(+v₄) − Z(−v₄) = 2·Z'·⟨sinh(push on α₄)⟩, and sinh is a series of odd\npowers with positive coefficients. With ferromagnetic couplings, positive\n"
        "pushes on α₁, α₂, α₃ and nonnegative moments, Griffiths' first inequality\nmakes every such average ≥ 0.  Picture: A, B, C push α₁, α₂, α₃ up;\n"
        "the parity rule α₄ = α₁α₂α₃ then leans α₄ up too. So Z(+v₄) ≥ Z(−v₄).  ∎", fontsize=10.5)
for yy in (3.6,): ax.plot([0.3, 15.7], [yy, yy], color="#dddddd", lw=1)
ax.plot([8.1, 8.1], [0.3, 7.2], color="#dddddd", lw=1)
fig.suptitle("Fig 16. Ours (with the Ellis–Monroe–Newman four-copy trick): why fuzzy min cuts can't break monogamy", fontsize=13.5, y=0.99)
fig.savefig(FIG + "fig16_four_copy_proof.png", bbox_inches="tight", dpi=130); plt.close(fig)

# ---------------------------------------------------------------- fig 17: five-party rules, numerically
fig, ax = plt.subplots(figsize=(10, 4.6))
names = list(f5.keys()); vals = [max(f5[n], 1e-16) for n in names]
ax.bar(range(len(names)), vals, color=TEAL); ax.set_yscale("log"); ax.set_ylim(1e-17, 1e-1)
ax.axhline(1e-9, color=GRAY, ls=":"); ax.text(len(names) - 0.5, 2e-9, "numerical zero", ha="right", fontsize=9)
ax.set_xticks(range(len(names))); ax.set_xticklabels(names, fontsize=10)
ax.set_ylabel("largest violation found\n(optimised weights, normalised)")
ax.set_title("Fig 17. Evidence, not proof: the five-region holographic rules also survive fuzzy cuts\n"
             "(12 optimised graphs per rule; bars at the floor = nothing positive found; the eighth rule was not retrieved)", fontsize=11)
fig.savefig(FIG + "fig17_five_party.png", bbox_inches="tight", dpi=130); plt.close(fig)
print("figures 13-17 written")
