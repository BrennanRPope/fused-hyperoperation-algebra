"""Pictures for the deeper pass (reads results2.json; run checks2.py first)."""
import json, math
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Circle, FancyBboxPatch, Polygon
from figures import geodesic_path, disk, interval, chord, BLUE, ORANGE, GREEN, RED, PURPLE, GRAY

R = json.load(open("results2.json")); FIG = "figures/"
TEAL = "#2a9d8f"

# ---------------------------------------------------------------- fig 7: gluing and the ladder of operations
fig = plt.figure(figsize=(15, 6.2))
ax = fig.add_axes([0.0, 0.02, 0.36, 0.86])
disk(ax, "(a) A, then B₁ and B₂ side by side.\nGluing B₁B₂ composes cross-ratios by the union rule")
A, B1, B2 = (0.15, 1.2), (1.75, 2.6), (2.6, 3.6)
interval(ax, "A", *A, BLUE); interval(ax, "B₁", *B1, ORANGE); interval(ax, "B₂", *B2, GREEN)
chord(ax, A[0], B2[1], "black"); chord(ax, A[1], B1[0], "black")
ax.text(0, -0.35, "x(A, B₁B₂) = x(A,B₁) ∘₋₁ x(A,B₂)\n= x₁ + x₂ − x₁x₂", ha="center", fontsize=12,
        bbox=dict(boxstyle="round", fc="#f4f4f4", ec="#cccccc"))
ax2 = fig.add_axes([0.39, 0.02, 0.6, 0.86]); ax2.axis("off"); ax2.set_xlim(0, 10); ax2.set_ylim(0, 6.3)
ax2.set_title("(b) Three readings of the same number, each with its own composition law", fontsize=12)
rows = [(5.0, "cross-ratio  x  (between 0 and 1)", "∘₋₁  (union rule, 'at least one')", BLUE),
        (3.0, "odds  X = φ₋₁(x) = x/(1−x)", "∘₁   (the height-flow law: E(a+b) = E(a) ∘₁ E(b))", ORANGE),
        (1.0, "length  t = −log(1−x),   X = E(t) = eᵗ − 1", "+   (plain addition)", GREEN)]
for yy, lab, law, col in rows:
    ax2.add_patch(FancyBboxPatch((0.2, yy - 0.55), 4.5, 1.1, boxstyle="round,pad=0.05", fc=col, alpha=0.15, ec=col))
    ax2.text(2.45, yy, lab, ha="center", va="center", fontsize=11)
    ax2.text(5.0, yy, "glue B₁, B₂  →  " + law, va="center", fontsize=11, color=col)
for y0, y1, lab in [(4.42, 3.58, "flip φ₋₁"), (2.42, 1.58, "X = E(t)")]:
    ax2.annotate("", xy=(2.45, y1), xytext=(2.45, y0), arrowprops=dict(arrowstyle="<->", lw=1.6))
    ax2.text(2.6, (y0 + y1) / 2, lab, va="center", fontsize=10)
ax2.text(5.0, 5.9, "Holographic mutual information  I(A:B) = (c/3) · max(0, log X)", fontsize=12, weight="bold")
fig.suptitle("Fig 7. Ours (a reading of known formulas): in AdS₃ the project's union rule, flip and height flow "
             "are exactly how entanglement adds up", fontsize=13, y=0.995)
fig.savefig(FIG + "fig7_gluing_ladder.png", bbox_inches="tight"); plt.close(fig)

# ---------------------------------------------------------------- fig 8: MMI is "o_1 beats multiplication"
t = np.linspace(0.02, 3.0, 300); T1, T2 = np.meshgrid(t, t)
E = lambda s: np.expm1(s)
mi = lambda s: np.maximum(0, np.log(E(s)))
I3 = mi(T1) + mi(T2) - mi(T1 + T2)
fig, axs = plt.subplots(1, 2, figsize=(13.5, 5.6))
im = axs[0].contourf(T1, T2, I3, levels=np.linspace(-1.2, 0, 13), cmap="Blues_r")
axs[0].contour(T1, T2, I3, levels=[-1e-9], colors="k", linewidths=0.8)
axs[0].axvline(math.log(2), color=RED, ls="--"); axs[0].axhline(math.log(2), color=RED, ls="--")
axs[0].text(math.log(2) + 0.04, 2.85, "B₁ connects to A\n(x₁ > ½, X₁ > 1)", color=RED, fontsize=9, va="top")
axs[0].text(2.85, math.log(2) + 0.06, "B₂ connects (x₂ > ½)", color=RED, fontsize=9, ha="right")
plt.colorbar(im, ax=axs[0], label="I₃(A:B₁:B₂)  in units of c/3  (never above 0)")
axs[0].set_xlabel("t₁ = −log(1 − x(A,B₁))"); axs[0].set_ylabel("t₂ = −log(1 − x(A,B₂))")
axs[0].set_title("(a) I₃ = log E(t₁)⁺ + log E(t₂)⁺ − log E(t₁+t₂)⁺ ≤ 0 everywhere", fontsize=11)
X = np.linspace(1.01, 12, 200)
for Y, col in [(1.5, ORANGE), (3, GREEN), (8, PURPLE)]:
    axs[1].plot(X, X * Y + X + Y, color=col, lw=2.5, label=f"X ∘₁ Y = XY + X + Y   (Y = {Y})")
    axs[1].plot(X, X * Y, color=col, lw=1.5, ls="--", label=f"X·Y   (Y = {Y})")
axs[1].set_xlabel("X = odds of A sharing with B₁"); axs[1].set_ylabel("odds")
axs[1].set_title("(b) Gluing gives X ∘₁ Y; MMI only needs X ∘₁ Y ≥ X·Y.\n"
                 "The slack is the '+X+Y' of the circle operation: I₃ = −log(1 + 1/X + 1/Y)", fontsize=11)
axs[1].legend(fontsize=8.5, ncol=1)
fig.suptitle("Fig 8. For three neighbouring regions in AdS₃, monogamy of mutual information is the circle "
             "operation ∘₁ beating plain multiplication", fontsize=13)
fig.tight_layout(); fig.savefig(FIG + "fig8_mmi_from_circle.png"); plt.close(fig)

# ---------------------------------------------------------------- fig 9: dial sweeps
fig, axs = plt.subplots(1, 3, figsize=(16, 5))
s = R["I_ksweep"]
axs[0].plot(s["k"], s["graph_neg"], color=PURPLE, lw=2.5, label="min-cut graphs (600)")
axs[0].plot(s["k"], s["haar_neg"], color=GRAY, lw=2, label="Haar random states")
axs[0].plot(s["k"], s["sparse_neg"], color=ORANGE, lw=2, label="sparse random states")
axs[0].axvline(-1, color=RED, ls="--"); axs[0].text(-0.97, 0.45, "k = −1", color=RED)
axs[0].set_xlabel("dial k"); axs[0].set_ylabel("fraction with corner sum ≤ 0")
axs[0].set_title("(a) Only at k = −1 is the cube sum ≤ 0 for every geometry", fontsize=11)
axs[0].legend(fontsize=9)
so = R["I_soft"]
axs[1].plot(so["h"], so["mmi_viol"], "o-", color=PURPLE, lw=2.5, label="positive weights (ferromagnet)")
axs[1].scatter([1.0], [R["I_soft_signed_viol"]], color=RED, s=80, zorder=3, label="mixed-sign weights, h = 1")
axs[1].set_xlabel("softness h   (h = 0: min cut, the tropical limit;  h = 1/log D)")
axs[1].set_ylabel("fraction breaking MMI"); axs[1].set_ylim(-0.03, 0.4)
axs[1].set_title("(b) Softening the min cut (log-sum-exp) never broke MMI\n"
                 "or the five-region rules, while weights stay positive", fontsize=11)
axs[1].legend(fontsize=9)
rn = R["I_renyi"]; al = [a if a < 1e5 else 8 for a in rn["alpha"]]
axs[2].plot(al, rn["sparse_viol"], "o-", color=ORANGE, lw=2, label="sparse random states")
axs[2].plot(al, rn["haar_viol"], "o-", color=GRAY, lw=2, label="Haar random states")
axs[2].set_xticks([0.25, 1, 2, 3, 5, 8]); axs[2].set_xticklabels(["¼", "1\n(von Neumann)", "2", "3", "5", "∞"])
axs[2].set_xlabel("Rényi index α"); axs[2].set_ylabel("fraction breaking MMI")
axs[2].set_title("(c) No Rényi index makes generic quantum states\npass; the test is about the state, not the entropy", fontsize=11)
axs[2].legend(fontsize=9)
fig.suptitle("Fig 9. Turning the dials away from k = −1 and away from the sharp min cut", fontsize=13)
fig.tight_layout(); fig.savefig(FIG + "fig9_dial_sweeps.png"); plt.close(fig)

# ---------------------------------------------------------------- fig 10: rule length vs k
g = R["G"]; fig, ax = plt.subplots(figsize=(9, 4.8))
for name, col in [("SA", ORANGE), ("MMI", BLUE), ("MMI uplift", GREEN), ("cyclic (5)", PURPLE)]:
    ax.plot(g["k"], g["counts"][name], "o-", ms=3, color=col, lw=2, label=name)
ax.axvline(-1, color=RED, ls="--"); ax.axvline(0, color=GRAY, ls=":")
ax.text(-0.97, 27, "k = −1\n(multipartite\ninformation)", color=RED, fontsize=9)
ax.text(0.03, 27, "k = 0\n(plain entropies)", color=GRAY, fontsize=9)
ax.set_xlabel("dial k of the corner basis"); ax.set_ylabel("number of terms in the rule")
ax.set_title("Fig 10. Rules get short only at k = −1 or k = 0; every other dial value spreads them out")
ax.legend(); fig.tight_layout(); fig.savefig(FIG + "fig10_rule_length.png"); plt.close(fig)

# ---------------------------------------------------------------- fig 11: first law and relative entropy
j = R["J"]; y = np.array(j["y"])
fig, axs = plt.subplots(1, 2, figsize=(14, 5.2))
axs[0].plot(y, j["dH"], color=BLUE, lw=2.5, label="Δ⟨H⟩: energy seen by the region's modular Hamiltonian")
axs[0].plot(y, j["dS"], color=GREEN, lw=2.5, label="ΔS: change of entanglement entropy (RT, BTZ)")
axs[0].fill_between(y, j["dS"], j["dH"], color=ORANGE, alpha=0.25, label="gap = relative entropy ≥ 0")
axs[0].set_xlim(0, 4); axs[0].set_ylim(0, 3)
axs[0].set_xlabel("region size × temperature  (y = πℓ/β)"); axs[0].set_ylabel("units of c/3")
axs[0].set_title("(a) First law: ΔS = Δ⟨H⟩ for small regions; the gap is\nrelative entropy, which gravity reads as positive energy", fontsize=11)
axs[0].legend(fontsize=9, loc="upper left")
ax = axs[1]; ax.axis("off"); ax.set_xlim(0, 10); ax.set_ylim(0, 6)
boxes = [(5.0, "Entanglement first law on every ball:  δS = δ⟨H⟩\n(checked here from RT in d = 2, 3, 4, 5, 6)", BLUE),
         (3.0, "Linearized Einstein equations in the bulk\n(Faulkner, Guica, Hartman, Myers, Van Raamsdonk 2013)", GREEN),
         (1.0, "Next order: relative entropy ≥ 0 (the shaded gap in (a))\n→ gravitational positive-energy conditions (Lashkari et al. 2016)", ORANGE)]
for yy, txt, col in boxes:
    ax.add_patch(FancyBboxPatch((0.4, yy - 0.7), 9.2, 1.4, boxstyle="round,pad=0.05", fc=col, alpha=0.13, ec=col))
    ax.text(5, yy, txt, ha="center", va="center", fontsize=11)
ax.annotate("", xy=(5, 3.75), xytext=(5, 4.28), arrowprops=dict(arrowstyle="<->", lw=2))
ax.text(5.15, 4.0, "equivalent (FGHMV)", fontsize=9.5, va="center")
ax.annotate("", xy=(5, 1.75), xytext=(5, 2.28), arrowprops=dict(arrowstyle="->", lw=2))
ax.text(5.15, 2.0, "next order", fontsize=9.5, va="center")
ax.set_title("(b) Where this meets Einstein's equations", fontsize=11)
fig.suptitle("Fig 11. The concrete North Star link: entanglement bookkeeping is Einstein's equations", fontsize=13)
fig.tight_layout(); fig.savefig(FIG + "fig11_first_law.png"); plt.close(fig)

# ---------------------------------------------------------------- fig 12: GHZ
k = np.linspace(-2, 1, 300); fig, ax = plt.subplots(figsize=(8, 5))
for n, col in zip(range(2, 7), [ORANGE, GREEN, BLUE, PURPLE, RED]):
    ax.plot(k, ((1 + k) ** n - 1) / np.where(np.abs(k) < 1e-9, 1e-9, k), color=col, lw=2.2,
            label=f"n = {n}:  1 ∘ₖ 1 ∘ₖ … ∘ₖ 1")
ax.scatter([-1], [1], color="black", s=60, zorder=3); ax.text(-0.95, 1.4, "every n passes through (−1, 1):\nGHZ breaks MMI by exactly 1 bit", fontsize=10)
ax.axhline(0, color="black", lw=0.8); ax.axvline(-1, color=RED, ls="--"); ax.set_ylim(-3, 8)
ax.set_xlabel("dial k"); ax.set_ylabel("cube sum of the GHZ entropies")
ax.set_title("Fig 12. GHZ entropies are the all-ones product, so its cube sum\nis the n-fold circle product of 1 (the q-number [n])")
ax.legend(fontsize=9); fig.tight_layout(); fig.savefig(FIG + "fig12_ghz.png"); plt.close(fig)
print("figures2 written")
