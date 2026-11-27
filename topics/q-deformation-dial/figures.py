"""Figures for the q-deformation vs plane-dial page. Run: python3 figures.py"""
import numpy as np, matplotlib
matplotlib.use("Agg"); import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch
BLUE, ORANGE, GREEN, RED, GRAY, PURPLE = "#2b6cb0", "#dd6b20", "#2f855a", "#c53030", "#4a5568", "#6b46c1"

# fig1: q-integers add by the line law
fig, axs = plt.subplots(1, 2, figsize=(13, 4.8), gridspec_kw={"width_ratios": [1.25, 1]})
ax = axs[0]; q = 1.35
vals = [sum(q**j for j in range(n)) for n in range(0, 7)]
for n, v in enumerate(vals):
    ax.scatter(v, 0, s=80, color=BLUE, zorder=3); ax.text(v, 0.12, f"[{n}]", ha="center", fontsize=11, color=BLUE)
m, n = 2, 3
ax.annotate("", xy=(vals[m+n], -0.35), xytext=(0, -0.35), arrowprops=dict(arrowstyle="-|>", color=ORANGE, lw=2))
ax.text(vals[m+n]/2, -0.55, f"[2] ∘ₖ [3] = [5]   with k = q − 1 = {q-1:.2f}", ha="center", fontsize=11, color=ORANGE)
ax.text(0, 0.55, "q-integers [n]_q = 1 + q + … + q^(n−1)  (here q = 1.35)", fontsize=11)
ax.text(0, -0.95, "[m + n] = [m] + [n] + (q − 1)[m][n]: the line dial, nothing else  (Q1, Lean qnum_add)", fontsize=10)
ax.set_ylim(-1.15, 0.8); ax.set_xlim(-0.3, vals[-1]+0.6); ax.axis("off")
ax = axs[1]
ax.text(0.02, 0.86, "Hoffman–Ihara's q letter product (§6.4):", fontsize=11.5)
ax.text(0.06, 0.70, "zᵢ ◇ zⱼ = z₍ᵢ₊ⱼ₎ + (1 − q)·z₍ᵢ₊ⱼ₋₁₎", fontsize=13, color=PURPLE)
ax.text(0.02, 0.52, "on the q-MZV letters fᵢ(n) = q^((i−1)n) / [n]^i the\ncorrection term is just", fontsize=10.5)
ax.text(0.06, 0.33, "qⁿ = 1 + (q − 1)·[n]_q", fontsize=13, color=ORANGE)
ax.text(0.02, 0.12, "so 1 − q = −k: Hoffman–Ihara's q is our k  (Q2, Q2b,\nLean q_letter_product)", fontsize=10.5)
ax.axis("off")
fig.suptitle("Fig 1. Hoffman–Ihara's q is the line dial's k (q = 1 + k), not d", fontsize=12.5)
fig.tight_layout(); fig.savefig("fig1-q-is-k.png", dpi=150); plt.close(fig)

# fig2: the dial plane, bases from eigenvalues
fig, ax = plt.subplots(figsize=(8.6, 6.6))
K = np.linspace(-3, 3, 400)
ax.fill_between(K, K**2/4, 3, color=PURPLE, alpha=0.08)
ax.plot(K, K**2/4, color=PURPLE, lw=2, label="k² = 4d: the two bases are equal")
ax.axhline(0, color=ORANGE, lw=4, alpha=0.8, label="d = 0: one base is 1 (plain q-numbers, Hoffman–Ihara)")
ax.text(-2.9, 2.5, "complex pair of bases\n(P, Q = conjugates: turning)", color=PURPLE, fontsize=10)
ax.text(-2.9, -1.6, "two real bases P ≠ Q, neither 1\n(p,q)-numbers, Chakrabarti–Jagannathan 1991", color=GREEN, fontsize=10)
pts = {"union (−1, 0)": (-1, 0), "parity (−2, 0)": (-2, 0), "velocity (0, −1)": (0, -1), "tangent (0, 1)": (0, 1)}
for lab, (a, b) in pts.items():
    ax.scatter(a, b, s=60, color="black", zorder=4); ax.text(a+0.08, b+0.12, lab, fontsize=9.5)
ax.text(0.6, -2.4, "N has eigenvalues α, β with k = α + β, d = αβ;\nbases P = 1 + x₀α, Q = 1 + x₀β  (Q4b, Q6)", fontsize=10)
ax.set_xlabel("k"); ax.set_ylabel("d"); ax.set_xlim(-3, 3); ax.set_ylim(-2.8, 3); ax.grid(alpha=0.3)
ax.legend(loc="upper right", fontsize=9)
ax.set_title("Fig 2. On the plane dial, d switches on the second base", fontsize=12)
fig.tight_layout(); fig.savefig("fig2-dial-bases.png", dpi=150); plt.close(fig)

# fig3: orbits, one base vs two bases
fig, axs = plt.subplots(1, 2, figsize=(13, 4.8))
x0 = 0.3; nn = np.arange(0, 9)
for (k, d, lab, c) in [(1.0, 0.0, "(k, d) = (1, 0): one base", ORANGE), (1.0, -0.5, "(k, d) = (1, −0.5)", GREEN), (1.0, 0.2, "(k, d) = (1, 0.2)", BLUE)]:
    disc = k*k - 4*d; al = (k + np.sqrt(disc))/2; be = (k - np.sqrt(disc))/2
    P, Q = 1 + x0*al, 1 + x0*be
    axs[0].plot(nn, P**nn, "-o", color=c, ms=4, label=f"{lab}: P = {P:.2f}")
    axs[0].plot(nn, Q**nn, "--s", color=c, ms=4, label=f"{' '*len(lab)}  Q = {Q:.2f}")
    A = (al*Q**nn - be*P**nn)/(al - be) if al != be else None
    B = (P**nn - Q**nn)/(al - be)
    axs[1].plot(nn, B/A, "-o", color=c, label=lab)
axs[0].set_title("The two bases Pⁿ and Qⁿ (dashed). d = 0 keeps Q = 1", fontsize=11)
axs[0].set_xlabel("n"); axs[0].legend(fontsize=8); axs[0].grid(alpha=0.3)
axs[1].set_title("The orbit xₙ = Bₙ/Aₙ of 1 + x₀N; x₍ₘ₊ₙ₎ = F(xₘ, xₙ)  (Q4, Lean orbit_add)", fontsize=11)
axs[1].set_xlabel("n"); axs[1].legend(fontsize=9); axs[1].grid(alpha=0.3)
fig.suptitle("Fig 3. One base gives q-numbers; a second base, set by d, gives the plane law", fontsize=12.5)
fig.tight_layout(); fig.savefig("fig3-one-base-two-bases.png", dpi=150); plt.close(fig)
print("wrote 3 figures")
