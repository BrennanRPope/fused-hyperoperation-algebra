"""Figures for the two-base letter product section. Run: python3 two_base_figures.py"""
import numpy as np, matplotlib
matplotlib.use("Agg"); import matplotlib.pyplot as plt
BLUE, ORANGE, GREEN, RED, GRAY, PURPLE = "#2b6cb0", "#dd6b20", "#2f855a", "#c53030", "#4a5568", "#6b46c1"

def brad(word, M, q):
    """truncated Bradley q-MZV, sum over M >= n1 > ... > nr >= 1 (float)"""
    import itertools
    tot = 0.0
    for ns in itertools.combinations(range(M, 0, -1), len(word)):
        t = 1.0
        for s, n in zip(word, ns): t *= q**((s-1)*n) / sum(q**j for j in range(n))**s
        tot += t
    return tot

# fig4: holding one base fixed draws a tangent line to the border
fig, ax = plt.subplots(figsize=(9, 6.8))
K = np.linspace(-3, 3, 400)
ax.fill_between(K, K**2/4, 3, color=PURPLE, alpha=0.07)
ax.plot(K, K**2/4, color=PURPLE, lw=2.5, label="border k² = 4d: the two bases are equal")
for c in [-0.5, 0, 0.25, 0.75, 1.5, 2.0]:
    b = c - 1
    ax.plot(K, b*K - b*b, color=GRAY, lw=1, alpha=0.55)
    ax.scatter(2*b, b*b, color=GRAY, s=14, zorder=3)
ax.plot(K, 0*K, color=ORANGE, lw=3.5, label="hold q₂ = 1: d = 0 (Hoffman–Ihara, Bradley)")
ax.plot(K, -K/2 - 0.25, color=GREEN, lw=3.5, label="hold q₂ = 1/2: d = −k/2 − 1/4")
ax.scatter([0], [0], color=ORANGE, s=70, zorder=5); ax.scatter([-1], [0.25], color=GREEN, s=70, zorder=5)
ax.annotate("touches the border at k = 0", (0, 0), (0.45, -0.75), fontsize=9.5, color=ORANGE,
            arrowprops=dict(arrowstyle="->", color=ORANGE))
ax.annotate("touches the border at (−1, 1/4)", (-1, 0.25), (-2.9, 1.3), fontsize=9.5, color=GREEN,
            arrowprops=dict(arrowstyle="->", color=GREEN))
ax.text(1.05, -2.2, "gray: other held values q₂ = c\n(line d = (c−1)k − (c−1)², T6)", fontsize=9.5, color=GRAY)
ax.text(-2.95, -2.55, "With the step fixed at 1, the shadows are q₁ = 1 + α, q₂ = 1 + β,\n"
        "k = q₁ + q₂ − 2 and d = (q₁ − 1)(q₂ − 1)  (T4, Lean dial_from_bases)", fontsize=9.5)
ax.set_xlim(-3, 3); ax.set_ylim(-2.8, 3); ax.set_xlabel("k"); ax.set_ylabel("d"); ax.grid(alpha=0.3)
ax.legend(loc="upper right", fontsize=9)
ax.set_title("Fig 4. Holding one base fixed walks a straight line tangent to the border", fontsize=12)
fig.tight_layout(); fig.savefig("fig4-hold-one-base.png", dpi=150); plt.close(fig)

# fig5: what the jet q-MZV is, region by region, for the weight-2 value zeta(2)
fig, axs = plt.subplots(1, 3, figsize=(15, 4.9))
M = 40
ax = axs[0]
al = np.linspace(-0.45, 0.45, 60)
ax.plot(1 + al, [brad((2,), M, 1 + a) for a in al], color=BLUE, lw=2, label="first shadow: Bradley at q₁ = 1 + α (sums to n = 40)")
ax.axhline(brad((2,), M, 0.5), color=GREEN, lw=2.5, label=f"second shadow, held: Bradley at q₂ = 1/2 = {brad((2,), 400, 0.5):.4f}")
ax.axhline(brad((2,), 4000, 1.0), color=ORANGE, lw=2, ls="--", label="held at q₂ = 1 instead: ordinary ζ(2) (truncated)")
ax.set_xlabel("q₁"); ax.set_title("Off the border: two separate q-MZVs (T4, T5, T7)", fontsize=10.5)
ax.legend(fontsize=8, loc="upper left"); ax.grid(alpha=0.3)
ax = axs[1]
qq = np.linspace(0.55, 0.9, 50); h = 1e-5; M2 = 400
val = np.array([brad((2,), M2, x) for x in qq]); der = np.array([(brad((2,), M2, x + h) - brad((2,), M2, x - h))/(2*h) for x in qq])
ax.plot(qq, val, color=PURPLE, lw=2, label="value ζ_q(2)")
ax.plot(qq, der, color=RED, lw=2, ls="--", label="derivative dζ_q(2)/dq")
ax.set_xlabel("q = q₁ = q₂"); ax.set_title("On the border: value and q-derivative, a literal jet (T8)", fontsize=10.5)
ax.legend(fontsize=9); ax.grid(alpha=0.3)
ax = axs[2]
ax.axis("off")
ax.text(0.0, 0.92, "The letter product with q a jet 1 + N:", fontsize=11.5)
ax.text(0.04, 0.78, "zᵢ ◇ zⱼ = z₍ᵢ₊ⱼ₎ − N·z₍ᵢ₊ⱼ₋₁₎", fontsize=13, color=PURPLE)
ax.text(0.0, 0.62, "Three letters at once:", fontsize=11)
ax.text(0.04, 0.49, "z₁◇z₁◇z₁ = z₃ − 2N·z₂ + (kN − d)·z₁", fontsize=12.5, color=GREEN)
ax.text(0.0, 0.30, "d enters as a plain number at third order,\nlike the plane law's x²y term (T2b, Q3;\nLean jet_letter_triple)", fontsize=10.5)
ax.text(0.0, 0.06, "H(2)·H(2) = 2H(2,2) + H(4) − N·H(3)  (T3)", fontsize=11, color=BLUE)
fig.suptitle("Fig 5. The jet q-MZV: two q-MZVs off the border, value plus derivative on it", fontsize=12.5)
fig.tight_layout(); fig.savefig("fig5-jet-letters.png", dpi=150); plt.close(fig)
print("wrote fig4, fig5")

# fig6: the border at the origin: jet slope vs true slope, and the residue that makes the 1/2
import mpmath as mp
fig, axs = plt.subplots(1, 2, figsize=(13.5, 5))
ax = axs[0]
S = [2, 3, 4, 5]
jet = [(s/2 - 1)*(float(mp.zeta(s-1)) if s > 2 else 0) + (s/2)*float(mp.zeta(s)) for s in S]
true = [j + (0.5 if s == 2 else 0) for s, j in zip(S, jet)]   # measured in border_checks.py B4 (s = 2..5)
ax.bar(np.array(S) - 0.18, jet, width=0.36, color=PURPLE, label="jet slope: (s/2 − 1)ζ(s−1) + (s/2)ζ(s)")
ax.bar(np.array(S) + 0.18, true, width=0.36, color=GREEN, label="true slope of Bradley's ζ_q(s) at q → 1")
ax.annotate("gap exactly 1/2", (2.18, true[0]), (2.6, true[0] + 0.6), fontsize=10.5, color=RED,
            arrowprops=dict(arrowstyle="->", color=RED))
ax.set_xticks(S); ax.set_xlabel("s"); ax.set_ylabel("dζ_q(s)/dq at q = 1")
ax.set_title("At the origin the jet gives the slope, except at s = 2 (B4)", fontsize=11)
ax.legend(fontsize=9, loc="upper left"); ax.grid(alpha=0.3, axis="y")
ax = axs[1]
ss = np.concatenate([np.linspace(1.6, 1.995, 80), np.linspace(2.005, 2.6, 80)])
ax.plot(ss, [float((s/2 - 1)*mp.zeta(s - 1)) for s in ss], color=BLUE, lw=2.2, label="(s/2 − 1)·ζ(s − 1)")
ax.scatter([2], [0.5], s=80, color=RED, zorder=5)
ax.axhline(0.5, color=RED, lw=1, ls=":")
ax.text(1.62, 0.72, "→ 1/2 at s = 2: a zero times a pole (B5a)", color=RED, fontsize=10)
ax.text(1.62, 0.05, "at s = 2 the jet's lowering coefficient is 0,\nbut ζ(1) is infinite; together they leave 1/2.\n"
        "In depth 2 a front 2 leaves half of the rest:\nslope ζ_q(2,2) = 2ζ(2,2) + ζ(2)/2  (B5b)", fontsize=9.5)
ax.set_xlabel("s"); ax.set_ylim(-0.05, 1.0); ax.grid(alpha=0.3); ax.legend(fontsize=9, loc="upper right")
ax.set_title("Where the 1/2 comes from", fontsize=11)
fig.suptitle("Fig 6. The border at the origin (both bases 1): the jet is the q → 1 slope, plus a residue at a front 2", fontsize=12.5)
fig.tight_layout(); fig.savefig("fig6-border-slope.png", dpi=150); plt.close(fig)
print("wrote fig6")

# fig7: the proof's picture: the 1/2 is an area, and the rule matches the measured slopes
import json as _json
fig, axs = plt.subplots(1, 2, figsize=(13.5, 5.2))
ax = axs[0]
X = np.linspace(0.05, 8, 600)
a1 = 1/(4*np.sinh(X/2)**2); a2 = 1/X**2
ax.plot(X, a2, color=GRAY, lw=2, label="1/x²: the letter with q = 1")
ax.plot(X, a1, color=BLUE, lw=2, label="1/(4 sinh²(x/2)): Bradley's letter, rescaled")
ax.fill_between(X, a1, a2, color=RED, alpha=0.25, label="area between them = 1/2 exactly")
ax.set_ylim(0, 3); ax.set_xlim(0, 8); ax.set_xlabel("x = n·t  (t = −log q)"); ax.grid(alpha=0.3)
ax.text(2.2, 1.6, "∫₀^∞ [1/x² − 1/(4 sinh²(x/2))] dx = 1/2\n(antiderivative coth(x/2)/2 − 1/x;  P4b)", fontsize=10, color=RED)
ax.legend(fontsize=9, loc="upper right")
ax.set_title("A front 2 sees this area: it is the 1/2", fontsize=11)
ax = axs[1]
rows = _json.load(open("slope_proof_checks.json"))["slopes"]
m = np.array([r["measured"] for r in rows]); p = np.array([r["rule"] for r in rows])
cols = [GREEN if r["word"][0] == 2 else PURPLE for r in rows]
ax.scatter(p, m, c=cols, s=60, zorder=3)
lim = [0, max(m.max(), p.max())*1.08]; ax.plot(lim, lim, color=GRAY, lw=1)
placed = []
for r in rows:
    kk = sum(1 for v in placed if abs(v - r["rule"]) < 0.12); placed.append(r["rule"])
    ax.annotate("(" + ",".join(map(str, r["word"])) + ")", (r["rule"], r["measured"]), (8, -13 - 11*kk), textcoords="offset points", fontsize=8.5)
ax.set_xlabel("slope from the rule"); ax.set_ylabel("slope measured from the q-series at q → 1")
ax.set_title(f"11 words, depth 1 to 3: largest gap {np.abs(m - p).max():.1e}  (P6)", fontsize=11)
ax.scatter([], [], c=GREEN, label="starts with 2 (uses the 1/2 term)"); ax.scatter([], [], c=PURPLE, label="starts with 3 or more")
ax.legend(fontsize=9, loc="upper left"); ax.grid(alpha=0.3)
fig.suptitle("Fig 7. The slope rule, proved: the 1/2 is an area, and the rule lands on every measured slope", fontsize=12.5)
fig.tight_layout(); fig.savefig("fig7-slope-proof.png", dpi=150); plt.close(fig)
print("wrote fig7")
