import numpy as np, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch, Circle, FancyBboxPatch
out = "figures/"
plt.rcParams.update({"font.size": 11})
DRAFT = dict(fontsize=9, color="0.5", ha="right", va="bottom")

# ---------- Fig 1: the regions of the odds dial ----------
fig, ax = plt.subplots(figsize=(9, 5.2))
for lo, hi, c in [(-4, -1, "#f6d5d5"), (-1, 0, "#dde7f7"), (0, 4, "#dcefdc")]:
    ax.axvspan(lo, hi, color=c, zorder=0)
t1 = np.linspace(-4, -1.02, 300); t2 = np.linspace(-0.98, 4, 400)
ax.plot(t1, t1/(1+t1), "k", lw=2); ax.plot(t2, t2/(1+t2), "k", lw=2)
ax.axvline(-1, color="crimson", lw=2.5, ls="--")
ax.axhline(0, color="0.4", lw=0.8); ax.axhline(1, color="0.4", lw=0.8, ls=":")
ax.set_ylim(-4, 5); ax.set_xlim(-4, 4)
ax.set_xlabel("odds t (chance over not-chance)"); ax.set_ylabel("chance p = t/(1+t)")
ax.text(2.4, -1.4, "real chances\n0 < p < 1", ha="center", fontsize=11)
ax.text(-0.42, -3.4, "rewinds\n(undo a chance)\np < 0", ha="center", fontsize=10)
ax.text(-2.5, 3.7, "past the wall\np > 1\n(sign flipped)", ha="center", fontsize=10)
ax.annotate("the wall t = −1\np → ±∞\nswallows everything", xy=(-1, -1.5), xytext=(-3.6, -2.6),
            color="crimson", arrowprops=dict(arrowstyle="->", color="crimson"))
for t, lab in [(0, "never (t=0)"), (1, "coin flip (t=1)"), (-2, "200% (t=−2)")]:
    ax.plot(t, t/(1+t), "o", color="navy"); ax.annotate(lab, (t, t/(1+t)), xytext=(6, 6), textcoords="offset points", color="navy")
ax.annotate("certain\n(t → ∞)", xy=(4, 0.8), xytext=(3.0, 2.2), arrowprops=dict(arrowstyle="->"))
ax.set_title("Fig 1. What each part of the odds dial means")
fig.text(0.99, 0.01, "DRAFT", **DRAFT); fig.tight_layout(); fig.savefig(out+"fig1_regions.png", dpi=140); plt.close()

# ---------- Fig 2: the triangle of special points and three reversals ----------
fig, ax = plt.subplots(figsize=(8.5, 7)); ax.set_aspect("equal"); ax.axis("off")
V = {"never": (0, 0), "certain": (4, 0), "wall": (2, 3.46)}
labs = {"never": "NEVER\nt = 0, p = 0", "certain": "CERTAIN\nt = ∞, p = 1", "wall": "THE WALL\nt = −1, p = ∞"}
cols = {"never": "#3a7d44", "certain": "#2b5fa8", "wall": "crimson"}
for k, (x, y) in V.items():
    ax.add_patch(Circle((x, y), 0.42, color=cols[k], alpha=0.9, zorder=3))
    off = {"never": (-0.2, -0.95), "certain": (0.2, -0.95), "wall": (0, 0.6)}[k]
    ax.text(x+off[0], y+off[1], labs[k], ha="center", va="center", fontsize=11, color=cols[k], weight="bold")
def side(a, b, txt, fixed, dx=0, dy=0):
    (x1, y1), (x2, y2) = V[a], V[b]
    ax.add_patch(FancyArrowPatch((x1, y1), (x2, y2), arrowstyle="<->", mutation_scale=20, lw=2, color="0.25",
                                 shrinkA=22, shrinkB=22, zorder=2))
    ax.text((x1+x2)/2+dx, (y1+y2)/2+dy, txt + f"\n(fixes {fixed})", ha="center", va="center", fontsize=10,
            bbox=dict(boxstyle="round", fc="white", ec="0.6"))
side("never", "certain", "complement  t ↦ 1/t,  p ↦ 1 − p\nthermodynamic time reversal\n(Crooks; Feynman–Stueckelberg E ↦ −E)", "the wall and the coin flip", 0, -0.05)
side("certain", "wall", "rewind  t ↦ −t/(1+t)\nclock run backward\n(undo; p ↦ −t)", "never and 200%", 1.15, 0.25)
side("never", "wall", "swap  t ↦ −1 − t,  p ↦ 1/p", "certain", -1.15, 0.25)
ax.text(2, 1.25, "rewind after complement\n= a turn of order 3\n(never → wall → certain)", ha="center", fontsize=9.5, style="italic", color="0.3")
ax.set_xlim(-1.6, 5.6); ax.set_ylim(-1.6, 4.6)
ax.set_title("Fig 2. Three special points, three mirrors: one six-map symmetry\n(the classical cross-ratio group of 0, 1, ∞ written in odds)")
fig.text(0.99, 0.01, "DRAFT", **DRAFT); fig.savefig(out+"fig2_triangle.png", dpi=140, bbox_inches="tight"); plt.close()

# ---------- Fig 3: the clock ----------
fig, (ax, bx) = plt.subplots(1, 2, figsize=(11, 4.6), gridspec_kw=dict(width_ratios=[1.6, 1]))
r = 1.0; T = np.linspace(-4, 2.2, 400)
ax.plot(T, np.expm1(r*T), "k", lw=2)
ax.axhline(-1, color="crimson", ls="--", lw=2); ax.axhline(0, color="0.5", lw=0.8); ax.axvline(0, color="0.5", lw=0.8)
ax.fill_between(T, np.expm1(r*T), 0, where=T < 0, color="#dde7f7"); ax.fill_between(T, np.expm1(r*T), 0, where=T > 0, color="#dcefdc")
ax.text(-3.9, -0.88, "wall = the infinite past", color="crimson")
ax.text(-3.9, 0.7, "negative time → odds in (−1, 0)\n(rewinds)", fontsize=10)
ax.text(-3.9, 6.3, "forward time → real chances\n(certain as T → ∞)", fontsize=10)
T1, T2 = 0.6, 0.9
ax.annotate("", xy=(T1+T2, np.expm1(T1+T2)), xytext=(T1, np.expm1(T1)), arrowprops=dict(arrowstyle="->", color="navy"))
ax.text(1.6, 2.6, "stacking items\n= adding time", color="navy", fontsize=10)
ax.set_xlabel("elapsed time T on a steady random clock (rate r = 1)"); ax.set_ylabel("odds t = e^{rT} − 1")
ax.set_ylim(-1.3, 8); ax.set_title("Fig 3a. Odds are a clock reading")
# complex time: t < -1 needs Im T = pi/r
bx.axhline(0, color="0.5", lw=0.8); bx.axhline(np.pi, color="crimson", lw=2)
bx.plot([-3, 2], [0, 0], color="#3a7d44", lw=4, alpha=0.6)
bx.text(-2.9, 0.25, "real time: odds in (−1, ∞)", color="#3a7d44")
bx.text(-2.9, np.pi+0.25, "half a turn of imaginary time:\nodds below −1 (sign flipped)", color="crimson")
bx.annotate("", xy=(0, np.pi), xytext=(0, 0), arrowprops=dict(arrowstyle="->", lw=1.5))
bx.text(0.1, 1.4, "+ iπ/r\nmultiplies 1+t by −1", fontsize=10)
bx.set_xlim(-3, 2); bx.set_ylim(-0.8, 4.6); bx.set_xlabel("Re T"); bx.set_ylabel("Im T")
bx.set_title("Fig 3b. No real time reaches t < −1")
fig.text(0.99, 0.01, "DRAFT", **DRAFT); fig.tight_layout(); fig.savefig(out+"fig3_clock.png", dpi=140); plt.close()

# ---------- Fig 4: log-odds plane: two lines, time reversal, fermion points ----------
fig, ax = plt.subplots(figsize=(9, 5))
ax.axhline(0, color="#3a7d44", lw=3); ax.axhline(np.pi, color="crimson", lw=3); ax.axhline(-np.pi, color="crimson", lw=3, alpha=0.35)
ax.text(1.6, 0.2, "positive odds (real chances)", color="#3a7d44")
ax.text(1.6, np.pi+0.2, "negative odds", color="crimson"); ax.text(1.6, -np.pi+0.2, "same line again (one full turn = 2π)", color="crimson", alpha=0.6)
ax.plot(0, 0, "o", color="navy", ms=9); ax.text(0.1, -0.5, "coin flip t = 1", color="navy")
ax.plot(0, np.pi, "o", color="crimson", ms=11); ax.text(0.1, np.pi-0.55, "the wall t = −1  (log = iπ)", color="crimson")
x0 = 1.2
ax.annotate("", xy=(-x0, 0.05), xytext=(x0, 0.05), arrowprops=dict(arrowstyle="<->", color="0.2"))
ax.annotate("", xy=(-x0, np.pi+0.05), xytext=(x0, np.pi+0.05), arrowprops=dict(arrowstyle="<->", color="0.2"))
ax.text(-3.9, 1.4, "time reversal λ ↦ −λ\n(t ↦ 1/t) flips each line\nend to end; its only fixed points\nare the coin flip and the wall", fontsize=10)
ax.set_xlabel("Re λ = log |t|   (with t = e^{−βE}: minus energy over temperature)"); ax.set_ylabel("Im λ")
ax.set_yticks([-np.pi, 0, np.pi]); ax.set_yticklabels(["−π", "0", "π"]); ax.set_xlim(-4, 4); ax.set_ylim(-4.2, 4.6)
ax.set_title("Fig 4. Log-odds plane: negative odds sit half a turn up, not mirrored in time.\nFor one fermion mode (t = e^{−βE}) the wall is hit at the Matsubara points E = i(2n+1)π/β")
fig.text(0.99, 0.01, "DRAFT", **DRAFT); fig.tight_layout(); fig.savefig(out+"fig4_log_odds.png", dpi=140); plt.close()

# ---------- Fig 5: game tokens ----------
fig, ax = plt.subplots(figsize=(10, 4.6)); ax.axis("off"); ax.set_xlim(0, 10); ax.set_ylim(0, 4.6)
toks = [("×4/3", "25% charm", "#dcefdc"), ("×2", "50% charm", "#dcefdc"), ("×3/4", "curse\n(−1/3 'chance')", "#dde7f7"),
        ("DOOM", "×0, kept as\none zero", "#f6d5d5"), ("DESTINY", "×∞, kept as\none anti-zero", "#fff2cc")]
for i, (a, b, c) in enumerate(toks):
    x = 0.4 + i*1.95
    ax.add_patch(FancyBboxPatch((x, 2.6), 1.6, 1.5, boxstyle="round,pad=0.05", fc=c, ec="0.3"))
    ax.text(x+0.8, 3.6, a, ha="center", fontsize=14, weight="bold"); ax.text(x+0.8, 2.95, b, ha="center", fontsize=9.5)
ax.text(0.4, 2.0, "Rule: multiply the tokens' multipliers m = 1 + t, then chance = 1 − 1/m.  Order never matters.", fontsize=11)
ax.text(0.4, 1.45, "Charms only: 4/3 × 2 = 8/3 → chance 1 − 3/8 = 5/8.     Add the curse: × 3/4 = 2 → chance 1/2.", fontsize=10.5)
ax.text(0.4, 0.9, "Doom and Destiny: count them.  More Doom → impossible.  More Destiny → certain.", fontsize=10.5)
ax.text(0.4, 0.4, "Equal → they cancel (0 × ∞ answered by keeping the count) and the charms decide.", fontsize=10.5)
ax.set_title("Fig 5. A token mechanic built from the odds dial (checked in checks.py)")
fig.text(0.99, 0.01, "DRAFT", **DRAFT); fig.savefig(out+"fig5_tokens.png", dpi=140, bbox_inches="tight"); plt.close()
print("figures done")
