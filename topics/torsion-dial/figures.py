"""DRAFT figures for torsion-dial. Run after checks.py (reads results.json). Writes figures/t1-t3."""
import json, numpy as np, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch
from scipy.linalg import expm, logm
R = json.load(open("results.json"))
BLUE, GREEN, ORANGE, RED, GREY, PURPLE = "#2b6cb0", "#2f855a", "#dd6b20", "#c53030", "#718096", "#6b46c1"
plt.rcParams.update({"font.size": 10})

# ---------------------------------------------------------------- t1: the torsion dial f
fig, ax = plt.subplots(1, 2, figsize=(12, 4.6))
f = np.linspace(-0.6, 2.6, 400)
a = ax[0]
a.axhline(0, color=GREY, lw=0.8)
a.plot(f, 2 * (f - 1), "--", color=GREY, lw=2, label="torsion 2(f − 1)")
a.plot(f, 2 * f - f**2, color=BLUE, lw=2.5, label="curvature 2f − f²")
for x0, lab, col in [(0, "f = 0\nleft product xy\n(the law's chains)", ORANGE), (1, "f = 1\nthe law's sphere sum\nx^½ y x^½", GREEN),
                     (2, "f = 2\nright product yx", ORANGE)]:
    a.plot([x0], [2 * x0 - x0**2], "o", color=col, ms=10, zorder=5)
    a.plot([x0], [2 * (x0 - 1)], "s", color=col, ms=7, zorder=5)
    a.annotate(lab, (x0, -2.25), ha="center", va="top", fontsize=8.5, color=col)
a.set_ylim(-3.6, 2.6); a.set_xlabel("f  (how a sideways direction is carried)")
a.set_title("A. One family, one number f: torsion is a line, curvature an arch")
a.legend(loc="upper left", fontsize=8.5)
top = a.secondary_xaxis("top", functions=(lambda x: 1 - x / 2, lambda l: 2 - 2 * l))
top.set_xlabel("order-dial page's λ = 1 − f/2  (mutation a + b + λab + (1 − λ)ba)")
b = ax[1]
fg = np.array(R["T5"]["f"]); ai = np.array(R["T5"]["auto_inverse_err"]); asx = np.array(R["T5"]["assoc_err"])
b.semilogy(fg, np.maximum(ai, 1e-16), color=GREEN, lw=2, label="breaks the law's rule −(a ⊕ b) = (−a) ⊕ (−b)")
b.semilogy(fg, np.maximum(asx, 1e-16), color=ORANGE, lw=2, label="breaks grouping (x ⊕ y) ⊕ z = x ⊕ (y ⊕ z)")
b.axvline(1, color=GREEN, ls=":", lw=1); b.axvline(0, color=ORANGE, ls=":", lw=1); b.axvline(2, color=ORANGE, ls=":", lw=1)
b.set_xlabel("f"); b.set_ylabel("largest error over random x, y, z")
b.set_title("B. The law's own rule picks f = 1; grouping picks 0 or 2")
b.legend(loc="lower right", fontsize=8.5); b.set_ylim(1e-17, 10)
fig.tight_layout(); fig.savefig("figures/t1_torsion_dial.png", dpi=150); plt.close(fig)

# ---------------------------------------------------------------- t2: where x is split around y
PA = [np.array([[0, 1], [1, 0]], complex), np.array([[0, -1j], [1j, 0]]), np.array([[1, 0], [0, -1]], complex)]
E = [-1j * m for m in PA]
mat = lambda v: sum(v[i] * E[i] for i in range(3))
def quat(M):  # back to (w, x, y, z)
    w = np.real(np.trace(M)) / 2; v = [np.real(np.trace(M @ np.conj(e).T)) / 2 for e in E]; return np.array([w] + v)
mp = lambda M, s: expm(s * logm(M))
def stereo(p): return p[1:] / (1 + p[0])
x, y = expm(mat(np.array([1.1, 0.2, 0.0]))), expm(mat(np.array([0.1, 1.0, 0.3])))
fs = np.linspace(-0.5, 2.5, 121); P3 = np.array([stereo(quat(mp(x, 1 - ff / 2) @ y @ mp(x, ff / 2))) for ff in fs])
cen = P3.mean(0); _, _, Vt = np.linalg.svd(P3 - cen); proj = lambda p: (p - cen) @ Vt[:2].T   # the plane the path lies in
pts = proj(P3)
fig = plt.figure(figsize=(12, 4.6))
a = fig.add_subplot(1, 2, 1); a.set_xlim(0, 10); a.set_ylim(0, 6); a.axis("off")
a.set_title("A. Every member is 'part of x, then y, then the rest of x'.\nf says how much of x comes after y")
for row, (ff, lab, col) in enumerate([(0, "f = 0:  x · y   (chains; torsion −2)", ORANGE), (1, "f = 1:  x^½ · y · x^½   (the law; torsion 0)", GREEN),
                                      (2, "f = 2:  y · x   (torsion +2)", ORANGE)]):
    yy = 4.6 - row * 1.7; w1 = 4 * (1 - ff / 2); w2 = 4 * ff / 2
    if w1 > 0: a.add_patch(FancyBboxPatch((0.5, yy), w1, 0.6, boxstyle="round,pad=0.02", color=BLUE, alpha=0.75))
    a.add_patch(FancyBboxPatch((0.5 + w1 + 0.1, yy), 1.6, 0.6, boxstyle="round,pad=0.02", color=PURPLE, alpha=0.75))
    if w2 > 0: a.add_patch(FancyBboxPatch((0.5 + w1 + 1.8, yy), w2, 0.6, boxstyle="round,pad=0.02", color=BLUE, alpha=0.75))
    a.text(0.5, yy + 0.85, lab, color=col, fontsize=9.5)
    if w1 > 0: a.text(0.5 + w1 / 2, yy + 0.3, f"x^{1 - ff/2:g}", ha="center", va="center", color="white", fontsize=9)
    a.text(0.5 + w1 + 0.9, yy + 0.3, "y", ha="center", va="center", color="white", fontsize=10)
    if w2 > 0: a.text(0.5 + w1 + 1.8 + w2 / 2, yy + 0.3, f"x^{ff/2:g}", ha="center", va="center", color="white", fontsize=9)
a.text(0.5, 0.1, "Torsion = 2(f − 1): how lopsided the split is.  If x and y lie on one line, every split gives\nthe same answer, the scalar law a ∘ₖ b. That is why one number can never show torsion.", fontsize=8.5, color=GREY)
b = fig.add_subplot(1, 2, 2)
b.plot(pts[:, 0], pts[:, 1], color=GREY, lw=1.5)
sc = b.scatter(pts[:, 0], pts[:, 1], c=fs, cmap="coolwarm", s=12, zorder=3)
for ff, lab, col in [(0, "f = 0: xy", ORANGE), (1, "f = 1: the law's x ⊕ y", GREEN), (2, "f = 2: yx", ORANGE)]:
    p = proj(stereo(quat(mp(x, 1 - ff / 2) @ y @ mp(x, ff / 2)))); b.plot(p[0], p[1], "o", color=col, ms=11, zorder=4)
    b.annotate(lab, p[:2], xytext=(8, 8), textcoords="offset points", color=col, fontsize=9)
plt.colorbar(sc, ax=b, label="f")
b.set_xlabel("law's chart, seen face-on to the path"); b.set_ylabel("")
b.set_title("B. One pair x, y: turning f swings the answer round an arc from xy to yx.\nThe law sits exactly halfway (checks T2, T3)")
b.set_aspect("equal", adjustable="datalim")
fig.tight_layout(); fig.savefig("figures/t2_split.png", dpi=150); plt.close(fig)

# ---------------------------------------------------------------- t3: where torsion lives, and how big it could be
fig, ax = plt.subplots(1, 2, figsize=(12, 4.6), gridspec_kw={"width_ratios": [1.25, 1]})
a = ax[0]; a.set_xlim(0, 10); a.set_ylim(0, 7); a.axis("off")
a.set_title("A. Where the law carries torsion, and what sets its size")
boxes = [(0.3, 4.4, "One number: a ∘ₖ b,\nalso (a + b + k·ab)/(1 − d·ab)", "torsion 0\n(inputs commute)", GREY),
         (0.3, 2.3, "Chains: A ∘ₖ B = A + B + k·AB\n(lifts I + kA multiply; f = 0)", "torsion = k × bracket\nflat (associative)\nκ-Poincaré's momentum group", ORANGE),
         (0.3, 0.2, "The sum on the 3-sphere\n(round 18; x^½ y x^½; f = 1)", "torsion 0, forced\ncurvature q", GREEN)]
for x0, y0, t1, t2, col in boxes:
    a.add_patch(FancyBboxPatch((x0, y0), 4.6, 1.6, boxstyle="round,pad=0.08", fc=col, alpha=0.18, ec=col))
    a.text(x0 + 0.15, y0 + 0.8, t1, va="center", fontsize=9)
    a.add_patch(FancyBboxPatch((5.4, y0), 4.3, 1.6, boxstyle="round,pad=0.08", fc="white", ec=col, lw=1.5))
    a.text(5.55, y0 + 0.8, t2, va="center", fontsize=9, color=col)
    a.annotate("", (5.35, y0 + 0.8), (4.95, y0 + 0.8), arrowprops=dict(arrowstyle="->", color=col))
b = ax[1]
labs = ["lab bounds\n(Kostelecký et al. 2008)", "Einstein–Cartan,\nall nuclear spins aligned", "flat members f = 0, 2\non today's cosmic slices", "the law (f = 1)"]
vals = [R["T8"]["bound_GeV"], R["T8"]["EC_polarized_nuclear_GeV"], R["T8"]["flat_members_cosmic_GeV"], None]
cols = [RED, PURPLE, ORANGE, GREEN]
for i, (v, c) in enumerate(zip(vals, cols)):
    if v: b.barh(i, np.log10(v) + 50, left=-50, color=c, alpha=0.8); b.text(np.log10(v) + 0.5, i, f"{v:.1e} GeV", va="center", fontsize=9)
    else: b.text(-49, i, "exactly zero: no bar", va="center", fontsize=9, color=c)
b.set_yticks(range(4)); b.set_yticklabels(labs, fontsize=8.5); b.set_ylim(3.6, -0.6)
b.set_xlim(-50, -25); b.set_xlabel("log₁₀ (torsion as an energy, GeV)")
b.set_title("B. Sizes, order of magnitude (check T8)")
fig.tight_layout(); fig.savefig("figures/t3_where_torsion.png", dpi=150); plt.close(fig)
print("figures written")
