"""DRAFT round 32 picture: f39_mirror.png, from results.json (panel A recomputes the mirror test on a finer grid of f)."""
import json, os
import numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch
HERE = os.path.dirname(os.path.abspath(__file__))
R = json.load(open(os.path.join(HERE, "results.json")))
GREY, BLUE, ORANGE, RED, GREEN = "#777777", "#1f6fc5", "#e08a1e", "#c8372d", "#2e9e4f"
fig = plt.figure(figsize=(20, 13.5))
gs = fig.add_gridspec(2, 2, hspace=0.28, wspace=0.16)
fig.suptitle("Round 32 (DRAFT): is there a reason inside the law to keep one hand of matter?", fontsize=16)

# A: mirror test across the family
def qmul(p, r):
    a, b, c, d = p; e, f_, g, k = r
    return np.array([a*e - b*f_ - c*g - d*k, a*f_ + b*e + c*k - d*g, a*g - b*k + c*e + d*f_, a*k + b*g - c*f_ + d*e])
def qpow(p, t):
    ang = np.arccos(np.clip(p[0], -1, 1)); v = p[1:]; nv = np.linalg.norm(v)
    return np.concatenate([[np.cos(t * ang)], np.sin(t * ang) * v / nv])
def qc(p): return p * np.array([1, -1, -1, -1])
def member(x, y, f): return qmul(qmul(qpow(x, 1 - f / 2), y), qpow(x, f / 2))
rng = np.random.default_rng(32); pairs = []
for _ in range(60):
    x = rng.normal(size=4); y = rng.normal(size=4); pairs.append((x / np.linalg.norm(x), y / np.linalg.norm(y)))
fg = np.linspace(0, 2, 201)
err = [max(np.abs(qc(member(x, y, f)) - member(qc(x), qc(y), f)).max() for x, y in pairs) for f in fg]
ax = fig.add_subplot(gs[0, 0])
ax.semilogy(fg, np.maximum(err, 1e-16), color=GREY, lw=2)
ax.plot([1], [max(err[100], 1e-16)], "o", ms=12, color=GREEN)
ax.set_xlabel("member f   (x ⊕_f y = x^(1 − f/2) · y · x^(f/2))"); ax.set_ylabel("how far the mirror changes the rule")
ax.set_title("A. The law's sum reads the same backwards: only f = 1 is its own mirror (K1)")
ax.annotate("the law, x^½ · y · x^½\nits own mirror (error 10⁻¹⁶)", (1, 3e-16), (1.15, 1e-10), color=GREEN, fontsize=11,
            arrowprops=dict(arrowstyle="->", color=GREEN))
ax.annotate("", (1.9, 8.0), (0.1, 8.0), arrowprops=dict(arrowstyle="<->", color=ORANGE, lw=2, connectionstyle="arc3,rad=-0.12"))
ax.text(0.45, 1e-5, "the mirror swaps f and 2 − f:\nxy (f = 0) ↔ yx (f = 2)", ha="center", color=ORANGE, fontsize=11)
ax.set_ylim(1e-17, 1e5)

# B: Dirac spectrum vs f
ax = fig.add_subplot(gs[0, 1])
ff = np.linspace(-0.3, 2.3, 300)
for m in range(4):
    for sg in (1, -1):
        lam = sg * (m + 1.5) + 1.5 * (1 - ff)
        ax.plot(ff, lam, color=BLUE, lw=1 + 0.6 * m, alpha=0.9)
    ax.text(2.32, (m + 1.5) - 1.5 * 1.3, f"{(m+1)*(m+2)}", color=BLUE, fontsize=9, va="center")
ax.axhline(0, color=GREY, lw=1); ax.axvline(1, color=GREEN, lw=2, ls="--")
ax.plot([0, 2], [0, 0], "o", ms=12, color=ORANGE)
ax.text(1.03, 5.6, "the law: ±3/2, ±5/2, …\nsymmetric, nothing massless", color=GREEN, fontsize=11)
ax.text(0.05, -5.6, "f = 0 and f = 2: 2 massless modes each\n(parallel spinors of the flat twists)", color=ORANGE, fontsize=11, bbox=dict(fc="white", ec="none", alpha=0.9))
ax.set_ylim(-6.5, 7); ax.set_xlim(-0.3, 2.5)
ax.set_xlabel("member f"); ax.set_ylabel("Dirac eigenvalue on the hidden 3-sphere (unit radius)")
ax.set_title("B. The whole spectrum slides by 3(1 − f)/2; line width and right-hand numbers = how many modes (K2)")

# C: eta invariant
ax = fig.add_subplot(gs[1, 0])
fe = np.linspace(0.001, 1.999, 400)
eta = 0.75 * (fe - 1) * (3 * fe**2 - 6 * fe + 2)
ax.plot(fe, eta, color=ORANGE, lw=2.5)
ax.axhline(0, color=GREY, lw=1)
z = [1 - 1 / np.sqrt(3), 1 + 1 / np.sqrt(3)]
ax.plot([1], [0], "o", ms=12, color=GREEN); ax.plot(z, [0, 0], "o", ms=9, mfc="white", mec=ORANGE, mew=2)
ax.text(1.03, 0.08, "the law: 0", color=GREEN, fontsize=11)
ax.text(z[0] - 0.05, 0.3, "accidental zero\n1 − 1/√3", ha="center", color=ORANGE, fontsize=10, bbox=dict(fc="white", ec="none", alpha=0.9))
ax.text(z[1] + 0.05, -0.45, "accidental zero\n1 + 1/√3", ha="center", color=ORANGE, fontsize=10, bbox=dict(fc="white", ec="none", alpha=0.9))
ax.text(1.0, -1.35, "η(0) = 3(f − 1)(3f² − 6f + 2)/4\nodd under f → 2 − f: a parity-odd number,\nnot massless matter of one hand",
        color=ORANGE, fontsize=11)
ax.set_xlabel("member f"); ax.set_ylabel("spectral asymmetry η(0)")
ax.set_title("C. The lopsided members break the mirror, as a number (K3)")

# D: chain
ax = fig.add_subplot(gs[1, 1]); ax.axis("off"); ax.set_xlim(0, 10); ax.set_ylim(0, 10)
ax.set_title("D. Why the law can't keep one hand")
boxes = [(5, 9.0, "the law's sum x^½ · y · x^½ reads the same backwards:\nno preferred hand (K1)", GREEN),
         (5, 6.95, "its Dirac spectrum on the hidden 3-sphere is exactly symmetric;\nthe other members slide it and break the mirror only as a number (K2, K3)", BLUE),
         (5, 4.9, "our 4 directions + 3 hidden = 7: no one-hand condition in odd dimensions;\nevery hidden massless mode is a full left + right pair, for every f (K4)", RED),
         (5, 2.85, "one hand needs a choice from outside the law:\na one-hand rule in 4 + 2 (round 31's 2-sphere), or one sign of weight", ORANGE),
         (5, 0.85, "assumption (b) stays; no lead", GREY)]
for xx, yy, txt, col in boxes:
    ax.add_patch(FancyBboxPatch((xx - 4.7, yy - 0.75), 9.4, 1.5, boxstyle="round,pad=0.08", fc="white", ec=col, lw=2))
    ax.text(xx, yy, txt, ha="center", va="center", fontsize=11, color=col)
for y0, y1 in [(8.25, 7.7), (6.2, 5.65), (4.15, 3.6), (2.1, 1.6)]:
    ax.annotate("", (5, y1), (5, y0), arrowprops=dict(arrowstyle="->", color=GREY, lw=2))
os.makedirs(os.path.join(HERE, "figures"), exist_ok=True)
fig.savefig(os.path.join(HERE, "figures", "f39_mirror.png"), dpi=120, bbox_inches="tight")
print("saved")
