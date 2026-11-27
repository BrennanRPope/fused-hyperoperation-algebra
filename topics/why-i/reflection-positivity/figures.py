"""Pictures for topics/why-i/reflection-positivity (Nov 27, 2026). Run: python3 figures.py"""
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from rp import slice_multiplier, gram, rp_violation  # noqa: E402

OUT = os.path.join(HERE, "figures")
os.makedirs(OUT, exist_ok=True)
BLUE, ORANGE, GREEN, RED, GREY, PURPLE = "#2b6cb0", "#dd6b20", "#2f855a", "#c53030", "#718096", "#6b46c1"
plt.rcParams.update({"font.size": 10, "axes.titlesize": 11, "figure.dpi": 130})
rng = np.random.default_rng(20261127)


def save(fig, name):
    fig.tight_layout()
    fig.savefig(os.path.join(OUT, name + ".png"), bbox_inches="tight")
    plt.close(fig)


def arrow(ax, x0, y0, x1, y1, col, lw=2, style="-|>", rad=0.0):
    ax.add_patch(FancyArrowPatch((x0, y0), (x1, y1), arrowstyle=style, mutation_scale=14,
                                 color=col, lw=lw, connectionstyle=f"arc3,rad={rad}"))


Q, _ = np.linalg.qr(rng.normal(size=(3, 3)))
EIGS = [-1.0, 0.5, 2.0]
H = Q @ np.diag(EIGS) @ Q.T
eps = 0.25

# ---------------------------------------------------------------- f1 the dial plane
xs = np.linspace(-4, 6, 81)
ys = np.linspace(-2.5, 2.5, 41)
V = np.zeros((len(ys), len(xs)))
for i, y in enumerate(ys):
    for j, x in enumerate(xs):
        V[i, j] = rp_violation(gram(slice_multiplier(x + 1j * y, eps, H), 3))
fig, ax = plt.subplots(figsize=(9, 5.2))
im = ax.imshow(np.log10(V + 1e-6), origin="lower", extent=[xs[0], xs[-1], ys[0], ys[-1]],
               aspect="auto", cmap="Greys", vmin=-6, vmax=0.5)
cb = fig.colorbar(im, ax=ax, shrink=0.8)
cb.set_label("how badly <θF, F> ≥ 0 fails (log scale; white = holds)")
ax.plot([-2, 4], [0, 0], color=GREEN, lw=6, solid_capstyle="butt", label="reflection positive (link): 1 + kεE ≥ 0 for every E")
ax.plot([xs[0], -2], [0, 0], color=GREEN, lw=1.5, ls=":")
ax.plot([4, xs[-1]], [0, 0], color=GREEN, lw=1.5, ls=":", label="site reflection only (slice has a negative eigenvalue)")
ax.axvline(0, color=PURPLE, lw=1.5, ls="--")
ax.text(0.1, 2.15, "imaginary axis: size kept,\nFeynman's dial k = −i/ħ", color=PURPLE, fontsize=9)
ax.plot([-1], [0], "o", color=BLUE, ms=8)
ax.text(-1.05, -0.45, "heat dial\nk = −1/ħ", color=BLUE, ha="center", fontsize=9)
ax.plot([0], [-1], "o", color=PURPLE, ms=8)
arrow(ax, -1, 0.05, -0.05, -0.95, ORANGE, lw=2.2, rad=0.45)
ax.text(-1.95, -1.25, "Osterwalder–Schrader:\npositivity here becomes\nunitarity there", color=ORANGE, fontsize=9)
ax.set_xlabel("Re k (dial, real part)")
ax.set_ylabel("Im k")
ax.set_title("1. Reflection positivity lives exactly on the real dial axis (three states, ε = 0.25, energies −1, 0.5, 2)")
ax.legend(loc="lower right", fontsize=8)
save(fig, "f1-the-dial-plane")

# ---------------------------------------------------------------- f2 the real axis sweep
ks = np.linspace(-3.5, 5.5, 181)
gmin, tmin = [], []
for k in ks:
    T = slice_multiplier(k, eps, H)
    G = gram(T, 5).real
    gmin.append(np.linalg.eigvalsh((G + G.T) / 2).min() / np.abs(G).max())
    tmin.append(np.linalg.eigvalsh(T).min())
fig, ax = plt.subplots(1, 2, figsize=(12, 4.3))
a = ax[0]
a.axhline(0, color=GREY, lw=1)
a.axvspan(-2, 4, color=GREEN, alpha=0.12)
a.plot(ks, tmin, color=BLUE, lw=2, label="lowest eigenvalue of the slice 1 + kεH")
a.set_xlabel("dial k (real)")
a.set_title("2a. the slice multiplier on each energy: 1 + kεE")
for E, c in zip(EIGS, [RED, GREY, ORANGE]):
    a.plot(ks, 1 + ks * eps * E, color=c, lw=1, ls="--", label=f"E = {E:g}")
a.text(1, 2.3, "every slice ≥ 0", color=GREEN, ha="center")
a.legend(fontsize=8, loc="upper left")
a.set_ylim(-1.5, 3)
a = ax[1]
a.axhline(0, color=GREY, lw=1)
a.axvspan(-2, 4, color=GREEN, alpha=0.12)
a.plot(ks, gmin, color=PURPLE, lw=2)
a.set_xlabel("dial k (real)")
a.set_title("2b. lowest eigenvalue of the brute-force Gram matrix <θFᵢ, Fⱼ>")
a.text(1, -0.3, "never negative:\nreflection positive", color=GREEN, ha="center")
a.text(-3, -0.05, "negative:\nfails", color=RED, ha="center", va="top")
a.text(5, -0.05, "negative:\nfails", color=RED, ha="center", va="top")
a.annotate("walls at k = −1/(ε·E_max) = −2\nand k = 1/(ε·|E_min|) = 4", xy=(4, 0), xytext=(0.2, -0.12),
           fontsize=8, arrowprops=dict(arrowstyle="->", color=GREY))
save(fig, "f2-the-real-axis-window")

# ---------------------------------------------------------------- f3 why: a mirror makes a square
fig, ax = plt.subplots(1, 2, figsize=(12, 4.2))
for a, kind in zip(ax, ["site", "link"]):
    a.axis("off")
    a.set_xlim(-0.5, 8.5)
    a.set_ylim(-2.2, 2.2)
    N = 6 if kind == "site" else 7
    xs_ = np.arange(N + 1) * (8 / N)
    path = [0.3, -0.6, 0.8, 0.2, 0.8, -0.6, 0.3, 0.9][:N + 1]
    if kind == "link":
        path = [0.3, -0.6, 0.8, 0.1, 0.1, 0.8, -0.6, 0.3]
    mid = 4.0
    a.axvline(mid, color=ORANGE, lw=2, ls="--")
    a.plot(xs_, path, "-o", color=GREY, lw=1.5)
    for x, y in zip(xs_, path):
        a.plot(x, y, "o", color=BLUE if x > mid + 1e-9 else (PURPLE if x < mid - 1e-9 else ORANGE), ms=9)
    a.add_patch(FancyBboxPatch((mid + 0.15, -1.25), 3.9, 2.4, boxstyle="round,pad=0.05", fill=False, ec=BLUE, lw=1.5))
    a.add_patch(FancyBboxPatch((-0.2, -1.25), 3.9, 2.4, boxstyle="round,pad=0.05", fill=False, ec=PURPLE, lw=1.5))
    a.text(6.1, 1.45, "F looks at the later half", color=BLUE, ha="center")
    a.text(1.75, 1.45, "θF: the same F, mirrored", color=PURPLE, ha="center")
    a.text(mid, -1.75, "mirror", color=ORANGE, ha="center")
    if kind == "site":
        a.set_title("3a. site mirror: both halves meet at one site, and if the slice is\n"
                    "symmetric they are equal, so <θF, F> = Σ (same number)² ≥ 0")
    else:
        a.text(mid, 0.75, "one slice T\nacross the mirror", ha="center", fontsize=8, color=ORANGE)
        a.set_title("3b. link mirror: the halves meet across one slice T, so\n"
                    "<θF, F> = gᵀ T g, never negative exactly when T ≥ 0")
save(fig, "f3-a-mirror-makes-a-square")

# ---------------------------------------------------------------- f4 rewind slices and the spectrum condition
fig, ax = plt.subplots(1, 2, figsize=(12, 4.4))
e = np.linspace(0.02, 4, 400)
a = ax[0]
a.fill_between(e, -1 / e, 3, color=BLUE, alpha=0.18)
a.plot(e, -1 / e, color=BLUE, lw=2)
a.set_ylim(-3, 3)
a.set_xlabel("step size ε")
a.set_ylabel("k · E  (dial times energy)")
a.axhline(0, color=GREY, lw=1)
a.text(2.3, 1.4, "explicit slice 1 + kεE ≥ 0", color=BLUE, ha="center")
a.text(2.3, -2.4, "slice negative: positivity fails", color=RED, ha="center")
a.set_title("4a. explicit slices: a window that shrinks as steps grow")
a = ax[1]
a.fill_between(e, -3, 1 / e, color=GREEN, alpha=0.18)
a.plot(e, 1 / e, color=GREEN, lw=2)
a.axhspan(-3, 0, color=GREEN, alpha=0.18)
a.set_ylim(-3, 3)
a.axhline(0, color=GREY, lw=1)
a.set_xlabel("step size ε")
a.text(2.3, -1.6, "k·E ≤ 0: positive at EVERY step size\n(at k = −1/ħ this is E ≥ 0, the spectrum condition)",
       color=GREEN, ha="center", fontsize=9)
a.text(2.3, 1.9, "k·E > 0: fails once ε > 1/(kE)", color=RED, ha="center")
a.set_title("4b. rewind slices (1 − kεE)⁻¹, the ∘ₖ inverse of a backward step")
save(fig, "f4-rewind-slices-and-positive-energy")

# ---------------------------------------------------------------- f5 the turn rotates the pairing
th = np.linspace(0, np.pi / 2, 200)
c = 1.0
kpath = -np.exp(1j * th)                                  # from k = -1 to k = -i
pairs = (kpath * c) ** 2
fig, ax = plt.subplots(1, 2, figsize=(11, 4.6))
a = ax[0]
a.plot(kpath.real, kpath.imag, color=ORANGE, lw=2.5)
a.plot([-1], [0], "o", color=BLUE, ms=9)
a.plot([0], [-1], "o", color=PURPLE, ms=9)
a.text(-1, 0.12, "k = −1 (heat)", color=BLUE, ha="center")
a.text(0.05, -1.15, "k = −i (Feynman)", color=PURPLE)
a.axhline(0, color=GREY, lw=0.8)
a.axvline(0, color=GREY, lw=0.8)
a.set_aspect("equal")
a.set_xlim(-1.4, 0.6)
a.set_ylim(-1.4, 0.6)
a.set_title("5a. turn the dial a quarter turn")
a = ax[1]
sc = a.scatter(pairs.real, pairs.imag, c=th, cmap="plasma", s=8)
a.plot([1], [0], "o", color=BLUE, ms=9)
a.plot([-1], [0], "o", color=PURPLE, ms=9)
a.text(1, -0.18, "+c² : a positive number", color=BLUE, ha="center")
a.text(-1, -0.18, "−c² : negative", color=PURPLE, ha="center")
a.text(0, 0.45, "in between: not even\na real number", ha="center", color=RED)
a.axhline(0, color=GREY, lw=0.8)
a.axvline(0, color=GREY, lw=0.8)
a.set_aspect("equal")
a.set_xlim(-1.5, 1.5)
a.set_ylim(-0.4, 1.3)
a.set_title("5b. the simplest pairing <θF, F> = (k·c)² turns a half turn")
save(fig, "f5-the-turn-rotates-the-pairing")

# ---------------------------------------------------------------- f6 the turn delivers one J
fig, ax = plt.subplots(figsize=(11, 4.6))
ax.axis("off")
ax.set_xlim(0, 11)
ax.set_ylim(0, 5)


def box(x, y, w, h, text, col):
    ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.08", fc=col, ec=col, alpha=0.15))
    ax.text(x + w / 2, y + h / 2, text, ha="center", va="center", fontsize=9)


box(0.2, 3.0, 3.2, 1.6, "Euclidean side (real dial)\nshift energy H → H + c\nweights × e^(−ct)\nONE positive number for every H", BLUE)
box(7.4, 3.0, 3.4, 1.6, "after the turn t → is\nU → e^(−ics) U\nONE rotation J for every H\n= the paper's (A4)", PURPLE)
arrow(ax, 3.5, 3.8, 7.3, 3.8, ORANGE, lw=2.5)
ax.text(5.4, 4.05, "reflection positivity + Osterwalder–Schrader", ha="center", color=ORANGE, fontsize=9)
ax.text(5.4, 3.45, "i enters here, as one scalar", ha="center", color=ORANGE, fontsize=9)
box(0.2, 0.4, 3.2, 1.6, "real QM, all symmetric observables:\neach dynamics reads its own J_H\nfrom its own polar part", RED)
box(7.4, 0.4, 3.4, 1.6, "no single J\n(Corollary 1 of the paper)\nso it is NOT the turn of any\nreflection-positive family", RED)
arrow(ax, 3.5, 1.2, 7.3, 1.2, RED, lw=1.5, style="-|>")
ax.text(5.4, 1.45, "no Euclidean family turns into this", ha="center", color=RED, fontsize=9)
ax.set_title("6. What the turn hands the paper: (A4) is the Wick image of a trivial Euclidean fact")
save(fig, "f6-the-turn-delivers-one-J")

# ---------------------------------------------------------------- f7 scorecard
rows = [
    ("Real dial, symmetric slices: site reflection positive", "holds, every real k", "known (transfer matrix)", GREEN),
    ("Real dial: link reflection positive", "holds iff 1 + kεE ≥ 0", "known; window is ours to draw", GREEN),
    ("Rewind slices at every step size", "holds iff k·E ≤ 0 (E ≥ 0)", "dial reading ours; Chernoff step known", GREEN),
    ("Off the real axis (any Im k ≠ 0)", "fails: pairing not real", "short proof, ours", RED),
    ("Feynman's dial k = −i/ħ", "fails: pairing negative", "short proof, ours", RED),
    ("Reversible Markov chains at k = −1", "holds; one-way cycles fail", "known (Klein 1978)", GREEN),
    ("Turn of an RP family satisfies (A1)–(A4)", "yes, one shared J", "ours (short)", GREEN),
    ("Real QM as the turn of an RP family", "impossible", "ours, via Corollary 1", RED),
    ("Infinite dimensions, field theory", "not done here", "open for the framework", GREY),
]
fig, ax = plt.subplots(figsize=(12, 4.6))
ax.axis("off")
for i, (q, a_, src, col) in enumerate(rows):
    y = len(rows) - i
    ax.add_patch(FancyBboxPatch((0, y - 0.42), 12, 0.84, boxstyle="round,pad=0.02", fc=col, alpha=0.12, ec="none"))
    ax.text(0.15, y, q, va="center", fontsize=9)
    ax.text(6.4, y, a_, va="center", fontsize=9, color=col, fontweight="bold")
    ax.text(9.3, y, src, va="center", fontsize=8.5, color=GREY)
ax.set_xlim(0, 12)
ax.set_ylim(0.4, len(rows) + 0.6)
ax.set_title("7. Scorecard: where reflection positivity holds, where it fails, and who proved what")
save(fig, "f7-scorecard")
print("figures written to", OUT)
