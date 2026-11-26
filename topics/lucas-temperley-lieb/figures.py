"""Pictures for lucas-temperley-lieb/README.md.

Reads checks.json (Gram eigenvalues, written by checks.py) and
../fractional-steps/grid_output.json (the fractional-step chances) for the same-shape panel.
"""
import json
import math
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "figures")
os.makedirs(OUT, exist_ok=True)
G, R, GR, B = "#2b8a3e", "#c92a2a", "#868e96", "#1c7ed6"
CM = plt.get_cmap("viridis")


def save(fig, name):
    fig.tight_layout()
    fig.savefig(os.path.join(OUT, name), dpi=140)
    plt.close(fig)
    print("wrote", name)


def lucas(k, d, mmax):
    U = [0.0, 1.0]
    for _ in range(mmax - 1):
        U.append(k * U[-1] - d * U[-2])
    return U


# Beraha numbers 4cos^2(pi j/m), with the first layer m at which each appears
beraha = {}
for m in range(2, 8):
    for j in range(1, m):
        b = round(4 * math.cos(math.pi * j / m) ** 2, 9)
        if b not in beraha and b < 3.999:
            beraha[b] = (m, j)

# ---------------------------------------------------------------- Fig 1: the dial plane
fig, ax = plt.subplots(figsize=(8.2, 5.6))
d = np.linspace(0, 2.6, 400)
ax.fill_between(d, -2 * np.sqrt(d), 2 * np.sqrt(d), color="#f1f3f5", zorder=0)
ax.plot(d, 2 * np.sqrt(d), color="k", lw=2)
ax.plot(d, -2 * np.sqrt(d), color="k", lw=2)
ax.axvspan(-1.2, 0, color="#e7f5ff", zorder=0)
ylab_prev = -9
for b, (m, j) in sorted(beraha.items()):
    col = CM((m - 2) / 6)
    ax.plot(d, np.sqrt(b) * np.sqrt(d), color=col, lw=1.4, ls="--")
    lab = f"layer {m} silent: k²/d = {b:.3g}" + (f"  (j = {j})" if j > 1 else "")
    yl = max(np.sqrt(b * 2.6), ylab_prev + 0.2); ylab_prev = yl
    ax.text(2.62, yl, lab, fontsize=7.5, color=col, va="center")
ax.axvline(1, color=GR, lw=0.8, ls=":")
for m in range(2, 8):
    k = 2 * math.cos(math.pi / m)
    ax.plot(1, k, "o", color=CM((m - 2) / 6), ms=6, zorder=5)
ax.annotate("d = 1: silent at k = 2cos(π/m)\n0, 1, √2, φ, √3, ...", (1, 1.0), (0.25, 2.55), fontsize=8,
            arrowprops=dict(arrowstyle="->", color=GR))
ax.text(1.7, -1.0, "inside the parabola k² < 4d\nloop value δ = k/√d below 2:\nsilent layers at the\nBeraha points (Jones' discrete series)",
        fontsize=8.5, ha="center")
ax.text(0.15, 3.05, "on and outside the parabola: δ ≥ 2, or d < 0\nno layer is silent, except k = 0 when d < 0\n(Jones' continuous range)", fontsize=8.5)
ax.plot([-1.2, 0], [0, 0], color=CM(0), lw=1.4, ls="--")
ax.text(-1.15, -2.9, "d < 0: δ imaginary; only the line\nk = 0 is silent (every even layer)", fontsize=8, color=B)
ax.set_xlim(-1.2, 3.9)
ax.set_ylim(-3.3, 4.0)
ax.set_xlabel("plane dial d")
ax.set_ylabel("plane dial k")
ax.set_title("Fig 1. The plane dial (k, d) read as a Temperley–Lieb loop value δ = k/√d\n"
             "each dashed curve is one silent layer of the Lucas weights", fontsize=10)
save(fig, "fig1-dial-plane.png")

# ---------------------------------------------------------------- Fig 2: weights at and near a Beraha point
fig, axs = plt.subplots(1, 3, figsize=(10.5, 3.6), sharey=True)
ms = np.arange(1, 9)
cases = [(1.55, "k = 1.55, just below φ:\nweight 5 goes negative"),
         (2 * math.cos(math.pi / 5), "k = φ = 2cos(π/5), a Beraha point:\nfirst failure is exactly zero (layer 5)"),
         (1.70, "k = 1.70, just above φ:\nlayers 1 to 5 positive, layer 6 negative")]
for ax, (k, title) in zip(axs, cases):
    U = lucas(k, 1, 8)[1:9]
    cols = []
    seen_bad = False
    for u in U:
        if abs(u) < 1e-12:
            cols.append("k"); seen_bad = True
        elif u < 0:
            cols.append(R); seen_bad = True
        else:
            cols.append(GR if seen_bad else G)
    ax.bar(ms, U, color=cols)
    ax.axhline(0, color="k", lw=0.8)
    for m, u in zip(ms, U):
        if abs(u) < 1e-12:
            ax.text(m, 0.08, "0", ha="center", fontsize=11, fontweight="bold")
    ax.set_title(title, fontsize=9)
    ax.set_xlabel("layer m")
    ax.set_xticks(ms)
axs[0].set_ylabel("Lucas weight U_m(k, 1)\n= Jones–Wenzl trace Δ_{m−1}(δ = k)")
fig.suptitle("Fig 2. Green = positive, red = negative, black zero = silent, grey = after the first failure", fontsize=10)
save(fig, "fig2-weights-at-beraha.png")

# ---------------------------------------------------------------- Fig 3: positivity, and the same shape
psd = json.load(open(os.path.join(HERE, "checks.json")))["psd"]
fig, axs = plt.subplots(1, 2, figsize=(11.5, 4.4))
ax = axs[0]
for n in ["2", "3", "4", "5", "6"]:
    g = np.array(psd[n]["grid"]); e = np.array(psd[n]["min_eig"])
    keep = g >= 0.3
    col = CM((int(n) - 2) / 5)
    ax.plot(g[keep], e[keep], color=col, lw=1.3, label=f"n = {n}: positive from 2cos(π/{int(n)+1}) = {psd[n]['threshold']:.3f}")
    ax.axvline(psd[n]["threshold"], color=col, lw=0.7, ls=":")
ax.set_yscale("symlog", linthresh=1e-3)
ax.axhline(0, color="k", lw=0.8)
for m in range(3, 8):
    ax.plot(2 * math.cos(math.pi / m), 0, "o", color="k", ms=4, zorder=6)
ax.set_xlabel("loop value δ")
ax.set_ylabel("smallest eigenvalue of the\nTemperley–Lieb Gram matrix (symlog)")
ax.set_title("(a) Diagrams on n strands: the inner product is positive only above a\nthreshold, or exactly at δ = 2cos(π/m) (black dots: zero, never negative)", fontsize=9)
ax.legend(fontsize=7, loc="lower right")

ax = axs[1]
grid = json.load(open(os.path.join(HERE, "..", "fractional-steps", "grid_output.json")))
ts, vals = [], []
for t, v in grid.items():
    try:
        tt = float(t)
    except ValueError:
        continue
    ts.append(tt); vals.append(v["min_rel"])
ts = np.array(ts); vals = np.array(vals)
o = np.argsort(ts); ts, vals = ts[o], vals[o]
depth = np.log10(np.maximum(-vals, 1e-30))
whole = np.isclose(ts, np.round(ts))
ax.plot(ts, depth, color=R, lw=1.2)
ax.plot(ts[~whole], depth[~whole], "o", color=R, ms=3)
ax.plot(ts[whole], depth[whole], "o", color=G, ms=7, zorder=5, label="whole steps: no negative chance (noise floor)")
ax.legend(fontsize=7.5, loc="lower left")
ax.axhspan(-31, -20, color="#f1f3f5", zorder=0)
ax.text(1.5, -24, "numerical noise", ha="center", fontsize=8, color=GR)
ax.set_xlabel("step size t of the rung step (fractional-steps/)")
ax.set_ylabel("log10 of the most negative chance\n(relative, first 56 counts)")
ax.set_title("(b) Fractional steps of the rung step: a negative chance at every\nnon-whole t, none at whole t", fontsize=9)
fig.suptitle("Fig 3. The same shape: positivity survives only on a discrete set (plus, for loops, everything from δ = 2 up)", fontsize=10)
save(fig, "fig3-same-shape.png")

# ---------------------------------------------------------------- Fig 4: 4-ended vs 6-ended tangles
rng = np.random.default_rng(1)
fig, axs = plt.subplots(1, 2, figsize=(11, 4.6), sharey=True)
for ax, ends in zip(axs, [4, 6]):
    d = np.linspace(0, 3, 300)
    ax.fill_between(d, -2 * np.sqrt(d), 2 * np.sqrt(d), color="#f1f3f5", zorder=0)
    ax.plot(d, 2 * np.sqrt(d), "k", lw=1.6); ax.plot(d, -2 * np.sqrt(d), "k", lw=1.6)
    delta = 1.5
    pts = []
    for _ in range(4000):
        if ends == 4:
            x, y = rng.normal(size=2)
            # x*1 + y*e on TL_2: its square law N^2 = kN - D has k = 2x + delta*y, D = x^2 + delta*x*y
            k, D = 2 * x + delta * y, x * x + delta * x * y
        else:
            a, b, c, f, g = rng.normal(size=5)
            e1 = np.array([[delta, 1], [0, 0]]); e2 = np.array([[0, 0], [1, delta]])
            M = a * np.eye(2) + b * e1 + c * e2 + f * e1 @ e2 + g * e2 @ e1
            k, D = np.trace(M), np.linalg.det(M)
        # rescale N -> lam N so the point lands in the frame; (k, D) -> (lam k, lam^2 D) keeps k^2/D
        lam = min(1.0, 2.8 / max(abs(k), 1e-9), math.sqrt(2.8 / max(abs(D), 1e-9)))
        pts.append((lam * k, lam * lam * D))
    pts = np.array(pts)
    inside = pts[:, 0] ** 2 < 4 * pts[:, 1] - 1e-12
    ax.plot(pts[~inside, 1], pts[~inside, 0], ".", color=B, ms=2, alpha=0.4)
    ax.plot(pts[inside, 1], pts[inside, 0], ".", color=R, ms=2, alpha=0.6)
    ax.set_xlim(-3, 3); ax.set_ylim(-3.6, 3.6)
    ax.set_xlabel("d")
    if ends == 4:
        ax.set_title(f"(a) 4-ended tangles (TL₂, δ = {delta}): 4000 random elements,\n"
                     f"{inside.sum()} land inside the parabola — always k² − 4d = (δy)² ≥ 0", fontsize=9)
    else:
        ax.plot([1.5 ** 2 - 1], [0], "*", color="k", ms=14, zorder=6)
        ax.annotate("e₁e₂ − e₂e₁: k = 0, d = δ² − 1", (1.25, 0), (-2.6, 3.25), fontsize=8.5,
                    arrowprops=dict(arrowstyle="->"), bbox=dict(fc="white", ec="k", lw=0.5))
        ax.set_title(f"(b) 6-ended tangles (TL₃, δ = {delta}), acting on its 2-dim block:\n"
                     f"{inside.sum()} of 4000 random elements land inside (red)", fontsize=9)
axs[0].set_ylabel("k")
fig.suptitle("Fig 4. The knot page's open question: 6-ended tangles do reach inside the parabola; 4-ended ones never do\n"
             "(each element N is shrunk to λN to fit the frame; that keeps k²/d, so inside stays inside)", fontsize=10)
save(fig, "fig4-tangles.png")
