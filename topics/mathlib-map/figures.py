"""Figures for the Mathlib map page. Run: python3 figures.py"""
import math, os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Circle, FancyArrowPatch

OUT = os.path.join(os.path.dirname(__file__), "figures")
os.makedirs(OUT, exist_ok=True)
INK, GREY, LIGHT, ORANGE, BLUE = "#141517", "#7a7d83", "#c9cbcf", "#ff4f00", "#2a5bd7"
plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 12, "axes.edgecolor": GREY,
                     "axes.labelcolor": INK, "xtick.color": INK, "ytick.color": INK,
                     "axes.spines.top": False, "axes.spines.right": False})

def save(fig, name):
    fig.savefig(os.path.join(OUT, name), dpi=160, bbox_inches="tight", facecolor="white")
    plt.close(fig)

# ---------- F2: what is built on the circle operation ----------
clusters = [  # name, files, theorems, x, y
    ("Circle operation\nPreQuasiregular,\nquasispectrum", 2, 79, 0.0, 0.0),
    ("Unitization", 8, 227, 2.4, 2.2),
    ("Neumann series", 1, 112, 2.4, -2.3),
    ("Spectrum and order", 9, 245, 2.4, 0.0),
    ("Non-unital functional calculus", 16, 459, 4.7, 0.6),
    ("Operator functions\nrpow, abs, positive part", 7, 268, 7.4, 1.9),
    ("C*-structure", 10, 159, 7.4, -1.3),
]
edges = [(0, 1), (0, 2), (0, 3), (1, 4), (3, 4), (4, 5), (4, 6)]
ours = [  # label, attaches to cluster, x, y
    ("Every dial k ≠ 0 is this monoid", 0, -1.3, 2.3),
    ("Walls = quasispectrum, flipped", 0, -2.3, 1.0),
    ("No wall while ‖k‖‖a‖ < 1", 2, 5.0, -2.9),
]
fig, ax = plt.subplots(figsize=(12, 6.2))
def r(n): return 0.22 + 0.034 * math.sqrt(n)
def trim(p, q, ra, rb):
    v = np.array(q, float) - np.array(p, float); L = np.hypot(*v); u = v / L
    return tuple(np.array(p) + u * (ra + 0.06)), tuple(np.array(q) - u * (rb + 0.06))
for i, j in edges:
    a, b = clusters[i], clusters[j]
    P, Q = trim((a[3], a[4]), (b[3], b[4]), r(a[2]), r(b[2]))
    ax.add_patch(FancyArrowPatch(P, Q, arrowstyle="-|>", mutation_scale=16, color=GREY, lw=1.6))
for k, (name, f, t, x, y) in enumerate(clusters):
    fc = INK if k == 0 else "#e9eaec"
    ax.add_patch(Circle((x, y), r(t), fc=fc, ec=INK, lw=1.2, zorder=3))
    ax.text(x, y, f"{t}", ha="center", va="center", fontsize=13, fontweight="bold",
            color="white" if k == 0 else INK, zorder=4)
    ax.text(x, y - r(t) - 0.12, f"{name}\n{f} file{'s' if f > 1 else ''}", ha="center", va="top",
            fontsize=10.5, color=INK, zorder=4)
for lab, c, x, y in ours:
    cx, cy = clusters[c][3], clusters[c][4]
    P, Q = trim((x, y), (cx, cy), 0.25, r(clusters[c][2]))
    ax.add_patch(FancyArrowPatch(P, Q, arrowstyle="-|>", mutation_scale=14, color=ORANGE, lw=2, linestyle="--"))
    ax.text(x, y, lab, ha="center", va="center", fontsize=10.5, color=ORANGE, fontweight="bold",
            bbox=dict(fc="white", ec=ORANGE, boxstyle="round,pad=0.35"))
ax.text(-3.4, 3.45, "About 1,550 theorems in 53 files use the circle operation directly.\n"
        "1,259 of Mathlib's 8,736 files import it.  Orange: our new links.",
        fontsize=11, color=GREY, va="top")
ax.set_xlim(-3.6, 8.8); ax.set_ylim(-3.6, 3.5); ax.set_aspect("equal"); ax.axis("off")
save(fig, "f2-circle-dependents.png")

# ---------- F3: walls along the dial ----------
A = np.array([[2.0, 1.5, 0.0], [0.0, 0.5, 1.2], [0.0, 0.0, -1.0]])  # eigenvalues 2, 0.5, -1
lam = np.linalg.eigvals(A).real
nrm = np.linalg.norm(A, 2); rho = max(abs(lam))
fig, (a1, a2) = plt.subplots(1, 2, figsize=(12, 4.8), gridspec_kw={"width_ratios": [1.35, 1]})
k = np.linspace(-3, 3, 4001)
d = np.abs(np.prod([1 + k * l for l in lam], axis=0))
a1.semilogy(k, np.maximum(d, 1e-4), color=INK, lw=2)
for l in lam:
    w = -1 / l
    a1.axvline(w, color=ORANGE, lw=1.5, ls="--")
    a1.text(w, 30, f"wall\nk = {w:g}\n(λ = {l:g})", color=ORANGE, ha="center", va="bottom", fontsize=10,
            bbox=dict(fc="white", ec="none", pad=1))
a1.axvspan(-1 / nrm, 1 / nrm, color=BLUE, alpha=0.10)
a1.text(0, 2.2e-4, "safe by\nNeumann", color=BLUE, ha="center", va="bottom", fontsize=10)
a1.set_ylim(1e-4, 1e3); a1.set_xlabel("dial k"); a1.set_ylabel("|det(1 + kA)|")
a1.set_title("A matrix A loses its way back exactly at k = −1/λ", loc="left", fontsize=12.5)
t = np.linspace(0, 2 * np.pi, 400)
a2.fill(np.cos(t) / rho, np.sin(t) / rho, color=LIGHT, alpha=0.6, lw=0)
a2.fill(np.cos(t) / nrm, np.sin(t) / nrm, color=BLUE, alpha=0.25, lw=0)
for l in lam:
    a2.plot(-1 / l, 0, "o", color=ORANGE, ms=9, zorder=5)
a2.text(0, 1 / nrm * 0.15, f"Neumann disc\nradius 1/‖A‖ = {1/nrm:.2f}", ha="center", color=BLUE, fontsize=10)
a2.text(0, -1 / rho * 0.85, f"wall-free disc  radius 1/ρ = {1/rho:.2f}", ha="center", color=GREY, fontsize=10)
a2.axhline(0, color=LIGHT, lw=0.8); a2.axvline(0, color=LIGHT, lw=0.8)
a2.set_xlim(-2.4, 1.4); a2.set_ylim(-1.2, 1.2); a2.set_aspect("equal")
a2.set_xlabel("Re k"); a2.set_ylabel("Im k")
a2.set_title("Complex dial: walls (orange) sit outside both discs", loc="left", fontsize=12.5)
save(fig, "f3-walls.png")

# ---------- F4: Jacobson grid for Z/12 ----------
n = 12
grid = np.array([[math.gcd((1 + kk * x) % n, n) == 1 for kk in range(n)] for x in range(n)])
fig, ax = plt.subplots(figsize=(7.4, 6.6))
from matplotlib.colors import ListedColormap
ax.imshow(grid, cmap=ListedColormap([ORANGE, "#e9eaec"]), aspect="equal")
for x in range(n):
    if grid[x].all():
        ax.add_patch(plt.Rectangle((-0.5, x - 0.5), n, 1, fill=False, ec=INK, lw=2.4))
ax.set_xticks(range(n)); ax.set_yticks(range(n))
ax.set_xlabel("dial k"); ax.set_ylabel("element x of ℤ/12")
ax.set_title("ℤ/12: orange = x hits a wall at dial k (1 + kx not a unit)\n"
             "Boxed rows never hit one: {0, 6} = the Jacobson radical", loc="left", fontsize=12)
for s in ax.spines.values(): s.set_visible(False)
save(fig, "f4-jacobson-grid.png")

# ---------- F5: open problems map ----------
probs = [  # code, label, close (0-10), matters (0-10), ours?
    ("O1", "Lee–Yang p-ball band", 6.5, 5.5, True),
    ("O2", "Soft min-cut keeps\nthe entropy rules", 5.0, 4.5, True),
    ("K3", "Group algebras are\nJacobson semisimple", 3.6, 6.5, False),
    ("K4", "Tree independence\npolynomials unimodal", 3.0, 4.8, False),
    ("K5", "Holographic entropy\ncone, 6+ parties", 3.2, 5.6, False),
    ("K2", "Köthe conjecture", 2.0, 7.4, False),
    ("J", "Jacobson conjecture\n(⋂ Jⁿ = 0)", 1.4, 6.0, False),
    ("K1", "Riemann hypothesis via\nPólya / de Bruijn–Newman", 1.0, 9.8, False),
    ("K6", "Navier–Stokes", 0.5, 9.3, False),
]
fig, ax = plt.subplots(figsize=(11, 6.4))
off = {"O1": (0.25, 0.25), "O2": (0.25, -0.55), "K3": (0.25, 0.1), "K4": (0.25, -0.5), "K5": (0.25, 0.15),
       "K2": (0.25, 0.1), "J": (0.25, -0.55), "K1": (0.25, -0.15), "K6": (0.25, -0.12)}
for c, lab, x, y, mine in probs:
    ax.plot(x, y, "o", ms=16, color=ORANGE if mine else INK, zorder=3)
    ax.text(x, y, c, ha="center", va="center", color="white", fontsize=8.5, fontweight="bold", zorder=4)
    dx, dy = off[c]
    ax.text(x + dx, y + dy, lab, fontsize=10.5, color=ORANGE if mine else INK, va="bottom")
ax.set_xlim(0, 8.6); ax.set_ylim(3.5, 10.6)
ax.set_xlabel("how close our tools look  →"); ax.set_ylabel("how much a proof would matter  →")
ax.set_xticks([]); ax.set_yticks([])
ax.set_title("Open problems near our dials (placement is judgment, not measurement). Orange = ours",
             loc="left", fontsize=12.5)
save(fig, "f5-open-problems.png")
print("ok", lam, nrm, rho, [x for x in range(n) if grid[x].all()])
