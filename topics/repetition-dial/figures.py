"""Figures for the repetition dial (run: python3 figures.py; writes fig*.png here)."""
import numpy as np, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle, Circle
import os
HERE = os.path.dirname(os.path.abspath(__file__))
plt.rcParams.update({"font.size": 11, "axes.spines.top": False, "axes.spines.right": False})
COL = {1: "#c0392b", 0.5: "#e67e22", 0: "#7f8c8d", -0.5: "#2980b9", -1: "#1a5276"}

def coeff(k, n):
    p = 1.0
    for j in range(n): p *= (1 - j*k)/(j + 1)
    return p
def Z(k, x): return np.exp(x) if k == 0 else (1 + k*x)**(1/k)

# ---------- Figure 1: the dial as cells, and the chance of each number of copies
fig = plt.figure(figsize=(12, 7.2))
fig.suptitle("Figure 1. One dial k, from \"at most once\" to \"any number of times\"", fontsize=14, weight="bold")
ax = fig.add_axes([0.03, 0.55, 0.94, 0.36]); ax.axis("off"); ax.set_xlim(-1.25, 1.25); ax.set_ylim(-0.75, 1.0)
ax.annotate("", xy=(-1.18, -0.45), xytext=(1.18, -0.45), arrowprops=dict(arrowstyle="<->", lw=1.5))
stops = [(1, "k = 1\none exclusive cell\n(subsets)"), (0.5, "k = 1/2\ntwo exclusive cells\n(up to 2 copies)"),
(0, "k = 0\ninfinitely many tiny cells\n(classical, e^x)"),
(-0.5, "k = -1/2\ntwo deep bins"), (-1, "k = -1\none deep bin\n(multisets)")]
for k, lab in stops:
    X = -k  # draw "at most once" on the left
    ax.plot([X], [-0.45], "o", color="k", ms=5)
    ax.text(X, -0.72, lab, ha="center", va="bottom" if False else "top", fontsize=8.5)
    if k > 0:
        m = round(1/k); w = 0.24/m
        for i in range(m):
            ax.add_patch(Rectangle((X - 0.12 + i*w, 0.0), w*0.9, 0.18, fc="#fadbd8", ec="#c0392b"))
        ax.text(X, 0.28, f"each cell holds 0 or 1\nso at most {m} cop{'y' if m==1 else 'ies'}", ha="center", fontsize=8)
    elif k < 0:
        m = round(-1/k); w = 0.24/m
        for i in range(m):
            ax.add_patch(Rectangle((X - 0.12 + i*w, 0.0), w*0.9, 0.62, fc="#d6eaf8", ec="#1a5276"))
        ax.text(X, 0.72, "each bin holds\nany number", ha="center", fontsize=8)
    else:
        rng = np.random.default_rng(3)
        ax.scatter(X + rng.uniform(-0.12, 0.12, 120), rng.uniform(0, 0.4, 120), s=2, color="#7f8c8d")
        ax.text(X, 0.5, "cells so fine that\ncopies never collide", ha="center", fontsize=8)
ax.text(-1.2, 0.95, "number of cells per element = 1/k   (positive k: cells that hold one copy; negative k: bins that hold any number)", fontsize=9.5, style="italic")
nmean = 0.6
for i, k in enumerate([1, 0.5, 0, -0.5, -1]):
    a2 = fig.add_axes([0.05 + i*0.19, 0.08, 0.16, 0.3])
    x = nmean/(1 - k*nmean)
    js = np.arange(0, 6); P = np.array([coeff(k, j)*x**j for j in js])/Z(k, x)
    a2.bar(js, P, color=COL[k]); a2.set_ylim(0, 0.65); a2.set_xticks(js)
    a2.set_title(f"k = {k:g}", color=COL[k], fontsize=11)
    a2.set_xlabel("copies of the element")
    if i == 0: a2.set_ylabel("chance")
fig.text(0.5, 0.43, f"Same average ({nmean} copies) at every k: turning k down lets copies pile up. The change is smooth.", ha="center", fontsize=10.5)
fig.savefig(os.path.join(HERE, "fig1-repetition-dial.png"), dpi=130); plt.close(fig)

# ---------- Figure 2: average copies against energy (Fermi, classical, Bose)
fig, ax = plt.subplots(figsize=(8.5, 5.2))
eps = np.linspace(-2.5, 4, 400)
for k in [1, 0.5, 0, -0.5, -1]:
    e = eps if k >= 0 else eps[eps > 0.02]
    n = 1/(np.exp(e) + k)
    lab = {1: "k = 1  Fermi-Dirac (at most once)", 0: "k = 0  classical (Boltzmann)", -1: "k = -1  Bose-Einstein (any number)"}.get(k, f"k = {k:g}")
    ax.plot(e, n, color=COL[k], lw=2.4 if k in (1, 0, -1) else 1.4, label=lab)
ax.set_ylim(0, 3); ax.set_xlabel("cost of one copy  (energy over temperature)"); ax.set_ylabel("average number of copies n")
ax.set_title("Figure 2. Average copies: n = 1 / (e^cost + k), one formula for all k", weight="bold")
ax.axhline(1, ls=":", color="#c0392b", lw=1); ax.text(-2.4, 1.05, "k = 1 can never pass 1 copy", color="#c0392b", fontsize=9)
ax.text(0.65, 2.05, "k = -1 piles up near zero cost\n(the wall at -1/k)", color="#1a5276", fontsize=9)
ax.legend(frameon=False, fontsize=9, loc="upper right")
fig.tight_layout(); fig.savefig(os.path.join(HERE, "fig2-average-copies.png"), dpi=130); plt.close(fig)

# ---------- Figure 3: the one-dial test: effective k = 1/n - 1/x against crowding
def wu_n(g, x):
    lo, hi = -60.0, 60.0
    for _ in range(200):
        mid = (lo + hi)/2
        f = g*mid + (1 - g)*np.log1p(np.exp(mid)) + np.log(x)
        hi, lo = (mid, lo) if f > 0 else (hi, mid)
    w = np.exp((lo + hi)/2); return 1/(w + g)
fig, ax = plt.subplots(figsize=(9, 5.6))
xs = np.logspace(-4, 4, 300)
for k in [1, 0.5, 0, -0.5, -1]:
    ax.plot([0, 4.2], [k, k], color=COL[k], lw=2, alpha=0.8)
    ax.text(4.22, k, f"k = {k:g}", color=COL[k], va="center", fontsize=9)
for g, c in [(0.25, "#8e44ad"), (0.5, "#16a085"), (0.75, "#d35400")]:
    n = np.array([wu_n(g, x) for x in xs]); ke = 1/n - 1/xs
    ax.plot(n, ke, "--", color=c, lw=2.2, label=f"Haldane exclusion g = {g}: slides from {2*g-1:g} to {g:g}")
m = 2
n = np.array([sum(j*x**j for j in range(m+1))/sum(x**j for j in range(m+1)) for x in xs])
ax.plot(n, 1/n - 1/xs, ":", color="k", lw=2.2, label="hard cap of 2 copies: slides from -1 to 1/2")
ax.set_xlim(0, 4.2); ax.set_ylim(-1.15, 1.15)
ax.set_xlabel("crowding: average copies n"); ax.set_ylabel("measured k  =  1/n - 1/x")
ax.set_title("Figure 3. The one-dial test: a true k is a flat line at every crowding", weight="bold")
ax.legend(frameon=False, fontsize=9, loc="lower right")
ax.text(0.05, 1.07, "flat lines: the k family (one number describes the element)", fontsize=9)
fig.tight_layout(); fig.savefig(os.path.join(HERE, "fig3-one-dial-test.png"), dpi=130); plt.close(fig)

# ---------- Figure 4: which k give honest counts
ks = np.linspace(-1.2, 1.2, 1201); ns = np.arange(0, 13)
S = np.zeros((len(ns), len(ks)))
for i, n in enumerate(ns):
    for j, k in enumerate(ks):
        f = [1 - jj*k for jj in range(n)]          # sign of prod(1 - jk); n! is positive
        S[i, j] = 0 if any(abs(v) < 1e-12 for v in f) else (-1)**sum(v < 0 for v in f)
fig, ax = plt.subplots(figsize=(9.5, 5))
from matplotlib.colors import ListedColormap
ax.imshow(S, origin="lower", aspect="auto", extent=[ks[0], ks[-1], -0.5, 12.5],
          cmap=ListedColormap(["#e74c3c", "#ecf0f1", "#5dade2"]), vmin=-1, vmax=1, interpolation="nearest")
for m in range(1, 6):
    nn = np.arange(0, 13)
    ax.scatter([1/m]*(m+1), nn[:m+1], s=28, color="#1b4f72", zorder=5)
    ax.scatter([1/m]*(12-m), nn[m+1:], s=28, facecolor="white", edgecolor="#555", zorder=5)
    ax.text(1/m, -1.3, f"1/{m}" if m > 1 else "1", ha="center", fontsize=8.5, color="#1b4f72")
ax.set_ylim(-1.8, 12.6)
ax.set_xlabel("k"); ax.set_ylabel("number of copies n")
ax.set_title("Figure 4. Which k give honest counts? Blue: positive. Red: negative. Dots: the lines k = 1/m", weight="bold", fontsize=11)
ax.text(-1.15, 10.8, "k ≤ 0: every number of copies\nhas a positive count", fontsize=9, color="#1b4f72")
ax.text(-0.75, 6.5, "dots on k = 1/m: filled = a real count (up to m copies),\nopen = zero (more than m copies is impossible).\nBetween those lines, some counts turn negative (red).", fontsize=9, color="#1b4f72", bbox=dict(fc="white", ec="none", alpha=0.85))
fig.tight_layout(); fig.savefig(os.path.join(HERE, "fig4-honest-counts.png"), dpi=130); plt.close(fig)
print("figures written")
