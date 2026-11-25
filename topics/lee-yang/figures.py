import pickle, math, cmath, numpy as np, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from dialgas import cycle_counts, z_zeros
BLUE, ORANGE, AQUA, VIOLET, INK, GREY = "#257de0", "#f47230", "#00b894", "#6a4ce2", "#191919", "#8a8a8a"
plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 10, "axes.edgecolor": "#bbbbbb",
                     "axes.labelcolor": INK, "xtick.color": "#555", "ytick.color": "#555"})

def squash(z):  # whole plane into the unit disk; |zeta|=1 lands on radius 1/2, 0 and infinity at centre and rim
    r = abs(z); return 0j if r == 0 else (2 / math.pi) * math.atan(r) * z / r

def disk(ax, title):
    t = np.linspace(0, 2 * np.pi, 400)
    ax.plot(np.cos(t), np.sin(t), color="#dddddd", lw=1)
    ax.plot(0.5 * np.cos(t), 0.5 * np.sin(t), color="#bbbbbb", lw=1, ls="--")
    ax.axhline(0, color="#e6e6e6", lw=0.8, zorder=0)
    ax.plot([-0.5], [0], marker="o", ms=11, mfc="none", mec=ORANGE, mew=2, zorder=5)
    ax.plot([0.5], [0], marker="s", ms=9, mfc="none", mec=VIOLET, mew=2, zorder=5)
    ax.set_xlim(-1.08, 1.08); ax.set_ylim(-1.08, 1.08); ax.set_aspect("equal"); ax.axis("off")
    ax.set_title(title, fontsize=10, color=INK)

# Figure 1: the 2D grid (6x6 torus) zeros in zeta = k^2 z across the dial
Z = pickle.load(open("torus_zeros.pkl", "rb"))
ks = [(4, "k = 4  (q = 5)"), (2, "k = 2  (q = 3)"), (1.0001, "k = 1  (q = 2)"), (0.5, "k = 0.5  (q = 1.5)"),
      (0.1, "k = 0.1  (q = 1.1)"), (-1, "k = −1  (q = 0)"), (-2, "k = −2  (q = −1)")]
fig, axs = plt.subplots(1, 7, figsize=(16, 2.9))
for ax, (k, lab) in zip(axs, ks):
    disk(ax, lab)
    for L, col, ms in ((5, AQUA, 3.5), (6, BLUE, 4.5)):
        zs = [squash(z * k * k) for z in Z[(L, k)]]
        ax.scatter([w.real for w in zs], [w.imag for w in zs], s=ms ** 2, color=col, zorder=4, lw=0)
fig.text(0.5, 0.01, "Each disk is the whole complex plane of ζ = k²·z, squeezed so the unit circle (dashed) sits halfway out. "
         "Orange ring: the wall ζ = −1.  Violet square: zero field ζ = +1, where a phase change shows up.  "
         "Blue: 6×6 grid.  Aqua: 5×5 grid.", ha="center", fontsize=9, color="#444")
plt.subplots_adjust(left=0.01, right=0.99, top=0.86, bottom=0.12, wspace=0.05)
plt.savefig("figures/f1_grid_zeros.png", dpi=160); plt.close()

# Figure 2: Potts colour counts in the dial plane k = Q - 1
P = pickle.load(open("potts_zeros.pkl", "rb"))
fig, ax = plt.subplots(figsize=(8.6, 4.6))
ax.axvspan(-3.2, -1, color="#f2f2f2", zorder=0)
ax.text(-2.6, 1.45, "no zeros, ever\n(no graph has a negative\ncolour count root)", ha="center", fontsize=9, color="#555")
sty = {1: (GREY, "o", 12, "ring, 40 sites"), 2: (AQUA, "o", 16, "width-2 strip, 2×20"),
       3: (BLUE, "o", 18, "width-3 strip, 3×14"), 4: (VIOLET, "D", 18, "width-4 strip, 4×10")}
for L, n in ((1, 40), (2, 20), (3, 14), (4, 10)):
    zs = [z - 1 for z in P[(L, n)]]
    c, m, s, lab = sty[L]
    ax.scatter([z.real for z in zs], [z.imag for z in zs], s=s, color=c, marker=m, label=lab, zorder=4, lw=0)
for kv, name in ((-2, "k = −2\nparity"), (-1, "k = −1\nor"), (0, "k = 0\naddition"), (1, "k = 1\ncircle operation")):
    ax.axvline(kv, color=ORANGE if kv > -2 else "#bbbbbb", lw=1, ls=":", zorder=1)
    ax.text(kv, -2.25, name, ha="center", va="top", fontsize=9, color=INK)
ax.set_xlim(-3.2, 3.6); ax.set_ylim(-2.0, 2.0); ax.axhline(0, color="#cccccc", lw=0.8)
ax.set_xlabel("dial k = Q − 1   (Q = number of colours)", labelpad=34); ax.set_ylabel("imaginary part")
ax.legend(loc="upper right", fontsize=8.5, frameon=False)
for sp in ("top", "right"): ax.spines[sp].set_visible(False)
plt.tight_layout(); plt.savefig("figures/f2_colour_zeros.png", dpi=160); plt.close()

# Figure 3: 2D Fisher zeros in the dial k at zero field (full view and close-up near k = 0)
Fz = pickle.load(open("fisher_zeros.pkl", "rb"))
kc = 3 + 2 * math.sqrt(2)
fig, (ax, bx) = plt.subplots(1, 2, figsize=(11, 4.2), gridspec_kw={"width_ratios": [1.6, 1]})
for a_, lim in ((ax, (-3, 8, -5, 5)), (bx, (-0.8, 0.8, -0.65, 0.65))):
    for L, c in ((4, GREY), (5, AQUA), (6, BLUE)):
        a_.scatter([z.real for z in Fz[L]], [z.imag for z in Fz[L]], s=18, color=c, label=f"{L}×{L} grid", zorder=4, lw=0)
    for kv in (-2, -1, 0, 1):
        a_.plot([kv], [0], marker="o", ms=9, mfc="none", mec=ORANGE, mew=2, zorder=5)
    for kv in (kc, 1 / kc):
        a_.plot([kv], [0], marker="s", ms=9, mfc="none", mec=VIOLET, mew=2, zorder=6)
    a_.axhline(0, color="#cccccc", lw=0.8); a_.set_xlim(*lim[:2]); a_.set_ylim(*lim[2:])
    for sp in ("top", "right"): a_.spines[sp].set_visible(False)
ax.text(kc, -0.6, "k = 3+2√2\nmagnet orders", ha="center", va="top", fontsize=9)
ax.text(-0.5, -4.4, "orange rings: k = −2, −1, 0, 1", ha="center", va="top", fontsize=9, color=ORANGE)
ax.add_patch(plt.Rectangle((-0.8, -0.65), 1.6, 1.3, fill=False, ec="#999", lw=0.8, ls="--"))
bx.text(0.5, -0.05, "k = 3−2√2\nchessboard\norders", ha="center", va="top", fontsize=9)
bx.set_title("close-up of the dashed box", fontsize=9, color="#555")
ax.set_xlabel("dial k (pair weight)"); ax.set_ylabel("imaginary part"); bx.set_xlabel("dial k")
ax.legend(loc="upper right", fontsize=8.5, frameon=False)
plt.tight_layout(); plt.savefig("figures/f3_fisher_zeros.png", dpi=160); plt.close()
print("figures written")
