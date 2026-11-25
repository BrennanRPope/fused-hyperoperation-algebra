"""Figures for the zoom dial (writes fig*.png next to this file)."""
import numpy as np, matplotlib, os
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle, FancyBboxPatch, FancyArrowPatch
np.seterr(all="ignore")
HERE = os.path.dirname(os.path.abspath(__file__))
plt.rcParams.update({"font.size": 11, "axes.spines.top": False, "axes.spines.right": False})

# ---------- Fig 1: zoom becomes a weight
fig, axes = plt.subplots(1, 2, figsize=(12, 4.6))
def draw(ax, rs, gaps, title, d):
    ax.axis("off"); ax.set_xlim(-0.05, 1.05); ax.set_ylim(-0.5, 3.6); ax.set_title(title, weight="bold")
    level = [(0.0, 1.0, 1.0)]
    for depth in range(4):
        y = 3 - depth
        for (a, L, w) in level:
            ax.add_patch(Rectangle((a, y), L, 0.32, fc=plt.cm.Blues(0.25 + 0.75*w), ec="#1b4f72", lw=0.6))
            if depth <= 2 and L > 0.05:
                ax.text(a + L/2, y + 0.42, f"weight {w:.3g}" if depth else "weight 1", ha="center", fontsize=8)
        new = []
        for (a, L, w) in level:
            pos = a
            for r, g in zip(rs, gaps):
                pos += g*L
                new.append((pos, r*L, w*r**d)); pos += r*L
        level = new
        ax.text(-0.04, y + 0.12, f"depth {depth}", ha="right", fontsize=8)
    ax.text(0.5, -0.4, "every depth's weights add to 1, which is what makes D the dimension", ha="center", fontsize=9, style="italic")
draw(axes[0], [0.25, 0.25], [0, 0.5], "Half box A: two pieces at zoom 1/4\neach carries weight (1/4)^(1/2) = 1/2", 0.5)
draw(axes[1], [0.09, 0.49], [0, 0.42], "Uneven half box D: zooms 9/100 and 49/100\nweights 3/10 and 7/10", 0.5)
fig.suptitle("Figure 1. A zoom becomes a weight: a piece shrunk by r carries r^D", fontsize=13, weight="bold", y=1.0)
fig.tight_layout(); fig.savefig(os.path.join(HERE, "fig1-zoom-is-weight.png"), dpi=130, bbox_inches="tight"); plt.close(fig)

# ---------- Fig 2: dimension is where the zoom sum hits the wall
SETS = [("A: two at 1/4", [0.25]*2, "#c0392b"), ("B: three at 1/9", [1/9]*3, "#e67e22"),
        ("D: 9/100 and 49/100", [0.09, 0.49], "#16a085"), ("E: four at 1/16", [1/16]*4, "#8e44ad"),
        ("middle-third Cantor", [1/3]*2, "#2c3e50"), ("gasket: three at 1/2", [0.5]*3, "#7f8c8d")]
fig, ax = plt.subplots(figsize=(9, 5.4))
s = np.linspace(0.05, 2.2, 400)
ax.axhspan(1, 4, color="#fdedec", zorder=0); ax.axhline(1, color="#c0392b", lw=1.5)
ax.text(0.62, 3.15, "x(s) ≥ 1: past the wall\nthe word series 1/(1 − x) blows up", color="#922b21", fontsize=9.5)
ax.text(0.95, 0.72, "x(s) < 1: words add up to 1/(1 − x)", color="#1b4f72", fontsize=9.5)
for name, rs, c in SETS:
    xs = np.sum(np.array(rs)[:, None]**s[None, :], axis=0)
    ax.plot(s, xs, color=c, lw=2, label=name)
    lo, hi = 0.0, 5.0
    for _ in range(100):
        mid = (lo + hi)/2; lo, hi = (mid, hi) if sum(r**mid for r in rs) > 1 else (lo, mid)
    ax.plot([lo], [1], "o", color=c, ms=7, zorder=5)
ax.annotate("four different half boxes,\none crossing: D = 1/2", xy=(0.5, 1), xytext=(0.12, 2.3), arrowprops=dict(arrowstyle="->"), fontsize=9.5)
ax.set_ylim(0, 3.5); ax.set_xlabel("exponent s"); ax.set_ylabel("zoom sum  x(s) = Σ rᵢˢ")
ax.set_title("Figure 2. Dimension is an output: where the zoom sum hits the wall x = 1", weight="bold")
ax.legend(frameon=False, fontsize=9, loc="upper left", bbox_to_anchor=(0.66, 0.78))
fig.tight_layout(); fig.savefig(os.path.join(HERE, "fig2-dimension-at-wall.png"), dpi=130); plt.close(fig)

# ---------- Fig 3: complex dimensions (the wobble) are the other wall points
fig, ax = plt.subplots(figsize=(9, 5.6))
for name, r, N, c, off in [("A: spacing 4.53", 0.25, 2, "#c0392b", -0.03), ("B: spacing 2.86", 1/9, 3, "#e67e22", 0.0), ("E: spacing 2.27", 1/16, 4, "#8e44ad", 0.03)]:
    d = np.log(N)/np.log(1/r); m = np.arange(0, 20); im = 2*np.pi*m/np.log(1/r); im = im[im < 40]
    ax.plot(np.full_like(im, d + off), im, "o", color=c, label=f"{name} (evenly spaced: Minkowski wobble)", ms=6)
L1, L2 = np.log(0.09), np.log(0.49); found = set()
for im0 in np.arange(0.5, 40, 0.25):
    for re0 in [0.2, 0.4, 0.5]:
        z = complex(re0, im0)
        for _ in range(60): z = z - (np.exp(L1*z) + np.exp(L2*z) - 1)/(L1*np.exp(L1*z) + L2*np.exp(L2*z))
        if abs(np.exp(L1*z) + np.exp(L2*z) - 1) < 1e-12 and 0 <= z.imag < 40 and -2 < z.real < 1: found.add((round(z.real, 5), round(z.imag, 5)))
found = np.array(sorted(found) + [(0.5, 0.0)])
ax.plot(found[:, 0], found[:, 1], "D", color="#16a085", ms=6, label="D: uneven zooms (scattered: no wobble)")
ax.axvline(0.5, color="#555", ls=":", lw=1)
ax.set_xlabel("real part of s"); ax.set_ylabel("imaginary part of s")
ax.set_title("Figure 3. Complex dimensions: the other points where x(s) = 1", weight="bold")
ax.legend(frameon=False, fontsize=9, loc="upper center", bbox_to_anchor=(0.5, -0.12), ncol=2); ax.set_xlim(-0.6, 0.85); ax.set_ylim(-1, 41)
ax.text(0.52, 38, "real dimension 1/2", fontsize=9, color="#555")
fig.tight_layout(); fig.savefig(os.path.join(HERE, "fig3-complex-dimensions.png"), dpi=130); plt.close(fig)

# ---------- Fig 4: how the zoom dial plugs into what the framework already has
fig, ax = plt.subplots(figsize=(12, 3.9)); ax.axis("off"); ax.set_xlim(0, 12); ax.set_ylim(0, 4)
boxes = [(0.2, "zoom r on a piece", "new input"), (2.6, "weight r^s\n(a zoom is a twist,\nk = r^D)", "weight = twist (Base)"),
         (5.0, "pieces side by side\nadd: x(s) = Σ rᵢˢ", "k = 0 (plain add)"), (7.4, "zoom again and again:\nall words = 1/(1 − x)", "k = −1 (any number\nof times)"),
         (9.8, "dimension D: the wall\nx(D) = 1; wobble:\ncomplex wall points", "output")]
for xb, txt, tag in boxes:
    ax.add_patch(FancyBboxPatch((xb, 1.3), 2.0, 1.7, boxstyle="round,pad=0.05", fc="#eaf2f8", ec="#1b4f72"))
    ax.text(xb + 1.0, 2.15, txt, ha="center", va="center", fontsize=9.5)
    ax.text(xb + 1.0, 0.95, tag, ha="center", va="top", fontsize=8.5, color="#7b241c", style="italic")
for xb in [2.2, 4.6, 7.0, 9.4]:
    ax.add_patch(FancyArrowPatch((xb, 2.15), (xb + 0.4, 2.15), arrowstyle="-|>", mutation_scale=14, color="#1b4f72"))
ax.text(6, 3.6, "Figure 4. The zoom dial is one new input; everything after it is a dial the framework already has", ha="center", fontsize=12, weight="bold")
fig.savefig(os.path.join(HERE, "fig4-how-it-plugs-in.png"), dpi=130, bbox_inches="tight"); plt.close(fig)
print("figures written")
