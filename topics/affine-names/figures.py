"""Figures for the affine names of the dials. Run: python3 figures.py"""
import os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch, Polygon

OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "figures")
os.makedirs(OUT, exist_ok=True)
NAVY, MAG, BLUE, ORANGE, GREY, GREEN = "#0b1f4b", "#e0115f", "#1f6fd1", "#e8892b", "#8a93a8", "#2a9d6f"
plt.rcParams.update({"font.size": 11, "axes.titleweight": "bold", "axes.titlesize": 13,
                     "axes.edgecolor": "#c9cfdd", "axes.labelcolor": NAVY, "text.color": NAVY,
                     "xtick.color": NAVY, "ytick.color": NAVY})


def lp_curve(p, n=800):
    th = np.linspace(0, 2 * np.pi, n)
    c, s = np.cos(th), np.sin(th)
    r = (np.abs(c) ** p + np.abs(s) ** p) ** (-1 / p) if np.isfinite(p) else 1 / np.maximum(np.abs(c), np.abs(s))
    return r * c, r * s


def conic(ax, k, d, color, lw=2.2, label=None):
    s = np.linspace(-60, 60, 40001)
    q = 1 + k * s + d * s**2
    x, y = s / q, s**2 / q
    bad = np.abs(q) < 2e-3
    x[bad] = np.nan; y[bad] = np.nan
    # break at sign changes of q (asymptotes)
    jump = np.r_[False, np.sign(q[1:]) != np.sign(q[:-1])]
    x[jump] = np.nan
    ax.plot(x, y, color=color, lw=lw, label=label)
    ax.plot(0, 0, "o", color=NAVY, ms=5)


# ---------- Figure 1: the four dials at a glance ----------
fig, axs = plt.subplots(1, 4, figsize=(17, 4.7))

ax = axs[0]
ax.set_title("k  ·  homothety")
cen = np.array([-1.0, 0.0])
tri = np.array([[0.0, 0.0], [0.35, 0.45], [0.55, -0.1]])
ax.add_patch(Polygon(tri, closed=True, fc=BLUE, alpha=0.25, ec=BLUE, lw=1.8))
for lam, col in [(1.6, ORANGE), (2.3, MAG)]:
    img = cen + lam * (tri - cen)
    ax.add_patch(Polygon(img, closed=True, fc=col, alpha=0.15, ec=col, lw=1.8))
for v in tri:
    far = cen + 2.45 * (v - cen)
    ax.plot([cen[0], far[0]], [cen[1], far[1]], color=GREY, lw=0.8, ls="--")
ax.plot(*cen, "s", color=NAVY, ms=9)
ax.text(cen[0], cen[1] - 0.22, "centre = the wall −1/k", ha="center", fontsize=9.5, color=NAVY)
ax.text(0.15, 0.62, "1 + ka = 1\n(identity)", color=BLUE, fontsize=9, ha="center")
ax.text(1.1, 1.1, "1 + ka = 1.6", color=ORANGE, fontsize=9)
ax.text(1.85, 1.55, "2.3", color=MAG, fontsize=9)
ax.text(0.6, -1.15, "a ∘ₖ b = doing scaling a, then scaling b.\nk = 0 moves the centre to infinity:\nevery scaling becomes a slide (addition)",
        ha="center", fontsize=9.5)
ax.set_xlim(-1.5, 2.6); ax.set_ylim(-1.55, 1.85); ax.set_aspect("equal"); ax.axis("off")

ax = axs[1]
ax.set_title("d  ·  curvature dial (−Δ/4)")
conic(ax, 0, 1.0, BLUE, label="Δ < 0  ellipse")
conic(ax, 0, 0.0, MAG, label="Δ = 0  parabola")
conic(ax, 0, -1.0, ORANGE, label="Δ > 0  hyperbola")
ax.set_xlim(-2.2, 2.2); ax.set_ylim(-1.6, 2.6); ax.set_aspect("equal"); ax.axis("off")
ax.legend(loc="upper center", fontsize=9, frameon=False, ncol=1, bbox_to_anchor=(0.5, -0.02))

ax = axs[2]
ax.set_title("p  ·  ℓp exponent")
cols = [ORANGE, "#f0a35c", BLUE, "#7aa7e0", GREEN]
for p, col in zip([1, 1.4, 2, 4, np.inf], cols):
    X, Y = lp_curve(p)
    ax.plot(X, Y, color=col, lw=2, label=f"p = {'∞' if not np.isfinite(p) else p:g}" if np.isfinite(p) else "p = ∞")
ax.set_aspect("equal"); ax.axis("off")
ax.legend(loc="center left", bbox_to_anchor=(0.98, 0.5), fontsize=9, frameon=False)
ax.text(0, -1.45, "cross-polytope → ball → cube", ha="center", fontsize=9.5)
ax.set_ylim(-1.6, 1.15)

ax = axs[3]
ax.set_title("C  ·  cubic form")
th = np.linspace(0, 2 * np.pi, 600)
for t, col in zip([0.0, 0.5, 1.0], [BLUE, "#8b6fd1", MAG]):
    tri = np.array([[np.cos(a), np.sin(a)] for a in np.pi / 2 + np.arange(3) * 2 * np.pi / 3])
    # support-function blend of disk and triangle
    u = np.c_[np.cos(th), np.sin(th)]
    h_disk = 0.5 * np.ones_like(th)
    h_tri = (u @ tri.T).max(axis=1)
    h = (1 - t) * h_disk + t * h_tri
    # boundary from support function: x = h u + h' u_perp
    hp = np.gradient(h, th)
    X = h * u[:, 0] - hp * u[:, 1]
    Y = h * u[:, 1] + hp * u[:, 0]
    ax.plot(X, Y, color=col, lw=2, label=["C = 0  round", "in between", "triangle, C largest"][int(t * 2)])
ax.set_aspect("equal"); ax.axis("off")
ax.legend(loc="upper center", fontsize=9, frameon=False, bbox_to_anchor=(0.5, -0.02))
fig.suptitle("The four dials, named from affine geometry", fontsize=15, fontweight="bold", color=NAVY, y=1.02)
fig.tight_layout()
fig.savefig(os.path.join(OUT, "fig1-four-dials.png"), dpi=150, bbox_inches="tight", facecolor="white")
plt.close(fig)

# ---------- Figure 2: k and kappa are the two halves of one move ----------
fig, axs = plt.subplots(1, 2, figsize=(13, 5.4))
ax = axs[0]
K = np.linspace(-3, 3, 400)
for kap, col in [(-1.5, ORANGE), (-0.5, ORANGE), (0, MAG), (0.5, BLUE), (1.5, BLUE)]:
    ax.plot(K, kap + K**2 / 4, color=col, lw=2 if kap == 0 else 1.2, ls="-" if kap == 0 else "--")
    ax.text(2.2, kap + 2.2**2 / 4 + 0.12, f"−Δ/4 = {kap:g}", color=col, fontsize=9, bbox=dict(fc="white", ec="none", pad=0.6, alpha=0.9))
ax.fill_between(K, K**2 / 4, 4, color=BLUE, alpha=0.06)
ax.fill_between(K, -3, K**2 / 4, color=ORANGE, alpha=0.06)
ax.axhline(0, color=GREY, lw=0.8); ax.axvline(0, color=GREY, lw=0.8)
ax.plot([0, 1.5, -1.0], [0, 0, 0], "o", color=NAVY)
ax.annotate("addition", (0, 0), (0.15, -0.55), fontsize=9)
ax.annotate("∘₁ (k = 1, d = 0)\nΔ = 1, a hyperbola", (1, 0), (0.6, -1.9), fontsize=9, arrowprops=dict(arrowstyle="-", color=GREY))
ax.plot([0, 0], [1, -1], "s", color=NAVY, ms=6)
ax.annotate("circle  d = 1", (0, 1), (-2.9, 1.45), fontsize=9, arrowprops=dict(arrowstyle="-", color=GREY))
ax.annotate("relativity  d = −1/c²", (0, -1), (-2.9, -2.2), fontsize=9, arrowprops=dict(arrowstyle="-", color=GREY))
ax.set_xlabel("k  (homothety)"); ax.set_ylabel("d  (Coleman–McGuinness's second number)")
ax.set_xlim(-3, 3); ax.set_ylim(-2.6, 3.6)
ax.set_title("Equal curvature −Δ/4 = d − k²/4: parabolas")

ax = axs[1]
ax.set_title("Multiplying by ε splits into a scaling and an area-keeping turn")
th = np.linspace(0, 2 * np.pi, 400)
circle = np.c_[np.cos(th), np.sin(th)]
ax.plot(circle[:, 0], circle[:, 1], color=GREY, lw=1, ls=":")
for kk, col, lab in [(0.0, BLUE, "k = 0: pure turn (d = 1)"), (0.6, MAG, "k = 0.6: turn + homothety")]:
    t = np.linspace(0, 2 * np.pi, 600)
    r = np.exp(kk / 2 * t / (2 * np.pi) * 2)
    ax.plot(r * np.cos(t), r * np.sin(t), color=col, lw=2, label=lab)
ax.set_aspect("equal"); ax.axis("off")
ax.text(0, -2.55, "ε = (k/2)·I  +  ε₀,   ε₀² = Δ/4\nhomothety part   +   curvature part", ha="center", fontsize=11, color=NAVY)
ax.legend(loc="upper right", fontsize=9, frameon=False)
fig.tight_layout()
fig.savefig(os.path.join(OUT, "fig2-homothety-plus-curvature.png"), dpi=150, bbox_inches="tight", facecolor="white")
plt.close(fig)

# ---------- Figure 3: three-lobed versus four-lobed ----------
fig, axs = plt.subplots(1, 3, figsize=(15, 5))
th = np.linspace(0, 2 * np.pi, 800)
ax = axs[0]
r = 1 + 0.25 * np.cos(3 * th)
ax.plot(r * np.cos(th), r * np.sin(th), color=MAG, lw=2.4)
ax.plot(np.cos(th), np.sin(th), color=GREY, lw=1, ls=":")
ax.set_title("Cubic form C: three lobes (cos 3θ)")
ax.set_aspect("equal"); ax.axis("off")
ax.text(0, -1.5, "the only pattern a trace-free cubic form\ncan make in the plane (check P1)\n= the trefoil of optics", ha="center", fontsize=9.5)
ax = axs[1]
p = 3.0
X, Y = lp_curve(p)
ax.plot(X, Y, color=BLUE, lw=2.4)
ax.plot(np.cos(th), np.sin(th), color=GREY, lw=1, ls=":")
ax.set_title("ℓp ball: four lobes (cos 4θ)")
ax.set_aspect("equal"); ax.axis("off")
ax.text(0, -1.5, "first non-round harmonic is 4, never 3\n(check L3); its cubic form is 0 at the centre", ha="center", fontsize=9.5)
ax = axs[2]
ps = np.r_[np.linspace(1.02, 2, 30), np.linspace(2, 30, 60)]
import importlib.util
spec = importlib.util.spec_from_file_location("chk_lp", os.path.join(os.path.dirname(os.path.abspath(__file__)), "lp_ratio.py"))
mod = importlib.util.module_from_spec(spec); spec.loader.exec_module(mod)
vals = [mod.lp_ratio(q) for q in ps]
ax.plot(ps, vals, color=NAVY, lw=2.2)
ax.set_xscale("log")
ax.axvline(2, color=BLUE, ls="--", lw=1)
ax.set_xlabel("p  (log scale)"); ax.set_ylabel("Ω³ / (8π²·area)")
ax.set_title("Blaschke's affine roundness along p")
ax.text(2.1, 0.08, "p = 2: ellipse, ratio 1", color=BLUE, fontsize=9.5)
for a, b in [(1.5, 3.0), (1.25, 5.0)]:
    ax.plot([a, b], [mod.lp_ratio(a), mod.lp_ratio(b)], "o", color=ORANGE)
ax.text(4.5, 0.62, "orange pairs: a ball and its polar\npartner p/(p−1) have different\naffine perimeter", fontsize=8.5, color=ORANGE)
ax.set_ylim(0, 1.05)
fig.tight_layout()
fig.savefig(os.path.join(OUT, "fig3-three-lobes-four-lobes.png"), dpi=150, bbox_inches="tight", facecolor="white")
plt.close(fig)

# ---------- Figure 4: lineage timeline ----------
rows = [
    ("k  homothety", BLUE, [(1827, "Chasles:\nhomothetic conics"), (1957, "Artin:\ndilatations"), (1974, "Sugeno λ"),
                            (1988, "Tsallis\nentropy"), (1991, "Coleman–\nMcGuinness c"), (2004, "q-sum, k = 1 − q:\nBorges; Nivanen et al.")]),
    ("d  curvature dial", MAG, [(1748, "Euler:\naffine curves"), (1841, "Transon:\naffine normal"), (1859, "Cayley"),
                                  (1871, "Klein"), (1923, "Blaschke:\nconics = constant\ncurvature"), (1991, "Coleman–\nMcGuinness d")]),
    ("p  ℓp exponent", GREEN, [(1818, "Lamé\ncurves"), (1896, "Minkowski:\nnorm = convex body"), (1910, "Riesz:\nLp spaces"),
                               (1923, "Blaschke:\naffine perimeter"), (1962, "Firey:\np-sums")]),
    ("C  cubic form", ORANGE, [(1908, "Țițeica:\naffine spheres"), (1916, "Fubini–Pick\n(1910s)"), (1923, "Blaschke:\napolarity"),
                               (1934, "Zernike:\ntrefoil (optics)"), (1960, "Benzécri:\nshape space"), (1972, "Calabi"),
                               (1986, "Cheng–Yau"), (2004, "Loftin 2001,\nLabourie 2007")]),
]
fig, ax = plt.subplots(figsize=(17, 6.6))
for i, (name, col, ev) in enumerate(rows):
    y = len(rows) - i
    ax.plot([1740, 2012], [y, y], color=col, lw=2, alpha=0.4)
    ax.text(1736, y, name, ha="right", va="center", fontsize=12, fontweight="bold", color=col)
    for j, (yr, lab) in enumerate(ev):
        ax.plot(yr, y, "o", color=col, ms=8)
        up = j % 2 == 0
        ha = "center"
        if lab.startswith("q-sum"):
            up, ha = True, "left"
        if lab.startswith("Blaschke:\napolarity"):
            ha = "left"
        if lab.startswith("Tsallis") or lab.startswith("Țițeica"):
            ha = "right"
        ax.text(yr + (-1 if ha == "right" else 0), y + (0.18 if up else -0.18), lab, ha=ha, va="bottom" if up else "top", fontsize=8.5)
ax.set_xlim(1655, 2045); ax.set_ylim(0.4, len(rows) + 0.7)
ax.set_yticks([]); ax.set_xticks(range(1750, 2021, 25))
for s in ["left", "right", "top"]:
    ax.spines[s].set_visible(False)
ax.set_title("Who found each idea, and when", fontsize=14)
fig.tight_layout()
fig.savefig(os.path.join(OUT, "fig4-lineage.png"), dpi=150, bbox_inches="tight", facecolor="white")
plt.close(fig)
print("ok")
