"""Figures for the second pass (other affine transformations). Run: python3 figures_more.py"""
import os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon, Circle
from scipy.linalg import expm

OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "figures")
NAVY, MAG, BLUE, ORANGE, GREY, GREEN = "#0b1f4b", "#e0115f", "#1f6fd1", "#e8892b", "#8a93a8", "#2a9d6f"
plt.rcParams.update({"font.size": 11, "axes.titleweight": "bold", "axes.titlesize": 13,
                     "axes.edgecolor": "#c9cfdd", "axes.labelcolor": NAVY, "text.color": NAVY,
                     "xtick.color": NAVY, "ytick.color": NAVY})


def E0(k, d):
    return np.array([[0, -d], [1, k]], float) - k / 2 * np.eye(2)


def unit_circle(k, d, tmax=3.2, n=600):
    kap = d - k * k / 4
    ts = np.linspace(-tmax, tmax, n) if kap <= 0 else np.linspace(0, 2 * np.pi / np.sqrt(kap), n)
    pts = np.array([expm(t * E0(k, d)) @ [1, 0] for t in ts])
    return pts, ts


def style(ax, lim):
    ax.set_xlim(-lim, lim); ax.set_ylim(-lim, lim); ax.set_aspect("equal")
    ax.axhline(0, color="#e3e6ee", lw=0.8, zorder=0); ax.axvline(0, color="#e3e6ee", lw=0.8, zorder=0)
    ax.set_xticks([]); ax.set_yticks([])


# ---------- Figure 5: centro-affine, rotation/shear/squeeze, projective, reflection ----------
fig, axs = plt.subplots(1, 4, figsize=(18, 4.9))

# (a) centro-affine: unit circles N = 1, with the swept wedge det(z, z') from the centre
ax = axs[0]
for (k, d), col, lab in [((0.6, 1.0), BLUE, "Δ < 0: ellipse"), ((0.6, -0.5), MAG, "Δ > 0: hyperbola")]:
    pts, ts = unit_circle(k, d, tmax=1.6)
    ax.plot(pts[:, 0], pts[:, 1], color=col, lw=2.4, label=lab)
    if d < 0:  # the other branch
        ax.plot(-pts[:, 0], -pts[:, 1], color=col, lw=2.4)
    sel = pts[(ts >= 0) & (ts <= 0.55)]
    ax.add_patch(Polygon(np.vstack([[0, 0], sel]), closed=True, color=col, alpha=0.18, lw=0))
# κ = 0: N = (x + k y / 2)^2 = 1 is two parallel lines
k0 = 0.6
yy = np.linspace(-3, 3, 2)
for sgn in (1, -1):
    ax.plot(sgn - k0 / 2 * yy, yy, color=GREEN, lw=1.8, ls="--", label="Δ = 0: two lines" if sgn == 1 else None)
ax.plot(0, 0, "o", color=NAVY, ms=6)
style(ax, 2.6)
ax.legend(loc="lower left", fontsize=9, frameon=False)
ax.set_title("Centro-affine: |z| = 1 has curvature −Δ/4")
ax.set_xlabel("shaded: equal centro-affine arc, swept from the centre", fontsize=9, color=GREY)

# (b) the curvature dial picks the rotation: rotation, shear, squeeze on a square
ax = axs[1]
sq = np.array([[0.3, 0.3], [1.1, 0.3], [1.1, 1.1], [0.3, 1.1]])
offs = [(-2.7, 0), (0, 0), (2.7, 0)]
for (k, d, col, name), (ox, oy) in zip([(0, 1, BLUE, "Δ < 0\nrotation"), (2, 1, GREEN, "Δ = 0\nshear"),
                                        (0, -1, MAG, "Δ > 0\nsqueeze")], offs):
    for j, t in enumerate(np.linspace(0, 0.6, 4)):
        g = expm(t * E0(k, d))
        P = sq @ g.T
        ax.add_patch(Polygon(P * 0.8 + [ox - 0.5, oy - 0.6], closed=True, fill=j == 3,
                             facecolor=col, alpha=0.25 if j == 3 else 1, edgecolor=col, lw=1.6 if j else 1))
    ax.text(ox, -1.75, name, ha="center", va="top", fontsize=10, color=col, fontweight="bold")
ax.set_xlim(-4.0, 4.0); ax.set_ylim(-2.7, 2.9); ax.set_aspect("equal"); ax.set_xticks([]); ax.set_yticks([])
ax.set_title("Unit multiplication: three kinds of turn")
ax.set_xlabel("every step keeps area (det = 1)", fontsize=9, color=GREY)

# (c) projective: chord conics meet the line at infinity in 0, 1, 2 points (plane squeezed into a disk)
ax = axs[2]


def squeeze(x, y):
    r = np.hypot(x, y)
    f = 1 / (1 + r)
    return x * f, y * f


def chord(k, d, lo=-400, hi=400, n=200001):
    s = np.linspace(lo, hi, n)
    q = 1 + k * s + d * s**2
    x, y = s / q, s**2 / q
    jump = np.r_[False, np.sign(q[1:]) != np.sign(q[:-1])] | (np.abs(q) < 1e-6)
    x[jump] = np.nan; y[jump] = np.nan
    return x, y


ax.add_patch(Circle((0, 0), 1, fill=False, ec=NAVY, lw=1.4, ls="--"))
for (k, d), col, lab in [((0, 1.4), BLUE, "ellipse: 0 points"), ((0, 0), GREEN, "parabola: 1 point"),
                         ((0, -0.5), MAG, "hyperbola: 2 points")]:
    x, y = chord(k, d)
    X, Y = squeeze(x, y)
    ax.plot(X, Y, color=col, lw=2.3, label=lab)
ax.plot(0, 0, "o", color=NAVY, ms=5)
style(ax, 1.12)
ax.set_xlim(-2.45, 1.25)
ax.legend(loc="lower left", fontsize=8.5, frameon=False, bbox_to_anchor=(-0.03, -0.03))
ax.set_title("Projective: where the conic meets infinity")
ax.set_xlabel("dashed circle = line at infinity (opposite points are one point)", fontsize=9, color=GREY)

# (d) reflection: conjugation eps -> k - eps is a slanted mirror; half-turn about the wall on the line
ax = axs[3]
k, d = 1.0, 1.0
C = np.array([[1, k], [0, -1]])
shape = np.array([[0.4, 0.5], [1.3, 0.6], [0.9, 1.4]])
ax.add_patch(Polygon(shape, closed=True, facecolor=BLUE, alpha=0.3, edgecolor=BLUE, lw=1.8))
ax.add_patch(Polygon(shape @ C.T, closed=True, facecolor=ORANGE, alpha=0.3, edgecolor=ORANGE, lw=1.8))
for p in shape:
    q = C @ p
    ax.plot([p[0], q[0]], [p[1], q[1]], color=GREY, lw=0.9, ls=":")
pts, _ = unit_circle(k, d)
ax.plot(pts[:, 0], pts[:, 1], color=NAVY, lw=1.2, alpha=0.6)
ax.axhline(0, color=NAVY, lw=1.4)
ax.text(1.9, 0.08, "mirror", fontsize=9, color=NAVY, ha="right")
ax.set_xlabel("z → z̄ keeps |z|; the mirror's lines slant by k", fontsize=9, color=GREY)
ax.set_xlim(-2.4, 2.4); ax.set_ylim(-1.9, 1.9); ax.set_aspect("equal"); ax.set_xticks([]); ax.set_yticks([])
ax.set_title("Reflection: conjugation")

fig.tight_layout()
fig.savefig(os.path.join(OUT, "fig5-other-transformations.png"), dpi=150)
print("wrote fig5")
