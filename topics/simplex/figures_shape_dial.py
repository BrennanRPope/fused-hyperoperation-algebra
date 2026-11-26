"""Pictures for the shape-dial and temperature tests.  Run: python3 figures_shape_dial.py (writes figures/fig10..13)."""
import os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon, Circle, RegularPolygon
from shape_dial import *

OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "figures")
BLUE, ORANGE, GREEN, RED, GRAY, PURPLE = "#2a6fdb", "#e8892b", "#2e9e5b", "#c8423b", "#888888", "#7a4fc4"
plt.rcParams.update({"font.size": 11, "axes.titlesize": 12, "figure.dpi": 130})
CMAP = "RdBu"  # -1 red-ish (round), 0 white-blue? use custom below
from matplotlib.colors import LinearSegmentedColormap
CM = LinearSegmentedColormap.from_list("turn", [ORANGE, "#f6dcc0", "#ffffff", ], N=256)
CM = LinearSegmentedColormap.from_list("turn", [ORANGE, "#f7e3cf", "#dce8fa", BLUE], N=256)


def save(fig, name):
    fig.savefig(os.path.join(OUT, name), bbox_inches="tight")
    plt.close(fig)


def heat(ax, B, khat_fn, n=72, margin_frac=0.06, title=""):
    xs = np.linspace(-1.05, 1.05, n)
    X, Y = np.meshgrid(xs, xs)
    P = np.column_stack([X.ravel(), Y.ravel()])
    width = (B.h + np.roll(B.h, len(B.h) // 2)).min()
    ok = B.inside(P, margin_frac * width)
    K = np.full(len(P), np.nan)
    for s in range(0, ok.sum(), 400):
        idx = np.where(ok)[0][s:s + 400]
        K[idx] = khat_fn(P[idx])
    K = np.clip(K, -1, 0)
    im = ax.pcolormesh(X, Y, K.reshape(X.shape), cmap=CM, vmin=-1, vmax=0, shading="auto")
    bd = B.boundary
    ang = np.arctan2(bd[:, 1], bd[:, 0])
    o = np.argsort(ang)
    ax.add_patch(Polygon(bd[o], closed=True, fill=False, lw=1.4, color="k"))
    ax.set_aspect("equal"); ax.axis("off"); ax.set_xlim(-1.1, 1.1); ax.set_ylim(-0.62, 1.08)
    ax.set_title(title, fontsize=11)
    return im


# fig10: heat maps along the Minkowski path
fig, axs = plt.subplots(1, 5, figsize=(15, 3.6))
ts = [0, 0.25, 0.5, 0.75, 1.0]
for ax, t in zip(axs, ts):
    if t == 1.0:
        B = PolyBody(TRI_V)
        fn = lambda P: 3 * gauss_curvature(B.metric, P, 1.5e-3)
    else:
        B = minkowski(t)
        fn = lambda P, B=B: khat(B, P)
    kc = fn(np.array([[0.0, 0.0]]))[0]
    lab = {0: "round (relativity)", 1.0: "triangle (Bayes)"}.get(t, f"t = {t}")
    im = heat(ax, B, fn, title=f"{lab}\ncentre turn rate {kc:+.2f}")
cb = fig.colorbar(im, ax=axs, fraction=0.02, pad=0.01)
cb.set_label("turn per unit area  (K̂)")
cb.set_ticks([-1, -0.5, 0]); cb.set_ticklabels(["−1 full\nThomas turn", "−0.5", "0 no turn"])
fig.suptitle("Plane dial held fixed, only the shape slides from round to triangle: the leftover turn fades, and in between it is a pattern, not one number",
             fontsize=12, y=1.04)
save(fig, "fig10-shape-dial-turn-maps.png")

# fig11: centre value along three paths, temperature on a log axis
fig, (a1, a2) = plt.subplots(1, 2, figsize=(12.5, 4.4))
tt = np.linspace(0, 1, 21)
km = [khat(minkowski(t), np.array([[0.0, 0.0]]))[0] if t < 1 else 0.0 for t in tt]
ss = np.linspace(0.5, 1, 21)
kcut = [khat(cut_triangle(s), np.array([[0.0, 0.0]]))[0] if s < 1 else 0.0 for s in ss]
a1.plot(tt, km, "-o", ms=3.5, color=BLUE, label="blend: (1−t)·disk + t·triangle")
a1.plot((ss - 0.5) / 0.5, kcut, "-s", ms=3.5, color=PURPLE, label="trim: triangle cut by a disk, radius ½ → 1")
a1.axhline(-1, color=ORANGE, lw=1, ls=":"); a1.axhline(0, color=GRAY, lw=1, ls=":")
a1.set_xlabel("shape dial (0 = round, 1 = triangle)"); a1.set_ylabel("centre turn per unit area  K̂")
a1.set_ylim(-1.08, 0.08); a1.legend(loc="lower right", fontsize=9, frameon=False)
a1.set_title("Two different routes, same two ends")
bb = np.logspace(-3, 3, 31)
kt = np.array([khat(free_energy_body(b), np.array([[0.0, 0.0]]))[0] for b in bb])
a2.semilogx(bb, kt, "-o", ms=3.5, color=GREEN, label="cooling the free-energy body")
a2.semilogx(bb[bb > 20], -1.45 / bb[bb > 20], "--", color="k", lw=1, label="−1.45/β  (cold tail)")
a2.semilogx(bb[bb < 0.3], -1 + 0.35 * bb[bb < 0.3], ":", color="k", lw=1.2, label="−1 + 0.35β  (hot start)")
a2.axhline(-1, color=ORANGE, lw=1, ls=":"); a2.axhline(0, color=GRAY, lw=1, ls=":")
a2.set_xlabel("inverse temperature β  (hot → cold)"); a2.set_ylim(-1.08, 0.08)
a2.set_title("Temperature: the turn fades smoothly to zero, no jump")
a2.legend(loc="center left", fontsize=9, frameon=False)
save(fig, "fig11-turn-fades-along-three-paths.png")

# fig12: two-dial map
fig, ax = plt.subplots(figsize=(9.5, 6.2))
S = np.linspace(0, 1, 120); Dd = np.linspace(0, 1, 120)
kint = np.interp(S, tt, km)
Z = Dd[:, None] * np.abs(kint)[None, :]
im = ax.pcolormesh(S, Dd, Z, cmap="Oranges", shading="auto", vmin=0, vmax=1)
cb = fig.colorbar(im, ax=ax, fraction=0.04)
cb.set_label("turn per unit velocity-area at the centre  = |d| × |K̂|")
ax.set_xlabel("shape dial: round  →  triangle   (pattern of turning)")
ax.set_ylabel("plane dial  |d| = 1/c²  (relativity side, d < 0)")
ax.set_xlim(-0.02, 1.02); ax.set_ylim(-0.02, 1.02)
for (x, y, txt, ha) in [(0.02, 0.93, "Special relativity\nround ball, turns", "left"),
                        (0.02, 0.06, "Galileo\nround ball, no turn", "left"),
                        (0.98, 0.93, "Bayes / simplex\nhexagonal ball, no turn\nat any d", "right"),
                        (0.98, 0.06, "Bayes again\n(d does nothing here)", "right")]:
    ax.text(x, y, txt, ha=ha, va="center", fontsize=10, bbox=dict(fc="white", ec="none", alpha=0.8))
ax.annotate("", xy=(0.5, 0.12), xytext=(0.5, 0.62), arrowprops=dict(arrowstyle="<->", color="k"))
ax.text(0.52, 0.37, "|d| scales the turn\nuniformly, never\nchanges its pattern", fontsize=9, va="center")
ax.annotate("", xy=(0.85, 0.75), xytext=(0.15, 0.75), arrowprops=dict(arrowstyle="->", color=BLUE, lw=1.6))
ax.text(0.5, 0.69, "shape changes the pattern and ends where d cannot reach", color=BLUE, fontsize=9, ha="center")
# unit balls
iax = ax.inset_axes([0.13, 0.30, 0.13, 0.2]); iax.add_patch(Circle((0, 0), 1, fill=False, color=ORANGE, lw=2))
iax.set_xlim(-1.3, 1.3); iax.set_ylim(-1.3, 1.3); iax.set_aspect("equal"); iax.axis("off"); iax.set_title("unit ball", fontsize=8)
iax2 = ax.inset_axes([0.76, 0.30, 0.13, 0.2]); iax2.add_patch(RegularPolygon((0, 0), 6, radius=1, fill=False, color=BLUE, lw=2))
iax2.set_xlim(-1.3, 1.3); iax2.set_ylim(-1.3, 1.3); iax2.set_aspect("equal"); iax2.axis("off"); iax2.set_title("unit ball", fontsize=8)
ax.set_title("Two dials: |d| = how strongly it turns, shape = whether and where it turns")
save(fig, "fig12-two-dial-map.png")

# fig13: regular polygons and the square's off-centre pattern
fig, (a1, a2) = plt.subplots(1, 2, figsize=(11.5, 4.2), gridspec_kw={"width_ratios": [1.2, 1]})
ns = list(range(3, 13))
kc = []
for n in ns:
    a = np.radians(90) + TWO_PI * np.arange(n) / n
    V = np.column_stack([np.cos(a), np.sin(a)])
    kc.append(3 * gauss_curvature(PolyBody(V).metric, np.array([[0.0, 0.0]]), 1e-3)[0])
a1.bar(ns, kc, color=[BLUE] + [ORANGE] * (len(ns) - 1))
a1.axhline(0, color="k", lw=0.8)
a1.plot([3], [kc[0]], "o", color=BLUE, ms=9)
a1.annotate("triangle: 0", (3, 0), (3.3, -0.25), color=BLUE, arrowprops=dict(arrowstyle="->", color=BLUE))
a1.set_xticks(ns); a1.set_xlabel("number of sides"); a1.set_ylabel("turn rate at the centre  K̂")
a1.set_ylim(-1.1, 0.1)
a1.set_title("Only the triangle is flat at its centre;\nany more symmetry makes the centre look round")
Sq = PolyBody(np.array([[1, 1], [-1, 1], [-1, -1], [1, -1]], float) / np.sqrt(2))
heat(a2, Sq, lambda P: 3 * gauss_curvature(Sq.metric, P, 1.5e-3), title="square: round in the middle,\nflatter toward the edges")
a2.set_ylim(-1.05, 1.05)
save(fig, "fig13-polygons-only-triangle-flat.png")
print("done")
