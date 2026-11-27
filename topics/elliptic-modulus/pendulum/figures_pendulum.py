"""Pictures for the pendulum section (Nov 27): fig6 phase plane and energy line, fig7 three motions, fig8 complex time.
Run from this folder: python3 figures_pendulum.py -> ../figures/fig6..fig8. omega = 1 throughout."""
import numpy as np, math, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.special import ellipj, ellipk
plt.rcParams.update({"font.size": 10, "axes.spines.top": False, "axes.spines.right": False})
BLUE, ORANGE, GREEN, GREY, RED, PURPLE = "#2f6fb5", "#e07b1f", "#2a9d55", "#9a9a9a", "#c8382f", "#7a4fb5"
cmap = matplotlib.colors.LinearSegmentedColormap.from_list("bo", [BLUE, "#7a7a7a", ORANGE])
def col(m): return cmap(min(m, 1.0)) if m <= 1 else PURPLE
def save(fig, name): fig.savefig(f"../figures/{name}.png", dpi=150, bbox_inches="tight"); plt.close(fig)

def G(u, m):
    """Slope reading sn dn / cn, real u, any real m (m > 1 by the reciprocal-modulus transformation)."""
    u = np.asarray(u, float)
    if m == 1: return np.tanh(u)
    if m < 1 and m >= 0:
        s, c, d, _ = ellipj(u, m); return s*d/c
    if m > 1:
        r = math.sqrt(m); s, c, d, _ = ellipj(r*u, 1/m); return (s/r)*c/d
    raise ValueError

# ---------------- fig6: the phase plane with energy (kappa) and the slope reading drawn on it
fig = plt.figure(figsize=(14.5, 7.4)); gs = fig.add_gridspec(2, 1, height_ratios=[4.2, 1.0], hspace=0.32)
ax = fig.add_subplot(gs[0])
th = np.linspace(-np.pi, np.pi, 1200)
# slope-reading curves t = sin th / th'  ->  th' = sin th / t
for t in [0.25, 0.5, 2, 4]:
    for sgn in (1, -1):
        ax.plot(th, sgn*np.sin(th)/t, color=GREEN, lw=0.8, ls="--", alpha=0.55)
for sgn in (1, -1): ax.plot(th, sgn*np.sin(th), color=GREEN, lw=2.4)
ax.axhline(0, color=RED, lw=2.0)
# energy contours m = th'^2/4 + sin^2(th/2)
for m, lab in [(0.067, "30°"), (0.25, "60°"), (0.5, "90°"), (0.75, "120°"), (0.933, "150°"), (1.0, "180°"), (1.35, ""), (2.2, "")]:
    v = 2*np.sqrt(np.clip(m - np.sin(th/2)**2, 0, None)); ok = m - np.sin(th/2)**2 >= 0
    lw = 3 if m in (0.5, 1.0) else 1.8
    if m <= 1:
        ax.plot(th[ok], v[ok], color=col(m), lw=lw); ax.plot(th[ok], -v[ok], color=col(m), lw=lw)
    else:
        ax.plot(th, v, color=PURPLE, lw=lw); ax.plot(th, -v, color=PURPLE, lw=lw)
# halfway points: cos th* = 2 cos(th0/2) - 1
for m in [0.067, 0.25, 0.5, 0.75, 0.933]:
    th0 = 2*math.asin(math.sqrt(m)); ts = math.acos(2*math.cos(th0/2) - 1)
    for sx in (1, -1):
        for sy in (1, -1):
            ax.plot(sx*ts, sy*sx*math.sin(ts), "o", color="white", mec=GREEN, mew=2, ms=7, zorder=5)
ax.set_xlim(-np.pi, np.pi); ax.set_ylim(-3.05, 3.05)
ax.set_xticks([-np.pi, -np.pi/2, 0, np.pi/2, np.pi]); ax.set_xticklabels(["−180°", "−90°", "0", "90°", "180° (top)"])
ax.set_xlabel("angle θ"); ax.set_ylabel("angular speed θ′  (ω = 1)")
ax.text(-0.9, 2.62, "over the top: the reading never blows up, it is capped below 1", color=PURPLE, fontsize=9.5, bbox=dict(fc="white", ec="none", alpha=0.85, pad=1.5))
ax.text(-2.75, 1.75, "separatrix (κ = −1): relativity's wall", color=ORANGE, fontsize=9.5, bbox=dict(fc="white", ec="none", alpha=0.85, pad=1.5))
ax.text(-3.0, -0.2, "the wall: θ′ = 0, the pendulum turns back, t = ∞", color=RED, fontsize=9.5, va="top", bbox=dict(fc="white", ec="none", alpha=0.85, pad=1.5))
ax.annotate("t = 1, reached at exactly half the time\nfrom the bottom to the turning point", xy=(1.571, 1.0), xytext=(1.75, 1.62),
            color=GREEN, fontsize=9.5, bbox=dict(fc="white", ec="none", alpha=0.85, pad=1.5), arrowprops=dict(arrowstyle="->", color=GREEN))
ax.text(-1.55, -2.62, "dashed: other values of the reading t = sin θ / θ′ (0.25, 0.5, 2, 4)", color=GREEN, fontsize=9)
ax.text(0.06, 1.47, "90° swing: κ = 0", color="#555", fontsize=9, bbox=dict(fc="white", ec="none", alpha=0.85, pad=1.5))
ax.set_title("(a) Every pendulum on one sheet. Colour is the energy, which is our curvature reading κ = 1 − 2E/E_top.\n"
             "Green is the slope reading t = sin θ/θ′. White dots: t = 1, at half the time to the turning point", fontsize=10.5)
# energy line
ax = fig.add_subplot(gs[1]); ax.set_xlim(-2.6, 3.2); ax.set_ylim(-1.2, 1.3); ax.axis("off")
ax.plot([-2.5, 3.1], [0, 0], color="k", lw=1)
ax.add_patch(plt.Rectangle((-2.5, -0.18), 1.5, 0.36, color=PURPLE, alpha=0.35))
for x0 in np.linspace(-1, 1, 41): ax.add_patch(plt.Rectangle((x0-0.025, -0.18), 0.05, 0.36, color=cmap((1-x0)/2), alpha=0.6, lw=0))
ax.add_patch(plt.Rectangle((1, -0.18), 2.1, 0.36, color=GREY, alpha=0.25))
for x, lab in [(-2, "−2"), (-1, "−1"), (0, "0"), (1, "+1"), (2, "+2")]:
    ax.plot([x, x], [-0.25, 0.25], color="k", lw=1); ax.text(x, -0.48, "κ = " + lab, ha="center", fontsize=9)
ax.text(-1.75, 0.38, "over the top: capped, no wall", color=PURPLE, ha="center", fontsize=9.5)
ax.text(0, 0.38, "swinging: a wall at every turning point", color="#555", ha="center", fontsize=9.5)
ax.text(2.05, 0.38, "the Wick image of going over the top\n(imaginary time); walls again", color="#666", ha="center", fontsize=9)
ax.text(-1, -0.95, "separatrix, tanh", color=ORANGE, ha="center", fontsize=9); ax.text(1, -0.95, "tiny swing, tan", color=BLUE, ha="center", fontsize=9)
ax.text(0, -0.95, "90° swing, square torus", color="#555", ha="center", fontsize=9)
ax.text(-2.55, 1.05, "(b) The whole energy line. Energy runs right to left: κ = 1 − 2E/E_top.", fontsize=10.5)
save(fig, "fig6-phase-plane")

# ---------------- fig7: three motions, the reach, the mirror
fig, axs = plt.subplots(1, 3, figsize=(16, 4.6), gridspec_kw=dict(width_ratios=[1.45, 1, 1], wspace=0.45))
ax = axs[0]; T = 9
m = 0.5; K = ellipk(m)
for j in range(-1, 5):
    a, b = (2*j-1)*K, (2*j+1)*K; uu = np.linspace(max(a, 0)+1e-3, min(b, T)-1e-3, 400)
    if len(uu) and uu[0] < uu[-1]: ax.plot(uu, G(uu, m), color="#7a7a7a", lw=1.8, label="swing 90° (m = 0.5)" if j == 0 else None)
    if 0 < b < T: ax.axvline(b, color=RED, ls=":", lw=0.9)
uu = np.linspace(0, T, 800)
ax.plot(uu, np.tanh(uu), color=ORANGE, lw=2.2, label="separatrix (m = 1): tanh")
for m, ls in [(1.3, "-"), (2.5, "--")]:
    ax.plot(uu, G(uu, m), color=PURPLE, lw=1.8, ls=ls, label=f"over the top (m = {m})")
ax.axhline(1, color=GREEN, lw=1.2); ax.text(T-0.05, 1.08, "t = 1", color=GREEN, ha="right")
ax.plot([K/2], [1], "o", color="white", mec=GREEN, mew=2, ms=7)
ax.set_ylim(-3, 3); ax.set_xlim(0, T); ax.set_xlabel("time (ω = 1)"); ax.set_ylabel("slope reading t")
ax.legend(loc="lower right", fontsize=8.5, frameon=False)
ax.set_title("(a) Three motions. Swinging hits a wall (red dotted) at each turn;\nthe separatrix creeps up to 1; over the top it is capped below 1", fontsize=10)
ax = axs[1]
ml = np.linspace(0, 0.9995, 500); ax.plot(ml, ellipk(ml), color=GREEN, lw=2, label="time to the wall, K(m)")
ax.plot(ml, ellipk(ml)/2, color=GREEN, lw=1.2, ls="--", label="time to t = 1, K(m)/2")
ax.set_ylim(0, 6); ax.set_xlabel("energy m = E/E_top"); ax.set_ylabel("time (ω = 1)", color=GREEN)
ax2 = ax.twinx(); ax2.spines["right"].set_visible(True); mr = np.linspace(1.0001, 4, 400)
ax2.plot(mr, np.sqrt(mr) - np.sqrt(mr-1), color=PURPLE, lw=2); ax2.set_ylim(0, 1.2); ax2.set_ylabel("cap on the reading, √m − √(m−1)", color=PURPLE)
ax.axvline(1, color=ORANGE, lw=1.5, ls=":"); ax.text(1.07, 3.3, "separatrix", color=ORANGE, fontsize=9)
ax.set_xlim(0, 4); ax.legend(loc="upper left", fontsize=8.5, frameon=False)
ax.set_title("(b) Below the separatrix the wall is reached\nin finite time; above it the reading is capped", fontsize=10)
ax = axs[2]
for m in [0.05, 0.3, 0.5, 0.7, 0.9, 0.99]:
    K = ellipk(m); x = np.linspace(0.01, 0.99, 300); ax.plot(x, np.log10(G(x*K, m)), color=col(m), lw=1.6)
ax.axhline(0, color=GREEN, lw=1); ax.axvline(0.5, color=GREEN, lw=1, ls="--"); ax.plot([0.5], [0], "o", color="white", mec=GREEN, mew=2, ms=7)
ax.set_xlabel("time as a fraction of the time to the wall, u/K"); ax.set_ylabel("log₁₀ t")
ax.set_ylim(-2.2, 2.2)
ax.text(0.53, 1.6, "G(u) · G(K − u) = 1:\neach curve is point-symmetric\nabout (½, t = 1)", color=GREEN, fontsize=9)
ax.set_title("(c) The mirror: time from the bottom and time\nleft to the wall give readings t and 1/t", fontsize=10)
save(fig, "fig7-three-motions")

# ---------------- fig8: walls and zeros in complex time
def lattice(m, box):
    """zeros and poles of G in a box, from the closed forms."""
    xr, yr = box; Z, P = [], []
    if m < 1:
        K, Kp = ellipk(m), ellipk(1-m)
        for j in range(-6, 7):
            for k in range(-6, 7):
                Z += [(2*j*K, 2*k*Kp), ((2*j+1)*K, (2*k+1)*Kp)]
                P += [((2*j+1)*K, 2*k*Kp), (2*j*K, (2*k+1)*Kp)]
    else:
        r = math.sqrt(m); K, Kp = ellipk(1/m)/r, ellipk(1-1/m)/r
        for j in range(-8, 9):
            for k in range(-6, 7):
                Z += [(2*j*K, 2*k*Kp), ((2*j+1)*K, 2*k*Kp)]
                P += [((2*j+1)*K, (2*k+1)*Kp), (2*j*K, (2*k+1)*Kp)]
    inb = lambda p: xr[0] <= p[0] <= xr[1] and yr[0] <= p[1] <= yr[1]
    return [p for p in Z if inb(p)], [p for p in P if inb(p)]
fig, axs = plt.subplots(1, 4, figsize=(16, 4.5), sharey=True, gridspec_kw=dict(wspace=0.08))
box = ((-4.2, 4.2), (-4.2, 4.2))
for ax, m, ttl in zip(axs, [0.3, 0.9, 1.1, 2.5], ["swing 66°", "swing 143°", "over the top", "over the top"]):
    Z, P = lattice(m, box)
    ax.axhline(0, color=col(m), lw=3, alpha=0.8)
    if Z: ax.plot(*zip(*Z), "o", color=GREEN, ms=5, label="t = 0")
    if P: ax.plot(*zip(*P), "x", color=RED, ms=7, mew=2, label="wall, t = ∞")
    if m < 1:
        K = ellipk(m)
        for c in [K/2, -K/2, 3*K/2, -3*K/2]: ax.axvline(c, color=GREEN, ls="--", lw=1)
        ax.text(K/2+0.08, 3.6, "|t| = 1", color=GREEN, fontsize=9)
    ax.set_xlim(*box[0]); ax.set_ylim(*box[1]); ax.set_aspect("equal")
    ax.set_title(f"{ttl} (m = {m})", fontsize=10, color=col(m)); ax.set_xlabel("real time")
axs[0].set_ylabel("imaginary time"); axs[0].legend(loc="lower left", fontsize=8.5, frameon=True)
fig.suptitle("Walls and zeros of the slope reading in complex time. Swinging: the real time line (coloured) runs straight into walls. "
             "Over the top: it slips between two rows of walls.\nDashed green: the halfway lines, where the reading has size exactly 1",
             fontsize=10.5, y=1.04)
save(fig, "fig8-complex-time")
print("done")
