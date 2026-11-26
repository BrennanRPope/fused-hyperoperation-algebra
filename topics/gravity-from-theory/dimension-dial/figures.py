"""DRAFT figures for dimension-dial/README.md. Run: python3 figures.py  (writes figures/*.png)"""
import math
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
from matplotlib.patches import FancyBboxPatch, Rectangle, Circle  # noqa: E402
from scipy.integrate import solve_ivp  # noqa: E402
from scipy.special import gamma as G, ive  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
FIG = os.path.join(HERE, "figures")
os.makedirs(FIG, exist_ok=True)
RED, GREEN, ORANGE, BLUE, PURPLE, GREY = "#d62728", "#2ca02c", "#ff7f0e", "#1f77b4", "#9467bd", "#7f7f7f"


def S(D):
    return 2 * np.pi**(D / 2) / G(D / 2)


# ---------------------------------------------------------------- fig1 the three pictures across the dial
def fig1():
    fig, ax = plt.subplots(figsize=(9, 5.5))
    Dc = np.linspace(1, 4.5, 300)
    ax.plot(Dc, Dc - 1, color=BLUE, lw=2.5, label="flux (Gauss): pull ~ 1/r^(D−1), any D")
    ax.plot(Dc, 0 * Dc + 2, color=RED, lw=2.5, label="velocity-line rule: pull ~ 1/r² at every D")
    ax.plot([3, 4], [2, 3], "s", color=GREEN, ms=12, label="Einstein's river (whole dimensions)")
    ax.plot([2], [0], "s", color=GREEN, ms=12, mfc="white", mew=2)
    ax.annotate("2+1: no pull outside matter", (2, 0), (2.15, 0.35), fontsize=9, color=GREEN,
                arrowprops=dict(arrowstyle="->", color=GREEN))
    ax.plot([3], [2], "o", ms=26, mfc="none", mec="black", mew=2)
    ax.text(3.08, 1.55, "all three agree:\nonly at D = 3", fontsize=10)
    ax.axvspan(4, 4.5, color=GREY, alpha=0.15)
    ax.text(4.03, 0.2, "no stable\norbits\n(D ≥ 4)", fontsize=8, color=GREY)
    ax.set_xlabel("dimension dial D (number of space dimensions)")
    ax.set_ylabel("power of the pull  (pull ~ 1/r^power)")
    ax.set_xlim(1, 4.5)
    ax.set_ylim(-0.3, 3.6)
    ax.set_xticks([1, 1.5, 2, 2.5, 3, 3.5, 4, 4.5])
    ax.legend(loc="upper left", fontsize=9)
    ax.set_title("Fig 1. Three pictures of a point source's pull, read across the dimension dial")
    fig.savefig(os.path.join(FIG, "fig1_three_pictures.png"), dpi=130, bbox_inches="tight")
    plt.close(fig)


# ---------------------------------------------------------------- fig2 right shape for the field
def fig2():
    fig, axs = plt.subplots(1, 3, figsize=(14, 4.8))
    # (a) field lines with sphere vs cube surface
    ax = axs[0]
    for ang in np.linspace(0, 2 * np.pi, 24, endpoint=False):
        ax.plot([0.08 * np.cos(ang), 1.6 * np.cos(ang)], [0.08 * np.sin(ang), 1.6 * np.sin(ang)], color=GREY, lw=0.8)
    ax.add_patch(Circle((0, 0), 1.0, fill=False, ec=GREEN, lw=2.5, label="sphere (p = 2)"))
    ax.add_patch(Rectangle((-1, -1), 2, 2, fill=False, ec=ORANGE, lw=2, ls="--", label="cube (p = ∞)"))
    ax.plot(0, 0, "ko")
    ax.set_aspect("equal")
    ax.set_xlim(-1.7, 1.7)
    ax.set_ylim(-1.7, 1.7)
    ax.axis("off")
    ax.legend(loc="lower center", fontsize=8, ncol=2, bbox_to_anchor=(0.5, -0.08))
    ax.set_title("Field lines cross the sphere at right angles,\nevenly; they cross the cube at a slant", fontsize=10)
    # (b) flux density on one face of the cube
    ax = axs[1]
    x = np.linspace(-1, 1, 201)
    X, Y = np.meshgrid(x, x)
    En = 1 / (4 * np.pi) / (X**2 + Y**2 + 1)**1.5
    im = ax.imshow(En / En.max(), extent=[-1, 1, -1, 1], origin="lower", cmap="viridis")
    fig.colorbar(im, ax=ax, fraction=0.046, label="flux density ÷ face-centre value")
    ax.set_title("Cube: flux through one face (3D)\ncentre is 3√3 ≈ 5.2× the corner (D4a)", fontsize=10)
    ax.set_xlabel("face x")
    ax.set_ylabel("face y")
    # (c) p dial: turned p-circles
    ax = axs[2]
    th = np.linspace(0, 2 * np.pi, 800)
    for p, col in [(1, PURPLE), (2, GREEN), (6, ORANGE)]:
        ux, uy = np.cos(th), np.sin(th)
        s = (np.abs(ux)**p + np.abs(uy)**p)**(1 / p)
        ux, uy = ux / s, uy / s
        ax.plot(ux, uy, color=col, lw=2, label=f"p = {p}" + (" (sphere)" if p == 2 else " (cross-polytope)" if p == 1 else " (toward cube)"))
        c, si = math.cos(0.4), math.sin(0.4)
        ax.plot(c * ux - si * uy, si * ux + c * uy, color=col, lw=1, ls=":")
    ax.set_aspect("equal")
    ax.axis("off")
    ax.legend(fontsize=8, loc="lower center", bbox_to_anchor=(0.5, -0.12), ncol=1)
    ax.set_title("p dial: solid = shape, dotted = turned by 0.4 rad\nonly the sphere lands on itself (D4b)", fontsize=10)
    fig.suptitle("Fig 2. The right shape for a point source's field is the sphere: its symmetry is turning", fontsize=12, y=1.03)
    fig.savefig(os.path.join(FIG, "fig2_right_shape_field.png"), dpi=130, bbox_inches="tight")
    plt.close(fig)


# ---------------------------------------------------------------- fig3 cube cells -> round field
def lattice_a2_grid(N=14):
    tt = np.geomspace(1e-3, 1e5, 6000)
    w = np.gradient(tt)
    I0 = ive(0, 2 * tt)
    xs = np.arange(-N, N + 1)
    cache = {k: ive(abs(k), 2 * tt) for k in xs}
    A = np.zeros((len(xs), len(xs)))
    for i, xi in enumerate(xs):
        for j, yj in enumerate(xs):
            A[j, i] = np.sum((I0 * I0 - cache[xi] * cache[yj]) * w)
    return xs, A


def fig3():
    xs, A = lattice_a2_grid()
    fig, axs = plt.subplots(1, 2, figsize=(12, 5.4))
    ax = axs[0]
    X, Y = np.meshgrid(xs, xs)
    ax.contour(X, Y, A, levels=[0.25, 0.32, 0.38, 0.45, 0.52, 0.6, 0.68], colors=BLUE)
    for k in xs:
        ax.axhline(k, color=GREY, lw=0.2)
        ax.axvline(k, color=GREY, lw=0.2)
    ax.plot(0, 0, "ko")
    ax.set_aspect("equal")
    ax.set_title("2D grid of square cells (the theory's cube in 2D)\npotential of one point: squarish near, round far", fontsize=10)
    ax.set_xlabel("cells")
    ax = axs[1]
    rr = np.hypot(X, Y).ravel()
    m = rr > 0
    ax.plot(rr[m], A.ravel()[m], ".", color=BLUE, ms=3, label="grid potential a(x) (all directions)")
    rc = np.linspace(0.8, 20, 200)
    ax.plot(rc, (np.log(rc) + 0.5772 + 1.5 * np.log(2)) / (2 * np.pi), color=GREEN, lw=2,
            label="round answer (ln r + γ + 1.5 ln 2)/(2π)")
    ax.set_xlabel("distance r (cells)")
    ax.set_ylabel("potential")
    ax.legend(fontsize=9)
    ax.set_title("Spread at each r is the cube's shape; it dies out with r\n(3D: 4πrG = 1.08 vs 0.95 at r ≈ 1.5, 1.002 vs 0.999 at r ≈ 12, D5a)", fontsize=10)
    fig.suptitle("Fig 3. The theory's own cells are cubes; the field they make turns into the sphere far away", fontsize=12, y=1.03)
    fig.savefig(os.path.join(FIG, "fig3_cube_cells_round_field.png"), dpi=130, bbox_inches="tight")
    plt.close(fig)


# ---------------------------------------------------------------- fig4 the Gaussian bridge and S(D)
def fig4():
    fig, axs = plt.subplots(1, 2, figsize=(12, 5))
    ax = axs[0]
    x = np.linspace(-2.5, 2.5, 300)
    X, Y = np.meshgrid(x, x)
    g1 = np.exp(-x**2)
    ax.contour(X, Y, np.outer(g1, g1), levels=6, colors=GREEN)
    ax.plot(x, -2.5 + 0.8 * g1, color=ORANGE, lw=2)
    ax.plot(-2.5 + 0.8 * g1, x, color=ORANGE, lw=2)
    ax.text(0.2, -2.3, "one axis: e^(−x²)", color=ORANGE, fontsize=9)
    ax.text(-1.0, 1.9, "product over axes\n= round contours", color=GREEN, fontsize=9)
    ax.set_aspect("equal")
    ax.set_title("A product over axes (like the box) that is also round:\nthe Gaussian, so it can be raised to any power D", fontsize=10)
    ax.axis("off")
    ax = axs[1]
    Dc = np.linspace(0.3, 6, 400)
    ax.plot(Dc, S(Dc), color=GREEN, lw=2.5, label="sphere area S(D) = 2π^(D/2)/Γ(D/2)")
    ax.plot(Dc, 2 * Dc * 2**(Dc - 1), color=ORANGE, lw=1.5, ls="--", label="cube of half-side 1: 2D·2^(D−1)")
    for k in [1, 2, 3, 4, 5]:
        ax.plot(k, S(k), "o", color=GREEN)
    ax.annotate("2π (circle)", (2, S(2)), (2.2, 3.8), fontsize=8)
    ax.annotate("4π (sphere)", (3, S(3)), (3.1, 10.5), fontsize=8)
    ax.set_xlabel("dimension dial D")
    ax.set_ylabel("surface of the unit shape")
    ax.set_ylim(0, 40)
    ax.legend(fontsize=9)
    ax.set_title("Both surfaces exist at every D (power r^(D−1) the same);\nthe pull's constant is 1/S(D) for the round field (D2b)", fontsize=10)
    fig.suptitle("Fig 4. Non-integer settings of the dial: the Gaussian carries the sphere off the whole numbers", fontsize=12, y=1.03)
    fig.savefig(os.path.join(FIG, "fig4_gaussian_bridge.png"), dpi=130, bbox_inches="tight")
    plt.close(fig)


# ---------------------------------------------------------------- fig5 orbits and velocity paths across the dial
def orbit(power, e=0.3, T=14.0):
    def f(tt, s):
        x, y, vx, vy = s
        r = math.hypot(x, y)
        a = -1 / r**power
        return [vx, vy, a * x / r, a * y / r]
    sol = solve_ivp(f, [0, T], [1, 0, 0, 1 - e], rtol=1e-10, atol=1e-12, max_step=0.01)
    return sol.y


def fig5():
    Ds = [2.0, 2.5, 3.0, 3.5]
    fig, axs = plt.subplots(2, 4, figsize=(15, 7.6))
    for j, Dv in enumerate(Ds):
        y = orbit(Dv - 1, T=14.0 if Dv != 3.5 else 9.0)
        col = GREEN if Dv == 3 else BLUE
        axs[0, j].plot(y[0], y[1], color=col, lw=1)
        axs[0, j].plot(0, 0, "o", color=ORANGE)
        axs[0, j].set_aspect("equal")
        axs[0, j].axis("off")
        axs[0, j].set_title(f"D = {Dv:g}: pull 1/r^{Dv - 1:g}\nnear-circular apsidal angle {180 / math.sqrt(4 - Dv):.0f}°", fontsize=10)
        axs[1, j].plot(y[2], y[3], color=col, lw=1)
        axs[1, j].plot(0, 0, "k+")
        axs[1, j].set_aspect("equal")
        axs[1, j].set_title("velocity path: " + ("one circle (straight on the dial)" if Dv == 3 else "not a circle"), fontsize=9)
        axs[1, j].set_xticks([])
        axs[1, j].set_yticks([])
    fig.suptitle("Fig 5. Orbits with the flux pull at four settings of the dimension dial (top) and their velocity paths (bottom).\n"
                 "Only D = 3 closes and keeps its velocity on one circle (D6a, D7a).", fontsize=12, y=1.02)
    fig.savefig(os.path.join(FIG, "fig5_orbits_across_dial.png"), dpi=130, bbox_inches="tight")
    plt.close(fig)


# ---------------------------------------------------------------- fig6 the 2D case: three answers
def fig6():
    fig, axs = plt.subplots(1, 3, figsize=(15, 5.2))
    y = orbit(1.0, T=30)
    axs[0].plot(y[0], y[1], color=BLUE, lw=0.8)
    axs[0].plot(0, 0, "o", color=ORANGE)
    axs[0].set_title("Flux in flatland: pull 1/r (log potential)\norbits rosette, never escape", fontsize=10)
    y = orbit(2.0, T=12)
    axs[1].plot(y[0], y[1], color=RED, lw=1.5)
    axs[1].plot(0, 0, "o", color=ORANGE)
    axs[1].set_title("Velocity-line rule in flatland: pull 1/r²\nclosed ellipses; velocity sphere S² (2D hydrogen)", fontsize=10)
    # cone: plane with a wedge removed, straight lines
    ax = axs[2]
    defi = np.deg2rad(70)
    th = np.linspace(defi / 2, 2 * np.pi - defi / 2, 200)
    ax.fill(np.r_[0, 1.6 * np.cos(th), 0], np.r_[0, 1.6 * np.sin(th), 0], color=GREEN, alpha=0.12)
    ax.plot([0, 1.6 * np.cos(defi / 2)], [0, 1.6 * np.sin(defi / 2)], color=GREEN, lw=1.5)
    ax.plot([0, 1.6 * np.cos(-defi / 2)], [0, 1.6 * np.sin(-defi / 2)], color=GREEN, lw=1.5)
    ax.text(1.05, -0.08, "wedge cut out;\nedges glued", fontsize=8, color=GREEN)
    hh = 0.45
    for sgn in (1, -1):
        xe = hh / math.tan(defi / 2)                       # where the path meets the wedge edge
        ax.plot([-1.5, xe], [sgn * hh, sgn * hh], color=PURPLE, lw=2)
        ax.plot(xe, sgn * hh, '*', color=RED, ms=16)
    ax.plot(0, 0, "o", color=ORANGE)
    ax.text(-1.5, 0.7, "two paths stay straight (no pull);\nthe two red stars are one point after\ngluing, so the paths cross there", fontsize=8, color=PURPLE)
    ax.set_title("Einstein in 2+1: river speed constant outside matter\nno pull; a mass only cuts a wedge (a cone)", fontsize=10)
    for a_ in axs:
        a_.set_aspect("equal")
        a_.axis("off")
    fig.suptitle("Fig 6. The 2D setting of the dial gives three different answers, so the pictures conflict there", fontsize=12, y=1.03)
    fig.savefig(os.path.join(FIG, "fig6_flatland_three_answers.png"), dpi=130, bbox_inches="tight")
    plt.close(fig)


# ---------------------------------------------------------------- fig7 shapes map
def box(ax, x, y, w, h, text, color, fs=9):
    ax.add_patch(FancyBboxPatch((x - w / 2, y - h / 2), w, h, boxstyle="round,pad=0.02", fc=color, ec="black", alpha=0.25))
    ax.text(x, y, text, ha="center", va="center", fontsize=fs)


def fig7():
    fig, ax = plt.subplots(figsize=(11, 6))
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.axis("off")
    rows = [
        ("Point source's field", "sphere (p = 2)", "its symmetry is turning; only the sphere is a level surface (D4)", GREEN),
        ("The theory's cells (lift, lattice)", "cube (p = ∞)", "the operation acts axis by axis; field turns round far away (D5)", ORANGE),
        ("Non-integer D", "Gaussian", "a product over axes that is round, so it takes any power D (D2)", BLUE),
        ("Velocity-line rule", "sphere one dimension up (S^D)", "the dial sphere in velocity space (D6c)", RED),
        ("Simplex, cross-polytope", "not used", "not invariant under turns; no point source has their symmetry", GREY),
    ]
    ax.text(0.15, 0.95, "job", weight="bold", fontsize=11, ha="center")
    ax.text(0.42, 0.95, "right shape", weight="bold", fontsize=11, ha="center")
    ax.text(0.76, 0.95, "why", weight="bold", fontsize=11, ha="center")
    for i, (job, shape, why, col) in enumerate(rows):
        yy = 0.82 - i * 0.16
        box(ax, 0.15, yy, 0.26, 0.11, job, col)
        box(ax, 0.42, yy, 0.22, 0.11, shape, col)
        ax.text(0.55, yy, why, va="center", fontsize=9)
    ax.set_title("Fig 7. The right shape for each job", fontsize=12)
    fig.savefig(os.path.join(FIG, "fig7_right_shapes.png"), dpi=130, bbox_inches="tight")
    plt.close(fig)


if __name__ == "__main__":
    for fn in (fig1, fig2, fig3, fig4, fig5, fig6, fig7):
        fn()
        print("wrote", fn.__name__)
