"""DRAFT figures for river-rebuild/README.md. Run: python3 figures.py  (writes figures/*.png)"""
import math
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
from matplotlib.patches import FancyBboxPatch  # noqa: E402
from scipy.integrate import solve_ivp  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
FIG = os.path.join(HERE, "figures")
os.makedirs(FIG, exist_ok=True)
RED, GREEN, ORANGE, BLUE, GREY = "#d62728", "#2ca02c", "#ff7f0e", "#1f77b4", "#7f7f7f"


def box(ax, x, y, w, h, text, color, fs=10):
    ax.add_patch(FancyBboxPatch((x - w / 2, y - h / 2), w, h, boxstyle="round,pad=0.02",
                                fc=color, ec="black", alpha=0.25, lw=1.2))
    ax.text(x, y, text, ha="center", va="center", fontsize=fs, wrap=True)


def arrow(ax, x0, y0, x1, y1, label=None):
    ax.annotate("", xy=(x1, y1), xytext=(x0, y0), arrowprops=dict(arrowstyle="->", lw=1.5))
    if label:
        ax.text((x0 + x1) / 2 + 0.02, (y0 + y1) / 2, label, fontsize=8, color=GREY, va="center")


# ---------------------------------------------------------------- fig1 chain
def fig1():
    fig, ax = plt.subplots(figsize=(9, 9.5))
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.axis("off")
    rows = [
        (0.93, "Velocity-line rule (ours, new input)\nfree-fall velocities, measured from the start with the dial addition,\nrun along one straight line", RED),
        (0.79, "Pull = K / r²   (Hamilton's converse; velocity-lines V7)\nK = G·M  (G is the one constant put in)", ORANGE),
        (0.66, "Energy E = v²/2 − K/r is conserved (R1)\nand sets the dial: velocity-space curvature = −2E (R2b)", ORANGE),
        (0.52, "The river = everything falling from rest far away: E = 0 (R2)\nflat middle of the dial, speed v² = 2K/r = escape speed", ORANGE),
        (0.37, "+ the theory's local relativity (dial at −1/c²)\nmetric  −dt² + (dr + v dt)² + r² dΩ²", GREEN),
        (0.23, "Vacuum (Einstein tensor = 0, R3) and exactly Schwarzschild (R4)\nriver clocks are free-fall proper time; γ = β = 1", ORANGE),
        (0.08, "Mercury 42.98\"/century · light bending 1.751\" · frame dragging (first write-up)", BLUE),
    ]
    for y, txt, col in rows:
        box(ax, 0.5, y, 0.86, 0.082, txt, col)
    for (y0, *_), (y1, *_) in zip(rows[:-1], rows[1:]):
        arrow(ax, 0.5, y0 - 0.045, 0.5, y1 + 0.045)
    ax.set_title("Fig 1. Gravity rebuilt from the velocity-line rule\n"
                 "red = new assumption · green = theory already had · orange = derived · blue = data", fontsize=11)
    fig.savefig(os.path.join(FIG, "fig1_chain.png"), dpi=130, bbox_inches="tight")
    plt.close(fig)


# ---------------------------------------------------------------- fig2 river = E=0
def fig2():
    fig, (a1, a2) = plt.subplots(1, 2, figsize=(12, 5.4))
    K = 1.0
    th = np.linspace(0, 2 * np.pi, 400)
    # same angular momentum L, three energies: velocity circles centre (0, e K/L) radius K/L in perifocal frame
    L = 1.0
    for e, col, lab in [(0.5, BLUE, "bound, E < 0"), (1.0, GREEN, "river family, E = 0"), (1.6, RED, "escaping, E > 0")]:
        R = K / L
        cx, cy = 0.0, e * K / L
        if e <= 1:
            a1.plot(cx + R * np.cos(th), cy + R * np.sin(th), color=col, lw=2, label=lab)
        else:
            phimax = math.acos(-1 / e)
            ph = np.linspace(-phimax + 1e-3, phimax - 1e-3, 300)
            vx, vy = -R * np.sin(ph), R * (e + np.cos(ph))
            a1.plot(vx, vy, color=col, lw=2, label=lab)
    a1.plot(0, 0, "ko")
    a1.annotate("v = 0\n(at rest far away)", (0, 0), (0.6, -0.55), arrowprops=dict(arrowstyle="->"), fontsize=9)
    a1.set_aspect("equal")
    a1.axhline(0, color=GREY, lw=0.5)
    a1.axvline(0, color=GREY, lw=0.5)
    a1.set_xlabel("velocity, x part")
    a1.set_ylabel("velocity, y part")
    a1.legend(loc="upper right", fontsize=9)
    a1.set_title("Velocity circles (same L): only the E = 0 circle\npasses through v = 0, the dial's flat middle")
    # river speed
    r = np.linspace(2.0, 30, 300)
    a2.plot(r, np.sqrt(2 * K / r), color=GREEN, lw=2.5, label="river speed √(2K/r) (from the rule)")
    a2.plot(r, np.sqrt(2 * K / r**1.5), color=GREY, lw=1, ls="--", label="other laws v² ∝ 1/r^p (p = 1.5, 0.5)")
    a2.plot(r, np.sqrt(2 * K / r**0.5) * 0.5, color=GREY, lw=1, ls="--")
    for rr in [4, 8, 14, 22]:
        v = math.sqrt(2 * K / rr)
        a2.annotate("", xy=(rr - 2.2 * v, 0.06), xytext=(rr, 0.06), arrowprops=dict(arrowstyle="->", color=BLUE))
    a2.text(12, 0.11, "space flows inward at the escape speed", color=BLUE, fontsize=9)
    a2.set_xlabel("areal radius r (units K/c²)")
    a2.set_ylabel("inflow speed (units c)")
    a2.set_ylim(0, 1.05)
    a2.legend(fontsize=9)
    a2.set_title("The river is the E = 0 family\n(its speed was Newton's input before; now an output)")
    fig.suptitle("Fig 2. Where the river comes from", fontsize=12, y=1.04)
    fig.savefig(os.path.join(FIG, "fig2_river_is_E0.png"), dpi=130, bbox_inches="tight")
    plt.close(fig)


# ---------------------------------------------------------------- fig3 Einstein for v^2 = 2M/r^p
def fig3():
    fig, ax = plt.subplots(figsize=(7, 4.6))
    r = np.linspace(3, 30, 300)
    M = 1.0
    for p_, col in [(0.5, ORANGE), (1.0, GREEN), (1.5, RED), (2.0, BLUE)]:
        Gtt = 2 * M * (p_ - 1) * r**(-p_ - 2)
        ax.plot(r, Gtt, color=col, lw=2.5 if p_ == 1 else 1.5, label=f"p = {p_:g}" + ("  (the rule's river)" if p_ == 1 else ""))
    ax.axhline(0, color="k", lw=0.6)
    ax.set_xlabel("r")
    ax.set_ylabel("Einstein Gᵗₜ  (= −8π × energy density needed)")
    ax.legend()
    ax.set_title("Fig 3. A river v² = 2M/rᵖ is empty space only for p = 1\n"
                 "Gᵗₜ = −(r v²)′ / r²  (R5; this is the shell theorem, known)")
    fig.savefig(os.path.join(FIG, "fig3_einstein_p.png"), dpi=130, bbox_inches="tight")
    plt.close(fig)


# ---------------------------------------------------------------- fig4 dimensions
def orbit_hodo(power, L, r0, tmax):
    def f(t, s):
        x, y, vx, vy = s
        r = math.hypot(x, y)
        a = -1 / r**power
        return [vx, vy, a * x / r, a * y / r]
    v_c = math.sqrt(1 / r0**(power - 1))
    s0 = [r0, 0, 0, (0.8 if power == 2 else 0.97) * v_c]
    stop = lambda t, s: math.hypot(s[0], s[1]) - (0.5 if power == 3 else 0.2)  # noqa: E731
    stop.terminal = True
    sol = solve_ivp(f, [0, tmax], s0, rtol=1e-10, atol=1e-12, max_step=tmax / 2000, events=stop)
    return sol.y


def fig4():
    fig, axs = plt.subplots(1, 2, figsize=(11, 5))
    for ax, power, title, col in [(axs[0], 2, "3+1: Einstein vacuum gives v² ∝ 1/r\npull 1/r² → velocity path is a circle (straight on the dial)", GREEN),
                                  (axs[1], 3, "4+1: Einstein vacuum gives v² ∝ 1/r²\npull 1/r³ → velocity path bends off any circle\n(shown until the orbit falls to half its start radius)", RED)]:
        y = orbit_hodo(power, 1, 1.0, 6.0 if power == 2 else 20.0)
        vx, vy = y[2], y[3]
        ax.plot(vx, vy, color=col, lw=2, label="velocity path")
        # best fit circle
        A = np.c_[2 * vx, 2 * vy, np.ones_like(vx)]
        sol, *_ = np.linalg.lstsq(A, vx**2 + vy**2, rcond=None)
        cx, cy = sol[0], sol[1]
        R = math.sqrt(sol[2] + cx**2 + cy**2)
        th = np.linspace(0, 2 * np.pi, 300)
        ax.plot(cx + R * np.cos(th), cy + R * np.sin(th), color=GREY, ls="--", lw=1, label="best-fit circle")
        resid = np.max(np.abs(np.hypot(vx - cx, vy - cy) - R)) / R
        ax.text(0.03, 0.03, f"misfit {resid:.1e}", transform=ax.transAxes, fontsize=9)
        ax.plot(0, 0, "k+")
        ax.set_aspect("equal")
        ax.set_title(title, fontsize=10)
        ax.set_xlabel("velocity x")
        ax.set_ylabel("velocity y")
        ax.legend(fontsize=8, loc="upper right")
    fig.suptitle("Fig 4. The rule and Einstein's equation agree only in 3 space dimensions (R6)", fontsize=12, y=1.06)
    fig.savefig(os.path.join(FIG, "fig4_dimensions.png"), dpi=130, bbox_inches="tight")
    plt.close(fig)


# ---------------------------------------------------------------- fig5 tests
def fig5():
    fig, (a1, a2) = plt.subplots(1, 2, figsize=(12, 5.2))
    # exaggerated rosette from Schwarzschild orbit equation u'' + u = M/L^2 + 3 M u^2 (M = 1)
    M, p, e = 1.0, 20.0, 0.5
    L2 = M * p / (1 - 3 * M / p * 0 - 0)  # approximate start; trace one equation
    u0 = (1 + e) / p

    def f(phi, s):
        u, du = s
        return [du, M / L2 - u + 3 * M * u**2]
    sol = solve_ivp(f, [0, 6 * 2 * np.pi], [u0, 0], rtol=1e-10, atol=1e-12, max_step=0.01)
    phi, u = sol.t, sol.y[0]
    a1.plot(np.cos(phi) / u, np.sin(phi) / u, color=BLUE, lw=1)
    a1.plot(0, 0, "o", color=ORANGE, ms=10)
    a1.set_aspect("equal")
    a1.set_title("Orbit in the rebuilt river (exaggerated, p = 20 GM/c²)\nthe ellipse turns: perihelion advance", fontsize=10)
    a1.axis("off")
    # data ratios
    names = ["Mercury perihelion\n42.98\"/cy", "Sun light bending\n1.751\" (via γ)", "Frame dragging GP-B\n39.1 mas/yr"]
    obs = [1.0, 1.0, 37.2 / 39.1]
    err = [0.04 / 42.98, 1.2e-4, 7.2 / 39.1]
    a2.errorbar(range(3), obs, yerr=err, fmt="o", color=BLUE, capsize=6, ms=8, label="observed / predicted")
    a2.axhline(1, color=GREEN, lw=2, label="prediction of the rebuilt river")
    a2.set_xticks(range(3))
    a2.set_xticklabels(names, fontsize=9)
    a2.set_ylim(0.6, 1.3)
    a2.set_ylabel("observed ÷ predicted")
    a2.legend(fontsize=9)
    a2.set_title("Data (Mercury: planetary ephemerides; light: VLBI γ − 1 = (−0.8 ± 1.2)×10⁻⁴;\nGP-B 37.2 ± 7.2 mas/yr)", fontsize=9)
    fig.suptitle("Fig 5. The rebuilt river passes the classic tests (R7)", fontsize=12, y=1.05)
    fig.savefig(os.path.join(FIG, "fig5_tests.png"), dpi=130, bbox_inches="tight")
    plt.close(fig)


# ---------------------------------------------------------------- fig6 inputs before / after
def fig6():
    fig, ax = plt.subplots(figsize=(11, 6.2))
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.axis("off")
    ax.text(0.25, 0.95, "Before (general-relativity page)", ha="center", fontsize=12, weight="bold")
    ax.text(0.75, 0.95, "After (this page)", ha="center", fontsize=12, weight="bold")
    box(ax, 0.25, 0.78, 0.42, 0.11, "Local relativity in free fall\n(the theory's dial at −1/c²)", GREEN)
    box(ax, 0.25, 0.58, 0.42, 0.13, "River speed v² = 2Gm/r\n(Newton's law, put in by hand)", RED)
    box(ax, 0.25, 0.30, 0.42, 0.11, "Schwarzschild, tests pass", BLUE)
    arrow(ax, 0.25, 0.51, 0.25, 0.36)
    box(ax, 0.75, 0.82, 0.42, 0.06, "Local relativity (theory's dial at −1/c²)", GREEN)
    box(ax, 0.75, 0.69, 0.42, 0.08, "Velocity-line rule (new, stated in the theory)\n+ G (sets K = G·M)", RED)
    box(ax, 0.75, 0.54, 0.42, 0.08, "Bridge to exact: rule's river holds in proper time\nand areal radius, OR vacuum Einstein", RED, fs=9)
    box(ax, 0.75, 0.30, 0.42, 0.11, "1/r², river speed, Schwarzschild, tests:\nall derived", ORANGE)
    arrow(ax, 0.75, 0.48, 0.75, 0.36)
    ax.text(0.5, 0.12, "Newton's law is no longer an input. The rule replaces it, picks 3 space dimensions with Einstein,\n"
            "and fixes the river's speed. The bridge from the weak-field rule to the exact river is still a choice:\n"
            "the rule alone is Newtonian (bound orbits in GR precess, so their velocity lines slowly turn).",
            ha="center", fontsize=10)
    ax.set_title("Fig 6. What gravity puts in, before and after", fontsize=12)
    fig.savefig(os.path.join(FIG, "fig6_inputs.png"), dpi=130, bbox_inches="tight")
    plt.close(fig)


if __name__ == "__main__":
    for fn in (fig1, fig2, fig3, fig4, fig5, fig6):
        fn()
        print("wrote", fn.__name__)
