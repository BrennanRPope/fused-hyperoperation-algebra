"""DRAFT figures for spinning-river/README.md. Run: python3 figures.py  (writes figures/*.png)"""
import json
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
from matplotlib.patches import Circle, FancyBboxPatch, Polygon  # noqa: E402
from scipy.integrate import solve_ivp  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
FIG = os.path.join(HERE, "figures")
os.makedirs(FIG, exist_ok=True)
RED, GREEN, ORANGE, BLUE, PURPLE, GREY = "#d62728", "#2ca02c", "#ff7f0e", "#1f77b4", "#9467bd", "#7f7f7f"


def save(fig, name):
    fig.savefig(os.path.join(FIG, name), dpi=130, bbox_inches="tight")
    plt.close(fig)


# ---------------------------------------------------------------- fig1 the cell's legs: one shear vs two
def sheared_square(ax, x0, y0, sx, sy, color, label=None):
    pts = np.array([[0, 0], [1, 0], [1, 1], [0, 1]], float)
    new = pts.copy()
    new[:, 1] = pts[:, 1] + sx * pts[:, 0]          # vertical (radial) leg picks up sx * horizontal leg
    ax.add_patch(Polygon(pts + [x0, y0], closed=True, fc="none", ec=GREY, ls=":", lw=1.4))
    ax.add_patch(Polygon(new + [x0, y0], closed=True, fc=color, alpha=0.35, ec=color, lw=2))
    if label:
        ax.text(x0 + 0.5, y0 - 0.32, label, ha="center", fontsize=9.5)


def fig1():
    fig, axs = plt.subplots(1, 2, figsize=(12, 5))
    ax = axs[0]
    sheared_square(ax, 0, 0, 0.55, None, BLUE, "time leg  →  radial leg\npushed by v (a Galilean boost)")
    ax.text(1.6, 0.6, "that is all a\nnon-spinning cell does", fontsize=10)
    ax.set_title("non-spinning hole: one shear per cell (Schwarzschild)")
    ax.set_xlim(-0.4, 2.8); ax.set_ylim(-0.8, 1.9); ax.set_aspect("equal"); ax.axis("off")
    ax = axs[1]
    sheared_square(ax, 0, 0, 0.55, None, BLUE, "time leg → radial leg\nby v (boost)")
    sheared_square(ax, 1.8, 0.35, -0.35, None, PURPLE, None)
    ax.text(2.3, -0.45, "azimuth leg → radial leg\nby −v·α (the turn)", ha="center", fontsize=9.5)
    ax.text(1.4, 1.75, "α = a sin θ / √(r² + a²)", ha="center", fontsize=10, color=PURPLE)
    ax.set_title("spinning hole: two shears per cell (Kerr), and they commute")
    ax.set_xlim(-0.4, 3.2); ax.set_ylim(-0.8, 1.9); ax.set_aspect("equal"); ax.axis("off")
    fig.text(0.5, 0.08, "dotted: the flat lattice's legs; colored: the cell's legs. Both shears live at the dial's flat "
             "middle (q = 0), where boosts and turns both become shears: (1 + εv)(1 + εw) = 1 + ε(v + w).",
             ha="center", fontsize=9.5)
    save(fig, "fig1_cell_legs.png")


# ---------------------------------------------------------------- fig2 the spheroidal lattice and the turn
def fig2():
    a, M = 0.7, 1.0
    fig, ax = plt.subplots(figsize=(8.5, 7))
    X, Z = np.meshgrid(np.linspace(-4.6, 4.6, 560), np.linspace(-3.8, 3.8, 460))
    # oblate spheroidal r from (rho_cyl, z): r^4 - (s^2 - a^2) r^2 - a^2 z^2 = 0, s^2 = x^2 + z^2
    s2 = X**2 + Z**2
    r2 = 0.5 * ((s2 - a**2) + np.sqrt((s2 - a**2)**2 + 4 * a**2 * Z**2))
    r = np.sqrt(r2)
    sinth2 = np.clip(X**2 / (r2 + a**2), 0, 1)
    alpha = a * np.sqrt(sinth2) / np.sqrt(r2 + a**2)
    v = np.sqrt(2 * M * r / (r2 + a**2 * (1 - sinth2)))
    turn = v * alpha
    im = ax.pcolormesh(X, Z, turn, cmap="Purples", shading="auto", vmin=0, vmax=0.8)
    plt.colorbar(im, ax=ax, label="size of the cell's turn  v·α")
    th = np.linspace(0, 2 * np.pi, 400)
    for rr in [0.35, 0.8, 1.4, 2.0, 2.7, 3.4]:
        ax.plot(np.sqrt(rr**2 + a**2) * np.sin(th), rr * np.cos(th), color="k", lw=0.6)
    for t0 in np.linspace(0.15, np.pi - 0.15, 9):
        rs = np.linspace(0, 4, 100)
        ax.plot(np.sqrt(rs**2 + a**2) * np.sin(t0), rs * np.cos(t0), color="k", lw=0.4, ls=":")
        ax.plot(-np.sqrt(rs**2 + a**2) * np.sin(t0), rs * np.cos(t0), color="k", lw=0.4, ls=":")
    rp = M + np.sqrt(M**2 - a**2)
    ax.plot(np.sqrt(rp**2 + a**2) * np.sin(th), rp * np.cos(th), color=RED, lw=2.5, label=f"horizon r₊ = {rp:.2f}")
    ax.plot([-a, a], [0, 0], "o", color=ORANGE, ms=8, label="the ring r = 0 (radius a)")
    ax.annotate("", (0, 3.3), (0, -3.3), arrowprops=dict(arrowstyle="->", color=GREY))
    ax.text(0.08, 3.15, "spin axis", color=GREY, fontsize=9)
    ax.set_xlim(-4.4, 4.4); ax.set_ylim(-3.6, 3.6)
    ax.set_aspect("equal"); ax.set_xlabel("distance from the axis"); ax.set_ylabel("height along the axis")
    ax.legend(loc="lower right", fontsize=9)
    ax.set_title("the flat lattice in the spinning body's shape (squashed spheres)\nshading: how much each cell turns "
                 f"(M = 1, a = {a})")
    save(fig, "fig2_spheroid_lattice.png")


# ---------------------------------------------------------------- fig3 the turn is forced
def fig3():
    res = json.load(open(os.path.join(HERE, "results.json")))

    def num(key, sub=None):
        d = res[key]["detail"]
        if sub is None:
            return float(d.split("=")[1].split()[0])
        return float(d.split(f"'{sub}': '")[1].split("'")[0])
    labels = ["velocity + turn\n(exact Kerr)", "near-max spin\na = 0.99", "no turn\n(velocity only)", "½ the turn",
              "1.5× the turn", "turn ∝ sin θ\n(not sin²θ)", "speed √(2Mr)", "speed √(2Mr)/ρ²"]
    k2a = [k for k in res if k.startswith("K2a")][0]; k2b = [k for k in res if k.startswith("K2b")][0]
    k3a = [k for k in res if k.startswith("K3a")][0]; k3b = [k for k in res if k.startswith("K3b")][0]
    k3c = [k for k in res if k.startswith("K3c")][0]
    vals = [num(k2a), num(k2b), num(k3a), num(k3b, "half the turn"), num(k3b, "1.5 x the turn"),
            num(k3b, "turn ~ sin(theta) not sin^2"), num(k3c, "speed sqrt(2Mr) (rho^0)"), num(k3c, "speed sqrt(2Mr)/rho^2")]
    fig, ax = plt.subplots(figsize=(11, 5))
    cols = [GREEN, GREEN] + [RED] * 6
    ax.bar(range(len(vals)), vals, color=cols)
    ax.set_yscale("log"); ax.set_ylim(1e-55, 1e3)
    ax.axhline(1e-30, color="k", ls="--", lw=0.8); ax.text(-0.4, 2e-29, "zero line (1e-30)", ha="left", fontsize=9)
    ax.set_xticks(range(len(vals))); ax.set_xticklabels(labels, fontsize=9)
    ax.set_ylabel("largest field-equation error  max |Ricci|  (log scale)")
    ax.set_title("only one turn and one speed law work: everything else breaks the field equations")
    save(fig, "fig3_turn_is_forced.png")


# ---------------------------------------------------------------- fig4 the clock gap
def fig4():
    fig, axs = plt.subplots(1, 2, figsize=(12, 5.2))
    ax = axs[0]
    a = 0.7
    ph = np.linspace(0, 2 * np.pi, 13)
    ax.add_patch(Circle((0, 0), 1, fc="none", ec="k"))
    for k, p in enumerate(ph[:-1]):
        x, y = np.cos(p), np.sin(p)
        ax.add_patch(Circle((1.22 * x, 1.22 * y), 0.13, fc="white", ec=PURPLE))
        ang = np.pi / 2 - k / 12 * 2 * np.pi * 0.85
        ax.plot([1.22 * x, 1.22 * x + 0.1 * np.cos(ang)], [1.22 * y, 1.22 * y + 0.1 * np.sin(ang)], color=PURPLE, lw=1.6)
    ax.annotate("", (1.05, -0.35), (1.05, 0.35), arrowprops=dict(arrowstyle="<->", color=RED))
    ax.text(1.5, 0.0, "after one lap the\nclocks disagree by\n2π a sin²θ", color=RED, fontsize=9.5, va="center")
    ax.text(0, 0, "spin\naxis\n(view from above)", ha="center", va="center", fontsize=9, color=GREY)
    ax.set_xlim(-1.6, 2.6); ax.set_ylim(-1.6, 1.6); ax.set_aspect("equal"); ax.axis("off")
    ax.set_title("the cell's time leg dt − a sin²θ dφ does not close")
    ax = axs[1]
    r = np.linspace(4, 30, 200)
    Tp, Tm = 2 * np.pi * (r**1.5 + a), 2 * np.pi * (r**1.5 - a)
    ax.plot(r, Tp - 2 * np.pi * r**1.5, color=BLUE, lw=2.5, label="prograde period − Kepler's")
    ax.plot(r, Tm - 2 * np.pi * r**1.5, color=ORANGE, lw=2.5, label="retrograde period − Kepler's")
    ax.axhline(0, color="k", lw=0.5)
    ax.set_ylim(-6, 6)
    ax.text(17, 3.0, "+2πa at every radius", color=BLUE); ax.text(17, -3.6, "−2πa at every radius", color=ORANGE)
    ax.set_xlabel("orbit radius r (units of M)"); ax.set_ylabel("period shift (units of GM/c³)")
    ax.legend(loc="center right", fontsize=9)
    ax.set_title("Kerr's clock effect: each shift is the clock gap\nof one lap at the equator (T₊ − T₋ = 4πa)")
    save(fig, "fig4_clock_gap.png")


# ---------------------------------------------------------------- fig5 the family of frames
def fig5():
    fig, ax = plt.subplots(figsize=(13, 7.5)); ax.axis("off")
    cols = {+1: "#e6f2ff", 0: "#f2f2f2", -1: "#ffe9e6"}
    heads = {+1: "d = +1: turning (circle)", 0: "d = 0: the flat middle (shear)", -1: "d = −1: boosting (relativity)"}
    xs = {+1: 0.02, 0: 0.355, -1: 0.69}
    for d in (+1, 0, -1):
        ax.add_patch(FancyBboxPatch((xs[d], 0.08), 0.3, 0.8, boxstyle="round,pad=0.01", fc=cols[d], ec="k"))
        ax.text(xs[d] + 0.15, 0.91, heads[d], ha="center", weight="bold")
    items = {
        +1: [("Foucault's pendulum", "on our pages (C5–C7)"), ("Coriolis = the swap rule", "on our pages"),
             ("tumbling body's phase", "on our pages (D1–D3)"), ("Larmor precession, Berry phase", "on our pages (F1–F4)"),
             ("Larmor's theorem: B = turning frame", "new here (F1)"), ("coil's twist (Fuller)", "on our pages")],
        0: [("Galilean boosts: nothing precesses", "on our pages"), ("river cell: boost as a shear", "jet-river"),
            ("Kerr cell: the turn as a shear", "new here (K1–K5)"), ("acoustic metric: sound in a flow", "new here (F2)"),
            ("Bargmann phase m·area/ħ", "prediction on our pages")],
        -1: [("Thomas–Wigner rotation", "on our pages"), ("geodetic (de Sitter)", "GP-B, on our pages"),
             ("Lense–Thirring frame dragging", "Part F, on our pages"), ("rotating BTZ drag", "precession §8"),
             ("Kerr clock effect = clock gap", "new here (K7)")],
    }
    for d, lst in items.items():
        for k, (name, where) in enumerate(lst):
            y = 0.8 - 0.115 * k
            new = "new here" in where
            ax.text(xs[d] + 0.015, y, "• " + name, fontsize=9.6, weight="bold" if new else "normal",
                    color=PURPLE if new else "k")
            ax.text(xs[d] + 0.035, y - 0.04, where, fontsize=8.3, color=GREY)
    ax.text(0.5, 0.02, "one law runs through all three columns: turn per loop = d × area (the swap rule over a loop).  "
            "Purple: added on this page.", ha="center", fontsize=10)
    ax.set_xlim(0, 1); ax.set_ylim(0, 1)
    ax.set_title("every idea of a frame in our work, sorted by the dial", fontsize=12)
    save(fig, "fig5_frame_family.png")


# ---------------------------------------------------------------- fig6 Larmor's theorem and Foucault
def fig6():
    fig, axs = plt.subplots(1, 2, figsize=(12, 5.4))
    ax = axs[0]
    B, tmax = 0.06, 40.0

    def run(Bv, s0):
        def rhs(t, s):
            x, y, vx, vy = s
            r3 = (x * x + y * y)**1.5
            return [vx, vy, -x / r3 + Bv * vy, -y / r3 - Bv * vx]
        ts = np.linspace(0, tmax, 3000)
        sol = solve_ivp(rhs, (0, tmax), s0, t_eval=ts, rtol=1e-10, atol=1e-12)
        return ts, sol.y
    ts, yB = run(B, [1.0, 0, 0, 1.1])
    w = -B / 2
    ts, y0 = run(0.0, [1.0, 0, 0, 1.1 - w])
    c, s = np.cos(w * ts), np.sin(w * ts)
    xl, yl = c * y0[0] - s * y0[1], s * y0[0] + c * y0[1]
    ax.plot(y0[0], y0[1], color=GREY, lw=0.8, label="no field (orbit stays put)")
    ax.plot(yB[0], yB[1], color=BLUE, lw=2.4, label="charge in a magnetic field")
    ax.plot(xl, yl, color=ORANGE, lw=1.0, ls="--", label="no field, watched from a frame\nturning at Larmor's rate −qB/2m")
    ax.plot(0, 0, "k*", ms=10)
    ax.set_aspect("equal"); ax.legend(fontsize=8.5, loc="upper left", bbox_to_anchor=(1.0, 1.0))
    ax.set_title("Larmor's theorem: a magnetic field is a turning frame\n(blue and dashed orange lie on top of each other)")
    ax = axs[1]
    lat = np.linspace(0, 90, 200)
    ax.plot(lat, 360 * np.sin(np.radians(lat)) / 0.99727, color=GREEN, lw=2.5, label="Foucault turn (deg per day)")
    ax.set_xlabel("latitude (degrees)"); ax.set_ylabel("degrees per day", color=GREEN)
    ax2 = ax.twinx()
    drag = 2 * 0.3307 * 6.674e-11 * 5.972e24 / (6.371e6 * 2.998e8**2) * 7.2921e-5 * np.sin(np.radians(lat))
    ax2.plot(lat, drag * 3.156e7 * 206264806, color=PURPLE, lw=2, ls="--", label="Earth's own frame dragging (mas/yr)")
    ax2.set_ylabel("milliarcseconds per year", color=PURPLE)
    ax.set_title("Foucault's pendulum, plus Earth's frame dragging\n(4.6 × 10⁻¹⁰ of Foucault's turn, at every latitude)")
    ax.legend(loc="upper left", fontsize=9); ax2.legend(loc="lower right", fontsize=9)
    save(fig, "fig6_larmor_foucault.png")


# ---------------------------------------------------------------- fig7 the bathtub vortex
def fig7():
    fig, axs = plt.subplots(1, 2, figsize=(12.5, 5.6))
    ax = axs[0]
    A, B = 1.0, 1.6
    Y, X = np.mgrid[-4:4:200j, -4:4:200j]
    r2 = X**2 + Y**2 + 1e-9
    U, V = (-A * X - B * Y) / r2, (-A * Y + B * X) / r2
    ax.streamplot(X, Y, U, V, color=BLUE, density=1.2, linewidth=0.8)
    th = np.linspace(0, 2 * np.pi, 300)
    ax.fill(A * np.cos(th), A * np.sin(th), color="k")
    ax.plot(np.hypot(A, B) * np.cos(th), np.hypot(A, B) * np.sin(th), color=RED, lw=2, ls="--")
    ax.text(0, -0.15, "sound\nhorizon", color="w", ha="center", fontsize=8)
    ax.text(-3.9, 3.35, "red dashed: ergoregion edge (waves must co-rotate)", color=RED, fontsize=9,
            bbox=dict(fc="white", ec="none", alpha=0.85))
    ax.set_aspect("equal"); ax.set_xlim(-4, 4); ax.set_ylim(-4, 4)
    ax.set_title("a draining, spinning bathtub drags sound\n(the acoustic metric is a river, Unruh 1981)")
    ax = axs[1]
    r = np.linspace(1.0, 6, 200)
    ax.plot(r, B / r**2, color=BLUE, lw=2.4, label="vortex: swirl rate B/r²")
    ax.plot(r, 2 * 1.6 / r**3, color=PURPLE, lw=2.4, label="spinning mass: swirl rate 2J/r³")
    ax.plot(r, -(A**2 + B**2) / r**4, color=RED, lw=2, ls="--", label="vortex's Q layer = −(A²+B²)/r⁴ < 0")
    ax.axhline(0, color="k", lw=0.6)
    ax.text(3.6, -0.25, "gravity's empty space: Q = 0 (black line)", fontsize=9)
    ax.set_xlabel("radius r"); ax.set_ylabel("rate, or Q (arbitrary units)"); ax.legend(fontsize=9)
    ax.set_title("same dragging kinematics, different field law:\nthe vortex would need negative energy")
    save(fig, "fig7_bathtub.png")


# ---------------------------------------------------------------- fig8 the derivation chain (round 7)
def fig8():
    fig, ax = plt.subplots(figsize=(13, 5.8))
    W, H = 0.21, 0.3
    boxes = {
        "s1": (0.0, 0.62, "1. Shells turn rigidly\n(assumption)", "each lattice spheroid r turns\nas a whole at its own rate Ω(r)", BLUE),
        "s2": (0.26, 0.62, "2. Area law (known: Sagnac)", "clock gap around a latitude circle\n= 2 Ω × enclosed area\n→ turn f = Ω (r²+b²) sin²θ", BLUE),
        "s3": (0.52, 0.62, "3. Circulation law\n(assumption)", "a falling ring of cells keeps its gap\n(cells fall on cones, D2)\n→ Ω(r) = κ / (r²+b²)", PURPLE),
        "s4": (0.78, 0.62, "4. Frame dragging far away", "match Part F's −2J sin²θ / r\n→ κ = J / M", GREEN),
        "s5": (0.26, 0.12, "5. Rim at light speed\n(a reading of the dial's edge)", "innermost shell = disk of radius b,\nrim speed κ / b = 1  →  b = κ", ORANGE),
        "out": (0.52, 0.12, "Output: Kerr, a = J/M", "every empty-space field\nequation holds; clock gap 2πJ/M;\nhorizon turns at Ω_H = a/(r₊²+a²),\nthe shell rate", RED),
    }
    for x, y, head, body, col in boxes.values():
        ax.add_patch(FancyBboxPatch((x, y), W, H, boxstyle="round,pad=0.01", fc=col, alpha=0.13, ec=col, lw=2))
        ax.text(x + W / 2, y + 0.26, head, ha="center", va="top", fontsize=10.5, weight="bold", color=col)
        ax.text(x + W / 2, y + 0.1, body, ha="center", va="center", fontsize=9.2)
    arr = dict(arrowstyle="->", lw=1.8, color="k")
    for x0 in (0.0, 0.26, 0.52):
        ax.annotate("", xy=(x0 + 0.26 - 0.015, 0.77), xytext=(x0 + W + 0.015, 0.77), arrowprops=arr)
    ax.annotate("", xy=(0.52 + W - 0.02, 0.445), xytext=(0.80, 0.595), arrowprops=arr)      # step 4 -> output
    ax.annotate("", xy=(0.505, 0.27), xytext=(0.485, 0.27), arrowprops=arr)               # step 5 -> output
    ax.annotate("", xy=(0.56, 0.595), xytext=(0.44, 0.445),
                arrowprops=dict(arrowstyle="->", lw=1.3, color=ORANGE, ls="--"))          # b = kappa fills in b
    ax.text(0.505, 0.52, "b = κ", fontsize=9, color=ORANGE, ha="center")
    ax.text(0.0, 0.42, "still put in by hand:\n• flat spheroidal lattice\n• speed v² = 2Mr/ρ²\n  (from the last page,\n   matched, not derived)",
            fontsize=9.5, va="top", color=GREY)
    ax.set_xlim(-0.02, 1.02); ax.set_ylim(0.05, 1.0); ax.axis("off")
    ax.set_title("Deriving the turn: two assumptions, one known law, one far-field match and one reading;\n"
                 "the field equations then only confirm", fontsize=12)
    save(fig, "fig8_turn_derivation.png")


# ---------------------------------------------------------------- fig9 the turning shells
def fig9():
    M, a = 1.0, 0.7
    rp = M + np.sqrt(M**2 - a**2)
    fig, axs = plt.subplots(1, 2, figsize=(13, 5.4))
    ax = axs[0]
    th = np.linspace(0, 2 * np.pi, 300)
    for r in [0.0, 0.5, 1.0, 1.6, 2.4]:
        A_, C_ = np.sqrt(r**2 + a**2), r
        col = RED if r == 0 else BLUE
        ax.plot(A_ * np.cos(th), C_ * np.sin(th), color=col, lw=2.6 if r == 0 else 1.4)
        speed = a / A_
        ax.annotate("", xy=(A_, 0.18 + 1.2 * speed), xytext=(A_, 0.18),
                    arrowprops=dict(arrowstyle="->", color=col, lw=2))
        ax.text(A_ + 0.04, 0.18 + 1.2 * speed, f"{speed:.2f}c", fontsize=8.5, color=col, va="center")
    ax.annotate("r = 0 shell: a flat disk\nradius a, rim at c", xy=(-0.35, 0), xytext=(-2.7, -2.6), ha="left",
                va="bottom", fontsize=9, color=RED, arrowprops=dict(arrowstyle="-", color=RED, lw=0.8))
    ax.text(-2.7, 3.2, "each spheroid turns rigidly at Ω(r) = a/(r²+a²)\n"
            "arrow length = speed at the equator, a/√(r²+a²)\n(schematic: the real motion is into the page)",
            fontsize=9.5, va="top")
    ax.set_aspect("equal"); ax.set_xlim(-2.8, 3.2); ax.set_ylim(-2.7, 3.3)
    ax.set_xlabel("distance from axis"); ax.set_ylabel("height along the spin axis")
    ax.set_title("the flat lattice, shell by shell, turning")
    ax = axs[1]
    r = np.linspace(0, 6, 400)
    ax.plot(r, a / (r**2 + a**2), color=PURPLE, lw=2.6, label="shell rate Ω(r) = a/(r²+a²)  (derived)")
    ax.plot(r[r > 1.3], 2 * M * a / r[r > 1.3]**3, color=GREEN, lw=1.8, ls="--", label="far-field drag 2J/r³ (Part F)")
    ax.axvline(rp, color=GREY, ls=":", lw=1.2)
    OmH = a / (rp**2 + a**2)
    ax.plot([rp], [OmH], "o", color=RED, ms=8)
    ax.text(rp + 0.1, OmH + 0.08, f"horizon r₊ = {rp:.2f}:\nΩ_H = a/(r₊²+a²) = {OmH:.3f}\n(Kerr's known horizon spin)",
            fontsize=9, color=RED)
    ax.plot([0], [1 / a], "s", color=ORANGE, ms=8)
    ax.text(0.15, 1.45, "r = 0: rate 1/a, rim speed (1/a)·a = c", fontsize=9, color=ORANGE, va="center")
    ax.set_ylim(0, 1.6); ax.set_xlabel("spheroid label r (M = 1, a = 0.7)"); ax.set_ylabel("turning rate")
    ax.legend(fontsize=9, loc="center right"); ax.set_title("one curve: light-speed rim, horizon spin, far-field drag")
    save(fig, "fig9_turning_shells.png")


# ---------------------------------------------------------------- fig10 each assumption matters
def fig10():
    import re
    res = json.load(open(os.path.join(HERE, "results.json")))
    bad = res["D5c_bad"]
    good = float(re.search(r"= ([0-9.e+-]+)", res["D5b and it solves every empty-space field equation (M = 1, kappa = b = 0.7)"]["detail"]).group(1))
    names = {"ring b = 0.5 kappa (rim faster than light)": "ring b = 0.5κ\n(rim faster\nthan light)",
             "ring b = 1.5 kappa (rim slower than light)": "ring b = 1.5κ\n(rim slower\nthan light)",
             "gap not kept: Omega ~ 1/r^2": "gap not kept:\nΩ ~ 1/r²",
             "gap not kept: one rigid rate for all space": "gap not kept:\none rigid rate\nfor all space",
             "gap drifts with r: f = kappa sin^2 (1 + M/r)": "gap drifts with r:\nf = κ sin²θ (1 + M/r)"}
    labs = ["derived chain\n(Kerr)"] + [names.get(k, k) for k in bad]
    vals = [good] + list(bad.values())
    fig, ax = plt.subplots(figsize=(12, 4.8))
    ax.bar(range(len(vals)), vals, color=[GREEN] + [RED] * len(bad), alpha=0.8)
    for i, v in enumerate(vals[1:], 1):
        ax.text(i, v * 3, f"{v:.2g}", ha="center", fontsize=8.5)
    ax.text(0, good * 3, f"{good:.1e}", ha="center", fontsize=8.5)
    ax.set_yscale("log"); ax.set_ylim(1e-52, 1e4)
    ax.set_xticks(range(len(vals))); ax.set_xticklabels(labs, fontsize=8.5)
    ax.axhline(1e-30, color=GREY, ls=":")
    ax.text(-0.45, 3e-30, "zero line (1e-30)", fontsize=8.5, color=GREY, ha="left", va="bottom")
    ax.set_ylabel("largest empty-space field-equation error |Ricci|")
    ax.set_title("change any step of the derivation and empty space breaks (M = 1, κ = 0.7)")
    save(fig, "fig10_each_step_matters.png")


# ---------------------------------------------------------------- fig11 the complex shift (round 8)
def fig11():
    a, M = 0.7, 1.0
    fig, axs = plt.subplots(1, 2, figsize=(13, 5.6), gridspec_kw=dict(width_ratios=[1, 1.6]))
    ax = axs[0]
    ax.axhline(0, color="k", lw=0.8); ax.axvline(0, color="k", lw=0.8)
    ax.plot([0], [0], "o", color=BLUE, ms=11)
    ax.plot([0], [a], "o", color=RED, ms=11)
    ax.annotate("", xy=(0, a - 0.06), xytext=(0, 0.07), arrowprops=dict(arrowstyle="->", lw=2, color=RED))
    ax.plot([a], [0], "o", mfc="none", mec=GREY, ms=11)
    ax.annotate("", xy=(a - 0.06, 0), xytext=(0.07, 0), arrowprops=dict(arrowstyle="->", lw=1.4, color=GREY, ls="--"))
    ax.text(0.08, a + 0.08, "spinning: source moved\nby i·a (a quarter turn)", color=RED, fontsize=10)
    ax.text(0.08, -0.22, "not spinning:\nsource at 0", color=BLUE, fontsize=10)
    ax.text(a + 0.05, 0.07, "a real shift by a\n(fails, S4b)", color=GREY, fontsize=9)
    ax.set_xlim(-0.6, 1.5); ax.set_ylim(-0.5, 1.2); ax.set_aspect("equal")
    ax.set_xlabel("real part of the source's height on the axis"); ax.set_ylabel("imaginary part")
    ax.set_title("one move: turn the source's height\nby the dial's rotation number i")
    ax = axs[1]
    X, Z = np.meshgrid(np.linspace(-3.2, 3.2, 640), np.linspace(-2.4, 2.4, 480))
    s2 = X**2 + Z**2
    r = np.sqrt(((s2 - a**2) + np.sqrt((s2 - a**2)**2 + 4 * a**2 * Z**2)) / 2)
    cth = np.where(r > 1e-9, Z / np.maximum(r, 1e-9), 0.0)
    v2 = 2 * M * r / (r**2 + a**2 * cth**2)
    cf = ax.contourf(X, Z, np.minimum(v2, 4), levels=np.linspace(0, 4, 17), cmap="Blues", alpha=0.8)
    ax.contour(X, Z, r, levels=[0.4, 0.8, 1.2, 1.6, 2.0, 2.4], colors=PURPLE, linewidths=1.0)
    ax.plot([-a, a], [0, 0], color=GREY, lw=1.5, ls="--")
    ax.plot([-a, a], [0, 0], "o", color=RED, ms=7)
    ax.annotate("D = 0: the ring of radius a\n(the same ring as the turn's)", xy=(a, 0), xytext=(1.7, -1.75),
                ha="center", va="center", fontsize=9, color=RED, bbox=dict(fc="white", ec="none", alpha=0.85),
                arrowprops=dict(arrowstyle="->", color=RED, lw=1.2))
    ax.text(-3.1, 2.15, "purple: Re D = r, the lattice's squashed spheres\nblue: speed² = 2 Re(M/D) = 2Mr/ρ²\ngrey dashes: the flat disk r = 0 (the ring is its edge)",
            fontsize=9.5, va="top", bbox=dict(fc="white", ec="none", alpha=0.85))
    fig.colorbar(cf, ax=ax, label="cell speed² (capped at 4)")
    ax.set_aspect("equal"); ax.set_xlabel("distance from axis"); ax.set_ylabel("height along the spin axis")
    ax.set_title("what the shifted source does in real space:\nlattice, ring and speed law at once")
    save(fig, "fig11_complex_shift.png")


# ---------------------------------------------------------------- fig12 harmonic is not enough
def fig12():
    import re
    res = json.load(open(os.path.join(HERE, "results.json")))
    bad = res["S4b_bad"]
    good = float(re.search(r"= ([0-9.e+-]+)", res["D5b and it solves every empty-space field equation (M = 1, kappa = b = 0.7)"]["detail"]).group(1))
    names = {"point at the centre: M/|x|": "point at the\ncentre: M/|x|",
             "point moved by a real a up the axis": "point moved\nby a real a",
             "point moved by i a/2 (ring a/2, not the turn's ring)": "point moved by i·a/2\n(ring a/2, not the\nturn's ring)"}
    labs = ["moved by i·a\n(Kerr)"] + [names.get(k, k) for k in bad]
    vals = [good] + list(bad.values())
    fig, ax = plt.subplots(figsize=(8, 4.8))
    ax.bar(range(len(vals)), vals, color=[GREEN] + [RED] * len(bad), alpha=0.8)
    for i, v in enumerate(vals[1:], 1):
        ax.text(i, v * 3, f"{v:.2g}", ha="center", fontsize=9)
    ax.text(0, good * 3, f"{good:.1e}", ha="center", fontsize=9)
    ax.set_yscale("log"); ax.set_ylim(1e-52, 1e4)
    ax.set_yticks([10.0**k for k in range(-50, 5, 10)])
    ax.set_xticks(range(len(vals))); ax.set_xticklabels(labs, fontsize=9.5)
    ax.axhline(1e-30, color=GREY, ls=":")
    ax.text(-0.45, 3e-30, "below this = zero to rounding (1e-30)", fontsize=8.5, color=GREY, ha="left", va="bottom")
    ax.set_ylabel("largest empty-space field-equation error |Ricci|")
    ax.set_title("all four speed laws are flat-space potentials of mass M;\nonly the i·a shift gives empty space")
    save(fig, "fig12_only_the_i_shift.png")


# ---------------------------------------------------------------- fig13 the shifted source's rays (round 9)
def fig13():
    fig, axs = plt.subplots(1, 2, figsize=(12.5, 6.2))
    for ax, a, title in [(axs[0], 0.0, "no spin: source at 0\nrays run straight out from the centre"),
                         (axs[1], 0.7, "source moved by i·a: every ray leans sideways,\ntangent to the ring of radius a")]:
        L = 3.4
        for ang in np.linspace(0, 2 * np.pi, 24, endpoint=False):
            u = np.array([np.cos(ang), np.sin(ang)])
            perp = np.array([np.sin(ang), -np.cos(ang)])     # clockwise tangent: sideways part -a sin^2 dphi
            start = a * perp
            ax.plot([start[0], start[0] + L * u[0]], [start[1], start[1] + L * u[1]], color=BLUE, lw=1.0, alpha=0.8)
            ax.annotate("", xy=start + 2.2 * u, xytext=start + 1.9 * u, arrowprops=dict(arrowstyle="->", color=BLUE, lw=1.2))
        th = np.linspace(0, 2 * np.pi, 300)
        for rr in [1.0, 2.0]:
            A = np.sqrt(rr**2 + a**2)
            ax.plot(A * np.cos(th), A * np.sin(th), color=PURPLE, lw=1.0, ls=":")
            ax.text(0.35, A + 0.08, f"r = {rr:g}", ha="center", fontsize=8, color=PURPLE)
        if a > 0:
            ax.plot(a * np.cos(th), a * np.sin(th), color=RED, lw=2.4)
            ax.plot([0, 2.0], [0, -0.7], color=GREY, ls="--", lw=1)
            ax.annotate("lean = a/√(r²+a²)\n= round 7's α (shell speed)", xy=(2.0, -0.7),
                        xytext=(3.55, -3.55), ha="right", va="bottom", fontsize=8.5, color=PURPLE,
                        arrowprops=dict(arrowstyle="->", color=PURPLE, lw=1))
            ax.text(0, 0, "ring\na = 0.7", ha="center", va="center", fontsize=9, color=RED)
        else:
            ax.plot([0], [0], "o", color=BLUE, ms=8)
        ax.set_xlim(-3.6, 3.6); ax.set_ylim(-3.6, 3.6); ax.set_aspect("equal")
        ax.set_xlabel("x (equator's plane)"); ax.set_ylabel("y")
        ax.set_title(title)
    fig.text(0.5, -0.02, "The rays are the real light rays of the shifted source: unit k with k·E = 1, E the complex unit "
             "direction to the source.\nTheir sideways part, −a sin²θ dφ, is the cell's turn.", ha="center", fontsize=9.5)
    save(fig, "fig13_shifted_rays.png")


# ---------------------------------------------------------------- fig14 one move, everything
def fig14():
    fig, ax = plt.subplots(figsize=(13, 5.4))

    def box(x, y, w, h, head, body, col):
        ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.01", fc=matplotlib.colors.to_rgba(col, 0.13),
                                    ec=col, lw=2))
        ax.text(x + w / 2, y + h - 0.04, head, ha="center", va="top", fontsize=10.5, weight="bold", color=col)
        ax.text(x + w / 2, y + 0.38 * h, body, ha="center", va="center", fontsize=9.2)

    box(0.0, 0.36, 0.2, 0.3, "One move", "source's height on the axis\nturned by i: 0 → i·a\n(the dial's rotation number)", RED)
    box(0.28, 0.36, 0.18, 0.3, "Complex distance", "D = r − i a cos θ\nE = grad D, a complex\nunit direction", PURPLE)
    box(0.56, 0.7, 0.2, 0.27, "Lattice (Part 4)", "Re D = r: the squashed\nspheres; D = 0: the ring", BLUE)
    box(0.56, 0.37, 0.2, 0.27, "Speed (Part 4)", "v² = 2 Re(M/D) = 2Mr/ρ²\n(jet river: v²/2 = −Φ)", BLUE)
    box(0.56, 0.04, 0.2, 0.27, "Turn (Part 5)", "light rays k·E = 1 lean\nsideways: f = a sin²θ", BLUE)
    box(0.83, 0.36, 0.17, 0.3, "Kerr", "every empty-space\nfield equation holds;\nfar field: a = J/M", GREEN)
    arr = dict(arrowstyle="->", lw=1.8, color="k")
    ax.annotate("", xy=(0.262, 0.51), xytext=(0.218, 0.51), arrowprops=arr)
    for yy in (0.83, 0.5, 0.17):
        ax.annotate("", xy=(0.542, yy), xytext=(0.478, {0.83: 0.58, 0.5: 0.51, 0.17: 0.44}[yy]), arrowprops=arr)
        ax.annotate("", xy=(0.812, {0.83: 0.60, 0.5: 0.51, 0.17: 0.42}[yy]), xytext=(0.778, yy), arrowprops=arr)
    ax.text(0.0, 0.25, "Still put in by hand:\nthe mass M, the river's rules,\nand the reading that spin is the shift",
            fontsize=9.5, color=GREY, va="top")
    ax.set_xlim(-0.02, 1.02); ax.set_ylim(0, 1.0); ax.axis("off")
    ax.set_title("One move gives the whole spinning black hole", fontsize=12)
    save(fig, "fig14_one_move.png")


# ---------------------------------------------------------------- fig15 charge on the shifted source (round 10)
def fig15():
    import re
    a, Q = 0.7, 0.5
    res = json.load(open(os.path.join(HERE, "results.json")))
    fig, axs = plt.subplots(1, 2, figsize=(13.5, 5.6), gridspec_kw=dict(width_ratios=[1.25, 1]))
    ax = axs[0]
    X, Z = np.meshgrid(np.linspace(-3.0, 3.0, 600), np.linspace(-2.4, 2.4, 480))
    s2 = X**2 + Z**2
    r = np.sqrt(((s2 - a**2) + np.sqrt((s2 - a**2)**2 + 4 * a**2 * Z**2)) / 2)
    cth = np.where(r > 1e-9, Z / np.maximum(r, 1e-9), 0.0)
    rho2 = r**2 + a**2 * cth**2
    re_ = Q * r / rho2
    im_ = Q * a * cth / rho2
    c1 = ax.contour(X, Z, re_, levels=[0.12, 0.18, 0.25, 0.35, 0.5], colors=BLUE, linewidths=1.3)
    c2 = ax.contour(X, Z, im_, levels=[-0.2, -0.12, -0.06, -0.03, 0.03, 0.06, 0.12, 0.2], colors=ORANGE, linewidths=1.1,
                    linestyles="--")
    ax.plot([-a, a], [0, 0], "o", color=RED, ms=6, zorder=5)
    ax.text(-2.9, 3.3, "Q/D, one complex potential of the shifted charge:\nblue: real part = electric potential Qr/ρ²\n"
            "orange dashes: imaginary part = magnetic potential\n   Qa cos θ/ρ², a magnet of moment Qa (g = 2)\n"
            "red dots: the ring of radius a, crossing this plane",
            fontsize=9, va="top", bbox=dict(fc="white", ec="none", alpha=0.85))
    ax.set_aspect("equal"); ax.set_xlim(-3.0, 3.0); ax.set_ylim(-2.4, 3.4)
    ax.set_xlabel("distance from axis"); ax.set_ylabel("height along the spin axis")
    ax.set_title("electric and magnetic fields from one shifted charge")
    ax = axs[1]
    bad = res["C5_bad"]
    names = {"energy from Re(Q/D)^2, not |Q/D|^2": "field energy\nfrom Re(Q/D)²\nnot |Q/D|²",
             "field without the turn: A = Re(Q/D) dt": "field without\nthe turn",
             "charge left at the centre (not shifted)": "charge left\nat the centre"}
    labs = ["shifted charge,\nenergy |Q/D|²\n(Kerr–Newman)"] + [names.get(k, k) for k in bad]
    vals = [res["C4_good"]] + list(bad.values())
    ax.bar(range(len(vals)), vals, color=[GREEN] + [RED] * len(bad), alpha=0.8)
    for i, v in enumerate(vals):
        ax.text(i, v * 3, f"{v:.1e}" if i == 0 else f"{v:.2g}", ha="center", fontsize=8.5)
    ax.set_yscale("log"); ax.set_ylim(1e-52, 1e2)
    ax.set_yticks([10.0**k for k in range(-50, 1, 10)])
    ax.axhline(1e-30, color=GREY, ls=":")
    ax.text(-0.45, 3e-30, "below this line:\nzero to rounding", fontsize=7.5, color=GREY, ha="left", va="bottom")
    ax.set_xticks(range(len(vals))); ax.set_xticklabels(labs, fontsize=8.5)
    ax.set_ylabel("largest Einstein–Maxwell error")
    ax.set_title("charge works only the shifted way\n(M = 1, a = 0.7, Q = 0.5)")
    save(fig, "fig15_charge.png")


# ---------------------------------------------------------------- fig16 two bodies at rest: what adds, and the strut
def fig16():
    M1, c1, M2, c2 = 1.0, 3.0, 0.5, -3.0
    fig, axs = plt.subplots(1, 2, figsize=(13.5, 6.0), gridspec_kw=dict(width_ratios=[1.0, 1.15]))
    ax = axs[0]
    P, Z = np.meshgrid(np.linspace(-5.5, 5.5, 560), np.linspace(-6.5, 6.5, 660))
    rho = np.abs(P)

    def rodpsi(M, c):
        Rp = np.sqrt(rho**2 + (Z - c - M)**2)
        Rm = np.sqrt(rho**2 + (Z - c + M)**2)
        return 0.5 * np.log((Rp + Rm - 2 * M) / (Rp + Rm + 2 * M))
    psi = rodpsi(M1, c1) + rodpsi(M2, c2)
    ax.contour(P, Z, psi, levels=[-1.2, -0.8, -0.55, -0.4, -0.3, -0.22, -0.16, -0.11], colors=BLUE, linewidths=1.1, linestyles="solid")
    ax.plot([0, 0], [c1 - M1, c1 + M1], color="k", lw=6, solid_capstyle="butt")
    ax.plot([0, 0], [c2 - M2, c2 + M2], color="k", lw=6, solid_capstyle="butt")
    ax.plot([0, 0], [c2 + M2, c1 - M1], color=ORANGE, lw=3, ls=(0, (2, 1.2)))
    ax.text(0.35, c1, "body 1, M = 1\n(a rod of length 2M\nin Weyl's coordinates)", fontsize=9, va="center",
            bbox=dict(fc="white", ec="none", alpha=0.85, pad=1))
    ax.text(0.35, c2, "body 2, M = ½", fontsize=9, va="center",
            bbox=dict(fc="white", ec="none", alpha=0.85, pad=1))
    ax.text(-0.3, 0.0, "strut", fontsize=10, color=ORANGE, ha="right", va="center")
    ax.set_aspect("equal"); ax.set_xlim(-5.5, 5.5); ax.set_ylim(-6.5, 8.6)
    ax.text(-5.3, 8.4, "blue: ψ = ½ ln(1 − v²), the log of the clock rate.\nIt is the two bodies' ψ added: exact, and empty\n"
            "space everywhere off the axis (check B1b).\nOn the axis the two bodies combine by Einstein's\n"
            "velocity addition (u₁ + u₂)/(1 + u₁u₂) (check B3).", fontsize=8.8, va="top",
            bbox=dict(fc="white", ec="none", alpha=0.9))
    ax.set_xlabel("distance from the axis"); ax.set_ylabel("height along the axis")
    ax.set_title("two black holes held at rest on one axis")
    ax = axs[1]
    ell = np.linspace(2.0, 12, 400)
    for (m1, m2), col in [((1.0, 0.5), BLUE), ((1.0, 1.0), PURPLE)]:
        ell_ = ell[ell > m1 + m2 + 0.15]
        ax.plot(ell_, m1 * m2 / (ell_**2 - (m1 + m2)**2), color=col, lw=2.2,
                label=f"strut push from the exact field, M₁ = {m1:g}, M₂ = {m2:g}")
        ax.plot(ell_, m1 * m2 / ell_**2, color=col, lw=1.4, ls="--", label=f"Newton M₁M₂/ℓ², M₁ = {m1:g}, M₂ = {m2:g}")
    res = json.load(open(os.path.join(HERE, "results.json")))
    ax.plot([6.0], [res["B4_force"][0]], "o", color=RED, ms=7, zorder=5)
    ax.annotate("measured from the metric in check B4\n(M₁ = 1, M₂ = ½, ℓ = 6)", (6.0, res["B4_force"][0]), (6.8, 0.09),
                fontsize=8.8, arrowprops=dict(arrowstyle="->", color=RED))
    ax.set_ylim(0, 0.3); ax.set_xlim(2, 12)
    ax.set_xlabel("distance ℓ between the centres"); ax.set_ylabel("force")
    ax.legend(fontsize=8.3, loc="upper right")
    ax.set_title("the strut's push is M₁M₂/(ℓ² − (M₁+M₂)²):\nNewton's inverse-square law once they are far apart")
    save(fig, "fig16_two_bodies.png")


# ---------------------------------------------------------------- fig17 why one shear per cell can't hold two bodies
def fig17():
    res = json.load(open(os.path.join(HERE, "results.json")))
    fig, axs = plt.subplots(1, 2, figsize=(13.5, 5.4))
    ax = axs[0]
    spec = res["B1_spec"]
    labs = ["one body\n(any point)", "two bodies:\nbetween them", "two bodies:\noff to the side", "two bodies:\nclose to body 1"]
    vals = [1e-50] + spec
    ax.bar(range(4), vals, color=[GREEN, RED, RED, ORANGE], alpha=0.8)
    for i, v in enumerate(vals):
        ax.text(i, v * 2.5, "0 to rounding" if i == 0 else f"{v:.3g}", ha="center", fontsize=8.5)
    ax.set_yscale("log"); ax.set_ylim(1e-52, 1e3)
    ax.set_yticks([10.0**k for k in range(-50, 1, 10)])
    ax.set_xticks(range(4)); ax.set_xticklabels(labs, fontsize=8.8)
    ax.set_ylabel("how far the curvature is from special\n(0 = two of its three stretch rates equal)")
    ax.set_title("one body's curvature is special, two bodies' is not\n(a single shear per cell only makes special fields)")
    ax = axs[1]
    left = res["B2_left"]
    names = {"add the two metrics": "add the\ntwo metrics", "add the two shears": "add the\ntwo shears",
             "compose the two cell maps": "compose the\ntwo cell maps"}
    labs = ["one body\nalone"] + [names.get(k, k) for k in left]
    vals = [1e-50] + list(left.values())
    ax.bar(range(len(vals)), vals, color=[GREEN] + [RED] * len(left), alpha=0.8)
    for i, v in enumerate(vals):
        ax.text(i, v * 2.5, "0 to rounding" if i == 0 else f"{v:.2g}", ha="center", fontsize=8.5)
    ax.set_yscale("log"); ax.set_ylim(1e-52, 1e-2)
    ax.set_yticks([10.0**k for k in range(-50, -1, 10)])
    ax.axhline(1e-30, color=GREY, ls=":")
    ax.set_xticks(range(len(vals))); ax.set_xticklabels(labs, fontsize=8.8)
    ax.set_ylabel("largest leftover in the empty-space equations")
    ax.set_title("two shears, one per body (M₁ = M₂ = 0.01): every way of\ncombining leaves a leftover that grows as M₁·M₂")
    save(fig, "fig17_one_shear_cannot.png")


def fig18():
    res = json.load(open(os.path.join(HERE, "results.json")))
    fig, axs = plt.subplots(1, 2, figsize=(13.5, 5.6), gridspec_kw={"width_ratios": [1.15, 1]})
    ax = axs[0]
    ax.set_xlim(0, 10); ax.set_ylim(0, 10); ax.axis("off")
    ax.add_patch(Circle((2.2, 5.0), 0.55, fc=BLUE, alpha=0.8)); ax.text(2.2, 5.0, "m₁", ha="center", va="center", color="w")
    ax.add_patch(Circle((7.8, 5.0), 0.4, fc=ORANGE, alpha=0.8)); ax.text(7.8, 5.0, "m₂", ha="center", va="center", color="w")
    ax.annotate("", (2.2, 6.4), (2.2, 5.6), arrowprops=dict(arrowstyle="->", color=BLUE, lw=2))
    ax.annotate("", (7.8, 3.6), (7.8, 4.4), arrowprops=dict(arrowstyle="->", color=ORANGE, lw=2))
    ax.text(2.2, 6.6, "v₁", ha="center", color=BLUE); ax.text(7.8, 3.2, "v₂", ha="center", color=ORANGE)
    for r_ in (1.2, 1.9, 2.6):
        ax.add_patch(Circle((2.2, 5.0), r_, fc="none", ec=BLUE, ls="-", lw=0.8, alpha=0.4))
    boxes = [(0.2, 8.3, "A  each body's own far field, sent out at light speed\n    (where it was a moment ago, not where it is now)", BLUE),
             (0.2, 1.4, "B  clock rates multiply:  rate = (rate from m₁) × (rate from m₂)\n    (the same rule as two bodies at rest, Part 7)", GREEN),
             (0.2, 0.0, "C  each mass is seen through the other's clock:\n    m₁ looks like m₁ (1 − m₂/r₁₂)", PURPLE)]
    for x_, y_, t_, c_ in boxes:
        ax.text(x_, y_, t_, fontsize=9.6, color=c_, va="bottom",
                bbox=dict(boxstyle="round,pad=0.35", fc="white", ec=c_, alpha=0.95))
    ax.set_title("three rules for two moving bodies")
    ax = axs[1]
    drop = res["O1_drop"]
    labs = ["all three\nrules"] + ["drop A\n(no light-speed lag)", "drop B\n(rates added)", "drop C\n(masses not dressed)"]
    vals = [1e-40] + [drop["no light-speed lag"], drop["clock rates added, not multiplied"],
                     drop["masses not seen through the other's clock"]]
    ax.bar(range(4), vals, color=[GREEN, RED, RED, RED], alpha=0.8)
    for i, v in enumerate(vals):
        ax.text(i, v * 3, "0 to rounding\n(40 digits)" if i == 0 else f"{v:.2g}", ha="center", fontsize=8.8)
    ax.set_yscale("log"); ax.set_ylim(1e-42, 1e2)
    ax.set_xticks(range(4)); ax.set_xticklabels(labs, fontsize=8.8)
    ax.set_ylabel("largest mismatch with gravity's published\nslow-orbit field (Blanchet-Faye-Ponsot 1998)")
    ax.set_title("the three rules give the published field exactly;\nleaving any one out breaks it")
    save(fig, "fig18_three_rules.png")


def fig19():
    from orbits import orbit
    res = json.load(open(os.path.join(HERE, "results.json")))
    fig, axs = plt.subplots(1, 3, figsize=(16, 5.4), gridspec_kw={"width_ratios": [1.2, 1, 0.9]}, constrained_layout=True)
    s = orbit(0.25, 50.0, 0.5, 1e6, rstop=6.0)
    t = np.linspace(0, s.t[-1], 60000)
    Y = s.sol(t)
    r = np.hypot(Y[0], Y[1])
    idx = [i for i in range(1, len(r) - 1) if r[i] < r[i - 1] and r[i] <= r[i + 1]]
    ax = axs[0]
    k6 = idx[6]
    ax.plot(Y[0][:k6], Y[1][:k6], color=BLUE, lw=0.8)
    ax.plot(0, 0, "k+")
    ang = np.degrees(np.unwrap(np.arctan2(Y[1][idx[:7]], Y[0][idx[:7]])))
    for i in [0] + idx[:6]:
        ax.plot([0, Y[0][i]], [0, Y[1][i]], color=RED, lw=1.2)
    ax.text(0.02, 0.02, f"red: closest points, turning about {np.diff(ang).mean():.0f}° per lap", color=RED,
            transform=ax.transAxes, fontsize=9)
    ax.set_aspect("equal"); ax.set_xlabel("x / M"); ax.set_ylabel("y / M")
    ax.set_title(f"the first six laps (equal masses, swinging between\n{r[0]:.0f} M and {r[:k6].max():.0f} M): the oval turns round")
    ax = axs[1]
    ax.plot(t, r, color=BLUE, lw=0.6)
    ax.set_xlabel("time / M"); ax.set_ylabel("separation r / M")
    ax.set_title(f"the waves carry energy away: {len(idx)} laps later the orbit\nhas shrunk and rounded off, and the bodies merge")
    ax = axs[2]
    obs = [(0.9983, 0.0016, "Hulse-Taylor pulsar:\norbit shrinking\n(Weisberg & Huang 2016)"),
           (4.77 / res["O4"]["OmB_deg"], 0.66 / res["O4"]["OmB_deg"], "double pulsar B:\nspin wobble\n(Breton et al. 2008)")]
    for i, (v, e, lab) in enumerate(obs):
        ax.errorbar(i, v, yerr=e, fmt="o", color=BLUE, capsize=6)
        ax.text(i + 0.08, v + e + 0.02, lab, fontsize=8.3, ha="center", va="bottom")
    ax.axhline(1, color=GREEN, lw=1.5); ax.text(1.45, 1.01, "prediction", color=GREEN, fontsize=9, ha="right")
    ax.set_xlim(-0.6, 1.6); ax.set_ylim(0.7, 1.45); ax.set_xticks([])
    ax.set_ylabel("observed / predicted")
    ax.set_title("against real binary pulsars\n(predicted: dP/dt = −2.402e-12; 5.07°/yr)")
    save(fig, "fig19_orbit_decays.png")


def fig20():
    from orbits import inspiral_spins, spin_angles, unit
    data = json.load(open(os.path.join(HERE, "orbits.json")))
    fig, axs = plt.subplots(1, 3, figsize=(15.5, 5.0))
    ax = axs[0]
    s = inspiral_spins(0.7, 0.8, 0.8, unit(30, 0), unit(120, 90), r0=200.0, r1=10.0, npts=1500)
    A = np.array([spin_angles(s.y[:, k]) for k in range(s.y.shape[1])])
    for j, (lab, c) in enumerate([("between the two spins", PURPLE), ("heavier spin to orbit axis", BLUE),
                                  ("lighter spin to orbit axis", ORANGE), ("between their shadows on the orbit plane", GREY)]):
        ax.plot(s.t, A[:, j], color=c, lw=1.4 if j < 3 else 0.6, label=lab)
    ax.set_xscale("log"); ax.invert_xaxis(); ax.set_ylim(0, 185)
    ax.set_xlabel("separation r / M  (orbit shrinking →)"); ax.set_ylabel("angle (degrees)")
    ax.legend(fontsize=7.8, loc="upper left")
    ax.set_title("one pair (mass ratio 0.7, spins 0.8, from 200 M):\nthe spins wobble round the orbit's axis as it shrinks")
    st, en = np.array(data["start"]), np.array(data["end"])
    ax = axs[1]
    bins = np.linspace(-1, 1, 11)
    c0, c1 = np.cos(np.radians(st[:, 0])), np.cos(np.radians(en[:, 0]))
    ax.hist(c0, bins, color=GREY, alpha=0.5, label=f"start (r = 1000 M), mean {c0.mean():+.2f}")
    ax.hist(c1, bins, histtype="step", color=PURPLE, lw=2, label=f"end (r = 10 M), mean {c1.mean():+.2f}")
    ax.set_xlabel("cosine of the angle between the two spins\n(+1 = lined up, −1 = opposite)"); ax.set_ylabel("pairs (of 80)")
    ax.legend(fontsize=8.2, loc="upper center"); ax.set_ylim(0, max(np.histogram(c0, bins)[0].max(), np.histogram(c1, bins)[0].max()) + 6)
    ax.set_title("80 pairs started pointing every which way:\non average they do not line up")
    ax = axs[2]
    lo = st[:, 1] < st[:, 2]
    b2 = np.linspace(0, 180, 10)
    ax.hist(en[lo, 3], b2, color=BLUE, alpha=0.6,
            label=f"heavier spin starts nearer the axis ({lo.sum()} pairs; median {np.median(en[lo, 3]):.0f}°)")
    ax.hist(en[~lo, 3], b2, color=ORANGE, alpha=0.6,
            label=f"heavier spin starts farther ({(~lo).sum()} pairs; median {np.median(en[~lo, 3]):.0f}°)")
    for m_, c_ in ((np.median(en[lo, 3]), BLUE), (np.median(en[~lo, 3]), ORANGE)):
        ax.axvline(m_, color=c_, ls="--")
    ax.set_xlabel("angle between the spins' shadows on the orbit plane at r = 10 M (degrees)")
    ax.set_ylabel("pairs"); ax.legend(fontsize=7.8, loc="upper center")
    ax.set_ylim(0, ax.get_ylim()[1] * 1.35)
    ax.set_title("but they sort into two families\n(near 0°, or near 180°: spin-orbit resonances)")
    save(fig, "fig20_spins.png")


def fig21():
    Tsun = 4.925490947e-6
    Mc = 30 * Tsun
    f0, f1 = 35.0, 150.0
    tau = lambda f: 5 / 256 * Mc**(-5 / 3) * (np.pi * f)**(-8 / 3)  # noqa: E731
    T = tau(f0) - tau(f1)
    tt = np.linspace(0, T, 20000)
    tc_left = tau(f1) + T - tt                    # time left to the formal coalescence
    f = (5 / (256 * tc_left))**(3 / 8) * Mc**(-5 / 8) / np.pi
    phase = 2 * np.pi * np.cumsum(f) * (tt[1] - tt[0])
    h = (f / f1)**(2 / 3) * np.cos(phase)
    Mf, af = 62 * 1.09 * Tsun, 0.67
    fr = (1.5251 - 1.1568 * (1 - af)**0.1292) / (2 * np.pi * Mf)
    Q = 0.7 + 1.4187 * (1 - af)**(-0.4990)
    td = Q / (np.pi * fr)
    gap = 0.006
    tr = np.linspace(0, 6 * td, 2000)
    hr = np.exp(-tr / td) * np.cos(2 * np.pi * fr * tr)
    fig, axs = plt.subplots(2, 1, figsize=(13, 7.2), sharex=True, gridspec_kw={"height_ratios": [1.3, 1]})
    ax = axs[0]
    ax.plot(tt, h, color=BLUE, lw=0.9, label="run-up: two bodies spiralling in (leading-order quadrupole formula)")
    ax.axvspan(T, T + gap, color=GREY, alpha=0.25)
    ax.text(T + gap / 2, 1.95, "merger: needs a full supercomputer\nsimulation (not done here)", ha="center", va="top", fontsize=8.5,
            bbox=dict(boxstyle="round,pad=0.3", fc="white", ec=GREY))
    ax.plot(T + gap + tr, hr, color=RED, lw=1.2, label=f"ring-down: the merged spinning hole's own tone ({fr:.0f} Hz)")
    ax.set_ylim(-1.4, 2.05); ax.set_xlim(T - 0.08, T + gap + 6 * td + 0.004)
    ax.set_ylabel("wave strain (shape only;\nheights matched by hand)")
    ax.legend(fontsize=8.5, loc="lower left")
    ax.set_title("the wave from a GW150914-like collision (30 solar-mass chirp mass; final hole 62 solar masses, spin 0.67)")
    ax = axs[1]
    ax.plot(tt, f, color=BLUE, label=f"run-up frequency: 35 → 150 Hz in {T:.2f} s, about 9 cycles")
    ax.axvspan(T, T + gap, color=GREY, alpha=0.25)
    ax.plot([T + gap, T + gap + 6 * td], [fr, fr], color=RED, lw=2, label=f"ring-down {fr:.0f} Hz, fades in {1e3 * td:.0f} ms")
    ax.set_ylim(20, 300); ax.set_xlabel("time (s) since the wave passed 35 Hz (last 0.08 s of the run-up shown)")
    ax.set_ylabel("frequency (Hz)"); ax.legend(fontsize=8.5, loc="upper left")
    save(fig, "fig21_wave.png")


def fig22():
    from fractions import Fraction as Fr
    nj = json.load(open(os.path.join(HERE, "next_order.json")))

    def truth3(a):
        b = 1 - a
        return -b * (54*a**5 + 224*a**4*b + 356*a**3*b**2 + 284*a**2*b**3 + 124*a*b**4 + 27*b**5) / 16

    def miss(a):
        return a * (1 - a) * (23*a**2 - 50*a + 36) / 3
    fig, axs = plt.subplots(1, 3, figsize=(16, 5.4), gridspec_kw={"width_ratios": [0.85, 1.1, 1.25]})
    # left: level by level, m1 = 1/3
    ax = axs[0]
    row = nj["two_pn"]["1/3"]
    rules = [float(Fr(c)) for c in row["rules"]][1:]
    tru = [float(Fr(c)) for c in row["truth"]][1:]
    xs_ = np.arange(3)
    ax.bar(xs_ - 0.18, tru, 0.36, color=GREY, label="general relativity (first law)")
    ax.bar(xs_ + 0.18, rules, 0.36, color=[GREEN, GREEN, RED], label="three rules")
    for i in range(3):
        ok = abs(rules[i] - tru[i]) < 1e-12
        ax.text(i, max(rules[i], tru[i], 0) + 0.12, "same" if ok else "off by\n" + f"{rules[i] - tru[i]:.2f}",
                ha="center", fontsize=9, color=GREEN if ok else RED)
    ax.axhline(0, color="k", lw=0.6)
    ax.set_xticks(xs_); ax.set_xticklabels(["x\n(Newton)", "x²\n(1PN)", "x³\n(2PN)"])
    ax.set_ylim(-1.6, 1.0); ax.set_ylabel("coefficient in body 1's clock rate z₁")
    ax.legend(fontsize=8.5, loc="lower left")
    ax.set_title("one level at a time (body 1 = ⅓ of the mass)", fontsize=10.5)
    # middle: the x^3 coefficient across mass ratio
    ax = axs[1]
    a = np.linspace(0, 1, 301)
    ax.plot(a, truth3(a), color=GREY, lw=2.2, label="general relativity")
    ax.plot(a, truth3(a) + miss(a), color=RED, lw=1.6, ls="--", label="three rules: fit through the dots")
    for k, v in nj["two_pn"].items():
        a_ = float(Fr(k)); y_ = float(Fr(v["rules"][3]))
        held = k == "2/5"
        ax.plot(a_, y_, "o", ms=7, mfc="white" if held else RED, mec=RED, zorder=5)
    ax.plot([], [], "o", mfc=RED, mec=RED, label="three rules, computed")
    ax.plot([], [], "o", mfc="white", mec=RED, label="computed after the fit (m₁ = 2/5): on the curve")
    ax.annotate("a body too light to pull back\n(test mass): exact", (0, truth3(0.0)), (0.22, -0.8), fontsize=8.5,
                arrowprops=dict(arrowstyle="->", color=GREEN), color=GREEN)
    ax.set_xlabel("m₁ / m,  with m = m₁ + m₂   (ν = m₁m₂/m²,  Δ = (m₂ − m₁)/m)"); ax.set_ylabel("x³ coefficient of z₁ (2PN)")
    ax.set_ylim(-1.9, 0.6)
    ax.legend(fontsize=8.2, loc="lower right")
    ax.set_title("at 2PN the rules miss by  ν(15/2 + 9Δ/2 − 23ν/3)", fontsize=10.5)
    # right: what each body feels, and what is missing
    ax = axs[2]
    ax.set_xlim(0, 10); ax.set_ylim(0, 10); ax.axis("off")
    ax.add_patch(Circle((2.0, 6.0), 0.55, fc=BLUE, alpha=0.85)); ax.text(2.0, 6.0, "m₁", ha="center", va="center", color="w")
    ax.add_patch(Circle((8.0, 6.0), 0.4, fc=ORANGE, alpha=0.85)); ax.text(8.0, 6.0, "m₂", ha="center", va="center", color="w")
    ax.add_patch(FancyBboxPatch((3.3, 4.9), 3.4, 2.2, boxstyle="round,pad=0.1", fc=PURPLE, alpha=0.12, ec=PURPLE, ls="--"))
    ax.text(5.0, 6.0, "field energy\nbetween them\n(missing?)", ha="center", va="center", fontsize=9, color=PURPLE)
    ax.annotate("", (2.6, 7.3), (7.6, 7.3), arrowprops=dict(arrowstyle="->", color=GREEN, lw=1.8))
    ax.text(5.0, 7.6, "m₁ feels m₂'s own field, lagged and dressed: kept", ha="center", fontsize=8.8, color=GREEN)
    ax.text(5.0, 8.9, "rule B (rates multiply) shapes the field between the bodies,\nbut drops out of what each body feels",
            ha="center", fontsize=8.8, color=GREY)
    boxes = [(0.1, 2.35, "1PN: rules = general relativity, and the test-mass limit is exact at 2PN", GREEN),
             (0.1, 1.15, "2PN: missing piece, one of (or both):\n  • gravity's own field energy between the bodies acting as a source\n  • the extra field of a body that is accelerating", RED)]
    for x_, y_, t_, c_ in boxes:
        ax.text(x_, y_, t_, fontsize=8.8, color=c_, va="bottom",
                bbox=dict(boxstyle="round,pad=0.35", fc="white", ec=c_, alpha=0.95))
    ax.set_title("why: the bodies only feel each other's own fields", fontsize=10.5)
    save(fig, "fig22_next_order.png")


# ---------------------------------------------------------------- fig23 a fourth rule: the field's own energy (round 14)
def fig23():
    from fractions import Fraction as Fr
    fj = json.load(open(os.path.join(HERE, "fourth_rule.json")))
    pr = fj["per_ratio"]
    fig, axs = plt.subplots(1, 3, figsize=(16.5, 5.4), gridspec_kw={"width_ratios": [0.95, 1.15, 1.2]})
    # left: the stages at m1 = 1/3
    ax = axs[0]
    row = pr["1/3"]
    vals = [float(Fr(row[k])) for k in ("miss_rules", "miss_closed", "miss_with_f", "miss_with_own_square")]
    labels = ["three\nrules", "+ field\nenergy,\nclosed form", "+ field on\nfield\n(4 integrals)",
              "+ own field\nsquared\n(round 15)"]
    cols = [RED, ORANGE, BLUE, GREEN]
    ax.bar(range(4), vals, 0.6, color=cols)
    for i, v in enumerate(vals):
        ax.text(i, v + 0.04, f"{v:.2f}" if i < 3 else "0", ha="center", fontsize=10, color=cols[i])
    ax.set_xticks(range(4)); ax.set_xticklabels(labels, fontsize=8.2)
    ax.set_ylabel("miss in the x³ (2PN) coefficient of z₁"); ax.set_ylim(0, 1.85)
    ax.text(2, 0.3, "92%\nclosed", ha="center", fontsize=9.5, color=BLUE)
    ax.set_title("what each step closes (body 1 = ⅓ of the mass)", fontsize=10.5)
    # middle: the stages across mass ratio
    ax = axs[1]
    a = np.linspace(0, 1, 301)
    curves = [(a*(1 - a)*(23*a**2 - 50*a + 36)/3, RED, "three rules: m₁m₂(23m₁² − 50m₁ + 36)/3", "miss_rules"),
              (a*(1 - a)*(5 - 10*a + 3*a**2), ORANGE, "+ closed-form field energy: m₁m₂(5 − 10m₁ + 3m₁²)", "miss_closed"),
              (a*(1 - a)**2*(3 - a)/3, BLUE, "+ field on field: m₁m₂²(3 − m₁)/3", "miss_with_f")]
    for y, c, lab, key in curves:
        ax.plot(a, y, color=c, lw=1.8, label=lab)
        for k, v in pr.items():
            held = k == "2/5"
            ax.plot(float(Fr(k)), float(Fr(v[key])), "o", ms=6.5, mfc="white" if held else c, mec=c, zorder=5)
    ax.plot(a, 0*a, color=GREEN, lw=1.8, label="+ own field squared (derived): 0 at every ratio")
    for k, v in pr.items():
        ax.plot(float(Fr(k)), float(Fr(v["miss_with_own_square"])), "o", ms=6.5,
                mfc="white" if k == "2/5" else GREEN, mec=GREEN, zorder=5)
    ax.plot([], [], "o", mfc="white", mec="k", label="m₁ = 2/5: never used in a fit")
    ax.axhline(0, color="k", lw=0.5)
    ax.set_xlabel("m₁ / m"); ax.set_ylabel("miss in the x³ coefficient of z₁")
    ax.set_ylim(-0.25, 2.55); ax.legend(fontsize=8, loc="upper right")
    ax.set_title("the same steps at every mass ratio", fontsize=10.5)
    # right: the picture
    ax = axs[2]
    ax.set_xlim(0, 10); ax.set_ylim(0, 10); ax.axis("off")
    ax.add_patch(Circle((2.0, 6.3), 0.55, fc=BLUE, alpha=0.85)); ax.text(2.0, 6.3, "m₁", ha="center", va="center", color="w")
    ax.add_patch(Circle((8.0, 6.3), 0.4, fc=ORANGE, alpha=0.85)); ax.text(8.0, 6.3, "m₂", ha="center", va="center", color="w")
    for k in range(9):
        xx = 3.0 + 0.5*k; yy = 6.3 + 0.9*np.sin(k*0.9)
        ax.add_patch(Circle((xx, yy), 0.12 + 0.05*(k % 3), fc=PURPLE, alpha=0.35))
    ax.text(5.0, 8.15, "the field's energy between the bodies is itself a source", ha="center", fontsize=9, color=PURPLE)
    ax.add_patch(Circle((8.0, 6.3), 0.75, fc="none", ec=GREEN, ls="--", lw=1.6))
    ax.add_patch(Circle((2.0, 6.3), 0.9, fc="none", ec=GREEN, ls="--", lw=1.6))
    ax.annotate("", xy=(6.9, 6.3), xytext=(7.55, 6.3), arrowprops=dict(arrowstyle="-|>", color=GREEN, lw=1.6))
    ax.annotate("", xy=(3.15, 6.3), xytext=(2.5, 6.3), arrowprops=dict(arrowstyle="-|>", color=GREEN, lw=1.6))
    ax.text(5.0, 4.2, "each body speeds up toward the other (arrows); its field shifts with it (dashed),\n"
            "and that shift times its own field is the last piece", ha="center", fontsize=8.3, color=GREEN)
    boxes = [(0.1, 2.6, "known: general relativity's 2PN field (Blanchet, Faye & Ponsot 1998);\nFock's function ln(r₁ + r₂ + r₁₂); Hadamard's finite part", GREY),
             (0.1, 1.2, "ours: rules + field energy rebuild the 2PN field piece by piece;\nthe last piece is a body's own field squared: the gap is 0 at\nevery ratio, and the orbit size matches the known 2PN law", BLUE)]
    for x_, y_, t_, c_ in boxes:
        ax.text(x_, y_, t_, fontsize=8.6, color=c_, va="bottom",
                bbox=dict(boxstyle="round,pad=0.35", fc="white", ec=c_, alpha=0.95))
    ax.set_title("a fourth rule: gravity's own energy gravitates, gap closed", fontsize=10.5)
    save(fig, "fig23_fourth_rule.png")

# ---------------------------------------------------------------- fig24 Kerr's orbit factor split, two ways (round 16)
def fig24():
    import mpmath as mp_
    import kerr_split as ks
    kj = json.load(open(os.path.join(HERE, "kerr_split.json")))
    fig, axs = plt.subplots(2, 2, figsize=(14.5, 10.5))
    # (a) the double copy: which piece of (k.p)^2 makes each number, and the single copy beside it
    ax = axs[0, 0]
    groups = ["u  (6 in Schwarzschild)", r"$\chi\, u^{3/2}$  (spin)", "χ² u²  (shape)"]
    x0 = np.arange(3)
    dbl = [[("P·P", -2, ORANGE), ("2E·P", -4, RED)], [("2E·S", 8, BLUE)], [("E·E", -3, PURPLE)]]
    sgl = [[("P", -1, RED)], [("S", 4, BLUE)], [("E", -3, PURPLE)]]
    for i in range(3):
        b = 0
        for lab, v, c in dbl[i]:
            ax.bar(x0[i] + 0.2, v, 0.36, bottom=b, color=c, alpha=0.9)
            ax.text(x0[i] + 0.2, b + v/2, f"{lab}\n{v:+d}", ha="center", va="center", color="w", fontsize=9)
            b += v
        for lab, v, c in sgl[i]:
            ax.bar(x0[i] - 0.2, v, 0.36, color="white", ec=c, hatch="//", lw=1.5)
            ax.text(x0[i] - 0.2, v + (0.35 if v > 0 else -0.35), f"{lab}: {v:+d}", ha="center",
                    va="bottom" if v > 0 else "top", color=c, fontsize=9)
    ax.axhline(0, color="k", lw=0.6)
    ax.set_xticks(x0); ax.set_xticklabels(groups, fontsize=9.5); ax.set_ylim(-7.5, 9.8)
    ax.set_ylabel(r"number in $T = 1 - 6u + 8\chi u^{3/2} - 3\chi^2u^2$")
    from matplotlib.patches import Patch
    ax.legend(handles=[Patch(fc="white", ec="k", hatch="//", label="single copy: −(2M/r)(E + P + S)"),
                       Patch(fc=GREY, label="double copy (gravity): −(2M/r)(E + P + S)²")], fontsize=8.5, loc="upper left")
    ax.text(1.45, -7.2, "spin: 8 = 2 × 4 (first order, doubled)\nSchwarzschild's 4 = 2² × 1 (second order, quadrupled)\n"
            "shape: the same −3 in both: Newton's pull of the shifted source", ha="center", va="bottom", fontsize=8.3)
    ax.set_title("(a) way 1: each number comes from one piece of the field", fontsize=10.5)
    # (b) the rules: dial the push and the shape
    ax = axs[0, 1]
    kap = np.linspace(0, 1.25, 50)
    ax.plot(kap, 8*kap, color=BLUE, lw=2)
    ax.plot(kap, -3*kap, color=PURPLE, lw=2)
    pn = kj["pn"]
    pts = [(0, float(sp_frac(pn["no push (kappa=0)"]["chi_u15"])), BLUE), (0.5, float(sp_frac(pn["half push (kappa=1/2)"]["chi_u15"])), BLUE),
           (1, float(sp_frac(pn["Kerr (beta=gamma=kappa=q=1)"]["chi_u15"])), BLUE),
           (0, float(sp_frac(pn["no shape (q=0)"]["chi2_u2"])), PURPLE), (1, float(sp_frac(pn["Kerr (beta=gamma=kappa=q=1)"]["chi2_u2"])), PURPLE),
           ]
    for xx, yy, c in pts:
        ax.plot(xx, yy, "o", ms=8, color=c, zorder=5)
    ax.plot(1.25, -3*1.25, alpha=0)
    ax.axvline(1, color=GREY, ls=":", lw=1)
    ax.annotate("rule A: a moving body's field is its rest field boosted:\n"
                "h₀ᵢ = −(2U [time part] + 2U [space part]) v\n→ κ = 1: spin number 4 + 4 = 8",
                xy=(1, 8), xytext=(0.03, 8.4), fontsize=8.6, color=BLUE, arrowprops=dict(arrowstyle="->", color=BLUE))
    ax.annotate("Kerr's shape (quadrupole Q = −Ma²): q = 1\nthe rules don't fix q; the shift by i·a does",
                xy=(1, -3), xytext=(0.03, -4.6), fontsize=8.6, color=PURPLE, arrowprops=dict(arrowstyle="->", color=PURPLE))
    ax.axhline(0, color="k", lw=0.6); ax.set_xlabel("dial value (κ or q)"); ax.set_ylabel("number in T")
    ax.text(0.55, 3.0, "spin number = 8κ", color=BLUE, fontsize=9.5)
    ax.text(0.55, -1.0, "shape number = −3q", color=PURPLE, fontsize=9.5)
    ax.set_ylim(-5, 10.5)
    ax.set_title("(b) way 2: the push rule makes the spin number, the shape makes the −3", fontsize=10.5)
    # (c) the picture
    ax = axs[1, 0]
    ax.set_xlim(-5, 5); ax.set_ylim(-4.2, 4.2); ax.set_aspect("equal"); ax.axis("off")
    th = np.linspace(0, 2*np.pi, 300)
    ax.plot(1.0*np.cos(th), 0.35*np.sin(th), color=GREY, lw=1.5)
    ax.text(0, -1.45, "shifted source:\nring of radius a", ha="center", fontsize=8.3, color=GREY)
    ax.plot(1.414*np.cos(th), 0.5*np.sin(th), color=PURPLE, lw=1.6, ls="--")
    ax.text(0, 0.75, "same mass, spin and shape:\na ring of radius √2·a\nturning at c/√2", ha="center", va="bottom", fontsize=8.3, color=PURPLE)
    ax.annotate("", xy=(-0.5, -0.48), xytext=(0.5, -0.48), arrowprops=dict(arrowstyle="-|>", color=PURPLE))
    rr_ = 3.4 + 0.25*np.cos(0.93*th*3)
    ax.plot(rr_*np.cos(th*3)[:300], 0.75*rr_*np.sin(th*3)[:300], color=BLUE, lw=1.3)
    ax.text(-4.8, 3.3, "orbit: the closest point turns each lap", fontsize=8.6, color=BLUE)
    ax.text(-4.9, -3.3, "E·E: Newton's pull of the shifted source → −3χ²u²\n2E·S: energy × the rays' twist → +8χu³ᐟ²\n"
            "P·P + 2E·P: energy × radial push, and push² → −6u", fontsize=8.6)
    ax.set_title("(c) where each number lives", fontsize=10.5)
    # (d) what doesn't carry over: the ISCO landmark 2v = c
    ax = axs[1, 1]
    chis = np.linspace(-0.98, 0.9999, 240)
    vs = [abs(float(ks.zamo_speed(ks.isco(mp_.mpf(c)), mp_.mpf(c)))) for c in chis]
    ax.plot(chis, vs, color=RED, lw=2, label="speed at the innermost stable orbit,\nseen by the local frame turning with the hole")
    ax.axhline(0.5, color=GREY, ls="--", lw=1.2, label="c/2: Schwarzschild's landmark (2v = c)")
    ax.set_xlabel("spin χ = a/M (negative: orbit against the spin)"); ax.set_ylabel("v / c")
    ax.legend(fontsize=8.5, loc="lower right"); ax.set_ylim(0.38, 0.68)
    ax.text(-0.95, 0.64, "the product form 1 − 6u = (1 − 2u)(1 − (2v)²)\nalso fails: off by 0.08 at r = 10M, χ = ½",
            fontsize=8.6, color=RED)
    ax.set_title("(d) what doesn't carry over to a spinning hole", fontsize=10.5)
    fig.tight_layout()
    save(fig, "fig24_kerr_split.png")

# ---------------------------------------------------------------- fig25 why the tower of shapes leaves one -3 (round 17)
def fig25():
    from fractions import Fraction as Fr
    tj = json.load(open(os.path.join(HERE, "tower.json")))
    fig, axs = plt.subplots(2, 2, figsize=(14.5, 10.5))
    # (a) the staircase, as a table of numbers
    ax = axs[0, 0]
    rows = [("shifted source, cut after its quadrupole", tj["staircase_in_r"]["up to (a/R)^2"]),
            ("... after its 2⁴-pole (hexadecapole)", tj["staircase_in_r"]["up to (a/R)^4"]),
            ("... after its 2⁶-pole", tj["staircase_in_r"]["up to (a/R)^6"]),
            ("... after its 2⁸-pole", tj["staircase_in_r"]["up to (a/R)^8"]),
            ("the whole shifted source", tj["staircase_in_r"]["up to (a/R)^10"]),
            ("rival: quadrupole only", tj["rivals_in_r"]["quadrupole only"]),
            ("rival: ring of radius √2·a", tj["rivals_in_r"]["ring of radius sqrt(2) a"])]
    ax.set_xlim(-3.6, 5); ax.set_ylim(-1.7, len(rows) + 0.4); ax.axis("off")
    for j, lab in enumerate(["(a/r)²", "(a/r)⁴", "(a/r)⁶", "(a/r)⁸", "(a/r)¹⁰"]):
        ax.text(j + 0.5, len(rows) - 0.1, lab, ha="center", fontsize=9.5, weight="bold")
    for i, (lab, vals) in enumerate(rows):
        y = len(rows) - 1 - i - 0.5
        if i == 5:
            ax.axhline(y + 0.5, xmin=0.02, xmax=0.98, color=GREY, lw=0.8)
        ax.text(-0.1, y, lab, ha="right", va="center", fontsize=8.8, color=BLUE if i < 5 else RED)
        for j, v in enumerate(vals):
            f = Fr(v)
            c = GREEN if f == 0 else (PURPLE if j == 0 else (ORANGE if i < 5 else RED))
            ax.add_patch(FancyBboxPatch((j + 0.08, y - 0.38), 0.84, 0.76, boxstyle="round,pad=0.02",
                                        fc=c, alpha=0.18 if f != 0 else 0.35, ec=c))
            ax.text(j + 0.5, y, "0" if f == 0 else (f"{float(f):+.4g}" if f.denominator > 16 else str(f)),
                    ha="center", va="center", fontsize=8.6, color="k")
    ax.text(-3.4, -1.62, "orbit factor T − 1 in the shift's radius r: each new shape wipes out the leftover\n"
            "of the ones below it (orange → green); a body with only Kerr's quadrupole never stops", fontsize=8.6)
    ax.set_title("(a) the staircase: Kerr's shapes cancel one another", fontsize=10.5)
    # (b) the same T in flat radius R and in the shift's radius r
    ax = axs[0, 1]
    q = np.linspace(0, 0.9, 200)            # a/R
    Tflat = (1 - 4*q**2)/(1 - q**2)
    ax.plot(q, Tflat, color=BLUE, lw=2.4, label="shifted source: exact T")
    part = np.ones_like(q)
    for n in range(1, 6):
        part = part - 3*q**(2*n)
        ax.plot(q, part, color=BLUE, lw=0.9, ls=":", alpha=0.8)
    ax.plot([], [], color=BLUE, lw=0.9, ls=":", label="1 − 3(a/R)² − 3(a/R)⁴ − …, cut after 1 to 5 terms")
    ax.set_xlabel("a / R   (R: flat radius)"); ax.set_ylabel("orbit factor T")
    ax.axhline(0, color="k", lw=0.5); ax.set_ylim(-1.6, 1.05)
    ax2 = ax.twiny()
    rq = q/np.sqrt(1 - q**2)
    ax2.set_xlim(ax.get_xlim()); tick_r = np.array([0, 0.25, 0.5, 0.75, 1.0, 1.5])
    ax2.set_xticks(tick_r/np.sqrt(1 + tick_r**2)); ax2.set_xticklabels([f"{t:g}" for t in tick_r])
    ax2.set_xlabel("a / r   (r = √(R² − a²): distance to the shifted source)", fontsize=9)
    ax.text(0.04, -0.9, "in flat radius every shape adds −3,\nin the shift's radius it is one term:\n"
            "T = 1 − 3a²/r²  and  Ω² = M/r³ (Kepler)", fontsize=9.2, color=BLUE)
    ax.legend(fontsize=8.3, loc="lower left")
    ax.set_title("(b) one curve, two radii", fontsize=10.5, pad=34)
    # (c) the picture
    ax = axs[1, 0]
    ax.set_xlim(-5, 5); ax.set_ylim(-4.3, 4.3); ax.set_aspect("equal"); ax.axis("off")
    th = np.linspace(0, 2*np.pi, 300)
    A_, R_ = 1.1, 3.4
    cx, cy = -1.6, 0.4
    ax.plot(cx + A_*np.cos(th), cy + A_*np.sin(th), color=GREY, lw=1.6)
    ax.text(cx, cy - A_ - 0.45, "Kerr's ring, radius a (seen from above)", ha="center", fontsize=8.3, color=GREY)
    tx, ty = cx + A_**2/R_, cy + A_*np.sqrt(R_**2 - A_**2)/R_
    ax.plot([cx, cx + R_], [cy, cy], color=PURPLE, lw=1.3); ax.text(cx + R_/2, cy - 0.35, "R (flat radius)", color=PURPLE, fontsize=9, ha="center")
    ax.plot([cx, tx], [cy, ty], color=GREY, lw=1.0, ls=":"); ax.text((cx + tx)/2 - 0.25, (cy + ty)/2, "a", color=GREY)
    ax.plot([tx, cx + R_], [ty, cy], color=BLUE, lw=2.0); ax.text((tx + cx + R_)/2 + 0.1, (ty + cy)/2 + 0.15, "r", color=BLUE, fontsize=12)
    ax.plot(cx + R_, cy, "o", color="k"); ax.text(cx + R_ + 0.15, cy - 0.1, "orbiting body", fontsize=8.3)
    ax.text(-4.8, 3.6, "the shift's radius r = √(R² − a²): the distance along a tangent to Kerr's ring", fontsize=8.6)
    ax.text(-4.8, 2.9, "in r, Newton's orbit is Kepler's: Ω² = M/r³ exactly,", fontsize=8.8, color=BLUE)
    ax.text(-4.8, 2.35, "and the only shape number left is −3a²/r²", fontsize=8.8, color=PURPLE)
    ax.text(-4.8, -2.6, "gravity, in the same r:\n• radial motion is Carter's quartic: only a² (shape) and L − aE (twist)\n"
            "• period = Kepler's in r + the clock gap: 1/Ω = r³ᐟ²/√M + a\n• so T stops at −3χ²u²", fontsize=8.6)
    ax.set_title("(c) why: Kepler's law holds in the shift's radius", fontsize=10.5)
    # (d) the measurable version for Kerr: against the period
    ax = axs[1, 1]
    chi = 0.7
    xs = np.linspace(0.005, 0.12, 200)
    def T_of_u(u):
        return 1 - 6*u + 8*chi*u**1.5 - 3*chi**2*u**2
    u_from_raw = xs*(1 - chi*xs**1.5)**(-2/3)
    for K, ls in ((4, "--"), (6, "-.")):
        coeffs = {2: -6, 3: 8*chi, 4: -3*chi**2, 5: -4*chi, 6: 8*chi**2, 7: -4*chi**3, 8: -10*chi**2/3}
        ser = 1 + sum(c*xs**(k/2) for k, c in coeffs.items() if k <= K)
        ax.plot(xs, ser - T_of_u(u_from_raw), ls, color=RED, lw=1.4,
                label=f"raw frequency x: the finite form misses" if K == 4 else "… even with 2 more terms")
    ax.axhline(0, color=GREEN, lw=2.2, label="clock gap taken off (1/Ω̃ = 1/Ω − a): the finite form is exact")
    ax.set_xlabel("x = (MΩ)²ᐟ³, a measurable orbital frequency (χ = 0.7)"); ax.set_ylabel("finite form − true T")
    ax.legend(fontsize=8.3, loc="lower left")
    ax.set_title("(d) measurable: subtract Part 1's clock gap 2πa from the period", fontsize=10.5)
    fig.tight_layout()
    save(fig, "fig25_tower.png")


def fig26():
    tj = json.load(open(os.path.join(HERE, "tilt.json")))
    fig, axs = plt.subplots(2, 2, figsize=(14.5, 10.5))
    av, rv = 0.8, 3.0
    # (a) the picture: one tilted orbit flown around the shifted source, seen from the side
    ax = axs[0, 0]
    ax.set_aspect("equal"); ax.set_xlim(-4.6, 4.6); ax.set_ylim(-4.3, 4.9); ax.axis("off")
    th = np.linspace(0, 2*np.pi, 400)
    ax.plot(np.sqrt(rv**2 + av**2)*np.sin(th), rv*np.cos(th), color=GREY, lw=1.0, ls="--")
    rows = tj["spheres"]["a=0.8, r=3.0"]["rows"]
    row = [x for x in rows if x["L_frac"] == "0.7"][0]
    c = np.sin(np.radians(row["tilt_deg"]))
    L = float(row["L_frac"])*tj["spheres"]["a=0.8, r=3.0"]["L_equator"]
    rho = np.sqrt(rv**2 + av**2)*np.sqrt(1 - c**2); zz = rv*c

    def f(t, y):
        x, yy, z, vx, vy, vz = y
        d = np.sqrt(complex(x*x + yy*yy + (z - 1j*av)**2))
        g = -np.array([x, yy, z - 1j*av])/d**3
        return [vx, vy, vz, g[0].real, g[1].real, g[2].real]
    sol = solve_ivp(f, [0, 900], [rho, 0, zz, 0, L/rho, 0], rtol=1e-10, atol=1e-10, max_step=0.05)
    ax.plot(sol.y[0], sol.y[2], color=BLUE, lw=0.5, alpha=0.45)
    ax.plot([-av, av], [0, 0], color="k", lw=3)
    ax.text(0, -0.45, "Kerr's ring\n(edge on)", ha="center", va="top", fontsize=8.3)
    ax.text(-4.5, 3.75, "an orbit tilted 43°, flown in ordinary space around the shifted source:\n"
            "it stays on one oblate shell (dashed), a fixed distance r = 3 from the shifted source,\n"
            "while its tilted plane slowly turns, so the path fills a band",
            fontsize=8.7)
    ax.text(-4.5, -4.1, "on that shell, only the equator keeps Kepler's period (Ω² = M/r³);\n"
            "tilted, the orbit feels the whole tower of shapes again", fontsize=8.7, color=BLUE)
    ax.set_title("(a) a tilted orbit (side view, a = 0.8, r = 3)", fontsize=10.5, pad=12)
    # (b) the rates against tilt
    ax = axs[0, 1]
    for key, ls in (("a=0.8, r=3.0", "-"), ("a=0.8, r=5.0", "--")):
        d = tj["spheres"][key]; rr = d["rows"]
        tilt = [x["tilt_deg"] for x in rr]
        ax.plot(tilt, [x["Om_th2"] for x in rr], ls, color=PURPLE, lw=2,
                label=f"up and down, Ω_θ² r³/M  ({key.split(', ')[1]})")
        ax.plot(tilt, [x["Om_ph2"] for x in rr], ls, color=ORANGE, lw=2,
                label=f"around the axis, Ω_φ² r³/M  ({key.split(', ')[1]})")
        ax.plot([0], [d["equator_vertical"]], "o", color=PURPLE, ms=6)
        ax.plot([90], [d["polar_exact"]], "s", color=PURPLE, ms=6)
    ax.axhline(1, color=GREEN, lw=2.2, label="Kepler's law in r")
    ax.text(2, 1.225, "1 + 3a²/r²  (exact)", fontsize=8.6, color=PURPLE)
    ax.text(58, 0.835, "pole: (π/2)² / (2K − E)², k = a/r\n= 1 − 3k²/2 + 33k⁴/32 − …", fontsize=8.6, color=PURPLE)
    ax.set_xlabel("tilt: highest latitude the orbit reaches (degrees)"); ax.set_ylabel("rate² × r³/M")
    ax.set_xlim(0, 92); ax.set_ylim(0.8, 1.26)
    ax.legend(fontsize=7.8, loc="upper right")
    ax.set_title("(b) Kepler's period law in r holds only in the equator", fontsize=10.5)
    # (c) Laplace's rule: radial + vertical, Newton and Kerr
    ax = axs[1, 0]
    ax.set_xlim(-3.4, 4.2); ax.set_ylim(-0.6, 7.8); ax.axis("off")
    cols = ["1", "u = M/r\n(mass)", "χ u³ᐟ²\n(spin)", "χ² u² = a²/r²\n(shape)"]
    for j, lab in enumerate(cols):
        ax.text(j + 0.5, 7.05, lab, ha="center", fontsize=8.6, weight="bold")
    table = [("Newton, shifted source", None),
             ("radial  (Ω_r/Ω)²", [1, 0, 0, -3]), ("up and down  (Ω_θ/Ω)²", [1, 0, 0, 3]), ("sum", [2, 0, 0, 0]),
             ("Kerr, exact", None),
             ("radial  (Ω_r/Ω_φ)²", [1, -6, 8, -3]), ("up and down  (Ω_θ/Ω_φ)²", [1, 0, -4, 3]),
             ("sum", [2, -6, 4, 0])]
    for i, (lab, vals) in enumerate(table):
        y = 6.6 - i*0.85
        if vals is None:
            ax.text(-3.3, y, lab, fontsize=9.2, weight="bold", va="center", color=BLUE)
            continue
        ax.text(-0.1, y, lab, ha="right", va="center", fontsize=8.8)
        for j, v in enumerate(vals):
            colr = GREEN if (lab == "sum" and j == 3) else (GREY if v == 0 else (PURPLE if j == 3 else ORANGE))
            ax.add_patch(FancyBboxPatch((j + 0.08, y - 0.33), 0.84, 0.66, boxstyle="round,pad=0.02",
                                        fc=colr, alpha=0.35 if (lab == "sum" and j == 3) else 0.15, ec=colr))
            ax.text(j + 0.5, y, f"{v:+d}" if v else "0", ha="center", va="center", fontsize=9.5,
                    weight="bold" if lab == "sum" else None)
    ax.text(-3.3, -0.45, "the shape term cancels in radial + vertical, in Newton and in Kerr alike (green):\n"
            "Laplace's rule for empty space; Kerr's two shape terms are the shifted source's, term for term",
            fontsize=8.6)
    ax.set_title("(c) what survives exactly in the equator: radial + vertical", fontsize=10.5, pad=12)
    # (d) Kepler's energy law survives for orbits with no spin about the axis
    ax = axs[1, 1]
    for key, ls in (("a=0.8, r=3.0", "-"), ("a=0.8, r=5.0", "--")):
        rr = tj["spheres"][key]["rows"]
        ax.plot([x["tilt_deg"] for x in rr], [x["energy"] for x in rr], ls, color=BLUE, lw=2,
                label=f"energy −2rE/M  ({key.split(', ')[1]})")
    ax.axhline(1, color=GREEN, lw=2.2, label="Kepler's energy law in r: E = −M/(2r)")
    ax.plot([90], [1], "s", color=GREEN, ms=7)
    ax.text(40, 0.94, "polar orbits (no spin about the axis):\nE = −M/(2r) exactly, radial motion is Kepler's,\n"
            "yet the period is not (panel b)", fontsize=8.8, color=BLUE)
    ax.set_xlabel("tilt: highest latitude the orbit reaches (degrees)"); ax.set_ylabel("−2 r E / M")
    ax.set_xlim(0, 92); ax.legend(fontsize=8.3, loc="lower right")
    ax.set_title("(d) what survives off the equator: Kepler's energy law, at the pole", fontsize=10.5)
    fig.tight_layout()
    save(fig, "fig26_tilt.png")


def fig27():
    import mpmath as mpm
    import polar as pl
    pj = json.load(open(os.path.join(HERE, "polar.json")))
    fig, axs = plt.subplots(2, 2, figsize=(14.5, 10.5))
    # (a) a polar orbit flown in Kerr, side view and top view
    ax = axs[0, 0]
    fl = pj["flights"][0]
    x, y, z = (np.array(v) for v in fl["path"])
    ax.plot(x, z, color=BLUE, lw=0.6, alpha=0.6)
    ax.plot(y + 15.5, x, color=PURPLE, lw=0.6, alpha=0.6)
    rp = 1 + np.sqrt(1 - 0.81)
    th = np.linspace(0, 2*np.pi, 200)
    for cx, cy in ((0, 0), (15.5, 0)):
        ax.fill(cx + np.sqrt(rp**2 + 0.81)*np.cos(th), cy + rp*np.sin(th) if cx == 0 else cy + np.sqrt(rp**2 + 0.81)*np.sin(th),
                color="k", alpha=0.85)
    ax.set_aspect("equal"); ax.axis("off"); ax.set_xlim(-8, 23); ax.set_ylim(-9.5, 9.5)
    ax.text(0, -7.6, "side view", ha="center", fontsize=9, color=BLUE)
    ax.text(15.5, -7.6, "from above the pole", ha="center", fontsize=9, color=PURPLE)
    ax.text(-7.8, 8.6, "a polar orbit around a black hole spinning at a = 0.9M, r = 6M, flown directly in Kerr's field:\n"
            "it passes over both poles; the spin drags its plane around (right)", fontsize=8.6)
    ax.set_title("(a) an orbit over the poles", fontsize=10.5, pad=12)
    # (b) the up-and-down rate against a/r: Kerr approaches the shifted source's law as the field weakens
    ax = axs[0, 1]
    kk = np.array(pj["curves"]["k"])
    ax.plot(kk, pj["curves"]["shifted source (Newton)"], "k--", lw=2.2, label="shifted source (Newton): (π/2)²/(2K − E)²")
    for uv, col in ((0.01, GREEN), (0.03, ORANGE), (0.1, RED)):
        ax.plot(kk, pj["curves"][f"Kerr, u = {uv}"], color=col, lw=2, label=f"Kerr, M/r = {uv}")
        kbh = uv                                    # a = M: the fastest a black hole can spin
        ax.plot([kbh], [float(pl.omega_theta_exact(uv, kbh, 1.0)**2/uv)], "o", color=col, ms=6)
    ax.text(0.02, 0.56, "dots: a = M, the fastest spin a black hole can have;\npast them the curves use a > M only to "
            "make the shape large", fontsize=8.4)
    ax.set_xlabel("k = a / r   (r: Boyer–Lindquist radius = the shift's distance)"); ax.set_ylabel("Ω_θ² r³ / M")
    ax.legend(fontsize=8.3, loc="lower left"); ax.set_ylim(0.3, 1.05)
    ax.set_title("(b) Kerr's polar rate approaches the shifted source's law", fontsize=10.5)
    # (c) the double series: Newtonian column = shifted source; other columns = mass and a² mixing
    ax = axs[1, 0]
    ax.set_xlim(-2.2, 4.3); ax.set_ylim(-2.3, 6.4); ax.axis("off")
    om = pj["series"]["omega_theta2"]
    heads = ["Newton\n(u⁰)", "× M/r\n(u¹)", "× (M/r)²\n(u²)", "× (M/r)³\n(u³)"]
    for j, h in enumerate(heads):
        ax.text(j + 0.5, 5.75, h, ha="center", fontsize=8.6, weight="bold")
    rowlab = ["1", "(a/r)²", "(a/r)⁴", "(a/r)⁶", "(a/r)⁸"]
    from fractions import Fraction as Fr
    for i, row in enumerate(om):
        yv = 4.9 - i*0.95
        ax.text(-0.15, yv, rowlab[i], ha="right", va="center", fontsize=9.5)
        for j, v in enumerate(row):
            fv = Fr(v)
            colr = GREEN if j == 0 else (GREY if fv == 0 else ORANGE)
            ax.add_patch(FancyBboxPatch((j + 0.06, yv - 0.38), 0.88, 0.76, boxstyle="round,pad=0.02", fc=colr,
                                        alpha=0.25 if fv != 0 else 0.12, ec=colr))
            ax.text(j + 0.5, yv, "0" if fv == 0 else (f"{float(fv):+.4g}" if fv.denominator > 64 else str(fv)),
                    ha="center", va="center", fontsize=8.8)
    ax.text(-2.1, -0.85, "Kerr's polar Ω_θ² r³/M, term by term. Green: exactly the shifted source's polar law.\n"
            "Orange: the a² terms mixing with the mass, which Newton's shifted source doesn't have.\n"
            "In the equator, the same radius gives Kerr's shape as −3a²/r² and nothing more.", fontsize=8.6)
    ax.set_title("(c) Newton's column is the shifted source; the rest is new", fontsize=10.5, pad=12)
    # (d) for a real black hole: how big each part is
    ax = axs[1, 1]
    chi = 0.9
    rs = np.logspace(np.log10(6), 3, 60)
    shape, left = [], []
    for rv in rs:
        F = float(pl.omega_theta_exact(1.0, chi, rv)*rv**1.5)**2
        Nw = float(pl.newton_polar(chi/rv))
        shape.append(abs(Nw - 1)); left.append(abs(F - Nw))
    ax.loglog(rs, shape, color=GREEN, lw=2.2, label="shape, as the shifted source gives it: |law − 1| ≈ (3/2) χ² (M/r)²")
    ax.loglog(rs, left, color=ORANGE, lw=2.2, label="what Kerr adds beyond it: ≈ 4 χ² (M/r)³")
    ax.loglog(rs, 1.5*chi**2/rs**2, ":", color=GREEN, lw=1); ax.loglog(rs, 4*chi**2/rs**3, ":", color=ORANGE, lw=1)
    ax.set_xlabel("orbit radius r / M   (black hole spinning at χ = 0.9)"); ax.set_ylabel("size of the term in Ω_θ² r³/M")
    ax.legend(fontsize=8.3, loc="lower left")
    ax.set_title("(d) for a black hole the extra terms fall one power of r faster", fontsize=10.5)
    fig.tight_layout()
    save(fig, "fig27_polar.png")


def fig28():
    import mpmath as mpm
    import sympy as sy
    import polar_split as ps
    sj = json.load(open(os.path.join(HERE, "polar_split.json")))
    pj = json.load(open(os.path.join(HERE, "polar.json")))
    fig, axs = plt.subplots(2, 2, figsize=(14.5, 10.5))
    # (a) the shape-meeting-mass term, piece by piece, over the poles and in the equator
    ax = axs[0, 0]
    pieces = ["E·E\n(time)", "E·P\n(push × time)", "P·P\n(push)", "S·S\n(twist)", "E·S\n(twist × time)", "P·S\n(twist × push)"]
    pole = [0, 16, -1, 0, 0, 0]
    equ = [0, 20, 8, -16, -12, 0]
    xs = np.arange(len(pieces))
    cols = [GREY, ORANGE, ORANGE, BLUE, BLUE, BLUE]
    ax.bar(xs - 0.2, pole, 0.38, color=cols, alpha=0.95, edgecolor="k", lw=0.6, label="over the poles")
    ax.bar(xs + 0.2, equ, 0.38, color=cols, alpha=0.4, edgecolor="k", lw=0.6, hatch="//", label="in the equator")
    for xv, v in zip(xs - 0.2, pole):
        ax.text(xv, v + (0.8 if v >= 0 else -2.2), f"{v:+d}" if v else "0", ha="center", fontsize=8.6)
    for xv, v in zip(xs + 0.2, equ):
        ax.text(xv, v + (0.8 if v >= 0 else -2.2), f"{v:+d}" if v else "0", ha="center", fontsize=8.6)
    for xv in (3, 4, 5):
        ax.text(xv - 0.2, 2.6, "gone:\nL = 0", ha="center", fontsize=7.6, color=BLUE)
    ax.axhline(0, color="k", lw=0.8)
    ax.set_xticks(xs); ax.set_xticklabels(pieces, fontsize=8.4)
    ax.set_ylabel("coefficient of (a/r)²(M/r) in (radial rate ÷ orbit rate)²")
    ax.set_ylim(-20, 26)
    ax.text(-0.45, 22.3, "totals: over the poles 16 − 1 = +15 (Kerr's term);  in the equator 20 + 8 − 16 − 12 = 0",
            fontsize=8.7, weight="bold")
    ax.legend(fontsize=8.4, loc="lower left")
    ax.set_title("(a) the extra polar term is the river's push, not the twist", fontsize=10.5)
    # (b) the plane turning: top view of the flown polar orbit, with its equator crossings
    ax = axs[0, 1]
    fl = pj["flights"][0]
    x, y, z = (np.array(v) for v in fl["path"])
    n = len(x)
    for i in range(n - 1):
        ax.plot(x[i:i + 2], y[i:i + 2], color=plt.cm.viridis(i/n), lw=0.9)
    rp = 1 + np.sqrt(1 - 0.81)
    th = np.linspace(0, 2*np.pi, 200)
    ax.fill(np.sqrt(rp**2 + 0.81)*np.cos(th), np.sqrt(rp**2 + 0.81)*np.sin(th), color="k", alpha=0.85)
    ii = np.where((z[:-1] < 0) & (z[1:] >= 0))[0]
    for j, i in enumerate(ii[:8]):
        ax.plot([0, 1.25*x[i]], [0, 1.25*y[i]], color=plt.cm.viridis(i/n), lw=1.4, ls="--")
        ax.plot(x[i], y[i], "o", color=plt.cm.viridis(i/n), ms=5)
    ax.annotate("", xy=(5.3, 5.3), xytext=(7.4, 1.6),
                arrowprops=dict(arrowstyle="->", color=RED, lw=1.8, connectionstyle="arc3,rad=0.25"))
    ax.text(7.6, 3.4, "spin and\nplane turn\nthe same way", fontsize=8.5, color=RED)
    ax.set_aspect("equal"); ax.axis("off"); ax.set_xlim(-8.5, 10.5); ax.set_ylim(-11.2, 8.8)
    ax.text(-8.4, -11.0, "seen from above the pole (a = 0.9M, r = 6M); dashed lines mark where the orbit\n"
            "crosses the equator going up, early (dark) to late (yellow). The crossing moves\n"
            "about 44° per lap: the twist turns the plane but leaves r, θ and the clock alone", fontsize=8.4)
    ax.set_title("(b) the twist only turns the polar orbit's plane", fontsize=10.5)
    # (c) the turning rate: exact, flown, and Lense-Thirring far out
    ax = axs[1, 0]
    rs = np.logspace(np.log10(6), np.log10(300), 50)
    ex = [float(ps.node_rate_exact(1.0, 0.9, rv)[0]) for rv in rs]
    ax.loglog(rs, ex, color=BLUE, lw=2.2, label="exact (from dφ/dλ = 2MraE/Δ in Mino time), all from E·S and P·S")
    ax.loglog(rs, 2*0.9/rs**3, "k--", lw=1.4, label="Lense–Thirring (1918): 2J/r³")
    nf = sj["node_flights"]
    ax.plot([6.0], [nf[0][0]], "o", color=RED, ms=8, label="flown in the exact field (a = 0.9, r = 6)")
    ax.set_xlabel("orbit radius r / M   (spin a = 0.9 M)"); ax.set_ylabel("turning rate of the plane (radians per unit of time t, M = 1)")
    ax.legend(fontsize=8.4, loc="lower left")
    ins = ax.inset_axes([0.6, 0.6, 0.36, 0.33])
    ins.semilogx(rs, np.array(ex)/(2*0.9/rs**3), color=BLUE, lw=1.8)
    ins.axhline(1, color="k", ls="--", lw=1); ins.set_ylim(0.95, 1.005)
    ins.set_title("exact ÷ Lense–Thirring", fontsize=8); ins.tick_params(labelsize=7)
    ax.set_title("(c) how fast the plane turns: Lense–Thirring far out, exact everywhere", fontsize=10.5)
    # (d) the push also moves the elliptic modulus; the push alone doesn't
    ax = axs[1, 1]
    rr, MM, aa = ps.r, ps.M, ps.a
    mEP = ps.polar_constants({"EE": 1, "EP": 1, "PP": 0})["m"]
    mFull = ps.polar_constants({"EE": 1, "EP": 1, "PP": 1})["m"]
    fEP = sy.lambdify((rr, MM, aa), mEP/(aa**2/rr**2)); fF = sy.lambdify((rr, MM, aa), mFull/(aa**2/rr**2))
    rs2 = np.linspace(5, 40, 200)
    ax.axhline(1, color=GREY, lw=2.2, label="E·E alone, or with P·P: m = a²/r², the shifted source's")
    ax.plot(rs2, fEP(rs2, 1.0, 0.9), color=ORANGE, lw=2.2, label="with E·P: m moves")
    ax.plot(rs2, fF(rs2, 1.0, 0.9), color=RED, lw=1.6, ls="--", label="the full field (Kerr)")
    ax.set_xlabel("orbit radius r / M   (spin a = 0.9 M)"); ax.set_ylabel("modulus m ÷ (a/r)²")
    ax.legend(fontsize=8.4, loc="lower right")
    ax.text(10, 0.5, "m sets the up-and-down period through K(m).\nThe push meeting the river's time part\n"
            "lowers it by a fraction of about 4M/r", fontsize=8.5)
    ax.set_title("(d) the push meeting time also reshapes the period itself", fontsize=10.5)
    fig.tight_layout()
    save(fig, "fig28_polar_split.png")


def fig29():
    import polar_split as ps
    jj = json.load(open(os.path.join(HERE, "jet_boost.json")))
    fig, axs = plt.subplots(2, 2, figsize=(14.5, 10.5))
    # (a) the one jet, squared: six pieces; the cross piece P.S exists only if push and twist are one jet
    ax = axs[0, 0]
    ax.set_xlim(-1.6, 4.6); ax.set_ylim(-1.9, 4.3); ax.set_aspect("equal"); ax.axis("off")
    labs = ["E\n(time)", "P\n(push)", "S\n(twist)"]
    cols = {"EE": GREY, "EP": ORANGE, "PP": ORANGE, "ES": BLUE, "SS": BLUE, "PS": RED}
    for i in range(3):
        ax.text(i + 0.5, 3.45, labs[i], ha="center", va="center", fontsize=9.5, weight="bold")
        ax.text(-0.35, 2.5 - i, labs[i], ha="center", va="center", fontsize=9.5, weight="bold")
    nm = ["E", "P", "S"]
    for i in range(3):
        for j in range(3):
            key = "".join(sorted([nm[i], nm[j]], key=lambda c: "EPS".index(c)))
            c = cols[key]
            ax.add_patch(FancyBboxPatch((j + 0.05, 2.05 - i), 0.9, 0.9, boxstyle="round,pad=0.02", fc=c,
                                        alpha=0.55 if key == "PS" else 0.25, ec=c, lw=2.2 if key == "PS" else 1))
            ax.text(j + 0.5, 2.5 - i, f"{nm[i]}·{nm[j]}", ha="center", va="center", fontsize=10,
                    weight="bold" if key == "PS" else "normal")
    ax.text(3.2, 2.95, "one jet:\n(E + P + S)²\nall nine tiles", fontsize=9.2, va="top")
    ax.text(3.2, 1.6, "two velocities:\n(E + P)² + (E + S)² − E²\nthe red tiles missing", fontsize=9.2, va="top",
            color=RED)
    ax.text(-1.5, -0.35, "The push and the twist are two shears into the same leg of each cell, so as jets they\n"
            "compose into one shear, a single ray k = time + push + twist. The field is that ray squared.\n"
            "Its cross piece P·S, push meeting twist, is the only thing the jet rule adds; without it\n"
            "the field is no longer empty space (Ricci ~ 10⁻² near the hole instead of 10⁻⁵⁰).", fontsize=8.6,
            va="top")
    ax.set_title("(a) one jet, squared: the cross piece P·S is the jet rule's", fontsize=10.5, pad=10)
    # (b) equator: the spin term and the clock gap, one jet vs two velocities
    ax = axs[0, 1]
    rr = np.linspace(2.0, 40, 400); u = 1/rr
    ax.axhline(1, color=GREY, lw=2.6, label=r"one jet (Kerr): spin term $8\chi u^{3/2}$ and clock gap $2\pi a$, at every radius")
    ax.plot(rr, (1 - 2*u)*(1 - 3*u), color=RED, lw=2.2, label="two velocities, spin term: × (1 − 2u)(1 − 3u)")
    ax.plot(rr, 1 - 4*u, color=ORANGE, lw=2.2, ls="--", label="two velocities, clock gap: × (1 − 4M/r)")
    gp = jj["clock_gap"]["two velocities"]["per unit a, small spin"]
    ax.plot([g[0] for g in gp if 4 < g[0] <= 40], [g[1] for g in gp if 4 < g[0] <= 40], "o", color=ORANGE, ms=6)
    ax.axhline(0, color="k", lw=0.6)
    for x0, lab, yv in ((3, "r = 3M: spin term 0", -0.3), (4, "r = 4M: clock gap 0", -0.42)):
        ax.axvline(x0, color="k", lw=0.6, ls=":"); ax.text(x0 + 0.4, yv, lab, fontsize=8)
    ax.set_xlim(2, 40)
    ax.set_xlabel("orbit radius r / M   (equator, small spin)"); ax.set_ylabel("÷ Kerr's value")
    ax.set_ylim(-0.5, 1.25); ax.legend(fontsize=8.3, loc="lower right")
    ax.set_title("(b) in the equator, two velocities spoil the spin term and the clock gap", fontsize=10.5)
    # (c) the spin tower piece by piece
    ax = axs[1, 0]
    tw = jj["spin_tower_chi_u15_u25_u35"]
    groups = [("push + E·S\n(twist meets time)", tw["EE+PP+EP+ES"]), ("push + P·S\n(twist meets push)", tw["EE+PP+EP+PS"]),
              ("all six pieces\n(one jet)", tw["EE+PP+EP+SS+ES+PS"]), ("all but P·S\n(two velocities)", tw["EE+PP+EP+SS+ES"])]
    xs = np.arange(len(groups)); wbar = 0.26
    for k, (lab, col) in enumerate(((r"$\chi u^{3/2}$", GREY), (r"$\chi u^{5/2}$", RED), (r"$\chi u^{7/2}$", PURPLE))):
        vals = [g[1][k] for g in groups]
        ax.bar(xs + (k - 1)*wbar, vals, wbar, color=col, alpha=0.8, label=lab)
        for xv, v in zip(xs + (k - 1)*wbar, vals):
            ax.text(xv, v + (2 if v >= 0 else -6), f"{v:+.0f}" if abs(v) > 1e-6 else "0", ha="center", fontsize=8.2)
    ax.axhline(0, color="k", lw=0.8)
    ax.set_xticks(xs); ax.set_xticklabels([g[0] for g in groups], fontsize=8.6)
    ax.set_ylabel("coefficient in (radial rate ÷ orbit rate)²"); ax.set_ylim(-60, 78)
    ax.legend(fontsize=8.4, loc="upper center", ncol=3, title="spin term", title_fontsize=8.4)
    ax.set_title("(c) the twist through E·S leaves a tower; through P·S it takes it out", fontsize=10.5)
    # (d) over the poles: the plane's turning, one jet vs two velocities
    ax = axs[1, 1]
    rs = np.logspace(np.log10(3), np.log10(300), 80); av = 0.9
    one = np.array([float(ps.node_rate_exact(1.0, av, rv)[0]) for rv in rs])
    ratio = (rs**2 - 2*rs + av**2)/(rs**2 + av**2)
    ax.loglog(rs, one, color=GREY, lw=2.6, label="one jet (Kerr): from 2MraE/Δ")
    ax.loglog(rs, one*ratio, color=RED, lw=2.2, label="two velocities: from 2MraE/(r² + a²)")
    ax.loglog(rs, 2*av/rs**3, "k--", lw=1.2, label="Lense–Thirring 2J/r³")
    ax.set_xlabel("orbit radius r / M   (polar orbit, spin a = 0.9 M)")
    ax.set_ylabel("turning rate of the plane (radians per unit of time t, M = 1)")
    ax.legend(fontsize=8.4, loc="lower left")
    ins = ax.inset_axes([0.58, 0.58, 0.38, 0.35])
    ins.semilogx(rs, ratio, color=RED, lw=1.8); ins.axhline(1, color=GREY, lw=1.6)
    ins.set_title("two velocities ÷ one jet = Δ/(r² + a²)", fontsize=8); ins.tick_params(labelsize=7)
    ins.plot([6], [(36 - 12 + 0.81)/36.81], "o", color=RED, ms=4); ins.text(7, 0.6, "0.67 at r = 6M", fontsize=7)
    ax.set_title("(d) over the poles, the cross piece is the strong-field part of the plane's turning", fontsize=10.5)
    fig.tight_layout()
    save(fig, "fig29_jet_boost.png")


def sp_frac(s):
    from fractions import Fraction as Fr
    return Fr(s)


# ---------------------------------------------------------------- fig30 round 22: one law for both wobbles
def fig30():
    wl = json.load(open(os.path.join(HERE, "wobble_law.json")))
    fig, axs = plt.subplots(2, 2, figsize=(14.5, 10.5))
    # (a) the two wobbles as functions of the spin, at several strengths: every curve has the same bend
    ax = axs[0, 0]
    qs = np.linspace(-0.45, 0.45, 300)
    cols = [GREY, BLUE, PURPLE, RED]
    for c, uu in zip(cols, (0.0, 0.04, 0.08, 0.12)):
        vv = np.sqrt(uu); Dd = 1 - 3*uu + 2*qs*vv
        ax.plot(qs, Dd - 3*(vv - qs)**2, color=c, lw=2.2, label=("M/r → 0 (Newton)" if uu == 0 else f"M/r = {uu:g}"))
        ax.plot(qs, Dd + 3*(vv - qs)**2, color=c, lw=1.6, ls="--")
    ax.text(0.05, 1.55, "up-and-down:  D + 3(v − q)²", fontsize=9.5, color="k")
    ax.text(0.05, -0.7, "in-and-out (solid):  D − 3(v − q)²", fontsize=9.5, color="k")
    ax.text(0.05, 1.42, "(dashed)", fontsize=9, color="k")
    ax.axhline(0, color="k", lw=0.6)
    ax.set_xlabel("spin seen from the orbit,  q = a/r   (negative: orbit runs against the spin)")
    ax.set_ylabel("(wobble rate ÷ orbit rate)²")
    ax.legend(fontsize=8.5, loc="lower center", ncol=4)
    ax.set_title("(a) one law for both wobbles: the bend in q is ∓3 at every strength\n"
                 "(the curves for stronger fields are the same parabolas, slid over)", fontsize=10.5)
    # (b) Mino-time stiffness: a seesaw that always adds to 2Mr
    ax = axs[0, 1]
    rr = np.linspace(1.6, 20, 400)
    for aa, ls, nm in ((0.9, "-", "with the spin (a = 0.9M)"), (-0.9, ":", "against the spin (a = −0.9M)")):
        vv = np.sqrt(1/rr); qq = aa/rr
        x2r2 = (vv - qq)**2/(1 - 3*vv**2 + 2*qq*vv)
        ok = (1 - 3*vv**2 + 2*qq*vv) > 0.02
        ax.plot(rr[ok], 1 - 3*x2r2[ok], color=ORANGE, ls=ls, lw=2.2, label=f"in-and-out, {nm}")
        ax.plot(rr[ok], 1 + 3*x2r2[ok], color=BLUE, ls=ls, lw=2.2, label=f"up-and-down, {nm}")
    ax.axhline(1, color=GREY, lw=1.2); ax.axhline(0, color="k", lw=0.6)
    ax.text(11.5, 0.9, "average: exactly 1 (the two add to 2Mr)", fontsize=9, color=GREY)
    ax.annotate("innermost stable orbit\n2.32M", (2.32, 0), (3.2, -0.45), fontsize=8, arrowprops=dict(arrowstyle="->"))
    ax.annotate("8.72M", (8.72, 0), (9.6, -0.4), fontsize=8, arrowprops=dict(arrowstyle="->"))
    ax.set_ylim(-0.6, 2.6); ax.set_xlim(1.6, 20)
    ax.set_xlabel("orbit radius r / M")
    ax.set_ylabel("stiffness in Mino time ÷ Mr  =  1 ∓ 3 x²/r²")
    ax.legend(fontsize=8, loc="upper right")
    ax.set_title("(b) the mass meets the orbit's spin-relative turning x = L − aE only as 2Mr·x²,\n"
                 "one power of r, so it moves stiffness from in-and-out to up-and-down and adds none", fontsize=10.5)
    # (c) the river's jet energy on the orbit: energy + push + twist = 1/sqrt(D), spin to the first power only
    ax = axs[1, 0]
    radii = np.array([3.0, 4.0, 6.0, 10.0, 20.0]); aa = 0.9
    vv = np.sqrt(1/radii); qq = aa/radii; Dd = 1 - 3*vv**2 + 2*qq*vv
    Ek = (1 - 2*vv**2 + qq*vv)/np.sqrt(Dd); Lk = radii*vv*(1 - 2*qq*vv + qq**2)/np.sqrt(Dd)
    Sk = -aa*Lk/(radii**2 + aa**2); Pk = 2*radii*(Ek + Sk)/(radii**2 - 2*radii + aa**2)
    xs = np.arange(len(radii)); wdt = 0.6
    ax.bar(xs, Ek + Sk, wdt, color=GREY, label="E (energy), less the twist")
    ax.bar(xs + wdt/2 + 0.07, -Sk, 0.12, bottom=Ek + Sk, color=BLUE,
           label="S (the rays' twist): takes this much off E")
    ax.bar(xs, Pk, wdt, bottom=Ek + Sk, color=ORANGE, alpha=0.9, label="P (the river's push, on the orbit)")
    ax.plot(xs, 1/np.sqrt(Dd), "k_", ms=34, mew=2.5, label="1/√D,  D = 1 − 3v² + 2qv")
    for i in range(len(radii)):
        ax.text(xs[i], Ek[i] + Sk[i] + Pk[i] + 0.04, f"{Ek[i] + Pk[i] + Sk[i]:.3f}", ha="center", fontsize=8.5)
    ax.set_xticks(xs); ax.set_xticklabels([f"r = {int(r_)}M" for r_ in radii])
    ax.set_ylabel("jet energy k·p = E + P + S on the circular orbit")
    ax.set_ylim(0, 2.25); ax.legend(fontsize=8.3, loc="upper right")
    ax.set_title("(c) on the orbit, energy + push + twist lands exactly on 1/√D (a = 0.9M, prograde)\n"
                 "x = L − aE comes out as r(v − q)·k·p, the turning rate as r v·k·p", fontsize=10.5)
    # (d) two ways to cut the a^2 M term; both sum to zero, the law says why
    ax = axs[1, 1]
    cut1 = [("2E·P", 20, ORANGE), ("P·P", 8, ORANGE), ("S·S", -16, BLUE), ("2E·S", -12, BLUE), ("2P·S", 0, RED)]
    jp = wl["jet_pieces"]
    e0 = float(sp_frac(jp["E alone"]["a2M"])); e1 = float(sp_frac(jp["E + push of E"]["a2M"]))
    e2 = float(sp_frac(jp["E + twist"]["a2M"]))
    cut2 = [("E alone", e0, GREY), ("push of E", e1 - e0, ORANGE), ("twist", e2 - e0, BLUE), ("push of twist", 0, RED)]
    xs1 = np.arange(len(cut1))
    ax.bar(xs1, [c[1] for c in cut1], 0.7, color=[c[2] for c in cut1])
    for xx, c in zip(xs1, cut1):
        ax.text(xx, c[1] + (1.0 if c[1] >= 0 else -2.4), f"{c[1]:+g}", ha="center", fontsize=9)
    ax.set_xticks(xs1); ax.set_xticklabels([c[0] for c in cut1], fontsize=8.5)
    ax.set_xlim(-0.6, 9.4); ax.axhline(0, color="k", lw=0.8); ax.set_ylim(-22, 27)
    ax.text(2, 23, "round 20: cut the field\n(a/r)²(M/r) in T, sum 0", ha="center", fontsize=9)
    ins = ax.inset_axes([0.56, 0.1, 0.42, 0.62])
    xs2 = np.arange(len(cut2))
    ins.bar(xs2, [c[1] for c in cut2], 0.7, color=[c[2] for c in cut2])
    for xx, c in zip(xs2, cut2):
        ins.text(xx, c[1] + (0.2 if c[1] >= 0 else -0.55), f"{c[1]:+g}", ha="center", fontsize=8.5)
    ins.set_xticks(xs2); ins.set_xticklabels([c[0] for c in cut2], fontsize=7.5); ins.tick_params(labelsize=7)
    ins.axhline(0, color="k", lw=0.8); ins.set_ylim(-4.5, 5.5)
    ins.set_title("cut the jet energy on the orbit:\n(a/r)²(M/r) in 1/(k·p)², sum 0", fontsize=8.5)
    ax.set_ylabel("share of the (a/r)²(M/r) term")
    ax.set_title("(d) different cuts give different shares and the same zero: the shares are bookkeeping;\n"
                 "the reason is that spin squared can only enter as (v − q)²", fontsize=10.5)
    fig.tight_layout()
    save(fig, "fig30_wobble_law.png")

# ---------------------------------------------------------------- fig31 Kerr = Newton's ring + one mass term (round 23)
def fig31():
    pc = json.load(open(os.path.join(HERE, "polar_coupling.json")))
    cv = pc["curves"]
    fig, axs = plt.subplots(2, 2, figsize=(14.5, 10.5))
    # (a) the picture: Newton's field of the ring-shifted mass, and a polar orbit around it
    ax = axs[0, 0]
    aa, r0 = 0.6, 1.0
    th = np.linspace(0, 2*np.pi, 400)
    for rr, al in ((0.75, 0.25), (1.0, 1.0), (1.4, 0.25)):
        ax.plot(np.sqrt(rr**2 + aa**2)*np.sin(th), rr*np.cos(th), color=(BLUE if al == 1 else GREY),
                lw=(2.6 if al == 1 else 0.9), ls=("-" if al == 1 else ":"))
    ax.plot([-aa, aa], [0, 0], color="k", lw=4, solid_capstyle="round")
    ax.plot([-aa, aa], [0, 0], "o", color="k", ms=7)
    ax.annotate("the ring: Newton's mass shifted to z = ±ia\n(seen edge-on; Euler's two-centre problem)",
                xy=(aa*0.6, 0), xytext=(-0.3, -0.5), fontsize=9, arrowprops=dict(arrowstyle="->", color="k"))
    ax.annotate("a polar orbit at fixed r (over both poles)", xy=(0.55, 0.85), xytext=(0.65, 1.15), fontsize=9,
                color=BLUE, arrowprops=dict(arrowstyle="->", color=BLUE))
    eqs = ("Kerr, in Mino time, term by term:\n"
           "  in-and-out:   Newton + 2Mr (K − a²)\n"
           "  up-and-down:  Newton, nothing added\n"
           "  clock:        E·Newton + 2Mr P/Δ\n"
           "  turning:      Newton + (a/(r²+a²))·2Mr P/Δ")
    ax.text(-1.65, -1.62, eqs, fontsize=9.2, family="monospace", va="top",
            bbox=dict(boxstyle="round", fc="white", ec=GREY))
    ax.set_xlim(-1.75, 1.85); ax.set_ylim(-2.75, 1.5); ax.set_aspect("equal"); ax.axis("off")
    ax.set_title("(a) Kerr's orbits are Newton's orbits around the spinning ring,\n"
                 "plus one term that carries the mass (K = Carter's constant)", fontsize=10.5)
    # (b) one Kerr polar orbit's radial function: Newton's part alone has no double root; the one term supplies it
    ax = axs[0, 1]
    Rs = cv["R_split"]; xs = np.array(Rs["r"])
    sc = xs**4
    ax.plot(xs, np.array(Rs["newton"])/sc, color=GREY, lw=2, ls="--", label="Newton's part (ring field), same constants")
    ax.plot(xs, np.array(Rs["extra"])/sc, color=ORANGE, lw=2, label="the one term  2Mr (K − a²)")
    ax.plot(xs, np.array(Rs["kerr"])/sc, color=BLUE, lw=2.6, label="Kerr = the sum")
    ax.axhline(0, color="k", lw=0.7); ax.axvline(Rs["r0"], color=BLUE, lw=0.8, ls=":")
    ax.text(Rs["r0"] + 0.15, 0.03, "orbit at r = 6M:\nthe sum just touches zero", fontsize=9, color=BLUE)
    ax.set_xlabel("r / M"); ax.set_ylabel("R(r) / r⁴   (radial speed² in Mino time, scaled)")
    ax.set_ylim(-0.13, 0.1)
    ax.legend(fontsize=8.5, loc="lower right")
    ax.set_title("(b) one Kerr orbit over the poles (a = 0.9M, r = 6M): Newton's part alone stays below zero;\n"
                 "the mass term lifts the sum until it just touches zero at the orbit", fontsize=10.5)
    # (c) the (a/r)^2 part of (Omega_r/Omega_theta)^2, exact in M/r, with pieces of the one term switched off
    ax = axs[1, 0]
    us = np.array(cv["u"])
    sty = [("Kerr", BLUE, "-", "Kerr:  (3/2)(1 + 10u − 8u²)"),
           ("2MrQ off", RED, "--", "only the binding part of x²:  (3/2)(1 + 4u)"),
           ("binding part of x^2 off", PURPLE, "-.", "only 2MrQ:  (3/2)(1 + 2u)(1 + 4u)"),
           ("whole term 2Mr(K-a^2) off", GREY, ":", "whole term off:  3/2  (Newton's ring)")]
    for key, c, ls, lab in sty:
        ax.plot(us, cv["ratio_k2: " + key], color=c, ls=ls, lw=2.3, label=lab)
    ax.set_xlabel("field strength  u = M/r"); ax.set_ylabel("coefficient of (a/r)² in (Ω_r/Ω_θ)²")
    ax.legend(fontsize=8.5, loc="upper left")
    ax.set_title("(c) over the poles, x² = a²E² is only part of the source: its binding part\n"
                 "gives 6 of round 19's 15; Schwarzschild's half, 2MrQ, gives the other 9", fontsize=10.5)
    # (d) bookkeeping of the mixed terms
    ax = axs[1, 1]
    groups = [("(a/r)²(M/r)\nin (Ω_r/Ω_θ)²", [("2MrQ", 9, PURPLE), ("binding of x²", 6, RED), ("total", 15, BLUE)]),
              ("(a/r)²(M/r)²\nin (Ω_r/Ω_θ)²", [("2MrQ", 12, PURPLE), ("binding of x²", 0, RED), ("together", -24, ORANGE),
                                             ("total", -12, BLUE)]),
              ("(a/r)²(M/r)\nin Ω_θ² r³/M", [("radial term", -10.5, PURPLE), ("clock terms", 6.5, GREEN),
                                            ("total", -4, BLUE)])]
    x0 = 0; xt = []; xl = []
    for gname, bars in groups:
        for i, (nm, val, c) in enumerate(bars):
            ax.bar(x0 + i, val, 0.75, color=c, alpha=(1.0 if nm in ("total",) else 0.7))
            ax.text(x0 + i, val + (0.6 if val >= 0 else -2.0), f"{val:+g}", ha="center", fontsize=8.8)
            xt.append(x0 + i); xl.append(nm)
        ax.text(x0 + (len(bars) - 1)/2, 19, gname, ha="center", fontsize=8.8)
        x0 += len(bars) + 1.2
    ax.axhline(0, color="k", lw=0.8)
    ax.set_xticks(xt); ax.set_xticklabels(xl, fontsize=8, rotation=35, ha="right")
    ax.set_ylim(-28, 24.5); ax.set_ylabel("coefficient")
    ax.set_title("(d) where round 19's mixed terms come from: all from the one term (and, for\n"
                 "the coordinate-time rate, the clock); switch both off and Newton's law returns", fontsize=10.5)
    fig.tight_layout()
    save(fig, "fig31_newton_plus_one.png")


for f in (fig1, fig2, fig3, fig4, fig5, fig6, fig7, fig8, fig9, fig10, fig11, fig12, fig13, fig14, fig15, fig16, fig17, fig18, fig19, fig20, fig21, fig22, fig23, fig24, fig25, fig26, fig27, fig28, fig29, fig30, fig31):
    f()
print("figures:", sorted(os.listdir(FIG)))
