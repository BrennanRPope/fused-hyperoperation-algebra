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


for f in (fig1, fig2, fig3, fig4, fig5, fig6, fig7, fig8, fig9, fig10, fig11, fig12, fig13, fig14, fig15, fig16, fig17, fig18, fig19, fig20, fig21, fig22, fig23):
    f()
print("figures:", sorted(os.listdir(FIG)))
