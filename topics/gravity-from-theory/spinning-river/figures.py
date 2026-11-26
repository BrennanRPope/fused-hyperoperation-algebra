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

for f in (fig1, fig2, fig3, fig4, fig5, fig6, fig7, fig8, fig9, fig10):
    f()
print("figures:", sorted(os.listdir(FIG)))
