"""DRAFT figures for jet-river/README.md. Run: python3 figures.py  (writes figures/*.png)"""
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
from matplotlib.patches import Ellipse, FancyBboxPatch, Circle  # noqa: E402
from mpl_toolkits.mplot3d.art3d import Poly3DCollection  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
FIG = os.path.join(HERE, "figures")
os.makedirs(FIG, exist_ok=True)
RED, GREEN, ORANGE, BLUE, PURPLE, GREY = "#d62728", "#2ca02c", "#ff7f0e", "#1f77b4", "#9467bd", "#7f7f7f"


def save(fig, name):
    fig.savefig(os.path.join(FIG, name), dpi=130, bbox_inches="tight")
    plt.close(fig)


# ---------------------------------------------------------------- fig1 the cell's one-step box
def box(ax, x0, y0, z0, dx, dy, dz, color, alpha=0.55):
    X = [x0, x0 + dx]; Y = [y0, y0 + dy]; Z = [z0, z0 + dz]
    v = np.array([[X[i], Y[j], Z[k]] for i in (0, 1) for j in (0, 1) for k in (0, 1)])
    faces = [[0, 1, 3, 2], [4, 5, 7, 6], [0, 1, 5, 4], [2, 3, 7, 6], [0, 2, 6, 4], [1, 3, 7, 5]]
    ax.add_collection3d(Poly3DCollection([v[f] for f in faces], facecolor=color, edgecolor="k", lw=0.4, alpha=alpha))


def fig1():
    fig = plt.figure(figsize=(12, 5.2))
    ax = fig.add_subplot(1, 2, 1, projection="3d")
    x, y, z = 0.55, 0.35, 0.35                       # h a, h b, h b
    box(ax, 0, 0, 0, 1, 1, 1, "#dddddd", 0.35)
    for (o, d, c) in [((1, 0, 0), (x, 1, 1), BLUE), ((0, 1, 0), (1, y, 1), BLUE), ((0, 0, 1), (1, 1, z), BLUE)]:
        box(ax, *o, *d, c)
    for (o, d) in [((1, 1, 0), (x, y, 1)), ((1, 0, 1), (x, 1, z)), ((0, 1, 1), (1, y, z))]:
        box(ax, *o, *d, ORANGE)
    box(ax, 1, 1, 1, x, y, z, RED, 0.85)
    ax.set_xlim(0, 1.6); ax.set_ylim(0, 1.6); ax.set_zlim(0, 1.6)
    ax.set_box_aspect((1, 1, 1)); ax.view_init(22, 35); ax.set_axis_off()
    ax.text2D(0.5, -0.02, "edges: radial 1 + h·v',  two sideways 1 + h·v/r\n(grey: the old cell; colored: what one step adds)",
              transform=ax.transAxes, ha="center", fontsize=9.5)
    ax.set_title("one time step of one river cell = a circle-operation cube\n(1 + ha)(1 + hb)(1 + hb)", fontsize=11)
    ax2 = fig.add_subplot(1, 2, 2); ax2.axis("off")
    rows = [("#dddddd", "1", "the cell itself"),
            (BLUE, "h·P  (3 slabs)", "P = trace = how fast the cell's volume changes"),
            (ORANGE, "h²·Q  (3 beams)", "Q = pairs of rates = the field equations of gravity"),
            (RED, "h³·R  (1 corner)", "R = all three rates at once")]
    for k, (c, a, b) in enumerate(rows):
        yy = 0.85 - 0.2 * k
        ax2.add_patch(FancyBboxPatch((0.0, yy - 0.06), 0.06, 0.1, boxstyle="round,pad=0.01", fc=c, ec="k"))
        ax2.text(0.09, yy, a, fontsize=11, weight="bold", va="center")
        ax2.text(0.09, yy - 0.07, b, fontsize=9.5, va="center", color="#333333")
    ax2.text(0.0, 0.05, "river cell: a = v' (radial rate), b = v/r (sideways rate)\n"
             "Q = 2ab + b² = (r v²)' / r²     (Lean: river_e2, river_e2_gauss)", fontsize=10,
             bbox=dict(fc="#fff7e6", ec=ORANGE))
    ax2.set_xlim(0, 1); ax2.set_ylim(0, 1)
    save(fig, "fig1_cell_box.png")


# ---------------------------------------------------------------- fig2 Q is the field equation
def fig2():
    fig, axs = plt.subplots(2, 1, figsize=(8.5, 7), sharex=True)
    r = np.linspace(0.05, 4, 800)
    Rb, M = 1.0, 1.0
    H = np.sqrt(2 * M / Rb**3)
    vS = np.sqrt(2 * M / r)
    vball = np.where(r < Rb, H * r, vS)
    vconst = np.full_like(r, 1.0)
    Qs = [np.zeros_like(r), np.where(r < Rb, 3 * H**2, 0.0), 1 / r**2]
    for vv, Q, c, lab in [(vS, Qs[0], BLUE, "point mass: v² = 2M/r (Schwarzschild)"),
                          (vball, Qs[1], GREEN, "uniform ball, radius 1, then empty space"),
                          (vconst, Qs[2], ORANGE, "constant speed v = 1 (not a vacuum)")]:
        axs[0].plot(r, vv, color=c, lw=2.2, label=lab)
        axs[1].plot(r, Q, color=c, lw=2.2)
    axs[0].set_ylim(0, 3); axs[0].set_ylabel("inflow speed v(r)")
    axs[0].legend(fontsize=9); axs[0].set_title("the river's speed law (top) and its Q layer (bottom)")
    axs[1].axhline(0, color="k", lw=0.6)
    axs[1].set_ylim(-0.5, 7.5); axs[1].set_ylabel("Q = (r v²)'/r² = 8πGρ")
    axs[1].set_xlabel("radius r")
    axs[1].annotate("inside the ball: Q = 3H²\n(Friedmann's equation)", (0.5, 3 * H**2), (1.4, 5.6), fontsize=9,
                    color=GREEN, arrowprops=dict(arrowstyle="->", color=GREEN))
    axs[1].annotate("Q = 0: empty space\nexactly Schwarzschild", (2.5, 0), (2.3, 1.6), fontsize=9, color=BLUE,
                    arrowprops=dict(arrowstyle="->", color=BLUE))
    axs[1].annotate("Q = 1/r² ≠ 0: needs matter", (3.0, 1 / 9), (2.8, 3.2), fontsize=9, color=ORANGE,
                    arrowprops=dict(arrowstyle="->", color=ORANGE))
    save(fig, "fig2_q_is_the_field_equation.png")


# ---------------------------------------------------------------- fig3 null Lagrangian: Gauss law
def fig3():
    fig, axs = plt.subplots(1, 2, figsize=(12, 5.3))
    ax = axs[0]
    th = np.linspace(0, 2 * np.pi, 200)
    for rr, a in [(1.0, 0.25)]:
        ax.add_patch(Circle((0, 0), rr, fc=GREEN, alpha=a, ec=GREEN, lw=2))
    ax.add_patch(Circle((0, 0), 1.7, fc="none", ec="k", ls="--", lw=1.2))
    for t in np.linspace(0, 2 * np.pi, 16, endpoint=False):
        ax.annotate("", (1.7 * np.cos(t) * 0.78, 1.7 * np.sin(t) * 0.78), (1.7 * np.cos(t) * 1.12, 1.7 * np.sin(t) * 1.12),
                    arrowprops=dict(arrowstyle="->", color=BLUE, lw=1.4))
    ax.text(0, 0, "Q > 0\nmatter", ha="center", va="center", fontsize=10, color="#145a14")
    ax.text(0, -2.25, "sum of Q over the ball  =  flow through its skin\n∫ Q dV = 4π r²·(v²/r) = 4π r v²",
            ha="center", fontsize=10)
    ax.text(1.25, 1.55, "skin at radius r\n(arrows: the inflow)", fontsize=9)
    ax.set_xlim(-2.6, 2.6); ax.set_ylim(-2.7, 2.3); ax.set_aspect("equal"); ax.axis("off")
    ax.set_title("Q is a pure divergence (a null Lagrangian)")
    ax = axs[1]
    r = np.linspace(0.01, 3, 600)
    rho = np.where(r < 1, 1.0, 0.0) + np.where((r > 1.5) & (r < 2), 0.4, 0.0)       # a ball and a shell
    m = np.cumsum(4 * np.pi * r**2 * rho) * (r[1] - r[0])
    v = np.sqrt(2 * m / r)
    dv = np.gradient(v, r)
    Q = (2 * r * dv * v + v**2) / r**2
    lhs = np.cumsum(4 * np.pi * r**2 * Q) * (r[1] - r[0])
    ax.plot(r, lhs, color=GREEN, lw=4, alpha=0.5, label="∫ Q dV over the ball of radius r")
    ax.plot(r, 4 * np.pi * r * v**2, "k--", lw=1.5, label="4π r v²  (flow through the skin)")
    ax.plot(r, 8 * np.pi * m, color=PURPLE, lw=1, ls=":", label="8π G m(r)  (enclosed mass)")
    ax.axvspan(0, 1, color=GREEN, alpha=0.08); ax.axvspan(1.5, 2, color=GREEN, alpha=0.08)
    ax.set_xlabel("radius r"); ax.legend(fontsize=9)
    ax.set_title("a ball and a shell: all three curves agree\n(empty gaps: flat, Q = 0)")
    save(fig, "fig3_gauss_from_the_box.png")


# ---------------------------------------------------------------- fig4 tides are the second jet
def fig4():
    fig, axs = plt.subplots(1, 2, figsize=(13, 5.2), gridspec_kw=dict(width_ratios=[1.25, 1]))
    ax = axs[0]
    M = 1.0
    ax.add_patch(Circle((0, 0), 1.2, fc="k"))
    ax.text(0, 0, "M", color="w", ha="center", va="center")
    for rr in [13, 10, 7.5, 5.2]:
        k = 1.0 / rr**3 * 40
        ax.add_patch(Ellipse((rr, 0), 1.2 * (1 + 2 * k), 1.2 * (1 - k), fc=BLUE, alpha=0.35, ec=BLUE))
        ax.add_patch(Circle((rr, 0), 0.6, fc="none", ec=GREY, ls=":"))
    ax.annotate("", (4.0, -1.5), (13.8, -1.5), arrowprops=dict(arrowstyle="->", color=GREY))
    ax.text(9, -1.95, "falling in with the river (dotted: starting shape)", ha="center", fontsize=9, color=GREY)
    ax.text(9, 1.35, "stretched along r by 2M/r³, squeezed sideways by M/r³", ha="center", fontsize=9.5)
    ax.set_xlim(-1.6, 14.5); ax.set_ylim(-2.4, 2.0); ax.set_aspect("equal"); ax.axis("off")
    ax.set_title("a ball of cells, carried by the river: the tide (stretch exaggerated)")
    ax = axs[1]
    r = np.linspace(2.2, 10, 300)
    v = np.sqrt(2 * M / r)
    dv = -v / (2 * r); d2v = 3 * v / (4 * r**2)
    jet_r = -(v * d2v + dv**2)
    jet_s = -v * dv / r
    ax.plot(r, jet_r, color=RED, lw=2.5, label="second jet, radial: −(v v'' + v'²)")
    ax.plot(r, jet_s, color=BLUE, lw=2.5, label="second jet, sideways: −v v'/r")
    rs = r[::25]
    ax.plot(rs, -2 * M / rs**3, "ko", mfc="none", label="Riemann tide −2M/r³ (radial)")
    ax.plot(rs, M / rs**3, "ks", mfc="none", label="Riemann tide +M/r³ (sideways)")
    ax.axhline(0, color="k", lw=0.5)
    ax.set_xlabel("radius r (units of M)"); ax.set_ylabel("tidal rate")
    ax.legend(fontsize=9); ax.set_title("curvature = the jet's second derivative\n(trace: Newton's ∇²Φ with Φ = −v²/2)")
    save(fig, "fig4_tides_second_jet.png")


# ---------------------------------------------------------------- fig5 spin: velocity works to first order
def fig5():
    fig, axs = plt.subplots(1, 2, figsize=(12, 5.4))
    ax = axs[0]
    M, J = 1.0, 12.0                                   # swirl exaggerated to be visible
    ax.add_patch(Circle((0, 0), 2 * M, fc="k"))
    for r0 in [9.5]:
        for ph0 in np.linspace(0, 2 * np.pi, 10, endpoint=False):
            rr, ph = [r0], [ph0]
            dt = 0.02
            while rr[-1] > 2.05 * M:
                v = np.sqrt(2 * M / rr[-1]); w = 2 * J / rr[-1]**3
                rr.append(rr[-1] - v * dt); ph.append(ph[-1] + w * dt)
            rr, ph = np.array(rr), np.array(ph)
            ax.plot(rr * np.cos(ph), rr * np.sin(ph), color=BLUE, lw=1)
    ax.set_aspect("equal"); ax.axis("off"); ax.set_xlim(-10, 10); ax.set_ylim(-10, 10)
    ax.set_title("spinning mass, first order in spin:\ninflow v = √(2M/r) plus swirl w = 2J/r³ (velocity only; swirl exaggerated)")
    ax.text(0, -10.6, "passes every field equation at first order (J7c)", ha="center", fontsize=9.5, color=GREEN)
    ax = axs[1]
    th = np.linspace(0, 2 * np.pi, 400)
    f = 3 * np.sin(th)**2 - 2
    rad = 1 + 0.3 * f
    ax.plot(rad * np.sin(th), rad * np.cos(th), color=PURPLE, lw=2.2)
    ax.plot(np.sin(th), np.cos(th), "k:", lw=1)
    ax.annotate("", (0, 1.45), (0, -1.45), arrowprops=dict(arrowstyle="<-", color=GREY))
    ax.text(0.05, 1.5, "spin axis", fontsize=9, color=GREY)
    ax.text(1.35, 0.0, "equator:\n3R > 0", fontsize=9, color=PURPLE)
    ax.text(0.12, 0.22, "poles:\n3R < 0", fontsize=9, color=PURPLE)
    ax.set_aspect("equal"); ax.axis("off"); ax.set_xlim(-1.8, 1.8); ax.set_ylim(-1.6, 1.7)
    ax.set_title("second order in spin: the lattice itself curves\n3R = 2Ma²(3 sin²θ − 2)/r⁵ on Doran's Kerr river (J7e)")
    save(fig, "fig5_spin.png")


# ---------------------------------------------------------------- fig6 inputs before and after
def fig6():
    fig, ax = plt.subplots(figsize=(12, 5.6)); ax.axis("off")

    def b(x, y, w, h, t, fc):
        ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.02", fc=fc, ec="k"))
        ax.text(x + w / 2, y + h / 2, t, ha="center", va="center", fontsize=9.3, wrap=True)
    ax.text(0.17, 0.95, "before (round 3: river rebuilt on the velocity rule)", ha="center", weight="bold")
    b(0.0, 0.66, 0.34, 0.2, "velocity rule for ONE free faller\n(ours: dial addition, 1/r²)", "#e6f2ff")
    b(0.0, 0.38, 0.34, 0.2, "bridge: Newton's t is proper time,\nNewton's r is the areal radius", "#ffe6e6")
    b(0.0, 0.10, 0.34, 0.2, "Einstein's equation to confirm\nthe result is exact", "#ffe6e6")
    ax.text(0.68, 0.95, "now (this round: the river as a jet field)", ha="center", weight="bold")
    b(0.44, 0.66, 0.46, 0.2, "the rule for EVERY cell at once: each cell a zero-energy\nfree faller, v² = 2Gm(r)/r  ⇔  Q layer = 8πGρ", "#e6f2ff")
    b(0.44, 0.38, 0.46, 0.2, "still an input, now one structural sentence: the cells form\nthe theory's flat cube lattice, each with its own clock", "#fff7e6")
    b(0.44, 0.10, 0.46, 0.2, "Einstein's equation not borrowed: for the river every\ncomponent of it IS the Q layer or its slope (checked, J2)", "#e8f8e8")
    for y in (0.76, 0.48, 0.20):
        ax.annotate("", (0.43, y), (0.35, y), arrowprops=dict(arrowstyle="->", lw=1.5))
    ax.text(0.95, 0.76, "gives: Schwarzschild,\nGauss law, Friedmann", fontsize=8.5, va="center")
    ax.text(0.95, 0.48, "fails for full spin\nand for waves", fontsize=8.5, va="center", color=RED)
    ax.text(0.95, 0.20, "tides and frame\ndragging come out", fontsize=8.5, va="center", color=GREEN)
    ax.set_xlim(-0.03, 1.12); ax.set_ylim(0.05, 1.0)
    save(fig, "fig6_inputs.png")


for f in (fig1, fig2, fig3, fig4, fig5, fig6):
    f()
print("figures written:", sorted(os.listdir(FIG)))
