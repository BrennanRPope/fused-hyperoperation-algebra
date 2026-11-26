"""Labeled pictures for the gravity-from-theory DRAFT. Writes figures/*.png."""
import json, math, os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle, Circle, FancyBboxPatch
from scipy.integrate import solve_ivp

HERE = os.path.dirname(os.path.abspath(__file__))
R = json.load(open(os.path.join(HERE, "results.json")))
OUT = os.path.join(HERE, "figures"); os.makedirs(OUT, exist_ok=True)
BLUE, ORANGE, RED, GREEN, GRAY, PURPLE = "#2a6fdb", "#e8833a", "#d23c3c", "#2e9e5b", "#777777", "#7b4fc9"
plt.rcParams.update({"font.size": 11, "axes.spines.top": False, "axes.spines.right": False})

def save(fig, name):
    fig.savefig(os.path.join(OUT, name), dpi=150, bbox_inches="tight"); plt.close(fig); print("wrote", name)

def orbit(p, x0, v0, T, n=3000):
    def f(t, y):
        r = math.hypot(y[0], y[1]); a = -1/r**p/r
        return [y[2], y[3], a*y[0], a*y[1]]
    s = solve_ivp(f, (0, T), list(x0) + list(v0), method="DOP853", rtol=1e-11, atol=1e-12, t_eval=np.linspace(0, T, n))
    return s.y[:2].T, s.y[2:].T

# ---------------------------------------------------------------- fig 1: flux through a growing box
fig, axs = plt.subplots(1, 2, figsize=(12, 5))
ax = axs[0]; ax.set_aspect("equal"); ax.axis("off")
for s, col in ((1, BLUE), (2, ORANGE), (3, GREEN)):
    ax.add_patch(Rectangle((-s, -s), 2*s, 2*s, fill=False, lw=2, ec=col))
    ax.text(s + 0.05, s - 0.3, f"side {2*s}: face area x{s*s}", color=col)
ax.plot(0, 0, "o", color=RED, ms=10); ax.text(0.12, 0.12, "mass", color=RED)
for t in np.linspace(0, 2*np.pi, 16, endpoint=False):
    ax.annotate("", xy=(3.3*np.cos(t), 3.3*np.sin(t)), xytext=(0.25*np.cos(t), 0.25*np.sin(t)),
                arrowprops=dict(arrowstyle="->", color=GRAY, lw=0.8))
ax.set_xlim(-3.6, 5.6); ax.set_ylim(-3.6, 3.6)
ax.set_title("Same number of field lines cross every box.\nFaces grow like r², so strength falls like 1/r²")
ax = axs[1]
rr = np.linspace(0.5, 4, 200)
for n, col in ((2, BLUE), (3, RED), (4, PURPLE)):
    ax.plot(rr, 1/rr**(n - 1), color=col, lw=2.5, label=f"{n} dimensions: faces ∝ r^{n-1}, pull ∝ 1/r^{n-1}")
ax.set_xlabel("distance r"); ax.set_ylabel("pull"); ax.set_ylim(0, 2.2)
ax.legend(frameon=False)
ax.set_title("The face layer of the cube's lift sets the falloff.\n3D box: (1+ka)(1+kb)(1+kc) → k² layer = ab+bc+ca = half the surface")
save(fig, "fig1_flux_box.png")

# ---------------------------------------------------------------- fig 2: hodographs and the dial
fig = plt.figure(figsize=(15, 5.2))
ax = fig.add_subplot(1, 3, 1); ax.set_aspect("equal")
cases = [((1, 0), (0, 1.15), BLUE, "ellipse (E < 0)"), ((1, 0), (0, math.sqrt(2)), GRAY, "parabola (E = 0)"),
         ((1, 0), (0, 1.7), ORANGE, "hyperbola (E > 0)")]
for x0, v0, col, lab in cases:
    X, V = orbit(2, x0, v0, 14 if col != BLUE else 20)
    ax.plot(X[:, 0], X[:, 1], color=col, lw=2, label=lab)
ax.plot(0, 0, "o", color=RED); ax.set_xlim(-3, 2); ax.set_ylim(-2.5, 2.5)
ax.legend(frameon=False, loc="lower left"); ax.set_title("Orbits under 1/r² (positions)")
ax = fig.add_subplot(1, 3, 2); ax.set_aspect("equal")
for x0, v0, col, lab in cases:
    X, V = orbit(2, x0, v0, 60 if col == BLUE else 200, n=8000)
    ax.plot(V[:, 0], V[:, 1], color=col, lw=2)
ax.plot(0, 0, "+", color="k"); ax.text(0.05, -0.15, "v = 0")
ax.set_xlabel("vₓ"); ax.set_ylabel("v_y")
ax.set_title("Their velocities trace circles\n(hodographs, Hamilton 1846)")
ax3 = fig.add_subplot(1, 3, 3, projection="3d")
X, V = orbit(2, (1, 0), (0, 1.15), 60, n=4000)
E = 0.5*np.sum(V[0]**2) - 1/np.linalg.norm(X[0]); c = -2*E
w = V/math.sqrt(c); w2 = np.sum(w**2, 1)[:, None]
S = np.hstack([2*w, w2 - 1])/(w2 + 1)
uu, vv = np.mgrid[0:2*np.pi:40j, 0:np.pi:20j]
ax3.plot_wireframe(np.cos(uu)*np.sin(vv), np.sin(uu)*np.sin(vv), np.cos(vv), color=GRAY, lw=0.3, alpha=0.5)
ax3.plot(S[:, 0], S[:, 1], S[:, 2], color=BLUE, lw=3, label="ellipse's velocities, dial κ = −1/(2E)")
wb = V/math.sqrt(2.5*c); w2b = np.sum(wb**2, 1)[:, None]
Sb = np.hstack([2*wb, w2b - 1])/(w2b + 1)
ax3.plot(Sb[:, 0], Sb[:, 1], Sb[:, 2], color=RED, lw=1.5, ls="--", label="same velocities, wrong dial")
ax3.set_box_aspect((1, 1, 1)); ax3.set_xticks([]); ax3.set_yticks([]); ax3.set_zticks([])
ax3.legend(frameon=False, loc="upper left", fontsize=9)
ax3.set_title("On the dial's sphere the circle is a great circle:\na straight line of the d = +1 geometry")
save(fig, "fig2_hodograph_dial.png")

# ---------------------------------------------------------------- fig 3: only 1/r² closes and has circular hodographs
fig, axs = plt.subplots(2, 3, figsize=(13, 8.5))
for j, (p, col, lab) in enumerate(((1.5, PURPLE, "pull 1/r^1.5  (2.5 dimensions)"), (2.0, BLUE, "pull 1/r²  (3 dimensions)"),
                                   (2.5, ORANGE, "pull 1/r^2.5  (3.5 dimensions)"))):
    X, V = orbit(p, (1, 0), (0, 0.8), 80, n=8000)
    ax = axs[0, j]; ax.set_aspect("equal"); ax.plot(X[:, 0], X[:, 1], color=col, lw=1)
    ax.plot(0, 0, "o", color=RED); ax.set_title(lab); ax.set_xticks([]); ax.set_yticks([])
    ax = axs[1, j]; ax.set_aspect("equal"); ax.plot(V[:, 0], V[:, 1], color=col, lw=1)
    res = R["N7"][str(p)]
    ax.set_title(f"velocities: circle misfit {res:.0e}"); ax.set_xticks([]); ax.set_yticks([])
axs[0, 0].set_ylabel("orbit (positions)"); axs[1, 0].set_ylabel("hodograph (velocities)")
fig.suptitle("Only the inverse-square pull closes its orbits and draws circles in velocity space (Bertrand 1873, Hamilton 1846)")
save(fig, "fig3_only_inverse_square.png")

# ---------------------------------------------------------------- fig 4: boosting Newton's field gives the drag term
fig, axs = plt.subplots(1, 2, figsize=(13, 5))
ax = axs[0]; ax.axis("off"); ax.set_xlim(0, 10); ax.set_ylim(0, 10)
ax.add_patch(FancyBboxPatch((0.2, 6.3), 3.6, 3, boxstyle="round,pad=0.1", fc="#eef3fb", ec=BLUE))
ax.text(0.4, 8.7, "mass at rest (Newton)", color=BLUE, weight="bold")
ax.text(0.4, 7.9, "clocks:   h_tt = 2U")
ax.text(0.4, 7.2, "rulers:   h_xx = 2γU/c²")
ax.text(0.4, 6.6, "no dt·dx term", color=GRAY)
ax.annotate("", xy=(6, 7.8), xytext=(4.1, 7.8), arrowprops=dict(arrowstyle="->", lw=2.5, color=RED))
ax.text(4.2, 8.2, "dial boost\nspeed u", color=RED)
ax.add_patch(FancyBboxPatch((6.2, 6.3), 3.6, 3, boxstyle="round,pad=0.1", fc="#fdf0e6", ec=ORANGE))
ax.text(6.4, 8.7, "same mass, moving", color=ORANGE, weight="bold")
ax.text(6.4, 7.9, "new dt·dx term:")
ax.text(6.4, 7.2, "h_tx = 2(1+γ) U u/c²", weight="bold")
ax.text(6.4, 6.6, "= the frame-drag field", color=GRAY)
ax.text(0.3, 4.8, "Where the pieces come from", weight="bold")
ax.text(0.3, 4.0, "U = GM/r   Newton's pull, put in (and G)")
ax.text(0.3, 3.3, "γ = 1        river page: one boost slows clocks and")
ax.text(0.3, 2.75, "                 shrinks rulers equally (mult-add inverses)")
ax.text(0.3, 2.1, "1/c² = −q   the plane dial: the boost's corner −q·u")
ax.text(0.3, 1.4, "boosts compose by (u+v)/(1 − quv): the theory's law")
ax.text(0.3, 0.5, "Galilean dial q = 0: the term vanishes, no frame dragging", color=RED)
ax = axs[1]
ax.bar([0], [2], color=BLUE, width=0.5, label="from clocks (h_tt)")
ax.bar([0], [2], bottom=[2], color=GREEN, width=0.5, label="from rulers (h_xx), γ = 1")
ax.bar([1], [2], color=BLUE, width=0.5)
ax.bar([1], [0.02], bottom=[2], color=GREEN, width=0.5)
ax.set_xticks([0, 1]); ax.set_xticklabels(["river (γ = 1)\nGR value", "dial-only theory\n(γ = 0)"])
ax.set_ylabel("size of h_tx in units of U·u/c²"); ax.set_ylim(0, 5)
ax.text(0, 4.15, "4", ha="center", weight="bold"); ax.text(1, 2.15, "2", ha="center", weight="bold")
ax.legend(frameon=False, loc="upper right")
ax.set_title("Frame-drag weight = clocks + rulers:\nthe same split that doubles light bending")
save(fig, "fig4_boost_newton.png")

# ---------------------------------------------------------------- fig 5: spinning ball, gyroscope turn field
fig, axs = plt.subplots(1, 2, figsize=(13, 5.6))
ax = axs[0]; ax.set_aspect("equal")
xg, zg = np.meshgrid(np.linspace(-3, 3, 300), np.linspace(-3, 3, 300))
r = np.sqrt(xg**2 + zg**2); r[r < 1] = np.nan
Ox = 3*zg*xg/r**5; Oz = (3*zg*zg/r**2 - 1)/r**3
ax.streamplot(xg, zg, Ox, Oz, color=np.log(np.hypot(Ox, Oz)), cmap="Blues", density=1.4, linewidth=1)
ax.add_patch(Circle((0, 0), 1, color=ORANGE, alpha=0.7)); ax.text(-0.55, -0.1, "spinning\nmass", fontsize=10)
ax.annotate("", xy=(0, 1.6), xytext=(0, 0.2), arrowprops=dict(arrowstyle="->", lw=3, color=RED)); ax.text(0.1, 1.4, "J", color=RED)
ax.plot([0], [2.2], "o", color=GREEN); ax.text(0.15, 2.25, "on the axis: turns with J at 2GJ/(c²r³)", color=GREEN, fontsize=9)
ax.plot([2.2], [0], "o", color=PURPLE); ax.text(1.2, -0.45, "at the equator: turns\nagainst J at GJ/(c²r³)", color=PURPLE, fontsize=9)
ax.set_xlim(-3, 3); ax.set_ylim(-3, 3); ax.set_xlabel("x"); ax.set_ylabel("z (spin axis)")
ax.set_title("Gyroscope turn rate = half the circulation of the\nboosted field per enclosed area (curl), a dipole")
ax = axs[1]
rs = np.linspace(5, 40, 50)
ax.plot(rs, np.ones_like(rs)*2, color=GREEN, lw=2.5, label="ours, axis: 2J/r³")
ax.plot(rs, np.ones_like(rs)*1, color=PURPLE, lw=2.5, label="ours, equator: J/r³")
import sympy as sp
def kerr(rv, thv, Mv=1.0, av=0.5):
    Mk, ak, r_, t_ = sp.symbols("Mk ak r_ t_", positive=True)
    D = r_**2 - 2*Mk*r_ + ak**2; Sg = r_**2 + ak**2*sp.cos(t_)**2
    gtt = -(1 - 2*Mk*r_/Sg); gtp = -2*Mk*ak*r_*sp.sin(t_)**2/Sg
    gpp = (r_**2 + ak**2 + 2*Mk*ak**2*r_*sp.sin(t_)**2/Sg)*sp.sin(t_)**2
    N2 = -gtt; A = gtp/N2; hpp = gpp + N2*A**2
    e = sp.sqrt(N2)/2*sp.sqrt(sp.diff(A, r_)**2*D/(Sg*hpp) + sp.diff(A, t_)**2/(Sg*hpp))
    return sp.lambdify((r_, t_), e.subs({Mk: Mv, ak: av}), "numpy")
fk = kerr(0, 0)
rk = np.array([6, 8, 12, 18, 25, 35, 40.0])
ax.plot(rk, [fk(x, 1e-4)*x**3/0.5 for x in rk], "*", color=GREEN, ms=12, label="exact Kerr (GR), axis")
ax.plot(rk, [fk(x, np.pi/2)*x**3/0.5 for x in rk], "*", color=PURPLE, ms=12, label="exact Kerr (GR), equator")
ax.set_xlabel("distance r (in units of GM/c²)"); ax.set_ylabel("turn rate × r³ / J"); ax.set_ylim(0, 3)
ax.legend(frameon=False); ax.set_title("Check against the exact spinning black hole:\nKerr approaches our field far out")
save(fig, "fig5_spinning_ball.png")

# ---------------------------------------------------------------- fig 6: data
f7 = R["F7"]
fig, axs = plt.subplots(1, 3, figsize=(14, 4.6))
items = [("Gravity Probe B\nframe dragging (mas/yr)", f7["gpb_proj"], f7["gpb_gamma0"], 37.2, 7.2, 7.2),
         ("LAGEOS node drag\n(% of GR value)", 100, 50, 99, 5, 5),
         ("double pulsar B\nspin precession (°/yr)", f7["dp1"], f7["dp0"], 4.77, 0.65, 0.66)]
for ax, (lab, g1, g0, m, el, eh) in zip(axs, items):
    ax.bar([0], [g1], color=BLUE, width=0.5); ax.bar([1], [g0], color=GRAY, width=0.5)
    ax.errorbar([2], [m], yerr=[[el], [eh]], fmt="o", color="k", capsize=6, ms=8)
    ax.set_xticks([0, 1, 2]); ax.set_xticklabels(["ours\n(γ = 1)", "dial only\n(γ = 0)", "measured"])
    ax.set_title(lab)
fig.suptitle("Same field, three measurements: the clocks-plus-rulers weight is needed")
save(fig, "fig6_data.png")

# ---------------------------------------------------------------- fig 7: map of derived vs put in
fig, ax = plt.subplots(figsize=(12, 6)); ax.axis("off"); ax.set_xlim(0, 12); ax.set_ylim(0, 7)
def box(x, y, w, h, txt, fc, ec):
    ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.1", fc=fc, ec=ec, lw=1.8))
    ax.text(x + w/2, y + h/2, txt, ha="center", va="center", fontsize=10)
box(0.2, 5.2, 3.3, 1.3, "PUT IN\nNewton's pull U = GM/r\n(strength G, and 3 dimensions)", "#fde8e8", RED)
box(4.3, 5.2, 3.3, 1.3, "FROM THE THEORY\nplane dial q = −1/c²\nboosts compose by ∘", "#e6f4ea", GREEN)
box(8.4, 5.2, 3.3, 1.3, "FROM AN EARLIER PAGE\nriver: γ = 1\n(clocks and rulers as inverses)", "#eef3fb", BLUE)
box(4.3, 2.9, 3.3, 1.3, "boosted field of moving mass\nh_tx = 2(1+γ)Uu/c²", "#fff6e0", ORANGE)
box(0.2, 0.5, 3.3, 1.3, "gyroscope turn = ½ curl\n(swap rule: turn ∝ area)\nLense–Thirring, GP-B", "#fff6e0", ORANGE)
box(4.3, 0.5, 3.3, 1.3, "orbit node drag 2J/a³\nLAGEOS", "#fff6e0", ORANGE)
box(8.4, 0.5, 3.3, 1.3, "drag weight 2m₁/M\ndouble pulsar", "#fff6e0", ORANGE)
for x0 in (1.85, 5.95, 10.05):
    ax.annotate("", xy=(5.95, 4.3), xytext=(x0, 5.1), arrowprops=dict(arrowstyle="->", lw=1.6, color=GRAY))
for x1 in (1.85, 5.95, 10.05):
    ax.annotate("", xy=(x1, 1.9), xytext=(5.95, 2.8), arrowprops=dict(arrowstyle="->", lw=1.6, color=GRAY))
ax.set_title("Frame dragging now needs nothing from gravity beyond Newton's static pull")
save(fig, "fig7_map.png")
