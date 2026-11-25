"""Labeled pictures for the precession DRAFT. Reads results.json from checks.py, writes figures/*.png."""
import json, math, os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch, Circle, Polygon
from mpl_toolkits.mplot3d import proj3d
from scipy.integrate import solve_ivp

HERE = os.path.dirname(os.path.abspath(__file__))
R = json.load(open(os.path.join(HERE, "results.json")))
OUT = os.path.join(HERE, "figures")
os.makedirs(OUT, exist_ok=True)
BLUE, ORANGE, RED, GREEN, GRAY, PURPLE = "#2a6fdb", "#e8833a", "#d23c3c", "#2e9e5b", "#777777", "#7b4fc9"
plt.rcParams.update({"font.size": 11, "axes.spines.top": False, "axes.spines.right": False})


class Arrow3D(FancyArrowPatch):
    def __init__(self, xs, ys, zs, *args, **kwargs):
        super().__init__((0, 0), (0, 0), *args, **kwargs)
        self._verts3d = xs, ys, zs

    def do_3d_projection(self, renderer=None):
        xs3d, ys3d, zs3d = self._verts3d
        xs, ys, zs = proj3d.proj_transform(xs3d, ys3d, zs3d, self.axes.M)
        self.set_positions((xs[0], ys[0]), (xs[1], ys[1]))
        return np.min(zs)


def arrow3(ax, p, q, color, lw=2.5):
    ax.add_artist(Arrow3D([p[0], q[0]], [p[1], q[1]], [p[2], q[2]], mutation_scale=16, lw=lw, arrowstyle="-|>", color=color))


def save(fig, name):
    fig.savefig(os.path.join(OUT, name), dpi=150, bbox_inches="tight")
    plt.close(fig)
    print("wrote", name)


# ---------------------------------------------------------------- fig 1: swap = turn
fig = plt.figure(figsize=(12, 5.2))
ax = fig.add_subplot(1, 2, 1, projection="3d")
u = np.array([1.0, 0.15, 0]); v = np.array([0.35, 1.0, 0]); c = np.cross(u, v)
O = np.zeros(3)
arrow3(ax, O, u, BLUE); arrow3(ax, O, v, ORANGE); arrow3(ax, O, c, RED, 3)
ax.text(*(u * 1.08), "u (first input)", color=BLUE)
ax.text(*(v * 1.08), "v (second input)", color=ORANGE)
ax.text(*(c * 1.05 + np.array([0.05, 0, 0])), "a∘b − b∘a = 2k (u × v)\nthe swap points here,\n90° to both inputs", color=RED)
xx, yy = np.meshgrid(np.linspace(-0.2, 1.3, 2), np.linspace(-0.2, 1.3, 2))
ax.plot_surface(xx, yy, 0 * xx, alpha=0.08, color=GRAY)
ax.text(1.2, -0.2, 0, "plane of u and v", color=GRAY)
ax.set_xlim(-0.2, 1.3); ax.set_ylim(-0.2, 1.3); ax.set_zlim(0, 1.0)
ax.set_xticks([]); ax.set_yticks([]); ax.set_zticks([])
ax.view_init(22, -60)
ax.set_title("Two directions: swapping the order leaves a turn\nat right angles to both (the corner's cross part)")
ax2 = fig.add_subplot(1, 2, 2)
ax2.set_xlim(0, 10); ax2.set_ylim(0, 10); ax2.axis("off")
ax2.set_title("Swap = turn, for any inputs")
ax2.text(0.2, 8.7, "a ∘ₖ b = a + b + k·ab", fontsize=14)
ax2.text(0.2, 7.6, "1 + k(a ∘ b) = (1 + ka)(1 + kb)   (still true when ab ≠ ba)", fontsize=11)
ax2.text(0.2, 6.3, "b ∘ a  =  X⁻¹ (a ∘ b) X,   X = 1 + ka", fontsize=14, color=RED)
ax2.text(0.2, 5.4, "swapping the inputs = turning the answer by the first input", fontsize=11, color=RED)
ax2.text(0.2, 4.1, "same size: |1 + k(a∘b)| = |1 + k(b∘a)|  (Euler's four squares)", fontsize=11)
ax2.annotate("", xy=(6.2, 1.6), xytext=(0.6, 1.6), arrowprops=dict(arrowstyle="-|>", color=BLUE, lw=2.5))
ax2.annotate("", xy=(9.2, 1.6), xytext=(6.2, 1.6), arrowprops=dict(arrowstyle="-|>", color=ORANGE, lw=2.5))
ax2.text(0.4, 2.2, "inputs on one line commute: no turn, so a one-dimensional", fontsize=10.5, color=GRAY)
ax2.text(0.4, 0.6, "theory can never precess. Precession needs two directions.", fontsize=10.5, color=GRAY)
save(fig, "fig1_swap_is_turn.png")


# ---------------------------------------------------------------- fig 2: one law, two ends
def slerp(p, q, n=60):
    om = math.acos(np.clip(p @ q, -1, 1))
    s = np.linspace(0, 1, n)[:, None]
    return (np.sin((1 - s) * om) * p + np.sin(s * om) * q) / math.sin(om)


fig = plt.figure(figsize=(12.5, 5.6))
ax = fig.add_subplot(1, 2, 1, projection="3d")
uu, vv = np.mgrid[0:2 * np.pi:40j, 0:np.pi:20j]
ax.plot_wireframe(np.cos(uu) * np.sin(vv), np.sin(uu) * np.sin(vv), np.cos(vv), color=GRAY, lw=0.3, alpha=0.4)
A = np.array([0, 0, 1.0]); B = np.array([math.sin(1.0), 0, math.cos(1.0)])
C = np.array([math.sin(1.1) * math.cos(1.2), math.sin(1.1) * math.sin(1.2), math.cos(1.1)])
for p, q in ((A, B), (B, C), (C, A)):
    seg = slerp(p, q); ax.plot(seg[:, 0], seg[:, 1], seg[:, 2], color=BLUE, lw=2.5)
pts = np.vstack([slerp(A, B), slerp(B, C), slerp(C, A)])
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
ax.add_collection3d(Poly3DCollection([pts * 1.001], alpha=0.35, color=BLUE))
area = 2 * math.atan2(A @ np.cross(B, C), 1 + A @ B + B @ C + C @ A)
ax.text(0.1, 0.1, 1.25, f"area = {area:.3f}\ntwist = +1 × area", color=BLUE)
ax.text(*(A * 1.08), "rest", color="k")
ax.set_box_aspect((1, 1, 1)); ax.axis("off"); ax.view_init(25, -40)
ax.set_title("Turns: dial d = +1 (sphere)\nFoucault, rigid-body phase, Berry (half of it)")
ax = fig.add_subplot(1, 2, 2)
ax.add_patch(Circle((0, 0), 1, fill=False, color="k", lw=1.5))
ax.text(0.62, 0.86, "light speed", fontsize=10)


def mob(p, z):
    return (z - p) / (1 - np.conj(p) * z)


def geo(p, q, n=80):
    w = mob(p, q)
    s = np.linspace(0, 1, n)
    zz = s * w
    return (zz + p) / (1 + np.conj(p) * zz)


P0, P1, P2 = 0j, 0.62 + 0j, 0.38 + 0.55j
segs = [geo(P0, P1), geo(P1, P2), geo(P2, P0)]
poly = np.concatenate(segs)
ax.add_patch(Polygon(np.c_[poly.real, poly.imag], alpha=0.3, color=ORANGE))
for s in segs:
    ax.plot(s.real, s.imag, color=ORANGE, lw=2.5)
ax.text(-0.15, -0.1, "rest", fontsize=10)
ax.text(0.63, -0.1, "u", fontsize=11); ax.text(0.4, 0.6, "u ⊕ v", fontsize=11)
ax.text(-0.95, -1.25, "twist = −1 × area  (Thomas–Wigner rotation;\nworked out on the twisted-boosts page, checked against GP-B and a pulsar)", fontsize=9.5)
ax.set_aspect("equal"); ax.set_xlim(-1.1, 1.1); ax.set_ylim(-1.35, 1.1); ax.axis("off")
ax.set_title("Boosts: dial d = −1 (velocity disk)\nThomas, geodetic (GP-B), double pulsar")
fig.suptitle("One composition law, two ends of the plane dial: precession = d × enclosed area", fontsize=13, y=1.06)
save(fig, "fig2_one_law_two_ends.png")

# ---------------------------------------------------------------- fig 3: holonomy family
fig, ax = plt.subplots(figsize=(9, 5.5))
rho = np.linspace(0, 1.6, 300)
cols = {1.0: BLUE, 0.5: "#4f8be6", 0.25: "#8fb6f2", 0.0: "k", -0.25: "#f5be93", -0.5: "#f0a46c", -1.0: ORANGE}
for d, col in cols.items():
    if d > 0:
        y = 2 * np.pi * (1 - np.cos(np.sqrt(d) * rho))
    elif d < 0:
        y = 2 * np.pi * (1 - np.cosh(np.sqrt(-d) * rho))
    else:
        y = 0 * rho
    ax.plot(rho, np.degrees(y), color=col, lw=2.2, label=f"d = {d:+g}" + ("  (Galilean: never precesses)" if d == 0 else ""))
for c in R["C_family"]:
    ax.plot(c["rho"], math.degrees(c["holonomy"]), "o", color="k", ms=5)
ax.plot([], [], "o", color="k", label="summed from ∘ twists (checks C4)")
pr = math.radians(90 - 48.85)
ax.plot(pr, math.degrees(2 * np.pi * (1 - math.cos(pr))), "*", color=BLUE, ms=15)
ax.annotate("Paris latitude: cap holonomy 89°;\nseen turn per day = 360° − 89° = 271°", (pr, math.degrees(2 * np.pi * (1 - math.cos(pr)))),
            (0.95, 330), arrowprops=dict(arrowstyle="->"), fontsize=9.5)
et = math.atanh(0.5)
ax.plot(et, math.degrees(2 * np.pi * (1 - math.cosh(et))), "*", color=ORANGE, ms=15)
ax.annotate("orbit at half light speed:\nThomas turn −55.7° per lap", (et, math.degrees(2 * np.pi * (1 - math.cosh(et)))), (0.7, -230),
            arrowprops=dict(arrowstyle="->"), fontsize=9.5)
ax.axhline(0, color=GRAY, lw=0.5)
ax.set_xlabel("radius ρ of the loop (angle on the sphere, or rapidity)")
ax.set_ylabel("turn per loop (degrees)")
ax.set_title("Turn per loop = 2π(1 − cos(√d ρ)): one curve family as the plane dial d moves")
ax.legend(loc="lower left", fontsize=9)
save(fig, "fig3_holonomy_family.png")

# ---------------------------------------------------------------- fig 4: Foucault
fig = plt.figure(figsize=(14, 5.6))
fig.subplots_adjust(wspace=0.3)
ax = fig.add_subplot(1, 2, 1, projection="3d")
ax.plot_wireframe(np.cos(uu) * np.sin(vv), np.sin(uu) * np.sin(vv), np.cos(vv), color=GRAY, lw=0.3, alpha=0.35)
lat = math.radians(48.85)
phs = np.linspace(0, 2 * np.pi, 200)
ax.plot(np.cos(lat) * np.cos(phs), np.cos(lat) * np.sin(phs), np.sin(lat) * np.ones_like(phs), color=BLUE, lw=2.5)
capv = np.linspace(0, np.pi / 2 - lat, 12)
for cv in capv:
    ax.plot(np.sin(cv) * np.cos(phs), np.sin(cv) * np.sin(phs), np.cos(cv) * np.ones_like(phs), color=BLUE, lw=0.4, alpha=0.5)
for k_ in range(9):
    p = 2 * np.pi * k_ / 8
    pos = np.array([math.cos(lat) * math.cos(p), math.cos(lat) * math.sin(p), math.sin(lat)])
    north = np.array([-math.sin(lat) * math.cos(p), -math.sin(lat) * math.sin(p), math.cos(lat)])
    east = np.array([-math.sin(p), math.cos(p), 0])
    al = -p * math.sin(lat)
    dvec = math.cos(al) * north + math.sin(al) * east
    s = 0.22
    ax.plot(*np.c_[pos - s * dvec, pos + s * dvec], color=RED, lw=2.5)
ax.quiver(0, 0, 0, 0, 0, 1.4, color="k", lw=1.5, arrow_length_ratio=0.08)
ax.text(0, 0, 1.45, "Earth's axis")
ax.text(-0.9, -0.9, 1.25, "shaded cap = enclosed area", color=BLUE)
ax.text(0.9, 0.9, -0.3, "red bars: swing plane\ncarried round one day\n(it lags the ground)", color=RED)
ax.set_box_aspect((1, 1, 1)); ax.axis("off"); ax.view_init(28, -55)
ax.set_title("Foucault's pendulum at Paris (48.85° N)")
ax = fig.add_subplot(1, 2, 2)
la = np.linspace(0, 90, 200)
ax.plot(la, 360 * np.sin(np.radians(la)), color=BLUE, lw=2.2, label="360° × sin(latitude)")
F = R["C_foucault"]
ax.plot([f["lat"] for f in F], [f["per_day_deg"] for f in F], "o", color="k", label="360° − (∘-twist sum around the latitude), checks C5")
ax.plot(48.85, 360 * math.sin(lat), "*", color=RED, ms=16, label="Paris: 271° per day = 11.3° per hour")
ax.set_xlabel("latitude (degrees)"); ax.set_ylabel("swing plane turn per sidereal day (degrees)")
ax.set_title("Turn per day from the d = +1 twist law alone")
ax.legend(fontsize=9, loc="lower right")
save(fig, "fig4_foucault.png")

# ---------------------------------------------------------------- fig 5: gyroscope
E = R["E"]
fig = plt.figure(figsize=(15, 6.2))
fig.subplots_adjust(wspace=0.25)
ax = fig.add_subplot(1, 2, 1, projection="3d")
ax.plot([0, 0], [0, 0], [-0.9, 0], color="k", lw=4)  # stand
ax.text(-0.1, -0.3, -1.0, "stand (pivot on top)")
ax.plot([0, 1.1], [0, 0], [0, 0], color="k", lw=3)
th = np.linspace(0, 2 * np.pi, 100)
ax.plot(np.ones_like(th), 0.45 * np.cos(th), 0.45 * np.sin(th), color=BLUE, lw=2.5)
ax.add_collection3d(Poly3DCollection([np.c_[np.ones_like(th), 0.45 * np.cos(th), 0.45 * np.sin(th)]], alpha=0.2, color=BLUE))
ax.text(0.9, -0.75, 0.55, "spin plane\n(the wheel)", color=BLUE)
arrow3(ax, [1, 0, 0], [1.8, 0, 0], BLUE)
ax.text(1.85, -0.15, -0.15, "L (spin)", color=BLUE)
arrow3(ax, [1, 0, 0], [1, 0, -0.7], GRAY)
ax.text(1.05, 0, -0.85, "gravity m g", color=GRAY)
ax.add_collection3d(Poly3DCollection([np.array([[0, 0, 0], [1.1, 0, 0], [1.1, 0, -0.7], [0, 0, -0.7]])], alpha=0.15, color=GRAY))
ax.text(0.05, 0, -0.45, "torque plane\n(axle and gravity)", color="#444444")
arrow3(ax, [1.8, 0, 0], [1.8, 0.7, 0], RED, 3)
ax.text(1.75, 0.8, 0.15, "dL = torque × dt:\nhorizontal, 90° to L", color=RED)
ang = np.linspace(-0.3, 1.2, 60)
ax.plot(1.3 * np.cos(ang), 1.3 * np.sin(ang), 0.75 * np.ones_like(ang), color=GREEN, lw=3)
arrow3(ax, [1.3*np.cos(1.1), 1.3*np.sin(1.1), 0.75], [1.3*np.cos(1.2), 1.3*np.sin(1.2), 0.75], GREEN)
ax.plot([0, 0], [0, 0], [0, 0.9], color=GREEN, lw=1, ls=":")
ax.text(0.1, 0.9, 0.9, "precession: the axle swings round the vertical", color=GREEN)
ax.set_xlim(-0.2, 2); ax.set_ylim(-0.8, 1); ax.set_zlim(-0.9, 0.9)
ax.axis("off"); ax.view_init(20, -65)
ax.set_title("Why the push lands at 90°: the spin plane and the torque plane\nshare the vertical line, so the wheel's plane swings about it")
ax = fig.add_subplot(1, 2, 2)
x = np.array(E["trace_x"]); y = np.array(E["trace_y"]); z = np.array(E["trace_z"])
az = np.degrees(np.unwrap(np.arctan2(y, x)))
ax.plot(az, np.degrees(np.arcsin(-z)), color=PURPLE, lw=1.3, label="released with no push sideways: loops (nutation)")
ax.axhline(math.degrees(E["pred_dip"]), color=RED, ls="--", lw=1, label=f"predicted depth 2·I₁·m·g·l/(I₃ω₃)² = {math.degrees(E['pred_dip']):.3f}°")
ax.axhline(0, color=GREEN, lw=2, label="started at the steady rate: stays level, no loops")
ax.invert_yaxis()
ax.set_xlabel("how far round the axle has precessed (degrees)")
ax.set_ylabel("axle dip below level (degrees)")
ax.set_title(f"Simulated by composing small turns with ∘ (checks E2–E6)\nprecession rate m g l/(I₃ ω₃) = {E['Wp']:.3f} rad/s, matched")
ax.legend(fontsize=8.5, loc="center right", framealpha=0.95)
save(fig, "fig5_gyroscope.png")

# ---------------------------------------------------------------- fig 6: rigid body (Montgomery)
I = np.array([1.0, 2.0, 3.0])


def euler_rhs(t, w):
    return [(I[1] - I[2]) / I[0] * w[1] * w[2], (I[2] - I[0]) / I[1] * w[2] * w[0], (I[0] - I[1]) / I[2] * w[0] * w[1]]


fig = plt.figure(figsize=(13, 5.6))
ax = fig.add_subplot(1, 2, 1, projection="3d")
ax.plot_wireframe(np.cos(uu) * np.sin(vv), np.sin(uu) * np.sin(vv), np.cos(vv), color=GRAY, lw=0.3, alpha=0.3)
cols = [BLUE, "#6aa0ee", ORANGE, "#f0a46c"]
for m, col in zip(R["D_montgomery"], cols):
    sol = solve_ivp(euler_rhs, (0, m["T"]), m["w0"], rtol=1e-11, atol=1e-12, dense_output=True)
    ts = np.linspace(0, m["T"], 400)
    L = (I[:, None] * sol.sol(ts)); L = L / np.linalg.norm(L, axis=0)
    ax.plot(L[0], L[1], L[2], color=col, lw=2.5)
ax.text(1.15, 0, 0, "axis 1 (smallest I)"); ax.text(0, 0, 1.15, "axis 3 (largest I)")
ax.set_box_aspect((1, 1, 1)); ax.axis("off"); ax.view_init(20, 35)
ax.set_title("A tumbling free body: its angular momentum,\nseen from the body, draws loops on the sphere")
ax = fig.add_subplot(1, 2, 2)
labels = [f"run {i+1}\n(axis {m['axis']+1})" for i, m in enumerate(R["D_montgomery"])]
xs = np.arange(4)
ax.bar(xs - 0.25, [m["twoET_over_L"] for m in R["D_montgomery"]], 0.25, color=GRAY, label="dynamic part 2ET/|L|")
ax.bar(xs, [-m["Lambda"] for m in R["D_montgomery"]], 0.25, color=BLUE, label="geometric part −Λ (enclosed area, from ∘ twists)")
tot = [((m["twoET_over_L"] - m["Lambda"] + math.pi) % (2 * math.pi)) - math.pi for m in R["D_montgomery"]]
meas = [m["dtheta"] for m in R["D_montgomery"]]
ax.bar(xs + 0.25, meas, 0.25, color=RED, label="measured turn, by composing turns (mod 2π)")
ax.plot(xs + 0.25, tot, "k_", ms=22, mew=2, label="sum of the two parts (mod 2π)")
ax.set_xticks(xs); ax.set_xticklabels(labels)
ax.set_ylabel("radians"); ax.axhline(0, color="k", lw=0.5)
ax.set_title("Montgomery's phase: turn per loop = 2ET/|L| − Λ\n(all four runs match to 2×10⁻⁸ rad; checks D2)")
ax.legend(fontsize=8.5)
save(fig, "fig6_rigid_body_phase.png")

# ---------------------------------------------------------------- fig 7: Larmor and Berry
fig = plt.figure(figsize=(12.5, 5.4))
ax = fig.add_subplot(1, 2, 1, projection="3d")
ax.plot_wireframe(np.cos(uu) * np.sin(vv), np.sin(uu) * np.sin(vv), np.cos(vv), color=GRAY, lw=0.3, alpha=0.3)
arrow3(ax, [0, 0, 0], [0, 0, 1.35], GREEN, 3)
ax.text(0.05, 0, 1.4, "magnetic field B", color=GREEN)
cth = 0.8
ax.plot(math.sin(cth) * np.cos(th), math.sin(cth) * np.sin(th), math.cos(cth) * np.ones_like(th), color=RED, lw=2)
arrow3(ax, [0, 0, 0], [math.sin(cth), 0, math.cos(cth)], BLUE, 3)
ax.text(math.sin(cth) + 0.05, 0, math.cos(cth), "spin S", color=BLUE)
arrow3(ax, [math.sin(cth), 0, math.cos(cth)], [math.sin(cth), 0.45, math.cos(cth)], RED)
ax.text(math.sin(cth), 0.5, math.cos(cth) + 0.1, "dS ∝ S × B:\n90° to S and to B", color=RED)
ax.set_box_aspect((1, 1, 1)); ax.axis("off"); ax.view_init(18, -60)
ax.set_title("Larmor precession: the spin circles the field at γB\n(length and tilt kept; checks F1, F2)")
ax = fig.add_subplot(1, 2, 2)
t0 = np.linspace(0, np.pi / 2, 200)
ax.plot(np.degrees(t0), -np.pi * (1 - np.cos(t0)), color=BLUE, lw=2.2, label="−½ × solid angle = −½ × (d = +1 twist)")
ax.plot(np.degrees(t0), -2 * np.pi * (1 - np.cos(t0)), color=GRAY, ls="--", lw=1.5, label="full twist (what a whole turn would give)")
for b in R["F"]:
    ax.plot(math.degrees(b["theta0"]), b["berry"], "o", color=RED, ms=8)
ax.plot([], [], "o", color=RED, label="spinor steps composed with ∘ (checks F3)")
ax.set_xlabel("cone angle of the field's slow loop (degrees)"); ax.set_ylabel("Berry phase (radians)")
ax.set_title("Berry's phase is half the twist:\nthe spinor lives on the half-angle chart")
ax.legend(fontsize=9, loc="lower left")
save(fig, "fig7_larmor_berry.png")

# ---------------------------------------------------------------- fig 8: BTZ frame dragging
G = R["G"]
fig, axs = plt.subplots(1, 2, figsize=(13, 5.2))
ax = axs[0]
r = np.array(G["r"])
for J, col in zip(G["J"], [BLUE, PURPLE, ORANGE]):
    ax.plot(r, G["w_static"][str(J)], color=col, lw=2.2, label=f"gyroscope held still, J = {J}")
    ax.plot(r, G["zamo"][str(J)], color=col, lw=1.2, ls="--")
ax.plot([], [], color="k", ls="--", lw=1.2, label="drag angular velocity J/(2r²)")
ax.axvline(1.0, color=GRAY, lw=1)
ax.text(1.12, 1.45, "r² = M: inside this,\nnothing can stand still", fontsize=9, color=GRAY, bbox=dict(fc="white", ec="none"))
fw = G["fw_check"]
ax.plot(fw["r0"], abs(fw["numeric"]) if fw["J"] in G["J"] else abs(fw["numeric"]), "k*", ms=13)
ax.annotate(f"Fermi–Walker check at r = {fw['r0']}, J = {fw['J']}", (fw["r0"], abs(fw["numeric"])), (2.0, 1.0),
            arrowprops=dict(arrowstyle="->"), fontsize=9)
ax.set_ylim(0, 2.2); ax.set_xlabel("radius r (M = 1, ℓ = 1)"); ax.set_ylabel("turn rate (per unit proper time)")
ax.set_title("Spinning BTZ black hole: a gyroscope held still turns\nat w = J / (2(r² − M))  (checks G2, G3)")
ax.legend(fontsize=8.5)
ax = axs[1]
ratio = np.logspace(-1.5, 1.5, 300)
ax.semilogx(ratio, (ratio - 1) / (ratio + 1), color=BLUE, lw=2.2)
ax.plot(0.25, -0.6, "*", color=RED, ms=16)
ax.annotate("left movers 4× hotter\n(the spinning bath tested\nin the entanglement thread):\ndrag = 3/5 of light speed", (0.25, -0.6), (0.6, -0.85),
            arrowprops=dict(arrowstyle="->"), fontsize=9)
ax.axhline(0, color=GRAY, lw=0.5); ax.axvline(1, color=GRAY, lw=0.5)
ax.set_xlabel("T_R / T_L (ratio of the two temperatures)"); ax.set_ylabel("horizon drag speed Ω_H ℓ")
ax.set_title("The temperature ratio is the lift: T_R/T_L − 1 composes by ∘₁,\nso stacking drags multiplies the ratios (checks G6–G9)")
save(fig, "fig8_btz_frame_dragging.png")

# ---------------------------------------------------------------- fig 9: map
fig, ax = plt.subplots(figsize=(13, 5.8))
ax.set_xlim(-1.35, 1.35); ax.set_ylim(-1.3, 1.25); ax.axis("off")
ax.annotate("", xy=(1.25, 0), xytext=(-1.25, 0), arrowprops=dict(arrowstyle="<|-|>", lw=2.5, color="k"))
for x, t in ((-1, "d = −1\nboosts"), (0, "d = 0\nGalilean"), (1, "d = +1\nturns")):
    ax.plot(x, 0, "o", color="k", ms=10); ax.text(x, -0.22, t, ha="center", va="top", fontsize=11)
left = ["Thomas precession: 2π(γ − 1) per lap",
        "Gravity Probe B geodetic: 6603.9 mas/yr",
        "   (measured 6601.8 ± 18.3)",
        "Double pulsar B: 5.07°/yr (4.77 ± 0.66)",
        "[twisted-boosts page, earlier today]"]
right = ["Foucault: 360° × sin(latitude) per day",
         "Free tumbling body: Montgomery phase",
         "Spin in a field: Berry phase = ½ twist",
         "Gyroscope under torque: the 90° is the",
         "   corner's cross part; rate m g l/(I₃ω₃)",
         "[this page]"]
for i, t in enumerate(left):
    ax.text(-1.3, 1.1 - 0.17 * i, t, fontsize=10.5, color=ORANGE if i < 4 else GRAY)
for i, t in enumerate(right):
    ax.text(0.28, 1.1 - 0.17 * i, t, fontsize=10.5, color=BLUE if i < 5 else GRAY)
ax.text(-0.35, 0.35, "no corner,\nno precession", fontsize=10.5, color="k", ha="center")
ax.text(-1.3, -0.75, "Frame dragging (a spinning mass adds its own twist):", fontsize=11, color=PURPLE)
ax.text(-1.3, -0.92, "Gravity Probe B 40.8 mas/yr (measured 37.2 ± 7.2); spinning BTZ: a gyroscope held still turns at J/(2(r² − M));", fontsize=10.5, color=PURPLE)
ax.text(-1.3, -1.07, "its horizon drag speed is the boost whose Doppler factor squared is T_R/T_L. Computed with GR; the theory supplies the bookkeeping, not the field.", fontsize=10.5, color=PURPLE)
ax.set_title("Where each precession sits: one rule, turn per loop = d × enclosed area, read off the plane dial", fontsize=13)
save(fig, "fig9_map.png")
