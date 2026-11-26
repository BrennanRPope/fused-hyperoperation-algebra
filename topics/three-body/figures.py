"""DRAFT figures for three-body/README.md. Run after checks.py: python3 figures.py  (numpy, scipy, matplotlib)."""
import sys, os, json
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from tbpy import *
from chore import Chore
import checks_lib as C

OUT = os.path.join(HERE, "figures"); os.makedirs(OUT, exist_ok=True)
KNOWN, OURS, GREY, RED, GREEN = "#d9822b", "#2b6cb0", "#888888", "#c53030", "#2f855a"
BODY = ["#2b6cb0", "#c53030", "#222222"]
plt.rcParams.update({"font.size": 11, "axes.titlesize": 12})
RES = json.load(open(os.path.join(HERE, "results.json")))
SP = special_points()

def save(fig, name):
    fig.savefig(os.path.join(OUT, name), dpi=130, bbox_inches="tight"); plt.close(fig); print("wrote", name)

def sphere(ax, alpha=0.10, marks=True, view=(22, -60)):
    u, v = np.mgrid[0:2 * np.pi:36j, 0:np.pi:18j]
    ax.plot_wireframe(np.cos(u) * np.sin(v), np.sin(u) * np.sin(v), np.cos(v), color=GREY, alpha=alpha, lw=0.5)
    th = np.linspace(0, 2 * np.pi, 200); ax.plot(np.cos(th), np.sin(th), 0 * th, color=GREY, lw=0.8)
    if marks:
        for k, p in SP.items():
            if k.startswith("collision"): ax.scatter(*p, color=RED, marker="x", s=50, depthshade=False)
            elif k.startswith("Lagrange"): ax.scatter(*p, color=KNOWN, s=40, depthshade=False)
            else: ax.scatter(*p, color=GREEN, s=25, depthshade=False)
    ax.set_box_aspect((1, 1, 1)); ax.set_axis_off(); ax.view_init(*view)

def orbit(name, n=3000):
    a, b, T = C.SUV[name]; q, r = refine(a, b, T); return traj(suvakov(*q[:2]), q[2], n)

def stereo(u):  # S^3 -> R^3 seen from the point (0,0,0,1)
    return u[:, :3] / (1 - u[:, 3])[:, None]

# ------------------------------------------------------------ fig1 two spheres
fig = plt.figure(figsize=(16, 5.2))
ax = fig.add_subplot(1, 3, 1)
r = np.array([[-1, 0], [1, 0], [0.35, 1.1]])
ax.add_patch(plt.Polygon(r, fill=False, color=GREY, ls=":"))
for k in range(3): ax.plot(*r[k], "o", color=BODY[k], ms=12); ax.text(r[k, 0] + 0.07, r[k, 1] + 0.07, f"body {k+1}")
ax.annotate("", xy=r[1], xytext=r[0], arrowprops=dict(arrowstyle="-|>", color=OURS, lw=2)); ax.text(-0.3, -0.18, "z₁ (1 to 2)", color=OURS)
ax.annotate("", xy=r[2], xytext=(0, 0), arrowprops=dict(arrowstyle="-|>", color=OURS, lw=2)); ax.text(0.25, 0.45, "z₂ (midpoint to 3)", color=OURS)
ax.text(-1.25, -0.65, "Two complex numbers (z₁, z₂) = a point of C² = R⁴.\nDivide by the size: a point of the 3-sphere S³.\nForget the overall turn: a point of the shape sphere S².", fontsize=10.5)
ax.set_aspect("equal"); ax.set_xlim(-1.4, 1.6); ax.set_ylim(-1.0, 1.5); ax.axis("off")
ax.set_title("A. A triangle of three bodies is a point on S³")
ax = fig.add_subplot(1, 3, 2, projection="3d"); sphere(ax)
for k, p in SP.items():
    lab = {"collision 12": "1–2 collide", "Euler 3 mid": "Euler: 3 in middle", "Lagrange +": "Lagrange ▲"}.get(k)
    if lab: ax.text(*(1.18 * p), lab, fontsize=9)
ax.set_title("B. The shape sphere (Hopf's map S³ → S²)\npoles = Lagrange triangles, equator = in a line\nred × = collisions, green = Euler lines")
ax = fig.add_subplot(1, 3, 3, projection="3d")
u, v = np.mgrid[0:2 * np.pi:36j, 0:np.pi:18j]
ax.plot_wireframe(np.cos(u) * np.sin(v), np.sin(u) * np.sin(v), np.cos(v), color=GREY, alpha=0.1, lw=0.5)
X0, Tk, p0 = C.homographic(0.6); t, Xs = traj(X0, Tk, 600)
F = C.fock(Xs[:, 6:8], p0); ax.plot(F[:, 0], F[:, 1], F[:, 3], color=KNOWN, lw=3)
ax.set_box_aspect((1, 1, 1)); ax.set_axis_off(); ax.view_init(20, -50)
ax.set_title("C. A different S³: Fock/Moser's momentum sphere.\nA Kepler orbit's momenta map to one great circle\n(planar orbits use its 2-sphere slice)")
save(fig, "fig1-two-three-spheres.png")

# ------------------------------------------------------------ fig2 Lagrange and Euler keep both spheres
fig = plt.figure(figsize=(16, 9))
for row, (name, euler) in enumerate((("Lagrange", False), ("Euler", True))):
    X0, Tk, p0 = C.homographic(0.5, euler); t, Xs = traj(X0, Tk, 800)
    ax = fig.add_subplot(2, 4, 4 * row + 1)
    for k in range(3): ax.plot(Xs[:, 2 * k], Xs[:, 2 * k + 1], color=BODY[k]); ax.plot(*Xs[0, 2 * k:2 * k + 2], "o", color=BODY[k])
    for i in (0, 200, 400):
        ax.add_patch(plt.Polygon(Xs[i, :6].reshape(3, 2), fill=False, color=GREY, ls=":"))
    ax.set_aspect("equal"); ax.axis("off"); ax.set_title(f"{name} (eccentricity 0.5): each body a Kepler ellipse,\nthe triangle only grows, shrinks and turns")
    ax = fig.add_subplot(2, 4, 4 * row + 2)
    for k in range(3): ax.plot(Xs[:, 6 + 2 * k], Xs[:, 7 + 2 * k], color=BODY[k])
    ax.plot(0, 0, "k+"); ax.set_aspect("equal"); ax.axis("off"); ax.set_title("velocities: exact circles\n→ great circles on Fock's sphere")
    ax = fig.add_subplot(2, 4, 4 * row + 3, projection="3d"); sphere(ax)
    w = shape(Xs); ax.scatter(*w[0], color=OURS, s=120, depthshade=False)
    ax.set_title("shape sphere: one fixed point\n(" + ("a pole" if not euler else "an equator point") + ")")
    ax = fig.add_subplot(2, 4, 4 * row + 4, projection="3d")
    P = stereo(s3point(Xs)); ax.plot(*P.T, color=OURS, lw=3)
    ax.set_box_aspect((1, 1, 1)); ax.set_axis_off()
    ax.set_title("configuration S³ (stereographic view):\none Hopf circle = a great circle")
save(fig, "fig2-lagrange-euler-keep-both.png")

# ------------------------------------------------------------ fig3 hodographs
fig = plt.figure(figsize=(16, 9.5))
items = [("Lagrange", None), ("figure-eight", "figure-eight"), ("butterfly I", "butterfly I"), ("moth I", "moth I")]
for col, (lab, name) in enumerate(items):
    if name is None:
        X0, Tk, p0 = C.homographic(0.5); t, Xs = traj(X0, Tk, 800)
    else:
        t, Xs = orbit(name, 40000)
    V = Xs[:, 6:8]; c, rho, res = C.circfit(V)
    ax = fig.add_subplot(2, 4, col + 1)
    ax.plot(V[:, 0], V[:, 1], color=KNOWN if name is None else (OURS if name == "figure-eight" else RED), lw=2)
    th = np.linspace(0, 2 * np.pi, 300); ax.plot(c[0] + rho * np.cos(th), c[1] + rho * np.sin(th), color=GREY, ls="--", lw=1)
    ax.plot(0, 0, "k+", ms=10); ax.set_aspect("equal"); ax.axis("off")
    ax.set_title(f"{lab}: body 1's velocities\nbest circle (dashed) misses by {100*res:.1f}%" if res < 1 else f"{lab}: body 1's velocities\nbest circle misses by {res:.1f}× its radius")
    ax = fig.add_subplot(2, 4, col + 5, projection="3d")
    u_, v_ = np.mgrid[0:2 * np.pi:36j, 0:np.pi:18j]
    ax.plot_wireframe(np.cos(u_) * np.sin(v_), np.sin(u_) * np.sin(v_), np.cos(v_), color=GREY, alpha=0.08, lw=0.5)
    p0 = np.sqrt(max(rho ** 2 - c @ c, 1e-3)) if res < 1 else 1.0
    F = C.fock(V, p0); ax.plot(F[:, 0], F[:, 1], F[:, 3], color=KNOWN if name is None else (OURS if name == "figure-eight" else RED), lw=2)
    if name == "figure-eight":
        Gc = C.fock(np.column_stack([c[0] + rho * np.cos(th), c[1] + rho * np.sin(th)]), p0); ax.plot(Gc[:, 0], Gc[:, 1], Gc[:, 3], color=GREY, ls="--", lw=1)
    ax.set_box_aspect((1, 1, 1)); ax.set_axis_off(); ax.view_init(20, -50)
    ax.set_title({0: "Fock's sphere: a full great circle", 1: "Fock's sphere: rocks along ~220° of\nnearly one great circle", 2: "Fock's sphere: no great circle", 3: "Fock's sphere: no great circle"}[col])
save(fig, "fig3-momenta-on-focks-sphere.png")

# ------------------------------------------------------------ fig4 shape sphere gallery
names = ["figure-eight", "butterfly I", "moth I", "yin-yang Ia", "goggles", "dragonfly"]
fig = plt.figure(figsize=(16, 10))
for k, name in enumerate(names):
    t, Xs = orbit(name); w = shape(Xs)
    ax = fig.add_subplot(2, 3, k + 1, projection="3d"); sphere(ax)
    ax.plot(*w.T, color=OURS, lw=1.6)
    ax.set_title(f"{name} on the shape sphere")
fig.suptitle("Every zero-angular-momentum orbit is a closed loop on the shape sphere. Its type is how it winds round the three collisions (red ×).", y=0.95)
save(fig, "fig4-shape-sphere-gallery.png")

# ------------------------------------------------------------ fig5 configuration S^3: fibre vs figure-eight
fig = plt.figure(figsize=(14, 6))
ax = fig.add_subplot(1, 2, 1, projection="3d")
X0, Tk, p0 = C.homographic(0.5); t, Xs = traj(X0, Tk, 600); P = stereo(s3point(Xs)); ax.plot(*P.T, color=KNOWN, lw=3, label="Lagrange: one Hopf circle")
X0, Tk, p0 = C.homographic(0.5, True); t, Xs = traj(X0, Tk, 600); P = stereo(s3point(Xs)); ax.plot(*P.T, color=GREEN, lw=3, label="Euler: another Hopf circle")
ax.legend(loc="upper left", fontsize=9); ax.set_box_aspect((1, 1, 1)); ax.set_axis_off()
ax.set_title("A. Homographic orbits ride Hopf circles of S³\n(great circles; the shape never changes)")
ax = fig.add_subplot(1, 2, 2, projection="3d")
t, Xs = orbit("figure-eight"); P = stereo(s3point(Xs)); ax.plot(*P.T, color=OURS, lw=2.5)
ax.set_box_aspect((1, 1, 1)); ax.set_axis_off()
ax.set_title("B. The figure-eight on S³: a closed curve that\nmoves across Hopf circles, never along one")
save(fig, "fig5-configuration-sphere.png")

# ------------------------------------------------------------ fig6 half-area law
rows = np.array(RES["B"]["rows"])
fig, axs = plt.subplots(1, 2, figsize=(14, 5.6))
ax = axs[0]
ax.plot(rows[:, 1], rows[:, 0], "o", ms=4, color=OURS, alpha=0.7)
lim = np.pi; ax.plot([-lim, lim], [-lim, lim], color=KNOWN, lw=1)
ax.set_xlabel("−½ × area enclosed on the shape sphere (radians)"); ax.set_ylabel("how far the triangle actually turned (radians)")
ax.set_title(f"A. {len(rows)} closed shape loops of one chaotic orbit:\nturn = −½ × enclosed area, every time (the swap rule)")
ax.set_aspect("equal")
ax = fig.add_subplot(1, 2, 2, projection="3d"); axs[1].remove(); sphere(ax)
X0 = suvakov(0.40, 0.30); t, Xs = traj(X0, 80.0, 160001, 1e-12); W = shape(Xs)
i0, i1 = int(0.333 / 80 * 160000), int(4.674 / 80 * 160000)
ax.plot(*W[i0:i1 + 1].T, color=OURS, lw=2)
ax.set_title("B. One such loop (t = 0.33 to 4.67): it encloses\nsolid angle −1.15, and the triangle turns by +0.58")
save(fig, "fig6-swap-rule-on-the-shape-sphere.png")

# ------------------------------------------------------------ fig7 dial = minus energy
fig, axs = plt.subplots(1, 3, figsize=(15, 5))
for ax, E, lab, col in zip(axs, (-0.5, 0.0, 0.5), ("E < 0: bound ellipse\nFock's space is a sphere (d > 0)", "E = 0: parabola\nflat space (d = 0)", "E > 0: hyperbola\nhyperbolic space (d < 0)"), (KNOWN, GREY, OURS)):
    h = np.sqrt(2 * (E + 1)); th = np.linspace(0, 2 * np.pi, 400)
    ax.plot(-np.sin(th) / h, (np.cos(th) + h * h - 1) / h, color=GREY, ls=":", lw=1)
    thr = np.linspace(0, 2 * np.pi, 400) if E < 0 else np.linspace(-np.arccos(-1 / (h * h - 1)) + 0.02 if E > 0 else -np.pi + 0.05, (np.arccos(-1 / (h * h - 1)) - 0.02) if E > 0 else np.pi - 0.05, 400)
    ax.plot(-np.sin(thr) / h, (np.cos(thr) + h * h - 1) / h, color=col, lw=3)
    ax.plot(0, 0, "k+", ms=12); ax.set_aspect("equal"); ax.axis("off"); ax.set_title(lab)
    ax.text(0.02, 0.02, "radius² − centre² = −2E", transform=ax.transAxes)
fig.suptitle("Velocity circles of Kepler orbits. Where the origin sits (inside, on, outside) is the sign of the energy, and that sign is the plane dial of Fock's space.", y=1.02)
save(fig, "fig7-dial-is-minus-energy.png")

# ------------------------------------------------------------ fig8 action minimizer
D = RES["D"]; As = np.array([d["A"] for d in D])
fig = plt.figure(figsize=(16, 5.4))
ax = fig.add_subplot(1, 3, 1)
ax.hist(As, bins=60, color=OURS); ax.set_xlabel("action at the stopping point"); ax.set_ylabel("random starts (of 40)")
ax.set_title("A. 40 random closed curves, each slid downhill\nin action: they pile up at two values")
ch = Chore(N=30, M=900)
best = {}
for d in D:
    k = round(d["A"], 3)
    if np.isfinite(d["ret"]) and d["ret"] < 1e-2: best.setdefault(k, d)
for ax_i, (k, d) in enumerate(sorted(best.items())[:2]):
    ax = fig.add_subplot(1, 3, 2 + ax_i)
    q = ch.curve(np.array(d["x"]))
    ax.plot(q[:, 0], q[:, 1], color=KNOWN, lw=2.5)
    for s_ in range(3):
        p = q[int(s_ * len(q) / 3)]; ax.plot(*p, "o", color=BODY[s_], ms=10)
    ax.set_aspect("equal"); ax.axis("off")
    ax.set_title(("B. Lowest: Lagrange's circle" if ax_i == 0 else "C. Next: the figure-eight") + f"\naction {k:.3f}, {sum(abs(x['A']-k)<1e-3 for x in D)} of 40 starts")
save(fig, "fig8-action-minimums.png")

# ------------------------------------------------------------ fig9 landscape
L = np.load(os.path.join(HERE, "landscape.npz")); E = RES["E"]
fig, axs = plt.subplots(1, 2, figsize=(17, 7.6), layout="constrained")
for ax, key, title in ((axs[0], "dabs", "A. How close each start comes back (dark = closer)"), (axs[1], "tabs", "B. When it comes closest (time)")):
    F = L[key]
    im = ax.imshow(np.log10(F) if key == "dabs" else F, origin="lower", extent=(0, 0.8, 0, 0.8), cmap="magma" if key == "dabs" else "viridis",
                   vmin=-2 if key == "dabs" else 0, vmax=0.3 if key == "dabs" else 40)
    plt.colorbar(im, ax=ax, shrink=0.85, label="log₁₀ closest return" if key == "dabs" else "time of closest return")
    for f in E:
        if f["kind"] == "absolute" or abs(f["alpha"]) < 1e-6:
            ax.plot(f["p1"], f["p2"], "o", mfc="none", mec=KNOWN if f["match"] else "#00e5ff", mew=1.6, ms=7 if f["match"] else 11)
    for name in ("figure-eight", "butterfly I", "moth I", "yin-yang Ia", "goggles", "dragonfly", "bumblebee"):
        a, b, _ = C.SUV[name]; ax.annotate(name, (a, b), xytext=(5, 5), textcoords="offset points", color="white", fontsize=8.5)
    ax.set_xlabel("p₁ (sideways push)"); ax.set_ylabel("p₂ (upward push)"); ax.set_title(title)
axs[0].text(0.01, 0.79, "orange rings: orbits whose loop type is in Li–Liao's 2017 list\ncyan rings: periodic orbits whose loop type is not in it", color="white", fontsize=9, va="top")
fig.suptitle("The minimum search. Each pixel starts the three bodies in a line (−1, 0, +1) with pushes (p₁, p₂), (p₁, p₂), (−2p₁, −2p₂); 449,773 starts, each run to time 40.")
save(fig, "fig9-return-landscape.png")

# ------------------------------------------------------------ fig10 the unlisted orbits
new = [f for f in E if not f["match"]]
seen = set(); show = []
for f in sorted(new, key=lambda f: f["T"]):
    c_ = canon(f["word"])
    if c_ in seen: continue
    seen.add(c_); show.append(f)
fig = plt.figure(figsize=(16, 4.2 * len(show) / 2 + 1))
for k, f in enumerate(show):
    t, Xs = traj(suvakov(f["p1"], f["p2"]), f["T"], int(f["T"] * 600))
    ax = fig.add_subplot(len(show) // 2 + len(show) % 2, 4, 2 * k + 1)
    for b in range(3): ax.plot(Xs[:, 2 * b], Xs[:, 2 * b + 1], color=BODY[b], lw=0.8)
    ax.set_aspect("equal"); ax.axis("off")
    ax.set_title(f"T = {f['T']:.3f}, closest pass {f['dmin']:.4f}\np = ({f['p1']:.6f}, {f['p2']:.6f})", fontsize=9.5)
    ax = fig.add_subplot(len(show) // 2 + len(show) % 2, 4, 2 * k + 2, projection="3d"); sphere(ax, marks=True)
    ax.plot(*shape(Xs).T, color="#0891b2", lw=0.9); ax.set_title(f"loop type, {len(f['word'])} letters", fontsize=9.5)
fig.suptitle("The 8 loop types the search found that are not in Li–Liao's 695-family list (each a real periodic orbit; real space left, shape sphere right)", y=1.0)
save(fig, "fig10-unlisted-orbits.png")

# ------------------------------------------------------------ fig11 sorting board
fig, ax = plt.subplots(figsize=(15, 6.4)); ax.axis("off")
rows = [("Question", "Answer", "Known or ours"),
        ("Lagrange and Euler on the 3-sphere", "Keep it exactly: one Hopf circle in shape space,\none great circle per body on Fock's sphere", "known (homographic = Kepler);\npictures ours"),
        ("Figure-eight on Fock's sphere", "Keeps a piece: momenta rock along ~220° of\nnearly one great circle (7.8% off)", "observation ours; not seen stated"),
        ("Šuvakov orbits on Fock's sphere", "Break it: best circles miss by 3 to 8 radii", "ours (check)"),
        ("Turn of the triangle (zero spin)", "−½ × area on the shape sphere: the swap rule,\nsame law as a spin-½ Berry phase", "known (Montgomery 1996);\nswap-rule reading ours"),
        ("Energy sign", "is the plane dial of Fock's space\n(sphere, flat, hyperbolic)", "known (Moser, Osipov, Belbruno, Milnor)"),
        ("Action minimum", "finds Lagrange's circle and the figure-eight, nothing else", "known (Moore 1993, Chenciner–Montgomery 2000)"),
        ("Return-distance minimum", "66 periodic orbits; 8 loop types not in\nLi–Liao's 2017 list, all with very close passes", "method known (Šuvakov 2013);\nthe 8 unlisted types ours, pending later catalogs")]
cols = [0.0, 0.34, 0.72]
for r_, row in enumerate(rows):
    y = 1 - r_ * 0.125
    for c_, txt in enumerate(row):
        col = "black"
        if r_ and c_ == 2: col = KNOWN if txt.startswith("known") or txt.startswith("method") else OURS
        ax.text(cols[c_], y, txt, fontsize=11 if r_ else 12, fontweight="bold" if r_ == 0 else None, va="top", color=col)
save(fig, "fig11-sorting-board.png")
