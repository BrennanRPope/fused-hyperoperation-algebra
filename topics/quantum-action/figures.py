"""Figures for the quantum action page. Run after checks.py:  python3 figures.py  -> figures/*.png"""
import json, math, os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon, FancyArrowPatch, Circle
from mpl_toolkits.mplot3d import Axes3D  # noqa
from mpl_toolkits.mplot3d.art3d import Poly3DCollection

HERE = os.path.dirname(os.path.abspath(__file__))
FIG = os.path.join(HERE, "figures"); os.makedirs(FIG, exist_ok=True)
R = json.load(open(os.path.join(HERE, "results.json")))
BLUE, ORANGE, GREEN, PURPLE, GRAY, RED = "#2a6fdb", "#e8833a", "#2e9e5b", "#8a4fcf", "#888888", "#cc3d3d"
plt.rcParams.update({"font.size": 11, "axes.spines.top": False, "axes.spines.right": False})

def save(fig, name):
    fig.savefig(os.path.join(FIG, name), dpi=140, bbox_inches="tight"); plt.close(fig)

# ---------------------------------------------------------------- fig1: action is area
fig, ax = plt.subplots(1, 3, figsize=(15, 4.8))
a = ax[0]
loop = np.array([(0, 0), (0.5, -0.3), (1, 0.4), (0.2, 0.9), (-0.4, 0.3)])
a.add_patch(Polygon(loop, closed=True, fc=BLUE, alpha=0.25, ec=BLUE, lw=2))
for p_, q_ in zip(loop, np.roll(loop, -1, axis=0)):
    a.annotate("", q_, p_, arrowprops=dict(arrowstyle="->", color=BLUE, lw=2))
for i in range(-6, 12):
    for k in range(-6, 12):
        a.add_patch(plt.Rectangle((i*0.1, k*0.1), 0.1, 0.1, fill=False, ec="#dddddd", lw=0.5))
a.text(0.3, 0.3, "area = action\n= 0.87", ha="center", fontsize=12, color=BLUE)
a.set_xlim(-0.6, 1.2); a.set_ylim(-0.5, 1.1); a.set_aspect("equal")
a.set_xlabel("position q"); a.set_ylabel("momentum p")
a.set_title("A path in the plane of position and momentum.\nAction is the area it encloses.")
a = ax[1]
u, v = np.array([0.7, 0.2]), np.array([-0.3, 0.9])
par = np.array([[0, 0], u, u+v, v])
a.add_patch(Polygon(par, closed=True, fc=ORANGE, alpha=0.25, ec=ORANGE, lw=1.5))
a.annotate("", u, (0, 0), arrowprops=dict(arrowstyle="->", color=BLUE, lw=2.5)); a.text(*(u/2 + [0.03, -0.08]), "move a", color=BLUE)
a.annotate("", u+v, u, arrowprops=dict(arrowstyle="->", color=PURPLE, lw=2.5)); a.text(*(u + v/2 + [0.06, -0.12]), "then b", color=PURPLE)
a.annotate("", v, (0, 0), arrowprops=dict(arrowstyle="->", color=PURPLE, lw=2.5, ls="--")); a.text(*(v/2 + [-0.32, 0]), "b first", color=PURPLE)
a.annotate("", u+v, v, arrowprops=dict(arrowstyle="->", color=BLUE, lw=2.5, ls="--")); a.text(*(v + u/2 + [-0.1, 0.06]), "then a", color=BLUE)
a.text(0.05, 0.55, "swap cost =\narea / ħ", color=ORANGE, fontsize=12)
a.set_xlim(-0.5, 0.9); a.set_ylim(-0.1, 1.25); a.set_aspect("equal")
a.set_xlabel("position q"); a.set_ylabel("momentum p")
a.set_title("The swap rule in phase space: doing two moves\nin the other order turns the state by area/ħ")
a = ax[2]
A1 = R["A1"]
xs = [r["area_over_hbar"] for r in A1]; ys = [r["phase"] for r in A1]
a.plot([-1.3, 2], [-1.3, 2], color=GRAY, lw=1)
a.scatter(xs, ys, s=70, color=BLUE, zorder=3, label="measured turn of the wavefunction")
a.set_xlabel("enclosed area / ħ"); a.set_ylabel("quantum phase (radians)")
a.set_title("Measured on a real wavefunction (check A1):\nfive loops, error 0")
a.legend(loc="upper left", frameon=False)
save(fig, "fig1-action-is-area.png")

# ---------------------------------------------------------------- fig2: sphere -> plane, phase = j x turn
fig = plt.figure(figsize=(17, 4.8))
a = fig.add_subplot(1, 3, 1, projection="3d")
uu, vv = np.mgrid[0:2*np.pi:40j, 0:np.pi:20j]
a.plot_wireframe(np.cos(uu)*np.sin(vv), np.sin(uu)*np.sin(vv), np.cos(vv), color="#cccccc", lw=0.5)
t = np.linspace(0, 2*np.pi, 200); th0 = 0.6
a.plot(np.sin(th0)*np.cos(t), np.sin(th0)*np.sin(t), np.cos(th0)*np.ones_like(t), color=BLUE, lw=3)
cap_v = np.linspace(0, th0, 10)
for cv in cap_v:
    a.plot(np.sin(cv)*np.cos(t), np.sin(cv)*np.sin(t), np.cos(cv)*np.ones_like(t), color=BLUE, alpha=0.15)
a.set_box_aspect((1, 1, 1)); a.set_axis_off()
a.set_title("Spin j: the loop's cap has solid angle Ω.\nA frame turns by Ω (precession);\nthe quantum state turns by j·Ω (Berry)")
a = fig.add_subplot(1, 3, 2)
A3 = R["A3"]
js = [r["j"] for r in A3]
a.plot(js, [-r["phase"] for r in A3], "o-", color=BLUE, lw=2, label="quantum phase (stays)")
a.plot(js, [r["turn"] for r in A3], "s-", color=ORANGE, lw=2, label="frame turn = d × area (fades)")
a.axhline(A3[0]["flat_area_over_hbar"], color=GRAY, ls="--"); a.text(20, A3[0]["flat_area_over_hbar"]-0.09, "flat phase space: area/ħ", color=GRAY)
a.set_xticks(js); a.set_xlabel("spin j  (sphere radius² = j·ħ, so d = 1/(jħ))"); a.set_ylabel("radians")
a.set_title("The same small square on bigger and bigger spheres\n(check A3)")
a.legend(frameon=False, loc="center right")
a = fig.add_subplot(1, 3, 3); a.set_axis_off()
a.text(0.0, 0.92, "One swap rule, three readings", fontsize=14, weight="bold")
lines = [("turn of a frame", "= d × area", "Thomas, Foucault, Berry (precession page)", ORANGE),
         ("turn of a quantum state", "= j × (turn of a frame)", "spin j on its sphere (A3a)", BLUE),
         ("flat phase space", "= area / ħ  with  ħ = 1/(j·d)", "the sphere flattens, ħ is what's left (A3b, A4)", PURPLE)]
for i, (h, f_, s_, c) in enumerate(lines):
    y = 0.72 - 0.27*i
    a.text(0.0, y, h, fontsize=12, color=c, weight="bold"); a.text(0.0, y-0.08, f_, fontsize=12); a.text(0.0, y-0.15, s_, fontsize=10, color=GRAY)
save(fig, "fig2-one-swap-rule-three-readings.png")

# ---------------------------------------------------------------- fig3: one dial: Fermat, Newton, Einstein
fig, ax = plt.subplots(1, 4, figsize=(18, 4.6))
T = 2.0; tt = np.linspace(0, T, 400)
titles = {0.5: "d > 0: Fermat (optics)\nlength of the path, least", 0.0: "d = 0: Newton\n∫ ½mv² dt, least",
          -0.5: "d < 0: Einstein\n−mc²·(proper time), clock time most"}
for a, d in zip(ax[:3], [0.5, 0.0, -0.5]):
    for e, c in [(0, BLUE), (0.25, ORANGE), (-0.25, ORANGE)]:
        a.plot(0.3*tt + e*np.sin(np.pi*tt/T), tt, color=c, lw=2.5 if e == 0 else 1.5)
    # unit "circles" of the dial geometry centred at the start event
    th = np.linspace(-np.pi, np.pi, 400)
    if d > 0:
        a.plot(0.5*np.sin(th)/math.sqrt(d), 0.5*np.cos(th), color=GRAY, ls=":")
    elif d < 0:
        sh = np.linspace(-1.2, 1.2, 200); a.plot(0.5*np.sinh(sh)/math.sqrt(-d), 0.5*np.cosh(sh), color=GRAY, ls=":")
    else:
        xx = np.linspace(-0.8, 0.8, 100); a.plot(xx, 0.5 + 0*xx, color=GRAY, ls=":")
    a.set_xlim(-0.6, 1.0); a.set_ylim(-0.05, 2.05); a.set_xlabel("place x"); a.set_ylabel("time t")
    a.set_title(titles[d]); a.text(-0.55, 1.85, "blue: straight path\norange: wiggles\ndotted: one unit of\nthis plane's length", fontsize=9, color=GRAY)
a = ax[3]
eps = np.array([-0.1, -0.05, 0, 0.05, 0.1])
for d, c in [(0.5, GREEN), (0.0, BLUE), (-0.5, PURPLE)]:
    dS = R["B1"][str(d)]["dS"]; vals = [dS[0], dS[1], 0, dS[2], dS[3]]
    a.plot(eps, vals, "o-", color=c, label=f"d = {d:+.1f}")
a.set_xlabel("wiggle size"); a.set_ylabel("extra action"); a.legend(frameon=False)
a.set_title("S_d = (m/d) ∫ √(1 + d v²) dt\nleast on the straight path at every d (B1a)")
save(fig, "fig3-one-dial-fermat-newton-einstein.png")

# ---------------------------------------------------------------- fig4: Kepler + fish-eye
fig = plt.figure(figsize=(16, 5))
a = fig.add_subplot(1, 3, 1)
cols = [BLUE, GREEN, ORANGE, PURPLE]
for e, c in zip([0.0, 0.3, 0.6, 0.9], cols):
    ang = np.linspace(0, 2*np.pi, 300); rr = (1-e*e)/(1+e*np.cos(ang))
    a.plot(rr*np.cos(ang), rr*np.sin(ang), color=c, lw=2, label=f"e = {e}")
a.plot(0, 0, "*", color="k", ms=12); a.set_aspect("equal"); a.set_ylim(-1.1, 1.6); a.legend(frameon=False, fontsize=9, loc="upper left", ncol=2)
a.set_title("Four orbits with the same energy")
a = fig.add_subplot(1, 3, 2)
B2 = R["B2"]
labels = [f"e={r['e']}" for r in B2]
a.bar(labels, [r["angular_part"] for r in B2], color=BLUE, label="turning part 2πL")
a.bar(labels, [r["radial_part"] for r in B2], bottom=[r["angular_part"] for r in B2], color=ORANGE, label="in-and-out part")
a.axhline(2*np.pi, color="k", lw=1); a.text(-0.4, 2*np.pi+0.15, "total action 2π, same for all (B2)")
a.set_ylim(0, 8.4); a.legend(frameon=False, loc="upper right", ncol=2, fontsize=9); a.set_ylabel("action ∮ p·dq")
a.set_title("Every one has the same action.\nThat is why hydrogen's levels are shared.")
a = fig.add_subplot(1, 3, 3, projection="3d")
uu, vv = np.mgrid[0:2*np.pi:40j, 0:np.pi:20j]
a.plot_wireframe(np.cos(uu)*np.sin(vv), np.sin(uu)*np.sin(vv), np.cos(vv), color="#dddddd", lw=0.4)
for (X, Y), c in zip(R["B3_rays"], cols):
    X, Y = np.array(X), np.array(Y); rr = X**2 + Y**2
    a.plot(2*X/(1+rr), 2*Y/(1+rr), (1-rr)/(1+rr), color=c, lw=2.5)
    m = rr < 9
    a.plot(X[m]/3, Y[m]/3, -1.3 + 0*X[m], color=c, lw=1.2, alpha=0.8)
a.text(0, 0, -1.6, "flat lens: rays are circles", ha="center", fontsize=9)
a.set_axis_off(); a.set_box_aspect((1, 1, 1.1))
a.set_title("Maxwell's fish-eye lens (below): its light rays\nlift to great circles of a sphere (B3)")
save(fig, "fig4-kepler-and-fisheye.png")

# ---------------------------------------------------------------- fig5: action shapes: cube and simplex
rng = np.random.default_rng(3)
fig = plt.figure(figsize=(15, 5.2))
a = fig.add_subplot(1, 3, 1, projection="3d")
v = rng.normal(size=(1500, 3, 3)); v /= np.linalg.norm(v, axis=2, keepdims=True)
z = v[:, :, 2]
a.scatter(z[:, 0], z[:, 1], z[:, 2], s=2, color=BLUE, alpha=0.4)
g = np.array([-1, -1/3, 1/3, 1])
G = np.array(np.meshgrid(g, g, g)).reshape(3, -1).T
a.scatter(G[:, 0], G[:, 1], G[:, 2], s=40, color=ORANGE, depthshade=False)
for s_ in [-1, 1]:
    for t_ in [-1, 1]:
        a.plot([-1, 1], [s_, s_], [t_, t_], color="k", lw=0.8); a.plot([s_, s_], [-1, 1], [t_, t_], color="k", lw=0.8); a.plot([s_, s_], [t_, t_], [-1, 1], color="k", lw=0.8)
a.set_title("Three spins: their actions (heights) fill\nthe CUBE evenly. Orange: the 4×4×4 quantum\nstates of three spin-3/2's (C1)")
a.set_xlabel("J₁"); a.set_ylabel("J₂"); a.set_zlabel("J₃")
a = fig.add_subplot(1, 3, 2)
psi = rng.normal(size=(3000, 3)) + 1j*rng.normal(size=(3000, 3)); psi /= np.linalg.norm(psi, axis=1, keepdims=True)
w = abs(psi)**2
P2 = w[:, 1] + 0.5*w[:, 2]; Q2 = (math.sqrt(3)/2)*w[:, 2]
a.scatter(P2, Q2, s=2, color=BLUE, alpha=0.4)
a.plot([0, 1, 0.5, 0], [0, 0, math.sqrt(3)/2, 0], color="k")
n = 6
for i in range(n+1):
    for k in range(n+1-i):
        l = n - i - k; ww = np.array([i, k, l])/n
        a.plot(ww[1] + 0.5*ww[2], math.sqrt(3)/2*ww[2], "o", color=ORANGE, ms=7)
a.text(-0.06, -0.06, "|ψ₁|²=1"); a.text(0.93, -0.06, "|ψ₂|²=1"); a.text(0.42, 0.9, "|ψ₃|²=1")
a.set_aspect("equal"); a.set_axis_off()
a.set_title("Random quantum states of a 3-level system:\nBorn weights fill the SIMPLEX evenly (C2a).\nOrange: the 28 states of 3 oscillators with 6 quanta", pad=18)
a = fig.add_subplot(1, 3, 3); a.set_axis_off()
a.text(0, 0.92, "Shapes of action", fontsize=14, weight="bold")
rows = [("independent spins", "cube", "each J between −j and +j"),
        ("oscillators sharing energy", "simplex", "quanta add up to n"),
        ("one quantum system's weights", "simplex", "probabilities add up to 1"),
        ("count quantum states", "lattice points", "spacing ħ, shifted by ½ (C2c)")]
for i, (s1, s2, s3) in enumerate(rows):
    y = 0.75 - 0.2*i
    a.text(0, y, s1, fontsize=11); a.text(0.78, y, s2, fontsize=12, color=ORANGE, weight="bold"); a.text(0, y-0.07, s3, fontsize=9, color=GRAY)
save(fig, "fig5-shapes-of-action.png")

# ---------------------------------------------------------------- fig6: Bohr-Sommerfeld with the half
C3 = R["C3"]
fig, ax = plt.subplots(1, 2, figsize=(13, 4.8))
a = ax[0]
for n, c in zip(range(5), [BLUE, GREEN, ORANGE, PURPLE, RED]):
    E = C3["bs_half"][n]; xt = E**0.25; X = np.linspace(-xt, xt, 300); Pp = np.sqrt(np.maximum(E - X**4, 0))
    a.fill_between(X, -Pp, Pp, color=c, alpha=0.08); a.plot(np.r_[X, X[::-1]], np.r_[Pp, -Pp[::-1]], color=c, lw=1.8, label=f"n = {n}: area = 2πħ({n}+½)")
    a.plot([xt, -xt], [0, 0], "|", color=c, ms=14, mew=2)
a.set_aspect("equal"); a.set_xlabel("position"); a.set_ylabel("momentum"); a.legend(frameon=False, fontsize=8, loc="upper right")
a.set_title("Allowed loops in a quartic well: area comes in steps of 2πħ,\nplus a half. Bars: the two places the loop turns vertical")
a = ax[1]
nn = np.arange(12)
eh = np.abs(np.array(C3["bs_half"]) - C3["exact"])/np.array(C3["exact"])
ez = np.abs(np.array(C3["bs_zero"], dtype=float) - C3["exact"])/np.array(C3["exact"])
a.semilogy(nn, eh, "o-", color=BLUE, label="with the half (n + ½)")
a.semilogy(nn[1:], ez[1:], "s-", color=ORANGE, label="without it (n)")
a.set_xlabel("level n"); a.set_ylabel("relative error vs exact level"); a.legend(frameon=False)
a.set_title("The half quantum is the loop's own turning:\ntwo quarter turns per lap (Maslov), checks C3a–c")
save(fig, "fig6-bohr-sommerfeld-and-the-half.png")

# ---------------------------------------------------------------- fig7: the path dial, with gravity placed on it
fig, ax = plt.subplots(figsize=(8.5, 7.5))
a = ax
a.axhline(0, color="k", lw=0.8); a.axvline(0, color="k", lw=0.8)
a.add_patch(Circle((0, 0), 3.2, fill=False, ls="--", color=GRAY, lw=1.5))
a.plot(0, 0, "o", color="k"); a.text(0.12, -0.35, "k = 0: count the paths\n(every path weighs 1)", fontsize=10)
a.annotate("", (0, 2.6), (0, 0.3), arrowprops=dict(arrowstyle="->", color=BLUE, lw=2.5))
a.text(0.12, 1.5, "imaginary k = i/ħ: QUANTUM\nsum over all paths of e^{iS/ħ}\nFeynman; Ponzano–Regge 3D gravity", color=BLUE, fontsize=10)
a.annotate("", (-2.6, 0), (-0.3, 0), arrowprops=dict(arrowstyle="->", color=ORANGE, lw=2.5))
a.text(-2.9, 0.25, "negative real k: STATISTICS\nsum of e^{−S/T}\nCDT and lattice gravity live here", color=ORANGE, fontsize=10)
a.annotate("", (-1.2, 1.2), (-0.15, 2.0), arrowprops=dict(arrowstyle="->", color=PURPLE, lw=1.8, connectionstyle="arc3,rad=0.4"))
a.text(-2.1, 1.75, "Wick rotation:\na quarter turn of the dial", color=PURPLE, fontsize=9)
a.text(1.0, 3.25, "ring at |k| = ∞: CLASSICAL\nonly the least / stationary action survives\n(Newton; Einstein's equations)", color=GREEN, fontsize=10)
a.set_xlim(-3.6, 3.6); a.set_ylim(-1.2, 4.1); a.set_aspect("equal"); a.set_axis_off()
a.set_title("Quantum vs least action is one dial (from the quantum-to-classical page),\nwith quantum gravity's attempts placed on it")
save(fig, "fig7-path-sum-vs-least-action-dial.png")

# ---------------------------------------------------------------- fig8: Regge
def icosphere(nsub):
    t = (1+5**0.5)/2
    V = [(-1,t,0),(1,t,0),(-1,-t,0),(1,-t,0),(0,-1,t),(0,1,t),(0,-1,-t),(0,1,-t),(t,0,-1),(t,0,1),(-t,0,-1),(-t,0,1)]
    V = [np.array(v)/np.linalg.norm(v) for v in V]
    F = [(0,11,5),(0,5,1),(0,1,7),(0,7,10),(0,10,11),(1,5,9),(5,11,4),(11,10,2),(10,7,6),(7,1,8),(3,9,4),(3,4,2),(3,2,6),(3,6,8),(3,8,9),(4,9,5),(2,4,11),(6,2,10),(8,6,7),(9,8,1)]
    for _ in range(nsub):
        cache = {}; NF = []
        def mid(i, j):
            key = (min(i, j), max(i, j))
            if key not in cache:
                m = V[i]+V[j]; V.append(m/np.linalg.norm(m)); cache[key] = len(V)-1
            return cache[key]
        for a_, b_, c_ in F:
            ab, bc, ca = mid(a_, b_), mid(b_, c_), mid(c_, a_)
            NF += [(a_, ab, ca), (b_, bc, ab), (c_, ca, bc), (ab, bc, ca)]
        F = NF
    return np.array(V), F
fig = plt.figure(figsize=(16, 5))
a = fig.add_subplot(1, 3, 1, projection="3d")
V, F = icosphere(1)
ang = np.zeros(len(V))
for tri in F:
    P = V[list(tri)]
    for i in range(3):
        u_ = P[(i+1) % 3]-P[i]; w_ = P[(i+2) % 3]-P[i]
        ang[tri[i]] += math.acos(u_@w_/np.linalg.norm(u_)/np.linalg.norm(w_))
defi = 2*np.pi - ang
pc = Poly3DCollection([V[list(t)] for t in F], fc="#e8eef9", ec="#7a8fb8", lw=0.6, alpha=0.9); a.add_collection3d(pc)
sc = a.scatter(V[:, 0], V[:, 1], V[:, 2], c=defi, cmap="Oranges", s=60, depthshade=False)
a.set_box_aspect((1, 1, 1)); a.set_axis_off(); a.set_xlim(-1, 1); a.set_ylim(-1, 1); a.set_zlim(-1, 1)
a.set_title("Regge: space built from flat triangles.\nAll the bending sits at the corners (deficit angles)")
a = fig.add_subplot(1, 3, 2)
E1 = R["E1"]
a.plot([r["nsub"] for r in E1 if r["R"] == 1.0], [r["ratio"] for r in E1 if r["R"] == 1.0], "o-", color=BLUE, label="R = 1")
a.plot([r["nsub"] for r in E1 if r["R"] == 2.0], [r["ratio"] for r in E1 if r["R"] == 2.0], "x--", color=ORANGE, label="R = 2", ms=10)
a.axhline(1, color=GRAY, ls=":"); a.set_xlabel("finer triangulation (subdivisions)"); a.set_ylabel("deficit ÷ (area/R²)")
a.set_title("Each corner's deficit = d × its area, d = 1/R²:\nthe swap rule's turn law (E1b). Sum is always 4π (E1a)")
a.legend(frameon=False)
a = fig.add_subplot(1, 3, 3); a.set_axis_off()
a.text(0, 0.92, "Einstein's action, Regge's way", fontsize=14, weight="bold")
rows = [("S = (1/8πG) Σ (hinge size) × (deficit turn)", ""),
        ("2D: hinges are points", "S = Σ turns = 4π, never moves (E1c)"),
        ("3D: hinges are edges", "length × turn"),
        ("4D: hinges are triangles", "area × turn"),
        ("hinge dimension = D − 2", "= G's units, the dial on the asymptotic-safety page")]
for i, (s1, s2) in enumerate(rows):
    y = 0.75 - 0.16*i
    a.text(0, y, s1, fontsize=11, weight="bold" if i == 0 else None); a.text(0.02, y-0.065, s2, fontsize=10, color=GRAY)
save(fig, "fig8-regge-turns-times-sizes.png")

# ---------------------------------------------------------------- fig9: Ponzano-Regge
sc_ = np.array(R["E2_scan"])
fig, ax = plt.subplots(1, 2, figsize=(14, 4.8), gridspec_kw={"width_ratios": [2.2, 1]})
a = ax[0]
a.plot(sc_[:, 0], sc_[:, 1], "o", color="k", ms=6, label="exact 6j symbol (quantum recoupling of spins)")
a.plot(sc_[:, 0], sc_[:, 2], "-", color=BLUE, lw=2, label="cos(Regge action + π/4) / √(12πV)")
a.plot(sc_[:, 0], sc_[:, 3], "--", color=ORANGE, lw=1.5, label="same without the quarter turn")
a.set_xlabel("sixth spin j₆ (others 30, 31, 32, 31, 30)"); a.set_ylabel("value"); a.legend(frameon=False, fontsize=9)
a.set_title("3D quantum gravity: summing over spins gives back Einstein's (Regge's) action,\nwith a Maslov quarter turn (E2a, E2b)")
a = ax[1]
Vt = np.array([[0, 0, 0], [1, 0, 0], [0.45, 0.85, 0], [0.5, 0.3, 0.8]])
for i in range(4):
    for k in range(i+1, 4):
        a.plot([Vt[i, 0], Vt[k, 0]], [Vt[i, 1]+0.35*Vt[i, 2], Vt[k, 1]+0.35*Vt[k, 2]], color=BLUE, lw=2)
a.text(0.4, -0.15, "edge length = spin + ½  (in Planck lengths)", ha="center", fontsize=9)
a.text(0.4, 1.25, "S = Σ length × (π − dihedral angle)", ha="center", fontsize=10, color=BLUE)
a.set_xlim(-0.3, 1.2); a.set_ylim(-0.3, 1.4); a.set_axis_off(); a.set_aspect("equal")
save(fig, "fig9-ponzano-regge.png")

# ---------------------------------------------------------------- fig10: q turns on the cosmological constant
fig, ax = plt.subplots(1, 2, figsize=(15, 4.8))
for a, r in zip(ax, ["160", "400"]):
    rows = R["E3"][r]["rows"]
    j6 = [q["j6"] for q in rows]; ex = np.array([q["exact"] for q in rows])
    wv = np.array([q["with_volume"] for q in rows]); fl = np.array([q["flat_q1"] for q in rows])
    s1 = np.dot(ex, wv)/np.dot(wv, wv); s2 = np.nanmax(np.abs(ex))/np.nanmax(np.abs(fl))
    a.plot(j6, ex, "o", color="k", ms=6, label="exact quantum-group 6j")
    a.plot(j6, s1*wv, "-", color=GREEN, lw=2, label="curved (spherical) Regge action + volume term")
    a.plot(j6, s2*fl, "--", color=ORANGE, lw=1.3, label="flat formula (q = 1)")
    a.set_xlabel("sixth spin j₆"); a.legend(frameon=False, fontsize=9)
    a.set_title(f"q = e^(iπ/{r}): correlation {R['E3'][r]['corr_with_volume']:.4f}, signs {R['E3'][r]['sign_agree']:.0%}\n"
                f"(cosmological constant ∝ (2π/{r})²)")
fig.suptitle("Turning on the q dial curves every tetrahedron: 3D quantum gravity with a cosmological constant (check E3)", y=1.03)
save(fig, "fig10-q-is-the-cosmological-constant.png")

# ---------------------------------------------------------------- fig11: scorecard
items = [
    ("Action is phase-space area; ħ converts area to turn", "known", "Weyl 1927; Heisenberg group"),
    ("Same swap rule as precession; phase = j × frame turn; ħ = 1/(j·d)", "ours (framing)", "Berry 1984, Inönü–Wigner known"),
    ("One action S_d: Fermat (d>0), Newton (0), Einstein (d<0)", "framing ours", "kinematic groups: Bacry–Lévy-Leblond 1968"),
    ("Our k is only a change of velocity ruler: d' = d − k²/4", "ours", "B1d"),
    ("Kepler: all orbits of one energy share one action", "known", "Bohr–Sommerfeld; Fock 1935"),
    ("Cube = spins' actions, simplex = Born weights; states = lattice points", "known; link ours", "Atiyah, Guillemin–Sternberg, Delzant"),
    ("The ½ quantum = the loop turning vertical twice", "known", "Maslov 1965; Keller 1958"),
    ("Regge deficit = swap-rule turn = d × area", "known; link ours", "Regge 1961"),
    ("3D quantum gravity: spins sum to cos(Regge + π/4)", "known", "Ponzano–Regge 1968; Roberts 1999"),
    ("q dial = cosmological constant in 3D quantum gravity", "known; link ours", "Turaev–Viro 1992; Mizoguchi–Tada 1992"),
    ("Why 4D gravity's action fails at small scales", "nothing new", "G has units (D−2); see asymptotic safety"),
]
fig, a = plt.subplots(figsize=(14, 6.2)); a.set_axis_off()
a.text(0, 1.0, "Scorecard", fontsize=14, weight="bold")
for i, (t1, t2, t3) in enumerate(items):
    y = 0.92 - 0.083*i
    col = {"known": GRAY, "nothing new": RED}.get(t2, GREEN if "ours" in t2 else GRAY)
    a.text(0, y, t1, fontsize=11); a.text(0.64, y, t2, fontsize=11, color=col, weight="bold"); a.text(0.78, y, t3, fontsize=9, color=GRAY)
save(fig, "fig11-scorecard.png")
print("figures written:", sorted(os.listdir(FIG)))
