"""DRAFT figures for qsharp-functors/README.md. Orange = known, blue = ours (a reading), red = does not work."""
import itertools, math
import numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Circle, FancyArrowPatch, Rectangle

OR, BL, RD, GR, INK = "#e07b00", "#1f5fbf", "#c62828", "#8a8f99", "#0b1f4b"
plt.rcParams.update({"font.size": 12, "axes.titlesize": 13, "axes.titleweight": "bold", "figure.dpi": 110})
def save(fig, name):
    fig.savefig(f"figures/{name}.png", bbox_inches="tight", facecolor="white"); plt.close(fig)
def arrow(ax, a, b, c=INK, lw=1.8, cs="arc3,rad=0", ms=14):
    ax.add_patch(FancyArrowPatch(a, b, arrowstyle="-|>", mutation_scale=ms, color=c, lw=lw, connectionstyle=cs))
def qmul(a, b):
    w1, x1, y1, z1 = a; w2, x2, y2, z2 = b
    return np.array([w1*w2-x1*x2-y1*y2-z1*z2, w1*x2+x1*w2+y1*z2-z1*y2, w1*y2-x1*z2+y1*w2+z1*x2, w1*z2+x1*y2-y1*x2+z1*w2])

# ---------------- F1: the gate sphere in our coordinates
fig, ax = plt.subplots(figsize=(9, 6.2))
for k, c in ((2, "#9ec1f0"), (1, BL), (0.5, "#4f86d6"), (0.25, "#7aa6e6")):
    ax.add_patch(Circle((-1/k, 0), 1/k, fill=False, color=c, lw=2.6 if k == 1 else 1.6))
    lab = {2: (-0.5, -0.75), 1: (-1.0, -1.3), 0.5: (-2.0, -2.3), 0.25: (-7.2, 2.7)}[k]
    ax.text(*lab, f"k = {k}", ha="center", color=c, fontsize=11, weight="bold")
ax.axvline(0, color=OR, lw=3)
ax.text(0.12, 3.3, "k → 0: the flat plane Re x = 0\n= Hamiltonians (generators).\nHere ∘₀ is plain adding.", color=OR, fontsize=11, va="top")
ax.plot(0, 0, "o", color=INK, ms=9); ax.text(0.1, -0.35, "identity (x = 0)\non every sphere", fontsize=10)
ax.plot(-2, 0, "o", color=INK, ms=7); ax.text(-2.1, 0.15, "−1: the 2π turn (k = 1)", ha="right", fontsize=10)
ax.plot(-1, 1, "s", color=BL, ms=8); ax.text(-0.95, 1.12, "X  (x = −1 + i at k = 1)", color=BL, fontsize=10)
ax.plot(-1+np.cos(np.pi/4)-0, np.sin(np.pi/4), "s", color=BL, ms=8)
ax.text(0.15, 0.55, "S (a quarter turn\non the Bloch sphere)", color=BL, fontsize=10)
ax.set_xlim(-8.6, 3.6); ax.set_ylim(-3.6, 3.6); ax.set_aspect("equal")
ax.set_xlabel("real part of x"); ax.set_ylabel("one imaginary direction of x (i, j or k)")
ax.set_title("Gates are the sphere |1 + kx| = 1: centre −1/k, radius 1/|k|, curvature d = k²")
ax.text(-8.4, -3.3, "Gates multiply on the round spheres (a ∘ₖ b = gate product).\nShrinking k flattens the sphere onto the generator plane.", fontsize=10.5, color=INK)
save(fig, "fig1-gates-are-the-sphere-1-plus-kx")

# ---------------- F2: Clifford corners
def cliff():
    H = np.array([0, 1, 0, 1])/np.sqrt(2); S = np.array([1, 0, 0, 1])/np.sqrt(2)
    G = {tuple(np.round([1, 0, 0, 0], 9))}; fr = [np.array([1., 0, 0, 0])]
    while fr:
        nx = []
        for a in fr:
            for g in (H, S):
                b = qmul(a, g); c = tuple(np.round(b, 9) + 0.)
                if c not in G: G.add(c); nx.append(b)
        fr = nx
    return np.array(sorted(G))
pts = cliff()
def kind(p):
    a = np.sort(np.abs(p))[::-1]
    return 0 if np.allclose(a, [1, 0, 0, 0]) else (1 if np.allclose(a, [.5]*4) else 2)
# a symmetric projection of R^4 to the plane (Coxeter-like)
th = np.array([0, 1, 2, 3])*np.pi/4 + 0.13
Pm = np.array([np.cos(th), np.sin(th)])
def proj(p): return Pm @ p
fig, axs = plt.subplots(1, 3, figsize=(16, 5.6))
titles = ["Paulis: the 8 corners of the\ncross-polytope (±1, ±i, ±j, ±k)",
          "Axis-cycling Cliffords: the 16\ncorners of the tesseract (±1±i±j±k)/2",
          "H, S, √X…: the 24 corners of the\ndual 24-cell (two of ±1,±i,±j,±k)/√2"]
labels = {(0, 1, 0, 0): "X", (0, 0, 1, 0): "Y", (0, 0, 0, 1): "Z", (1, 0, 0, 0): "1", (-1, 0, 0, 0): "−1"}
for kk, ax in enumerate(axs):
    sub = [p for p in pts if kind(p) == kk]
    other = [p for p in pts if kind(p) != kk]
    for p in other:
        q = proj(p); ax.plot(*q, "o", color="#d5d9e2", ms=4)
    dmin = min(np.linalg.norm(a-b) for a, b in itertools.combinations(sub, 2))
    for a, b in itertools.combinations(sub, 2):
        if abs(np.linalg.norm(a-b) - dmin) < 1e-6:
            pa, pb = proj(a), proj(b); ax.plot([pa[0], pb[0]], [pa[1], pb[1]], color=OR, lw=0.8, alpha=0.6)
    for p in sub:
        q = proj(p); ax.plot(*q, "o", color=OR, ms=8)
        key = tuple(int(round(v)) for v in p) if kk == 0 else None
        if key in labels: ax.text(q[0]+0.05, q[1]+0.05, labels[key], fontsize=12, color=INK, weight="bold")
    if kk == 2:
        for nm, p in (("H", np.array([0, 1, 0, 1])/np.sqrt(2)), ("S", np.array([1, 0, 0, 1])/np.sqrt(2))):
            q = proj(p); ax.text(q[0]+0.05, q[1]+0.05, nm, fontsize=12, weight="bold")
        t = np.array([np.cos(np.pi/8), 0, 0, np.sin(np.pi/8)]); q = proj(t)
        ax.plot(*q, "*", color=RD, ms=18); ax.text(0, -1.25, "red star T: on no corner (0.39 away in 4D;\nit only looks close in this view)", color=RD, fontsize=11, ha="center")
    ax.add_patch(Circle((0, 0), max(np.linalg.norm(proj(p)) for p in pts)+0.02, fill=False, color=GR, ls=":"))
    ax.set_aspect("equal"); ax.axis("off"); ax.set_title(titles[kk])
fig.suptitle("The 48 single-qubit Clifford gates on S³ are the corners of our shapes (one fixed 4D → 2D view; grey = the other two sets)", fontsize=13)
save(fig, "fig2-clifford-gates-are-corners-of-our-shapes")

# ---------------- F3: Adjoint is the mirror
fig, axs = plt.subplots(1, 3, figsize=(16.5, 5.6))
ax = axs[0]
ax.add_patch(Circle((-1, 0), 1, fill=False, color=BL, lw=2.4))
ax.axhline(0, color=GR, lw=1, ls="--"); ax.text(-2.75, 0.08, "real axis\n(the mirror)", fontsize=10, color=GR)
ang = 1.0; u = (-1+np.cos(ang), np.sin(ang)); ub = (u[0], -u[1])
ax.plot(*u, "o", color=BL, ms=10); ax.text(u[0]+0.06, u[1]+0.05, "U (x)", color=BL, fontsize=12, weight="bold")
ax.plot(*ub, "o", color=BL, ms=10); ax.text(ub[0]+0.06, ub[1]-0.15, "Adjoint U = x̄\n= ∘ₖ inverse", color=BL, fontsize=11, weight="bold")
arrow(ax, (u[0], u[1]-0.05), (ub[0], ub[1]+0.05), BL, cs="arc3,rad=-0.3")
anti = (-2-u[0], -u[1]); ax.plot(*anti, "X", color=RD, ms=11)
ax.text(-2.75, -1.35, "antipode −U: the same gate up to phase −1,\nNOT the adjoint", color=RD, fontsize=10, ha="left")
ax.plot(0, 0, "o", color=INK); ax.text(0.05, 0.05, "1", fontsize=11); ax.plot(-2, 0, "o", color=INK); ax.text(-2.2, -0.2, "−1", fontsize=11)
ax.set_xlim(-2.8, 0.9); ax.set_ylim(-1.6, 1.4); ax.set_aspect("equal"); ax.axis("off")
ax.set_title("A. Adjoint is the mirror through the real axis\n(k = 1 shown; the same for every k)")
ax = axs[1]
steps = ["A", "B", "C"]
for i, s in enumerate(steps):
    ax.add_patch(Rectangle((i*1.2, 1.0), 0.9, 0.7, color=BL, alpha=0.18)); ax.text(i*1.2+0.45, 1.35, s, ha="center", va="center", fontsize=16, weight="bold")
    j = 2 - i
    ax.add_patch(Rectangle((i*1.2, -0.6), 0.9, 0.7, color=BL, alpha=0.18)); ax.text(i*1.2+0.45, -0.25, steps[j]+"†", ha="center", va="center", fontsize=16, weight="bold")
arrow(ax, (-0.1, 2.0), (3.4, 2.0), INK); ax.text(1.65, 2.15, "body: A then B then C", ha="center", fontsize=11)
arrow(ax, (-0.1, -0.85), (3.4, -0.85), INK); ax.text(1.65, -1.2, "Adjoint: C† then B† then A†", ha="center", fontsize=11)
ax.text(1.65, 0.45, "conj(x ∘ₖ y) = ȳ ∘ₖ x̄", ha="center", fontsize=14, color=BL, weight="bold")
ax.set_xlim(-0.4, 3.8); ax.set_ylim(-1.5, 2.5); ax.axis("off")
ax.set_title("B. The mirror reverses order:\nthis is Q#'s 'adjoint invert'")
ax = axs[2]
t = np.linspace(0, 2*np.pi, 200)
ax.plot(np.cos(t), np.sin(t), color=GR, lw=1)
ax.plot(np.cos(t), 0.32*np.sin(t), color=OR, lw=3)
ax.text(0, -0.55, "the sphere of i's (Re = 0): π-turns\nhere the mirror and the antipode agree", ha="center", color=OR, fontsize=10.5)
for nm, a in (("X", 0.3), ("Y", 1.6), ("Z", 2.7), ("H", 4.4)):
    ax.plot(np.cos(a), 0.32*np.sin(a), "o", color=OR, ms=9); ax.text(np.cos(a)+0.04, 0.32*np.sin(a)+0.07, nm, fontsize=12, weight="bold")
ax.plot(0, 1, "o", color=INK); ax.text(0.06, 1.02, "1 (fixed)", fontsize=10)
ax.plot(0, -1, "o", color=INK); ax.text(0.06, -1.1, "−1 (fixed)", fontsize=10)
for nm, (xx, yy) in (("S", (0.5, 0.75)), ("T", (0.25, 0.93))):
    ax.plot(xx, yy, "s", color=BL, ms=8); ax.text(xx+0.05, yy, nm+" (not self-adjoint)", color=BL, fontsize=10)
ax.set_aspect("equal"); ax.axis("off"); ax.set_xlim(-1.3, 1.6)
ax.set_title("C. 'adjoint self' gates (X, Y, Z, H)\nlive exactly on the sphere of i's")
save(fig, "fig3-adjoint-is-the-mirror")

# ---------------- F4: Controlled is a cube corner
fig, axs = plt.subplots(1, 3, figsize=(16.5, 5.6))
ax = axs[0]
V = np.array(list(itertools.product((0, 1), repeat=3)))
def p3(v): return np.array([v[0] + 0.45*v[2], v[1] + 0.32*v[2]])
for a, b in itertools.combinations(V, 2):
    if np.sum(np.abs(a-b)) == 1:
        pa, pb = p3(a), p3(b); ax.plot([pa[0], pb[0]], [pa[1], pb[1]], color=GR, lw=1.4)
for v in V:
    q = p3(v); top = tuple(v) == (1, 1, 1)
    ax.plot(*q, "o", color=BL if top else "#c7ccd6", ms=20 if top else 12)
    ax.text(q[0]+0.06, q[1]-0.13, "".join(map(str, v)), fontsize=10, color=INK)
q = p3((1, 1, 1)); ax.text(q[0]-0.1, q[1]+0.12, "U acts only here", color=BL, fontsize=12, weight="bold", ha="center")
ax.text(0.7, -0.45, "three controls: Controlled = 1 + P1⊗P1⊗P1⊗(U − 1)\nevery other corner: identity", ha="center", fontsize=10.5)
ax.set_aspect("equal"); ax.axis("off"); ax.set_xlim(-0.4, 1.9); ax.set_ylim(-0.7, 1.75)
ax.set_title("A. Controlled lights one corner\nof the control cube")
ax = axs[1]
ax.add_patch(Rectangle((0, 1), 1, 1, color="#e9edf5")); ax.text(0.5, 1.5, "1", ha="center", va="center", fontsize=22)
ax.add_patch(Rectangle((1, 0), 1, 1, color=BL, alpha=0.25)); ax.text(1.5, 0.5, "U", ha="center", va="center", fontsize=22, weight="bold", color=BL)
ax.text(0.5, 0.5, "0", ha="center", va="center", fontsize=18, color=GR); ax.text(1.5, 1.5, "0", ha="center", va="center", fontsize=18, color=GR)
ax.plot([0, 2, 2, 0, 0], [0, 0, 2, 2, 0], color=INK, lw=1.5)
ax.text(1, -0.35, "in our coordinates:  x  ↦  P1 ⊗ x   (for every k)", ha="center", fontsize=12, color=BL, weight="bold")
ax.text(1, -0.75, "keeps ∘ₖ because P1·P1 = P1:  P1x ∘ₖ P1y = P1(x ∘ₖ y)\nso Controlled can be pushed into every step ('controlled distribute')", ha="center", fontsize=10.5)
ax.set_xlim(-0.6, 2.6); ax.set_ylim(-1.3, 2.3); ax.set_aspect("equal"); ax.axis("off")
ax.set_title("B. The controlled gate: identity on the\ncontrol-0 block, U on the control-1 block")
ax = axs[2]
for (a, b) in itertools.product((0, 1), repeat=2):
    lit = a or b
    ax.plot(a, b, "o", color=OR if lit else "#c7ccd6", ms=24 if lit else 14)
    ax.text(a+0.08, b-0.14, f"{a}{b}", fontsize=11)
ax.plot([0, 1, 1, 0, 0], [0, 0, 1, 1, 0], color=GR)
ax.text(0.5, -0.45, "'control on a OR b':  P_a ∘₋₁ P_b = P_a + P_b − P_a P_b\nour union rule, applied to projectors", ha="center", fontsize=10.5, color=INK)
ax.set_xlim(-0.5, 1.5); ax.set_ylim(-0.8, 1.4); ax.set_aspect("equal"); ax.axis("off")
ax.set_title("C. OR of two controls is our union rule:\nthree corners lit, AND is one")
save(fig, "fig4-controlled-lights-one-corner")

# ---------------- F5: commuting square, and what breaks it
fig, axs = plt.subplots(1, 2, figsize=(15, 5.6))
ax = axs[0]
P = {"U": (0, 1), "U†": (1.6, 1), "C(U)": (0, 0), "C(U)† = C(U†)": (1.6, 0)}
for nm, (x, y) in P.items():
    ax.text(x, y, nm, ha="center", va="center", fontsize=15, weight="bold", bbox=dict(boxstyle="round,pad=0.4", fc="#eef3fb", ec=BL))
arrow(ax, (0.25, 1), (1.35, 1), BL); ax.text(0.8, 1.1, "mirror (Adjoint)", ha="center", fontsize=11, color=BL)
arrow(ax, (0.25, 0), (1.05, 0), BL); ax.text(0.65, 0.1, "mirror", ha="center", fontsize=11, color=BL)
arrow(ax, (0, 0.85), (0, 0.15), BL); ax.text(-0.07, 0.5, "corner\n(Controlled)", ha="right", fontsize=11, color=BL)
arrow(ax, (1.6, 0.85), (1.6, 0.15), BL); ax.text(1.67, 0.5, "corner", ha="left", fontsize=11, color=BL)
ax.text(0.8, -0.45, "(1 + P1⊗(U − 1))† = 1 + P1⊗(U† − 1)  because  P1† = P1", ha="center", fontsize=12)
ax.set_xlim(-0.7, 2.4); ax.set_ylim(-0.7, 1.35); ax.axis("off")
ax.set_title("A. Controlled Adjoint = Adjoint Controlled: the two moves commute")
ax = axs[1]
props = ["distributes over a body\n(needs e·e = e)", "commutes with Adjoint\n(needs e† = e)", "stays a gate (unitary)\n(needs both)"]
rows = [("P1 = |1⟩⟨1|\n(Q#'s control)", [1, 1, 1]), ("oblique E = [[0,0],[1,1]]\n(E·E = E, E† ≠ E)", [1, 0, 0]), ("½·P1\n(e·e ≠ e, e† = e)", [0, 1, 0])]
for j, p in enumerate(props): ax.text(j*1.9+1.3, 3.2, p, ha="center", fontsize=10.5)
for i, (nm, vals) in enumerate(rows):
    y = 2.3 - i*0.95; ax.text(-0.1, y, nm, ha="right", va="center", fontsize=11)
    for j, v in enumerate(vals):
        ax.add_patch(Rectangle((j*1.9+0.85, y-0.35), 0.9, 0.7, color=BL if v else RD, alpha=0.85))
        ax.text(j*1.9+1.3, y, "yes" if v else "no", ha="center", va="center", color="white", fontsize=13, weight="bold")
ax.set_xlim(-2.4, 5.6); ax.set_ylim(-0.3, 3.7); ax.axis("off")
ax.set_title("B. Which property of the corner does which job")
save(fig, "fig5-adjoint-and-controlled-commute")

# ---------------- F6: Hopf circle and kickback
fig, axs = plt.subplots(1, 3, figsize=(16.5, 5.6))
def sphere_outline(ax, c=GR):
    t = np.linspace(0, 2*np.pi, 200); ax.plot(np.cos(t), np.sin(t), color=c, lw=1.2); ax.plot(np.cos(t), 0.3*np.sin(t), color=c, lw=0.8, ls=":")
ax = axs[0]; sphere_outline(ax)
pt = np.array([0.45, 0.55]); ax.plot(*pt, "o", color=INK, ms=9)
t = np.linspace(0, 2*np.pi, 100); ax.plot(pt[0] + 0.22*np.cos(t), pt[1] + 0.55 + 0.08*np.sin(t), color=OR, lw=2.5)
ax.plot([pt[0], pt[0]], [pt[1], pt[1]+0.47], color=OR, lw=1, ls="--")
ax.text(pt[0]+0.27, pt[1]+0.6, "Hopf circle: the state's phase\n(invisible on its own)", color=OR, fontsize=10.5)
ax.text(0, -1.3, "target's Bloch sphere; an eigenstate of U only\nslides round its circle by the eigenvalue angle φ", ha="center", fontsize=10.5)
ax.set_aspect("equal"); ax.axis("off"); ax.set_xlim(-1.3, 1.9); ax.set_ylim(-1.6, 1.5)
ax.set_title("A. Target in an eigenstate: it moves\nonly along its Hopf circle")
ax = axs[1]; sphere_outline(ax)
ax.plot(1, 0, "o", color=GR, ms=8); ax.text(1.04, -0.12, "|+⟩", fontsize=11)
phi = 1.1; e = (np.cos(phi), 0.3*np.sin(phi)); ax.plot(*e, "o", color=OR, ms=11)
ax.text(e[0]+0.05, e[1]+0.08, "turned by φ", color=OR, fontsize=11, weight="bold")
arrow(ax, (1, 0.02), (e[0]+0.03, e[1]), OR, cs="arc3,rad=0.25")
ax.text(0, -1.3, "the control turns about Z by the same φ:\nControlled(global phase) = phase gate on the control", ha="center", fontsize=10.5)
ax.set_aspect("equal"); ax.axis("off"); ax.set_xlim(-1.3, 1.6); ax.set_ylim(-1.6, 1.5)
ax.set_title("B. ...and the control turns by that angle\n(phase kickback, standard)")
ax = axs[2]
rng = np.random.default_rng(3); lens = []
SX = np.array([[0, 1], [1, 0]], complex); SY = np.array([[0, -1j], [1j, 0]]); SZ = np.diag([1, -1]).astype(complex)
for _ in range(4000):
    q = rng.normal(size=4); q /= np.linalg.norm(q); w, x, y, z = q
    U = w*np.eye(2) - 1j*(x*SX + y*SY + z*SZ)
    tt = rng.normal(size=2) + 1j*rng.normal(size=2); tt /= np.linalg.norm(tt)
    st = np.vstack([tt, U@tt])/np.sqrt(2); rho = st@st.conj().T
    lens.append(np.linalg.norm(np.real([np.trace(rho@m) for m in (SX, SY, SZ)])))
ax.hist(lens, bins=40, color=RD, alpha=0.8)
ax.axvline(1, color=OR, lw=2.5); ax.text(0.97, ax.get_ylim()[1]*1.02, "eigenstate (length 1)", color=OR, ha="right", fontsize=10.5)
ax.set_xlabel("length of the control's Bloch vector after Controlled-U"); ax.set_ylabel("random targets and U")
ax.set_title("C. Typical target: NOT one fiber. The target's\npoint moves and the control gets entangled")
save(fig, "fig6-controlled-and-the-hopf-circle")

# ---------------- F7: quaternionic Hopf picture for two qubits
fig, ax = plt.subplots(figsize=(9.5, 7))
sphere_outline(ax)
for h in (-0.6, 0, 0.6):
    r = math.sqrt(1-h*h); t = np.linspace(0, 2*np.pi, 200); ax.plot(r*np.cos(t), h + 0.25*r*np.sin(t), color=BL, lw=1.2 if h else 2)
ax.plot([0, 0], [-1.15, 1.15], color=INK, lw=1, ls="--")
ax.plot(0, 1, "o", color=INK, ms=9); ax.text(0.07, 1.06, "control surely 0: Controlled does nothing", fontsize=11)
ax.plot(0, -1, "o", color=OR, ms=11); ax.text(0.08, -1.15, "control surely 1: U acts inside one S³ fiber\n(the voice-call hypothesis holds here)", fontsize=11, color=OR)
arrow(ax, (0.75, 0.12), (-0.3, 0.24), BL, cs="arc3,rad=0.35", lw=2.4, ms=18)
ax.text(0.9, 0.3, "in between: Controlled-U turns S⁴\nabout the control axis by left-multiplying\nthe 'relative frame' of the two branches by U;\nlatitude (control probabilities) is kept", fontsize=10.5, color=BL)
ax.text(-1.25, 1.35, "height = P(control 0) − P(control 1)", fontsize=11)
ax.set_aspect("equal"); ax.axis("off"); ax.set_xlim(-1.4, 2.8); ax.set_ylim(-1.5, 1.5)
ax.set_title("Two qubits on S⁷, seen on the quaternionic Hopf base S⁴ (fibers are S³)")
save(fig, "fig7-two-qubits-on-the-quaternionic-hopf-base")

# ---------------- F8: swap rule
def rot(n, th):
    n = np.asarray(n, float)/np.linalg.norm(n); return math.cos(th/2)*np.eye(2) - 1j*math.sin(th/2)*(n[0]*SX + n[1]*SY + n[2]*SZ)
fig, axs = plt.subplots(1, 2, figsize=(14, 5.4))
ax = axs[0]
ax.add_patch(Rectangle((0, 0), 1, 0.7, fill=True, color=OR, alpha=0.2))
for (a, b, txt) in (((0, 0), (1, 0), "Rx(a)"), ((1, 0), (1, 0.7), "Ry(b)"), ((1, 0.7), (0, 0.7), "Rx(−a)"), ((0, 0.7), (0, 0.04), "Ry(−b)")):
    arrow(ax, a, b, INK, lw=2.2); ax.text((a[0]+b[0])/2 + (0.07 if a[0] == b[0] else 0), (a[1]+b[1])/2 + (0.05 if a[1] == b[1] else 0), txt, fontsize=12)
ax.text(0.5, 0.35, "area a·b", ha="center", fontsize=14, color=OR, weight="bold")
ax.text(0.5, -0.25, "what is left: a turn about Z by a·b\n(turn = dial × area; the dial is 1 on the Bloch sphere)", ha="center", fontsize=11)
ax.set_xlim(-0.3, 1.4); ax.set_ylim(-0.5, 0.95); ax.set_aspect("equal"); ax.axis("off")
ax.set_title("A. Swap two small gates and undo them")
ax = axs[1]
sizes = np.geomspace(0.003, 0.6, 40); r = []
for s in sizes:
    Q = rot([1, 0, 0], s)@rot([0, 1, 0], s)@rot([1, 0, 0], -s)@rot([0, 1, 0], -s)
    ang = 2*math.acos(np.clip((np.trace(Q)/2).real, -1, 1)); r.append(ang/(s*s))
ax.semilogx(sizes, r, color=OR, lw=2.5); ax.axhline(1, color=INK, ls="--", lw=1)
ax.set_xlabel("size of each turn (radians)"); ax.set_ylabel("leftover turn ÷ area")
ax.set_title("B. The ratio goes to exactly 1 for small turns\n(known: the Lie bracket; used in Solovay–Kitaev)")
save(fig, "fig8-swap-rule-for-gates")

# ---------------- F9: sorting board
rows = [("Single-qubit gates = points of S³ (unit quaternions)", "known", OR, "proved"),
        ("Gate product = our ∘ₖ on the sphere |1+kx| = 1, any k ≠ 0", "ours (reading)", BL, "proved"),
        ("Dial k = the gate sphere's curvature; k → 0 gives Hamiltonians that add", "known math, our reading", BL, "proved"),
        ("Clifford gates = corners of cross-polytope, tesseract, dual 24-cell", "known", OR, "proved"),
        ("Adjoint = our ∘ₖ inverse = mirror x → x̄", "ours (reading)", BL, "proved"),
        ("Adjoint = antipode (the 'opposite point')", "—", RD, "no"),
        ("Adjoint reverses order ('adjoint invert')", "standard", OR, "proved"),
        ("'adjoint self' gates = the sphere of i's", "standard, our picture", BL, "proved"),
        ("Controlled = x ↦ P1⊗x, a corner of the control cube", "ours (reading)", BL, "proved"),
        ("OR of controls = our union rule ∘₋₁", "ours (reading)", BL, "proved"),
        ("Controlled acts on one Hopf fiber", "only at the control-1 pole", RD, "no (as stated)"),
        ("Controlled exposes the phase circle (kickback)", "standard", OR, "proved"),
        ("Controlled Adjoint = Adjoint Controlled, from P1·P1 = P1, P1† = P1", "standard, our split", BL, "proved"),
        ("Gate commutator = swap rule (turn = area)", "known", OR, "proved")]
fig, ax = plt.subplots(figsize=(14, 0.55*len(rows) + 1.2))
for i, (txt, src, c, verdict) in enumerate(rows):
    y = len(rows) - i
    ax.add_patch(Rectangle((0, y-0.42), 0.18, 0.84, color=c))
    ax.text(0.3, y, txt, va="center", fontsize=12)
    ax.text(9.6, y, src, va="center", fontsize=11, color=c)
    ax.text(12.6, y, verdict, va="center", fontsize=12, weight="bold", color=RD if verdict.startswith("no") else INK)
ax.text(9.6, len(rows)+0.9, "where it comes from", fontsize=11, weight="bold"); ax.text(12.6, len(rows)+0.9, "status", fontsize=11, weight="bold")
ax.set_xlim(0, 14.5); ax.set_ylim(0.3, len(rows)+1.4); ax.axis("off")
ax.set_title("Q#'s three forms on our shapes: orange = known, blue = our reading, red = does not hold", loc="left")
save(fig, "fig9-sorting-board")
print("ok")
