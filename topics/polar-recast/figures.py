"""Pictures for the polar recast (DRAFT). Colour key: orange = known math, blue = ours (bridge or reading)."""
import itertools, math, cmath
import numpy as np, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon, Circle
from scipy.linalg import expm
OR, BL, GR, RD, GY = "#e07b00", "#1f6fd1", "#2a9d55", "#c0392b", "#888888"
plt.rcParams.update({"font.size": 10, "axes.titlesize": 11})
def save(fig, name): fig.tight_layout(); fig.savefig(f"figures/{name}.png", dpi=130); plt.close(fig)
rng = np.random.default_rng(3)

# ---------- fig1: the right shape for each family
fig, axs = plt.subplots(3, 3, figsize=(11, 10.5))
def clean(ax, t, sub, col):
    ax.set_aspect("equal"); ax.axis("off"); ax.set_title(t, fontsize=11, fontweight="bold")
    ax.text(0.5, -0.06, sub, transform=ax.transAxes, ha="center", va="top", fontsize=9, color=col, wrap=True)
ax = axs[0, 0]  # cube
c = np.array(list(itertools.product([0, 1], repeat=3)), float); P = c[:, :2] + 0.45*c[:, 2:3]*np.array([1, 0.6])
for i, j in itertools.combinations(range(8), 2):
    if np.abs(c[i] - c[j]).sum() == 1: ax.plot(*zip(P[i], P[j]), color=BL, lw=2)
ax.scatter(*P.T, color=BL, zorder=3)
clean(ax, "1 Core law: CUBE", "1+(a∘b∘c) has one term per corner.\nThe cube really is right here:\nindependent axes, one switch each.", BL)
ax = axs[0, 1]  # nested rulers
for row, (f, lab) in enumerate([(lambda x: x, "x"), (np.log1p, "log(1+x)"), (lambda x: np.log1p(np.log1p(x)), "log(1+log(1+x))")]):
    y = -row
    ax.plot([0, 6], [y, y], color="k", lw=1)
    for v in ([100, 300, 600, 1000] if row == 0 else [1, 2, 5, 10, 20, 50, 100, 1000]):
        p = f(v)*6/f(1000); ax.plot([p, p], [y - 0.08, y + 0.08], color=BL); ax.text(p, y + 0.13, str(v), ha="center", fontsize=6)
    ax.text(-0.2, y, lab, ha="right", va="center", fontsize=8)
ax.set_xlim(-3, 6.4); ax.set_ylim(-2.6, 0.6)
clean(ax, "2 Rungs: HALF-LINE, read again and again", "Each rung up = one more radius taken.\nRung 3 becomes plain adding\nafter two radii.", BL)
ax = axs[0, 2]  # conics
t = np.linspace(0, 2*np.pi, 200); ax.plot(np.cos(t), np.sin(t), color=OR, lw=2, label="circle: angle")
s = np.linspace(-1.2, 1.2, 100); ax.plot(np.cosh(s), np.sinh(s), color=BL, lw=2, label="hyperbola: radius"); ax.plot(-np.cosh(s), np.sinh(s), color=BL, lw=2)
ax.plot([0.0, 0.0], [-1.6, 1.6], color=GR, lw=2, label="line: straight")
ax.legend(fontsize=7, loc="lower center", bbox_to_anchor=(0.5, -0.02)); ax.set_xlim(-2, 2); ax.set_ylim(-1.9, 1.9)
clean(ax, "3 Plane dial: CONICS", "\"Radius 1\" for each sign of k²−4d.\nCircle = angle, hyperbola = radius.", OR)
ax = axs[1, 0]  # sphere, plane, saddle schematic
ax.add_patch(Circle((-1.6, 0), 0.75, fill=False, color=OR, lw=2)); ax.text(-1.6, -1.05, "sphere\nd=+1", ha="center", fontsize=8)
ax.plot([-0.55, 0.55], [0, 0], color=GR, lw=2); ax.text(0, -1.05, "plane\nd=0", ha="center", fontsize=8)
xx = np.linspace(-0.7, 0.7, 50); ax.plot(1.6 + xx, 0.5*xx**2*2 - 0.2, color=BL, lw=2); ax.plot(1.6 + xx, -0.5*xx**2*2 + 0.2, color=BL, lw=2)
ax.text(1.6, -1.05, "saddle\nd=−1", ha="center", fontsize=8); ax.set_xlim(-2.6, 2.6); ax.set_ylim(-1.5, 1.2)
clean(ax, "4 Swap rule: SPHERE / PLANE / SADDLE", "Swap order around a loop:\nturn = d × area (Gauss–Bonnet).", OR)
ax = axs[1, 1]  # generalized permutohedra hexagon slider
def hexa(y12, y123, off):
    pts = []
    for a in np.linspace(0, 2*np.pi, 7)[:-1]: pts.append([math.cos(a), math.sin(a)])
    return np.array(pts)
for i, (sc, col) in enumerate([(0.0, BL), (0.5, BL), (1.0, BL)]):
    tri = np.array([[0, 1], [-0.866, -0.5], [0.866, -0.5]])*(1 - sc)*0.8
    hexv = np.array([[math.cos(a), math.sin(a)] for a in np.linspace(0, 2*np.pi, 7)[:-1] + np.pi/6])*sc*0.8
    pts = np.array([a + b for a in tri for b in hexv]) if sc > 0 else tri
    from scipy.spatial import ConvexHull
    pts = pts[ConvexHull(pts).vertices] if len(pts) > 3 else pts
    ax.add_patch(Polygon(pts + [2.1*i - 2.1, 0], closed=True, fill=False, color=col, lw=2))
ax.set_xlim(-3.2, 3.2); ax.set_ylim(-1.3, 1.3)
clean(ax, "5 Corners: GENERALIZED PERMUTOHEDRA", "Slide corner weights: triangle → hexagon.\nCube, simplex, associahedron,\npermutohedron are one family.", OR)
ax = axs[1, 2]  # simplex and elliptope
tri = np.array([[0, 0], [1, 0], [0.5, 0.866]]); ax.add_patch(Polygon(tri - [1.3, 0.3], fill=False, color=OR, lw=2))
ax.text(-0.8, -0.55, "simplex\n(probabilities)", ha="center", fontsize=8)
sq = np.array([[-1, -1], [1, -1], [1, 1], [-1, 1]])*0.5 + [0.9, 0.15]; ax.add_patch(Polygon(sq, fill=False, color=GY, lw=1.5))
tt = np.linspace(0, 2*np.pi, 200); r = 0.5/(np.abs(np.cos(tt))**3 + np.abs(np.sin(tt))**3)**(1/3)
ax.plot(0.9 + r*np.cos(tt), 0.15 + r*np.sin(tt), color=OR, lw=2)
dm = np.array([[0, 1], [1, 0], [0, -1], [-1, 0]])*0.5 + [0.9, 0.15]; ax.add_patch(Polygon(dm, fill=False, color=GY, lw=1.5, ls="--"))
ax.text(0.9, -0.55, "elliptope between\nlocal and no-signal", ha="center", fontsize=8); ax.set_xlim(-1.6, 1.6); ax.set_ylim(-0.9, 0.9)
clean(ax, "6 Probability: SIMPLEX; quantum: ELLIPTOPE", "Independent events keep the cube.\nCorrelations: cut polytope ⊂ elliptope ⊂ cube.", OR)
ax = axs[2, 0]  # Lumley triangle-ish: traceless strain plane
h = np.array([[0, 1], [-0.866, -0.5], [0.866, -0.5]]); ax.add_patch(Polygon(h, fill=False, color=OR, lw=2))
ax.scatter([0.866], [-0.5], color=RD, s=60, zorder=3); ax.text(0.9, -0.75, "two stretch,\none squeeze", fontsize=8, color=RD)
ax.set_xlim(-1.3, 1.6); ax.set_ylim(-1.1, 1.2)
clean(ax, "7 Flows: STRAIN TRIANGLE", "Stretch rates sum to 0 (a plane).\nPolar adds nothing new here.", GY)
ax = axs[2, 1]  # torus
u = np.linspace(0, 2*np.pi, 60); v = np.linspace(0, 2*np.pi, 13)
for vv in v: ax.plot((1 + 0.35*np.cos(vv))*np.cos(u), 0.45*(1 + 0.35*np.cos(vv))*np.sin(u) + 0.35*np.sin(vv), color=OR, lw=0.8)
ax.set_xlim(-1.6, 1.6); ax.set_ylim(-1.0, 1.0)
clean(ax, "8 Number theory: TORUS OF ANGLES", "Roots of unity live on the angles.\nRadii carry nothing about n.", GY)
ax = axs[2, 2]  # log axis
ax.plot([0, 6], [0, 0], color="k")
for i, v in enumerate([1e-7, 1e-5, 1e-3, 1e-1, 10]):
    ax.plot([i*1.4 + 0.2]*2, [-0.1, 0.1], color=OR); ax.text(i*1.4 + 0.2, 0.2, f"{v:g}", ha="center", fontsize=7)
ax.text(3, -0.45, "pH, log K, log rate, log return", ha="center", fontsize=8); ax.set_xlim(-0.3, 6.3); ax.set_ylim(-0.8, 0.6)
clean(ax, "9 Applied sciences: LOG AXIS", "Already written in radius\ncoordinates all along.", GY)
fig.suptitle("The right shape for each family (orange = known shape, blue = our framing)", fontsize=13)
save(fig, "fig1-right-shape-per-family")

# ---------- fig2: core law in radius coordinates, flip disappears
fig, axs = plt.subplots(1, 2, figsize=(10, 4))
k = 1.0; a = 0.15; xs = [0.0]
for _ in range(14): xs.append(xs[-1] + a + k*xs[-1]*a)
axs[0].plot(range(15), xs, "o-", color=OR, label="Cartesian: x, x∘a, x∘a∘a, …  (spreads out)")
axs[0].plot(range(15), [math.log1p(x) for x in xs], "o-", color=BL, label="radius log(1+x): even steps")
axs[0].set_xlabel("number of compositions"); axs[0].legend(fontsize=8); axs[0].set_title("Composing = adding radii (A1, A3)")
x = np.linspace(-0.45, 2, 300)
for kk, col in [(1, BL), (0.5, OR)]:
    axs[1].plot(x, np.log1p(kk*x)/kk, color=col, lw=3, alpha=0.5, label=f"radius of x, k={kk}")
    ph = x/(1 + kk*x); axs[1].plot(x, np.log1p(-kk*ph)/(-kk), "--", color="k", lw=1)
axs[1].plot([], [], "--", color="k", label="radius of flip φ(x), read with −k")
axs[1].legend(fontsize=8); axs[1].set_title("The flip changes nothing radial (A2)"); axs[1].set_xlabel("x")
save(fig, "fig2-core-law-radius")

# ---------- fig3: elliptic degeneration = count of angles
fig, axs = plt.subplots(1, 3, figsize=(11, 3.6))
ax = axs[0]
for vv in np.linspace(0, 2*np.pi, 14): ax.plot((1 + 0.4*np.cos(vv))*np.cos(u), 0.45*(1 + 0.4*np.cos(vv))*np.sin(u) + 0.4*np.sin(vv), color=OR, lw=0.8)
ax.set_title("torus: 2 angles\n(no wall, elliptic law)"); ax.set_aspect("equal"); ax.axis("off")
ax = axs[1]
for yy in np.linspace(-1, 1, 7): ax.plot(np.cos(u)*0.6, 0.25*np.sin(u) + yy, color=BL, lw=0.8)
ax.plot([-0.6, -0.6], [-1, 1], color=BL); ax.plot([0.6, 0.6], [-1, 1], color=BL)
ax.set_title("cylinder: 1 angle + 1 radius\n(mult-add, ∘ₖ)"); ax.set_aspect("equal"); ax.axis("off")
ax = axs[2]
for yy in np.linspace(-1, 1, 7): ax.plot([-1, 1], [yy, yy], color=GR, lw=0.8); ax.plot([yy, yy], [-1, 1], color=GR, lw=0.8)
ax.set_title("plane: 0 angles\n(addition)"); ax.set_aspect("equal"); ax.axis("off")
fig.suptitle("Pinch the torus once, then again: the shape-dependency line counts angles 2 → 1 → 0", fontsize=11)
save(fig, "fig3-count-of-angles")

# ---------- fig4: swap rule
fig, axs = plt.subplots(1, 3, figsize=(12, 3.8))
def Rx(a): c, s = math.cos(a), math.sin(a); return np.array([[1, 0, 0], [0, c, -s], [0, s, c]])
def Ry(a): c, s = math.cos(a), math.sin(a); return np.array([[c, 0, s], [0, 1, 0], [-s, 0, c]])
def Bx(r): M = np.eye(4); M[0, 0] = M[1, 1] = math.cosh(r); M[0, 1] = M[1, 0] = -math.sinh(r); return M
def By(r): M = np.eye(4); M[0, 0] = M[2, 2] = math.cosh(r); M[0, 2] = M[2, 0] = -math.sinh(r); return M
areas = np.linspace(0.0005, 0.09, 25); rt, bt = [], []
for A in areas:
    s = math.sqrt(A); C = Rx(s) @ Ry(s) @ Rx(-s) @ Ry(-s); rt.append(math.atan2(C[1, 0], C[0, 0]))
    Cb = Bx(s) @ By(s) @ Bx(-s) @ By(-s); bt.append(math.atan2(Cb[2, 1], Cb[1, 1]))
axs[0].plot(areas, np.abs(rt), color=OR, lw=2, label="turns (sphere): +area")
axs[0].plot(areas, -np.abs(bt), color=BL, lw=2, label="boosts (saddle): −area"); axs[0].plot(areas, 0*areas, color=GR, lw=2, label="Galilean (flat): 0")
axs[0].set_xlabel("area of the swap loop"); axs[0].set_ylabel("turn"); axs[0].legend(fontsize=8); axs[0].set_title("Precession: turn = d × area (D1)")
ax = axs[1]
for off, pos, col in [(0, (0, 1, 2), GR), (4.2, (1, 3, 4), OR), (8.4, (3, 4, 5), BL)]:
    path = [(0, 0)]; xq = yq = 0
    for st in range(6):
        if st in pos: yq += 1
        else: xq += 1
        path.append((xq, yq))
    P = np.array(path, float) + [off, 0]; area = sum(p - i for i, p in enumerate(pos))
    ax.fill_between(P[:, 0], P[:, 1], 3, color=col, alpha=0.25)
    ax.plot(P[:, 0], P[:, 1], color=col, lw=2); ax.add_patch(Polygon([[off, 0], [off+3, 0], [off+3, 3], [off, 3]], fill=False, color=GY, lw=0.6))
    ax.text(off + 1.5, -0.7, f"pays q^{area}", ha="center", color=col)
ax.set_xlim(-0.5, 12); ax.set_ylim(-1.2, 3.5); ax.set_aspect("equal"); ax.axis("off")
ax.set_title("Cube face weights: a path pays q^(boxes above it)\nsum over paths = q-binomial (D2, D3)")
js = np.array([2, 5, 10, 20, 40, 80, 160, 320]); ph = []
def spin(j):
    m = np.arange(j, -j - 1, -1); dim = len(m); Jp = np.zeros((dim, dim))
    for i in range(1, dim): Jp[i - 1, i] = math.sqrt(j*(j + 1) - m[i]*(m[i] + 1))
    return (Jp + Jp.T)/2, (Jp - Jp.T)/(2j)
for j in js:
    Jx, Jy = spin(j); X, P = Jx/math.sqrt(j), Jy/math.sqrt(j); a = 0.4
    Uc = expm(-1j*a*X) @ expm(-1j*a*P) @ expm(1j*a*X) @ expm(1j*a*P); ph.append(abs(cmath.phase(Uc[0, 0]))/a**2)
axs[2].semilogx(js, ph, "o-", color=BL); axs[2].axhline(1, color=OR, ls="--", label="flat quantum plane: phase = area/ħ")
axs[2].set_xlabel("spin j (sphere grows)"); axs[2].set_ylabel("phase ÷ area"); axs[2].legend(fontsize=8)
axs[2].set_title("Spin sphere flattens into\nquantum phase space (D4)")
save(fig, "fig4-swap-costs-a-turn")

# ---------- fig5: statistical models and their shapes (moment maps)
fig, axs = plt.subplots(1, 3, figsize=(12, 4))
th = rng.normal(scale=2, size=(3000, 2)); m = 1/(1 + np.exp(-th))
axs[0].scatter(*m.T, s=2, color=BL, alpha=0.4); axs[0].add_patch(Polygon([[0, 0], [1, 0], [1, 1], [0, 1]], fill=False, color=OR, lw=2))
axs[0].set_title("Two independent yes/no:\nlogistic fills the SQUARE"); axs[0].set_aspect("equal")
th = rng.normal(scale=2, size=(3000, 3)); p = np.exp(th)/np.exp(th).sum(1, keepdims=True)
Pm = p[:, 1:2]*np.array([1, 0]) + p[:, 2:3]*np.array([0.5, 0.866])
axs[1].scatter(*Pm.T, s=2, color=BL, alpha=0.4); axs[1].add_patch(Polygon([[0, 0], [1, 0], [0.5, 0.866]], fill=False, color=OR, lw=2))
axs[1].set_title("One choice of three:\nsoftmax fills the TRIANGLE"); axs[1].set_aspect("equal")
pts = []
for _ in range(1500):
    w = np.exp(rng.normal(scale=1.5, size=3)); E = np.zeros(3)
    for perm in itertools.permutations(range(3)):
        pr, rest = 1, w.sum()
        for it in perm: pr *= w[it]/rest; rest -= w[it]
        for pos, it in enumerate(perm): E[it] += pr*pos
    pts.append(E)
pts = np.array(pts); B = np.array([[1, 0], [-0.5, 0.866], [-0.5, -0.866]]); Q = pts @ B
hexv = np.array([np.array(perm) @ B for perm in itertools.permutations(range(3))])
from scipy.spatial import ConvexHull
hv = hexv[ConvexHull(hexv).vertices]
axs[2].scatter(*Q.T, s=2, color=BL, alpha=0.4); axs[2].add_patch(Polygon(hv, fill=False, color=OR, lw=2))
axs[2].set_title("Ranking three items:\nexpected ranks fill the HEXAGON (E1)"); axs[2].set_aspect("equal")
for a in axs: a.set_xticks([]); a.set_yticks([])
fig.suptitle("Statistics: each model has its own shape (moment maps of exponential families)", fontsize=11)
save(fig, "fig5-statistics-shapes")

# ---------- fig6: quantum correlations straight in angle coordinates
fig, axs = plt.subplots(1, 2, figsize=(10, 4.6))
S = np.linspace(-1, 1, 601); s, t = np.meshgrid(S, S)
Ang = np.linspace(-np.pi/2, np.pi/2, 601); Ag, Bg = np.meshgrid(Ang, Ang)
loc = np.ones_like(s, bool); qua = np.ones_like(s, bool)
combos = [(3*s - t), (s + t), (s + t), (s + t)]
for sg in itertools.product([1, -1], repeat=1): pass
E = [s, s, s, t]
def regions(E):
    L = np.ones_like(s, bool); Qr = np.ones_like(s, bool)
    for kpos in range(4):
        sg = np.ones(4); sg[kpos] = -1
        L &= np.abs(sum(g*e for g, e in zip(sg, E))) <= 2 + 1e-9
        Qr &= np.abs(sum(g*np.arcsin(e) for g, e in zip(sg, E))) <= np.pi + 1e-9
    return L, Qr
L, Qr = regions(E)
axs[0].contourf(s, t, Qr.astype(float), levels=[0.5, 1.5], colors=[OR], alpha=0.35)
axs[0].contour(s, t, L.astype(float), levels=[0.5], colors=[GY]); axs[0].contour(s, t, Qr.astype(float), levels=[0.5], colors=[OR])
axs[0].set_title("Correlations E (slice E11=E12=E21=s, E22=t)\nquantum (orange) is ROUND; local (grey) straight")
axs[0].set_xlabel("s"); axs[0].set_ylabel("t"); axs[0].set_aspect("equal")
s, t = np.sin(Ag), np.sin(Bg); L, Qr = regions([s, s, s, t]); A, Bv = Ag, Bg
axs[1].contourf(A, Bv, Qr.astype(float), levels=[0.5, 1.5], colors=[OR], alpha=0.35)
axs[1].contour(A, Bv, L.astype(float), levels=[0.5], colors=[GY]); axs[1].contour(A, Bv, Qr.astype(float), levels=[0.5], colors=[OR])
axs[1].set_title("Same region in angle coordinates θ = arcsin E\nquantum is now STRAIGHT; local bends (E4)")
axs[1].set_xlabel("arcsin s"); axs[1].set_ylabel("arcsin t"); axs[1].set_aspect("equal")
save(fig, "fig6-quantum-correlations-angles")

# ---------- fig7: entropy as a shape
def ent_vec(S1, S2, S3):
    return {(): 0, (0,): S1, (1,): S1, (2,): S1, (0, 1): S2, (0, 2): S2, (1, 2): S2, (0, 1, 2): S3}
def verts(Sf):
    V = []
    for perm in [(0, 1, 2), (1, 0, 2), (1, 2, 0), (2, 1, 0), (2, 0, 1), (0, 2, 1)]:
        x = np.zeros(3); prev = ()
        for i in perm: cur = tuple(sorted(prev + (i,))); x[i] = Sf[cur] - Sf[prev]; prev = cur
        V.append(x)
    return np.array(V)
cases = [("4-qubit GHZ\nI3 = +1: triangle ADDED\n(breaks holography)", ent_vec(1, 1, 1), RD),
         ("classical secret (XOR)\nI3 = −1: triangle CARVED", ent_vec(1, 2, 2), BL),
         ("perfect tensor AME(4,3)\nI3 = −2: triangle CARVED\n(holographic)", ent_vec(1, 2, 1), BL),
         ("random 4-qubit pure state\n(computed)", None, BL)]
def rand_state_ent():
    psi = rng.normal(size=16) + 1j*rng.normal(size=16); psi /= np.linalg.norm(psi); T = psi.reshape(2, 2, 2, 2)
    def S(keep):
        if not keep: return 0.0
        rest = [i for i in range(4) if i not in keep]; M = np.transpose(T, list(keep) + rest).reshape(2**len(keep), -1)
        sv = np.linalg.svd(M, compute_uv=False)**2; sv = sv[sv > 1e-15]; return float(-(sv*np.log2(sv)).sum())
    return {I: S(I) for r in range(4) for I in itertools.combinations(range(3), r)}
fig, axs = plt.subplots(1, 4, figsize=(14, 4))
Bm = np.array([[1, 0], [-0.5, 0.866], [-0.5, -0.866]])
for ax, (title, Sf, col) in zip(axs, cases):
    if Sf is None:
        Sf = rand_state_ent()
        I3 = Sf[(0,)] + Sf[(1,)] + Sf[(2,)] - Sf[(0, 1)] - Sf[(0, 2)] - Sf[(1, 2)] + Sf[(0, 1, 2)]
        title = f"random 4-qubit pure state\nI3 = {I3:+.2f}" + (": carved" if I3 < 0 else ": added")
    V = verts(Sf) @ Bm; V = V - V.mean(0)
    ax.add_patch(Polygon(V, closed=True, fill=True, alpha=0.25, color=col)); ax.add_patch(Polygon(V, closed=True, fill=False, color=col, lw=2))
    ax.scatter(*V.T, color="k", s=10)
    lim = np.abs(V).max()*1.4 + 0.3; ax.set_xlim(-lim, lim); ax.set_ylim(-lim, lim); ax.set_aspect("equal"); ax.axis("off"); ax.set_title(title, fontsize=9)
fig.suptitle("Entropy of three regions as a shape: the top corner weight I3 adds a triangle (points right) or carves one (points left) (E5)", fontsize=11, y=1.02)
save(fig, "fig7-entropy-shapes")

# ---------- fig8: the 5-party cyclic inequality on the corners
fig, ax = plt.subplots(figsize=(8, 4.2))
groups = {1: [], 2: [], 3: [], 4: [], 5: []}
coef = {}
cyc = lambda i, m: tuple(sorted((i + j) % 5 for j in range(m)))
lhs = [cyc(i, 3) for i in range(5)]; rhs = [cyc(i, 2) for i in range(5)] + [tuple(range(5))]
for r in range(1, 6):
    for J in itertools.combinations(range(5), r):
        c = sum(set(J) <= set(I) for I in lhs) - sum(set(J) <= set(I) for I in rhs); groups[r].append((J, c))
for r, items in groups.items():
    for i, (J, c) in enumerate(items):
        x = i - (len(items) - 1)/2; col = RD if c < 0 else ("#dddddd" if c == 0 else GR)
        ax.scatter(x, r, s=260, color=col, edgecolor="k", zorder=3); ax.text(x, r - 0.33, "".join("ABCDE"[q] for q in J), ha="center", fontsize=6)
ax.set_yticks(range(1, 6)); ax.set_yticklabels(["single", "pairs", "triples", "fours", "all five"])
ax.set_xticks([]); ax.set_ylim(0.4, 5.5)
ax.set_title("5-region holographic inequality on the cube's corners: red = weight −1, grey = untouched (E6, E7)\n"
             "Singles and pairs (segments) never appear; only simplices of 3+ regions. Known: 'superbalance'.", fontsize=9)
save(fig, "fig8-holographic-inequality-on-corners")

# ---------- fig9: KdV 2-soliton with the polar (tropical) grid
k1, k2 = 2.0, 1.0; Aw = ((k1 - k2)/(k1 + k2))**2
x = np.linspace(-25, 25, 500); tt = np.linspace(-6, 6, 300); X, T = np.meshgrid(x, tt)
e1 = k1*X - k1**3*T; e2 = k2*X - k2**3*T
m = np.maximum.reduce([np.zeros_like(X), e1, e2, e1 + e2 + math.log(Aw)])
fig, ax = plt.subplots(figsize=(8, 4.5))
ax.contourf(X, T, m, levels=40, cmap="Blues", alpha=0.6)
lab = np.argmax(np.stack([np.zeros_like(X), e1, e2, e1 + e2 + math.log(Aw)]), axis=0)
ax.contour(X, T, lab, levels=[0.5, 1.5, 2.5], colors=[OR], linewidths=2)
ax.set_xlabel("x"); ax.set_ylabel("t")
ax.set_title("Two solitons = the polar grid of the cube 1+X1+X2+A·X1X2 (orange lines)\n"
             "The pair weight A is a height: log A is the jog in the slow wave's line (E2)", fontsize=10)
save(fig, "fig9-solitons-on-polar-grid")

# ---------- fig10: Lee-Yang and Fermat
fig, axs = plt.subplots(1, 2, figsize=(10, 4.4))
edges = [(0, 1), (1, 2), (2, 3), (3, 0), (0, 2), (1, 3)]
for aa, col in [(0.3, BL), (0.7, OR)]:
    coeffs = np.zeros(5)
    for Sx in itertools.product([0, 1], repeat=4):
        cut = sum(Sx[i] != Sx[j] for i, j in edges); coeffs[sum(Sx)] += aa**cut
    r = np.roots(coeffs[::-1]); axs[0].scatter(r.real, r.imag, color=col, s=50, label=f"coupling a={aa}")
axs[0].add_patch(Circle((0, 0), 1, fill=False, color="k", ls="--")); axs[0].set_aspect("equal"); axs[0].legend(fontsize=8)
axs[0].set_title("Lee–Yang: every zero has log-radius 0\n(pure angle) (E3)")
for nexp, col in [(1, "k"), (3, OR), (5, BL), (9, GR)]:
    Xs = rng.normal(size=4000)*3 + 1j*rng.normal(size=4000)*3; Ys = 1 - Xs
    axs[1].scatter(np.log(np.abs(Xs))/nexp, np.log(np.abs(Ys))/nexp, s=1, color=col, alpha=0.4, label=f"n={nexp}")
axs[1].set_xlim(-3, 3); axs[1].set_ylim(-3, 3); axs[1].set_aspect("equal"); axs[1].legend(fontsize=8, markerscale=6)
axs[1].set_xlabel("log|x|"); axs[1].set_ylabel("log|y|")
axs[1].set_title("Fermat xⁿ+yⁿ=1 in radius coordinates:\nthe line's picture shrunk by n (F1)")
save(fig, "fig10-lee-yang-and-fermat")

# ---------- fig11: the sorting board
cols = [("Gets simpler\nor unifies", GR, ["1 Core law:\nradius adds,\nflip vanishes", "2 Rungs:\nrung = depth\nof radius",
          "3 Plane dial:\ncircle vs hyperbola;\nangles 2→1→0", "4 Swap rule:\none area law", "6 Statistics:\none shape per model",
          "6 Quantum correlations:\nstraight in angles", "6 Lee–Yang:\nzeros at radius 0"]),
        ("Already\npolar", OR, ["6 Information:\nRényi = radius\nof Tsallis", "6 Repetition dial:\noccupation = flip",
          "6 Heat:\ntemperature =\nroundness", "9 Sciences:\nlog axes"]),
        ("Unchanged", GY, ["7 Navier–Stokes", "Newton's G lead", "Curved space\nby counting", "8 Hopf one,\nKummer–Wilson"]),
        ("Harder /\nno help", RD, ["8 Fermat:\nall content\nin the angles"]),
        ("Opens a\nlead", BL, ["5+6 Entropy\nas a shape", "4+5 Swap = area:\nprecession, solitons, ħ", "5 Pair weights\n= heights"])]
fig, axs = plt.subplots(1, 5, figsize=(14, 8))
for ax, (h, c, items) in zip(axs, cols):
    ax.axis("off"); ax.set_xlim(0, 1); ax.set_ylim(0, 1); ax.set_title(h, color=c, fontweight="bold")
    for j, it in enumerate(items):
        ax.text(0.5, 0.94 - j*0.135, it, ha="center", va="top", fontsize=9, bbox=dict(boxstyle="round", fc=c, alpha=0.15, ec=c))
fig.suptitle("Where every family lands in the polar view", fontsize=13)
save(fig, "fig11-sorting-board")
print("figures written")
