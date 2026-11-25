"""Figures for 'Entropy and quantum noise on the cube'. Run: python3 figures.py"""
import itertools, math, os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d.art3d import Poly3DCollection

OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "figures")
os.makedirs(OUT, exist_ok=True)
INK = "#1d2430"; OCT = "#d9473b"; SPH = "#2f6fb0"; CUBE = "#7d8a96"; GOOD = "#2f9e6e"; GOLD = "#c98a12"
plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 11, "axes.edgecolor": INK,
                     "axes.labelcolor": INK, "xtick.color": INK, "ytick.color": INK, "text.color": INK})


def h2(p):
    p = np.clip(p, 1e-300, 1)
    q = np.clip(1 - p, 1e-300, 1)
    return -(p * np.log2(p) + q * np.log2(q))


def Hvec(P):
    P = np.clip(P, 1e-300, None)
    return -(P * np.log2(P)).sum(0)


# --- e1: three norm balls inside the cube, kissing at the tetrahedron's edge midpoints
fig = plt.figure(figsize=(8, 7.2))
ax = fig.add_subplot(111, projection="3d")
for s, e in itertools.combinations(list(itertools.product([-1, 1], repeat=3)), 2):
    if sum(a != b for a, b in zip(s, e)) == 1:
        ax.plot(*zip(s, e), color=CUBE, lw=1, alpha=0.7)
T = {"I": (1, 1, 1), "X": (1, -1, -1), "Y": (-1, 1, -1), "Z": (-1, -1, 1)}
for a, b in itertools.combinations(T.values(), 2):
    ax.plot(*zip(a, b), color=INK, lw=1.8)
for n, v in T.items():
    ax.scatter(*v, color=INK, s=30)
    ax.text(v[0] * 1.13, v[1] * 1.13, v[2] * 1.13, {"I": "no fault", "X": "X", "Y": "Y", "Z": "Z"}[n], fontsize=11)
V = [(1, 0, 0), (-1, 0, 0), (0, 1, 0), (0, -1, 0), (0, 0, 1), (0, 0, -1)]
faces = [[(sx, 0, 0), (0, sy, 0), (0, 0, sz)] for sx in (1, -1) for sy in (1, -1) for sz in (1, -1)]
ax.add_collection3d(Poly3DCollection(faces, facecolor=OCT, alpha=0.28, edgecolor=OCT, lw=0.8))
u, v = np.mgrid[0:2 * np.pi:40j, 0:np.pi:20j]
ax.plot_wireframe(np.cos(u) * np.sin(v), np.sin(u) * np.sin(v), np.cos(v), color=SPH, lw=0.4, alpha=0.45)
for p in V:
    ax.scatter(*p, color=GOLD, s=60, edgecolor=INK, zorder=10)
ax.scatter(0, 0, 0, color=INK, s=25)
ax.text(0.05, 0.05, -0.18, "total noise", fontsize=9)
ax.set_box_aspect((1, 1, 1)); ax.set_axis_off(); ax.view_init(elev=18, azim=32)
ax.set_title("Cube (classical capacity), sphere (collision entropy), octahedron (entanglement breaking)\n"
             "all touch at the six gold points: the tetrahedron's edge midpoints", fontsize=11)
fig.tight_layout(); fig.savefig(os.path.join(OUT, "e1_three_balls.png"), dpi=150); plt.close(fig)

# --- e2: slice through the plane lambda = (u, v/sqrt2, v/sqrt2): the hashing sandwich
r2 = math.sqrt(2)
U, Vv = np.meshgrid(np.linspace(-1.05, 1.05, 700), np.linspace(-1.5, 1.5, 900))
lx, ly, lz = U, Vv / r2, Vv / r2
P = np.stack([1 + lx + ly + lz, 1 + lx - ly - lz, 1 - lx + ly - lz, 1 - lx - ly + lz]) / 4
inside = (P >= -1e-12).all(0)
Hm = np.where(inside, Hvec(np.clip(P, 0, None)), np.nan)
fig, ax = plt.subplots(figsize=(8, 7))
ax.contourf(U, Vv, np.where(inside & (Hm < 1), 1, np.nan), levels=[0.5, 1.5], colors=[GOOD], alpha=0.25)
ax.contour(U, Vv, Hm, levels=[1], colors=[GOOD], linewidths=2.2)
for c in [0.4, 0.7]:
    ax.contour(U, Vv, np.where(inside, np.maximum(np.abs(lx), np.abs(ly)), np.nan), levels=[c], colors=[CUBE], linewidths=1, linestyles=":")
tri = np.array([[1, r2], [1, -r2], [-1, 0], [1, r2]])
ax.plot(tri[:, 0], tri[:, 1], color=INK, lw=2)
t = np.linspace(0, 2 * np.pi, 400)
ax.plot(np.cos(t), np.sin(t), color=SPH, lw=2)
d = np.array([[1, 0], [0, 1 / r2], [-1, 0], [0, -1 / r2], [1, 0]])
ax.fill(d[:, 0], d[:, 1], color=OCT, alpha=0.25); ax.plot(d[:, 0], d[:, 1], color=OCT, lw=2)
pst = 0.8107103
uu = np.linspace(-1, 1, 10)
ax.plot(uu, (4 * pst - 1 - uu) / r2, color=GOLD, lw=1.6, ls="--")
ax.plot(uu, -(4 * pst - 1 - uu) / r2, color=GOLD, lw=1.6, ls="--")
ax.plot([0, 1], [0, r2], color=INK, lw=0.8, ls="-.")
lam = 1 - 4 * 0.18929 / 3
ax.scatter([lam], [r2 * lam], color=GOOD, s=40, zorder=5, edgecolor=INK)
ax.scatter([1, -1], [0, 0], color=GOLD, s=70, zorder=6, edgecolor=INK)
ax.text(1.03, r2 + 0.02, "no fault", fontsize=10); ax.text(1.03, -r2 - 0.08, "X", fontsize=10)
ax.text(-1.05, 0.1, "Y-Z edge\nmidpoint", fontsize=9, ha="right")
ax.text(lam + 0.04, r2 * lam - 0.1, "depolarizing\nhashing point\np = 0.1893", fontsize=9)
ax.text(-1.4, 0.82, "midsphere |λ| = 1", color=SPH, fontsize=10)
ax.text(-0.45, 0.05, "octahedron\n(entanglement\nbreaking)", color=OCT, fontsize=9)
ax.annotate("hashing rate > 0", (0.93, 0.75), (0.05, 1.42), color=GOOD, fontsize=10, arrowprops=dict(arrowstyle="->", color=GOOD))
ax.text(-1.4, 1.2, "gold dashes: some fault chance\nabove 0.8107 (enough)", color=GOLD, fontsize=9)
ax.text(-0.95, -1.2, "grey dots: capacity levels\n(squares of max|λ|)", color=CUBE, fontsize=9)
ax.set_aspect("equal"); ax.set_xlim(-1.45, 1.4); ax.set_ylim(-1.55, 1.55)
ax.set_xlabel("λx"); ax.set_ylabel("√2 · λy  (with λz = λy)")
ax.set_title("A slice of the tetrahedron: positive hashing rate (green) lies outside the sphere\n"
             "and contains the gold-capped corners; it meets the sphere only at the edge midpoint", fontsize=11)
fig.tight_layout(); fig.savefig(os.path.join(OUT, "e2_hashing_slice.png"), dpi=150); plt.close(fig)

# --- e3: the weight is a temperature
fig, ax = plt.subplots(figsize=(8.5, 4.6))
tt = np.logspace(-0.7, 2, 400)
ax.plot(tt, 1 - 1 / tt, color=INK, lw=2.2)
ax.axhline(1, color=CUBE, ls=":", lw=1); ax.axhline(0, color=CUBE, ls=":", lw=1); ax.axhline(-1, color=CUBE, ls=":", lw=1)
marks = [(1, 0, "T = T₀: Shannon, k = 0 (addition)"), (0.5, -1, "T = T₀/2: Tsallis q_T = 2, k = −1 (union rule;\nlinear entropy, the sphere)"),
         (60, 1 - 1 / 60, "T → ∞: q_T → 0, k → 1 (the plain circle operation;\ncounting states, Hartley)")]
for x0, y0, lab in marks:
    ax.scatter([x0], [y0], color=GOLD, s=60, edgecolor=INK, zorder=5)
ax.annotate(marks[0][2], (1, 0), (2.2, -0.45), arrowprops=dict(arrowstyle="-", color=INK), fontsize=10)
ax.annotate(marks[1][2], (0.5, -1), (1.2, -1.35), arrowprops=dict(arrowstyle="-", color=INK), fontsize=10)
ax.annotate(marks[2][2], (60, 1 - 1 / 60), (4, 0.55), arrowprops=dict(arrowstyle="-", color=INK), fontsize=10)
ax.set_xscale("log"); ax.set_ylim(-1.6, 1.25)
ax.set_xlabel("temperature T, in units of the state's own temperature T₀  (Tsallis q_T = T₀/T)")
ax.set_ylabel("weight k = 1 − T₀/T")
ax.set_title("The Tsallis weight is Carnot's efficiency between T and T₀, so it cascades by the union rule", fontsize=11)
fig.tight_layout(); fig.savefig(os.path.join(OUT, "e3_temperature_dial.png"), dpi=150); plt.close(fig)

# --- e4: the 3-bit code, syndromes as body diagonals, entropy split
fig = plt.figure(figsize=(11, 4.8))
ax = fig.add_subplot(121, projection="3d")
for s, e in itertools.combinations(list(itertools.product([0, 1], repeat=3)), 2):
    if sum(a != b for a, b in zip(s, e)) == 1:
        ax.plot(*zip(s, e), color=CUBE, lw=1)
cols = {(0, 0): GOOD, (1, 0): OCT, (1, 1): SPH, (0, 1): GOLD}
for c in itertools.product([0, 1], repeat=3):
    s = (c[0] ^ c[1], c[1] ^ c[2])
    ax.scatter(*c, color=cols[s], s=90, edgecolor=INK, zorder=5)
    ax.text(c[0] + 0.05, c[1] + 0.05, c[2] + 0.07, "".join(map(str, c)), fontsize=9)
    if sum(c) <= 1:
        f = tuple(1 - x for x in c)
        ax.plot(*zip(c, f), color=cols[s], lw=3 if s == (0, 0) else 2, alpha=0.9)
ax.set_axis_off(); ax.set_box_aspect((1, 1, 1)); ax.view_init(elev=20, azim=35)
ax.set_title("Each syndrome is one body diagonal.\nGreen: no alarm. The vote picks the lighter end.", fontsize=10)
ax2 = fig.add_subplot(122)
ps = np.linspace(0.001, 0.5, 300)
HS, HL, HE = [], [], []
for p in ps:
    w = lambda e: p ** sum(e) * (1 - p) ** (3 - sum(e))
    groups = {}
    for e in itertools.product([0, 1], repeat=3):
        groups.setdefault((e[0] ^ e[1], e[1] ^ e[2]), []).append(w(e))
    Ps = np.array([sum(g) for g in groups.values()])
    HS.append(Hvec(Ps[:, None])[0])
    HL.append(sum(sum(g) * h2(g[0] / sum(g)) for g in groups.values()))
    HE.append(3 * h2(p))
ax2.fill_between(ps, 0, HS, color=SPH, alpha=0.35, label="removed by the syndrome")
ax2.fill_between(ps, HS, np.array(HS) + np.array(HL), color=OCT, alpha=0.35, label="left on the diagonal")
ax2.plot(ps, HE, color=INK, lw=1.5, label="all fault entropy, 3·h(p)")
ax2.set_xlabel("flip chance p"); ax2.set_ylabel("bits"); ax2.legend(loc="upper left", fontsize=9)
ax2.set_title("Error correction as entropy removal", fontsize=10)
fig.tight_layout(); fig.savefig(os.path.join(OUT, "e4_three_bit_code.png"), dpi=150); plt.close(fig)

# --- e5: what stacking hides, and the half-bit threshold
fig, (a1, a2) = plt.subplots(1, 2, figsize=(11, 4.4))
n = np.arange(1, 17)
for lam, col in [(0.9, SPH), (0.6, OCT)]:
    p = (1 - lam) / 2
    a1.plot(n, n * h2(p), color=col, lw=1.5, ls="--")
    a1.plot(n, h2((1 - lam ** n) / 2), color=col, lw=2.2, label=f"λ = {lam}")
a1.axhline(1, color=CUBE, ls=":")
a1.set_xlabel("channels stacked, n"); a1.set_ylabel("bits")
a1.set_title("Dashed: entropy of the fault pattern (adds).\nSolid: entropy of the net flip (stays under 1 bit).", fontsize=10)
a1.legend(fontsize=9)
pp = np.linspace(0.001, 0.25, 400)
a2.plot(pp, h2(pp), color=INK, lw=2, label="h(p), one qubit's flip entropy")
a2.plot(pp, 1 - 2 * h2(pp), color=GOOD, lw=1.6, label="1 − 2h(p), random CSS codes' rate")
a2.axhline(0.5, color=CUBE, ls=":"); a2.axhline(0, color=CUBE, lw=0.8)
for x0, lab, col in [(0.103, "matching ≈ 0.103", SPH), (0.1093, "best decoder ≈ 0.1093", OCT), (0.110028, "h(p) = ½: 0.1100", GOLD)]:
    a2.axvline(x0, color=col, lw=1.2)
a2.text(0.114, 0.52, "half a bit", fontsize=9)
a2.text(0.0, -0.95, "blue: matching ≈ 0.103   red: best ≈ 0.1093   gold: h(p) = ½ at 0.1100", fontsize=8.5)
a2.set_ylim(-1.05, 1.05); a2.set_xlabel("flip chance p")
a2.set_title("Surface code: the threshold sits just under half a bit", fontsize=10)
a2.legend(fontsize=8.5, loc="upper right")
fig.tight_layout(); fig.savefig(os.path.join(OUT, "e5_stack_and_threshold.png"), dpi=150); plt.close(fig)
print("figures written to", OUT)
