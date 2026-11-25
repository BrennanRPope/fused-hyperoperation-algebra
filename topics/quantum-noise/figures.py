"""Pictures for the quantum-noise page. Run after checks.py (uses surface_results.npy)."""
import numpy as np, itertools, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from math import comb
B, O, A, Yl, V, G = "#257de0", "#f47230", "#00b894", "#f2a600", "#6a4ce2", "#8a8a8a"
TXT, TXT2 = "#191919", "#666666"
plt.rcParams.update({"font.size": 11, "axes.edgecolor": "#bbbbbb", "axes.labelcolor": TXT, "xtick.color": TXT2,
    "ytick.color": TXT2, "axes.spines.top": False, "axes.spines.right": False, "axes.grid": True,
    "grid.color": "#eeeeee", "lines.linewidth": 2, "figure.dpi": 130})
F = "figures/"

# q1: three rules, one shape
fig, ax = plt.subplots(1, 2, figsize=(11, 4.2))
n = np.arange(0, 41); p = 0.05
rules = [("amplitude damping (T1), k = −1", -1, A), ("depolarizing, k = −4/3", -4/3, B), ("bit-flip / phase-flip, k = −2", -2, O)]
for lab, k, c in rules:
    tot = (1 - (1 + k*p)**n) / (-k)
    ax[0].plot(n, tot, color=c, label=lab); ax[0].axhline(-1/k, color=c, ls=":", lw=1.2)
    ax[0].text(40.5, -1/k, f" wall {-1/k:g}", color=TXT2, va="center", fontsize=9)
    ax[1].plot(n, (1 + k*p)**n, color=c)
ax[0].set_xlabel("noise channels stacked (each 5%)"); ax[0].set_ylabel("combined error chance")
ax[0].set_title("Stacked noise climbs to a wall at −1/k", color=TXT, loc="left"); ax[0].legend(frameon=False, fontsize=9, loc="lower right")
ax[1].set_yscale("log"); ax[1].set_xlabel("noise channels stacked"); ax[1].set_ylabel("1 + k·(error)   (log scale)")
ax[1].set_title("In the (1 + k·a) coordinate it is a straight line", color=TXT, loc="left")
fig.tight_layout(); fig.savefig(F+"q1_three_rules.png"); plt.close(fig)

# q2: one cube, two colorings
fig = plt.figure(figsize=(11, 4.8))
corners = list(itertools.product([0,1], repeat=3))
def proj(c): x, y, z = c; return np.array([x + 0.45*z, y + 0.3*z])
edges = [(a, b) for a in corners for b in corners if sum(abs(i-j) for i, j in zip(a, b)) == 1 and a < b]
for k, (title, rule) in enumerate([("Amplitude damping: any fault counts (k = −1)", lambda c: sum(c) > 0),
                                   ("Bit-flip: only an odd number of faults flips (k = −2)", lambda c: sum(c) % 2 == 1)]):
    ax = fig.add_subplot(1, 2, k+1); ax.set_axis_off(); ax.set_aspect("equal")
    for a, b in edges:
        pa, pb = proj(a), proj(b); ax.plot([pa[0], pb[0]], [pa[1], pb[1]], color="#cccccc", lw=1.5, zorder=1)
    for c in corners:
        q = proj(c); bad = rule(c)
        ax.scatter(*q, s=900, color=(O if bad else "#ffffff"), edgecolor=(O if bad else G), lw=2, zorder=2)
        ax.text(q[0], q[1], "".join(map(str, c)), ha="center", va="center", fontsize=9, color=("white" if bad else TXT), zorder=3)
        w = "".join(["p" if b else "(1−p)" for b in c])
        ax.text(q[0], q[1]-0.2, w, ha="center", va="top", fontsize=7, color=TXT2)
    ax.set_title(title, color=TXT, fontsize=11)
    ax.set_xlim(-0.4, 1.9); ax.set_ylim(-0.45, 1.55)
fig.text(0.5, 0.02, "Three stacked channels = the 3-cube of fault patterns (1 = that channel fired). Orange corners are the ones that hurt.\n"
         "Left: 1 − (1−p)³  (layers 3, 3, 1).   Right: 3p(1−p)² + p³ = (1 − (1−2p)³)/2  (layers 3 and 1, the even layer cancels).",
         ha="center", fontsize=9.5, color=TXT2)
fig.subplots_adjust(bottom=0.18); fig.savefig(F+"q2_fault_cube.png"); plt.close(fig)

# q3: tetrahedron inside the cube
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
fig = plt.figure(figsize=(7.5, 7)); ax = fig.add_subplot(projection="3d")
cub = np.array(list(itertools.product([-1,1], repeat=3)))
for a in cub:
    for b in cub:
        if np.sum(a != b) == 1 and tuple(a) < tuple(b): ax.plot(*zip(a, b), color="#cccccc", lw=1)
T = np.array([[1,1,1],[1,-1,-1],[-1,1,-1],[-1,-1,1]])
faces = [[T[i] for i in f] for f in [(0,1,2),(0,1,3),(0,2,3),(1,2,3)]]
ax.add_collection3d(Poly3DCollection(faces, alpha=0.12, facecolor=B, edgecolor=B, lw=1.5))
for t, lab in zip(T, ["identity (no noise)", "X  (certain bit-flip)", "Y", "Z  (certain phase-flip)"]):
    ax.scatter(*t, s=50, color=B); ax.text(*(t*1.12 + (np.array([0,0,0.12]) if t[0]==1 and t[1]==1 else 0)), lab, fontsize=9, color=TXT)
d = np.linspace(-1/3, 1, 50); ax.plot(d, d, d, color=O, lw=2.5); ax.text(-0.45, -0.45, -0.6, "depolarizing line", color=TXT, fontsize=9)
ax.scatter(0, 0, 0, s=60, color=G); ax.text(0.05, 0.05, -0.15, "total noise (the −1 of mult-add)", fontsize=8.5, color=TXT2)
lam = np.array([0.95, 0.80, 0.85]); path = np.array([lam**m for m in range(0, 25)])
ax.plot(*path.T, color=A, lw=2, marker="o", ms=3); ax.text(0.55, 0.2, -0.35, "same channel\napplied again and again", fontsize=8.5, color=TXT)
ax.set_xlabel("λx"); ax.set_ylabel("λy"); ax.set_zlabel("λz"); ax.view_init(elev=20, azim=35)
ax.set_title("Every Pauli noise channel is a point in this tetrahedron.\nComposing channels multiplies the coordinates.", color=TXT, fontsize=11)
fig.tight_layout(); fig.savefig(F+"q3_tetrahedron.png"); plt.close(fig)

# q4: helices and the echo
rng = np.random.default_rng(3)
fig = plt.figure(figsize=(12, 4.6))
t = np.linspace(0, 1, 400); dets = np.array([-2.1, -1.2, -0.5, 0.3, 0.9, 1.6, 2.4])
for k, echo in enumerate([False, True]):
    ax = fig.add_subplot(1, 3, k+1, projection="3d")
    for dlt in dets:
        ph = dlt * 2*np.pi * t if not echo else dlt * 2*np.pi * np.where(t < 0.5, t, 1 - t)
        ax.plot(t, np.cos(ph + 4*np.pi*t*0), np.sin(ph), lw=1.2, color=B, alpha=0.8)
    ax.scatter([1]*len(dets), *[f(dets*2*np.pi*(0 if echo else 1)) for f in (np.cos, np.sin)], color=O, s=25, zorder=5)
    ax.set_xlabel("time"); ax.set_yticks([]); ax.set_zticks([])
    ax.set_title("No echo: helices of different pitch fan out" if not echo else "Echo at the halfway mark: mirror, unwind, refocus", fontsize=10, color=TXT)
    ax.view_init(elev=18, azim=-60)
ax = fig.add_subplot(1, 3, 3)
T2 = np.linspace(0, 3, 200); sig = 1.5*2*np.pi
ax.plot(T2, np.exp(-(sig*T2)**2/2), color=O, label="no echo")
# slow (Ornstein–Uhlenbeck) noise, echo: simulate
N, steps = 3000, 900; tau = 1.0; dt = 3.0/steps
dl = rng.normal(0, sig, N); phF = np.zeros(N); hist = []
d_hist = np.zeros((steps, N))
for s in range(steps):
    dl = dl - dl/tau*dt + sig*np.sqrt(2*dt/tau)*rng.normal(size=N); d_hist[s] = dl
cum = np.cumsum(d_hist, 0)*dt
echo_c = []
for s in range(1, steps, 6):
    h = s//2; ph = cum[h-1] - (cum[s-1] - cum[h-1]) if h > 0 else -cum[s-1]
    echo_c.append((s*dt, abs(np.mean(np.exp(1j*ph)))))
echo_c = np.array(echo_c)
ax.plot(T2, np.ones_like(T2), color=A, label="echo, frozen noise")
ax.plot(echo_c[:,0], echo_c[:,1], color=B, label="echo, slowly drifting noise")
ax.set_xlabel("time"); ax.set_ylabel("coherence left"); ax.legend(frameon=False, fontsize=9)
ax.set_title("Average over many qubits", fontsize=10, color=TXT, loc="left")
fig.tight_layout(); fig.savefig(F+"q4_echo_helix.png"); plt.close(fig)

# q5: error correction as a flow on the square
f = lambda e: (3*e - e**3)/2
fig, ax = plt.subplots(1, 3, figsize=(13, 4.4))
g = np.linspace(0.0, 1, 60); Xg, Zg = np.meshgrid(g, g)
for a_, (title, mp) in zip(ax[:2], [("3-qubit bit-flip code", lambda x, z: (f(x), z**3)),
                                     ("3-qubit phase-flip code", lambda x, z: (x**3, f(z)))]):
    U, W = mp(Xg, Zg); a_.streamplot(Xg, Zg, U - Xg, W - Zg, color=B, density=1.1, linewidth=1, arrowsize=1)
    a_.set_title(title, color=TXT, loc="left"); a_.set_xlabel("ηx = 1 − 2·(bit-flip chance)"); a_.set_ylabel("ηz = 1 − 2·(phase-flip chance)")
    a_.set_xlim(0, 1); a_.set_ylim(0, 1); a_.scatter([1], [1], color=A, s=70, zorder=5, clip_on=False)
def shor(x, z): return f(z**3), f(x)**3
gg = np.linspace(0.0, 1, 301); Xb, Zb = np.meshgrid(gg, gg); x, z = Xb.copy(), Zb.copy()
for _ in range(80): x, z = shor(x, z)
ok = (x > 0.99) & (z > 0.99)
ax[2].contourf(Xb, Zb, ok.astype(float), levels=[-0.5, 0.5, 1.5], colors=["#f6f6f6", "#cfe6fb"])
ax[2].contour(Xb, Zb, ok.astype(float), levels=[0.5], colors=[B], linewidths=2)
ax[2].text(0.905, 0.93, "noise crushed", ha="center", fontsize=9.5, color=TXT)
ax[2].text(0.66, 0.74, "noise wins", ha="center", fontsize=9.5, color=TXT2)
ax[2].text(0.807, 0.605, " bit-flips up to 9.67%", fontsize=9, color=TXT, rotation=90, va="bottom")
ax[2].text(0.605, 0.853, "phase-flips up to 7.47%", fontsize=9, color=TXT, va="bottom")
ax[2].set_title("Concatenated Shor code: the safe region is a box", color=TXT, loc="left", fontsize=10.5)
ax[2].set_xlabel("ηx"); ax[2].set_ylabel("ηz"); ax[2].set_xlim(0.6, 1); ax[2].set_ylim(0.6, 1); ax[2].grid(False)
fig.tight_layout(); fig.savefig(F+"q5_code_flow.png"); plt.close(fig)

# q6: repetition code as cube layers
fig, ax = plt.subplots(1, 2, figsize=(11, 4.2))
nn = 7; layers = [comb(nn, j) for j in range(nn+1)]
ax[0].bar(range(nn+1), layers, color=[O if j > 3 else "#d9d9d9" for j in range(nn+1)], width=0.7)
for j, v in enumerate(layers): ax[0].text(j, v + 0.8, str(v), ha="center", fontsize=9, color=TXT2)
ax[0].set_xlabel("layer of the 7-cube (how many qubits flipped)"); ax[0].set_ylabel("corners in the layer")
ax[0].set_title("Majority vote on 7 copies: orange layers fail", color=TXT, loc="left"); ax[0].grid(False)
pp = np.linspace(0, 1, 300)
for nn, c in zip([1, 3, 5, 9, 21], [G, A, B, V, O]):
    t = (nn-1)//2; pl = sum(comb(nn, j)*pp**j*(1-pp)**(nn-j) for j in range(t+1, nn+1))
    ax[1].plot(pp, pl, color=c, label=(f"{nn} copies" if nn > 1 else "1 copy"))
ax[1].axvline(0.5, color=G, ls=":", lw=1); ax[1].set_xlabel("flip chance per copy"); ax[1].set_ylabel("chance the vote is wrong")
ax[1].set_title("Every curve crosses at 1/2, the wall for k = −2", color=TXT, loc="left"); ax[1].legend(frameon=False, fontsize=9)
fig.tight_layout(); fig.savefig(F+"q6_repetition_layers.png"); plt.close(fig)

# q7: surface code
surf = np.load("surface_results.npy", allow_pickle=True).item()
fig, ax = plt.subplots(figsize=(6.5, 4.3))
for (L, row), c in zip(surf.items(), [G, A, B, O]):
    r = np.array(row); ax.plot(r[:,0], r[:,1], marker="o", ms=4, color=c, label=f"{L}×{L} grid")
ax.axvline(0.103, color=G, ls=":", lw=1); ax.text(0.1035, 0.03, " known: about 10.3%", fontsize=9, color=TXT2)
ax.set_xlabel("bit-flip chance per qubit"); ax.set_ylabel("chance the code fails")
ax.set_title("Surface code: bigger grids win below the crossing", color=TXT, loc="left"); ax.legend(frameon=False, fontsize=9)
fig.tight_layout(); fig.savefig(F+"q7_surface.png"); plt.close(fig)
print("done")
