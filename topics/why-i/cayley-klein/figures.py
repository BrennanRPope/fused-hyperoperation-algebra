"""Pictures for topics/why-i/cayley-klein (Nov 25, 2026). Run: python3 figures.py"""
import numpy as np, matplotlib, os, itertools
matplotlib.use("Agg")
import matplotlib.pyplot as plt
HERE = os.path.dirname(os.path.abspath(__file__)); OUT = os.path.join(HERE, "figures"); os.makedirs(OUT, exist_ok=True)
BLUE, ORANGE, GREEN, RED, GREY, PURPLE = "#2b6cb0", "#dd6b20", "#2f855a", "#c53030", "#718096", "#6b46c1"
plt.rcParams.update({"font.size": 10, "axes.titlesize": 11, "figure.dpi": 130})
def save(fig, name):
    fig.tight_layout(); fig.savefig(os.path.join(OUT, name + ".png"), bbox_inches="tight"); plt.close(fig)

# f1. the nine Cayley-Klein geometries of the plane, and where one o_k dial can sit
names = {(1, 1): "sphere\n(elliptic)", (0, 1): "Euclidean\nplane", (-1, 1): "hyperbolic\nplane",
         (1, 0): "oscillating\nNewton-Hooke", (0, 0): "Galilean", (-1, 0): "expanding\nNewton-Hooke",
         (1, -1): "anti-de Sitter", (0, -1): "Minkowski", (-1, -1): "de Sitter"}
rowcol = {1: "#ebf8ff", 0: "#f7fafc", -1: "#fffaf0"}
fig, ax = plt.subplots(1, 2, figsize=(12.5, 5.6), gridspec_kw={"width_ratios": [1.35, 1]})
a0 = ax[0]
for (k1, k2), nm in names.items():
    a0.add_patch(plt.Rectangle((k1 - .45, k2 - .45), .9, .9, color=rowcol[k2], ec="#cbd5e0"))
    a0.text(k1, k2 + .12, nm, ha="center", va="center", fontsize=9)
a0.axhline(0, color="#a0aec0", lw=.8); a0.axvline(0, color="#a0aec0", lw=.8)
a0.set_xticks([-1, 0, 1]); a0.set_xticklabels(["kappa1 < 0", "kappa1 = 0", "kappa1 > 0"])
a0.set_yticks([-1, 0, 1]); a0.set_yticklabels(["kappa2 < 0\n(Lorentzian)", "kappa2 = 0\n(Newtonian)", "kappa2 > 0\n(Riemannian)"])
a0.set_xlim(-1.6, 1.6); a0.set_ylim(-1.6, 1.85)
a0.set_xlabel("first dial kappa1: curvature (time-translation label)"); a0.set_ylabel("second dial kappa2: signature (boost label)")
a0.annotate("", xy=(0, -1.35), xytext=(0, 1.35), arrowprops=dict(arrowstyle="<->", color=BLUE, lw=2.5, alpha=.6))
a0.plot(0, -1, "o", color=RED, ms=11, mfc="none", mew=2.5)
a0.text(0.12, -1.38, "repo's relativity:\n(k, d) = (0, -1/c^2)\nsets kappa2 = -1/c^2", color=RED, fontsize=8)
a0.annotate("", xy=(1.35, -1.55), xytext=(-1.35, -1.55), arrowprops=dict(arrowstyle="<->", color=GREEN, lw=2.5, alpha=.6), annotation_clip=False)
a0.text(-1.55, 1.62, "a single o_k chain is ONE line through this plane:\nblue if it composes boosts, green if time shifts", fontsize=8, color=GREY)
a0.set_title("Cayley-Klein's nine plane geometries (Yaglom; Herranz-Ortega-Santander)")
a1 = ax[1]; a1.set_xlim(-1, 1); a1.set_ylim(-0.5, 1.6); a1.axis("off")
a1.annotate("", xy=(0.95, 0.9), xytext=(-0.95, 0.9), arrowprops=dict(arrowstyle="->", color="k"))
a1.text(0, 1.05, "one label kappa = -k^2/4  (ours: K2)", ha="center", fontsize=10)
for kv, kp, col, lab in [(-1.0, -0.25, GREEN, "k = -1\n(random walk)"), (0, 0, GREY, "k = 0\n(plain +)"), (-1j, 0.25, BLUE, "k = -i\n(quantum)")]:
    x = kp*3
    a1.plot(x, 0.9, "o", color=col, ms=10); a1.text(x, 0.62, lab, ha="center", fontsize=9, color=col)
    a1.text(x, 1.22 if kv != 0 else 1.22, f"kappa = {kp:+.2f}", ha="center", fontsize=8, color=col)
a1.text(-0.95, 0.3, "kappa < 0: hyperbolic, boost-like,\nlift 1 + k a is a positive number", fontsize=8, color=GREEN)
a1.text(0.15, 0.3, "kappa > 0: elliptic, rotation-like,\nlift 1 + k a is a unit phase", fontsize=8, color=BLUE)
a1.text(0, -0.25, "k -> i k flips the sign: the Cayley-Klein analytic continuation\n(Gromov's j = 1 -> i), a 'Wick rotation' of that one plane", ha="center", fontsize=8.5, color=PURPLE)
a1.set_title("where our dial settings land on one label")
save(fig, "f1-nine-geometries")

# f2. the (k, d) plane is one dial: kappa = d - k^2/4
fig, ax = plt.subplots(figsize=(7.2, 5.8))
X = np.linspace(-1.2, 1.2, 400); D = np.linspace(-1.4, 1.2, 400); XX, DD = np.meshgrid(X, D); KAP = DD - XX
ax.contourf(XX, DD, KAP, levels=[-5, -1e-9, 1e-9, 5], colors=["#fffaf0", "#f7fafc", "#ebf8ff"])
cs = ax.contour(XX, DD, KAP, levels=[-1.2, -0.8, -0.4, 0, 0.4, 0.8], colors=GREY, linewidths=.8)
ax.clabel(cs, fmt=lambda v: f"kappa={v:+.1f}", fontsize=7)
ax.plot(X, X, color="k", lw=1.5)
pts = [(0.25, 0, "k = -1 (d = 0)", GREEN), (-0.25, 0, "k = -i (d = 0)", BLUE), (0, 0, "", GREY), (0, -1, "velocity addition\n(k, d) = (0, -1/c^2), c = 1", RED)]
for x, d, lab, col in pts:
    ax.plot(x, d, "o", color=col, ms=9); ax.text(x + .04, d + .05, lab, color=col, fontsize=8)
ax.plot(X, X - 1, color=RED, ls="--", lw=1.2); ax.text(0.55, -0.35, "repo's line d = k^2/4 - 1:\nall the same geometry", color=RED, fontsize=8)
ax.text(-1.15, 1.0, "kappa > 0: elliptic\n(no walls)", color=BLUE, fontsize=9); ax.text(0.5, -1.3, "kappa < 0: hyperbolic (two walls)", color="#975a16", fontsize=9)
ax.set_ylim(-1.4, 1.2); ax.set_xlabel("k^2/4   (negative for imaginary k)"); ax.set_ylabel("d  (Coleman-McGuinness's second dial; v0.1.0's q)")
ax.set_title("law (a + b + k ab)/(1 - d ab): equal-kappa lines are parallel,\nso k and d are two coordinates on one Cayley-Klein dial (ours: D1)")
save(fig, "f2-k-d-plane-is-one-dial")

# f3. the three one-dial worlds: hyperbola, line, circle, and the lift 1 + k a
fig, ax = plt.subplots(1, 3, figsize=(13, 4.2))
th = np.linspace(-1.4, 1.4, 200)
ax[0].plot(np.cosh(th), np.sinh(th), color=GREEN, lw=2); ax[0].set_title("kappa < 0 (real k): hyperbola\nboosts compose like velocities", color=GREEN)
ax[1].plot(1 + 0*th, th, color=GREY, lw=2); ax[1].set_title("kappa = 0 (k = 0): line\nGalilean, plain addition", color=GREY)
t2 = np.linspace(0, 2*np.pi, 200); ax[2].plot(np.cos(t2), np.sin(t2), color=BLUE, lw=2); ax[2].set_title("kappa > 0 (imaginary k): circle\nrotations, phases", color=BLUE)
for a_, (t1, t2_), f in [(ax[0], (0.5, 0.4), lambda t: (np.cosh(t), np.sinh(t))), (ax[2], (0.6, 0.9), lambda t: (np.cos(t), np.sin(t)))]:
    for t, col in [(t1, ORANGE), (t1 + t2_, RED)]:
        x, y = f(t); a_.plot([0, x], [0, y], color=col, lw=1.5); a_.plot(x, y, "o", color=col)
    a_.text(0.05, -0.9 if a_ is ax[2] else -1.6, "orange: theta1; red: theta1 + theta2\nlog(1 + k a)/k = theta (K3)", fontsize=8)
ax[1].plot([1, 1], [0, 0.5], color=ORANGE, lw=3); ax[1].plot([1, 1], [0.5, 0.9], color=RED, lw=3)
for a_ in ax: a_.set_aspect("equal"); a_.axhline(0, color="#e2e8f0"); a_.axvline(0, color="#e2e8f0")
ax[0].set_xlim(-0.2, 2.4); ax[1].set_xlim(0, 2); ax[2].set_xlim(-1.3, 1.3)
save(fig, "f3-three-one-dial-worlds")

# f4. higher dimensions: one label per pair, but N free dials
fig, ax = plt.subplots(1, 2, figsize=(12, 4.6))
a0 = ax[0]; a0.axis("off"); a0.set_xlim(-0.5, 3.5); a0.set_ylim(-0.6, 2.2)
for i in range(4): a0.plot(i, 0, "o", ms=16, color=GREY); a0.text(i, -0.35, f"axis {i}", ha="center", fontsize=8)
for i, lab in enumerate(["kappa1", "kappa2", "kappa3"]):
    a0.plot([i, i + 1], [0, 0], color=RED, lw=3); a0.text(i + .5, 0.1, lab, ha="center", color=RED, fontsize=10)
from matplotlib.patches import Arc
for (a, b), lab, h in [((0, 2), "kappa1 kappa2", .9), ((1, 3), "kappa2 kappa3", .9), ((0, 3), "kappa1 kappa2 kappa3", 1.7)]:
    a0.add_patch(Arc(((a + b)/2, 0), b - a, h*2, theta1=0, theta2=180, color=BLUE))
    a0.text((a + b)/2, h + .05, lab, ha="center", color=BLUE, fontsize=8)
a0.set_title("a 3-dimensional geometry (axes 0 to 3): 6 pairs, 6 labels,\nbut each label is a product of 3 dials (3+1 spacetime: 10 labels, 4 dials)")
import importlib.util, sys
spec = importlib.util.spec_from_file_location("ck", os.path.join(HERE, "checks.py"))
from scipy.linalg import expm
def gens(labels, N):
    G = {}
    for (a, b), kab in labels.items():
        M = np.zeros((N + 1, N + 1)); M[b, a] = 1; M[a, b] = -kab; G[(a, b)] = M
    return G
def resid(G):
    basis = np.array([M.ravel() for M in G.values()]).T; worst = 0.0
    for p, q in itertools.combinations(G.values(), 2):
        Cm = (p @ q - q @ p).ravel(); coef, *_ = np.linalg.lstsq(basis, Cm, rcond=None); worst = max(worst, np.abs(basis @ coef - Cm).max())
    return worst
rng = np.random.default_rng(3)
rows = []
for N in [2, 3, 4]:
    kaps = rng.choice([-1.3, -0.4, 0.6, 1.1], N)
    L = {(a, b): float(np.prod(kaps[a:b])) for a in range(N + 1) for b in range(a + 1, N + 1)}
    free = {p: rng.uniform(-1.5, 1.5) for p in L}
    rows.append((N, resid(gens(L, N)), resid(gens(free, N))))
xs = np.arange(len(rows))
ax[1].bar(xs - .18, [max(r[1], 1e-17) for r in rows], .34, color=GREEN, label="labels = products of N dials")
ax[1].bar(xs + .18, [r[2] for r in rows], .34, color=RED, label="a free label for every pair")
ax[1].set_yscale("log"); ax[1].set_xticks(xs); ax[1].set_xticklabels([f"N = {r[0]}" for r in rows])
ax[1].set_ylabel("how far the motions fail to close (0 = a geometry)"); ax[1].legend(fontsize=8, loc="center right")
ax[1].set_title("one free dial per pair doesn't give a geometry (N1)")
save(fig, "f4-higher-dimensions")

# f5. scorecard
rows = [("Is k the Cayley-Klein kappa?", "yes, for one label: kappa = -k^2/4 after a Cayley map [ours]"),
        ("Which kappa?", "whichever plane the chain composes; velocity addition is kappa2 = -1/c^2"),
        ("Our two dial settings", "k = -1 is hyperbolic (boost), k = -i is elliptic (rotation)"),
        ("Is d the second dial kappa2?", "no: (k, d) collapse to one kappa = d - k^2/4 [ours]"),
        ("Can a scalar law carry two dials?", "no: all nine plane groups are non-commutative"),
        ("Higher dimensions", "one label per pair, but only N free dials (products)"),
        ("What does carry every dial", "o_k on matrices: the dials sit in the matrices, not in k")]
fig, ax = plt.subplots(figsize=(11.5, 3.9)); ax.axis("off"); ax.set_xlim(0, 1); ax.set_ylim(0, 1)
for i, (q, a_) in enumerate(rows):
    y = 0.92 - i*0.13; col = GREEN if i in (0, 1, 2, 6) else ORANGE
    ax.add_patch(plt.Rectangle((0.0, y - 0.05), 0.012, 0.1, color=col))
    ax.text(0.025, y, q, va="center", fontsize=10, weight="bold"); ax.text(0.38, y, a_, va="center", fontsize=10, color=col)
ax.set_title("does the weighted circle operation carry both Cayley-Klein dials?")
save(fig, "f5-scorecard")
print("figures written")
