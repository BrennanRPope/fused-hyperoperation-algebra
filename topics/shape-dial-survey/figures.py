"""Pictures for the shape-dial survey. Run after checks.py and laws.py."""
import json, itertools
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
from scipy.spatial import ConvexHull
from scipy import special

R = json.load(open("results.json")); LAWS = json.load(open("laws.json"))
BLUE, ORANGE, GREY, RED, GREEN = "#2b6cb0", "#dd6b20", "#888888", "#c53030", "#2f855a"
plt.rcParams.update({"font.size": 11, "axes.spines.top": False, "axes.spines.right": False})
def save(fig, name):
    fig.savefig(f"figures/{name}.png", dpi=130, bbox_inches="tight"); plt.close(fig); print("wrote", name)
def pball(p, m=400):
    t = np.linspace(0, 2*np.pi, m); D = np.column_stack([np.cos(t), np.sin(t)])
    g = np.max(np.abs(D), 1) if np.isinf(p) else np.sum(np.abs(D)**p, 1)**(1/p)
    return D/g[:, None]
TRI = np.array([[np.cos(a), np.sin(a)] for a in np.pi/2 + 2*np.pi*np.arange(3)/3])

# ---------------------------------------------------------------- fig0: the shape dials
fig, ax = plt.subplots(3, 5, figsize=(14, 8.6))
for j, p in enumerate([1, 1.5, 2, 4, np.inf]):
    B = pball(p); a = ax[0, j]; a.fill(B[:, 0], B[:, 1], color=BLUE, alpha=.25); a.plot(B[:, 0], B[:, 1], color=BLUE)
    a.set_title(["p = 1\ncross-polytope", "p = 1.5", "p = 2\nsphere", "p = 4", "p = ∞\ncube"][j]); a.set_aspect("equal"); a.axis("off")
for j, t in enumerate([0, .25, .5, .75, 1]):
    th = np.linspace(0, 2*np.pi, 600); D = np.column_stack([np.cos(th), np.sin(th)])
    nrm = np.array([[np.cos(a), np.sin(a)] for a in -np.pi/2 + 2*np.pi*np.arange(3)/3])
    g = (1 - t)*np.max(D @ nrm.T, 1)/0.5 + t*1.0
    B = D/g[:, None]; a = ax[1, j]; a.fill(B[:, 0], B[:, 1], color=ORANGE, alpha=.25); a.plot(B[:, 0], B[:, 1], color=ORANGE)
    a.set_title(["triangle\n(simplex)", "", "halfway", "", "round"][j]); a.set_aspect("equal"); a.axis("off")
shapes = [("simplex", [[0, 0], [1, 0], [0, 1]]), ("4 corners", [[0, 0], [1, 0], [.45, 1], [0, 1]]),
          ("hexagon\n(permutohedron)", [[.3, 0], [1, 0], [1.3, .5], [1, 1], [.3, 1], [0, .5]]),
          ("pentagon\n(associahedron)", [[0, 0], [1, 0], [1, .5], [.5, 1], [0, 1]]), ("square\n(cube)", [[0, 0], [1, 0], [1, 1], [0, 1]])]
for j, (nm, P) in enumerate(shapes):
    P = np.array(P, float); a = ax[2, j]; a.add_patch(Polygon(P, color=GREEN, alpha=.25)); a.add_patch(Polygon(P, fill=False, color=GREEN))
    a.set_xlim(-.1, 1.4); a.set_ylim(-.1, 1.1); a.set_aspect("equal"); a.axis("off"); a.set_title(nm)
for i, lab in enumerate(["Dial A: roundness p\n(norm balls)", "Dial B: round ↔ triangle\n(the simplex thread's route)", "Dial C: corner weights\n(generalized permutohedra)"]):
    fig.text(0.005, 0.80 - 0.31*i, lab, fontsize=12, weight="bold", va="center")
fig.suptitle("Three ways to slide the shape while the plane dial (k, d) stays put", fontsize=14, x=0.55)
fig.subplots_adjust(left=0.17); save(fig, "fig0-three-shape-dials")

# ---------------------------------------------------------------- fig1: relativity <-> Bayes (homogeneity + causal order)
fig, ax = plt.subplots(1, 3, figsize=(17, 5), gridspec_kw={"width_ratios": [1, 1.5, 1]})
h = R["homog"]; names = ["triangle", "p=1.5", "p=2", "p=3", "p=4", "square"]
lab = ["triangle", "p = 1.5", "disk", "p = 3", "p = 4", "square"]
vals = [max(h[n], 1e-16) for n in names]
cols = [BLUE if v < 1e-6 else GREY for v in vals]
ax[0].bar(lab, vals, color=cols); ax[0].set_yscale("log"); ax[0].set_ylim(1e-17, 1)
ax[0].set_title("Can a projective move carry the centre\nanywhere and keep the shape?"); ax[0].set_ylabel("best boundary mismatch (0 = yes)")
for x_, v in enumerate(vals): ax[0].text(x_, v*2, "yes" if v < 1e-6 else "no", ha="center", fontsize=10)
ax[0].text(0.5, -0.27, "Only the triangle (Bayes moves) and the disk (Lorentz boosts) can.\nEverything between has no addition law.", transform=ax[0].transAxes, ha="center", fontsize=10)
img = plt.imread("../simplex/figures/fig10-shape-dial-turn-maps.png")
ax[1].imshow(img); ax[1].axis("off"); ax[1].set_title("From the simplex thread (section 7): the leftover turn,\nplane dial fixed, shape sliding round → triangle")
td = np.array(R["ordering_tri_disk"]); ax[2].plot(td[:, 0], td[:, 1], "o-", color=BLUE, label="ours: light cone with this section")
ax[2].axhline(0.25, color=ORANGE, ls="--"); ax[2].axhline(0.2286, color=ORANGE, ls=":")
ax[2].text(0.02, 0.2505, "cube's product order 1/4 (known)", color=ORANGE, fontsize=9, va="bottom")
ax[2].text(0.35, 0.2275, "round light cone 0.229 (known)", color=ORANGE, fontsize=9, va="top")
ax[2].set_xticks([0, .5, 1]); ax[2].set_xticklabels(["triangle", "halfway", "round"]); ax[2].set_ylim(0.222, 0.254)
ax[2].set_title("Causal order in 2+1 dimensions:\nchance two random events are related"); ax[2].legend(fontsize=9, loc="center left")
fig.suptitle("C1 · Relativity (ball) ↔ Bayes (triangle) ↔ the cube's causal order", fontsize=14, y=1.04); save(fig, "fig1-relativity-bayes-causal-order")

# ---------------------------------------------------------------- fig2: Bell bound
fig, ax = plt.subplots(1, 3, figsize=(15, 4.4), gridspec_kw={"width_ratios": [1.4, 1, 1]})
pp = np.concatenate([np.linspace(1, 2, 60), np.linspace(2, 12, 120)])
ax[0].plot(pp, 2**(1 + np.maximum(1/pp, 1 - 1/pp)), color=BLUE, lw=2.5, label="ours: 2^(1 + max(1/p, 1 − 1/p))")
for k_, v in R["chsh_p"].items():
    x_ = 13 if k_ == "inf" else float(k_); ax[0].plot(x_, v, "o", color="k", ms=5)
ax[0].plot(13, 4, "o", color="k"); ax[0].axhline(2, color=GREY, ls="--"); ax[0].axhline(2*np.sqrt(2), color=ORANGE, ls=":"); ax[0].axhline(4, color=ORANGE, ls=":")
ax[0].text(1.05, 2.05, "classical bit (segment): 2", color=GREY, fontsize=9); ax[0].text(6, 2.86, "Tsirelson 2√2 (known)", color=ORANGE, fontsize=9)
ax[0].text(6, 3.93, "PR box 4 (known)", color=ORANGE, fontsize=9, va="top")
ax[0].set_xticks([1, 2, 4, 8, 13]); ax[0].set_xticklabels(["1\n(cross)", "\n2\n(sphere)", "4", "8", "∞\n(cube)"])
ax[0].set_xlabel("roundness p of the shape the observables live in"); ax[0].set_ylabel("largest CHSH value"); ax[0].legend(fontsize=9, loc="lower right")
ax[0].set_title("Bell's bound slides with the shape")
ns = sorted(int(k_) for k_ in R["chsh_poly"]); ax[1].plot(ns, [R["chsh_poly"][str(n)] for n in ns], "o-", color=BLUE)
ax[1].axhline(2*np.sqrt(2), color=ORANGE, ls=":"); ax[1].set_xscale("log"); ax[1].set_xlabel("sides of the polygon"); ax[1].set_title("Even polygons: square 4, hexagon 3;\nmultiples of 8 exactly 2√2")
B = pball(np.inf); ax[2].plot(B[:, 0], B[:, 1], color=BLUE, lw=2); B = pball(2); ax[2].plot(B[:, 0], B[:, 1], color=ORANGE, lw=2)
ax[2].plot([-1, 1], [0, 0], color=GREY, lw=3); ax[2].text(0, 0.06, "segment (bit): 2", ha="center", color=GREY, fontsize=9)
ax[2].text(0, 1.06, "square (cube): 4", ha="center", color=BLUE, fontsize=9); ax[2].text(0, 0.72, "circle (sphere): 2√2", ha="center", color=ORANGE, fontsize=9)
ax[2].set_aspect("equal"); ax[2].axis("off"); ax[2].set_title("same law, three shapes")
fig.suptitle("C2 · Quantum correlations: one bound, read on any shape", fontsize=14, y=1.04); save(fig, "fig2-bell-bound-slides")

# ---------------------------------------------------------------- fig3: Born step and norm-keeping dynamics
fig, ax = plt.subplots(1, 3, figsize=(15, 4.4))
th = np.linspace(0, np.pi/2, 100); ax[0].plot(np.cos(th), np.sin(th), color=ORANGE, lw=2.5, label="amplitudes (p = 2)")
ax[0].plot([1, 0], [0, 1], color=BLUE, lw=2.5, label="probabilities (p = 1)")
for t in np.linspace(0.15, 1.4, 6):
    ax[0].annotate("", xy=(np.cos(t)**2, np.sin(t)**2), xytext=(np.cos(t), np.sin(t)), arrowprops=dict(arrowstyle="->", color=GREY))
ax[0].set_aspect("equal"); ax[0].legend(fontsize=9, loc="upper right"); ax[0].set_title("Born rule = squaring:\nthe step from p = 2 to p = 1")
for d, c, nm in [(R["iso_signed"], ORANGE, "all vectors (signed maps)"), (R["iso_pos"], BLUE, "positive vectors (positive maps)")]:
    ks = sorted(d, key=float); ax[1].plot([float(k_) for k_ in ks], [max(d[k_], 1e-16) for k_ in ks], "o-", color=c, label=nm)
ax[1].set_yscale("log"); ax[1].set_xlabel("p"); ax[1].set_ylabel("best mismatch for a map near the identity")
ax[1].set_title("Is there a smooth family of maps\nthat keeps the p-norm?"); ax[1].legend(fontsize=9, loc="center right")
ax[1].annotate("p = 1: stochastic\n(classical chance)", (1, 1e-16), (0.8, 1e-6), color=BLUE, fontsize=9, arrowprops=dict(arrowstyle="->", color=BLUE))
ax[1].annotate("p = 2: rotations\n(quantum, unitary)", (2, 1e-16), (2.35, 1e-13), color=ORANGE, fontsize=9, arrowprops=dict(arrowstyle="->", color=ORANGE))
ax[2].axis("off")
ax[2].text(0, .95, "What slides and what doesn't", weight="bold", fontsize=12, va="top")
ax[2].text(0, .8, "• The state law slides: squaring carries the\n  sphere of amplitudes onto the simplex.\n\n"
           "• The motion law does not slide, it jumps:\n  a continuous family of norm-keeping moves\n  exists only at p = 1 (Markov chains) and\n  p = 2 (unitary quantum). Every other p\n  has only finitely many.\n\n"
           "Known: Aaronson 2004 ('quantum mechanics\nis the 2-norm'); Banach–Lamperti.", fontsize=10, va="top")
fig.suptitle("C3 · Quantum (sphere) ↔ classical probability (simplex): two stops, nothing between", fontsize=14, y=1.04); save(fig, "fig3-born-step-two-stops")

# ---------------------------------------------------------------- fig4: replicator <-> imaginary-time Schrodinger
rp = R["replicator"]; P = np.array(rp["P"]); fit = np.array(rp["fit"])
fig = plt.figure(figsize=(15, 4.6))
a0 = fig.add_subplot(1, 3, 1)
p3 = P[:, :3]/P[:, :3].sum(1, keepdims=True); C2 = np.array([[0, 0], [1, 0], [.5, np.sqrt(3)/2]])
xy = p3 @ C2; a0.add_patch(Polygon(C2, fill=False, color=GREY)); a0.plot(xy[:, 0], xy[:, 1], color=BLUE, lw=2.5); a0.plot(*xy[0], "o", color=BLUE)
a0.set_aspect("equal"); a0.axis("off"); a0.set_title("Replicator flow on the simplex\n(three of five types shown)")
a1 = fig.add_subplot(1, 3, 2, projection="3d")
u = np.linspace(0, np.pi/2, 30); v = np.linspace(0, np.pi/2, 30); U, Vv = np.meshgrid(u, v)
a1.plot_surface(np.cos(U)*np.sin(Vv), np.sin(U)*np.sin(Vv), np.cos(Vv), color=ORANGE, alpha=.15)
S3 = np.sqrt(p3); a1.plot(S3[:, 0], S3[:, 1], S3[:, 2], color=ORANGE, lw=2.5); a1.set_title("Same flow read as √p on the sphere:\nimaginary-time Schrödinger, H = −fitness/2")
a1.set_xticks([]); a1.set_yticks([]); a1.set_zticks([])
a2 = fig.add_subplot(1, 3, 3); t = np.array(rp["t"]); mean = P @ fit; var = P @ fit**2 - mean**2
a2.plot(t, np.gradient(mean, t), color=BLUE, lw=2.5, label="d(mean fitness)/dt"); a2.plot(t, var, "--", color=ORANGE, lw=2, label="variance of fitness")
a2.set_xlabel("time"); a2.legend(fontsize=9); a2.set_title("Fisher's theorem = 'energy falls by\ntwice its variance' (exact)")
fig.suptitle("C4 · Evolution (simplex) ↔ quantum ground-state seeking (sphere), through the Born step", fontsize=14, y=1.04); save(fig, "fig4-replicator-imaginary-time")

# ---------------------------------------------------------------- fig5: genus slide
def loday(n):
    from functools import lru_cache
    @lru_cache(None)
    def trees(i, j):
        if i > j: return [()]
        return [(L, r, Rr) for r in range(i, j + 1) for L in trees(i, r - 1) for Rr in trees(r + 1, j)]
    def size(t): return 0 if t == () else size(t[0]) + 1 + size(t[2])
    def coords(t, acc):
        if t == (): return
        acc[t[1]] = (size(t[0]) + 1)*(size(t[2]) + 1); coords(t[0], acc); coords(t[2], acc)
    out = []
    for t in trees(0, n - 1):
        acc = [0]*n; coords(t, acc); out.append(acc)
    return np.array(out, float)
Q = np.linalg.qr(np.array([[1, -1, 0, 0], [0, 1, -1, 0], [0, 0, 1, -1]], float).T)[0]
polys = {"simplex": np.array([[0, 0, 0], [1, 0, 0], [0, 1, 0], [0, 0, 1]], float),
         "cube": np.array(list(itertools.product([0, 1], repeat=3)), float),
         "permutohedron": np.array(list(itertools.permutations([1, 2, 3, 4])), float) @ Q,
         "associahedron": (loday(4) - loday(4).mean(0)) @ Q}
poly_txt = {"simplex": "1 + q + q² + q³\n= [4]_q  (q-counting)", "cube": "1 + 3q + 3q² + q³\n= (1 + q)³  (powers)",
            "permutohedron": "1 + 11q + 11q² + q³\n(Eulerian)", "associahedron": "1 + 6q + 6q² + q³\n(Narayana)"}
fig = plt.figure(figsize=(15, 7.6))
for j, (nm, Pt) in enumerate(polys.items()):
    a = fig.add_subplot(2, 4, j + 1, projection="3d"); hull = ConvexHull(Pt); c = Pt.mean(0)
    a.add_collection3d(Poly3DCollection([Pt[s] for s in hull.simplices], facecolor=GREEN, alpha=.25, edgecolor="none"))
    eq = np.round(hull.equations, 6); groups = {}
    for i_, e_ in enumerate(map(tuple, eq)): groups.setdefault(e_, set()).update(hull.simplices[i_])
    Fs = list(groups.values())
    for u_, v_ in itertools.combinations(range(len(Pt)), 2):
        if sum(1 for F_ in Fs if u_ in F_ and v_ in F_) >= 2: a.plot(*np.array([Pt[u_], Pt[v_]]).T, color="#22543d", lw=1)
    a.scatter(*Pt.T, color="k", s=10); lim = np.abs(Pt - c).max()
    for d_ in range(3): getattr(a, "set_" + "xyz"[d_] + "lim")(c[d_] - lim, c[d_] + lim)
    a.set_axis_off(); a.set_title(f"{nm}\n{poly_txt[nm]}", fontsize=10)
g = R["genus"]; keys = list(g); shp = ["simplex", "cube", "permutohedron", "associahedron"]
a = fig.add_subplot(2, 1, 2); w = 0.25
for i, kk in enumerate(keys):
    a.bar(np.arange(4) + (i - 1)*w, [g[kk][s] for s in shp], w, label=kk, color=[GREY, ORANGE, BLUE][i])
    for j, s in enumerate(shp): a.text(j + (i - 1)*w, g[kk][s] + .4, f"{g[kk][s]:.0f}", ha="center", fontsize=9)
a.set_xticks(range(4)); a.set_xticklabels(shp); a.legend(fontsize=9, loc="upper left"); a.set_ylabel("genus value")
a.set_title("Hold the plane dial, slide the shape. At our core law every shape gives 1 (the slide does nothing);\n"
            "at the signature point every shape gives 0 (odd dimension); at the Euler point it counts corners: 4, 8, 24, 14.", fontsize=10)
fig.suptitle("C5 · One formula (Hirzebruch's genus = our plane dial), the shape picks the polynomial", fontsize=14); save(fig, "fig5-genus-slide")

# ---------------------------------------------------------------- fig6: volume constants
fig, ax = plt.subplots(1, 2, figsize=(14, 4.6))
vs = R["vol_sum_log10"]; ps = [float(k_) for k_ in vs]
ax[0].plot(ps, [10**vs[k_] for k_ in vs], "o-", color=BLUE, lw=2); ax[0].set_yscale("log")
ax[0].plot([1, 2], [10**vs["1"], 10**vs["2"]], "o", color=ORANGE, ms=9, zorder=5)
ax[0].text(1.1, 1e30, "orange stops, already known:\n  p = 1 (cross-polytope): e² = 7.39\n  p = 2 (sphere): e^π (1 + erf √π) = 46.0\n\ntoward the cube (p = ∞) the sum diverges", color=ORANGE, fontsize=10, va="bottom")
ax[0].set_xlabel("roundness p"); ax[0].set_ylabel("sum of unit-ball volumes over all dimensions")
ax[0].set_title(r"One Mittag-Leffler function $E_{1/p}(2\Gamma(1 + 1/p))$ along the dial")
t = np.linspace(0, .45, 200)
from scipy.special import erf
ax[1].plot(t, np.exp(t), color=GREY, lw=2, label="corner simplex: eᵗ (exponential)")
ax[1].plot(t, np.exp(2*t), color=GREEN, lw=2, label="cross-polytope: e²ᵗ (exponential)")
ax[1].plot(t, np.exp(np.pi*t**2)*(1 + erf(np.sqrt(np.pi)*t)), color=ORANGE, lw=2, label="sphere: e^{πt²}(1 + erf √π t) (Gaussian)")
ax[1].plot(t, 1/(1 - 2*t), color=BLUE, lw=2, label="cube: 1/(1 − 2t) (geometric series)")
ax[1].set_xlabel("t"); ax[1].set_ylabel("Σ volume(n) tⁿ"); ax[1].legend(fontsize=9); ax[1].set_ylim(0, 12)
ax[1].set_title("The same sum as a function of t:\nexponential → Gaussian → geometric")
fig.suptitle("C6 · e (simplex) and π (sphere) are two readings of one function; the cube's reading is the geometric series", fontsize=13, y=1.04)
save(fig, "fig6-volume-constants-slide")

# ---------------------------------------------------------------- fig7: Gauss-Bonnet
def surf(p, n=120):
    th = np.linspace(1e-3, np.pi - 1e-3, n); ph = np.linspace(0, 2*np.pi, 2*n); TH, PH = np.meshgrid(th, ph, indexing="ij")
    D = np.stack([np.sin(TH)*np.cos(PH), np.sin(TH)*np.sin(PH), np.cos(TH)], -1)
    Pt = D/(np.sum(np.abs(D)**p, -1)**(1/p))[..., None]
    x, y, z = [Pt[..., i] for i in range(3)]
    ax_ = np.abs(Pt) + 1e-12
    # Gaussian curvature of |x|^p+|y|^p+|z|^p=1 (implicit surface formula)
    gx, gy, gz = [np.sign(Pt[..., i])*ax_[..., i]**(p - 1) for i in range(3)]
    hx, hy, hz = [(p - 1)*ax_[..., i]**(p - 2) for i in range(3)]
    K = (hy*hz*gx**2 + hx*hz*gy**2 + hx*hy*gz**2)/(gx**2 + gy**2 + gz**2)**2
    return x, y, z, K
fig = plt.figure(figsize=(16, 7.6))
for j, p in enumerate([1.3, 2, 4, 12]):
    a = fig.add_subplot(2, 4, j + 1, projection="3d"); x, y, z, K = surf(p)
    Kn = np.clip(np.log10(K + 1e-9), -2, 2); colors = plt.cm.viridis((Kn + 2)/4)
    a.plot_surface(x, y, z, facecolors=colors, rstride=2, cstride=2, linewidth=0, antialiased=False, shade=False)
    a.set_box_aspect((1, 1, 1)); a.set_axis_off(); a.set_title(f"p = {p}" + {1.3: "  (toward octahedron)", 2: "  (sphere)", 12: "  (toward cube)"}.get(p, "") + ("\ncolour: curvature, log scale (yellow = high)" if p == 2 else ""))
gb = R["gauss_bonnet"]; ks = sorted(gb, key=float)
a = fig.add_subplot(2, 2, 3)
a.plot([float(k_) for k_ in ks], [100*gb[k_]["curv_frac_near_corners"] for k_ in ks], "o-", color=BLUE, lw=2, label="share of curvature near the corners")
a.plot([float(k_) for k_ in ks], [100*gb[k_]["area_frac_near_corners"] for k_ in ks], "s--", color=GREY, label="share of area near the corners")
a.set_xscale("log"); a.set_xlabel("p"); a.set_ylabel("%"); a.legend(fontsize=9); a.set_title("Total stays 4π at every p (checked); where it sits slides")
a = fig.add_subplot(2, 2, 4); a.axis("off")
a.text(0, .95, "Reading", weight="bold", fontsize=12, va="top")
a.text(0, .82, "Our swap rule: the turn around a loop = curvature × area.\n\n"
       "• Sphere: curvature spread evenly, so every loop turns in\n  proportion to its area (Thomas, Foucault, Berry).\n"
       "• Cube: all curvature sits in 8 corners, π/2 each. A loop turns\n  only if it goes round a corner, and then by exactly π/2.\n"
       "• Octahedron: 6 corners, 2π/3 each.\n\n"
       "So sliding p carries the smooth swap rule into 'curvature by\ncounting corners' (Descartes' defect law, known), with the\ntotal 4π fixed throughout (Gauss–Bonnet, known).", fontsize=10, va="top")
fig.suptitle("C7 · The swap rule (sphere) slides into corner-counting curvature (cube, octahedron)", fontsize=14); save(fig, "fig7-gauss-bonnet-slide")

# ---------------------------------------------------------------- fig8: Fermat on the p dial
from fractions import Fraction
def rat_points(p, Nmax=24):
    pts = set()
    for c in range(1, Nmax + 1):
        for a_ in range(-c, c + 1):
            for b_ in range(-c, c + 1):
                if (p == "inf" and max(abs(a_), abs(b_)) == c) or (p != "inf" and abs(a_)**p + abs(b_)**p == c**p):
                    pts.add((Fraction(a_, c), Fraction(b_, c)))
    return np.array([[float(x), float(y)] for x, y in pts])
fig, ax = plt.subplots(1, 6, figsize=(18, 3.8), gridspec_kw={"width_ratios": [1, 1, 1, 1, 1, 1.7], "wspace": 0.35})
for j, p in enumerate([1, 2, 3, 4, "inf"]):
    B = pball(np.inf if p == "inf" else p); a = ax[j]; a.plot(B[:, 0], B[:, 1], color=GREY)
    pts = rat_points(p); a.plot(pts[:, 0], pts[:, 1], ".", color=BLUE if p in (1, "inf") else (ORANGE if p == 2 else RED), ms=3 if len(pts) > 50 else 7)
    a.set_aspect("equal"); a.axis("off"); a.set_title({1: "p = 1: every\nrational point", 2: "p = 2: Pythagorean\ntriples (dense)", 3: "p = 3: only 4\n(Fermat, our Lean)", 4: "p = 4: only 4\n(Fermat, our Lean)", "inf": "p = ∞: every\nrational point"}[p], fontsize=10)
fc = R["fermat_counts"]; N = np.arange(1, len(fc["1"]) + 1)
for k_, c in [("1", BLUE), ("inf", GREEN), ("2", ORANGE), ("3", RED), ("4", "#7b341e")]:
    ax[5].loglog(N, fc[k_], color=c, lw=2, label=f"p = {'∞' if k_ == 'inf' else k_}")
ax[5].set_xlabel("largest denominator"); ax[5].set_ylabel("points", labelpad=2); ax[5].legend(fontsize=8); ax[5].set_title("counts: N², N, then flat 4")
fig.suptitle("C8 · Fermat's theorem is the roundness dial at whole-number p: all → dense → none → all", fontsize=14, y=1.06); save(fig, "fig8-fermat-on-the-dial")

# ---------------------------------------------------------------- fig9: Lee-Yang
scan = R["lee_yang_scan"]; prow = list(scan); acol = list(next(iter(scan.values())))
M = np.zeros((len(prow), len(acol)))
for i, p in enumerate(prow):
    for j, a_ in enumerate(acol):
        # corrected with exact sums (new_laws.py, N3): on the line iff a <= p, or p an even whole number
        pf, af = float(p), float(a_); M[i, j] = 1 if (af <= pf + 1e-12 or abs(pf/2 - round(pf/2)) < 1e-12) else 0
fig, ax = plt.subplots(1, 2, figsize=(14, 4.8))
cm = matplotlib.colors.ListedColormap([RED, BLUE]); ax[0].imshow(M, cmap=cm, vmin=0, vmax=1, aspect="auto")
for i in range(len(prow)):
    for j in range(len(acol)):
        ax[0].text(j, i, "?" if np.isnan(M[i, j]) else ("on" if M[i, j] else "off"), ha="center", va="center", color="w" if not np.isnan(M[i, j]) else "k", fontsize=8)
ax[0].set_xticks(range(len(acol))); ax[0].set_xticklabels(acol); ax[0].set_yticks(range(len(prow))); ax[0].set_yticklabels(prow)
ax[0].set_xlabel("a = (n − 1)/p   (a spin uniform in an n-dimensional p-ball has a = (n − 1)/p)"); ax[0].set_ylabel("roundness p")
ax[0].set_title("Lee–Yang for one spin read off a p-shape (exact sums):\nblue = every zero on the line, red = some zeros off it")
def Gfun(p, a, s, n=200):
    tj, wj = special.roots_jacobi(n, a, 0); x = (tj + 1)/2; w = wj*.5**(a + 1)
    with np.errstate(divide="ignore", invalid="ignore"): r = np.where(x < 1 - 1e-14, (1 - x**p)/(1 - x), p)
    return np.cos(np.outer(s, x)) @ (r**a*w)
X, Y = np.meshgrid(np.linspace(0.2, 26, 500), np.linspace(-6, 6, 260)); Z = X + 1j*Y
for p, a_, c, lab in [(2, 2, ORANGE, "ball, n = 3, p = 2: zeros on the line"), (1, 2, RED, "octahedron, n = 3, p = 1: zeros off it")]:
    G = Gfun(p, a_, Z.ravel()).reshape(Z.shape)
    ax[1].contour(X, Y, G.real, [0], colors=c, linewidths=1.2); ax[1].contour(X, Y, G.imag, [0], colors=c, linewidths=1.2, linestyles=":")
    ax[1].plot([], [], color=c, label=lab)
ax[1].axhline(0, color=GREY, lw=.8); ax[1].legend(fontsize=9, loc="upper right"); ax[1].set_xlabel("Re (imaginary field)"); ax[1].set_ylabel("Im")
ax[1].set_title("Zeros sit where solid and dotted curves cross")
fig.suptitle("C9 · Lee–Yang survives on the round and cube side, and breaks near the cross-polytope", fontsize=14, y=1.03); save(fig, "fig9-lee-yang-survives-partly")

# ---------------------------------------------------------------- fig10: k dial in disguise
fig, ax = plt.subplots(1, 5, figsize=(17, 3.8), gridspec_kw={"width_ratios": [1, 1, 1, 1, 1.3]})
U1, U2 = np.meshgrid(np.linspace(-4, 4, 400), np.linspace(-4, 4, 400))
for j, kk in enumerate(R["K_tropical_k"]):
    terms = [np.zeros_like(U1), U1, U2] + ([np.log(kk) + U1 + U2] if kk > 0 else [])
    T = np.stack(terms); top = np.sort(T, 0); edge = (top[-1] - top[-2]) < 0.05
    a = ax[j]; a.imshow(edge, extent=[-4, 4, -4, 4], origin="lower", cmap=matplotlib.colors.ListedColormap(["white", BLUE]))
    a.set_title(f"k = {kk}" + ("  simplex\n1 + x + y" if kk == 0 else ("  cube\n(1 + x)(1 + y)" if kk == 1 else "\n1 + x + y + kxy")), fontsize=10)
    a.set_xticks([]); a.set_yticks([])
P3 = np.array([[0, 0], [1, 0], [.5, np.sqrt(3)/2]]); g = np.linspace(0, 1, 300); A, Bb = np.meshgrid(g, g); Cc = 1 - A - Bb
mask = Cc > 0; pts = np.stack([A, Bb, np.where(mask, Cc, np.nan)], -1)
xy = pts @ P3
for q, c in [(0.5, GREEN), (1.0001, GREY), (2, ORANGE), (8, BLUE)]:
    H = np.log(np.nansum(pts**q, -1))/(1 - q)
    ax[4].contour(xy[..., 0], xy[..., 1], np.where(mask, H, np.nan), [0.6], colors=c)
    lab = {0.5: "order 1/2  (k = 1/2)", 1.0001: "order 1, Shannon  (k = 0)", 2: "order 2  (k = −1)", 8: "order 8  (k = −7)"}[q]
    ax[4].plot([], [], color=c, label=lab)
ax[4].add_patch(Polygon(P3, fill=False, color="k")); ax[4].set_aspect("equal"); ax[4].axis("off")
ax[4].legend(fontsize=8, loc="center left", bbox_to_anchor=(1.0, 0.5), frameon=False)
ax[4].set_title("Rényi level sets on the simplex:\nthe norm's p is the order, and k = 1 − p", fontsize=10)
fig.suptitle("C11 · Where the law is a sum over outcomes, sliding the shape and turning k are the same move", fontsize=14, y=1.07); save(fig, "fig10-k-dial-in-disguise")

# ---------------------------------------------------------------- fig11: stick-breaking
fig, ax = plt.subplots(1, 2, figsize=(11, 4.6))
for v in np.linspace(0, 1, 11):
    s = np.linspace(0, 1, 100)
    ax[0].plot(s, v + 0*s, color=GREY, lw=.8); ax[0].plot(v + 0*s, s, color=GREY, lw=.8)
    for (x1, y1) in [(s, v + 0*s), (v + 0*s, s)]:
        p1 = x1; p2 = (1 - x1)*y1; xy = np.column_stack([p1, p2])
        ax[1].plot(xy[:, 0], xy[:, 1], color=BLUE, lw=.8)
ax[0].set_aspect("equal"); ax[0].set_title("cube of cuts (v₁, v₂): independent"); ax[0].set_xlabel("v₁"); ax[0].set_ylabel("v₂")
ax[1].set_aspect("equal"); ax[1].set_title("simplex of shares (p₁, p₂):\np₁ + p₂ = v₁ ∘₋₁ v₂"); ax[1].set_xlabel("p₁"); ax[1].set_ylabel("p₂")
fig.suptitle("C10 · Union rule (cube) → shares of a whole (simplex), at the same k = −1", fontsize=13, y=1.03); save(fig, "fig11-stick-breaking")

# ---------------------------------------------------------------- fig12: sorting board
rows = LAWS["rows"]; order = ["slides", "plausible", "partly", "kdial", "self", "no"]
colr = {"slides": BLUE, "plausible": "#63b3ed", "partly": "#9f7aea", "kdial": ORANGE, "self": GREEN, "no": "#cbd5e0"}
title = {"slides": "Slides into another law", "plausible": "Slides (plausible)", "partly": "Survives partly",
         "kdial": "k dial in disguise", "self": "Is the dial", "no": "No slide"}
fig, ax = plt.subplots(figsize=(17, 10.5)); ax.axis("off"); xcol = 0
cols_layout = [["slides"], ["plausible", "partly", "kdial", "self"], ["no"]]
for ci, group in enumerate(cols_layout):
    y = 1.0
    for vk in group:
        rs = [r for r in rows if r["verdict"] == vk]
        ax.text(ci*0.34, y, f"{title[vk]} ({len(rs)})", weight="bold", fontsize=12, color=colr[vk] if vk != "no" else "#4a5568", transform=ax.transAxes)
        y -= 0.035
        for r in rs:
            tag = f"[{r['candidate']}] " if r["candidate"] else ""
            ax.text(ci*0.34, y, f"{tag}{r['law']}", fontsize=8.6, transform=ax.transAxes,
                    bbox=dict(facecolor=colr[vk], alpha=.25 if vk != "no" else .5, edgecolor="none", pad=1.5))
            ax.text(ci*0.34 + 0.005, y - 0.017, f"lives on: {r['shape']}", fontsize=7, color="#4a5568", transform=ax.transAxes)
            y -= 0.04
        y -= 0.02
fig.suptitle(f"All {len(rows)} proved laws (grouped by Lean module) and what the shape dial does to each", fontsize=15, y=0.99)
save(fig, "fig12-sorting-board")
