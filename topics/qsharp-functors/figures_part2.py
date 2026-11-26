"""DRAFT figures for sections 8 (Trotter terms) and 9 (ZX-calculus). Orange = known, blue = ours, red = does not hold."""
import itertools, math, sys
import numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Circle, FancyArrowPatch, Rectangle, FancyBboxPatch
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
OR, BL, RD, GR, INK = "#e07b00", "#1f5fbf", "#c62828", "#8a8f99", "#0b1f4b"
plt.rcParams.update({"font.size": 12, "axes.titlesize": 13, "axes.titleweight": "bold", "figure.dpi": 110})
def save(fig, name): fig.savefig(f"figures/{name}.png", bbox_inches="tight", facecolor="white"); plt.close(fig)
def arrow(ax, a, b, c=INK, lw=1.8, cs="arc3,rad=0"):
    ax.add_patch(FancyArrowPatch(a, b, arrowstyle="-|>", mutation_scale=14, color=c, lw=lw, connectionstyle=cs))

# F10: the 3-cube cut into slabs, weights by descents
fig = plt.figure(figsize=(16, 6))
ax = fig.add_subplot(1, 3, 1, projection="3d")
def slab_poly(s):  # polygon x+y+z = s inside unit cube
    pts = []
    for (i, j) in itertools.combinations(range(3), 2):
        pass
    cand = []
    for a in (0, 1):
        for b in (0, 1):
            for free in range(3):
                others = [q for q in range(3) if q != free]
                v = [0, 0, 0]; v[others[0]] = a; v[others[1]] = b; v[free] = s - a - b
                if 0 <= v[free] <= 1: cand.append(tuple(v))
    cand = np.unique(np.round(cand, 9), axis=0)
    c = cand.mean(0); n = np.ones(3)/np.sqrt(3); u = np.cross(n, [1, 0, 0]); u /= np.linalg.norm(u); w = np.cross(n, u)
    ang = np.arctan2((cand-c)@w, (cand-c)@u); return cand[np.argsort(ang)]
for s, col in ((1, BL), (2, BL)):
    ax.add_collection3d(Poly3DCollection([slab_poly(s)], alpha=0.35, facecolor=col, edgecolor=INK))
for a, b in itertools.combinations(list(itertools.product((0, 1), repeat=3)), 2):
    if sum(abs(np.array(a)-np.array(b))) == 1: ax.plot(*zip(a, b), color=GR, lw=1)
ax.text(0, 0, 0, "  000", fontsize=10); ax.text(1, 1, 1, "  111", fontsize=10)
ax.text2D(0.02, -0.02, "slab 0 (corner 000): volume 1/6, weight +1/3\nslab 1 (middle): volume 4/6, weight −1/6\nslab 2 (corner 111): volume 1/6, weight +1/3", transform=ax.transAxes, fontsize=10.5)
ax.set_axis_off(); ax.view_init(28, -35)
ax.set_title("A. The 3-cube cut across its long diagonal\n(planes x+y+z = 1 and 2)")
ax = fig.add_subplot(1, 3, 2)
perms = list(itertools.permutations((1, 2, 3)))
def des(p): return sum(p[i] > p[i+1] for i in range(len(p)-1))
w = {0: "+1/3", 1: "−1/6", 2: "+1/3"}
for i, p in enumerate(sorted(perms, key=des)):
    d = des(p); y = 5 - i
    ax.add_patch(Rectangle((0, y-0.4), 2.2, 0.8, color=OR if d != 1 else BL, alpha=0.25))
    ax.text(0.15, y, "".join("x" + "₀₁₂₃"[q] for q in p), va="center", fontsize=14)
    ax.text(1.35, y, f"{d} descent{'s' if d != 1 else ''}", va="center", fontsize=11)
    ax.text(2.4, y, w[d], va="center", fontsize=14, weight="bold")
ax.text(1.2, -0.9, "the six orderings of three gates:\nweights add to 1/3 − 4/6 + 1/3 = 0", ha="center", fontsize=11)
ax.set_xlim(-0.2, 3.2); ax.set_ylim(-1.5, 5.8); ax.axis("off")
ax.set_title("B. The leftover of x₁∘x₂∘x₃ (top corner):\none weight per ordering, set by its descents")
ax = fig.add_subplot(1, 3, 3)
for row, n in enumerate((3, 4, 5, 6)):
    cnt = [sum(1 for s in itertools.permutations(range(n)) if des(s) == d) for d in range(n)]
    wt = [(-1)**d/(n*math.comb(n-1, d)) for d in range(n)]
    tot = [c*x for c, x in zip(cnt, wt)]
    x0 = 0
    for d, v in enumerate(tot):
        ax.barh(row, v, left=0 if v > 0 else 0, color=OR if v > 0 else BL, alpha=0.5 + 0.1*d, height=0.6)
    ax.text(1.75, row, f"n = {n}: counts {cnt}", va="center", fontsize=10)
pos = [sum(c*x for c, x in zip([sum(1 for s in itertools.permutations(range(n)) if des(s) == d) for d in range(n)], [(-1)**d/(n*math.comb(n-1, d)) for d in range(n)]) if c*x > 0) for n in (3, 4, 5, 6)]
ax.cla()
for row, n in enumerate((3, 4, 5, 6)):
    cnt = [sum(1 for s in itertools.permutations(range(n)) if des(s) == d) for d in range(n)]
    wt = [(-1)**d/(n*math.comb(n-1, d)) for d in range(n)]
    p = sum(c*x for c, x in zip(cnt, wt) if x > 0); m = sum(c*x for c, x in zip(cnt, wt) if x < 0)
    ax.barh(row, p, color=OR, height=0.55); ax.barh(row, m, color=BL, height=0.55)
    ax.text(7.8, row, f"n = {n}: orderings per slab\n{cnt}", va="center", fontsize=10)
ax.axvline(0, color=INK)
ax.set_yticks([]); ax.set_xlim(-7.8, 15.5); ax.set_xlabel("total weight: positive slabs (orange) vs negative (blue)")
ax.set_title("C. They cancel at every level, so gates\nthat commute leave no error")
save(fig, "fig10-trotter-weights-live-on-the-cubes-slabs")

# F11: the k^2, k^3, k^4 terms + numerical convergence
sys.path.insert(0, "trotter")
fig, axs = plt.subplots(1, 2, figsize=(15, 5.6))
ax = axs[0]
rows = [("k²", "½ [A,B]", "2 orderings"), ("k³", "1/12 ([A,[A,B]] + [B,[B,A]])", "6 orderings"),
        ("k⁴", "−1/24 [B,[A,[A,B]]]", "24 orderings"), ("k⁵", "single-B part −1/720 ad_A⁴B …", "120 orderings")]
for i, (o, t, c) in enumerate(rows):
    y = 3 - i
    ax.text(0, y, o, fontsize=18, weight="bold", va="center", color=OR)
    ax.text(0.7, y, t, fontsize=14, va="center"); ax.text(4.9, y, c, fontsize=11, va="center", color=BL)
ax.text(0, -1.0, "each line = the corner weights of section 8 summed over\nthe time slices of 'A then B' (exact, word for word)", fontsize=11)
ax.set_xlim(-0.2, 6.2); ax.set_ylim(-1.6, 3.6); ax.axis("off")
ax.set_title("A. The Trotter error terms, read off the cube's orderings")
ax = axs[1]
def qmul(a, b):
    w1, x1, y1, z1 = a; w2, x2, y2, z2 = b
    return np.array([w1*w2-x1*x2-y1*y2-z1*z2, w1*x2+x1*w2+y1*z2-z1*y2, w1*y2-x1*z2+y1*w2+z1*x2, w1*z2+x1*y2-y1*x2+z1*w2])
def qexp(h):
    t = np.linalg.norm(h[1:]); return np.r_[math.cos(t), h[1:]/t*math.sin(t)]
Aq = np.array([0, 0.9, -0.3, 0.4]); Bq = np.array([0, 0.2, 0.7, -0.5])
ks = np.geomspace(0.01, 1, 30)
e1 = [np.linalg.norm(qmul(qexp(k*Aq), qexp(k*Bq)) - qexp(k*(Aq+Bq))) for k in ks]
e2 = [np.linalg.norm(qmul(qmul(qexp(k*Aq/2), qexp(k*Bq)), qexp(k*Aq/2)) - qexp(k*(Aq+Bq))) for k in ks]
ax.loglog(ks, e1, color=OR, lw=2.5, label="A then B: error ∝ k² (½[A,B])")
ax.loglog(ks, e2, color=BL, lw=2.5, label="½A, B, ½A: error ∝ k³ (Strang)")
ax.loglog(ks, 0.5*ks**2*np.linalg.norm(qmul(Aq, Bq)-qmul(Bq, Aq)), "--", color=OR, lw=1)
ax.set_xlabel("step size k"); ax.set_ylabel("distance from the true gate"); ax.legend(fontsize=10)
ax.set_title("B. Same weights, three time blocks: the k² term cancels")
save(fig, "fig11-trotter-terms-from-the-corner-weights")

# F12: rungs no; Bernoulli = reciprocal of simplex volumes
fig, axs = plt.subplots(1, 2, figsize=(15, 5.2))
ax = axs[0]
n = np.arange(7)
simp = [1/math.factorial(i+1) for i in n]
bern = [1, -1/2, 1/12, 0, -1/720, 0, 1/30240]
ax.bar(n-0.18, simp, width=0.36, color=OR, label="corner-simplex volumes 1/(n+1)!  → e")
ax.bar(n+0.18, np.abs(bern), width=0.36, color=BL, label="|Bernoulli weights| B_n/n! in the error")
ax.set_yscale("log"); ax.set_xlabel("n"); ax.legend(fontsize=10)
ax.set_title("A. The error's single-B weights are the reciprocal\nseries of the simplex volumes (product = 1)")
ax = axs[1]
labels = ["rung 1\nadd", "rung 2\nmultiply", "rung 3\npower"]
for i, l in enumerate(labels):
    ax.add_patch(Rectangle((i*1.6, 0), 1.2, 0.8, color=GR, alpha=0.2)); ax.text(i*1.6+0.6, 0.4, l, ha="center", va="center", fontsize=12)
arrow(ax, (1.55, 0.95), (0.65, 0.95), OR, cs="arc3,rad=0.4"); ax.text(1.1, 1.35, "one log: gates → generator", ha="center", color=OR, fontsize=11)
arrow(ax, (2.75, 0.95), (3.65, 0.95), RD, cs="arc3,rad=-0.4"); ax.text(3.2, 1.35, "no further climb", ha="center", color=RD, fontsize=11)
ax.text(1.9, -0.7, "k³, k⁴, … never climb again: every term is a plain\npolynomial in the gates. Not rungs; corner terms.", ha="center", fontsize=11)
ax.set_xlim(-0.3, 4.6); ax.set_ylim(-1.2, 1.8); ax.axis("off")
ax.set_title("B. Rungs: one step down at the start, none after")
save(fig, "fig12-not-rungs-bernoulli-is-the-simplex-reciprocal")

# F13: ZX dictionary
fig, ax = plt.subplots(figsize=(15, 7.4))
rows = [("ZX spider fusion: phases add", "our core law on the circle (lifts multiply)", "same", OR),
        ("ZX adjoint: flip, negate phases", "our ∘ₖ inverse = mirror + order reversal", "same", OR),
        ("X-spider = parity (XOR)", "our ∘₋₂ on projectors (logic-sets page)", "same", OR),
        ("phase gadget / phase polynomial", "the cube's corner expansion ∏(1 − Zᵢ)/2", "same", OR),
        ("phase folding (gadgets on one parity fuse)", "one corner, its weights add", "same", OR),
        ("Euler-angle rule (merge Z-X-Z runs)", "one ∘ₖ product on S³", "ours simpler to state", BL),
        ("AND / OR / XOR controls", "one dial: product, k = −1, k = −2", "ours unifies", BL),
        ("no error terms for non-commuting gadgets", "corner weights give them (= BCH)", "ours adds (known math)", BL),
        ("bialgebra, local complementation, pivoting", "nothing", "ZX only", RD)]
for i, (zx, ours, verdict, c) in enumerate(rows):
    y = len(rows) - i
    ax.add_patch(Rectangle((0, y-0.42), 0.15, 0.84, color=c))
    ax.text(0.3, y, zx, va="center", fontsize=12); ax.text(6.1, y, ours, va="center", fontsize=12); ax.text(11.6, y, verdict, va="center", fontsize=12, weight="bold", color=c)
for x, h in ((0.3, "ZX-calculus"), (6.1, "Hypercube Composition Theory"), (11.6, "verdict")):
    ax.text(x, len(rows)+0.9, h, fontsize=12, weight="bold")
ax.set_xlim(0, 14.6); ax.set_ylim(0.3, len(rows)+1.4); ax.axis("off")
ax.set_title("Our framework next to ZX: orange = the same idea, blue = ours says more, red = only ZX has it", loc="left")
save(fig, "fig13-zx-dictionary")

# F14: before / after circuits with T-counts
fig, axs = plt.subplots(1, 3, figsize=(17, 6))
ax = axs[0]
qs = ["a", "b", "c", "d"]
def wires(ax, n, L):
    for i in range(n): ax.plot([0, L], [n-1-i]*2, color=INK, lw=1); ax.text(-0.3, n-1-i, qs[i], va="center", fontsize=12)
def ccz(ax, x, qubits, n, col=INK):
    ys = [n-1-q for q in qubits]; ax.plot([x, x], [min(ys), max(ys)], color=col, lw=2)
    for y in ys: ax.plot(x, y, "o", color=col, ms=10)
wires(ax, 4, 4); ccz(ax, 1.2, [0, 1, 2], 4); ccz(ax, 2.6, [0, 1, 3], 4)
ax.text(2, -1.1, "before: two CCZ, T-count 7 + 7 = 14", ha="center", fontsize=11)
ax.text(2, -1.75, "after: corners a, b, ab are shared with the same sign,\ntheir π/8 phases add to Clifford π/4: T-count 8", ha="center", fontsize=11, color=OR)
ax.set_xlim(-0.6, 4.4); ax.set_ylim(-2.3, 3.6); ax.axis("off"); ax.set_title("A. Phase folding: 14 → 8 T gates\n(known; ZX and ours agree)")
ax = axs[1]
def box(ax, x, y, t, c=BL, w=0.55):
    ax.add_patch(FancyBboxPatch((x-w/2, y-0.28), w, 0.56, boxstyle="round,pad=0.02", fc="white", ec=c, lw=1.6, zorder=3)); ax.texts and None; ax.text(x, y, t, ha="center", va="center", fontsize=10, zorder=4)
ax.plot([0, 8], [1, 1], color=INK, lw=1)
for i, g in enumerate(["H", "T", "H", "S", "Rz", "Rx", "T†", "H", "Ry", "S", "H", "T"]): box(ax, 0.4 + i*0.65, 1, g, GR, 0.5)
ax.plot([0, 8], [-0.6, -0.6], color=INK, lw=1)
for i, g in enumerate(["Rz", "Ry", "Rz"]): box(ax, 2.6 + i*1.2, -0.6, g, OR, 0.7)
arrow(ax, (4, 0.6), (4, -0.2), INK)
ax.text(4, -1.5, "12 gates → 3 rotations: one quaternion product on S³.\nStandard optimizers do this already; ZX needs its Euler rule.", ha="center", fontsize=11)
ax.set_xlim(-0.3, 8.3); ax.set_ylim(-2.2, 1.8); ax.axis("off"); ax.set_title("B. A single-qubit run: 12 → 3\n(ours simpler to state, not cheaper)")
ax = axs[2]
ax.plot([0, 6], [1, 1], color=INK); ax.plot([0, 6], [0, 0], color=INK)
for x, g in ((1, "H"), (3, "S"), (5, "H")):
    ax.plot([x, x], [0, 1], color=INK); ax.plot(x, 1, "o", color=INK, ms=8); box(ax, x, 0, g, GR, 0.6)
ax.text(3, 1.5, "C(H) C(S) C(H): T-count 2 + 3 + 2 = 7", ha="center", fontsize=11)
ax.plot([0, 6], [-1.6, -1.6], color=INK); ax.plot([0, 6], [-2.6, -2.6], color=INK)
box(ax, 1, -2.6, "H", OR, 0.6); box(ax, 5, -2.6, "H", OR, 0.6)
ax.plot([3, 3], [-2.6, -1.6], color=INK); ax.plot(3, -1.6, "o", color=INK, ms=8); box(ax, 3, -2.6, "S", OR, 0.6)
ax.text(3, -3.4, "H C(S) H: T-count 3 (the corner commutes with conjugation;\nknown: Q#'s within/apply)", ha="center", fontsize=11)
ax.set_xlim(-0.3, 6.3); ax.set_ylim(-4.0, 2.0); ax.axis("off"); ax.set_title("C. Controlled conjugation: 7 → 3 T gates\n(known; ZX and ours agree)")
save(fig, "fig14-zx-circuits-before-and-after")
print("ok")
