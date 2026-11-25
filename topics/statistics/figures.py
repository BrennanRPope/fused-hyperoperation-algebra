"""Figures for the statistics page. Run: python3 figures.py  (writes figures/*.png)"""
import math, numpy as np, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle, Polygon
from scipy import stats

INK = "#15131a"; MUTE = "#5d5a66"; BLUE = "#1d2bd8"; RED = "#cf3519"; GOLD = "#d99a00"; WASH = "#eef0ff"; GREEN = "#1f8a5b"
plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 11, "axes.edgecolor": INK, "axes.labelcolor": INK,
                     "xtick.color": MUTE, "ytick.color": MUTE, "savefig.dpi": 160})

def save(fig, name):
    fig.savefig(f"figures/{name}", bbox_inches="tight", facecolor="white"); plt.close(fig)

# f1: dependence is the weight. Unit square of outcomes, A and B as regions.
fig, axs = plt.subplots(1, 3, figsize=(12, 4.3))
p, q = 0.45, 0.35
cases = [
    ("Disjoint: k = 0", [(0, 0, p, 1)], [(p, 0, q, 1)], None, f"P(A or B) = p + q = {p+q:.2f}"),
    ("Independent: k = −1", [(0, 0, p, 1)], [(0, 0, 1, q)], (0, 0, p, q), f"P(A or B) = p + q − pq = {p+q-p*q:.4f}"),
    ("B inside A: k = −1/p", [(0, 0, p, 1)], [(0, 0, p, q/p)], (0, 0, p, q/p), f"P(A or B) = p = {p:.2f}"),
]
for ax, (title, A, B, C, txt) in zip(axs, cases):
    ax.add_patch(Rectangle((0, 0), 1, 1, fill=False, lw=1.5, ec=INK))
    for r in A: ax.add_patch(Rectangle(r[:2], r[2], r[3], fc=BLUE, alpha=.22, ec=BLUE, lw=1.2))
    for r in B: ax.add_patch(Rectangle(r[:2], r[2], r[3], fc=GOLD, alpha=.30, ec=GOLD, lw=1.2))
    if C: ax.add_patch(Rectangle(C[:2], C[2], C[3], fc="none", ec=RED, lw=2.4, hatch="///"))
    ax.text(p/2, .93, "A", color=BLUE, ha="center", fontsize=14, weight="bold")
    bx, by = (p + q/2, .5) if title.startswith("Dis") else (0.75, q/2) if title.startswith("Ind") else (p/2, (q/p)/2)
    ax.text(bx, by, "B", color="#8a6100", ha="center", va="center", fontsize=14, weight="bold", bbox=dict(fc="white", ec="none", pad=1.5))
    ax.set_xlim(-.02, 1.02); ax.set_ylim(-.2, 1.02); ax.set_aspect("equal"); ax.axis("off")
    ax.set_title(title, color=INK, fontsize=13)
    ax.text(.5, -.1, txt, ha="center", fontsize=10.5, color=INK)
fig.suptitle("The overlap (red) is the corner cell. Its size, as a multiple of p·q, is the weight with its sign flipped.", y=1.02, fontsize=12, color=MUTE)
save(fig, "f1_dependence_weight.png")

# f2: the generating-function dial, t = 1 + k
ps = [0.2, 0.5, 0.7]
G = lambda t: np.prod([(1 - pp + pp*t) for pp in ps], axis=0)
t = np.linspace(-1.3, 2.3, 400)
fig, ax = plt.subplots(figsize=(9, 5.2))
ax.axhline(0, color=MUTE, lw=.8); ax.axvline(0, color=MUTE, lw=.6, ls=":")
ax.plot(t, G(t), color=BLUE, lw=2.6)
marks = [(-1, RED, "t = −1, k = −2\nP(even) − P(odd)\nparity, piling-up"),
         (0, GOLD, "t = 0, k = −1\nP(no heads)\nthe union rule's survival"),
         (1, GREEN, "t = 1, k = 0\ntotal = 1, slope = mean\naddition"),
         (2, INK, "t = 2, k = 1\nE[2^X]\nplain the circle operation")]
for tt, c, lab in marks:
    y = float(G(np.array(tt)))
    ax.plot([tt], [y], "o", color=c, ms=9, zorder=5)
    off = {-1: (0.12, .9), 0: (0.12, 1.1), 1: (0.3, -0.85), 2: (-0.85, 0.55)}[tt]
    ax.annotate(lab, (tt, y), (tt + off[0], y + off[1]), fontsize=9.5, color=c,
                arrowprops=dict(arrowstyle="-", color=c, lw=.8))
tt = np.linspace(0.6, 1.4, 2); ax.plot(tt, 1 + sum(ps)*(tt - 1), color=GREEN, lw=1.2, ls="--")
ax.set_xlabel("t   (the weight is k = t − 1)"); ax.set_ylabel("G(t) = E[t^X]")
ax.set_title("Three coins (0.2, 0.5, 0.7): one curve, every weight is a point on it", color=INK)
ax.set_ylim(-0.6, 4.6); ax.set_xlim(-1.35, 2.35)
for s in ["top", "right"]: ax.spines[s].set_visible(False)
save(fig, "f2_generating_dial.png")

# f3: binomial as Pascal layers of the cube
fig = plt.figure(figsize=(11, 4.8))
ax = fig.add_axes([0.0, 0.05, 0.5, 0.9]); ax.set_aspect("equal"); ax.axis("off")
P = lambda x, y, z: (x + 0.45*z, y + 0.32*z)
corners = [(a, b, c) for a in (0, 1) for b in (0, 1) for c in (0, 1)]
cols = {0: MUTE, 1: BLUE, 2: GOLD, 3: RED}
for u in corners:
    for v in corners:
        if sum(abs(i - j) for i, j in zip(u, v)) == 1 and u < v:
            (x1, y1), (x2, y2) = P(*u), P(*v); ax.plot([x1, x2], [y1, y2], color="#b9b7c4", lw=1.4, zorder=1)
for u in corners:
    x, y = P(*u); j = sum(u)
    lab = {0: "q³", 1: "pq²", 2: "p²q", 3: "p³"}[j]
    ax.plot(x, y, "o", ms=26, color=cols[j], alpha=.9, zorder=3)
    ax.text(x, y, lab, ha="center", va="center", color="white", fontsize=10, weight="bold", zorder=4)
ax.text(0.72, -0.28, "corner = which of 3 coins came up heads\nweight = p per head, q per tail", ha="center", color=MUTE, fontsize=10)
ax.set_xlim(-.3, 1.75); ax.set_ylim(-.45, 1.55)
ax2 = fig.add_axes([0.58, 0.16, 0.38, 0.72])
pp = 0.3; qq = 1 - pp
vals = [qq**3, 3*pp*qq**2, 3*pp**2*qq, pp**3]
ax2.bar(range(4), vals, color=[cols[j] for j in range(4)], width=.62)
for j, v in enumerate(vals):
    ax2.text(j, v + .015, f"{[1,3,3,1][j]} × {['q³','pq²','p²q','p³'][j]}\n= {v:.3f}", ha="center", fontsize=9.5, color=INK)
ax2.set_xticks(range(4)); ax2.set_xticklabels(["0 heads", "1", "2", "3"]); ax2.set_ylim(0, .56)
ax2.set_title("Layer sums: the binomial, p = 0.3", color=INK, fontsize=12)
for s in ["top", "right"]: ax2.spines[s].set_visible(False)
save(fig, "f3_binomial_layers.png")

# f4: CLT, cube layers standardized for growing n
fig, ax = plt.subplots(figsize=(9, 5))
x = np.linspace(-4, 4, 400); ax.plot(x, stats.norm.pdf(x), color=INK, lw=2.2, label="normal")
for n, c in [(1, RED), (4, GOLD), (16, GREEN), (64, BLUE)]:
    pp = 0.3; mu = n*pp; sd = math.sqrt(n*pp*(1 - pp)); k = np.arange(n + 1)
    z = (k - mu)/sd; ax.plot(z, stats.binom.pmf(k, n, pp)*sd, "o-", color=c, ms=4 if n > 10 else 7, lw=1.2,
                             label=f"{n} coin{'s' if n > 1 else ''}: kurtosis {(1 - 6*pp*(1 - pp))/(n*pp*(1 - pp)):+.3f}")
ax.set_xlim(-3.5, 3.8); ax.set_ylim(0, .5)
ax.set_xlabel("standardized number of heads"); ax.set_ylabel("height × spread")
ax.set_title("The cube's layers settle onto the bell; the shape numbers fall as 1/n", color=INK)
ax.legend(frameon=False, fontsize=9.5)
for s in ["top", "right"]: ax.spines[s].set_visible(False)
save(fig, "f4_clt_layers.png")

# f5: the landscape plane with statistics placed
fig, ax = plt.subplots(figsize=(9.5, 6.2))
K = np.linspace(-3, 1.8, 400)
ax.fill_between(K, K**2/4, 2.0, color="#f3eefe", zorder=0); ax.fill_between(K, -1.6, K**2/4, color="#eef6f1", zorder=0)
ax.plot(K, K**2/4, color=MUTE, lw=1.4)
ax.axhline(0, color=INK, lw=1.4)
th = np.linspace(-1, 1, 50); ax.plot(-(1 + th), th, color=RED, lw=3, zorder=3)
pts = [(1, 0, BLUE, "k = 1: generating functions − 1,\nodds − 1 (Bayes), moments", (0.2, 0.32)),
       (0, 0, GREEN, "k = 0: means, variances, cumulants,\nentropy, disjoint events", (-0.2, -0.62)),
       (-1, 0, GOLD, "k = −1: independent union,\nsurvival multiplies", (-1.55, -0.55)),
       (-2, 0, INK, "k = −2: parity, XOR,\npiling-up", (-2.95, -0.42)),
       (0, -1, RED, "velocity addition: pooling\nopinions, ρ = 2p − 1", (0.1, -1.35)),
       (-2, 1, RED, "cusp: double wall\nat certainty", (-2.95, 1.25))]
for k, qv, c, lab, (tx, ty) in pts:
    ax.plot(k, qv, "o", color=c, ms=10, zorder=5, mec="white")
    ax.annotate(lab, (k, qv), (tx, ty), fontsize=9.3, color=c, arrowprops=dict(arrowstyle="-", color=c, lw=.7))
ax.text(-1.5, 0.82, "dependence line k + q = −1\n(AMH copulas, θ from −1 to 1)", color=RED, fontsize=9.5, rotation=-28, ha="center")
ax.text(1.1, 1.7, "no real wall\n(turning)", color="#6b4fb3", fontsize=9.5, ha="center")
ax.text(1.25, -1.15, "two real walls\n(saturating)", color=GREEN, fontsize=9.5, ha="center")
ax.set_xlim(-3, 1.8); ax.set_ylim(-1.6, 2.0)
ax.set_xlabel("k: how full the corner cell is"); ax.set_ylabel("q: how much the next layer bends")
ax.set_title("Statistics on the landscape plane: the laws (a + b + k·ab)/(1 − q·ab)", color=INK)
for s in ["top", "right"]: ax.spines[s].set_visible(False)
save(fig, "f5_landscape_statistics.png")
print("figures written")
