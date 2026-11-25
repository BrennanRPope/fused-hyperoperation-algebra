"""Figures for the actuarial page. Run: python3 figures.py  (writes figures/*.png)"""
import os, math, itertools, warnings
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon, FancyArrowPatch
warnings.filterwarnings("ignore")
HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "figures"); os.makedirs(OUT, exist_ok=True)
INK, MUTE, GRID = "#15131a", "#5d5a66", "#e4e3ea"
BLUE, RED, GREEN, AMBER = "#1d2bd8", "#cf3519", "#1b8a5a", "#b07a00"
plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 11, "axes.edgecolor": MUTE,
    "axes.labelcolor": INK, "xtick.color": MUTE, "ytick.color": MUTE, "axes.spines.top": False,
    "axes.spines.right": False, "axes.grid": True, "grid.color": GRID, "grid.linewidth": 0.8,
    "lines.linewidth": 2, "savefig.facecolor": "white", "figure.facecolor": "white"})
rng = np.random.default_rng(7)
def load(name, key):
    import rdata
    return rdata.read_rda(os.path.join(HERE, "data", name))[key]
def save(fig, name):
    fig.savefig(os.path.join(OUT, name), dpi=170, bbox_inches="tight"); plt.close(fig)

def gen(n, k):
    g = 1.0
    for t in range(k):
        g *= (n - t) / (t + 1)
    return g
def coeffs(n, c, K):
    x = c / (1 - c)
    return np.array([(1 - c) ** n * gen(n, k) * x ** k for k in range(K)])

# ---------- a1: the dimension line for claim-count models ----------
fig = plt.figure(figsize=(11, 6.4))
gs = fig.add_gridspec(2, 4, height_ratios=[1.05, 1], hspace=0.45, wspace=0.35)
# top row: cube (binomial m=3), orthant (neg. binomial r=2), forbidden n=2.5, log
ax = fig.add_subplot(gs[0, 0]); ax.set_title("n = 3: corners of a cube\n(binomial, 3 policies)", fontsize=10.5)
P = {}
for c in itertools.product([0, 1], repeat=3):
    P[c] = (c[0] + 0.45 * c[2], c[1] + 0.35 * c[2])
for a_, b_ in itertools.combinations(P, 2):
    if sum(abs(np.array(a_) - np.array(b_))) == 1:
        ax.plot(*zip(P[a_], P[b_]), color=MUTE, lw=1)
cols = {0: INK, 1: BLUE, 2: GREEN, 3: RED}
for c, (x, y) in P.items():
    ax.scatter(x, y, s=70, color=cols[sum(c)], zorder=3, edgecolor="white", linewidth=1.5)
ax.text(0.72, -0.32, "layers 1-3-3-1 = P(0,1,2,3 claims)", ha="center", fontsize=9, color=MUTE)
ax.set_xlim(-0.3, 1.75); ax.set_ylim(-0.45, 1.55); ax.axis("off")
ax = fig.add_subplot(gs[0, 1]); ax.set_title("n = -2: lattice points of an orthant\n(negative binomial, r = 2)", fontsize=10.5)
for i in range(6):
    for j in range(6):
        if i + j <= 6:
            ax.scatter(i, j, s=55 * 0.8 ** (i + j), color=cols.get(min(i + j, 3), RED) if i + j < 4 else MUTE, zorder=3)
for t in range(1, 5):
    ax.plot([0, t], [t, 0], color=GRID, lw=1, zorder=1)
ax.text(3.0, -1.75, "layer k has k + 1 points, with no end:\nthe series blows up at z = 1 + 1/β", fontsize=9, color=MUTE, ha="center")
ax.set_xlim(-0.5, 6.5); ax.set_ylim(-2.4, 6.5); ax.axis("off")
ax = fig.add_subplot(gs[0, 2]); ax.set_title("n = 2.5: no such model\n(a negative probability)", fontsize=10.5)
p = coeffs(2.5, 0.3, 7)
ax.bar(range(7), p, color=[BLUE if v >= 0 else RED for v in p], width=0.7)
ax.axhline(0, color=INK, lw=1)
ax.set_xticks(range(7)); ax.set_xlabel("claims k", fontsize=9)
ax.annotate("p₄ < 0", (4, p[4]), (4.6, 0.2), color=RED, fontsize=10, arrowprops=dict(arrowstyle="-", color=RED))
ax = fig.add_subplot(gs[0, 3]); ax.set_title("n = 0: the circle-operation log\n(logarithmic, zero removed)", fontsize=10.5)
b = 3.0; k = np.arange(1, 9); pl = (b / (1 + b)) ** k / k / np.log(1 + b)
ax.bar(k, pl, color=GREEN, width=0.7); ax.set_xlabel("claims k", fontsize=9)
ax.text(2.2, 0.32, "P(z) ∝ log(1 − β(z−1))", fontsize=9, color=MUTE)
# bottom: the line
ax = fig.add_subplot(gs[1, :])
ax.set_xlim(-6.6, 6.6); ax.set_ylim(-1.3, 1.5); ax.axis("off")
ax.plot([-6.3, 0], [0, 0], color=GREEN, lw=7, solid_capstyle="butt")
ax.plot([0, 1], [0, 0], color=GREEN, lw=7, alpha=0.35, solid_capstyle="butt")
ax.plot([1, 6.3], [0, 0], color=RED, lw=2, ls=(0, (3, 3)))
for m in range(1, 7):
    ax.scatter(m, 0, s=130, color=BLUE, zorder=4, edgecolor="white", linewidth=2)
    ax.text(m, -0.38, str(m), ha="center", color=INK, fontsize=11)
for m in range(-6, 1):
    ax.text(m, -0.38, str(m), ha="center", color=INK, fontsize=11)
ax.scatter(0, 0, s=130, facecolor="white", edgecolor=GREEN, linewidth=2.5, zorder=5)
ax.text(-3.2, 0.32, "negative binomial: any negative n\n(pole of order |n| in the claim-count generating function)", ha="center", fontsize=10, color=INK)
ax.text(0.5, 0.75, "0 < n < 1 and n = 0 (log):\nonly with the zero count removed", ha="center", fontsize=9.5, color=INK)
ax.text(3.7, 0.32, "binomial: whole numbers only\n(n = number of policies)", ha="center", fontsize=10, color=INK)
ax.text(3.5, -0.75, "between the whole numbers: forbidden", ha="center", fontsize=9.5, color=RED)
ax.text(-6.5, -1.15, "Poisson ← as n → −∞", fontsize=10, color=INK)
ax.text(6.6, -1.15, "as n → +∞ → Poisson", fontsize=10, color=INK, ha="right")
ax.text(0, -1.15, "dimension n in P(z) = (1 + c(z − 1))ⁿ", ha="center", fontsize=11, color=INK, weight="bold")
save(fig, "a1_dimension_line.png")

# ---------- a2: the infinity test (moments vs dimension) ----------
fig, axs = plt.subplots(1, 2, figsize=(11, 4.3))
ax = axs[0]
al = np.linspace(0.02, 6, 4000)
for kk, col in [(1, BLUE), (2, RED), (3, GREEN)]:
    den = np.prod([al - t for t in range(1, kk + 1)], axis=0)
    y = math.factorial(kk) / den
    y[np.abs(den) < 2e-3] = np.nan
    y[al <= kk] = np.nan          # only the finite branch is a real moment
    ax.plot(al, y, color=col, label=f"E[X^{kk}]")
    ax.axvline(kk, color=col, lw=1, ls=":")
ax.set_ylim(0, 12); ax.set_xlim(0, 6)
ax.set_xlabel("tail dimension α  (survival (1 + x)^−α)"); ax.set_ylabel("moment of a Lomax loss")
ax.legend(frameon=False, loc="upper right")
for kk in [1, 2, 3]:
    ax.text(kk + 0.05, 11.3, f"∞ at α = {kk}", fontsize=9, color=MUTE)
ax.set_title("The k-th moment exists only past dimension k", fontsize=11.5, loc="left")
ax = axs[1]
dims = [1, 2, 3, 4, 5]
est = []
for d in dims:
    R = rng.random(300000) ** (1 / d)
    xs = np.sort(1 / R)[::-1]
    est.append(1 / np.mean(np.log(xs[:3000] / xs[3000])))
ax.plot([0.5, 5.5], [0.5, 5.5], color=GRID, lw=1.5, zorder=1)
ax.scatter(dims, est, s=80, color=BLUE, zorder=3, edgecolor="white", linewidth=1.5)
for d, e in zip(dims, est):
    ax.text(d + 0.12, e - 0.25, f"{e:.2f}", fontsize=9, color=MUTE)
ax.set_xlabel("dimension d of the ball"); ax.set_ylabel("tail index measured from data")
ax.set_xticks(dims)
ax.set_title("1 / distance to a random point in a d-ball:\nthe tail index comes out as d", fontsize=11.5, loc="left")
fig.tight_layout(); save(fig, "a2_infinity_test.png")

# ---------- a3: Danish fire data ----------
x = np.asarray(load("danish.rda", "danish"), float)
fig, axs = plt.subplots(1, 2, figsize=(11, 4.3))
ax = axs[0]
xs = np.sort(x)[::-1]
ks = np.arange(20, 700)
hill = np.array([1 / np.mean(np.log(xs[:k] / xs[k])) for k in ks])
ax.axhspan(1, 2, color="#eef0ff", zorder=0)
ax.plot(ks, hill, color=BLUE)
ax.text(690, 1.06, "mean exists (α > 1)", ha="right", fontsize=9.5, color=INK)
ax.text(690, 2.07, "variance would need α > 2", ha="right", fontsize=9.5, color=INK)
ax.set_ylim(0.8, 2.8); ax.set_xlabel("largest claims used, k"); ax.set_ylabel("Hill tail index α")
ax.set_title("Danish fire losses: about 1.4 dimensions", fontsize=11.5, loc="left")
ax = axs[1]
ns = np.arange(50, len(x) + 1, 10)
for pw, col in [(1, BLUE), (2, RED), (3, GREEN)]:
    r = [np.max(x[:m] ** pw) / np.sum(x[:m] ** pw) for m in ns]
    ax.plot(ns, r, color=col)
    ax.text(ns[-1] + 25, r[-1], f"x^{pw}", va="center", fontsize=10, color=INK)
ax.set_ylim(0, 1); ax.set_xlim(0, len(x) + 200)
ax.set_xlabel("claims seen so far (time order)"); ax.set_ylabel("largest / sum")
ax.set_title("The infinity test on each moment:\nshare of the total held by one claim", fontsize=11.5, loc="left")
fig.tight_layout(); save(fig, "a3_danish.png")

# ---------- a4: ruin ----------
fig, axs = plt.subplots(1, 2, figsize=(11, 4.4))
ax = axs[0]
theta = 0.2
s = np.linspace(0, 1.4, 800)
line = 1 + (1 + theta) * s
ax.plot(s, line, color=INK, lw=1.5, ls="--")
ax.text(1.05, 1 + 1.2 * 1.05 - 0.45, "1 + (1 + θ) μ s\n(premium line)", fontsize=9, color=INK)
se = s[s < 0.985]; ax.plot(se, 1 / (1 - se), color=BLUE)
sg = s[s < 2.9]; ax.plot(sg, (1 - sg / 3) ** -3, color=GREEN)
ax.axvline(1, color=BLUE, lw=1, ls=":"); ax.text(1.01, 5.6, "pole at s = 1\n(order 1)", fontsize=9, color=MUTE)
ax.text(0.18, 5.3, "exponential", fontsize=10, color=BLUE)
ax.text(1.22, 4.25, "gamma 3", fontsize=10, color=GREEN)
for R, col in [(1 / 6, BLUE), (0.257705, GREEN)]:
    ax.scatter(R, 1 + 1.2 * R, s=60, color=col, zorder=4, edgecolor="white", linewidth=1.5)
ax.text(0.2, 1.05, "R", fontsize=10, color=INK)
ax.plot([0, 0], [1, 6], color=RED, lw=3)
ax.text(0.03, 2.2, "Lomax: ∞ for every s > 0,\nno crossing, no R", fontsize=9, color=RED)
ax.set_ylim(0.8, 6.2); ax.set_xlim(-0.03, 1.4)
ax.set_xlabel("s"); ax.set_ylabel("moment generating function M(s)")
ax.set_title("Lundberg's root must come before the first infinity", fontsize=11.5, loc="left")
ax = axs[1]
u = np.linspace(0.1, 300, 600)
ax.semilogy(u, np.exp(-u / 6) / 1.2, color=BLUE)
al = 2.5; mu = 1 / (al - 1)
ax.semilogy(u, (1 / (theta * mu)) * (1 + u) ** (-(al - 1)) / (al - 1), color=RED)
us = [0, 5, 10, 20, 40]; sim_e = [0.83447, 0.36394, 0.15859, 0.02976, 0.00111]
ax.scatter(us, sim_e, color=BLUE, s=40, zorder=4, edgecolor="white")
us2 = [10, 30, 100, 300]; sim_l = [0.26594, 0.06766, 0.00693, 0.00106]
ax.scatter(us2, sim_l, color=RED, s=40, zorder=4, edgecolor="white")
ax.set_ylim(1e-6, 1.5); ax.set_xlim(0, 300)
ax.text(90, 4e-5, "light tail: ruin falls like e^(−Ru)", fontsize=9.5, color=BLUE)
ax.text(120, 0.012, "heavy tail α = 2.5: ruin falls like u^(−1.5)", fontsize=9.5, color=RED)
ax.set_xlabel("starting capital u (in mean claims)"); ax.set_ylabel("probability of ruin")
ax.set_title("Dots: simulated. Lines: formula or asymptote", fontsize=11.5, loc="left")
fig.tight_layout(); save(fig, "a4_ruin.png")

# ---------- a5: decrement paths ----------
q1, q2 = 0.3, 0.5
fig, axs = plt.subplots(1, 2, figsize=(11, 4.6), gridspec_kw={"width_ratios": [1.15, 1]})
ax = axs[0]
ax.add_patch(Polygon([[0, 0], [q1, 0], [q1, q2], [0, q2]], closed=True, fill=False, edgecolor=INK, lw=1.5))
t = np.linspace(0, 1, 200)
ax.plot(q1 * t, q2 * t, color=BLUE, label="uniform within the year (UDD)")
m1, m2 = -np.log(1 - q1), -np.log(1 - q2)
ax.plot(1 - np.exp(-m1 * t), 1 - np.exp(-m2 * t), color=RED, ls="--", label="constant force")
ax.plot(np.concatenate([q1 * np.minimum(2 * t, 1)]), q2 * t, color=GREEN, ls=":", lw=2.5, label="deaths early in the year")
ax.scatter([q1], [q2], s=60, color=INK, zorder=4)
ax.text(q1 + 0.005, q2 + 0.01, "year end", fontsize=9, color=INK)
ax.set_xlabel("P(death would have happened by t)"); ax.set_ylabel("P(lapse would have happened by t)")
ax.set_xlim(-0.01, 0.34); ax.set_ylim(-0.01, 0.56); ax.set_aspect("equal")
ax.legend(frameon=False, fontsize=9, loc="upper center", bbox_to_anchor=(0.5, -0.17), ncol=1)
ax.set_title("Three timing assumptions = three paths\nthrough the same square", fontsize=11.5, loc="left")
ax = axs[1]
def rates(Q):
    ts = np.linspace(0, 1, 100001); P = np.array([Q(tt) for tt in ts]); dP = np.diff(P, axis=0); mid = (P[1:] + P[:-1]) / 2
    return np.array([np.sum((1 - mid[:, 1]) * dP[:, 0]), np.sum((1 - mid[:, 0]) * dP[:, 1])])
r_u = rates(lambda tt: np.array([q1 * tt, q2 * tt]))
r_c = rates(lambda tt: np.array([1 - np.exp(-m1 * tt), 1 - np.exp(-m2 * tt)]))
r_e = rates(lambda tt: np.array([q1 * min(2 * tt, 1), q2 * tt]))
tot = 1 - (1 - q1) * (1 - q2)
labels = ["UDD", "constant force", "deaths early"]
for i_, (r, col) in enumerate([(r_u, BLUE), (r_c, RED), (r_e, GREEN)]):
    ax.barh(i_, r[0], color=col, height=0.55)
    ax.barh(i_, r[1], left=r[0] + 0.004, color=col, alpha=0.45, height=0.55)
    ax.text(r[0] / 2, i_, f"death {r[0]:.4f}", ha="center", va="center", color="white", fontsize=9)
    ax.text(r[0] + r[1] / 2, i_, f"lapse {r[1]:.4f}", ha="center", va="center", color=INK, fontsize=9)
ax.set_yticks(range(3)); ax.set_yticklabels(labels); ax.invert_yaxis(); ax.grid(False)
ax.axvline(tot + 0.004, color=INK, lw=1, ls=":"); ax.text(tot + 0.006, -0.42, f"total {tot:.2f} on every path", fontsize=9, color=INK)
ax.set_xlim(0, 0.75); ax.set_xlabel("dependent rate (probability the year ends this way)")
ax.set_title("The total is the same; the split depends on the path\n(illustrative rates 30% and 50%)", fontsize=11.5, loc="left")
fig.tight_layout(); save(fig, "a5_decrement_paths.png")

# ---------- a6: interest cube ----------
fig, axs = plt.subplots(1, 2, figsize=(11, 4.2), gridspec_kw={"width_ratios": [1, 1.2]})
ax = axs[0]
P = {c: (c[0] + 0.45 * c[2], c[1] + 0.35 * c[2]) for c in itertools.product([0, 1], repeat=3)}
for a_, b_ in itertools.combinations(P, 2):
    if sum(abs(np.array(a_) - np.array(b_))) == 1:
        ax.plot(*zip(P[a_], P[b_]), color=MUTE, lw=1)
lab = {0: "1", 1: "i", 2: "i²", 3: "i³"}
for c, (x_, y_) in P.items():
    ax.scatter(x_, y_, s=90, color=cols[sum(c)], zorder=3, edgecolor="white", linewidth=1.5)
    ax.text(x_ + 0.05, y_ + 0.06, lab[sum(c)], fontsize=11, color=INK)
ax.text(0.7, -0.35, "(1 + i)³ = 1 + 3i + 3i² + i³", ha="center", fontsize=10.5, color=INK)
ax.set_xlim(-0.2, 1.75); ax.set_ylim(-0.5, 1.55); ax.axis("off")
ax.set_title("Three years of interest: each corner is a\nset of years whose interest earns interest", fontsize=11.5, loc="left")
ax = axs[1]
i = 0.05
for n_, col in [(10, BLUE), (30, RED)]:
    k = np.arange(1, n_ + 1)
    sh = np.array([math.comb(n_, kk) * i ** kk for kk in k]) / ((1 + i) ** n_ - 1)
    ax.plot(k[:12], sh[:12], marker="o", color=col, ms=6, label=f"{n_} years")
ax.set_xlabel("cube layer (how many years compound together)"); ax.set_ylabel("share of total interest at 5%")
ax.legend(frameon=False); ax.set_xticks(range(1, 13))
ax.text(1.3, 0.72, "layer 1 = simple interest", fontsize=9.5, color=MUTE)
ax.set_title("Over 30 years most interest lives in the higher layers", fontsize=11.5, loc="left")
fig.tight_layout(); save(fig, "a6_interest_cube.png")

# ---------- a7: mortality ----------
us = load("demoUsa.rda", "demoUsa"); ilt = load("soa08Act.rda", "soa08Act")
lm = us["USSS2007M"].to_numpy(float); age = us["age"].to_numpy()
g = np.isfinite(lm) & (lm > 0); lm, age = lm[g], age[g]
mum = -np.log(lm[1:] / lm[:-1])
lx = np.array(ilt.lx, float); ax_ = np.array(ilt.x)
mui = -np.log(lx[1:] / lx[:-1])
fig, axs = plt.subplots(1, 2, figsize=(11, 4.4))
ax = axs[0]
ax.semilogy(age[:-1], mum, color=BLUE, label="US 2007 men (Social Security)")
sel = ax_[:-1] <= 115
ax.semilogy(ax_[:-1][sel], mui[sel], color=GREEN, ls="--", label="SOA illustrative table (Makeham)")
xx = np.linspace(60, 109.5, 400)
for kk, ls_ in [(3, ":")]:
    ax.semilogy(xx, kk / (110 - xx), color=RED, ls=ls_, label=f"cube of dimension {kk}, end age 110")
ax.set_xlim(0, 115); ax.set_ylim(1e-4, 30)
ax.set_xlabel("age"); ax.set_ylabel("force of mortality (per year)")
ax.legend(frameon=False, fontsize=9, loc="lower right")
ax.annotate("table flattens at about 0.7:\nan assumption, not a pole", (108, 0.69), (35, 4), fontsize=9, color=MUTE, arrowprops=dict(arrowstyle="-", color=MUTE))
ax.set_title("Gompertz line vs a pole at a fixed end age", fontsize=11.5, loc="left")
ax = axs[1]
s50 = (age[:-1] >= 50)
ax.plot(age[:-1][s50], 1 / mum[s50], color=BLUE, label="US 2007 men")
xx = np.linspace(50, 110, 200)
for kk, col, ls_ in [(1, RED, ":"), (3, AMBER, "-.")]:
    ax.plot(xx, (110 - xx) / kk, color=col, ls=ls_, label=f"dimension {kk}: hits zero at 110")
ax.set_ylim(0, 45); ax.set_xlim(50, 112)
ax.set_xlabel("age"); ax.set_ylabel("1 / force of mortality (years)")
ax.legend(frameon=False, fontsize=9)
ax.set_title("A pole shows as a straight line to zero;\nslope = 1 / dimension", fontsize=11.5, loc="left")
fig.tight_layout(); save(fig, "a7_mortality.png")
print("figures written to", OUT)
