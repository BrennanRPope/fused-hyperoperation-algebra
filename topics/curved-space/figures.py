"""Pictures for the curved-space write-up (DRAFT). Run after checks.py (reads results.json).
Writes figures/f1 ... f8. Colour key: orange = known, blue = ours (reading or result), grey = context."""
import os, json, math
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon, Arc, Rectangle

HERE = os.path.dirname(os.path.abspath(__file__))
FIG = os.path.join(HERE, "figures"); os.makedirs(FIG, exist_ok=True)
RES = json.load(open(os.path.join(HERE, "results.json")))

BLUE, ORANGE, AQUA, YELLOW, MAGENTA, GREEN, VIOLET, RED = "#257de0", "#f47230", "#00b894", "#f2a600", "#ee7fac", "#008a17", "#6a4ce2", "#ec4647"
INK, INK2, GRID = "#191919", "#666666", "#e6e6e6"
plt.rcParams.update({"font.size": 10, "axes.edgecolor": INK2, "axes.labelcolor": INK, "xtick.color": INK2,
                     "ytick.color": INK2, "axes.grid": True, "grid.color": GRID, "grid.linewidth": 0.8,
                     "axes.spines.top": False, "axes.spines.right": False, "lines.linewidth": 2,
                     "legend.frameon": False, "axes.titleweight": "bold", "axes.titlesize": 11})


def save(fig, name):
    fig.tight_layout(rect=[0, 0, 1, 0.9])
    fig.savefig(os.path.join(FIG, name), dpi=150, bbox_inches="tight"); plt.close(fig); print("wrote", name)


def sprinkle(K, rho, lo, hi, seed=2):
    r = np.random.default_rng(seed)
    dmax = max(2 / (1 - K * x * y)**2 for x in (lo, hi) for y in (lo, hi))
    n = r.poisson(rho * dmax * (hi - lo)**2)
    p = r.uniform(lo, hi, (n, 2))
    keep = r.uniform(0, dmax, n) < 2 / (1 - K * p[:, 0] * p[:, 1])**2
    return p[keep]


def to_tx(u, v):          # draw light-cone coordinates as a diamond: time up, space across
    return (u - v), (u + v)


def corner_vol(c, K):
    return 2 * c if K == 0 else -(2 / K) * math.log(1 - K * c)


# ---------- f1: cut a causal diamond with a light ray ----------
fig, axs = plt.subplots(1, 3, figsize=(13.5, 4.8))
a, b, m = 0.9, 0.9, 0.42
for ax, K, name, rule in zip(axs, (0.0, 1.0, -1.0),
                             ("Flat (K = 0)", "Curved like a sphere (de Sitter, K = +1)", "Curved like a saddle (anti-de Sitter, K = −1)"),
                             ("c = c₁ + c₂", "c = c₁ + c₂ − c₁c₂   (∘₋₁)", "c = c₁ + c₂ + c₁c₂   (∘₁)")):
    pts = sprinkle(K, 900, 0, 0.97, seed=5)
    inside = (pts[:, 0] < a) & (pts[:, 1] < b)
    left = inside & (pts[:, 0] < m)
    right = inside & (pts[:, 0] >= m)
    for sel, col in ((left, BLUE), (right, AQUA)):
        x, t = to_tx(pts[sel, 0], pts[sel, 1])
        ax.scatter(x, t, s=4, color=col, lw=0)
    for (u0, v0, u1, v1) in ((0, 0, a, 0), (a, 0, a, b), (a, b, 0, b), (0, b, 0, 0)):
        x0, t0 = to_tx(u0, v0); x1, t1 = to_tx(u1, v1)
        ax.plot([x0, x1], [t0, t1], color=INK, lw=1.4)
    x0, t0 = to_tx(m, 0); x1, t1 = to_tx(m, b)
    ax.plot([x0, x1], [t0, t1], color=RED, lw=2.2, ls="--")
    c1, c2_ = m * b, (a * b - m * b) / (1 - K * m * b) if K else (a - m) * b
    n1, n2 = left.sum(), right.sum()
    ax.text(-0.8, 1.45, f"left piece\nc₁ = {c1:.3f}\n{n1} points", color=BLUE, ha="center", fontsize=9)
    ax.text(0.8, 1.45, f"right piece\nc₂ = {c2_:.3f}\n{n2} points", color="#00806a", ha="center", fontsize=9)
    ax.set_title(name, fontsize=10.5)
    ax.text(0, -0.2, "whole corner c = ab = " + f"{a*b:.3f}\n" + rule, ha="center", fontsize=10, color=BLUE if K else ORANGE,
            fontweight="bold")
    ax.set_xlim(-1.05, 1.05); ax.set_ylim(-0.32, 2.0); ax.set_aspect("equal"); ax.grid(False)
    ax.set_xticks([]); ax.set_yticks([])
    ax.text(0.5, 0.2, "cut along\na light ray", color=RED, fontsize=8.5)
    for sp_ in ax.spines.values():
        sp_.set_visible(False)
fig.suptitle("f1. Cut a causal diamond along a light ray. Points are sprinkled at random, the same number per unit of spacetime volume everywhere.\n"
             "Flat: the two corner areas add (orange, known). Curved: corners combine by our operation with dial k = −K (blue, ours).",
             fontsize=10.5)
save(fig, "f1_cut_diamond.png")

# ---------- f2: counting = the ruler of the corner ----------
fig, ax = plt.subplots(1, 2, figsize=(12, 4.4))
cs = np.linspace(0, 0.95, 300)
for K, col in ((-1.0, VIOLET), (-0.5, "#9a86ee"), (0.0, INK2), (0.5, "#6fb0f0"), (1.0, BLUE)):
    V = [corner_vol(c, K) / 2 for c in cs]
    ax[0].plot(cs, V, color=col, label=f"K = {K:+g}" + ("  (flat: V/2 = c)" if K == 0 else ""))
# Monte Carlo dots
rng = np.random.default_rng(11)
for K, col in ((1.0, BLUE), (-1.0, VIOLET)):
    pts = sprinkle(K, 3000, 0, 0.97, seed=9)
    for c in (0.15, 0.35, 0.55, 0.75, 0.9):
        side = math.sqrt(c)
        n = np.sum((pts[:, 0] < side) & (pts[:, 1] < side))
        ax[0].scatter([c], [n / 3000 / 2], color=col, edgecolor="white", zorder=5, s=40)
ax[0].set_xlabel("corner c = (light-cone width) × (light-cone width)")
ax[0].set_ylabel("counted volume / 2  =  L₋K(c) = −log(1 − Kc)/K")
ax[0].set_title("A. Points counted in a diamond (dots) lie on the ruler of its corner")
ax[0].legend(fontsize=8.5, loc="upper left")
taus = np.linspace(0.02, 2.0, 300)
for K, col in ((1.0, BLUE), (0.0, INK2), (-1.0, VIOLET)):
    if K > 0:
        V = (4 / K) * np.log(np.cosh(np.sqrt(K) * taus / 2))
    elif K < 0:
        V = (4 / K) * np.log(np.cos(np.sqrt(-K) * taus / 2))
    else:
        V = taus**2 / 2
    ax[1].plot(taus, V / (taus**2 / 2 + 1e-30), color=col, label=f"K = {K:+g}")
ax[1].plot(taus, 1 - taus**2 / 24, color=ORANGE, ls="--", lw=1.4, label="1 − Kτ²/24 (known small-diamond law, K = +1)")
ax[1].set_xlabel("proper time τ across the diamond"); ax[1].set_ylabel("volume ÷ flat volume τ²/2")
ax[1].set_title("B. Curvature shows up as a shortfall (K > 0) or surplus (K < 0)")
ax[1].legend(fontsize=8.5)
fig.suptitle("f2. In a curved 2D spacetime the number of points in a causal diamond is the ruler L_k of its corner, with dial k = −K",
             fontsize=10.5)
save(fig, "f2_count_is_ruler.png")

# ---------- f3: the same shape as the AdS3 gluing ----------
fig, ax = plt.subplots(1, 3, figsize=(14.5, 4.6))
a1_, a2_, b1_, mm_, b2_ = 0.0, 1.0, 1.6, 2.4, 3.6
A = ax[0]
A.axhline(0, color=INK, lw=1.2)
for (p, q, col, lab) in ((a1_, a2_, ORANGE, "A"), (b1_, mm_, BLUE, "B₁"), (mm_, b2_, AQUA, "B₂")):
    A.plot([p, q], [0, 0], color=col, lw=6, solid_capstyle="butt")
    A.text((p + q) / 2, -0.18, lab, ha="center", color=col, fontweight="bold")
for (p, q, col) in ((a1_, b2_, INK2), (a2_, b1_, INK2), (a1_, b1_, "#bbbbbb"), (a2_, b2_, "#bbbbbb"), (a2_, mm_, "#bbbbbb"),
                    (a1_, mm_, "#bbbbbb")):
    A.add_patch(Arc(((p + q) / 2, 0), q - p, q - p, theta1=0, theta2=180, color=col, lw=1.3))
A.set_xlim(-0.3, 3.9); A.set_ylim(-0.35, 2.0); A.set_aspect("equal"); A.grid(False); A.set_yticks([])
A.set_title("A. Space at one moment in AdS₃ (known)\nregions on the edge, shortest paths inside")
A.text(1.8, 1.85, "x(A, B₁B₂) = x(A,B₁) ∘₋₁ x(A,B₂)", ha="center", color=BLUE, fontsize=9.5)
B = ax[1]
# kinematic space: each interval [u, w] is a point; draw with w on the vertical axis
B.plot([-0.3, 3.9], [-0.3, 3.9], color="#bbbbbb", lw=1)
B.text(2.3, 1.95, "w = u (empty intervals)", color=INK2, fontsize=8, rotation=45)
P, Q = (a1_, b2_), (a2_, b1_)
corners = [(a1_, b2_), (a2_, b2_), (a2_, b1_), (a1_, b1_)]
B.add_patch(Polygon(corners, closed=True, facecolor="#dceafb", edgecolor=INK, lw=1.3))
B.plot([a1_, a2_], [mm_, mm_], color=RED, ls="--", lw=2)
for pt, lab in ((P, "[a₁,b₂] = A ∪ gap ∪ B"), (Q, "[a₂,b₁] = the gap"), ((a1_, mm_), "[a₁,m]"), ((a2_, mm_), "[a₂,m]")):
    B.scatter(*pt, color=INK, s=20, zorder=5)
    B.text(pt[0] + 0.06, pt[1] + 0.06, lab, fontsize=8)
B.text(0.5, 3.05, "piece 2:\nx(A,B₂)", ha="center", color="#00806a", fontsize=8.5)
B.text(0.5, 1.95, "piece 1:\nx(A,B₁)", ha="center", color=BLUE, fontsize=8.5)
B.set_xlabel("u (left end of an interval)"); B.set_ylabel("w (right end)")
B.set_xlim(-0.3, 3.9); B.set_ylim(-0.3, 3.9); B.set_aspect("equal")
B.set_title("B. The plane of intervals is a curved 2D spacetime\n(kinematic space, known); gluing = cutting a diamond")
C = ax[2]
xx = np.linspace(0, 0.97, 300)
C.plot(xx, -np.log(1 - xx), color=BLUE)
C.set_xlabel("cross-ratio x  (= K × corner)"); C.set_ylabel("t = −log(1 − x)")
C.set_title("C. One law for both: 1 − x = e^(−t)\nt adds when pieces join, x joins by ∘₋₁")
C.text(0.05, 2.6, "AdS₃: t ∝ conditional mutual information\n2D causal set: t = K V/2 (points counted)",
       fontsize=9, color=INK)
fig.suptitle("f3. Same shape, and in fact the same space: the curved plane whose points are intervals (de Sitter 2D) turns the AdS₃\n"
             "gluing rule found earlier into the diamond-cutting rule of f1. The identification is known (Czech et al. 2015); the ∘₋₁ reading is ours.",
             fontsize=10)
save(fig, "f3_same_shape_ads3.png")

# ---------- f4: curvature from counting ----------
fig, ax = plt.subplots(1, 2, figsize=(12.5, 4.8))
A = ax[0]
r = np.random.default_rng(4)
n = 220
pts = r.uniform(0, 1, (n, 2))
top = np.array([1.0, 1.0])


def between(p):
    return np.sum((pts[:, 0] > p[0]) & (pts[:, 1] > p[1]) & (pts[:, 0] < 1) & (pts[:, 1] < 1))


cols = {0: BLUE, 1: RED, 2: AQUA}
for p in pts:
    k = between(p)
    x, t = to_tx(p[0], p[1])
    A.scatter(x, t, s=16 if k <= 2 else 5, color=cols.get(k, "#cccccc"), lw=0, zorder=3 if k <= 2 else 1)
xt = to_tx(1, 1)
A.scatter(*xt, s=60, color=INK, zorder=6, marker="*")
A.text(xt[0] + 0.05, xt[1], "x", fontsize=12)
for (u0, v0, u1, v1) in ((1, 1, 0, 1), (1, 1, 1, 0)):
    x0, t0 = to_tx(u0, v0); x1, t1 = to_tx(u1, v1)
    A.plot([x0, x1], [t0, t1], color=INK2, lw=1)
A.text(-2.3, 1.55, "weight +1: nothing in between (links)", color=BLUE, fontsize=9)
A.text(-2.3, 1.4, "weight −2: one point in between", color=RED, fontsize=9)
A.text(-2.3, 1.25, "weight +1: two points in between", color="#00806a", fontsize=9)
A.text(-2.3, 0.95, "(1, −2, 1) = Möbius signs of the 2-cube\n            grouped by size: + / − − / +", color=BLUE, fontsize=9,
       fontweight="bold")
A.set_aspect("equal"); A.grid(False); A.set_xticks([]); A.set_yticks([])
A.set_xlim(-2.35, 1.05); A.set_ylim(-0.05, 2.1)
for sp_ in A.spines.values():
    sp_.set_visible(False)
A.set_title("A. Benincasa–Dowker: weigh the past of x by how many points sit in between")
Bx = ax[1]
vals = RES["B2_values"]
rhos = [1e2, 1e3, 1e4]
for key, col, lab in (("1.0", BLUE, "K = +1 → −1"), ("-1.0", VIOLET, "K = −1 → +1"), ("0.0", INK2, "flat → 0")):
    Bx.plot(rhos, vals[key], "o-", color=col, label=lab)
Bx.axhline(-1, color=BLUE, ls=":", lw=1); Bx.axhline(1, color=VIOLET, ls=":", lw=1)
Bx.set_xscale("log"); Bx.set_xlabel("points per unit volume ρ"); Bx.set_ylabel("expected weighted count  (→ −R/2 = −K)")
Bx.set_title("B. The weighted count converges to the curvature (exact expectation)")
Bx.legend(fontsize=9, loc="center left", bbox_to_anchor=(0.02, 0.6))
fig.suptitle("f4. Curvature from counting (known: Benincasa–Dowker 2010). Ours: the weights are the cube's corner weights at k = −1,\n"
             "and the curvature enters the count only through the inverse lift 1/(1 − K·corner).", fontsize=10.5)
save(fig, "f4_counting_curvature.png")

# ---------- f5: cube order vs light-cone order ----------
fig, ax = plt.subplots(figsize=(7.2, 4.4))
ns = np.arange(2, 9)
cube = 2.0**(1 - ns)
mm = [math.gamma(d + 1) * math.gamma(d / 2) / (2 * math.gamma(3 * d / 2)) for d in ns]
ax.plot(ns, cube, "o-", color=BLUE, label="n-cube product order: 2^(1−n)")
ax.plot(ns, mm, "s-", color=ORANGE, label="round light cone in n dimensions (Myrheim–Meyer, known)")
ax.set_yscale("log"); ax.set_xlabel("dimension n"); ax.set_ylabel("chance two random points are causally related")
ax.annotate("equal only at n = 2", xy=(2, 0.5), xytext=(3.2, 0.42), arrowprops=dict(arrowstyle="->", color=INK2))
ax.set_title("f5. The cube's order is spacetime's order only in 1 + 1 dimensions")
ax.legend(fontsize=9)
save(fig, "f5_cube_vs_lightcone.png")

# ---------- f6: river profiles ----------
fig, ax = plt.subplots(1, 2, figsize=(12, 4.4))
T = np.linspace(0, 2.0, 300)
M_, H_ = 1.0, 0.6
r_L = 1.0 * np.exp(H_ * T)
r_M = (1.5 * np.sqrt(2 * M_) * T + 1.0)**(2 / 3)
r_S = (np.sqrt(2 * M_) / H_ * np.sinh(1.5 * H_ * T + np.arcsinh(H_ / np.sqrt(2 * M_))))**(2 / 3)
ax[0].plot(T, r_L, color=BLUE, label="dark energy only: r = r₀ e^(HT)  →  r/r₀ − 1 = E(HT)")
ax[0].plot(T, r_M, color=ORANGE, label="mass only: r^(3/2) grows at a steady rate")
ax[0].plot(T, r_S, color=GREEN, label="both: r^(3/2) ∝ sinh(3HT/2)  (the ΛCDM law)")
ax[0].set_xlabel("river time T (a falling clock's own time)"); ax[0].set_ylabel("distance r of a river parcel")
ax[0].set_title("A. Three rivers, run outward (inflow is the same played backward)")
ax[0].legend(fontsize=8.5)
ax[1].plot(T, np.log(r_L), color=BLUE, label="log r  (dark energy): straight")
ax[1].plot(T, r_M**1.5, color=ORANGE, label="r^(3/2)  (mass): straight")
ax[1].plot(T, np.arcsinh(H_ * r_S**1.5 / np.sqrt(2 * M_)), color=GREEN, label="asinh(H r^(3/2)/√2M)  (both): straight")
ax[1].set_xlabel("river time T"); ax[1].set_ylabel("the variable in which the river is a steady drift")
ax[1].set_title("B. Only the dark-energy river drifts in log r, i.e. is our flow E")
ax[1].legend(fontsize=8.5)
fig.suptitle("f6. Does the height flow E(x) = eˣ − 1 describe space falling into a mass? Only the dark-energy part. The mass part\n"
             "is a steady drift in r^(3/2) (known: Painlevé, Gullstrand; the two together give the ΛCDM scale-factor law).", fontsize=10.5)
save(fig, "f6_river_profiles.png")

# ---------- f7: drop speeds and orbits ----------
fig, ax = plt.subplots(1, 2, figsize=(12.5, 4.5))
rr = np.linspace(2.05, 40, 400)
ax[0].plot(rr, 1 - 2 / rr, color=INK2, label="clock rate f(r) = 1 − 2M/r")
for (ra, rb, col) in ((30, 12, BLUE), (12, 4, AQUA), (30, 4, VIOLET)):
    w2 = 1 - (1 - 2 / rb) / (1 - 2 / ra)
    ax[0].annotate("", xy=(rb, 1 - 2 / rb), xytext=(ra, 1 - 2 / ra), arrowprops=dict(arrowstyle="->", color=col, lw=2))
    ax[0].text((ra + rb) / 2, 1 - 2 / ((ra + rb) / 2) - (0.07 if col != VIOLET else 0.2), f"w² = {w2:.3f}", color=col, fontsize=9)
_w = lambda ra, rb: 1 - (1 - 2 / rb) / (1 - 2 / ra)
_p, _q, _s = _w(30, 12), _w(12, 4), _w(30, 4)
ax[0].text(13, 0.22, f"{_p:.3f} ∘₋₁ {_q:.3f} = {_p:.3f} + {_q:.3f} − {_p:.3f}·{_q:.3f} = {_s:.3f}\n(drop 30M → 12M, then 12M → 4M  =  drop 30M → 4M)",
           fontsize=8.5, color=BLUE)
ax[0].set_xlabel("r / M"); ax[0].set_ylabel("f(r)")
ax[0].set_title("A. Squared speeds of things dropped from rest stack by ∘₋₁")
ax[0].legend(fontsize=8.5, loc="lower right")
xm = np.linspace(0, 0.49, 400)
ax[1].plot(xm, xm / (1 - 2 * xm), color=BLUE, label="relativity: φ₋₂(M/r) = (M/r)/(1 − 2M/r)")
ax[1].plot(xm, xm, color=ORANGE, ls="--", label="Newton: v² = M/r")
for xv, lab, yv in ((1 / 6, "ISCO, v = c/2", 0.25), (1 / 3, "photon sphere, v = c", 1.0)):
    ax[1].scatter([xv], [yv], color=INK, zorder=5)
    ax[1].text(xv + 0.01, yv + 0.08, lab, fontsize=9)
ax[1].axvline(0.5, color=RED, ls=":", lw=1.2); ax[1].text(0.47, 2.5, "flip's wall = horizon", rotation=90, color=RED, fontsize=9)
ax[1].set_ylim(0, 3); ax[1].set_xlabel("M/r"); ax[1].set_ylabel("orbital speed² seen by a hovering observer")
ax[1].set_title("B. Circular orbits: the flip φ₋₂ of Newton's speed²")
ax[1].legend(fontsize=8.5)
fig.suptitle("f7. The river's speeds in our operations (exact; the physics is known, the readings are ours)", fontsize=10.5)
save(fig, "f7_drop_and_orbits.png")

# ---------- f8: the horizon ----------
fig, ax = plt.subplots(1, 2, figsize=(12, 4.6))
A = ax[0]
A.plot([-2, 2], [-2, 2], color=INK, lw=1.2); A.plot([-2, 2], [2, -2], color=INK, lw=1.2)
for tv in np.linspace(-1.5, 1.5, 7):
    th = np.tanh(tv / 2)          # slope of t = const lines in Kruskal (T/X = tanh(t/4M))
    X = np.linspace(0, 2, 10)
    A.plot(X, th * X, color=BLUE, lw=1)
for rv in (2.4, 3.0, 4.0):
    X0 = math.sqrt((rv / 2 - 1) * math.exp(rv / 2))
    eta = np.linspace(-1.2, 1.2, 100)
    A.plot(X0 * np.cosh(eta), X0 * np.sinh(eta), color=ORANGE, lw=1.2)
A.text(1.0, 1.75, "horizon r = 2M", rotation=45, fontsize=8.5)
A.text(1.15, -0.55, "hovering observers\n(r = const, orange)", fontsize=8.5, color=ORANGE)
A.text(1.2, 0.95, "t = const\n(blue)", fontsize=8.5, color=BLUE)
A.set_xlim(-0.3, 2.1); A.set_ylim(-1.6, 2.1); A.set_aspect("equal"); A.grid(False)
A.set_xlabel("Kruskal X"); A.set_ylabel("Kruskal T")
A.set_title("A. Waiting Δt is a boost: V → e^(κΔt) V, U → e^(−κΔt) U")
Bx = ax[1]
c7 = RES["C7_curve"]
t = np.array(c7["t"]); y = np.array(c7["log1pz"])
o = np.argsort(t); t, y = t[o], y[o]
Bx.plot(t - t.min(), y, color=BLUE, label="falling source, computed")
sel = slice(-12, None)
line = y[-1] + 0.25 * (t - t[-1])
Bx.plot(t - t.min(), line, color=ORANGE, ls="--", lw=1.4, label="straight line of slope κ = 1/(4M)")
Bx.set_xlabel("arrival time far away (units of M)"); Bx.set_ylabel("log(1 + z)")
Bx.set_title("B. A falling source fades as 1 + z ∝ e^(κt): the flow E")
Bx.legend(fontsize=9)
fig.suptitle("f8. Where E does describe the river: near the horizon, time translation is an exact boost (Kruskal 1960, known), so\n"
             "Doppler factors compose as 1 + E(κt) and redshifts by ∘₁. This is the same step that gives the Hawking temperature.",
             fontsize=10.5)
save(fig, "f8_horizon_boost.png")
