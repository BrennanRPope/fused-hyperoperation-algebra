"""Pictures for the homothety-dial survey. Run: python3 figures.py -> figures/*.png
Colour key: blue = our reading, orange = already known, grey = endpoints/background."""
import json, math, itertools, os
import numpy as np, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy import stats
import mpmath as mp

BLUE, ORANGE, GREY, RED = "#1f5fbf", "#e07b00", "#8a8a8a", "#c0392b"
plt.rcParams.update({"font.size": 14, "axes.titlesize": 15, "axes.labelsize": 14, "legend.fontsize": 12,
                     "xtick.labelsize": 12, "ytick.labelsize": 12, "figure.dpi": 130})
os.makedirs("figures", exist_ok=True)
def circ(k, a, b): return a + b + k * a * b
def save(fig, name):
    fig.tight_layout(); fig.savefig("figures/" + name, dpi=130); plt.close(fig); print("figures/" + name)

# ---------- fig0 hero ----------
fig, axs = plt.subplots(2, 1, figsize=(12, 8.2), gridspec_kw={"height_ratios": [1, 1.5]})
ax = axs[0]
ax.axhline(0, color="k", lw=1.5); ax.set_xlim(-2.6, 2.6); ax.set_ylim(-1.2, 1.4); ax.set_yticks([])
ax.axvspan(-2.6, -0.04, color=GREY, alpha=0.12); ax.axvspan(0.04, 2.6, color=GREY, alpha=0.12)
ax.plot([0], [0], "o", color=GREY, ms=14)
ax.text(-1.3, 0.75, "k < 0\none law (a ↦ λa, λ > 0)", ha="center", fontsize=14)
ax.text(1.3, 0.75, "k > 0\none law (a ↦ λa, λ > 0)", ha="center", fontsize=14)
ax.text(0, 0.75, "k = 0\naddition", ha="center", fontsize=14, color="k")
for s in (0.5, 1.0, 2.0):
    ax.annotate("", xy=(1.0 / s, -0.45), xytext=(1.0, -0.45), arrowprops=dict(arrowstyle="->", color=GREY))
    ax.annotate("", xy=(-1.0 / s, -0.45), xytext=(-1.0, -0.45), arrowprops=dict(arrowstyle="->", color=GREY))
ax.annotate("", xy=(1.0, -0.95), xytext=(-1.0, -0.95), arrowprops=dict(arrowstyle="<->", color=GREY))
ax.text(0, -0.85, "a ↦ −a swaps the two sides", ha="center", fontsize=12, color=GREY)
ax.set_title("On the bare line, rescaling sends k to k/λ: only three stops survive")
ax.set_xlabel("homothety k")
ax = axs[1]
ax.axhline(0, color="k", lw=1.2); ax.set_xlim(-2.6, 4.6); ax.set_ylim(-0.5, 6.4); ax.set_yticks([])
rows = [(5.6, "copies per mode −1/m (Katz/Panjer, Shirai–Takahashi)", [-1 / m for m in range(1, 12)], ORANGE),
        (4.6, "alphabet of q letters −q/(q−1) (symmetric channels)", [-q / (q - 1) for q in range(2, 14)], ORANGE),
        (3.6, "bits stay bits: OR −1, XOR −2", [-1, -2], ORANGE),
        (2.6, "q-numbers k = q − 1 (subspace counts at prime powers)", [1, 2, 3, 4], ORANGE),
        (1.6, "Lee–Yang on graphs holds for k ≥ 1", None, ORANGE),
        (0.6, "corner sum sign-definite: k = −1; positive for k ≥ √6 − 3", [-1, math.sqrt(6) - 3], BLUE)]
for y, lab, pts, col in rows:
    if pts is None:
        ax.plot([1, 4.5], [y, y], color=col, lw=6, alpha=0.7)
    else:
        ax.plot(pts, [y] * len(pts), "o", color=col, ms=9)
    ax.text(-2.55, y + 0.33, lab, fontsize=12.5, color=col if col == BLUE else "k")
ax.axvline(math.sqrt(6) - 3, color=BLUE, ls=":", lw=1.2)
for xv, lab in ((-2, "−2"), (-1, "−1"), (0, "0"), (1, "1"), (2, "2")):
    ax.axvline(xv, color=GREY, lw=0.6, alpha=0.5)
ax.set_xticks([-2, -1, 0, 1, 2, 3, 4]); ax.set_xlabel("homothety k (scale pinned: a top value 1, whole numbers, a graph, a spectrum)")
ax.set_title("Once something pins the scale, k becomes a real dial with whole-number stops")
save(fig, "fig0-hero-k-dial.png")

# ---------- fig12 sorting board ----------
import sys
sys.path.insert(0, "../dial-catalog"); from catalog import GROUPS
L = json.load(open("laws.json"))["rows"]
short = [g[0].split(":")[0].split(" (")[0].split(",")[0][:30] for g in GROUPS]
OVR = {6: "Plane law (k, d)", 36: "Fermat n = 3, 4", 42: "i, Euler, complex analysis", 14: "Curvature by counting",
       22: "Logistic and softmax", 30: "Holographic entropy rules", 53: "Jet fields R[e]/(e² − ke)", 52: "Number planes u² = ku − q",
       47: "Logic, sets, algebra", 16: "Bell, Catalan, partitions", 28: "Quantum noise channels", 45: "Biology (genetics)"}
for i_, v_ in OVR.items(): short[i_ - 1] = v_
order = ["slides", "partly", "self", "disguise", "plausible", "no"]
cols = {"slides": BLUE, "partly": "#5a8fd8", "self": GREY, "disguise": "#b0b0b0", "plausible": "#9fbde8", "no": "#dddddd"}
fig, ax = plt.subplots(figsize=(13, 8.6)); ax.axis("off")
xpos = {"slides": 0.0, "partly": 0.25, "selfdis": 0.25, "no": 0.5}
def put(items, x0, y0, title, col, ncol=1, w=0.24):
    ax.text(x0, y0 + 0.02, title, fontsize=15, weight="bold")
    for i, (lab, cand) in enumerate(items):
        cx = x0 + (i % ncol) * w; cy = y0 - 0.032 * (i // ncol + 1)
        ax.add_patch(plt.Rectangle((cx, cy - 0.012), w - 0.01, 0.026, color=col, alpha=0.85))
        ax.text(cx + 0.005, cy + 0.001, (cand + " " if cand else "") + lab, fontsize=9.5, va="center",
                color="white" if col in (BLUE, "#5a8fd8") else "k")
grp = {v: [(short[i], r["candidate"]) for i, r in enumerate(L) if r["verdict"] == v] for v in order}
put(grp["slides"], 0.0, 0.95, f"slides ({len(grp['slides'])})", BLUE)
put(grp["partly"], 0.25, 0.95, f"partly ({len(grp['partly'])})", "#5a8fd8")
put(grp["self"], 0.25, 0.70, f"self ({len(grp['self'])})", GREY)
put(grp["disguise"], 0.25, 0.56, f"t dial in disguise ({len(grp['disguise'])})", "#b0b0b0")
put(grp["no"], 0.5, 0.95, f"no slide ({len(grp['no'])})", "#dddddd", ncol=2)
ax.set_xlim(0, 1); ax.set_ylim(0.25, 1)
ax.set_title("All 65 law groups, sorted by what turning k (d = 0) does", fontsize=16)
save(fig, "fig12-sorting-board.png")

# ---------- fig1 C1 corner expansion ----------
ks = np.linspace(-2, 2, 401)
fig, axs = plt.subplots(1, 2, figsize=(12, 4.8))
p = [0.3, 0.4, 0.5]
fold = lambda k: circ(k, circ(k, p[0], p[1]), p[2])
ax = axs[0]; ax.plot(ks, [fold(k) for k in ks], color=BLUE, lw=2.5)
for kv, lab in ((-1, "k = −1\nP(at least one)\n= 0.79"), (0, "k = 0\nsum 1.2"), (1, "k = 1\n(1.3·1.4·1.5) − 1\n= 1.73")):
    ax.plot(kv, fold(kv), "o", color=ORANGE, ms=10); ax.annotate(lab, (kv, fold(kv)), (kv - 0.35, fold(kv) + 0.35), fontsize=11)
ax.set_xlabel("homothety k"); ax.set_ylabel("0.3 ∘ₖ 0.4 ∘ₖ 0.5"); ax.set_title("One fold, three readings")
ax.set_ylim(0, 2.6)
ax = axs[1]
for n, ls in ((2, ":"), (5, "-"), (10, "--")):
    ax.plot(ks, [((1 + k * 0.3) ** n - 1) / k if abs(k) > 1e-9 else n * 0.3 for k in ks], ls, color=BLUE, lw=2, label=f"n = {n}")
for kv in (-1, 0, 1): ax.axvline(kv, color=GREY, lw=0.8)
ax.text(-1, 6.3, "1 − 0.7ⁿ", ha="center", color=ORANGE); ax.text(0, 6.3, "0.3 n", ha="center", color=ORANGE); ax.text(1, 6.3, "1.3ⁿ − 1", ha="center", color=ORANGE)
ax.set_ylim(0, 6); ax.set_xlabel("homothety k"); ax.set_ylabel("0.3 ∘ₖ-power n"); ax.legend(loc="upper left"); ax.set_title("Powers: repeated union → repeated addition → compounding")
save(fig, "fig1-corner-expansion.png")

# ---------- fig2 C2 Panjer ----------
def dialpmf(kv, z, emax=12):
    w = [1.0]
    for e in range(1, emax + 1): w.append(w[-1] * z * (1 + (e - 1) * kv) / e)
    w = np.array(w); w[w < 0] = 0; return w / w.sum()
fig, axs = plt.subplots(1, 2, figsize=(12.5, 4.8))
ax = axs[0]; mean = 0.5
cases = [(-1, "k = −1 Bernoulli (fermion)"), (-1 / 3, "k = −1/3 binomial, ≤ 3 copies"), (0, "k = 0 Poisson"), (1, "k = 1 geometric (boson)")]
for j, (kv, lab) in enumerate(cases):
    z = mean / (1 + kv * mean)  # mean = z/(1 - k z) -> z = mean/(1 + k mean)
    pm = dialpmf(kv, z, 6)
    ax.bar(np.arange(7) + (j - 1.5) * 0.2, pm, 0.2, color=[GREY, ORANGE, "#f2b366", BLUE][j], label=lab)
ax.set_xlabel("copies in one mode"); ax.set_ylabel("chance"); ax.set_title("Same mean 0.5, four settings of k"); ax.legend(fontsize=10.5)
ax = axs[1]; ax.axhline(0, color="k")
ax.plot([0, 2.2], [0, 0], color=ORANGE, lw=7, alpha=0.7)
ax.plot([-1 / m for m in range(1, 15)], [0] * 14, "o", color=ORANGE, ms=9)
ax.plot(np.linspace(-0.98, -0.02, 300), [0] * 300, color=RED, lw=1, ls=":")
ax.text(-0.95, 0.35, "k < 0: only −1, −1/2, −1/3, …\nbinomial (Katz 1965; Sundt–Jewell 1981)", fontsize=11.5)
ax.text(0.15, -0.55, "k ≥ 0: every value\nPoisson, then negative binomial r = 1/k", fontsize=11.5)
ax.text(-0.95, -0.25, "in between: a negative 'count' (k = −2/3: −1/9)", fontsize=10.5, color=RED)
ax.set_xlim(-1.1, 2.2); ax.set_ylim(-0.9, 0.9); ax.set_yticks([]); ax.set_xlabel("homothety k (probability convention)")
ax.set_title("Which k give a real count")
save(fig, "fig2-panjer-class.png")

# ---------- fig3 C3 Tsallis/Renyi ----------
rng = np.random.default_rng(3); P = np.array([0.5, 0.25, 0.15, 0.1])
ks3 = np.linspace(-1.5, 0.95, 300)
T = [(1 - np.sum(P ** (1 - k))) / (-k) if abs(k) > 1e-9 else -np.sum(P * np.log(P)) for k in ks3]
Rn = [np.log(np.sum(P ** (1 - k))) / k if abs(k) > 1e-9 else -np.sum(P * np.log(P)) for k in ks3]
fig, axs = plt.subplots(1, 2, figsize=(12, 4.6))
ax = axs[0]; ax.plot(ks3, T, color=ORANGE, lw=2.5, label="Tsallis T (composes by ∘ₖ)"); ax.plot(ks3, Rn, color=BLUE, lw=2.5, label="Rényi = log(1 + kT)/k (adds)")
ax.axvline(0, color=GREY); ax.text(0.03, 2.3, "k = 0: Shannon", color=GREY)
ax.set_xlabel("homothety k = 1 − q (Tsallis index q)"); ax.set_ylabel("entropy (nats)"); ax.legend(); ax.set_title("One distribution along the dial")
ax = axs[1]
P2 = np.array([0.6, 0.3, 0.1]); PJ = np.outer(P, P2).ravel()
gap = []; pred = []
for k in ks3:
    if abs(k) < 1e-9: gap.append(0); pred.append(0); continue
    q = 1 - k; t1 = (1 - np.sum(P ** q)) / (q - 1); t2 = (1 - np.sum(P2 ** q)) / (q - 1); tj = (1 - np.sum(PJ ** q)) / (q - 1)
    gap.append(tj - t1 - t2); pred.append(k * t1 * t2)
ax.plot(ks3, gap, color=ORANGE, lw=4, alpha=0.5, label="T(A,B) − T(A) − T(B)"); ax.plot(ks3, pred, "--", color=BLUE, lw=2, label="k·T(A)·T(B)")
ax.set_xlabel("homothety k"); ax.set_ylabel("excess over adding"); ax.legend(); ax.set_title("Pseudo-additivity is the corner term")
save(fig, "fig3-tsallis-renyi.png")

# ---------- fig4 C4 Sugeno ----------
from scipy import optimize
def lam_of(g):
    f = lambda l: np.prod(1 + l * g) - (1 + l); s = g.sum()
    if abs(s - 1) < 1e-12: return 0.0
    return optimize.brentq(f, 1e-12, 1e7) if s < 1 else optimize.brentq(f, -1 + 1e-13, -1e-13)
ss = np.linspace(0.15, 2.9, 200); lams = [lam_of(np.full(3, s / 3)) for s in ss]
fig, axs = plt.subplots(1, 2, figsize=(12.5, 4.8))
ax = axs[0]; ax.plot(ss, lams, color=ORANGE, lw=2.5); ax.axhline(0, color=GREY); ax.axhline(-1, color=RED, ls=":")
ax.axvline(1, color=GREY, ls=":")
ax.text(0.2, 6, "Σg < 1: k > 0\nbelief side", fontsize=12); ax.text(1.5, 3, "Σg > 1: −1 < k < 0\nplausibility side", fontsize=12)
ax.text(1.05, -0.8, "never reaches −1", color=RED, fontsize=11)
ax.set_yscale("symlog", linthresh=1); ax.set_xlabel("sum of the three point weights g"); ax.set_ylabel("k = λ (pinned)"); ax.set_title("Normalisation pins Sugeno's k")
ax = axs[1]; ax.axhline(0, color="k"); ax.set_xlim(-3.2, 2.2); ax.set_ylim(-1, 1); ax.set_yticks([])
ax.axvspan(-3.2, -1, color=ORANGE, alpha=0.15); ax.axvspan(-1, 2.2, color=BLUE, alpha=0.10)
ax.axvline(-1, color="k", lw=2)
ax.text(-3.1, 0.45, "k ≤ −1: an 'or' of chances\nwall −1/k = cure share π ≤ 1\n(cure model; alphabets at −q/(q−1))", fontsize=11.5)
ax.text(-0.85, 0.45, "k > −1: Sugeno measures\nwall outside [0, 1]\n(probability at k = 0)", fontsize=11.5)
ax.text(-1, -0.55, "k = −1: both\n(union rule)", ha="center", fontsize=12, weight="bold")
ax.set_xlabel("homothety k, for values in [0, 1]"); ax.set_title("On [0, 1] the line splits at k = −1")
save(fig, "fig4-sugeno-split.png")

# ---------- fig5 C5 calculus ----------
fig, axs = plt.subplots(1, 2, figsize=(12, 4.6))
u = np.linspace(0, 2, 200)
ax = axs[0]
for kv, col, lab in ((-1, GREY, "k = −1: 1 − e^(−u), a CDF from its hazard"), (0, ORANGE, "k = 0: u, the ordinary integral"), (1, BLUE, "k = 1: e^u − 1, the product integral − 1")):
    ax.plot(u, (np.exp(kv * u) - 1) / kv if kv else u, color=col, lw=2.5, label=lab)
ax.set_xlabel("ordinary integral u = ∫ f"); ax.set_ylabel("fused integral E_k(u)"); ax.legend(fontsize=10.5, loc="upper left"); ax.set_title("One integrand, three calculi")
ax = axs[1]; kk = np.linspace(-0.9, 2, 300)
ax.plot(kk, (1 + kk) ** (1 / kk), color=BLUE, lw=2.5); ax.axhline(math.e, color=ORANGE, ls="--"); ax.text(1.2, 2.8, "e", color=ORANGE, fontsize=16)
ax.plot([1], [2], "o", color=GREY, ms=9); ax.text(1.05, 1.85, "k = 1: one step, 2", fontsize=11)
ax.set_xlabel("grain k (step length)"); ax.set_ylabel("(1 + k)^(1/k)"); ax.set_title("Shrinking the grain: compound → continuous")
save(fig, "fig5-calculus-grain.png")

# ---------- fig6 C6 primes ----------
PR = np.array([p for p in range(2, 200000) if all(p % d for d in range(2, int(p ** 0.5) + 1))]) if False else None
def primes_upto(N):
    s = np.ones(N + 1, bool); s[:2] = False
    for i in range(2, int(N ** 0.5) + 1):
        if s[i]: s[i * i::i] = False
    return np.nonzero(s)[0].astype(float)
PR = primes_upto(300000)
kk6 = np.linspace(-1.5, 3.5, 200)
Z2 = [np.exp(np.sum(-np.log1p(-k * PR ** -2.0) / k)) if abs(k) > 1e-9 else np.exp(np.sum(PR ** -2.0)) for k in kk6]
def mob(n):
    r, m, p = 1, n, 2
    while p * p <= m:
        if m % p == 0:
            m //= p
            if m % p == 0: return 0
            r = -r
        p += 1
    return -r if m > 1 else r
M = lambda k, N: sum(mob(d) * k ** (N // d) for d in range(1, N + 1) if N % d == 0) / N
fig, axs = plt.subplots(1, 2, figsize=(12.5, 4.6))
ax = axs[0]; ax.plot(kk6, Z2, color=BLUE, lw=2.5)
for kv, lab in ((-1, "ζ(2)/ζ(4) = 15/π²"), (1, "ζ(2) = π²/6"), (0, "exp P(2)")):
    zv = np.interp(kv, kk6, Z2); ax.plot(kv, zv, "o", color=ORANGE, ms=9); ax.annotate(lab, (kv, zv), (kv + 0.15, zv - 0.25), fontsize=11)
ax.set_xlabel("homothety k"); ax.set_ylabel("Z_k(2)"); ax.set_title("Prime gas at s = 2 along the dial")
ax = axs[1]; Ns = np.arange(1, 9)
for j, (kv, col) in enumerate(((-1, ORANGE), (-0.5, GREY), (0.5, "#bbbbbb"), (1, BLUE))):
    ax.bar(Ns + (j - 1.5) * 0.2, [M(kv, N) / kv for N in Ns], 0.2, color=col, label=f"k = {kv}")
ax.axhline(0, color="k", lw=0.8); ax.set_xlabel("N (zeta factor ζ(Ns))"); ax.set_ylabel("exponent M(k,N)/k"); ax.legend()
ax.set_title("Finitely many zeta factors only at k = ±1")
save(fig, "fig6-primes-wall.png")

# ---------- fig7 C7 tunnelling ----------
fig, axs = plt.subplots(1, 2, figsize=(12, 4.6))
ax = axs[0]; Kk = np.linspace(-2, 2, 200)
ax.plot(Kk, -2 * Kk, color=BLUE, lw=2.5)
for K, lab in ((1, "sphere: mixers, bit flips\npiling-up lemma (k = −2)"), (0, "flat: variances add (k = 0)"), (-1, "hyperbolic: barriers\nmean resistance (k = +2)")):
    ax.plot(K, -2 * K, "o", color=ORANGE, ms=10); ax.annotate(lab, (K, -2 * K), (K - 0.9 if K > 0 else K + 0.1, -2 * K + (0.6 if K <= 0 else 0.8)), fontsize=11)
ax.set_xlabel("curvature K of the background"); ax.set_ylabel("homothety k of the averaged law"); ax.set_title("Averaging over the angle: k = −2K")
ax = axs[1]
D = np.linspace(0, 3, 200)
for K, col in ((1, ORANGE), (0, GREY), (-1, BLUE)):
    t = (1 - np.cos(np.sqrt(K) * D)) / (2 * K) if K > 0 else ((1 - np.cosh(np.sqrt(-K) * D)) / (2 * K) if K < 0 else D ** 2 / 4)
    ax.plot(D, t, color=col, lw=2.5, label=f"K = {K}: wall {'at t = 1/2' if K > 0 else 'none'}")
ax.axhline(0.5, color=ORANGE, ls=":"); ax.set_ylim(0, 3); ax.set_xlabel("step length D"); ax.set_ylabel("squared half-chord t"); ax.legend()
ax.set_title("The k < 0 side has a wall, the k > 0 side has none")
save(fig, "fig7-tunnelling-k.png")

# ---------- fig8 C8 Lee-Yang on graphs ----------
ring = [(i, (i + 1) % 10) for i in range(10)]
def gas(kv, edges=ring, nv=10):
    cf = np.zeros(nv + 1)
    for S in range(1 << nv):
        m = bin(S).count("1"); e = sum(1 for (i, j) in edges if (S >> i) & 1 and (S >> j) & 1); cf[m] += kv ** e
    return np.roots(cf[::-1])
fig, axs = plt.subplots(1, 2, figsize=(12, 5.2))
ax = axs[0]; th = np.linspace(0, 2 * np.pi, 300); ax.plot(np.cos(th), np.sin(th), color=GREY, lw=1)
for kv, col, mk in ((0.0, RED, "x"), (0.6, "#e67e22", "s"), (2.0, BLUE, "o"), (4.0, ORANGE, "^")):
    r = gas(kv) * kv if kv > 0 else gas(kv)
    ax.plot(r.real, r.imag, mk, color=col, ms=8, label=f"k = {kv:g}" + (" (hard core)" if kv == 0 else ""))
ax.set_aspect("equal"); ax.set_xlim(-3.2, 1.6); ax.legend(fontsize=10.5, loc="lower left"); ax.set_title("Zeros × k on a ring of 10 (circle = Lee–Yang)")
ax.set_xlabel("Re (k z)"); ax.set_ylabel("Im (k z)")
ax = axs[1]; kk8 = np.concatenate([np.linspace(0.05, 0.98, 25), np.linspace(1.02, 4, 25)])
dev = [np.max(np.abs(np.abs(gas(kv)) * kv - 1)) for kv in kk8]
ax.semilogy(kk8, np.maximum(dev, 1e-16), "o-", color=BLUE); ax.axvline(1, color=ORANGE, lw=2); ax.text(1.1, 1, "k ≥ 1: on the circle\n(Lee & Yang 1952)", color=ORANGE, fontsize=12)
ax.set_xlabel("edge weight k"); ax.set_ylabel("largest distance off the circle"); ax.set_title("On a graph, k is not a unit")
save(fig, "fig8-lee-yang-graph.png")

# ---------- fig9 C9 counting ----------
fig, axs = plt.subplots(1, 2, figsize=(12.5, 4.8))
ax = axs[0]; kvals = [-2, -1, 0, 1, 2, 3]; nvals = range(1, 8)
tab = np.array([[sum((1 + kv) ** j for j in range(n)) for n in nvals] for kv in kvals], float)
im = ax.imshow(np.log10(np.abs(tab) + 1), cmap="Blues", aspect="auto")
for i in range(len(kvals)):
    for j in range(len(nvals)): ax.text(j, i, int(tab[i, j]), ha="center", va="center", fontsize=10.5)
ax.set_yticks(range(len(kvals))); ax.set_yticklabels([f"k={kv}  " + s for kv, s in zip(kvals, ["parity", "OR", "count", "subsets", "q=3", "q=4"])])
ax.set_xticks(range(7)); ax.set_xticklabels(list(nvals)); ax.set_xlabel("n trues folded: [n] at q = 1 + k"); ax.set_title("One fold, a counting class at each stop")
ax = axs[1]; kk9 = np.linspace(-3, 1, 300); ax.plot(kk9, (kk9 + 1) * (kk9 + 2), color=BLUE, lw=2.5); ax.axhline(0, color="k", lw=0.8)
ax.plot([-2, -1], [0, 0], "o", color=ORANGE, ms=11); ax.text(-2.05, 0.4, "XOR", ha="center", color=ORANGE); ax.text(-0.95, 0.4, "OR", ha="center", color=ORANGE)
ax.set_xlabel("homothety k"); ax.set_ylabel("(k+1)(k+2)"); ax.set_title("Bits stay bits only where this is 0")
save(fig, "fig9-counting-stops.png")

# ---------- fig10 C10 alphabet ladder ----------
fig, axs = plt.subplots(1, 2, figsize=(12.5, 4.8))
ax = axs[0]; ax.axhline(0, color="k")
qsv = list(range(2, 12)); kq = [-q / (q - 1) for q in qsv]
ax.plot(kq, [0] * len(kq), "o", color=ORANGE, ms=10)
for q, kv in zip(qsv[:4], kq[:4]):
    ax.text(kv, 0.12, f"q={q}", ha="center", fontsize=11)
ax.text(-2, -0.3, "bits: Haldane map,\nbit flips, piling-up", ha="center", fontsize=10.5)
ax.text(-4 / 3, -0.3, "DNA (Jukes–Cantor),\nqubit depolarizing", ha="center", fontsize=10.5)
ax.text(-1.02, 0.25, "q→∞: erasure,\nunion (−1)", ha="center", fontsize=10.5)
ax.plot(np.linspace(-2.2, -1, 100), [-0.6] * 100, color=BLUE, lw=5, alpha=0.6); ax.text(-2.2, -0.75, "cure model k = −1/π fills the gaps (any share π)", color=BLUE, fontsize=10.5)
ax.set_xlim(-2.3, -0.85); ax.set_ylim(-0.95, 0.5); ax.set_yticks([]); ax.set_xlabel("homothety k"); ax.set_title("Alphabet stops on [−2, −1]")
ax = axs[1]
for q, col in ((2, ORANGE), (4, BLUE), (20, GREY)):
    kv = -q / (q - 1); e = 0.1; ser = [0]
    for _ in range(40): ser.append(circ(kv, ser[-1], e))
    ax.plot(ser, color=col, lw=2.5, label=f"q = {q}: wall {(q - 1) / q:.2f}"); ax.axhline((q - 1) / q, color=col, ls=":")
ax.set_xlabel("channels in a row (each 10% error)"); ax.set_ylabel("error rate"); ax.legend(); ax.set_title("Errors pile up to the wall (q−1)/q")
save(fig, "fig10-alphabet-ladder.png")

# ---------- fig11 C11 holographic ----------
ks11 = np.linspace(-1.6, 0.6, 300); kst = math.sqrt(6) - 3
fig, ax = plt.subplots(figsize=(10, 5))
ax.plot(ks11, 3 + 6 * ks11 + ks11 ** 2, color=BLUE, lw=3, label="perfect tensor (star graph): 3 + 6k + k²")
ax.plot(ks11, 3 * (1 + ks11) ** 2, color=GREY, lw=2, label="three unentangled regions: 3(1+k)²")
ax.plot(ks11, 6 * (1 + ks11), color=GREY, lw=2, ls="--", label="three shared Bell pairs: 6(1+k)")
ax.axhline(0, color="k", lw=0.8); ax.axvline(-1, color=ORANGE, lw=2); ax.axvline(kst, color=BLUE, ls=":", lw=2); ax.axvline(-0.5, color=GREY, ls=":", lw=1)
ax.text(-1.0, 7.5, "k = −1: ≤ 0\n(MMI, known)", color=ORANGE, ha="center", fontsize=12)
ax.text(kst + 0.03, 7.5, "k* = √6 − 3 ≈ −0.551:\n≥ 0 from here up", color=BLUE, fontsize=12)
ax.text(-0.48, -2.5, "−1/2 (old bound)", color=GREY, fontsize=10)
ax.set_ylim(-3.5, 9.5); ax.set_xlabel("homothety k (corner weight k^(rank−1))"); ax.set_ylabel("corner sum of entropies"); ax.legend(loc="lower right", fontsize=10.5)
ax.set_title("Holographic entropies: where the corner sum has a fixed sign")
save(fig, "fig11-holographic-threshold.png")

# ---------- fig13 N: three whole-number ladders, one order-3 map ----------
fig, ax = plt.subplots(figsize=(12, 4.6))
ax.axhline(0, color="k"); ax.set_xlim(-2.3, 4.5); ax.set_ylim(-1.3, 1.3); ax.set_yticks([])
for m, col in ((2, BLUE), (3, ORANGE), (4, "#6aa84f")):
    pts = [-1 / m, -m / (m - 1), m - 1]
    ax.plot(pts, [0, 0, 0], "o", color=col, ms=11)
    for (x1, x2), rad in zip(((pts[0], pts[1]), (pts[1], pts[2]), (pts[2], pts[0])), (0.4, 0.3, 0.3)):
        ax.annotate("", xy=(x2, 0.02), xytext=(x1, 0.02), arrowprops=dict(arrowstyle="->", color=col, connectionstyle=f"arc3,rad={-rad if x2 > x1 else rad}"))
    ax.text(m - 1, -0.3, f"m = {m}", ha="center", color=col)
ax.text(-2.25, 1.05, "k ↦ −1/(1+k): copies −1/m → alphabet −m/(m−1) → q-numbers m − 1 → copies", fontsize=12.5)
ax.text(-2.25, -1.15, "Checked algebra; no shared meaning found (N1). A curiosity, not a lead.", fontsize=12, color=GREY)
ax.set_xlabel("homothety k"); ax.set_title("Three whole-number ladders are one orbit of an order-3 map")
save(fig, "fig13-three-ladders-orbit.png")
