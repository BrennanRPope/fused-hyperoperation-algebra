"""Figures for the order-dial survey. Run: python3 figures.py -> figures/*.png
Colour key: blue = our reading, orange = already known, grey = endpoints/background."""
import itertools, json, math, cmath
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon, Circle, FancyArrowPatch
import sympy as sp
from qlib import q, gauss, qint, qfact, inv, words

BLUE, ORANGE, GREY, DARK, RED = "#1f6fd1", "#e8840c", "#8c8c8c", "#222222", "#c0392b"
plt.rcParams.update({"font.size": 14, "axes.titlesize": 15, "axes.labelsize": 14, "legend.fontsize": 12,
                     "xtick.labelsize": 12, "ytick.labelsize": 12, "figure.dpi": 130})


def save(fig, name):
    fig.tight_layout()
    fig.savefig("figures/" + name, dpi=130)
    plt.close(fig)


def gnum(n, m, z):
    return complex(sp.Poly(gauss(n, m), q).as_expr().subs(q, z)) if not isinstance(z, (int, float)) else float(gauss(n, m).subs(q, z))


def gpoly(n, m):
    return [int(c) for c in sp.Poly(gauss(n, m), q).all_coeffs()[::-1]]


def geval(coeffs, z):
    return sum(c * z ** j for j, c in enumerate(coeffs))


# ---------------- fig0 hero: the q plane with its proved stops ----------------
fig, ax = plt.subplots(figsize=(10, 8))
th = np.linspace(0, 2 * np.pi, 400)
ax.plot(np.cos(th), np.sin(th), color=BLUE, lw=2.5)
ax.axhline(0, color=GREY, lw=1); ax.axvline(0, color=GREY, lw=0.6)
ax.plot([0, 4.6], [0, 0], color=ORANGE, lw=3, alpha=0.5)
for m in range(3, 9):
    for j in range(1, m):
        if math.gcd(j, m) == 1:
            z = cmath.exp(2j * math.pi * j / m)
            ax.plot(z.real, z.imag, "o", color=ORANGE, ms=6)
stops = [(0, 0, "q = 0\nonly the sorted word\n(one bag; the minimum)", (-2.95, -1.25)),
         (1, 0, "q = 1\ncommutative:\nbinomial counts", (1.35, 0.25)),
         (-1, 0, "q = -1 = e^{iπ}\nanticommute: Pauli,\nClifford, at most once", (-2.95, 0.35))]
for x, y, lab, (tx, ty) in stops:
    ax.plot(x, y, "o", color=DARK, ms=12, zorder=5)
    ax.annotate(lab, (x, y), (tx, ty), fontsize=12.5, color=DARK, arrowprops=dict(arrowstyle="-", color=GREY))
for n_ in (2, 3, 4, 5):
    ax.plot(n_, 0, "s", color=ORANGE, ms=8)
ax.text(2.0, -0.42, "q = 2, 3, 4, 5: subspace counts over\nthe field with q elements", color=ORANGE, fontsize=12)
ax.annotate("q → ∞: only the reversed word\n(the maximum area)", (4.6, 0), (2.9, 0.75), fontsize=12.5,
            arrowprops=dict(arrowstyle="->", color=GREY))
ax.text(-1.2, 1.15, "unit circle q = e^{iħ}: Weyl and Moyal,\nphase = area / ħ", color=BLUE, fontsize=12.5)
ax.annotate("roots of unity (orange dots):\nclock-and-shift matrices,\n(x+y)^m = x^m + y^m, q-Lucas,\nTuraev–Viro truncation",
            (math.cos(2 * math.pi / 5), -math.sin(2 * math.pi / 5)), (0.6, -1.85), fontsize=12, color=ORANGE,
            arrowprops=dict(arrowstyle="-", color=ORANGE))
ax.set_xlim(-3.1, 4.9); ax.set_ylim(-2.05, 1.65); ax.set_aspect("equal")
ax.set_xlabel("real part of q"); ax.set_ylabel("imaginary part of q")
ax.set_title("The order dial q: yx = q·xy. Black = proved stops")
save(fig, "fig0-order-dial-stops.png")

# ---------------- sorting board ----------------
rows = json.load(open("laws.json"))["rows"]
order = ["slides", "plausible", "partly", "disguise", "self", "no"]
cols = {"slides": BLUE, "plausible": BLUE, "partly": ORANGE, "disguise": GREY, "self": DARK, "no": GREY}
fig, ax = plt.subplots(figsize=(11, 11))
y = 0.0
short = lambda s: (s[:44] + "…") if len(s) > 45 else s
for v in order:
    rs = [r for r in rows if r["verdict"] == v]
    ax.text(0.0, y, f"{v} ({len(rs)})", fontsize=15, fontweight="bold", color=cols[v], va="top")
    y -= 0.045
    if v == "no":
        line = ""
        nums = [str(r["n"]) for r in rs]
        for i in range(0, len(nums), 16):
            ax.text(0.02, y, "  ".join(nums[i:i + 16]), fontsize=12.5, color=GREY, va="top")
            y -= 0.04
    else:
        for i, r in enumerate(rs):
            cx = 0.02 if i % 2 == 0 else 0.52
            ax.text(cx, y, f"{r['n']}. {short(r['law'])}  [{r['candidate']}]", fontsize=11.5, color=cols[v], va="top")
            if i % 2 == 1 or i == len(rs) - 1:
                y -= 0.04
    y -= 0.025
ax.set_xlim(0, 1); ax.set_ylim(y, 0.03); ax.axis("off")
ax.set_title("All 65 law groups, sorted by what turning q does\n(the 48 'no' groups by number)")
save(fig, "fig13-sorting-board.png")

# ---------------- C1 ----------------
fig, axs = plt.subplots(1, 2, figsize=(11, 4.6))
qq = np.linspace(-1.5, 2.5, 300)
axs[0].plot(qq, 1 - qq, color=BLUE, lw=2.5)
axs[0].axhline(0, color=GREY)
for x, lab in [(1, "q = 1: order free"), (-1, "q = -1: gap doubles,\n(x+y)² = x² + y²")]:
    axs[0].plot(x, 1 - x, "o", color=DARK, ms=10)
    axs[0].annotate(lab, (x, 1 - x), (x + 0.1, 1 - x + 0.5), fontsize=12)
axs[0].set_xlabel("q (real)"); axs[0].set_ylabel("order gap / (k·αβ·xy)")
axs[0].set_title("a ∘ₖ b − b ∘ₖ a = k(1 − q)·ab")
hh = np.linspace(0, 2 * np.pi, 300)
axs[1].plot(hh, 2 * np.abs(np.sin(hh / 2)), color=BLUE, lw=2.5)
axs[1].set_xticks([0, np.pi / 2, np.pi, 3 * np.pi / 2, 2 * np.pi]); axs[1].set_xticklabels(["0", "π/2", "π", "3π/2", "2π"])
axs[1].plot([0, np.pi], [0, 2], "o", color=DARK, ms=10)
axs[1].set_xlabel("ħ, with q = e^{iħ}"); axs[1].set_ylabel("size of the gap |1 − q|")
axs[1].set_title("On the circle: zero at ħ = 0, largest at ħ = π")
save(fig, "fig1-core-law-order-gap.png")

# ---------------- C2: corner weights in reversed order ----------------
fig, ax = plt.subplots(figsize=(8.5, 6.5))
P = {(0, 0, 0): (0, 0), (1, 0, 0): (3, 0), (0, 1, 0): (1.2, 1.1), (0, 0, 1): (0, 3)}
def pt(c):
    return (3 * c[0] + 1.2 * c[1], 1.1 * c[1] + 3 * c[2])
corners = list(itertools.product([0, 1], repeat=3))
for a in corners:
    for b in corners:
        if sum(abs(i - j) for i, j in zip(a, b)) == 1 and a < b:
            (x1, y1), (x2, y2) = pt(a), pt(b)
            ax.plot([x1, x2], [y1, y2], color=GREY, lw=1.5)
labs = {0: "1", 1: "k", 2: "k² · q", 3: "k³ · q³"}
for c in corners:
    x, y = pt(c); m = sum(c)
    ax.plot(x, y, "o", color=BLUE if m >= 2 else GREY, ms=14)
    ax.text(x + 0.15, y + 0.15, labs[m], fontsize=14, color=BLUE if m >= 2 else DARK)
ax.text(-0.3, -1.3, "Chain run in reversed order over x_j x_i = q·x_i x_j: a corner of rank m picks up\n"
        "q^C(m,2), one q per pair. That is the soliton corner A^C(m,2) with A = q (orange).",
        fontsize=12, color=DARK)
ax.text(-0.3, 4.5, "soliton tau function: 1 + Σ A^C(m,2)·(corner)", color=ORANGE, fontsize=13)
ax.set_xlim(-0.6, 5.6); ax.set_ylim(-1.7, 4.9); ax.axis("off")
ax.set_title("Corner expansion slides into soliton corners as q turns")
save(fig, "fig2-corners-to-solitons.png")

# ---------------- C3: q-binomial row and q-Lucas zeros ----------------
fig, axs = plt.subplots(1, 2, figsize=(12, 5))
n = 8
for qv, col, lab in [(0.0, GREY, "q = 0: every word once"), (0.5, BLUE, "q = 0.5"), (1.0, ORANGE, "q = 1: Pascal")]:
    axs[0].plot(range(n + 1), [float(geval(gpoly(n, m), qv)) for m in range(n + 1)], "o-", color=col, lw=2, label=lab, ms=7)
axs[0].set_yscale("log"); axs[0].set_xlabel("m (number of y's)"); axs[0].set_ylabel("[8, m]_q")
axs[0].set_title("Row 8 of the triangle as q turns"); axs[0].legend()
N = 21; w = cmath.exp(2j * math.pi / 3)
M = np.full((N, N), np.nan)
for nn in range(N):
    for m in range(nn + 1):
        M[nn, m] = abs(geval(gpoly(nn, m), w))
im = axs[1].imshow(np.where(M < 1e-9, -1, M), cmap="Blues", origin="upper", vmin=-1, vmax=np.nanmax(M))
zs = np.argwhere(M < 1e-9)
axs[1].plot(zs[:, 1], zs[:, 0], "s", color=ORANGE, ms=3.2, label="exactly 0")
axs[1].set_xlabel("m"); axs[1].set_ylabel("n")
axs[1].set_title("|[n, m]_q| at q = e^{2πi/3}\n(q-Lucas: blocks of 3; orange = 0)"); axs[1].legend(loc="lower left")
save(fig, "fig3-binomial-to-q-binomial.png")

# ---------------- C4: q-integers ----------------
fig, axs = plt.subplots(1, 2, figsize=(12, 4.8))
qq = np.linspace(-1.3, 2.05, 300)
for nn in range(1, 7):
    axs[0].plot(qq, sum(qq ** i for i in range(nn)), lw=2, color=plt.cm.Blues(0.35 + 0.1 * nn), label=f"[{nn}]_q")
for x, lab in [(0, "q=0: 1"), (1, "q=1: n"), (-1, "q=-1: parity"), (2, "q=2: 2ⁿ−1")]:
    axs[0].axvline(x, color=GREY, ls=":", lw=1)
    axs[0].text(x + 0.03, 33, lab, fontsize=11, color=DARK)
axs[0].set_ylim(-1.5, 36); axs[0].set_xlabel("q (real)"); axs[0].set_ylabel("[n]_q"); axs[0].legend(ncol=2, fontsize=10)
axs[0].set_title("q-counting: 1 + q + … + q^(n−1)")
hh = np.linspace(0, 2 * np.pi, 600)
for nn, col in [(4, BLUE), (6, ORANGE)]:
    axs[1].plot(hh, np.abs(sum(np.exp(1j * hh * i) for i in range(nn))), color=col, lw=2.2, label=f"|[{nn}]_q|, q = e^{{iħ}}")
axs[1].set_xticks([0, np.pi / 2, np.pi, 3 * np.pi / 2, 2 * np.pi]); axs[1].set_xticklabels(["0", "π/2", "π", "3π/2", "2π"])
axs[1].set_xlabel("ħ"); axs[1].set_title("On the circle: zeros at ħ = 2πj/n"); axs[1].legend()
save(fig, "fig4-q-counting.png")

# ---------------- C5: face weights ----------------
fig, axs = plt.subplots(1, 2, figsize=(12, 5))
word = (0, 1, 1, 0, 1, 0, 0)  # 0 = x (right), 1 = y (up); 3 y's, 4 x's
x, y = 0, 0; path = [(0, 0)]
for s in word:
    if s == 0: x += 1
    else: y += 1
    path.append((x, y))
ax = axs[0]
for i in range(5): ax.plot([i, i], [0, 3], color="#dddddd")
for j in range(4): ax.plot([0, 4], [j, j], color="#dddddd")
# boxes above the path and below the top: area = inversions
ys = {}
for (x1, y1), (x2, y2) in zip(path, path[1:]):
    if y2 == y1: ys[x1] = y1
for col in range(4):
    for row in range(ys[col], 3):
        ax.add_patch(Polygon([(col, row), (col + 1, row), (col + 1, row + 1), (col, row + 1)], color=BLUE, alpha=0.35))
px, py = zip(*path)
ax.plot(px, py, color=DARK, lw=3)
ax.set_title(f"word x y y x y x x: {inv(word)} out-of-order pairs\n= {inv(word)} shaded boxes (weight q^{inv(word)})")
ax.set_aspect("equal"); ax.axis("off")
n, m = 10, 5
cnt = {}
for w_ in words(n, m):
    cnt[inv(w_)] = cnt.get(inv(w_), 0) + 1
A = np.array(sorted(cnt)); c = np.array([cnt[a] for a in A], dtype=float)
for qv, col, lab in [(0.4, GREY, "q = 0.4: piles at the minimum"), (1.0, ORANGE, "q = 1: plain count"), (2.5, BLUE, "q = 2.5: piles at the maximum")]:
    wts = c * qv ** A; axs[1].plot(A, wts / wts.sum(), "o-", color=col, label=lab, lw=2, ms=5)
axs[1].set_xlabel("area of the path (0 to 25)"); axs[1].set_ylabel("share of the total weight")
axs[1].set_title("Paths in a 5 × 5 box, weighted q^area"); axs[1].legend(fontsize=11)
save(fig, "fig5-face-weights-q-area.png")

# ---------------- C6: routes on the permutohedron ----------------
fig, axs = plt.subplots(1, 2, figsize=(12, 5))
ax = axs[0]
perms3 = list(itertools.permutations([1, 2, 3]))
ang = {(1, 2, 3): 270, (2, 1, 3): 210, (1, 3, 2): 330, (2, 3, 1): 150, (3, 1, 2): 30, (3, 2, 1): 90}
pos = {p: (math.cos(math.radians(a)), math.sin(math.radians(a))) for p, a in ang.items()}
for a in perms3:
    for b in perms3:
        if a < b and sum(1 for i in range(2) if (a[i], a[i + 1]) == (b[i + 1], b[i]) and a[:i] == b[:i] and a[i + 2:] == b[i + 2:]):
            ax.plot(*zip(pos[a], pos[b]), color=GREY, lw=2)
for p in perms3:
    x, y = pos[p]
    ax.plot(x, y, "o", ms=16, color=BLUE)
    ax.text(x * 1.42, y * 1.42 - (0.25 if y < 0 else 0.05), "".join(map(str, p)) + f"\nq^{inv(p)}", ha="center", fontsize=13)
ax.text(0, -0.15, "1 + 2q + 2q² + q³\n= [3]_q!", ha="center", fontsize=13, color=ORANGE)
ax.set_xlim(-1.7, 1.7); ax.set_ylim(-1.7, 1.7); ax.set_aspect("equal"); ax.axis("off")
ax.set_title("Orders of 3 things, weighted q^(out-of-order pairs)")
cnt = {}
for p in itertools.permutations(range(6)):
    cnt[inv(p)] = cnt.get(inv(p), 0) + 1
A = np.array(sorted(cnt)); c = np.array([cnt[a] for a in A], dtype=float)
for qv, col, lab in [(0.5, GREY, "q = 0.5"), (1.0, ORANGE, "q = 1 (all 720 alike)"), (2.0, BLUE, "q = 2")]:
    wts = c * qv ** A; axs[1].plot(A, wts / wts.sum(), "o-", color=col, label=lab, lw=2, ms=5)
axs[1].set_xlabel("out-of-order pairs (0 to 15)"); axs[1].set_ylabel("share of weight"); axs[1].legend()
axs[1].set_title("720 orders of 6 things: the weight slides\nfrom sorted (q→0) to reversed (q→∞)")
save(fig, "fig6-permutohedron-routes.png")

# ---------------- C7: Rulecraft with Mallows weights ----------------
def outcome(order, lam, tgt, x0=0.5):
    x = x0
    for i in order:
        x = (1 - lam[i]) * x + lam[i] * tgt[i]
    return x


rng = np.random.default_rng(3)
n = 5
perms = list(itertools.permutations(range(n))); invs = np.array([inv(p) for p in perms])
qg = np.logspace(-2, 2, 200)
fig, axs = plt.subplots(1, 2, figsize=(12, 5))
cols_ = [BLUE, ORANGE, DARK]
for j in range(3):
    lam = rng.uniform(0.05, 0.9, n); tgt = np.sort(rng.uniform(0, 1, n))
    outs = np.array([outcome(p, lam, tgt) for p in perms])
    means = [np.dot(qv ** invs, outs) / np.sum(qv ** invs) for qv in qg]
    axs[0].plot(qg, means, color=cols_[j], lw=2.3, label=f"package {j + 1}")
    axs[0].plot([qg[0], qg[-1]], [outs.max(), outs.min()], "o", color=cols_[j], ms=7)
    axs[0].plot(1, outs.mean(), "s", color=GREY, ms=8)
axs[0].set_xscale("log"); axs[0].axvline(1, color=GREY, ls=":")
axs[0].set_xlabel("q (weight per out-of-order pair)"); axs[0].set_ylabel("expected slider after 5 rules")
axs[0].set_title("Best order (q→0) → plain average (q=1, grey)\n→ worst order (q→∞)"); axs[0].legend(fontsize=11)
n4 = list(itertools.permutations(range(4))); iv4 = np.array([inv(p) for p in n4])
for qv, col, lab in [(0.3, BLUE, "q = 0.3"), (1.0, ORANGE, "q = 1"), (3.0, DARK, "q = 3")]:
    w_ = qv ** iv4; w_ = w_ / w_.sum()
    jit = np.array([sum(1 for kk in range(i) if iv4[kk] == iv4[i]) for i in range(24)])
    axs[1].scatter(iv4 + 0.12 * jit - 0.3, w_, color=col, s=40, label=lab)
axs[1].set_xlabel("swaps away from 'sorted by target'")
axs[1].set_ylabel("weight of each of the 24 orders"); axs[1].legend()
axs[1].set_title("Mallows weights on the ladder of orders")
save(fig, "fig7-rulecraft-mallows.png")

# ---------------- C8: Weyl / Moyal ----------------
fig, axs = plt.subplots(1, 2, figsize=(12, 5))
ax = axs[0]
a_ = np.array([2.0, 0.5]); b_ = np.array([0.6, 1.6])
ax.add_patch(Polygon([(0, 0), tuple(a_), tuple(a_ + b_), tuple(b_)], color=BLUE, alpha=0.3))
for v, lab, col in [(a_, "move a (shift)", DARK), (b_, "move b (boost)", DARK)]:
    ax.annotate("", v, (0, 0), arrowprops=dict(arrowstyle="->", lw=2.5, color=col))
    ax.text(v[0] * 0.5 + 0.05, v[1] * 0.5 - 0.25 if v is a_ else v[1] * 0.5, lab, fontsize=12)
ax.annotate("", a_ + b_, a_, arrowprops=dict(arrowstyle="->", lw=1.5, color=GREY))
ax.annotate("", a_ + b_, b_, arrowprops=dict(arrowstyle="->", lw=1.5, color=GREY))
area = a_[0] * b_[1] - a_[1] * b_[0]
ax.text(0.7, 0.95, f"area ω = {area:.2f}", fontsize=14, color=BLUE)
ax.text(0.0, 2.45, "two orders end at the same point;\nthe states differ by q^ω = e^{iħω}", fontsize=12)
ax.set_xlim(-0.2, 3.0); ax.set_ylim(-0.3, 2.9); ax.set_aspect("equal"); ax.set_xlabel("position"); ax.set_ylabel("momentum")
ax.set_title("Weyl: the swap costs a phase = area")
ax = axs[1]
ax.plot(np.cos(th), np.sin(th), color=BLUE, lw=2)
for m in (2, 3, 4, 6):
    z = cmath.exp(2j * math.pi / m)
    ax.plot(z.real, z.imag, "o", color=ORANGE, ms=9)
    ax.text(z.real * 1.18 - 0.1, z.imag * 1.18, f"ħ=2π/{m}\n{m}×{m}", fontsize=10.5, color=ORANGE)
ax.plot(1, 0, "o", color=DARK, ms=12); ax.text(1.08, -0.18, "ħ = 0:\ncommutative", fontsize=11)
ax.text(-0.75, -0.3, "other ħ: only\ninfinite matrices\n(Stone–von Neumann)", fontsize=11, color=BLUE)
ax.set_xlim(-1.6, 1.8); ax.set_ylim(-1.5, 1.5); ax.set_aspect("equal"); ax.axis("off")
ax.set_title("q = e^{iħ}: finite clock-and-shift\nonly at roots of unity")
save(fig, "fig8-weyl-moyal.png")

# ---------------- C9: Clifford and repetition ----------------
fig, axs = plt.subplots(1, 2, figsize=(12, 4.8))
ax = axs[0]
ax.plot(np.cos(th), np.sin(th), color="#dddddd", lw=1.5); ax.axhline(0, color="#dddddd")
ax.plot([-1, 1], [0, 0], "o", color=DARK, ms=14)
ax.text(-1.45, 0.15, "q = -1\nClifford,\nPauli", fontsize=12); ax.text(1.05, 0.15, "q = 1\ncommuting", fontsize=12)
ax.text(-0.6, -0.9, "everywhere else: no elements that\nsquare to a number (q² = 1 forced)", fontsize=11.5, color=RED)
ax.set_xlim(-1.7, 1.8); ax.set_ylim(-1.3, 1.3); ax.set_aspect("equal"); ax.axis("off")
ax.set_title("Things that square to numbers")
ax = axs[1]
hh = np.linspace(0, 2 * np.pi, 400)
ax.plot(hh, np.where(np.isclose(hh, 0) | np.isclose(hh, 2 * np.pi), 1, 0), color=BLUE, lw=3)
ax.plot([0, 2 * np.pi], [1, 1], "o", color=DARK, ms=10)
ax.text(0.35, 0.8, "q = 1: repeats allowed (Bose)", fontsize=12)
ax.text(1.2, 0.12, "any q ≠ 1: x·x = 0, at most once (Fermi at q = -1)", fontsize=12, color=BLUE)
ax.set_xticks([0, np.pi, 2 * np.pi]); ax.set_xticklabels(["0", "π", "2π"]); ax.set_yticks([0, 1]); ax.set_yticklabels(["no", "yes"])
ax.set_xlabel("ħ, with q = e^{iħ}"); ax.set_ylabel("can x repeat?")
ax.set_title("Repeats: a jump, not a slide")
save(fig, "fig9-clifford-and-repetition.png")

# ---------------- C10: knots ----------------
V = {"trefoil 3₁": lambda z: z ** -1 + z ** -3 - z ** -4,
     "figure-eight 4₁": lambda z: z ** 2 - z + 1 - z ** -1 + z ** -2,
     "cinquefoil 5₁": lambda z: z ** -2 + z ** -4 - z ** -5 + z ** -6 - z ** -7}
fig, ax = plt.subplots(figsize=(10, 5))
hh = np.linspace(0, np.pi, 500)
for (name, f), col in zip(V.items(), [BLUE, ORANGE, DARK]):
    ax.plot(hh, np.abs(f(np.exp(1j * hh))), color=col, lw=2.3, label=name)
for x, lab in [(0, "t = 1\ncounts parts: 1\n(Lean)"), (np.pi / 2, "t = i\nArf: |V| = 1"), (2 * np.pi / 3, "t = e^{2πi/3}\nV = 1"),
               (np.pi, "t = -1\ndeterminant\n3, 5, 5 (Lean)")]:
    ax.axvline(x, color=GREY, ls=":")
    ax.text(x + 0.03, 8.0, lab, fontsize=10.5)
ax.set_xticks([0, np.pi / 2, 2 * np.pi / 3, np.pi]); ax.set_xticklabels(["0", "π/2", "2π/3", "π"])
ax.set_xlabel("angle θ, t = e^{iθ}"); ax.set_ylabel("|V(t)|"); ax.set_ylim(0, 10.5); ax.legend(loc="center left", bbox_to_anchor=(0.0, 0.62), fontsize=11)
ax.set_title("Jones polynomial around the circle: Lean proves the two ends")
save(fig, "fig10-jones-on-the-circle.png")

# ---------------- C11: q-exponential ----------------
fig, ax = plt.subplots(figsize=(9, 5))
qq = np.linspace(0, 3, 300)
def eq(z, qv, N=60):
    s, f = 0.0, 1.0
    for nn in range(N):
        if nn > 0:
            f *= sum(qv ** i for i in range(nn))
        s += z ** nn / f
    return s
for z, col in [(0.5, BLUE), (0.8, ORANGE)]:
    ax.plot(qq, [eq(z, qv) for qv in qq], color=col, lw=2.3, label=f"e_q({z})")
    ax.plot([0, 1], [1 / (1 - z), math.exp(z)], "o", color=DARK, ms=8)
ax.text(0.05, 4.7, "q = 0: 1/(1 − z)\n(geometric)", fontsize=12)
ax.text(1.05, 2.3, "q = 1: e^z\n(Lean E: Σ 1/n! = e)", fontsize=12)
ax.set_xlabel("q"); ax.set_ylabel("e_q(z) = Σ zⁿ/[n]_q!"); ax.legend()
ax.set_title("e_q(x + y) = e_q(x)·e_q(y) when yx = q·xy")
save(fig, "fig11-q-exponential.png")

# ---------------- C12: disguises ----------------
fig, axs = plt.subplots(1, 2, figsize=(12, 5))
ax = axs[0]
K = np.linspace(-3, 3, 400); D = np.linspace(-1, 3, 400)
KK, DD = np.meshgrid(K, D)
Delta = KK ** 2 - 4 * DD
ax.contourf(KK, DD, (Delta < 0).astype(float), levels=[-0.5, 0.5, 1.5], colors=["#fbe3c7", "#d4e4fa"])
ax.plot(K, K ** 2 / 4, color=DARK, lw=2.5)
ax.text(-0.6, 2.2, "|q| = 1\nq = e^{iħ}\n(angle side)", fontsize=12, color=BLUE)
ax.text(1.6, 0.0, "q real\n(radius side)", fontsize=12, color=ORANGE)
ax.text(-2.9, 2.5, "parabola k² = 4d:\nq = 1", fontsize=12)
ax.set_xlabel("k (homothety)"); ax.set_ylabel("d (curvature dial)"); ax.set_ylim(-1, 3)
ax.set_title("q = β/α is a reading of (k, d)")
ax = axs[1]
qq = np.linspace(-1.2, 2.2, 100)
ax.plot(qq, qq - 1, color=GREY, lw=2.5)
for x, lab in [(0, "q=0 ↔ k=−1 (or)"), (1, "q=1 ↔ k=0 (adding)"), (2, "q=2 ↔ k=1 (∘₁)"), (-1, "q=−1 ↔ k=−2 (xor)")]:
    ax.plot(x, x - 1, "o", color=DARK, ms=9); ax.text(x + 0.07, x - 1 - 0.3, lab, fontsize=11)
ax.set_xlabel("q in q-counting, toric counts, Jackson"); ax.set_ylabel("k = q − 1")
ax.set_title("Where q is the k dial: k = q − 1")
save(fig, "fig12-disguises.png")

# ---------------- N: in-between laws ----------------
fig, axs = plt.subplots(1, 2, figsize=(12, 5))
hh = np.linspace(0.001, 2 * np.pi - 0.001, 1500)
for m_, col in zip(range(1, 6), plt.cm.Blues(np.linspace(0.4, 1, 5))):
    cf = gpoly(6, m_)
    axs[0].plot(hh, np.abs(geval(cf, np.exp(1j * hh))), color=col, lw=2, label=f"|[6,{m_}]|")
for x in (2 * np.pi / 6, 10 * np.pi / 6):
    axs[0].axvline(x, color=ORANGE, lw=2, ls="--")
axs[0].text(2 * np.pi / 6 + 0.05, 17, "primitive 6th roots:\nall five vanish,\n(x+y)⁶ = x⁶ + y⁶", fontsize=11, color=ORANGE)
axs[0].set_xticks([0, np.pi / 3, np.pi, 5 * np.pi / 3, 2 * np.pi]); axs[0].set_xticklabels(["0", "π/3", "π", "5π/3", "2π"])
axs[0].set_xlabel("ħ, q = e^{iħ}"); axs[0].set_ylabel("|[6, m]_q|"); axs[0].legend(fontsize=10)
axs[0].set_title("The binomial collapses to two terms\nonly at primitive roots")
ax = axs[1]
ms_ = np.arange(2, 13)
ax.bar(ms_, ms_, color=ORANGE, alpha=0.8)
ax.set_xlabel("order m of the root of unity q"); ax.set_ylabel("smallest matrix size")
ax.set_title("Finite matrices with VU = q·UV:\nonly at roots of unity, size = order")
ax.text(2.2, 10.2, "proof: det(VU) = qᵐ·det(UV)\nforces qᵐ = 1", fontsize=12)
save(fig, "figN-in-between-laws.png")
print("figures written")
