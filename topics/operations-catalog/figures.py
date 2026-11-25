"""Pictures for the operations catalog. One picture per operation. Run: python3 -I figures.py (writes figures/*.png)."""
import os, math, json
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch, Polygon
os.makedirs("figures", exist_ok=True)
plt.rcParams.update({"font.size": 11, "axes.spines.top": False, "axes.spines.right": False, "figure.dpi": 110})
BLUE, ORANGE, GREEN, RED, PURPLE, GRAY = "#2b6cb0", "#dd6b20", "#2f855a", "#c53030", "#6b46c1", "#718096"
def save(fig, name):
    fig.tight_layout(); fig.savefig(f"figures/{name}.png", bbox_inches="tight"); plt.close(fig); print("wrote", name)
def arrow(ax, p, q, text="", color=GRAY, rad=0.0, fs=10, off=(0, 0)):
    ax.add_patch(FancyArrowPatch(p, q, arrowstyle="-|>", mutation_scale=16, color=color, lw=1.8, connectionstyle=f"arc3,rad={rad}"))
    if text:
        m = ((p[0] + q[0]) / 2 + off[0], (p[1] + q[1]) / 2 + off[1])
        ax.text(*m, text, ha="center", va="center", fontsize=fs, color=color, bbox=dict(fc="white", ec="none", alpha=0.9))
def box(ax, xy, text, color, fs=11, w=None):
    ax.text(*xy, text, ha="center", va="center", fontsize=fs, bbox=dict(boxstyle="round,pad=0.5", fc="white", ec=color, lw=2))

# fractional iteration of E(x)=e^x-1 by the Levy-type limit (known method; see fractional-rungs/e-every-height)
def Et(t, x, n=60):
    x = float(x)
    if abs(x) < 1e-15: return x
    if x > 0:
        y = x
        for _ in range(n): y = math.log1p(y)
        y = y / (1 - t * y / 2)
        for _ in range(n):
            y = math.expm1(y)
            if y > 1e6: return float("inf")
        return y
    y = x
    for _ in range(n): y = math.expm1(y)
    y = y / (1 - t * y / 2)
    for _ in range(n):
        if y <= -1: return float("nan")
        y = math.log1p(y)
    return y
def Bs(s, a, b): return Et(s, Et(-s, a) + Et(-s, b))

# 01 transport hub
fig, ax = plt.subplots(figsize=(10, 7)); ax.set_xlim(-5, 5); ax.set_ylim(-4, 4); ax.axis("off")
box(ax, (0, 0), "plain addition\na + b", BLUE, 13)
spokes = [((0, 3.2), "a ∘ₖ b = a + b + k·ab", "Lₖ(x) = log(1 + kx)/k", GREEN),
          ((3.6, 1.6), "two-dial law\n(a+b+k·ab)/(1−d·ab)", "L_{k,d}(x) = ∫ dt/(1+kt+dt²)", GREEN),
          ((3.6, -1.9), "height-s operation\nB_s(a, b)", "E^(−s)  (fractional iterate)", PURPLE),
          ((-3.6, -1.9), "multiplication a·b\n(positive numbers)", "log", ORANGE),
          ((-3.6, 1.6), "a ∘₁ b, a+b+1, (1+a)^(1+b) − 1\n(fused rungs)", "shift x ↦ x + 1 on the rungs\n(the lift)", RED),
          ((0, -3.3), "max(a, b)", "limit h → 0 of h·log(e^(a/h) + e^(b/h))\n(Maslov dequantization)", GRAY)]
for (xy, name, f, c) in spokes:
    box(ax, xy, name, c, 10)
    arrow(ax, (xy[0] * 0.72, xy[1] * 0.72), (xy[0] * 0.2, xy[1] * 0.22), f, c, fs=9, off=(0, 0))
ax.set_title("One move does most of the work: transport an operation along a map\nf₊(∗)(a, b) = f( f⁻¹(a) ∗ f⁻¹(b) ).  Each arrow is labelled with the map that straightens that law into +",
             fontsize=12)
save(fig, "fig01-transport-hub")

# 02 lift
fig, axs = plt.subplots(1, 2, figsize=(11, 4.2))
k, a, b = 1.0, 0.5, 0.8; c = a + b + k * a * b
ax = axs[0]; ax.axhline(0, color="k", lw=0.8)
for v, lab, col in [(a, "a", BLUE), (b, "b", BLUE), (c, "a ∘ₖ b", RED), (-1 / k, "wall −1/k", GRAY), (0, "0 (identity)", GREEN)]:
    ax.plot(v, 0, "o", color=col, ms=9); ax.annotate(lab, (v, 0), (v, 0.35 if lab != "0 (identity)" else -0.45), ha="center", color=col)
ax.set_ylim(-1, 1); ax.set_xlim(-1.4, 2.6); ax.set_yticks([]); ax.set_title("on the line, combine with ∘ₖ  (k = 1)")
ax = axs[1]; ax.axhline(0, color="k", lw=0.8)
for v, lab, col in [(1 + k * a, "1 + ka", BLUE), (1 + k * b, "1 + kb", BLUE), ((1 + k * a) * (1 + k * b), "(1 + ka)(1 + kb)", RED), (0, "0 (wall lands here)", GRAY), (1, "1 (identity)", GREEN)]:
    ax.plot(v, 0, "o", color=col, ms=9); ax.annotate(lab, (v, 0), (v, 0.35 if v not in (1,) else -0.45), ha="center", color=col)
ax.set_ylim(-1, 1); ax.set_xlim(-0.4, 3.6); ax.set_yticks([]); ax.set_title("after the lift x ↦ 1 + kx, just multiply")
fig.suptitle("Lift ℓₖ(x) = 1 + kx:  1 + k(a ∘ₖ b) = (1 + ka)(1 + kb).  Known: Jacobson's adjoint-group map a ↦ 1 + a; the growth factor 1 + i in finance", fontsize=11)
save(fig, "fig02-lift")

# 03 wall and way back
fig, ax = plt.subplots(figsize=(7.5, 5))
k = 1.0; x = np.linspace(-0.98, 3, 600); y = -x / (1 + k * x)
ax.plot(x, y, color=BLUE, lw=2.5, label="way back  a′ = −a/(1 + ka)")
x2 = np.linspace(-3, -1.02, 300); ax.plot(x2, -x2 / (1 + k * x2), color=BLUE, lw=2.5)
ax.axvline(-1 / k, color=RED, ls="--", lw=2); ax.text(-1 / k + 0.05, 3.2, "wall −1/k:\nno way back,\nand it swallows\neverything:\nwall ∘ₖ b = wall", color=RED, fontsize=10, va="top")
ax.plot([0], [0], "o", color=GREEN, ms=8); ax.annotate("0 is its own way back", (0, 0), (0.6, 1.2), arrowprops=dict(arrowstyle="->", color=GREEN), color=GREEN)
ax.plot([-2 / k], [-2 / k], "o", color=ORANGE, ms=8); ax.annotate("−2/k: the other self-inverse point", (-2, -2), (-2.9, -3.2), color=ORANGE, arrowprops=dict(arrowstyle="->", color=ORANGE))
ax.set_ylim(-4, 4); ax.set_xlim(-3, 3); ax.axhline(0, color="k", lw=0.5); ax.axvline(0, color="k", lw=0.5)
ax.set_xlabel("a"); ax.set_ylabel("a′"); ax.legend(loc="upper right")
ax.set_title("Way back = quasi-inverse (Perlis, Jacobson).  Wall = absorbing element (it lifts to 0)")
save(fig, "fig03-wall-and-way-back")

# 04 circle power
fig, ax = plt.subplots(figsize=(7.5, 4.8))
k = 1.0; t = np.linspace(0, 5, 300)
for a, col in [(0.1, GREEN), (0.25, BLUE), (0.5, ORANGE)]:
    ax.plot(t, ((1 + k * a) ** t - 1) / k, color=col, lw=2, label=f"a = {a}")
    n = np.arange(0, 6); ax.plot(n, ((1 + k * a) ** n - 1) / k, "o", color=col)
ax.set_xlabel("n  (how many copies of a are combined; dots = whole numbers, curve = real powers)")
ax.set_ylabel("a^(∘n)"); ax.legend()
ax.set_title("Circle power a^(∘ₖ n) = ((1 + ka)ⁿ − 1)/k: compound interest at rate a for n periods")
save(fig, "fig04-circle-power")

# 05 n-fold combine on the cube
fig, axs = plt.subplots(1, 2, figsize=(11, 4.8))
ax = axs[0]; ax.axis("off"); ax.set_aspect("equal")
sq = {(0, 0): "∅ : 0", (1, 0): "a", (0, 1): "b", (1, 1): "k·ab"}
for (i, j), lab in sq.items():
    ax.plot(i, j, "o", ms=14, color=[GRAY, BLUE, BLUE, RED][i + j + (1 if i + j == 2 else 0) - (1 if i + j == 2 else 0)] if False else (GRAY if i + j == 0 else BLUE if i + j == 1 else RED))
    ax.text(i + (0.12 if i else -0.12), j + (0.12 if j else -0.12), lab, ha="left" if i else "right", va="bottom" if j else "top", fontsize=12)
for p, q in [((0, 0), (1, 0)), ((0, 0), (0, 1)), ((1, 0), (1, 1)), ((0, 1), (1, 1))]: ax.plot(*zip(p, q), color="k", lw=1)
ax.set_xlim(-0.8, 1.8); ax.set_ylim(-0.6, 1.6); ax.set_title("2 inputs: a ∘ₖ b = a + b + k·ab")
ax = axs[1]; ax.axis("off"); ax.set_aspect("equal")
proj = lambda v: (v[0] + 0.45 * v[2], v[1] + 0.35 * v[2])
names = "abc"
import itertools
for v in itertools.product([0, 1], repeat=3):
    P = proj(v); m = sum(v)
    col = [GRAY, BLUE, ORANGE, RED][m]
    lab = "0" if m == 0 else ("" if m == 1 else ("k·" if m == 2 else "k²·")) + "".join(n for n, s in zip(names, v) if s)
    ax.plot(*P, "o", ms=12, color=col); ax.text(P[0] + 0.07, P[1] + 0.07, lab, fontsize=11)
    for i in range(3):
        if v[i] == 0:
            w = list(v); w[i] = 1; ax.plot(*zip(P, proj(w)), color="k", lw=0.8)
ax.set_xlim(-0.4, 1.9); ax.set_ylim(-0.4, 1.7)
ax.set_title("3 inputs: sum of corners, rank-m corner weighted k^(m−1)")
fig.suptitle("n-fold combine ∘ₖ(a₁, …, aₙ) = (Π(1 + kaᵢ) − 1)/k = Σ k^(m−1)·e_m.  Matrix form: det(I + kA) = 1 + k·(e₁ + k e₂ + k² e₃ + …)", fontsize=11)
save(fig, "fig05-cube-combine")

# 06 Mobius inversion
fig, ax = plt.subplots(figsize=(7, 5)); ax.axis("off"); ax.set_aspect("equal")
for v in itertools.product([0, 1], repeat=3):
    P = proj(v); m = sum(v); sign = (-1) ** (3 - m)
    col = GREEN if sign > 0 else RED
    lab = ("+" if sign > 0 else "−") + " f(" + ("{" + ",".join(n for n, s in zip(names, v) if s) + "}" if m else "∅") + ")"
    ax.plot(*P, "o", ms=13, color=col); ax.text(P[0] + 0.07, P[1] + 0.07, lab, fontsize=11, color=col)
    for i in range(3):
        if v[i] == 0:
            w = list(v); w[i] = 1; ax.plot(*zip(P, proj(w)), color="k", lw=0.8)
ax.set_xlim(-0.4, 2.2); ax.set_ylim(-0.4, 1.7)
ax.set_title("Möbius inversion (Rota 1964): alternate signs down the cube\n= the top corner alone, k²·abc.  It undoes the n-fold combine")
save(fig, "fig06-mobius-inversion")

# 07 flip
fig, axs = plt.subplots(1, 2, figsize=(11, 4.5))
ax = axs[0]; x = np.linspace(0, 6, 400)
for k, col in [(1, BLUE), (0.5, GREEN), (2, ORANGE)]:
    ax.plot(x, x / (1 + k * x), color=col, lw=2, label=f"k = {k}:  x/(1 + kx), levels off at 1/k")
ax.set_xlabel("weight x (k = 1: the odds)"); ax.set_ylabel("flipped value μ (k = 1: the probability)"); ax.legend(fontsize=9)
ax.set_title("the flip squeezes weights into probabilities")
ax = axs[1]; ax.axis("off"); ax.set_xlim(0, 10); ax.set_ylim(0, 6)
box(ax, (2, 4.5), "weights x, y\ncombine by ∘ₖ", BLUE); box(ax, (8, 4.5), "probabilities μ, ν\ncombine by ∘₋ₖ", ORANGE)
arrow(ax, (3.6, 4.7), (6.3, 4.7), "flip x ↦ x/(1 + kx)", GRAY, fs=10, off=(0, 0.5))
arrow(ax, (6.3, 4.1), (3.6, 4.1), "flip back x ↦ x/(1 − kx)", GRAY, fs=10, off=(0, -0.5))
ax.text(5, 1.7, "k = 1: odds 1 and 1 join to 3;\nprobabilities 1/2 and 1/2 join to 3/4 by a + b − ab\n\nthe flip is minus the way back: x/(1+kx) = −a′", ha="center", fontsize=10)
fig.suptitle("Flip φₖ(x) = x/(1 + kx): carries ∘ₖ to ∘₋ₖ.  Known: odds ↔ probability (k = 1); Tsallis additive duality q ↔ 2 − q", fontsize=11)
save(fig, "fig07-flip")

# 08 mirror
fig, ax = plt.subplots(figsize=(7.5, 5))
x = np.linspace(-3, 3, 400); b = 0.6
for k, col in [(1, BLUE), (-1, ORANGE)]:
    y = x + b + k * x * b; ax.plot(x, y, color=col, lw=2, label=f"a ∘_{k:+d} b with b = {b if k > 0 else -b}" if False else f"k = {k:+d}")
    ax.axvline(-1 / k, color=col, ls=":", lw=1.5); ax.text(-1 / k, 3.6, f"wall {-1/k:+.0f}", color=col, ha="center")
ax.plot(x, x + 0.6 + 0.6 * x, color=BLUE, lw=2)
ax.plot(x, -((-x) + 0.6 + 0.6 * (-x)), color=ORANGE, lw=2, ls="--", label="−((−a) ∘₊₁ b): the same line as k = −1, b = −0.6")
ax.set_ylim(-4, 4); ax.axhline(0, color="k", lw=0.5); ax.axvline(0, color="k", lw=0.5); ax.legend(fontsize=9, loc="lower right")
ax.set_xlabel("a"); ax.set_title("Mirror x ↦ −x:  a ∘₋ₖ b = −((−a) ∘ₖ (−b)).  Swaps k ↔ −k and the walls ±1/k\nKnown: the two ring-theory sign conventions (Divinsky + ; Coleman–Easdown −)")
save(fig, "fig08-mirror")

# 09 logarithm
fig, axs = plt.subplots(1, 2, figsize=(12, 4.8))
ax = axs[0]; x = np.linspace(-0.6, 3, 400)
for k, col in [(0.001, GRAY), (0.5, GREEN), (1, BLUE), (2, ORANGE)]:
    xx = x[1 + k * x > 0]; ax.plot(xx, np.log1p(k * xx) / k, color=col, lw=2, label=("k → 0: Lₖ(x) = x" if k < 0.01 else f"k = {k}"))
ax.set_xlabel("x"); ax.set_ylabel("Lₖ(x) = log(1 + kx)/k"); ax.legend(fontsize=9); ax.set_title("the ruler: equal steps in Lₖ are equal ∘ₖ-steps")
for i in range(1, 6):  # ruler ticks for k=1
    v = math.expm1(i * 0.25); ax.plot([v, v], [-0.9, -0.75], color=BLUE)
ax.text(0.9, -0.62, "ticks: 0.25, 0.5, … on the k = 1 ruler", color=BLUE, fontsize=9)
ax = axs[1]; x = np.linspace(-0.9, 0.9, 400)
ax.plot(x, np.arctanh(x), color=BLUE, lw=2, label="(k, d) = (0, −1): artanh x  (two walls ±1)")
ax.plot(x, x / (1 + x / 2) if False else x / (1 + 0 * x), color=GRAY, lw=1, ls=":")
xx = x[x > -0.99]; ax.plot(xx, xx / (1 + xx), color=GREEN, lw=2, label="(k, d) = (2, 1): x/(1 + x)  (one double wall)")
ax.plot(x, np.arctan(x), color=ORANGE, lw=2, label="(k, d) = (0, 1): arctan x  (no real wall)")
ax.set_ylim(-2, 2); ax.legend(fontsize=9); ax.set_xlabel("x"); ax.set_title("the two-dial law's ruler L_{k,d} = ∫₀ˣ dt/(1 + kt + dt²)")
fig.suptitle("Logarithm Lₖ (project words: bend, ruler, flat coordinate): Lₖ(a ∘ₖ b) = Lₖ(a) + Lₖ(b).  Known: the logarithm of the formal group law; Aczél's generator", fontsize=11)
save(fig, "fig09-logarithm")

# 10 k-log k-exp
fig, axs = plt.subplots(1, 2, figsize=(11, 4.5))
ax = axs[0]; x = np.linspace(0.05, 4, 300)
for k, col in [(-1, PURPLE), (-0.5, GREEN), (1e-6, GRAY), (0.5, BLUE), (1, ORANGE)]:
    ax.plot(x, (x ** k - 1) / k, color=col, lw=2, label="k → 0: log x" if abs(k) < 1e-3 else f"k = {k}")
ax.set_ylim(-4, 4); ax.legend(fontsize=9); ax.set_xlabel("x"); ax.set_title("k-log (Box–Cox): (xᵏ − 1)/k\ncarries multiplication to ∘ₖ")
ax = axs[1]; x = np.linspace(-2, 2, 300)
for k, col in [(-0.5, GREEN), (1e-6, GRAY), (0.5, BLUE), (1, ORANGE)]:
    base = 1 + k * x; yy = np.where(base > 0, np.abs(base) ** (1 / k), np.nan); ax.plot(x, yy, color=col, lw=2, label="k → 0: eˣ" if abs(k) < 1e-3 else f"k = {k}")
ax.set_ylim(0, 6); ax.legend(fontsize=9); ax.set_xlabel("x"); ax.set_title("k-exp (Tsallis q-exponential, q = 1 − k): (1 + kx)^(1/k)\nand Lₖ = log ∘ k-exp")
fig.suptitle("k-log and k-exp: the bridges between × and ∘ₖ.  Known: Box & Cox 1964; Tsallis 1988", fontsize=11)
save(fig, "fig10-k-log-k-exp")

# 11 landscape moves
fig, ax = plt.subplots(figsize=(9, 7))
K = np.linspace(-3, 3, 400)
ax.fill_between(K, K ** 2 / 4, 4, color=ORANGE, alpha=0.12); ax.fill_between(K, -3, K ** 2 / 4, color=BLUE, alpha=0.10)
ax.plot(K, K ** 2 / 4, color="k", lw=2)
ax.text(-2.9, 3.5, "k² − 4d < 0: no real wall\n(class of tangent addition (a+b)/(1−ab))", color=ORANGE, fontsize=10)
ax.text(0.6, -2.7, "k² − 4d > 0: two walls\n(class of ∘₁, velocity addition (a+b)/(1+ab))", color=BLUE, fontsize=10)
ax.annotate("k² − 4d = 0: one double wall\n(class of plain +)", (2.3, 2.3 ** 2 / 4), (2.05, 0.3), fontsize=10, arrowprops=dict(arrowstyle="->"), bbox=dict(fc="white", ec="none", alpha=0.85))
# slide orbits: parabolas d = k^2/4 - D/4
for D, col in [(4, BLUE), (-4, ORANGE)]:
    ax.plot(K, K ** 2 / 4 - D / 4, color=col, lw=1.2, ls="--")
k0, d0 = -1.6, K[0] * 0
ax.annotate("", (1.0, 1.0 ** 2 / 4 + 1), (-1.0, 1.0 ** 2 / 4 + 1), arrowprops=dict(arrowstyle="-|>", color=GREEN, lw=2, connectionstyle="arc3,rad=0.25"))
ax.text(-0.9, 1.75, "slide y/(1+cy):  (k, d) → (k+2c, d+kc+c²)\nmoves along a dashed parabola (same discriminant)", color=GREEN, fontsize=9)
# rescale rays d = alpha k^2 through origin
for al in [-0.6, 0.6]:
    kk = np.linspace(0.2, 2.6, 50); ax.plot(kk, al * kk ** 2, color=PURPLE, lw=1, ls=":")
ax.annotate("", (1.9, -0.6 * 1.9 ** 2), (0.9, -0.6 * 0.9 ** 2), arrowprops=dict(arrowstyle="-|>", color=PURPLE, lw=2.5))
ax.text(2.35, -2.2, "rescale x ↦ λx:\n(k, d) → (λk, λ²d)", color=PURPLE, fontsize=9)
ax.annotate("", (-1.0, 0.0), (1.0, 0.0), arrowprops=dict(arrowstyle="<|-|>", color=RED, lw=1.5))
ax.text(-0.5, -0.35, "mirror k ↔ −k", color=RED, fontsize=9)
for (kk, dd, lab) in [(1, 0, "∘₁"), (-1, 0, "∘₋₁"), (0, 0, "+"), (0, 1, "tan add"), (0, -1, "velocity add")]:
    ax.plot(kk, dd, "ko", ms=7); ax.text(kk + 0.07, dd + 0.12, lab, fontsize=10)
ax.set_xlim(-3, 3.6); ax.set_ylim(-3, 4); ax.set_xlabel("k (corner weight)"); ax.set_ylabel("d (second dial)")
ax.set_title("Dial moves on the (k, d) plane of (a + b + k·ab)/(1 − d·ab).  Slide + rescale leave exactly 3 classes\nKnown: Coleman–McGuinness 1991 (rational formal group laws); the sign of a quadratic form's discriminant")
save(fig, "fig11-landscape-moves")

# 12 dequantization / contraction
fig, axs = plt.subplots(1, 2, figsize=(11, 4.4))
ax = axs[0]; a = np.linspace(-2, 2, 400)
for h, col in [(1, ORANGE), (0.4, BLUE), (0.1, GREEN), (0.02, GRAY)]:
    ax.plot(a, h * np.logaddexp(a / h, 0), color=col, lw=2, label=f"h = {h}")
ax.set_xlabel("a  (b = 0)"); ax.legend(fontsize=9); ax.set_title("rescale to zero: h·log(e^(a/h) + e^(b/h)) → max(a, b)\n(Maslov dequantization; the 'tropical limit')")
ax = axs[1]; a = np.linspace(-0.9, 2, 400); b = 0.5
for k, col in [(2, ORANGE), (1, BLUE), (0.3, GREEN), (0, GRAY)]:
    ax.plot(a, a + b + k * a * b, color=col, lw=2, label=f"k = {k}")
ax.set_xlabel("a  (b = 0.5)"); ax.set_ylabel("a ∘ₖ b"); ax.legend(fontsize=9)
ax.set_title("turning k → 0: ∘ₖ contracts to +\n(a deformation; Inönü–Wigner contraction in Lie terms)")
fig.suptitle("Turning a dial to its end: two known limits", fontsize=11)
save(fig, "fig12-limits")

# 13 height flow
fig, ax = plt.subplots(figsize=(7.5, 5))
xs = np.linspace(0.001, 2.5, 160)
for t, col in [(-1, PURPLE), (-0.5, GREEN), (0, GRAY), (0.5, BLUE), (1, ORANGE), (1.5, RED)]:
    ys = [Et(t, v) for v in xs]; ax.plot(xs, ys, color=col, lw=2, label=f"t = {t}" + ("  (log(1+x))" if t == -1 else "  (x)" if t == 0 else "  (eˣ − 1)" if t == 1 else ""))
ax.set_ylim(0, 6); ax.set_xlabel("x"); ax.set_ylabel("E^t(x)"); ax.legend(fontsize=9)
ax.set_title("Height flow E^t: 'do eˣ − 1 t times', for any real t\nKnown: regular iteration at a parabolic fixed point (Szekeres 1958, Écalle); Abel function = Fatou coordinate")
save(fig, "fig13-height-flow")

# 14 B_s
fig, ax = plt.subplots(figsize=(7.5, 5))
xs = np.linspace(0.0, 1.6, 120)
for s, col in [(-1, PURPLE), (0, GRAY), (0.5, GREEN), (1, BLUE), (1.5, ORANGE), (2, RED)]:
    ax.plot(xs, [Bs(s, v, 0.5) for v in xs], color=col, lw=2, label=f"s = {s}" + ("  (a + b)" if s == 0 else "  (a ∘₁ b)" if s == 1 else ""))
ax.set_ylim(0, 8); ax.set_xlabel("a  (b = 0.5)"); ax.set_ylabel("B_s(a, 0.5)"); ax.legend(fontsize=9)
ax.set_title("Height-s operation B_s(a, b) = E^s(E^(−s)a + E^(−s)b)\nthe transport of + by E^s; its corner weight at 0 is s (same as k)")
save(fig, "fig14-height-s-operation")

# 15 climb
fig, ax = plt.subplots(figsize=(11, 5)); ax.axis("off"); ax.set_xlim(0, 12); ax.set_ylim(0, 6)
ours = ["a + b", "a + b + ab", "(1+a)(1+b)e^(L(a)L(b)) − 1", "…"]
benn = ["a + b", "a · b", "a^(log b)", "…"]
for i, (o, bb) in enumerate(zip(ours, benn)):
    box(ax, (1.5 + 3.3 * i, 4.4), o, BLUE, 10); box(ax, (1.5 + 3.3 * i, 1.4), bb, ORANGE, 10)
    if i < 3:
        arrow(ax, (2.6 + 3.3 * i, 4.4), (3.9 + 3.3 * i, 4.4), "E₊", BLUE, off=(0, 0.45))
        arrow(ax, (2.6 + 3.3 * i, 1.4), (3.9 + 3.3 * i, 1.4), "exp₊", ORANGE, off=(0, 0.45))
ax.text(6, 5.6, "ours: climb one rung = transport by E(x) = eˣ − 1.  Identity stays 0 at every height", ha="center", color=BLUE)
ax.text(6, 2.6, "Bennett 1915, commutative hyperoperations: transport by exp.  Identity moves (0, then 1, then e)", ha="center", color=ORANGE)
ax.text(6, 0.3, "E = (x ↦ x − 1) ∘ exp, so our climb is Bennett's climb followed by a fuse (shift back by 1)", ha="center", fontsize=10)
ax.set_title("Climb a rung: B_(s+1) = E₊(B_s).  Half a climb is transport by E^(1/2)")
save(fig, "fig15-climb")

# 16 fuse
fig, axs = plt.subplots(1, 2, figsize=(11, 4.2))
for ax, title, f, idv in [(axs[0], "plain rung: a × b, identity 1, absorber 0", lambda a, b: a * b, 1),
                          (axs[1], "fused rung: (1+a)(1+b) − 1 = a ∘₁ b, identity 0, absorber −1", lambda a, b: (1 + a) * (1 + b) - 1, 0)]:
    g = np.arange(-1, 4) if idv == 1 else np.arange(-2, 3)
    M = np.array([[f(a, b) for a in g] for b in g])
    ax.imshow(M, cmap="Blues", origin="lower", extent=[g[0] - 0.5, g[-1] + 0.5, g[0] - 0.5, g[-1] + 0.5])
    for i, bb in enumerate(g):
        for j, a in enumerate(g): ax.text(a, bb, str(M[i, j]), ha="center", va="center", fontsize=9, color="k")
    ax.axvline(idv, color=GREEN, lw=2); ax.axhline(idv, color=GREEN, lw=2); ax.set_title(title, fontsize=10); ax.set_xlabel("a"); ax.set_ylabel("b")
fig.suptitle("Fuse = transport by x ↦ x − 1 (i.e. read every rung through the lift x ↦ 1 + x).  The table just slides one step down-left", fontsize=11)
save(fig, "fig16-fuse")

# 17 translation fan
fig, ax = plt.subplots(figsize=(7, 5.5))
k = 1.0; x = np.linspace(-3, 2, 100)
for a, col in zip([-0.75, -0.5, 0, 0.5, 1, 2], [PURPLE, GREEN, GRAY, BLUE, ORANGE, RED]):
    ax.plot(x, a + (1 + k * a) * x, color=col, lw=2, label=f"T_a, a = {a}")
ax.plot(-1 / k, -1 / k, "ko", ms=9); ax.annotate("every T_a passes through\nthe wall (−1/k, −1/k)", (-1, -1), (-2.9, 1.6), arrowprops=dict(arrowstyle="->"))
ax.set_ylim(-4, 4); ax.axhline(0, color="k", lw=0.5); ax.axvline(0, color="k", lw=0.5); ax.legend(fontsize=9, loc="lower right"); ax.set_xlabel("x")
ax.set_title("Translation T_a(x) = a ∘ₖ x = (1 + ka)x + a, and T_a ∘ T_b = T_(a∘ₖb)\nKnown: the regular representation; a stretch about the wall")
save(fig, "fig17-translation")

# 18 tilt
fig, ax = plt.subplots(figsize=(7.5, 4.6))
n = np.arange(0, 12); mu = 2.0
p = np.array([math.exp(-mu) * mu ** int(i) / math.factorial(int(i)) for i in n])
for t, col, off in [(0.5, GREEN, -0.27), (1, GRAY, 0), (2, ORANGE, 0.27)]:
    q = p * t ** n; q = q / q.sum(); ax.bar(n + off, q, width=0.27, color=col, label=f"t = {t}" + ("  (original)" if t == 1 else ""))
ax.set_xlabel("count N"); ax.set_ylabel("probability"); ax.legend()
ax.set_title("Exponential tilting: reweight each outcome by t^N and renormalise\nKnown: Esscher 1932.  Used on the point-process and statistics pages")
save(fig, "fig18-tilt")

# 19 zoom
fig, ax = plt.subplots(figsize=(7.5, 4.8))
s = np.linspace(0, 2.5, 300)
for (lab, rs, col) in [("Cantor set: 2 pieces at 1/3", [1 / 3] * 2, BLUE), ("half box: 2 pieces at 1/4", [1 / 4] * 2, GREEN), ("gasket: 3 pieces at 1/2", [1 / 2] * 3, ORANGE)]:
    xs = sum(np.array(r) ** s for r in rs); ax.plot(s, xs, color=col, lw=2, label=lab)
    D = math.log(len(rs)) / -math.log(rs[0]); ax.plot(D, 1, "o", color=col); ax.text(D, 1.12, f"D = {D:.3f}", color=col, ha="center", fontsize=9)
ax.axhline(1, color=RED, ls="--"); ax.text(2.0, 1.05, "wall: x(s) = 1", color=RED)
ax.set_ylim(0, 3.2); ax.set_xlabel("zoom exponent s"); ax.set_ylabel("x(s) = Σ rᵢˢ"); ax.legend(fontsize=9)
ax.set_title("Zoom dial: a piece shrunk by r weighs rˢ.  Dimension D is where x(D) = 1\nKnown: Moran 1946; complex dimensions (Lapidus–van Frankenhuijsen)")
save(fig, "fig19-zoom")

# 20 order gap
fig, ax = plt.subplots(figsize=(6.5, 5.5)); ax.set_aspect("equal")
u, v = np.array([0.3, 0.5]), np.array([-0.2, 0.9])
ax.add_patch(Polygon([[0, 0], u, u + v, v], closed=True, fc=ORANGE, alpha=0.3, ec=ORANGE))
for w, lab, col in [(u, "(s₁, c₁) = (0.3, 0.5)", BLUE), (v, "(s₂, c₂) = (−0.2, 0.9)", GREEN)]:
    ax.annotate("", w, (0, 0), arrowprops=dict(arrowstyle="-|>", color=col, lw=2)); ax.text(*(w * 1.05), lab, color=col)
ax.text(0.02, 0.62, f"area = s₁c₂ − s₂c₁ = {u[0]*v[1]-v[0]*u[1]:.2f}", color=ORANGE, fontsize=11)
ax.set_xlim(-0.5, 0.7); ax.set_ylim(-0.1, 1.6); ax.set_xlabel("stretch s"); ax.set_ylabel("shift c")
ax.set_title("Order gap: x ↦ (1+s)x + c done in the two orders\ndiffers by a shift equal to this oriented area (a wedge)\nKnown: the commutator in the affine group")
save(fig, "fig20-order-gap")

# 21 usage
cnt = json.load(open("counts.json", encoding="utf-8"))
items = sorted(cnt.items(), key=lambda kv: kv[1]["uses"])
fig, ax = plt.subplots(figsize=(9, 7.5))
ax.barh([k for k, _ in items], [v["uses"] for _, v in items], color=BLUE)
for i, (_, v) in enumerate(items): ax.text(v["uses"] + 15, i, f'{v["uses"]} in {v["files"]} files', va="center", fontsize=8)
ax.set_xlabel("matches of the move's words and formulas across math-theory (rough: words like 'wall' and 'lift' include some other senses)")
ax.set_title("How often each move shows up (946 files: pages, checks, Lean)")
save(fig, "fig21-usage")
