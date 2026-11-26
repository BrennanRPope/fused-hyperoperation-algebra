"""Pictures for the Yang–Mills gap DRAFT. Run after checks.py (reads results.json).
python3 -I figures.py -> figures/*.png"""
import json, math
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch, Circle, Rectangle
from scipy import special
import mpmath as mp

R = json.load(open("results.json"))
BLUE, ORANGE, GREY, RED, GREEN = "#2b6cb0", "#dd6b20", "#718096", "#c53030", "#2f855a"
plt.rcParams.update({"font.size": 10, "axes.spines.top": False, "axes.spines.right": False})


def save(fig, name):
    fig.tight_layout()
    fig.savefig(f"figures/{name}.png", dpi=140)
    plt.close(fig)

# f1: which laws close up -------------------------------------------------------------
fig, ax = plt.subplots(1, 3, figsize=(13, 4))
x = np.linspace(-3, 3, 600)
for (k, d, lab, c) in [(3, 1, "line class  Δ = 5 > 0  (two walls)", ORANGE), (2, 1, "border  Δ = 0  (double wall)", GREY),
                       (1, 1, "circle class  Δ = −3 < 0  (no wall)", BLUE)]:
    ax[0].plot(x, 1 + k * x + d * x * x, color=c, lw=2, label=lab)
ax[0].axhline(0, color="k", lw=0.6)
ax[0].set_ylim(-2, 6); ax[0].set_xlabel("x"); ax[0].set_ylabel("generator v(x) = 1 + kx + dx²")
ax[0].set_title("A. Walls are where v(x) = 0"); ax[0].legend(fontsize=8, loc="upper left")
# B: the ruler closes up into a circle only for Delta<0
t = np.linspace(0, 2 * np.pi, 300)
cx = 0
for D in [-8, -3, -1]:
    P = 2 * np.pi / np.sqrt(-D)
    r = P / (2 * np.pi)
    ax[1].plot(cx + r * np.cos(t), r * np.sin(t), color=BLUE, lw=1.8)
    ax[1].text(cx, -r - 0.25, f"Δ = {D}\nlength {P:.2f}", fontsize=8, color=BLUE, ha="center", va="top")
    cx += r + 1.3
ax[1].annotate("", xy=(4.6, 1.6), xytext=(-0.6, 1.6), arrowprops=dict(arrowstyle="<->", color=ORANGE, lw=2))
ax[1].text(2.0, 1.8, "line class / border: ruler runs off to infinity at a wall", fontsize=8, color=ORANGE, ha="center")
ax[1].set_xlim(-0.8, 4.8); ax[1].set_ylim(-2.2, 2.2)
ax[1].set_aspect("equal"); ax[1].axis("off")
ax[1].set_title("B. In ruler units the law is +; only Δ < 0 closes\ninto a circle, of length 2π/√(−Δ)")
# C: spectra
for i, D in enumerate([-3, -8]):
    for n in range(0, 4 if D == -3 else 3):
        ax[2].plot([i - 0.3, i + 0.3], [n * n * (-D)] * 2, color=BLUE, lw=2)
        ax[2].text(i + 0.33, n * n * (-D), f"rung {n}", fontsize=7, va="center", color=BLUE)
    ax[2].annotate("", xy=(i, -D), xytext=(i, 0), arrowprops=dict(arrowstyle="<->", color=RED))
    ax[2].text(i - 0.05, -D / 2, f"gap {-D}", color=RED, fontsize=8, ha="right")
ax[2].fill_between([1.7, 2.3], 0, 48, color=ORANGE, alpha=0.35)
ax[2].text(2, 49, "line class:\nunbroken band\nfrom 0, no gap", ha="center", fontsize=8, color=ORANGE)
ax[2].set_xticks([0, 1, 2]); ax[2].set_xticklabels(["circle Δ=−3", "circle Δ=−8", "line / border"])
ax[2].set_ylim(-2, 60); ax[2].set_ylabel("Laplacian level (in the dial's own units)")
ax[2].set_title("C. Gap = −Δ = 4 × curvature reading")
save(fig, "f1_which_laws_close_up")

# f2: the plaquette ----------------------------------------------------------------------
fig, ax = plt.subplots(1, 2, figsize=(11, 4.2))
a = ax[0]
pts = [(0, 0), (1, 0), (1, 1), (0, 1), (0, 0)]
labs = ["U₁ (go)", "U₂ (go)", "U₃⁻¹ (back)", "U₄⁻¹ (back)"]
for i in range(4):
    (x0, y0), (x1, y1) = pts[i], pts[i + 1]
    a.add_patch(FancyArrowPatch((x0, y0), (x1, y1), arrowstyle="-|>", mutation_scale=18, color=BLUE, lw=2.2, shrinkA=6, shrinkB=6))
    a.text((x0 + x1) / 2 + (0.12 if x0 == x1 else 0) * (1 if x0 > 0.5 else -1) - (0.12 if x0 == x1 and x0 < 0.5 else 0),
           (y0 + y1) / 2 + (0.08 if y0 == y1 and y0 > 0.5 else -0.1 if y0 == y1 else 0), labs[i], ha="center", fontsize=9)
a.text(0.5, 0.5, "holonomy\n= U₁ ∘ U₂ ∘ U₃⁻¹ ∘ U₄⁻¹\n(lift: I + kA per link)", ha="center", va="center", fontsize=9)
a.text(0.5, -0.38, "shared centre (abelian): order free, only the total ruler counts\n"
       "different centres (non-abelian): an extra order gap = [A_x, A_y]", ha="center", fontsize=8.5, color=GREY)
a.set_xlim(-0.5, 1.5); a.set_ylim(-0.55, 1.3); a.axis("off"); a.set_title("A. A plaquette is four links composed round a square")
F = np.array(R["B_errors"])
b = ax[1]
b.loglog([0.1, 0.05, 0.025], F, "o-", color=BLUE, label="|log(plaquette)/a² − F|")
b.loglog([0.1, 0.025], [F[0], F[0] / 4], "--", color=GREY, label="slope 1 (error ∝ a)")
b.set_xlabel("lattice spacing a"); b.set_ylabel("error"); b.legend(fontsize=8)
b.set_title("B. Shrinking the square gives the field strength\nF = ∂ₓA_y − ∂_yAₓ + [Aₓ, A_y]  (here |[A,A]| ≈ ½|F|)")
save(fig, "f2_plaquette")

# f3: one plaquette -----------------------------------------------------------------------
NMAX = 60


def u1(g):
    n = np.arange(-NMAX, NMAX + 1)
    H = np.diag(2 * g * g * n ** 2 + 1 / g ** 2) - (np.eye(len(n), k=1) + np.eye(len(n), k=-1)) / (2 * g * g)
    return np.linalg.eigvalsh(H)


def su2(g):
    js = np.arange(0, 121) / 2
    H = np.diag(2 * g * g * js * (js + 1) + 1 / g ** 2)
    for i in range(len(js) - 1):
        H[i, i + 1] = H[i + 1, i] = -1 / (2 * g * g)
    return np.linalg.eigvalsh(H)

gs = np.geomspace(0.15, 4, 120)
fig, ax = plt.subplots(1, 3, figsize=(14, 4))
lev = np.array([u1(g)[:7] for g in gs])
for i in range(7):
    ax[0].plot(gs, (lev[:, i] - lev[:, 0]) / 2, color=BLUE if i else "k", lw=1.4)
ax[0].set_xscale("log"); ax[0].set_ylim(0, 12)
ax[0].set_xlabel("coupling g"); ax[0].set_ylabel("level − ground   (units: 2)")
ax[0].text(0.17, 5.4, "weak coupling:\nevenly spaced\n(oscillator)", fontsize=8)
ax[0].text(1.6, 1.4, "strong coupling:\nrungs 0, ±1, ±2 ...\nspacing ∝ n²", fontsize=8)
ax[0].set_title("A. U(1) plaquette: levels from oscillator to rungs")
from matplotlib.ticker import NullFormatter
for a_ in ax[:2]:
    a_.set_xticks([0.2, 0.5, 1, 2, 4]); a_.set_xticklabels(["0.2", "0.5", "1", "2", "4"]); a_.xaxis.set_minor_formatter(NullFormatter())
gapU = np.array(R["C"]["gap_u1"]); gapS = np.array(R["C"]["gap_su2"]); gg = np.array(R["C"]["g"])
ax[1].plot(gg, gapU, color=BLUE, lw=2, label="U(1)")
ax[1].plot(gg, gapS, color=GREEN, lw=2, label="SU(2)")
ax[1].plot(gg, 2 * gg ** 2, "--", color=GREY, label="2g² (rung 1, U(1))")
ax[1].plot(gg, 1.5 * gg ** 2, ":", color=GREY, label="1.5g² (spin ½, SU(2))")
ax[1].axhline(2, color=ORANGE, lw=1, ls="--"); ax[1].text(0.16, 2.15, "oscillator 2", color=ORANGE, fontsize=8)
ax[1].set_xscale("log"); ax[1].set_yscale("log"); ax[1].set_ylim(1, 50)
ax[1].set_xlabel("coupling g"); ax[1].set_ylabel("gap"); ax[1].legend(fontsize=8)
ax[1].set_title("B. The gap dips in the crossover\nbut never closes (min 1.76 and 1.55)")
# C rung picture: wavefunction weights on rungs at strong vs weak
for g, c, lab in [(3.0, BLUE, "strong g = 3"), (0.3, ORANGE, "weak g = 0.3")]:
    n = np.arange(-NMAX, NMAX + 1)
    H = np.diag(2 * g * g * n ** 2 + 1 / g ** 2) - (np.eye(len(n), k=1) + np.eye(len(n), k=-1)) / (2 * g * g)
    w, v = np.linalg.eigh(H)
    sel = np.abs(n) <= 8
    ax[2].bar(n[sel] + (0.2 if c == ORANGE else -0.2), v[sel, 1] ** 2, width=0.4, color=c, label=lab)
ax[2].set_xlabel("rung n (electric flux)"); ax[2].set_ylabel("weight in first excited state")
ax[2].legend(fontsize=8); ax[2].set_title("C. Strong coupling: one rung (±1). Weak: many rungs")
save(fig, "f3_one_plaquette")

# f4: two dimensions -------------------------------------------------------------------------
fig, ax = plt.subplots(1, 2, figsize=(12, 4.6))
b = np.array(R["D_ratio2"]["beta"]); r = np.array(R["D_ratio2"]["ratio"])
ax[0].plot(b, r, color=BLUE, lw=2, label="σ₂/σ₁ (exact)")
bs = b[b < 0.5]
ax[0].plot(bs, 2 + np.log(2) / np.log(2 / bs), "--", color=GREY, label="2 + ln2/ln(2/β) (strong)")
ax[0].axhline(4, color=ORANGE, ls=":"); ax[0].text(80, 3.85, "4 = 2² (Casimir)", color=ORANGE, fontsize=8, ha="right")
ax[0].axhline(2, color=ORANGE, ls=":"); ax[0].text(0.06, 1.85, "2 (counting)", color=ORANGE, fontsize=8)
ax[0].set_xscale("log"); ax[0].set_ylim(1.6, 4.3); ax[0].set_xlabel("β = 1/g²"); ax[0].set_ylabel("string tension, rung 2 / rung 1")
ax[0].legend(fontsize=8, loc="center right"); ax[0].set_title("A. The rung ladder's spacing: linear at strong coupling,\nsquares (n²) at weak coupling")
# zeros
X, Y = np.meshgrid(np.linspace(0.02, 6, 300), np.linspace(-15, 15, 600))
rr = np.abs(special.iv(1, X + 1j * Y) / special.iv(0, X + 1j * Y))
ax[1].contourf(X, Y, rr >= 1, levels=[0.5, 1.5], colors=["#feebc8"])
ax[1].contour(X, Y, rr, levels=[1], colors=[ORANGE], linewidths=1)
for N, mk, c in [("8", "o", "#90cdf4"), ("16", "s", "#4299e1"), ("32", "^", "#2a4365")]:
    z = np.array(R["D_zeros"][N])
    ax[1].scatter(z[:, 0], z[:, 1], marker=mk, color=c, s=22, label=f"zeros, N = {N} plaquettes", zorder=3)
    ax[1].scatter(z[:, 0], -z[:, 1], marker=mk, color=c, s=22, zorder=3)
for j in R["D_j0"][:3]:
    ax[1].plot(0, j, "x", color=RED); ax[1].plot(0, -j, "x", color=RED)
ax[1].plot([0, 6], [0, 0], color=GREEN, lw=4, alpha=0.6)
ax[1].text(3.2, 0.5, "real β: rung 0 always heaviest, no zeros, gap open", color=GREEN, fontsize=8)
ax[1].text(2.6, 11, "shaded: rung 1 outweighs rung 0\n× : rung 0's weight is zero", fontsize=8, color="#975a16")
ax[1].set_xlabel("Re β"); ax[1].set_ylabel("Im β"); ax[1].legend(fontsize=7, loc="lower right")
ax[1].set_title("B. 2D U(1) torus: zeros gather where the top rungs tie,\nnear the imaginary axis, never on real β")
save(fig, "f4_two_dimensions")

# f5: twist ------------------------------------------------------------------------------------
fig, ax = plt.subplots(1, 2, figsize=(12, 4.2))
th = np.linspace(0, 2 * np.pi, 400)
for n in range(-2, 2):
    ax[0].plot(th, (n + th / (2 * np.pi)) ** 2, color=BLUE, lw=1.6)
for n, xpos in [(0, 0.3), (-1, 0.3), (1, 0.3)]:
    ax[0].text(xpos, (n + xpos / (2 * np.pi)) ** 2 + 0.07, f"rung {n}", fontsize=8, color=BLUE)
ax[0].axvline(np.pi, color=RED, ls="--"); ax[0].text(np.pi + 0.05, 2.0, "half a rung:\ntwo lowest tie,\ngap closes", color=RED, fontsize=8)
ax[0].annotate("", xy=(np.pi, 2.25), xytext=(np.pi, 0.25), arrowprops=dict(arrowstyle="<->", color=GREEN))
ax[0].text(np.pi - 0.1, 1.2, "each vacuum\nstill gapped", color=GREEN, fontsize=8, ha="right")
ax[0].set_ylim(0, 2.6); ax[0].set_xlabel("twist θ (holonomy round the circle)"); ax[0].set_ylabel("energy / (g²L/2)")
ax[0].set_xticks([0, np.pi, 2 * np.pi]); ax[0].set_xticklabels(["0", "π", "2π"])
ax[0].set_title("A. A twist slides the ladder by a fractional rung θ/2π")
p = np.array(R["E_pinch"]); Ns = [10, 20, 40, 80]
ax[1].loglog(Ns, p[:, 1], "o-", color=BLUE, label="nearest zero, Im θ (Re θ = π exactly)")
ax[1].loglog(Ns, p[0, 1] * 10 / np.array(Ns), "--", color=GREY, label="∝ 1/N")
ax[1].set_xlabel("size N (plaquettes)"); ax[1].set_ylabel("distance of zero from the real θ axis")
ax[1].legend(fontsize=8); ax[1].set_title("B. Zeros pinch onto θ = π: a genuine phase change,\nbut a first-order one")
save(fig, "f5_twist")

# f6: dial gas strips --------------------------------------------------------------------------
import importlib.util


def dg(k, W):
    z = k ** -2
    S = np.array([[(s >> i) & 1 for i in range(W)] for s in range(2 ** W)], float)
    D = z ** (S.sum(1) / 2) * k ** ((S * np.roll(S, -1, axis=1)).sum(1) / 2)
    T = D[:, None] * k ** (S @ S.T) * D[None, :]
    ev = np.sort(np.linalg.eigvalsh(T))[::-1]
    return math.log(ev[0] / ev[1]), math.log(ev[0] / ev[2])

kc = 3 + 2 * math.sqrt(2)
ks = np.linspace(2.5, 14, 47)
fig, ax = plt.subplots(1, 2, figsize=(12, 4.2))
cols = ["#bee3f8", "#63b3ed", "#3182ce", "#2c5282"]
for W, c in zip([4, 6, 8, 10], cols):
    g1 = [dg(k, W) for k in ks]
    ax[0].plot(ks, [x[0] for x in g1], color=c, lw=1.6, label=f"width {W}")
    ax[1].plot(ks, [W * x[0] for x in g1], color=c, lw=1.6, label=f"width {W}")
ax[0].axvline(kc, color=RED, ls="--"); ax[0].text(kc + 0.2, 1.5, "k_c = 3 + 2√2", color=RED, fontsize=8)
ax[0].text(3, 0.08, "disordered:\ngap settles", fontsize=8); ax[0].text(10.5, 0.4, "ordered: top two\nlevels merge\n(two states, not a particle)", fontsize=8)
ax[0].set_xlabel("dial k (pair weight)"); ax[0].set_ylabel("gap ln(λ₀/λ₁)"); ax[0].set_ylim(0, 2); ax[0].legend(fontsize=8)
ax[0].set_title("A. Our dial gas on strips: gap against the dial")
ax[1].axvline(kc, color=RED, ls="--"); ax[1].axhline(np.pi / 4, color=ORANGE, ls=":")
ax[1].text(3, np.pi / 4 + 0.1, "π/4", color=ORANGE)
ax[1].set_xlabel("dial k"); ax[1].set_ylabel("width × gap"); ax[1].set_ylim(0, 4); ax[1].legend(fontsize=8)
ax[1].set_title("B. Curves cross at k_c: gap ∝ 1/width there,\nso it closes in the infinite plane (continuous change)")
save(fig, "f6_dial_gas_gap")

# f7: height flow -------------------------------------------------------------------------------
mp.mp.dps = 40
C = [mp.mpf(-1) / 36, mp.mpf(1) / 540, mp.mpf(1) / 7776, mp.mpf(-71) / 435456, mp.mpf(8759) / 163296000, mp.mpf(31) / 20995200]


def abel(x, n=1500):
    y = mp.mpf(x)
    for _ in range(n):
        y = mp.log(1 + y)
    return -2 / y + mp.log(y) / 3 + sum(c * y ** (i + 1) for i, c in enumerate(C)) + n

A1 = abel(1)
x0 = np.geomspace(0.03, 1, 60)
m = [float(mp.e ** (-(A1 - abel(x)))) for x in x0]
fig, ax = plt.subplots(1, 2, figsize=(12, 4.2))
ax[0].semilogy(x0, m, color=BLUE, lw=2.5, label="gap m = e^(−heights to strong coupling)")
for pw, c in [(2, "#cbd5e0"), (6, "#a0aec0"), (12, "#718096")]:
    ax[0].semilogy(x0, x0 ** pw * m[-1], "--", color=c, label=f"x₀^{pw}")
ax[0].set_xlabel("starting coupling x₀"); ax[0].set_ylabel("gap in lattice units"); ax[0].set_ylim(1e-30, 2)
ax[0].legend(fontsize=8); ax[0].set_title("A. The gap from our height flow: ∝ x₀^(1/3) e^(−2/x₀)\nflatter at 0 than every power (invisible to every order)")
# B staircase: x_n = E^n(x0)
x = 0.05; xs = [x]
for _ in range(45):
    x = math.expm1(x); xs.append(x)
    if x > 1.5:
        break
ax[1].step(range(len(xs)), xs, where="post", color=BLUE, lw=2)
ax[1].axhline(1, color=RED, ls="--"); ax[1].text(1, 1.05, "strong coupling", color=RED, fontsize=8)
ax[1].set_xlabel("height s (each step = one rung of e^x − 1)"); ax[1].set_ylabel("coupling x")
n1 = next(i for i, v in enumerate(xs) if v > 1)
ax[1].text(2, 0.6, f"from x₀ = 0.05: {n1} rungs to pass 1; 2/x₀ − (1/3)ln x₀ = {2/0.05 - math.log(0.05)/3:.1f}\nzoom e^(rungs) is the size where the gap lives", fontsize=8)
ax[1].set_title("B. Weak coupling crawls for ~2/x₀ rungs, then shoots up")
save(fig, "f7_height_flow_gap")

# f8: four ways a gap closes --------------------------------------------------------------------
fig, ax = plt.subplots(figsize=(12, 3.8))
ax.axis("off")
boxes = [("1. Continuous change", "gap → 0 at one dial value\n(dial gas at k_c, f6)", GREEN),
         ("2. First-order change", "zeros pinch, two vacua tie,\neach still gapped (twist θ = π, f5)", ORANGE),
         ("3. Massless phase", "gap is 0 over a whole range\n(abelian, weak coupling;\nnot computed here)", GREY),
         ("4. End of the line", "gap > 0 for every coupling,\n→ 0 only as coupling → 0\n(height flow, f7)", BLUE)]
for i, (t, s, c) in enumerate(boxes):
    ax.add_patch(Rectangle((i * 3 + 0.1, 0.2), 2.7, 2.2, color=c, alpha=0.18))
    ax.text(i * 3 + 1.45, 2.05, t, ha="center", fontsize=11, weight="bold", color=c)
    ax.text(i * 3 + 1.45, 1.1, s, ha="center", fontsize=9)
ax.text(10.45, -0.2, "↑ where 4D Yang–Mills is expected to sit:\nthe Millennium problem is to prove the zoomed-in gap stays positive here",
        ha="center", fontsize=9, color=BLUE)
ax.set_xlim(0, 12); ax.set_ylim(-0.6, 2.6)
ax.set_title("Four ways a gap can close (1, 2 and 4 computed here with our own tools)")
save(fig, "f8_four_ways")
print("figures done")
