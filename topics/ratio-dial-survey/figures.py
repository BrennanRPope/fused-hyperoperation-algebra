"""Figures for the similarity-ratio survey. Run from this folder: python3 figures.py -> figures/*.png
Colour key: blue = our reading, orange = already known, grey = endpoints / background."""
import os, json, math, itertools
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle

BLUE, ORANGE, GREY, LGREY = "#1f5fbf", "#e07b00", "#7a7a7a", "#d0d0d0"
plt.rcParams.update({"font.size": 15, "axes.titlesize": 16, "axes.labelsize": 15, "legend.fontsize": 12,
                     "xtick.labelsize": 13, "ytick.labelsize": 13})
here = os.path.dirname(os.path.abspath(__file__)); out = os.path.join(here, "figures"); os.makedirs(out, exist_ok=True)
def save(fig, name):
    fig.tight_layout(); fig.savefig(os.path.join(out, name), dpi=130); plt.close(fig); print("wrote", name)
lg = math.log

def level_lefts(r, n, shifts=None):
    shifts = shifts if shifts is not None else [0, 1 - r]
    L = np.array([0.0])
    for _ in range(n):
        L = np.concatenate([r*L + s for s in shifts])
    return L

def draw_row(ax, r, y, n=7, h=0.6, color=BLUE, shifts=None):
    n = min(n, max(3, int(math.ceil(math.log(0.004)/math.log(r))))) if r < 0.5 else n
    L = level_lefts(r, n, shifts); w = r**n
    for x in L:
        ax.add_patch(Rectangle((x, y - h/2), w, h, color=color, alpha=0.35 if r > 0.5 else 0.9, lw=0))

# ---------------------------------------------------------------- hero
fig, (a, b) = plt.subplots(1, 2, figsize=(14, 6.2), gridspec_kw={"width_ratios": [1.25, 1]})
r = np.linspace(0.01, 0.95, 400)
for N, c, lab in ((8, LGREY, "8 pieces"), (4, LGREY, "4 pieces"), (3, LGREY, "3 pieces")):
    a.plot(r, np.log(N)/np.log(1/r), color=GREY, lw=1.5, alpha=0.7)
    xr = (1/N)**(1/3.3); a.text(xr - 0.012, 3.3, lab, color=GREY, fontsize=12, ha="right", va="center")
D2 = np.log(2)/np.log(1/r)
a.plot(r, D2, color=BLUE, lw=3, label="two pieces: D = log 2 / log(1/r)")
a.plot(r, np.minimum(1, D2), color=ORANGE, lw=2.2, ls="--", label="dimension of the set itself (Hutchinson)")
a.axvspan(0.5, 0.95, color=LGREY, alpha=0.35)
a.text(0.72, 0.15, "two pieces overlap\n(open set condition fails)", ha="center", fontsize=12, color=GREY)
stops = [(0.25, 0.5, "half box"), (1/3, lg(2)/lg(3), "Cantor"), (0.5, 1, "segment"), (2**-0.5, 2, "square\nedge-on"),
         ((math.sqrt(5) - 1)/2, lg(2)/lg(2/(math.sqrt(5) - 1)), "golden")]
for x, y, t in stops:
    a.plot(x, y, "o", color=BLUE if t not in ("golden",) else ORANGE, ms=10, mec="k")
    off = {"golden": (10, -18), "segment": (10, -20), "square\nedge-on": (-12, 10)}.get(t, (-12, 10))
    a.annotate(t, (x, y), xytext=off, textcoords="offset points", ha="left" if t in ("golden", "segment") else "right", fontsize=12)
for x, y, t in ((0.5, lg(3)/lg(2), "Sierpinski"), (1/3, lg(4)/lg(3), "Koch"), (0.5, 3, "cube")):
    a.plot(x, y, "s", color=GREY, ms=8, mec="k"); a.annotate(t, (x, y), xytext=(6, -4), textcoords="offset points", fontsize=12, color="k")
a.set_xlim(0, 0.95); a.set_ylim(0, 4); a.set_xlabel("similarity ratio r (the dial)"); a.set_ylabel("dimension D (a reading)")
a.set_title("Turning r: D comes out of Moran's equation"); a.legend(loc="upper left", frameon=False)
for i, rr in enumerate((0.2, 0.25, 1/3, 0.45, 0.5, 0.618)):
    draw_row(b, rr, -i, n=8 if rr < 0.5 else 7, color=BLUE if rr <= 0.5 else ORANGE)
    b.text(-0.03, -i, f"r = {rr:.3g}", ha="right", va="center", fontsize=13)
b.text(0.5, -6.1, "(each row drawn to the depth where pieces are still visible)", ha="center", fontsize=11, color=GREY)
b.set_xlim(-0.33, 1.02); b.set_ylim(-6.4, 0.6); b.axis("off")
b.set_title("Two pieces at ratio r\n(pale orange: pieces pile on top of each other)")
save(fig, "fig0-ratio-dial-hero.png")

# ---------------------------------------------------------------- sorting board
rows = json.load(open(os.path.join(here, "laws.json")))
order = ["slides", "partly", "plausible", "disguise", "self"]
cols = {"slides": BLUE, "partly": "#6fa0e6", "plausible": "#a9c4ee", "disguise": ORANGE, "self": "#2e2e2e", "no": LGREY}
fig, ax = plt.subplots(figsize=(14, 7.6))
y = 0
for v in order:
    items = [rw for rw in rows if rw["verdict"] == v]
    ax.text(0, y, f"{v} ({len(items)})", fontsize=15, weight="bold", color=cols[v], va="center"); y -= 0.75
    for rw in items:
        ax.add_patch(Rectangle((0.2, y - 0.32), 0.8, 0.64, color=cols[v]))
        ax.text(0.6, y, str(rw["n"]), ha="center", va="center", fontsize=12, color="white" if v in ("slides", "disguise", "self", "partly") else "k")
        lawtxt = rw["law"][:58] + ("  (ZoomDial lives here)" if rw["n"] == 4 else "")
        ax.text(1.15, y, lawtxt + (f"   [{rw['candidate']}]" if rw["candidate"] else ""), va="center", fontsize=12)
        y -= 0.75
    y -= 0.25
nos = [rw for rw in rows if rw["verdict"] == "no"]
ax.text(9.6, 0, f"no ({len(nos)})", fontsize=15, weight="bold", color=GREY, va="center")
for j, rw in enumerate(nos):
    cx, cy = 9.6 + (j % 6)*0.85, -0.75 - (j // 6)*0.85
    ax.add_patch(Rectangle((cx, cy - 0.38), 0.78, 0.76, color=LGREY))
    ax.text(cx + 0.39, cy, str(rw["n"]), ha="center", va="center", fontsize=12)
ax.text(9.6, -9.3, "no: a plane-dial or shape law, or no piece\nwith a size to turn", fontsize=12, color=GREY, va="top")
ax.set_xlim(-0.2, 14.9); ax.set_ylim(min(y, -10.6), 0.6); ax.axis("off")
ax.set_title("All 65 law groups, sorted by what turning the similarity ratio r does")
save(fig, "fig1-sorting-board.png")

# ---------------------------------------------------------------- C1
def boxdim(r, n=12):
    L = level_lefts(r, n); w = r**n
    def count(eps):
        a_ = np.floor(L/eps).astype(np.int64); b_ = np.floor((L + w)/eps).astype(np.int64)
        s = set()
        for lo, hi in zip(a_, b_): s.update(range(lo, hi + 1))
        return len(s)
    e1, e2 = r**(n - 4), r**n
    return (lg(count(e2)) - lg(count(e1)))/lg(e1/e2)
fig, (a, b) = plt.subplots(1, 2, figsize=(13, 5.6))
rg = np.linspace(0.12, 0.8, 18)
a.plot(r, D2, color=ORANGE, lw=2, label="similarity dimension log2/log(1/r) (Moran)")
a.plot(rg, [boxdim(x) for x in rg], "o", color=BLUE, ms=8, label="box count of the set (measured)")
a.axvline(0.5, color=GREY, ls=":"); a.text(0.51, 0.2, "r = 1/2:\nsegment", color=GREY, fontsize=12)
a.set_xlim(0.1, 0.82); a.set_ylim(0, 3.2); a.set_xlabel("similarity ratio r"); a.set_ylabel("dimension")
a.set_title("They agree up to r = 1/2, then part"); a.legend(frameon=False, loc="upper left")
for i, rr in enumerate((0.55, 0.6, 0.7)):
    L = level_lefts(rr, 12); xs = np.linspace(0, 1, 400)
    dens = np.histogram(L + rr**12/2, bins=xs)[0]
    b.plot(xs[:-1], dens/dens.max() + 1.2*i, color=ORANGE if i else BLUE, lw=1)
    b.text(1.02, 1.2*i + 0.3, f"r = {rr}", fontsize=13)
b.set_yticks([]); b.set_xlabel("position on the segment"); b.set_xlim(0, 1.25)
b.set_title("Past 1/2: how many depth-12 pieces\npile up at each point")
save(fig, "fig2-C1-two-pieces.png")

# ---------------------------------------------------------------- C2
from scipy.special import binom
fig, (a, b) = plt.subplots(1, 2, figsize=(13, 5.6))
ks = np.arange(0, 6); Ds = [0.5, 1, 1.5, 2, 3]; wbar = 0.16
for i, D in enumerate(Ds):
    c = GREY if float(D).is_integer() else BLUE
    a.bar(ks + (i - 2)*wbar, binom(D, ks), width=wbar, color=c, alpha=0.5 + 0.1*i, label=f"D = {D}")
a.axhline(0, color="k", lw=0.8); a.set_xlabel("corner layer k (pieces touching k sides)"); a.set_ylabel("layer size C(D, k)")
a.set_title("Box layers slide with D\n(grey: whole D, blue: fractional)"); a.legend(frameon=False, ncol=2)
for rr, c, lab in ((0.6, BLUE, "r = 0.6 (a = 2/3)"), (0.5, GREY, "r = 1/2 (a = 1)"), (1/3, ORANGE, "r = 1/3, Cantor (a = 2)")):
    D = lg(2)/lg(1/rr); aa = 1/rr - 1; kk = np.arange(0, 26)
    ps = np.cumsum(binom(D, kk)*aa**kk)
    b.plot(kk, ps, "o-", color=c, ms=4, label=lab)
b.axhline(2, color="k", ls=":"); b.text(25, 2.5, "N = 2 pieces", ha="right")
b.set_ylim(-3, 7); b.set_xlabel("layers added"); b.set_ylabel("running total")
b.set_title("Moran box: side 1/r holds N pieces,\nbut its layer sum converges only for r >= 1/2"); b.legend(frameon=False, loc="upper left")
save(fig, "fig3-C2-box-layers.png")

# ---------------------------------------------------------------- C3
fig, (a, b) = plt.subplots(1, 2, figsize=(13, 5.6))
rr = np.linspace(0.02, 0.5, 200)
for s, al in ((0.3, 0.5), (0.5, 0.75), (0.8, 1)):
    a.plot(rr, rr**s, color=ORANGE, alpha=al, lw=2, label=f"r^s at fixed s = {s} (k moves)")
a.plot(rr, np.full_like(rr, 0.5), color=BLUE, lw=3.5, label="r^D at the Moran D: stays 1/N = 1/2")
a.set_xlabel("similarity ratio r (two pieces, both turned)"); a.set_ylabel("weight of one piece")
a.set_title("Uniform turning: k is frozen at the true D"); a.legend(frameon=False, fontsize=11, loc="lower right")
from scipy.optimize import brentq
r1s = np.linspace(0.02, 0.58, 120); r2 = 0.4
Dm = [brentq(lambda s: x**s + r2**s - 1, 1e-6, 50) for x in r1s]
p1 = [x**d for x, d in zip(r1s, Dm)]
b.plot(r1s, p1, color=BLUE, lw=3, label="weight of the turned piece r1^D")
b.plot(r1s, 1 - np.array(p1), color=GREY, lw=2, label="weight of the other piece r2^D")
b.plot(r1s, Dm, color=ORANGE, lw=2, ls="--", label="D (reading)")
b.axvline(r2, color=GREY, ls=":"); b.text(r2 + 0.01, 0.3, "r1 = r2", color=GREY)
b.set_xlabel("ratio of one piece, r1 (other fixed at 0.4)"); b.set_ylabel("weight / D")
b.set_title("One piece alone: weights move on the simplex\n(the odds (r1/r2)^D, not k)"); b.legend(frameon=False, fontsize=11)
save(fig, "fig4-C3-weights.png")

# ---------------------------------------------------------------- C4
fig, (a, b) = plt.subplots(1, 2, figsize=(13, 5.8))
r1, r2 = 0.3, 0.6
D = brentq(lambda s: r1**s + r2**s - 1, 1e-6, 50)
tree = {"": None, "0": None, "1": None, "00": None, "01": None, "10": None, "11": None, "010": None, "011": None}
leaves = ["00", "010", "011", "10", "11"]
def pos(code):
    x = 0.5; dx = 0.25
    for ch in code: x += (-dx if ch == "0" else dx); dx /= 2
    return x, -len(code)
for code in tree:
    if code:
        x, y = pos(code); px, py = pos(code[:-1])
        a.plot([px, x], [py, y], color=GREY, lw=2)
        a.text((px + x)/2 + (-0.03 if code[-1] == "0" else 0.03), (py + y)/2, f"{r1 if code[-1]=='0' else r2}",
               fontsize=11, color=ORANGE, ha="center")
tot = 0
for lf in leaves:
    x, y = pos(lf); wgt = np.prod([(r1 if ch == "0" else r2)**D for ch in lf]); tot += wgt
    a.plot(x, y, "o", color=BLUE, ms=12); a.text(x, y - 0.35, f"{wgt:.3f}", ha="center", fontsize=12, color=BLUE)
a.text(0.5, -4.1, f"leaf weights add to {tot:.6f} at D = {D:.3f}\n(Kraft's 1 is the case r1 = r2 = 1/2, D = 1)", ha="center", fontsize=12)
a.set_xlim(0, 1); a.set_ylim(-4.6, 0.4); a.axis("off"); a.set_title("Every full tree weighs 1 at Moran's D")
q = np.linspace(-3, 6, 200); p = np.array([0.5, 0.3, 0.2])
def renyi(q):
    return np.where(np.abs(q - 1) < 1e-9, -np.sum(p*np.log(p)), np.log(np.sum(p[None, :]**q[:, None], axis=1))/(1 - q + 1e-300))
for rr, c in ((0.25, BLUE), (0.15, "#6fa0e6")):
    b.plot(q, renyi(q)/np.log(1/rr), color=c, lw=2.5, label=f"weights (0.5, 0.3, 0.2), r = {rr}: H_q / log(1/r)")
b.axhline(np.log(3)/np.log(4), color=ORANGE, lw=2, ls="--", label="natural weights 1/3 each, r = 1/4: flat at D")
b.set_xlabel("order q"); b.set_ylabel("Renyi dimension D_q"); b.legend(frameon=False, fontsize=11)
b.set_title("Renyi entropy, divided by log(1/r),\nis the Renyi dimension (uniform r)")
save(fig, "fig5-C4-kraft-renyi.png")

# ---------------------------------------------------------------- C5
fig, ax = plt.subplots(figsize=(11, 6))
ms = np.arange(1, 9); rstops = 1/(ms + 1)
rc = np.linspace(0.001, 0.5, 300)
for beta, al in ((1, 0.55), (2, 0.75), (3, 1.0)):
    ax.plot(rc, 1 - rc**beta, color=GREY, ls=":", lw=1.5)
    ax.plot(rstops, 1 - rstops**beta, "o", color=BLUE, ms=9, alpha=al, label=f"growth degree beta = {beta}")
ax.plot(0.5, 7/8, "o", color=ORANGE, ms=15, mfc="none", mew=3)
ax.annotate("7/8: Stefan-Boltzmann\n(fermions / bosons, beta = 3)", (0.5, 7/8), xytext=(0.3, 0.6), fontsize=13,
            arrowprops=dict(arrowstyle="->", color=ORANGE))
for m, x in zip(ms[:4], rstops[:4]):
    ax.text(x, 1.015, f"cap {m}", ha="center", fontsize=11)
ax.text(0.005, 1.015, "Bose", fontsize=11)
ax.set_xlabel("similarity ratio r = 1/(m+1) of the leftover copy"); ax.set_ylabel("share kept by the capped count  1 - r^beta")
ax.set_xlim(-0.01, 0.53); ax.set_ylim(0.45, 1.06); ax.legend(frameon=False, loc="lower left")
ax.set_title("Capping copies at m leaves the same set at size (m+1):\nshare r^beta; only the stops r = 1/(m+1) count anything")
save(fig, "fig6-C5-cap-ladder.png")

# ---------------------------------------------------------------- C6
fig, ax = plt.subplots(figsize=(11, 5.8))
lr = np.linspace(0, 2.5, 100)
ax.plot(lr, lr/3, color=BLUE, lw=3, label="clock r: rho = (ln r)/3 (newton-g Q1)")
for name, val in (("r = 2", lg(2)), ("r = e", 1), ("r = 6", lg(6)), ("r = e^(9/5)", 9/5)):
    ax.plot(val, val/3, "o", color=BLUE, ms=9)
    ax.annotate(name, (val, val/3), xytext=(-70, 4) if name == "r = 6" else (6, -18), textcoords="offset points", fontsize=12)
ax.axhline(0.6, color=ORANGE, ls="--"); ax.text(0.05, 0.62, "gravity's two-loop 3/5 (one scheme)", color=ORANGE, fontsize=12)
ax2 = ax.twiny(); ax2.set_xlim(ax.get_xlim()); ticks = [0, 0.5, 1, 1.5, 2, 2.5]
ax2.set_xticks(ticks); ax2.set_xticklabels([f"{(1 - t)/4:.2f}" for t in ticks]); ax2.set_xlabel("the same rho from the curvature dial: d / k^2")
ax.set_xlabel("ln r (e-folds per height step)"); ax.set_ylabel("two-loop number rho"); ax.legend(frameon=False, loc="lower right")
ax.set_title("The clock is the curvature dial in disguise (for rho)", pad=12)
save(fig, "fig7-C6-clock.png")

# ---------------------------------------------------------------- complex dimensions
fig, axs = plt.subplots(1, 2, figsize=(13, 6))
X, Y = np.meshgrid(np.linspace(-1.5, 1.2, 420), np.linspace(0, 60, 900))
S = X + 1j*Y
for ax, (ra, rb, title, c) in zip(axs, ((0.25, 0.5, "lattice: ratios 1/4 and 1/2\n(log ratios rationally related)", ORANGE),
                                       (0.3, 0.5, "non-lattice: ratios 0.3 and 0.5", BLUE))):
    F = np.abs(ra**S + rb**S - 1)
    ax.imshow(np.log10(F + 1e-6), extent=[-1.5, 1.2, 0, 60], origin="lower", aspect="auto", cmap="Greys_r", vmin=-3, vmax=0.5)
    D = brentq(lambda s: ra**s + rb**s - 1, 1e-6, 5)
    ax.axvline(D, color=c, lw=1.5, ls="--"); ax.text(D + 0.03, 57, f"Re s = D = {D:.3f}", color=c, fontsize=12)
    ax.set_xlabel("real part of s"); ax.set_ylabel("imaginary part of s"); ax.set_title(title)
fig.suptitle("Complex dimensions (dark dots: zeros of 1 - r1^s - r2^s). Turning one ratio alone flips the pattern.", fontsize=14)
save(fig, "fig8-complex-dimensions.png")

# ---------------------------------------------------------------- N: Bernoulli convolutions
def bern_hist(lam, n=20, bins=500):
    v = np.zeros(1)
    for i in range(n):
        v = np.concatenate([v, v + lam**i])
    tot = 1/(1 - lam)
    h, e = np.histogram(v, bins=bins, range=(0, tot)); return e[:-1]/tot, h/h.max()
fig, axs = plt.subplots(1, 3, figsize=(15, 5.2))
for i, (lam, c, lab) in enumerate((((math.sqrt(5) - 1)/2, ORANGE, "r = 1/phi = 0.618 (golden, Pisot):\nsingular, D = 0.9957"),
                                   (0.65, GREY, "r = 0.65: typical\n(absolutely continuous a.e., Solomyak)"),
                                   (2**-0.5, BLUE, "r = 1/sqrt 2 = 0.707:\nsquare seen edge-on, trapezoid"))):
    x, h = bern_hist(lam); axs[0].plot(x, h + 1.25*i, color=c, lw=0.9); axs[0].text(1.04, 1.25*i + 0.35, lab, fontsize=10.5, ha="left", color=c)
axs[0].set_yticks([]); axs[0].set_xlim(0, 2.05); axs[0].set_xticks([0, 0.5, 1]); axs[0].set_xlabel("position (scaled to [0, 1])")
axs[0].set_title("The natural measure at three ratios")
for sp_ in ("top", "right", "left"): axs[0].spines[sp_].set_visible(False)
rr = np.linspace(0.3, 0.95, 300)
axs[1].plot(rr, np.minimum(1, np.log(2)/np.log(1/rr)), color=GREY, lw=2.5, label="min(1, log2/log(1/r)): every r\nwithout exact overlaps (Hochman; Varju)")
for lam, d, t, off in (((math.sqrt(5) - 1)/2, 0.995713, "golden 0.9957", (10, -28)), (0.543689, 0.980409, "tribonacci 0.9804", (-40, -34))):
    axs[1].plot(lam, d, "o", color=BLUE, ms=10, mec="k"); axs[1].annotate(t, (lam, d), xytext=off, textcoords="offset points", fontsize=12)
axs[1].set_ylim(0.6, 1.03); axs[1].set_xlabel("similarity ratio r (two pieces)"); axs[1].set_ylabel("dimension of the natural measure")
axs[1].set_title("Exact overlaps drop the dimension\n(blue: recomputed here)"); axs[1].legend(frameon=False, fontsize=10.5, loc="lower right")
# square edge-on: the four corners of a square grid projected onto the line with slope sqrt2
pts = np.array(list(itertools.product(range(4), repeat=2)), float)
axs[2].scatter(pts[:, 0], pts[:, 1], color=BLUE, s=60)
d = np.array([1, 2**-0.5]); d /= np.linalg.norm(d)
for p_ in pts:
    proj = p_.dot(d)*d; axs[2].plot([p_[0], proj[0]], [p_[1], proj[1]], color=LGREY, lw=1)
axs[2].plot([0, 4.2*d[0]], [0, 4.2*d[1]], color=ORANGE, lw=2.5)
axs[2].set_aspect("equal"); axs[2].set_xlim(-0.5, 4); axs[2].set_ylim(-0.5, 3.6); axs[2].axis("off")
axs[2].set_title("r = 2^(-1/2): the square's\npieces seen along one line")
save(fig, "fig9-N-bernoulli.png")

# ---------------------------------------------------------------- N: sum sets and the four-corner set
fig, (a, b) = plt.subplots(1, 2, figsize=(14, 5.8))
for i, rr in enumerate((0.25, 0.3, 1/3, 0.4)):
    L = np.unique(np.round(level_lefts(rr, 6, [0, 1 - rr, 2*(1 - rr)]), 12)); w = 2*rr**6
    for x in L:
        a.add_patch(Rectangle((x, -i - 0.3), w, 0.6, color=BLUE if rr < 1/3 else ORANGE, lw=0))
    a.text(-0.05, -i, f"r = {rr:.3g}", ha="right", va="center", fontsize=13)
a.set_xlim(-0.45, 2.05); a.set_ylim(-3.7, 0.6); a.axis("off")
a.set_title("C_r + C_r: gaps until r = 1/3, then a full interval")
cq = [sum(3*e/4**(j + 1) for j, e in enumerate(eps)) for eps in itertools.product((0, 1), repeat=3)]
for x in cq:
    for y in cq:
        b.add_patch(Rectangle((x, y), 4**-3, 4**-3, color=BLUE))
for x in {round((2*a_ + c_)/3, 12) for a_ in cq for c_ in cq}:
    b.add_patch(Rectangle((x, -0.28), 4**-3, 0.12, color=GREY, lw=0))
for x in {round((a_ + c_)/2, 12) for a_ in cq for c_ in cq}:
    b.add_patch(Rectangle((x, -0.52), 4**-3, 0.12, color=ORANGE, lw=0))
b.text(1.08, -0.22, "shadow along (2, 1):\nthe whole segment", fontsize=12, color=GREY, va="center")
b.text(1.08, -0.46, "shadow along (1, 1):\na dust, D = log3/log4", fontsize=12, color=ORANGE, va="center")
b.text(1.08, 0.5, "four-corner set\n= half box x half box\n(depth 3)", fontsize=12, color=BLUE, va="center")
b.set_xlim(-0.05, 1.75); b.set_ylim(-0.6, 1.05); b.set_aspect("equal"); b.axis("off")
b.set_title("Two half boxes joined: the four-corner set.\nIts shadow depends on the direction")
save(fig, "fig10-N-sumsets.png")
