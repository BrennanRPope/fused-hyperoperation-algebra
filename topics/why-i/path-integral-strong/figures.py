"""Pictures for path-integral-strong. Run: python3 figures.py  (writes figures/*.png)"""
import numpy as np, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon, Circle, FancyArrowPatch, Rectangle
from scipy.linalg import expm, eigh
from scipy.special import fresnel
from scipy.integrate import quad
from scipy.optimize import brentq
import itertools, os
rng = np.random.default_rng(1125)
OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "figures"); os.makedirs(OUT, exist_ok=True)
BLUE, ORANGE, GREEN, RED, GREY, PURPLE = "#2b6cb0", "#dd6b20", "#2f855a", "#c53030", "#718096", "#6b46c1"
plt.rcParams.update({"font.size": 10, "axes.titlesize": 11, "figure.dpi": 130})
def save(fig, name):
    fig.tight_layout(); fig.savefig(os.path.join(OUT, name + ".png"), bbox_inches="tight"); plt.close(fig)
def arrow(ax, p, q, c="k", lw=1.5, ms=12):
    ax.add_patch(FancyArrowPatch(p, q, arrowstyle="-|>", mutation_scale=ms, color=c, lw=lw))
def hdist(z1, z2): return np.arccosh(1 + np.abs(z1 - z2)**2/(2*z1.imag*z2.imag))

# ---------------- GAP 1: DYNAMICS ----------------
fig, ax = plt.subplots(1, 2, figsize=(11, 4.3))
a_, v_ = 2.0, 1.4
ax[0].add_patch(Rectangle((0, 0), a_, v_, color=ORANGE, alpha=.25))
arrow(ax[0], (0, 0), (a_, 0), BLUE); arrow(ax[0], (a_, 0), (a_, v_), BLUE)
arrow(ax[0], (0, 0), (0, v_), GREEN); arrow(ax[0], (0, v_), (a_, v_), GREEN)
ax[0].text(a_/2, -.18, "shift by a", ha="center", color=BLUE); ax[0].text(a_ + .05, v_/2, "boost by v", color=BLUE)
ax[0].text(-.05, v_/2, "boost by v", ha="right", color=GREEN); ax[0].text(a_/2, v_ + .08, "shift by a", ha="center", color=GREEN)
ax[0].text(a_/2, v_/2, "the two orders end in the same place\nbut differ by the phase  e^(i m v a)\n(area m v a, in units of hbar)", ha="center", va="center")
ax[0].set_xlim(-.9, 2.9); ax[0].set_ylim(-.45, 1.85); ax[0].axis("off")
ax[0].set_title("Bargmann (1954): the boost/shift twist has a size, and the size is the mass")
pg = np.linspace(-2, 2, 200)
for mm, c in [(0.5, ORANGE), (1.0, BLUE), (2.0, GREEN)]:
    ax[1].plot(pg, pg**2/(2*mm), color=c, label=f"twist size m = {mm}")
ax[1].set_ylim(0, 4); ax[1].set_xlabel("momentum p"); ax[1].set_ylabel("energy H")
ax[1].legend(loc="upper center"); ax[1].set_title("shift symmetry + boost rule force H = p^2/(2m) (check K2)\nfree motion is picked out; a potential V(x) is not")
save(fig, "g1-mass-is-the-boost-shift-twist")

fig, ax = plt.subplots(1, 2, figsize=(11, 4.6))
for A in range(5):
    for Bx in range(5):
        inq = A + Bx <= 2
        ax[0].plot(A, Bx, "o", ms=11 if inq else 7, color=GREEN if inq else GREY)
ax[0].add_patch(Polygon([(-.4, -.4), (2.75, -.4), (-.4, 2.75)], color=GREEN, alpha=.12))
ax[0].text(.25, .25, "degree <= 2\ncloses", color=GREEN)
arrow(ax[0], (3, 0), (2.1, 1.9), RED); arrow(ax[0], (0, 3), (1.9, 2.1), RED)
ax[0].text(2.15, 2.2, "{q^3, p^3} = 9 q^2 p^2\nlands at degree 4", color=RED)
ax[0].set_xlabel("power of q  (a)"); ax[0].set_ylabel("power of p  (b)"); ax[0].set_xticks(range(5)); ax[0].set_yticks(range(5))
ax[0].set_title("which generators keep the chain closed (Groenewold-van Hove)")
ax[1].add_patch(Polygon([(0, 0), (2, 0), (2 + 1, 0 + 2), (1, 2)], color=ORANGE, alpha=.3))
arrow(ax[1], (0, 0), (2, 0), BLUE); arrow(ax[1], (0, 0), (1, 2), BLUE)
ax[1].text(1.0, -.25, "(a, b) = (2, 0): q^2", ha="center", color=BLUE); ax[1].text(-.2, 2.2, "(c, d) = (1, 2): q p^2", color=BLUE)
ax[1].text(1.5, 1.0, "area = a d - b c = 4", ha="center")
ax[1].text(1.5, 2.75, "{q^a p^b, q^c p^d} = (a d - b c) q^(a+c-1) p^(b+d-1)\nthe order dial's wedge, on exponents (check Q1)", ha="center")
ax[1].set_xlim(-.4, 3.4); ax[1].set_ylim(-.5, 3.3); ax[1].set_aspect("equal"); ax[1].axis("off")
ax[1].set_title("the bracket's coefficient is the same wedge as the order gap")
save(fig, "g2-quadratics-close")

w, T = 1.3, 2*np.pi/1.3
fig, ax = plt.subplots(1, 2, figsize=(11, 4.3))
th = np.linspace(0, 2*np.pi, 300)
ax[0].plot(np.cos(th), -w*np.sin(th), color=GREY, lw=3, alpha=.5, label="exact flow (one period)")
for N, c in [(8, RED), (40, BLUE)]:
    dt = T/N; half = np.array([[1, 0], [-w*w*dt/2, 1.]]); drift = np.array([[1, dt], [0, 1.]]); step = half @ drift @ half
    z = np.array([1., 0]); pts = [z]
    for _ in range(N): z = step @ z; pts.append(z)
    pts = np.array(pts); ax[0].plot(pts[:, 0], pts[:, 1], "o-", ms=3, color=c, label=f"chain of {N} kick-drift-kick steps")
ax[0].set_xlabel("x"); ax[0].set_ylabel("p"); ax[0].set_aspect("equal"); ax[0].legend(fontsize=8, loc="lower left")
ax[0].set_title("oscillator: the chain is a product of 2x2 area-keeping matrices")
Ns = np.array([10, 30, 100, 300, 1000]); T2 = 2.0
exact = np.array([[np.cos(w*T2), np.sin(w*T2)/w], [-w*np.sin(w*T2), np.cos(w*T2)]])
errs = []
for N in Ns:
    dt = T2/N; half = np.array([[1, 0], [-w*w*dt/2, 1.]]); drift = np.array([[1, dt], [0, 1.]])
    errs.append(np.linalg.norm(np.linalg.matrix_power(half @ drift @ half, N) - exact))
ax[1].loglog(Ns, errs, "o-", color=BLUE, label="chain error"); ax[1].loglog(Ns, errs[0]*(Ns[0]/Ns)**2, "--", color=GREY, label="slope -2")
ax[1].set_xlabel("number of slices N"); ax[1].set_ylabel("distance from exact flow"); ax[1].legend()
ax[1].set_title("quadratic world: slices converge, no measure needed (check Q5)")
save(fig, "g3-quadratic-chain-is-exact")

# ---------------- GAP 2: hbar ----------------
fig, ax = plt.subplots(1, 2, figsize=(11, 4.4))
hbar = 1.0
for n in range(5):
    E = hbar*w*(n + .5); r = np.sqrt(2*E)
    ax[0].plot(r/w*np.cos(th), r*np.sin(th), color=plt.cm.viridis(n/5))
    ax[0].text(r/w*np.cos(.7) + .05, r*np.sin(.7), f"n={n}", fontsize=8, color=plt.cm.viridis(n/5))
ax[0].set_aspect("equal"); ax[0].set_xlabel("x"); ax[0].set_ylabel("p")
ax[0].set_title("allowed orbits enclose area 2 pi hbar (n + 1/2)")
L, Ng = 24.0, 1200; xs = np.linspace(-L/2, L/2, Ng); dx = xs[1] - xs[0]
lap = (np.diag(np.ones(Ng-1), 1) - 2*np.eye(Ng) + np.diag(np.ones(Ng-1), -1))/dx**2
E4 = eigh(-0.5*lap + np.diag(xs**4), eigvals_only=True)[:6]
def area4(En):
    xt = En**0.25; return 2*quad(lambda X: np.sqrt(max(0, 2*(En - X**4))), -xt, xt)[0]
bs4 = np.array([brentq(lambda En: area4(En) - 2*np.pi*(n + .5), 1e-3, 50) for n in range(6)])
Eo = eigh(-0.5*lap + np.diag(0.5*w*w*xs**2), eigvals_only=True)[:6]
ax[1].semilogy(range(6), np.abs(Eo - w*(np.arange(6) + .5))/Eo + 1e-16, "o-", color=BLUE, label="oscillator: area rule exact (left at grid error)")
ax[1].semilogy(range(6), np.abs(bs4 - E4)/E4, "s-", color=ORANGE, label="x^4: area rule off 18% at n=0, improves")
ax[1].set_xlabel("level n"); ax[1].set_ylabel("relative error of the area rule"); ax[1].legend(fontsize=8)
ax[1].set_title("Bohr-Sommerfeld + Maslov half (checks B1, B2)")
save(fig, "h1-bohr-sommerfeld-rings")

fig = plt.figure(figsize=(11, 4.4))
ax3 = fig.add_subplot(1, 2, 1, projection="3d")
j = 2; nb = int(2*j + 1)
u = np.linspace(0, 2*np.pi, 60)
for k in range(nb):
    z = np.linspace(-1 + 2*k/nb, -1 + 2*(k + 1)/nb, 12); U, Z = np.meshgrid(u, z); R = np.sqrt(1 - Z**2)
    ax3.plot_surface(R*np.cos(U), R*np.sin(U), Z, color=plt.cm.Set2(k), alpha=.85, linewidth=0)
ax3.set_box_aspect((1, 1, 1)); ax3.axis("off")
ax3.set_title(f"spin j = {j}: sphere of area 4 pi hbar (j + 1/2)\nsplits into {nb} equal bands of area 2 pi hbar")
ax = fig.add_subplot(1, 2, 2)
sq = np.array([[0, 0], [2, 0], [2, np.pi], [0, np.pi], [0, 0]])
ax.fill(sq[:, 0], sq[:, 1], color=ORANGE, alpha=.2); ax.plot(sq[:, 0], sq[:, 1], color=ORANGE)
for i in range(4): arrow(ax, sq[i], sq[i] + .55*(sq[i + 1] - sq[i]), ORANGE)
ax.text(1, np.pi/2, "closed chain of steps\nenclosed area 2 pi", ha="center")
for cx, hb, lab in [(3.6, 1.0, "hbar = 1: phase = 1\n(loop is invisible)"), (5.6, 0.7, "hbar = 0.7: phase != 1")]:
    ax.add_patch(Circle((cx, 1.6), .7, fill=False, color=GREY)); ph = 2*np.pi/hb
    arrow(ax, (cx, 1.6), (cx + .7*np.cos(ph), 1.6 + .7*np.sin(ph)), GREEN if hb == 1 else RED)
    ax.text(cx, .55, lab, ha="center", fontsize=9)
ax.set_xlim(-.3, 6.6); ax.set_ylim(-.2, 3.6); ax.set_aspect("equal"); ax.axis("off")
ax.set_title("the rule fixes area in units of hbar, never hbar itself (B3-B5)")
save(fig, "h2-sphere-cells-and-loop-phase")

fig, ax = plt.subplots(1, 2, figsize=(11, 4.3))
bet = np.linspace(0.05, 3, 400)
for g, c in [(2, BLUE), (3, GREEN), (4, ORANGE)]:
    ax[0].plot(bet, (2*g - 2)*bet, color=c, label=f"genus {g}")
    bi = np.arange(1, int(3*(2*g - 2)) + 1)/(2*g - 2); bi = bi[bi <= 3]
    ax[0].plot(bi, (2*g - 2)*bi, "o", color=c, ms=4)
ax[0].set_xlabel("beta  (Klauder's dimensionless scale)"); ax[0].set_ylabel("number of hbar-cells = (2g - 2) beta")
ax[0].legend(); ax[0].set_title("[ours] closed half-plane surface: dots = allowed beta (check B6)")
ax[1].axis("off")
ax[1].text(0, .95, "Klauder's affine phase space (curvature -1/(beta hbar), check M6)", fontsize=10, weight="bold")
ax[1].text(0, .78, "area element = dp dq  (metric area = symplectic area)")
ax[1].text(0, .64, "Gauss-Bonnet, genus g:   area = 2 pi (2g - 2) beta hbar")
ax[1].text(0, .50, "Kostant-Souriau:   area / (2 pi hbar) must be a whole number")
ax[1].text(0, .36, "=>   (2g - 2) beta  in  Z      hbar cancels", color=RED)
ax[1].text(0, .18, "Reading: geometry can fix the ratio (curvature radius)^2 / hbar,\nnot hbar. Needs isometries beyond the order dial (SL(2,R)),\nand physical phase space is flat (beta -> infinity).", color=GREY)
ax[1].set_title("the curvature trade")
save(fig, "h3-the-curvature-trade")

# ---------------- GAP 3: MEASURE ----------------
fig, ax = plt.subplots(1, 2, figsize=(11, 4.3))
Np, Ns_ = 40, 200; dt = 1/Ns_; tg = np.linspace(0, 1, Ns_ + 1)
B = 0.4 + np.concatenate([np.zeros((Np, 1)), np.cumsum(rng.normal(scale=np.sqrt(dt), size=(Np, Ns_)), axis=1)], axis=1)
wts = np.exp(-dt*np.sum(0.5*B**2, axis=1)); wn = (wts - wts.min())/(wts.max() - wts.min())
for i in range(Np): ax[0].plot(tg, B[i], color=plt.cm.plasma(wn[i]), lw=.9)
ax[0].plot(tg, 0*tg, color=GREY, ls=":")
ax[0].set_xlabel("time"); ax[0].set_ylabel("position"); ax[0].set_title("real dial: Brownian paths, each weighted by exp(-int V dt)\n(bright = heavier; paths far from V's floor fade)")
ax[1].bar([0, 1], [0.3834, 0.3836], yerr=[0.0014, 0], color=[BLUE, GREEN], width=.5)
ax[1].set_xticks([0, 1]); ax[1].set_xticklabels(["40,000 weighted paths", "exact e^(-tH) f"]); ax[1].set_ylim(.37, .39)
ax[1].set_title("Feynman-Kac: a true measure on paths (Kac 1949; check M1)")
save(fig, "m1-feynman-kac-paths")

fig, ax = plt.subplots(1, 2, figsize=(11, 4.6))
s = np.linspace(-6, 6, 600)
for k in [-2, -1, -.5, .5, 1, 2]:
    lam = 1 + k*s; m_ = lam >= 0
    ax[0].plot(s[m_], lam[m_], color=GREY, lw=.8); ax[0].plot(-1/k, 0, "x", color=GREY, ms=5)
def hyp_bm(z0, T_, n_, N_):
    xx = np.full(N_, z0.real); yy = np.full(N_, z0.imag); dt_ = T_/n_; tr = [xx + 1j*yy]
    for _ in range(n_):
        dW = rng.normal(scale=np.sqrt(dt_), size=(2, N_)); xx = xx + yy*dW[0]; yy = yy*np.exp(dW[1] - .5*dt_); tr.append(xx + 1j*yy)
    return np.array(tr)
for z0, c in [(1j, BLUE), (3 + .25j, ORANGE)]:
    tr = hyp_bm(z0, 1.0, 400, 6)
    for i in range(6): ax[0].plot(tr[:, i].real, tr[:, i].imag, color=c, lw=.7, alpha=.8)
    ax[0].plot(z0.real, z0.imag, "o", color=c)
ax[0].text(0.08, 1.06, "i", fontsize=12); ax[0].set_xlim(-4, 5); ax[0].set_ylim(0, 3.2)
ax[0].set_xlabel("shift b"); ax[0].set_ylabel("stretch a  (x -> a x + b)")
ax[0].set_title("Brownian motion on the order dial's half-plane\n(grey: dial lines through i; x: their centres -1/k)")
Z1 = hyp_bm(1j, 1.0, 400, 20000)[-1]; Z2 = hyp_bm(3 + .25j, 1.0, 400, 20000)[-1]
ax[1].hist(hdist(Z1, 1j), bins=60, density=True, alpha=.5, color=BLUE, label="start at i")
ax[1].hist(hdist(Z2, 3 + .25j), bins=60, density=True, alpha=.5, color=ORANGE, label="start at 3 + 0.25 i")
ax[1].set_xlabel("hyperbolic distance travelled by time 1"); ax[1].legend()
ax[1].set_title("same law from every start: it is the group's own measure\n(Klauder metric = Poincare metric, checks M2-M5)")
save(fig, "m2-brownian-motion-on-the-dial-half-plane")

fig, ax = plt.subplots(figsize=(9.5, 5.2))
ax.axhline(0, color="k", lw=.8); ax.axvline(0, color="k", lw=.8)
ax.plot([0, 3.2], [0, 0], color=GREEN, lw=4); ax.plot([0, 0], [0, -3.2], color=RED, lw=4)
ax.text(1.6, .25, "k real (k = -1): Wiener measure exists\nFeynman-Kac, rigorous", ha="center", color=GREEN)
ax.text(-.15, -2.8, "k = -i/hbar: no countably additive\nmeasure on paths (Cameron 1960)", color=RED, ha="right")
th_ = np.linspace(0, -np.pi/2, 50); ax.plot(2.2*np.cos(th_), 2.2*np.sin(th_), "--", color=BLUE)
arrow(ax, (2.2*np.cos(-1.2), 2.2*np.sin(-1.2)), (2.2*np.cos(-1.35), 2.2*np.sin(-1.35)), BLUE)
ax.text(1.75, -1.65, "reach it by turning the dial:\nNelson 1964 (continue in mass),\nOsterwalder-Schrader (QFT)", color=BLUE, fontsize=9)
ax.text(-3.4, 2.4, "Fresnel integrals (Ito 1965; Albeverio-Hoegh-Krohn 1976):\ndefine the k = -i/hbar sum as a limit, not a measure", fontsize=9)
ax.text(-3.4, 1.2, "Daubechies-Klauder 1985 / DKP 1987: keep a real\nBrownian motion in phase space, speed nu -> infinity.\nFor the affine group it runs on the half-plane (m2).", fontsize=9, color=PURPLE)
ax.set_xlim(-3.6, 3.6); ax.set_ylim(-3.4, 3.2); ax.set_xlabel("Re k"); ax.set_ylabel("Im k")
ax.set_title("the measure question on the k-dial: where a real measure exists and how each route reaches k = -i/hbar")
save(fig, "m3-measure-routes-on-the-k-dial")

# ---------------- SHAPES ----------------
fig, ax = plt.subplots(1, 2, figsize=(11, 4.4))
steps = np.array([1, 1, -1, -1, -1, 1, 1, 1, -1, 1]); xp = np.concatenate([[0], np.cumsum(steps)]); tp = np.arange(len(xp))
ax[0].plot(xp, tp, "o-", color=BLUE)
for i in range(1, len(steps)):
    if steps[i] != steps[i - 1]: ax[0].plot(xp[i], tp[i], "s", ms=11, mfc="none", color=ORANGE)
ax[0].text(3.2, 1, "each turn (square) multiplies by\n  i m eps   (Feynman checkerboard)\n  a eps     (Kac telegraph, real)", color=ORANGE, fontsize=9)
ax[0].set_xlabel("position"); ax[0].set_ylabel("time step"); ax[0].set_xlim(-3, 7)
ax[0].set_title("elsewhere: two-way walk, weight per turn")
pmom, mm = 1.1, .8; sz = np.diag([1., -1.]); sx = np.array([[0, 1.], [1, 0]]); tt = np.linspace(0, 4, 120)
D = [np.linalg.norm(expm(t*(-1j*pmom*sz + 1j*mm*sx))[:, 0]) for t in tt]
Tl = [np.linalg.norm(expm(t*(-1j*pmom*sz + mm*sx))[:, 0]) for t in tt]
ax[1].plot(tt, D, color=BLUE, label="turn weight i m: Dirac, length kept (unitary)")
ax[1].plot(tt, Tl, color=ORANGE, label="turn weight m (real): telegraph, size not kept")
ax[1].set_xlabel("time"); ax[1].set_ylabel("size of the state"); ax[1].legend(fontsize=8)
ax[1].set_title("path-integral version: same chain, dial turned (check S2)")
save(fig, "s1-checkerboard-and-telegraph")

fig, ax = plt.subplots(1, 2, figsize=(11, 4.4))
NT = 12; ang = 2*np.pi*np.arange(NT)/NT; spins = rng.choice([1, -1], NT)
ax[0].add_patch(Circle((0, 0), 1, fill=False, color=GREY))
for a__, s__ in zip(ang, spins):
    arrow(ax[0], (np.cos(a__), np.sin(a__) - .15*s__), (np.cos(a__), np.sin(a__) + .15*s__), BLUE if s__ > 0 else RED, ms=9)
ax[0].text(0, 0, "classical Ising ring\none site per time slice\ncoupling K = -(1/2) ln tanh(eps h)", ha="center", va="center", fontsize=9)
ax[0].set_xlim(-1.5, 1.5); ax[0].set_ylim(-1.4, 1.4); ax[0].set_aspect("equal"); ax[0].axis("off")
ax[0].set_title("elsewhere: a classical chain one dimension up")
hT, gT, bT = .7, .4, 2.0; Zex = np.trace(expm(bT*(hT*sx + gT*sz)))
Nl = np.arange(2, 40); Zs = [np.trace(np.linalg.matrix_power(expm(bT/N*hT*sx) @ expm(bT/N*gT*sz), N)) for N in Nl]
ax[1].plot(Nl, Zs, "o-", ms=3, color=BLUE, label="time-sliced quantum trace (= ring sum, S3)"); ax[1].axhline(Zex, color=GREEN, ls="--", label="Tr exp(beta H)")
ax[1].set_xlabel("slices N"); ax[1].set_ylabel("partition function"); ax[1].legend(fontsize=8)
ax[1].set_title("path-integral version: Suzuki-Trotter slices")
save(fig, "s2-trotter-slices-are-a-classical-chain")

fig, ax = plt.subplots(1, 2, figsize=(11, 4.4))
for i in range(6):
    ax[0].plot([0, 5], [i, i], color=GREY, lw=.6); ax[0].plot([i, i], [0, 5], color=GREY, lw=.6)
for i in range(1, 4):
    for j_ in range(1, 4): ax[0].add_patch(Rectangle((i, j_), 1, 1, color=ORANGE, alpha=.15 + .1*((i + j_) % 2)))
loopp = [(1, 1), (4, 1), (4, 4), (1, 4), (1, 1)]
for p0, p1 in zip(loopp[:-1], loopp[1:]): arrow(ax[0], p0, p1, BLUE, lw=2)
ax[0].text(2.5, 2.5, "loop phase =\nproduct of the\nplaquette phases\n(check S4)", ha="center", va="center", fontsize=9)
ax[0].set_aspect("equal"); ax[0].axis("off"); ax[0].set_title("elsewhere: Wilson loop on a lattice (holonomy)")
ax[1].add_patch(Circle((0, 0), .35, color=PURPLE, alpha=.4)); ax[1].text(0, 0, "flux", ha="center", va="center")
xx_ = np.linspace(-2, 2, 100)
ax[1].plot(xx_, 1.1*np.cos(xx_*np.pi/4), color=BLUE); ax[1].plot(xx_, -1.1*np.cos(xx_*np.pi/4), color=GREEN)
ax[1].plot([-2, 2], [0, 0], "ko"); ax[1].text(-2, .2, "start"); ax[1].text(1.8, .2, "end")
ax[1].text(0, 1.35, "path above", ha="center", color=BLUE); ax[1].text(0, -1.5, "path below", ha="center", color=GREEN)
ax[1].text(0, -2.1, "phase difference = e^(i q flux / hbar): a closed loop sees its enclosed flux", ha="center", fontsize=9)
ax[1].set_xlim(-2.5, 2.5); ax[1].set_ylim(-2.4, 1.8); ax[1].set_aspect("equal"); ax[1].axis("off")
ax[1].set_title("path-integral version: Aharonov-Bohm, two paths around a flux")
save(fig, "s3-wilson-loop-and-aharonov-bohm")

fig, ax = plt.subplots(1, 2, figsize=(11, 4.6))
xx_ = np.linspace(-.5, .5, 100); ax[0].fill_between(xx_, np.sqrt(1 - xx_**2), 2.6, color=BLUE, alpha=.15)
ax[0].plot(xx_, np.sqrt(1 - xx_**2), color=BLUE); ax[0].plot([-.5, -.5], [np.sqrt(.75), 2.6], color=BLUE); ax[0].plot([.5, .5], [np.sqrt(.75), 2.6], color=BLUE)
fp = [(1 + np.sqrt(5))/2, (1 - np.sqrt(5))/2]; cen, rad = sum(fp)/2, abs(fp[0] - fp[1])/2
tha = np.linspace(0, np.pi, 200); ax[0].plot(cen + rad*np.cos(tha), rad*np.sin(tha), color=ORANGE)
z0 = cen + rad*np.exp(1.1j); gz = (2*z0 + 1)/(z0 + 1); ta = np.linspace(min(np.angle(z0 - cen), np.angle(gz - cen)), max(np.angle(z0 - cen), np.angle(gz - cen)), 50)
ax[0].plot(cen + rad*np.cos(ta), rad*np.sin(ta), color=RED, lw=3)
ax[0].text(cen, rad + .12, "axis of [[2,1],[1,1]]: closed geodesic\nof length 2 arccosh(3/2) = 1.925 (S5)", ha="center", fontsize=8, color=RED)
ax[0].text(0, 2.2, "fundamental domain\nof SL(2,Z)", ha="center", fontsize=8, color=BLUE)
ax[0].set_xlim(-1.6, 2.2); ax[0].set_ylim(0, 2.6); ax[0].set_aspect("equal"); ax[0].set_xlabel("Re z"); ax[0].set_ylabel("Im z")
ax[0].set_title("elsewhere: Selberg trace formula on the half-plane\n(spectrum = sum over closed geodesics, exact)")
for (aa, bb, c) in [(1, .6, BLUE), (1.6, 1.0, GREEN)]: ax[1].plot(aa*np.cos(th), bb*np.sin(th), color=c)
t8 = np.linspace(0, 2*np.pi, 300); ax[1].plot(2.2*np.sin(t8), 1.1*np.sin(2*t8), color=ORANGE)
ax[1].text(0, -1.45, "each closed orbit adds  A e^(i S/hbar)", ha="center")
ax[1].set_aspect("equal"); ax[1].axis("off")
ax[1].set_title("path-integral version: Gutzwiller trace formula\n(sum over periodic orbits, leading order only)")
save(fig, "s4-selberg-and-gutzwiller")

fig, ax = plt.subplots(1, 2, figsize=(11, 4.4))
r = np.linspace(-3, 3, 400)
ax[0].plot(r, np.cos(r**2), color=RED, label="along the real line: Re e^(i x^2) keeps oscillating")
ax[0].plot(r, np.exp(-r**2), color=GREEN, label="along the 45-degree line: e^(-s^2) decays")
ax[0].set_xlabel("distance along the line"); ax[0].legend(fontsize=8, loc="lower left")
ax[0].set_title("elsewhere: Lefschetz thimble (turn the contour)")
R = np.linspace(0, 6, 800); S_, C_ = fresnel(R*np.sqrt(2/np.pi)); part = 2*np.sqrt(np.pi/2)*(C_ + 1j*S_)
gpart = np.exp(1j*np.pi/4)*np.sqrt(np.pi)*np.array([__import__("math").erf(x) for x in R])
ax[1].plot(part.real, part.imag, color=RED, label="real-line partial sums (spiral, slow)")
ax[1].plot(gpart.real, gpart.imag, color=GREEN, lw=2, label="turned-line partial sums (straight, fast)")
ax[1].plot(np.sqrt(np.pi/2), np.sqrt(np.pi/2), "k*", ms=12); ax[1].text(1.35, .95, "sqrt(pi) e^(i pi/4)")
ax[1].set_aspect("equal"); ax[1].legend(fontsize=8, loc="upper left"); ax[1].set_xlabel("Re"); ax[1].set_ylabel("Im")
ax[1].set_title("path-integral version: Wick rotation = a dial turn (S6)")
save(fig, "s5-thimble-is-a-dial-turn")

fig, ax = plt.subplots(1, 2, figsize=(11, 4.4))
ax[0].add_patch(Polygon([(0, 0), (2, 0), (2, 1.5), (0, 1.5)], color=ORANGE, alpha=.25))
arrow(ax[0], (0, 0), (2, 0), BLUE); arrow(ax[0], (2, 0), (2, 1.5), BLUE)
arrow(ax[0], (0, 0), (0, 1.5), GREEN); arrow(ax[0], (0, 1.5), (2, 1.5), GREEN)
ax[0].text(1, -.25, "A then B", ha="center", color=BLUE); ax[0].text(-.1, .75, "B then A", ha="right", color=GREEN)
ax[0].text(1, .75, "area between = Levy area\nlog of the product picks up (1/2)[A,B]", ha="center", fontsize=9)
ax[0].set_xlim(-1.2, 3); ax[0].set_ylim(-.5, 2); ax[0].axis("off")
ax[0].set_title("elsewhere: signature / Magnus (log of an ordered product)")
Am = .05*rng.normal(size=(3, 3)); Bm = .05*rng.normal(size=(3, 3)); from scipy.linalg import logm
cm = lambda X, Y: X @ Y - Y @ X; Lg = logm(expm(Am) @ expm(Bm))
e1 = abs(Lg - (Am + Bm)).max(); e2 = abs(Lg - (Am + Bm + cm(Am, Bm)/2)).max()
e3 = abs(Lg - (Am + Bm + cm(Am, Bm)/2 + (cm(Am, cm(Am, Bm)) - cm(Bm, cm(Am, Bm)))/12)).max()
ax[1].bar(["A + B", "+ [A,B]/2", "+ nested /12"], [e1, e2, e3], color=[GREY, BLUE, GREEN]); ax[1].set_yscale("log")
ax[1].set_ylabel("error in log(e^A e^B)")
ax[1].set_title("path-integral version: Trotter error led by the twist (S7)\nbookkeeping only: no H, no hbar, no measure")
save(fig, "s6-signature-and-magnus")

fig = plt.figure(figsize=(11, 4.4))
ax3 = fig.add_subplot(1, 2, 1, projection="3d")
U, V = np.meshgrid(np.linspace(0, 2*np.pi, 60), np.linspace(0, np.pi, 30))
X, Y, Z = np.sin(V)*np.cos(U), np.sin(V)*np.sin(U), np.cos(V)
ax3.plot_surface(X, Y, Z, facecolors=plt.cm.coolwarm((Z + 1)/2), alpha=.6, linewidth=0)
ax3.scatter([0, 0], [0, 0], [1, -1], color="k", s=40); ax3.text(0, 0, 1.25, "north pole"); ax3.text(0, 0, -1.4, "south pole")
ax3.set_box_aspect((1, 1, 1)); ax3.axis("off"); ax3.set_title("elsewhere: Duistermaat-Heckman\nheight z generates a circle action; only the poles count")
ax = fig.add_subplot(1, 2, 2); tt = np.linspace(.2, 15, 300)
ex = [2*np.pi*(quad(lambda z: np.cos(t*z), -1, 1)[0]) for t in tt]
st = (2*np.pi*(np.exp(1j*tt) - np.exp(-1j*tt))/(1j*tt)).real
ax.plot(tt, ex, color=BLUE, lw=4, alpha=.4, label="exact integral over the sphere"); ax.plot(tt, st, "--", color=RED, label="stationary phase at the two poles")
ax.set_xlabel("t  (1/hbar)"); ax.legend(fontsize=8); ax.set_title("path-integral version: semiclassical = exact for circle actions (S1)")
save(fig, "s7-duistermaat-heckman")

# ---------------- SCORECARD ----------------
fig, ax = plt.subplots(figsize=(15, 7.5)); ax.axis("off"); ax.set_xlim(0, 1); ax.set_ylim(0, 1)
rows = [("Dynamics", "Bargmann twist = mass; quadratic closure (Jacobi group)", "established", "a fixed central twist size; generators of degree <= 2 for exactness", "free H picked out; V(x) is data", GREEN),
        ("hbar", "Bohr-Sommerfeld+Maslov, Kostant-Souriau, Dirac", "established", "closed loops' area in units of 2 pi hbar must be whole", "fixes area/hbar, never hbar", ORANGE),
        ("hbar (ours)", "curvature trade on Klauder's half-plane", "ours, arithmetic", "compact quotient needs SL(2,R), beyond the dial", "fixes beta, hbar cancels", ORANGE),
        ("Measure", "Feynman-Kac (real k); Cameron no-go at k = -i/hbar", "established", "real Brownian measure only at real k", "turning the dial loses the measure", ORANGE),
        ("Measure", "Daubechies-Klauder-Paul: Brownian motion on the half-plane", "established (rigorous)", "the dial's left-invariant metric, one scale beta hbar", "best fit: uses the dial's own group", GREEN)]
hdr = ["gap", "published piece", "status", "what the framework must add", "what it buys"]
xcol = [0, .085, .40, .51, .80]
for xc, hname in zip(xcol, hdr): ax.text(xc, .97, hname, weight="bold")
for i, (g_, piece, st_, need, buys, c) in enumerate(rows):
    y = .89 - i*.075
    for xc, txt in zip(xcol, [g_, piece, st_, need, buys]): ax.text(xc, y, txt, fontsize=8.5, color=c if xc == xcol[4] else "k")
ax.plot([0, 1], [.5, .5], color=GREY, lw=.6)
ax.text(0, .46, "shape matches", weight="bold")
shapes = [("checkerboard / telegraph", "the chain with a 2-state dial, real vs imaginary turn weight", "strong lead", GREEN),
          ("Selberg on the half-plane", "exact sum over closed geodesics on the dial's own space", "strong lead", GREEN),
          ("Duistermaat-Heckman", "exactness of stationary phase for circle actions", "lead", GREEN),
          ("Lefschetz thimbles", "contour turn = dial turn; converts oscillation to decay", "lead (finite dims)", ORANGE),
          ("Wilson loop / Aharonov-Bohm", "holonomy; already used in gap 2", "covered", GREY),
          ("Suzuki-Trotter", "quantum chain = classical chain one dimension up", "tool, not a lead", GREY),
          ("signature / Magnus / rough paths", "bookkeeping of ordered products", "none of the three gaps", GREY)]
for i, (nm, what, verdict, c) in enumerate(shapes):
    y = .39 - i*.055
    ax.text(0, y, nm, fontsize=8.5); ax.text(.28, y, what, fontsize=8.5); ax.text(.83, y, verdict, fontsize=8.5, color=c, weight="bold")
ax.set_title("scorecard: what each published piece supplies, and which shape matches are worth following", fontsize=12)
save(fig, "z-scorecard")
print("figures written:", sorted(os.listdir(OUT)))
