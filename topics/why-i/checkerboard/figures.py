"""Pictures for topics/why-i/checkerboard. Run: python3 figures.py"""
import numpy as np, matplotlib, itertools, os, sys
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch
from matplotlib.colors import LogNorm
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from cb import evolve, kernel
from scipy.special import iv
OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "figures"); os.makedirs(OUT, exist_ok=True)
BLUE, ORANGE, GREEN, RED, GREY, PURPLE = "#2b6cb0", "#dd6b20", "#2f855a", "#c53030", "#718096", "#6b46c1"
plt.rcParams.update({"font.size": 10, "axes.titlesize": 11, "figure.dpi": 130})
def save(fig, name):
    fig.tight_layout(); fig.savefig(os.path.join(OUT, name + ".png"), bbox_inches="tight"); plt.close(fig)
def path_from(S_):
    d, x, xs = 0, 0, [0]
    for turn in S_:
        if turn: d = 1 - d
        x += 1 if d == 0 else -1; xs.append(x)
    return np.array(xs), d

# 1. one path = one set of turn times = one ordered sub-chain
fig, ax = plt.subplots(1, 2, figsize=(11, 4.8), gridspec_kw={"width_ratios": [1.2, 1]})
S_ = [0, 0, 1, 0, 0, 0, 1, 0, 1, 0]
xs, _ = path_from(S_); ts = np.arange(len(xs))
for xg in range(-4, 6):
    for tg in range(11): ax[0].plot(xg, tg, ".", color="#cbd5e0", ms=4)
ax[0].plot(xs, ts, "-o", color=BLUE, ms=4)
for j, turn in enumerate(S_):
    if turn: ax[0].plot(xs[j], ts[j], "s", ms=13, mfc="none", color=ORANGE, mew=2)
ax[0].text(3.0, 4.5, "each orange square is a turn:\nk = -1: weight a eps (Kac)\nk = -i: weight i m eps (Feynman)", color=ORANGE, fontsize=9)
ax[0].set_xlabel("position (steps of eps)"); ax[0].set_ylabel("time step"); ax[0].set_xlim(-4.5, 7.5)
ax[0].set_title("a checkerboard path: move one step each tick, sometimes turn")
ax[1].axis("off"); ax[1].set_xlim(0, 10); ax[1].set_ylim(0, 10)
for j in range(10):
    c = ORANGE if S_[j] else "#e2e8f0"
    ax[1].add_patch(plt.Rectangle((j, 7.4), .9, .9, color=c)); ax[1].text(j + .45, 7.0, f"t{j+1}", ha="center", fontsize=8)
ax[1].text(5, 8.8, "the same path, as a choice of turn times", ha="center")
ax[1].text(5, 5.6, "chain:  t10 o_k t9 o_k ... o_k t1", ha="center", fontsize=11)
ax[1].text(5, 4.5, "expands into every ordered sub-chain once;\nthe sub-chain {t3, t7, t9} IS this path", ha="center", fontsize=9)
ax[1].text(5, 2.4, "weight of the path = (-k eps mu)^3\nk = -1, mu = a:  (a eps)^3\nk = -i, mu = m:  (i m eps)^3 = -i (m eps)^3", ha="center", fontsize=9, color=PURPLE)
ax[1].set_title("one path = one ordered sub-chain (check B1, B2)")
save(fig, "c1-a-path-is-a-sub-chain")

# 2. all 16 paths for N = 4, grouped by endpoint, with their weights at both dials
fig, axs = plt.subplots(2, 8, figsize=(14, 5))
allS = sorted(itertools.product([0, 1], repeat=4), key=lambda s: (path_from(s)[0][-1], sum(s)))
for axx, s in zip(axs.flat, allS):
    xs, d = path_from(s)
    axx.plot(xs, range(5), "-o", ms=3, color=BLUE); R = sum(s)
    for j, turn in enumerate(s):
        if turn: axx.plot(xs[j], j, "s", mfc="none", color=ORANGE, ms=7)
    axx.set_xlim(-4.5, 4.5); axx.set_ylim(-.5, 4.5); axx.set_xticks([]); axx.set_yticks([])
    wk = {0: "1", 1: "i", 2: "-1", 3: "-i", 4: "1"}[R]
    axx.set_title(f"end x={xs[-1]}, R={R}\nk=-1: (a e)^{R}  k=-i: {wk}(m e)^{R}", fontsize=7.5)
fig.suptitle("all 16 four-step paths (start moving right). The chain's 16 ordered sub-chains are exactly these paths.", y=1.02)
save(fig, "c2-all-paths-for-four-steps")

# 3. adding vs cancelling: weights at an endpoint, both dials
fig, ax = plt.subplots(1, 2, figsize=(11, 4.4))
N, eps, mu, xend = 12, 0.15, 1.3, 2
ws = {}
for s in itertools.product([0, 1], repeat=N):
    xs, d = path_from(s)
    if xs[-1] == xend and d == 1: ws.setdefault(sum(s), 0); ws[sum(s)] += 1
Rs = sorted(ws)
ax[0].bar(Rs, [ws[R]*(eps*mu)**R for R in Rs], color=GREEN)
ax[0].set_xlabel("number of turns R"); ax[0].set_ylabel("total weight N(R) (a eps)^R")
ax[0].set_title(f"k = -1: every path adds a positive weight\n(paths ending at x = {xend}, moving left, N = {N})")
vals = np.array([ws[R]*((1j*eps*mu)**R).imag for R in Rs])
ax[1].bar(Rs, vals, color=[BLUE if v > 0 else RED for v in vals])
ax[1].axhline(vals.sum(), color="k", ls="--", label=f"sum = {vals.sum():.3f} i  (vs {sum(ws[R]*(eps*mu)**R for R in Rs):.3f} at k = -1)")
ax[1].axhline(0, color=GREY, lw=.6)
ax[1].set_xlabel("number of turns R"); ax[1].set_ylabel("total weight N(R) (i m eps)^R, imaginary part"); ax[1].legend(fontsize=8)
ax[1].set_title("k = -i: the same groups, each turn multiplies by i,\nso the signs alternate and the groups cancel")
save(fig, "c3-adding-versus-cancelling")

# 4. space-time pictures at both dials
fig, ax = plt.subplots(1, 2, figsize=(11, 4.6))
N, eps, mu = 200, 0.02, 1.3
P = evolve(N, eps, mu, -1, normalized=True); A = evolve(N, eps, mu, -1j, normalized=True)
c0 = P.shape[1]//2; ext = [-(c0)*eps, (P.shape[1] - 1 - c0)*eps, 0, N*eps]
for i in range(N + 1):                                     # remove the never-turned atom on the right light-cone edge
    P[i, c0 + i, 0] = 0; A[i, c0 + i, 0] = 0
Pd = P.real.sum(axis=2); Ad = (np.abs(A)**2).sum(axis=2)
for i in range(N + 1):
    for arr in (Pd, Ad):                                   # merge the empty parity sites with their neighbour for display
        row = arr[i]; arr[i] = row + np.roll(row, 1)
ax[0].imshow(Pd + 1e-12, norm=LogNorm(vmin=Pd.max()*1e-3, vmax=Pd.max()), origin="lower", extent=ext, aspect="auto", cmap="Greens"); ax[0].set_title("k = -1: Kac's telegraph, probability\n(log colour, never-turned atom removed; smooth spread)")
ax[1].imshow(Ad + 1e-12, norm=LogNorm(vmin=Ad.max()*1e-3, vmax=Ad.max()), origin="lower", extent=ext, aspect="auto", cmap="Blues"); ax[1].set_title("k = -i: Feynman's checkerboard, |amplitude|^2\n(log colour, never-turned atom removed; dark bands = cancellation)")
for a_ in ax: a_.set_xlim(-4.2, 4.2); a_.set_xlabel("position x"); a_.set_ylabel("time t"); a_.plot([0, 4], [0, 4], "--", color=GREY, lw=.8); a_.plot([0, -4], [0, 4], "--", color=GREY, lw=.8)
save(fig, "c4-space-time-at-both-dials")

# 5. the continuum limit at k = -1
T, a = 2.0, 1.3
xx = np.linspace(-1.98, 1.98, 400)
GR, GL = kernel(xx, T, a)
fig, ax = plt.subplots(1, 2, figsize=(11, 4.4))
for N, c in [(20, "#9ae6b4"), (40, "#48bb78"), (160, "#22543d")]:
    e_ = T/N; ev = evolve(N, e_, a, -1, normalized=True)[-1]; c0 = ev.shape[0]//2
    idx = np.arange(-N + 2, N - 1, 2); x_ = idx*e_   # reachable sites have the parity of N
    ax[0].plot(x_, ev[c0 + idx, 0].real/(2*e_), "o", ms=3, color=c, label=f"chain, eps = {e_:.3f}")
    ax[1].plot(x_, ev[c0 + idx, 1].real/(2*e_), "o", ms=3, color=c)
ax[0].plot(xx, (np.exp(-a*T)*GR).real, color="k", label="Kac's kernel"); ax[1].plot(xx, (np.exp(-a*T)*GL).real, color="k")
ax[0].text(-1.9, .12, "plus an atom of probability e^(-aT) = 0.074\nat x = T: the path that never turned", fontsize=8)
ax[0].set_title("k = -1, moving right at time T = 2"); ax[1].set_title("k = -1, moving left at time T = 2")
for a_ in ax: a_.set_xlabel("position x"); a_.set_ylabel("probability per unit length")
ax[0].legend(fontsize=8); save(fig, "c5-continuum-limit-telegraph")

# 6. the continuum limit at k = -i
m_ = 1.3; GR, GL = kernel(xx, T, 1j*m_)
fig, ax = plt.subplots(1, 2, figsize=(11, 4.4))
for N, c in [(20, "#90cdf4"), (40, "#4299e1"), (160, "#1a365d")]:
    e_ = T/N; ev = evolve(N, e_, m_, -1j)[-1]; c0 = ev.shape[0]//2
    idx = np.arange(-N + 2, N - 1, 2); x_ = idx*e_   # reachable sites have the parity of N
    ax[0].plot(x_, ev[c0 + idx, 0].real/(2*e_), "o", ms=3, color=c, label=f"chain, eps = {e_:.3f}")
    ax[1].plot(x_, ev[c0 + idx, 1].imag/(2*e_), "o", ms=3, color=c)
ax[0].plot(xx, GR.real, color="k", label="Dirac kernel, -(m/2) sqrt((t+x)/(t-x)) J1"); ax[1].plot(xx, GL.imag, color="k", label="Dirac kernel, (m/2) J0")
ax[0].set_ylim(-3, 1.5)
ax[0].set_title("k = -i, moving right (real part)"); ax[1].set_title("k = -i, moving left (imaginary part)")
for a_ in ax: a_.set_xlabel("position x"); a_.set_ylabel("amplitude per unit length"); a_.legend(fontsize=8)
save(fig, "c6-continuum-limit-checkerboard")

# 7. one chain, one dial: what changes and what stays
fig, ax = plt.subplots(1, 2, figsize=(11, 4.6))
ax[0].axis("off"); ax[0].set_xlim(0, 10); ax[0].set_ylim(0, 10)
ax[0].text(5, 9.4, "the same chain of turn elements t = -eps mu X", ha="center", weight="bold")
ax[0].text(2.4, 7.9, "k = -1", ha="center", color=GREEN, fontsize=12); ax[0].text(7.6, 7.9, "k = -i", ha="center", color=BLUE, fontsize=12)
L = [("turn weight", "a eps", "i m eps"), ("normalized step", "stochastic (1 - a e, a e)", "unitary (cos m e, i sin m e)"),
     ("smooth limit", "telegraph equation", "1+1 Dirac equation"), ("kernel", "I0, I1 (grow)", "J0, J1 (oscillate)"),
     ("long times", "diffusion, D = 1/(2a)", "Schrodinger, 1/(2m)"), ("over paths", "a probability measure", "complex weights, no probability")]
for i, (lab, l, r) in enumerate(L):
    y = 6.9 - i*1.1
    ax[0].text(5, y, lab, ha="center", color=GREY, fontsize=8.5); ax[0].text(2.4, y - .45, l, ha="center", color=GREEN, fontsize=9); ax[0].text(7.6, y - .45, r, ha="center", color=BLUE, fontsize=9)
ax[0].set_title("Gaveau-Jacobson-Kac-Schulman (1984): one model, continued in the turn rate")
p = np.linspace(-6, 6, 400); a = m_ = 2.0
lam = (1 - 0*p)
rate_slow = a - np.sqrt(np.maximum(a**2 - p**2, 0) + 0j).real
ax[1].plot(p, rate_slow, color=GREEN, label="k = -1: decay rate of the slow mode")
ax[1].plot(p, p**2/(2*a), "--", color=GREEN, lw=.8, label="p^2/(2a): diffusion")
ax[1].plot(p, np.sqrt(m_**2 + p**2) - m_, color=BLUE, label="k = -i: frequency minus rest energy")
ax[1].plot(p, p**2/(2*m_), ":", color=BLUE, lw=1.2, label="p^2/(2m): Schrodinger")
ax[1].set_ylim(0, 4); ax[1].set_xlabel("momentum p"); ax[1].legend(fontsize=8)
ax[1].set_title("long-time limits share the parabola p^2/(2 mu) (checks D1, D2)")
save(fig, "c7-one-chain-one-dial")

# 8. the dial is on the turn weight, not on time
fig, ax = plt.subplots(figsize=(9, 5))
ax.axhline(0, color="k", lw=.8); ax.axvline(0, color="k", lw=.8)
ax.plot(-1, 0, "o", color=GREEN, ms=10); ax.text(-1, .18, "k = -1\ntelegraph", ha="center", color=GREEN)
ax.plot(0, -1, "o", color=BLUE, ms=10); ax.text(.08, -1.25, "k = -i: checkerboard", color=BLUE)
th = np.linspace(np.pi, 1.5*np.pi, 50); ax.plot(np.cos(th), np.sin(th), "--", color=PURPLE)
ax.text(-1.55, -1.05, "turn weight a eps -> i m eps\n(continue in the mass:\nGJKS 1984, Nelson 1964)", color=PURPLE, fontsize=9)
ax.text(.25, .55, "what the dial multiplies here: only the turn element.\nThe moves are fixed lattice shifts at both settings.\nSo this is NOT a Wick rotation of time: rotating t in the\nDirac equation does not give the telegraph process.", fontsize=9)
ax.set_xlim(-1.8, 2.2); ax.set_ylim(-1.6, 1.2); ax.set_xlabel("Re k"); ax.set_ylabel("Im k"); ax.set_aspect("equal")
ax.set_title("where the two models sit on the k-dial")
save(fig, "c8-the-dial-is-on-the-turns")
print(sorted(os.listdir(OUT)))
