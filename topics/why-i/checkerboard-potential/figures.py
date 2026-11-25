"""Pictures for topics/why-i/checkerboard-potential (Nov 25, 2026). Run: python3 figures.py"""
import numpy as np, matplotlib, os, sys, itertools
matplotlib.use("Agg")
import matplotlib.pyplot as plt
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from cbv import evolve_v, reference, schrodinger, heat_killing
OUT = os.path.join(HERE, "figures"); os.makedirs(OUT, exist_ok=True)
BLUE, ORANGE, GREEN, RED, GREY, PURPLE = "#2b6cb0", "#dd6b20", "#2f855a", "#c53030", "#718096", "#6b46c1"
plt.rcParams.update({"font.size": 10, "axes.titlesize": 11, "figure.dpi": 130})
def save(fig, name):
    fig.tight_layout(); fig.savefig(os.path.join(OUT, name + ".png"), bbox_inches="tight"); plt.close(fig)
def path_from(S_):
    d, x, xs = 0, 0, [0]
    for turn in S_:
        if turn: d = 1 - d
        x += 1 if d == 0 else -1; xs.append(x)
    return np.array(xs)

# p1. a path over a potential: one tick = turn element + potential element
eps = 0.25; S_ = [0, 0, 1, 0, 0, 0, 1, 0, 1, 0, 0, 1]
xs = path_from(S_); ts = np.arange(len(xs))
Vf = lambda x: 2.0*np.exp(-((x - 0.5)/0.6)**2)
fig, ax = plt.subplots(1, 2, figsize=(11.5, 5.2), gridspec_kw={"width_ratios": [1.15, 1]})
xx = np.linspace(-5, 6, 300)
ax[0].fill_between(xx, -1.6, -1.6 + 0.9*Vf(xx*eps), color="#fbd38d", alpha=.7)
ax[0].text(2.3, -1.35, "potential V(x)", color="#975a16", fontsize=9)
for xg in range(-4, 6):
    for tg in range(13): ax[0].plot(xg, tg, ".", color="#cbd5e0", ms=4)
ax[0].plot(xs, ts, "-", color=BLUE, lw=1.5)
sc = ax[0].scatter(xs[:-1], ts[:-1], c=Vf(xs[:-1]*eps), cmap="YlOrBr", s=46, edgecolor=GREY, zorder=3, vmin=0, vmax=2)
for j, turn in enumerate(S_):
    if turn: ax[0].plot(xs[j], ts[j], "s", ms=14, mfc="none", color=PURPLE, mew=2)
ax[0].set_xlim(-4.5, 6); ax[0].set_ylim(-1.7, 12.6)
ax[0].set_xlabel("position (steps of eps)"); ax[0].set_ylabel("tick")
ax[0].set_title("one path: purple squares are turns, dot colour is V where it stands")
fig.colorbar(sc, ax=ax[0], fraction=.04, label="V at that tick")
a1 = ax[1]; a1.axis("off"); a1.set_xlim(0, 10); a1.set_ylim(0, 10)
a1.text(5, 9.4, "each tick now holds two elements", ha="center", fontsize=11)
a1.text(5, 8.5, "turn  t = -eps m X        potential  v = eps V(x)", ha="center", color=PURPLE)
a1.text(5, 7.6, "they commute at a site, so  (1 + k v)(1 + k t) = 1 + k (v o_k t)", ha="center", fontsize=9)
a1.text(5, 6.4, "k = -i (quantum):   turn  i m eps,     potential  1 - i eps V  ->  e^(-i eps V)", ha="center", color=BLUE, fontsize=9)
a1.text(5, 5.7, "k = -1 (random):     turn  a eps,         potential  1 - eps V  ->  e^(-eps V)", ha="center", color=GREEN, fontsize=9)
a1.text(5, 4.4, "[ours]  N ticks = free moves x (1 + k (chain of all 2N elements))", ha="center", fontsize=10, color=RED)
a1.text(5, 3.3, "a sub-chain = (which turns, which potential kicks)\ngroup by the turns:  every path still counted once,\n"
        "with weight  (turns) x product of (1 + k eps V) along it", ha="center", fontsize=9)
a1.text(5, 1.3, "Feynman-Hibbs problem 2-6; Jacobson-Schulman 1984; Cedzich et al. 2013", ha="center", fontsize=8, color=GREY)
save(fig, "p1-potential-as-a-chain-element")

# p2. what V does to one path's weight at the two dials
N = 200; e_ = 0.01; m = 1.3
rng = np.random.default_rng(7)
d, x = 0, 0; xpath = []; turns = []
for j in range(N):
    turn = rng.random() < 3*m*e_
    if turn: d = 1 - d
    turns.append(turn); xpath.append(x); x += 1 if d == 0 else -1
xpath = np.array(xpath)*e_; Vp = 3*np.cos(3*xpath) + 2 + 4*xpath
cum = np.cumsum(e_*Vp)
fig, ax = plt.subplots(1, 3, figsize=(13, 4.2))
ax[0].plot(xpath, np.arange(N)*e_, color=BLUE); ax[0].set_xlabel("x"); ax[0].set_ylabel("time")
ax[0].set_title("one sample path (time up)")
th = np.linspace(0, 2*np.pi, 200); ax[1].plot(np.cos(th), np.sin(th), color="#e2e8f0")
z = np.exp(-1j*cum)
ax[1].scatter(z.real, z.imag, c=np.arange(N), cmap="Blues", s=10)
ax[1].set_aspect("equal"); ax[1].set_xlim(-1.3, 1.3); ax[1].set_ylim(-1.3, 1.3)
ax[1].set_title("k = -i: V only turns the weight's phase\n(size stays 1 on the circle)", color=BLUE)
ax[1].set_xlabel("real part"); ax[1].set_ylabel("imaginary part")
ax[2].plot(np.arange(N)*e_, np.exp(-cum), color=GREEN, lw=2, label="k = -1: e^(-sum eps V)  (killing)")
ax[2].plot(np.arange(N)*e_, np.abs(z), color=BLUE, lw=2, ls="--", label="k = -i: |e^(-i sum eps V)| = 1")
ax[2].set_xlabel("time"); ax[2].set_ylabel("size of the potential factor"); ax[2].legend(fontsize=8)
ax[2].set_title("same path, same V: shrinks at k = -1, keeps size at k = -i")
save(fig, "p2-one-path-two-dials")

# p3. smooth limit at both dials: chain vs an independent continuum solver
L = 8.0
Vfun = lambda x: 3.0*np.exp(-x**2) + 0.5*x
sig, p0, T = 0.4, 2.0, 1.5
fig, ax = plt.subplots(1, 3, figsize=(13.5, 4.2))
cols = {-1.0: GREEN, -1j: BLUE}
for k, lab, norm_tel in [(-1.0, "k = -1: telegraph with killing", True), (-1j, "k = -i: Dirac with V", False)]:
    errs, es = [], [0.04, 0.02, 0.01]
    for e_ in es:
        nn = int(round(L/e_)); x = (np.arange(nn) - nn//2)*e_
        phi = np.exp(-x**2/(2*sig**2))*(np.exp(1j*p0*x) if k == -1j else 1)
        psi0 = np.stack([phi, 0.4*phi], axis=1); Nt = int(round(T/e_))
        ch = evolve_v(psi0, Nt, e_, 1.3, k, Vfun(x), linear=False, normalized=True)
        ref = reference(psi0, Nt*e_, e_, 1.3, k, Vfun(x), normalized_telegraph=norm_tel, steps=4*Nt)
        errs.append(np.abs(ch - ref).max()/np.abs(ref).max())
    ax[0].loglog(es, errs, "o-", color=cols[k], label=lab)
    j = 1 if k == -1.0 else 2
    dens = lambda P: (np.abs(P)**2).sum(axis=1) if k == -1j else P.real.sum(axis=1)
    ax[j].plot(x, dens(ref), color="k", lw=2.5, alpha=.35, label="continuum solver")
    ax[j].plot(x, dens(ch), color=cols[k], lw=1.2, label="chain, eps = 0.01")
    ax[j].plot(x, Vfun(x)/Vfun(x).max()*dens(ref).max()*0.6, color="#d69e2e", ls=":", label="V (scaled)")
    ax[j].set_xlim(-4, 4); ax[j].set_xlabel("x"); ax[j].legend(fontsize=8)
ax[0].loglog([0.04, 0.01], [0.05, 0.0125], color=GREY, ls="--", label="slope 1 (error halves)")
ax[0].set_xlabel("eps"); ax[0].set_ylabel("max relative error at t = 1.5"); ax[0].legend(fontsize=8)
ax[0].set_title("the chain converges at both dials")
ax[1].set_title("k = -1: density at t = 1.5"); ax[2].set_title("k = -i: |psi|^2 at t = 1.5")
save(fig, "p3-smooth-limit-both-dials")

# p4. long times at k = -i: Schrodinger with V, and the relativistic edge
def centres(m, e_, x0, every, w=1.0, L2=6.0):
    nn = int(round(L2/e_)); x = (np.arange(nn) - nn//2)*e_
    Vh = 0.5*m*w**2*x**2; s0 = 1/np.sqrt(m*w)
    chi0 = np.exp(-(x - x0)**2/(2*s0**2)); chi0 /= np.sqrt((np.abs(chi0)**2).sum()*e_)
    psi = np.stack([chi0/np.sqrt(2), -chi0/np.sqrt(2)], axis=1).astype(complex)
    Nq = int(round(2*np.pi/w/e_)); tt, cc = [0.0], [x0]; done = 0
    while done < Nq:
        st = min(every, Nq - done); psi = evolve_v(psi, st, e_, m, -1j, Vh, linear=False, normalized=True); done += st
        dens = (np.abs(psi)**2).sum(axis=1); tt.append(done*e_); cc.append((x*dens).sum()/dens.sum())
    return np.array(tt), np.array(cc)
fig, ax = plt.subplots(1, 2, figsize=(11.5, 4.2))
tg = np.linspace(0, 2*np.pi, 300)
for x0, col, lab in [(0.1, BLUE, "slow swing, top speed 0.1"), (0.5, RED, "fast swing, top speed 0.5")]:
    tt, cc = centres(40.0, 0.00125, x0, 80)
    ax[0].plot(tt, cc/x0, "o", ms=3, color=col, label=f"chain, {lab}")
ax[0].plot(tg, np.cos(tg), color="k", alpha=.4, lw=2, label="Schrodinger with V: cos t")
ax[0].set_xlabel("time"); ax[0].set_ylabel("packet centre / starting offset"); ax[0].legend(fontsize=8)
ax[0].set_title("k = -i, m = 40, V = (m/2) x^2 (light speed = 1)")
ms = [10, 40]; slow = [0.9973, 0.9994]; fast = [0.9571, 0.9265]   # overlaps from checks.py D1 (fast at m = 10 from the same run with x0 = 0.5)
ax[1].bar(np.arange(2) - .18, slow, .34, color=BLUE, label="slow swing (0.1)")
ax[1].bar(np.arange(2) + .18, fast, .34, color=RED, label="fast swing (0.5)")
ax[1].set_xticks(range(2)); ax[1].set_xticklabels([f"m = {m_}" for m_ in ms]); ax[1].set_ylim(0.9, 1.012)
ax[1].set_ylabel("overlap with Schrodinger after half a period"); ax[1].legend(fontsize=8, loc="center right")
ax[1].set_title("slow packets become Schrodinger with V;\nfast ones stay Dirac (relativistic), as they should")
for i_, (s_, f_) in enumerate(zip(slow, fast)):
    ax[1].text(i_ - .18, s_ + .002, f"{s_:.4f}", ha="center", fontsize=8); ax[1].text(i_ + .18, f_ + .002, f"{f_:.4f}", ha="center", fontsize=8)
save(fig, "p4-schrodinger-with-V")

# p5. long times at k = -1: Feynman-Kac (diffusion with killing) and the oscillator ladder
a = 20.0; e_ = 0.002; nn = int(round(6.0/e_)); x = (np.arange(nn) - nn//2)*e_
c = 4.0; Vk = 0.5*c*x**2
u0 = np.exp(-x**2/(2*0.3**2)); u0 /= u0.sum()*e_
P0 = np.stack([u0/2, u0/2], axis=1)*e_
fig, ax = plt.subplots(1, 2, figsize=(11.5, 4.2))
for Tk, col in [(0.25, "#9ae6b4"), (1.0, GREEN), (3.0, "#22543d")]:
    Nk = int(round(Tk/e_))
    Pt = evolve_v(P0, Nk, e_, a, -1.0, Vk, linear=False, normalized=True).real
    ax[0].plot(x, Pt.sum(axis=1)/e_, color=col, lw=1.4, label=f"chain, t = {Tk}")
    ax[0].plot(x, heat_killing(u0, Nk*e_, e_, 1/(2*a), Vk, steps=2000), color="k", ls=":", lw=1)
ax[0].plot([], [], color="k", ls=":", label="u_t = (1/2a) u'' - V u")
ax[0].set_xlim(-2, 2); ax[0].set_xlabel("x"); ax[0].set_ylabel("surviving density"); ax[0].legend(fontsize=8)
ax[0].set_title("k = -1, a = 20, V = 2 x^2: Feynman-Kac with killing")
lim = np.sqrt(c/a)*(np.arange(3) + .5)
lat01 = [0.1998, 0.6022, 1.007]; lat004 = [0.2142, 0.6464, 1.0819]   # rates from checks.py D3
for i_ in range(3):
    ax[1].hlines(lim[i_], 0.6, 1.4, color="k", lw=2)
    ax[1].hlines(lat004[i_], 1.6, 2.4, color=GREEN, lw=2)
    ax[1].hlines(lat01[i_], 2.6, 3.4, color="#9ae6b4", lw=2)
    ax[1].text(0.45, lim[i_], f"n = {i_}", ha="right", va="center", fontsize=8)
ax[1].set_xticks([1, 2, 3]); ax[1].set_xticklabels(["limit sqrt(c/a)(n + 1/2)", "chain, eps = 0.004", "chain, eps = 0.01"], fontsize=8)
ax[1].set_ylabel("decay rate"); ax[1].set_xlim(0, 3.6)
ax[1].set_title("slowest decay rates: the oscillator ladder in imaginary time")
save(fig, "p5-feynman-kac-with-killing")

# p6. the measure at k = -i survives: each path's size equals its telegraph weight with no V
Nn = 12; e_ = 0.1; m = 1.3; Vs = lambda x: 1.5*np.cos(2*x) + 0.7*x
sz_q, sz_c = [], []
for S_ in itertools.product([0, 1], repeat=Nn):
    xs_ = path_from(S_)[:-1]*e_; R = sum(S_)
    wq = (1j*m*e_)**R*np.prod(np.exp(-1j*e_*Vs(xs_))); wc = (m*e_)**R
    sz_q.append(abs(wq)); sz_c.append(wc)
fig, ax = plt.subplots(1, 2, figsize=(11, 4.2))
ax[0].loglog(sz_c, sz_q, "o", ms=4, color=BLUE, alpha=.5)
ax[0].loglog([1e-12, 1], [1e-12, 1], color=GREY, ls="--")
ax[0].set_xlabel("telegraph weight of the path, no V  (m eps)^R"); ax[0].set_ylabel("|quantum weight with V|")
ax[0].set_title(f"all 2^{Nn} paths sit on the diagonal: V never changes a path's size")
N = 300; e_ = 0.01; n = 2*N + 3; xg = (np.arange(n) - n//2)*e_; Vb = 3*np.cos(3*xg) + 2*xg
psi0 = np.zeros((n, 2)); psi0[n//2, 0] = 1
Aq = evolve_v(psi0, N, e_, m, -1j, Vb, linear=False, keep=True)      # quantum amplitudes with V
Ac = evolve_v(psi0, N, e_, m, -1.0, np.zeros(n), keep=True).real      # path sizes summed: telegraph, no V
tt = np.arange(N + 1)*e_
ax[1].plot(tt, Ac.sum(axis=(1, 2)), color=BLUE, lw=2.5, label="sum of path sizes with V (= without V)")
ax[1].plot(tt, np.exp(m*tt), color="k", ls=":", lw=1.5, label="e^(m t)")
ax[1].plot(tt, np.abs(Aq).sum(axis=(1, 2)), color=ORANGE, lw=1.5, label="sum of |amplitude| over sites (paths interfere)")
ax[1].set_xlabel("time"); ax[1].set_ylabel("sum of |path weights|"); ax[1].legend(fontsize=8)
ax[1].set_title("bounded total size = a genuine (complex) path measure\n(Ichinose 1982, Ichinose-Tamura 1984)")
save(fig, "p6-measure-survives")

# p7. a vector potential = a phase on each move, and in 1+1 it is pure gauge
fig, ax = plt.subplots(1, 2, figsize=(11, 4.0))
a0 = ax[0]; a0.axis("off"); a0.set_xlim(-0.5, 6.5); a0.set_ylim(-1.5, 2.5)
for s in range(7): a0.plot(s, 0, "o", color=GREY, ms=9)
for s in range(6):
    a0.annotate("", xy=(s + .85, 0.15), xytext=(s + .15, 0.15), arrowprops=dict(arrowstyle="->", color=BLUE))
    a0.text(s + .5, 0.45, f"e^(i th{s})", ha="center", fontsize=8, color=BLUE)
    a0.annotate("", xy=(s + .15, -0.15), xytext=(s + .85, -0.15), arrowprops=dict(arrowstyle="->", color=ORANGE))
    a0.text(s + .5, -0.6, f"e^(-i th{s})", ha="center", fontsize=8, color=ORANGE)
a0.text(3, 1.6, "move right across a link: phase e^(i th);  left: e^(-i th)\nth = q A dx  (Feynman-Hibbs problem 2-6)", ha="center", fontsize=9)
a0.text(3, -1.3, "if th(s) = chi(s+1) - chi(s): same as multiplying by e^(i chi)", ha="center", fontsize=9, color=RED)
import cbv
n = 203; e_ = 0.05; x = (np.arange(n) - n//2)*e_
chi = np.sin(1.7*x) + 0.3*x**2; theta = np.roll(chi, -1) - chi; theta[-1] = 0
phi = np.exp(-x**2/0.2); psi0 = np.stack([phi, 0.5*phi], axis=1).astype(complex)
withA = evolve_v(psi0*np.exp(1j*chi)[:, None], 20, e_, 1.3, -1j, Vfun(x), linear=False, normalized=True, theta=theta)
gauge = evolve_v(psi0, 20, e_, 1.3, -1j, Vfun(x), linear=False, normalized=True)*np.exp(1j*chi)[:, None]
ax[1].plot(x, withA[:, 0].real, color=BLUE, lw=2.5, alpha=.5, label="with link phases (vector potential)")
ax[1].plot(x, gauge[:, 0].real, color="k", lw=1, ls="--", label="no vector potential, then times e^(i chi)")
ax[1].set_xlim(-3, 3); ax[1].set_xlabel("x"); ax[1].set_ylabel("real part of the right-mover")
ax[1].legend(fontsize=8); ax[1].set_title(f"identical (max difference {np.abs(withA - gauge).max():.0e})")
save(fig, "p7-vector-potential-is-gauge")

# p8. scorecard: what survived adding a potential
rows = [("V enters as one more element of the o_k chain", "yes, exact at finite eps [ours: the 2N-element chain]"),
        ("every path still counted once", "yes: group the sub-chains by their turns"),
        ("smooth limit, k = -i", "Dirac with V; slow packets give Schrodinger with V"),
        ("smooth limit, k = -1", "telegraph with killing; fast turning gives Feynman-Kac"),
        ("genuine path measure at k = -i", "survives any real V: sizes unchanged, total e^(m t)"),
        ("vector potential", "a phase per move; pure gauge in 1+1"),
        ("not covered", "a limit proof for V growing without bound (checked in a box), 3+1 dims, 2+1 magnetic fields")]
fig, ax = plt.subplots(figsize=(11, 3.8)); ax.axis("off"); ax.set_xlim(0, 1); ax.set_ylim(0, 1)
for i_, (q, a_) in enumerate(rows):
    y = 0.92 - i_*0.13; col = GREY if q == "not covered" else GREEN
    ax.add_patch(plt.Rectangle((0.0, y - 0.05), 0.012, 0.1, color=col))
    ax.text(0.025, y, q, va="center", fontsize=10, weight="bold"); ax.text(0.42, y, a_, va="center", fontsize=10, color=col)
ax.set_title("adding a potential: what holds")
save(fig, "p8-scorecard")
print("figures written to", OUT)
