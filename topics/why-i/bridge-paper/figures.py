"""Figures w1-w7 for the widened paper (DRAFT). Run: python3 figures.py -> figures/*.png"""
import os, numpy as np, matplotlib
matplotlib.use('Agg'); import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, Ellipse, Rectangle
D = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'figures'); os.makedirs(D, exist_ok=True)
C = dict(proved='#2e8b57', observed='#1f6fb4', assumed='#e08a00', open='#8c8c8c', ours='#c2185b')
plt.rcParams.update({'font.size': 10})
def box(ax, x, y, w, h, text, col, fs=9):
    ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle='round,pad=0.02', fc=col, ec='k', alpha=0.18))
    ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle='round,pad=0.02', fc='none', ec=col, lw=2))
    ax.text(x + w/2, y + h/2, text, ha='center', va='center', fontsize=fs, wrap=True)
def key(ax, y=-0.06):
    for i, (k, v) in enumerate([('proved', 'proved'), ('observed', 'observed'), ('assumed', 'assumed or put in'), ('ours', 'ours (a reading)'), ('open', 'open')]):
        ax.add_patch(Rectangle((0.02 + 0.19*i, y), 0.025, 0.03, color=C[k], transform=ax.transAxes, clip_on=False))
        ax.text(0.05 + 0.19*i, y + 0.015, v, transform=ax.transAxes, va='center', fontsize=8)

# w1 the map of the paper
fig, ax = plt.subplots(figsize=(11, 6.2)); ax.set_xlim(0, 11); ax.set_ylim(0, 6.4); ax.axis('off')
ax.text(5.5, 6.2, 'w1. The bridge, not a unification: three sides meet at one shift rule', ha='center', fontsize=12, weight='bold')
box(ax, 0.2, 4.3, 3.2, 1.4, 'Part 1. Why quantum needs i\nenergy zero is a convention\n-> one shared J, J^2 = -1', C['proved'])
box(ax, 0.2, 2.5, 3.2, 1.5, 'Part 2. Information side\nlinear mixing -> k = 0\nquantum interference stops at 2\nball 2.83, square 4', C['proved'])
box(ax, 0.2, 1.0, 3.2, 1.2, 'three-slit data:\n|k| < 0.07', C['observed'])
box(ax, 3.9, 2.6, 3.2, 2.2, 'Part 3. The hinge\na constant shift is a convention\nexactly when its generator\ncommutes with everything\n+ the action picture', C['ours'], 9.5)
box(ax, 7.6, 4.3, 3.2, 1.4, 'Part 4. Gravity weighs the\nenergy zero: rho * volume\n(cosmological constant)', C['proved'])
box(ax, 7.6, 2.5, 3.2, 1.5, 'spinning black hole from a\nshift of the source by i*a\n(Newman-Janis; checked,\nnot explained)', C['assumed'])
box(ax, 7.6, 1.0, 3.2, 1.2, 'spinning bubble: potential zero\nis physical -> spin needs charge', C['proved'])
box(ax, 3.9, 0.35, 3.2, 1.6, 'Part 5. Open\nquantum gravity not reached\nlead: fixed 4-volume\n(unimodular) restores the\nconvention and a clock', C['open'])
for (x0, y0, x1, y1) in [(3.4, 5.0, 3.9, 4.2), (3.4, 3.2, 3.9, 3.5), (7.1, 4.2, 7.6, 5.0), (7.1, 3.4, 7.6, 3.2), (7.1, 3.0, 7.6, 1.6), (5.5, 2.6, 5.5, 1.95)]:
    ax.annotate('', (x1, y1), (x0, y0), arrowprops=dict(arrowstyle='->', lw=1.5))
key(ax, 0.0); fig.savefig(f'{D}/w1-the-bridge-map.png', dpi=150, bbox_inches='tight'); plt.close()

# w2 one shift rule, three uses
fig, axs = plt.subplots(1, 3, figsize=(13, 3.9)); cs = np.linspace(0, 10, 200)
rng = np.random.default_rng(1)
from scipy.linalg import expm
def herm(n): a_ = rng.normal(size=(n, n)) + 1j*rng.normal(size=(n, n)); return (a_ + a_.conj().T)/2
Hm, Gm, Om = herm(4), herm(4), herm(4); ps = rng.normal(size=4) + 1j*rng.normal(size=4); ps /= np.linalg.norm(ps)
ev = lambda M: np.real(np.vdot(expm(-1j*M*1.3) @ ps, Om @ expm(-1j*M*1.3) @ ps))
axs[0].plot(cs, [ev(Hm + c_*np.eye(4)) for c_ in cs], color=C['proved'], lw=2.5, label='shift by c * 1: unchanged (computed)')
axs[0].plot(cs, [ev(Hm + c_*Gm) for c_ in cs], color=C['open'], lw=1.5, ls='--', label='shift by c * G, G not commuting (computed)')
axs[0].set(title='A. Quantum: energy zero is a convention', xlabel='constant added to every energy, c', ylabel='a measured average'); axs[0].legend(fontsize=7.5, loc='lower left')
axs[1].plot(cs, 0*cs, color=C['open'], lw=2, label='Abelian field: E = 0')
axs[1].plot(cs, 0.086*cs, color=C['proved'], lw=2.5, label='non-Abelian: |E| = |[A, phi]| grows')
axs[1].set(title='B. Spinning bubble: potential zero is physical', xlabel='constant potential, c0', ylabel='field from the constant (one random A)'); axs[1].legend(fontsize=8)
t = np.linspace(0, 2, 200)
for r0, ls in [(0.01, ':'), (0.04, '--'), (0.16, '-')]:
    axs[2].plot(t, np.exp(np.sqrt(8*np.pi*r0/3)*t), color=C['proved'], ls=ls, lw=2, label=f'rho = {r0}')
axs[2].set(title='C. Gravity: the energy zero is weighed', xlabel='time (G = 1)', ylabel='size of the universe a(t)'); axs[2].legend(fontsize=8)
fig.suptitle('w2. One rule, three uses: a constant shift is invisible only when it commutes with everything', weight='bold')
fig.tight_layout(); fig.savefig(f'{D}/w2-one-rule-three-uses.png', dpi=150, bbox_inches='tight'); plt.close()

# w3 the action picture
fig, axs = plt.subplots(1, 3, figsize=(13, 4.2))
tt = np.linspace(0, 1, 100)
for q, lab in [(np.sin(np.pi*tt), 'path 1'), (0.3*np.sin(2*np.pi*tt), 'path 2')]: axs[0].plot(tt, q, lw=2, label=lab)
axs[0].text(0.5, -0.92, 'each path: S -> S - c T\nsame T for every path\n-> one shared phase e^{-icT/hbar}', ha='center', fontsize=9, color=C['proved'])
axs[0].set(title='A. Quantum: the constant is a shared phase', xlabel='time (0 to T)', ylabel='position'); axs[0].set_ylim(-1.25, 1.2); axs[0].legend(fontsize=8)
axs[1].add_patch(Rectangle((0.1, 0.2), 0.3, 0.5, fc=C['ours'], alpha=0.25, ec='k')); axs[1].text(0.25, 0.45, 'volume V1', ha='center')
axs[1].add_patch(Rectangle((0.55, 0.1), 0.4, 0.75, fc=C['ours'], alpha=0.45, ec='k')); axs[1].text(0.75, 0.47, 'volume V2', ha='center')
axs[1].text(0.5, -0.08, 'each geometry: S -> S - rho V\nV differs -> the constant pushes:\nT_mn = -rho g_mn (a cosmological constant)', ha='center', va='top', fontsize=9, color=C['proved'])
axs[1].set(title='B. Gravity: the constant is weighed by volume', xlim=(0, 1), ylim=(-0.35, 1)); axs[1].axis('off')
for x0 in (0.1, 0.55): axs[2].add_patch(Rectangle((x0, 0.2), 0.35, 0.55, fc=C['open'], alpha=0.3, ec='k'))
axs[2].text(0.275, 0.47, 'V0', ha='center'); axs[2].text(0.725, 0.47, 'V0', ha='center')
axs[2].text(0.5, -0.08, 'fix the total volume (unimodular):\nS -> S - rho V0 for every geometry\n-> a shared factor again, and a clock', ha='center', va='top', fontsize=9, color=C['open'])
axs[2].set(title='C. The lead: fixed volume restores the convention', xlim=(0, 1), ylim=(-0.35, 1)); axs[2].axis('off')
fig.suptitle('w3. Is it action? Yes. A constant energy adds (constant) x (how much spacetime) to the action', weight='bold')
fig.tight_layout(); fig.savefig(f'{D}/w3-the-action-picture.png', dpi=150, bbox_inches='tight'); plt.close()

# w4 interference orders on the dial
fig, ax = plt.subplots(figsize=(7.5, 4.4)); k = np.linspace(-1, 1, 300); w = 1/3
ax.plot(k, k*w*w, color=C['ours'], lw=2, label='law: two-path term k w^2')
ax.plot(k, k**2*w**3, color=C['ours'], lw=2, ls='--', label='law: three-path term k^2 w^3')
ax.axhline(0, color=C['proved'], lw=3, alpha=0.6, label='quantum: three-path term is exactly 0 (Sorkin)')
ax.axvspan(-0.073, 0.073, color=C['observed'], alpha=0.25, label='three-slit data allow |k| < 0.073')
ax.set(xlabel='the law\'s dial k', ylabel='interference term (equal slits, w = 1/3)', title='w4. Interference order: the law at k != 0 has a three-path term, quantum does not')
ax.legend(fontsize=8); fig.tight_layout(); fig.savefig(f'{D}/w4-interference-orders.png', dpi=150, bbox_inches='tight'); plt.close()

# w5 ball vs square
fig, axs = plt.subplots(1, 2, figsize=(11, 4.4))
th = np.linspace(0, 2*np.pi, 300)
axs[0].add_patch(Rectangle((-1, -1), 2, 2, fc=C['open'], alpha=0.15, ec=C['open'], lw=2))
axs[0].plot(np.cos(th), np.sin(th), color=C['proved'], lw=2.5)
axs[0].plot([1, 0, -1, 0, 1], [0, 1, 0, -1, 0], color=C['assumed'], lw=2)
axs[0].text(0, 0, 'circle a = 2: quantum\nsquare a -> inf: PR box\ndiamond a = 1: score 2', ha='center', va='center', fontsize=8)
axs[0].set(aspect='equal', title='A. State shapes: the circle is quantum', xlim=(-1.3, 1.3), ylim=(-1.3, 1.3))
a = np.linspace(1, 8, 200); axs[1].plot(a, 2**(2 - 1/a), color=C['ours'], lw=2)
for val, lab, col in [(2, 'classical 2 (a = 1)', C['assumed']), (2*np.sqrt(2), 'quantum 2.83 (a = 2)', C['proved']), (4, 'square 4 (a -> inf)', C['open'])]:
    axs[1].axhline(val, color=col, ls='--'); axs[1].text(8.1, val, lab, va='center', fontsize=8, color=col)
axs[1].plot([2], [2*np.sqrt(2)], 'o', color=C['proved'])
axs[1].set(xlabel='shape exponent a (|x|^a + |y|^a <= 1)', ylabel='best CHSH score 2^(2 - 1/a)', title='B. Bit-flip symmetry plus smoothness picks a = 2'); axs[1].set_xlim(1, 8)
fig.suptitle('w5. Why a ball and not a square: the shape sets the Bell ceiling', weight='bold')
fig.tight_layout(); fig.savefig(f'{D}/w5-ball-vs-square.png', dpi=150, bbox_inches='tight'); plt.close()

# w6 Euclidean to real time
fig, ax = plt.subplots(figsize=(6.4, 5.4)); s = np.linspace(0, 2*np.pi, 300)
ax.plot(np.cos(s), np.sin(s), color=C['proved'], lw=2.5, label='imaginary k: e^{i c t} turns (a phase, invisible)')
x = np.linspace(0, 2.2, 50); ax.plot(np.exp(0.5*x), 0*x, color=C['ours'], lw=2.5, label='real k: e^{c tau} grows (Euclidean weight; de Sitter)')
ax.plot([0], [0], 'k+'); ax.annotate('Wick rotation\ntau -> i t', (0.75, 0.68), (1.3, 1.25), arrowprops=dict(arrowstyle='->'))
ax.set(aspect='equal', xlim=(-1.4, 3.2), ylim=(-1.4, 1.6), xlabel='real part', ylabel='imaginary part',
       title='w6. One fused shift 1 + k a: the sign of k^2 decides turn or grow')
ax.legend(fontsize=8, loc='lower right'); fig.tight_layout(); fig.savefig(f'{D}/w6-turn-or-grow.png', dpi=150, bbox_inches='tight'); plt.close()

# w7 scorecard
rows = [('A single shared i from the energy-zero convention', 'proved', 'theorem + 32 checks, parts in Lean'),
        ('Real formulations that match experiment carry that i', 'proved', 'three papers checked'),
        ('Linear mixing pins the entropy dial at k = 0', 'proved', 'Faddeev 1956; lift is ours'),
        ('Quantum interference stops at order 2; the law at k != 0 does not', 'proved', 'Sorkin 1994; Sugeno form'),
        ('Experiment caps |k| < 0.073', 'observed', 'Kauten 2017, forgiving reading'),
        ('One shift rule read in three places', 'ours', 'our reading of textbook facts'),
        ('Gravity weighs a constant energy (rho x volume)', 'proved', 'textbook; Weinberg 1989'),
        ('Spin from a shift by i a (Newman-Janis)', 'assumed', 'checked to 50 digits, unexplained'),
        ('Fixed 4-volume restores the convention and a clock', 'open', 'Unruh 1989 known; the tie to i is ours, untested'),
        ('Quantum gravity', 'open', 'not reached')]
fig, ax = plt.subplots(figsize=(11.5, 5.2)); ax.axis('off')
ax.text(0, 1.02, 'w7. Scorecard: what is proved, observed, put in, ours and open', weight='bold', fontsize=12, transform=ax.transAxes)
for i, (claim, st, src) in enumerate(rows):
    y = 0.92 - i*0.09
    ax.add_patch(Rectangle((0, y - 0.03), 0.02, 0.06, color=C[st], transform=ax.transAxes))
    ax.text(0.03, y, claim, va='center', transform=ax.transAxes, fontsize=9.5)
    ax.text(0.62, y, st, va='center', color=C[st], weight='bold', transform=ax.transAxes, fontsize=9)
    ax.text(0.72, y, src, va='center', transform=ax.transAxes, fontsize=8.5)
fig.savefig(f'{D}/w7-scorecard.png', dpi=150, bbox_inches='tight'); plt.close()
print('wrote', sorted(os.listdir(D)))
