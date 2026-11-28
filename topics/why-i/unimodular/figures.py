"""Figures m1-m5 for the unimodular page (DRAFT). Run: python3 figures.py"""
import os, numpy as np, matplotlib
matplotlib.use('Agg'); import matplotlib.pyplot as plt
from matplotlib.patches import Polygon, Rectangle
from scipy.linalg import expm, eigh
D = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'figures'); os.makedirs(D, exist_ok=True)
C = dict(proved='#2e8b57', observed='#1f6fb4', assumed='#e08a00', open='#8c8c8c', ours='#c2185b')

# m1 fixed cells, stretched space (drawn in 1+1 dimensions: cell area = a * dt; in 3+1 it is a^3 * dt)
fig, axs = plt.subplots(1, 2, figsize=(11, 5.2), sharey=True)
a = lambda t: 1 + 1.2*t
def draw(ax, slices, title, col):
    xs = np.arange(-2, 3)
    for t0, t1 in zip(slices[:-1], slices[1:]):
        for x0 in xs[:-1]:
            tt = np.linspace(t0, t1, 10)
            pts = [(x0*a(t0), t0), ((x0 + 1)*a(t0), t0), ((x0 + 1)*a(t1), t1), (x0*a(t1), t1)]
            area = (t1 - t0)*np.mean([a(t0), a(t1)])
            ax.add_patch(Polygon(pts, closed=True, fc=col, alpha=min(0.9, 0.25 + 0.9*area), ec='k', lw=0.6))
    ax.set(xlim=(-6.5, 6.5), ylim=(0, 2.05), title=title, xlabel='space (stretching outward)')
draw(axs[0], np.linspace(0, 2, 9), 'A. Equal clock ticks: cells grow\n(darker = more spacetime per cell)', C['open'])
sl = [0.0]
while sl[-1] < 2: sl.append(sl[-1] + 0.25/a(sl[-1]))
draw(axs[1], np.array(sl), 'B. Unimodular ticks: every cell the same size\n(ticks shorten as space stretches)', C['ours'])
axs[0].set_ylabel('time')
fig.suptitle('m1. "Fixed but stretched": each cell of spacetime keeps its volume, space stretches, the tick shortens', weight='bold')
fig.tight_layout(); fig.savefig(f'{D}/m1-fixed-cells-stretched-space.png', dpi=150, bbox_inches='tight'); plt.close()

# m2 unimodular time against cosmic time
fig, ax = plt.subplots(figsize=(6.8, 4.4)); t = np.linspace(0, 2, 200); Hh = 1.0
ax.plot(t, (np.exp(3*Hh*t) - 1)/(3*Hh), color=C['ours'], lw=2.5, label='de Sitter, a = e^t: T = (e^{3t} - 1)/3')
ax.plot(t, 0.5*((1 + t)**(3*2/3 + 1) - 1)/(3*2/3 + 1)*2, color=C['proved'], lw=2, ls='--', label='matter era, a = (1 + t)^{2/3}')
ax.plot(t, t, color=C['open'], lw=1, ls=':', label='T = t (no expansion)')
ax.set(xlabel='cosmic time t', ylabel='unimodular time T = four-volume to the past', yscale='log', ylim=(1e-2, 200),
       title='m2. The unimodular clock is the spacetime volume so far')
ax.legend(fontsize=8); fig.tight_layout(); fig.savefig(f'{D}/m2-unimodular-clock.png', dpi=150, bbox_inches='tight'); plt.close()

# m3 the test
Ngrid, L = 600, 30.0; h = L/Ngrid; z = h*np.arange(1, Ngrid + 1)
lap = (np.diag(-2*np.ones(Ngrid)) + np.diag(np.ones(Ngrid - 1), 1) + np.diag(np.ones(Ngrid - 1), -1))/h**2
Hu = lambda c: 3/8*(lap + np.diag(1/(4*z**2))) + np.diag(1/z) + c*np.eye(Ngrid)
u0 = np.exp(-(z - 8)**2/2)*np.exp(-1j*0.6*z); u0 /= np.linalg.norm(u0)
Ts = np.linspace(0, 6, 25)
fig, axs = plt.subplots(1, 3, figsize=(14, 4.2))
for cv, ls, lw in [(0, '-', 4), (5, '--', 2), (40, ':', 2)]:
    vals = []
    for TT in Ts: u = expm(-1j*Hu(cv)*TT) @ u0; vals.append(np.real(np.vdot(u, z*u)))
    axs[0].plot(Ts, vals, ls=ls, lw=lw, color=C['proved'] if cv == 0 else C['ours'], label=f'vacuum energy c = {cv}')
axs[0].set(xlabel='unimodular time T', ylabel='<a^3>, the size of space', title='A. Unimodular: the expansion ignores c'); axs[0].legend(fontsize=8)
w0 = eigh(Hu(0), eigvals_only=True)[-12:]; w5 = eigh(Hu(5), eigvals_only=True)[-12:]
for w in w0: axs[1].plot([0, 1], [w, w], color=C['proved'])
for w in w5: axs[1].plot([1.5, 2.5], [w, w], color=C['ours'])
for w_a, w_b in zip(w0[::3], w5[::3]): axs[1].annotate('', (1.5, w_b), (1, w_a), arrowprops=dict(arrowstyle='->', color='k', lw=0.6))
axs[1].set(xticks=[0.5, 2], xticklabels=['c = 0', 'c = 5'], ylabel='energy levels (top 12)', title='B. Every level moves by exactly c\n(same states): assumption A4')
Zm = np.diag(z)
for cv, col, ls in [(0, C['proved'], '-'), (0.05, C['open'], '--')]:
    Cm = Zm @ Hu(0); Cm = (Cm + Cm.T)/2 + cv*Zm; w, v = eigh(Cm); s = v[:, np.argmin(np.abs(w))]
    s = s*np.sign(s[np.argmax(np.abs(s))]); axs[2].plot(z, s, color=col, ls=ls, label=f'c = {cv}')
axs[2].set(xlabel='a^3', ylabel='timeless state (real)', title='C. Without fixed volume: c reshapes the state'); axs[2].legend(fontsize=8)
fig.suptitle('m3. The test: with fixed cell volume, vacuum energy is a shared phase, exactly why-i\'s assumption A4', weight='bold')
fig.tight_layout(); fig.savefig(f'{D}/m3-the-test.png', dpi=150, bbox_inches='tight'); plt.close()

# m4 the size of Lambda
fig, ax = plt.subplots(figsize=(6.8, 4.2))
names = ['naive vacuum\nenergy', 'Sorkin:\n1/sqrt(V)', 'observed']
vals = [1.0, 1.39e-122, 2.85e-122]; cols = [C['open'], C['assumed'], C['observed']]
ax.bar(names, np.log10(vals), color=cols)
for i, v in enumerate(vals): ax.text(i, np.log10(v) + (3 if v < 1 else -10), f'{v:.2g}', ha='center')
ax.set(ylabel='log10 of Lambda (Planck units)', title='m4. Lambda paired with four-volume gives the right order (heuristic)')
fig.tight_layout(); fig.savefig(f'{D}/m4-size-of-lambda.png', dpi=150, bbox_inches='tight'); plt.close()

# m5 scorecard
rows = [('Each cell of spacetime has a fixed volume; space still stretches', 'proved', 'definition; check M1'),
        ('Lambda is the conserved energy of unimodular time', 'proved', 'known (Unruh 1989); check M2'),
        ('Vacuum energy is a shared phase, so A4 holds and i is forced', 'ours', 'checks M3, M4, M7; theorem applies'),
        ('Same fact stops vacuum energy gravitating', 'proved', 'known (Smolin 2010; Weinberg 1989)'),
        ('Without fixed volume, vacuum energy reshapes the state', 'proved', 'check M6; WdW state is real'),
        ('A big-bang boundary condition is a free choice', 'open', 'check M5; Fabris et al. 2025'),
        ('Lambda ~ 1/sqrt(V) matches the observed order', 'assumed', 'Sorkin, heuristic; check M8'),
        ('Why Lambda is small is not solved (only moved)', 'open', 'Unruh 1989; Padilla-Saltas 2015'),
        ('Four-volume time fixes hypersurfaces only up to a class', 'open', 'Kuchar 1991'),
        ('Energy is unbounded below; observables may be scarce', 'open', 'Unruh-Wald 1989; Fabris et al.')]
fig, ax = plt.subplots(figsize=(11.5, 5.2)); ax.axis('off')
ax.text(0, 1.02, 'm5. Scorecard: what unimodular gravity gives the why-i story, and what it does not', weight='bold', fontsize=12, transform=ax.transAxes)
for i, (cl, st, src) in enumerate(rows):
    y = 0.92 - i*0.09
    ax.add_patch(Rectangle((0, y - 0.03), 0.02, 0.06, color=C[st], transform=ax.transAxes))
    ax.text(0.03, y, cl, va='center', transform=ax.transAxes, fontsize=9.5)
    ax.text(0.62, y, st, va='center', color=C[st], weight='bold', transform=ax.transAxes, fontsize=9)
    ax.text(0.72, y, src, va='center', transform=ax.transAxes, fontsize=8.5)
fig.savefig(f'{D}/m5-scorecard.png', dpi=150, bbox_inches='tight'); plt.close()
print(sorted(os.listdir(D)))
