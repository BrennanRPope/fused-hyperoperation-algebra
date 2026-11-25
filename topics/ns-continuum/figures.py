import json, numpy as np, matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
plt.rcParams.update({'font.size': 11, 'axes.spines.top': False, 'axes.spines.right': False})
C = dict(blue='#2b6cb0', orange='#dd6b20', green='#2f855a', red='#c53030', gray='#718096', purple='#6b46c1')

# ---------- F1: the fan ----------
fig, ax = plt.subplots(figsize=(7.6, 5.6))
n = np.linspace(-2, 6.2, 200)
ax.fill_between(n, (n + 2) / 4, 2.4, color=C['green'], alpha=0.10, lw=0)
for g, lab, col, ls in [(1/2, 'critical ray α=(n+2)/4  (Lions: smooth on and above)', C['green'], '-'),
                       (2/5, 'ray α=(n+2)/5  (dyadic model proved smooth; 3D NS sits on it)', C['blue'], '--'),
                       (1/3, 'ray α=(n+2)/6  (dyadic blow-up below, Cheskidov; numerics: the true dyadic edge)', C['red'], ':')]:
    ax.plot(n, g * (n + 2) / 2 * 1, color=col, ls=ls, lw=2, label=lab)
ax.plot([2, 2], [0, 2.4], color=C['purple'], lw=5, alpha=0.55, solid_capstyle='butt',
        label='2D column: smooth for every α ≥ 0 (stretching corner empty)')
pts = [(3, 1, '3D NS', C['red']), (2, 1, '2D NS', C['purple']), (3, 1.25, 'Lions 5/4', C['green']),
       (2, 0.8, '2D, α=4/5\n(same ray as 3D NS)', C['purple']), (4, 1, '4D', C['gray']), (5, 1, '5D', C['gray'])]
for x, y, t, c in pts:
    ax.plot(x, y, 'o', color=c, ms=8, zorder=5)
    ax.annotate(t, (x, y), xytext=(6, -14 if t.startswith('2D, ') else 6), textcoords='offset points', fontsize=9.5)
ax.plot([2, 3], [1, 1], color=C['orange'], lw=3, alpha=0.8)
ax.annotate('fractional dimension at α=1\n(slides 2D→3D across rays)', (2.5, 1), xytext=(2.15, 0.45), fontsize=9,
            arrowprops=dict(arrowstyle='->', color=C['orange']))
ax.plot(-2, 0, 'k*', ms=12); ax.annotate('pivot (−2, 0)', (-2, 0), xytext=(6, 6), textcoords='offset points', fontsize=9.5)
ax.set_xlim(-2.2, 6.2); ax.set_ylim(0, 2.2)
ax.set_xlabel('dimension n'); ax.set_ylabel('dissipation power α   (−Δ)^α')
ax.legend(loc='upper left', fontsize=8.3, frameon=False)
ax.set_title('Every scaling threshold is a ray from one pivot')
fig.tight_layout(); fig.savefig('figures/f1_fan.png', dpi=150); plt.close(fig)

# ---------- F2: rung height of the dissipation deficiency ----------
d = json.load(open('rung_height.json'))
L = np.array(d['L']); log10s = L / np.log(10)
fig, ax = plt.subplots(figsize=(7.6, 5.2))
ax.plot(log10s, np.log10(L ** 0.25), color=C['green'], lw=2.5, label="Tao 2009, proved: D = (log s)^¼")
ax.plot(log10s, np.log10(L ** 0.5), color=C['green'], lw=2, ls='--', label="Tao's conjecture (proved for dyadic): (log s)^½")
cols = {'0.25': '#90cdf4', '0.5': '#4299e1', '0.75': '#2b6cb0'}
for t in ['0.25', '0.5', '0.75']:
    ax.plot(log10s, d['fam'][t], color=cols[t], lw=2, label=f'rung {t}:  E^{t}(½ log s)')
ax.plot(log10s, d['fam']['1.0'], color=C['red'], lw=2.5, label='3D NS:  s^½  (rung 1)')
ax.plot(log10s, d['fam']['0.0'], color=C['gray'], lw=1.5, ls=':', label='rung 0:  ½ log s')
ax.set_xlabel('frequency  log₁₀ s'); ax.set_ylabel('deficiency  log₁₀ D(s)   (dissipation = s^{5/2} / D)')
ax.set_ylim(-0.3, 30); ax.legend(fontsize=8.5, frameon=False, loc='upper left')
ax.set_title('The gap is one full rung of the circle-operation flow wide')
fig.tight_layout(); fig.savefig('figures/f2_rung_height.png', dpi=150); plt.close(fig)

# ---------- F3: Osgood integrals ----------
fig, ax = plt.subplots(figsize=(7.0, 4.4))
Lo = L
for name, col, ls in [('Tao edge L^(1/4)', C['green'], '-'), ('conj edge L^(1/2)', C['green'], '--'),
                      ('rung t=0.25', '#90cdf4', '-'), ('rung t=0.5', '#4299e1', '-'), ('rung t=1.0', C['red'], '-')]:
    I = np.array(d['osg'][name]['I4']); ax.plot(Lo[:len(I)], I, color=col, ls=ls, lw=2, label=name.replace('L^', '(log s)^'))
ax.set_xscale('log'); ax.set_xlabel('log s'); ax.set_ylabel('∫ d(log s) / D⁴   (Tao\'s test)')
ax.legend(fontsize=9, frameon=False); ax.set_title('Only the bottom rung keeps the integral growing')
fig.tight_layout(); fig.savefig('figures/f3_osgood.png', dpi=150); plt.close(fig)

# ---------- F4: dyadic dimension continuum ----------
r = json.load(open('viscous_refine.json'))
fig, ax = plt.subplots(figsize=(7.0, 4.6))
ns = sorted(set(o['n'] for o in r)); cm = plt.cm.viridis(np.linspace(0, 0.9, len(ns)))
for c, nn in zip(cm, ns):
    rows = [o for o in r if o['n'] == nn]
    ax.plot([o['N'] for o in rows], [o['deepest'] for o in rows], 'o-', color=c, lw=2, label=f'n = {nn}  (γ = {rows[0]["gam"]:.3f})')
Ns = np.array([24, 56]); ax.plot(Ns, Ns, 'k:', lw=1); ax.text(50, 52, 'deepest = N\n(blow-up)', fontsize=8.5)
ax.set_xlabel('number of shells N'); ax.set_ylabel('deepest shell the cascade fills')
ax.legend(fontsize=8.5, frameon=False, loc='upper left'); ax.set_title('Dyadic model at α = 1: it breaks exactly at n = 4')
fig.tight_layout(); fig.savefig('figures/f4_dyadic.png', dpi=150); plt.close(fig)

# ---------- F5: annulus as cube layers ----------
fig, axs = plt.subplots(1, 2, figsize=(9.0, 4.2), gridspec_kw=dict(width_ratios=[1, 1.15]))
ax = axs[0]
for i in range(8):
    for j in range(8):
        hi = (i >= 4) + (j >= 4)
        col = ['#edf2f7', '#90cdf4', '#2b6cb0'][hi]
        ax.add_patch(Rectangle((i, j), 1, 1, fc=col, ec='white', lw=1))
ax.set_xlim(0, 8); ax.set_ylim(0, 8); ax.set_aspect('equal'); ax.axis('off')
ax.set_title('Doubling a 4×4 box: 48 new cells\n= 16 · 1^(⊞2) = 16 · (2 strips + 1 corner)', fontsize=10)
ax = axs[1]
dims = [1, 2, 3, 4, 5]
from math import comb
bottom = np.zeros(len(dims))
cols = ['#bee3f8', '#90cdf4', '#4299e1', '#2b6cb0', '#1a365d']
for k in range(1, 6):
    vals = np.array([comb(nn, k) if k <= nn else 0 for nn in dims])
    ax.bar(dims, vals, bottom=bottom, color=cols[k - 1], label=f'layer {k}')
    bottom += vals
for nn, b in zip(dims, bottom): ax.text(nn, b + 0.5, f'{int(b)}', ha='center', fontsize=9)
ax.set_xlabel('dimension n'); ax.set_ylabel('new cells per old cell = 1^(⊞n) = 2ⁿ − 1')
ax.legend(fontsize=8.5, frameon=False); ax.set_title('The shell is the circle-operation cube, graded', fontsize=10)
fig.tight_layout(); fig.savefig('figures/f5_annulus.png', dpi=150); plt.close(fig)
print('figures done')
