"""Figures for the Gram-test and prime-drift DRAFT (torus, Klein bottle, odd half)."""
import json, numpy as np, matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch, Rectangle

C = {'torus': '#257de0', 'klein': '#f47230', 'odd': '#00b894', 'zeta': '#f2a600'}
NAME = {'torus': 'torus', 'klein': 'Klein bottle', 'odd': 'odd half', 'zeta': 'Riemann zeta (matched spacing)'}
plt.rcParams.update({'font.size': 11, 'axes.spines.top': False, 'axes.spines.right': False,
                     'axes.grid': True, 'grid.color': '#e6e6e6', 'grid.linewidth': 0.6})
Z = json.load(open('zeros.json')); G = json.load(open('gram_2000.json'))
W = json.load(open('gram_windows.json')); D = json.load(open('drift.json'))
cz = lambda L: np.array([complex(a, b) for a, b in L])

# ---------- fig0: the set-up ----------
fig, ax = plt.subplots(figsize=(11, 4.2)); ax.axis('off'); ax.set_xlim(0, 11); ax.set_ylim(0, 4.2)
ax.add_patch(Rectangle((0.4, 0.9), 2.4, 2.4, fill=False, lw=2, ec=C['torus']))
for y in (1.5, 2.7):
    ax.annotate('', xy=(2.8, y), xytext=(0.4, y), arrowprops=dict(arrowstyle='->', color='#999', lw=1))
ax.text(1.6, 3.55, 'square torus', ha='center', color='#333')
ax.text(1.6, 0.45, 'opposite edges glued straight', ha='center', fontsize=9, color='#666')
ax.add_patch(Rectangle((4.3, 0.9), 2.4, 2.4, fill=False, lw=2, ec=C['klein']))
ax.annotate('', xy=(6.7, 1.3), xytext=(4.3, 2.9), arrowprops=dict(arrowstyle='->', color=C['klein'], lw=1.4))
ax.annotate('', xy=(6.7, 2.9), xytext=(4.3, 1.3), arrowprops=dict(arrowstyle='->', color=C['klein'], lw=1.4))
ax.text(5.5, 3.55, 'Klein bottle = torus with one twist', ha='center', color='#333')
ax.text(5.5, 0.45, 'glide g(x, y) = (x + 1/2, −y)', ha='center', fontsize=9, color='#666')
txt = ('torus     Z_T(s) = 4 ζ(s) L(s, χ₋₄)       (Euler product, zeros on the line)\n'
       'Klein     Z_K(s) = Z_T(s)/2 − η(2s)        (tones the twist keeps)\n'
       'odd half  Z_O(s) = Z_T(s)/2 + η(2s)        (tones the twist flips)\n'
       '          Z_K + Z_O = Z_T')
ax.text(7.2, 2.1, txt, family='monospace', fontsize=8.6, va='center', color='#222')
ax.text(5.5, 4.05, 'Three spectral zetas from one square (Gram test and prime drift, DRAFT)', ha='center', fontsize=13, weight='bold')
fig.savefig('figures/fig0-setup.png', dpi=150, bbox_inches='tight'); plt.close(fig)

# ---------- fig1: where the zeros are ----------
fig, axs = plt.subplots(1, 3, figsize=(13, 6), sharey=True)
for ax, k in zip(axs, ('torus', 'klein', 'odd')):
    z = cz(Z[k]['low'])
    ax.axvspan(-3, -1, color='#f3f3f3', zorder=0); ax.axvline(0.5, color='#888', lw=1, ls='--')
    ax.scatter(z.real, z.imag, s=6, color=C[k], lw=0)
    ax.set_xlim(-3, 2); ax.set_ylim(0, 600); ax.set_xlabel('real part')
    n_left = int((z.real < -1).sum()); n_right = int((z.real > 1).sum())
    ax.set_title(f"{NAME[k]}: {len(z)} zeros\n{n_left} far left, {n_right} right of Re = 1", fontsize=11)
axs[0].set_ylabel('height (imaginary part)')
axs[0].text(-2.9, 560, 'far-left\nband', fontsize=9, color='#777')
fig.suptitle('Where the zeros sit (height up to 600). Dashed line: Re = 1/2.', y=1.0)
fig.savefig('figures/fig1-zero-map.png', dpi=150, bbox_inches='tight'); plt.close(fig)

# ---------- fig2: near-line family drifts toward the line ----------
fig, axs = plt.subplots(1, 2, figsize=(12, 4.2), sharey=True)
bins = np.linspace(-0.45, 0.85, 53)
for ax, k in zip(axs, ('klein', 'odd')):
    lo = cz(Z[k]['low']); hi = cz(Z[k]['high'])
    a = lo.real[(lo.imag > 400) & (lo.real > -1) & (lo.real < 1.5)] - 0.5
    b = hi.real[hi.real < 1.5] - 0.5
    ax.hist(a, bins=bins, density=True, color=C[k], alpha=0.35, label=f'height 400–600 (mean |offset| {np.mean(np.abs(a[np.abs(a)<0.5])):.3f})')
    ax.hist(b, bins=bins, density=True, histtype='step', lw=2, color=C[k], label=f'height 1900–2000 (mean |offset| {np.mean(np.abs(b[np.abs(b)<0.5])):.3f})')
    ax.axvline(0, color='#888', ls='--', lw=1); ax.set_xlabel('real part − 1/2'); ax.legend(fontsize=9, frameon=False)
    ax.set_title(f'{NAME[k]}: near-line zeros creep toward the line, never onto it', fontsize=11)
axs[0].set_ylabel('share of zeros (density)')
fig.savefig('figures/fig2-offset-shrinks.png', dpi=150, bbox_inches='tight'); plt.close(fig)

# ---------- fig3: Gram rates ----------
fig, axs = plt.subplots(1, 2, figsize=(13, 4.6), gridspec_kw={'width_ratios': [1.15, 1]})
order = ['zeta', 'odd', 'torus', 'klein']
last = W['windows'][-1]
mid = G['strips']['middle']['counts']; rt = G['strips']['right']['counts']
allc = {k: np.array(mid[k]) + np.array(rt[k]) for k in ('torus', 'klein', 'odd')}
def st(c): return [np.mean(c == 1), np.mean(c == 0), np.mean(c >= 2)]
vals = {k: st(allc[k]) for k in allc}
zc = W['windows']; zeta_one = np.mean([w['zeta']['one'] for w in zc]); zeta_empty = np.mean([w['zeta']['empty'] for w in zc])
vals['zeta'] = [zeta_one, zeta_empty, 1 - zeta_one - zeta_empty]
y = np.arange(len(order))[::-1]
for i, k in enumerate(order):
    left = 0
    for j, (v, lab, alpha) in enumerate(zip(vals[k], ('exactly one zero', 'empty', 'two or more'), (1, 0.45, 0.2))):
        axs[0].barh(y[i], v, left=left, color=C[k], alpha=alpha, height=0.6, edgecolor='white', lw=2)
        if v > 0.07: axs[0].text(left + v/2, y[i], f'{100*v:.0f}%', ha='center', va='center', fontsize=9, color='#111')
        left += v
axs[0].set_yticks(y); axs[0].set_yticklabels([NAME[k] for k in order]); axs[0].set_xlim(0, 1)
axs[0].set_xlabel('share of Gram intervals (heights 10–2000; zeta at matched spacing)')
axs[0].set_title('Gram test: dark = exactly one zero, mid = empty, pale = two or more', fontsize=11); axs[0].grid(False)
mids = [np.mean(w['torus_window']) for w in zc]
for k in order:
    v = [w['zeta']['one'] if k == 'zeta' else w[k]['one'] for w in zc]
    axs[1].plot(mids, v, '-o', color=C[k], lw=2, ms=7, label=NAME[k])
axs[1].set_ylim(0.3, 0.9); axs[1].set_xlabel('height of the torus window'); axs[1].set_ylabel('share with exactly one zero')
axs[1].legend(fontsize=9, frameon=False, loc='lower left'); axs[1].set_title('By height window', fontsize=11)
fig.savefig('figures/fig3-gram-rates.png', dpi=150, bbox_inches='tight'); plt.close(fig)

# ---------- fig4: why: the Gram values ----------
fig, ax = plt.subplots(figsize=(10, 4.4))
gv = W['gram_values']; b = np.linspace(-15, 25, 81)
for k, scale, lab in (('torus', 0.5, 'torus ÷ 2'), ('klein', 1, 'Klein bottle'), ('odd', 1, 'odd half')):
    v = np.array(gv[k])*scale
    ax.hist(v, bins=b, histtype='step', lw=2, color=C[k], density=True, label=f'{lab}: mean {np.mean(v):.2f}')
    ax.axvline(np.mean(v), color=C[k], lw=1.5, ls=':')
ax.axvline(0, color='#555', lw=1)
ax.text(9, 0.12, 'right of 0: the Gram point\nsits on the expected side', fontsize=9, color='#555')
ax.set_xlabel('Gram value  (−1)ⁿ Re[e^{iθ(gₙ)} F(1/2 + i gₙ)]'); ax.set_ylabel('density')
ax.set_title('Why: the twist shifts every Gram value by ∓1 (the leading tone of η), not the zeros’ distance from the line', fontsize=11)
ax.legend(frameon=False, fontsize=9)
fig.savefig('figures/fig4-gram-values.png', dpi=150, bbox_inches='tight'); plt.close(fig)

# ---------- fig5: lattice count ----------
X = np.array(D['grid'], float)
fig, ax = plt.subplots(figsize=(10, 4.2))
NT = np.array(D['torus']['N']); NK = np.array(D['klein']['N'])
ax.plot(X, (NT - (np.pi*X - 1))/X**0.25, color=C['torus'], lw=1.5, label='torus: (N − πX + 1)/X^¼')
ax.plot(X, (2*NK - (np.pi*X - 1))/X**0.25, color=C['klein'], lw=1.2, ls='--', label='Klein bottle ×2, same scale')
ax.set_xscale('log'); ax.set_xlabel('lattice size X (eigenvalue bound)'); ax.set_ylabel('count minus main term, ÷ X^¼')
ax.set_title('Lattice count: the twist changes it by at most 1, so the off-line zeros leave no trace here', fontsize=11)
ax.legend(frameon=False, fontsize=9)
fig.savefig('figures/fig5-lattice-count.png', dpi=150, bbox_inches='tight'); plt.close(fig)

# ---------- fig6: prime drift ----------
fig, axs = plt.subplots(1, 2, figsize=(13, 4.6))
for k in ('torus', 'odd', 'klein'):
    p = np.array(D[k]['psi']); d = np.abs(p - X)
    axs[0].loglog(X, np.maximum(d, 1e-1), color=C[k], lw=1.4, label=NAME[k])
bK = max(cz(Z['klein']['low']).real); bO = max(cz(Z['odd']['low']).real)
for e, lab, col in ((0.5, 'X^½ (zeros on the line)', '#888'), (bK, f'X^{bK:.3f} (Klein’s rightmost zero)', C['klein']), (bO, f'X^{bO:.3f} (odd half’s rightmost zero)', C['odd'])):
    axs[0].loglog(X, 0.6*X**e if e > 1 else 1.5*X**e, ls=':', color=col, lw=1.2, label=lab)
axs[0].set_xlabel('size X'); axs[0].set_ylabel('|ψ(X) − X|'); axs[0].legend(fontsize=8.5, frameon=False)
axs[0].set_title('Prime count drift: Klein outruns the main term itself', fontsize=11)
zk = cz(Z['klein']['low']); zk = zk[zk.real > -1]
p1 = np.array(D['klein']['psi1'])
pred = -2*np.real(np.sum(X[:, None]**(zk[None, :] + 1)/(zk*(zk + 1))[None, :], axis=1))
sel = X > 1e3
axs[1].plot(X[sel], p1[sel]/X[sel]**2.2855, color=C['klein'], lw=2.2, label='counted from the lattice')
axs[1].plot(X[sel], pred[sel]/X[sel]**2.2855, color='#222', lw=1, ls='--', label=f'rebuilt from {len(zk)} zeros (height ≤ 600)')
axs[1].set_xscale('log'); axs[1].set_xlabel('size X'); axs[1].set_ylabel('∫(ψ − x)dx ÷ X^2.2855')
axs[1].set_title('Klein: the drift is exactly its zeros (explicit formula)', fontsize=11); axs[1].legend(fontsize=9, frameon=False)
fig.savefig('figures/fig6-prime-drift.png', dpi=150, bbox_inches='tight'); plt.close(fig)
print('figures written')
