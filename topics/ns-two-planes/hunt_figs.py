"""DRAFT. Blow-up hunt pictures: python3 hunt_figs.py -> figures/fig15, fig16 (and hunt_summary.json)"""
import json, numpy as np, matplotlib
matplotlib.use('Agg'); import matplotlib.pyplot as plt
import hunt_an as H
R = H.load_all()
FOLD = {'tg': 1.0, 'rings': 0.1105}                     # exact fold time of pressureless dust (1 / max squeeze rate at t = 0)
grids = [64, 128, 256]; col = {64: '#9ecae1', 128: '#3182bd', 256: '#08306b'}
summ = {f'{ic} N={N} nu={nu} theta={th}': dict(t_trust=H.t_trust(r), t_end=r['t_end'], stop=r['stop'],
        max_spin=max(q['wmax'] for q in r['rows']), max_grad=max(q['gmax'] for q in r['rows'])) for (ic, N, nu, th), r in sorted(R.items())}
json.dump(summ, open('hunt_summary.json', 'w'), indent=1, default=float)

fig, axs = plt.subplots(1, 3, figsize=(16, 4.8))
TH = sorted({th for (i, n, nu, th) in R if i == 'tg' and nu == 0})
ax = axs[0]
for N in grids:
    pts = [(th, H.t_trust(R[('tg', N, 0.0, th)])) for th in TH if ('tg', N, 0.0, th) in R]
    pts = [p for p in pts if p[1] and np.isfinite(p[1])]
    if pts: ax.plot([p[0] for p in pts], [p[1] for p in pts], 'o-', color=col[N], label=f'{N}³: how far it can follow the flow')
ext = [(th, H.grid_close(R, 'tg', th)) for th in TH]
cl = [(th, g['T']) for th, g in ext if 'ratio' in g and g['ratio'] < 0.75]   # clear closing-in only
if cl: ax.plot([p[0] for p in cl], [p[1] for p in cl], 's', color='#c0392b', ms=8, mfc='none', mew=2, label='blow-up time the grids close in on')
ax.plot([0], [1.0], marker='*', ms=14, color='#c0392b', ls='none', label='exact fold time of dust')
ax.set_xlabel('push-back θ   (0 = pressureless dust;  1 = real flow)'); ax.set_ylabel('time')
ax.set_title('Taylor–Green, no viscosity: where each grid has to stop', fontsize=10); ax.legend(fontsize=7.5, loc='upper left')
ax = axs[1]
rr = [(th, g['ratio']) for th, g in ext if 'ratio' in g]
if rr: ax.plot([p[0] for p in rr], [p[1] for p in rr], 'o-', color='k')
ax.axhline(1, color='#2e86c1', ls='--'); ax.text(0.02, 1.02, 'each finer grid goes the same step further: no blow-up in reach', fontsize=7.5, color='#2e86c1')
ax.axhline(2**(-2/3), color='#c0392b', ls=':'); ax.text(0.02, 2**(-2/3) - 0.06, 'what a fold like dust gives (0.63)', fontsize=7.5, color='#c0392b')
ax.set_ylim(0, 1.4); ax.set_xlim(-0.03, 1.03); ax.set_xlabel('push-back θ'); ax.set_ylabel('(256³ step) ÷ (128³ step)')
ax.set_title('Do the grids close in?  below 1 = yes, on a finite blow-up time', fontsize=10)
ax = axs[2]
for N in grids:
    pts = sorted((th, H.t_trust(r)) for (i, n, nu, th), r in R.items() if i == 'rings' and n == N and nu == 0 and H.t_trust(r))
    pts = [p for p in pts if np.isfinite(p[1])]
    if pts: ax.plot([p[0] for p in pts], [p[1] for p in pts], 'o-', color=col[N], label=f'{N}³')
ax.plot([0], [FOLD['rings']], marker='*', ms=14, color='#c0392b', ls='none', label='exact fold time of dust')
ax.set_xlabel('push-back θ'); ax.set_ylabel('how far the grid can follow the flow')
ax.set_title('colliding rings (start near grid scale: support only)', fontsize=10); ax.legend(fontsize=7.5, loc='upper left')
fig.text(0.01, 0.01, 'DRAFT · Euler (no viscosity). "How far" = last time the distance to the nearest complex-time singularity is at least 2 grid spacings (δ·kmax ≥ 2).', fontsize=7.5, color='0.35')
plt.tight_layout(rect=(0, 0.04, 1, 1)); plt.savefig('figures/fig15-push-back-dial.png', dpi=130); plt.close()

fig, axs = plt.subplots(2, 3, figsize=(13.5, 7))
for c, (ic, name) in enumerate([('tg', 'Taylor–Green'), ('klein', 'bent antiparallel tube pair'), ('rings', 'colliding rings')]):
    for N in grids:
        r = R.get((ic, N, 0.0, 1.0))
        if not r: continue
        rows = [q for q in r['rows'] if np.isfinite(q['delta']) and q['delta'] > 0 and q.get('nfit', 99) >= 0.5*(N/3 - 3)]; tt = H.t_trust(r)
        t = np.array([q['t'] for q in rows]); d = np.array([q['delta'] for q in rows]); w = np.array([q['wmax'] for q in rows])
        ok = t <= (tt if tt and np.isfinite(tt) else t[-1])
        axs[0, c].semilogy(t[ok], d[ok], color=col[N], label=f'{N}³'); axs[0, c].semilogy(t[~ok], d[~ok], color=col[N], ls=':', lw=0.8)
        axs[0, c].axhline(2/(N/3), color=col[N], lw=0.6, ls='--')
        axs[1, c].plot(t[ok], 1/w[ok], color=col[N], label=f'{N}³'); axs[1, c].plot(t[~ok], 1/w[~ok], color=col[N], ls=':', lw=0.8)
    axs[0, c].set_title(f'{name} (Euler)', fontsize=10); axs[0, c].set_ylabel('δ: distance to the nearest\ncomplex-time singularity')
    axs[1, c].set_ylabel('1 / peak spin'); axs[1, c].set_xlabel('time'); axs[1, c].set_ylim(bottom=0)
    axs[0, c].legend(fontsize=7.5); axs[1, c].legend(fontsize=7.5)
fig.suptitle('DRAFT · real flow, no viscosity: a blow-up needs δ and 1/peak spin to hit zero at one time on every grid.\n'
             'Solid = resolved (δ above the dashed line for that grid); dotted = past what that grid can follow.', fontsize=10)
plt.tight_layout(rect=(0, 0, 1, 0.92)); plt.savefig('figures/fig16-real-flow-candidates.png', dpi=120); plt.close()
print('ok')
