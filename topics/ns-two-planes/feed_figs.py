"""DRAFT. Round 3 pictures: python3 feed_figs.py -> figures/fig11, fig12 (and feed_summary.json)"""
import json, glob, os, numpy as np, matplotlib
matplotlib.use('Agg'); import matplotlib.pyplot as plt
L = lambda f: json.load(open(f))
rows = [L(f) for f in glob.glob('runs/feed/f_*.json')]
best = {}
for r in rows:
    key = (r['cls'], r['N'], r['nu'])
    if key not in best or r['R'] > best[key]['R']: best[key] = r
summ = {f'{c} N={N} nu={nu}': dict(R=b['R'], P=b['P'], k_rms=b['k_rms'], seeds=sum(1 for r in rows if (r['cls'], r['N'], r['nu']) == (c, N, nu)),
        spread=max(r['R'] for r in rows if (r['cls'], r['N'], r['nu']) == (c, N, nu)) - min(r['R'] for r in rows if (r['cls'], r['N'], r['nu']) == (c, N, nu)),
        M_eig=b.get('M_eig')) for (c, N, nu), b in sorted(best.items())}
json.dump(summ, open('feed_summary.json', 'w'), indent=1)
col = {'any': '#1f618d', 'two': '#c0392b', 'one': '#7d7d7d'}
lab = {'any': 'any 3D flow', 'two': 'two planes only (no spin along the shared line\'s normal... i.e. none along y)', 'one': 'one layered plane only'}
lab = {'any': 'any 3D flow', 'two': 'two planes sharing a box axis', 'twod': 'two planes sharing the box diagonal', 'one': 'one layered plane only'}
col['twod'] = '#e67e22'
mk = {32: 'o', 64: 's', 96: 'D'}
fig, axs = plt.subplots(1, 2, figsize=(11.5, 4.4))
for c in ['any', 'twod', 'two', 'one']:
    for N in [32, 64, 96]:
        pts = sorted((nu, b['R']) for (cc, NN, nu), b in best.items() if cc == c and NN == N)
        if not pts: continue
        x = [1/nu for nu, _ in pts]; y = [R for _, R in pts]
        axs[0].plot(x, y, marker=mk[N], color=col[c], ls='-' if N == 64 else ':', ms=7 if c == 'twod' else 4,
                    mfc='none' if c == 'twod' else col[c], zorder=3 if c == 'any' else 2,
                    label=f'{lab[c]}, {N}³' if c != 'one' or N == 32 else None)
axs[0].axhline(0, color='k', lw=0.6)
axs[0].set_xscale('log'); axs[0].set_xlabel('1 / viscosity  (total spin fixed at 1)')
axs[0].set_ylabel('fastest possible spin growth  (dZ/dt)/2')
axs[0].set_title('Fastest growth of total spin at one instant\n(below 0 = no flow of that kind can grow at all)', fontsize=10)
axs[0].legend(fontsize=7)
axs[0].text(0.98, 0.40, 'dotted 32³ points at small viscosity are too coarse;\nthe 96³ points show the rate itself is not settled at 64³ either.\nCompare classes on the same grid.', transform=axs[0].transAxes, ha='right', va='bottom', fontsize=7, bbox=dict(fc='w', ec='0.8'))
for N in [32, 64, 96]:
    pts = sorted((nu, best[('two', N, nu)]['R']/best[('any', N, nu)]['R']) for (cc, NN, nu) in best
                 if cc == 'any' and NN == N and ('two', N, nu) in best and best[('any', N, nu)]['R'] > 0)
    if pts:
        axs[1].plot([1/nu for nu, _ in pts], [q for _, q in pts], marker=mk[N], color='#c0392b', ls='-' if N == 64 else ':', label=f'sharing a box axis, {N}³')
    pts = sorted((nu, best[('twod', N, nu)]['R']/best[('any', N, nu)]['R']) for (cc, NN, nu) in best
                 if cc == 'any' and NN == N and ('twod', N, nu) in best and best[('any', N, nu)]['R'] > 0)
    if pts:
        axs[1].plot([1/nu for nu, _ in pts], [q for _, q in pts], marker=mk[N], color='#e67e22', ls='-' if N == 64 else ':', label=f'sharing the box diagonal, {N}³')
axs[1].axhline(1, color='k', lw=0.6, ls='--'); axs[1].set_xscale('log'); axs[1].set_ylim(-0.5, 1.1)
axs[1].text(0.02, 0.05, 'at viscosity 0.1 and 0.07 the ratio is below 0:\ntwo planes on a box axis cannot grow at all there', transform=axs[1].transAxes, fontsize=7.5)
axs[1].set_xlabel('1 / viscosity'); axs[1].set_ylabel('two planes ÷ any flow')
axs[1].set_title('How close two planes come to the fastest 3D growth\n(1 = no limit at all from having only two planes)', fontsize=10)
axs[1].legend(fontsize=8)
fig.suptitle('Round 3 (DRAFT): two planes feed each other as fast as any 3D flow can grow its spin; on the box diagonal they match it within 0.2%', fontsize=10.5)
fig.tight_layout(); fig.savefig('figures/fig11-fastest-feeding.png', dpi=140); plt.close(fig)
# fig12: the fastest-feeding flows themselves (saved fields), spin magnitude on the slice through its peak
fs = {c: f'runs/feed/f_n64_nu0.025_{c}_s0_field.npz' for c in ['any', 'twod', 'two']}
if all(os.path.exists(f) for f in fs.values()):
    fig, axs = plt.subplots(3, 3, figsize=(10.5, 10.2))
    short = {'any': 'any 3D flow', 'twod': 'two planes, diagonal', 'two': 'two planes, y axis'}
    for row, c in enumerate(['any', 'twod', 'two']):
        w = np.load(fs[c])['w']; m = np.sqrt((w**2).sum(0)); N = m.shape[0]
        i, j, k = np.unravel_index(np.argmax(m), m.shape)
        m = np.roll(np.roll(np.roll(m, N//2 - i, 0), N//2 - j, 1), N//2 - k, 2)   # centre the peak
        h = N//5; c0 = N//2; sl_ = slice(c0 - h, c0 + h); ext = [-h*2*np.pi/N, h*2*np.pi/N]*2
        for col_, (img, xl, yl) in enumerate([(m[c0, sl_, sl_].T, 'y', 'z'), (m[sl_, c0, sl_].T, 'x', 'z'), (m[sl_, sl_, c0].T, 'x', 'y')]):
            axs[row, col_].imshow(img, origin='lower', extent=ext, cmap='magma')
            axs[row, col_].set_title(f'{short[c]}: slice {"xyz"[col_]} = peak', fontsize=9)
            axs[row, col_].set_xlabel(xl); axs[row, col_].set_ylabel(yl)
        comp = (w**2).mean(axis=(1, 2, 3)); comp = comp/comp.sum()
        axs[row, 0].text(0.03, 0.96, f'spin share x,y,z\n{comp[0]:.2f}, {comp[1]:.2f}, {comp[2]:.2f}', transform=axs[row, 0].transAxes,
                         va='top', color='w', fontsize=8)
    axs[2, 1].text(0.03, 0.04, 'a ring around the y axis:\ntwo colliding vortex rings', transform=axs[2, 1].transAxes, color='w', fontsize=8)
    fig.suptitle('Size of the spin in the flows that grow total spin fastest (viscosity 0.025, total spin 1, 64³; zoomed on the peak).\n'
                 'Top: any flow. Middle: two planes sharing the box diagonal. Bottom: two planes sharing the y axis.\n'
                 'All three are colliding vortex rings. The top one has no spin along a diagonal, so it is the same kind as the middle one.', fontsize=9.5)
    fig.tight_layout(); fig.savefig('figures/fig12-fastest-flows.png', dpi=140); plt.close(fig)
print(json.dumps(summ, indent=0)[:3000])
