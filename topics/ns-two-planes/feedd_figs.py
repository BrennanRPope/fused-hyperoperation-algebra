"""DRAFT. Round 4 pictures: python3 feedd_figs.py -> figures/fig13, fig14 (and feedd_summary.json)"""
import json, glob, numpy as np, matplotlib
matplotlib.use('Agg'); import matplotlib.pyplot as plt
L = lambda f: json.load(open(f))
rows = [L(f) for f in glob.glob('runs/feedd/g_*.json')]
best = {}
for r in rows:
    k = (r['d'], r['N'], r['nu'])
    if k not in best or r['R'] > best[k]['R']: best[k] = r
json.dump({f"d={d} N={N} nu={nu}": dict(R=b['R'], iters=b['iters'], capped=b['iters'] >= 150 and b['N'] == 32 and d == 4,
           k_rms=b['k_rms'], flat_share=b['flat_share'], vary_eig=b.get('vary_eig')) for (d, N, nu), b in sorted(best.items())},
          open('feedd_summary.json', 'w'), indent=1)
D = L('runs/drain/drain_rings.json')

fig, axs = plt.subplots(1, 2, figsize=(12, 4.6))
ax = axs[0]
col = {0.1: '#1f618d', 0.05: '#c0392b'}
for d, ls, mk in [(3, '-', 'o'), (4, '--', 's')]:
    for nu in [0.1, 0.05]:
        pts = sorted((N, b['R'], b['iters']) for (dd, N, n), b in best.items() if dd == d and n == nu)
        x = [p[0] for p in pts]; y = [p[1] for p in pts]
        ax.plot(x, y, ls=ls, marker=mk, color=col[nu], label=f'{d}D, viscosity {nu}')
        for N, R, it in pts:
            if d == 4 and N == 32:
                ax.annotate('', xy=(N, R*1.35), xytext=(N, R), arrowprops=dict(arrowstyle='->', color=col[nu]))
                ax.text(N*1.03, R, 'climb stopped early:\nat least this', fontsize=7, va='center', color=col[nu])
ax.set_xscale('log', base=2); ax.set_yscale('log'); ax.set_ylim(0.03, 6)
ax.set_xticks([16, 24, 32, 48]); ax.set_xticklabels(['16', '24', '32', '48'])
ax.set_xlabel('grid points per side (finer grid = smaller swirls allowed)')
ax.set_ylabel('fastest possible spin growth (dZ/dt)/2, spin fixed at 1')
ax.set_title('3D settles on a finite worst case; 4D keeps climbing\n(4D: about doubles every refinement, no finite worst case)', fontsize=10)
ax.legend(fontsize=8, loc='upper left')
ax = axs[1]
s = D['scan']; h = [q['hidden_spin']/D['Z_view'] for q in s]
ax.plot(h, [q['view'] for q in s], 'o-', color='#1f618d', label='growth seen in the 3D view')
ax.plot(h, [q['hidden'] for q in s], 's-', color='#e67e22', label='growth of the hidden 4th-direction part')
ax.plot(h, [q['total'] for q in s], 'D-', color='k', label='total growth (4D)')
ax.axhline(0, color='k', lw=0.6); ax.axvline(D['hidden_spin_to_cancel']/D['Z_view'], color='0.6', ls=':')
ax.text(D['hidden_spin_to_cancel']/D['Z_view']*1.05, 0.55, 'hidden spin = 4% of visible:\n3D view stops growing', fontsize=7.5)
ax.set_xlabel('hidden spin (varies along the 4th direction) ÷ visible spin')
ax.set_ylabel('spin growth rate')
ax.set_title('Can hidden 4th-direction motion stop the 3D view from growing?\nYes, but the total grows faster: spin is moved out of sight, not limited', fontsize=10)
ax.legend(fontsize=8)
fig.text(0.01, 0.005, 'DRAFT · start: the fastest-growing 3D flow (colliding rings), 32³, viscosity 0.05; hidden motion = the strongest-draining pattern (exact eigenproblem)', fontsize=7, color='0.4')
plt.tight_layout(rect=(0, 0.03, 1, 1)); plt.savefig('figures/fig13-four-dimensions.png', dpi=130); plt.close()

# 4D winner pictures
u = np.load('runs/feedd/pic_d4_n16_nu0.05_s0_u.npz')['u']; N = u.shape[1]; M = 48
def up(f):  # Fourier interpolation to M^4 would be heavy; interpolate 2D slices instead
    return f
k = np.fft.fftfreq(N, 1/N); K = np.meshgrid(k, k, k, k, indexing='ij')
uh = np.fft.fftn(u, axes=(1, 2, 3, 4))
A = np.array([[np.real(np.fft.ifftn(1j*K[j]*uh[i])) for j in range(4)] for i in range(4)])  # A[i,j] = d_j u_i
W = np.sqrt(sum((A[i, j] - A[j, i])**2 for i in range(4) for j in range(i+1, 4)))  # spin size (4D spin is a 2-form)
def interp2(sl):
    F = np.fft.fft2(sl); P = np.zeros((M, M), complex); h = N//2
    P[:h, :h] = F[:h, :h]; P[:h, -h:] = F[:h, -h:]; P[-h:, :h] = F[-h:, :h]; P[-h:, -h:] = F[-h:, -h:]
    return np.real(np.fft.ifft2(P))*(M/N)**2
i0 = np.unravel_index(np.argmax(W), W.shape)
names = 'xyzw'
fig, axs = plt.subplots(2, 4, figsize=(13, 6.6))
vmax = W.max()
for c, wi in enumerate([(i0[3] + s) % N for s in [0, N//8, N//4, N//2]]):
    ax = axs[0, c]; ax.imshow(interp2(W[:, :, i0[2], wi]).T, origin='lower', cmap='magma', vmin=0, vmax=vmax, extent=(0, 2*np.pi, 0, 2*np.pi))
    ax.set_title(f'x–y slice, w = {2*np.pi*(wi - i0[3])/N:+.2f} from the peak', fontsize=9); ax.set_xticks([]); ax.set_yticks([])
for c, (a, b) in enumerate([(0, 3), (1, 3), (2, 3)]):
    idx = list(i0); sl = [slice(None) if q in (a, b) else idx[q] for q in range(4)]
    ax = axs[1, c]; ax.imshow(interp2(W[tuple(sl)]).T, origin='lower', cmap='magma', vmin=0, vmax=vmax, extent=(0, 2*np.pi, 0, 2*np.pi))
    ax.set_title(f'{names[a]}–w slice through the peak', fontsize=9); ax.set_xlabel(names[a]); ax.set_ylabel('w (4th direction)'); ax.set_xticks([]); ax.set_yticks([])
# 3D view: w-average of the velocity, then its spin
ub = u.mean(axis=4)[:3]; kk = np.meshgrid(k, k, k, indexing='ij'); ubh = np.fft.fftn(ub, axes=(1, 2, 3))
g = lambda i, j: np.real(np.fft.ifftn(1j*kk[j]*ubh[i]))
Wv = np.sqrt((g(2, 1) - g(1, 2))**2 + (g(0, 2) - g(2, 0))**2 + (g(1, 0) - g(0, 1))**2)
ax = axs[1, 3]; ax.imshow(interp2(Wv[:, :, i0[2]]).T, origin='lower', cmap='magma', vmin=0, vmax=vmax, extent=(0, 2*np.pi, 0, 2*np.pi))
ax.set_title('3D view (averaged over w), x–y slice', fontsize=9); ax.set_xticks([]); ax.set_yticks([])
fig.suptitle('DRAFT · the fastest-growing 4D flow (16⁴, viscosity 0.05): spin size on slices, same colour scale\n'
             'It changes along w as much as along x, y, z, so it is a genuinely 4D flow, and its 3D view keeps only a faint trace', fontsize=10)
plt.tight_layout(rect=(0, 0, 1, 0.93)); plt.savefig('figures/fig14-fastest-4d-flow.png', dpi=120); plt.close()
print('ok', best[(4, 32, 0.05)]['R'], Wv.max()/vmax)
