import json, numpy as np, matplotlib
matplotlib.use('Agg'); import matplotlib.pyplot as plt
from math import pi
G, B, O, R, K = '#2e8b57', '#1f77b4', '#ff7f0e', '#d62728', '#444444'
A = json.load(open('apart.json')); a = np.load('apart.npz')
C = json.load(open('closing.json')); c = np.load('closing.npz')
P = json.load(open('polys.json')); p = np.load('polys.npz')
# Fig 10: pushing the gaps apart
fig, ax = plt.subplots(1, 3, figsize=(16, 4.8))
z0, zb = a['z0_8'], a['zb_8']
lo, hi = z0[20], z0[45]
for zz, y, col, lab in ((z0, 1, B, "zeta's own slice"), (zb, 0, O, 'angles turned to push gaps apart')):
    m = (zz >= lo - 1) & (zz <= hi + 1); ax[0].plot(zz[m], np.full(m.sum(), y), '|', ms=30, mew=2, color=col)
    ax[0].text(lo - 1, y + 0.28, lab, fontsize=9, color=col)
g0 = np.diff(z0) * np.log(((z0[1:] + z0[:-1]) / 2) / (2 * pi)) / (2 * pi); i0 = np.argmin(g0)
ax[0].annotate(f'closest pair, {g0[i0]:.2f} mean gaps', ((z0[i0] + z0[i0 + 1]) / 2, 1.05), xytext=((z0[i0] + z0[i0 + 1]) / 2 - 3, 1.55),
               fontsize=8, arrowprops=dict(arrowstyle='->'))
ax[0].set_ylim(-0.6, 2.0); ax[0].set_yticks([]); ax[0].set_xlabel('height t'); ax[0].set_xlim(lo - 1.5, hi + 1)
ax[0].set_title(f"Zeros near t = {lo:.0f} (N = 8). A turn of {A['8']['best_distance_from_zeta']:.2f} radians\nwidens the closest pair; every zero stays on the line", fontsize=9.5)
Ns = [5, 6, 8, 10, 12]; x = np.arange(len(Ns))
ax[1].bar(x - 0.2, [A[str(n)]['zeta_min_gap'] for n in Ns], 0.4, color=B, label="zeta's own slice")
ax[1].bar(x + 0.2, [A[str(n)]['best_min_gap'] for n in Ns], 0.4, color=O, label='best turn found')
ax[1].set_xticks(x); ax[1].set_xticklabels([f'N = {n}' for n in Ns]); ax[1].set_ylabel('smallest gap in the block (mean gaps)')
ax[1].legend(fontsize=8); ax[1].set_title('Smallest gap, zeta vs the best turn: it can be\nraised 1.1 to 2 times, but only by leaving zeta', fontsize=9.5)
ax[2].bar(x, [100 * A[str(n)]['random_keep_all_zeros'] for n in Ns], 0.5, color=K)
for i, n in enumerate(Ns): ax[2].text(i, 100 * A[str(n)]['random_keep_all_zeros'] + 0.1, f"{100 * A[str(n)]['random_keep_all_zeros']:.1f}%", ha='center', fontsize=8)
ax[2].set_xticks(x); ax[2].set_xticklabels([f'N = {n}' for n in Ns]); ax[2].set_ylabel('random angle settings that keep every zero (%)')
ax[2].set_title('Almost no random setting keeps every zero on the line:\nzeta sits in a small pocket, near but not at its best spot', fontsize=9.5)
plt.tight_layout(); plt.savefig('../figures/fig10-pushing-zero-gaps-apart.png', dpi=130); plt.close()
# Fig 11: the fastest closing turn as a detector
fig, ax = plt.subplots(1, 2, figsize=(13.5, 4.8))
ax[0].semilogy(c['tB'], c['mbB'], 'o-', color=B, ms=4)
top = C['B']['smallest'][0]
ax[0].annotate(f"Lehmer's pair, t = 7005.06 / 7005.10\nturn of {top['m_best']:.4f} rad closes it", (top['t0'], top['m_best']),
               xytext=(top['t0'] + 2.5, 0.004), fontsize=9, arrowprops=dict(arrowstyle='->'))
ax[0].set_xlabel('height t (zeta: main sum plus first correction, N = 33)'); ax[0].set_ylabel('smallest turn that closes this gap (radians)')
ax[0].set_title('The fastest closing turn as a function of height:\none gap stands out by a factor of 60', fontsize=9.5)
ins = ax[0].inset_axes([0.06, 0.08, 0.3, 0.3]); m = (c['tfB'] > 7004.6) & (c['tfB'] < 7005.6)
ins.plot(c['tfB'][m], c['fB'][m], color=K, lw=1); ins.axhline(0, color=R, lw=0.6); ins.set_title('Z(t) there', fontsize=7); ins.tick_params(labelsize=6)
ax[1].loglog(c['s'], c['mp'], '.', color=K, alpha=0.25, ms=3, label='best single prime (section 4)')
ax[1].loglog(c['s'], c['mb'], '.', color=B, alpha=0.5, ms=3, label='all primes, steepest direction')
ok = np.isfinite(c['ex']); ax[1].loglog(c['s'][c['sel']][ok], c['ex'][ok], 'o', mfc='none', color=R, ms=5, label='exact, by bisection')
xx = np.geomspace(0.1, 3, 50); ax[1].plot(xx, C['A']['C_median'] * xx**2, '--', color=O, label=f"{C['A']['C_median']:.2f} x gap squared (median)")
ax[1].set_xlabel('gap between the two zeros (mean gaps), heights 57 to 1414'); ax[1].set_ylabel('turn that closes it (radians)')
ax[1].legend(fontsize=8); ax[1].set_title('Closing turn vs gap: about 0.85 times the gap squared;\nturning all primes together beats one prime by 1.2 times', fontsize=9.5)
plt.tight_layout(); plt.savefig('../figures/fig11-fastest-closing-turn.png', dpi=130); plt.close()
# Fig 12: polynomials whose values are primes with small gaps
fig, ax = plt.subplots(1, 3, figsize=(16.5, 4.8), gridspec_kw=dict(width_ratios=[1.1, 1, 1.2]))
cs = np.arange(1, len(p['C1']) + 1); odd = cs % 2 == 1
ax[0].plot(cs[odd], p['C1'][odd], '.', color=K, alpha=0.15, ms=1.5)
res = P['results']
for cc, col in (('41', R), ('55661', O), ('46127', G)):
    lab = {'41': "c = 41 (Euler)", '55661': 'c = 55661 (best here)', '46127': 'c = 46127 (best for twin values)'}[cc]
    ax[0].plot(int(cc), p['C1'][int(cc) - 1], 'o', color=col, mec='k', ms=8, zorder=5, label=lab)
ax[0].axhline(1, color=B, lw=1, label='the integers themselves (C = 1)')
ax[0].set_xscale('log'); ax[0].set_xlabel('c in n^2 + n + c'); ax[0].set_ylabel('prime-density constant C (primes up to 2000)')
ax[0].legend(fontsize=8, loc='center left', bbox_to_anchor=(0.0, 0.45)); ax[0].set_title('How much denser the prime values are than\nrandom numbers of the same size, for c up to 100,000', fontsize=9.5)
for i, cc in enumerate(['1', '41', '55661']):
    ok = np.unpackbits(p[f'ok_{cc}'])[:200]
    ax[1].bar(np.nonzero(ok)[0], 0.8, bottom=2 - i, width=0.8, color=[B, R, O][i])
    ax[1].text(201, 2.4 - i, f"c = {cc}\nC = {res[cc]['C']:.2f}, gap {res[cc]['mean_gap_in_n']:.1f}", fontsize=8, va='center')
ax[1].set_xlim(0, 260); ax[1].set_yticks([]); ax[1].set_xlabel('n (first 200); a bar where n^2 + n + c is prime')
ax[1].set_title('Barcodes: bigger C, more bars, smaller gaps\n(mean gap in n over the first million)', fontsize=9.5)
names = sorted(res, key=lambda k: res[k]['C'])
pr = np.array([res[k]['predicted'] for k in names]); ac = np.array([res[k]['primes'] for k in names])
tp = np.array([res[k]['twins_predicted'] for k in names if res[k]['C_twin'] > 1e-6]); ta = np.array([res[k]['twins'] for k in names if res[k]['C_twin'] > 1e-6])
ax[2].loglog(pr, ac, 'o', color=B, label='prime values f(n)')
ax[2].loglog(tp, ta, 's', color=G, label='twin values f(n), f(n)+2')
xx = np.array([3e3, 4e5]); ax[2].plot(xx, xx, ':', color=K)
ax[2].set_xlabel('predicted by Bateman-Horn, n < 1,000,000'); ax[2].set_ylabel('counted')
ax[2].legend(fontsize=8); ax[2].set_title('The prediction holds to within 0.3% (primes)\nand 1.6% (twins) for all 20 polynomials', fontsize=9.5)
plt.tight_layout(); plt.savefig('../figures/fig12-prime-values-with-small-gaps.png', dpi=130); plt.close()
print('ok')
