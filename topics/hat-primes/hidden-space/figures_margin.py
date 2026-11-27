import json, numpy as np, matplotlib
matplotlib.use('Agg'); import matplotlib.pyplot as plt
from math import pi, log
G, B, O, R, K = '#2e8b57', '#1f77b4', '#ff7f0e', '#d62728', '#444444'
W = json.load(open('margin.json'))['windows']; L = json.load(open('margin_local.json')); d = np.load('margin_local.npz')
Rr = d['R']; Z = np.load('../../primes/zeta_zeros_1000.npy')
def theta(t): return t / 2 * np.log(t / (2 * pi)) - t / 2 - pi / 8 + 1 / (48 * t) + 7 / (5760 * t**3)
# Fig 8: margins per window and per prime
fig, ax = plt.subplots(1, 2, figsize=(13, 4.8))
cols = {2: R, 3: O, 5: G, 7: B, 11: '#9467bd', 13: '#8c564b'}
for w in W:
    for p, v in w['margin_per_prime'].items():
        ax[0].plot(w['N'], v, 'o', color=cols[int(p)], ms=5)
for p, c in cols.items(): ax[0].plot([], [], 'o', color=c, label=f'turn prime {p} only')
ax[0].plot([w['N'] for w in W], [w['margin_random_dir'] for w in W], 'k_', ms=14, mew=2, label='turn all angles (worst of 8 directions)')
ax[0].set_yscale('log'); ax[0].set_xlabel('terms N (heights 2 pi N^2 to 2 pi (N+1)^2)'); ax[0].set_ylabel('smallest turn that changes the zeros (radians)')
ax[0].annotate('t = 564.3: the slice just misses a close pair\nthat zeta has; a 0.00003 turn makes it', (9, 3e-5), xytext=(9.6, 4e-4), fontsize=8, arrowprops=dict(arrowstyle='->'))
ax[0].legend(fontsize=7, loc='lower left', ncol=2); ax[0].set_ylim(2e-6, 2)
ax[0].set_title('Whole windows: a turn of a few hundredths of a radian already\nchanges the zeros. Prime 2 is almost always the most sensitive', fontsize=9.5)
Ns = sorted(int(k) for k in L['median_margin_by_N'])
for n in Ns:
    v = Rr[Rr[:, 0] == n, 4]
    ax[1].plot(np.full(len(v), n) + np.random.default_rng(n).uniform(-0.25, 0.25, len(v)), v, '.', color=B, alpha=0.25, ms=3)
ax[1].plot(Ns, [L['median_margin_by_N'][str(n)] for n in Ns], 's-', color=O, label='median over the gaps')
ax[1].set_yscale('log'); ax[1].set_xlabel('terms N'); ax[1].set_ylabel('margin of one gap (radians, first order)')
ax[1].legend(fontsize=8)
ax[1].set_title('Gap by gap: the typical margin is about 1 radian and drifts down\nwith height (1.44 to 0.81); the smallest few set the window\'s margin', fontsize=9.5)
plt.tight_layout(); plt.savefig('figures/fig8-stability-margin.png', dpi=130); plt.close()
# Fig 9: what sets it
fig, ax = plt.subplots(1, 3, figsize=(15.5, 4.8))
ax[0].loglog(Rr[:, 5], Rr[:, 4], '.', color=B, alpha=0.4, ms=4)
ax[0].set_xlabel('gap between the two zeros (in mean gaps)')
from matplotlib.ticker import FixedLocator, NullLocator, FormatStrFormatter
ax[0].xaxis.set_major_locator(FixedLocator([0.2, 0.5, 1, 2])); ax[0].xaxis.set_minor_locator(NullLocator()); ax[0].xaxis.set_major_formatter(FormatStrFormatter('%g')); ax[0].set_ylabel('margin of that gap (radians)')
ax[0].set_title(f"Margin follows the gap (correlation of logs {L['corr_logmargin_loggap']:.2f}):\nclose pairs of zeros are the fragile spots", fontsize=9.5)
pick, tz = d['pick'], d['trueZ']
ax[1].loglog(2 * Rr[pick[:40], 2], tz[:40], 'o', color=R, ms=4, label='40 smallest margins')
ax[1].loglog(2 * Rr[pick[40:], 2], tz[40:], 'o', color=G, ms=4, label='40 random gaps')
xx = np.array([0.03, 10]); ax[1].plot(xx, xx, ':', color=K)
ax[1].set_xlabel('peak of the slice between the two zeros'); ax[1].set_ylabel('peak of true zeta (|Z|) there')
ax[1].legend(fontsize=8)
ax[1].set_title(f"The slice's weak peaks are zeta's weak peaks (corr {L['corr_slice_true_peak']:.2f}):\nmedian true peak {L['median_true_peak_smallest40']:.2f} at the fragile gaps, {L['median_true_peak_random40']:.2f} at random ones.\nThese are Lehmer-type near misses", fontsize=9)
N = 9; t = np.arange(563.95, 564.75, 0.0005); ns = np.arange(1, N + 1)
e2 = np.array([1 if n % 2 == 0 else 0 for n in ns]) + np.array([1 if n % 4 == 0 else 0 for n in ns]) + np.array([1 if n % 8 == 0 else 0 for n in ns])
for dd, c, lab in [(0, B, "zeta's slice (no turn)"), (-0.03, O, 'prime 2 turned by -0.03')]:
    f = 2 * (ns[:, None]**-0.5 * np.cos(theta(t)[None, :] - t[None, :] * np.log(ns)[:, None] + dd * e2[:, None])).sum(0)
    ax[2].plot(t, f, color=c, label=lab)
import mpmath as mp
ax[2].plot(t[::20], [float(mp.siegelz(x)) for x in t[::20]], color=K, lw=1, ls='--', label='true Z(t)')
tz9 = Z[(Z > 563.95) & (Z < 564.75)]; ax[2].plot(tz9, 0 * tz9, 'x', color=R, ms=9, mew=2, label='true zeros')
ax[2].axhline(0, color=K, lw=0.6); ax[2].set_ylim(-0.35, 0.25); ax[2].set_xlabel('height t'); ax[2].legend(fontsize=7.5)
ax[2].set_title('At t = 564.3 the slice touches zero without crossing;\nzeta has a close pair there (corrections restore it).\nThe slightest turn of prime 2 makes the pair', fontsize=9)
plt.tight_layout(); plt.savefig('figures/fig9-what-sets-the-margin.png', dpi=130); plt.close()
print('ok')
