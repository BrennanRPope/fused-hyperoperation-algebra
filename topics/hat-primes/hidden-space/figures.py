import json, numpy as np, matplotlib
matplotlib.use('Agg'); import matplotlib.pyplot as plt
from math import gcd, pi, log, exp
import sympy as sp
G, B, O, R, K = '#2e8b57', '#1f77b4', '#ff7f0e', '#d62728', '#444444'
ad = json.load(open('additive.json')); mu = json.load(open('multiplicative.json')); cv = np.load('mult_curves.npz')
plt.rcParams.update({'font.size': 10})

# Fig 1: the hidden space seen through 2, 3, 5, 7 (residues mod 210 on a 15 x 14 CRT grid)
P = np.array(list(sp.primerange(11, 2 * 10**6)))
cnt = np.bincount(P % 210, minlength=210)
grid = np.zeros((14, 15))
for r in range(210): grid[r % 14, r % 15] = cnt[r]
fig, ax = plt.subplots(1, 2, figsize=(12, 5.2), gridspec_kw=dict(width_ratios=[1.15, 1]))
im = ax[0].imshow(grid, cmap='Greens', origin='lower')
for r in range(210):
    if gcd(r, 210) == 1: ax[0].text(r % 15, r % 14, str(r), ha='center', va='center', fontsize=6.3)
ax[0].set_xlabel('n mod 15  (that is, n mod 3 and n mod 5)'); ax[0].set_ylabel('n mod 14  (n mod 2 and n mod 7)')
ax[0].set_title('A finite piece of the hidden space: every n has one cell (mod 210).\nPrimes up to 2 million (shade) fall only in the 48 numbered cells,\nthe window "not divisible by 2, 3, 5 or 7", and fill them evenly', fontsize=9.5)
plt.colorbar(im, ax=ax[0], shrink=0.75, label='primes in the cell')
Ps = [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47]
w = np.cumprod([1 - 1 / p for p in Ps])
ax[1].bar(range(len(Ps)), w, color=G)
ax[1].set_xticks(range(len(Ps))); ax[1].set_xticklabels(Ps)
ax[1].set_xlabel('add the prime p to the hidden space'); ax[1].set_ylabel('share of the space still open to primes')
ax[1].set_title('The prime window shrinks with every prime added:\nshare = product of (1 - 1/p), which goes to 0 (Mertens).\nThe primes are a window of size zero', fontsize=9.5)
for i, v in enumerate(w[:6]): ax[1].text(i, v + 0.01, f'{v:.3f}', ha='center', fontsize=7.5)
plt.tight_layout(); plt.savefig('figures/fig1-hidden-space-of-all-primes.png', dpi=130); plt.close()

# Fig 2: the window overshoots the primes by 2e^-gamma (the sieve's parity gap)
rows = ad['mertens']['rows']; x = [r['x'] for r in rows]
fig, ax = plt.subplots(1, 2, figsize=(12, 4.6))
ax[0].loglog(x, [r['window'] for r in rows], 'o-', color=O, label='window size: product of (1 - 1/p), p up to the square root of x')
ax[0].loglog(x, [r['density'] for r in rows], 's-', color=G, label='true share of primes up to x')
ax[0].loglog(x, [1 / log(v) for v in x], ':', color=K, label='1 / log x')
ax[0].set_xlabel('x'); ax[0].legend(fontsize=8); ax[0].set_title('The window cut by all primes up to the square root of x\nnearly matches the primes, but stays above them', fontsize=9.5)
ax[1].semilogx(x, [r['ratio_loglaw'] for r in rows], 'o-', color=O, label='window size times log x')
ax[1].axhline(2 * exp(-0.5772156649), color=R, ls='--', label='2 e^(-gamma) = 1.1229 (Mertens 1874)')
ax[1].axhline(1, color=G, ls=':', label='what the primes need: 1')
ax[1].set_ylim(0.95, 1.16); ax[1].legend(fontsize=8); ax[1].set_xlabel('x')
ax[1].set_title('The gap never closes: the window is 12% too big.\nA window in the hidden space can\'t tell primes from numbers\nwith two big factors (the sieve\'s parity problem, Selberg)', fontsize=9.5)
plt.tight_layout(); plt.savefig('figures/fig2-window-overshoots-by-2e-gamma.png', dpi=130); plt.close()

# Fig 3: the spectrum the hidden space predicts, against the primes
fig, ax = plt.subplots(1, 3, figsize=(15, 4.6))
for q in range(1, 31):
    m = int(sp.mobius(q))
    if m == 0: continue
    h = (m / int(sp.totient(q)))**2
    for a in range(0, q):
        if gcd(a, q) == 1: ax[0].vlines(a / q, 0, h, color=B, lw=1.2 if q < 8 else 0.7)
meas = {b['q']: b['re']**2 + b['im']**2 for b in ad['bragg']}
for q in [1, 2, 3, 5, 6, 7, 10, 11, 13, 14, 15]:
    ax[0].plot(1 / q if q > 1 else 0, meas[q], 'o', color=O, ms=4)
ax[0].set_yscale('log'); ax[0].set_ylim(1e-4, 1.5); ax[0].set_xlabel('frequency (fraction of a full turn)')
ax[0].set_title('Bright spots at every fraction a/q, brightness (mu(q)/phi(q))^2\n(blue, q up to 30). Orange: measured on the log-weighted primes\nto 10^7. Proven: prime number theorem in progressions', fontsize=9)
pr = ad['pairs']; hs = [p['h'] for p in pr]
ax[1].bar(hs, [p['measured'] for p in pr], color=G, label='measured, primes to 10^7')
ax[1].plot(hs, [p['hl'] for p in pr], '_', color=R, ms=9, mew=2, label='Hardy-Littlewood singular series')
ax[1].set_xlabel('gap h'); ax[1].legend(fontsize=8)
ax[1].set_title('Pairs of primes h apart: the spots add up to the\nHardy-Littlewood prediction (a conjecture, open even for h = 2)', fontsize=9)
tot = ad['bragg_total']
ax[2].semilogx([t['Q'] for t in tot], [t['bragg_total'] for t in tot], 'o-', color=B, label='total brightness of spots with q up to Q')
QQ = np.logspace(1, 7, 50); ax[2].semilogx(QQ, np.log(QQ) + 1.3325, ':', color=K, label='log Q + 1.3325 (known asymptotic)')
lm = ad['lambda_sq_mean']
Qc = np.exp(lm['mean'] - 1.3325); ax[2].axvline(Qc, color=O, ls='--', lw=1); ax[2].axhline(lm['mean'], color=O, ls='--', lw=1)
ax[2].text(Qc * 1.3, 4.5, f'q up to {Qc:,.0f}\n(about N / 10)', color=O, fontsize=8)
ax[2].plot([lm['N']], [lm['mean']], 's', color=O, ms=8, label='all the brightness there is (primes to 10^7)')
ax[2].legend(fontsize=8); ax[2].set_xlabel('Q')
ax[2].set_title('Spots with q up to about N/10 already hold all the light\nthe primes up to N = 10^7 have: at each scale the spectrum\nlooks all spots (Torquato et al. 2018-19; conditional on H-L)', fontsize=9)
plt.tight_layout(); plt.savefig('figures/fig3-spectrum-from-the-hidden-space.png', dpi=130); plt.close()

# Fig 4: the zeta zeros as a slice of the prime torus (2-torus, t from 25 to 56)
fig, ax = plt.subplots(1, 2, figsize=(12, 5.4), gridspec_kw=dict(width_ratios=[1, 1.25]))
a = np.linspace(0, 2 * pi, 400); A_, B_ = np.meshgrid(a, a, indexing='ij')
F = np.cos(A_) + 2**-0.5 * np.cos(A_ - B_)
ax[0].contour(a, a, F.T, levels=[0], colors=R, linewidths=2)
A, Bc = cv['A'], cv['B']
br = np.nonzero((np.abs(np.diff(A)) > 3) | (np.abs(np.diff(Bc)) > 3))[0] + 1
for seg_a, seg_b in zip(np.split(A, br), np.split(Bc, br)): ax[0].plot(seg_a, seg_b, color=B, lw=0.9)
t2 = cv['t2']; th = lambda t: t / 2 * np.log(t / (2 * pi)) - t / 2 - pi / 8 + 1 / (48 * t)
for z in cv['z2']: ax[0].plot(th(z) % (2 * pi), (z * log(2)) % (2 * pi), 'o', color=O, ms=7, zorder=5)
ax[0].set_xlabel('first angle: theta(t)'); ax[0].set_ylabel('second angle: t log 2  (the prime 2)')
ax[0].set_title('Hidden torus for heights 25 to 56: red is where\ncos A + cos(A - B)/sqrt2 = 0; the blue path winds as t grows;\neach crossing (orange) is a zero of the slice', fontsize=9)
ax[0].set_xlim(0, 2 * pi); ax[0].set_ylim(0, 2 * pi)
f2 = np.cos(th(t2)) + 2**-0.5 * np.cos(th(t2) - t2 * log(2))
ax[1].plot(t2, f2, color=B, label='the slice along the path')
ax[1].axhline(0, color=K, lw=0.6)
ax[1].plot(cv['z2'], 0 * cv['z2'], 'o', color=O, label='slice zeros')
ax[1].plot(cv['true2'], 0 * cv['true2'], 'x', color=R, ms=9, mew=2, label='true zeta zeros')
ax[1].legend(fontsize=8); ax[1].set_xlabel('height t')
ax[1].set_title('Just two terms (n = 1, 2): 8 of the 9 zeros here land within a tenth\nof a gap; the 9th sits at the edge, where the n = 3 term switches on', fontsize=9.5)
plt.tight_layout(); plt.savefig('figures/fig4-zeros-as-a-slice-of-the-prime-torus.png', dpi=130); plt.close()

# Fig 5: how many primes (torus angles) each height needs, and how good the slice is
fig, ax = plt.subplots(1, 2, figsize=(12, 4.6))
T = np.logspace(1.2, 12, 300)
ax[0].semilogx(T, [sp.primepi(int(np.sqrt(t / (2 * pi)))) for t in T], color=B, drawstyle='steps-post')
for d in mu['dimension'][2:8]:
    ax[0].annotate(f"{d['primes']}", (d['t'], d['primes']), textcoords='offset points', xytext=(-6, 6), fontsize=8)
ax[0].set_xlabel('height t'); ax[0].set_ylabel('primes in play (torus angles, plus theta)')
ax[0].set_title('The torus the zeros need grows with height: it has one angle\nfor each prime up to sqrt(t / 2 pi). Never finite, so the zeros\ncan\'t be a quasicrystal of a fixed torus (Kurasov-Sarnak type)', fontsize=9)
ax[1].hist(cv['err'], bins=np.linspace(0, 0.6, 61), color=G)
m = mu['main_sum']
ax[1].set_xlabel('distance from a true zero to the nearest slice zero, in gaps')
ax[1].set_title(f"First 1,000 zeros against the slice (Riemann-Siegel main sum):\nmedian {m['median_err_gaps']:.3f} of a gap; {100*m['windows_exactly_one']:.1f}% get exactly one slice zero;\n{m['windows_zero']} slice zeros missing, put back by the correction terms", fontsize=9)
plt.tight_layout(); plt.savefig('figures/fig5-torus-size-and-slice-quality.png', dpi=130); plt.close()

# Fig 6: zeta's slice is special: shifted slices of the same torus lose zeros
sh = mu['shifts']
fig, ax = plt.subplots(figsize=(8.5, 4.8))
Ns = [s['N'] for s in sh]
ax.plot(Ns, [s['random_mean_lost'] for s in sh], 'o-', color=R, label='slice with random angles per prime (mean of 40)')
ax.plot(Ns, [s['cos_theta_zeros'] - s['zeta_slice_zeros'] for s in sh], 's-', color=G, label="zeta's own slice (every angle 0)")
ax.axhline(0, color=K, lw=0.6)
ax.set_xlabel('terms in the slice, N (heights 2 pi N^2 to 2 pi (N+1)^2)'); ax.set_ylabel('zeros lost (against the count theta / pi)')
ax.set_title("Same torus, other slices: they lose more and more zeros (they go complex).\nZeta's slice is the one tied to its functional equation, and keeps them.\nWith one or two terms nothing can be lost (the first term dominates)", fontsize=9.5)
ax.legend(fontsize=8.5)
plt.tight_layout(); plt.savefig('figures/fig6-zetas-slice-is-special.png', dpi=130); plt.close()

# Fig 7: the two hidden spaces and where they meet
fig, ax = plt.subplots(figsize=(12, 6)); ax.axis('off')
def box(x, y, w, h, txt, c):
    ax.add_patch(plt.Rectangle((x, y), w, h, fc=c, alpha=0.18, ec=c, lw=2)); ax.text(x + w / 2, y + h / 2, txt, ha='center', va='center', fontsize=9.2)
box(0.02, 0.45, 0.42, 0.5, 'ADDING: the profinite integers\nZ^ = one circle of residues for every prime power\n\nSees: arithmetic progressions, twin-prime pairs,\nsieves. Spectrum: spots at every fraction a/q.\nSquare-free numbers: a true window (pure spots, proven).\nPrimes: a window of size zero; the cut-out\nwindow is off by 2e^-gamma (parity).\nPeaks proven; "all spots" = Hardy-Littlewood (open)', B)
box(0.56, 0.45, 0.42, 0.5, 'MULTIPLYING: the prime torus prod S^1\n= one angle (p^-it) for every prime (Bohr)\n\nSees: zeta and its zeros. Spectrum: lines at log p^m\n(explicit formula). Zeros = where the path\n(theta(t), t log 2, t log 3, ...) cuts a zero surface.\nAngles needed grow like the primes up to sqrt(t/2pi).\nZeta\'s slice keeps all zeros; other slices lose them', O)
box(0.25, 0.04, 0.5, 0.3, 'WHERE THEY MEET: the adeles (Tate 1950)\nand Connes\' adele class space (1999).\nThe join gives the functional equation\n(repo, Riemann round 10), not the zeros:\nthat step is the open problem', R)
ax.annotate('', xy=(0.56, 0.7), xytext=(0.44, 0.7), arrowprops=dict(arrowstyle='<->', lw=2, color=K))
ax.text(0.5, 0.74, 'explicit formula\n(Riemann, von Mangoldt)', ha='center', fontsize=8.5)
ax.annotate('', xy=(0.42, 0.34), xytext=(0.25, 0.45), arrowprops=dict(arrowstyle='->', lw=1.5, color=K))
ax.annotate('', xy=(0.58, 0.34), xytext=(0.75, 0.45), arrowprops=dict(arrowstyle='->', lw=1.5, color=K))
ax.set_title('A hidden space built from all the primes comes in two halves, one for adding and one for multiplying', fontsize=11)
plt.savefig('figures/fig7-two-halves-of-the-hidden-space.png', dpi=130, bbox_inches='tight'); plt.close()
print('ok')
