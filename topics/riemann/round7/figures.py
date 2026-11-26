"""Pictures for Riemann round 7.  Colour key: blue = our reading, orange = known (literature),
grey = reference, red = where a road stops / the hypothesis fails."""
import json, math, numpy as np, matplotlib
import mpmath as mp
matplotlib.use('Agg'); import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, Rectangle
B, O, G_, R, GR = '#1f77b4', '#ff7f0e', '#888888', '#d62728', '#2ca02c'
plt.rcParams.update({'font.size': 10, 'axes.spines.top': False, 'axes.spines.right': False})
PG = json.load(open('prime_gas.json')); FI = json.load(open('fisher.json')); HP = json.load(open('hilbert_polya.json'))
RM = json.load(open('rmt.json')); TR = json.load(open('trace.json')); SS = json.load(open('solved_shapes.json'))

# ---------- fig0: atlas of solved shapes
rows = [
 ('Curves over a finite field\n(Weil 1948, Deligne 1974)', 'one prime scale q;\ncurvature dial, Δ ≤ 0', 'Frobenius acting on a\nfinite cohomology +\nHodge-index positivity', 'infinitely many unrelated\nscales log p; no surface\nSpec Z × Spec Z with positivity'),
 ('Ramanujan graphs\n(Ihara, Sunada)', 'curvature dial per eigenvalue:\n1 − λu + qu², Δ ≤ 0', 'finite symmetric matrix;\nRH = spectral gap |λ| ≤ 2√q', 'no finite matrix; and not\nevery graph is Ramanujan'),
 ('Hyperbolic surfaces\n(Selberg zeta)', 'periodic geodesics =\nprimes; lengths ≈ log p', 'self-adjoint Laplacian\n(Hilbert–Pólya realised)\n+ trace formula', 'no surface whose geodesics\nare the primes; RH there needs\nno eigenvalue below 1/4'),
 ('Lee–Yang circle theorem\n(Ising ferromagnet)', 'odds / fugacity dial t', 'finite multi-affine polynomial;\nferromagnetic couplings\n(Asano contraction)', 'ζ is not a finite polynomial;\nfugacity zeros are just z = p^s'),
 ('Heat flow (de Bruijn,\nNewman, Rodgers–Tao)', 'our rewind\n(backward heat flow)', 'forward heat makes all\nzeros real after time Λ;\n0 ≤ Λ ≤ 0.22 (Polymath 15)', 'RH = Λ = 0 exactly:\nno room to spare'),
 ('Local RH: Ramanujan–\nPetersson (Deligne)', 'curvature dial per prime,\nΔ_p ≤ 0 (s+d = Hasse)', 'Weil conjectures for\nKuga–Sato varieties', 'ζ’s own factor 1 − u has\nd = 0: the local test is empty'),
]
fig, ax = plt.subplots(figsize=(13.5, 9.2)); ax.axis('off'); ax.set_xlim(0, 4); ax.set_ylim(-0.2, len(rows)+0.9)
heads = ['Where RH is a theorem', 'Which dial it is', 'What made it provable', 'What ζ lacks']
cols = [O, B, GR, R]
for j, h in enumerate(heads):
    ax.text(j+0.5, len(rows)+0.45, h, ha='center', va='center', fontsize=11.5, weight='bold', color=cols[j])
for i, r in enumerate(rows):
    y = len(rows)-1-i
    for j, txt in enumerate(r):
        ax.add_patch(FancyBboxPatch((j+0.03, y+0.05), 0.94, 0.88, boxstyle='round,pad=0.01', fc='white', ec=cols[j], lw=1.3))
        ax.text(j+0.5, y+0.49, txt, ha='center', va='center', fontsize=8.6)
ax.set_title('Six shapes where a Riemann-type hypothesis is proved, read through our dials', fontsize=13)
plt.savefig('figures/fig0-solved-shapes-atlas.png', dpi=140, bbox_inches='tight'); plt.close()

# ---------- fig1: the order dial as statistics
fig, axs = plt.subplots(1, 2, figsize=(12.5, 4.9))
ax = axs[0]
q = np.linspace(-0.98, 0.95, 400)
ax.plot(q, 1/(1-q), color=B, lw=2, label='power γ = 1/(1−q): Z ~ (β−1)^(−γ)')
ex = PG['exponent']
ax.plot([float(k) for k in ex], [v['slope'] for v in ex.values()], 'o', color=B, ms=6, label='measured slope (β−1 = 10⁻⁶)')
ax.plot([-1], [0], 's', color=R, ms=8); ax.annotate('q = −1: (1+P(β))/(1−P(2β))\nlog divergence only', (-1, 0), (-0.95, 3.6), fontsize=8.5, color=R, arrowprops=dict(arrowstyle='->', color=R))
ax.plot([0], [1], 'o', color=O, ms=9, mfc='none', mew=2); ax.annotate('q = 0: Bose gas, ζ(β)', (0, 1), (0.12, 0.25), fontsize=8.5, color=O, arrowprops=dict(arrowstyle='->', color=O))
ax.axvline(1.0, color=R, lw=1, ls='--'); ax.text(-0.62, 2.3, 'q = 1 (dashed): Boltzmann gas 1/(1−P(β));\nthe wall jumps to Hagedorn β = 1.39943', fontsize=8.5, color=R, ha='left')
ax.set_xlabel('order dial q'); ax.set_ylabel('critical power at β = 1'); ax.set_ylim(-0.4, 8); ax.legend(fontsize=8.5, loc='upper left', bbox_to_anchor=(0.0, 0.82))
ax.set_title('Statistics dial: the wall stays at β = 1, its strength is 1/(1−q)', fontsize=10)
ax = axs[1]
th = PG['thermo']; bet = np.array(th['betas'])
for qq, c in zip(['0.0', '0.3', '0.6'], ['#08306b', '#2171b5', '#6baed6']):
    ax.plot(bet, th[qq], color=c, lw=2, label=f'order dial q = {qq}')
ax.plot(bet, th['fermi'], color=G_, lw=1.5, ls='--', label='Fermi: ζ(β)/ζ(2β)')
from sympy import mobius as _mu
mp.mp.dps = 15
def _P(b): return sum(_mu(n)/mp.mpf(n)*mp.log(mp.zeta(n*b)) for n in range(1, 60) if _mu(n) != 0)
_bb = PG['hagedorn'] + np.geomspace(2e-4, 1.6, 60)
ax.plot(_bb, [float(-mp.log(1-_P(b))) for b in _bb], color=R, lw=1.5, label='Boltzmann (q = 1)')
ax.axvline(1, color='k', lw=0.6, ls=':'); ax.axvline(PG['hagedorn'], color=R, lw=0.6, ls=':')
ax.text(1.02, 13, 'β = 1', fontsize=8.5); ax.text(PG['hagedorn']+0.02, 13, 'Hagedorn 1.39943', fontsize=8.5, color=R)
ax.set_xlim(0.95, 3.0); ax.set_ylim(0, 14); ax.set_xlabel('inverse temperature β'); ax.set_ylabel('log Z')
ax.legend(fontsize=8.5, loc='upper right'); ax.set_title('Free energy of the prime gas for each statistics', fontsize=10)
plt.tight_layout(); plt.savefig('figures/fig1-prime-gas-statistics.png', dpi=150, bbox_inches='tight'); plt.close()

# ---------- fig2: Fisher zeros
fig, ax = plt.subplots(figsize=(7.6, 7.2))
rz = FI['riemann']
ax.plot([0.5]*len(rz), rz, 'x', color='k', ms=8, mew=2, label='Riemann zeros: log branch points of Z_q (all q ≠ 0)')
for qq, c, m in [('-0.5', '#9e9ac8', 'o'), ('-0.8', '#54278f', 's')]:
    z = np.array([r[:2] for r in FI[qq]['zeros']])
    ax.plot(z[:, 0], z[:, 1], m, color=c, ms=7, label=f'Fisher zeros, q = {qq}')
ax.axvline(0.5, color='k', lw=0.8)
ax.fill_betweenx([0, 51], 0.505, 1.605, color=B, alpha=0.06)
ax.text(1.12, 3.5, 'searched: q = 0.2, 0.5, 0.8\nno Fisher zeros here', fontsize=9, color=B)
ax.set_xlim(0.4, 1.65); ax.set_ylim(0, 51); ax.set_xlabel('Re β'); ax.set_ylabel('Im β')
ax.legend(fontsize=8.5, loc='upper right'); ax.set_title('Zeros of the order-dial prime gas Z_q(β) for Re β > 1/2', fontsize=10.5)
plt.savefig('figures/fig2-fisher-zeros.png', dpi=150, bbox_inches='tight'); plt.close()

# ---------- fig3: Hilbert-Polya Euler blindness
fig, ax = plt.subplots(figsize=(9.5, 4.8))
import mpmath as mp
zz = [float(mp.im(mp.zetazero(n))) for n in range(1, 120)]
E = [-2*g for g in zz if 2*g < 380]
ax.plot(E, [0]*len(E), '|', color=O, ms=16, mew=1.4, label='ζ: E = i(2ρ − 1), all real')
dh = np.array([d['E'] for d in HP['dh']])
ax.plot(dh[:, 0], dh[:, 1], 'o', color=R, ms=8, label='Davenport–Heilbronn: same recipe, E not real')
e5 = HP['x2+5y2']['E']; ax.plot([e5[0]], [e5[1]], 's', color='#a50f15', ms=8, label='x² + 5y² (round 6 zero)')
for d in HP['dh']:
    ax.annotate(f"ρ = {d['zero'][0]:.4f} + {d['zero'][1]:.2f}i", (d['E'][0], d['E'][1]), (d['E'][0]+4, d['E'][1]+0.05), fontsize=7.8, color=R)
ax.set_xlabel('Re E'); ax.set_ylabel('Im E'); ax.set_ylim(-0.15, 1.05); ax.set_xlim(-380, 0)
ax.legend(fontsize=8.5, loc='upper left'); ax.set_title('Bender–Brody–Müller "eigenvalues": the construction never uses the Euler product', fontsize=10.5)
plt.savefig('figures/fig3-hilbert-polya-euler-blind.png', dpi=150, bbox_inches='tight'); plt.close()

# ---------- fig4: spacing statistics
sg = np.array(RM['sgrid'])
fig, axs = plt.subplots(2, 2, figsize=(12, 8.4))
panels = [('zeta', 'ζ (one Euler product), heights 10000–11600', False), ('x2+y2', 'x² + y² = 4ζ·L(χ₋₄) (two Euler products)', False),
          ('x2+5y2', 'x² + 5y² (no Euler product), zeros on the line', True), ('DH', 'Davenport–Heilbronn (no Euler product), on the line', True)]
for ax, (k, title, resc) in zip(axs.ravel(), panels):
    sp = np.array(RM[k]['spacings']); sp = sp/sp.mean() if resc else sp
    ax.hist(sp, bins=40, range=(0, 4), density=True, color='#c6dbef', ec='white')
    ax.plot(sg, RM['p_gue'], color=O, lw=2, label='GUE')
    ax.plot(sg, RM['p_sup'], color=GR, lw=2, ls='--', label='two independent GUE')
    ax.plot(sg, RM['p_poi'], color=G_, lw=1.2, ls=':', label='Poisson')
    ksr = RM[k]['ks_rescaled'] if resc else RM[k]['ks']
    ax.text(2.25, 0.62, f"{RM[k]['n']} zeros\nKS distance:\n GUE {ksr[0]:.3f}\n two GUE {ksr[1]:.3f}", fontsize=8.5)
    ax.set_title(title, fontsize=9.5); ax.set_xlim(0, 4); ax.set_ylim(0, 1.15); ax.set_xlabel('unfolded spacing' + (' (rescaled to mean 1)' if resc else ''))
    ax.legend(fontsize=8, loc='upper right')
fig.suptitle('Nearest-neighbour spacings: one Euler product → GUE; two → two GUE superposed; none → neither', fontsize=11.5)
plt.tight_layout(); plt.savefig('figures/fig4-spacing-statistics.png', dpi=150, bbox_inches='tight'); plt.close()

# ---------- fig5: off-line fraction
fig, ax = plt.subplots(figsize=(8.5, 4.2))
names = [('zeta', 'ζ'), ('L-4', 'L(s, χ₋₄)'), ('x2+y2', 'x² + y²'), ('DH', 'Davenport–Heilbronn'), ('x2+5y2', 'x² + 5y²')]
tot = [RM[k]['n'] + max(0, RM[k]['deficit']) for k, _ in names]
off = [max(0, RM[k]['deficit']) for k, _ in names]
yy = np.arange(len(names))
ax.barh(yy, [100*o/t for o, t in zip(off, tot)], color=[O, O, O, R, R])
for y, (k, _), o, t in zip(yy, names, off, tot):
    T = RM[k]['T']; ax.text(100*o/t + 0.6, y, f"{o:.0f} of {t:.0f} off the line (heights {T[0]}–{T[1]})", va='center', fontsize=8.5)
ax.set_yticks(yy); ax.set_yticklabels([n for _, n in names]); ax.set_xlim(0, 60); ax.set_xlabel('% of zeros off the critical line')
ax.set_title('Zeros off the line, counted as phase-function count minus sign changes', fontsize=10.5)
plt.savefig('figures/fig5-off-line-fraction.png', dpi=150, bbox_inches='tight'); plt.close()

# ---------- fig6: trace formula both ways
fig, axs = plt.subplots(2, 1, figsize=(11, 7.4))
zp = TR['zeros_from_primes']; Ee = np.array(zp['E'])
axs[0].plot(Ee, np.array(zp['smooth'])+np.array(zp['osc']), color=B, lw=1.3, label='smooth term + orbit sum over primes < 3000')
for g in TR['rz']: axs[0].axvline(g, color=O, lw=0.9, ls='--')
axs[0].set_xlim(5, 60); axs[0].set_xlabel('E'); axs[0].set_ylabel('zero density'); axs[0].legend(fontsize=8.5, loc='upper left')
axs[0].set_title('Primes → zeros: periodic orbits of period m·log p rebuild the zero positions (dashed: true zeros)', fontsize=10)
pz = TR['primes_from_zeros']; x = np.array(pz['x'])
axs[1].plot(x, pz['S'], color=B, lw=1.1, label=f"Σ cos(γx) over {pz['nzeros']} zeros")
lab = {2: '2', 3: '3', 4: '4', 5: '5', 7: '7', 8: '8', 9: '9', 11: '11', 13: '13', 16: '16', 17: '17', 19: '19', 23: '23'}
for n, l in lab.items():
    axs[1].axvline(math.log(n), color=O, lw=0.8, ls='--'); axs[1].text(math.log(n), max(pz['S'])*0.9, l, ha='center', fontsize=8, color=O)
axs[1].set_xlabel('x (dashed: log of prime powers)'); axs[1].legend(fontsize=8.5, loc='lower right')
axs[1].set_title('Zeros → primes: dips (absorption, opposite sign to Gutzwiller) at log 2, log 3, log 4, …', fontsize=10)
plt.tight_layout(); plt.savefig('figures/fig6-trace-formula.png', dpi=150, bbox_inches='tight'); plt.close()

# ---------- fig7: orbit weights -> multiplicities
fig, ax = plt.subplots(figsize=(10, 4.8))
N = np.arange(1, 9)
for t, c, m in [('2', B, 'o'), ('3', '#08306b', 's'), ('0.5', R, '^'), ('-1', G_, 'D')]:
    v = np.array(TR['mult_t'][t]); vv = np.sign(v)*np.log10(1+np.abs(v))
    ax.plot(N, vv, m+'-', color=c, label=f'odds t = {t}')
ax.axhline(0, color='k', lw=0.6)
ax.set_xlabel('N: the zeros ρ/N of ζ(Ns)'); ax.set_ylabel('multiplicity M(t, N)  (signed log scale)')
ax.text(4.55, 0.12, 't = 2, 3: 2, 1, 2, 3, 6, 9, 18, 30 / 3, 3, 8, 18, 48 …\nWitt numbers: Lyndon words, free Lie algebra', fontsize=8.5, color=B)
ax.set_ylim(-1.05, 3.1); ax.text(3.0, -0.95, 't = 1/2: 0.5, −0.125, … not a spectrum\nt = −1: zeros of ζ(2s) enter as poles (Fermi)', fontsize=8.5, color=R)
ax.legend(fontsize=8.5, loc='upper left'); ax.set_title('Weighting repetition m of each prime orbit by t^m dresses the spectrum: whole multiplicities only at whole t', fontsize=10)
plt.savefig('figures/fig7-orbit-weights.png', dpi=150, bbox_inches='tight'); plt.close()

# ---------- fig8: solved shapes, numerically
fig, axs = plt.subplots(1, 3, figsize=(15, 5))
th_ = np.linspace(0, 2*np.pi, 400)
ax = axs[0]
for name, c, m in [('petersen', GR, 'o'), ('necklace10', R, 'x')]:
    r = np.array(SS['graphs'][name]['roots']); ax.plot(r[:, 0], r[:, 1], m, color=c, ms=6, label={'petersen': 'Petersen (Ramanujan)', 'necklace10': 'necklace of 10 beads (not)'}[name])
ax.plot(np.cos(th_)/math.sqrt(2), np.sin(th_)/math.sqrt(2), color='k', lw=0.8)
ax.set_aspect('equal'); ax.set_title('Ihara zeta of 3-regular graphs:\nroots of 1 − λu + 2u²; circle |u| = 1/√2', fontsize=9.5); ax.legend(fontsize=8, loc='upper center', bbox_to_anchor=(0.5, -0.08))
ax = axs[1]
for key, c, m in [('petersen ferro J=0.5', GR, 'o'), ('petersen antiferro J=-0.5', R, 'x')]:
    z = np.array(SS['lee_yang'][key]); w = z[:, 0]+1j*z[:, 1]; w = np.where(np.abs(w) > 1, 1/np.conj(w), w)
    ax.plot(w.real, w.imag, m, color=c, ms=6, label=key.replace('petersen ', 'Petersen Ising, ').replace('J=', 'J = '))
ax.plot(np.cos(th_), np.sin(th_), color='k', lw=0.8); ax.set_aspect('equal')
ax.set_title('Lee–Yang zeros in the fugacity z\n(points outside the circle reflected in)', fontsize=9.5); ax.legend(fontsize=8, loc='upper center', bbox_to_anchor=(0.5, -0.08))
ax = axs[2]
ap = np.array(SS['elliptic']['ap']); ax.plot(ap[:, 0], ap[:, 1]/(2*np.sqrt(ap[:, 0])), '.', color=GR, ms=4, label='y² = x³ + x + 1: a_p / 2√p')
ax.axhline(1, color='k', lw=0.8); ax.axhline(-1, color='k', lw=0.8); ax.axhline(0, color=G_, lw=0.4)
ax.set_ylim(-1.6, 1.6); ax.set_xlabel('p'); ax.set_title('Hasse: every local factor 1 − a_p u + p u²\nhas its roots on |u| = 1/√p (Δ_p ≤ 0)', fontsize=9.5)
ax.text(100, 1.25, 'ζ’s own factor 1 − u: d = 0, Δ = 1 > 0\nroot at |u| = 1: the local test is empty', fontsize=8.5, color=R)
ax.legend(fontsize=8, loc='lower right')
plt.tight_layout(); plt.savefig('figures/fig8-solved-shapes-numbers.png', dpi=150, bbox_inches='tight'); plt.close()
print('figures done')
