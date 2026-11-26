"""Pictures for the elliptic-object DRAFT. Colour key: blue = our reading, orange = already known,
grey = reference, red = where a candidate breaks."""
import json, math, numpy as np, matplotlib
matplotlib.use('Agg'); import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, Rectangle
B, O, G, R, GR = '#1f77b4', '#ff7f0e', '#888888', '#d62728', '#2ca02c'
plt.rcParams.update({'font.size': 10, 'axes.spines.top': False, 'axes.spines.right': False})
d = json.load(open('dial.json')); li = json.load(open('li.json'))
gam = np.array([r['gamma'] for r in d['C']])
EG = 0.5772156649015329
def asym(n): n = np.asarray(n, float); return n / 2 * (np.log(n) - 1 + EG - np.log(2 * np.pi))

# ---------------- fig0: the map
fig, ax = plt.subplots(figsize=(11.5, 5.6)); ax.axis('off'); ax.set_xlim(0, 11.5); ax.set_ylim(0, 5.7)
ax.add_patch(FancyBboxPatch((4.25, 2.35), 3.0, 1.0, boxstyle='round,pad=0.1', fc='#eef', ec=B))
ax.text(5.75, 2.85, 'An elliptic object for zeta:\nfactors of degree 2, all on the\nelliptic side Δ = k² − 4d ≤ 0', ha='center', va='center', fontsize=10)
items = [((0.15, 4.3), 'A. Pair zeta\'s prime factors\n(any pairing, incl. s ↔ 1 − s)', 'never elliptic: Δ = (a − b)² ≥ 0', R),
         ((0.15, 2.4), 'B. The prime at infinity\n(Gamma factor, trivial zeros)', 'hyperbolic side, k < −2', R),
         ((0.15, 0.5), 'L(E, s): every prime factor\nelliptic (Hasse), yet its own\nzeros are still open', 'local elliptic ≠ global RH', O),
         ((8.2, 4.3), 'C. Cayley dial: zeros paired\nby s ↔ 1 − s, read at z = 1 − 1/ρ', 'elliptic ⇔ RH, but built\nfrom the zeros (circular)', B),
         ((8.2, 2.4), 'Positivity reading:\nλ_n = Σ 2(1 − T_n(−k/2)) ≥ 0', '= Li\'s criterion (known)', O),
         ((8.2, 0.5), 'Numerical check on\nthe known zeros', 'passes, but blind past\nthe verified height', G)]
for (x, y), txt, verdict, col in items:
    ax.add_patch(FancyBboxPatch((x, y), 3.15, 1.2, boxstyle='round,pad=0.08', fc='white', ec=col, lw=1.8))
    ax.text(x + 1.575, y + 0.78, txt, ha='center', va='center', fontsize=8.6)
    ax.text(x + 1.575, y + 0.2, verdict, ha='center', va='center', fontsize=8.4, color=col, weight='bold')
    ax.annotate('', xy=(4.25 if x < 4 else 7.25, 2.85), xytext=(x + 3.15 if x < 4 else x, y + 0.6),
                arrowprops=dict(arrowstyle='->', color=col, lw=1.2))
ax.set_title('Three candidate elliptic objects for zeta, a positivity reading, and where each one ends', fontsize=11)
plt.savefig('figures/fig0-map.png', dpi=150, bbox_inches='tight'); plt.close()

# ---------------- fig1: pairing prime factors never goes elliptic
fig, axs = plt.subplots(1, 2, figsize=(13, 4.8))
ax = axs[0]
kk = np.linspace(-7, 6.5, 400); ax.fill_between(kk, kk ** 2 / 4, 12, color='#e8f0fb'); ax.plot(kk, kk ** 2 / 4, color='k', lw=1)
rng = np.random.default_rng(1); a = rng.uniform(-3, 3, 400); b = rng.uniform(-3, 3, 400)
ax.scatter(-(a + b), a * b, s=6, color=R, label='(1 − a·u)(1 − b·u), a, b real: Δ = (a − b)²')
ax.scatter([-(1 + p) / 1 for p in (2, 3, 5)], [p for p in (2, 3, 5)], s=55, color='k', marker='s', zorder=3,
           label="zeta's factor with its mirror: (1 − u)(1 − p·u)")
for p in (2, 3, 5): ax.annotate(f'p = {p}', (-(1 + p), p), textcoords='offset points', xytext=(6, -4), fontsize=8)
ax.text(0, 7.5, 'elliptic side\n(Δ < 0: both roots on one circle)', ha='center', color=B, fontsize=9)
ax.text(2.2, -6, 'hyperbolic side (Δ > 0)', color=R, fontsize=9)
ax.set_xlim(-7, 6.5); ax.set_ylim(-9, 12); ax.set_xlabel('k'); ax.set_ylabel('d')
ax.set_title('(a) The dial plane: any two real degree-1 factors\nmultiply to a point on or below the parabola', fontsize=10)
ax.legend(fontsize=7.5, loc='lower left')
ax = axs[1]
ax.axvspan(-2, 2, color='#e8f0fb'); ax.axvline(-2, color='k', lw=0.8); ax.axvline(2, color='k', lw=0.8)
A1 = d['A1']; A2 = d['A2']; A3 = d['A3']
ax.scatter([r['k_hat'] for r in A1], np.full(len(A1), 3), s=14, color=R)
ax.text(-9.8, 3.2, "zeta, prime paired with mirror:\nk̂ = −(√p + 1/√p) ≤ −2.12", fontsize=8, color=R)
ax.scatter([r['k_hat'] for r in A2], 2 + 0.06 * rng.standard_normal(len(A2)), s=14, color=O)
ax.text(-9.8, 2.25, 'complex characters (mod 5, mod 7) with their conjugates:\nelliptic, but d = 1, so the roots sit on Re s = 0', fontsize=8, color=O)
ax.scatter([r['k_hat'] for r in A3], 1 + 0.03 * rng.standard_normal(len(A3)), s=14, color=B)
ax.text(-9.8, 1.25, 'elliptic curve y² = x³ + x + 1: elliptic with d = p (Hasse),\nroots on Re s = 1/2 for each prime factor', fontsize=8, color=B)
ax.set_xlim(-10, 3); ax.set_ylim(0.4, 3.8); ax.set_yticks([]); ax.set_xlabel('k̂ = k/√d   (elliptic segment |k̂| ≤ 2 shaded)')
ax.set_title('(b) Every prime factor read on one ruler\n(zeta and characters: primes below 200; curve: primes below 3000)', fontsize=10)
plt.savefig('figures/fig1-pairing-never-elliptic.png', dpi=150, bbox_inches='tight'); plt.close()

# ---------------- fig2: the Cayley map
fig, axs = plt.subplots(1, 3, figsize=(16, 5.2), gridspec_kw={'width_ratios': [1, 1.1, 0.8]})
ax = axs[0]
ax.axvspan(0, 1, color='#f2f2f2'); ax.axvline(0.5, color=B, lw=1.5)
g40 = gam[:12]
ax.plot(np.full(len(g40), 0.5), g40, 'o', color=B, ms=5, label='first zeros ρ = 1/2 + iγ')
ax.plot(np.full(len(g40), 0.5), -g40, 'o', color=B, ms=5)
ax.plot([-2, -4, -6, -8], [0] * 4, 's', color=R, ms=6, label='trivial zeros −2, −4, … (Gamma factor)')
ax.plot([1], [0], 'x', color='k', ms=9, mew=2, label='pole at s = 1')
ax.plot([0.5 + 0.3, 0.5 - 0.3], [20, 20], 'D', color=R, mfc='none', ms=7, label='a toy off-line pair ρ, 1 − ρ̄')
ax.set_xlim(-9, 3); ax.set_ylim(-45, 45); ax.set_xlabel('Re s'); ax.set_ylabel('Im s')
ax.set_title('(a) s-plane: the strip and the critical line', fontsize=10); ax.legend(fontsize=7.5, loc='lower left')
ax = axs[1]; ax.set_aspect('equal')
th = np.linspace(0, 2 * np.pi, 400); ax.plot(np.cos(th), np.sin(th), color=B, lw=1.5)
z = 1 - 1 / (0.5 + 1j * g40); ax.plot(z.real, z.imag, 'o', color=B, ms=5); ax.plot(z.real, -z.imag, 'o', color=B, ms=5)
zt = np.array([1 + 1 / (2 * m) for m in range(1, 5)]); ax.plot(zt, 0 * zt, 's', color=R, ms=6)
ax.plot(1 / zt, 0 * zt, 's', color=R, mfc='none', ms=6)
for rho in (0.8 + 20j, 0.2 + 20j):
    zz = 1 - 1 / rho; ax.plot(zz.real, zz.imag, 'D', color=R, mfc='none', ms=7)
ax.plot([0], [0], 'x', color='k', ms=9, mew=2)
ax.text(-0.95, 1.12, 'critical line → unit circle |z| = 1\nmirror s ↔ 1 − s → z ↔ 1/z', fontsize=8.5, color=B)
ax.set_xlim(-1.2, 1.7); ax.set_ylim(-1.25, 1.35); ax.set_xlabel('Re z'); ax.set_ylabel('Im z')
ax.set_title('(b) z = 1 − 1/s: zeros crowd toward z = 1; filled squares\ntrivial zeros, open squares their mirrors', fontsize=10)
tz = np.linspace(-0.01, 0.09, 600); ins = axs[2]; ins.plot(np.cos(tz), np.sin(tz), color=B, lw=1, label='unit circle')
zz = 1 - 1 / (0.5 + 1j * gam[:60]); ins.plot(zz.real, zz.imag, '.', color=B, ms=4, label='first 60 zeros')
for i, rho in enumerate((0.8 + 20j, 0.2 + 20j)):
    w_ = 1 - 1 / rho; ins.plot(w_.real, w_.imag, 'D', color=R, mfc='none', ms=7, label='toy off-line pair (one inside, one outside)' if i == 0 else None)
ins.set_xlim(0.994, 1.0015); ins.set_ylim(-0.003, 0.075); ins.set_xlabel('Re z'); ins.legend(fontsize=7.5, loc='upper right')
ins.set_title('(c) zoom near z = 1', fontsize=10)
plt.savefig('figures/fig2-cayley-map.png', dpi=150, bbox_inches='tight'); plt.close()

# ---------------- fig3: the Cayley dial line, off-line lift, curvature
fig, axs = plt.subplots(1, 3, figsize=(15, 4.6))
ax = axs[0]
kp2 = np.array([r['k_plus_2'] for r in d['C']])
ax.loglog(gam, kp2, '.', color=B, ms=2, label='5,000 zeros: k + 2')
gg = np.logspace(0.5, 4, 100); ax.loglog(gg, 1 / (0.25 + gg ** 2), color='k', lw=0.8, ls='--', label='1/|ρ|² = 1/(1/4 + γ²)')
ax.set_xlabel('height γ'); ax.set_ylabel('distance from the edge k = −2')
ax.set_title('(a) On the line d = 1 every zero sits inside\nthe elliptic segment, 1/|ρ|² from its edge', fontsize=10); ax.legend(fontsize=8)
ax = axs[1]
ax.plot([-2, 2], [0, 0], color='#9cc0ea', lw=8, solid_capstyle='butt', zorder=0, label='elliptic segment |k| ≤ 2')
ax.annotate('real zeros just inside −2,\ntrivial zeros just outside', (-2, 0), xytext=(0.3, 0.0022), fontsize=7.5, arrowprops=dict(arrowstyle='->', color=G))
ax.plot([r['k_re'] for r in d['C'][:40]], [0] * 40, 'o', color=B, ms=4, label='real zeros (all inside, near −2)')
ax.plot([r['k'] for r in d['B'][:12]], [0] * 12, 's', color=R, ms=5, label='trivial zeros: k < −2 (hyperbolic)')
off = [r for r in d['C_off'] if 'gamma0' in r and abs(r['gamma0'] - 5) < 1e-9]
ax.plot([r['k_re'] for r in off], [r['k_im'] for r in off], 'D', color=R, mfc='none', ms=6, label='toy off-line zeros (γ = 5): k leaves the real line')
for r in off: ax.annotate(f"δ={r['delta']}", (r['k_re'], r['k_im']), textcoords='offset points', xytext=(4, 2), fontsize=6.5)
roff = [r for r in d['C_off'] if 'real_rho' in r]
ax.plot([r['k_re'] for r in roff], [0] * len(roff), 'v', color=R, ms=6, label='toy real zeros in (1/2, 1): k > 2')
ax.set_xlim(-2.6, 7); ax.set_xlabel('Re k'); ax.set_ylabel('Im k')
ax.set_title('(b) What an off-line zero does to its dial point', fontsize=10); ax.legend(fontsize=7, loc='upper right')
ax = axs[2]
ax.loglog(gam, [r['curv'] for r in d['C']], '.', color=B, ms=2, label='curvature d − k²/4 of each zero pair')
ax.loglog(gg, gg ** 2 / (gg ** 2 + 0.25) ** 2, color='k', lw=0.8, ls='--', label='γ²/(γ² + 1/4)²')
ax.set_xlabel('height γ'); ax.set_ylabel('curvature d − k²/4')
ax.set_title('(c) Every pair reads as an ellipsoid\n(positive curvature); off-line → complex', fontsize=10); ax.legend(fontsize=8)
plt.savefig('figures/fig3-cayley-dial.png', dpi=150, bbox_inches='tight'); plt.close()

# ---------------- fig4: Li = Chebyshev traces
fig, axs = plt.subplots(1, 3, figsize=(15, 4.6))
ax = axs[0]; x = np.linspace(-1.15, 1.15, 800)
for n, c in ((3, G), (6, O), (12, B)):
    ax.plot(x, np.polynomial.chebyshev.chebval(x, [0] * n + [1]), color=c, lw=1.2, label=f'T_{n}(x)')
ax.axvspan(-1, 1, color='#e8f0fb'); ax.set_ylim(-2.5, 4); ax.axhline(1, color='k', lw=0.6, ls=':')
ax.set_xlabel('x = −k/2'); ax.set_title('(a) Elliptic (|x| ≤ 1) ⇔ every T_n(x) ≤ 1;\noutside, T_n blows up', fontsize=10); ax.legend(fontsize=8)
ax = axs[1]; lam = np.array([float(v) for v in li['lambda']]); n = np.arange(1, len(lam) + 1)
ax.plot(n, lam, color=B, lw=1.5, label='λ_n from zeta itself (Cauchy integral, 140 digits)')
ax.plot(n, asym(n), color=O, lw=1, ls='--', label='(n/2)(log n − 1 + γ − log 2π)')
ax.set_xlabel('n'); ax.set_ylabel('λ_n'); ax.set_title('(b) λ_n = Σ over zero pairs of 2(1 − T_n(−k/2)):\npositive for every n ≤ 1000', fontsize=10); ax.legend(fontsize=8)
ax = axs[2]; ax.plot(n, lam - asym(n), color=B, lw=1)
ax.axhline(0, color='k', lw=0.6); ax.set_xlabel('n'); ax.set_ylabel('λ_n − trend')
ax.set_title('(c) What is left after the trend: small, bounded wiggle', fontsize=10)
plt.savefig('figures/fig4-li-chebyshev.png', dpi=150, bbox_inches='tight'); plt.close()

# ---------------- fig5: detection
fig, axs = plt.subplots(1, 2, figsize=(13.5, 4.8))
ax = axs[0]
base = np.array([[r['n'], r['zero_sum'] + r['tail']] for r in d['D']])
ax.plot(base[base[:, 0] <= 9000, 0], base[base[:, 0] <= 9000, 1], color=B, lw=2, label='λ_n, all zeros on the line')
cols = [R, '#9467bd', O, '#8c564b']
for t, c in zip(d['toy'], cols):
    cv = np.array(t['curve']); ax.plot(cv[:, 0], cv[:, 1], color=c, lw=0.7,
        label=f"double zero at γ = {t['gamma0']:.2f} split by δ = {t['delta']}: first λ_n < 0 at n = {t['first_negative_n']}")
    ax.plot(t['first_negative_n'], 0, 'o', color=c, ms=6, zorder=4)
ax.axhline(0, color='k', lw=0.6); ax.set_yscale('symlog', linthresh=100); ax.set_xlim(0, 9000)
ax.set_xlabel('n'); ax.set_ylabel('λ_n (symmetric log scale)')
ax.set_title('(a) Toy: one pair moved off the line stays hidden until\nits wobble outgrows the trend, then λ_n turns negative', fontsize=10); ax.legend(fontsize=7, loc='lower left')
ax = axs[1]
def nstar(g0, dl):
    eps = 0.5 * math.log1p(2 * dl / ((0.5 - dl) ** 2 + g0 ** 2))
    n = 10.0
    for _ in range(200):
        tr = n / 2 * (math.log(n) - 1 + EG - math.log(2 * math.pi)) + 4
        n = max(10.0, math.log(max(tr, 1) / 2) / eps)
    return n
G0 = np.logspace(0.7, 13, 60)
for dl, c in ((0.45, R), (0.1, O), (0.001, G)):
    ax.loglog(G0, [nstar(g, dl) for g in G0], color=c, label=f'off by δ = {dl}')
for t, c in zip(d['toy'], cols): ax.loglog(t['gamma0'], t['first_negative_n'], 'o', color=c, ms=6)
ax.axvline(3e12, color='k', ls=':', lw=1); ax.text(3e12, 1e3, ' zeros checked\n on the line\n up to here', fontsize=8)
ax.axhline(1e5, color='k', ls='--', lw=0.8); ax.text(1e2, 2.5e5, 'largest n computed so far (10^5, Johansson)', fontsize=8)
ax.set_xlabel('height γ of the off-line pair'); ax.set_ylabel('first n with λ_n < 0')
ax.set_title('(b) Needed n grows like γ²/δ: dots are the toy runs', fontsize=10); ax.legend(fontsize=8, loc='upper left')
plt.savefig('figures/fig5-detection-reach.png', dpi=150, bbox_inches='tight'); plt.close()

# ---------------- fig6: standard names
fig, ax = plt.subplots(figsize=(11.5, 5.2)); ax.axis('off'); ax.set_xlim(0, 11.5); ax.set_ylim(0, 6.2)
rows = [('Dial point 1 + k·u + d·u², elliptic side Δ ≤ 0', 'Frobenius pair; Hasse (1933), Weil (1948) for curves; Deligne (1974)'),
        ('k̂ = k/√d on the elliptic segment |k̂| ≤ 2', 'Trace of SL₂(ℝ): elliptic / parabolic / hyperbolic; Sato–Tate angle'),
        ('Pairing degree-1 factors gives Δ = (a − b)² ≥ 0', 'Zeta has no local H¹: its prime factors are the "poles" (H⁰, H²) part'),
        ('Cayley dial z = 1 − 1/ρ, mirror z ↔ 1/z', 'Li (1997); Keiper (1992); the Cayley transform of the critical line'),
        ('λ_n = Σ 2(1 − T_n(−k/2)) ≥ 0 for all n', "Li's criterion; Bombieri–Lagarias (1999); Voros"),
        ('Curvature d − k²/4 > 0 for each pair (ellipsoid)', 'Blaschke affine spheres (our affine-sphere/ folder)'),
        ('The missing elliptic object (a "Frobenius" for Spec ℤ)', 'Hilbert–Pólya; Weil positivity; Deninger; Connes; Manin, F₁')]
for i, (a, b) in enumerate(rows):
    y = 5.6 - i * 0.8
    ax.add_patch(FancyBboxPatch((0.2, y - 0.3), 4.9, 0.6, boxstyle='round,pad=0.05', fc='white', ec=B))
    ax.text(2.65, y, a, ha='center', va='center', fontsize=8.6)
    ax.add_patch(FancyBboxPatch((6.3, y - 0.3), 5.0, 0.6, boxstyle='round,pad=0.05', fc='white', ec=O))
    ax.text(8.8, y, b, ha='center', va='center', fontsize=8.3)
    ax.annotate('', xy=(6.3, y), xytext=(5.1, y), arrowprops=dict(arrowstyle='->', color=G))
ax.text(2.65, 6.15, 'Our reading', ha='center', color=B, weight='bold'); ax.text(8.8, 6.15, 'Standard name', ha='center', color=O, weight='bold')
plt.savefig('figures/fig6-standard-names.png', dpi=150, bbox_inches='tight'); plt.close()

# ---------------- fig7: probe spin (precession) picture
pr = json.load(open('probe.json')); t = np.array(pr['t'])
fig, axs = plt.subplots(1, 2, figsize=(13.5, 4.8))
ax = axs[0]
for beta, c, ls in (('2.0', G, '-'), ('1.5', O, '-'), ('1.1', R, '-'), ('0.75', '#9467bd', '--'), ('0.5', B, '-')):
    ax.plot(t, pr['abs_zeta'][beta], color=c, ls=ls, lw=1.3 if beta == '0.5' else 1,
            label=f'β = {beta}' + ('  (bath exists)' if float(beta) > 1 else '  (no bath: sum diverges)'))
ax.plot(gam[gam < 50], 0 * gam[gam < 50], 'o', color=B, ms=5)
ax.set_xlabel('probe time t'); ax.set_ylabel('|ζ(β + i t)|  (probe coherence × ζ(β))'); ax.set_ylim(0, 4)
ax.set_title('(a) A probe spin on the primon gas: its coherence hits zero\nonly at β = 1/2 (dots), past the edge β = 1 where the gas exists', fontsize=10)
ax.legend(fontsize=7.5, loc='upper right')
ax = axs[1]
bs = np.linspace(1.02, 3, 100)
import mpmath as _mp
lb = [float(_mp.zeta(2 * b) / _mp.zeta(b) ** 2) for b in bs]
ax.plot(bs, lb, color=O, label='lower bound ζ(2β)/ζ(β)² (from the Möbius series)')
for beta, v in pr['min_coherence'].items(): ax.plot(float(beta), v, 'o', color=B, ms=7)
ax.plot([], [], 'o', color=B, label='smallest coherence seen, t in (0, 50]')
ax.axvline(1, color='k', ls=':', lw=1); ax.text(1.02, 0.45, ' edge: the gas\n stops existing', fontsize=8)
ax.set_xlabel('β'); ax.set_ylabel('coherence |L_β(t)|'); ax.set_xlim(0.9, 3); ax.set_ylim(0, 0.8)
ax.set_title('(b) Wherever the bath is real, coherence never vanishes', fontsize=10); ax.legend(fontsize=8, loc='lower right')
plt.savefig('figures/fig7-probe-spin.png', dpi=150, bbox_inches='tight'); plt.close()

# ---------------- fig8: Li split, prime at infinity vs finite primes
sp_ = json.load(open('li_split.json')); A_ = np.array([float(v) for v in sp_['arch']]); P_ = np.array([float(v) for v in sp_['prime']])
nn = np.arange(1, len(A_) + 1)
fig, axs = plt.subplots(1, 2, figsize=(13.5, 4.8))
ax = axs[0]
ax.plot(nn[:30], A_[:30], 'o-', color=O, ms=3, label='prime at infinity (Gamma factor)')
ax.plot(nn[:30], P_[:30], 's-', color=R, ms=3, label='finite primes ((s − 1)ζ(s))')
ax.plot(nn[:30], (A_ + P_)[:30], 'D-', color=B, ms=3, label='λ_n (their sum)')
ax.axhline(0, color='k', lw=0.6); ax.axvline(10.5, color=G, ls=':'); ax.text(11, 8, 'from n = 11 the trend\nalone outweighs the primes', fontsize=8)
ax.set_xlabel('n'); ax.set_title('(a) First 30: the primes carry the positivity for n ≤ 10', fontsize=10); ax.legend(fontsize=8)
ax = axs[1]
ax.plot(nn, A_, color=O, lw=1.5, label='prime at infinity: grows like (n/2) log n')
ax.plot(nn, P_, color=R, lw=1, label='finite primes: a bounded wiggle')
ax.set_yscale('symlog', linthresh=10); ax.axhline(0, color='k', lw=0.6)
ax.set_xlabel('n'); ax.set_title(f'(b) n ≤ 1000: the prime part stays between {P_.min():.1f} and {P_.max():.1f}', fontsize=10); ax.legend(fontsize=8, loc='upper left')
plt.savefig('figures/fig8-li-split.png', dpi=150, bbox_inches='tight'); plt.close()
print('figures done')
