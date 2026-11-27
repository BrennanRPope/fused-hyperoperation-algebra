"""Pictures for the hat-and-primes write-up. Colour key everywhere:
green = built and measured here, blue = our earlier pieces, orange = known results, grey = reference."""
import json, numpy as np, matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.collections import PolyCollection
import hatgen as h

GREEN, BLUE, ORANGE, GREY, RED, GOLD = '#2e8b57', '#2b6cb0', '#dd6b20', '#888888', '#c53030', '#d69e2e'
KC = {'H': '#90cdf4', 'H1': '#1a365d', 'T': '#cbd5e0', 'P': '#fefcbf', 'F': '#c6f6d5'}
sub = json.load(open('substitution.json')); lat = json.load(open('lattice_primes.json'))
spec = json.load(open('spectra.json')); curves = np.load('spectra_curves.npz')
plt.rcParams.update({'font.size': 11})

def draw_hats(ax, hats, colors, lw=0.3):
    polys = [h.hat_polygon(T) for T, k in hats]
    ax.add_collection(PolyCollection(polys, facecolors=colors, edgecolors='k', linewidths=lw))
    allp = np.vstack(polys); ax.set_xlim(allp[:, 0].min(), allp[:, 0].max()); ax.set_ylim(allp[:, 1].min(), allp[:, 1].max())
    ax.set_aspect('equal'); ax.axis('off')

L = h.levels(3)
# ---------------- Figure 1: a real hat patch ----------------
fig, ax = plt.subplots(figsize=(10, 9))
hs = h.hats_of(L[3][0])
draw_hats(ax, hs, [KC[k] for T, k in hs])
from matplotlib.patches import Patch
ax.legend(handles=[Patch(fc=KC['H'], ec='k', label='hat in an H cluster'),
                   Patch(fc=KC['H1'], ec='k', label='reflected hat (one per H)'),
                   Patch(fc=KC['T'], ec='k', label='hat in a T'), Patch(fc=KC['P'], ec='k', label='hat in a P'),
                   Patch(fc=KC['F'], ec='k', label='hat in an F')], loc='lower left', fontsize=9, framealpha=0.95)
ax.set_title('Figure 1. 1,156 hats (= 34², a level-3 H supertile), generated and checked here:\n'
             'no overlaps, every corner on one triangular lattice. It never repeats (Smith–Myers–Kaplan–Goodman-Strauss 2023).',
             fontsize=11)
fig.savefig('figures/fig1-hat-patch.png', dpi=130, bbox_inches='tight'); plt.close(fig)

# ---------------- Figure 2: how it avoids repeating ----------------
fig = plt.figure(figsize=(15, 5.2))
counts = [len(h.hats_of(L[n][0])) for n in range(4)]
fibs = [2, 5, 13, 34]
for n in range(4):
    ax = fig.add_axes([0.01 + n*0.18, 0.22, 0.17, 0.66])
    hs = h.hats_of(L[n][0])
    draw_hats(ax, hs, [KC[k] for T, k in hs], lw=0.5 if n < 3 else 0.2)
    ax.set_title(f'level {n} H supertile\n{counts[n]} hats = {fibs[n]}²', fontsize=11)
ax = fig.add_axes([0.74, 0.22, 0.25, 0.66]); ax.axis('off')
A = sub['A']
txt = ('Every supertile is built from the\nlevel below by the same rule:\n\n'
       '        old H  T  P  F\n'
       f'new H  [ {A[0][0]}  {A[0][1]}  {A[0][2]}  {A[0][3]} ]\n'
       f'new T  [ {A[1][0]}  {A[1][1]}  {A[1][2]}  {A[1][3]} ]\n'
       f'new P  [ {A[2][0]}  {A[2][1]}  {A[2][2]}  {A[2][3]} ]\n'
       f'new F  [ {A[3][0]}  {A[3][1]}  {A[3][2]}  {A[3][3]} ]\n\n'
       'Growth per level: φ⁴ = (7+3√5)/2 ≈ 6.854\n(area; lengths grow by φ²).\n'
       'φ⁴ is irrational, so no finite block\ncan ever repeat: that is the\nwhole reason it never repeats.')
ax.text(0, 1, txt, va='top', family='monospace', fontsize=10.5,
        bbox=dict(fc='#fffaf0', ec=ORANGE))
ax2 = fig.add_axes([0.02, 0.02, 0.7, 0.14]); ax2.axis('off')
ax2.text(0, 0.9, 'Hats per H supertile: 4, 25, 169, 1156, 7921, 54289 …  = (2, 5, 13, 34, 89, 233 …)²  — every other Fibonacci number, squared '
         '(green: checked here; the matrix is known).\n'
         'Unreflected : reflected hats → φ⁴ ≈ 6.854 : 1.   Metatile shares: H 1/3, T 7/6−√5/2, P √5−2, F 3/2−√5/2 — all in the golden field ℚ(√5).',
         va='top', fontsize=10.5, color=GREEN)
fig.suptitle('Figure 2. How the hat avoids repeating: one growth rule, golden-ratio scaling', fontsize=13, y=0.99)
fig.savefig('figures/fig2-how-it-never-repeats.png', dpi=120, bbox_inches='tight'); plt.close(fig)

# ---------------- Figure 3: counts modulo each prime ----------------
per = sub['periods']
ps = np.array([int(p) for p in per]); T_ = np.array([per[str(p)]['period'] for p in ps])
cls = [per[str(p)]['cls'] for p in ps]
fig, ax = plt.subplots(figsize=(11, 6))
xx = np.linspace(0, ps.max(), 10)
ax.plot(xx, xx - 1, color=GREY, lw=1, ls='--'); ax.plot(xx, xx + 1, color=GREY, lw=1, ls=':')
for c, col, lab in (('split', BLUE, 'p ends in 1 or 9 (p = ±1 mod 5): period divides p − 1'),
                    ('inert', ORANGE, 'p = ±2 mod 5: period divides p + 1'),
                    ('ramified', RED, 'p = 5: period 10')):
    m = np.array([q == c for q in cls])
    ax.scatter(ps[m], T_[m], color=col, s=28, label=lab, zorder=3)
ax.annotate('p = 2 and p = 3: period 6\n(an eigenvalue repeats mod 2 and mod 3,\nso these two sit off both lines)', (3, 6), xytext=(25, 150),
            textcoords='data', fontsize=9, arrowprops=dict(arrowstyle='->', color=GREY))
ax.text(250, 60, 'Every prime lands on the p − 1 line\nor the p + 1 line (or a divisor of it),\n'
        'and WHICH line is decided only by\np mod 5. The hat feels a prime only\nthrough its golden ratio √5.',
        fontsize=11, bbox=dict(fc='white', ec=GREEN))
ax.set_xlabel('prime p'); ax.set_ylabel('period of the hat counts mod p')
ax.set_title('Figure 3. Count the hats in each supertile and reduce mod p: the counts repeat, with this period\n'
             '(green: computed here for every prime below 400; the rule behind it is Frobenius in ℚ(√5), known)', fontsize=11)
ax.legend(loc='upper left', fontsize=9.5)
fig.savefig('figures/fig3-counts-mod-p.png', dpi=130, bbox_inches='tight'); plt.close(fig)

# ---------------- Figure 4: the hats on the Eisenstein lattice, sieve colours ----------------
d = np.load('patch_level4.npz'); Ts = d['T']; Aa = d['anchors']; kinds = d['kinds']
# window: a disc of radius 22 hexagons placed where the patch fully covers it
def nrm(P): return P[:, 0]**2 + P[:, 0]*P[:, 1] + P[:, 1]**2
cands = Aa[::97]
cnt = [int((nrm(Aa - c) <= 22**2).sum()) for c in cands]
C = cands[int(np.argmax(cnt))]
sel = nrm(Aa - C) <= 22**2
import lattice_primes_helpers as lph
sqf = lph.squarefree_mask(Aa[sel]); prm = lph.eis_prime_mask(Aa[sel])
fig, axs = plt.subplots(1, 2, figsize=(16, 8))
for ax, mode in zip(axs, ('prime', 'sieve')):
    hats = [(T, k) for T, k in zip(Ts[sel], kinds[sel])]
    if mode == 'prime':
        cols = [GOLD if p else '#edf2f7' for p in prm]
        title = 'Hats whose anchor is an Eisenstein prime (gold)'
    else:
        cols = ['#c6f6d5' if s else RED for s in sqf]
        title = 'Square-free sieve: anchor divisible by a prime squared (red)'
    draw_hats(ax, hats, cols, lw=0.3)
    # every hexagon centre in view, with an anchor or not
    pts = np.mgrid[-24:25, -24:25].reshape(2, -1).T; nn = pts[:, 0]**2 + pts[:, 0]*pts[:, 1] + pts[:, 1]**2
    pts = pts[nn <= 22**2] + C
    o = d['origin']
    # hex-centre coordinates (u,v) -> plane: z = 2(1+w6)(u + v w6) at half scale, plus origin
    def plane(P):
        u, v = P[:, 0], P[:, 1]
        # (1+w6)(u + v w6) = u + (u + v) w6 + v w6^2 = (u - v) + (u + 2v) w6   (w6^2 = w6 - 1)
        x, y = 2*(u - v), 2*(u + 2*v)
        return np.stack([0.5*(x + 0.5*y), 0.5*np.sqrt(3)/2*y], 1) + o
    Pp = plane(pts); occupied = set(map(tuple, Aa[sel].tolist()))
    emp = np.array([tuple(p) not in occupied for p in pts.tolist()])
    ax.scatter(Pp[emp, 0], Pp[emp, 1], s=9, facecolors='none', edgecolors=GREY, lw=0.6, zorder=4)
    ax.set_title(title + '\nhollow dots: hexagon centres with no hat anchor (1 in 4)', fontsize=11)
fig.suptitle('Figure 4. A surface that never repeats AND maps to primes: the hat lattice coloured by Eisenstein arithmetic.\n'
             'The primes come from the lattice colouring, not from the hat: anchors and empty centres get primes at the same rate (Figure 5).',
             fontsize=12)
fig.savefig('figures/fig4-hats-on-eisenstein-primes.png', dpi=120, bbox_inches='tight'); plt.close(fig)

# ---------------- Figure 5: the hat is blind to every prime ----------------
fig, axs = plt.subplots(1, 2, figsize=(15, 5.8), gridspec_kw={'width_ratios': [1.6, 1]})
ax = axs[0]
rows = lat['anchors']
labels = [f"N={r['norm']}" + ('' if r['kind'] == 'split' else f"\n({r['kind']})") for r in rows]
x = np.arange(len(rows))
lo = [r['min_share'] for r in rows]; hi = [r['max_share'] for r in rows]
ax.vlines(x, lo, hi, color=GREEN, lw=6)
ax.axhline(1, color=GREY, lw=1)
ax.set_xticks(x); ax.set_xticklabels(labels, fontsize=8.5)
ax.set_ylim(0.97, 1.03); ax.set_ylabel('share of each residue class\n(1 = perfectly even)')
ax.set_title('Hat anchors over the residue classes mod each Eisenstein prime\n(green bar = lowest to highest class; 129,314 hats)', fontsize=11)
ax = axs[1]
rc = lat['region_control']; sv = lat['sieve']
vals = [rc['prime_share_anchors'], rc['prime_share_empty_centres'], rc['squarefree_share_anchors'], rc['squarefree_share_empty'], sv['inv_zetaK2']]
names = ['prime\n(hats)', 'prime\n(empty)', 'sq-free\n(hats)', 'sq-free\n(empty)', 'sq-free\n1/ζ_K(2)']
ax.bar(range(5), vals, color=[GREEN, GREY, GREEN, GREY, ORANGE])
for i, v in enumerate(vals): ax.text(i, v + 0.01, f'{v:.4f}', ha='center', fontsize=9.5)
ax.set_xticks(range(5)); ax.set_xticklabels(names, fontsize=8.5); ax.set_ylim(0, 0.9)
ax.set_title('Where hats sit vs where they don\'t:\nsame prime rate, same square-free rate', fontsize=11)
fig.suptitle('Figure 5. The hat is blind to every prime of its lattice: no residue class is favoured, no prime is preferred', fontsize=12.5, y=1.04)
fig.savefig('figures/fig5-hat-is-blind-to-primes.png', dpi=130, bbox_inches='tight'); plt.close(fig)

# ---------------- Figure 6: three spectra ----------------
fig, axs = plt.subplots(1, 3, figsize=(19, 6))
ax = axs[0]
I4 = np.load('hat_diffraction_L4.npy'); M = I4.shape[0]
Ish = np.fft.fftshift(I4)
# reciprocal basis of the hex-centre lattice (a1=(1,0), a2=(1/2, sqrt3/2)): b1=(1,-1/sqrt3), b2=(0, 2/sqrt3)
hk = (np.arange(M) - M//2) / M
Hh, Kk = np.meshgrid(hk, hk, indexing='ij')
qx = Hh; qy = -Hh/np.sqrt(3) + 2*Kk/np.sqrt(3)
cvals = np.log10(Ish.ravel() + 1e-7)
for sh in (-1, 0, 1):
    for sk in (-1, 0, 1):
        H2, K2 = Hh + sh, Kk + sk
        X2_, Y2_ = H2, -H2/np.sqrt(3) + 2*K2/np.sqrt(3)
        m = (np.abs(X2_) <= 0.52) & (np.abs(Y2_) <= 0.52)
        ax.scatter(X2_[m], Y2_[m], c=cvals.reshape(M, M)[m], s=0.05, cmap='magma', vmin=-5, vmax=-1.5, rasterized=True)
for row in spec['pure_point_test'][:4]:
    hq, kq = row['q']; hq = hq - round(hq); kq = kq - round(kq)
    X_, Y_ = hq, -hq/np.sqrt(3) + 2*kq/np.sqrt(3)
    if abs(X_) < 0.5 and abs(Y_) < 0.5:
        ax.annotate('Bragg spot,\nsame height\nat 18.9k and\n129k hats', (X_, Y_), xytext=(X_*0.5, Y_*0.5 + 0.08), color='w', fontsize=7.5,
                    arrowprops=dict(arrowstyle='->', color='w', lw=0.6))
        break
ax.set_aspect('equal'); ax.set_xlim(-0.5, 0.5); ax.set_ylim(-0.5, 0.5)
ax.set_title('(a) the hat: sharp spots (pure point, Baake–Gähler–Sadun)\nspot positions = (a + b√5)/8: built from √5 only', fontsize=10.5)
ax.set_xlabel('frequency'); ax.set_ylabel('frequency')
ax = axs[1]
ax.semilogy(curves['ks'], curves['sqf_I'] + 1e-9, color=GREEN, lw=0.7)
for (a, q) in [(1, 4), (1, 2), (1, 9), (2, 9), (1, 3), (1, 12), (5, 12), (1, 25), (2, 25), (1, 5), (1, 7)]:
    if a/q > 0.5: continue
    ax.text(a/q, 0.12, f'{a}/{q}', fontsize=8.5, ha='center', color=ORANGE, rotation=90)
ax.set_xlim(0, 0.5); ax.set_ylim(1e-6, 0.4)
ax.set_ylabel('intensity (log scale)')
ax.set_title('(b) the square-free integers (sieved by every prime²):\nspots at fractions a/q, height ∏ over p | q of 1/(p²−1)²', fontsize=10.5)
ax.set_xlabel('frequency')
ax = axs[2]
ax.plot(curves['xs'], curves['S'], color=GREEN, lw=0.8)
for m in (2, 3, 4, 5, 7, 8, 9, 11, 13, 16, 17, 19, 23, 25, 27, 29, 31):
    ax.axvline(np.log(m), color=ORANGE, lw=0.6, alpha=0.5)
for m in (2, 3, 4, 5, 7, 9, 11, 13, 17, 19, 23, 29):
    ax.text(np.log(m), 70, str(m), ha='center', fontsize=8, color=ORANGE)
ax.set_title('(c) the first 1,000 zeta zeros: Σ cos(γ x) dips exactly at x = log(p^m)\n(Guinand–Weil explicit formula, known)', fontsize=10.5)
ax.set_xlabel('x'); ax.set_ylim(-180, 85)
fig.suptitle('Figure 6. Three ordered-but-never-repeating things, and what their spectra are built from: '
             '√5 (hat) · every prime at once (sieve) · log of every prime (zeros)', fontsize=12.5)
fig.savefig('figures/fig6-three-spectra.png', dpi=120, bbox_inches='tight'); plt.close(fig)

# ---------------- Figure 7: verdict map ----------------
fig, ax = plt.subplots(figsize=(14, 7.5)); ax.axis('off'); ax.set_xlim(0, 14); ax.set_ylim(0, 7.5)
def box(x, y, w, hh, text, ec, fc='white', fs=10.5):
    ax.add_patch(plt.Rectangle((x, y), w, hh, ec=ec, fc=fc, lw=2))
    ax.text(x + w/2, y + hh/2, text, ha='center', va='center', fontsize=fs, wrap=True)
box(0.3, 4.6, 4.2, 2.4, 'HAT / SPECTRE\nhidden space: a flat plane\n(cut and project from 4D)\nnumbers inside: √5 (hat), √15 (Spectre)\nsees a prime only as p mod 5', ORANGE, '#fffaf0')
box(4.9, 4.6, 4.2, 2.4, 'PRIME SIEVES\n(visible points, square-free sets)\nhidden space: every prime at once\n(Baake–Moody–Pleasants 2000)\nsees every prime individually', ORANGE, '#fffaf0')
box(9.5, 4.6, 4.2, 2.4, 'ZETA ZEROS\nspectrum at log p, log p², …\ninfinitely many independent\nfrequencies (log p are\nindependent over ℚ)', ORANGE, '#fffaf0')
box(0.3, 1.0, 4.2, 2.8, 'CLOSED as a prime map\nTested on 129,314 hats:\nflat over every Eisenstein prime\n(norm 3 to 43), same prime rate\non and off the hats. Its spectrum\nhas finite rank, so it can never\nbe "the zeros\' quasicrystal".', RED)
box(4.9, 1.0, 4.2, 2.8, 'OPEN, modest\nA sieve on the hat lattice is a\nnon-repeating surface that maps to\nprimes (Fig. 4). Its density is a\nk = −1 (union) composite over\nprimes: 1 − (∘₋₁ of 1/N(π)²) = 1/ζ_K(2)\nchecked 0.7780 vs 0.7781.', GREEN)
box(9.5, 1.0, 4.2, 2.8, 'WHERE THE REAL QUESTION IS\nDyson: is there a 1D quasicrystal\nwhose spectrum is the primes?\nIt would need a hidden space built\nfrom all primes (like the sieves),\nnot a flat one (like the hat).', BLUE)
for x in (2.4, 7.0, 11.6): ax.annotate('', (x, 3.85), (x, 4.55), arrowprops=dict(arrowstyle='->', lw=1.5))
ax.set_title('Figure 7. The verdict in one picture. Orange = known, green = built and checked here, red = closed, blue = the open question', fontsize=12)
fig.savefig('figures/fig7-verdict-map.png', dpi=120, bbox_inches='tight'); plt.close(fig)
print('figures written')
