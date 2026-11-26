"""Pictures for Riemann round 5. Colour key: blue = our reading, orange = already known, grey = reference,
red = where a road stops."""
import json, math, numpy as np, matplotlib
matplotlib.use('Agg'); import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, Rectangle
B, O, G_, R = '#1f77b4', '#ff7f0e', '#888888', '#d62728'
plt.rcParams.update({'font.size': 10, 'axes.spines.top': False, 'axes.spines.right': False})
p1 = json.load(open('pair1.json')); p2 = json.load(open('pair2.json')); lf = json.load(open('pass_light_flux.json'))
p3 = json.load(open('pair3.json'))
Z = np.array([float(x) for x in json.load(open('../elliptic-object/zeros.json'))['gamma']])

# ---------------- fig0: map
fig, ax = plt.subplots(figsize=(12, 6.4)); ax.axis('off'); ax.set_xlim(0, 12); ax.set_ylim(0, 6.4)
ax.add_patch(FancyBboxPatch((4.45, 2.7), 3.1, 1.0, boxstyle='round,pad=0.1', fc='#eef', ec=B))
ax.text(6.0, 3.2, 'The untried dial pairs on Riemann,\nmost fundamental first\n(C dropped: derived from k and d)', ha='center', va='center', fontsize=10)
items = [((0.15, 4.8), '1. Height s with homothety k\ntwo ways to bend the log', 'only the homothety bend keeps\nthe zeros; height moves the\nsingular value off 0', R),
         ((0.15, 2.65), '2. Ratio r with order q\nharmonic string ↔ prime fractal', 'order dial = a sum of homothety\ncubes; Moran dimension 1.3994', B),
         ((0.15, 0.5), '3. Light and flux tubes\nWick rotation, Berry–Keating', 'light turns each zero’s speed\n90°: off the line at once', R),
         ((8.65, 4.8), '4. Height s across roundness p\nLaguerre as a detector', 'fires exactly where kernels\nhave non-real zeros (p = 2.5, 3)', B),
         ((8.65, 2.65), 'RH equivalents found', 'generalised Mertens (k = 1, 2, 3);\norder-dial zeta analytic in\nRe s > 1/2: both known kinds', O),
         ((8.65, 0.5), 'Is any of it a lead?', 'No: each pair either restates\na known equivalent or loses\nthe zeros', G_)]
for (x, y), txt, verdict, col in items:
    ax.add_patch(FancyBboxPatch((x, y), 3.2, 1.4, boxstyle='round,pad=0.08', fc='white', ec=col, lw=1.8))
    ax.text(x + 1.6, y + 1.0, txt, ha='center', va='center', fontsize=8.4)
    ax.text(x + 1.6, y + 0.33, verdict, ha='center', va='center', fontsize=8.0, color=col, weight='bold')
    ax.annotate('', xy=(4.45 if x < 4 else 7.55, 3.2), xytext=(x + 3.2 if x < 4 else x, y + 0.7),
                arrowprops=dict(arrowstyle='->', color=col, lw=1.2))
ax.set_title('Round 5: the dial pairs we had not tried, plus light and flux tubes', fontsize=11)
plt.savefig('figures/fig0-map.png', dpi=150, bbox_inches='tight'); plt.close()

# ---------------- fig1: pair 1, the two bends
fig, axs = plt.subplots(1, 2, figsize=(13, 4.6))
ax = axs[0]
k = np.linspace(-1.6, 0.6, 400)
ax.plot(k, k*(k + 1)/12, color=B, lw=2)
ax.axhline(0, color=G_, lw=0.8)
for kk, lab in [(0, 'k = 0, height 1:\nboth are log'), (-1, 'k = −1, height 0:\nboth are X − 1')]:
    ax.plot([kk], [0], 'o', color=O, ms=8); ax.annotate(lab, (kk, 0), (kk - 0.05, 0.035), fontsize=8.5, ha='center')
ax.set_xlabel('homothety k   (height s = 1 + k gives the same bend)')
ax.set_ylabel('third-order gap between the two logs')
ax.set_title('The height bend and the homothety bend agree to second order;\nat third order they meet only at two settings (Lean bends_meet_twice)', fontsize=9.5)
ax = axs[1]
sv = p1['singular_value_X_s']; sa = p1['singular_value_above_1']
ss = [float(a) for a in sv] + [1.0] + [float(a) for a in sa]
xv = [sv[a] for a in sv] + [0.0] + [sa[a] for a in sa]
ax.plot(ss, xv, 'o-', color=B, lw=2, label='height dial: where ln_s(ζ) breaks')
ax.axhline(0, color=O, lw=2, ls='--', label='homothety dial: breaks at ζ = 0 for every k > 0')
ax.axhline(p1['zeta_half'], color=G_, lw=0.8, ls=':')
ax.text(0.45, p1['zeta_half'] - 0.18, 'ζ(1/2) = −1.4604: below this the break is a real point in (1/2, 1)', fontsize=7.5, color=G_)
ax.plot([1], [0], 'o', color=R, ms=9); ax.annotate('s = 1: the zeros of ζ\n(the only RH height)', (1, 0), (1.25, -0.9), fontsize=8.5,
                                                   arrowprops=dict(arrowstyle='->', color=R))
ax.set_xlabel('height s'); ax.set_ylabel('value of ζ where the log breaks'); ax.set_ylim(-2.8, 0.8)
ax.legend(fontsize=8, loc='lower right')
ax.set_title('Lowering the height slides the break from ζ = 0 to negative values;\nraising it slides it into (0, 1)', fontsize=9.5)
plt.tight_layout(); plt.savefig('figures/fig1-two-bends.png', dpi=150, bbox_inches='tight'); plt.close()

# ---------------- fig2: pair 1, a-points and generalised Mertens
fig, axs = plt.subplots(1, 3, figsize=(17, 4.8), gridspec_kw=dict(width_ratios=[1, 1.2, 0.8]))
ax = axs[0]
ax.axvline(0.5, color=G_, lw=1)
zz = Z[Z < 60]; ax.plot(0.5*np.ones_like(zz), zz, 'x', color='k', ms=6, label='zeros of ζ (s = 1)')
cols = {'0.5': '#9ecae1', '0.8': '#4292c6', '0.95': '#08519c'}
for s_, pts in p1['complex_a_points'].items():
    a = np.array(pts)
    if len(a): ax.plot(a[:, 0], a[:, 1], 'o', color=cols[s_], ms=5, label=f'height {s_}: ζ = {p1["singular_value_X_s"][s_]:.3f}')
for s_, r in p1['real_a_points_off_line'].items():
    ax.plot([r], [0.6], 's', color=R, ms=7)
ax.text(0.6, 2.0, 'heights 0.05, 0.1:\nreal breaks in (1/2, 1)', fontsize=8, color=R)
ax.set_xlim(0, 1.2); ax.set_ylim(0, 60); ax.set_xlabel('Re s'); ax.set_ylabel('Im s')
ax.legend(fontsize=7.5, loc='upper right')
ax.set_title('Where ln_s(ζ) breaks below height 1: a-points of ζ, off the line\n(left of it here; Bohr: infinitely many in Re s > 1 too)', fontsize=9.5)
ax = axs[1]
mk = p1['mertens_k']
for kk, c in [('1.0', 'k'), ('2.0', B), ('3.0', '#6baed6')]:
    x = np.array(mk[kk]['x']); ax.plot(x, mk[kk]['ratio'], color=c, lw=1.2,
                                        label=f'k = {kk[0]}: slope of running max {mk[kk]["slope_runmax"]:.2f}')
ax.set_xscale('log'); ax.axhline(0, color=G_, lw=0.6)
ax.set_xlabel('x'); ax.set_ylabel('M_k(x) / √x,   M_k(x) = Σ μ(n) k^ω(n)')
ax.legend(fontsize=8, loc='lower left')
ax.set_title('The homothety bend’s corner sums: RH ⇔ M_k(x) = O(x^(1/2+ε))\nfor k = 1, 2, 3 (k = 1 is Mertens/Littlewood)', fontsize=9.5)
ax = axs[2]
x = np.array(mk['0.5']['x']); ax.plot(x, mk['0.5']['main_ratio'], color=O, lw=1); ax.set_xscale('log'); ax.set_ylim(0.8, 2.2)
ax.axhline(1, color=G_, lw=0.8); ax.set_xlabel('x'); ax.set_ylabel('M_(1/2)(x) / main term')
ax.set_title('k = 1/2: no longer a test; it follows its\nSelberg–Delange main term (ratio 1.19 at 10⁷, still falling)', fontsize=9.5)
plt.tight_layout(); plt.savefig('figures/fig2-apoints-mertens.png', dpi=150, bbox_inches='tight'); plt.close()

# ---------------- fig3: pair 2, the prime fractal and its dimension
fig, axs = plt.subplots(1, 2, figsize=(13, 5.2))
ax = axs[0]; ax.set_aspect('equal'); ax.axis('off')
def primes(n):
    return [p for p in range(2, n) if all(p % d for d in range(2, int(p**0.5) + 1))]
PR = primes(140)
def pack(x0, y0, size, depth):
    # shelf packing of squares of side size/p inside a square of side size
    shelves = []; y = y0; x = x0; h = 0
    placed = []
    for p in PR:
        s = size/p
        if s*size < 1e-4 and depth > 0: break
        if x + s > x0 + size + 1e-12:
            y += h; x = x0; h = 0
        if y + s > y0 + size + 1e-12: break
        placed.append((x, y, s)); x += s; h = max(h, s)
    return placed
ax.add_patch(Rectangle((0, 0), 1, 1, fc='none', ec='k', lw=1.2))
lvl1 = pack(0, 0, 1.0, 0)
for (x, y, s) in lvl1:
    ax.add_patch(Rectangle((x, y), s, s, fc='#c6dbef', ec=B, lw=0.6))
    for (x2, y2, s2) in pack(x, y, s, 1)[:40]:
        ax.add_patch(Rectangle((x2, y2), s2, s2, fc='#2171b5', ec='white', lw=0.3))
ax.set_xlim(-0.02, 1.02); ax.set_ylim(-0.02, 1.02)
ax.set_title('The prime fractal: one copy shrunk by 1/p for every prime p\n(two levels; it fits in the plane because Σ 1/p² = 0.45 < 1)', fontsize=9.5)
ax = axs[1]
import mpmath as mp
Ds = np.linspace(1.1, 2.2, 200)
Pv = [float(mp.primezeta(d)) for d in Ds]
ax.plot(Ds, Pv, color=B, lw=2); ax.axhline(1, color=G_, lw=0.8)
ax.axvline(p2['moran_dimension'], color=R, lw=1.2, ls='--')
ax.text(p2['moran_dimension'] + 0.03, 1.6, f'D = {p2["moran_dimension"]:.5f}\n(round 4’s wall at q = 1)', fontsize=9, color=R)
ax.set_ylim(0, 3); ax.set_xlabel('D'); ax.set_ylabel('Σ_p p^(−D)   (prime zeta)')
ax.set_title('Moran’s equation Σ r_p^D = 1 with r_p = 1/p gives the dimension', fontsize=9.5)
plt.tight_layout(); plt.savefig('figures/fig3-prime-fractal.png', dpi=150, bbox_inches='tight'); plt.close()

# ---------------- fig4: pair 2, complex dimensions, the order dial, Z_q on the line
fig, axs = plt.subplots(1, 3, figsize=(15, 4.8))
ax = axs[0]
vs = json.load(open('pair2_valueset.json')); cdj = json.load(open('pair2_complexdims.json'))
th = np.linspace(0, 2*np.pi, 300)
ax.fill(vs['P_sigma']*np.cos(th), vs['P_sigma']*np.sin(th), color='#deebf7', label=f"every value P can take at σ = 1.2 (radius {vs['P_sigma']:.2f})")
ax.plot(vs['re'], vs['im'], color=B, lw=0.6, label='P(1.2 + it), 0 < t < 200')
ax.plot([1], [0], '*', color=R, ms=13, label=f"the value 1 (closest pass {vs['min_dist_t_gt_1']:.2f})")
ax.set_aspect('equal'); ax.set_xlim(-1.7, 1.7); ax.set_ylim(-1.7, 1.7); ax.set_xlabel('Re P'); ax.set_ylabel('Im P')
ax.legend(fontsize=7, loc='lower left')
ax.set_title('q = 1: complex dimensions are where P(s) = 1. They exist in\n1 < Re s < 1.3994 (1 is inside the disc) but none below t = 200', fontsize=9)
ax = axs[1]
qq = np.linspace(0, 0.8, 200)
ax.plot(qq, 1/(1 - qq), color=B, lw=2, label='our identity: leading power of ζ is 1/(1 − q)')
br = json.load(open('../round4/barrier.json'))['growth']
qs = [float(a) for a in br]; bf = [br[a]['beta_fit'] + 1 for a in br]
ax.plot(qs, bf, 'o', color=O, ms=8, label='round 4: fitted growth power + 1')
ax.set_xlabel('order dial q'); ax.set_ylabel('power of ζ at s = 1')
ax.legend(fontsize=8); ax.set_title('0 < q < 1: Z_q = (q;q)∞ Σ_m q^m/(q;q)_m · Π_j Z_(q^(m+j))^(q^(m+j)),\na sum of homothety cubes; the wall stays at 1 until q = 1', fontsize=9)
ax = axs[2]
for qk, c in [('0.5', B), ('0.3', '#6baed6')]:
    L = p2['critical_line'][qk]; ax.plot(L['t'], L['absZq'], color=c, lw=1.4, label=f'|Z_q(1/2 + it)|, q = {qk}')
ax.plot(L['t'], L['absZeta'], color=G_, lw=1, label='|ζ(1/2 + it)|')
for g in Z[(Z > 10) & (Z < 35)]: ax.axvline(g, color=R, lw=0.5, ls=':')
ax.set_xlabel('t'); ax.legend(fontsize=7.5); ax.set_yscale('log')
ax.set_title('Continued through the identity, Z_q vanishes at every zero of ζ\n(as a branch point: powers ζ^(q^m/(1−q)))', fontsize=9.5)
plt.tight_layout(); plt.savefig('figures/fig4-order-dial.png', dpi=150, bbox_inches='tight'); plt.close()

# ---------------- fig5: light, Wick rotation
fig, axs = plt.subplots(1, 2, figsize=(15, 4.8), gridspec_kw=dict(width_ratios=[1.5, 1]))
ax = axs[0]
v = np.array(lf['v_gamma_first40'])[:10]; g = Z[:10]
for n, (gi, vi) in enumerate(zip(g, v)):
    x0 = 3.0*n
    ax.annotate('', xy=(x0, vi), xytext=(x0, 0), arrowprops=dict(arrowstyle='-|>', color=O, lw=2))
    ax.annotate('', xy=(x0 + vi, 0), xytext=(x0, 0), arrowprops=dict(arrowstyle='-|>', color=B, lw=2))
    ax.plot([x0], [0], 'o', color='k', ms=4)
    ax.text(x0, 0.35, f'{gi:.1f}', ha='center', fontsize=7.5)
ax.set_aspect('equal'); ax.set_xlim(-2.6, 28.5); ax.set_ylim(-2.7, 3.4)
ax.set_xticks([]); ax.set_ylabel('speed (γ units)'); ax.set_xlabel('the first ten zeros (labelled by γ)')
ax.plot([], [], color=O, lw=3, label='heat (forward or backward heat flow): along the line (drawn up/down)')
ax.plot([], [], color=B, lw=3, label='light (t = iτ): the same speed, turned 90°, off the line (drawn sideways)')
ax.legend(fontsize=7.5, loc='upper center', frameon=False)
ax.set_title('Each zero’s velocity under heat, and under light', fontsize=9.5)
ax = axs[1]
w = lf['wick_check']; tau = [a['tau'] for a in w]
ax.plot(tau, [-a['zero'][1] for a in w], 'o', color=B, ms=8, label='zero of H_(iτ) from Φ (moved off the line)')
ax.plot(tau, [-a['predicted'][1] for a in w], '-', color=O, label='prediction τ·H″/H′ at the first zero')
h = lf['heat_check']
ax.plot([a['t'] for a in h], [lf['first_zero_z'] - a['zero'] for a in h], 's', color=G_, ms=7, label='heat: shift along the line (stays real)')
ax.set_xlabel('τ (light) or t (heat)'); ax.set_ylabel('displacement of the first zero (z units)')
ax.legend(fontsize=8)
ax.set_title(f'Checked from Riemann’s Φ: the same size, at a right angle\n(off-line part agrees to 1e-12 at τ = 0.002)', fontsize=9.5)
plt.tight_layout(); plt.savefig('figures/fig5-light-rotation.png', dpi=150, bbox_inches='tight'); plt.close()

# ---------------- fig6: Berry-Keating boost hyperbola and the flux-tube dimension
fig, axs = plt.subplots(1, 3, figsize=(15, 4.6))
ax = axs[0]
E = 6.0; lx = lp = 1.0
xx = np.linspace(lx, E/lp, 300)
ax.fill_between(xx, lp, E/xx, color='#c6dbef')
xs_ = np.linspace(0.4, 9, 400); ax.plot(xs_, E/xs_, color=B, lw=2, label='xp = E: a boost orbit')
ax.axvline(lx, color=G_, lw=0.8, ls='--'); ax.axhline(lp, color=G_, lw=0.8, ls='--')
ax.set_xlim(0, 9); ax.set_ylim(0, 9); ax.set_xlabel('x'); ax.set_ylabel('p')
ax.text(2.3, 1.5, 'area = E·η − E + l_x l_p\nη = log(E/(l_x l_p)), the rapidity span', fontsize=8.5)
ax.legend(fontsize=8); ax.set_title('Berry–Keating H = xp: states = energy × rapidity / 2π', fontsize=9.5)
ax = axs[1]
T = np.linspace(10, 120, 2000); ax.step(T, np.searchsorted(Z, T), where='post', color='k', lw=1, label='zeros counted')
ax.plot(T, T/(2*np.pi)*(np.log(T/(2*np.pi)) - 1) + 7/8, color=O, lw=2, label='Berry–Keating smooth count')
ax.set_xlabel('T'); ax.set_ylabel('N(T)'); ax.legend(fontsize=8)
ax.set_title(f'The smooth count fits all 5,000 zeros\n(mean residual {lf["bk_residual_mean"]:.3f}, max {lf["bk_residual_maxabs"]:.2f})', fontsize=9.5)
ax = axs[2]
D = lf['Deff']; Tc = np.array(D['T'])
ax.plot(Tc, D['measured'], 'o', color=B, ms=4, label='measured from the zeros')
ax.plot(Tc, D['smooth'], color=O, lw=2, label='log(T/2π) / (log(T/2π) − 1)')
ax.axhline(1, color=R, lw=1, ls='--'); ax.text(Tc[-1]*0.4, 1.03, 'tube: D = 1', color=R, fontsize=8.5)
ax.set_xscale('log'); ax.set_ylim(0.8, 2.5); ax.set_xlabel('height T'); ax.set_ylabel('effective dimension of the count')
ax.legend(fontsize=8); ax.set_title('Flux-tube reading: the zeros’ count is a tube (D → 1)\nwith one extra slowly closing log, the same log as the rapidity', fontsize=9.5)
plt.tight_layout(); plt.savefig('figures/fig6-boost-and-tube.png', dpi=150, bbox_inches='tight'); plt.close()

# ---------------- fig7: pair 3, the height detector across the roundness dial
fig, axs = plt.subplots(1, 3, figsize=(15, 4.8))
ax = axs[0]
u = np.linspace(0, 2, 400)
for p, c in [(2.5, '#fdae6b'), (3.0, O), (4.0, '#9ecae1'), (8.0, B)]:
    ax.plot(u, np.exp(-u**p), color=c, lw=1.6, label=f'p = {p:g}')
def Phi(x):
    tot = 0
    for n in range(1, 8):
        a = math.pi*n*n; tot += (2*a*a*math.exp(9*x) - 3*a*math.exp(5*x))*math.exp(-a*math.exp(4*x))
    return tot
ph = np.array([Phi(x) for x in u]); ax.plot(u, ph/ph[0], color='k', lw=1.4, ls='--', label='Riemann’s Φ (scaled)')
ax.set_xlabel('u'); ax.set_ylabel('kernel'); ax.legend(fontsize=8)
ax.set_title('Kernels on the roundness dial, e^(−|u|^p)', fontsize=9.5)
ax = axs[1]
for p, c in [(2.5, '#fdae6b'), (3.0, O), (4.0, '#9ecae1'), (8.0, B)]:
    K = p3['kernels'][str(p)]; ax.plot(K['xs'], K['h'], color=c, lw=1.2, label=f'p = {p:g}')
ax.axhline(1, color=R, lw=1.2); ax.text(0.5, 1.08, 's = 1: Laguerre (all zeros real)', color=R, fontsize=8)
for n in [2, 3, 4]: ax.axhline(1 - 1/n, color=G_, lw=0.6, ls=':')
for n in [2, 4]: ax.text(20.3, 1 - 1/n + 0.03, f'degree {n}', fontsize=7, color=G_)
ax.axhline(p3['xi_round4_hmax'], color='k', lw=0.8, ls='--'); ax.text(12.8, p3['xi_round4_hmax'] + 0.09, 'Ξ’s max 0.705 (round 4)', fontsize=7)
ax.set_ylim(-1.5, 2.2); ax.set_xlabel('x'); ax.set_ylabel('critical height h* = f f″/f′²')
ax.legend(fontsize=8, loc='lower left'); ax.set_title('The height detector: above 1 means a zero has left the line', fontsize=9.5)
ax = axs[2]
for p, c, mk_ in [(2.5, '#fdae6b', 'o'), (3.0, O, 's'), (4.0, '#9ecae1', 'D')]:
    K = p3['kernels'][str(p)]
    ax.plot(K['zeros'], np.zeros(len(K['zeros'])), mk_, color=c, ms=6, mfc='none')
    nz = p3['nonreal_zeros'].get(str(p), [])
    if nz:
        a = np.array(nz); ax.plot(a[:, 0], a[:, 1], mk_, color=c, ms=7, label=f'p = {p:g}: non-real zeros')
    for pk in p3['detector_peaks'].get(str(p), []):
        ax.axvline(pk, color=c, lw=0.8, ls=':')
ax.plot([], [], 'D', color='#9ecae1', mfc='none', label='p = 4: all real (Pólya)')
ax.set_xlim(0, 22); ax.set_ylim(-0.3, 6); ax.set_xlabel('Re x'); ax.set_ylabel('Im x')
ax.legend(fontsize=8, loc='upper left'); ax.set_title('Dotted: where the detector fires. Markers: the zeros\n(hollow on the axis = real, filled = non-real)', fontsize=9.5)
plt.tight_layout(); plt.savefig('figures/fig7-height-roundness.png', dpi=150, bbox_inches='tight'); plt.close()
print('figures done')
