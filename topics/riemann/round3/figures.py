"""Pictures for Riemann round 3. Colour key: blue = our reading, orange = already known, grey = reference,
red = where a road stops."""
import json, numpy as np, matplotlib
matplotlib.use('Agg'); import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, Circle
from scipy.special import ellipk
B, O, G, R, GR = '#1f77b4', '#ff7f0e', '#888888', '#d62728', '#2ca02c'
plt.rcParams.update({'font.size': 10, 'axes.spines.top': False, 'axes.spines.right': False})
rw = json.load(open('rewind.json')); rd = json.load(open('roads.json'))

# ---------------- fig0: map
fig, ax = plt.subplots(figsize=(12, 6.2)); ax.axis('off'); ax.set_xlim(0, 12); ax.set_ylim(0, 6.3)
ax.add_patch(FancyBboxPatch((4.55, 2.65), 2.9, 1.0, boxstyle='round,pad=0.1', fc='#eef', ec=B))
ax.text(6.0, 3.15, 'Riemann zeros\n(on the line Re s = 1/2?)', ha='center', va='center', fontsize=11)
items = [((0.15, 4.75), '1. Fused rewind (Brennan):\nrun the heat flow backward\nin fused (pair) coordinates', 'reversible, but RH is a sign,\nnot lost information', B),
         ((0.15, 2.6), '2. The Kerr turn: a shift\nturned by i (Newman–Janis)', 'turning by i swaps rewind\nand forward (Pólya, de Bruijn)', O),
         ((0.15, 0.45), '3. Barrier chain with\nrandom phases (tunnelling)', '= the Euler product\n(Bohr–Jessen, moments)', O),
         ((8.65, 4.75), '4. Lee–Yang far zone:\nthe margin that settled ℓp', 'would have to beat S(T),\nwhich is unbounded', R),
         ((8.65, 2.6), '5. Liebig knob (Clayton):\nproduct ↔ sum ↔ min', 'positive arithmetic,\ncannot make a zero', G),
         ((8.65, 0.45), 'Pointwise fused steps\n(rung flow E^t applied to Ξ)', 'fix every zero:\nE^t(0) = 0', G)]
for (x, y), txt, verdict, col in items:
    ax.add_patch(FancyBboxPatch((x, y), 3.2, 1.35, boxstyle='round,pad=0.08', fc='white', ec=col, lw=1.8))
    ax.text(x + 1.6, y + 0.9, txt, ha='center', va='center', fontsize=8.6)
    ax.text(x + 1.6, y + 0.25, verdict, ha='center', va='center', fontsize=8.4, color=col, weight='bold')
    ax.annotate('', xy=(4.55 if x < 4 else 7.45, 3.15), xytext=(x + 3.2 if x < 4 else x, y + 0.68),
                arrowprops=dict(arrowstyle='->', color=col, lw=1.2))
ax.set_title('Round 3: five new roads toward the Riemann zeros, and where each one ends', fontsize=11)
plt.savefig('figures/fig0-map.png', dpi=150, bbox_inches='tight'); plt.close()

# ---------------- fig1: fused rewind on a toy
T1 = rd['T1']; tr = T1['track']; t = np.array([x['t'] for x in tr]); tc = T1['t_collide']
fig, axs = plt.subplots(1, 3, figsize=(15, 4.6))
ax = axs[0]
allre = np.array([x['all_re'] for x in tr]); pos = allre[:, allre.shape[1] // 2:]
for j in range(pos.shape[1]): ax.plot(t, pos[:, j], color=G, lw=1)
re = np.array([x['re'] for x in tr]); im = np.array([x['im'] for x in tr])
real = im < 1e-7
ax.plot(t[real], re[real, 0], color=B, lw=2); ax.plot(t[real], re[real, 1], color=B, lw=2)
ax.plot(t[~real], re[~real, 0], color=R, lw=2, ls='--', label='the pair, now complex (real part)')
ax.axvline(tc, color=R, lw=0.8, ls=':'); ax.invert_xaxis()
ax.annotate('the two roots meet\n(speed → ∞ here)', (tc, 2.16), xytext=(0.03, 2.8), fontsize=8.5, arrowprops=dict(arrowstyle='->', color=R))
ax.set_xlabel('time t  (the rewind runs left to right →)'); ax.set_ylabel('root position z'); ax.set_ylim(0.5, 5.4)
ax.set_title('(a) Root coordinates: the rewind folds', fontsize=10); ax.legend(fontsize=8, loc='upper left')
ax = axs[1]
C = np.array([x['coef'] for x in tr]); C0 = C[np.argmin(abs(t))]
for j in range(0, C.shape[1] - 1, 2):
    if abs(C0[j]) > 0: ax.plot(t, C[:, j] / abs(C0[j]), lw=1.5, label=f'z^{j}')
ax.axvline(tc, color=R, lw=0.8, ls=':'); ax.invert_xaxis()
ax.text(tc, ax.get_ylim()[1] * 0.92, ' nothing happens here', color=R, fontsize=8.5)
ax.set_xlabel('time t'); ax.set_ylabel('coefficient / its value at t = 0')
ax.set_title('(b) Fused coordinates (coefficients): straight\nthrough, and back again to 1e-31', fontsize=10); ax.legend(fontsize=8)
ax = axs[2]
s = np.array([x['s'] for x in tr]); p = np.array([x['p'] for x in tr])
ss = np.linspace(4.04, 4.48, 100); ax.plot(ss, ss ** 2 / 4, color='k', lw=1.2, label='wall p = s²/4  (Δ = 0)')
ax.fill_between(ss, ss ** 2 / 4, 5.1, color='#fbe9e9', label='complex pair (Δ < 0)')
ax.set_ylim(4.2, 5.0)
ax.plot(s[t >= tc], p[t >= tc], color=B, lw=2.2); ax.plot(s[t <= tc], p[t <= tc], color=R, lw=2.2, ls='--')
for tv in [0.04, 0.0, tc, -0.06]:
    k = np.argmin(abs(t - tv)); ax.plot(s[k], p[k], 'o', color=B if tv > tc + 1e-4 else R)
    ax.annotate(f't = {tv:.3f}' + ('  (the roots meet)' if tv == tc else ''), (s[k], p[k]), textcoords='offset points', xytext=(8, -12), fontsize=8)
ax.set_xlabel('s = sum of the pair (the add part)'); ax.set_ylabel('p = product of the pair (the multiply part)')
ax.set_title('(c) The pair as one fused point (s, p):\nit slides through the wall, never breaks', fontsize=10); ax.legend(fontsize=8, loc='upper left')
fig.suptitle('Fusion keeps the rewind reversible. Toy: roots ±1, ±2, ±2.3, ±3.5, ±5 under the heat flow ∂H/∂t = −∂²H/∂z²', fontsize=11)
plt.tight_layout(); plt.savefig('figures/fig1-fused-rewind.png', dpi=150, bbox_inches='tight'); plt.close()

# ---------------- fig2: the rewind budget of the real zeros
pr = rw['pairs']; gm = np.array(pr['gamma_mid']); tau = np.array(pr['tau_isolated'])
fig, ax = plt.subplots(figsize=(11, 5))
ax.scatter(gm, tau, s=2, color=G, alpha=0.5, label='each neighbour pair alone: τ = Δγ²/2')
GLw = np.array([float(x) for x in json.load(open('zeros_lehmer.json'))['gamma']])
ax.scatter((GLw[1:] + GLw[:-1]) / 2, np.diff(GLw) ** 2 / 2, s=2, color=G, alpha=0.5)
ax.text(6200, 0.0035, 'zeros 6551–6870\n(window around\nLehmer’s pair)', fontsize=8, ha='center', color=G)
f = [x for x in rw['full'] if x['same_pair']]
ax.scatter([x['gamma'] for x in f], [-x['t_full'] for x in f], s=30, color=B, zorder=3, label='with every other zero pulling (full rewind)')
L = rw['lehmer']; ax.scatter([L['gamma']], [-L['t_full']], s=120, marker='*', color=O, zorder=4, label=f"Lehmer's pair (zeros 6709, 6710): {-L['t_full']:.1e}")
gg = np.linspace(20, 7500, 400); sp = 2 * np.pi / np.log(gg / (2 * np.pi))
ax.plot(gg, sp ** 2 / 2, color='k', ls='--', lw=1, label='typical pair: (mean spacing)²/2')
srt = np.argsort(gm); runmin = np.minimum.accumulate(tau[srt]); ax.step(gm[srt], runmin, color=B, lw=1.2, where='post', label='closest pair so far')
ax.set_yscale('log'); ax.set_ylim(4e-4, 40); ax.set_xlabel('height γ'); ax.set_ylabel('rewind time before this pair merges, |t|')
m = rw['full'][0]
ax.annotate(f"first merger among the first 5,000 zeros:\nt ≈ {m['t_full']:.5f}, at height {m['gamma']:.1f}", (m['gamma'], -m['t_full']), xytext=(2000, 6.5e-4),
            fontsize=8.5, arrowprops=dict(arrowstyle='->', color=B))
ax.set_title('How far back the backward heat flow can run before two real zeros merge. The budget shrinks with height,\nwhich is why it reaches 0 (Rodgers–Tao: Λ ≥ 0). RH says it is exactly 0, never below.', fontsize=10)
ax.legend(fontsize=8, loc='upper right'); ax.set_xlim(0, 7300)
plt.tight_layout(); plt.savefig('figures/fig2-rewind-budget.png', dpi=150, bbox_inches='tight'); plt.close()

# ---------------- fig3: the Kerr turn
T2 = rd['T2']
fig, axs = plt.subplots(1, 3, figsize=(15, 4.6))
ax = axs[0]; ax.set_aspect('equal'); ax.set_xlim(-1.6, 1.6); ax.set_ylim(-1.4, 1.4)
ax.axhline(0, color=G, lw=0.8); ax.axvline(0, color=G, lw=0.8)
ax.plot([-1, 1], [0.04, 0.04], color=R, lw=4, alpha=0.6); ax.plot([0.04, 0.04], [-1, 1], color=B, lw=4, alpha=0.6)
ax.add_patch(Circle((0, 0), 1, fill=False, ls=':', color='k'))
ax.text(0.55, -0.32, 'real smear\n= rewind', color=R, fontsize=9, ha='center')
ax.text(-1.05, 0.5, 'smear turned by i\n= forward\n(Kerr ring seen\nedge-on)', color=B, fontsize=9, ha='center')
ax.text(0, -1.3, 'z + a·cos θ  versus  z + i·a·cos θ', ha='center', fontsize=9)
ax.axis('off'); ax.set_title('(a) Same ring of shifts, turned by i', fontsize=10)
for ax, key, col, title in [(axs[1], 'A', B, '(b) Off-line pair 2 ± 0.3i: the i-turned ring\naverage puts it on the line'),
                            (axs[2], 'B', R, '(c) Close real pair 2, 2.3: the real ring\naverage pushes it off the line')]:
    for mode, c, ls, lab in [('imag', B, '-', 'ring turned by i: J₀(aD)'), ('real', R, '--', 'real ring: I₀(aD)')]:
        tr2 = [x for x in T2[f'{key}_{mode}_ring'] if x['a'] <= 0.6]; a = [x['a'] for x in tr2]; im2 = [x['im'] for x in tr2]
        ax.plot(a, im2, color=c, ls=ls, lw=2, label=lab)
    ax.set_xlabel('ring radius a'); ax.set_ylabel('|imaginary part| of the pair'); ax.legend(fontsize=8)
    ax.set_title(title, fontsize=10)
axs[1].axvline(T2['A_lands_at'], color=B, lw=0.8, ls=':'); axs[1].annotate(f"lands at a = {T2['A_lands_at']:.3f}\n(pair alone: √2·0.3 = 0.424)", (T2['A_lands_at'], 0.02), xytext=(0.36, 0.12), fontsize=8.5, arrowprops=dict(arrowstyle='->'))
axs[2].axvline(T2['B_leaves_at'], color=R, lw=0.8, ls=':'); axs[2].annotate(f"leaves at a = {T2['B_leaves_at']:.3f}\n(pair alone: 0.3/√2 = 0.212)", (T2['B_leaves_at'], 0.01), xytext=(0.3, 0.06), fontsize=8.5, arrowprops=dict(arrowstyle='->'))
fig.suptitle('The Kerr turn on zeros: averaging a function over a ring of shifts. Turned by i it is the forward heat flow\'s cousin; left real it is a rewind', fontsize=11)
plt.tight_layout(); plt.savefig('figures/fig3-kerr-turn.png', dpi=150, bbox_inches='tight'); plt.close()

# ---------------- fig4: barrier chain = Euler product
T3 = rd['T3']
def primes(n):
    s = np.ones(n + 1, bool); s[:2] = False
    for i in range(2, int(n ** 0.5) + 1):
        if s[i]: s[i * i::i] = False
    return np.nonzero(s)[0].astype(float)
PR = primes(10 ** 6)
fig, axs = plt.subplots(1, 3, figsize=(15, 4.6))
ax = axs[0]; sg = np.linspace(1.03, 3.2, 120)
ax.plot(sg, [np.exp(np.sum(np.log(2 / np.pi * ellipk(PR ** (-2 * x))))) for x in sg], color=O, lw=2, label='chain of barriers: ∏ (2/π)·K(p^(−σ))')
ms = sorted(T3['mean'], key=float)
ax.errorbar([float(k) for k in ms], [T3['mean'][k]['data'][0] for k in ms], yerr=[3 * T3['mean'][k]['data'][1] for k in ms], fmt='o', color=B, label='mean of |ζ(σ + it)| over 20,000 random t (bars: 3 s.e.)')
ax.set_xlabel('σ'); ax.set_ylabel('mean |ζ|'); ax.legend(fontsize=8); ax.set_title('(a) Legendre\'s K per barrier = mean size of ζ', fontsize=10)
ax = axs[1]; vc = T3['var_curve']
ax.plot(vc['sigma'], vc['chain'], color=O, lw=2, label='chain: ½ Σ Li₂(p^(−2σ))  (→ ∞ at σ = ½)')
vs = sorted(T3['var'], key=float)
ax.plot([float(k) for k in vs], [T3['var'][k]['data'] for k in vs], 'o', color=B, label='spread of log|ζ(σ + it)|')
ax.axvline(0.5, color=R, lw=1, ls=':'); ax.text(0.52, ax.get_ylim()[1] * 0.85 if ax.get_ylim()[1] < 5 else 1.4, 'critical\nline', color=R, fontsize=8.5)
ax.set_xlabel('σ'); ax.set_ylabel('variance of log|ζ|'); ax.legend(fontsize=8); ax.set_title('(b) The chain\'s spread blows up exactly at the line', fontsize=10)
ax = axs[2]
for k, c in [('1.0', G), ('0.75', B), ('0.6', R)]:
    ax.hist(T3['dist'][k], bins=70, density=True, histtype='step', color=c, lw=1.6, label=f'ζ data, σ = {k}')
    ax.hist(T3['mc'][k], bins=70, density=True, histtype='stepfilled', color=c, alpha=0.18)
ax.set_xlabel('log|ζ(σ + it)|'); ax.set_ylabel('density'); ax.legend(fontsize=8, title='shaded: random-phase chain', title_fontsize=8)
ax.set_title('(c) Same spread, data versus chain', fontsize=10)
fig.suptitle('A chain of barriers with random phases is the Euler product: the tunnelling law gives ζ\'s value distribution off the line (Bohr–Jessen)', fontsize=11)
plt.tight_layout(); plt.savefig('figures/fig4-barrier-chain.png', dpi=150, bbox_inches='tight'); plt.close()

# ---------------- fig5: far zone and S(T)
T4 = rd['T4']; g4 = np.array(T4['gamma']); S = np.array(T4['S'])
fig, ax = plt.subplots(figsize=(11, 4.4))
ax.scatter(g4, S, s=2, color=B, alpha=0.6, label='S at each zero: actual count minus the smooth count')
ax.axhline(0, color='k', lw=0.8)
srt = np.argsort(g4); ax.plot(g4[srt], np.maximum.accumulate(abs(S[srt])), color=R, lw=1.2, label='largest |S| so far')
ax.plot(g4[srt], -np.maximum.accumulate(abs(S[srt])), color=R, lw=1.2)
ax.set_xlabel('height γ'); ax.set_ylabel('S(γ)'); ax.legend(fontsize=8, loc='lower right')
ax.set_title('Where the ℓp far-zone argument would have to work for zeta: the zeros wander off the smooth count by S(T),\n'
             'which keeps growing without bound (Selberg). The ℓp margin had a fixed size; zeta\'s would need to beat this.', fontsize=10)
plt.tight_layout(); plt.savefig('figures/fig5-far-zone.png', dpi=150, bbox_inches='tight'); plt.close()

# ---------------- fig6: standard names
rows = [('Fused (pair) coordinates s, p of two zeros', 'Vieta / elementary symmetric functions; heat flow acts linearly on Taylor coefficients'),
        ('Rewind collision of two real zeros', 'Lehmer pair; lower bounds on Λ (Csordas–Smith–Varga 1994, Saouter–Gourdon–Demichel 2011)'),
        ('Zero motion under the flow', 'dz/dt = 2 Σ 1/(z_j − z_k) (Csordas–Smith–Varga; Rodgers–Tao 2020)'),
        ('Budget → 0 with height', 'de Bruijn–Newman constant Λ ≥ 0 (Rodgers–Tao); Λ ≤ 0.22 (Polymath 15)'),
        ('Smear turned by i = forward flow', 'Gaussian factor e^(λu²) as a universal factor (Pólya 1927, de Bruijn 1950, strip theorem)'),
        ('Ring average J₀(aD) / I₀(aD)', 'Laguerre–Pólya multiplier sequences; link to Newman–Janis not found in print'),
        ('Barrier chain with random phases', 'random Euler product (Bohr–Jessen 1930); Σ d_k(n)² n^(−2σ) = ∏ ₂F₁(k,k;1;p^(−2σ))'),
        ('Spread → ∞ at the line', 'Selberg\'s central limit theorem for log|ζ(½ + it)|'),
        ('Far-zone wobble', 'S(T) = (1/π) arg ζ(½ + iT); unbounded (Selberg)'),
        ('Pointwise fused step fixes zeros', 'composition with g, g(0) = 0, g invertible near 0')]
fig, ax = plt.subplots(figsize=(12.5, 4.6)); ax.axis('off')
tb = ax.table(cellText=rows, colLabels=['Our name in round 3', 'Standard name / where it lives'], loc='center', cellLoc='left', colWidths=[0.3, 0.7])
tb.auto_set_font_size(False); tb.set_fontsize(8.6); tb.scale(1, 1.45)
for (r_, c_), cell in tb.get_celld().items():
    if r_ == 0: cell.set_facecolor('#eef'); cell.set_text_props(weight='bold')
ax.set_title('Round 3 in standard terms', fontsize=11)
plt.savefig('figures/fig6-standard-names.png', dpi=150, bbox_inches='tight'); plt.close()
print('ok')
