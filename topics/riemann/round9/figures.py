import json, os
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from math import pi, log
from stringlib import String, theta, gram, true_S, mean_spacing

os.makedirs('figures', exist_ok=True)
gam = np.array(json.load(open('zeros.json'))['gamma'])
P = json.load(open('probe.json'))
L = json.load(open('learn.json'))
G = json.load(open('gram.json'))
D = json.load(open('drift.json'))
PR = json.load(open('predict.json'))
BLUE, ORANGE, GREY, RED, GREEN = '#2b6cb0', '#dd6b20', '#718096', '#c53030', '#2f855a'


def save(fig, name):
    fig.savefig(f'figures/{name}.png', dpi=130, bbox_inches='tight')
    plt.close(fig)


# fig0: the string and its parts
fig, ax = plt.subplots(figsize=(11, 5.2)); ax.axis('off')
def box(x, y, w, h, txt, c):
    ax.add_patch(plt.Rectangle((x, y), w, h, fc=c, ec='k', alpha=0.18))
    ax.text(x + w / 2, y + h / 2, txt, ha='center', va='center', fontsize=10, wrap=True)
box(0.02, 0.62, 0.30, 0.30, "Smooth law (archimedean place)\nN(E) = θ(E)/π + 1\nxp dynamics, extension angle α", BLUE)
box(0.35, 0.62, 0.30, 0.30, "Orbits (one p-adic place per prime)\nlength m·log p, weight p^(−m/2)/m\nsign −, no time reversal", ORANGE)
box(0.68, 0.62, 0.30, 0.30, "Levels E_n solve\nθ(E)/π + 1 + S_orbits(E) = n − ½\n(França–LeClair equation,\nLeClair–Mussardo scattering)", GREEN)
ax.annotate('', (0.35, 0.77), (0.32, 0.77), arrowprops=dict(arrowstyle='->'))
ax.annotate('', (0.68, 0.77), (0.65, 0.77), arrowprops=dict(arrowstyle='->'))
box(0.02, 0.08, 0.46, 0.44, "Dials on the orbits\nk: repeat weight k^(m−1) (primon-gas statistics)\nt: odds t^m   s: line σ the orbit sum is read on\nd: curvature 1 − u + d u²   q: free words\nr: orbit length r·log p   p: lattice (no orbits off p = 2)", GREY)
box(0.52, 0.08, 0.46, 0.44, "Settings (not dials)\ncutoff X and taper (a regulator)\nextension angle α, amplitude A, Maslov phase φ\ngrowth law of the dynamics\nfree weights or free frequencies (the trap)", GREY)
ax.set_title('Round 9: the string we score, and what can be turned', fontsize=13)
save(fig, 'fig0-round9-string')

# fig1: the seven dials
fig, axs = plt.subplots(2, 3, figsize=(13, 7))
names = dict(k='k homothety (repeat weight)', t='t odds', s='s height (line σ)', d='d curvature', q='q order', r='r ratio (orbit length)')
truth = dict(k=1.0, t=1.0, s=0.5, d=0.0, q=0.0, r=1.0)
for ax, dn in zip(axs.flat, ['k', 't', 's', 'd', 'q', 'r']):
    dd = P['dials'][dn]
    v = np.array([float(x) for x in dd]); o = np.argsort(v)
    tr = np.array([dd[x]['train']['rms'] for x in dd])[o]; te = np.array([dd[x]['test']['rms'] for x in dd])[o]
    ax.plot(v[o], tr, 'o-', c=BLUE, label='fit zeros 1–500'); ax.plot(v[o], te, 's-', c=ORANGE, label='held-out 501–5500')
    ax.axvline(truth[dn], c='k', ls=':', lw=1)
    ax.axhline(P['smooth']['test']['rms'], c=GREY, ls='--', lw=1, label='smooth law only (Gram)')
    ax.set_title(names[dn]); ax.set_ylabel('rms error (mean spacings)')
    if dn == 'r':
        ax.set_xscale('log')
axs[0, 0].legend(fontsize=8)
fig.suptitle('The seven dials on the string (p only acts on lattices: at p = 2 it adds L(χ₋₄)\'s levels and the error is off the scale)', fontsize=11)
fig.tight_layout(); save(fig, 'fig1-dials')

# fig2: settings
fig, axs = plt.subplots(1, 4, figsize=(16, 3.8))
grid = P['cutoff_grid']
for taper, c in [('sharp', BLUE), ('fejer', GREEN), ('gauss', ORANGE)]:
    ks = [k for k in grid if k.startswith(taper)]
    X = np.array([float(k.split('|')[1]) for k in ks])
    o = np.argsort(X)
    axs[0].plot(X[o], np.array([grid[k]['test']['rms'] for k in ks])[o], 'o-', c=c, label=taper)
axs[0].set_xscale('log'); axs[0].set_yscale('log'); axs[0].set_xlabel('cutoff X'); axs[0].legend(); axs[0].set_title('cutoff: more orbits is always better')
for ax, nm, tl in [(axs[1], 'alpha', 'extension angle α'), (axs[2], 'A', 'amplitude A (A = 2: time reversal)'), (axs[3], 'phi', 'Maslov phase φ (π: sign flip)')]:
    sw = P['setting_sweeps'][nm]
    v = np.array([float(x) for x in sw]); o = np.argsort(v)
    ax.plot(v[o], np.array([sw[x]['train']['rms'] for x in sw])[o], 'o-', c=BLUE, label='fit')
    ax.plot(v[o], np.array([sw[x]['test']['rms'] for x in sw])[o], 's-', c=ORANGE, label='held-out')
    ax.set_title(tl)
axs[1].legend(); axs[0].set_ylabel('held-out rms (mean spacings)')
fig.tight_layout(); save(fig, 'fig2-settings')

# fig3: regulator drift
fig, ax = plt.subplots(figsize=(7, 4.2))
for taper, c, mk in [('sharp', BLUE, 'o'), ('gauss', ORANGE, 's')]:
    rows = [r for r in D if r['taper'] == taper]
    X = [r['X'] for r in rows]
    ax.plot(X, [r['s']['best_train'] for r in rows], mk + '-', c=c, label=f'{taper} taper: best line σ (fit zeros)')
ax.axhline(0.5, c='k', ls=':', label='σ = 1/2')
ax.set_xscale('log'); ax.set_xlabel('cutoff X (orbits kept)'); ax.set_ylabel('best height dial s')
ax.set_title('The best line drifts back to 1/2 as orbits are added:\nthe offset is the regulator, not the string'); ax.legend(fontsize=8)
save(fig, 'fig3-regulator-drift')

# fig4: what the zeros teach
fig, axs = plt.subplots(1, 2, figsize=(13, 4.2))
fw = L['free_weights']
n = np.array(fw['n']); r = np.array(fw['fitted']) / np.array(fw['truth'])
axs[0].plot(n, r, 'o', c=BLUE); axs[0].axhline(1, c='k', ls=':')
axs[0].set_xscale('log'); axs[0].set_ylim(0, 2); axs[0].set_xlabel('prime power n = p^m')
axs[0].set_ylabel('fitted weight / (p^(−m/2)/m)'); axs[0].set_title('Weights learned from the first 500 zeros match p^(−m/2)/m')
ell = np.array(L['ell']); F = np.array(L['landau_all'])
axs[1].plot(ell, F, c=GREY, lw=0.8)
for v in [log(x) for x in [2, 3, 4, 5, 7, 8, 9, 11, 13, 16, 17, 19, 23, 25, 27, 29, 31, 32, 37, 41, 43, 47, 49, 53, 59, 61]]:
    axs[1].axvline(v, c=ORANGE, lw=0.6, alpha=0.6)
axs[1].set_xlabel('length ℓ'); axs[1].set_title('Σ over 5500 zeros of cos(γℓ): dips at ℓ = log p^m (Landau 1911)')
fig.tight_layout(); save(fig, 'fig4-zeros-teach-primes')

# fig5: the inverse-spectral trap
fig, ax = plt.subplots(figsize=(7.5, 4.5))
tp = L['trap']
K = [t['params'] for t in tp]
ax.plot(K, [t['train']['rms'] for t in tp], 'o-', c=BLUE, label='free frequencies: fit zeros')
ax.plot(K, [t['test']['rms'] for t in tp], 's-', c=ORANGE, label='free frequencies: held-out zeros')
ax.plot([2], [P['cutoff_grid']['gauss|200']['test']['rms']], '*', ms=16, c=GREEN, label='prime string (2 settings): held-out')
ax.set_xscale('log'); ax.set_yscale('log'); ax.set_xlabel('number of adjustable settings'); ax.set_ylabel('rms error (mean spacings)')
ax.set_title('The inverse-spectral trap: more settings fit the 500,\nand predict the next 5000 worse'); ax.legend(fontsize=8)
save(fig, 'fig5-inverse-spectral-trap')

# fig6: Gram's law
fig, axs = plt.subplots(1, 3, figsize=(16, 4.2))
t = np.linspace(255, 330, 6000)
axs[0].plot(t, true_S(t, gam), c=BLUE, lw=1, label='S(t)')
st = String(X=200.0, taper='gauss'); axs[0].plot(t, st.S(t), c=ORANGE, lw=1, label='prime string S (X = 200)')
ms = np.arange(int(np.floor(theta(255) / pi)), int(np.floor(theta(330) / pi)) + 1)
for m in ms:
    axs[0].axvline(gram(m), c=GREY, lw=0.5)
for z in gam[(gam > 255) & (gam < 330)]:
    axs[0].plot(z, -1.25, '|', c='k', ms=8)
axs[0].set_title("Gram's first failure (Hutchinson 1925): an empty\ninterval [g125, g126) next to one with two zeros")
axs[0].legend(fontsize=8); axs[0].set_xlabel('t')
hist = G['hist']; tot = sum(hist.values())
k = np.arange(4); obs = [hist.get(str(i), 0) / tot for i in k]
axs[1].bar(k - 0.2, obs, 0.4, color=BLUE, label='first 5500 zeros')
axs[1].bar(k[:3] + 0.2, [0.17022, 0.66143, 0.16649], 0.4, color=ORANGE, label="GUE limit (Fujii's conjecture)")
axs[1].bar(k + 0.6, [0.13708, 0.72763, 0.13351, 0.00178], 0.2, color=GREY, label='first 10⁹ (Brent, van de Lune et al.)')
axs[1].set_xticks(k); axs[1].set_xlabel('zeros in a Gram interval'); axs[1].legend(fontsize=7); axs[1].set_title('Where the law sits on its way to 66%')
b = G['law_by_height']
axs[2].plot([x['t_from'] for x in b], [x['rate'] for x in b], 'o-', c=BLUE)
axs[2].axhline(0.66143, c=ORANGE, ls='--', label='GUE limit 0.661'); axs[2].axhline(0.7276, c=GREY, ls=':', label='at 10⁹ zeros: 0.728')
axs[2].set_xlabel('height t'); axs[2].set_ylabel("share of intervals obeying Gram's law"); axs[2].legend(fontsize=8)
axs[2].set_title('It falls slowly with height')
fig.tight_layout(); save(fig, 'fig6-gram-law')

# fig7: next-zero prediction
fig, axs = plt.subplots(1, 2, figsize=(13, 4.2))
te = slice(500, len(gam))
ms_ = mean_spacing(gam[te])
for lab, c, nm in [('G', GREY, 'Gram midpoint (smooth law)'), ('L', RED, 'last zero + mean gap'),
                   ('PL', GREEN, 'prime string anchored on last zero'), ('P', ORANGE, 'prime string (X = 200)')]:
    e = (np.array(PR[lab], float)[te] - gam[te]) / ms_
    axs[0].hist(e[np.isfinite(e)], bins=np.linspace(-1, 1, 161), histtype='step', color=c, lw=1.4,
                label=f"{nm}: rms {PR[lab + '_stats']['rms'] if (lab + '_stats') in PR else np.sqrt(np.nanmean(e**2)):.3f}")
axs[0].set_yscale('log'); axs[0].set_xlabel('error (mean spacings), held-out zeros 501–5500'); axs[0].legend(fontsize=8)
axs[0].set_title('The prime string places every held-out zero; Gram\'s point is off by a quarter gap')
labs = ['G', 'P', 'PL']
x = np.arange(3)
axs[1].bar(x - 0.2, [PR[f'{l}_err_lawful'] for l in labs], 0.4, color=BLUE, label="intervals obeying Gram's law")
axs[1].bar(x + 0.2, [PR[f'{l}_err_lawless'] for l in labs], 0.4, color=RED, label='intervals breaking it')
axs[1].set_xticks(x); axs[1].set_xticklabels(['Gram midpoint', 'prime string', 'string, anchored']); axs[1].set_ylabel('mean |error| (spacings)')
axs[1].legend(fontsize=8); axs[1].set_title('Errors double where Gram\'s law breaks, for every rule')
fig.tight_layout(); save(fig, 'fig7-next-zero')

# fig8: string-theory checks
sto = P['string_theory']
labels, vals, cols = [], [], []
labels.append('truth (gauss X=200)'); vals.append(P['cutoff_grid']['gauss|200']['test']['rms']); cols.append(GREEN)
for p_, r_ in sto['drop_place'].items():
    labels.append(f'drop the {p_}-adic place'); vals.append(r_['test']['rms']); cols.append(ORANGE)
for lab, v in sto['statistics'].items():
    labels.append(f'primon gas: {lab}'); vals.append(v); cols.append(BLUE)
labels.append("sign flipped (φ = π)"); vals.append(sto['sign_flip']); cols.append(RED)
labels.append('time reversal (A = 2)'); vals.append(sto['time_reversal']); cols.append(RED)
labels.append('smooth law only (no orbits)'); vals.append(P['smooth']['test']['rms']); cols.append(GREY)
fig, ax = plt.subplots(figsize=(9, 6.5))
y = np.arange(len(labels))
ax.barh(y, vals, color=cols); ax.set_yticks(y); ax.set_yticklabels(labels, fontsize=9); ax.invert_yaxis()
ax.set_xscale('log'); ax.set_xlabel('held-out rms error (mean spacings)')
ax.set_title('String-theory readings as settings: each place, Bose statistics,\nthe absorption sign and broken time reversal are all needed')
save(fig, 'fig8-string-theory')

# fig9: the quantisation picture near Gram's first failure
fig, ax = plt.subplots(figsize=(9, 4.5))
t = np.linspace(270, 300, 4000)
N_true = np.searchsorted(gam, t, side='right')
ax.step(t, N_true, where='post', c='k', lw=1, label='true count N(t)')
ax.plot(t, theta(t) / pi + 1, c=GREY, label='smooth law θ/π + 1')
ax.plot(t, st.F(t), c=ORANGE, label='string F = smooth + prime orbits')
for n_ in range(int(N_true.min()), int(N_true.max()) + 2):
    ax.axhline(n_ - 0.5, c=BLUE, lw=0.4, alpha=0.5)
ax.set_xlabel('E'); ax.set_ylabel('count'); ax.legend(fontsize=8)
ax.set_title('Levels sit where F crosses n − ½: the orbits bend the smooth law onto the staircase')
save(fig, 'fig9-quantisation')
print('figures written')
