"""DRAFT. Pictures for the Conley-Zehnder follow-up (figures 8-11). Run: python3 figures_cz.py (needs matplotlib)."""
import json, numpy as np, matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import cz
from cz_summary import margin

OR, BL, GR, RD, GY = '#e08a1e', '#2a6fdb', '#2e9e4f', '#d43b3b', '#888888'
CLS_COL = {'E-E': GR, 'E-H': '#9b59b6', 'E-IH': '#c39bd3', 'H-H': RD, 'H-IH': '#e67e22', 'IH-IH': '#a04000', 'CU': '#5d6d7e'}
D = json.load(open('cz.json')); KR = json.load(open('krein.json'))
ROB = [r for r in D if margin(r) > 0.01 and r['rel_pal'] < 1e-6]

def save(fig, name):
    fig.savefig('figures/' + name, dpi=130, bbox_inches='tight'); plt.close(fig); print(name)

def L(r, M=30):
    mu = r['muhat']; return np.array([mu[m - 1] - m * mu[0] for m in range(1, M + 1)])

def first(r):
    l = L(r, 60); nz = np.nonzero(l)[0]
    return int(nz[0]) + 1 if len(nz) else None

def shortest(k):
    return min((r for r in ROB if r['klass'] == k), key=lambda r: r['T'])

# ---- Figure 8: the iterate test on four orbits ----
ex = [next(r for r in D if r['cls'] == 'I.A1'), shortest('E-H'), shortest('H-H'), shortest('CU')]
names = ['the figure-eight (stable, E-E)', None, None, None]
fig, ax = plt.subplots(1, 2, figsize=(15, 5.2))
m = np.arange(1, 31)
for r, nm, mk in zip(ex, names, ['o', 's', '^', 'v']):
    lab = nm or '%s, %s (%s)' % (r['cls'], r['word'][:10], r['klass'])
    ax[0].plot(m, L(r), mk + '-', color=CLS_COL[r['klass']], ms=5, lw=1.2, label=lab)
ax[0].axhline(0, color='k', lw=0.6)
ax[0].set_xlabel('iterate m (laps of the orbit)')
ax[0].set_ylabel('mu(m-th lap) - m * mu(1 lap)')
ax[0].set_title('Conley-Zehnder index of the m-th lap, minus m times the first')
ax[0].set_ylim(-5, 29)
ax[0].text(17, -3.9, 'H-H and complex-unstable: exactly 0 for every m\n(hyperbolic pairs only add a fixed amount per lap)',
           fontsize=9, color=RD, ha='center')
ax[0].text(18, 6.5, 'any elliptic pair breaks the straight line;\nthe eight first shows it at lap 4', fontsize=9, color=GR)
ax[0].legend(fontsize=8, loc='upper left', bbox_to_anchor=(0, 0.88))
# right: same data split into one sawtooth per pair
f8 = ex[0]; mm = np.arange(1, 41)
for p, c, nm in zip(f8['pairs'], [OR, BL], ['pair 1', 'pair 2']):
    saw = 2 * cz.E(mm * p['theta'] / (2 * np.pi)) - 1 - mm * p['theta'] / np.pi
    ax[1].step(mm, saw, where='mid', color=c, lw=1.6,
               label='%s: multiplier angle %.3f, Krein sign %s' % (nm, p['phi'], '+' if p['s'] > 0 else '-'))
ax[1].plot(mm, cz.deviation(f8['pairs'], mm), 'k.', ms=6, label='their sum = index minus mean growth')
ax[1].axhline(0, color='k', lw=0.6); ax[1].set_ylim(-2.3, 2.3)
ax[1].set_xlabel('iterate m'); ax[1].set_ylabel('index minus m * mean index')
ax[1].set_title('The eight: one sawtooth per elliptic pair (Bott, Long)')
ax[1].legend(fontsize=8, loc='upper center', bbox_to_anchor=(0.5, -0.12), ncol=2)
save(fig, 'fig8-cz-iterates.png')

# ---- Figure 9: does it separate stable from unstable? ----
fig, ax = plt.subplots(1, 2, figsize=(15, 5))
pp = [r for r in ROB if r['klass'] in ('E-E', 'H-H', 'IH-IH', 'CU')]
pp.sort(key=lambda r: (['E-E', 'H-H', 'IH-IH', 'CU'].index(r['klass']), r['T']))
for i, r in enumerate(pp):
    c = CLS_COL[r['klass']]
    ax[0].plot(i, 1, 's', color=c, ms=7)
    f = first(r)
    ax[0].plot(i, 0, 'o' if f else 'x', color=c, ms=7, mew=2)
    if f: ax[0].text(i, -0.17, str(f), ha='center', fontsize=7)
ax[0].set_yticks([0, 1]); ax[0].set_yticklabels(['Conley-Zehnder laps:\nfirst lap that bends\n(x = never)', 'round one:\nfixed-point index\n(+1 for all)'])
ax[0].set_ylim(-0.5, 1.4); ax[0].set_xlabel('the 52 safe orbits whose fixed-point index is +1 at every lap')
for k in ('E-E', 'H-H', 'IH-IH', 'CU'): ax[0].plot([], [], 's', color=CLS_COL[k], label=k)
ax[0].legend(fontsize=8, ncol=4, loc='upper center')
ax[0].set_title('Same index, now told apart: only the 4 stable ones bend')
cls = ['E-E', 'E-H', 'E-IH', 'H-H', 'H-IH', 'IH-IH', 'CU']
bins = np.arange(1, 32)
fs = [first(r) for r in D if r['klass'] == 'E-E']
ax[1].hist([f for f in fs], bins=bins, color=GR, alpha=0.8, label='71 stable (E-E): lap where the line first bends')
ax[1].bar([1], [0])
nH = sum(1 for r in D if r['klass'] in ('H-H', 'H-IH', 'IH-IH', 'CU'))
ax[1].text(9, 25, '%d orbits with no elliptic pair: never bends (60 laps)\n218 with one elliptic pair: bends at lap 2, every time' % nH, fontsize=9)
ax[1].set_xlabel('lap'); ax[1].set_ylabel('orbits'); ax[1].legend(fontsize=8)
ax[1].set_title('All 692 catalog orbits: the integers alone recover the class')
save(fig, 'fig9-separation.png')

# ---- Figure 10: the dial reading ----
fig, ax = plt.subplots(1, 2, figsize=(13, 6))
th = np.linspace(0, 2 * np.pi, 400)
for a in ax:
    a.plot(np.cos(th), np.sin(th), color=GY, lw=1); a.set_aspect('equal'); a.axis('off')
    a.plot(-1, 0, 'ko', ms=6); a.text(-1.0, -0.16, 'x = infinity', fontsize=9, ha='center')
    a.plot(1, 0, 'k*', ms=14); a.text(1.22, -0.18, 'x = 0\n(the identity)', fontsize=9, ha='center')
p = f8['pairs'][1]
for k in range(1, 12):
    ang = k * p['theta']
    ax[0].annotate('', (np.cos(ang), np.sin(ang)), (np.cos(ang - p['theta']), np.sin(ang - p['theta'])),
                   arrowprops=dict(arrowstyle='->', color=BL, lw=1.2, connectionstyle='arc3,rad=0.25'))
    ax[0].text(1.13 * np.cos(ang), 1.13 * np.sin(ang), str(k), fontsize=8, ha='center', va='center')
ax[0].set_title('Elliptic pair (tangent law, d = +1): x = tan(m theta / 2)\nwinds round RP1; each pass through x = 0 adds 2 to the index', fontsize=10)
ax[0].text(0, -1.45, "the eight's pair 2 (theta / 2 pi = 0.298) passes x = 0 at laps 4, 7, 11: the index jumps there", ha='center', fontsize=9)
mu_h = 0.6
xs = np.tanh(np.arange(0, 8) * mu_h / 2)
angs = 2 * np.arctan(xs)
ax[1].plot(np.cos(angs), np.sin(angs), 'o-', color=RD, ms=6)
for k, a_ in enumerate(angs): ax[1].text(1.12 * np.cos(a_), 1.12 * np.sin(a_), str(k), fontsize=8, ha='center')
ax[1].set_title('Hyperbolic pair (velocity addition, d = -1): x = tanh(m mu / 2)\nleaves x = 0 and never comes back: no jumps', fontsize=10)
ax[1].text(0, -1.45, 'the index then grows by the same amount every lap', ha='center', fontsize=9)
save(fig, 'fig10-dial-pole.png')

# ---- Figure 11: Krein signs of the confirmed stable orbits ----
fig, ax = plt.subplots(1, 2, figsize=(14, 6))
ax[0].plot(np.cos(th), np.sin(th), color=GY, lw=1); ax[0].set_aspect('equal'); ax[0].axis('off')
for x in KR['orbits']:
    for ph, s in zip(x['phi'], x['s']):
        rr = 1.0 + (0.06 if x['krein'] == 'mixed' else -0.06)
        z = rr * np.exp(1j * ph)
        ax[0].plot(z.real, z.imag, 'o' if s > 0 else 'D', color=BL if s > 0 else OR, ms=7,
                   mec='k' if x['cls'] == 'I.A1' else 'none', mew=1.5)
f = [x for x in KR['orbits'] if x['cls'] == 'I.A1'][0]
ax[0].annotate('the eight', (np.cos(f['phi'][1]) * 1.06, np.sin(f['phi'][1]) * 1.06), (0.2, 1.3), arrowprops=dict(arrowstyle='->'))
ax[0].plot([], [], 'o', color=BL, label='Krein sign +'); ax[0].plot([], [], 'D', color=OR, label='Krein sign -')
ax[0].text(-1.25, -1.2, 'outer ring: mixed orbits (signs differ)\ninner ring: definite orbits (same sign)', fontsize=8)
ax[0].legend(fontsize=8, loc='lower right'); ax[0].set_xlim(-1.3, 1.3); ax[0].set_ylim(-1.35, 1.4)
ax[0].set_title('Multipliers e^(i phi) of the 19 confirmed stable orbits, upper half')
km = [x for x in KR['orbits'] if x['krein'] == 'mixed']; kd = [x for x in KR['orbits'] if x['krein'] == 'definite']
ax[1].scatter([x['collision_gap'] for x in kd], [x['edge_gap'] for x in kd], color=GR, s=40, label='definite (%d): cannot turn complex-unstable' % len(kd))
ax[1].scatter([x['collision_gap'] for x in km], [x['edge_gap'] for x in km], color=RD, marker='D', s=40, label='mixed (%d): can, by a Krein collision' % len(km))
ax[1].set_xlabel('gap between the two angles (distance to a collision)')
ax[1].set_ylabel('distance of the nearer angle to 0 or pi')
ax[1].set_title('Two ways a stable orbit can be lost (Krein 1950)')
ax[1].legend(fontsize=8)
save(fig, 'fig11-krein.png')
