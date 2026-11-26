import json, math, collections, os
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, Circle
D = os.path.dirname(os.path.abspath(__file__)); F = D + '/figures/'
rows = [x for x in json.load(open(D + '/bridges.json')) if not x['dup']]
COL = {'KNOWN': '#2b6cb0', 'FRAMING': '#8fb3d9', 'POSSIBLY-NEW': '#e0115f', 'NEGATIVE': '#c4c9d4'}
PRED = '#f08c00'; SHAPE = '#7b2cbf'
plt.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 11})

# ---------- Figure 1: the bridge map ----------
hubs = ['Dial k (one law a+b+kab)', 'Plane of laws (k, d)', 'Cube corners & weights', 'Rungs & height flow',
        'Swap rule (turn = d x area)', 'Shapes: simplex, ball, sphere', 'Chains & path sums']
hub_short = {'Dial k (one law a+b+kab)': 'dial k\na+b+kab', 'Plane of laws (k, d)': 'plane of\nlaws (k,d)',
             'Cube corners & weights': 'cube\ncorners', 'Rungs & height flow': 'rungs &\nheight flow',
             'Swap rule (turn = d x area)': 'swap rule\nturn = d·area', 'Shapes: simplex, ball, sphere': 'simplex,\nball, sphere',
             'Chains & path sums': 'chains &\npath sums'}
fields = ['Combinatorics & games', 'Logic, sets & computation', 'Number theory', 'Formal groups & algebra', 'Knots & topology',
          'Geometry', 'Analysis & dynamics', 'Probability & statistics', 'Information & entropy', 'Statistical physics & thermo',
          'Solitons & integrable', 'Fluids', 'Classical mechanics', 'Relativity & EM', 'Gravity & GR', 'Quantum mechanics',
          'Quantum gravity & holography', 'Chemistry', 'Biology', 'Economics, actuarial & engineering']
fig, ax = plt.subplots(figsize=(15, 15)); ax.set_aspect('equal'); ax.axis('off')
hp = {h: (2.9 * math.cos(2 * math.pi * i / 7 + .3), 2.9 * math.sin(2 * math.pi * i / 7 + .3)) for i, h in enumerate(hubs)}
fp = {f: (6.3 * math.cos(math.pi / 2 - 2 * math.pi * i / 20), 6.3 * math.sin(math.pi / 2 - 2 * math.pi * i / 20)) for i, f in enumerate(fields)}
E_ = collections.defaultdict(collections.Counter)
for x in rows: E_[(x['hub'], x['node'])][x['v']] += 1
for (h, f), c in E_.items():
    n = sum(c.values()); (x1, y1), (x2, y2) = hp[h], fp[f]
    newish = c['POSSIBLY-NEW']
    ax.plot([x1, x2], [y1, y2], color=COL['KNOWN'] if c['KNOWN'] >= n / 2 else '#7aa6d6', lw=0.4 + 0.9 * n ** 0.75, alpha=.38, zorder=1, solid_capstyle='round')
    if newish:
        ax.plot([x1, x2], [y1, y2], color=COL['POSSIBLY-NEW'], lw=0.9 + 0.5 * newish, alpha=.75, zorder=2)
cnt_h = collections.Counter(x['hub'] for x in rows); cnt_f = collections.Counter(x['node'] for x in rows)
newf = collections.Counter(x['node'] for x in rows if x['v'] == 'POSSIBLY-NEW')
for h, (x, y) in hp.items():
    r = 0.62 + cnt_h[h] ** .5 / 60
    ax.add_patch(Circle((x, y), r, fc='#0b1f4b', ec='white', lw=2, zorder=4))
    ax.text(x, y, hub_short[h] + f'\n{cnt_h[h]}', color='white', ha='center', va='center', fontsize=10.5, fontweight='bold', zorder=5)
for f, (x, y) in fp.items():
    lab = f.replace(' & ', ' &\n').replace(', actuarial', ',\nactuarial')
    ax.text(x * 1.03, y * 1.03, f'{lab}\n{cnt_f[f]} bridges' + (f' · {newf[f]} new?' if newf[f] else ''), ha='center', va='center', fontsize=10.5, zorder=5,
            bbox=dict(boxstyle='round,pad=.4', fc='white', ec=COL['POSSIBLY-NEW'] if newf[f] >= 5 else '#9aa6c0', lw=1.6 if newf[f] >= 5 else 1))
ax.add_patch(Circle((0, 0), 0.95, fc='#e0115f', ec='white', lw=2, zorder=4))
ax.text(0, 0, 'a ∘ₖ b\n= a+b+kab', color='white', ha='center', va='center', fontsize=12, fontweight='bold', zorder=5)
for h, (x, y) in hp.items(): ax.plot([0, x], [0, y], color='#0b1f4b', lw=2, zorder=1)
from matplotlib.lines import Line2D
ax.legend(handles=[Line2D([], [], color=COL['KNOWN'], lw=5, alpha=.6, label='mostly known bridges (width = how many)'),
                   Line2D([], [], color='#7aa6d6', lw=5, alpha=.6, label='mostly our own reading of known pieces'),
                   Line2D([], [], color=COL['POSSIBLY-NEW'], lw=3, label='carries a possibly-new bridge'),
                   Line2D([], [], color='#0b1f4b', lw=2, label='the seven parts of the framework (numbers = bridges)')],
          loc='lower center', bbox_to_anchor=(.5, -0.04), ncol=2, frameon=False, fontsize=12)
ax.set_xlim(-8.2, 8.2); ax.set_ylim(-8.0, 7.8)
ax.set_title('The bridge map: 739 bridges from 7 parts of the framework to 20 fields', fontsize=17, fontweight='bold', pad=4)
fig.savefig(F + 'fig1-bridge-map.png', dpi=110, bbox_inches='tight'); plt.close()

# ---------- Figure 2: verdicts by field ----------
order = sorted(fields, key=lambda f: cnt_f[f])
fig, ax = plt.subplots(figsize=(12, 9))
left = [0] * len(order)
for v, lab in [('KNOWN', 'known (source named)'), ('FRAMING', 'our reading, both sides known'), ('POSSIBLY-NEW', 'possibly new'), ('NEGATIVE', 'tested and failed')]:
    vals = [sum(1 for x in rows if x['node'] == f and x['v'] == v) for f in order]
    ax.barh(order, vals, left=left, color=COL[v], label=lab, edgecolor='white', hatch='//' if v == 'NEGATIVE' else None)
    for i, (l, w) in enumerate(zip(left, vals)):
        if w >= 4: ax.text(l + w / 2, i, str(w), ha='center', va='center', fontsize=9.5, color='white' if v in ('KNOWN', 'POSSIBLY-NEW') else '#0b1f4b')
    left = [a + b for a, b in zip(left, vals)]
ax.legend(loc='lower right', frameon=False); ax.spines[['top', 'right']].set_visible(False)
ax.set_xlabel('bridges'); ax.set_title('Every field, split by verdict', fontsize=15, fontweight='bold', loc='left')
fig.savefig(F + 'fig2-verdicts.png', dpi=110, bbox_inches='tight'); plt.close()

# ---------- Figure 3: the dial line ----------
fig, ax = plt.subplots(figsize=(15, 7.6)); ax.axis('off')
ax.plot([-2.6, 2.6], [0, 0], color='#0b1f4b', lw=3)
marks = {-2: ('k = −2', ['XOR / parity', 'piling-up lemma', 'Haldane map function ★']),
         -1: ('k = −1', ['union rule p+q−pq', 'inclusion–exclusion', 'stick-breaking', 'fermions (DPP)', 'AdS₃ entanglement gluing', 'holographic entropy test', "Kac's telegraph", 'actuarial UDD = Shapley ★']),
         -.5: ('k = −1/m', ['at most m copies', 'sphere marginal (k = −½)']),
         0: ('k = 0', ['plain addition', 'Poisson', 'Galilean boosts', 'Lavoisier sums']),
         1: ('k = 1', ['circle op = (1+a)(1+b)−1', 'K-theory, Todd genus', 'bosons, ζ(s)', 'Bayes on odds − 1', 'only size-consistent CI']),
         2: ('k = 2', ['Euler characteristic', 'law of cosines (c = 1)'])}
XP = {-2: -2.35, -1: -1.35, -.5: -0.35, 0: 0.55, 1: 1.45, 2: 2.35}
for k0, (lab, items) in marks.items():
    k = XP[k0]
    ax.plot(k, 0, 'o', ms=14, color='#e0115f', zorder=3); ax.text(k, -0.18, lab, ha='center', va='top', fontsize=13, fontweight='bold')
    for j, it in enumerate(items):
        ax.text(k, 0.25 + 0.27 * j, it, ha='center', va='bottom', fontsize=10.5, color=SHAPE if '★' in it else '#0b1f4b',
                bbox=dict(boxstyle='round,pad=.25', fc='white', ec=SHAPE if '★' in it else '#c7cfe0'))
ax.text(-2.6, -0.75, 'Movable dial values:', fontsize=12, fontweight='bold')
movers = ['k = λ  Sugeno λ-measure', 'k = 1 − q  Tsallis entropy', 'k = α  α-determinantal processes', 'k = −1/E_Planck  κ-Poincaré momenta',
          'k = −K  curvature in causal diamonds', 'k = −g*  one-loop running of G', 'k = −i/ħ  path integral, k = −i Feynman checkerboard']
for j, m in enumerate(movers): ax.text(-2.6 + (j % 3) * 1.85, -1.05 - 0.3 * (j // 3), m, fontsize=10.5, color='#2b6cb0')
ax.text(0, 2.75, 'Where bridges land on the dial k  (★ = found in the shape grouping, not named before)', ha='center', fontsize=15, fontweight='bold')
ax.set_xlim(-2.8, 2.9); ax.set_ylim(-1.9, 2.9)
fig.savefig(F + 'fig3-dial-line.png', dpi=110, bbox_inches='tight'); plt.close()

# ---------- Figure 4: the (k, d) plane with grouping finds ----------
import numpy as np
fig, ax = plt.subplots(figsize=(11, 10))
kk = np.linspace(-4.5, 4.5, 400); ax.plot(kk, kk ** 2 / 4, color='#0b1f4b', lw=2.5); ax.text(3.3, 3.25, 'k² = 4d\n(one double wall)', fontsize=11)
ax.fill_between(kk, kk ** 2 / 4, 5, color='#fde8ef', alpha=.6); ax.fill_between(kk, -5, kk ** 2 / 4, color='#e6eef9', alpha=.6)
ax.text(-4.2, 4.4, 'turning laws (no real wall): tan addition, waves', fontsize=11, color='#a1103f')
ax.text(-4.2, -4.6, 'two real walls: relativity, logistic, Kalman, Bayes pooling', fontsize=11, color='#1d4f8f')
pts = [((0, 0), 'addition'), ((1, 0), 'circle op'), ((0, -1), 'velocity addition\n(Kosambi)'), ((0, 1), 'tan addition'),
       ((2, 1), 'Euler characteristic,\nLavoisier c = 1'), ((-2, 0), 'XOR / Haldane'), ((-1, 0), 'union rule')]
OFF = {'tan addition': (8, 14), 'XOR / Haldane': (-30, -22), 'union rule': (-6, -22)}
for (k, d), lab in pts:
    ax.plot(k, d, 'o', ms=9, color='#2b6cb0'); ax.annotate(lab, (k, d), xytext=OFF.get(lab, (8, 6)), textcoords='offset points', fontsize=10.5)
ax.plot(-2, 1, '*', ms=22, color=SHAPE); ax.annotate('second-order kinetics = one-compartment synapse\n= AMH/Hamacher θ = 1 copula  (grouping find 1)',
                                                      (-2, 1), xytext=(-4.4, 2.3), fontsize=10.5, color=SHAPE, arrowprops=dict(arrowstyle='->', color=SHAPE))
x = np.linspace(-4.5, 1.5, 10); ax.plot(x, -1 - x, color=SHAPE, lw=2, ls='--')
ax.text(-4.4, 2.75, '', fontsize=1); ax.annotate('copula line k + d = −1 = genetics map-function line\n(after x = 2r; grouping find 2)', (0.9, -1.9), xytext=(0.6, -3.4), fontsize=10.5, color=SHAPE,
                                                 arrowprops=dict(arrowstyle='->', color=SHAPE))
ax.annotate('', xy=(2, 1), xytext=(-2, 1), arrowprops=dict(arrowstyle='<->', color=SHAPE, ls=':', lw=1.5)); ax.text(-1.2, 0.8, 'x → −x', color=SHAPE, fontsize=10)
ax.set_xlim(-4.6, 4.6); ax.set_ylim(-5, 5); ax.set_xlabel('k (corner weight)', fontsize=12); ax.set_ylabel('d (plane dial)', fontsize=12)
ax.axhline(0, color='#b8c2d6', lw=.8); ax.axvline(0, color='#b8c2d6', lw=.8); ax.spines[['top', 'right']].set_visible(False)
ax.set_title('The plane of laws (a+b+kab)/(1−d·ab): two bridges the grouping turned up', fontsize=14, fontweight='bold', loc='left')
fig.savefig(F + 'fig4-plane.png', dpi=110, bbox_inches='tight'); plt.close()

# ---------- Figure 5: the shape map ----------
fig, ax = plt.subplots(figsize=(15, 10.5)); ax.axis('off'); ax.set_xlim(0, 15); ax.set_ylim(-0.3, 10.2)
S = {'cross-polytope\n(p = 1)': (1.8, 8.6), 'sphere / ball\n(p = 2)': (7.5, 8.6), 'cube\n(p = ∞)': (13.2, 8.6), 'simplex': (7.5, 5.0),
     'permutohedron': (13.2, 5.0), 'elliptope': (1.8, 5.0), 'associahedron': (13.2, 1.4), 'torus / circle\nof angles': (1.8, 1.4), 'log axis': (7.5, 1.4)}
def link(a, b, txt, at, rad=0.0, col='#0b1f4b'):
    (x1, y1), (x2, y2) = S[a], S[b]
    ax.annotate('', xy=(x2, y2), xytext=(x1, y1), arrowprops=dict(arrowstyle='-|>', color=col, lw=2, shrinkA=34, shrinkB=34, connectionstyle=f'arc3,rad={rad}'))
    ax.text(at[0], at[1], txt, ha='center', va='center', fontsize=9.5, color=col, bbox=dict(fc='white', ec='none', pad=1.5))
link('cross-polytope\n(p = 1)', 'sphere / ball\n(p = 2)', 'roundness dial p', (4.65, 8.85))
link('sphere / ball\n(p = 2)', 'cube\n(p = ∞)', 'roundness dial p (→ tropical max)', (10.35, 8.85))
link('sphere / ball\n(p = 2)', 'simplex', 'square coordinates,\ndrop phases (Born rule)', (6.2, 6.8))
link('sphere / ball\n(p = 2)', 'simplex', 'Hilbert metric: ball curved −1\n→ simplex flat 0\n(1D: log odds = 2 × rapidity)', (8.95, 6.8), rad=-0.0001)
link('cube\n(p = ∞)', 'simplex', 'stick-breaking k = −1', (10.9, 7.3))
link('simplex', 'permutohedron', 'corner weights grow', (10.35, 5.25))
link('permutohedron', 'cube\n(p = ∞)', 'its vertices = paths\nup the cube', (14.25, 6.8))
link('cube\n(p = ∞)', 'associahedron', 'cluster coupling\nbends the cube', (14.45, 3.0), rad=-0.25)
link('elliptope', 'cube\n(p = ∞)', 'cut polytope (cube\ncorners) sits inside', (3.1, 7.55), rad=-0.12)
fams = {'cross-polytope\n(p = 1)': 'no science laws', 'sphere / ball\n(p = 2)': 'relativity, Kepler/Fock S³,\nhydrogen, spin, Berry, Tsirelson',
        'cube\n(p = ∞)': 'probability, solitons, CI, risks,\nfactorial designs, holography', 'simplex': 'Bayes, Gibbs, Born rule,\nreplicator, Fisher metric',
        'permutohedron': 'orderings, Plackett–Luce,\ngame cores, Shapley', 'associahedron': 'φ³ amplitudes, cluster\nalgebras, AdS₃ linking', 'elliptope': 'quantum correlations only',
        'torus / circle\nof angles': 'roots of unity, elliptic\ncurves, Lee–Yang circle', 'log axis': 'most of chemistry, biology,\nactuarial, economics'}
for s_, (x, y) in S.items():
    ax.add_patch(FancyBboxPatch((x - 1.1, y - .42), 2.2, .84, boxstyle='round,pad=.08', fc='#0b1f4b', ec='white'))
    ax.text(x, y, s_, color='white', ha='center', va='center', fontsize=11, fontweight='bold')
    ax.text(x, y - .58, fams[s_], ha='center', va='top', fontsize=9, color='#4a5778')
def dash(p1, p2, txt, at):
    ax.annotate('', xy=p2, xytext=p1, arrowprops=dict(arrowstyle='-', color=SHAPE, lw=2.2, ls=(0, (4, 3))))
    ax.text(at[0], at[1], txt, fontsize=9.5, color=SHAPE, ha='center', va='center', bbox=dict(fc='white', ec=SHAPE, boxstyle='round,pad=.25'))
dash((6.4, 1.4), (2.9, 1.4), 'binding zeros on the Lee–Yang circle (5)\nHückel ring levels = roots of unity (4)', (4.65, 1.95))
dash((8.6, 1.4), (12.1, 1.4), 'synthesis routes =\nbracketings (10)', (10.35, 1.95))
dash((8.3, 1.85), (12.4, 4.55), 'decrement splits: waterfall = vertices,\nShapley = centroid (7)', (10.6, 3.4))
dash((6.6, 1.85), (2.5, 4.55), 'effective dimension\nn²/‖C‖² lives here (9)', (4.4, 3.4))
dash((7.5, 1.85), (7.5, 4.5), 'selection on odds =\nflat Bayes move (11)', (7.5, 3.25))
ax.text(7.5, 10.0, 'The shape families and how they slide into each other', ha='center', fontsize=15, fontweight='bold')
ax.text(7.5, 9.6, 'purple dashes = bridges the grouping implies but no page named (numbers match the list on the page)', ha='center', fontsize=11.5, color=SHAPE)
fig.savefig(F + 'fig5-shape-map.png', dpi=110, bbox_inches='tight'); plt.close()

# ---------- Figure 6: predicted bridges ----------
fig, ax = plt.subplots(figsize=(15, 9.5)); ax.axis('off'); ax.set_xlim(0, 15); ax.set_ylim(0, 9.5)
P = [('P1', 'cube corners at k = −1', 'Hirota soliton corner weights\n(Cauchy matrix, det(I+tM))', 'holographic entropy test\n(MMI: corner sum at k = −1 ≤ 0)', 'log of soliton corner sums obeys\nmonogamy of mutual information'),
     ('P2', 'river picture + velocity disc', 'Schwarzschild ISCO factor\n1−6M/r = (1−2M/r)(1−(2v)²)', 'Kerr ISCO\n(spinning black hole)', 'Kerr factor splits the same way,\nwith v the frame-dragged speed'),
     ('P3', 'swap rule on the ball', 'Thomas–Wigner turn on the\nvelocity ball (d = −1)', 'qubit filtering measurements\non the Bloch ball', 'two filters leave a Bloch-vector\nturn = −(triangle area)'),
     ('P4', 'dial k = −1', 'fermions as k = −1\npoint process', 'AdS₃ / CFT entanglement\ngluing by ∘₋₁', 'free-fermion entanglement\nglues by ∘₋₁ at every cross-ratio'),
     ('P5', 'wall −1/k of one law', 'κ-Poincaré: k = −1/E_Planck,\nwall = maximum momentum', 'dimension-flow curves compose\nby velocity addition (ours)', 'κ-Minkowski spectral dimension\nis a two-wall curve, wall at E_Planck')]
for i, (lab, hubtxt, a, b, claim) in enumerate(P):
    y = 8.4 - i * 1.75
    for x, t, c in [(2.2, a, '#2b6cb0'), (8.2, b, '#2b6cb0')]:
        ax.add_patch(FancyBboxPatch((x - 1.9, y - .5), 3.8, 1.0, boxstyle='round,pad=.06', fc='white', ec=c, lw=1.6)); ax.text(x, y, t, ha='center', va='center', fontsize=10)
    ax.annotate('', xy=(6.25, y), xytext=(4.15, y), arrowprops=dict(arrowstyle='<->', color=PRED, lw=3, ls=(0, (4, 3))))
    ax.text(5.2, y + .62, hubtxt, ha='center', fontsize=9.5, color='#0b1f4b', fontweight='bold')
    ax.text(5.2, y - .05, lab, ha='center', va='top', fontsize=11, color='white', fontweight='bold', bbox=dict(boxstyle='round,pad=.25', fc=PRED, ec='none'))
    ax.text(10.5, y, 'test: ' + claim, va='center', fontsize=10.5, color='#7a4600')
ax.text(7.5, 9.35, 'Five bridges the map predicts, written down before anyone checks or looks them up', ha='center', fontsize=14, fontweight='bold')
fig.savefig(F + 'fig6-predictions.png', dpi=110, bbox_inches='tight'); plt.close()
print('ok6')
