"""DRAFT. Pictures for the three-body Conley-index test. Run: python3 figures.py (needs matplotlib)."""
import json, numpy as np, matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import cb3, section as S

OR, BL, GR, RD, GY = '#e08a1e', '#2a6fdb', '#2e9e4f', '#d43b3b', '#888888'
CLS_COL = {'E-E': GR, 'E-H': '#9b59b6', 'E-IH': '#c39bd3', 'H-H': RD, 'H-IH': '#e67e22', 'IH-IH': '#a04000', 'CU': '#5d6d7e'}
FL = [o for o in json.load(open('floquet.json')) if o['ok']]
CU = [o for o in json.load(open('curated.json')) if o['ok']]
TW = {(o['name'], o['which']): o for o in json.load(open('twins.json'))}
AN = json.load(open('analysis.json'))
v1, v2 = 0.347116888120318, 0.532724945387253
X0 = cb3.suvakov(v1, v2); E8 = cb3.energy(X0); T8 = 6.325913982936183
s8 = np.array([2, 0, v1, v1])

def save(fig, name):
    fig.savefig('figures/' + name, dpi=130, bbox_inches='tight'); plt.close(fig); print(name)

# ---- Figure 1: the section and the return map near the eight ----
fig, ax = plt.subplots(1, 3, figsize=(16, 5))
ts = np.linspace(0, T8, 1500); P = np.array([cb3.flow(X0, t) for t in ts])
ax[0].plot(P[:, 4], P[:, 5], color=GY, lw=1.4, label='the one curve all three bodies follow')
ax[0].plot([-1, 1, 0], [0, 0, 0], 'k--', lw=0.8)
for i, c in enumerate([BL, OR, GY]): ax[0].plot(*X0[2 * i:2 * i + 2], 'o', color=c, ms=9, mec='k', label='body %d at the start' % (i + 1))
ax[0].annotate('section: a line-up with body 3 between 1 and 2,\ntriangle area shrinking through zero', (0, 0), (-1.1, -0.75), arrowprops=dict(arrowstyle='->'), annotation_clip=False)
ax[0].set_aspect('equal'); ax[0].set_ylim(-0.85, 0.75); ax[0].legend(loc='upper right', fontsize=7); ax[0].set_title('The figure-eight (Moore 1993) and where the map looks')
# coordinates sketch
ax[1].plot([-1, 1], [0, 0], 'k-', lw=1)
for x, c, lab in [(-1, BL, '1'), (1, OR, '2'), (0.25, GY, '3')]:
    ax[1].plot(x, 0, 'o', color=c, ms=14, mec='k'); ax[1].text(x, 0.08, lab, ha='center')
ax[1].annotate('', (1, -0.15), (-1, -0.15), arrowprops=dict(arrowstyle='<->')); ax[1].text(0, -0.25, 'p = distance 1 to 2', ha='center')
ax[1].annotate('', (0.25, 0.25), (0, 0.25), arrowprops=dict(arrowstyle='<->')); ax[1].text(0.12, 0.3, 'q', ha='center')
ax[1].arrow(-1, 0, 0.35, 0, head_width=0.04, color=BL); ax[1].text(-0.8, -0.07, 'u1', color=BL)
ax[1].arrow(1, 0, 0.35, 0, head_width=0.04, color=OR); ax[1].text(1.15, -0.07, 'u2', color=OR)
ax[1].text(-1.2, -0.55, 'Four numbers (p, q, u1, u2) fix a point.\nThe sideways speeds follow from total momentum 0,\nangular momentum 0 and the energy.\nF sends a point to the next line-up of the same kind.', fontsize=9)
ax[1].set_xlim(-1.4, 1.6); ax[1].set_ylim(-0.75, 0.5); ax[1].axis('off'); ax[1].set_title('The four coordinates on the section')
# invariant tori slices near the eight
J8 = S.jac(s8, E8, 1); w, V = np.linalg.eig(J8)
for amp, c in [(2e-4, GR), (6e-4, BL), (1.2e-3, OR)]:
    x = s8 + amp * np.real(V[:, 0] + V[:, 2]); pts = []
    for _ in range(1500):
        x = S.F(x, E8, 1)[0]; pts.append(x - s8)
    pts = np.array(pts); ax[2].plot(pts[:, 1], pts[:, 2] - pts[:, 3], '.', ms=1.2, color=c)
ax[2].plot(0, 0, '*', color=GR, ms=14, mec='k', label='the eight: fixed point of F')
ax[2].set_xlabel('q (offset of body 3)'); ax[2].set_ylabel('u1 - u2'); ax[2].legend(loc='upper right')
ax[2].set_title('1,500 returns from three nearby starts:\nthey stay on rings (slices of invariant tori)')
save(fig, 'fig1-section-and-map.png')

# ---- Figure 2: Floquet classes and their index pairs ----
fig, ax = plt.subplots(1, 2, figsize=(15, 5.5), gridspec_kw=dict(width_ratios=[1.25, 1]))
classes = ['E-E', 'H-H', 'IH-IH', 'CU', 'E-H', 'E-IH', 'H-IH']
desc = {'E-E': 'both pairs on the unit circle\n(linearly stable)', 'H-H': 'both pairs real, positive', 'IH-IH': 'both pairs real, negative',
        'CU': 'a quartet off the circle\n(Krein collision)', 'E-H': 'one circle pair, one real pair', 'E-IH': 'circle pair + negative pair', 'H-IH': 'positive pair + negative pair'}
signs = {'E-E': (1, 1), 'H-H': (1, 1), 'IH-IH': (1, 1), 'CU': (1, 1), 'E-H': (-1, -1), 'E-IH': (1, -1), 'H-IH': (-1, 1)}
for i, c in enumerate(classes):
    y = len(classes) - i
    ax[0].add_patch(plt.Rectangle((0, y - 0.45), 1, 0.9, color=CLS_COL[c], alpha=0.25))
    ax[0].text(0.03, y, c, fontsize=12, weight='bold', va='center')
    ax[0].text(0.2, y, desc[c], fontsize=8.5, va='center')
    a, b = signs[c]
    ax[0].text(0.78, y, '%+d' % a, fontsize=14, va='center', color=GR if a > 0 else RD)
    ax[0].text(0.9, y, '%+d' % b, fontsize=14, va='center', color=GR if b > 0 else RD)
ax[0].add_patch(plt.Rectangle((0.72, 3.5), 0.27, 4.0, fill=False, ec=BL, lw=2.5))
ax[0].text(0.73, 7.65, 'i(F)   i(F²)', fontsize=10)
ax[0].text(0.0, 0.3, 'Blue box: four classes, one index pair. Only the top one is stable.', color=BL, fontsize=10)
ax[0].set_xlim(0, 1); ax[0].set_ylim(0, 8.1); ax[0].axis('off'); ax[0].set_title('Index of every iterate sees only the parity of hyperbolic pairs')
# eigenvalues of curated orbits in the complex plane (log-polar so big ones fit)
th = np.linspace(0, 2 * np.pi, 300); ax[1].plot(np.cos(th), np.sin(th), color=GY, lw=0.8)
for o in CU:
    lam = np.array(o['eig_abs']) * np.exp(1j * np.array(o['eig_arg']))
    z = (1 + np.log(np.abs(lam)) / 4) * np.exp(1j * np.angle(lam))
    ax[1].plot(z.real, z.imag, 'o', color=CLS_COL[o['klass']], ms=5, alpha=0.8)
for c in classes: ax[1].plot([], [], 'o', color=CLS_COL[c], label=c)
ax[1].set_aspect('equal'); ax[1].legend(fontsize=8, loc='upper left', ncol=2)
ax[1].set_title('Multipliers of 21 refined orbits\n(radius 1 + log|λ|/4; unit circle grey)')
save(fig, 'fig2-classes-and-index-pairs.png')

# ---- Figure 3: catalog tally ----
fig, ax = plt.subplots(1, 2, figsize=(14, 4.8))
rob = AN['classes_robust']; allc = AN['classes_all']
xs = np.arange(len(classes))
ax[0].bar(xs - 0.2, [allc.get(c, 0) for c in classes], 0.4, color=GY, label='all 692 (catalog precision)')
ax[0].bar(xs + 0.2, [rob.get(c, 0) for c in classes], 0.4, color=[CLS_COL[c] for c in classes], label='108 with a safe margin')
ax[0].set_xticks(xs); ax[0].set_xticklabels(classes); ax[0].legend(); ax[0].set_title('Floquet classes of the Li–Liao orbits')
pp = AN['plusplus']; lab = ['E-E', 'H-H', 'CU', 'IH-IH']
ax[1].pie([pp[c] for c in lab], labels=['%s: %d' % (c, pp[c]) for c in lab], colors=[CLS_COL[c] for c in lab], startangle=90,
          wedgeprops=dict(edgecolor='w'))
ax[1].set_title('The 52 safe-margin orbits whose index is +1 for every iterate:\nonly 4 are linearly stable')
save(fig, 'fig3-catalog-tally.png')

# ---- Figure 4: value-only count against the Jacobian ----
fig, ax = plt.subplots(1, 2, figsize=(14, 5))
good = [o for o in CU if o['residual'] < 1e-8]
for i, o in enumerate(good):
    v = o['deg_F'] if o['deg_F'] not in (0, None) else o['degrees']['F_cube']
    ax[0].plot(i, o['idx_jac'][0], 's', color=OR, ms=12, alpha=0.5)
    if v in (0, None): ax[0].plot(i, 0, 'x', color=RD, ms=10, mew=2)
    else: ax[0].plot(i, v, 'o', color=BL, ms=6)
ax[0].set_xticks(range(len(good))); ax[0].set_xticklabels(['%s\n%s' % (o['cls'], o['klass']) for o in good], rotation=90, fontsize=7)
ax[0].set_yticks([-1, 0, 1]); ax[0].set_ylabel('index of F^k')
ax[0].plot([], [], 's', color=OR, alpha=0.5, label='sign det(I − DF^k) (Jacobian)'); ax[0].plot([], [], 'o', color=BL, label='count from map values only')
ax[0].plot([], [], 'x', color=RD, label='count did not settle'); ax[0].legend(fontsize=8, loc='center right')
ax[0].set_title('Value-only index agrees wherever it settles (14 of 15)')
# triangulated boundary sketch (a 3-cube face cut into Kuhn simplices, shown in 2D analogue)
for k in range(5):
    ax[1].plot([k, k], [0, 4], color=GY, lw=0.6); ax[1].plot([0, 4], [k, k], color=GY, lw=0.6)
    ax[1].plot([0, 4 - k], [k, 4], color=GY, lw=0.4); ax[1].plot([k, 4], [0, 4 - k], color=GY, lw=0.4)
ax[1].plot([0, 4, 4, 0, 0], [0, 0, 4, 4, 0], 'k-', lw=2)
ax[1].text(2, -0.6, 'the box boundary is cut into small simplices; F^k − id is evaluated only at their corners', ha='center', fontsize=9)
ax[1].text(2, 4.4, 'count the simplices whose image surrounds one fixed direction, with orientation signs', ha='center', fontsize=9)
ax[1].set_xlim(-0.5, 4.5); ax[1].set_ylim(-1, 5); ax[1].axis('off'); ax[1].set_title('How the value-only count works (2D picture of the 4D method)')
save(fig, 'fig4-value-only-index.png')

# ---- Figure 5: twins ----
fig, ax = plt.subplots(1, 2, figsize=(15, 5))
names = ['I.A1', 'II.C73', 'II.B8', 'II.C192', 'II.C246', 'II.C66', 'II.C75']
for i, nm in enumerate(names):
    for j, w in enumerate(('cat', 'twin')):
        o = TW[(nm, w)]
        ax[0].scatter(i + 0.18 * (2 * j - 1), o['idx_jac'], s=40 + 25 * o['logmaxlam'], color=CLS_COL[o['klass']], edgecolor='k', zorder=3)
        ax[0].text(i + 0.18 * (2 * j - 1), o['idx_jac'] + (0.25 if j == 0 else -0.35), o['klass'], ha='center', fontsize=7)
    a, b = TW[(nm, 'cat')]['idx_jac'], TW[(nm, 'twin')]['idx_jac']
    ax[0].plot([i - 0.18, i + 0.18], [a, b], color=GY if a == b else BL, lw=1.5)
ax[0].axvline(4.5, color=GY, ls=':'); ax[0].text(5.5, 1.6, 'same index', ha='center', color=GY)
ax[0].text(2, 1.6, 'opposite index', ha='center', color=BL)
ax[0].set_xticks(range(len(names))); ax[0].set_xticklabels(names); ax[0].set_yticks([-1, 1]); ax[0].set_ylim(-1.9, 1.9)
ax[0].set_ylabel('index of F^k'); ax[0].set_title('Same-word twins (dot size ~ instability)')
a, b = TW[('I.A1', 'cat')], TW[('I.A1', 'twin')]
ta = [x[0] for x in a['tau'] if abs(x[0] - 2) < 0.1][0]; tb = [x[0] for x in b['tau'] if abs(x[0] - 2) < 0.1][0]
ax[1].axvline(2, color='k'); ax[1].axvspan(1.98, 2, color=GR, alpha=0.12); ax[1].axvspan(2, 2.02, color=RD, alpha=0.12)
ax[1].plot(ta, 0, 'o', color=GR, ms=14, mec='k'); ax[1].plot(tb, 0, 'o', color='#9b59b6', ms=14, mec='k')
ax[1].text(ta - 0.0015, 0.12, 'figure-eight\ntrace %.4f\nelliptic, index +1' % ta, ha='right'); ax[1].text(tb, -0.32, "Simó's H orbit\ntrace %.4f\nhyperbolic, index −1" % tb, ha='center')
ax[1].text(1.981, 0.42, 'stable side (2 − trace > 0)', color=GR); ax[1].text(2.001, 0.42, 'unstable side', color=RD)
ax[1].set_ylim(-0.5, 0.55); ax[1].set_xlim(1.98, 2.02); ax[1].set_yticks([]); ax[1].set_xlabel('trace of the fold pair, λ + 1/λ')
ax[1].set_title('The eight and its twin straddle the wall trace = 2\n(their other pair: −0.595 and −0.589, both elliptic)')
save(fig, 'fig5-twins.png')

# ---- Figure 6: isolation ----
fig, ax = plt.subplots(figsize=(8, 5))
cols = {'figure-eight I.A1 (E-E)': GR, "its twin, Simo's H orbit (E-H)": '#9b59b6', 'I.A4 (H-H)': RD, 'II.C24 (CU)': '#5d6d7e'}
for nm, v in AN['isolation'].items():
    ax.semilogx(v['radii'], v['frac_stay'], 'o-', color=cols[nm], label=nm.replace('Simo', 'Simó'), lw=2)
ax.set_xlabel('starting distance from the fixed point (section units)'); ax.set_ylabel('fraction of 24 starts that stay (300 returns)')
ax.set_title('Stable points hold their neighbours; fully unstable points lose them all.\nStaying sets that reach any box edge are why an elliptic point has no Conley index of its own.', fontsize=10)
ax.legend(fontsize=9); ax.set_ylim(-0.05, 1.1)
save(fig, 'fig6-isolation.png')

# ---- Figure 7: dial reading of one Floquet pair ----
fig, ax = plt.subplots(1, 3, figsize=(16, 4.8))
R = np.linspace(-0.6, 1.6, 400)
ax[0].plot(R, 4 * R * (1 - R), color=OR, lw=2, label='R ↦ 4R(1 − R)'); ax[0].plot(R, R, color=GY, ls='--')
ax[0].axvspan(0, 1, color=GR, alpha=0.1); ax[0].text(0.5, -1.6, 'elliptic\n0 < R < 1', ha='center', color=GR)
ax[0].text(-0.45, -1.6, 'hyperbolic\nR < 0', color=RD); ax[0].text(1.15, -1.6, 'inverse\nhyperbolic R > 1', color='#a04000')
for o in CU:
    for tr, ti in o['tau']:
        if abs(ti) < 1e-9 and -0.6 < (2 - tr) / 4 < 1.6:
            r_ = (2 - tr) / 4; ax[0].plot(r_, 4 * r_ * (1 - r_), 'o', color=BL, ms=4)
ax[0].set_ylim(-2.2, 1.3); ax[0].set_xlabel("Greene's residue R of one pair"); ax[0].set_ylabel('residue of the pair for F²')
ax[0].legend(); ax[0].set_title('Squaring the map: each pair\'s residue follows the logistic rule')
th = np.linspace(-np.pi + 0.05, np.pi - 0.05, 300)
ax[1].plot(th, np.tan(th / 2), color=GR, lw=2, label='elliptic: x = tan(θ/2), composes by\nx ∘ y = (x + y)/(1 − xy)  (k = 0, d = +1)')
mu = np.linspace(-3, 3, 300)
ax[1].plot(mu, np.tanh(mu / 2), color=RD, lw=2, label='hyperbolic: x = tanh(μ/2), composes by\nx ∘ y = (x + y)/(1 + xy)  (k = 0, d = −1)')
ax[1].set_xlabel('angle θ  or  log-multiplier μ of one pair'); ax[1].set_ylabel('dial coordinate x'); ax[1].set_ylim(-3, 3)
ax[1].legend(fontsize=8, loc='upper left'); ax[1].set_title('Iterating a pair = composing in a dial law')
ax[2].axis('off')
txt = ('Each Floquet pair carries its own curvature dial:\n\n'
       '   elliptic pair   →  d = +1  →  2 − trace > 0  →  factor +1\n'
       '   hyperbolic pair →  d = −1  →  2 − trace < 0  →  factor −1\n\n'
       'The fixed-point index is the PRODUCT of the factors.\n\n'
       'In the plane (one pair) the product is the pair:\n'
       'the index tells stable from unstable.\n\n'
       'In the reduced three-body map (two pairs)\n'
       '(+1)(+1) = (−1)(−1): the product forgets\n'
       'which pairs are stable.\n\n'
       'To keep the information you need the dials\n'
       'pair by pair, not their product.')
ax[2].text(0, 0.98, txt, va='top', fontsize=10.5, family='monospace', color=BL)
save(fig, 'fig7-dial-reading.png')
