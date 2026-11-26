"""Pictures for the dial audit (DRAFT). python3 figures.py -> figures/*.png"""
import os, math, itertools
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch, Polygon

os.makedirs('figures', exist_ok=True)
GREEN, AMBER, BLUE, GREY, RED = '#2e7d32', '#ef8f00', '#1f5fae', '#8a8a8a', '#b23b3b'
plt.rcParams.update({'font.size': 11, 'axes.spines.top': False, 'axes.spines.right': False})


def save(fig, name):
    fig.savefig('figures/' + name, dpi=150, bbox_inches='tight')
    plt.close(fig)


# ---- fig 1: the verdict board
fig, ax = plt.subplots(figsize=(14, 7.4))
ax.set_xlim(0, 14.4); ax.set_ylim(-0.9, 7.9); ax.axis('off')
boxes = {
    'k': (1.5, 5.5, GREEN, 'k  homothety', 'fundamental\n(its size is a unit)'),
    'd': (5.5, 5.5, GREEN, 'd  curvature dial', 'fundamental\n(d/k² is unit-free)'),
    't': (9.5, 5.5, GREEN, 't  odds', 'fundamental input\n(one per element)'),
    'r': (13.0, 5.5, AMBER, 'r  similarity ratio', 'partly derived\n= odds + one zoom'),
    's': (1.5, 2.3, GREEN, 's  height', 'fundamental\n(along an axis)'),
    'p': (5.5, 2.3, GREEN, 'p  roundness', 'fundamental\n(across axes)'),
    'q': (9.5, 2.3, AMBER, 'q  order', 'partly derived: d × area,\nunless put in by hand'),
}
for key, (cx, cy, col, title, sub) in boxes.items():
    ax.add_patch(FancyBboxPatch((cx - 1.3, cy - 0.8), 2.6, 1.6, boxstyle='round,pad=0.05', fc=col, alpha=0.13, ec=col, lw=2))
    ax.text(cx, cy + 0.4, title, ha='center', va='center', fontsize=12.5, weight='bold', color=col)
    ax.text(cx, cy - 0.25, sub, ha='center', va='center', fontsize=9.6)


def link(k1, k2, label, lx, ly, col=GREY, rad=0.0, style='<->'):
    x1, y1 = boxes[k1][:2]; x2, y2 = boxes[k2][:2]
    ar = FancyArrowPatch((x1, y1), (x2, y2), arrowstyle=style, mutation_scale=14, color=col, lw=1.6,
                         connectionstyle='arc3,rad=%s' % rad, shrinkA=62, shrinkB=62)
    ax.add_patch(ar)
    ax.text(lx, ly, label, ha='center', va='center', fontsize=9, color=col, bbox=dict(fc='white', ec='none', pad=1))


link('k', 'd', 'unit change: k→λk, d→λ²d;\none unit-free number d/k²', 3.5, 4.5, BLUE)
link('k', 's', 'rung s of ∘ₖ has\nhomothety reading s·k', 2.6, 3.45, BLUE)
link('s', 'p', 'same operation:\np = 2ˢ on the Pythagorean sum', 3.5, 1.05, BLUE)
link('q', 'd', 'moves built from the law:\norder turn = d × area', 7.45, 3.9, AMBER, style='->')
link('r', 't', 'its piece-to-piece\npart is odds', 11.25, 4.4, AMBER, style='->')
link('t', 'k', 'with odds on every element, k is a common odds scale', 7.0, 7.3, BLUE, rad=0.3)
ax.text(7.2, 0.0, 'Rank tests, fixed ruler:  (k, d) rank 2   ·   (k, d, s) rank 3   ·   (s, p) rank 2   ·   r ↔ (D, n−1 weights) full rank',
        ha='center', fontsize=10.5)
ax.text(7.2, -0.45, 'Up to a change of unit: k and d keep one number between them (d/k²) and their signs; s, p, t, r and q are unit-free.',
        ha='center', fontsize=10.5)
ax.text(7.2, -0.85, 'Green: no other dial or combination reaches what it reaches.   Amber: part of what it does is another dial.',
        ha='center', fontsize=10, color=GREY)
ax.set_title('Figure 1. The seven dials after the audit', fontsize=14, weight='bold', loc='left')
save(fig, 'fig1-verdict-board.png')

# ---- fig 2: the law space in generator coordinates
sv = np.linspace(-1.2, 2.4, 400)
c1, c2, c3 = sv, sv * (sv - 1) / 4, -sv * (sv - 1) / 12
fig, axs = plt.subplots(1, 2, figsize=(12, 5))
ax = axs[0]
K, Dg = np.meshgrid(np.linspace(-1.5, 2.6, 2), np.linspace(-0.6, 1.2, 2))
ax.fill_between([-1.5, 2.6], -0.6, 1.2, color=GREEN, alpha=0.07)
ax.text(-1.4, 1.08, 'every (k, d) law is a point here\n(generator 1 + kx + dx², cubic term 0)', fontsize=9.5, color=GREEN, va='top')
ax.plot(c1, c2, color=BLUE, lw=2.4, label='rung s of ∘₁ (shadow on the plane)')
for sp_, lab in [(0, 's = 0  addition'), (1, 's = 1  ∘₁'), (2, 's = 2'), (0.5, 's = ½'), (-1, 's = −1')]:
    ax.plot(sp_, sp_ * (sp_ - 1) / 4, 'o', color=BLUE if sp_ not in (0, 1) else 'k', ms=7)
    ax.annotate(lab, (sp_, sp_ * (sp_ - 1) / 4), xytext=(-8, -18) if sp_ == 0.5 else (6, 8), textcoords='offset points', fontsize=9.5)
ax.axhline(0, color=GREY, lw=0.8); ax.axvline(0, color=GREY, lw=0.8)
ax.set_xlabel('c₁ = homothety reading k'); ax.set_ylabel('c₂ = curvature dial d')
ax.set_title('(a) shadow on the (k, d) plane', loc='left')
ax.legend(loc='lower right', fontsize=9)
ax = axs[1]
ax.plot(c1, c3, color=BLUE, lw=2.4)
ax.axhline(0, color=GREEN, lw=3, alpha=0.5)
ax.text(1.6, 0.012, 'the whole (k, d) plane\nlives on this line (c₃ = 0)', color=GREEN, fontsize=9.5)
for sp_ in (0, 1, 2, 0.5, -1):
    ax.plot(sp_, -sp_ * (sp_ - 1) / 12, 'o', color='k' if sp_ in (0, 1) else BLUE, ms=7)
ax.annotate('touches the plane only\nat s = 0 and s = 1', (0.5, 0.0), xytext=(-0.9, -0.28), fontsize=9.5,
            arrowprops=dict(arrowstyle='->', color=GREY))
ax.set_xlabel('c₁ = s·k   (here k = 1)'); ax.set_ylabel('c₃ = −k³ s(s−1)/12')
ax.set_title('(b) the cubic coordinate: the height leaves the plane', loc='left')
fig.suptitle('Figure 2. A law is its generator P(x) = 1 + c₁x + c₂x² + c₃x³ + …; the height s draws a curve the plane cannot reach',
             fontsize=12, weight='bold', x=0.02, ha='left')
fig.tight_layout()
save(fig, 'fig2-law-space.png')

# ---- fig 3: units
fig, ax = plt.subplots(figsize=(8.5, 6.4))
kk = np.linspace(-2.2, 2.2, 300)
ax.plot(kk, kk**2 / 4, color='k', lw=1.4)
ax.text(1.55, 0.8, 'Δ = 0\n(d = k²/4)', fontsize=9.5)
for ratio, col in [(-0.5, BLUE), (0.0, BLUE), (0.15, BLUE), (0.25, 'k'), (0.6, BLUE)]:
    for sgn in (1, -1):
        kr = np.linspace(0, 2.2, 100) * sgn
        ax.plot(kr, ratio * kr**2, color=col, lw=1.2, alpha=0.6)
    ax.text(2.22, ratio * 2.2**2, 'd/k² = %g' % ratio, fontsize=8.5, va='center', color=col)
lams = [0.5, 1, 1.6, 2.0]
for L in lams:
    ax.plot(L * 1.0, 0.6 * L**2, 's', color=AMBER, ms=8)
ax.annotate('one law, four units:\n(k, d) → (λk, λ²d)', (1.6, 0.6 * 1.6**2), xytext=(-2.1, 2.3), fontsize=10, color=AMBER,
            arrowprops=dict(arrowstyle='->', color=AMBER))
for dd, lab in [(-1, 'Einstein  d = −1/c²'), (0, 'Galileo  d = 0'), (1, 'rotations  d > 0')]:
    ax.plot(0, dd * 0.9, 'o', color=RED, ms=8)
    ax.text(0.08, dd * 0.9 - 0.12, lab, color=RED, fontsize=9.5)
ax.set_xlim(-2.3, 2.9); ax.set_ylim(-1.3, 3.0)
ax.set_xlabel('k  (homothety)'); ax.set_ylabel('d  (curvature dial)')
ax.set_title('Figure 3. Change the unit and a law slides along its curve d/k² = const.\nOn the k = 0 line only the sign of d is left (Lévy-Leblond 1976).',
             loc='left', fontsize=11.5, weight='bold')
save(fig, 'fig3-units.png')

# ---- fig 4: one operation, two places
fig, axs = plt.subplots(1, 2, figsize=(12, 5.2))
ax = axs[0]
xx = np.linspace(-0.6, 1.6, 300)
gen = {0: 1 + 0 * xx, 0.5: None, 1: 1 + xx, 2: (1 + xx) * (1 + np.log1p(xx))}
# s = 1/2 from the series (first six terms, k = 1)
gen[0.5] = 1 + 0.5 * xx + (0.5 * -0.5 / 4) * xx**2 - (0.5 * -0.5 / 12) * xx**3 + (0.5 * -0.5 * 4.5 / 144) * xx**4
cols = {0: GREY, 0.5: '#7aa6d8', 1: BLUE, 2: '#0b2f63'}
for sv_, g in gen.items():
    ax.plot(xx, g, color=cols[sv_], lw=2.3, label='s = %g' % sv_ + ('  (addition)' if sv_ == 0 else '  (∘₁)' if sv_ == 1 else ''))
ax.set_xlabel('x'); ax.set_ylabel('generator P_s(x) = (Eˢ)′(E⁻ˢ x)')
ax.set_title('(a) along one axis: rung s of the fused law ∘₁', loc='left')
ax.legend(fontsize=9.5)
ax = axs[1]
th = np.linspace(0, 2 * np.pi, 600)
for sv_, col in [(0, GREY), (0.5, '#7aa6d8'), (1, BLUE), (2, '#0b2f63'), (4, 'k')]:
    pp = 2**sv_
    c, s_ = np.cos(th), np.sin(th)
    n = (np.abs(c)**pp + np.abs(s_)**pp)**(1 / pp)
    ax.plot(c / n, s_ / n, color=col, lw=2, label='s = %g  →  p = %g' % (sv_, round(pp, 3)))
ax.set_aspect('equal'); ax.set_xlim(-1.5, 2.6)
ax.legend(fontsize=9.5, loc='center right')
ax.set_title('(b) across the axes: rung s of (a² + b²)^½ is ℓp, p = 2ˢ', loc='left')
fig.suptitle('Figure 4. Height and roundness are one operation (iterate the generator s times) in two places; a shape can set both',
             fontsize=12, weight='bold', x=0.02, ha='left')
fig.tight_layout()
save(fig, 'fig4-height-two-places.png')

# ---- fig 5: r = odds + zoom
fig, axs = plt.subplots(1, 2, figsize=(12, 5.4), gridspec_kw={'width_ratios': [1.1, 1]})
ax = axs[0]
V = np.array([[0, 0], [1, 0], [0.5, math.sqrt(3) / 2]])
ax.add_patch(Polygon(V, closed=True, fc=AMBER, alpha=0.08, ec=AMBER))
for i, lab in enumerate(['piece 1', 'piece 2', 'piece 3']):
    ax.text(*(V[i] + np.array([[-0.08, -0.06], [0.02, -0.06], [-0.04, 0.04]])[i]), lab, fontsize=9.5)


def moran(rs):
    lo, hi = 1e-9, 60
    for _ in range(200):
        mid = (lo + hi) / 2
        if sum(r**mid for r in rs) > 1:
            lo = mid
        else:
            hi = mid
    return (lo + hi) / 2


examples = [([0.8, 0.2, 0.2], 'r = 0.8, 0.2, 0.2', (8, -4)), ([0.2, 0.7, 0.2], 'r = 0.2, 0.7, 0.2', (-40, -30)),
            ([0.3, 0.3, 0.75], 'r = 0.3, 0.3, 0.75', (10, -2)), ([0.6, 0.5, 0.15], 'r = 0.6, 0.5, 0.15', (6, -26))]
for rs, lab, off in examples:
    D = moran(rs); w = np.array([r**D for r in rs])
    P = w @ V
    ax.plot(*P, 'o', color=BLUE, ms=8)
    ax.annotate('%s\nD = %.3f' % (lab, D), P, xytext=off, textcoords='offset points', fontsize=8.5)
C = np.ones(3) / 3 @ V
ax.plot(*C, 'o', color=AMBER, ms=10)
ax.annotate('all equal (½,½,½ Sierpinski D = 1.585;\n0.3,0.3,0.3 D = 0.912): the centre', C, xytext=(-60, 22), textcoords='offset points', fontsize=8.5, color=AMBER)
ax.set_aspect('equal'); ax.axis('off')
ax.set_title('(a) where the natural measure puts its weight w = rᴰ\n(equal ratios always sit at the centre, whatever their size)', loc='left', fontsize=10.5, pad=22)
ax = axs[1]
rr = np.linspace(0.05, 0.95, 200)
for Nn, col in [(2, GREY), (3, BLUE), (4, '#0b2f63')]:
    ax.plot(rr, np.log(Nn) / np.log(1 / rr), color=col, lw=2, label='N = %d equal pieces' % Nn)
ax.set_ylim(0, 6); ax.set_xlabel('common ratio r'); ax.set_ylabel('dimension D = log N / log(1/r)')
ax.legend(fontsize=9.5)
ax.set_title('(b) the one number of r that is its own: the zoom ↔ D', loc='left', fontsize=10.5)
fig.suptitle('Figure 5. n ratios = (D, a point of the simplex): the point is relative odds between pieces, the D is the zoom',
             fontsize=12, weight='bold', x=0.02, ha='left')
fig.tight_layout()
save(fig, 'fig5-ratio-is-odds-plus-zoom.png')

# ---- fig 6: odds absorb k
t = np.array([0.4, 1.0, 2.5])
subsets = [S for r_ in range(4) for S in itertools.combinations(range(3), r_)]
labels = ['∅' if not S else '{' + ','.join(str(i + 1) for i in S) + '}' for S in subsets]


def weights(kv, tv):
    return np.array([kv**(len(S) - 1) * np.prod(tv[list(S)]) if S else 1 / kv for S in subsets])


wA, wB = weights(2.0, t), weights(1.0, 2.0 * t)
fig, axs = plt.subplots(1, 2, figsize=(12, 4.4))
X = np.arange(len(subsets))
axs[0].bar(X - 0.2, wA, 0.4, color=BLUE, label='k = 2, odds t')
axs[0].bar(X + 0.2, wB, 0.4, color=AMBER, label='k = 1, odds 2t')
axs[0].set_xticks(X); axs[0].set_xticklabels(labels); axs[0].set_title('(a) corner weights k^(|S|−1)·∏t: they differ (by the factor k)', loc='left', fontsize=10.5)
axs[0].legend(fontsize=9.5)
axs[1].bar(X - 0.2, wA / wA.sum(), 0.4, color=BLUE, label='k = 2, odds t')
axs[1].bar(X + 0.2, wB / wB.sum(), 0.4, color=AMBER, label='k = 1, odds 2t')
axs[1].set_xticks(X); axs[1].set_xticklabels(labels); axs[1].set_title('(b) chance of each subset: identical, k is absorbed', loc='left', fontsize=10.5)
axs[1].legend(fontsize=9.5)
fig.suptitle('Figure 6. With odds on every element, the size of k is just a common odds scale (rank 5 → 4)', fontsize=12, weight='bold', x=0.02, ha='left')
fig.tight_layout()
save(fig, 'fig6-odds-absorb-k.png')

# ---- fig 7: q
fig, axs = plt.subplots(1, 2, figsize=(12, 4.8), gridspec_kw={'width_ratios': [1, 1.2]})
ax = axs[0]
dd = np.array([-1, -0.4, 0, 0.3, 1.0])
ax.plot(dd, dd, 'o', color=AMBER, ms=9, label='measured: turn / area of two small moves')
ax.plot([-1.2, 1.2], [-1.2, 1.2], color=GREY, lw=1, label='turn = d × area')
ax.set_xlabel('curvature dial d'); ax.set_ylabel('order turn per unit area')
ax.legend(fontsize=9.5, loc='upper left')
ax.set_title('(a) moves built from the law: the order cost is d', loc='left', fontsize=10.5)
ax = axs[1]; ax.axis('off'); ax.set_xlim(0, 10); ax.set_ylim(0, 6)
rows = [(AMBER, 'Boosts, rotations, translations\nof the (k, d) geometry', 'q = e^(i·d·area): a reading of d\n(Gauss–Bonnet)'),
        (AMBER, 'Counting (q-binomials, q-integers)', '[n]_q is the ∘ₖ power of 1 at k = q − 1'),
        (GREEN, 'Weyl / Moyal / quantum plane\n(moves on functions)', 'new constant put in by hand; continuous moves:\nħ is a unit; lattice moves: θ is a true invariant')]
for i, (col, a_, b_) in enumerate(rows):
    y = 5 - 1.8 * i
    ax.add_patch(FancyBboxPatch((0.1, y - 0.7), 9.7, 1.4, boxstyle='round,pad=0.04', fc=col, alpha=0.1, ec=col))
    ax.text(0.3, y, a_, va='center', fontsize=9.8, weight='bold')
    ax.text(5.2, y, b_, va='center', fontsize=9.4)
ax.set_title('(b) the three places q shows up', loc='left', fontsize=10.5)
fig.suptitle('Figure 7. The order dial is partly derived: only the third row is a number no other dial supplies', fontsize=12, weight='bold', x=0.02, ha='left')
fig.tight_layout()
save(fig, 'fig7-order-dial.png')
print('wrote', sorted(os.listdir('figures')))

# ---- fig 8: the mixing chart, primaries and their secondaries
# Every derived dial found is built from one or two primaries, so the picture is a table of pairs:
# cell (row, column) holds the dials mixed from those two primaries, coloured as their mix; the diagonal
# holds one-primary readings (tints). A wheel stacked different pairs at the same angle, and a graph
# tangled because k mixes with five other dials; a table shows every pair once, empty pairs included.
from derived_dials import PRIMARIES, PRIMARY_NAMES, DERIVED
import matplotlib.colors as mcolors
order = ['k', 'd', 's', 'p', 't', 'r', 'q']
pcol = {'k': '#d62728', 'd': '#ff7f0e', 'q': '#bcbd22', 'r': '#2ca02c', 't': '#17becf', 'p': '#1f5fae', 's': '#9467bd'}
cells = {}
for name, sym, par, *_ in DERIVED:
    a, b = (par[0], par[0]) if len(par) == 1 else sorted(par, key=order.index)
    cells.setdefault((a, b), []).append((name, sym))
n = len(order); W, H = 2.35, 1.75
fig, ax = plt.subplots(figsize=(16.5, 13.2))
ax.set_aspect('equal'); ax.axis('off')
ax.set_xlim(-2.6, n * W + 0.2); ax.set_ylim(-n * H - 1.6, 2.0)
for j, key in enumerate(order):   # column heads
    pale = key in ('r', 'q')
    ax.add_patch(FancyBboxPatch((j * W + 0.08, 0.15), W - 0.16, 1.2, boxstyle='round,pad=0.02,rounding_size=0.12',
                                fc=pcol[key], alpha=0.55 if pale else 0.92, ec='white'))
    ax.text(j * W + W / 2, 0.92, key, ha='center', va='center', fontsize=22, weight='bold', color='white')
    ax.text(j * W + W / 2, 0.42, PRIMARY_NAMES[key], ha='center', va='center', fontsize=9, color='white', weight='bold')
for i, key in enumerate(order):   # row heads
    pale = key in ('r', 'q')
    y = -(i + 1) * H
    ax.add_patch(FancyBboxPatch((-2.45, y + 0.08), 2.25, H - 0.16, boxstyle='round,pad=0.02,rounding_size=0.12',
                                fc=pcol[key], alpha=0.55 if pale else 0.92, ec='white'))
    ax.text(-1.32, y + H / 2 + 0.22, key, ha='center', va='center', fontsize=22, weight='bold', color='white')
    ax.text(-1.32, y + H / 2 - 0.42, PRIMARY_NAMES[key], ha='center', va='center', fontsize=9, color='white', weight='bold')
for i, a in enumerate(order):
    for j, b in enumerate(order):
        if j < i:
            continue
        x, y = j * W, -(i + 1) * H
        mix = np.mean([mcolors.to_rgb(pcol[a]), mcolors.to_rgb(pcol[b])], axis=0)
        got = cells.get((a, b), [])
        if got:
            ax.add_patch(FancyBboxPatch((x + 0.08, y + 0.08), W - 0.16, H - 0.16, boxstyle='round,pad=0.02,rounding_size=0.12',
                                        fc=mix, alpha=0.22, ec=mix, lw=2.4))
            ax.add_patch(FancyBboxPatch((x + 0.08, y + H - 0.44), W - 0.16, 0.36, boxstyle='round,pad=0.0,rounding_size=0.1',
                                        fc=mix, alpha=0.95, ec='none'))
            ax.text(x + W / 2, y + H - 0.26, (a if a == b else a + ' + ' + b) + ('  (tint)' if a == b else ''),
                    ha='center', va='center', fontsize=9.5, weight='bold', color='white')
            step = min(0.32, (H - 0.62) / max(len(got), 1))
            for m, (name, sym) in enumerate(got):
                ax.text(x + W / 2, y + H - 0.66 - m * step, name + '  ' + sym, ha='center', va='center',
                        fontsize=8.6 if len(got) < 4 else 8.0)
        else:
            ax.add_patch(FancyBboxPatch((x + 0.08, y + 0.08), W - 0.16, H - 0.16, boxstyle='round,pad=0.02,rounding_size=0.12',
                                        fc='#f3f3f3', ec='#dddddd', lw=1))
            ax.text(x + W / 2, y + H / 2, 'none found', ha='center', va='center', fontsize=8, color='#aaaaaa')
ax.text(-2.45, -n * H - 0.55, 'Diagonal: tints, one primary read another way.  Off the diagonal: secondaries, two primaries mixed, '
        'coloured as their mix.\nGrey: no derived dial found for that pair.  r and q are paler because they are partly derived '
        'themselves (README sections 6 and 7).', fontsize=10, va='top')
ax.set_title('Figure 8. The mixing chart: every derived dial sits in the cell of the one or two primaries it is built from',
             fontsize=14, weight='bold', loc='left')
save(fig, 'fig8-mixing-chart.png')
print('mixing chart done')
