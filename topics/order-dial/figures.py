"""Figures for the order dial write-up. Run: python3 figures.py (needs numpy, matplotlib)."""
import itertools
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle, FancyArrowPatch, Polygon
OUT = 'figures/'
plt.rcParams.update({'font.size': 11, 'axes.titlesize': 13, 'figure.dpi': 130})
BLUE, ORANGE, GREEN, RED, GREY, PURPLE = '#2f6db5', '#e08a1e', '#3a9a4a', '#c53b3b', '#888888', '#7a4fb0'

def arrow(ax, p, q, color='k', lw=1.5, style='-|>', ms=14):
    ax.add_patch(FancyArrowPatch(p, q, arrowstyle=style, mutation_scale=ms, color=color, lw=lw))

# ---------- fig1: the cube of laws ----------
fig = plt.figure(figsize=(9.5, 7.5)); ax = fig.add_subplot(111); ax.set_axis_off()
ox, oy = 0.08, 0.10; dx = (0.55, 0.0); dy = (0.0, 0.50); dz = (0.25, 0.20)
def P(g, o, r):  # g=grouping law, o=order law, r=repetition law (0/1)
    return (ox + g*dx[0] + o*dz[0] + r*dy[0], oy + g*dx[1] + o*dz[1] + r*dy[1])
edges = [((0,0,0),(1,0,0)),((0,1,0),(1,1,0)),((0,0,1),(1,0,1)),((0,1,1),(1,1,1)),
         ((0,0,0),(0,1,0)),((1,0,0),(1,1,0)),((0,0,1),(0,1,1)),((1,0,1),(1,1,1)),
         ((0,0,0),(0,0,1)),((1,0,0),(1,0,1)),((0,1,0),(0,1,1)),((1,1,0),(1,1,1))]
for u, v in edges:
    a_, b_ = P(*u), P(*v); ax.plot([a_[0], b_[0]], [a_[1], b_[1]], color='#bbbbbb', lw=1.5, transform=ax.transAxes)
labels = {
 (0,0,0): ('Trees', 'no laws: grouping, order and\nrepeats all count'),
 (1,0,0): ('Lists / arrays', '+ associative (grouping\nwashes out): free monoid'),
 (1,1,0): ('Bags / multisets', '+ commutative (order\nwashes out): counts only'),
 (1,1,1): ('Finite sets', '+ idempotent (x∪x = x):\nfree semilattice; CRDTs'),
 (1,0,1): ('Ordered sets', 'associative + idempotent, not\ncommutative: LinkedHashSet,\nfree (left regular) band'),
 (0,1,0): ('Mobiles', 'commutative but not\nassociative (unordered trees)'),
 (0,0,1): ('', ''), (0,1,1): ('', ''),
}
for key, (name, sub) in labels.items():
    x_, y_ = P(*key)
    col = BLUE if name else '#cccccc'
    ax.plot(x_, y_, 'o', ms=11, color=col, transform=ax.transAxes, zorder=3)
    if name:
        ax.text(x_ + 0.012, y_ + 0.018, name, fontsize=12, weight='bold', transform=ax.transAxes)
        ax.text(x_ + 0.012, y_ - 0.012, sub, fontsize=8.5, va='top', color='#333', transform=ax.transAxes)
mx, my = P(1,1,0)
ax.plot(mx, my, 'o', ms=22, mfc='none', mec=ORANGE, mew=3, transform=ax.transAxes, zorder=4)
ax.annotate('The circle operation sits here:\nassociative and commutative,\nbut 1 ✳ 1 = 3, not 1.\nIts only idempotents are\n0 and the wall −1/k',
            xy=(mx, my), xytext=(0.30, 0.86), textcoords='axes fraction', xycoords='axes fraction',
            fontsize=10, color=ORANGE, arrowprops=dict(arrowstyle='->', color=ORANGE, lw=1.8))
a0 = P(0,0,0)
arrow(ax, (a0[0]-0.02, a0[1]-0.11), (a0[0]+0.25, a0[1]-0.11)); ax.text(a0[0]+0.27, a0[1]-0.12, 'grouping stops mattering\n(associative)', fontsize=9, transform=ax.transAxes)
arrow(ax, (a0[0]-0.04, a0[1]-0.02), (a0[0]-0.04, a0[1]+0.22)); ax.text(a0[0]-0.075, a0[1]+0.24, 'repeats stop\ncounting\n(idempotent)', fontsize=9, transform=ax.transAxes)
ax.text(P(0,1,0)[0]-0.06, P(0,1,0)[1]+0.03, '', transform=ax.transAxes)
ax.text(0.20, 0.20, 'diagonal edges:\norder stops mattering\n(commutative)', fontsize=9, color=GREY, transform=ax.transAxes)
ax.set_title('The cube of laws: every collection is a choice of three yes/no rules\n(the Boom hierarchy is the path tree → list → bag → set; known)', fontsize=12)
fig.savefig(OUT + 'fig1-cube-of-laws.png', bbox_inches='tight'); plt.close(fig)

# ---------- fig2: everything is a map ----------
fig, ax = plt.subplots(figsize=(10.5, 5.6)); ax.set_axis_off(); ax.set_xlim(0, 10.5); ax.set_ylim(0, 6)
rows = [
 ('Set', 'key → yes/no', ['apple: 1', 'pear: 0', 'fig: 1'], 'join = OR', 'k = 1 corner\n(each at most once)', BLUE),
 ('Bag', 'key → count', ['apple: 2', 'pear: 0', 'fig: 1'], 'join = add counts', 'k = −1\n(any number)', GREEN),
 ('Random set', 'key → chance', ['apple: .3', 'pear: .5', 'fig: .9'], 'join = p + q − pq', 'this IS the circle operation\nat k = −1', ORANGE),
 ('Dictionary', 'key → value', ['apple: red', 'pear: green', 'fig: ?'], 'merge = last write wins', 'order matters only\non shared keys', PURPLE),
 ('List / array', 'position → value', ['0: fig', '1: apple', '2: apple'], 'join = append', 'order lives in\nthe keys 0, 1, 2', RED),
]
ax.text(0.2, 5.6, 'structure', weight='bold'); ax.text(1.9, 5.6, 'what each key holds', weight='bold')
ax.text(4.1, 5.6, 'example', weight='bold'); ax.text(6.9, 5.6, 'how two combine', weight='bold'); ax.text(8.9, 5.6, 'in the framework', weight='bold')
for i, (nm, holds, ex, join, fw, col) in enumerate(rows):
    y = 4.8 - i*1.05
    ax.add_patch(Rectangle((0.1, y-0.35), 10.3, 0.85, color=col, alpha=0.08))
    ax.text(0.2, y, nm, fontsize=12, weight='bold', color=col, va='center')
    ax.text(1.9, y, holds, fontsize=10, va='center')
    for j, e_ in enumerate(ex):
        ax.add_patch(Rectangle((4.0 + j*0.95, y-0.25), 0.9, 0.5, fc='white', ec=col))
        ax.text(4.45 + j*0.95, y, e_, fontsize=7.5, ha='center', va='center')
    ax.text(6.9, y, join, fontsize=10, va='center')
    ax.text(8.9, y, fw, fontsize=8.5, va='center', color=col)
ax.set_title('Every collection is a map from keys to "how present". The value rule picks the structure.\nSet theory builds lists the same way: a sequence is a function from {0, 1, ..., n−1}. (Known.)', fontsize=11.5)
fig.savefig(OUT + 'fig2-everything-is-a-map.png', bbox_inches='tight'); plt.close(fig)

# ---------- fig3: the order dial q ----------
def inv(w): return sum(1 for i in range(len(w)) for j in range(i+1, len(w)) if w[i] > w[j])
words = sorted(set(itertools.permutations('xxyy')), key=lambda w: (inv(w), w))
fig, (a1, a2) = plt.subplots(1, 2, figsize=(11.5, 4.8), gridspec_kw={'width_ratios': [1.05, 1]})
a1.set_axis_off(); a1.set_xlim(0, 6); a1.set_ylim(-0.3, 6.6)
a1.text(0.05, 6.3, 'The 6 orders of x, x, y, y', weight='bold')
a1.text(3.9, 6.3, 'out-of-order\npairs', weight='bold', fontsize=9, ha='center', va='center')
a1.text(5.3, 6.3, 'weight', weight='bold', fontsize=9, ha='center')
for i, w in enumerate(words):
    y = 5.4 - i*0.95
    for j, ch in enumerate(w):
        a1.add_patch(Rectangle((0.2 + j*0.62, y-0.27), 0.55, 0.55, color=BLUE if ch == 'x' else ORANGE, alpha=0.85))
        a1.text(0.475 + j*0.62, y, ch, ha='center', va='center', color='white', weight='bold')
    a1.text(3.9, y, str(inv(w)), ha='center', va='center', fontsize=12)
    a1.text(5.3, y, ['1', 'q', 'q²', 'q²', 'q³', 'q⁴'][i], ha='center', va='center', fontsize=12)
a1.text(0.05, -0.2, 'Each y that sits before an x is one swap away from sorted.', fontsize=9, color='#444')
qs = np.linspace(0, 1, 200); val = 1 + qs + 2*qs**2 + qs**3 + qs**4
a2.plot(qs, val, color=PURPLE, lw=3)
a2.plot([0], [1], 'o', color=GREEN, ms=10); a2.plot([1], [6], 'o', color=RED, ms=10)
a2.annotate('q = 0: only the sorted order counts.\nOrder is ignored: one bag.', (0, 1), (0.05, 3.4), fontsize=10, color=GREEN, arrowprops=dict(arrowstyle='->', color=GREEN))
a2.annotate('q = 1: every order counts.\nSix lists.', (1, 6), (0.45, 5.4), fontsize=10, color=RED, arrowprops=dict(arrowstyle='->', color=RED))
a2.set_xlabel('order dial q (weight per out-of-order pair)'); a2.set_ylabel('how many arrangements count')
a2.set_title('1 + q + 2q² + q³ + q⁴, the Gaussian binomial [4, 2]_q')
a2.grid(alpha=0.3)
fig.suptitle('A smooth order dial: weigh each arrangement by q^(out-of-order pairs). Bags at q = 0, lists at q = 1. (Known: Gauss, MacMahon.)', fontsize=11.5, y=1.04)
fig.savefig(OUT + 'fig3-order-dial-q.png', bbox_inches='tight'); plt.close(fig)

# ---------- fig4: order in the formula costs grouping ----------
rng = np.random.default_rng(1)
Mx = [rng.standard_normal((3, 3)) for _ in range(3)]
def mm(X, Y, l): return X + Y + l*X@Y + (1-l)*Y@X
lams = np.linspace(0, 1, 201)
X, Y, Z = Mx
assoc = [np.abs(mm(mm(X,Y,l),Z,l) - mm(X,mm(Y,Z,l),l)).max() for l in lams]
flex = [np.abs(mm(mm(X,Y,l),X,l) - mm(X,mm(Y,X,l),l)).max() for l in lams]
comm = [np.abs(mm(X,Y,l) - mm(Y,X,l)).max() for l in lams]
fig, ax = plt.subplots(figsize=(9.5, 5.2))
ax.plot(lams, np.array(assoc)/max(assoc), color=RED, lw=3, label='grouping error (a✳b)✳c − a✳(b✳c)')
ax.plot(lams, np.array(comm)/max(comm), color=BLUE, lw=3, label='order error a✳b − b✳a')
ax.plot(lams, np.array(flex)/max(assoc), color=GREEN, lw=3, ls='--', label='flexibility error (a✳b)✳a − a✳(b✳a): always 0')
for xv, txt in [(0, 'λ = 0\nb·a only\n(associative)'), (0.5, 'λ = ½: Jordan\n(commutative,\nnot associative)'), (1, 'λ = 1\na·b only\n(associative)')]:
    ax.axvline(xv, color=GREY, lw=1, ls=':'); ax.text(xv, 1.08, txt, ha='center', fontsize=9)
ax.set_xlabel('λ in  a ✳ b = a + b + λ·ab + (1 − λ)·ba   (3×3 matrices)'); ax.set_ylabel('size of error (scaled)')
ax.set_ylim(-0.05, 1.32); ax.legend(loc='center left', bbox_to_anchor=(0.0, 0.55), fontsize=9); ax.grid(alpha=0.3)
ax.set_title('Putting an order dial into the formula breaks grouping in between.\nThe grouping error is λ(1 − λ) × fixed matrices, so only the two ends are associative. (Known: quasi-associative / mutation algebras.)', fontsize=10.5)
fig.savefig(OUT + 'fig4-order-in-the-formula.png', bbox_inches='tight'); plt.close(fig)

# ---------- fig5: Cayley-Dickson ladder ----------
algs = ['Reals', 'Complex', 'Quaternions', 'Octonions', 'Sedenions']
laws = ['ordered (<)', 'commutative', 'associative', 'alternative', 'no zero\ndivisors']
has = np.array([[1,1,1,1,1],[0,1,1,1,1],[0,0,1,1,1],[0,0,0,1,1],[0,0,0,0,0]])
fig, ax = plt.subplots(figsize=(10.5, 4.6))
for i in range(5):
    for j in range(5):
        ax.add_patch(Rectangle((j, 4-i), 0.95, 0.9, color=GREEN if has[i, j] else '#e8e8e8'))
        ax.text(j+0.475, 4-i+0.45, 'yes' if has[i, j] else 'no', ha='center', va='center', color='white' if has[i, j] else '#999', weight='bold')
    ax.text(-0.1, 4-i+0.45, f'{algs[i]} (dim {2**i})', ha='right', va='center', fontsize=11)
for j in range(5): ax.text(j+0.475, 5.05, laws[j], ha='center', va='bottom', fontsize=10)
ax.set_xlim(-2.6, 7.6); ax.set_ylim(-0.1, 5.8); ax.set_axis_off()
ax.text(5.25, 3.0, 'The circle operation over each one,\na ✳ b = a + b + ab,\nhas exactly the same\nrow of yes and no.\n\nThe shift x ↦ 1 + x turns\n✳ into the product, so\nevery law moves across.\nSedenion zero divisors\nbecome extra walls.', fontsize=10, va='center', color=ORANGE)
ax.set_title('Each step up the Cayley–Dickson ladder drops one law. The circle operation copies its inputs\' laws exactly. (Checked numerically; the shift argument is elementary.)', fontsize=11)
fig.savefig(OUT + 'fig5-cayley-dickson-ladder.png', bbox_inches='tight'); plt.close(fig)

# ---------- fig6: compounding order ----------
fig, (a1, a2) = plt.subplots(1, 2, figsize=(11, 4.8))
Bv, d, r = 1.0, 0.5, 0.4
for ax, title, dep_first in [(a1, 'Grow, then deposit\n(end-of-period deposit)', False), (a2, 'Deposit, then grow\n(start-of-period deposit)', True)]:
    ax.add_patch(Rectangle((0, 0), 1, Bv, color=BLUE, alpha=0.8)); ax.text(0.5, Bv/2, 'balance B', ha='center', va='center', color='white', weight='bold')
    ax.add_patch(Rectangle((1, 0), r, Bv, color=BLUE, alpha=0.35)); ax.text(1+r/2, Bv/2, 'interest\nr·B', ha='center', va='center', fontsize=9)
    ax.add_patch(Rectangle((0, Bv), 1, d, color=GREEN, alpha=0.8)); ax.text(0.5, Bv+d/2, 'deposit d', ha='center', va='center', color='white', weight='bold')
    if dep_first:
        ax.add_patch(Rectangle((1, Bv), r, d, color=ORANGE, alpha=0.9)); ax.text(1+r/2, Bv+d/2, 'r·d', ha='center', va='center', weight='bold')
        ax.text(0.7, Bv+d+0.12, '(B + d)(1 + r) = B + rB + d + r·d', ha='center', fontsize=10)
    else:
        ax.add_patch(Rectangle((1, Bv), r, d, fill=False, ec=ORANGE, ls='--', lw=2)); ax.text(1+r/2, Bv+d/2, 'missing', ha='center', va='center', color=ORANGE, fontsize=8)
        ax.text(0.7, Bv+d+0.12, 'B(1 + r) + d = B + rB + d', ha='center', fontsize=10)
    ax.set_xlim(-0.1, 1.6); ax.set_ylim(-0.1, 1.9); ax.set_aspect('equal'); ax.set_axis_off(); ax.set_title(title)
    ax.annotate('', (1+r, -0.05), (0, -0.05), arrowprops=dict(arrowstyle='<->')); ax.text((1+r)/2, -0.09, 'stretch 1 + r', ha='center', va='top', fontsize=9)
fig.suptitle('Order changes a compounding result by exactly the corner r·d: the circle-operation box\'s corner piece is the price of swapping "add" and "multiply".', fontsize=11, y=1.04)
fig.savefig(OUT + 'fig6-compounding-order.png', bbox_inches='tight'); plt.close(fig)

# ---------- fig7: stretches about centres (the fan) ----------
fig, (a1, a2) = plt.subplots(1, 2, figsize=(12, 5.4))
k = 1.0; c = -1/k; xs = np.linspace(-3, 2, 50)
for s, col in zip([-0.8, -0.4, 0.3, 0.8, 1.5], ['#9ecae1', '#6baed6', '#4292c6', '#2171b5', '#08519c']):
    a1.plot(xs, (1+k*s)*xs + s, color=col, lw=2, label=f'a = {s}: x ↦ (1 + {s})x + {s}')
a1.plot(xs, xs, color=GREY, ls=':', lw=1)
a1.plot([c], [c], 'o', color=ORANGE, ms=12, zorder=5)
a1.annotate('every map in the family fixes\nthe wall −1/k (here −1).\nOne shared centre ⇒ any order\ngives the same result: a ✳ b = b ✳ a', (c, c), (-2.9, 1.3), fontsize=10, color=ORANGE, arrowprops=dict(arrowstyle='->', color=ORANGE))
a1.set_xlim(-3, 2); a1.set_ylim(-4, 3.5); a1.grid(alpha=0.3); a1.legend(fontsize=7.5, loc='lower right')
a1.set_title('The circle operation = composing stretches about one centre\n(T_comp in Lean: T(a) after T(b) = T(a ✳ b))'); a1.set_xlabel('x'); a1.set_ylabel('image of x')
# right: two centres
f = lambda x: 2*(x - 0) + 0      # stretch x2 about 0 (grow by r=1)
g = lambda x: 2*(x - 1) + 1      # stretch x2 about 1
a2.plot(xs, f(g(xs)), color=BLUE, lw=3, label='f after g')
a2.plot(xs, g(f(xs)), color=RED, lw=3, label='g after f')
a2.plot(xs, xs, color=GREY, ls=':', lw=1)
a2.plot([0], [0], 'o', color=BLUE, ms=9); a2.text(0.08, -0.5, 'centre of f', color=BLUE, fontsize=9)
a2.plot([1], [1], 'o', color=RED, ms=9); a2.text(1.05, 0.55, 'centre of g', color=RED, fontsize=9)
a2.annotate('', (-1.5, f(g(-1.5))), (-1.5, g(f(-1.5))), arrowprops=dict(arrowstyle='<->', color='k', lw=1.5))
a2.text(-1.45, (f(g(-1.5))+g(f(-1.5)))/2, '  gap = a·d − b·c = 1·(−1) − 1·0\n  the same everywhere: a pure shift', fontsize=9, va='center')
a2.set_xlim(-3, 2); a2.set_ylim(-9, 6); a2.grid(alpha=0.3); a2.legend(fontsize=9, loc='lower right'); a2.set_xlabel('x')
a2.set_title('Different centres ⇒ order matters.\nThe gap is the wedge of (stretch, shift) pairs.')
fig.suptitle('Composition is the home: the commuting slice through the wall is the circle operation; the full stretch-and-shift (affine) group is where order begins.\n(Known group; circle-operation reading ours.)', fontsize=11, y=1.06)
fig.savefig(OUT + 'fig7-centres-and-order.png', bbox_inches='tight'); plt.close(fig)

# ---------- fig8: route to quantum ----------
fig, (a1, a2) = plt.subplots(1, 2, figsize=(12, 5.2))
a_, b_ = 1.6, 1.0
a1.add_patch(Rectangle((0, 0), 1, 1, color=BLUE, alpha=0.75)); a1.text(0.5, 0.5, '1', ha='center', va='center', color='white', fontsize=14)
a1.add_patch(Rectangle((1, 0), a_-1+0.0, 1, color=BLUE, alpha=0.35)); a1.text(1+(a_-1)/2, 0.5, 'a', ha='center', va='center', fontsize=13)
a1.add_patch(Rectangle((0, 1), 1, b_*0.6, color=GREEN, alpha=0.4)); a1.text(0.5, 1+0.3, 'b', ha='center', va='center', fontsize=13)
a1.add_patch(Rectangle((1, 1), a_-1, b_*0.6, color=ORANGE, alpha=0.9)); a1.text(1+(a_-1)/2, 1.3, 'ab', ha='center', va='center', fontsize=13, weight='bold')
a1.set_xlim(-0.1, 2.0); a1.set_ylim(-0.1, 2.0); a1.set_aspect('equal'); a1.set_axis_off()
a1.set_title('Commutative the circle operation: the corner ab\nis kept as extra area (a ✳ b = a + b + ab)')
# phase plane loop
a2.add_patch(Rectangle((0, 0), a_, b_, color=ORANGE, alpha=0.35))
arrow(a2, (0, 0), (a_, 0), color=BLUE, lw=2.5); a2.text(a_/2, -0.15, 'shift position by a', ha='center', va='top', color=BLUE)
arrow(a2, (a_, 0), (a_, b_), color=GREEN, lw=2.5); a2.text(a_+0.07, b_/2, 'boost momentum\nby b', va='center', color=GREEN)
arrow(a2, (0, 0), (0, b_), color=GREEN, lw=2.5, style='-|>'); arrow(a2, (0, b_), (a_, b_), color=BLUE, lw=2.5)
a2.text(a_/2, b_/2, 'area ab', ha='center', va='center', fontsize=14, weight='bold')
a2.text(a_/2, b_+0.12, 'other order', ha='center', fontsize=9)
a2.set_xlim(-0.3, 2.6); a2.set_ylim(-0.5, 1.6); a2.set_aspect('equal'); a2.set_axis_off()
a2.set_xlabel('position'); a2.set_title('Quantum (Weyl relation): the two orders differ\nby the phase e^(i·ab/ħ), set by the same area')
fig.suptitle('Compose shifts (additions) with multiplications by e^(ibx) on functions, not points, and the corner becomes the cost of order.\nThat is the Heisenberg group of quantum mechanics, unique by Stone–von Neumann. The i and the size ħ are inputs, not forced. (Known: Weyl 1927.)', fontsize=10.5, y=1.10)
fig.savefig(OUT + 'fig8-route-to-quantum.png', bbox_inches='tight'); plt.close(fig)

# ---------- fig9: one wedge, three places ----------
fig, axs = plt.subplots(1, 3, figsize=(12.5, 4.3))
specs = [('Clifford (vectors)', 'x ✳ y − y ✳ x = 2 x∧y', 'x', 'y', BLUE),
         ('Affine maps (stretch, shift)', 'f∘g − g∘f shift = a·d − b·c', '(a, c)', '(b, d)', PURPLE),
         ('Heisenberg (shift, boost)', 'phase = (a₁b₂ − a₂b₁)/2ħ', '(a₁, b₁)', '(a₂, b₂)', ORANGE)]
for ax, (t, form, l1, l2, col) in zip(axs, specs):
    u = np.array([1.6, 0.4]); v = np.array([0.5, 1.3])
    ax.add_patch(Polygon([[0, 0], u, u+v, v], color=col, alpha=0.3))
    arrow(ax, (0, 0), tuple(u), color='k'); arrow(ax, (0, 0), tuple(v), color='k')
    ax.text(*(u + [0.05, -0.15]), l1, fontsize=10); ax.text(*(v + [-0.35, 0.08]), l2, fontsize=10)
    ax.text(*((u+v)/2), 'oriented\narea', ha='center', va='center', fontsize=10, weight='bold')
    ax.set_xlim(-0.4, 2.5); ax.set_ylim(-0.4, 2.0); ax.set_aspect('equal'); ax.set_axis_off()
    ax.set_title(t + '\n' + form, fontsize=10.5)
fig.suptitle('Wherever order matters, the cost of swapping is the same thing: an oriented area (a wedge). Zero exactly when the two arrows line up.', fontsize=11, y=1.04)
fig.savefig(OUT + 'fig9-one-wedge-three-places.png', bbox_inches='tight'); plt.close(fig)

# ---------- fig10: order x repetition grid ----------
fig, ax = plt.subplots(figsize=(11, 5.4)); ax.set_axis_off(); ax.set_xlim(-2.3, 8.2); ax.set_ylim(-0.6, 4.4)
cols = ['k = 1\neach at most once', 'k = 1/m\nat most m', 'k = 0\nclassical', 'k = −1\nany number']
rowsq = ['q = 0\norder ignored', 'q = 1\nevery order distinct', 'q = −1\nswap costs a sign']
cells = [['set\n(HashSet, TreeSet)', 'capped bag', 'Poisson count', 'bag / multiset\n(C++ multiset)'],
         ['arrangement\n(list, no repeats)', 'list, ≤ m repeats', 'list (Boltzmann\ncount, ÷ r!)', 'list / array\n(words)'],
         ['fermions,\nGrassmann numbers', 'forbidden\nfor m ≥ 2', 'forbidden', 'forbidden: a swap sign\nforces x² = 0']]
okc = [[1, 1, 1, 1], [1, 1, 1, 1], [1, 0, 0, 0]]
for j, cn in enumerate(cols): ax.text(j*2 + 0.9, 4.1, cn, ha='center', fontsize=10, weight='bold')
for i, rn in enumerate(rowsq):
    y = 2.8 - i*1.3
    ax.text(-0.25, y + 0.5, rn, ha='right', va='center', fontsize=10, weight='bold')
    for j in range(4):
        ok = okc[i][j]
        ax.add_patch(Rectangle((j*2, y), 1.85, 1.1, color=(GREEN if ok else RED), alpha=0.18 if ok else 0.12))
        ax.text(j*2 + 0.925, y + 0.55, cells[i][j], ha='center', va='center', fontsize=8.8, color='k' if ok else RED)
ax.set_title('Order dial q (rows) against the repetition dial k (columns). Anticommuting kills repetition:\nx·x = −x·x forces x² = 0, the Pauli rule. So the two dials are linked, not independent. (Known physics.)', fontsize=11)
fig.savefig(OUT + 'fig10-order-vs-repetition.png', bbox_inches='tight'); plt.close(fig)
print('figures written')
