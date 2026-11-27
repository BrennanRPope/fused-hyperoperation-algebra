"""Figures for topics/chaos-composition. Run: python3 figures.py (numpy, matplotlib, mpmath)."""
import math
from fractions import Fraction
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch
import mpmath as mp

OUT = 'figures/'
plt.rcParams.update({'font.size': 11, 'axes.titlesize': 13, 'figure.dpi': 130})
BLUE, ORANGE, GREEN, RED, GREY, PURPLE = '#2f6db5', '#e08a1e', '#3a9a4a', '#c53b3b', '#888888', '#7a4fb0'
PHI = (1 + 5 ** .5) / 2


def arrow(ax, p, q, color='k', lw=1.5, style='-|>', ms=14, **kw):
    ax.add_patch(FancyArrowPatch(p, q, arrowstyle=style, mutation_scale=ms, color=color, lw=lw, **kw))


def box(ax, xy, w, h, text, fc, ec='k', fs=10, **kw):
    ax.add_patch(FancyBboxPatch(xy, w, h, boxstyle='round,pad=0.01', fc=fc, ec=ec, lw=1.0, transform=ax.transAxes))
    ax.text(xy[0] + w / 2, xy[1] + h / 2, text, ha='center', va='center', fontsize=fs, transform=ax.transAxes, **kw)


# ---------- fig1: the dictionary ----------
rows = [
    ('a dynamical system (one map f, iterated)', 'a functor from (ℕ, +) to sets: n ↦ fⁿ', 'known', 'Lawvere–Schanuel'),
    ('orbit of a point', 'the image of that functor at the point', 'known', ''),
    ('conjugacy / factor map', 'natural isomorphism / natural transformation', 'known', 'iter_semiconj'),
    ('every dial square x ∘ₖ x', 'one system: x ∘₁ x in the units x ↦ kx', 'ours', 'dial_iter_units'),
    ('logistic map at r = 4', 'x ∘₁ x on the circle, folded by the circle inverse', 'known + reading', 'logistic_iter_fold'),
    ('Lyapunov exponent', 'growth rate of the composed slope (chain rule)', 'known', 'iter_slope'),
    ('f ∘ g against g ∘ f', 'the same system seen through g (shift equivalence)', 'known', 'order_invisible'),
    ('the order cost, inside iteration', 'only moves centres; the commutator is a shift by the wedge', 'known + reading', 'order_is_a_shift'),
    ('topological entropy of branches', 'growth of the branch monoid: lists log 2, bags 0', 'known + reading', 'runB_free, runT_count'),
    ('some branches share a centre', 'trace monoid: entropy from the centre pattern', 'not found in print', "checks F'"),
    ('symbolic dynamics (itineraries)', 'streams; doubling is the list tail', 'known', 'doubling_shift'),
    ('finite orbits (Pollard rho)', 'a list that loops back on itself: tail + cycle', 'known', 'checks H'),
    ('period doubling', 'compose with yourself, then change units', 'known', 'checks D'),
]
fig = plt.figure(figsize=(13, 8.2)); ax = fig.add_subplot(111); ax.set_axis_off()
ax.text(0.0, 0.985, 'Chaos words, read as composition words', fontsize=15, weight='bold', transform=ax.transAxes)
ax.text(0.0, 0.952, 'orange = known; blue = known pieces, our reading; green = not found in print. Right column: where it is checked or proved.',
        fontsize=9.5, color='#444', transform=ax.transAxes)
colmap = {'known': ORANGE, 'known + reading': BLUE, 'ours': BLUE, 'not found in print': GREEN}
H = 0.066
for i, (a, b, tag, where) in enumerate(rows):
    y0 = 0.865 - i * H
    c = colmap[tag]
    box(ax, (0.0, y0), 0.33, H * 0.82, a, '#f4f4f4', fs=10)
    arrow(ax, (0.335, y0 + H * 0.41), (0.375, y0 + H * 0.41), transform=ax.transAxes, color='#555')
    box(ax, (0.38, y0), 0.45, H * 0.82, b, c + '22', ec=c, fs=10)
    ax.text(0.84, y0 + H * 0.41, where, fontsize=8.5, va='center', color='#333', transform=ax.transAxes, family='monospace')
fig.savefig(OUT + 'fig1-dictionary.png', bbox_inches='tight'); plt.close(fig)

# ---------- fig2: one arrow, composed with itself ----------
fig, axs = plt.subplots(1, 2, figsize=(12.5, 5.6))
ax = axs[0]
t = np.linspace(0, 2 * np.pi, 400)
ax.plot(-1 + np.cos(t), np.sin(t), color='#bbb')
ax.plot([-1], [0], 'o', color=RED); ax.text(-1.05, -0.13, 'wall −1\n(superattracting)', ha='center', va='top', fontsize=9, color=RED)
ax.plot([0], [0], 'o', color='k'); ax.text(0.06, 0.05, '0 (repels, slope 2)', fontsize=9)
ang = 0.7
pts = []
for n in range(7):
    pts.append(ang); ang = (2 * ang) % (2 * np.pi)
for n, a in enumerate(pts):
    z = complex(-1 + np.cos(a), np.sin(a))
    ax.plot(z.real, z.imag, 'o', color=BLUE, ms=7)
    ax.text(z.real * 1.0 + 0.06 * np.cos(a), z.imag + 0.08 * np.sin(a) + 0.02, f'x{n}', fontsize=9, color=BLUE)
for a, b in zip(pts[:-1], pts[1:]):
    za = (-1 + np.cos(a), np.sin(a)); zb = (-1 + np.cos(b), np.sin(b))
    arrow(ax, za, zb, color=BLUE, lw=0.9, ms=9, connectionstyle='arc3,rad=0.25', alpha=0.7)
ax.set_aspect('equal'); ax.set_xlim(-2.4, 0.9); ax.set_ylim(-1.35, 1.35)
ax.set_title('x ↦ x ∘₁ x on the circle |1 + x| = 1:\nthe angle doubles each step (Lyapunov log 2)')
ax.set_xlabel('Re x'); ax.set_ylabel('Im x')
ax = axs[1]
xs = np.linspace(-3, 1.5, 400)
for kk, col in [(1, BLUE), (0.5, GREEN), (2, PURPLE), (-2, ORANGE)]:
    ax.plot(xs, xs + xs + kk * xs * xs, color=col, label=f'x ∘ₖ x, k = {kk}')
ax.plot(xs, xs, color='#999', ls='--', lw=1)
ax.set_ylim(-3, 3); ax.axhline(0, color='#ddd', lw=0.8); ax.axvline(0, color='#ddd', lw=0.8)
ax.legend(fontsize=9, loc='lower right')
ax.set_title('Every dial square is one map in other units:\nk·(x ∘ₖ x) = (kx) ∘₁ (kx), orbit for orbit')
ax.text(-2.95, -2.75, 'k = −2 is the logistic map at r = 2\n(dashed: y = x; fixed points sit on it)', fontsize=9, color=ORANGE)
ax.set_xlabel('x')
fig.tight_layout(); fig.savefig(OUT + 'fig2-one-arrow.png', bbox_inches='tight'); plt.close(fig)

# ---------- fig3: folding the circle ----------
fig = plt.figure(figsize=(12.5, 6))
ax = fig.add_axes([0.03, 0.08, 0.40, 0.84]); ax.set_axis_off()
pos = {'TL': (0.05, 0.72), 'TR': (0.62, 0.72), 'BL': (0.05, 0.18), 'BR': (0.62, 0.18)}
txt = {'TL': 'circle |1+x| = 1', 'TR': 'circle |1+x| = 1', 'BL': 'interval [0, 1]', 'BR': 'interval [0, 1]'}
for k_, (x_, y_) in pos.items():
    box(ax, (x_, y_), 0.33, 0.12, txt[k_], '#f4f4f4', fs=11)
arrow(ax, (0.39, 0.78), (0.61, 0.78), transform=ax.transAxes, color=BLUE)
ax.text(0.50, 0.82, 'x ↦ x ∘₁ x', ha='center', color=BLUE, transform=ax.transAxes)
arrow(ax, (0.39, 0.24), (0.61, 0.24), transform=ax.transAxes, color=ORANGE)
ax.text(0.50, 0.28, 'logistic, r = 4', ha='center', color=ORANGE, transform=ax.transAxes)
for xx_ in (0.21, 0.78):
    arrow(ax, (xx_, 0.71), (xx_, 0.31), transform=ax.transAxes, color='k')
ax.text(0.23, 0.50, 'fold\nw ↦ w·ι(w)/4', transform=ax.transAxes, fontsize=10)
ax.text(0.80, 0.50, 'fold', transform=ax.transAxes, fontsize=10)
ax.text(0.0, 0.02, 'The square commutes at every step (logistic_iter_fold):\nn logistic steps = n squarings, folded. ι = circle inverse; '
        'the fold\nforgets ι, which is why two circle points land on one number.', fontsize=9.5, transform=ax.transAxes)
ax.set_title('Chaos at r = 4 is a quotient of a group law', fontsize=13)
ax2 = fig.add_axes([0.52, 0.12, 0.45, 0.78])
t = np.linspace(0, 2 * np.pi, 400)
ax2.plot(t, np.sin(t / 2) ** 2, color='k', lw=1.2, label='fold of 1 + x = e^{iθ}: sin²(θ/2)')
for th0, col in [(1.1, BLUE), (2 * np.pi - 1.1, PURPLE)]:
    ax2.plot([th0], [np.sin(th0 / 2) ** 2], 'o', color=col, ms=8)
ax2.annotate('θ and 2π − θ (x and ι(x))\nfold to the same number', (2 * np.pi - 1.1, np.sin((2 * np.pi - 1.1) / 2) ** 2),
             (3.9, 0.12), fontsize=9, arrowprops=dict(arrowstyle='->', color='#555'))
ax2.set_xlabel('angle θ on the circle'); ax2.set_ylabel('folded value')
ax2.set_xticks([0, np.pi, 2 * np.pi]); ax2.set_xticklabels(['0', 'π', '2π'])
ax2.legend(fontsize=9, loc='upper left'); ax2.set_title('The fold, angle by angle')
fig.savefig(OUT + 'fig3-fold.png', bbox_inches='tight'); plt.close(fig)

# ---------- fig4: bifurcation and Lyapunov ----------
rs = np.linspace(1.8, 4.0, 2200)
X = np.full_like(rs, 0.3)
for _ in range(600):
    X = rs * X * (1 - X)
fig, (a1, a2) = plt.subplots(2, 1, figsize=(12, 8), sharex=True, gridspec_kw={'height_ratios': [3, 2]})
lyap = np.zeros_like(rs)
for i in range(300):
    X = rs * X * (1 - X)
    a1.plot(rs, X, ',', color='k', alpha=0.25)
    lyap += np.log(np.abs(rs * (1 - 2 * X)) + 1e-300)
lyap /= 300
supers = [2, 1 + 5 ** .5, 3.4985616993, 3.5546408628, 3.5666673799]
for j, r in enumerate(supers):
    a1.axvline(r, color=BLUE if j else RED, lw=0.8, ls=':')
a1.annotate('r = 2: every circle-dial square\n(x ∘ₖ x is logistic at r = 2 in other units)', (2.0, 0.5), (2.05, 0.83),
            fontsize=9.5, color=RED, arrowprops=dict(arrowstyle='->', color=RED))
a1.annotate('r = 4: x ∘₁ x on the circle, folded', (4.0, 0.98), (3.25, 0.995), fontsize=9.5, color=ORANGE,
            arrowprops=dict(arrowstyle='->', color=ORANGE))
a1.text(3.0, 0.06, 'blue dotted: superstable r₀, r₁, r₂, … ;\ngaps shrink by δ = 4.6692 (check D3)', fontsize=9, color=BLUE)
a1.set_ylabel('long-run x'); a1.set_title('Iterating rx(1 − x): where the dial sits, and where chaos starts')
a2.plot(rs, np.clip(lyap, -3, 1), color='k', lw=0.6)
a2.axhline(0, color=GREY, lw=0.8); a2.axhline(np.log(2), color=ORANGE, lw=0.8, ls='--')
a2.text(1.82, np.log(2) + 0.05, 'log 2', color=ORANGE, fontsize=9)
a2.text(1.85, -2.6, 'below 0: orbits settle; above 0: nearby orbits spread (chaos)', fontsize=9)
a2.set_ylabel('Lyapunov exponent\n(growth of composed slope)'); a2.set_xlabel('r'); a2.set_ylim(-3, 1)
fig.tight_layout(); fig.savefig(OUT + 'fig4-bifurcation-lyapunov.png', bbox_inches='tight'); plt.close(fig)

# ---------- fig5: order is invisible to iteration ----------
fig, axs = plt.subplots(1, 2, figsize=(12.5, 5.4))
ax = axs[0]; ax.set_axis_off()
P = {'a0': (0.08, 0.70), 'a1': (0.40, 0.70), 'a2': (0.72, 0.70), 'b0': (0.08, 0.25), 'b1': (0.40, 0.25), 'b2': (0.72, 0.25)}
for k_, (x_, y_) in P.items():
    ax.plot(x_ + 0.05, y_ + 0.03, 'o', color=BLUE if k_[0] == 'a' else PURPLE, ms=11, transform=ax.transAxes)
for i in range(2):
    arrow(ax, (P[f'a{i}'][0] + 0.09, 0.73), (P[f'a{i + 1}'][0] + 0.01, 0.73), transform=ax.transAxes, color=BLUE)
    ax.text((P[f'a{i}'][0] + P[f'a{i + 1}'][0]) / 2 + 0.05, 0.77, 'f ∘ g', ha='center', color=BLUE, transform=ax.transAxes)
    arrow(ax, (P[f'b{i}'][0] + 0.09, 0.28), (P[f'b{i + 1}'][0] + 0.01, 0.28), transform=ax.transAxes, color=PURPLE)
    ax.text((P[f'b{i}'][0] + P[f'b{i + 1}'][0]) / 2 + 0.05, 0.20, 'g ∘ f', ha='center', color=PURPLE, transform=ax.transAxes)
for i in range(3):
    arrow(ax, (P[f'a{i}'][0] + 0.05, 0.69), (P[f'b{i}'][0] + 0.05, 0.33), transform=ax.transAxes, color='k', lw=1)
    ax.text(P[f'a{i}'][0] + 0.065, 0.50, 'g', transform=ax.transAxes)
ax.text(0.0, 0.95, 'g carries every orbit of f ∘ g onto an orbit of g ∘ f', fontsize=11, weight='bold', transform=ax.transAxes)
ax.text(0.0, 0.03, 'Same periods, same multipliers, same Lyapunov exponent\n(order_invisible, order_periodic; known: shift equivalence)',
        fontsize=9.5, transform=ax.transAxes)
ax = axs[1]
sa, ta, sb, tb = 0.5, 1.0, 0.7, -0.6
xs = np.linspace(-3, 4, 50)
AB = lambda z: sa * (sb * z + tb) + ta
BA = lambda z: sb * (sa * z + ta) + tb
ax.plot(xs, AB(xs), color=BLUE, label='A ∘ B'); ax.plot(xs, BA(xs), color=PURPLE, label='B ∘ A')
ax.plot(xs, xs, color='#999', ls='--', lw=1, label='y = x')
s = sa * sb; cAB = (ta + sa * tb) / (1 - s); cBA = (tb + sb * ta) / (1 - s)
ax.plot(cAB, cAB, 'o', color=BLUE); ax.text(cAB - 0.2, cAB + 0.35, 'centre of A∘B', color=BLUE, fontsize=9, ha='right')
ax.plot(cBA, cBA, 'o', color=PURPLE); ax.text(cBA + 0.25, cBA - 0.35, 'centre of B∘A', color=PURPLE, fontsize=9)
arrow(ax, (cAB, cAB), (cBA, cBA), color='k', lw=1.4, connectionstyle='arc3,rad=0.5')
ax.text(1.4, -1.9, 'B carries one centre to the other;\nthe gap is wedge / (1 − s_A s_B)', fontsize=9)
ax.set_title('Two affine maps: same slope both ways,\nthe order only moves the centre')
ax.legend(fontsize=9, loc='upper left'); ax.set_xlabel('x')
fig.tight_layout(); fig.savefig(OUT + 'fig5-order-invisible.png', bbox_inches='tight'); plt.close(fig)


# ---------- fig6: order makes chaos ----------
def tree(maps, n, x0):
    levels = [[x0]]
    for _ in range(n):
        levels.append([s_ * v + t_ for v in levels[-1] for (s_, t_) in maps])
    return levels


fig, axs = plt.subplots(1, 3, figsize=(14, 5.4), sharey=True)
cases = [
    ('Two shifts x, x + 1 (k = 0 corner)\nbag: only the count of 1s matters', [(1, 0), (1, 1)], 0.0, ORANGE),
    ('Two stretches about one centre\n(slopes 1/2 and 1/3): bag again', [(0.5, 0.25), (1 / 3, 1 / 3)], 1.0, ORANGE),
    ('Two stretches, different centres 0 and 1\n(slope 1/3): list, every word its own point', [(1 / 3, 0), (1 / 3, 2 / 3)], 0.5, BLUE),
]
for ax, (title, maps, x0, col) in zip(axs, cases):
    lv = tree(maps, 6, x0)
    for d in range(len(lv) - 1):
        for i, v in enumerate(lv[d]):
            for j in range(2):
                w = lv[d + 1][2 * i + j]
                ax.plot([v, w], [d, d + 1], color=col, lw=0.6, alpha=0.6)
    for d, L in enumerate(lv):
        u = np.unique(np.round(L, 9))
        ax.plot(u, [d] * len(u), 'o', color='k', ms=3)
        ax.annotate(f'{len(u)}', (1.02, d), xycoords=('axes fraction', 'data'), fontsize=9, va='center')
    ax.set_title(title, fontsize=10.5)
    ax.set_xlabel('where the composite sends the start point')
axs[0].set_ylabel('word length n'); axs[0].invert_yaxis()
fig.text(0.5, -0.04, 'Numbers on the right: distinct composites at each length. Shared centre (wedge 0): n + 1, entropy 0 '
         '(runT_count, stretchRun_perm). Different centres: 2ⁿ, entropy log 2 (runB_free).', ha='center', fontsize=10)
fig.suptitle('Order is what makes chaos: the inverse branches must not commute', y=1.02, fontsize=13)
fig.tight_layout(); fig.savefig(OUT + 'fig6-order-makes-chaos.png', bbox_inches='tight'); plt.close(fig)

# ---------- fig7: entropy from the centre pattern ----------
fig, axs = plt.subplots(1, 2, figsize=(13, 5.4))
ax = axs[0]
pats = [('2 maps,\none centre', (2,), 1), ('2 maps,\n2 centres', (1, 1), 2), ('3 maps,\none centre', (3,), 1),
        ('3 maps:\n2 share', (2, 1), PHI ** 2), ('3 maps,\n3 centres', (1, 1, 1), 3),
        ('4 maps:\n2 + 2', (2, 2), 2 + 2 ** .5), ('4 maps:\n2 + 1 + 1', (2, 1, 1), 2 + 3 ** .5)]
xs_ = np.arange(len(pats))
hs = [math.log(p[2]) for p in pats]
cols = [ORANGE if h_ == 0 else (GREEN if len(p[1]) > 1 and max(p[1]) > 1 else BLUE) for p, h_ in zip(pats, hs)]
ax.bar(xs_, hs, color=cols)
for i, (nm, pat, rate) in enumerate(pats):
    lab = {1: '1', 2: '2', 3: '3'}.get(rate, None) or {PHI ** 2: 'φ²', 2 + 2 ** .5: '2+√2', 2 + 3 ** .5: '2+√3'}[rate]
    ax.text(i, hs[i] + 0.03, f'growth {lab}', ha='center', fontsize=8.5)
ax.set_xticks(xs_); ax.set_xticklabels([p[0] for p in pats], fontsize=8.5)
ax.set_ylabel('entropy = log(growth of distinct composites)')
ax.set_title('Which branches share a centre decides the entropy\n(trace monoid; growth = 1/root of the Cartier–Foata polynomial)')
ax.text(0.02, 0.97, 'green: partly commuting (not found in print as an affine statement)\nblue: free (list), orange: bag',
        transform=ax.transAxes, fontsize=8.5, va='top')
ax.set_ylim(0, 1.55)
ax = axs[1]


def count_exact(s_frac, n):
    vals = {Fraction(0)}; out = []
    for i in range(n):
        vals = vals | {v + s_frac ** i for v in vals}; out.append(len(vals))
    return out


def golden_counts(n):
    A, B = [1], [0]
    for i in range(n):
        A.append(B[-1]); B.append(A[-2] - B[-1])
    vals = {(0, 0)}; out = []
    for i in range(n):
        vals = vals | {(a + A[i], b + B[i]) for a, b in vals}; out.append(len(vals))
    return out


nmax = 16
ns = np.arange(1, nmax + 1)
ax.semilogy(ns, count_exact(Fraction(1, 2), nmax), 'o-', color=BLUE, ms=4, label='s = 1/2: free, 2ⁿ')
ax.semilogy(ns, count_exact(Fraction(3, 5), nmax), 's--', color=BLUE, ms=4, mfc='none', label='s = 3/5: free, 2ⁿ')
ax.semilogy(ns, golden_counts(nmax), 'o-', color=GREEN, ms=4, label='s = 1/φ (s + s² = 1): about φⁿ')
ax.semilogy(ns, ns + 1, 'o-', color=ORANGE, ms=4, label='s = 1 (shifts): n + 1, a bag')
ax.set_xlabel('word length n'); ax.set_ylabel('distinct composites (exact count)')
ax.set_title('Two branches x ↦ sx, x ↦ sx + 1 (centres 0 and 1/(1 − s)):\na relation like s + s² = 1 cuts the growth, but it stays exponential')
ax.legend(fontsize=9, loc='upper left')
fig.tight_layout(); fig.savefig(OUT + 'fig7-centre-pattern-entropy.png', bbox_inches='tight'); plt.close(fig)

# ---------- fig8: streams ----------
fig, axs = plt.subplots(1, 2, figsize=(13, 5.2))
ax = axs[0]; ax.set_axis_off()
bits = [1, 0, 1, 1, 0, 0, 1, 0, 1, 1, 1, 0, 1, 0, 0]
for row in range(5):
    y_ = 0.86 - row * 0.15
    ax.text(0.0, y_, f'Dⁿ x, n = {row}:', fontsize=10, transform=ax.transAxes, va='center')
    ax.text(0.22, y_, '0.', fontsize=12, family='monospace', transform=ax.transAxes, va='center')
    for j, b in enumerate(bits[row:row + 9]):
        ax.text(0.27 + j * 0.055, y_, str(b), fontsize=13, family='monospace', transform=ax.transAxes, va='center',
                color=RED if j == 0 else 'k', weight='bold' if j == 0 else 'normal')
    if row < 4:
        arrow(ax, (0.27, y_ - 0.03), (0.27, y_ - 0.12), transform=ax.transAxes, color=GREY, lw=1)
ax.text(0.0, -0.02, 'Doubling mod 1 drops the leading binary digit (red): it is the list tail\non an infinite list (a stream). '
        'Two points that share n leading digits\nstay together for n steps only (doubling_shift; checks I1, I2).',
        fontsize=9.5, transform=ax.transAxes)
ax.set_title('Chaos is iterating tail on a stream', fontsize=13)
ax = axs[1]
mp.mp.dps = 60
bb = [int(c) for c in '1011001011100101001110100101011001']
th = mp.mpf(int(''.join(map(str, bb)), 2)) / mp.mpf(2) ** len(bb)
N = 26
xs = [float(mp.sin(mp.pi * th * 2 ** n) ** 2) for n in range(N)]
ax.plot(range(N), xs, 'o-', color=ORANGE, ms=4, lw=0.8, label='logistic r = 4 orbit')
ax.axhline(0.5, color=GREY, lw=0.8, ls='--')
for n in range(N):
    g = bb[n] ^ bb[n + 1]
    ax.text(n, 1.08, str(bb[n]), ha='center', fontsize=8, family='monospace', color=BLUE)
    ax.text(n, 1.17, str(g), ha='center', fontsize=8, family='monospace', color=ORANGE)
ax.text(-1.6, 1.08, 'angle digits', ha='right', fontsize=8.5, color=BLUE)
ax.text(-1.6, 1.17, 'above ½?', ha='right', fontsize=8.5, color=ORANGE)
ax.set_ylim(-0.05, 1.25); ax.set_xlim(-6.5, N); ax.set_xticks(range(0, N, 5))
ax.set_xlabel('step n'); ax.set_ylabel('x')
ax.set_title('After the fold, "above ½ or not" is the Gray code\nb(n) xor b(n+1) of the angle digits (check I3)')
fig.tight_layout(); fig.savefig(OUT + 'fig8-streams.png', bbox_inches='tight'); plt.close(fig)

# ---------- fig9: Pollard rho ----------
fig, axs = plt.subplots(1, 2, figsize=(12.5, 5.6))
for ax, p, x0 in [(axs[0], 83, 2), (axs[1], 1009, 3)]:
    seen = {}; orb = []; v = x0
    while v not in seen:
        seen[v] = len(orb); orb.append(v); v = (v * v + 1) % p
    mu = seen[v]; lam = len(orb) - mu
    coords = []
    for i in range(len(orb)):
        if i < mu:
            coords.append((0, -(mu - i)))
        else:
            a = 2 * np.pi * (i - mu) / lam
            R_ = max(lam / (2 * np.pi), 1.0)
            coords.append((R_ * np.sin(a), R_ * (1 - np.cos(a))))
    cx, cy = zip(*coords)
    for i in range(len(orb)):
        j = i + 1 if i + 1 < len(orb) else mu
        arrow(ax, coords[i], coords[j], color=BLUE if i >= mu else PURPLE, lw=0.8, ms=7)
    ax.plot(cx, cy, 'o', color='k', ms=3)
    if len(orb) < 30:
        for i, (a_, b_) in enumerate(coords):
            ax.text(a_ + 0.12, b_ + 0.1, str(orb[i]), fontsize=8)
    ax.set_aspect('equal'); ax.set_axis_off()
    ax.set_title(f'x ↦ x² + 1 mod {p} from {x0}: tail {mu}, loop {lam}', fontsize=11)
fig.suptitle('On a finite set every orbit is a ρ: a linked list whose last node points back into itself.\n'
             'Pollard rho finds factors by noticing the loop (Floyd: tortoise and hare). In lift coordinates x² + 1 is y ∘₁ y + 1.',
             fontsize=11, y=1.03)
fig.savefig(OUT + 'fig9-rho.png', bbox_inches='tight'); plt.close(fig)

# ---------- fig10: renormalization ----------
supers = [2.0, 1 + 5 ** .5, 3.4985616993277016, 3.5546408627688025, 3.5666673798562816]
fig, ax = plt.subplots(figsize=(8.5, 6.2))
alpha = -2.502907875
u = np.linspace(-1, 1, 600)
cols = [RED, ORANGE, GREEN, BLUE, PURPLE]
for n, (r, col) in enumerate(zip(supers, cols)):
    scale = alpha ** n
    x_ = 0.5 + u / scale * 0.5
    y_ = x_.copy()
    for _ in range(2 ** n):
        y_ = r * y_ * (1 - y_)
    ax.plot(u, (y_ - 0.5) / 0.5 * scale, color=col, lw=1.2, label=f'f^(2^{n}) at r{n}, rescaled by α^{n}')
ax.set_ylim(-1.25, 0.45); ax.axhline(0, color='#ddd'); ax.axvline(0, color='#ddd')
ax.set_title('Period doubling = compose with yourself, then change units:\nthe rescaled curves settle on one shape (Feigenbaum)')
ax.legend(fontsize=9, loc='lower center', bbox_to_anchor=(0.5, 0.25)); ax.set_xlabel('rescaled x − ½'); ax.set_ylabel('rescaled value')
ax.text(-0.98, 0.25, 'red is the circle-dial square (r = 2); it is the first curve,\nnot the limit: the dial gives no Feigenbaum link',
        fontsize=9, color=RED)
fig.tight_layout(); fig.savefig(OUT + 'fig10-renormalization.png', bbox_inches='tight'); plt.close(fig)
print('figures written')


# ---------- fig11: turning the order dial q partway ----------
def W_rs(n, qq):
    a, b = mp.mpf(1), mp.mpf(2)
    for m in range(1, n):
        a, b = b, 2 * b - (1 - mp.mpf(qq) ** m) * a
    return b


mp.mp.dps = 30
fig, axs = plt.subplots(1, 2, figsize=(13, 5.2))
ax = axs[0]
ns = np.unique(np.round(np.logspace(0, np.log10(5000), 60)).astype(int))
for qq, col, lab in [(1, BLUE, 'q = 1 (list): rate log 2 forever'), (0.999, PURPLE, 'q = 0.999'), (0.99, GREEN, 'q = 0.99'),
                     (0.9, ORANGE, 'q = 0.9'), (0.5, RED, 'q = 0.5'), (0, GREY, 'q = 0 (bag): n + 1')]:
    ys = [math.log(n + 1) / n if qq == 0 else (math.log(2) if qq == 1 else float(mp.log(W_rs(n, qq)) / n)) for n in ns]
    ax.semilogx(ns, ys, '-', color=col, lw=2, label=lab)
ax.axhline(math.log(2), color=BLUE, ls=':', lw=1)
ax.set_xlabel('word length n'); ax.set_ylabel('(1/n) · log(weighted count of words)')
ax.set_title('Partway along the dial, the count looks like a list for a while,\nthen falls back: every q below 1 ends at rate 0')
ax.legend(fontsize=9, loc='lower left'); ax.set_ylim(0, 0.8)
ax.text(1.5, 0.74, 'log 2', color=BLUE, fontsize=9)
ax = axs[1]
ax.plot([0, 1], [0, 0], '-', color=ORANGE, lw=4)
ax.plot([1], [math.log(2)], 'o', color=BLUE, ms=11)
ax.plot([1], [0], 'o', color='white', mec=ORANGE, mew=2, ms=9)
ax.annotate('', xy=(1.45, 1.45), xytext=(1.3, 0.95), arrowprops=dict(arrowstyle='-|>', color=RED, lw=2))
ax.text(1.04, 0.62, 'q > 1: grows like q^(n²/4),\nfaster than any exponential', color=RED, fontsize=9.5)
ax.text(0.05, 0.08, 'q < 1: entropy 0 (bag-like)', color=ORANGE, fontsize=10)
ax.text(0.62, math.log(2) + 0.05, 'q = 1: log 2 (list)', color=BLUE, fontsize=10)
for val, lab in [(math.log(PHI ** 2), 'φ²: 3 maps, 2 share'), (math.log(2 + 2 ** .5), '2+√2'), (math.log(2 + 3 ** .5), '2+√3')]:
    ax.axhline(val, color=GREEN, ls='--', lw=1, xmin=0.0, xmax=0.33)
    ax.text(0.02, val + 0.02, lab, color=GREEN, fontsize=8.5)
ax.text(0.53, 1.30, 'green: in-between amounts come\nfrom which maps share a centre\n(fig 7), not from q', color=GREEN, fontsize=8.5, va='top')
ax.set_xlim(0, 1.5); ax.set_ylim(-0.05, 1.5)
ax.set_xlabel('order dial q (weight per out-of-order pair)'); ax.set_ylabel('long-run entropy')
ax.set_title('For chaos, the order dial is a switch at q = 1')
fig.tight_layout(); fig.savefig(OUT + 'fig11-order-dial-switch.png', bbox_inches='tight'); plt.close(fig)
