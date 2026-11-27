"""Figures for topics/union-find. Run: python3 figures.py (numpy, matplotlib, mpmath)."""
import math, random
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch, Circle
from ufcore import (alpha_k, ack, UF, workload, petersen, cycle, grid, merge_stages, gas_roots, components)

OUT = 'figures/'
plt.rcParams.update({'font.size': 11, 'axes.titlesize': 13, 'figure.dpi': 130})
BLUE, ORANGE, GREEN, RED, GREY, PURPLE = '#2f6db5', '#e08a1e', '#3a9a4a', '#c53b3b', '#888888', '#7a4fb0'


def arrow(ax, p, q, color='k', lw=1.4, ms=12, **kw):
    ax.add_patch(FancyArrowPatch(p, q, arrowstyle='-|>', mutation_scale=ms, color=color, lw=lw, **kw))


# ---------- fig1: climbing the ladder up, and coming down it ----------
fig, axs = plt.subplots(1, 2, figsize=(13.5, 5.6))
ax = axs[0]
rows = [('rung 1: + 2', [ack(1, m) for m in range(5)]), ('rung 2: × 2', [ack(2, m) for m in range(5)]),
        ('rung 3: 2^', [ack(3, m) for m in range(5)]), ('rung 4: tower', [ack(4, 0), ack(4, 1), '2↑↑5−3', '2↑↑6−3', '2↑↑7−3'])]
for i, (nm, vals) in enumerate(rows):
    y = 3 - i
    ax.text(-0.3, y, nm, ha='right', va='center', fontsize=11, color=[GREY, ORANGE, BLUE, PURPLE][i])
    for m, v in enumerate(vals):
        ax.text(1.25 * m, y, str(v), ha='center', va='center', fontsize=10 if i < 3 else 9,
                bbox=dict(boxstyle='round', fc='#f4f4f4', ec=[GREY, ORANGE, BLUE, PURPLE][i]))
for m in range(5):
    ax.text(1.25 * m, 4.0, f'm = {m}', ha='center', fontsize=10, color=GREY)
ax.annotate('', xy=(5.7, -0.25), xytext=(5.7, 3.25), arrowprops=dict(arrowstyle='-|>', color=RED, lw=2))
ax.text(5.85, 1.5, 'each row\niterates\nthe one\nabove', color=RED, fontsize=9.5, va='center')
ax.set_xlim(-2.6, 6.8); ax.set_ylim(-0.7, 4.4); ax.axis('off')
ax.set_title('Ackermann A(n, m) = 2 [n] (m + 3) − 3:\nthe hyperoperation ladder (our fused ladder at a = 1)')
ax = axs[1]
N80 = 10 ** 80
labs = ['undo + : n − 2', 'undo × 2 : n/2', 'undo 2^ : log₂ n', 'undo tower : log* n', 'undo rung 5', 'undo rung 6']
vals = [N80 - 2, alpha_k(1, N80), alpha_k(2, N80), alpha_k(3, N80), alpha_k(4, N80), alpha_k(5, N80)]
ylog = [math.log10(max(v, 1)) for v in vals]
cols = [GREY, ORANGE, BLUE, PURPLE, GREEN, GREEN]
ax.barh(range(len(vals))[::-1], [max(y, 0.05) for y in ylog], color=cols)
for i, (v, y) in enumerate(zip(vals, ylog)):
    txt = '10⁸⁰ − 2' if i == 0 else ('5 × 10⁷⁹' if i == 1 else str(v))
    ax.text(y + 0.8, len(vals) - 1 - i, txt, va='center', fontsize=10.5)
ax.set_yticks(range(len(vals))[::-1]); ax.set_yticklabels(labs)
ax.set_xlabel('log₁₀ of the result (bar length)'); ax.set_xlim(0, 95)
ax.set_title('Coming back down from n = 10⁸⁰: each rung\'s inverse\n(counting how often the one below must be applied)')
ax.text(30, 0.6, 'inverse Ackermann α(10⁸⁰) = 4:\nunion-find\'s cost per step', color=RED, fontsize=10.5)
fig.tight_layout(); fig.savefig(OUT + 'fig1-ladder-up-and-down.png', bbox_inches='tight'); plt.close(fig)


# ---------- fig2: union-find trees, path compression, order ----------
def draw_forest(ax, parent, title, hl=None, pos=None):
    n = len(parent)
    kids = {i: [] for i in range(n)}
    roots_ = []
    for v, p in enumerate(parent):
        (roots_ if p == v else kids[p]).append(v)
    P = {}
    xcur = [0]

    def place(v, d):
        if not kids[v]:
            P[v] = (xcur[0], -d); xcur[0] += 1
        else:
            for c in kids[v]:
                place(c, d + 1)
            P[v] = (np.mean([P[c][0] for c in kids[v]]), -d)
    for r in roots_:
        place(r, 0); xcur[0] += 0.6
    for v, p in enumerate(parent):
        if p != v:
            arrow(ax, P[v], P[p], color=GREY, lw=1.1, ms=9, shrinkA=9, shrinkB=9)
    for v, (x, y) in P.items():
        fc = '#ffe2b3' if hl and v in hl else ('#cfe0f5' if parent[v] == v else 'white')
        ax.add_patch(Circle((x, y), 0.28, fc=fc, ec='k', zorder=3))
        ax.text(x, y, str(v), ha='center', va='center', fontsize=9, zorder=4)
    ax.set_title(title, fontsize=11.5); ax.set_aspect('equal'); ax.axis('off')
    xs = [p[0] for p in P.values()]; ys = [p[1] for p in P.values()]
    ax.set_xlim(min(xs) - .6, max(xs) + .6); ax.set_ylim(min(ys) - .6, .6)


fig, axs = plt.subplots(1, 4, figsize=(15, 4.6), gridspec_kw=dict(width_ratios=[1, 1, 1.4, 1.4]))
chain = UF(7, rank=False, compress=False)
for i in range(1, 7):
    chain.union(i - 1, i)
draw_forest(axs[0], chain.p[:], 'unions in one order:\na chain (depth 6)', hl={0})
uc = UF(7, rank=False, compress=True)
uc.p = chain.p[:]
uc.find(0)
draw_forest(axs[1], uc.p, 'after find(0) with path\ncompression: depth 1', hl={0})
U = [(0, 1), (2, 3), (4, 5), (0, 2), (4, 6), (0, 4), (7, 8), (9, 8)]
a = UF(10, rank=True, compress=False)
for x, y in U:
    a.union(x, y)
draw_forest(axs[2], a.p[:], 'the same unions, by rank:\nshallow trees')
b = UF(10, rank=False, compress=False)
for x, y in U[::-1]:
    b.union(y, x)
draw_forest(axs[3], b.p[:], 'the same unions, reversed:\nother trees, same groups')
fig.suptitle('Union-find keeps groups as trees. The groups never depend on the order of the unions; the trees, and so the cost, do',
             fontsize=12.5, y=1.03)
fig.tight_layout(); fig.savefig(OUT + 'fig2-trees-and-order.png', bbox_inches='tight'); plt.close(fig)

# ---------- fig3: cost per operation ----------
fig, ax = plt.subplots(figsize=(8.5, 5.2))
ns = [2 ** e for e in range(5, 18)]
for rank, comp, adv, col, lab, nmax in [(False, False, True, RED, 'naive linking, worst order: grows like n', 2 ** 11),
                                        (True, False, False, ORANGE, 'union by rank: at most log₂ n deep', 2 ** 17),
                                        (False, True, False, PURPLE, 'path compression alone', 2 ** 17),
                                        (True, True, False, BLUE, 'rank + compression: flat (α(n) ≤ 4)', 2 ** 17)]:
    xs = [n for n in ns if n <= nmax]
    ys = [workload(n, rank, comp, 1, adversarial=adv) for n in xs]
    ax.loglog(xs, ys, 'o-', color=col, label=lab, ms=4)
ax.set_xlabel('number of elements n'); ax.set_ylabel('pointer steps per operation (average)')
ax.set_title('What union-find costs: the inverse of the hyperoperation tower')
ax.legend(fontsize=9.5); ax.grid(alpha=.3, which='both')
fig.tight_layout(); fig.savefig(OUT + 'fig3-cost.png', bbox_inches='tight'); plt.close(fig)

# ---------- fig4: one Asano merge is one union ----------
fig, axs = plt.subplots(1, 3, figsize=(14, 4.4))
ax = axs[0]
for (x1, x2, y, lab) in [(0, 1.4, 1, 'edge 1'), (2.2, 3.6, 1, 'edge 2')]:
    ax.plot([x1, x2], [y, y], color='k', lw=2)
    ax.text((x1 + x2) / 2, y + 0.25, lab, ha='center', fontsize=10)
for x, nm, col in [(0, 'x', 'white'), (1.4, 'y₁', '#ffe2b3'), (2.2, 'y₂', '#ffe2b3'), (3.6, 'u', 'white')]:
    ax.add_patch(Circle((x, 1), 0.2, fc=col, ec='k', zorder=3)); ax.text(x, 1, nm, ha='center', va='center', zorder=4)
ax.annotate('', xy=(1.75, 0.45), xytext=(1.75, 0.8), arrowprops=dict(arrowstyle='-|>', lw=2, color=RED))
ax.text(1.95, 0.6, 'union(y₁, y₂)', color=RED, fontsize=10)
ax.plot([0.4, 1.8, 3.2], [0, 0, 0], color='k', lw=2)
for x, nm, col in [(0.4, 'x', 'white'), (1.8, 'y', '#ffe2b3'), (3.2, 'u', 'white')]:
    ax.add_patch(Circle((x, 0), 0.2, fc=col, ec='k', zorder=3)); ax.text(x, 0, nm, ha='center', va='center', zorder=4)
ax.set_xlim(-.4, 4); ax.set_ylim(-.5, 1.6); ax.set_aspect('equal'); ax.axis('off')
ax.set_title('Two loose edges; one union glues\nthe two copies of y into one site')
ax = axs[1]; ax.axis('off')
ax.text(0.02, 0.92, 'Each edge is the dial factor 1 + p ⊕ₖ q = 1 + p + q + k·pq.', fontsize=10.5, transform=ax.transAxes)
ax.text(0.02, 0.76, 'Multiply the two edges, then group by y₁, y₂:', fontsize=10.5, transform=ax.transAxes)
ax.text(0.05, 0.62, '(no y)       (1 + x)(1 + u)', fontsize=10.5, transform=ax.transAxes, color=GREEN)
ax.text(0.05, 0.52, '(y₁ only)    y₁ (1 + kx)(1 + u)', fontsize=10.5, transform=ax.transAxes, color=GREY)
ax.text(0.05, 0.42, '(y₂ only)    y₂ (1 + x)(1 + ku)', fontsize=10.5, transform=ax.transAxes, color=GREY)
ax.text(0.05, 0.32, '(both)       y₁y₂ (1 + kx)(1 + ku)', fontsize=10.5, transform=ax.transAxes, color=GREEN)
ax.text(0.02, 0.16, 'Asano\'s merge: drop the grey (the copies disagree),\nturn y₁y₂ into y. That is the gas on the path x–y–u.',
        fontsize=10.5, transform=ax.transAxes, color=RED)
ax.set_title('The union, in algebra (Asano 1970)')
ax = axs[2]
th = np.linspace(0, 2 * np.pi, 400)
ax.plot(np.cos(th), np.sin(th), color=GREY, lw=1)
a = 3 ** -0.5
ws = 1.6 * np.exp(1j * th)
img = -(1 + a * ws) / (ws + a)
ax.plot(ws.real, ws.imag, color=BLUE, lw=2, label='|w₁| = 1.6 (outside)')
ax.plot(img.real, img.imag, color=ORANGE, lw=2, label='where w₂ must sit for a zero')
ax.set_aspect('equal'); ax.set_xlim(-2, 2); ax.set_ylim(-2, 2)
ax.legend(fontsize=8.5, loc='upper left')
ax.set_title('One edge (k = 3): if w₁ is outside the circle,\na zero needs w₂ inside, so both-outside is zero-free')
ax.text(-1.95, -1.15, '|w+a|² − |1+aw|² = (1 − a²)(|w|² − 1),\na = k^(−1/2) ≤ 1 exactly when k ≥ 1', fontsize=8.5)
fig.tight_layout(); fig.savefig(OUT + 'fig4-one-merge.png', bbox_inches='tight'); plt.close(fig)

# ---------- fig5: zeros stay on the circle through the whole union-find build ----------
stages = list(merge_stages(*petersen(), order_seed=7))
pick = [0, 5, 10, 15, 20]
fig, axs = plt.subplots(2, 5, figsize=(16, 6.8), gridspec_kw=dict(height_ratios=[1, 1.15]))
for j, si in enumerate(pick):
    (n_, E_), i = stages[si]
    ax = axs[0, j]
    comps = components(n_, E_)
    # draw stage graph with a spring-ish circle layout per component
    rnd = random.Random(1)
    ang = np.linspace(0, 2 * np.pi, n_, endpoint=False)
    P = {v: (np.cos(t), np.sin(t)) for v, t in enumerate(ang)}
    for u, v in E_:
        ax.plot([P[u][0], P[v][0]], [P[u][1], P[v][1]], color='k', lw=0.8)
    for v in range(n_):
        ax.add_patch(Circle(P[v], 0.06, fc='#cfe0f5', ec='k', zorder=3))
    ax.set_aspect('equal'); ax.axis('off'); ax.set_xlim(-1.2, 1.2); ax.set_ylim(-1.2, 1.2)
    pc = len(comps)
    ax.set_title(f'after {i} unions:\n{n_} sites, {pc} piece' + ('s' if pc != 1 else ''), fontsize=10.5)
    ax = axs[1, j]
    rs = np.array(gas_roots(n_, E_, 3.0))
    ax.plot(np.cos(th), np.sin(th), color=GREY, lw=1)
    ax.plot(rs.real, rs.imag, 'o', color=BLUE, ms=5)
    ax.set_aspect('equal'); ax.set_xlim(-1.3, 1.3); ax.set_ylim(-1.3, 1.3)
    distinct = len({(round(r.real, 6), round(r.imag, 6)) for r in rs})
    ax.set_title(f'{len(rs)} zeros, all on |w| = 1' + (f'\n({distinct} distinct points, rest repeat)' if distinct < len(rs) else ''),
                 fontsize=10.5)
    ax.set_xticks([]); ax.set_yticks([])
fig.suptitle('Building the Petersen graph by union-find from 15 loose edges (dial k = 3). Every union removes one zero; '
             'the rest slide along the circle and never leave it', fontsize=12.5)
fig.tight_layout(); fig.savefig(OUT + 'fig5-build-on-the-circle.png', bbox_inches='tight'); plt.close(fig)

# ---------- fig6: the dial k decides ----------
fig, axs = plt.subplots(1, 4, figsize=(15.5, 4.3))
for ax, k, col, ttl in zip(axs, [0.3, 0.7, 1.0, 4.0], [RED, RED, GREY, BLUE],
                           ['k = 0.3: repulsive', 'k = 0.7: repulsive', 'k = 1: free gas', 'k = 4: attractive']):
    ax.plot(np.cos(th), np.sin(th), color=GREY, lw=1)
    xmin = -1.8
    for (N, E), mk in [(petersen(), 'o'), (cycle(9), 's'), (grid(3, 3), '^')]:
        if k == 1.0:
            rs = np.array([-1.0 + 0j])
        else:
            rs = np.array(gas_roots(N, E, k))
        xmin = min(xmin, rs.real.min() - 1.5)
        ax.plot(rs.real, rs.imag, mk, color=col, ms=5, mfc='none' if mk != 'o' else col)
    W = 1.5 - xmin
    ax.set_aspect('equal'); ax.set_xlim(xmin, 1.5); ax.set_ylim(-0.5 * W, 0.5 * W)
    ax.set_title(ttl + ('\nall N zeros at w = −1' if k == 1 else ''), fontsize=11)
    if k < 1:
        ax.text(xmin + 0.04 * W, 0.3 * W, 'zeros leave the circle\nand spread along the\nnegative real axis',
                fontsize=9, color=RED)
        ax.text(xmin + 0.04 * W, -0.45 * W, 'o Petersen  □ ring C9  △ grid 3×3', fontsize=8)
fig.suptitle('The dial k is the switch: the Lee–Yang circle holds for k ≥ 1 (ferromagnetic) and fails below (w = scaled fugacity)',
             fontsize=12.5)
fig.tight_layout(); fig.savefig(OUT + 'fig6-dial-switch.png', bbox_inches='tight'); plt.close(fig)

# ---------- fig7: which moves keep Lee-Yang, and the RH picture ----------
fig, ax = plt.subplots(figsize=(14, 5.4)); ax.axis('off')
cells = [('product\n(two systems side by side)', 'keeps the circle\n(zeros: both lists)', GREEN, 'no union',
          'finite primon gas: products only,\nno unions, so no zeros at all', ORANGE),
         ('union with an attractive link\n(Asano merge, k ≥ 1)', 'keeps the circle', GREEN, "union-find's union",
          'the missing bridge would be built here:\nattractive unions, field on the total,\nconverging at σ = ½ (nobody has it)', BLUE),
         ('union with a repulsive link\n(k < 1)', 'breaks it', RED, 'union, wrong sign', '', GREY),
         ('sum / mixture of two systems', 'breaks it\n(e.g. (1 + w)⁴ + (1 − w)⁴)', RED, 'no union-find move at all',
          "Newman's spin Φ is a sum of\nshifted bumps, so it starts here", ORANGE)]
for i, (move, eff, col, uf, rh, rcol) in enumerate(cells):
    y = 0.8 - i * 0.21
    ax.text(0.0, y, move, fontsize=11, va='center', transform=ax.transAxes)
    ax.text(0.27, y, eff, fontsize=11, va='center', color=col, transform=ax.transAxes)
    ax.text(0.46, y, uf, fontsize=10.5, va='center', color=GREY, transform=ax.transAxes)
    ax.text(0.66, y, rh, fontsize=10, va='center', color=rcol, transform=ax.transAxes)
for x, h in [(0.0, 'move'), (0.27, 'Lee–Yang circle'), (0.46, 'in union-find terms'), (0.66, 'the Riemann route (spin-build page)')]:
    ax.text(x, 0.97, h, fontsize=11.5, weight='bold', transform=ax.transAxes)
fig.suptitle('Of these four moves, only products and attractive unions keep the zeros on the circle', fontsize=13)
fig.savefig(OUT + 'fig7-moves.png', bbox_inches='tight'); plt.close(fig)
print('figures written')
