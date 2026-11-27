"""Figures for topics/entropy-knots. Run after checks.py (reads windows.json, zeros.json, posbraid.json).
Colour key: orange = known, blue = known pieces in our reading, green = not found in print."""
import json, math
import numpy as np, mpmath as mp, sympy as sp
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Circle
from windows import F, orbit, kneading_entropy

R = json.load(open('windows.json')); Z = json.load(open('zeros.json')); PB = json.load(open('posbraid.json'))
phi = (1 + 5 ** .5) / 2
ORANGE, BLUE, GREEN, GREY = '#d9822b', '#2b6cb0', '#2f855a', '#718096'
pcol = {2: '#999999', 3: '#e53e3e', 4: '#dd6b20', 5: '#d69e2e', 6: '#38a169', 7: '#319795', 8: '#3182ce', 9: '#805ad5'}
plt.rcParams.update({'font.size': 10, 'axes.spines.top': False, 'axes.spines.right': False})
def save(fig, name):
    fig.tight_layout(); fig.savefig('figures/' + name, dpi=140); plt.close(fig)

# ---------- fig1: entropy staircase with every stable window on it ----------
mp.mp.dps = 40
rs = np.linspace(3.4, 4.0, 1200)
def h_curve(r):
    # a truncated kneading series can have a spurious zero near t = 1; keep a zero only if it is
    # stable between two truncation lengths
    a = kneading_entropy(mp.mpf(r), N=110)[0]; b = kneading_entropy(mp.mpf(r), N=170)[0]
    return a if abs(a - b) < 2e-3 else 0.0
hs = [h_curve(r) for r in rs]
fig, axes = plt.subplots(1, 2, figsize=(12, 4.6), gridspec_kw={'width_ratios': [2.2, 1]})
ax = axes[0]
ax.plot(rs, hs, color=GREY, lw=1.2, label='map entropy h(r) (kneading determinant)')
for n in range(2, 10):
    pts = [w for w in R if w['n'] == n]
    ax.scatter([w['r'] for w in pts], [w['h_pattern'] for w in pts], s=22, color=pcol[n], zorder=3, label='period %d window' % n)
for w in R:
    if w['n'] <= 5:
        ax.annotate(''.join(map(str, w['perm'])), (w['r'], w['h_pattern']), xytext=(-8, 8), textcoords='offset points', fontsize=7, color=pcol[w['n']])
ax.set_xlabel('r in r·x(1 − x)'); ax.set_ylabel('entropy'); ax.set_title('Each stable window sits where map entropy = the entropy its orbit pattern forces')
ax.legend(fontsize=7, ncol=2, loc='upper left', frameon=False)
ax = axes[1]
m = (rs > 3.82) & (rs < 3.865)
ax.plot(rs[m], np.array(hs)[m], color=GREY, lw=1.5)
for w in R:
    if 3.82 < w['r'] < 3.865:
        ax.scatter(w['r'], w['h_pattern'], color=pcol[w['n']], s=40, zorder=3)
        ax.annotate('period %d' % w['n'], (w['r'], w['h_pattern']), xytext=(0, 10 + 8 * (w['n'] % 2)), textcoords='offset points', ha='center', fontsize=8)
ax.axhline(math.log(phi), color=ORANGE, ls=':', lw=1)
ax.text(3.8215, math.log(phi) + 0.004, 'log φ plateau', color=ORANGE, fontsize=8)
ax.axvspan(1 + 8 ** .5, 3.8568, color=ORANGE, alpha=0.08)
ax.set_title('Zoom: the period-3 window is one plateau'); ax.set_xlabel('r'); ax.set_ylim(0.44, 0.53)
save(fig, 'fig1-entropy-staircase.png')

# ---------- fig2: entropy is the nearest zero ----------
fig, axes = plt.subplots(1, 3, figsize=(12, 4))
for ax, perm in zip(axes, [[1, 2, 0], [1, 3, 4, 2, 0], [2, 4, 6, 7, 5, 3, 1, 0]]):
    w = [x for x in R if x['perm'] == perm][0]
    th = w['kneading']
    roots = np.polynomial.Polynomial([float(v) for v in th]).roots()
    ax.add_patch(Circle((0, 0), 1, fill=False, color=GREY, ls='--', lw=0.8))
    ax.scatter(roots.real, roots.imag, color=BLUE, s=30, zorder=3, label='zeros of the kneading polynomial')
    z = math.exp(-w['h_pattern'])
    ax.scatter([z], [0], s=160, facecolor='none', edgecolor='red', lw=2, zorder=4, label='nearest zero = e^(−h)')
    ax.set_aspect('equal'); ax.set_xlim(-1.6, 1.6); ax.set_ylim(-1.6, 1.6)
    ax.axhline(0, color='#ddd', lw=0.5); ax.axvline(0, color='#ddd', lw=0.5)
    ax.set_title('period %d, pattern %s\nh = %.4f, e^(−h) = %.4f' % (w['n'], ''.join(map(str, perm)), w['h_pattern'], z), fontsize=9)
axes[0].legend(fontsize=7, loc='lower left', frameon=False)
fig.suptitle('Entropy is read off the zero nearest the origin (Milnor–Thurston). Dashed: unit circle', fontsize=10)
save(fig, 'fig2-entropy-is-a-zero.png')

# ---------- fig3: line entropy vs plane (braid) entropy ----------
fig, ax = plt.subplots(figsize=(7, 6))
ax.plot([0, 0.72], [0, 0.72], color=GREY, lw=0.8, ls='--')
for w in R:
    hb = w['h_braid'] if w['h_braid'] > 0.01 else 0.0
    ax.scatter(w['h_pattern'], hb, color=pcol[w['n']], s=28, zorder=3)
for n in range(2, 10): ax.scatter([], [], color=pcol[n], label='period %d' % n)
w3 = [w for w in R if w['n'] == 3][0]
ax.annotate('period 3: a rotation braid,\nforces nothing in the plane', (w3['h_pattern'], 0), xytext=(0.2, 0.2), arrowprops=dict(arrowstyle='->', color=GREY), fontsize=8)
w5 = [w for w in R if w['perm'] == [1, 3, 4, 2, 0]][0]
ax.annotate('period 5: smallest 5-strand\ndilatation 1.7221, no loss', (w5['h_pattern'], w5['h_braid']), xytext=(0.18, 0.6), arrowprops=dict(arrowstyle='->', color=GREY), fontsize=8)
ax.text(0.45, 0.66, 'on the diagonal: the braid\nforces all the line entropy', fontsize=8, color=GREY)
ax.set_xlabel('entropy forced on a line (orbit pattern)'); ax.set_ylabel('entropy forced in the plane (orbit braid)')
ax.set_title('65 stable-window orbits: the braid never forces more, often less')
ax.legend(fontsize=7, frameon=False, loc='upper left')
save(fig, 'fig3-line-vs-plane.png')

# ---------- fig4: two orbit braids drawn ----------
def draw_braid(ax, word, n, title):
    xs = list(range(n)); y = 0
    pos = list(range(n))
    for k, s in enumerate(word):
        i = abs(s) - 1
        for j in range(n):
            if j not in (i, i + 1): ax.plot([j, j], [y, y - 1], color='k', lw=1.5)
        over = (i, i + 1) if s > 0 else (i + 1, i)
        ax.plot([over[1], over[0]], [y, y - 1], color='k', lw=1.5, zorder=2)
        tt = np.linspace(0, 1, 30); xa = over[0] + (over[1] - over[0]) * tt; ya = y - tt
        msk = np.abs(tt - 0.5) > 0.12
        ax.plot(np.where(msk, xa, np.nan), ya, color='k', lw=1.5, zorder=1)
        y -= 1
    ax.set_xlim(-0.5, n - 0.5); ax.set_ylim(y - 0.3, 0.3); ax.axis('off'); ax.set_title(title, fontsize=9)
fig, axes = plt.subplots(1, 3, figsize=(11, 5.5), gridspec_kw={'width_ratios': [1, 1.4, 1.6]})
draw_braid(axes[0], w3['braid'], 3, 'period 3 window orbit\nbraid %s: entropy 0' % w3['braid'])
draw_braid(axes[1], w5['braid'], 5, 'period 5 window [1,3,4,2,0]\nentropy %.4f = log 1.7221' % w5['h_braid'])
w7 = [w for w in R if w['perm'] == [1, 3, 4, 5, 6, 2, 0]][0]
draw_braid(axes[2], w7['braid'], 7, 'period 7 window [1,3,4,5,6,2,0]\nentropy %.4f = log 1.4656 (line: %.4f)' % (w7['h_braid'], w7['h_pattern']))
fig.suptitle('Braids traced by stable orbits of the folded (horseshoe) logistic map; time runs down', fontsize=10)
save(fig, 'fig4-orbit-braids.png')

# ---------- fig5: a stable and an unstable orbit with the same knot type ----------
mp.mp.dps = 30
r = 3.835
xs = np.linspace(0, 1, 2000)
f3 = xs.copy()
for _ in range(3): f3 = r * f3 * (1 - f3)
fig, ax = plt.subplots(figsize=(7, 6))
ax.plot(xs, f3, color=GREY, lw=1, label='f³ at r = 3.835')
ax.plot([0, 1], [0, 1], color='#bbb', lw=0.8)
def per3(x0):
    x = mp.findroot(lambda y: F(mp.mpf(r), F(mp.mpf(r), F(mp.mpf(r), y))) - y, mp.mpf(x0))
    o = [x, F(mp.mpf(r), x), F(mp.mpf(r), F(mp.mpf(r), x))]
    m_ = 1
    for y in o: m_ *= r * (1 - 2 * y)
    return [float(v) for v in o], float(m_)
seen = []
for x0 in [0.15, 0.17, 0.5, 0.52, 0.95, 0.96]:
    o, m_ = per3(x0)
    if any(abs(min(o) - min(s[0])) < 1e-9 for s in seen) or abs(o[0] - o[1]) < 1e-6: continue
    seen.append((o, m_))
for (o, m_), c in zip(seen, [GREEN, 'red']):
    lab = ('stable' if abs(m_) < 1 else 'unstable') + ' period-3 orbit, multiplier %.2f' % m_
    ax.scatter(o, o, s=60, color=c, zorder=3, label=lab)
    xx = o + [o[0]]
    for a, b in zip(xx, xx[1:]): ax.annotate('', (b, b), (a, a), arrowprops=dict(arrowstyle='->', color=c, lw=0.8, connectionstyle='arc3,rad=0.3'))
ax.set_title('Same pattern 1-2-0 (same knot type), opposite stability:\nborn together in one saddle-node at r = 1 + √8')
ax.legend(fontsize=8, frameon=False, loc='upper left'); ax.set_xlabel('x'); ax.set_ylabel('f³(x)')
save(fig, 'fig5-same-knot-opposite-stability.png')

# ---------- fig6: hard-core zeros at the Jones index ----------
fig, axes = plt.subplots(1, 2, figsize=(13, 4.8), gridspec_kw={'width_ratios': [1.3, 1]})
ax = axes[0]
for m in range(1, 9):
    zs = [z for z in Z['P%d' % m]['zeros']]
    ax.scatter(zs, [m] * len(zs), color=BLUE, s=25)
    ax.scatter([min(zs)], [m], s=110, facecolor='none', edgecolor='red', lw=1.6)
    ax.text(1.02, (m - 0.5) / 8.4 + 0.02, 'growth %.4f = 4cos²(π/%d)' % (Z['P%d' % m]['growth'], m + 2), va='center', fontsize=8, transform=ax.transAxes)
ax.set_xlim(0, 2.8); ax.set_ylim(0.5, 8.9); ax.set_yticks(range(1, 9)); ax.set_ylabel('m letters in a row (path P_m)')
ax.set_xlabel('t, i.e. minus the fugacity of a hard-core gas on the path')
ax.set_title('Zeros of the trace-monoid polynomial μ(t): all real;\nthe nearest (red) is e^(−entropy)', fontsize=10)
ax = axes[1]
dd = np.linspace(0, 2.05, 800)
for m, c in zip([2, 3, 4, 6], ['#a0aec0', ORANGE, BLUE, GREEN]):
    U = [float(sp.chebyshevu(m + 1, v / 2)) for v in dd]
    ax.plot(dd, U, color=c, lw=1.4, label='Jones–Wenzl trace Δ_%d(d), path P_%d' % (m + 1, m))
ax.axhline(0, color='k', lw=0.6)
for n in range(3, 13):
    ax.axvline(2 * math.cos(math.pi / n), color='#ddd', lw=0.7, zorder=0)
ax.text(2 * math.cos(math.pi / 5), 2.75, 'φ', ha='center', fontsize=9, color=ORANGE)
ax.text(2 * math.cos(math.pi / 6), 2.75, '√3', ha='center', fontsize=9, color=BLUE)
ax.text(2 * math.cos(math.pi / 4), 2.75, '√2', ha='center', fontsize=9, color=GREY)
ax.set_ylim(-2.2, 3); ax.set_xlabel('loop value d = −A² − A⁻² (the knot page\'s dial weight)')
ax.set_title('Same zeros in the loop value: t = 1/d² at d = 2cos(jπ/(m+2)),\nwhere Jones–Wenzl projectors fail (grey lines: 2cos(π/n))', fontsize=10)
ax.legend(fontsize=7, frameon=False, loc='upper left')
save(fig, 'fig6-zeros-at-jones-index.png')

# ---------- fig7: affine centre patterns, real vs complex zeros ----------
fig, ax = plt.subplots(figsize=(7.5, 5.5))
marks = ['o', 's', '^', 'D', 'v', 'P', 'X']
for p, mk in zip(Z['multipartite'], marks):
    zs = np.array(p['zeros'])
    lab = 'classes %s: growth %.4f %s' % (p['sizes'], p['growth'], ('= B_%d' % p['beraha_n']) if p['is_beraha'] else '(not Beraha)')
    ax.scatter(zs[:, 0], zs[:, 1], marker=mk, s=45, label=lab, color=GREEN if not p['real_zeros'] else BLUE)
ax.axhline(0, color='#ccc', lw=0.6); ax.set_xlabel('Re t'); ax.set_ylabel('Im t')
ax.set_title('Affine maps grouped by shared centre: zeros of μ(t)\nblue all real; green: the claw [3,1] and [4,1] leave the axis', fontsize=10)
ax.legend(fontsize=7, frameon=False)
save(fig, 'fig7-centre-patterns-zeros.png')

# ---------- fig8: growth ladder ----------
fig, ax = plt.subplots(figsize=(8, 5))
ms = np.arange(1, 9)
ax.plot(ms + 1, [math.log(m) if m > 1 else 0 for m in ms], 'o-', color=GREY, label='free: no letters commute (m letters)')
ax.plot(ms + 1, [math.log(Z['P%d' % m]['growth']) for m in ms], 'o-', color=BLUE, label='trace monoid, neighbours only clash (strand maps)')
ax.plot([3, 4, 5], [math.log(PB['B%d' % n]['growth']) for n in [3, 4, 5]], 's-', color=ORANGE, label='positive braids B_n⁺ (add the braid relation)')
ax.axhline(math.log(4), color='#ccc', ls=':'); ax.text(8.6, math.log(4) + 0.02, 'log 4 (Jones index 4)', fontsize=8, ha='right')
ax.set_xlabel('strands (m letters act on m + 1 strands)'); ax.set_ylabel('entropy (log growth of distinct composites)')
ax.set_title('How much each commuting rule cuts entropy')
ax.legend(fontsize=8, frameon=False)
save(fig, 'fig8-growth-ladder.png')
print('figures written')
