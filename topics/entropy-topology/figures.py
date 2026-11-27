"""Figures for topics/entropy-topology (DRAFT). Run: python3 figures.py  (about 1 minute).
Colour key: orange = known, blue = known pieces in our reading, green = not found in print."""
import json, math
from math import gcd
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from dialtop import (dial_power_h, rp1_degree_alg, newton_parts, distinct_real_roots,
                     newton_preimage_counts)
from orbits import (logistic_iter, orbits_of_exact_period, interval_index, birkhoff_pair,
                    linear_type, stdmap, stdmap_n)
from morse_graph import morse_graph, set_centre

OR, BL, GR, RD, GY = '#e08214', '#2166ac', '#1a9850', '#b2182b', '#777777'
plt.rcParams.update({'font.size': 10, 'axes.spines.top': False, 'axes.spines.right': False})

def save(fig, name):
    fig.tight_layout(); fig.savefig('figures/' + name, dpi=140); plt.close(fig)

# ---- fig 1: the dial plane sorted by what squaring does ----
fig, ax = plt.subplots(1, 3, figsize=(13, 4.2))
kk = np.linspace(-3, 3, 400)
a = ax[0]
a.fill_between(kk, kk ** 2 / 4, 3, color=BL, alpha=.15)
a.fill_between(kk, -3, kk ** 2 / 4, color=OR, alpha=.15)
a.plot(kk, kk ** 2 / 4, 'k', lw=1.5)
a.text(-1.2, 1.9, 'Δ < 0: zeros complex\nno stable orbit\ndegree 2, entropy log 2', ha='center', color=BL)
a.text(-1.9, -2.3, 'Δ > 0: two real zeros\n= two stable orbits\ndegree 0, entropy 0', ha='left', color=OR)
a.annotate('wall Δ = 0 (degree 1)', (-2.4, 1.44), (-2.9, 2.7), fontsize=9, arrowprops=dict(arrowstyle='->', lw=.8))
a.plot([0], [1], 'o', color=BL); a.annotate('tangent law', (0, 1), (0.25, 0.6), fontsize=9)
a.plot([0], [-1], 's', color=OR); a.annotate('velocity addition\n(c = 1)', (0, -1), (0.25, -1.4), fontsize=9)
a.plot([1], [0], '^', color=OR); a.annotate('circle op ∘₁', (1, 0), (1.2, 0.2), fontsize=9)
a.set_xlabel('k'); a.set_ylabel('d'); a.set_xlim(-3, 3); a.set_ylim(-3, 3)
a.set_title('Squaring sorts the dial plane by Δ = k² − 4d')

def cobweb(a, f, x0, n, col):
    x = x0
    for _ in range(n):
        y = f(x)
        a.plot([x, x], [x, y], color=col, lw=.8); a.plot([x, y], [y, y], color=col, lw=.8)
        x = y
xs = np.linspace(-4, 4, 2000)
a = ax[1]
fh = lambda t: (2 * t + t * t) / (1 + t * t)          # k = 1, d = -1
a.plot(xs, fh(xs), color=OR); a.plot(xs, xs, color=GY, lw=.8)
for x0 in (0.15, -0.4, 3.5):
    cobweb(a, fh, x0, 8, BL)
z1, z2 = (1 + 5 ** .5) / 2, (1 - 5 ** .5) / 2
a.plot([z1, z2], [z1, z2], 'o', color=GR, ms=8, zorder=5); a.plot([0], [0], 'x', color=RD, ms=9, mew=2, zorder=5)
a.text(z1 + .3, z1 + .5, 'zero 1.618\nstable', color=GR); a.text(z2 - 1.5, z2 + .4, 'zero −0.618\nstable', color=GR)
a.text(0.2, -0.6, '0: repeller', color=RD)
a.set_xlim(-4, 4); a.set_ylim(-2, 3); a.set_xlabel('x'); a.set_ylabel('x ∘ x')
a.set_title('k = 1, d = −1: orbits fall onto the\nzeros of 1 + kx + dx² (Newton\'s method)')
a = ax[2]
th = np.linspace(0, np.pi, 4000, endpoint=False)
fe = lambda t: 2 * t / (1 - t * t)                    # k = 0, d = 1
x = 0.3; orb = []
for _ in range(400):
    orb.append(math.atan(x) % math.pi); x = fe(x)
a.plot(th, np.mod(2 * th, np.pi), '.', ms=.6, color=BL)
a.plot(orb[:-1], orb[1:], '.', color=RD, ms=3)
a.set_xlabel('angle θ on RP¹ (x = tan θ)'); a.set_ylabel('angle of x ∘ x')
a.set_title('k = 0, d = 1 (tangent law): angle doubles,\nno stable orbit (red: 400 steps of one\norbit spread over the whole circle)')
save(fig, 'fig1-dial-plane-newton.png')

# ---- fig 2: circle degree and Lefschetz ----
fig, ax = plt.subplots(1, 3, figsize=(13, 4))
for a, (k, d, title) in zip(ax, [(0, 1, 'Δ < 0: degree 2, L = 1 − 2 = −1'), (2, 1, 'Δ = 0: degree 1, L = 0'), (0, -1, 'Δ > 0: degree 0, L = 1 − 0 = 1')]):
    th = np.linspace(0, np.pi, 20002)   # even spacing that skips pi/4 and 3pi/4, so no fixed point sits on a grid point
    U, V = dial_power_h(np.sin(th), np.cos(th), k, d, 2)
    phi = np.unwrap(2 * np.arctan2(U, V)) / 2
    phi -= phi[0] - (np.arctan2(U[0], V[0]) % np.pi)
    a.plot(th, phi, color=BL, lw=2, label='image angle, followed round')
    for j in range(-1, 3):
        a.plot(th, th + j * np.pi, color=GY, lw=.6)
    g = phi - th
    s = np.sign(np.sin(g)); idx = np.nonzero(s[:-1] * s[1:] < 0)[0]
    for i in idx:
        slope = (phi[i + 1] - phi[i]) / (th[i + 1] - th[i])
        ind = int(np.sign(1 - slope))
        a.plot(th[i], phi[i], 'o', color=GR if ind > 0 else RD, ms=8)
        a.annotate('%+d' % ind, (th[i], phi[i]), (th[i] + .05, phi[i] + .25), color=GR if ind > 0 else RD)
    a.plot(0, 0, 'o', color=RD, ms=8)
    a.set_title(title); a.set_xlabel('θ'); a.set_ylabel('image angle (lifted)')
ax[0].legend(loc='upper left', fontsize=8)
fig.suptitle('The circle degree is how many half-turns the image makes; fixed points (where the curve meets a grey diagonal) have index +1 (green, attracting) or −1 (red, repelling), and they add to L = 1 − degree', fontsize=10)
save(fig, 'fig2-degree-lefschetz.png')

# ---- fig 3: Conley-Morse graph along the curvature dial ----
fig, ax = plt.subplots(figsize=(9, 4.6))
for dd in np.linspace(-1, 1, 81):
    S, R, K, edges = morse_graph(lambda u, v: dial_power_h(u, v, 1, dd, 2), ncell=900, samples=12)
    for s, kd in zip(S, K):
        if len(s) > 400:
            ax.plot([dd, dd], [0, np.pi], color=BL, lw=3, alpha=.35, solid_capstyle='butt')
        else:
            c = math.atan(set_centre(edges, s)) % math.pi
            ax.plot(dd, c, 'o' if kd == 'attractor' else 'x', color=GR if kd == 'attractor' else RD, ms=4)
ax.axvline(.25, color='k', ls='--', lw=1); ax.text(.26, 2.9, 'wall d = 1/4 (Δ = 0)')
ax.text(-0.95, 2.85, 'green: attractors = real zeros of 1 + x + dx²\nred: repeller at 0', fontsize=9)
ax.text(.45, 1.4, 'one Morse set:\nthe whole circle;\nits Conley index is\nH₁(S¹) with the map ×2', color=BL, backgroundcolor='white')
ax.set_xlabel('curvature dial d   (k = 1)'); ax.set_ylabel('point of RP¹ (angle θ, x = tan θ)')
ax.set_title('Conley–Morse graph of x ∘ x on a 2000-cell grid (k = 1):\nthe real zeros are the stable Morse sets until they leave the line at d = 1/4')
save(fig, 'fig3-morse-graph-dial.png')

# ---- fig 4: Nielsen counts and entropy ----
fig, ax = plt.subplots(1, 2, figsize=(12, 4.2))
j = np.arange(1, 11)
ax[0].semilogy(j, 2 ** j - 1, 'o-', color=BL, label='Δ < 0: 2^j − 1 (= Nielsen number |L(f^j)|)')
ax[0].semilogy(j, 3 * np.ones_like(j), 's-', color=OR, label='Δ > 0: 3 (0 and the two zeros)')
ax[0].set_xlabel('j'); ax[0].set_ylabel('points of period dividing j'); ax[0].legend(fontsize=9)
ax[0].set_title('Periodic points of x ∘ x on RP¹ (check B3)')
cases = [('y² + 1, r = 0', [1, 0, 1], BL), ('quartic, r = 4', np.poly([-2, -.5, .7, 1.9]), GR),
         ('quintic, r = 5', np.poly([-2, -1, .2, 1.1, 2.3]), OR), ('sextic, r = 6', np.poly([-2.5, -1.4, -.3, .6, 1.5, 2.7]), RD)]
for name, cf, col in cases:
    cnt = newton_preimage_counts(cf, 12, cap=60000)
    r = len(distinct_real_roots(cf))
    ax[1].semilogy(range(len(cnt)), cnt, 'o-', color=col, label='%s: slope log %d' % (name, abs(2 - r)))
ax[1].set_xlabel('backward steps j'); ax[1].set_ylabel('|N^{-j}(∞)|')
ax[1].legend(fontsize=9); ax[1].set_title('Real Newton maps: the backward orbit grows at log|2 − r| (check C3)')
save(fig, 'fig4-nielsen-entropy.png')

# ---- fig 5: degree = 2 - r ----
rng = np.random.default_rng(3)
pts = []
for n in range(2, 11):
    for _ in range(120):
        r = int(rng.integers(0, n + 1));  r -= (n - r) % 2; r = max(r, 0)
        roots = list(rng.normal(size=r) * 2)
        for _ in range((n - r) // 2):
            a_, b_ = rng.normal(), abs(rng.normal()) + .2; roots += [a_ + 1j * b_, a_ - 1j * b_]
        cf = np.real(np.poly(roots))
        pts.append((len(distinct_real_roots(cf)), rp1_degree_alg(*newton_parts(cf)[:2]), n))
pts = np.array(pts)
fig, ax = plt.subplots(1, 2, figsize=(12, 4.3))
ax[0].scatter(pts[:, 0] + rng.normal(0, .08, len(pts)), pts[:, 1] + rng.normal(0, .08, len(pts)), c=pts[:, 2], cmap='viridis', s=6)
rr = np.arange(0, 11); ax[0].plot(rr, 2 - rr, color=GR, lw=1.5, label='degree = 2 − r')
ax[0].set_xlabel('distinct real roots r'); ax[0].set_ylabel('circle degree of real Newton map')
ax[0].legend(); ax[0].set_title('1,080 random polynomials (colour = degree n), every one on the line')
fl = np.log(np.maximum(np.abs(2 - rr), 1)); ax[1].plot(rr, fl, 'o-', color=OR, label='Lefschetz floor log|2 − r|')
ax[1].plot(rr[2:], np.log(rr[2:]), '--', color=GY, label='ceiling log n when n = r (all roots real)')
for n_, h_ in [(4, .6931), (5, 1.0986), (6, 1.3863)]:
    ax[1].plot(n_, h_, '*', color=GR, ms=12)
ax[1].plot(1, .609, '*', color=RD, ms=12); ax[1].annotate('y³ − 2y + 2: 0.61 > floor 0\n(folds add entropy)', (1, .609), (0.3, 1.25), color=RD, fontsize=9)
ax[1].text(4.2, .35, 'stars: measured; all-real\npolynomials sit on the floor', color=GR, fontsize=9)
ax[1].set_xlabel('distinct real roots r'); ax[1].set_ylabel('entropy'); ax[1].legend(fontsize=8, loc='upper left')
ax[1].set_title('Entropy against the floor')
save(fig, 'fig5-newton-degree.png')

# ---- fig 6: same knot type, opposite index ----
fig, ax = plt.subplots(1, 2, figsize=(12, 4.3))
r = 3.835
xs = np.linspace(0.44, 0.58, 2000)
g = xs - logistic_iter(xs, r, 3)
ax[0].plot(xs, g, color=BL); ax[0].axhline(0, color=GY, lw=.8)
for orb, m in orbits_of_exact_period(r, 3):
    p = [z for z in orb if .44 < z < .58][0]
    ind = int(np.sign(1 - m))
    ax[0].plot(p, 0, 'o', color=GR if ind > 0 else RD, ms=9)
    ax[0].annotate('multiplier %.2f\nindex %+d' % (m, ind), (p, 0), (p - .005 if ind > 0 else p + .003, -.03 if ind > 0 else -.06), color=GR if ind > 0 else RD)
ax[0].axvspan(.44, .4415, color=OR, alpha=.5); ax[0].axvspan(.5785, .58, color=OR, alpha=.5)
ax[0].text(.445, -.16, 'ends of the isolating interval:\nsame sign, so total index 0', color=OR, fontsize=9)
ax[0].set_xlabel('x'); ax[0].set_ylabel('x − f³(x)')
ax[0].set_title('r = 3.835: both period-3 orbits visit middle, right, left\n(one knot type), but their indices are +1 and −1')
rs = np.linspace(3.80, 3.86, 121)
tot = [interval_index(lambda z: logistic_iter(z, rr_, 3), .44, .58) for rr_ in rs]
plus = []; minus = []
for rr_ in rs:
    ob = orbits_of_exact_period(rr_, 3)
    plus.append(sum(1 for o, m in ob if m < 1)); minus.append(-sum(1 for o, m in ob if m > 1))
ax[1].step(rs, plus, color=GR, label='stable orbit (+1)'); ax[1].step(rs, minus, color=RD, label='unstable twin (−1)')
ax[1].plot(rs, tot, 'k', lw=2, label='Conley/fixed-point index of the interval')
ax[1].axvline(1 + 8 ** .5, color=GY, ls='--'); ax[1].text(1 + 8 ** .5 + .001, .45, 'birth at 1 + √8')
ax[1].set_xlabel('r'); ax[1].set_ylabel('index'); ax[1].legend(fontsize=8, loc='lower right')
ax[1].set_title('The pair is born from nothing because their indices cancel')
save(fig, 'fig6-same-knot-opposite-index.png')

# ---- fig 7: standard map, Birkhoff pairs ----
K = 0.9
fig, ax = plt.subplots(1, 2, figsize=(12.5, 5))
a = ax[0]
rng = np.random.default_rng(7)
for _ in range(120):
    z = np.array([rng.uniform(0, 2 * np.pi), rng.uniform(-np.pi, np.pi)])
    pts_ = []
    for _ in range(300):
        z = stdmap(z, K); pts_.append([z[0] % (2 * np.pi), (z[1] + np.pi) % (2 * np.pi) - np.pi])
    pts_ = np.array(pts_); a.plot(pts_[:, 0], pts_[:, 1], ',', color=GY, alpha=.6)
for (m, n) in [(0, 1), (1, 3), (2, 5), (1, 2)]:
    for o in birkhoff_pair(K, m, n):
        O = o['orbit']; col = {'elliptic': GR, 'hyperbolic': RD, 'inverse hyperbolic': OR}[linear_type(o['trace'])]
        a.plot(O[:, 0] % (2 * np.pi), (O[:, 1] + np.pi) % (2 * np.pi) - np.pi, 'o' if col == GR else 'X', color=col, ms=7)
a.set_xlabel('x'); a.set_ylabel('p'); a.set_title('Standard map K = 0.9: Birkhoff pairs 0/1, 1/3, 2/5, 1/2\ngreen o elliptic (+1, +1), red X hyperbolic (−1, −1)')
a = ax[1]
pr = birkhoff_pair(K, 1, 3)
for o, off in zip(pr, (0, 1)):
    z0 = o['orbit'][0]
    Fn = lambda w: stdmap_n(w, K, 3) - np.array([2 * np.pi, 0])
    t = np.linspace(0, 2 * np.pi, 25)[:-1]; rho = .03
    circ = z0 + rho * np.stack([np.cos(t), np.sin(t)], 1)
    v = Fn(circ) - circ; v = v / np.linalg.norm(v, axis=1)[:, None]
    cx = np.cos(t) + 3 * off; cy = np.sin(t)
    a.plot(np.cos(np.linspace(0, 2 * np.pi, 100)) + 3 * off, np.sin(np.linspace(0, 2 * np.pi, 100)), color=GY, lw=.6)
    a.quiver(cx, cy, v[:, 0], v[:, 1], color=GR if o['i1'] > 0 else RD, scale=12)
    a.text(3 * off, -1.6, '%s orbit, %s\nF³(z) − z turns %+d time%s' % (o['kind'], linear_type(o['trace']), o['i1'], '' if abs(o['i1']) == 1 else 's'),
           ha='center', color=GR if o['i1'] > 0 else RD)
a.set_xlim(-1.6, 4.6); a.set_ylim(-2.1, 1.5); a.set_aspect('equal'); a.axis('off')
a.set_title('Period-3 pair, same cyclic order: arrows of F³(z) − z round a small circle')
save(fig, 'fig7-standard-map-pairs.png')

# ---- fig 8: index pair classifies linear stability ----
fig, ax = plt.subplots(figsize=(10, 3.6))
recs = []
for n in range(1, 14):
    for m in range(0, n):
        if gcd(m, n) != 1: continue
        pr = birkhoff_pair(K, m, n)
        if all(o['ok'] for o in pr): recs += [(o['trace'], o['i1'], o['i2']) for o in pr]
recs = np.array(recs)
y = {(1, 1): 0, (-1, -1): 1, (1, -1): 2}
yy = [y[(int(a_), int(b_))] for _, a_, b_ in recs]
ax.scatter(np.sign(recs[:, 0]) * np.log10(1 + np.abs(recs[:, 0])), yy, color=BL, s=18)
for xv in (-2, 2):
    ax.axvline(np.sign(xv) * np.log10(1 + abs(xv)), color='k', ls='--', lw=1)
ax.set_yticks([0, 1, 2]); ax.set_yticklabels(['(+1, +1)', '(−1, −1)', '(+1, −1)'])
ax.set_xlabel('trace of D(Fⁿ), signed log scale (dashed: trace = ±2)'); ax.set_ylabel('(index Fⁿ, index F²ⁿ)')
ax.text(0, .25, 'elliptic', ha='center'); ax.text(1.0, 1.25, 'hyperbolic'); ax.text(-1.25, 1.75, 'inverse hyperbolic', ha='center')
ax.set_title('100 standard-map orbits: two winding numbers, taken without derivatives, sort every orbit by stability')
save(fig, 'fig8-index-pair-chart.png')

# ---- fig 9: the Riemann door ----
J = json.load(open('jensen.json'))
fig, ax = plt.subplots(1, 2, figsize=(12.5, 4.2))
M = np.full((8, 31), np.nan)
for r_ in J['rows']:
    if r_['d'] >= 2: M[r_['d'] - 1, r_['n']] = r_['newton_degree']
im = ax[0].imshow(M[1:], aspect='auto', origin='lower', extent=[-.5, 30.5, 1.5, 8.5], cmap='RdBu')
plt.colorbar(im, ax=ax[0], label='circle degree of real Newton map')
ax[0].set_xlabel('shift n'); ax[0].set_ylabel('Jensen degree d')
ax[0].set_title('Jensen polynomials of Riemann\'s Ξ: every row has degree 2 − d\n(all roots real, as RH requires; 248 checked)')
sv = np.linspace(-1.5, 4, 440) + 1e-3
dg = [rp1_degree_alg(*newton_parts(np.poly([-2, -.3, .5, 1.2]) + np.array([0, 0, 0, 0, s]))[:2]) for s in sv]
ax[1].step(sv, dg, color=BL, lw=2)
for w in (-0.372, 0.419, 3.27):
    ax[1].axvline(w, color=RD, ls='--', lw=1)
ax[1].text(-.37, -1.85, '4 real', fontsize=9)
ax[1].text(-1.45, .2, '2 real', fontsize=9); ax[1].text(.5, .2, '2 real', fontsize=9); ax[1].text(3.3, 1.7, '0 real', fontsize=9)
ax[1].set_xlabel('shift s in q(y) + s, q = (y + 2)(y + .3)(y − .5)(y − 1.2)'); ax[1].set_ylabel('circle degree')
ax[1].set_title('Real roots are lost only at a double-root wall (red),\nwhere the degree jumps by 2')
save(fig, 'fig9-riemann-door.png')

# ---- fig 10: Shub floor, measured entropy, ceiling ----
fig, ax = plt.subplots(figsize=(12, 4.4))
rows = [('one element x∘a\n(Möbius)', 0, 0, 0), ('x∘x, Δ > 0', 0, 0, math.log(2)), ('x∘x, Δ < 0', math.log(2), math.log(2), math.log(2)),
        ('x∘x∘x, Δ < 0', math.log(3), math.log(3), math.log(3)), ('Newton, n = r = 5\n(Hurley–Martin)', math.log(3), math.log(3), math.log(5)),
        ('Newton, n = 6, r = 0\nno folds', math.log(2), math.log(2), math.log(6)),
        ('Newton, n = 6, r = 0\nwith folds (check C6)', math.log(2), .904, math.log(6)),
        ('Newton y³ − 2y + 2\n(check C5)', 0, .609, math.log(3))]
xpos = np.arange(len(rows))
ax.bar(xpos - .25, [r_[1] for r_ in rows], .25, color=OR, label='Shub / Lefschetz floor log|degree|')
ax.bar(xpos, [r_[2] for r_ in rows], .25, color=GR, label='real entropy (exact where a theorem gives it, else lap growth measured)')
ax.bar(xpos + .25, [r_[3] for r_ in rows], .25, color=GY, label='ceiling: log(algebraic degree) (Gromov–Lyubich)')
ax.set_xticks(xpos); ax.set_xticklabels([r_[0] for r_ in rows], fontsize=8.5); ax.set_ylabel('entropy')
ax.legend(fontsize=8.5); ax.set_title('Where the real entropy sits between the topological floor and the complex ceiling')
save(fig, 'fig10-shub-ladder.png')
print('figures written')
