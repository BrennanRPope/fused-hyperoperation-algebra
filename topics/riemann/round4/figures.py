"""Pictures for Riemann round 4. Colour key: blue = our reading, orange = already known, grey = reference,
red = where a road stops."""
import json, numpy as np, matplotlib
matplotlib.use('Agg'); import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch
B, O, G, R = '#1f77b4', '#ff7f0e', '#888888', '#d62728'
plt.rcParams.update({'font.size': 10, 'axes.spines.top': False, 'axes.spines.right': False})
pm = json.load(open('polymath.json')); ex = json.load(open('extrema.json')); rw = json.load(open('rewind4.json'))
br = json.load(open('barrier.json')); bd = json.load(open('braid_deep.json')); r3 = json.load(open('../round3/rewind.json'))

# ---------------- fig0: map
fig, ax = plt.subplots(figsize=(12, 6.4)); ax.axis('off'); ax.set_xlim(0, 12); ax.set_ylim(0, 6.4)
ax.add_patch(FancyBboxPatch((4.45, 2.7), 3.1, 1.0, boxstyle='round,pad=0.1', fc='#eef', ec=B))
ax.text(6.0, 3.2, 'What has to stay positive\nwhile the heat flow is rewound\nfrom t = 0.22 to t = 0?', ha='center', va='center', fontsize=10.5)
items = [((0.15, 4.8), 'Polymath 15: a dominance margin\n(first Dirichlet term beats the rest)', 'the forward flow builds it;\nit is gone at t = 0', O),
         ((0.15, 2.65), 'Ours: the local fused discriminant\nD = −8Ξ/Ξ″ at every extremum\n(= Laguerre’s expression, pair by pair)', 'read straight off Ξ;\nall 1,318 positive', B),
         ((0.15, 0.5), 'Its exact rewind rate\ndD/dt = 8 + D·(log Ξ″)″', 'one point of Ξ predicts\nLehmer’s merger to 1e-6', B),
         ((8.65, 4.8), 'The Laguerre ladder:\neach rung’s sink is the next rung', 'never closes: no finite\nset of rungs controls itself', R),
         ((8.65, 2.65), 'Barrier chain as composed letters\nz ↦ p^(−s)·z + c, with dials\n(shift c, order q, link k, d)', 'order counts once centres\nsplit; no dial tames the line', R),
         ((8.65, 0.5), 'Topological word of the rewind\n(three-body style syzygies)', 'empty: 70 mergers,\nno zero changes its slot', G)]
for (x, y), txt, verdict, col in items:
    ax.add_patch(FancyBboxPatch((x, y), 3.2, 1.4, boxstyle='round,pad=0.08', fc='white', ec=col, lw=1.8))
    ax.text(x + 1.6, y + 0.93, txt, ha='center', va='center', fontsize=8.4)
    ax.text(x + 1.6, y + 0.27, verdict, ha='center', va='center', fontsize=8.3, color=col, weight='bold')
    ax.annotate('', xy=(4.45 if x < 4 else 7.55, 3.2), xytext=(x + 3.2 if x < 4 else x, y + 0.7),
                arrowprops=dict(arrowstyle='->', color=col, lw=1.2))
ax.set_title('Round 4: the quantity that has to stay positive along the rewind, and the barrier chain deformed by our dials', fontsize=11)
plt.savefig('figures/fig0-map.png', dpi=150, bbox_inches='tight'); plt.close()

# ---------------- fig1: Polymath dominance margin
fig, axs = plt.subplots(1, 2, figsize=(13, 4.6))
ts = np.array(pm['ts']); ax = axs[0]
cols = plt.cm.Greys(np.linspace(0.35, 0.95, len(pm['xs'])))
for c, x in zip(cols, pm['xs']):
    ax.plot(ts, pm['m'][str(x)], color=c, lw=1.6, label=f'x = {x:.0e}')
ax.axhline(0, color=R, lw=1); ax.set_ylim(-6, 1.05)
ax.set_xlabel('heat-flow time t'); ax.set_ylabel('toy margin  1 − Σ (terms n ≥ 2)')
ax.set_title('Bare first-term dominance (no mollifier): positive only for t ≳ 0.25 at x = 6e10', fontsize=9.5)
ax.legend(fontsize=7.5, loc='lower right')
ax = axs[1]
mt = np.array(pm['moll_ts'])
ax.plot(ts, pm['m']['60000000000.0'], color=G, lw=1.4, label='bare, x = 6e10')
ax.plot(mt, pm['moll']['60000000000.0'], color=O, lw=2.2, marker='o', ms=3, label='Euler mollifier p = 2, 3, 5, 7, x = 6e10')
ax.axhline(0, color=R, lw=1); ax.axvline(0.2, color=O, ls=':', lw=1)
ax.text(0.205, -3.5, 'Polymath 15\nt0 = 0.2', color=O, fontsize=8.5)
ax.text(0.27, -4.6, 'at t = 0 the bare margin is\n−212 at x = 6e10:\nthe rewind takes it away', fontsize=8.5, color=R)
ax.set_ylim(-6, 1.05); ax.set_xlabel('heat-flow time t'); ax.legend(fontsize=8, loc='lower right')
ax.set_title('With the mollifier the toy turns positive at t = 0.2, matching Polymath’s choice', fontsize=9.5)
plt.savefig('figures/fig1-polymath-margin.png', dpi=150, bbox_inches='tight'); plt.close()

# ---------------- fig2: local fused discriminant from Xi vs gap^2
fig, axs = plt.subplots(1, 2, figsize=(13, 4.6))
for ax, key, lab in [(axs[0], 'main', 'zeros 1 to 1,000'), (axs[1], 'lehmer', 'zeros 6551 to 6870 (Lehmer window)')]:
    S = ex[key]; ge = np.array([r['ge'] for r in S]); De = np.array([r['De'] for r in S]); g2 = np.array([r['gap'] for r in S])**2
    ax.scatter(ge, De/g2, s=6, color=B, label='D / gap², D = −8Ξ/Ξ″ at the extremum')
    ax.axhline(1, color=G, ls='--', lw=1); ax.text(ge[0], 1.02, 'lone pair: D = gap² exactly', color=G, fontsize=8)
    ax.axhline(0, color=R, lw=1); ax.text(ge[0], 0.03, 'D < 0 would be a complex pair here', color=R, fontsize=8)
    ax.set_ylim(-0.08, 1.1); ax.set_xlabel('height γ of the extremum'); ax.set_ylabel('D / gap²')
    i = int(np.argmin(De)); ax.scatter([ge[i]], [De[i]/g2[i]], s=40, facecolor='none', edgecolor=O)
    ax.set_title(f'{lab}: all {len(S)} positive; smallest D = {De[i]:.5f} at γ = {ge[i]:.2f}', fontsize=9.5)
    ax.legend(fontsize=8, loc='lower right')
plt.savefig('figures/fig2-local-disc.png', dpi=150, bbox_inches='tight'); plt.close()

# ---------------- fig3: through the rewind
fig, axs = plt.subplots(1, 2, figsize=(13, 4.6))
ty = rw['toy']; tt = np.array([r['t'] for r in ty['traj']]); DD = np.array([r['De'] for r in ty['traj']])
ax = axs[0]
ax.plot(tt, DD, color=B, lw=2.2, label='D(t), exact heat flow of the toy')
ax.plot(tt, 0.09 + 8*tt, color=G, ls='--', label='lone pair: gap² + 8t')
k = (ty['rate0']-8)/ty['De0']; pred = (ty['De0'] + 8/k)*np.exp(k*tt) - 8/k
ax.plot(tt, pred, color=O, ls=':', lw=2, label='one-point prediction from D, rate at t = 0')
ax.axhline(0, color=R, lw=1); ax.axvline(ty['t_merge'], color=R, ls=':', lw=1)
ax.text(ty['t_merge'], 0.06, f" merger t = {ty['t_merge']:.5f}", color=R, fontsize=8.5)
ax.invert_xaxis(); ax.set_xlabel('t (rewind runs to the right)'); ax.set_ylabel('local fused discriminant D')
ax.set_title('Toy roots ±{1, 2, 2.3, 3.5, 5}: D crosses 0 exactly at the merger and keeps going', fontsize=9.5)
ax.legend(fontsize=8)
ax = axs[1]
names, full, onept, lone = [], [], [], []
Lh = ex['lehmer']; il = min(range(len(Lh)), key=lambda i: Lh[i]['gap']); e = Lh[il]
kk = e['kappa']; onept.append((4/kk)*np.log(8/(8+kk*e['De']))); full.append(r3['lehmer']['t_full']); lone.append(-e['gap']**2/2); names.append('Lehmer 6709')
for f_ in r3['full']:
    if f_['j']-1 < len(ex['main']):
        e = ex['main'][f_['j']-1]; kk = e['kappa']
        onept.append((4/kk)*np.log(8/(8+kk*e['De']))); full.append(f_['t_full']); lone.append(-f_['tau_isolated']); names.append(f"zero {f_['j']}")
full, onept, lone = map(np.array, (full, onept, lone))
xx = np.arange(len(names)); w = 0.35
ax.bar(xx-w/2, np.abs(lone/full-1), w, color=G, label='lone pair, gap²/2')
ax.bar(xx+w/2, np.abs(onept/full-1), w, color=B, label='one point of Ξ: D and (log Ξ″)″')
ax.set_yscale('log'); ax.set_xticks(xx); ax.set_xticklabels(names)
ax.set_ylabel('relative error vs the full rewind (round 3)')
ax.set_title('Predicting the first merger time without moving any zero', fontsize=9.5); ax.legend(fontsize=8)
plt.savefig('figures/fig3-rewind-disc.png', dpi=150, bbox_inches='tight'); plt.close()

# ---------------- fig4: rate field D|kappa| vs 8
fig, ax = plt.subplots(figsize=(7.5, 5))
for key, col, lab in [('main', B, 'zeros 1 to 1,000'), ('lehmer', O, 'Lehmer window')]:
    S = ex[key]; De = np.array([r['De'] for r in S]); kap = np.abs([r['kappa'] for r in S])
    ax.scatter(kap, De, s=6, color=col, alpha=0.7, label=lab)
kx = np.logspace(-1, 3, 100); ax.plot(kx, 8/kx, color=R, lw=1.5, label='D·|(log Ξ″)″| = 8: steady (like sin)')
ax.set_xscale('log'); ax.set_yscale('log'); ax.set_xlabel('|(log Ξ″)″| at the extremum'); ax.set_ylabel('D')
ax.text(0.15, 30, 'above the line: the gap opens\nunder the rewind', fontsize=8.5, color=G)
ax.text(30, 0.002, 'below: the gap closes\n(Lehmer pair is far below)', fontsize=8.5, color=G)
ax.set_title('Rewind rate of D: −8 + D·|(log Ξ″)″|; it needs the next rung (Ξ″) to be log-concave', fontsize=9.5)
ax.legend(fontsize=8, loc='lower left')
plt.savefig('figures/fig4-rate-field.png', dpi=150, bbox_inches='tight'); plt.close()

# ---------------- fig5: letters and the order dial
fig, axs = plt.subplots(1, 2, figsize=(13, 4.6))
ax = axs[0]
s = 0.5 + 14.134725j; u2, u3 = 2**(-s), 3**(-s)
z23 = u2*(u3*0+1)+1; z32 = u3*(u2*0+1)+1
c2, c3 = 1/(1-u2), 1/(1-u3)
ax.scatter([c2.real, c3.real], [c2.imag, c3.imag], color=O, s=50, zorder=3)
ax.text(c2.real-0.02, c2.imag+0.03, 'centre of f₂ =\nEuler factor at 2', color=O, fontsize=8.5, ha='right')
ax.text(c3.real-0.02, c3.imag-0.07, 'centre of f₃ =\nEuler factor at 3', color=O, fontsize=8.5, ha='right')
ax.scatter([z23.real, z32.real], [z23.imag, z32.imag], color=B, s=50, zorder=3)
ax.annotate('', xy=(z32.real, z32.imag), xytext=(z23.real, z23.imag), arrowprops=dict(arrowstyle='<->', color=R))
ax.text(z23.real+0.02, z23.imag+0.05, 'f₂∘f₃(0)', color=B); ax.text(z32.real+0.02, z32.imag-0.08, 'f₃∘f₂(0)', color=B)
ax.text(0.04, 0.17, 'order gap = 2^(−s) − 3^(−s)\n(a wedge; zero only on Re s = 0)', color=R, fontsize=8.5)
ax.scatter([0], [0], color=G, s=30); ax.text(0.02, -0.08, 'c = 0: one shared centre, the letters commute', color=G, fontsize=8)
ax.set_aspect('equal', 'datalim'); ax.set_xlabel('Re'); ax.set_ylabel('Im')
ax.set_title('Letters z ↦ p^(−s)·z + 1 at s = 1/2 + 14.1347i: order counts', fontsize=9.5)
ax = axs[1]
qs = [r[0] for r in br['Zq_two_prime_check']]; zq = [r[2] for r in br['Zq_two_prime_check']]
qq = np.linspace(0, 1, 41); x_, y_ = 0.25, 1/9
def closed(q):
    tot = 0; 
    for a in range(200):
        den = np.prod([1-q**j*y_ for j in range(a+1)]); tot += x_**a/den
    return tot
ax.plot(qq, [closed(q) for q in qq], color=B, lw=2)
ax.scatter([0, 1], [br['euler_23'], br['free_23']], color=O, zorder=3)
ax.text(0.04, br['euler_23']-0.001, 'q = 0: Euler product (bag)', color=O, fontsize=8.5)
ax.text(0.55, br['free_23']-0.006, 'q = 1: free words\n1/(1 − 2^(−s) − 3^(−s))', color=O, fontsize=8.5)
ax.set_xlabel('order dial q (each out-of-order pair weighs q)'); ax.set_ylabel('Z_q(s = 2), primes 2 and 3')
ax.set_title(f"Order dial: with every prime the wall moves from σ = 1 (q < 1) to σ = {br['sigma0_all_primes']:.3f} (q = 1)", fontsize=9.5)
plt.savefig('figures/fig5-letters-order-dial.png', dpi=150, bbox_inches='tight'); plt.close()

# ---------------- fig6: link dial on the line
fig, ax = plt.subplots(figsize=(8, 4.8))
ld = br['link_dial']; lx = ld['loglog']
cmap = {'zeta (k=-1, d=0)': O, 'elliptic, fixed angle pi/3 (k=-1, d=1)': B, 'elliptic, Sato-Tate angles (d=1)': B,
        'k = 0, d = 1 (alpha = i, beta = -i)': G, 'hyperbolic (k=-1.2, d=0.2)': R}
ls = {'elliptic, Sato-Tate angles (d=1)': '--'}
for kname, v in ld['var'].items():
    ax.plot(lx, v, color=cmap[kname], ls=ls.get(kname, '-'), marker='o', ms=3, label=f"{kname}; mean k² = {ld['k2_mean'][kname]:.2f}")
ax.set_xlabel('log log X  (primes up to X = 1e2 … 1e6)'); ax.set_ylabel('spread of log|L| on σ = 1/2 (chain variance)')
ax.set_title('Link dial (k, d): the spread grows like ½·mean(k²)·log log X; only k = 0 stays bounded', fontsize=9.5)
ax.legend(fontsize=7.5)
plt.savefig('figures/fig6-link-dial.png', dpi=150, bbox_inches='tight'); plt.close()

# ---------------- fig7: topological word of the rewind
fig, ax = plt.subplots(figsize=(11, 4.8))
cen = np.array(bd['centres'])/2; hh = np.array(bd['heights'])/2
real = hh == 0
ax.scatter(cen[real], np.zeros(real.sum()), s=8, color=G, label='still real at t = −0.4 (γ units)')
ax.scatter(cen[~real], hh[~real], s=14, color=B, label='complex pair (upper member)')
ax.scatter(cen[~real], -hh[~real], s=14, color=B, alpha=0.4)
GL = [float(x) for x in json.load(open('../round3/zeros_lehmer.json'))['gamma']]
ma = [m['at'] for m in bd['mergers']]
ax.scatter(ma, [-0.78]*len(ma), marker='^', s=14, color=O, label='where an inner merger happened')
ax.axvspan(GL[0], GL[40], color=G, alpha=0.12); ax.axvspan(GL[-41], GL[-1], color=G, alpha=0.12)
ax.text(GL[2], -0.7, 'edge (frozen\nbackground)', fontsize=7.5, color=G)
ax.set_xlabel('real part (γ units)'); ax.set_ylabel('imaginary part (γ units)')
ax.set_title(f"Rewind of the Lehmer window to t = −0.4: {len(bd['mergers'])} mergers inside, {len(bd['events'])} syzygies; every pair stays in its slot", fontsize=9.5)
ax.legend(fontsize=8, loc='upper right')
plt.savefig('figures/fig7-topological-word.png', dpi=150, bbox_inches='tight'); plt.close()
print('ok')

# ---------------- fig8: hypernatural log as a Laguerre dial
hl = json.load(open('hyperlog.json')); pi_ = np.array(hl['per_interval'])
fig, ax = plt.subplots(figsize=(11, 4.6))
ax.scatter(pi_[:, 0], pi_[:, 2], s=4, color=B, label='largest h* = ΞΞ″/Ξ′² in each gap (every 7th gap shown)')
ax.scatter([5], [hl['hstar_below_first_zero']], s=30, color=B, marker='D')
ax.text(120, hl['hstar_below_first_zero'] - 0.035, f"below the first zero: {hl['hstar_below_first_zero']:.3f}", fontsize=8.5, color=B)
for sv, n in [(0.5, 2), (2/3, 3), (0.75, 4)]:
    ax.axhline(sv, color=G, ls='--', lw=0.9); ax.text(4500, sv + 0.01, f's = {sv:.2f}: degree {n}', fontsize=8, color=G, va='bottom')
ax.axhline(1, color=R, lw=1.2); ax.text(4300, 1.01, 's = 1: Laguerre (RH here)', fontsize=8, color=R, va='bottom')
ax.set_ylim(-0.1, 1.1); ax.set_xlabel('height γ'); ax.set_ylabel('critical height h* of the hypernatural log')
ax.set_title('Hypernatural log height s as a Laguerre dial: Ξ passes the height-3/4 (degree-4) test on 0 < γ < 5355', fontsize=9.5)
ax.legend(fontsize=8, loc='upper center')
plt.savefig('figures/fig8-hyperlog-dial.png', dpi=150, bbox_inches='tight'); plt.close()
print('fig8 ok')
