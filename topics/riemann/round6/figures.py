"""Pictures for Riemann round 6. Colour key: blue = our reading, orange = already known, grey = reused or reference,
red = where a road stops."""
import json, math, numpy as np, matplotlib
matplotlib.use('Agg'); import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from lattice_engine import gauge, square_setup, hex_setup
B, O, G_, R = '#1f77b4', '#ff7f0e', '#888888', '#d62728'
plt.rcParams.update({'font.size': 10, 'axes.spines.top': False, 'axes.spines.right': False})
LZ = json.load(open('lattice_zeros.json')); EP = json.load(open('euler_pairs.json'))
LP = json.load(open('lattice_pairs.json')); RB = json.load(open('robustness.json'))

# ---------------- fig0: the pair map
D = ['k\nhomothety', 't\nodds', 's\nheight', 'd\ncurvature', 'q\norder', 'r\nratio', 'p\nroundness', 'f\nface power\n(p on faces)']
code = {}   # (i, j) with i > j : (label, colour)
reuse = {('d', 'k'): 'R1', ('t', 'k'): 'R1/R3', ('q', 't'): 'R4', ('s', 'k'): 'R5', ('r', 'q'): 'R5', ('p', 's'): 'R5'}
verdict = {('s', 't'): ('collapses', R), ('d', 't'): ('lines split\nper prime', R), ('q', 'k'): ('power\n1/(1−q)', B),
           ('r', 'k'): ('zeros kept,\nwall moves', B), ('d', 's'): ('= Hasse\nbound', O), ('p', 'k'): ('no Euler\nproduct', R),
           ('r', 't'): ('some zeros\nblow up', R), ('q', 's'): ('collapses', R), ('f', 'k'): ('Euler only\np = 2, 4, ∞', B),
           ('p', 't'): ('no Euler\nproduct', R), ('r', 's'): ('collapses', R), ('q', 'd'): ('field\ndimensions', O),
           ('f', 't'): ('Euler only\np = 2, 4, ∞', B), ('r', 'd'): ('lines split\nper prime', R), ('f', 's'): ('collapses', R),
           ('p', 'd'): ('zeros off,\nsome Re > 1', R), ('f', 'd'): ('zeros off,\nsome Re > 1', R), ('p', 'q'): ('no Euler\nproduct', R),
           ('f', 'q'): ('Euler only\np = 2, 4, ∞', B), ('p', 'r'): ('scale only:\nzeros fixed', R), ('f', 'r'): ('scale only:\nzeros fixed', R),
           ('f', 'p'): ('line only\nwhen round', B)}
keys = ['k', 't', 's', 'd', 'q', 'r', 'p', 'f']
fig, ax = plt.subplots(figsize=(11, 9.6)); ax.set_xlim(0, 8); ax.set_ylim(0, 8.6); ax.axis('off')
for i in range(8):
    for j in range(i):
        a, b = keys[i], keys[j]
        x, y = j, 7 - i
        if (a, b) in reuse:
            lab, col, fc = 'reused\n' + reuse[(a, b)], G_, '#eeeeee'
        else:
            lab, col = verdict[(a, b)]; fc = {R: '#fde8e8', B: '#e6effa', O: '#fff1e0', G_: '#f2f2f2'}[col]
        ax.add_patch(Rectangle((x + 0.04, y + 0.04), 0.92, 0.92, fc=fc, ec=col, lw=1.4))
        ax.text(x + 0.5, y + 0.5, lab, ha='center', va='center', fontsize=7.6, color='black')
for i in range(8):
    ax.text(i + 0.5, 7 - i + 0.5, D[i], ha='center', va='center', fontsize=9, weight='bold')
for col, lab, yy in [(G_, 'reused from rounds 1, 3, 4, 5', 8.35), (B, 'our reading: structure, nothing on RH', 8.1),
                     (O, 'restates something known', 7.85), (R, 'road stops (collapses, splits, zeros leave)', 7.6)]:
    ax.add_patch(Rectangle((4.3, yy - 0.08), 0.25, 0.18, fc='white', ec=col, lw=1.6)); ax.text(4.65, yy, lab, va='center', fontsize=9)
ax.set_title('All 28 pairs on Riemann: the seven dials plus face power (p on a polytope) as its own column', fontsize=11)
plt.savefig('figures/fig0-pair-map.png', dpi=150, bbox_inches='tight'); plt.close()

# ---------------- fig1: shapes of the two power dials
th = np.linspace(0, 2*np.pi, 2000); U = np.stack([np.cos(th), np.sin(th)], 1)
fig, axs = plt.subplots(1, 2, figsize=(11.5, 5.4))
for ax, (name, (Bm, N)), plist in [(axs[0], ('square: roundness p', square_setup()), [2, 3, 4, 8, np.inf]),
                                  (axs[1], ('hexagon: face power p', hex_setup()), [2, 3, 4, 6, np.inf])]:
    cols = plt.cm.viridis(np.linspace(0, 0.9, len(plist)))
    for p, c in zip(plist, cols):
        g = gauge(U, N, p); r = 1/g; r = r/r.max()
        lab = f"p = {'∞' if np.isinf(p) else int(p)}"
        if name.startswith('hexagon') and p in (2, 4): lab += '  (a circle)'
        ax.plot(r*U[:, 0], r*U[:, 1], color=c, lw=2 if p in (2, 4) else 1.3, label=lab)
    ax.set_aspect('equal'); ax.set_title(name); ax.legend(fontsize=8.5, loc='lower left', bbox_to_anchor=(0.62, 0.0))
    ax.set_xlim(-1.15, 1.9); ax.set_ylim(-1.15, 1.15); ax.set_xticks([]); ax.set_yticks([])
fig.suptitle('The unit shapes {g ≤ 1}, scaled to touch the same circle. The hexagon is exactly round at p = 2 and p = 4', fontsize=11)
plt.savefig('figures/fig1-power-dial-shapes.png', dpi=150, bbox_inches='tight'); plt.close()

# ---------------- fig2/3: zero trajectories
def traj(name, fname, title, mark4):
    data = LZ[name]; fig, ax = plt.subplots(figsize=(11, 5.6))
    for key, v in data.items():
        p = np.inf if key == 'inf' else float(key); xpos = 1 - 2/p
        z = np.array(v['zeros']); sc = ax.scatter(np.full(len(z), xpos), z[:, 0], c=z[:, 1], cmap='plasma', vmin=2, vmax=30, s=22, zorder=3)
    ax.axhline(0.5, color=B, lw=1.2, ls='--'); ax.text(1.035, 0.51, 'Re s = 1/2', color=B, va='bottom', fontsize=9)
    ax.axhline(0.75, color=O, lw=1.2, ls='--'); ax.text(1.035, 0.76, 'Re s = 3/4 (ζ(2s−1))', color=O, va='bottom', fontsize=9)
    ax.axhspan(1.0, 1.35, color='#fde8e8', zorder=0); ax.text(0.02, 1.33, 'Re s > 1: where the series converges absolutely', color=R, fontsize=9, va='top')
    tick_p = [2, 2.5, 3, 4, 6, 8, 20, np.inf]; ax.set_xticks([1 - 2/p for p in tick_p])
    ax.set_xticklabels(['2', '2.5', '3', '4', '6', '8', '20', '∞'])
    ax.set_xlabel('power p (axis spaced as 1 − 2/p)'); ax.set_ylabel('real part of each zero, 2 < t < 30')
    ax.set_ylim(0.0, 1.35); ax.set_xlim(-0.03, 1.25)
    if mark4:
        ax.annotate('p = 4: the face-power hexagon\nis round again, all zeros\nback on the line', (0.5, 0.5), (0.27, 0.95), fontsize=9,
                    arrowprops=dict(arrowstyle='->', color=B), color=B)
    ax.annotate('p = 2: ζ(s)L(s,χ), all on the line', (0, 0.5), (0.03, 0.04), fontsize=9, color=B, arrowprops=dict(arrowstyle='->', color=B))
    cb = plt.colorbar(sc, ax=ax, pad=0.09); cb.set_label('height t of the zero')
    ax.set_title(title, fontsize=11); plt.savefig(f'figures/{fname}', dpi=150, bbox_inches='tight'); plt.close()
traj('square', 'fig2-square-zeros.png', 'Roundness p on the square lattice: the zeros leave the line at once and settle on Re s = 3/4', False)
traj('hexagon', 'fig3-hexagon-zeros.png', 'Face power p on the hexagonal lattice: on the line at p = 2 and p = 4 only, then on to 3/4', True)

# ---------------- fig4: curvature (lattice shape) with the power dials
fig, axs = plt.subplots(1, 2, figsize=(12.5, 5.2), sharey=True)
for ax, name, ttl in [(axs[0], 'disc-20_square', 'Rectangle lattice Z × √5 Z, square normals\n(p = 2 is x² + 5y², class number 2)'),
                      (axs[1], 'disc-23_voronoi', 'Skew lattice of x² + xy + 6y², its own Voronoi hexagon\n(p = 2 is NOT the lattice norm: irregular normals)')]:
    for key, mk, col in [('2.0', 'o', B), ('3.0', 's', O), ('4.0', '^', G_)]:
        z = np.array(LZ[name][key]['zeros']); ax.scatter(z[:, 0], z[:, 1], marker=mk, color=col, s=30, label=f'p = {key[0]}', zorder=3)
    ax.axvline(0.5, color=B, ls='--', lw=1); ax.axvspan(1.0, 1.35, color='#fde8e8', zorder=0); ax.set_ylim(2, 31)
    ax.set_xlabel('Re s'); ax.set_title(ttl, fontsize=10); ax.set_xlim(0.1, 1.35); ax.legend(loc='lower right', fontsize=8.5)
ex20 = [r['genus_formula'] for r in RB['p2_curvature_exact_check']['disc-20_square'] if abs(r['genus_formula'][0] - 0.5) > 1e-6]
ex23 = [r['exact'] for r in RB['p2_curvature_exact_check']['disc-23_voronoi']['rows'] if abs(r['exact'][0] - 0.5) > 1e-6]
for z in ex20: axs[0].annotate(f'{z[0]:.4f} + {z[1]:.4f}i\n(checked on the exact\ngenus formula)', (z[0], z[1]), (0.12, 19.0), fontsize=8.3, arrowprops=dict(arrowstyle='->', color=B))
for z, ty in zip(ex23, [16.5, 26.0]): axs[1].annotate(f'{z[0]:.4f} + {z[1]:.4f}i', (z[0], z[1]), (0.985, ty), fontsize=8.3, arrowprops=dict(arrowstyle='->', color=B))
axs[0].set_ylabel('height t')
fig.suptitle('Curvature with the power dials: off the square and hexagon, even p = 2 has zeros off the line (Davenport–Heilbronn)', fontsize=11, y=1.05)
plt.savefig('figures/fig4-curvature.png', dpi=150, bbox_inches='tight'); plt.close()

# ---------------- fig5: Euler-side pairs
fig, axs = plt.subplots(2, 2, figsize=(13, 9.4))
ax = axs[0, 0]; P = np.array([p for p in range(2, 3000) if all(p % q for q in range(2, int(p**0.5) + 1)) and p != 37])
for t, col in [(0.8, O), (1.25, R)]: ax.plot(P, 0.5 + np.log(t)/np.log(P), color=col, lw=1.5, label=f'odds t = {t}')
for rho, col, ls in [(0.8, B, '--'), (1.25, G_, '--')]: ax.plot(P, 0.5*np.log(P)/np.log(P/rho), color=col, ls=ls, lw=1.5, label=f'ratio ρ = {rho}')
ax.axhline(0.5, color='k', lw=0.8); ax.set_xscale('log'); ax.set_xlabel('prime p'); ax.set_ylabel('line where the local factor vanishes')
ax.set_title('(a) t + d and d + r on curve 37a: one line per prime, no common line', fontsize=10); ax.legend(fontsize=8.3)
ax = axs[0, 1]
for key, col in [('k=0.6,q=0.3', B), ('k=1.3,q=0.6', O)]:
    rows = np.array(EP['kq_leading_power'][key]['rows']); x = 1/np.log(1/rows[:, 0]); y = rows[:, 1]
    fit = np.polyfit(x, y, 1); ax.plot(x, y, 'o', color=col, label=f'{key.replace(",", ", ")}: limit {fit[1]:.3f}, predicted {EP["kq_leading_power"][key]["predicted"]:.3f}')
    xx = np.linspace(0, x.max(), 50); ax.plot(xx, np.polyval(fit, xx), color=col, lw=1, ls=':')
ax.set_xlabel('1 / log(1/ε)'); ax.set_ylabel('log Z / log ζ at s = 1 + ε'); ax.set_xlim(left=0)
ax.set_title('(b) k + q: the leading power of ζ is 1/(1 − q), whatever k is', fontsize=10); ax.legend(fontsize=8)
ax = axs[1, 0]
hs = [0.8, 0.95, 1.0]
hc = EP['height_collapse_on_kr']; off_kr = [abs(hc['0.8']['off_line']), abs(hc['0.95']['off_line']), 0]
ax.plot(hs, off_kr, 'o-', color=B, label='ζ · H (homothety at k ρ^s): first zero')
for key, col in [('square_p2.0', O), ('hexagon_p4.0', G_)]:
    v = LP['sf_height_a_points'][key]; ax.plot(hs, [v[h]['max_off_line'] for h in ['0.8', '0.95', '1.0']], 's-', color=col, label=key.replace('_p', ' lattice, p = ').replace('.0', '') + ': first 6 zeros')
ax.set_xlabel('height s'); ax.set_ylabel('distance of the singular points from Re = 1/2')
ax.set_title('(c) s with t, q, r, f: the height sees the zeros only at s = 1', fontsize=10); ax.legend(fontsize=8)
ax = axs[1, 1]; rhos = np.linspace(0.3, 1.9, 200)
for k, col in [(1.2, B), (1.5, O), (2.0, R)]: ax.plot(rhos, np.log(k)/np.log(2/rhos), color=col, label=f'k = {k}')
ax.axhline(0.5, color='k', lw=0.8); ax.text(1.5, 0.42, 'critical line', fontsize=8.5)
ax.set_ylim(0, 2); ax.set_xlabel('ratio ρ'); ax.set_ylabel('wall σ_w = log k / log(2/ρ)')
ax.set_title('(d) k + r: zeros kept exactly while the wall stays left of 1/2', fontsize=10); ax.legend(fontsize=8.3, loc='upper left')
plt.tight_layout(); plt.savefig('figures/fig5-euler-pairs.png', dpi=150, bbox_inches='tight'); plt.close()

# ---------------- fig6: light and flux on face power
fig, axs = plt.subplots(1, 2, figsize=(12.5, 4.8))
ax = axs[0]; fl = LP['flux_discrepancy']
names = ['square_p2.0', 'square_p3.0', 'square_p4.0', 'square_pinf', 'hexagon_p3.0', 'hexagon_p4.0', 'hexagon_p6.0', 'hexagon_pinf']
pred = [0.5, 2/3, 3/4, 1.0, None, 0.5, None, 1.0]
labs = ['sq 2', 'sq 3', 'sq 4', 'sq ∞', 'hex 3', 'hex 4', 'hex 6', 'hex ∞']
ax.bar(range(8), [fl[n]['rms_exponent'] for n in names], color=[B]*8, alpha=0.8, label='measured, R = 40 to 900')
for i, pv in enumerate(pred):
    if pv is not None: ax.plot([i - 0.4, i + 0.4], [pv, pv], color=O, lw=2.2, label='known size (1/2, 1 − 1/p, 1)' if i == 0 else None)
ax.set_xticks(range(8)); ax.set_xticklabels(labs); ax.set_ylabel('growth exponent of N(R) − area·R²')
ax.set_title('(a) Flux tubes: boundary flux grows as the shape flattens', fontsize=10); ax.legend(fontsize=8.3, loc='upper left'); ax.set_ylim(0, 1.25)
ax = axs[1]; ax.axis('off')
lw = EP['light_wick_hexagon']
lines = ['(b) Light: Wick-rotate the face-power hexagon (y → i y)', '']
for p in ['2', '4', '6', '8']:
    lines.append(f"p = {p}: real form {'= power of x² − y²  → indefinite quadratic form (Hecke, Q(√3))' if lw[p]['equals_power_of_x2_minus_y2'] else 'not a norm form  → no Euler product'}")
lines += ['', 'In lattice terms a² + ac + c² becomes 2a² + 2ac − c²', '(discriminant 12). Only the round points survive light.']
for i, l in enumerate(lines): ax.text(0.0, 0.95 - 0.09*i, l, fontsize=9.6 if i else 10.5, transform=ax.transAxes, weight='bold' if i == 0 else None)
plt.tight_layout(); plt.savefig('figures/fig6-light-flux.png', dpi=150, bbox_inches='tight'); plt.close()

# ---------------- fig7: order dial on quadratic fields
dq = EP['dq_prime_fractal_dimension']; fig, ax = plt.subplots(figsize=(8, 4.4))
ks = list(dq); ax.bar(ks, [dq[k] for k in ks], color=[G_, B, B, B])
for i, k in enumerate(ks): ax.text(i, dq[k] + 0.02, f'{dq[k]:.4f}', ha='center', fontsize=9.5)
ax.set_ylim(1, 1.65); ax.set_ylabel('D with Σ over prime ideals of N(P)^−D = 1')
ax.set_title('d + q: where the order dial’s wall sits for each quadratic field (Q is round 4’s 1.3994)', fontsize=10)
plt.savefig('figures/fig7-field-dimensions.png', dpi=150, bbox_inches='tight'); plt.close()
print('ok')
