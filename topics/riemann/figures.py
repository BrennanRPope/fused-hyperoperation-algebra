"""Pictures for the Riemann DRAFT. Colour key: blue = our reading, orange = already known, grey = reference."""
import json, math, numpy as np, matplotlib
matplotlib.use('Agg'); import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch
B, O, G, R = '#1f77b4', '#ff7f0e', '#888888', '#d62728'
plt.rcParams.update({'font.size': 10, 'axes.spines.top': False, 'axes.spines.right': False})
cube = json.load(open('cube.json')); weil = json.load(open('weil.json')); up = json.load(open('up_step.json'))
dbn = json.load(open('dbn.json')); lad = json.load(open('ladder.json'))
gam = np.load('../primes/zeta_zeros_1000.npy')

# ---- fig0: the map
fig, ax = plt.subplots(figsize=(11, 5.2)); ax.axis('off'); ax.set_xlim(0, 11); ax.set_ylim(0, 5.45)
ax.add_patch(FancyBboxPatch((4.1, 2.2), 2.8, 0.9, boxstyle='round,pad=0.1', fc='#eef', ec=B))
ax.text(5.5, 2.65, 'Riemann zeros\n(on the line Re s = 1/2?)', ha='center', va='center', fontsize=11)
items = [((0.2, 4.1), '1. Prime cube fold\nk dial changes the order of a zero,\nnever its place', 'not a lead', G),
         ((0.2, 2.3), '2. Corner signs (k = -1)\nRH = the signs balance\nlike coin flips', 'picture only', G),
         ((0.2, 0.4), '3. Local factor = dial point\nRH over finite fields =\nevery point elliptic', 'known (Weil); our picture', O),
         ((7.7, 4.1), '4. Lee–Yang ladder, one step up\nfrom Riemann: every zero\nleaves the line', 'route closed (proved)', B),
         ((7.7, 2.3), '5. Our p = 4 step, repeated\n= de Bruijn–Newman heat flow\nzeros meet and leave on rewind', 'known; same wall', O),
         ((7.7, 0.4), '6. Kernel ladder exp(-|u|^p)\nkink push decides; Riemann\'s\nkernel has no kink', 'consistent, proves nothing', G)]
for (x, y), txt, verdict, col in items:
    ax.add_patch(FancyBboxPatch((x, y), 3.0, 1.1, boxstyle='round,pad=0.08', fc='white', ec=col, lw=1.8))
    ax.text(x+1.5, y+0.66, txt, ha='center', va='center', fontsize=8.6)
    ax.text(x+1.5, y+0.1, verdict, ha='center', va='center', fontsize=8.6, color=col, weight='bold')
    ax.annotate('', xy=(4.1 if x < 4 else 6.9, 2.65), xytext=(x+3.0 if x < 4 else x, y+0.55), arrowprops=dict(arrowstyle='->', color=col, lw=1.2))
ax.set_title('Six roads from Affine Composition Theory toward the Riemann zeros, with where each one ends', fontsize=11)
plt.savefig('figures/fig0-map.png', dpi=150, bbox_inches='tight'); plt.close()

# ---- fig1: prime cube and local-zero lines
fig, axs = plt.subplots(1, 2, figsize=(12, 4.8))
ax = axs[0]; ax.axis('off'); ax.set_aspect('equal')
P = [2, 3, 5]; pos = {}
for m in range(8):
    S = [P[i] for i in range(3) if m >> i & 1]
    x = (m & 1)*2 + (m >> 2 & 1)*0.9; y = (m >> 1 & 1)*2 + (m >> 2 & 1)*0.9
    pos[m] = (x, y); n = int(np.prod(S)) if S else 1
    ax.plot(x, y, 'o', ms=26, mfc='white', mec=B, mew=2)
    ax.text(x, y, str(n), ha='center', va='center', fontsize=11, weight='bold')
    ax.text(x+0.28, y-0.38, '$k^{%d}$' % len(S), fontsize=9, color=B)
for a in range(8):
    for i in range(3):
        b = a | (1 << i)
        if b != a: ax.plot(*zip(pos[a], pos[b]), color=G, lw=1, zorder=0)
ax.set_xlim(-0.6, 3.6); ax.set_ylim(-0.7, 3.4)
ax.set_title('The prime cube on {2, 3, 5}: corners are squarefree n,\ncorner weight $k^{\\omega(n)} n^{-s}$; folding $\\circ_k$ sums them all', fontsize=10)
ax = axs[1]
ax.axvspan(0, 1, color='#f4f4f4'); ax.axvline(0.5, color=O, lw=2, label='zeta zeros, Re s = 1/2 (order k)')
for k, c in ((0.5, '#9ecae1'), (1, B), (2, R), (3, '#8c2d04')):
    for p, sig in cube['local_zero_lines'][str(k)]:
        if -0.6 <= sig <= 1.7:
            for m in range(-25, 26):
                ax.plot(sig, (2*m+1)*math.pi/math.log(p), 's', ms=3.5, color=c)
    ax.plot([], [], 's', color=c, label='k = %g: cube factors vanish at Re s = log k / log p' % k)
g6 = gam[gam < 40]; ax.plot(np.full(2*len(g6), 0.5), np.concatenate([g6, -g6]), 'o', color=O, ms=4)
ax.set_xlim(-0.6, 1.7); ax.set_ylim(-40, 40); ax.set_xlabel('Re s'); ax.set_ylabel('Im s')
ax.text(0.5, 37, 'critical strip', ha='center', color=G)
ax.legend(fontsize=7.5, loc='lower left'); ax.set_title('Where the dial puts zeros: for |k| <= 1 no cube factor vanishes in Re s > 0;\nk > 1 scatters zeros off the line (k = 2: p = 3 at Re s = 0.63)', fontsize=10)
plt.savefig('figures/fig1-prime-cube.png', dpi=150, bbox_inches='tight'); plt.close()

# ---- fig2: Mertens walk
fig, ax = plt.subplots(figsize=(10, 4.2))
mc = np.array(cube['mertens_curve']); cc = np.array(cube['coin_curve']); x = mc[:, 0]
ax.fill_between(x, -np.sqrt(x), np.sqrt(x), color='#eee', label='$\\pm\\sqrt{x}$ (coin-flip size)')
ax.plot(cc[:, 0], cc[:, 1], color=G, lw=0.8, label='random signs on the same corners')
ax.plot(x, mc[:, 1], color=B, lw=1.2, label='Möbius signs, M(x) = sum of $\\mu(n)$ (k = -1 corners)')
ax.set_xscale('log'); ax.set_xlabel('x'); ax.set_ylabel('running total of corner signs')
ax.set_title('RH says the prime cube\'s corner signs cancel like coin flips: $M(x) = O(x^{1/2+\\epsilon})$ (checked to $10^7$)', fontsize=10)
ax.legend(fontsize=8, loc='upper left'); plt.savefig('figures/fig2-corner-sign-walk.png', dpi=150, bbox_inches='tight'); plt.close()

# ---- fig3: dial plane
fig, axs = plt.subplots(1, 3, figsize=(15, 4.6))
ax = axs[0]; kk = np.linspace(-110, 110, 400)
ax.fill_between(kk, kk**2/4, 3200, color='#e8f1fb'); ax.plot(kk, kk**2/4, color=G, lw=1.5, label='$\\Delta = k^2 - 4d = 0$')
g1 = weil['genus1']; ax.plot([-r['a'] for r in g1], [r['p'] for r in g1], '.', color=B, ms=3, label='genus 1: $(-a_p, p)$, %d primes' % len(g1))
g2 = weil['genus2']; ax.plot([-r['A'] for r in g2] + [-r['B'] for r in g2], [r['p'] for r in g2]*2, '.', color=O, ms=4, label='genus 2: both factors, %d primes' % len(g2))
ax.set_xlabel('k'); ax.set_ylabel('d'); ax.set_ylim(0, 3100); ax.legend(fontsize=7.5, loc='upper right')
ax.text(0, 2200, 'elliptic side\n(positive curvature)', ha='center', color=B); ax.text(-95, 300, 'hyperbolic', color=G)
ax.set_title('Every local factor $1 + ku + du^2$ lands\non the elliptic side of our dial plane', fontsize=10)
ax = axs[1]; th = np.linspace(0, 2*np.pi, 300); p = 101
ax.plot(np.cos(th)/math.sqrt(p), np.sin(th)/math.sqrt(p), color=O, lw=1.5, label='$|u| = p^{-1/2}$  (Re s = 1/2)')
for r in g1:
    if r['p'] in (101, 103, 107, 109, 113):
        roots = np.roots([r['p'], r['k'], 1])*math.sqrt(r['p']/101); ax.plot(roots.real, roots.imag, 'o', color=B, ms=5)
ax.plot([], [], 'o', color=B, label='roots for p = 101..113, rescaled to p = 101')
ax.set_aspect('equal'); ax.axhline(0, color='#ddd'); ax.axvline(0, color='#ddd'); ax.legend(fontsize=7.5, loc='lower left')
ax.set_title('Elliptic = both roots on one circle:\nthe local Frobenius is a rotation', fontsize=10)
ax = axs[2]; thv = np.array([r['theta'] for r in g1])
ax.hist(thv, bins=16, range=(0, np.pi), density=True, color=B, alpha=0.6, label='rotation angles, genus 1')
t = np.linspace(0, np.pi, 200); ax.plot(t, 2/np.pi*np.sin(t)**2, color=O, lw=2, label='Sato–Tate $(2/\\pi)\\sin^2\\theta$')
ax.set_xlabel('rotation angle $\\theta_p$'); ax.legend(fontsize=7.5); ax.set_title('How the rotation angles spread', fontsize=10)
plt.savefig('figures/fig3-dial-plane-weil.png', dpi=150, bbox_inches='tight'); plt.close()

# ---- fig4: one step up
fig, axs = plt.subplots(1, 2, figsize=(13, 5))
ax = axs[0]
d1 = up['1']; xs = np.array(d1['grid_x']); ys = np.array(d1['grid_y']); L = np.array(d1['logabs'])
# grid_up stored first-quadrant rows; mirror to show both half-planes
im = ax.pcolormesh(xs, np.concatenate([-ys[::-1], ys]), np.vstack([L.T[::-1], L.T]), cmap='Greys_r', shading='auto', vmin=-25, vmax=5)
for z in d1['zeros']: ax.plot([z[0]]*2, [z[1], -z[1]], 'o', color=B, ms=6)
for z in up['2']['zeros']: ax.plot([z[0]]*2, [z[1], -z[1]], 'd', color='#6baed6', ms=5)
ax.plot(gam[gam < 60], 0*gam[gam < 60], 'o', color=O, ms=6)
ax.plot([], [], 'o', color=O, label='zeros of Riemann $\\Xi$ (all on the line, as far as checked)')
ax.plot([], [], 'o', color=B, label='one step up, $F_1(z) = \\frac{1}{z}\\int_0^z \\Xi$')
ax.plot([], [], 'd', color='#6baed6', label='$F_2$ ($\\kappa = 2$)')
ax.set_xlabel('Re z'); ax.set_ylabel('Im z'); ax.legend(fontsize=7.5, loc='lower left', framealpha=0.9)
ax.set_title('Step identity run upward from $\\Xi$: every zero leaves the line\n(shading: $\\log|F_1|$)', fontsize=10)
ax = axs[1]
import mpmath as mp
mp.mp.dps = 20
xi = lambda s: s*(s-1)/2*mp.pi**(-s/2)*mp.gamma(s/2)*mp.zeta(s)
tt = np.linspace(0.5, 60, 300); vals = [abs(float(mp.re(xi(0.5+1j*t)))) for t in tt]
ax.semilogy(tt, vals, color=O, lw=1, label='the wave: $|\\Xi(x)|$, dies like $e^{-\\pi x/4}$')
ax.axhline(2.8067, color=B, lw=2, label='the push: $M_1 = \\int_0^\\infty \\Xi = \\pi\\Phi(0) = 2.807$')
ax.set_xlabel('x (on the line)'); ax.legend(fontsize=8); ax.set_ylim(1e-25, 10)
ax.set_title('Why: the same push-against-wave as our Lee–Yang break.\nThe push is a constant, the wave is exponentially small', fontsize=10)
plt.savefig('figures/fig4-one-step-up.png', dpi=150, bbox_inches='tight'); plt.close()

# ---- fig5: heat flow zero paths
fig, ax = plt.subplots(figsize=(10, 5.2))
for t, v in dbn.items():
    ax.plot(v['changes'], [float(t)]*len(v['changes']), '.', color=B, ms=3)
ax.axhspan(0, 0.2, color='#fde0c5', label='$0 \\leq \\Lambda \\leq 0.2$: known range of the constant')
ax.axhline(0, color=O, lw=1.5, label='t = 0: Riemann\'s own zeros ($x = 2\\gamma$)')
ax.set_xlabel('x (zeros of $H_t$, real ones only)'); ax.set_ylabel('t  (heat time; our p = 4 step repeated)')
ax.set_xlim(0, 150); ax.legend(fontsize=8, loc='lower right')
ax.set_title('de Bruijn–Newman heat flow: on rewind (backward heat flow), zeros meet in pairs and leave the line.\nRH = no meeting has happened yet at t = 0; running forward never makes a meeting', fontsize=10)
plt.savefig('figures/fig5-heat-flow-zero-paths.png', dpi=150, bbox_inches='tight'); plt.close()

# ---- fig6: kernel ladder
fig, axs = plt.subplots(1, 2, figsize=(13, 4.6))
ax = axs[0]; u = np.linspace(-2.5, 2.5, 600)
for p, c in ((2.5, '#9ecae1'), (3, B), (4, O), (6, '#e6550d')):
    ax.plot(u, np.exp(-np.abs(u)**p), color=c, label='$e^{-|u|^{%g}}$' % p)
PhiT = lambda uu: 2*sum((2*np.pi**2*n**4*np.exp(4.5*uu) - 3*np.pi*n**2*np.exp(2.5*uu))*np.exp(-np.pi*n**2*np.exp(2*uu)) for n in range(1, 6))
ax.plot(u, PhiT(u)/PhiT(0), 'k--', lw=1.5, label='Riemann\'s $\\Phi$ (scaled): smooth centre,\ndouble-exponential fall')
ax.legend(fontsize=8); ax.set_xlabel('u'); ax.set_title('Kernels on the roundness dial, and Riemann\'s', fontsize=10)
ax = axs[1]
for p, c in (('3', B), ('5', '#9ecae1'), ('4', O), ('6', '#e6550d')):
    x = np.array(lad[p]['x']); v = np.array(lad[p]['val'])
    ax.plot(x[1:], np.sign(v[1:])*np.log10(np.abs(v[1:]) + 1e-60), color=c, label='p = %s' % p)
    if lad[p]['push'][1] is not None and p in ('3', '5'):
        pu = np.array(lad[p]['push'][1:], dtype=float); ax.plot(x[1:], np.sign(pu)*np.log10(np.abs(pu)), ':', color=c)
ax.set_ylim(-40, 3); ax.set_xlabel('x'); ax.set_ylabel('sign · log10 |transform|')
ax.legend(fontsize=8); ax.set_title('Non-even p: the kink push (dotted) takes over and the zeros stop.\nEven p: no kink, zeros keep coming', fontsize=10)
plt.savefig('figures/fig6-kernel-ladder.png', dpi=150, bbox_inches='tight'); plt.close()

# ---- fig7: what our tools are, in standard names
fig, ax = plt.subplots(figsize=(11, 3.4)); ax.axis('off'); ax.set_xlim(0, 11); ax.set_ylim(0, 3.4)
rows = [('Lee–Yang step identity  $(a+1)p\\,G_a = (\\kappa + s\\,d/ds)\\,G_{a+1}$', 'Laguerre\'s operator; $z\\,d/dz$ is the dilation generator (Berry–Keating "xp")'),
        ('our p = 4 trick  $(1 - D^2)$, applied many times', 'de Bruijn\'s heat flow $e^{-tD^2}$ (1950); Newman\'s constant $\\Lambda$ (1976)'),
        ('the k dial on the prime cube,  $\\prod (1 + k p^{-s})$', 'fugacity of the free primon gas; k = -1 is $(-1)^F$ (Spector 1990, Julia 1990)'),
        ('local factor as dial point, elliptic side', 'Hasse 1933, Weil 1948: Frobenius roots of size $\\sqrt{p}$')]
for i, (a, b) in enumerate(rows):
    y = 2.9 - i*0.8
    ax.text(0.1, y, a, color=B, fontsize=9.5, va='center'); ax.annotate('', xy=(6.0, y), xytext=(5.4, y), arrowprops=dict(arrowstyle='->', color=G))
    ax.text(6.1, y, b, color=O, fontsize=9.5, va='center')
ax.set_title('Our tools and their standard names (blue: ours; orange: in print)', fontsize=11)
plt.savefig('figures/fig7-standard-names.png', dpi=150, bbox_inches='tight'); plt.close()
print('figures done')
