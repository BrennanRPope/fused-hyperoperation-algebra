"""Pictures for riemann-neighbours/README.md. Reads results.json (run compute.py first).
Colour key: blue = our reading, orange = already known, grey = reference, red = where a road stops or fails."""
import json, os, numpy as np, matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt, networkx as nx, mpmath as mp
from matplotlib.patches import FancyBboxPatch, Ellipse
HERE = os.path.dirname(os.path.abspath(__file__)); FIG = os.path.join(HERE, 'figures')
R = json.load(open(os.path.join(HERE, 'results.json')))
BLUE, ORANGE, GREY, RED, GREEN = '#1f6fb4', '#e8790c', '#8a8a8a', '#c62828', '#2e7d32'
plt.rcParams.update({'font.size': 10, 'axes.titlesize': 11})
def save(fig, name):
    fig.savefig(os.path.join(FIG, name), dpi=140, bbox_inches='tight'); plt.close(fig)

# ---------- fig0: Riemann's layer and its neighbours ----------
fig, ax = plt.subplots(figsize=(13.4, 7.6)); ax.set_xlim(0, 13.3); ax.set_ylim(0, 7.6); ax.axis('off')
layers = [(5.9, 'Layer 1: the k dial alone on a cube of "primes"\n(the fold  ∏(1 + k·x_p))'),
          (3.55, 'Layer 2: k with d, the conic plane\n(quadratic factors 1 + k·u + d·u²)'),
          (1.2, 'Layer 3: order turned on (q dial)\n(affine maps x ↦ n·x + r)')]
for y, lab in layers:
    ax.add_patch(FancyBboxPatch((0.15, y - 0.95), 13.0, 1.9, boxstyle='round,pad=0.02', fc='#f4f4f4', ec='#cccccc'))
    ax.text(0.3, y + 0.62, lab, fontsize=9.5, color='#444', va='center')
def node(x, y, title, status, sub, col):
    ax.add_patch(FancyBboxPatch((x - 1.08, y - 0.5), 2.16, 1.0, boxstyle='round,pad=0.03', fc='white', ec=col, lw=2))
    ax.text(x, y + 0.22, title, ha='center', fontsize=9.5, weight='bold')
    ax.text(x, y - 0.02, status, ha='center', fontsize=8.5, color=col)
    ax.text(x, y - 0.27, sub, ha='center', fontsize=7.6, color='#333')
node(4.25, 5.75, 'Riemann ζ', 'OPEN', 'primes 2, 3, 5, …', 'black')
node(6.75, 5.75, 'Ihara (graphs)', 'proved for Ramanujan graphs', 'primes = prime cycles', ORANGE)
node(9.3, 5.75, 'Lee–Yang cube', 'proved (1952)', 'cube + pair coupling', ORANGE)
node(11.85, 5.75, 'Davenport–Heilbronn', 'FALSE (zero off line)', 'no fold, same symmetry', RED)
node(4.25, 3.4, 'Weil (curves / F_q)', 'proved (Weil 1948)', 'Frobenius pairs', ORANGE)
node(6.75, 3.4, 'Selberg (surfaces)', 'proved, λ ≥ 1/4 part', 'Laplace eigenvalues', ORANGE)
node(9.3, 3.4, 'Dirichlet, Dedekind,', 'OPEN (GRH)', 'Hecke, Selberg class', 'black')
node(11.85, 3.4, 'Ramanujan–Petersson', 'proved (Deligne 1974)', 'modular forms', ORANGE)
node(6.75, 1.05, 'Bost–Connes', 'no RH content', 'ax + b Hecke algebra', GREY)
node(9.3, 1.05, 'Connes trace formula', 'RH ⇔ a positivity', 'adele class space', GREY)
ax.annotate('', xy=(6.75, 3.92), xytext=(6.75, 5.23), arrowprops=dict(arrowstyle='->', color=BLUE, lw=2))
ax.text(6.85, 4.62, 'the fold collapses to quadratic factors\n(the same happens for curves and surfaces)', color=BLUE, fontsize=8.5)
ax.annotate('', xy=(4.25, 3.92), xytext=(4.25, 5.23), arrowprops=dict(arrowstyle='->', color=RED, lw=2, ls=':'))
ax.text(2.0, 4.62, 'for ζ no such collapse\nis known (Hilbert–Pólya)', color=RED, fontsize=8.5)
ax.text(1.0, 0.12, 'Orange = proved neighbour, red = fails or stops, black = open, grey = reference. Blue = our reading.', fontsize=8.5, color='#444')
ax.set_title("Figure 0. Riemann's layer of the dial map and its neighbours, most fundamental at the top", loc='left')
save(fig, 'fig0-neighbours-map.png')

# ---------- fig1: one ruler for every neighbour ----------
fig, ax = plt.subplots(figsize=(12, 6.2))
rows = ['ζ (pairs of\nprime factors)', 'Lee–Yang,\ntwo spins', 'Selberg\n(surfaces)', 'Ihara\n(graphs)', 'Weil\n(curve over F_p)']
ax.axvspan(-2, 2, color='#e3f0e3'); ax.axvline(-2, color=GREEN); ax.axvline(2, color=GREEN)
ax.text(0, 4.62, 'elliptic segment |k̂| ≤ 2: roots on the circle, the "Riemann hypothesis holds" side', ha='center', color=GREEN, fontsize=9)
kh = np.array(R['W4_khat']); ax.scatter(kh, 4 + 0.25*(np.random.default_rng(1).random(len(kh)) - .5), s=6, color=ORANGE, label='known theorem')
ax.text(2.15, 4.0, f'all {len(kh)} primes\n(Hasse 1933)', fontsize=8, va='center')
g = R['G4_graphs']
for name, y0, col in [('random60', 3.12, ORANGE), ('twocluster60', 2.88, RED)]:
    ev = np.array(g[name]['eig'][1:]); kk = -ev/np.sqrt(2)
    ax.scatter(kk, np.full(len(kk), y0), s=10, color=col)
ax.text(2.15, 3.12, 'random 3-regular graph: Ramanujan', fontsize=8, va='center')
ax.text(2.15, 2.88, 'two-cluster graph: one factor outside', fontsize=8, va='center', color=RED)
lam1 = 0.25 + 9.533695261353557**2
ax.scatter([-1/np.sqrt(lam1)], [2], s=40, color=ORANGE, marker='D'); ax.text(-1/np.sqrt(lam1), 2.17, 'PSL(2,ℤ), λ₁ ≈ 91.14 (quoted)', fontsize=8, ha='center')
ax.scatter([-1/np.sqrt(975/4096)], [1.88], s=30, color=GREY, marker='|'); ax.text(-1/np.sqrt(975/4096) - .05, 1.72, 'Kim–Sarnak bound\nλ₁ ≥ 975/4096', fontsize=7.5, ha='right')
ax.scatter([-1/np.sqrt(0.05)], [2], s=30, color=RED, marker='x'); ax.text(-1/np.sqrt(0.05), 2.17, 'small eigenvalue\n(allowed on some surfaces)', fontsize=7.5, ha='center', color=RED)
aa = np.linspace(-1.6, 1.6, 33); ax.scatter(2*aa, np.full(len(aa), 1), s=10, c=[ORANGE if abs(x) <= 1 else RED for x in aa])
ax.text(2.15, 1.18, 'k̂ = 2a: ferromagnetic |a| ≤ 1 inside', fontsize=8)
ps = np.array([2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47])
ax.scatter(-(np.sqrt(ps) + 1/np.sqrt(ps)), np.zeros(len(ps)), s=12, color=RED)
ax.text(-2.15, -0.22, '(1 − u)(1 − p·u): always |k̂| ≥ 2.12 (round 2)', fontsize=8, ha='right', color=RED)
ax.set_yticks(range(5)); ax.set_yticklabels(rows); ax.set_xlim(-7.5, 7.5); ax.set_ylim(-0.5, 4.85)
ax.set_xlabel('k̂ = k / √d  for each quadratic factor 1 + k·u + d·u²')
ax.set_title('Figure 1. One ruler for every neighbour: each spectral factor is a dial point, and "RH" = every point on the green segment', loc='left')
save(fig, 'fig1-one-ruler.png')

# ---------- fig2: graphs ----------
fig, axs = plt.subplots(1, 3, figsize=(15, 4.9))
P = nx.petersen_graph(); pos = nx.shell_layout(P, [range(5, 10), range(5)])
nx.draw(P, pos, ax=axs[0], node_size=120, node_color=BLUE, edge_color='#555')
axs[0].set_title('(a) Petersen graph: 3-regular, q = 2.\nIts "primes" are the prime cycles.', loc='left')
th = np.linspace(0, 2*np.pi, 300)
for name, col, lab in [('random60', ORANGE, 'random 3-regular (Ramanujan)'), ('twocluster60', RED, 'two-cluster (not Ramanujan)')]:
    p = np.array(g[name]['poles']); axs[1].scatter(p[:, 0], p[:, 1], s=8, color=col, label=lab, alpha=.8)
axs[1].plot(np.cos(th)/np.sqrt(2), np.sin(th)/np.sqrt(2), color=GREEN, lw=1, label='|u| = 1/√q')
axs[1].set_aspect('equal'); axs[1].legend(fontsize=7.5, loc='lower left'); axs[1].set_xlim(-1.2, 1.2); axs[1].set_ylim(-1.2, 1.2)
axs[1].annotate('off the circle', xy=(0.53, 0.0), xytext=(0.2, 0.35), color=RED, fontsize=8, arrowprops=dict(arrowstyle='->', color=RED))
axs[1].text(-1.15, 0.12, 'trivial poles at u = ±1, ±1/2', fontsize=7.5, color=GREY)
axs[1].set_title('(b) Poles of the Ihara zeta in u.\nRH for graphs = every non-trivial pole on the circle.', loc='left')
pc = R['G2_prime_cycles']; ls = np.arange(5, 21); cnt = np.array([pc[str(l)] for l in ls])
axs[2].semilogy(ls, cnt, 'o', color=ORANGE, label='prime cycles of length ℓ (counted)')
axs[2].semilogy(ls, 2.0**ls/ls, '-', color=GREY, label='2^ℓ / ℓ (graph prime number theorem)')
axs[2].legend(fontsize=7.5); axs[2].set_xlabel('length ℓ')
axs[2].text(5.2, 2e4, 'the fold ∏(1 − u^ℓ) over all these\ncycles is exactly det(I − uB):\na polynomial of degree 30', color=BLUE, fontsize=8.5)
axs[2].set_title('(c) Infinitely many primes, one finite\ndeterminant (check G2)', loc='left')
fig.suptitle('Figure 2. Ihara zeta: Riemann\'s own fold, on a graph', x=0.02, y=1.06, ha='left')
save(fig, 'fig2-graphs.png')

# ---------- fig3: curves ----------
fig, axs = plt.subplots(1, 2, figsize=(12.5, 5))
p, a = 499, R['W3_ap']['499']
X, Y = np.meshgrid(np.linspace(-60, 60, 400), np.linspace(-3, 3, 400))
Q = X**2 - a*X*Y + p*Y**2
axs[0].contour(X, Y, Q, levels=[100, 400, 900, 1600], colors=ORANGE)
mm, nn = np.meshgrid(np.arange(-60, 61), np.arange(-3, 4)); axs[0].scatter(mm, nn, s=2, color=GREY)
axs[0].scatter([1], [1], s=40, color=BLUE, zorder=5); axs[0].text(3, 1.15, f'deg(1 − F) = #E(F_p) = {1 - a + p}', color=BLUE, fontsize=8.5)
axs[0].set_xlabel('m'); axs[0].set_ylabel('n')
axs[0].set_title(f'(a) Degree form m² − a·mn + p·n² for p = {p}, a = {a}:\nellipses, so positive definite. That is Hasse\'s proof.', loc='left')
good = np.array(sorted(int(p_) for p_ in R['W3_ap'].keys()))
kh = np.array(R['W4_khat']); pp = np.array(R.get('W4_primes', []))
ks = np.linspace(-95, 95, 400); axs[1].plot(ks, ks**2/4, color=GREY, label='Δ = k² − 4d = 0')
axs[1].fill_between(ks, ks**2/4, 2100, color='#e3f0e3')
pp = np.array(R['W4_primes']); kk = kh*np.sqrt(pp)
axs[1].scatter(kk, pp, s=5, color=ORANGE, label='Frobenius dial point (k, d) = (−a_p, p)')
axs[1].text(-30, 1700, 'elliptic side:\nnorm form positive definite', color=GREEN, fontsize=9, ha='center', bbox=dict(fc='white', ec='none', alpha=.85))
axs[1].text(60, 300, 'hyperbolic side:\nnorm form indefinite', color=RED, fontsize=9, ha='center', bbox=dict(fc='white', ec='none', alpha=.85))
axs[1].set_xlabel('k'); axs[1].set_ylabel('d'); axs[1].set_ylim(0, 2100); axs[1].legend(fontsize=7.5, loc='lower left')
axs[1].set_title('(b) The (k, d) plane: our law\'s norm a² + k·ab + d·b² is the degree form,\nand Weil\'s RH is "every Frobenius point on the elliptic side".', loc='left')
fig.suptitle('Figure 3. Curves over finite fields: the proof ingredient is a positive-definite norm, which our dial plane already has', x=0.02, y=1.06, ha='left')
save(fig, 'fig3-curves.png')

# ---------- fig4: Lee-Yang on the cube ----------
fig, axs = plt.subplots(1, 3, figsize=(15, 4.9))
rf = np.array(R['L1_roots']['ferro']); ra = np.array(R['L1_roots']['anti'])
axs[0].plot(np.cos(th), np.sin(th), color=GREY, lw=.8)
axs[0].scatter(rf[:, 0], rf[:, 1], color=ORANGE, s=25, label='ferromagnetic couplings')
ra2 = ra[np.abs(ra[:, 0] + 1j*ra[:, 1]) < 3]
axs[0].scatter(ra2[:, 0], ra2[:, 1], color=RED, s=25, marker='x', label='antiferromagnetic (some far off)')
axs[0].set_aspect('equal'); axs[0].set_xlim(-2.5, 2.5); axs[0].set_ylim(-2.5, 2.5); axs[0].legend(fontsize=7.5, loc='lower left')
axs[0].set_title('(a) Nine spins: zeros in z on |z| = 1\n(Lee–Yang 1952)', loc='left')
zf = np.array(R['L4_ferro_zeros']); za = np.array(R['L4_anti_zeros'])
axs[1].scatter(zf[:, 0], zf[:, 1], s=10, color=BLUE, label=f'ferromagnetic: {len(zf)} zeros, all on Re s = 1/2')
axs[1].scatter(za[:, 0], za[:, 1], s=10, color=RED, marker='x', label='antiferromagnetic: off the line')
axs[1].axvline(0.5, color=GREY, lw=.8); axs[1].set_ylim(0, 37); axs[1].legend(fontsize=7.5, loc='upper right'); axs[1].set_xlim(-1.6, 2.6)
axs[1].set_xlabel('Re s'); axs[1].set_ylabel('Im s')
axs[1].set_title('(b) A coupled prime cube, primes 2 to 19,\nz_p = p^(1/2 − s): the coupling keeps the line', loc='left')
axs[2].axis('off')
axs[2].text(0, 0.95, '(c) What the coupling is, in our terms', weight='bold', fontsize=10, transform=axs[2].transAxes)
txt = ['Two directions:  1 + a(z₁ + z₂) + z₁z₂ = 1 + a·(z₁ ∘ₖ z₂)',
       'with k = 1/a. So two coupled spins are our fused law,',
       'and Lee–Yang\'s condition |a| ≤ 1 is our elliptic side.', '',
       'Three or more directions: the weights a^(|S|(n−|S|)) are',
       'not a product form (n = 3) and not an order-2 recurrence',
       '(n ≥ 4, Hankel test): no single k, d or t fold gives them.',
       'But they are k placed on each edge (k_ij = a_ij^−2), our',
       '"dial gas" in topics/lee-yang. Not a new dial (check L6).', '',
       "ζ's cube has no coupling at all (the primes are",
       'independent), and the weights, not the coupling, set',
       'which line: p^(−s) gives Re s = 0, p^(1/2−s) gives 1/2.']
for i, l in enumerate(txt):
    axs[2].text(0, 0.84 - i*0.07, l, fontsize=8.6, transform=axs[2].transAxes, color=BLUE if i < 3 else '#222')
fig.suptitle('Figure 4. Lee–Yang lives on the cube itself: a pair coupling between directions', x=0.02, y=1.06, ha='left')
save(fig, 'fig4-lee-yang-cube.png')

# ---------- fig5: Davenport-Heilbronn ----------
mp.mp.dps = 15
chi = {1: 1, 2: 1j, 3: -1j, 4: -1}; chib = {k_: complex(v).conjugate() for k_, v in chi.items()}
L5 = lambda s, c: mp.power(5, -s)*sum(c[k_]*mp.zeta(s, mp.mpf(k_)/5) for k_ in range(1, 5))
kap = (mp.sqrt(10 - 2*mp.sqrt(5)) - 2)/(mp.sqrt(5) - 1)
f = lambda s: (1 - 1j*kap)/2*L5(s, chi) + (1 + 1j*kap)/2*L5(s, chib)
S = np.linspace(0.1, 1.3, 49); T = np.linspace(83.5, 88.5, 81)
V = np.array([[float(abs(f(mp.mpc(s, t)))) for s in S] for t in T])
fig, ax = plt.subplots(figsize=(7.5, 5.6))
cs = ax.contourf(S, T, np.log10(V), levels=24, cmap='Greys_r'); fig.colorbar(cs, label='log₁₀ |f(s)|')
ax.axvline(0.5, color=BLUE, ls='--', label='critical line'); ax.axvline(1.0, color=GREY, ls=':', label='edge of the strip')
zr = R['H1_DH_zero']; ax.scatter([zr[0]], [zr[1]], color=RED, s=80, marker='x', label=f'zero at {zr[0]:.6f} + {zr[1]:.6f} i')
ax.annotate('its mirror at 1 − s̄\n(the functional equation)', xy=(1 - zr[0], zr[1]), xytext=(0.12, 84.4), color=RED, fontsize=8, arrowprops=dict(arrowstyle='->', color=RED))
ax.legend(fontsize=8, loc='upper left'); ax.set_xlabel('Re s'); ax.set_ylabel('Im s')
ax.set_title('Figure 5. Davenport–Heilbronn: same symmetry as ζ, no fold\n(no Euler product), and a zero off the line. The fold is needed.', loc='left')
save(fig, 'fig5-davenport-heilbronn.png')

# ---------- fig6: ingredients ----------
fig, ax = plt.subplots(figsize=(13.5, 5.6)); ax.axis('off')
cols = ['Neighbour', 'Fold becomes…', 'Ingredient that pins the circle', 'Can our dials say it?', 'For ζ?']
rows = [['Ihara (graphs)', 'det(I − uB), finite', 'Ramanujan bound: Eichler–Shimura + Weil (LPS, Margulis 1988)\nor interlacing (Marcus–Spielman–Srivastava 2015)', 'Yes: |k̂| ≤ 2 on the ruler', 'No operator B known'],
        ['Weil (curves)', 'det(1 − uF | H¹), finite', 'Positivity: Castelnuovo–Severi / Hodge index\n(Weil 1948); Hasse: degree form > 0', 'Yes: norm a² + kab + db²\npositive definite', 'No Frobenius, no norm form'],
        ['Selberg (surfaces)', 'regularized det(Δ)', 'Self-adjoint Laplacian: λ real ≥ 0\n(Selberg 1956 trace formula)', 'Yes: factor s² − s + λ,\nelliptic iff λ ≥ 1/4', 'No Laplacian known'],
        ['Lee–Yang', 'polynomial on the cube', 'Ferromagnetic coupling + Asano contraction', 'Yes: n = 2 is ∘ₖ; n ≥ 3 is\nk on each edge (dial gas)', 'ζ\'s cube is uncoupled'],
        ['Davenport–Heilbronn', '(no fold)', 'none: RH fails', '—', 'Shows the fold is needed'],
        ['Bost–Connes', 'partition function ζ(β)', 'none for zeros; phase transition at β = 1', 'Yes: q dial on ax + b maps', 'Sees the pole, not the zeros']]
tb = ax.table(cellText=rows, colLabels=cols, loc='center', cellLoc='left', colWidths=[.13, .16, .33, .2, .18])
tb.auto_set_font_size(False); tb.set_fontsize(8.3); tb.scale(1, 3.0)
for (r_, c_), cell in tb.get_celld().items():
    if r_ == 0: cell.set_facecolor('#e8e8e8'); cell.set_text_props(weight='bold')
    elif c_ == 3: cell.set_text_props(color=BLUE)
    elif c_ == 4: cell.set_text_props(color=RED)
ax.set_title('Figure 6. What each proof uses, whether our dials can express it, and where it stops for ζ', loc='left')
save(fig, 'fig6-ingredients.png')

# ---------- fig7: Mathlib hooks ----------
fig, ax = plt.subplots(figsize=(13.5, 6.4)); ax.set_xlim(0, 13.5); ax.set_ylim(0, 6.4); ax.axis('off')
boxes = [(0.2, 3.4, 'Riemann ζ', ['riemannZeta_eulerProduct_tprod', 'LSeries_zeta_mul_Lseries_moebius  (1/ζ = k = −1 fold)', 'LSeriesSummable_moebius_iff  (fold converges ⇔ Re s > 1)', 'completedRiemannZeta_one_sub, RiemannHypothesis', 'riemannZeta_ne_zero_of_one_le_re, hardyZ_eq_zero_iff']),
         (6.9, 3.4, 'Dirichlet, Dedekind, Hecke, modular L', ['DirichletCharacter.LSeries_eulerProduct_tprod', 'completedLFunction_one_sub, LFunction_ne_zero_of_one_le_re', 'norm_gaussSum_stdAddChar  (|Gauss sum| = √N)', 'NumberField.dedekindZeta (class number formula)', 'ModularForm weakFEPair, hasSum_L']),
         (0.2, 0.3, 'Curves and graphs', ['WeierstrassCurve.localPolynomial = 1 − a·T + q·T²', '   (exactly our dial point; Hasse bound itself: not found)', 'Matrix.discr_fin_two: discr = trace² − 4·det = our Δ', 'SimpleGraph.adjMatrix_pow_apply_eq_card_walk', 'isHermitian_adjMatrix, posSemidef_lapMatrix']),
         (6.9, 0.3, 'Not found in Mathlib (gaps)', ['Hasse bound, Weil\'s theorem for curves', 'Ihara zeta, non-backtracking matrix, Ramanujan graphs', 'Lee–Yang circle theorem, Asano contraction', 'Selberg zeta and trace formula', 'Bost–Connes system, de Bruijn–Newman flow'])]
for x, y, t, items in boxes:
    col = RED if 'Not found' in t else ORANGE
    ax.add_patch(FancyBboxPatch((x, y), 6.4, 2.85, boxstyle='round,pad=0.03', fc='white', ec=col, lw=1.8))
    ax.text(x + 0.15, y + 2.55, t, weight='bold', fontsize=10)
    for i, it in enumerate(items):
        ax.text(x + 0.2, y + 2.1 - i*0.45, it, fontsize=8.4, family='monospace' if not it.startswith('   ') else None)
ax.set_title('Figure 7. Mathlib (v4.37.0-rc1, the version the repo pins) pieces each neighbour could depend on, and the gaps', loc='left')
save(fig, 'fig7-mathlib-hooks.png')

# ---------- fig8: new dial lookout ----------
fig, axs = plt.subplots(1, 2, figsize=(13, 4.8))
ns = np.arange(4, 9); hv = []
import sympy as sp
for nn in ns:
    c = [sp.Rational(1, 2)**(j*(nn-j)) for j in range(nn+1)]
    hv.append(abs(float(sp.Matrix(3, 3, lambda i, j: c[i+j]).det())))
axs[0].semilogy(ns, hv, 'o-', color=BLUE)
axs[0].axhline(1e-16, color=GREY, ls=':'); axs[0].text(4, 3e-16, 'exact rational values, all nonzero; zero would mean "a (k, d) fold in disguise"', fontsize=7.5, color=GREY)
axs[0].set_xlabel('number of coupled directions n'); axs[0].set_ylabel('|3×3 Hankel determinant| of the weights')
axs[0].set_title('(a) Lee–Yang coupling a = 1/2: not one (k, d) fold over the cube for n ≥ 4,\nbut it is k placed on each edge (k = a^−2): the dial gas, check L6.', loc='left')
tt = np.linspace(0, 1, 50); deff = (1 - tt)*(2 + tt)
for lamv, col, lab in [(1.0, ORANGE, 'Petersen λ = 1'), (-2.0, ORANGE, 'Petersen λ = −2'), (2.9, RED, 'λ = 2.9 (outside already at t = 0)')]:
    axs[1].plot(-lamv/np.sqrt(deff[:-1]), tt[:-1], color=col, label=lab, ls='-' if lamv != -2.0 else '--')
axs[1].legend(fontsize=7.5, loc='lower left')
axs[1].axvspan(-2, 2, color='#e3f0e3'); axs[1].set_xlim(-8, 8)
axs[1].set_xlabel('k̂ = −λ/√d_eff'); axs[1].set_ylabel("Bartholdi's bump weight t")
axs[1].set_title("(b) Bartholdi's bump weight t only changes d: d_eff = (1 − t)(q + t).\nThe d dial in disguise (check B2).", loc='left')
fig.suptitle('Figure 8. Looking for new dials in the neighbours', x=0.02, y=1.06, ha='left')
save(fig, 'fig8-new-dials.png')
print('figures done')
