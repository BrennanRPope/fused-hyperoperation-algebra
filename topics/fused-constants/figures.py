"""Figures for the fused-constants page.  Fused coordinates: origin at the unit corner (x -> x-1)."""
import json, math, numpy as np, matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
C1, C2, C3, C4, C5, C6 = '#257de0', '#f47230', '#00b894', '#f2a600', '#ee7fac', '#6a4ce2'
INK, INK2, GRID = '#191919', '#666666', '#d9d9d9'
plt.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 12, 'axes.edgecolor': INK2, 'axes.labelcolor': INK,
                     'xtick.color': INK2, 'ytick.color': INK2, 'savefig.dpi': 160, 'savefig.bbox': 'tight',
                     'savefig.facecolor': 'white', 'lines.linewidth': 2})
E, ETA = math.e, math.e ** (1 / math.e)
data = json.load(open('constants.json'))
lin = {int(k): float(v['c']) for k, v in data['linear'].items()}
linx = {int(k): float(v['x']) for k, v in data['linear'].items()}
NU = 1.63532449671527639934534
def clean(ax):
    ax.spines[['top', 'right']].set_visible(False); ax.grid(color=GRID, lw=0.6); ax.set_axisbelow(True)

def Hf(k, b, x):
    if x > 200: return float('inf')
    if k == 3: return b ** x
    if x <= 0: return 1 + x
    return Hf(k - 1, b, Hf(k, b, x - 1))

# ---------- C1: the kiss on rung 3, and the same kiss is e's slope test ----------
fig, axs = plt.subplots(1, 2, figsize=(13, 5.6))
ax = axs[0]
x = np.linspace(-1, 3.2, 600)
for c, col, lab in [(ETA - 1 - 0.07, C3, 'c = 0.375: crosses, tower settles'),
                    (ETA - 1, C1, 'c = e^(1/e) - 1 = 0.4447: just kisses'),
                    (ETA - 1 + 0.08, C2, 'c = 0.525: clears, tower runs off')]:
    ax.plot(x, (1 + c) ** (1 + x) - 1, color=col, label=lab)
ax.plot(x, x, color=INK, lw=1.2, ls='--', label='diagonal  y = x')
ax.plot([E - 1], [E - 1], 'o', ms=10, color=C1, mec='white', mew=2, zorder=5)
ax.annotate('kiss at x = e - 1 = 1.718\nslope exactly 1', (E - 1, E - 1), (2.0, 0.55), fontsize=11,
            arrowprops=dict(arrowstyle='->', color=INK2))
ax.set_xlim(-1, 3.2); ax.set_ylim(-1, 3.4); ax.set_xlabel('x  (fused: unit corner at 0)'); ax.set_ylabel('c (3) x  =  (1+c)^(1+x) - 1')
ax.set_title('Fused rung 3: the tower test', fontsize=12); ax.legend(frameon=False, fontsize=9.5, loc='upper left'); clean(ax)
ax = axs[1]
u = np.linspace(-1.2, 1.3, 400)
ax.plot(u, np.exp(u) - 1, color=C1, label='e^u - 1  (the rung-3 kiss, shrunk by e and slid)')
ax.plot(u, np.log1p(np.clip(u, -0.999, None)), color=C6, label='ln(1+u)  (mult-add\'s own log)')
ax.plot(u, u, color=INK, lw=1.2, ls='--', label='diagonal')
ax.plot([0], [0], 'o', ms=10, color=C1, mec='white', mew=2, zorder=5)
ax.annotate("e's slope-one test:\nslope 1 at the unit corner", (0, 0), (0.25, -0.9), fontsize=11, arrowprops=dict(arrowstyle='->', color=INK2))
ax.set_xlim(-1.2, 1.3); ax.set_ylim(-1.2, 1.6); ax.set_xlabel('u'); ax.set_title('The same kiss, seen at the corner', fontsize=12)
ax.legend(frameon=False, fontsize=9.5, loc='upper left'); clean(ax)
plt.savefig('figures/c1_kiss.png'); plt.close()

# ---------- C2: the cube whose volume equals its dimension ----------
fig = plt.figure(figsize=(13, 5.4))
ax = fig.add_subplot(121)
d = np.linspace(0.3, 5, 500)
for s, col, lab in [(1.35, C3, 'side 1.35: volume catches its dimension twice'),
                    (ETA, C1, 'side e^(1/e) = 1.4447: touches once, at d = e'),
                    (1.55, C2, 'side 1.55: volume always ahead')]:
    ax.plot(d, s ** d, color=col, label=lab)
ax.plot(d, d, color=INK, lw=1.2, ls='--', label='volume = number of dimensions')
ax.plot([E], [E], 'o', ms=10, color=C1, mec='white', mew=2, zorder=5)
for k in (1, 2, 3):
    ax.plot([k], [ETA ** k], 's', ms=7, color=C1)
    ax.annotate(f'd={k}: {ETA**k:.3f}', (k, ETA ** k), (k + 0.12, ETA ** k - 0.55), fontsize=10, color=INK2)
ax.set_xlabel('dimension d'); ax.set_ylabel('volume of a cube of side s'); ax.set_ylim(0, 6); ax.set_xlim(0.3, 5)
ax.set_title('e^(1/e) is the largest side a cube can have\nand still have volume equal to its dimension', fontsize=12)
ax.legend(frameon=False, fontsize=9.5, loc='upper left'); clean(ax)
ax = fig.add_subplot(122, projection='3d')
s = ETA
def box(ax, o, sz, fc, ec, a):
    X, Y, Z = o; dx, dy, dz = sz
    v = np.array([[X, Y, Z], [X+dx, Y, Z], [X+dx, Y+dy, Z], [X, Y+dy, Z], [X, Y, Z+dz], [X+dx, Y, Z+dz], [X+dx, Y+dy, Z+dz], [X, Y+dy, Z+dz]])
    f = [[0,1,2,3],[4,5,6,7],[0,1,5,4],[2,3,7,6],[1,2,6,5],[0,3,7,4]]
    ax.add_collection3d(Poly3DCollection([v[i] for i in f], facecolors=fc, edgecolors=ec, alpha=a, linewidths=0.8))
box(ax, (0, 0, 0), (s, s, s), C1, C1, 0.10)
box(ax, (0, 0, 0), (1, 1, 1), '#bbbbbb', INK2, 0.35)
ax.text(s/2, -0.35, -0.15, f'{s:.4f}', color=INK, fontsize=10)
ax.text2D(0.02, 0.02, f'side e^(1/e) = 1.4447, volume {s**3:.4f} (just above 3)\nunit corner cube in gray', transform=ax.transAxes, fontsize=10, color=INK)
ax.set_xlim(0, 1.6); ax.set_ylim(0, 1.6); ax.set_zlim(0, 1.6); ax.set_box_aspect((1, 1, 1)); ax.view_init(22, -55); ax.set_axis_off()
ax.set_title('Fused: extra side 0.4447 beyond the unit corner;\nextra volume equals extra dimension at d - 1 = e - 1', fontsize=11)
plt.savefig('figures/c2_cube.png'); plt.close()

# ---------- C3: the constants climb to 1 (fused) ----------
fig, axs = plt.subplots(1, 2, figsize=(13, 5.2), gridspec_kw={'width_ratios': [1.5, 1]})
ax = axs[0]
ks = list(range(1, 21))
vals = [-1, 0, ETA - 1] + [lin[k] - 1 for k in range(4, 21)]
ax.axhline(1, color=INK2, lw=1.2, ls='--'); ax.text(20.3, 1.0, 'limit 1\n(ordinary 2)', va='center', fontsize=10, color=INK2)
ax.plot(ks, vals, '-', color=C1, lw=1.5, alpha=0.5)
ax.plot(ks[:3], vals[:3], 'o', ms=9, color=C1, mec='white', mew=1.5, label='exact: -1, 0, e^(1/e) - 1')
ax.plot(ks[3:], vals[3:], 'o', ms=8, color=C2, mec='white', mew=1.5, label='linear real-height extension')
ax.plot([4], [NU - 1], 'D', ms=9, color=C6, mec='white', mew=1.5, label='rung 4, analytic extension: 0.6353')
for k, lab in [(1, 'add: -1'), (2, 'mult-add: 0'), (3, '0.4447'), (5, '0.7570'), (8, '0.9199')]:
    ax.annotate(lab, (k, vals[k - 1]), (k + 0.4, vals[k - 1] - 0.14), fontsize=10, color=INK2)
ax.set_xticks(ks); ax.set_xlabel('rung k  (the map x -> c (k) x; its tower lives on rung k+1)')
ax.set_ylabel('fused natural constant c_k'); ax.set_ylim(-1.15, 1.15)
ax.set_title('Each rung\'s natural constant, fused coordinates', fontsize=12); ax.legend(frameon=False, fontsize=10, loc='lower right'); clean(ax)
ax = axs[1]
k2 = list(range(4, 21)); gap = [2 - lin[k] for k in k2]
ax.semilogy(k2, gap, 'o-', color=C2, ms=7, lw=1.5)
ax.semilogy(k2, [gap[0] * math.log(2) ** (k - 4) for k in k2], color=INK2, ls='--', lw=1.2, label='shrinks by ln 2 = 0.693 per rung')
ax.set_xlabel('rung k'); ax.set_ylabel('distance to the limit, 1 - c_k'); ax.set_title('The gap closes about 30% per rung', fontsize=12)
ax.legend(frameon=False, fontsize=10); clean(ax)
plt.savefig('figures/c3_ladder.png'); plt.close()

# ---------- C4: why the limit is 1: the curves pinch onto the corner (1,1) ----------
fig, axs = plt.subplots(1, 2, figsize=(13, 5.4))
xs = np.linspace(-0.999, 2.6, 900)
cols = [C1, C3, C4, C2, C6]
ax = axs[0]
for k, col in zip((3, 5, 8, 12, 20), cols):
    ax.plot(xs, [Hf(k, 1.8, t + 1) - 1 for t in xs], color=col, label=f'rung {k}')
ax.plot(xs, xs, color=INK, lw=1.2, ls='--', label='diagonal')
ax.axhline(0.8, color=INK2, lw=0.8, ls=':')
ax.plot([0.8], [0.8], 'o', ms=9, color=INK, mec='white', mew=1.5)
ax.annotate('every rung high enough\ncrosses here: x = c = 0.8', (0.8, 0.8), (1.25, 0.15), fontsize=10, arrowprops=dict(arrowstyle='->', color=INK2))
ax.set_ylim(-1, 3); ax.set_xlim(-1, 2.6); ax.set_xlabel('x (fused)'); ax.set_ylabel('0.8 (k) x')
ax.set_title('Fused base 0.8: high rungs flatten onto the step\n(1.8)^x - 1, then flat at 0.8', fontsize=12); ax.legend(frameon=False, fontsize=9.5, loc='upper left'); clean(ax)
ax = axs[1]
xs2 = np.linspace(-0.999, 1.6, 900)
for k, col in zip((4, 8, 16, 30, 60), cols):
    ax.plot(xs2, [Hf(k, 2.0, t + 1) - 1 for t in xs2], color=col, label=f'rung {k}')
ax.plot(xs2, xs2, color=INK, lw=1.2, ls='--', label='diagonal')
ax.plot([1], [1], 'o', ms=10, color=INK, mec='white', mew=1.5)
ax.annotate('(1, 1): 1 (k) 1 = 3 on every rung,\nand the curves pinch toward this corner', (1, 1), (-0.9, 1.75), fontsize=10, arrowprops=dict(arrowstyle='->', color=INK2))
ax.set_ylim(-1, 3); ax.set_xlim(-1, 1.6); ax.set_xlabel('x (fused)'); ax.set_ylabel('1 (k) x')
ax.set_title('Fused base 1: never touches, but the gap\nshrinks to 0 at the corner (1,1)', fontsize=12); ax.legend(frameon=False, fontsize=9.5, loc='lower right'); clean(ax)
plt.savefig('figures/c4_pinch.png'); plt.close()

# ---------- C5: the unit square doubled ----------
fig, ax = plt.subplots(figsize=(6, 6))
ax.add_patch(plt.Rectangle((-1, -1), 2, 2, fc=C1, alpha=0.12, ec=C1, lw=2))
ax.add_patch(plt.Rectangle((-1, -1), 1, 1, fc='#bbbbbb', alpha=0.6, ec=INK2, lw=1.5))
ax.text(-0.5, -0.5, 'unit\ncorner', ha='center', va='center', fontsize=11)
for (px, py) in ((0.5, -0.5), (-0.5, 0.5), (0.3, 0.65)): ax.text(px, py, 'new', ha='center', va='center', fontsize=11, color=INK2)
ax.text(0.15, -1.15, '1 (+) 1 = 3: the three new cells', ha='center', fontsize=11)
ax.plot([1], [1], 'o', ms=12, color=INK); ax.text(1.05, 1.08, '(1, 1)', fontsize=12)
ax.plot([-1, 1.3], [-1, 1.3], color=INK, ls='--', lw=1.2)
for t in (0, 1): ax.plot([t, t], [-1, 1], color=C1, lw=1); ax.plot([-1, 1], [t, t], color=C1, lw=1)
ax.set_xlim(-1.2, 1.5); ax.set_ylim(-1.2, 1.5); ax.set_aspect('equal'); ax.axis('off')
ax.set_title('Side 2 = fused 1: 2+2 = 2·2 = 2^2 = 2^^2 = ... = 4\nThe one square every rung agrees on', fontsize=12)
plt.savefig('figures/c5_square.png'); plt.close()
print('ok')
