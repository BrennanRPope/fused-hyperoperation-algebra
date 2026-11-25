"""Figures for the closing-rate section."""
import json, math, numpy as np, matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
C1, C2, C3, C4, C5, C6 = '#257de0', '#f47230', '#00b894', '#f2a600', '#ee7fac', '#6a4ce2'
INK, INK2, GRID = '#191919', '#666666', '#d9d9d9'
plt.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 12, 'axes.edgecolor': INK2, 'axes.labelcolor': INK,
                     'xtick.color': INK2, 'ytick.color': INK2, 'savefig.dpi': 160, 'savefig.bbox': 'tight',
                     'savefig.facecolor': 'white', 'lines.linewidth': 2})
L2, L4 = math.log(2), math.log(4)
def clean(ax):
    ax.spines[['top', 'right']].set_visible(False); ax.grid(color=GRID, lw=0.6); ax.set_axisbelow(True)

# ---------- C6: the corner (1,1) ----------
fig, axs = plt.subplots(1, 2, figsize=(13, 5.8))
ax = axs[0]
y = np.linspace(-0.15, 1.35, 400)
ax.add_patch(plt.Rectangle((0, 0), 1, 1, fc='#eef4fc', ec=C1, lw=1, zorder=0))
ax.plot(y, 2 ** y - 1, color=C1, label='y -> 2^y - 1   (1 mult-added to itself y times)')
ax.plot(y, y, color=INK, lw=1.2, ls='--', label='diagonal')
t = np.linspace(-0.15, 0.5, 10); ax.plot(t, L2 * t, color=C3, lw=1.6, ls=':')
t = np.linspace(0.75, 1.3, 10); ax.plot(t, 1 + L4 * (t - 1), color=C2, lw=1.6, ls=':')
# cobweb at c = 0.97 from the corner y = 1
c = 0.97; z = 1.0; xs, ys = [z], [z]
for _ in range(9):
    w = (1 + c) ** z - 1; xs += [z, w]; ys += [w, w]; z = w
ax.plot(xs, ys, color=C6, lw=1.1, label='steps at c = 0.97, starting at the corner')
ax.plot([0, 1], [0, 1], 'o', ms=9, color=INK, mec='white', mew=2, zorder=5)
ax.annotate('slope ln 2 = 0.69\nsteps pull in', (0.02, 0.01), (0.15, -0.12), fontsize=10.5, color=C3)
ax.annotate('slope ln 4 = ln(2*2) = 1.39\nsteps push out', (1, 1), (1.03, 0.62), fontsize=10.5, color=C2, ha='left')
ax.set_xlim(-0.15, 1.75); ax.set_ylim(-0.15, 1.35); ax.set_aspect('equal')
ax.set_xlabel('y'); ax.set_ylabel('(1+c)^y - 1')
ax.set_title('The two corners of the unit square', fontsize=12); ax.legend(frameon=False, fontsize=9, loc='upper left'); clean(ax)
ax = axs[1]
for eps, col in [(1e-3, C1), (1e-4, C2), (1e-5, C3), (1e-6, C6)]:
    b = 2 - eps; o = [2.0]
    for _ in range(80):
        o.append(b ** (o[-1] - 1))
    # repelling fixed point x2 > 2
    lo, hi = 2.0, 2.5
    for _ in range(200):
        mid = (lo + hi) / 2
        if b ** (mid - 1) - mid > 0: hi = mid
        else: lo = mid
    d = np.array([lo - v for v in o])
    j = np.arange(len(d)); keep = d < 0.6
    ax.semilogy(j[keep], d[keep], 'o-', ms=3, lw=1.2, color=col, label=f'1 - c = {eps:g}')
jj = np.arange(0, 60); ax.semilogy(jj, 2e-4 * L4 ** jj, color=INK2, lw=1, ls='--', label='growth by ln 4 each step')
ax.set_ylim(1e-7, 1); ax.set_xlim(0, 44)
ax.set_xlabel('step j'); ax.set_ylabel('distance from the push-out point near (1,1)')
ax.set_title('Escape from the corner: 7 steps per factor of 10', fontsize=12)
ax.legend(frameon=False, fontsize=9.5, loc='lower right'); clean(ax)
plt.tight_layout(); plt.savefig('figures/c6_corner.png'); plt.close()

# ---------- C7: the measured rate ----------
R = json.load(open('rate_scan.json'))
ks = np.array([r['k'] for r in R]); g = np.array([2 - float(r['c']) for r in R])
fig, axs = plt.subplots(1, 2, figsize=(13, 5))
ax = axs[0]
ax.plot(ks[1:], g[1:] / g[:-1], 'o-', color=C1, ms=5, label='gap at rung k / gap at rung k-1')
ax.axhline(1 / L4, color=C2, lw=1.4, ls='--', label='1 / ln 4 = 0.7213')
ax.axvline(19, color=GRID, lw=1); ax.text(19.3, 0.69, 'kiss point jumps\n(3.49 to 3.24)', fontsize=9.5, color=INK2)
ax.set_xlabel('rung k'); ax.set_ylabel('shrink factor'); ax.set_title('Shrink per rung', fontsize=12)
ax.legend(frameon=False, fontsize=10); clean(ax)
ax = axs[1]
ax.plot(ks, g * L4 ** ks, 'o-', color=C6, ms=5)
ax.set_ylim(0, 1.45); ax.set_xlabel('rung k'); ax.set_ylabel('gap x (ln 4)^k')
ax.set_title('Gap times (ln 4)^k stays between two constants', fontsize=12); clean(ax)
plt.tight_layout(); plt.savefig('figures/c7_rate.png'); plt.close()
print('ok')
