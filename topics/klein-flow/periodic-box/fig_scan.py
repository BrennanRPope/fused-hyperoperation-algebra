import json, glob, numpy as np, matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from analyse import R, hs, As, seeds, get, wall
cK, cT = 'tab:purple', 'tab:orange'

# ---- figure C: enstrophy gain vs 3D amount, one panel per height, and the wall vs height ----
fig, axs = plt.subplots(1, len(hs) + 1, figsize=(16, 4.2))
for ax, h in zip(axs, hs):
    for s, c, lab in (('K', cK, 'Klein bottle x circle'), ('TR', cT, 'torus x circle (same size)')):
        G = np.array([get(s, h, a, 'G') for a in As])
        for j in range(G.shape[1]):
            ax.plot(As, G[:, j], color=c, alpha=0.25, lw=1)
        ax.plot(As, G.mean(1), color=c, lw=2.5, marker='o', label=lab)
    ax.axhline(1.0, color='0.5', lw=0.8, ls=':')
    ax.set_xscale('log'); ax.set_ylim(0.95, 1.95)
    ax.set_title(f'height h = {h:g}' + ('  (thin slab)' if h < 0.1 else '  (cube)' if h == 1 else ''), fontsize=10)
    ax.set_xlabel('3D energy / 2D energy at start')
axs[0].set_ylabel('enstrophy gain  max Z(t) / Z(0)')
axs[0].legend(fontsize=8, loc='upper left')
ax = axs[-1]
for s, c, off in (('K', cK, -0.03), ('TR', cT, 0.03)):
    for h in hs:
        w = np.array([wall(s, h, sd) for sd in seeds])
        never = np.isnan(w)
        ax.scatter(np.full((~never).sum(), h*(1+off)), w[~never], color=c, s=25)
        if never.any():
            ax.scatter([h*(1+off)], [2.0], marker='^', color=c, s=40)
    med = []
    for h in hs:
        w = np.array([wall(s, h, sd) for sd in seeds])
        med.append(np.nanmedian(w) if np.isnan(w).sum() < len(w)/2 else np.nan)
    ax.plot(hs, med, color=c, lw=2)
ax.set_xscale('log'); ax.set_yscale('log'); ax.set_ylim(0.01, 3)
ax.set_xlabel('box height h'); ax.set_ylabel('wall: 3D amount where enstrophy first rises')
ax.set_title('the wall (triangles: never, up to 1)', fontsize=10)
fig.suptitle('Figure C. Turning the height dial from a thin slab to the full cube, with the same strong 2D flow plus a 3D part. '
             'Below the wall the flow stays in control (enstrophy only falls, as in 2D).\nThin slabs never lose control even with as much 3D energy as 2D. '
             'The Klein bottle (purple) and the same-size torus (orange) hit the wall at the same place at every height.\n'
             '(The amount is measured in energy; a 3D part squeezed into a slab carries more spin per unit of energy, which is why the wall is lowest at h = 1/4.)',
             fontsize=10, x=0.01, ha='left')
fig.tight_layout(rect=(0, 0, 1, 0.83)); fig.savefig('figures/figC-height-dial-wall.png', dpi=120); plt.close(fig)

# ---- figure D: how much the 3D part grows (small seeds), by height ----
fig, axs = plt.subplots(1, 2, figsize=(12, 4.4))
ax = axs[0]
for s, c, lab in (('K', cK, 'Klein bottle x circle'), ('TR', cT, 'torus x circle')):
    for h in hs:
        rs = [r for r in R if r['surf'] == s and r['h'] == h and r['a'] == 0.001]
        g = [r['E3'][-1]/r['E3'][0] for r in rs]
        ax.scatter(np.full(len(g), h*(0.96 if s == 'K' else 1.04)), g, color=c, s=25, label=lab if h == hs[0] else None)
ax.set_xscale('log'); ax.set_yscale('log'); ax.axhline(1, color='0.5', ls=':', lw=0.8)
ax.set_xlabel('box height h'); ax.set_ylabel('3D share of the energy: end / start')
ax.set_title('(a) a tiny 3D seed (0.1% of the energy): does the third direction get in?', fontsize=9.5, loc='left')
ax.legend(fontsize=8)
ax = axs[1]
for s, c in (('K', cK), ('TR', cT)):
    for r in [r for r in R if r['surf'] == s and r['h'] == 1.0 and r['a'] == 0.001]:
        ax.plot(r['t'], np.array(r['E3'])/r['E3'][0], color=c, alpha=0.8)
ax.set_yscale('log'); ax.set_xlabel('time'); ax.set_ylabel('3D share / its start')
ax.set_title('(b) the full cube, four seeds each', fontsize=9.5, loc='left')
fig.suptitle('Figure D. In slabs up to height 1/8 the 3D part dies (left), as the thin-domain theorems lead you to expect. From quarter height up it grows, '
             'and it grows 1.4 to 6 times more on the Klein bottle (12 of 12 matched pairs).\nBut each surface here starts from its own 2D flow, so part of the gap may come from the 2D flow, '
             'not the twist (Figure E separates the two).', fontsize=10, x=0.01, ha='left')
fig.tight_layout(rect=(0, 0, 1, 0.86)); fig.savefig('figures/figD-third-direction-gets-in.png', dpi=120); plt.close(fig)
print('ok')
