import json, glob, numpy as np, matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
J = sum((json.load(open(f)) for f in sorted(glob.glob('runs/joint*_s*.json'))), [])
C = sum((json.load(open(f)) for f in sorted(glob.glob('runs/classes_s*.json'))), [])
col = {'K': 'tab:purple', 'TR': 'tab:orange', 'T3': '0.4'}
fig, axs = plt.subplots(1, 2, figsize=(14, 4.8))
ax = axs[0]
xs = []; labels = []
for i, (h, a) in enumerate([(0.25, 0.001), (0.25, 0.01), (0.5, 0.001), (0.5, 0.01), (1.0, 0.001), (1.0, 0.01)]):
    for j, s in enumerate(('TR', 'K', 'T3')):
        g = [r['E3abs'][-1]/r['E3abs'][0] for r in J if r['surf'] == s and r['h'] == h and r['a'] == a]
        if not g: continue
        x = i*4 + j
        ax.bar(x, np.mean(g), color=col[s], alpha=0.6, label={'TR': 'torus x circle (straight gluing)', 'K': 'Klein bottle x circle (twist)', 'T3': 'doubled box (all disturbances, mixed)'}[s] if i == 0 or (s == 'T3' and i == 2) else None)
        ax.scatter([x]*len(g), g, color='k', s=10, zorder=3)
    xs.append(i*4 + 1); labels.append(f'h={h:g}\nseed {a:g}')
ax.set_xticks(xs); ax.set_xticklabels(labels, fontsize=8.5); ax.set_yscale('log')
ax.set_ylabel('3D energy at t = 5 / at start')
ax.legend(fontsize=8, loc='upper left')
ax.set_title('(a) Same 2D flow on both, only the gluing differs: the 3D part\ngrows more on the Klein bottle in 24 of 24 matched runs', fontsize=9.5, loc='left')
ax = axs[1]
names = {(1, 1): 'same under shift and glide\n(mirror-symmetric; both keep it)', (-1, 1): 'flips under shift, same under glide\n(period doubled; only the twist keeps it)',
         (1, -1): 'same under shift, flips under glide\n(only the straight gluing keeps it)', (-1, -1): 'flips under both\n(period doubled; neither keeps it)'}
cc = {(1, 1): 'tab:green', (-1, 1): 'tab:purple', (1, -1): 'tab:orange', (-1, -1): '0.5'}
for k, ((ts, gs), nm) in enumerate(names.items()):
    g = [r['growth'] for r in C if r['tau'] == ts and r['g'] == gs]
    if not g: continue
    ax.barh(k, np.mean(g), color=cc[(ts, gs)], alpha=0.6); ax.scatter(g, [k]*len(g), color='k', s=10, zorder=3)
ax.set_yticks(range(4)); ax.set_yticklabels(list(names.values()), fontsize=8.5); ax.set_xscale('log')
from matplotlib.ticker import NullFormatter, LogLocator
ax.xaxis.set_minor_formatter(NullFormatter()); ax.set_xlabel('3D energy at t = 5 / at start (cube, seed 0.1%)')
ax.set_title('(b) Why: sort disturbances by the two gluings (bars: means). The class only the\ntwist keeps (period doubled) beats the class only the straight gluing keeps, every seed', fontsize=9.5, loc='left')
fig.suptitle('Figure E. The twist lets the third direction in faster because it keeps disturbances twice as long along the glide '
             'as a straight gluing of the same box allows.', fontsize=10.5, x=0.01, ha='left')
fig.tight_layout(rect=(0, 0, 1, 0.92)); fig.savefig('figures/figE-why-the-twist-lets-3D-in.png', dpi=125); plt.close(fig)
print('ok')
