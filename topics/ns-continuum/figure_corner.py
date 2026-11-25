import re, glob, json, numpy as np, matplotlib
matplotlib.use('Agg'); import matplotlib.pyplot as plt
plt.rcParams.update({'font.size': 11, 'axes.spines.top': False, 'axes.spines.right': False})
cells = {}
for f in ['corner_b.log', 'cf1.log', 'cf2.log', 'cf3.log', 'cf4.log']:
    try: lines = open(f).read().splitlines()
    except FileNotFoundError: continue
    for ln in lines:
        m = re.match(r'eps=(\S+) alpha=(\S+): (.*)', ln)
        if not m: continue
        eps, a = m.group(1), float(m.group(2))
        runs = re.findall(r'N=(\d+) (REACH|stop) t=([\d.]+)', m.group(3))
        cells.setdefault((eps, a), []).extend((int(N), r == 'REACH', float(t)) for N, r, t in runs)
json.dump({f'{k[0]}|{k[1]}': v for k, v in cells.items()}, open('corner_cells.json', 'w'), indent=0)
epss = ['inf', '4.0', '1.0', '0.25', '0.1', '0.05', '0.0']
labels = {'inf': 'only cross-handed', '4.0': '4', '1.0': '1', '0.25': '0.25', '0.1': '0.1', '0.05': '0.05', '0.0': '0 (none)'}
alphas = [0.2, 0.25, 0.3, 0.4, 0.5, 0.75, 1.0]
fig, ax = plt.subplots(figsize=(7.6, 4.6))
for i, e in enumerate(epss):
    for j, a in enumerate(alphas):
        v = cells.get((e, a))
        if not v: ax.text(j, i, '·', ha='center', va='center', color='#a0aec0'); continue
        allr = all(r for _, r, _ in v); none = not any(r for _, r, _ in v)
        col = '#c53030' if allr else ('#2b6cb0' if none else '#dd6b20')
        ax.add_patch(plt.Rectangle((j - .42, i - .38), .84, .76, color=col, alpha=.85))
        tmin = min(t for _, r, t in v if r) if allr else None
        ax.text(j, i, f'reach\nt≈{tmin:.1f}' if allr else ('held' if none else 'mixed'), ha='center', va='center', color='white', fontsize=8.5)
ax.set_xticks(range(len(alphas))); ax.set_xticklabels([str(a) for a in alphas])
ax.set_yticks(range(len(epss))); ax.set_yticklabels([labels[e] for e in epss])
ax.set_xlim(-.6, len(alphas) - .4); ax.set_ylim(len(epss) - .5, -.5)
ax.set_xlabel('dissipation power α  (ν k^{2α}, ν = 0.001)'); ax.set_ylabel('cross-handed corner weight ε')
ax.set_title('Corner dial: the coupling sets when, not whether')
fig.tight_layout(); fig.savefig('figures/f6_corner.png', dpi=150); plt.close(fig)
for k in sorted(cells): print(k, cells[k])
