"""figG: peak spin and enstrophy for the antiparallel pair, base alone vs a small seed in each class.
usage: python3 fig_classes.py TAG seed out.png"""
import sys, json, glob, numpy as np, matplotlib
matplotlib.use('Agg'); import matplotlib.pyplot as plt
tag, seed, out = sys.argv[1], sys.argv[2], sys.argv[3]
lab = {'none': 'base pair alone', 't-g+': 'Klein only (period doubled)', 't+g-': 'torus only (mirror broken)',
       't+g+': 'both (Kerr period, mirror)', 't-g-': 'neither (plain box only)'}
col = {'none': '#555555', 't-g+': '#c0392b', 't+g-': '#2471a3', 't+g+': '#7d3c98', 't-g-': '#b7950b'}
fig, ax = plt.subplots(1, 2, figsize=(11, 4.2))
for c in ['none', 't+g+', 't-g+', 't+g-', 't-g-']:
    f = 'runs/%s_%s_s%s.json' % (tag, c, seed)
    try: r = json.load(open(f))
    except FileNotFoundError: continue
    t = np.array(r['t'])
    ax[0].plot(t, np.array(r['wmax'])/r['W0'], color=col[c], lw=2, label=lab[c])
    ax[1].plot(t, np.array(r['Z'])/r['Z0'], color=col[c], lw=2, label=lab[c])
ax[0].set_xlabel('time'); ax[0].set_ylabel('peak spin  max|ω| / start')
ax[1].set_xlabel('time'); ax[1].set_ylabel('total spin  enstrophy / start')
ax[0].set_title('How hard the tubes pinch'); ax[1].set_title('How much spin the flow makes')
ax[1].legend(frameon=False, fontsize=9)
for a in ax: a.grid(alpha=0.25)
fig.tight_layout(); fig.savefig(out, dpi=140)
