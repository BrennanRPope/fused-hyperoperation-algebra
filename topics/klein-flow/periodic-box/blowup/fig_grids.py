"""figI: peak spin and total spin for the pair alone at viscosity 0.001 on three grids (128, 192, 256).
usage: python3 fig_grids.py out.png"""
import sys, json, numpy as np, matplotlib
matplotlib.use('Agg'); import matplotlib.pyplot as plt
runs = [(128, 'nu0.001_none_s0.json', '#f0b27a'), (192, 'res192_none_s0.json', '#cb4335'), (256, 'res256_none_s0.json', '#1b2631')]
fig, ax = plt.subplots(1, 3, figsize=(14, 4.2))
peaks = []
for N, f, c in runs:
    try: r = json.load(open('runs/' + f))
    except FileNotFoundError: continue
    t = np.array(r['t']); m = t <= 1.9
    ax[0].plot(t[m], np.array(r['wmax'])[m], color=c, lw=2.2, label='%d³ grid' % N)
    ax[1].plot(t[m], np.array(r['Z'])[m]/r['Z0'], color=c, lw=2.2, label='%d³ grid' % N)
    peaks.append((N, max(np.array(r['wmax'])[m])))
ax[0].set_title('Peak spin, viscosity 0.001'); ax[0].set_xlabel('time'); ax[0].set_ylabel('max|ω| on the grid')
ax[1].set_title('Total spin, same runs'); ax[1].set_xlabel('time'); ax[1].set_ylabel('enstrophy / start')
if peaks:
    Ns, P = zip(*peaks)
    ax[2].plot(Ns, P, 'o-', color='#1b2631', lw=2, ms=8)
    for n, p in peaks: ax[2].annotate('%.0f' % p, (n, p), textcoords='offset points', xytext=(6, -12))
    ax[2].set_xticks(Ns); ax[2].set_xlabel('grid points per side'); ax[2].set_ylabel('highest peak spin')
    ax[2].set_title('Does the peak settle as the grid gets finer?')
    ax[2].set_ylim(0, max(P)*1.25)
for a in ax: a.grid(alpha=0.25)
ax[0].legend(frameon=False)
fig.tight_layout(); fig.savefig(sys.argv[1], dpi=140)
