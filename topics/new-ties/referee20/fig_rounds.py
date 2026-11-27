import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
mine = {'18.1 (3 v 3)': [6.438e45, 422, 42, 33],
        '19.1 (4 v 2)': [6.438e45, 558, 50, 38],
        '19.2 (4 v 3)': [2.775e62, 765, 59, 47, 45],
        '19.3 (4 v 4)': [1.550e79, 1.570e19, 263, 61, 51, 49]}
claimed = {'19.1 (4 v 2)': [2.1e29, 369, 61, 50, 50],
           '19.2 (4 v 3)': [2.1e62, 4.0e18, 269, 68, 59, 55, 55],
           '19.3 (4 v 4)': [4.8e112, 4.0e40, 1.0e23, 371, 74, 66, 61, 59, 59]}
fig, axs = plt.subplots(1, 4, figsize=(15, 3.8), sharey=True)
for ax, k in zip(axs, mine):
    ax.semilogy(range(len(mine[k])), mine[k], 'o-', color='#2a6fb0', label='referee 20 (gap-sum bound)')
    if k in claimed:
        ax.semilogy(range(len(claimed[k])), claimed[k], 's--', color='#c0612b', label='claimed')
    ax.set_title(k); ax.set_xlabel('round (0 = Matveev start)'); ax.grid(alpha=.3)
axs[0].set_ylabel('bound'); axs[3].legend(fontsize=8)
plt.suptitle('Bound per round: start -> lattice/2-adic rounds -> final box for the exact check')
plt.tight_layout(); plt.savefig('fig_rounds.png', dpi=120)
