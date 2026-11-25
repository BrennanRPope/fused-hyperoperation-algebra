import csv, math, matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from collections import Counter
BLUE, ORANGE, INK, MUTE = '#257de0', '#f47230', '#1d2330', '#8a93a3'
rows = [r for r in csv.DictReader(open('cluster/landscape_census.csv')) if r['k'] != 'family']
cnt = Counter(round(float(r['k']), 6) for r in rows)
fig, ax = plt.subplots(figsize=(9, 3.6), dpi=200)
for k, c in cnt.items():
    ax.bar(k, c, width=0.07, color=BLUE, zorder=3)
    ax.text(k, c + 0.8, str(c), ha='center', va='bottom', fontsize=8, color=INK)
for n in range(5, 40):
    b = 1 + 2 * math.cos(2 * math.pi / n)
    ax.axvline(b, color=ORANGE, lw=1.0 if n <= 12 else 0.4, alpha=0.9 if n <= 12 else 0.4, zorder=1)
for n, lab in [(5, '1.618'), (6, '2'), (7, '2.247'), (8, '2.414')]:
    ax.text(1 + 2 * math.cos(2 * math.pi / n), 38.5, lab, rotation=90, ha='right', va='top', fontsize=7, color=ORANGE)
for b in (-1, 0, 1):
    ax.axvline(b, color=ORANGE, lw=1, ls=(0, (2, 2)), alpha=0.6, zorder=1)
ax.set_xlim(-2.5, 3.15); ax.set_ylim(0, 40)
ax.set_xlabel('dial value k  (placed laws, fixed points only)', fontsize=9, color=INK)
ax.set_ylabel('laws placed', fontsize=9, color=INK)
for s in ('top', 'right'): ax.spines[s].set_visible(False)
for s in ('left', 'bottom'): ax.spines[s].set_color(MUTE)
ax.tick_params(colors=MUTE, labelsize=8)
ax.grid(axis='y', color='#e6e8ec', lw=0.6, zorder=0)
ax.annotate('k = −2: 4 laws, colour count −1,\nwhere no ring ever vanishes', xy=(-2, 4.3), xytext=(-2.4, 16),
            fontsize=8, color=INK, arrowprops=dict(arrowstyle='-', color=MUTE, lw=0.8))
ax.annotate('k = 2: 3 laws, all\nk = 1 in other units', xy=(2, 3.3), xytext=(0.12, 16),
            fontsize=8, color=INK, arrowprops=dict(arrowstyle='-', color=MUTE, lw=0.8))
fig.text(0.99, 0.01, 'blue: laws per dial value   orange: predicted values 1 + 2cos(2π/n), n ≥ 5 (dashed: n = 2, 3, 4)',
         ha='right', fontsize=7, color=MUTE)
plt.tight_layout(rect=(0, 0.04, 1, 1))
plt.savefig('figures/f4_cluster.png'); plt.savefig('page/figures/f4_cluster.png')
