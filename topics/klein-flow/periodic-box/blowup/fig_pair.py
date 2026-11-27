"""figF: the antiparallel pair as one tube plus its glide image, and what each kind of disturbance does.
Top view (x along the tubes, y across). usage: python3 fig_pair.py out.png"""
import sys, numpy as np, matplotlib
matplotlib.use('Agg'); import matplotlib.pyplot as plt
x = np.linspace(0, 2*np.pi, 400); d, a, e = 0.5, 0.25, 0.12
yA = d - a*np.cos(2*x)
fig, ax = plt.subplots(1, 3, figsize=(13, 3.6), sharey=True)
def pair(axi, yA, yB, title, note):
    axi.plot(x, yA, color='#c0392b', lw=3); axi.plot(x, yB, color='#2471a3', lw=3)
    for xx in (0.6, 3.7):
        i = np.searchsorted(x, xx)
        axi.annotate('', (x[i+25], yA[i+25]), (x[i], yA[i]), arrowprops=dict(arrowstyle='->', color='#c0392b', lw=2))
        axi.annotate('', (x[i], yB[i]), (x[i+25], yB[i+25]), arrowprops=dict(arrowstyle='->', color='#2471a3', lw=2))
    axi.axhline(0, color='#aaa', lw=0.8, ls='--'); axi.axvline(np.pi, color='#ddd', lw=0.8)
    axi.set_title(title, fontsize=11); axi.text(0.1, -1.05, note, fontsize=9, color='#444')
    axi.set_xticks([0, np.pi, 2*np.pi]); axi.set_xticklabels(['0', 'π', '2π']); axi.set_xlabel('along the tubes')
    axi.set_ylim(-1.15, 1.0)
# glide image of a centreline y(x) is -y(x - pi)
img = lambda f: -np.interp((x - np.pi) % (2*np.pi), x, f)
pair(ax[0], yA, img(yA), 'The pair: red tube + its glide image', 'spin arrows opposite: the twist makes\nthe partner antiparallel. Kerr bend repeats every π')
yA1 = yA + e*np.cos(x)
pair(ax[1], yA1, img(yA1), 'Klein-only nudge (period doubled)', 'still mirror symmetric across the twist,\nbut no longer repeats every π')
yA2 = yA + e*np.cos(2*x + 1.0); yB2 = img(yA) + e*np.cos(2*x + 1.0)
pair(ax[2], yA2, yB2, 'Torus-only nudge (mirror broken)', 'still repeats every π,\nbut the two tubes no longer mirror')
ax[0].set_ylabel('across (y)')
fig.tight_layout(); fig.savefig(sys.argv[1], dpi=140)
