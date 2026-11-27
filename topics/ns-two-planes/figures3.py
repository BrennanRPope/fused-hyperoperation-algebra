"""DRAFT. Round 2 pictures (three perpendicular planes). python3 figures3.py -> figures/fig8..fig10"""
import json, numpy as np, matplotlib
matplotlib.use('Agg'); import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
L = lambda f: json.load(open(f))
col = {0: '#922b21', 1: '#117a65', 2: '#1f618d'}; nm = {0: 'plane x (yz)', 1: 'plane y (xz)', 2: 'plane z (xy)'}
lc = {0: '#8e44ad', 1: '#e67e22', 2: '#2c3e50'}
def lines_of(T):
    T = np.array(T); return [T[1, :, 2].sum()+T[2, :, 1].sum(), T[0, :, 2].sum()+T[2, :, 0].sum(), T[0, :, 1].sum()+T[1, :, 0].sum()]
# fig8 concept
fig = plt.figure(figsize=(11, 5.2)); ax = fig.add_subplot(1, 2, 1, projection='3d'); s = 1
P = {2: [[-s, -s, 0], [s, -s, 0], [s, s, 0], [-s, s, 0]], 0: [[0, -s, -s], [0, s, -s], [0, s, s], [0, -s, s]], 1: [[-s, 0, -s], [s, 0, -s], [s, 0, s], [-s, 0, s]]}
for m, q in P.items(): ax.add_collection3d(Poly3DCollection([q], alpha=0.18, facecolor=col[m]))
for a, v in enumerate(np.eye(3)):
    ax.plot(*[[-1.25*v[i], 1.25*v[i]] for i in range(3)], color=lc[a], lw=3.5)
    ax.text(*(1.35*v), f'line {"xyz"[a]}', color=lc[a], fontsize=9)
for m, v in enumerate(np.eye(3)):
    ax.quiver(*(0.55*(np.ones(3)-v)*np.array([1, -1, 1])), *(0.45*v), color=col[m], lw=2, arrow_length_ratio=0.3)
ax.set_box_aspect((1, 1, 1)); ax.view_init(22, -55); ax.set_xlabel('x'); ax.set_ylabel('y'); ax.set_zlabel('z')
ax.set_title('Three perpendicular planes, three crossing lines\n(arrow = each plane\'s own spin, along its normal)', fontsize=10)
ax2 = fig.add_subplot(1, 2, 2); ax2.axis('off')
ax2.text(0, 0.98, "Exact (algebra3.py, checks T1-T4):\n\n"
 "1. Any 3D flow = three flat flows, one per plane,\n   with no plane singled out (Biot-Savart applied\n   to one spin component at a time).\n\n"
 "2. Plane m's own spin is the m-th spin component\n   (plus a gradient that cancels across planes;\n   zero for strict flat flows).\n\n"
 "3. A plane never stretches its own spin. Every\n   stretching term pairs two different planes, so\n   all stretching sits on the three crossing lines.\n\n"
 "4. No flat twist of the cube cycles the three planes\n   (0 of 144 closing screws act freely). The quarter\n   turn swaps two planes; the half-turn space\n   (Hantzsche-Wendt) keeps all three.", va='top', family='monospace', fontsize=8.8)
fig.tight_layout(); fig.savefig('figures/fig8-three-planes-three-lines.png', dpi=140); plt.close(fig)
# fig9 two-plane start: the third plane fills in; the lines share out
fig, axs = plt.subplots(1, 2, figsize=(11, 4.2))
for s in range(3):
    r = L(f'runs/p3_two-un_s{s}.json')['rows']; t = [x['t'] for x in r]
    for i in range(3):
        axs[0].plot(t, [x['comp_Z'][i]/sum(x['comp_Z']) for x in r], color=col[i], lw=1.3, alpha=0.8, label=nm[i] if s == 0 else None)
    Ls = np.array([lines_of(x['T']) for x in r[1:]]); tot = Ls.sum(1)
    for a in range(3):
        axs[1].plot(t[1:], Ls[:, a]/tot, color=lc[a], lw=1.3, alpha=0.8, label=f'line {"xyz"[a]}' if s == 0 else None)
axs[0].axhline(1/3, color='k', lw=0.5, ls=':'); axs[1].axhline(1/3, color='k', lw=0.5, ls=':')
axs[0].set_title('share of total spin held by each plane\n(start: only planes x and z have spin)', fontsize=10)
axs[1].set_title('share of stretching on each crossing line\n(line y = where planes x and z meet)', fontsize=10)
for a in axs: a.set_xlabel('time'); a.legend(fontsize=8)
axs[1].set_ylim(-0.1, 1)
fig.suptitle('Two-plane start (untwisted, 3 seeds, viscosity 0.005, 64³): the third plane fills in to its full third, and the stretching spreads from one line to all three', fontsize=10)
fig.tight_layout(); fig.savefig('figures/fig9-third-plane-fills-in.png', dpi=140); plt.close(fig)
# fig10 two vs three planes, twisted vs untwisted
fig, axs = plt.subplots(1, 2, figsize=(11, 4.2)); modes = ['two-tw', 'two-un', 'three-tw', 'three-un']
lab = {'two-tw': 'two planes\ntwisted', 'two-un': 'two planes\nuntwisted', 'three-tw': 'three planes\ntwisted', 'three-un': 'three planes\nuntwisted'}
cm = {'two-tw': '#c0392b', 'two-un': '#2471a3', 'three-tw': '#e74c3c', 'three-un': '#5dade2'}
for k, m in enumerate(modes):
    zz = []
    for s in range(3):
        r = L(f'runs/p3_{m}_s{s}.json')['rows']; zz.append(max(x['Z'] for x in r)/r[0]['Z'])
        axs[1].plot([x['t'] for x in r], [x['Z']/r[0]['Z'] for x in r], color=cm[m], lw=1.2, alpha=0.8, label=lab[m].replace('\n', ' ') if s == 0 else None)
    axs[0].bar(k, np.mean(zz), color=cm[m], alpha=0.8); axs[0].plot([k]*3, zz, 'ko', ms=4)
axs[0].set_xticks(range(4)); axs[0].set_xticklabels([lab[m] for m in modes], fontsize=8); axs[0].set_ylabel('peak total spin / start')
axs[0].set_title('same starting energy; dots = seeds', fontsize=10); axs[1].legend(fontsize=7); axs[1].set_xlabel('time'); axs[1].set_ylabel('Z(t)/Z(0)')
axs[1].set_title('total spin over time', fontsize=10)
fig.suptitle('Three planes cross more than two: more growth at the same energy. The twist does not help with three planes either', fontsize=10)
fig.tight_layout(); fig.savefig('figures/fig10-two-vs-three-planes.png', dpi=140); plt.close(fig)
print('figures3 done')
