"""DRAFT. Pictures for ns-two-planes. Run after batch1.sh, snapshot run and klein_planes.py.
python3 figures.py  -> figures/fig1..fig7 (png)"""
import json, glob, os, numpy as np, matplotlib
matplotlib.use('Agg'); import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
os.makedirs('figures', exist_ok=True)
L = lambda f: json.load(open(f))
C = {'twisted': '#c0392b', 'untwisted': '#2471a3', 'flat': '#7f8c8d', 'line': '#e67e22'}

# ---------- fig1: the idea and the exact crossing term
fig = plt.figure(figsize=(11, 5.2))
ax = fig.add_subplot(1, 2, 1, projection='3d')
s = 1.0
P1 = [[-s, -s, 0], [s, -s, 0], [s, s, 0], [-s, s, 0]]       # xy-plane
P2 = [[0, -s, -s], [0, s, -s], [0, s, s], [0, -s, s]]       # yz-plane
ax.add_collection3d(Poly3DCollection([P1], alpha=0.25, facecolor='#5dade2'))
ax.add_collection3d(Poly3DCollection([P2], alpha=0.25, facecolor='#f1948a'))
th = np.linspace(0, 2*np.pi, 60)
for cx, cy, sg in [(-0.5, -0.5, 1), (0.5, 0.5, 1), (-0.5, 0.5, -1), (0.5, -0.5, -1)]:
    ax.plot(cx + 0.3*np.cos(th), cy + 0.3*np.sin(th), 0*th, color='#1f618d', lw=1.5)
    ax.quiver(cx, cy, 0, 0, 0, 0.35*sg, color='#1f618d', lw=1.2, arrow_length_ratio=0.3)
for cy, cz, sg in [(-0.5, -0.5, 1), (0.5, 0.5, 1), (-0.5, 0.5, -1), (0.5, -0.5, -1)]:
    ax.plot(0*th, cy + 0.3*np.cos(th), cz + 0.3*np.sin(th), color='#922b21', lw=1.5)
    ax.quiver(0, cy, cz, 0.35*sg, 0, 0, color='#922b21', lw=1.2, arrow_length_ratio=0.3)
ax.plot([0, 0], [-1.2, 1.2], [0, 0], color=C['line'], lw=4)
ax.quiver(0, 0.0, 0, 0, 0.6, 0, color=C['line'], lw=2.5, arrow_length_ratio=0.25)
ax.text(0.05, 1.25, 0.05, 'shared line (y)\nnew spin is made here', color=C['line'], fontsize=9)
ax.text(1.0, -1.0, 0.05, 'flat flow A\non the xy-plane\nspin points along z', color='#1f618d', fontsize=8)
ax.text(0.05, -1.1, 1.05, 'flat flow B\non the yz-plane\nspin points along x', color='#922b21', fontsize=8)
ax.set_xlabel('x'); ax.set_ylabel('y'); ax.set_zlabel('z'); ax.set_box_aspect((1, 1, 1))
ax.set_xlim(-1.2, 1.2); ax.set_ylim(-1.2, 1.2); ax.set_zlim(-1.2, 1.2); ax.view_init(22, -55)
ax.set_title('Two flat flows on crossing planes', fontsize=11)
ax2 = fig.add_subplot(1, 2, 2); ax2.axis('off')
txt = ("Exact (sympy, checks A2-A4):\n\n"
       "1. Each flat flow alone: no stretching at all.\n\n"
       "2. Together, the only new term is the crossing term\n"
       "      (shared-line speed of A) x (change of B along the line)\n"
       "    + (shared-line speed of B) x (change of A along the line).\n"
       "   Nothing else from either plane enters.\n\n"
       "3. All vortex stretching is crossing:\n"
       "   B's spin is stretched by A's strain, A's spin by B's,\n"
       "   and both are tilted onto the shared line, where\n"
       "   neither plane had any spin at the start.\n\n"
       "4. The crossing is harmless (a pure pressure push) only if\n"
       "   both flows are flat along the shared line, i.e. both are\n"
       "   plain shears in the same direction. Otherwise it spins\n"
       "   the fluid out of both planes.")
ax2.text(0.0, 0.98, txt, va='top', family='monospace', fontsize=8.6)
fig.tight_layout(); fig.savefig('figures/fig1-two-planes-and-the-crossing-term.png', dpi=140); plt.close(fig)

# ---------- fig2: Fourier picture: two planes in wavenumber space, a crossing triad leaves both
fig = plt.figure(figsize=(6.2, 5.6)); ax = fig.add_subplot(projection='3d')
s = 3
ax.add_collection3d(Poly3DCollection([[[-s, -s, 0], [s, -s, 0], [s, s, 0], [-s, s, 0]]], alpha=0.22, facecolor='#5dade2'))
ax.add_collection3d(Poly3DCollection([[[0, -s, -s], [0, s, -s], [0, s, s], [0, -s, s]]], alpha=0.22, facecolor='#f1948a'))
k = np.array([2, 1, 0]); p = np.array([0, 1, 2]); q = -(k + p)
for v, col, lab in [(k, '#1f618d', 'k in plane A (kz=0)'), (p, '#922b21', 'p in plane B (kx=0)'), (-q, '#117a65', 'k+p: off both planes')]:
    ax.quiver(0, 0, 0, *v, color=col, lw=2.5, arrow_length_ratio=0.12); ax.text(*(v*1.08), lab, color=col, fontsize=8.5)
ax.plot([0, 0], [-s, s], [0, 0], color=C['line'], lw=3); ax.text(0, s, 0.2, 'ky axis (shared)', color=C['line'], fontsize=8)
ax.set_xlabel('kx'); ax.set_ylabel('ky'); ax.set_zlabel('kz'); ax.set_box_aspect((1, 1, 1)); ax.view_init(20, -60)
ax.set_title('In wavenumbers: two flat flows = two planes.\nTriads inside a plane stay in it (2D);\n'
             'a triad with one leg in each plane lands off both (3D)', fontsize=9.5)
fig.tight_layout(); fig.savefig('figures/fig2-two-planes-in-wavenumber-space.png', dpi=140); plt.close(fig)

# ---------- fig3: does the crossing stay tame? enstrophy and planarity over time
fig, axs = plt.subplots(2, 2, figsize=(11, 7.4), sharex=True)
for j, nu in enumerate(['0.01', '0.005']):
    for s in range(3):
        rA = L(f'runs/n64_nu{nu}_A_k1_s{s}.json')['rows']
        for m in ['twisted', 'untwisted']:
            f = f'runs/n64_nu{nu}_{m}_k1_s{s}.json'
            if not os.path.exists(f): continue
            r = L(f)['rows']; t = [x['t'] for x in r]; Z0 = r[0]['Z']
            axs[0, j].plot(t, [x['Z']/Z0 for x in r], color=C[m], lw=1.4, alpha=0.85, label=m + ' (3 seeds)' if s == 0 else None)
            axs[1, j].plot(t, [x['off_frac'] for x in r], color=C[m], lw=1.4, alpha=0.85)
        # 2D-only reference: plane A alone (spin of a flat flow never grows)
        axs[0, j].plot([x['t'] for x in rA], [x['Z']/rA[0]['Z'] for x in rA], color=C['flat'], ls='--', lw=1.2,
                       label='one flat flow alone (2D)' if s == 0 else None)
    axs[0, j].axhline(1, color='k', lw=0.5)
    axs[0, j].set_title(f'viscosity {nu}: total spin (enstrophy), relative to start', fontsize=10)
    axs[1, j].set_title('share of energy that is no longer in two flat flows', fontsize=10)
    axs[1, j].set_xlabel('time'); axs[1, j].set_ylim(0, 1)
axs[0, 0].legend(fontsize=8); axs[0, 0].set_ylabel('Z(t)/Z(0)'); axs[1, 0].set_ylabel('off-plane energy share')
fig.suptitle('3D Navier-Stokes started as two crossing flat flows (64³): the crossing grows spin that 2D never does,\n'
             'and within one or two turnovers most of the energy has left the two planes. Twisted vs untwisted: no clear difference.', fontsize=10)
fig.tight_layout(); fig.savefig('figures/fig3-does-the-crossing-stay-tame.png', dpi=140); plt.close(fig)

# ---------- fig4: the k dial on the crossing
ks = [0, 0.1, 0.25, 0.5, 1]
fig, axs = plt.subplots(1, 2, figsize=(11, 4.2))
for m in ['twisted', 'untwisted']:
    zz, oo, kk = [], [], []
    for kv in ks:
        f = f'runs/n64_nu0.005_{m}_k{kv if kv != 1 else 1}_s0.json'
        if not os.path.exists(f): continue
        r = L(f)['rows']; kk.append(kv); zz.append(max(x['Z'] for x in r)/r[0]['Z']); oo.append(r[-1]['off_frac'])
        axs[1].plot([x['t'] for x in r], [x['Z']/r[0]['Z'] for x in r], color=C[m], alpha=0.25 + 0.75*kv, lw=1.3)
    axs[0].plot(kk, zz, 'o-', color=C[m], label=m)
axs[0].set_xlabel('k (weight on the crossing term)'); axs[0].set_ylabel('peak enstrophy / start')
axs[0].set_title('a o_k b = a + b + k ab applied to the two flows', fontsize=10); axs[0].legend(fontsize=8); axs[0].axhline(1, color='k', lw=0.5)
axs[1].set_xlabel('time'); axs[1].set_ylabel('Z(t)/Z(0)'); axs[1].set_title('curves for k = 0, 0.1, 0.25, 0.5, 1 (darker = larger k)', fontsize=10)
fig.suptitle('Turning the crossing down: at k = 0 the planes stay flat and spin never grows; any k > 0 lets it grow (viscosity 0.005, 64³)', fontsize=10)
fig.tight_layout(); fig.savefig('figures/fig4-k-dial-on-the-crossing.png', dpi=140); plt.close(fig)

# ---------- fig5: where the planes interact (snapshots of shared-line spin)
f = 'runs/snap_twisted_snap.npz'
if os.path.exists(f):
    S = np.load(f); keys = sorted(S.files, key=lambda k: float(k[1:]))[:4]
    fig, axs = plt.subplots(2, len(keys), figsize=(3.2*len(keys), 6.4))
    for i, kk in enumerate(keys):
        w = S[kk]; N = w.shape[1]
        wy = w[1]; wn = np.sqrt((w**2).sum(0))
        im = axs[0, i].imshow(wn[:, N//3, :].T, origin='lower', cmap='magma', extent=[0, 2*np.pi, 0, 2*np.pi])
        axs[0, i].set_title(f't = {float(kk[1:]):.1f}: all spin |w|', fontsize=9); axs[0, i].set_xlabel('x'); axs[0, i].set_ylabel('z')
        vm = max(1e-6, np.abs(wy).max())
        axs[1, i].imshow(wy[:, N//3, :].T, origin='lower', cmap='RdBu_r', vmin=-vm, vmax=vm, extent=[0, 2*np.pi, 0, 2*np.pi])
        axs[1, i].set_title(f'spin along the shared line w_y\n(max {vm:.1f}; zero at start)', fontsize=9); axs[1, i].set_xlabel('x'); axs[1, i].set_ylabel('z')
    fig.suptitle('Where the two planes interact: a slice across the shared line (y fixed). Top: total spin. Bottom: the new spin\n'
                 'along the shared line, which only the crossing can make (twisted run, viscosity 0.005, 64³)', fontsize=10)
    fig.tight_layout(); fig.savefig('figures/fig5-where-the-planes-interact.png', dpi=140); plt.close(fig)

# ---------- fig6: Klein antiparallel pair through the two-plane lens
f = 'runs/klein_planes.json'
if os.path.exists(f):
    r = L(f)['rows']; t = [x['t'] for x in r]
    fig, axs = plt.subplots(1, 3, figsize=(13, 4))
    axs[0].plot(t, [x['Z'] for x in r], color='k'); axs[0].set_title('total spin (enstrophy)', fontsize=10)
    axs[1].plot(t, [x['plane'] for x in r], color='#922b21'); axs[1].set_ylim(0, 1)
    axs[1].set_title('share of spin in the flat part\n(straight pair = flat flow on the yz-plane)', fontsize=10)
    keys = list(r[0]['P'].keys())
    lab = {'F': 'flat', 'D': 'bend'}
    tot = np.array([sum(x['P'].values()) for x in r])
    axs[2].plot(t, tot, color='k', lw=2, label='all stretching')
    axs[2].plot(t, [x['P']['DDF'] for x in r], color='#c0392b', lw=1.6, label="crossing feed: the flat pair's strain stretches the bend")
    axs[2].plot(t, [x['P']['DDD'] for x in r], color='#2471a3', lw=1.6, label='the bend stretching itself')
    axs[2].plot(t, [sum(v for k, v in x['P'].items() if k not in ('DDF', 'DDD')) for x in r], color='#7f8c8d', lw=1.2, label='everything else (flat part alone: exactly 0)')
    tw = t[int(np.argmax([x['wmax'] for x in r]))]
    for a in axs: a.axvline(tw, color='#e67e22', ls=':', lw=1)
    axs[0].text(tw, axs[0].get_ylim()[0], ' peak spin', color='#e67e22', fontsize=8, va='bottom')
    axs[2].axhline(0, color='k', lw=0.5); axs[2].legend(fontsize=7); axs[2].set_title('stretching split by flat part / bend', fontsize=10)
    for a in axs: a.set_xlabel('time')
    fig.suptitle(f"Klein-class antiparallel pair (N={L(f)['N']}, viscosity {L(f)['nu']}): the crossing feed peaks just before peak spin, then runs down as the flat pair is used up", fontsize=10)
    fig.tight_layout(); fig.savefig('figures/fig6-klein-pair-through-two-planes.png', dpi=140); plt.close(fig)

# ---------- fig7: grid check
fig, ax = plt.subplots(figsize=(6.5, 4))
for m in ['twisted', 'untwisted']:
    for N, ls in [(64, '-'), (96, '--')]:
        f = f'runs/n{N}_nu0.005_{m}_k1_s0.json'
        if os.path.exists(f):
            r = L(f)['rows']; ax.plot([x['t'] for x in r], [x['Z'] for x in r], ls, color=C[m], label=f'{m}, {N}³')
ax.set_xlabel('time'); ax.set_ylabel('enstrophy'); ax.legend(fontsize=8)
ax.set_title('Grid check, viscosity 0.005, twisted seed 0: 64³ vs 96³', fontsize=10)
fig.tight_layout(); fig.savefig('figures/fig7-grid-check.png', dpi=140); plt.close(fig)
print('figures done')
