"""DRAFT. Pictures for the hyperspherical and glued-spheres sections. Orange = known, blue = ours, red = does not work."""
import json, os, sys, numpy as np
H = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, H); sys.argv = sys.argv[:1] + sys.argv[1:]
import figures as F
from figures import plt, OR, BL, RD, GR, save, traj, shape
from tbw import raw
from hyper import rf, runs
Z = json.load(open(os.path.join(H, 'hyper.json'))); FIT = json.load(open(os.path.join(H, 'hyper_fit.json')))
def coll():
    out = []
    for X in ([[0, 0], [1e-7, 0], [1, 0]], [[1, 0], [0, 0], [1e-7, 0]], [[0, 0], [1, 0], [1e-7, 0]]):
        Q = np.zeros((1, 12)); Q[0, :6] = np.ravel(X); out.append(shape(Q)[0])
    return out
def fig9():
    O = Z['ours']; want = ['I.A1', 'I.A2', 'I.B1', 'I.B2', 'I.A10', None]; picks = []
    for m in want:
        for i, o in enumerate(O):
            if o.get('match') == m and i not in picks and len(o['word']) < 40: picks.append(i); break
    f = plt.figure(figsize=(18, 7.4))
    for col, i in enumerate(picks):
        o = O[i]; Xs = traj(o['v1'], o['v2'], o['T'], 3000); S = shape(Xs); r = Xs[:, :6]
        R = np.sqrt((r**2).sum(1)); ts = np.linspace(0, o['T'], len(R))
        a = f.add_subplot(2, 6, col + 1, projection='3d')
        u, vv = np.mgrid[0:2 * np.pi:40j, 0:np.pi:20j]
        a.plot_wireframe(np.cos(u) * np.sin(vv), np.sin(u) * np.sin(vv), np.cos(vv), color='#e6e6e6', lw=.4)
        a.scatter(S[:, 1], S[:, 2], S[:, 0], c=R, cmap='Blues', s=1.2, vmin=R.min() * .7, vmax=R.max())
        for q in coll(): a.scatter(q[1], q[2], q[0], color='k', marker='x', s=50)
        a.set_axis_off(); a.view_init(25, 35)
        w, tl, _ = raw(o['v1'], o['v2'], o['T'] * 0.999999, maxev=400)
        nm = o.get('match') or 'not in Li-Liao list'
        a.set_title(f"{nm}\n{len(w)} crossings, period {o['T']*abs(o['E'])**1.5:.1f}", fontsize=10)
        b = f.add_subplot(2, 6, col + 7); b.plot(ts, R, color=BL, lw=1)
        for t in tl: b.axvline(t, color=GR, lw=.5, alpha=.6)
        b.set_xlabel('time'); b.set_ylabel('size R' if col == 0 else ''); b.tick_params(labelsize=8)
        b.set_title(f"size breathes {o['nRmax']} time{'s' if o['nRmax'] != 1 else ''}, {len(runs(w))} runs", fontsize=9)
    f.suptitle('Six of our polished orbits in hyperspherical coordinates: direction on the shape sphere (top, darker = bigger), size R over one period (bottom, grey = line-ups)', y=1.0)
    save(f, 'fig9-hyperspherical-orbits.png')
def fig10():
    f, ax = plt.subplots(1, 3, figsize=(17, 4.8))
    a = ax[0]; T, tau, arc = FIT['Ts'], FIT['taus'], FIT['arc']; br = FIT['R breaths']['Ts_per_breath']['sd']
    bars = [('ordinary time,\nby letters (2015)', T['letters']['sd'], OR), ('ordinary time,\nby run prices', T['run prices']['sd'], BL),
            ('hyper-angle time\n∫dt/R², best fit', min(tau['run prices']['sd'], tau['number of runs']['sd']), RD),
            ('length of the path\non the shape sphere', arc['run prices']['sd'], RD), ('breaths of\nthe size R', br, RD)]
    a.bar(range(5), [100 * b[1] for b in bars], color=[b[2] for b in bars]); a.set_yscale('log')
    for i, b in enumerate(bars): a.text(i, 100 * b[1] * 1.15, f'{100*b[1]:.2g}%', ha='center', fontsize=9)
    a.set_xticks(range(5)); a.set_xticklabels([b[0] for b in bars], fontsize=8, rotation=30, ha='right'); a.set_ylabel('typical held-out miss (%)')
    a.set_title('A. Which clock does the word predict best?\n(693 Li-Liao orbits)', fontsize=10)
    a = ax[1]; O = Z['ours']; c = np.array(T['run prices']['coef']); xs, y1, y2 = [], [], []
    for o in O:
        Ts = o['T'] * abs(o['E'])**1.5; w = raw(o['v1'], o['v2'], o['T'] * 0.999999, maxev=4000)[0]
        xs.append(Ts); y1.append(np.array(rf(o['word'])) @ c); y2.append(np.array(rf(w)) @ c)
    xs, y1, y2 = map(np.array, (xs, y1, y2)); m = np.max(xs) * 1.05
    a.plot([0, m], [0, m], 'k--', lw=.8)
    a.scatter(xs, y1, s=26, facecolor='none', edgecolor=RD, label=f'reduced word: {np.sum(abs(y1-xs)/xs<.01)} of 66 within 1%')
    a.scatter(xs, y2, s=12, color=BL, label=f'every line-up kept: {np.sum(abs(y2-xs)/xs<.01)} of 66 within 1%')
    a.set_xlabel('true period (common energy)'); a.set_ylabel('run-price period'); a.legend(fontsize=8)
    a.set_title('B. Our 66 polished orbits: count every line-up,\nnot the shortened word', fontsize=10)
    a = ax[2]; per = FIT['stability']['per']; ks = list(per)
    a.bar(range(len(ks)), [per[k]['cv'] for k in ks], color=RD); a.axhline(0, color='k', lw=.5)
    a.set_xticks(range(len(ks))); a.set_xticklabels(ks, fontsize=8, rotation=30, ha='right')
    a.set_ylabel('spread of the instability rate\n(std / mean)'); a.set_ylim(0, 1.4)
    a.set_title(f"C. Instability per unit of each clock: no clock\nmakes it a constant ({per[ks[0]]['n']} unstable orbits)", fontsize=10)
    save(f, 'fig10-hyperspherical-scores.png')
def fig11():
    from matplotlib.patches import Circle, FancyArrowPatch
    f, ax = plt.subplots(1, 4, figsize=(19, 5.2))
    def sph(a, x, y, r, col, al=.25): a.add_patch(Circle((x, y), r, fc=col, ec=col, alpha=al, lw=2)); a.add_patch(Circle((x, y), r, fc='none', ec=col, lw=1.5))
    a = ax[0]
    for k in range(3):
        th = np.pi / 2 + 2 * np.pi * k / 3; sph(a, np.cos(th), np.sin(th), .866, RD)
    a.plot(*zip(*[(1.6 * np.cos(np.pi / 2 + 2 * np.pi * k / 3 + .0), 1.6 * np.sin(np.pi / 2 + 2 * np.pi * k / 3)) for k in range(4)]), alpha=0)
    th = np.linspace(0, 2 * np.pi, 200); a.plot(.55 * np.cos(th), .55 * np.sin(th), color='k', lw=1.5, ls='--')
    a.set_title('1. One sphere per pair, touching in a ring\nloops: just one (Z). Two letters need two.', fontsize=10)
    a = ax[1]
    for k in range(3): sph(a, -1.6 + 1.6 * k, 0, .8, RD)
    a.set_title('2. Three spheres in a chain, or glued on discs\nloops: none. No letters at all.', fontsize=10)
    a = ax[2]; sph(a, 0, 0, 1.6, OR, .12)
    P = [(0, 1.0), (-.8, -.35), (.8, -.35)]; lab = ['bodies 1,2 meet', 'bodies 2,3 meet', 'bodies 1,3 meet']
    for (x, y), l in zip(P, lab):
        a.add_patch(Circle((x, y), .3, fc=OR, ec=OR, alpha=.35)); a.plot(x, y, 'kx', ms=10, mew=2); a.text(x, y - .5 if y < 0 else y + .4, l, ha='center', fontsize=8)
    t = np.linspace(0, 2 * np.pi, 100)
    a.plot(-.8 + .5 * np.cos(t), -.35 + .5 * np.sin(t), color='k', lw=1.4); a.text(-1.5, .2, 'a', fontsize=14, weight='bold')
    a.plot(.8 + .5 * np.cos(t), -.35 + .5 * np.sin(t), color='k', lw=1.4, ls='--'); a.text(1.35, .2, 'b', fontsize=14, weight='bold')
    a.set_title('3. Three collision caps glued into one sphere,\ncollision points removed: two loops a, b, the words', fontsize=10)
    a = ax[3]; sph(a, 0, 0, 1.6, OR, .12)
    for k, c in enumerate(['k', '#555', '#999']):
        a.add_patch(Circle((.35 * np.cos(2 * np.pi * k / 3), .35 * np.sin(2 * np.pi * k / 3)), .8, fc='none', ec=c, lw=2))
    a.set_title('4. The same in 4D: the 3-sphere minus three\nlinked circles (one per collision): words × one turn', fontsize=10)
    for a in ax: a.set_xlim(-2.6, 2.6); a.set_ylim(-2.2, 2.4); a.set_aspect('equal'); a.axis('off')
    save(f, 'fig11-gluing-spheres.png')
if __name__ == '__main__':
    for nm in sys.argv[1:]: globals()[nm]()
