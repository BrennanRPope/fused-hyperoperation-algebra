"""DRAFT. Pictures for the twin-orbit and period-rule follow-up. Orange = known, blue = ours, red = does not work."""
import json, os, sys, numpy as np
H = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, H); sys.argv = sys.argv[:1] + sys.argv[1:]
import figures as F
from figures import plt, OR, BL, RD, GR, save
import twins as W
PU = '#8e44ad'
def fig12():
    T = json.load(open(os.path.join(H, 'twins.json'))); Fd = json.load(open(os.path.join(H, 'fold.json')))['I.A1']['down']
    picks = ['I.A1', 'II.C192', 'II.C66']; names = {'I.A1': 'figure-eight', 'II.C192': 'II.C192', 'II.C66': 'II.C66'}
    f = plt.figure(figsize=(18, 9))
    for col, c in enumerate(picks):
        r = T[c]
        for row, key in enumerate(('cat', 'twin')):
            p = r[key]; P = W.pos(*p, 6000); a = f.add_subplot(2, 4, col + 1 + 4 * row)
            for k, cc in enumerate(['k', OR, BL] if row == 0 else ['k', PU, BL]): a.plot(P[:, k, 0], P[:, k, 1], color=cc, lw=.8)
            a.set_aspect('equal'); a.axis('off')
            lab = ('Li-Liao ' + names[c]) if row == 0 else ('same word, second orbit' + (' (Simó\'s H orbit, known)' if c == 'I.A1' else ''))
            a.set_title(f"{lab}\nperiod {r['Ts'][row]:.6f}, worst multiplier {r['maxlam'][row]:.3g}", fontsize=10)
    a = f.add_subplot(1, 4, 4)
    A = np.array(Fd['A']); B = np.array(Fd['B'])
    a.plot(A[:, 0], A[:, 1], color=OR, lw=2, label='figure-eight branch'); a.plot(B[:, 0], B[:, 1], color=PU, lw=2, label='second-orbit branch')
    mstar = 1.000045  # measured: both close at m3 = 1.00004, neither at 1.00005 (fine polish sweep)
    a.axvline(mstar, color='k', ls=':', lw=1); a.text(mstar, a.get_ylim()[0] if False else np.min(B[:, 1]), ' they meet and vanish\n at m3 ≈ 1.00004', fontsize=9, va='bottom')
    a.set_xlabel('mass of the middle body, m3 (others 1)'); a.set_ylabel('starting push v1'); a.legend(fontsize=9, loc='upper left')
    a.set_title('Changing one mass: the two orbits\nrun into each other just above equal masses', fontsize=10)
    f.suptitle('One word, two orbits: each pair differs under every symmetry (path, size swings, stability)', y=1.0)
    save(f, 'fig12-one-word-two-orbits.png'); return mstar
def fig13():
    R = json.load(open(os.path.join(H, 'runcost.json'))); pr = json.load(open(os.path.join(H, 'period_runs.json')))['cost per run length']['coeffs']
    allr = [r for o in R for r in o['runs']]
    f, ax = plt.subplots(1, 3, figsize=(18, 4.8))
    a = ax[0]
    for k, c in zip((2, 3, 4, 5), ('#9cc3e6', '#5b9bd5', BL, '#0b3d91')):
        du = [r['dur'] for r in allr if r['k'] == k]; a.hist(du, bins=np.linspace(2, 14, 121), color=c, alpha=.85, label=f'run of {k}: {np.mean(du):.2f} ± {100*np.std(du)/np.mean(du):.1f}%')
        a.axvline(pr[k - 1], color='k', lw=.8, ls='--')
    a.set_yscale('log'); a.set_xlabel('time the run lasts (common energy)'); a.set_ylabel('runs'); a.legend(fontsize=8)
    a.set_title('A. Measured run durations (dashed = fitted prices)', fontsize=10)
    a = ax[1]; ks = (1, 2, 3, 4, 5); bottom = np.zeros(5)
    for nc, c in zip(range(0, 6), ('#eeeeee', '#cfe0f3', '#9cc3e6', '#5b9bd5', BL, '#0b3d91')):
        v = np.array([np.mean([len(r['close']) == nc for r in allr if r['k'] == k]) for k in ks])
        a.bar(ks, v, bottom=bottom, color=c, label=f'{nc} close passes'); bottom += v
    a.set_xlabel('letters in the run'); a.set_ylabel('share of runs'); a.legend(fontsize=8, loc='lower right'); a.set_xticks(ks)
    a.set_title('B. A run of k letters holds k close passes\n(every run of 4, 99.8% of runs of 3)', fontsize=10)
    a = ax[2]; o = [x for x in F.d if x['cls'] == 'II.C1'][0]
    P = W.pos(o['v1'], o['v2'], o['T'], 6000); E = abs(-2.5 + 3 * (o['v1']**2 + o['v2']**2)); ts = np.linspace(0, o['T'], len(P)) * E**1.5
    d = np.stack([np.hypot(*(P[:, i] - P[:, j]).T) for i, j in ((0, 1), (1, 2), (0, 2))], 1)
    for k, c in enumerate(['k', OR, BL]): a.plot(ts, d[:, k] * E, color=c, lw=.9, label=['bodies 1,2', 'bodies 2,3', 'bodies 1,3'][k])
    from tbw import raw
    w, tl, _ = raw(o['v1'], o['v2'], o['T'] * 0.999999, maxev=100)
    for t, ch in zip(tl, w): a.text(t * E**1.5, a.get_ylim()[1] if False else 3.3, ch, ha='center', fontsize=9, color=BL if ch.islower() else RD)
    a.set_ylim(0, 3.6); a.set_xlabel('time (common energy)'); a.set_ylabel('pair distance (common energy)'); a.legend(fontsize=8, loc='upper right')
    a.set_title(f'C. Yin-yang orbit {w}: each letter comes with\none close pass; capitals red, small letters blue', fontsize=10)
    save(f, 'fig13-why-runs-cost-fixed-time.png')
def fig14():
    X = json.load(open(os.path.join(H, 'period_ext.json'))); pr = json.load(open(os.path.join(H, 'period_runs.json')))
    rows = [('m3 = 0.5', X['unequal']['0.5']), ('m3 = 0.75', X['unequal']['0.75']), ('equal masses\n(Li-Liao 695)', None), ('m3 = 2', X['unequal']['2.0']), ('m3 = 4', X['unequal']['4.0']),
            ('free fall,\nequal masses', X['free_fall']['(1.0, 1.0, 1.0)']), ('free fall,\n1, 0.8, 0.6', X['free_fall']['(1.0, 0.8, 0.6)'])]
    f, ax = plt.subplots(figsize=(13, 4.8)); w = .38
    for i, (lab, r) in enumerate(rows):
        if r is None: L, Rr = pr['letters only']['heldout_rel_sd'], pr['cost per run length']['heldout_rel_sd']; n = 695
        else: L, Rr = r['letters']['sd'], r['refit']['sd']; n = r['n']
        ax.bar(i - w / 2, 100 * L, w, color=OR); ax.bar(i + w / 2, 100 * Rr, w, color=BL if Rr < L else RD)
        ax.text(i, max(L, Rr) * 100 * 1.25, f'n={n}', ha='center', fontsize=8)
    ax.set_yscale('log'); ax.set_xticks(range(len(rows))); ax.set_xticklabels([r[0] for r in rows], fontsize=9); ax.set_ylabel('typical held-out miss in the period (%)')
    from matplotlib.patches import Patch
    ax.legend(handles=[Patch(color=OR, label='counting letters (known law, refitted per group)'), Patch(color=BL, label='run prices, refitted per group: better'), Patch(color=RD, label='run prices: no better')], fontsize=9, loc='upper left')
    ax.set_title('The run rule outside its fit: sharp at equal masses, still ahead at m3 = 0.75, no help elsewhere')
    save(f, 'fig14-run-rule-outside-the-fit.png')
if __name__ == '__main__':
    for nm in sys.argv[1:]: print(globals()[nm]())
