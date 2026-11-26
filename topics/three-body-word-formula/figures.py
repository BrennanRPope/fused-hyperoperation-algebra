"""DRAFT. Pictures for the word -> orbit page. Colour key: orange = known, blue = ours, red = does not work."""
import json, os, sys, numpy as np
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
H = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, H)
from tbw import ic, flow, raw
FIG = os.path.join(H, 'figures'); os.makedirs(FIG, exist_ok=True)
OR, BL, RD, GR = '#e07b00', '#1f6fd1', '#d0312d', '#888888'
plt.rcParams.update({'font.size': 11, 'axes.titlesize': 12, 'figure.dpi': 110})
d = json.load(open(os.path.join(H, 'data', 'li-liao-695.json')))
Lf = np.array([o['Lf'] for o in d], float); Ts = np.array([o['Tstar'] for o in d])
v1 = np.array([o['v1'] for o in d]); v2 = np.array([o['v2'] for o in d])
def inverse(w): return w[::-1].swapcase()
def pairs(w):
    c = w + w[0]; a = b = e = 0
    for x, y in zip(c, c[1:]):
        if x == y: a += 1
        elif x.isupper() == y.isupper(): b += 1
        else: e += 1
    return a, b, e
P = np.array([pairs(o['word']) for o in d], float)
def save(f, name): f.savefig(os.path.join(FIG, name), bbox_inches='tight'); plt.close(f); print('wrote', name)
def traj(p1, p2, T, n=1500):
    X = ic(p1, p2); out = [X]; dt = T / n
    for _ in range(n): X = flow(X, dt); out.append(X)
    return np.array(out)
def shape(Xs):
    r = Xs[:, :6:2] + 1j * Xs[:, 1:6:2]
    z1 = (r[:, 1] - r[:, 0]) / np.sqrt(2); z2 = np.sqrt(2 / 3) * (r[:, 2] - 0.5 * (r[:, 0] + r[:, 1]))
    R2 = abs(z1)**2 + abs(z2)**2; w = np.conj(z1) * z2
    return np.stack([(abs(z1)**2 - abs(z2)**2) / R2, 2 * w.real / R2, 2 * w.imag / R2], 1)

def fig1():
    pr = json.load(open(os.path.join(H, 'period_runs.json')))
    c1 = pr['letters only']['coeffs'][0]; cr = pr['cost per run length']['coeffs']
    f, ax = plt.subplots(1, 3, figsize=(16.5, 4.9))
    ax[0].scatter(Lf, Ts, s=8, c=OR, alpha=.7, label='695 orbits (Li & Liao 2017)')
    x = np.linspace(0, Lf.max(), 10); ax[0].plot(x, c1 * x, 'k--', lw=1, label=f'period = {c1:.2f} × letters')
    ax[0].set_xlabel('letters in the word'); ax[0].set_ylabel('period at a common energy, T* = T|E|^1.5')
    ax[0].set_title('A. Longer word, longer period (known, 2015)'); ax[0].legend(loc='upper left', fontsize=9)
    ks = np.arange(1, 6); cnt = [pr['run_length_counts'][str(k)] for k in ks]
    ax[1].bar(ks, cr, color=BL, width=.6)
    for k, c, n in zip(ks, cr, cnt): ax[1].text(k, c + .25, f'{c:.2f}\n({c/k:.2f} per letter)\nseen {n}×', ha='center', fontsize=8)
    ax[1].set_xticks(ks); ax[1].set_xticklabels(['a', 'ab', 'abA…\n3', 'abAB…\n4', '5'][:0] or ['1', '2', '3', '4', '5'])
    ax[1].set_xlabel('length of a run (letters turning the same way before a reversal)'); ax[1].set_ylabel('time the run costs')
    ax[1].set_ylim(0, 15.5); ax[1].set_title('B. Each run of letters has a fixed price (ours)')
    bins = np.linspace(-6, 6, 97)
    for nm, col, lab in [('letters only', GR, 'letters only'), ('letters + reversals', OR, 'letters + reversals'), ('cost per run length', BL, 'price per run')]:
        rel = np.array(pr[nm]['rel']) * 100
        ax[2].hist(np.clip(rel, -6, 6), bins=bins, color=col, alpha=.6 if col != BL else .85, label=f"{lab}: ±{np.std(rel):.2f}%")
    ax[2].set_xlabel('miss in period, held-out orbits (%)'); ax[2].set_ylabel('orbits'); ax[2].set_yscale('log')
    ax[2].set_title('C. Misses shrink about 12 times'); ax[2].legend(fontsize=9)
    save(f, 'fig1-period-from-the-word.png')

def run_of_each_letter(w):
    """length of the cyclic same-case run each letter belongs to"""
    n = len(w); lab = [0] * n
    i0 = next((i for i in range(n) if w[i].isupper() != w[i - 1].isupper()), 0)
    i = 0
    while i < n:
        j = i
        while j + 1 < n and w[(i0 + j + 1) % n].isupper() == w[(i0 + i) % n].isupper(): j += 1
        for k in range(i, j + 1): lab[(i0 + k) % n] = j - i + 1
        i = j + 1
    return lab
RUNCOL = {1: '#999999', 2: BL, 3: '#8e44ad', 4: OR, 5: '#16a085'}
def fig2():
    pr = json.load(open(os.path.join(H, 'period_runs.json'))); cr = pr['cost per run length']['coeffs']
    picks = ['I.A1', 'I.A2', 'II.C1']; names = {'I.A1': 'figure-eight', 'I.A2': 'butterfly I', 'II.C1': 'yin-yang (II.C1)'}
    f = plt.figure(figsize=(16, 9.6))
    for col, cls in enumerate(picks):
        o = [x for x in d if x['cls'] == cls][0]
        Xs = traj(o['v1'], o['v2'], o['T'], 3000); S = shape(Xs); ts = np.linspace(0, o['T'], len(S))
        w, t, _ = raw(o['v1'], o['v2'], o['T'] * 0.99999, maxev=200); lab = run_of_each_letter(w)
        a1 = f.add_subplot(2, 3, col + 1); r = Xs[:, :6].reshape(-1, 3, 2)
        for k, c in enumerate(['k', '#777', '#bbb']): a1.plot(r[:, k, 0], r[:, k, 1], color=c, lw=1)
        a1.set_aspect('equal'); a1.axis('off'); a1.set_title(f"{names[cls]}: word {w}", fontsize=10)
        a2 = f.add_subplot(2, 3, col + 4, projection='3d')
        u, vv = np.mgrid[0:2 * np.pi:40j, 0:np.pi:20j]
        a2.plot_wireframe(np.cos(u) * np.sin(vv), np.sin(u) * np.sin(vv), np.cos(vv), color='#e4e4e4', lw=.4)
        edges = np.concatenate([[0], t, [o['T']]])
        for i in range(len(edges) - 1):
            m = (ts >= edges[i]) & (ts <= edges[i + 1]); li = min(i, len(lab) - 1)
            a2.plot(S[m, 1], S[m, 2], S[m, 0], color=RUNCOL[lab[li]], lw=2)
        for X in ([[0, 0], [1e-7, 0], [1, 0]], [[1, 0], [0, 0], [1e-7, 0]], [[0, 0], [1, 0], [1e-7, 0]]):
            Z = np.zeros((1, 12)); Z[0, :6] = np.ravel(X); q = shape(Z)[0]; a2.scatter(q[1], q[2], q[0], color='k', marker='x', s=60)
        a2.set_axis_off(); a2.view_init(25, 35)
        from collections import Counter
        rr = Counter(); i = 0
        rl = [int(x) for x in __import__('predict2').runs(o['word'])]
        price = sum(cr[k - 1] for k in rl)
        a2.set_title('runs: ' + ' + '.join(str(k) for k in rl) + f'\nprice {price:.2f}   true period {o["Tstar"]:.2f}', fontsize=10)
    from matplotlib.lines import Line2D
    f.legend(handles=[Line2D([0], [0], color=RUNCOL[k], lw=3, label=f'run of {k} (price {cr[k-1]:.2f})') for k in (2, 3, 4)],
             loc='lower center', ncol=3, fontsize=10)
    f.suptitle('On the shape sphere (black × = collisions), each run of same-way letters gets one colour', fontsize=12)
    save(f, 'fig2-runs-on-the-shape-sphere.png')

def fig3():
    pl = json.load(open(os.path.join(H, 'period_law.json')))
    f, ax = plt.subplots(1, 2, figsize=(13, 5.2))
    ax[0].scatter(v1, v2, s=5, c=OR); th = np.linspace(0, np.pi / 2, 100); ax[0].plot(np.sqrt(5 / 6) * np.cos(th), np.sqrt(5 / 6) * np.sin(th), 'k:', lw=1)
    ax[0].set_xlabel('starting push v₁'); ax[0].set_ylabel('starting push v₂'); ax[0].set_aspect('equal')
    ax[0].set_title('A. Where the 695 orbits start (dotted: escape)')
    ins = ax[0].inset_axes([0.62, 0.62, 0.36, 0.36]); m = (abs(v1 - 0.4) < 0.03) & (abs(v2 - 0.25) < 0.03)
    ins.scatter(v1[m], v2[m], s=10, c=OR); ins.set_xticks([]); ins.set_yticks([]); ins.set_title('zoom ×15', fontsize=8)
    cv = pl['feature_cv']
    names = ['guess the average', 'straight-line fit\non letter counts', 'nearest look-alike\nwords']
    vals = [cv['mean']['median_dist'], cv['linear']['median_dist'], cv['knn5']['median_dist']]
    ax[1].bar(names, vals, color=[GR, RD, RD]); ax[1].axhline(pl['nearest_neighbour_dist_median'], color=OR, lw=2, ls='--')
    ax[1].text(2.4, pl['nearest_neighbour_dist_median'] * 1.6, 'typical gap to the\nnext orbit', ha='right', color=OR, fontsize=9)
    ax[1].set_yscale('log'); ax[1].set_ylabel('typical miss in the start (held-out orbits)')
    ax[1].set_title('B. Counting letters cannot place an orbit'); save(f, 'fig3-counting-letters-misses.png')

def fig4():
    Z = np.load(os.path.join(H, 'grid800.npz')); g = Z['g']; L = Z['L']; C = Z['C']
    o = [x for x in d if x['cls'] == 'I.A2'][0]; s = inverse(o['word']); code = {'a': 1, 'A': -1, 'b': 2, 'B': -2}
    tgt = np.array([code[c] for c in (s * 6)[:32]], np.int8)
    pref = np.cumprod(L == tgt[None, None, :], axis=2); ext = [0, g[-1], 0, g[-1]]
    f, ax = plt.subplots(1, 5, figsize=(21, 4.8))
    from matplotlib.colors import ListedColormap
    for i, k in enumerate((2, 4, 8, 16, 32)):
        m = np.where(C >= 0, pref[..., k - 1], np.nan)
        ax[i].imshow(m.T, origin='lower', extent=ext, cmap=ListedColormap(['#f2f2f2', BL]), vmin=0, vmax=1)
        ax[i].set_title(f'first {k} letters match\n{int(np.nansum(m)):,} of 502,041 starts', fontsize=11)
    for a in ax:
        a.plot(o['v1'], o['v2'], marker='o', mfc='none', mec='k', ms=12, mew=2); a.set_xlabel('v₁'); a.set_xticks([0, .4, .8]); a.set_yticks([0, .4, .8])
    ax[0].set_ylabel('v₂'); f.suptitle('Butterfly I (word abABABab, black ring): every start in the quarter disc, kept while its letters keep agreeing', y=1.03)
    save(f, 'fig4-letters-as-digits.png')

def fig5():
    cy = json.load(open(os.path.join(H, 'cylinders.json'))); cnt = np.array([r['counts'] for r in cy], float)
    k = np.arange(1, cnt.shape[1] + 1); med = np.median(cnt, 0); lo = np.percentile(cnt, 25, 0); hi = np.percentile(cnt, 75, 0)
    f, ax = plt.subplots(figsize=(7.5, 4.8))
    ax.fill_between(k, lo, hi, color=BL, alpha=.2); ax.semilogy(k, med, color=BL, lw=2.5, label='median over 695 words')
    fit = np.polyfit(k[7:], np.log(med[7:]), 1); ax.semilogy(k, np.exp(np.polyval(fit, k)), 'k--', lw=1, label=f'each letter keeps {np.exp(fit[0])*100:.0f}% of the starts')
    ax.set_xlabel('letters read'); ax.set_ylabel('starts (of 502,041) that still match'); ax.legend()
    ax.set_title('Each letter narrows the start, like a digit'); save(f, 'fig5-each-letter-narrows.png'); return float(np.exp(fit[0]))

def fig6(cls='I.A9'):
    import predict2 as Q
    from tbw import retdist, polish
    i = [j for j, o in enumerate(Q.d) if o['cls'] == cls][0]; o = Q.d[i]; w = Q.inverse(o['word'])
    f_ = (w * (Q.K0 // len(w) + 2))[:Q.K0]; tgt = np.array([Q.code[c] for c in f_], np.int8)
    pref = np.cumprod(Q.L0 == tgt[None, None, :], axis=2).sum(2); pref[Q.C0 < 0] = -1
    k = Q.K0
    while k > 4 and (pref >= k).sum() == 0: k -= 1
    ii, jj = np.nonzero(pref >= k); P1, P2 = Q.g[ii], Q.g[jj]
    Tst = float(np.array(Q.rfeat(o['word'])) @ Q.period_law_without(i))
    E = -2.5 + 3 * (P1**2 + P2**2); T = Tst / np.abs(E)**1.5
    dm, tm = retdist(P1, P2, T * .97, T * 1.03)
    cands, _, _ = Q.predict(w, Tst); best, c1, c2, tg = cands[0][:4]
    q, r = polish(c1, c2, tg)
    f, ax = plt.subplots(1, 3, figsize=(15, 4.8))
    a = ax[0]; a.scatter(v1, v2, s=3, color=GR, alpha=.4, label='all 695 catalog starts')
    a.scatter(P1, P2, s=4, color=BL, label=f'starts whose first {k} letters match the word')
    a.scatter([o['v1']], [o['v2']], s=120, facecolor='none', edgecolor='k', lw=2, label='catalog start (answer key, not used)')
    a.set_xlabel('v1'); a.set_ylabel('v2'); a.set_aspect('equal'); a.legend(fontsize=8, loc='lower left')
    a.set_title(f'1. The letters pick out thin strips\n({len(P1)} grid starts of 502,041)', fontsize=10)
    a = ax[1]; m = np.argsort(-dm); sc = a.scatter(P1[m], P2[m], c=np.log10(dm[m]), cmap='Blues_r', s=14)
    cb = plt.colorbar(sc, ax=a, shrink=.8); cb.set_label('log10 closest return to the start')
    a.scatter([c1], [c2], marker='x', s=90, color=BL, lw=2.5, label='best guess (closest return)')
    a.scatter([o['v1']], [o['v2']], s=120, facecolor='none', edgecolor='k', lw=2, label='catalog start')
    hw = 0.006; a.set_xlim(o['v1'] - hw, o['v1'] + hw); a.set_ylim(o['v2'] - hw, o['v2'] + hw); a.set_aspect('equal')
    a.legend(fontsize=8, loc='lower left'); a.ticklabel_format(useOffset=False); a.tick_params(labelsize=8)
    a.set_title(f'2. The run-law period ({Tst:.2f}) asks: which strip point\ncomes back to itself? Guess is {np.hypot(c1-o["v1"], c2-o["v2"]):.0e} off', fontsize=10)
    a = ax[2]; Xs = traj(q[0], q[1], q[2], 2500); rr = Xs[:, :6].reshape(-1, 3, 2)
    for kk, c in enumerate(['k', OR, BL]): a.plot(rr[:, kk, 0], rr[:, kk, 1], color=c, lw=.9)
    a.set_aspect('equal'); a.axis('off')
    a.set_title(f"3. After the polish: orbit {cls} ({len(w)} letters)\nlands on the catalog start, miss {np.hypot(q[0]-o['v1'], q[1]-o['v2']):.0e}", fontsize=10)
    save(f, 'fig6-word-to-start.png')

def tally():
    """One verdict per held-out orbit: hit (polish lands on the catalog start: residual < 1e-7, within 1e-6),
    twin (a periodic orbit with exactly the same word but another start), other (a periodic orbit with another word), none."""
    R = json.load(open(os.path.join(H, 'predict2_sample.json')))
    rp = os.path.join(H, 'retry.json'); Q = json.load(open(rp)) if os.path.exists(rp) else {}
    out = {}
    for c, r in R.items():
        fs = [f for rr in r['rotations'] for f in rr['finals']] + Q.get(c, {}).get('finals', [])
        if any(f['res'] < 1e-7 and f['dist'] < 1e-6 for f in fs): v = 'hit'
        elif any(f['res'] < 1e-7 and f.get('word_ok') for f in fs): v = 'twin'
        elif any(f['res'] < 1e-7 for f in fs): v = 'other'
        else: v = 'none'
        out[c] = dict(Lf=r['Lf'], v=v, retry=(c in Q and Q[c]['hit']))
    return out

def fig7():
    T = tally(); rows = sorted(T.items(), key=lambda x: x[1]['Lf'])
    col = dict(hit=BL, twin='#8e44ad', other=GR, none=RD)
    f, ax = plt.subplots(figsize=(12, 4.4))
    for i, (c, r) in enumerate(rows): ax.bar(i, r['Lf'], color=col[r['v']], hatch='//' if r['retry'] else None, edgecolor='white' if r['retry'] else None)
    ax.set_xticks(range(len(rows))); ax.set_xticklabels([c for c, _ in rows], rotation=90, fontsize=7)
    ax.set_ylabel('letters in the word')
    from matplotlib.patches import Patch
    n = {k: sum(r['v'] == k for _, r in rows) for k in col}
    ax.legend(handles=[Patch(color=BL, label=f"word alone → the catalog orbit ({n['hit']})"), Patch(color='#8e44ad', label=f"→ another orbit with the same word ({n['twin']})"),
                       Patch(color=GR, label=f"→ an orbit with a different word ({n['other']})"), Patch(color=RD, label=f"no orbit found ({n['none']})")], fontsize=9, loc='upper left')
    ax.set_title(f"From the word alone, {len(rows)} held-out orbits: {n['hit']} land on the catalog orbit (hatched = found on the second, finer pass)")
    save(f, 'fig7-word-alone-results.png')

def fig8():
    a = (0.3471168881, 0.5327249454, 6.3259139829); b = (0.3394074119878319, 0.5361849781525929, 6.290818765065431)
    f, ax = plt.subplots(1, 2, figsize=(11, 4.2))
    for A, p, ttl in [(ax[0], a, 'figure-eight (Moore 1993): all three share one curve'), (ax[1], b, 'same word aBAb, same period to 6 digits,\nbut each body has its own curve')]:
        Xs = traj(*p, 2000); r = Xs[:, :6].reshape(-1, 3, 2)
        for k, c in enumerate(['k', OR, BL]): A.plot(r[:, k, 0], r[:, k, 1], color=c, lw=1)
        A.set_aspect('equal'); A.axis('off'); A.set_title(ttl, fontsize=10)
    save(f, 'fig8-same-word-two-orbits.png')

if __name__ == '__main__':
    for nm in sys.argv[1:]: globals()[nm]()
