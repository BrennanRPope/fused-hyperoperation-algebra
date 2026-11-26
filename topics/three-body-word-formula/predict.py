"""DRAFT. Word -> start, digit by digit. Uses no catalog starts: only the word.
1. The start must be a symmetric point of the word: s = reverse(s) with a<->b (every catalog word has exactly two).
2. On the 800x800 letter map, keep the starts whose first letters agree with s s s ...; group them into islands.
3. Zoom into each island (24x24 sub-grid), read more letters, keep the best-matching cells; repeat (beam of 6 islands).
4. Period guess from the time of letter |s| (period per letter ~2.43 at scale-free energy); polish with LM shooting.
Success = polished orbit lands on the catalog start (to 1e-6) with the same word."""
import json, os, sys, time, numpy as np
from scipy import ndimage
H = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, H)
from tbw import grid, raw, polish, creduce
Z = np.load(os.path.join(H, 'grid800.npz')); g = Z['g']; L0 = Z['L']; C0 = Z['C']; K0 = int(Z['K']); dv0 = g[1] - g[0]
code = {'a': 1, 'A': -1, 'b': 2, 'B': -2}
def inverse(w): return w[::-1].swapcase()
def sym(s): return s == s[::-1].translate(str.maketrans('abAB', 'baBA'))
def matchlen(L, tgt):
    eq = (L[..., :len(tgt)] == tgt)
    return np.cumprod(eq, axis=-1).sum(-1)
def islands(M, k, maxn=6):
    lab, n = ndimage.label(M >= k)
    out = []
    for sl, idx in zip(ndimage.find_objects(lab), range(1, n + 1)):
        out.append((int(M[sl][lab[sl] == idx].max()), sl))
    out.sort(key=lambda x: -x[0]); return out[:maxn]
def predict(s, beam=4, levels=12, sub=24, ntmax=480, verbose=False, hist=None):
    Lw = len(s); nt = min(max(2 * Lw + 2, 40), ntmax)
    full = (s * (ntmax // Lw + 2))[:ntmax]; TGT = np.array([code[c] for c in full], np.int8)
    M = np.where(C0 >= 0, matchlen(L0, TGT[:K0]), -1)
    best = M.max(); boxes = []
    for m, (si, sj) in islands(M, best - 1, beam):
        boxes.append((m, g[max(si.start - 1, 0)] - dv0 / 2, g[min(si.stop, len(g) - 1)] + dv0 / 2,
                      g[max(sj.start - 1, 0)] - dv0 / 2, g[min(sj.stop, len(g) - 1)] + dv0 / 2))
    for lev in range(levels):
        tgt = TGT[:nt]; new = []
        for (m, a0, a1, b0, b1) in boxes:
            x = np.linspace(a0, a1, sub); y = np.linspace(b0, b1, sub); X, Y = np.meshgrid(x, y, indexing='ij')
            E = np.minimum(-2.5 + 3 * (X**2 + Y**2), -0.05)
            Tmax = float(np.max(1.6 * 2.6 * nt / np.abs(E)**1.5)) + 5
            Lg, cg = grid(X, Y, Tmax, nt, tol=1e-11)
            Mg = matchlen(Lg.reshape(sub, sub, nt), tgt).astype(int)
            hx, hy = x[1] - x[0], y[1] - y[0]
            for mm, (si, sj) in islands(Mg, Mg.max(), 3):
                new.append((mm, x[si.start] - hx, x[si.stop - 1] + hx, y[sj.start] - hy, y[sj.stop - 1] + hy))
        new.sort(key=lambda b: (-b[0], b[2] - b[1])); boxes = new[:beam]
        if hist is not None: hist.append((nt, [tuple(map(float, b)) for b in boxes]))
        if verbose: print(lev, nt, [(b[0], float('%.3g' % (b[2] - b[1]))) for b in boxes])
        if not boxes: break
        if boxes[0][0] >= nt:
            if (boxes[0][2] - boxes[0][1]) < 1e-6: break
            if nt < ntmax and (boxes[0][2] - boxes[0][1]) > 0.1 * (a1 - a0 if False else 1e-30): nt = min(nt + Lw, ntmax)
    return boxes, nt
def finish(s, boxes, nt):
    """polish the top boxes; return list of (v1, v2, T, residual, word_ok)"""
    out = []
    for (m, a0, a1, b0, b1) in boxes[:3]:
        p1, p2 = (a0 + a1) / 2, (b0 + b1) / 2
        w, t, ok = raw(p1, p2, 600.0, maxev=len(s) + 1)
        if len(t) < len(s) + 1: continue
        Tg = 0.5 * (t[len(s) - 1] + t[len(s)])
        q, r = polish(p1, p2, Tg)
        w2, _, _ = raw(q[0], q[1], q[2] * 0.999999, maxev=4 * len(s))
        out.append(dict(match=int(m), guess=[p1, p2, Tg], v1=q[0], v2=q[1], T=q[2], res=r, word_ok=creduce(w2) == creduce(s) or w2 == s))
    return out
if __name__ == '__main__':
    d = json.load(open(os.path.join(H, 'data', 'li-liao-695.json')))
    which = sys.argv[1] if len(sys.argv) > 1 else 'sample'
    rng = np.random.default_rng(7)
    if which == 'sample':
        order = np.argsort([o['Lf'] for o in d]); pick = order[np.linspace(0, len(d) - 1, int(sys.argv[2]) if len(sys.argv) > 2 else 60).astype(int)]
    else:
        pick = [i for i, o in enumerate(d) if o['cls'] in which.split(',')]
    outp = os.path.join(H, f'predict_{which}.json'); results = json.load(open(outp)) if os.path.exists(outp) else {}
    for i in pick:
        o = d[i]
        if o['cls'] in results: continue
        t0 = time.time(); w = inverse(o['word'])
        rots = [w[r:] + w[:r] for r in range(len(w)) if sym(w[r:] + w[:r])]
        rec = dict(cls=o['cls'], Lf=o['Lf'], true=[o['v1'], o['v2'], o['T']], rotations=[])
        for s in rots:
            boxes, nt = predict(s)
            fin = finish(s, boxes, nt)
            hit = [f for f in fin if f['res'] < 1e-8 and np.hypot(f['v1'] - o['v1'], f['v2'] - o['v2']) < 1e-6]
            same_word = [f for f in fin if f['res'] < 1e-8 and f['word_ok']]
            rec['rotations'].append(dict(is_t0=(s == w), nt=nt, top_match=boxes[0][0] if boxes else 0, finals=fin, hit=bool(hit), same_word_orbit=bool(same_word)))
        rec['hit'] = any(r['hit'] for r in rec['rotations']); rec['same_word_orbit'] = any(r['same_word_orbit'] for r in rec['rotations']); rec['secs'] = time.time() - t0
        results[o['cls']] = rec
        print(o['cls'], o['Lf'], 'hit' if rec['hit'] else ('same-word orbit' if rec['same_word_orbit'] else 'MISS'), [(r['is_t0'], r['top_match'], r['nt'], [round(f['res'], 12) for f in r['finals']]) for r in rec['rotations']], round(rec['secs'], 1), flush=True)
        json.dump(results, open(outp, 'w'), default=float)
