"""DRAFT. Word -> start, using both things the word tells us. No catalog start is used, and the period law is refitted
with the target orbit left out.
 1. Letters as digits: on the 800x800 letter map keep the starts whose first 32 letters agree with s s s ...
    (s = the word read forward from its symmetric starting letter).
 2. Period from the runs: T* = sum over runs of the price of that run length; T = T*/|E|^1.5 at each start.
 3. Among the kept starts, the one that comes back closest to itself within 3% of that period.
 4. Zoom twice around the best few, then the light polishing step (Levenberg-Marquardt shooting).
Success = the polished orbit is the catalog orbit (start within 1e-6, residual < 1e-8)."""
import json, os, sys, time, numpy as np
H = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, H)
from tbw import retdist, polish, raw, creduce
Z = np.load(os.path.join(H, 'grid800.npz')); g = Z['g']; L0 = Z['L']; C0 = Z['C']; K0 = int(Z['K']); dv = g[1] - g[0]
d = json.load(open(os.path.join(H, 'data', 'li-liao-695.json')))
code = {'a': 1, 'A': -1, 'b': 2, 'B': -2}
def inverse(w): return w[::-1].swapcase()
def sym(s): return s == s[::-1].translate(str.maketrans('abAB', 'baBA'))
def runs(w):
    n = len(w); i0 = next((i for i in range(n) if w[i].isupper() != w[i - 1].isupper()), 0); u = w[i0:] + w[:i0]
    out, cur = [], 1
    for x, y in zip(u, u[1:]):
        if x.isupper() == y.isupper(): cur += 1
        else: out.append(cur); cur = 1
    out.append(cur); return out
def rfeat(w): r = runs(w); return [sum(1 for x in r if x == k) for k in (1, 2, 3, 4, 5)]
RF = np.array([rfeat(o['word']) for o in d], float); TS = np.array([o['Tstar'] for o in d])
def period_law_without(i):
    m = np.ones(len(d), bool); m[i] = False
    return np.linalg.lstsq(RF[m], TS[m], rcond=None)[0]
def predict(s, Tstar, nbest=5, win=0.03, levels=3, shrink=4):
    f = (s * (K0 // len(s) + 2))[:K0]; tgt = np.array([code[c] for c in f], np.int8)
    pref = np.cumprod(L0 == tgt[None, None, :], axis=2).sum(2); pref[C0 < 0] = -1
    k = K0
    while k > 4 and (pref >= k).sum() == 0: k -= 1
    ii, jj = np.nonzero(pref >= k); P1, P2 = g[ii], g[jj]
    if len(P1) > 20000: sel = np.random.default_rng(0).choice(len(P1), 20000, replace=False); P1, P2 = P1[sel], P2[sel]
    E = -2.5 + 3 * (P1**2 + P2**2); T = Tstar / np.abs(E)**1.5
    dm, tm = retdist(P1, P2, T * (1 - win), T * (1 + win))
    order = np.argsort(dm); picks = []
    for q in order:
        if all(np.hypot(P1[q] - P1[p], P2[q] - P2[p]) > 3 * dv for p in picks): picks.append(q)
        if len(picks) >= nbest: break
    cands = []
    for q in picks:
        c1, c2, tg, best, hw = P1[q], P2[q], tm[q], dm[q], 1.5 * dv
        for lev in range(levels):
            x = np.linspace(c1 - hw, c1 + hw, 11); y = np.linspace(c2 - hw, c2 + hw, 11); X, Y = np.meshgrid(x, y, indexing='ij')
            ok = X**2 + Y**2 < 5 / 6 - 1e-3; X, Y = X[ok], Y[ok]
            dd, tt = retdist(X, Y, np.full(X.size, tg * (1 - win / 2)), np.full(X.size, tg * (1 + win / 2)))
            b = np.argmin(dd); c1, c2, tg, best = X[b], Y[b], tt[b], dd[b]; hw = hw / shrink
        cands.append((best, c1, c2, tg, int(k), len(P1)))
    return cands, k, len(P1)
if __name__ == '__main__':
    n = int(sys.argv[1]) if len(sys.argv) > 1 else 40
    order = np.argsort([o['Lf'] for o in d]); pick = order[np.linspace(0, len(d) - 1, n).astype(int)]
    outp = os.path.join(H, 'predict2_sample.json'); R = json.load(open(outp)) if os.path.exists(outp) else {}
    for i in pick:
        o = d[i]
        if o['cls'] in R: continue
        t0 = time.time(); w = inverse(o['word']); coef = period_law_without(i); Tst = float(np.array(rfeat(o['word'])) @ coef)
        rots = [w[r:] + w[:r] for r in range(len(w)) if sym(w[r:] + w[:r])]
        rec = dict(cls=o['cls'], Lf=o['Lf'], true=[o['v1'], o['v2'], o['T']], Tstar_true=o['Tstar'], Tstar_pred=Tst, rotations=[])
        for s in rots:
            cands, k, nreg = predict(s, Tst)
            fins = []
            for (best, c1, c2, tg, _, _) in cands[:4]:
                q, r = polish(c1, c2, tg)
                dist = float(np.hypot(q[0] - o['v1'], q[1] - o['v2']))
                wq, _, okq = raw(q[0], q[1], q[2] * 0.999999, maxev=4 * len(s) + 8) if r < 1e-6 else ('', None, False)
                fins.append(dict(guess=[float(c1), float(c2), float(tg)], guess_return=float(best), guess_dist=float(np.hypot(c1 - o['v1'], c2 - o['v2'])),
                                 v1=float(q[0]), v2=float(q[1]), T=float(q[2]), res=r, dist=dist, word_ok=(wq == s)))
                if r < 1e-8 and dist < 1e-6: break
            rec['rotations'].append(dict(is_t0=(s == w), k=k, region_cells=nreg, finals=fins,
                                         hit=any(f['res'] < 1e-8 and f['dist'] < 1e-6 for f in fins),
                                         same_word=any(f['res'] < 1e-8 and f['word_ok'] and f['dist'] >= 1e-6 for f in fins),
                                         other_orbit=any(f['res'] < 1e-8 and not f['word_ok'] for f in fins)))
        rec['hit'] = any(r['hit'] for r in rec['rotations']); rec['same_word'] = any(r['same_word'] for r in rec['rotations'])
        rec['other_orbit'] = any(r['other_orbit'] for r in rec['rotations']); rec['secs'] = time.time() - t0
        R[o['cls']] = rec; json.dump(R, open(outp, 'w'))
        r0 = rec['rotations'][0]
        print(o['cls'], o['Lf'], 'HIT' if rec['hit'] else ('twin' if rec['same_word'] else ('other' if rec['other_orbit'] else 'miss')),
              'T*err %.4f' % ((Tst - o['Tstar']) / o['Tstar']), 'k', r0['k'], 'cells', r0['region_cells'],
              'best guess dist %.2e' % min(f['guess_dist'] for f in r0['finals']), 'secs %.0f' % rec['secs'], flush=True)
