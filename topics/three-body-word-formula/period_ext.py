"""DRAFT. The run rule on orbits outside its fit: Li-Jing-Liao unequal masses (1,1,m3) and Li-Liao free fall.
For each mass group: (0) equal-mass prices used as they are, (1) the same prices times one scale for the group,
(2) prices refitted for the group, against (L) the straight letters law refitted for the group. All held out (leave one out)."""
import json, os, sys, numpy as np
H = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, H)
from hyper import runs
def rf(w, K=40): r = runs(w); return [sum(1 for x in r if x == k) for k in range(1, K + 1)]
base = np.array(json.load(open(os.path.join(H, 'period_runs.json')))['cost per run length']['coeffs'])
def loo(X, y):
    X = np.asarray(X, float); y = np.asarray(y, float); A = X.T @ X; b = X.T @ y; r = []
    for i in range(len(y)):
        Ai = A - np.outer(X[i], X[i]); bi = b - X[i] * y[i]
        try: r.append((X[i] @ np.linalg.lstsq(Ai, bi, rcond=None)[0] - y[i]) / y[i])
        except Exception: r.append(np.nan)
    return np.array(r)
def stats(r): r = r[np.isfinite(r)]; return dict(sd=float(np.sqrt(np.mean(r**2))), max=float(np.max(np.abs(r))), within1=float(np.mean(abs(r) < .01)))
def group(rows, key):
    out = {}
    for g in sorted(set(key(o) for o in rows), key=str):
        R = [o for o in rows if key(o) == g]; y = np.array([o['Tstar'] for o in R]); F = np.array([rf(o['word']) for o in R], float); nr = np.array([len(runs(o['word'])) for o in R], float)
        n = np.array([[len(o['word'])] for o in R], float)
        used = [k for k in range(F.shape[1]) if F[:, k].any()]
        res = dict(n=len(R), run_lengths={str(k + 1): int(F[:, k].sum()) for k in range(F.shape[1]) if F[:, k].any()},
                   direct=stats((F[:, :5] @ base - y) / y) if not F[:, 5:].any() else None, one_scale=stats(loo((F[:, :5] @ base)[:, None], y)) if not F[:, 5:].any() else None, letters_plus_reversals=stats(loo(np.c_[n[:, 0], nr], y)),
                   refit=stats(loo(F[:, used], y)) if len(R) > len(used) + 2 else None, letters=stats(loo(n, y)))
        if len(R) > len(used) + 2:
            c = np.linalg.lstsq(F[:, used], y, rcond=None)[0]; res['refit_prices'] = {str(k + 1): float(v) for k, v in zip(used, c)}
            res['scale'] = float(np.linalg.lstsq((F[:, :5] @ base)[:, None], y, rcond=None)[0][0])
        big = np.argsort(-abs(loo(F[:, used], y)))[:3] if len(R) > len(used) + 2 else []
        res['worst'] = [(R[i]['cls'], R[i]['word'][:40]) for i in big]
        out[str(g)] = res
    return out
if __name__ == '__main__':
    U = json.load(open(os.path.join(H, 'data', 'unequal-mass.json'))); Fo = json.load(open(os.path.join(H, 'data', 'free-fall.json')))
    out = dict(unequal=group(U, lambda o: o['m3']), free_fall=group(Fo, lambda o: tuple(o['masses'])))
    json.dump(out, open(os.path.join(H, 'period_ext.json'), 'w'), indent=1)
    for k, G in out.items():
        print(k)
        for g, r in G.items():
            f = lambda s: '%.2f%%' % (100 * s['sd']) if s else '-'
            print(' ', g, 'n', r['n'], 'runs', r['run_lengths'], 'direct', f(r['direct']), 'one scale', f(r['one_scale']), 'refit', f(r['refit']), 'letters', f(r['letters']), 'let+rev', f(r['letters_plus_reversals']), 'prices', {k: round(v, 2) for k, v in (r.get('refit_prices') or {}).items()})
