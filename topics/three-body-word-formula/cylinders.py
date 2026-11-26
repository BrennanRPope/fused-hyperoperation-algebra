"""DRAFT. How fast does each letter of an orbit's word narrow down where it starts?
For every catalog orbit: the forward letters from its start are the inverse of Li & Liao's word, repeated.
On the 800x800 letter map, count the starts whose first k letters agree, for k = 1..32."""
import json, os, numpy as np
H = os.path.dirname(os.path.abspath(__file__))
Z = np.load(os.path.join(H, 'grid800.npz')); g = Z['g']; L = Z['L']; C = Z['C']; K = int(Z['K'])
N = len(g); dv = g[1] - g[0]
d = json.load(open(os.path.join(H, 'data', 'li-liao-695.json')))
code = {'a': 1, 'A': -1, 'b': 2, 'B': -2}
def inverse(w): return w[::-1].swapcase()
def sym(s): return s == s[::-1].translate(str.maketrans('abAB', 'baBA'))
valid = C >= 0
res = []
for o in d:
    s = inverse(o['word']); f = (s * (K // len(s) + 2))[:K]
    tgt = np.array([code[c] for c in f], np.int8)
    i = min(N - 1, int(o['v1'] / dv)); j = min(N - 1, int(o['v2'] / dv))
    eq = (L == tgt[None, None, :]) & valid[..., None]
    pref = np.cumprod(eq, axis=2)              # starts agreeing on the first k letters
    counts = pref.sum((0, 1))
    # cells matching, and how far the farthest matching cell is from the true start
    spread = []
    for k in (4, 8, 12, 16, 20, 24, 32):
        m = pref[..., k - 1].astype(bool)
        if m.any():
            ii, jj = np.nonzero(m); dist = np.hypot(g[ii] - o['v1'], g[jj] - o['v2'])
            spread.append((k, int(m.sum()), float(dist.min()), float(np.median(dist))))
        else: spread.append((k, 0, None, None))
    own = int(pref[i, j].sum())  # how many letters the true start's own grid cell gets right
    nsym = sum(sym(s[r:] + s[:r]) for r in range(len(s)))
    res.append(dict(cls=o['cls'], Lf=o['Lf'], counts=counts.tolist(), spread=spread, own_cell_match=own, sym_t0=sym(s), n_sym_rotations=nsym))
json.dump(res, open(os.path.join(H, 'cylinders.json'), 'w'))
cnt = np.array([r['counts'] for r in res], float)
print('sym at t0:', np.mean([r['sym_t0'] for r in res]), ' sym rotations per word:', np.bincount([r['n_sym_rotations'] for r in res]))
med = np.median(cnt, 0); print('median matching cells by k:', [int(x) for x in med])
print('fraction of words with >=1 matching cell at k=8,16,24,32:', [(cnt[:, k-1] > 0).mean() for k in (8,16,24,32)])
print('own-cell match median/mean:', np.median([r['own_cell_match'] for r in res]), np.mean([r['own_cell_match'] for r in res]))
