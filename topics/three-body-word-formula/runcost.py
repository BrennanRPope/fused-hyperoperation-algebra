"""DRAFT. Why does each run cost a fixed time? Cut every catalog orbit into its runs (boundaries halfway between the
last letter of one run and the first letter of the next), and measure each run's duration (at the common energy),
the line-ups with body 3 in the middle (not written in the word), and the close passes (a pair distance at a local
minimum) inside it. Writes runcost.json."""
import json, os, sys, ctypes, numpy as np
from multiprocessing import Pool
H = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, H)
from tbw import raw
if not os.path.exists(os.path.join(H, "libhyp.so")):  # the compiled file is not kept in the repo; build it from hyp.c
    import subprocess
    subprocess.run(["gcc", "-O3", "-fopenmp", "-shared", "-fPIC", "-o", os.path.join(H, "libhyp.so"), os.path.join(H, "hyp.c"), "-lm"], check=True)
lib = ctypes.CDLL(os.path.join(H, 'libhyp.so')); Dd = ctypes.c_double
lib.sample.argtypes = [Dd, Dd, Dd, ctypes.c_int, Dd, ctypes.POINTER(Dd)]
def one(o):
    p1, p2, T = o['v1'], o['v2'], o['T']; E = abs(-2.5 + 3 * (p1**2 + p2**2)); s = E**1.5
    w, t, ok = raw(p1, p2, T * 0.999999, maxev=4000)
    if not ok or len(w) < 2: return None
    n = len(w); up = [c.isupper() for c in w]
    if all(u == up[0] for u in up): return None
    k = next(i for i in range(n) if up[i] != up[i - 1])   # first letter of a run
    order = [(k + j) % n for j in range(n)]
    tt = np.array([t[i] if i >= k else t[i] + T for i in order])  # times in run order, increasing
    runs = []; cur = [0]
    for j in range(1, n):
        if up[order[j]] == up[order[j - 1]]: cur.append(j)
        else: runs.append(cur); cur = [j]
    runs.append(cur)
    # boundaries halfway between runs (cyclic)
    b = []
    for i, r in enumerate(runs):
        prev_last = tt[runs[i - 1][-1]] - (T if i == 0 else 0); b.append(0.5 * (prev_last + tt[r[0]]))
    b.append(b[0] + T)
    N = int(max(4000, 60 * T)); out = (Dd * (6 * (N + 1)))(); lib.sample(p1, p2, 2 * T, 2 * N, 1e-11, out) if False else None
    out = (Dd * (6 * (2 * N + 1)))(); lib.sample(p1, p2, 2 * T, 2 * N, 1e-11, out)
    r = np.frombuffer(out, float).reshape(-1, 3, 2); ts = np.linspace(0, 2 * T, 2 * N + 1)
    d = np.stack([np.hypot(*(r[:, i] - r[:, j]).T) for i, j in ((0, 1), (1, 2), (0, 2))], 1)
    area = (r[:, 1, 0] - r[:, 0, 0]) * (r[:, 2, 1] - r[:, 0, 1]) - (r[:, 1, 1] - r[:, 0, 1]) * (r[:, 2, 0] - r[:, 0, 0])
    sc = np.nonzero(np.sign(area[1:]) != np.sign(area[:-1]))[0]
    mid3 = [ts[i] for i in sc if (r[i, 2, 0] - r[i, 0, 0]) * (r[i, 2, 0] - r[i, 1, 0]) + (r[i, 2, 1] - r[i, 0, 1]) * (r[i, 2, 1] - r[i, 1, 1]) < 0]
    dm = d.min(1); i = np.arange(1, len(dm) - 1); cp = i[(dm[i] < dm[i - 1]) & (dm[i] <= dm[i + 1])]
    R = np.sqrt((r**2).sum((1, 2)))
    res = []
    for j, rr in enumerate(runs):
        a, c = b[j], b[j + 1]
        if a < 0: a += T; c += T
        res.append(dict(k=len(rr), dur=(c - a) * s, hidden=int(sum(a <= x < c for x in mid3)),
                        close=[float(dm[q] * E) for q in cp if a <= ts[q] < c], letters_gap=[float((tt[rr[q + 1]] - tt[rr[q]]) * s) for q in range(len(rr) - 1)]))
    return dict(cls=o['cls'], runs=res, Tstar=T * s)
if __name__ == '__main__':
    d = json.load(open(os.path.join(H, 'data', 'li-liao-695.json')))
    with Pool(4) as P: R = [x for x in P.map(one, d) if x]
    json.dump(R, open(os.path.join(H, 'runcost.json'), 'w'))
    allr = [r for o in R for r in o['runs']]
    for k in range(1, 6):
        x = [r for r in allr if r['k'] == k]
        if not x: continue
        du = np.array([r['dur'] for r in x]); hd = np.array([r['hidden'] for r in x]); nc = np.array([len(r['close']) for r in x])
        print(k, len(x), 'duration %.3f sd %.3f (%.1f%%)' % (du.mean(), du.std(), 100 * du.std() / du.mean()), 'hidden line-ups', np.bincount(hd)[:6], 'close passes', np.bincount(nc)[:8])
