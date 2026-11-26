"""DRAFT. The orbits in hyperspherical coordinates (hyperradius R and the shape sphere), and whether period,
word length and stability come out simpler there. Uses Li & Liao's 695 orbits and the 66 polished orbits
(read only from ../three-body/polish/polished.json). Writes hyper.json."""
import json, os, sys, ctypes, numpy as np
from multiprocessing import Pool
H = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, H)
if not os.path.exists(os.path.join(H, "libhyp.so")):  # the compiled file is not kept in the repo; build it from hyp.c
    import subprocess
    subprocess.run(["gcc", "-O3", "-fopenmp", "-shared", "-fPIC", "-o", os.path.join(H, "libhyp.so"), os.path.join(H, "hyp.c"), "-lm"], check=True)
lib = ctypes.CDLL(os.path.join(H, 'libhyp.so')); D = ctypes.c_double
lib.hyper_stats.argtypes = [D, D, D, D, D, ctypes.POINTER(D)]; lib.hyper_stats.restype = ctypes.c_int
def stats(a):
    p1, p2, T = a; out = (D * 7)(); ns = lib.hyper_stats(p1, p2, T, 1e-12, 0.005, out); return list(out) + [ns]
def runs(w):
    c = [x.isupper() for x in w]; n = len(c); k = next((i for i in range(n) if c[i] != c[i - 1]), None)
    if k is None: return [n]
    c = c[k:] + c[:k]; out = []; cur = 1
    for i in range(1, n):
        if c[i] == c[i - 1]: cur += 1
        else: out.append(cur); cur = 1
    out.append(cur); return out
def rf(w): r = runs(w); return [sum(1 for x in r if x == k) for k in (1, 2, 3, 4, 5)]
def floquet(p1, p2, T, eps=1e-7):
    from tbw import ic, flow
    X0 = ic(p1, p2); M = np.zeros((12, 12))
    for j in range(12):
        e = np.zeros(12); e[j] = eps; M[:, j] = (flow(X0 + e, T) - flow(X0 - e, T)) / (2 * eps)
    return float(np.max(np.abs(np.linalg.eigvals(M))))
if __name__ == '__main__':
    d = json.load(open(os.path.join(H, 'data', 'li-liao-695.json')))
    pol = json.load(open(os.path.join(H, '..', 'three-body', 'polish', 'polished.json')))
    cat = [(o['v1'], o['v2'], o['T']) for o in d]; ours = [(float(o['p1']), float(o['p2']), float(o['T'])) for o in pol]
    with Pool(2) as P: sc = P.map(stats, cat); so = P.map(stats, ours)
    with Pool(2) as P: fl = P.starmap(floquet, ours)
    rec = lambda src, st, w: [dict(v1=a[0], v2=a[1], T=a[2], E=-2.5 + 3 * (a[0]**2 + a[1]**2), tau=s[0], arc=s[1], nRmax=int(s[2]), Rmin=s[3], Rmax=s[4], Imean=s[5], ok=s[7] > 0, word=ww) for a, s, ww in zip(src, st, w)]
    C = rec(cat, sc, [o['word'] for o in d]); O = rec(ours, so, [o['word'] for o in pol])
    for o, c in zip(O, pol): o['match'] = c.get('match')
    for o, f in zip(O, fl): o['maxlam'] = f
    json.dump(dict(catalog=C, ours=O), open(os.path.join(H, 'hyper.json'), 'w'))
    print('done', sum(c['ok'] for c in C), sum(o['ok'] for o in O))
