"""Exact final check: all reduced words with n items (n = 0..3) and gap sum <= HB[n]; no two different
words share a fingerprint, and only the word N has fingerprint 1.

Fingerprints are hashed: K mod p1, K mod p2 (p1, p2 primes < 2^31) packed in 64 bits, where
K = T / (2^e2 5^e5 19^g) and e2, e5, g are read off T mod 32, T mod 25, T mod 19^4.  This is exact
because e2 <= 4, e5 <= 1, g <= 3 for words with <= 4 letters (item4.py, v19dp.py).
Every repeated hash (and every hash equal to that of 1) is re-checked with exact big integers.
"""
import numpy as np, sys, time, itertools
from engine import Tval, strip, items_upto, gap

P1, P2 = 2147483647, 2147483629
MODS = [32, 25, 19 ** 4, P1, P2]

def item_arrays(Hmax):
    its = items_upto(Hmax)
    L = np.array([i[0] for i in its], dtype=np.int64)
    G = np.array([gap(i) for i in its], dtype=np.int64)
    rp = {m: np.array([((19 if y == 'N' else 10) * pow(200, l, m) * pow(20, mm, m)) % m for (l, mm, y) in its],
                      dtype=np.int64) for m in MODS}
    return its, L, G, rp

def hashes_from_T(Tm):
    """Tm: dict mod -> residue array. returns uint64 hash of K."""
    t32, t25, t19 = Tm[32], Tm[25], Tm[19 ** 4]
    e2 = np.zeros_like(t32)
    for k in range(1, 5):
        e2 += (t32 % (2 ** k) == 0)
    assert not np.any(t32 % 32 == 0), 'v2 >= 5 found'
    e5 = (t25 % 5 == 0).astype(np.int64)
    assert not np.any(t25 % 25 == 0), 'v5 >= 2 found'
    g = np.zeros_like(t19)
    for k in range(1, 4):
        g += (t19 % (19 ** k) == 0)
    assert not np.any(t19 % 19 ** 4 == 0), 'v19 >= 4 found'
    out = []
    for p in (P1, P2):
        inv2, inv5, inv19 = pow(2, -1, p), pow(5, -1, p), pow(19, -1, p)
        K = Tm[p].copy()
        for (e, inv) in ((e2, inv2), (e5, inv5), (g, inv19)):
            for k in range(1, int(e.max()) + 1 if e.size else 1):
                K = np.where(e >= k, K * inv % p, K)
        out.append(K)
    return (out[0].astype(np.uint64) << np.uint64(32)) | out[1].astype(np.uint64)

def enumerate_hashes(n, H, its, L, G, rp, record_targets=None):
    """yield (hash array, decoder) chunks for all words with n items, gap sum <= H."""
    pw = {m: [pow(19, a, m) for a in range(H + 1)] for m in MODS}
    nits = len(its)
    if n == 0:
        Tm = {m: np.array([1], dtype=np.int64) for m in MODS}
        yield hashes_from_T(Tm), lambda k: ()
        return
    first = [i for i in range(nits) if its[i] != (0, 0, 'N') and G[i] <= H]
    if n == 1:
        idx = np.array(first)
        Tm = {m: (np.array([pw[m][int(a)] for a in L[idx]]) + rp[m][idx] * 1) % m for m in MODS}
        yield hashes_from_T(Tm), lambda k, idx=idx: (its[idx[k]],)
        return
    for i1 in first:
        A1 = int(L[i1]); h1 = int(G[i1])
        T1 = {m: (pw[m][A1] + int(rp[m][i1])) % m for m in MODS}
        rem = H - h1
        sel2 = np.nonzero(G <= rem)[0]
        if n == 2:
            A2 = A1 + L[sel2]
            Tm = {m: (np.array(pw[m], dtype=np.int64)[A2] + rp[m][sel2] * T1[m]) % m for m in MODS}
            yield hashes_from_T(Tm), lambda k, i1=i1, sel2=sel2: (its[i1], its[sel2[k]])
            continue
        # n == 3: pairs (i2, i3) with G2 + G3 <= rem
        I2, I3 = [], []
        for i2 in sel2:
            s3 = np.nonzero(G <= rem - G[i2])[0]
            I2.append(np.full(len(s3), i2)); I3.append(s3)
        I2 = np.concatenate(I2); I3 = np.concatenate(I3)
        A2 = A1 + L[I2]; A3 = A2 + L[I3]
        Tm = {}
        for m in MODS:
            pwa = np.array(pw[m], dtype=np.int64)
            T2 = (pwa[A2] + rp[m][I2] * T1[m]) % m
            Tm[m] = (pwa[A3] + rp[m][I3] * T2) % m
        yield hashes_from_T(Tm), lambda k, i1=i1, I2=I2, I3=I3: (its[i1], its[I2[k]], its[I3[k]])

def run(HB):
    t0 = time.time()
    Hmax = max(HB.values())
    its, L, G, rp = item_arrays(Hmax)
    allh = []
    counts = {}
    for n, H in sorted(HB.items()):
        c = 0
        for h, dec in enumerate_hashes(n, H, its, L, G, rp):
            allh.append(h); c += len(h)
        counts[n] = c
        print('n=%d items (%d letters), gap sum <= %d: %d words  (%.0fs)' % (n, n + 1, H, c, time.time() - t0), flush=True)
    allh = np.concatenate(allh)
    s = np.sort(allh)
    dup = s[1:][s[1:] == s[:-1]]
    one = hashes_from_T({m: np.array([1], dtype=np.int64) for m in MODS})[0]
    targets = set(int(x) for x in np.unique(dup)) | {int(one)}
    print('total %d words; repeated hashes: %d; re-checking %d target hashes exactly' % (len(allh), len(np.unique(dup)), len(targets)), flush=True)
    del s, allh
    found = {}
    tarr = np.array(sorted(targets), dtype=np.uint64)
    for n, H in sorted(HB.items()):
        for h, dec in enumerate_hashes(n, H, its, L, G, rp):
            ks = np.nonzero(np.isin(h, tarr))[0]
            for k in ks:
                w = dec(int(k))
                found.setdefault(int(h[k]), []).append(w)
    bad = []
    for hv, ws in found.items():
        exact = {}
        for w in ws:
            K = strip(Tval(list(w)))
            exact.setdefault(K, []).append(w)
        for K, wl in exact.items():
            if K == 1 and wl != [()]:
                bad.append(('kappa 1', wl))
            elif K != 1 and len(wl) > 1:
                bad.append(('collision', K, wl))
    print('exact re-check: %d hash classes, problems: %s  (%.0fs)' % (len(found), bad if bad else 'none', time.time() - t0), flush=True)
    return counts, bad

if __name__ == '__main__':
    HB = {int(a): int(b) for a, b in (x.split(':') for x in sys.argv[1:])}
    run(HB)
