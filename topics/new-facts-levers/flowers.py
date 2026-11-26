"""Holographic candidates for n named regions: flower graphs (Czech-Shuai extreme rays of the cone of average
entropies): centre node joined to every leg; one leg (the stem) has weight w, the others weight 1.
Legs = the n named regions plus O, with possibly some legs trivial (absent). Entropy of a leg set I = min cut:
min(sum of weights in I, sum of weights not in I) for a single centre node."""
import numpy as np, itertools, json
def corner_poly(n, S):
    # S(maskset) over named subsets
    co = np.zeros(n)
    for m in range(1, 1 << n):
        j = bin(m).count('1'); co[j - 1] += S(m)
    return co
def root(co):
    r = np.roots(co[::-1]); r = [x.real for x in r if abs(x.imag) < 1e-9 and -1 < x.real < 0]
    return max(r) if r else None
def best_flower(n):
    best = (None, None)
    for legs in range(2, n + 2):                 # number of non-trivial legs among the n+1 parties
        for useO in (False, True):
            if useO and legs > n + 1: continue
            if not useO and legs > n: continue
            for w in range(1, 2 * legs + 1):
                for stem_is_O in ((True, False) if useO else (False,)):
                    # legs: named regions 0..L-1 nontrivial (L = legs - useO), O nontrivial if useO
                    L = legs - (1 if useO else 0)
                    wts = [1.0] * L
                    wO = 0.0
                    if useO: wO = w if stem_is_O else 1.0
                    if not stem_is_O and L > 0: wts[0] = float(w)
                    tot = sum(wts) + wO
                    def S(m):
                        a = sum(wts[i] for i in range(L) if m >> i & 1)
                        return min(a, tot - a)
                    r = root(corner_poly(n, S))
                    if r is not None and (best[0] is None or r > best[0]):
                        best = (r, dict(L=L, useO=useO, w=w, stem_is_O=stem_is_O))
    return best
out = {}
for n in range(3, 13):
    out[n] = best_flower(n); print(n, out[n])
json.dump({str(k): v for k, v in out.items()}, open('flowers.json', 'w'), indent=1)
