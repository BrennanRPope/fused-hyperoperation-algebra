"""Thresholds of explicit state families, n named regions, purifier O trivial or included.
PT_N: perfect tensor / star graph on N legs, S = min(j, N-j)  (holographic)
MDS(N, r): GHZ-type state from an [N, r] MDS code, S = min(j, r, N-j)  (quantum; holographic only when r >= N/2, = PT)
"""
import numpy as np, json, itertools
from math import comb
def root(coeffs):
    r = np.roots(coeffs[::-1]); r = [x.real for x in r if abs(x.imag) < 1e-9 and -1 < x.real < 0]
    return max(r) if r else None
def poly_named(n, g, m):
    """state on m named regions (others trivial) + O trivial: f = (1+k)^(n-m) * f_m ; root of f_m"""
    return [comb(m, j) * g(j) for j in range(1, m + 1)]
out = {'PT': {}, 'MDS': {}}
for N in range(4, 41, 2):
    out['PT'][N] = root(poly_named(N, lambda j: min(j, N - j), N))
for N in range(4, 41):
    best = None
    for r in range(1, N // 2 + 1):
        x = root(poly_named(N, lambda j: min(j, r, N - j), N))
        if x is not None and (best is None or x > best[0]): best = (x, r)
    out['MDS'][N] = best
for N in range(4, 21):
    print(N, 'PT', out['PT'].get(N), 'MDS best', out['MDS'][N])
json.dump(out, open('families.json', 'w'), indent=1)
