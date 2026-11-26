"""Sliding-route check for the open band of the Lee-Yang rule.
Zeros of G_a(s) = int_0^1 (1 - x^p)^a cos(s x) dx on 0 < s <= S, tracked as a slides up.
(1) Interlacing: do the zeros of G_a and G_(a+1) alternate (as for Bessel J_nu, J_(nu+1))?
(2) Collisions: does the gap between neighbouring zeros stay open while a <= p?
Polish: every bracketed zero refined with mpmath findroot, then rechecked.
"""
import json, time, sys
import numpy as np
import mpmath as mp
S = 30.0; DS = 0.05
def coeffs(p, a, K=130, dps=45):
    mp.mp.dps = dps; p = mp.mpf(p); a = mp.mpf(a)
    return [(-1)**k/mp.factorial(2*k)*mp.beta((2*k + 1)/p, a + 1)/p for k in range(K)]
def ev(c, s):
    s2 = mp.mpf(s)**2; acc = mp.mpf(0)
    for ck in reversed(c): acc = acc*s2 + ck
    return acc
def zeros(p, a, polish=False):
    c = coeffs(p, a); grid = np.arange(DS, S + 1e-9, DS)
    v = [ev(c, s) for s in grid]
    br = [(grid[i], grid[i + 1]) for i in range(len(v) - 1) if v[i]*v[i + 1] < 0]
    if not polish: return [float((l + r)/2) for l, r in br], c
    return [float(mp.findroot(lambda s: ev(c, s), (mp.mpf(l), mp.mpf(r)), solver="anderson")) for l, r in br], c
def interlace(z1, z2):
    """z1 zeros of G_a, z2 zeros of G_(a+1): require z1[0] < z2[0] < z1[1] < z2[1] < ... over the common range."""
    m = min(len(z1), len(z2)); seq = []
    if m == 0: return False
    for k in range(m):
        seq += [(z1[k], 0), (z2[k], 1)]
    seq.sort()
    return all(seq[i][1] != seq[i + 1][1] for i in range(len(seq) - 1)) and seq[0][1] == 0
if __name__ == "__main__":
    t0 = time.perf_counter(); R = {"S": S}
    # (1) interlacing in the band and around it
    il = {}
    for p in [1.5, 2.0, 3.0, 4.0, 5.0, 6.0]:
        for a in sorted(set([1.0, 1.25, 1.5] + list(np.arange(1.0, p + 1e-9, 0.5)))):
            za, _ = zeros(p, a, polish=True); zb, _ = zeros(p, a + 1, polish=True)
            il[f"p={p:g}, a={a:g}"] = dict(ok=interlace(za, zb), n_a=len(za), n_b=len(zb),
                                           in_band=bool(1 < a <= p), next_in_band=bool(a + 1 <= p))
    R["interlace"] = il
    # (2) zero tracks as a slides from 0 to p + 1.5 (for the picture and the collision test)
    tracks = {}
    for p in [1.5, 3.0]:
        A = np.round(np.arange(0.05, p + 1.5 + 1e-9, 0.05), 3); T = []
        for a in A:
            z, _ = zeros(p, a); T.append(z)
        tracks[str(p)] = dict(a=A.tolist(), zeros=T)
    R["tracks"] = tracks
    R["seconds"] = round(time.perf_counter() - t0, 1)
    json.dump(R, open("interlace_check.json", "w"), indent=1)
    print("seconds", R["seconds"])
