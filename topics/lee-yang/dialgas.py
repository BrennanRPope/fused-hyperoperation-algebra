"""Dial gas: the lattice gas whose two-site partition function is 1 + (x (+)_k y).

Z_G(z, k) = sum over occupied sets S of z^|S| * k^(edges inside S).
For one edge: 1 + x + y + k x y = 1 + wAdd(k, x, y).
counts(G)[m][e] = number of sets with m sites and e internal edges (exact integers).
"""
import numpy as np
from functools import lru_cache

def pc(x): return bin(x).count("1")

def cycle_counts(n):
    # transfer over a ring, exact bivariate polynomials as int arrays [m][e]
    import itertools
    C = np.zeros((n + 1, n + 1), dtype=object); C[:, :] = 0
    for s0 in (0, 1):
        v = {0: np.zeros((n + 1, n + 1), dtype=object), 1: np.zeros((n + 1, n + 1), dtype=object)}
        v[0][:, :] = 0; v[1][:, :] = 0
        v[s0][s0, 0] = 1
        for _ in range(n - 1):
            w = {0: v[0] * 0, 1: v[1] * 0}
            for s in (0, 1):
                for t in (0, 1):
                    sh = v[s]
                    if t == 0:
                        w[0] = w[0] + sh
                    else:
                        add = np.zeros_like(sh); add[:, :] = 0
                        de = s * t
                        add[1:, de:] = sh[:-1, :sh.shape[1] - de]
                        w[1] = w[1] + add
            v = w
        for s in (0, 1):
            de = s * s0
            add = np.zeros_like(v[s]); add[:, :] = 0
            add[:, de:] = v[s][:, :v[s].shape[1] - de]
            C = C + add
    return C.astype(np.int64) if n < 60 else C

def torus_counts(L, M=None):
    """L x M torus (L, M >= 3), sites N = L*M, edges 2N."""
    M = M or L
    S = 1 << L; N = L * M; E = 2 * N
    rot = lambda s: ((s << 1) | (s >> (L - 1))) & (S - 1)
    intra = np.array([pc(s & rot(s)) for s in range(S)])
    size = np.array([pc(s) for s in range(S)])
    band = np.array([[pc(s & t) for t in range(S)] for s in range(S)])  # between-column bonds
    total = np.zeros((N + 1, E + 1), dtype=np.int64)
    for s0 in range(S):
        v = np.zeros((S, N + 1, E + 1), dtype=np.int64)
        v[s0, size[s0], intra[s0]] = 1
        for _ in range(M - 1):
            w = np.zeros_like(v)
            for t in range(S):
                for b in range(L + 1):
                    mask = band[:, t] == b
                    if not mask.any(): continue
                    g = v[mask].sum(axis=0)
                    dm, de = size[t], intra[t] + b
                    if dm or de:
                        w[t, dm:, de:] += g[:N + 1 - dm, :E + 1 - de]
                    else:
                        w[t] += g
            v = w
        for s in range(S):
            de = band[s, s0]
            total[:, de:] += v[s][:, :E + 1 - de]
    return total

def brute_counts(nverts, edges):
    N = nverts
    C = np.zeros((N + 1, len(edges) + 1), dtype=np.int64)
    for S_ in range(1 << N):
        m = pc(S_); e = sum(1 for a, b in edges if (S_ >> a) & 1 and (S_ >> b) & 1)
        C[m, e] += 1
    return C

def torus_edges(L, M):
    idx = lambda i, j: (i % L) * M + (j % M)
    ed = []
    for i in range(L):
        for j in range(M):
            ed.append((idx(i, j), idx(i + 1, j))); ed.append((idx(i, j), idx(i, j + 1)))
    return ed

def z_coeffs(C, k):
    """coefficients a_m (m = 0..N) of Z as a polynomial in z, for numeric k (mpmath or float)."""
    import mpmath as mp
    N, E = C.shape[0] - 1, C.shape[1] - 1
    kp = [mp.mpf(1)] + [None] * E
    for e in range(1, E + 1): kp[e] = kp[e - 1] * k
    return [sum(int(C[m, e]) * kp[e] for e in range(E + 1) if C[m, e]) for m in range(N + 1)]

def z_zeros(C, k, dps=60):
    import mpmath as mp
    mp.mp.dps = dps
    a = z_coeffs(C, mp.mpf(k) if not isinstance(k, mp.mpc) else k)
    while a and a[-1] == 0: a.pop()
    if len(a) <= 1: return []
    for steps, extra in ((400, 400), (2000, 1500)):
        try:
            return mp.polyroots(a[::-1], maxsteps=steps, extraprec=extra)
        except Exception:
            pass
    import numpy as np  # repeated roots (k = 1): fall back to numpy
    return [mp.mpc(complex(r)) for r in np.roots([complex(x) for x in a[::-1]])]

def k_coeffs(C, z):
    import mpmath as mp
    N, E = C.shape[0] - 1, C.shape[1] - 1
    zp = [mp.mpf(1)] + [None] * N
    for m in range(1, N + 1): zp[m] = zp[m - 1] * z
    return [sum(int(C[m, e]) * zp[m] for m in range(N + 1) if C[m, e]) for e in range(E + 1)]

def evalZ(C, z, k):
    N, E = C.shape[0] - 1, C.shape[1] - 1
    return sum(int(C[m, e]) * z ** m * k ** e for m in range(N + 1) for e in range(E + 1) if C[m, e])
