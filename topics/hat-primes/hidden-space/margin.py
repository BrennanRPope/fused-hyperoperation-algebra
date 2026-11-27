"""Stability margin of zeta's slice of the prime torus.
Slice: Z_N(t) = sum_{n<=N} n^{-1/2} cos(theta(t) - t log n + sum_p v_p(n) d_p), N = floor(sqrt(t/2pi)).
Zeta's own slice has every turn d_p = 0. For each height window [2 pi N^2, 2 pi (N+1)^2] and for
blocks of 10 consecutive zeros inside it, find the smallest turn |d| of one prime's angle (either sign)
that changes the zero count, and the smallest turn along random directions of all angles at once.
Compare with how close the nearest pair of zeros is (Lehmer pairs). Writes margin.json, margin_blocks.npz."""
import json, numpy as np, sympy as sp
from math import pi, log, sqrt
def theta(t): return t / 2 * np.log(t / (2 * pi)) - t / 2 - pi / 8 + 1 / (48 * t) + 7 / (5760 * t**3)
DT = 0.0005
def setup(N, T0, T1):
    t = np.arange(T0, T1, DT); th = theta(t)
    ns = np.arange(1, N + 1); amp = ns**-0.5
    ph = th[None, :] - t[None, :] * np.log(ns)[:, None]
    C, S = amp[:, None] * np.cos(ph), amp[:, None] * np.sin(ph)
    ps = list(sp.primerange(2, N + 1))
    E = np.array([[sp.multiplicity(p, int(n)) for p in ps] for n in ns], float).reshape(N, len(ps))
    return t, C, S, E, ps
def zeros(t, C, S, E, d):
    phase = E @ d if E.shape[1] else np.zeros(C.shape[0])
    f = np.cos(phase) @ C - np.sin(phase) @ S
    sc = np.nonzero(np.sign(f[1:]) != np.sign(f[:-1]))[0]
    return (t[sc] + t[sc + 1]) / 2
def threshold(t, C, S, E, direction, base, dmax=np.pi):
    """smallest s > 0 where a pair of zeros is lost or made: |count - base| >= 2 (one zero sliding across
    the window edge changes the count by 1 and is not a loss)."""
    grid = np.concatenate([np.geomspace(1e-4, dmax, 60)])
    prev = 0.0
    for s in grid:
        if abs(len(zeros(t, C, S, E, s * direction)) - base) >= 2:
            lo, hi = prev, s
            for _ in range(30):
                mid = (lo + hi) / 2
                if abs(len(zeros(t, C, S, E, mid * direction)) - base) >= 2: hi = mid
                else: lo = mid
            return hi
        prev = s
    return np.inf
rng = np.random.default_rng(3)
win, blocks = [], []
for N in range(3, 15):
    T0, T1 = 2 * pi * N**2 + 0.3, 2 * pi * (N + 1)**2 - 0.3
    t, C, S, E, ps = setup(N, T0, T1)
    z0 = zeros(t, C, S, E, np.zeros(len(ps)))
    gaps = np.diff(z0) * np.log(z0[1:] / (2 * pi)) / (2 * pi)    # in mean gaps
    # whole-window margins, one prime at a time
    per_p = {}
    for i, p in enumerate(ps):
        e = np.zeros(len(ps)); e[i] = 1
        per_p[p] = min(threshold(t, C, S, E, e, len(z0)), threshold(t, C, S, E, -e, len(z0)))
    rnd = []
    for _ in range(8):
        v = rng.normal(size=len(ps)); v /= np.linalg.norm(v)
        rnd.append(threshold(t, C, S, E, v, len(z0)))
    win.append(dict(N=N, primes=ps, t=[T0, T1], zeros=len(z0), min_gap=float(gaps.min()),
                    margin_per_prime={str(k): float(v) for k, v in per_p.items()},
                    margin_one_prime=float(min(per_p.values())), margin_random_dir=float(np.min(rnd)),
                    margin_random_median=float(np.median(rnd))))
    print(win[-1], flush=True)
    # blocks of 10 zeros
    mids = np.concatenate([[T0], (z0[1:] + z0[:-1]) / 2, [T1]])
    for b in range(0, len(z0) - 9, 10):
        a0, a1 = mids[b], mids[b + 10]
        m = (t >= a0) & (t < a1)
        tb, Cb, Sb = t[m], C[:, m], S[:, m]
        base = len(zeros(tb, Cb, Sb, E, np.zeros(len(ps))))
        if base != 10: continue
        mm = np.inf
        for i, p in enumerate(ps):
            e = np.zeros(len(ps)); e[i] = 1
            mm = min(mm, threshold(tb, Cb, Sb, E, e, base), threshold(tb, Cb, Sb, E, -e, base))
        zb = z0[b:b + 10]
        g = np.diff(zb) * np.log(zb[1:] / (2 * pi)) / (2 * pi)
        # slope of Z at its zeros, normalised: small slope = near-miss
        blocks.append((N, float(zb[0]), float(mm), float(g.min())))
json.dump(dict(windows=win), open('margin.json', 'w'), indent=1)
np.savez('margin_blocks.npz', blocks=np.array(blocks))
B = np.array(blocks)
print('blocks', len(B), 'corr log margin vs log min gap', np.corrcoef(np.log(B[:, 2]), np.log(B[:, 3]))[0, 1])
