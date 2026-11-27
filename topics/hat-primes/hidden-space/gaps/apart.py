"""Approach 1: push zeta's zero gaps apart by turning the prime angles.
The slice is F(t; d) = sum_{n<=N} n^-1/2 cos(theta(t) - t log n + sum_p v_p(n) d_p), N = floor(sqrt(t/2pi)).
d = 0 is zeta's own main sum. For each block of fixed N we search the torus of angles d (one per prime <= N)
for the point whose smallest zero gap (in mean-gap units s = g L / 2pi) is as large as possible, while every
zero stays real (the zero count in the block must not drop). We compare with zeta's own point and with random points.
Writes apart.json, apart.npz."""
import json, numpy as np, sympy as sp
from math import pi
from scipy.optimize import minimize
def theta(t): return t / 2 * np.log(t / (2 * pi)) - t / 2 - pi / 8 + 1 / (48 * t) + 7 / (5760 * t**3)
DT = 0.001
rng = np.random.default_rng(5)
def block(N):
    T0, T1 = 2 * pi * N**2 + 0.3, 2 * pi * (N + 1)**2 - 0.3
    t = np.arange(T0, T1, DT); ns = np.arange(1, N + 1)
    ph = theta(t)[None, :] - t[None, :] * np.log(ns)[:, None]
    A = ns[:, None]**-0.5 * np.cos(ph); B = ns[:, None]**-0.5 * np.sin(ph)
    ps = list(sp.primerange(2, N + 1))
    E = np.array([[sp.multiplicity(p, int(n)) for p in ps] for n in ns], float)
    return t, A, B, E, ps
def zeros(t, A, B, E, d):
    ph = E @ d; F = np.cos(ph) @ A - np.sin(ph) @ B
    i = np.nonzero(np.sign(F[1:]) != np.sign(F[:-1]))[0]
    return t[i] - F[i] * DT / (F[i + 1] - F[i])
def gaps(z):
    g = np.diff(z); m = (z[1:] + z[:-1]) / 2
    return g * np.log(m / (2 * pi)) / (2 * pi)
out = {}; keep = {}
for N in (5, 6, 8, 10, 12):
    t, A, B, E, ps = block(N)
    z0 = zeros(t, A, B, E, np.zeros(len(ps))); base = len(z0); s0 = gaps(z0)
    def score(d, tau=None):
        z = zeros(t, A, B, E, d)
        if len(z) < base: return -1.0 - (base - len(z)), None
        s = gaps(z)
        if tau is None: return s.min(), s
        return -tau * np.log(np.sum(np.exp(-s / tau))), s     # soft minimum
    # random points: how many keep every zero, and their smallest gap
    R = rng.uniform(-pi, pi, (400, len(ps)))
    rs = np.array([score(d)[0] for d in R]); okr = rs > 0
    # search: best random starts, then Nelder-Mead on the soft minimum
    starts = list(R[np.argsort(-rs)[:6]]) + [np.zeros(len(ps))]
    best_d, best = None, -np.inf
    for d0 in starts:
        res = minimize(lambda d: -score(d, 0.03)[0], d0, method='Nelder-Mead',
                       options=dict(maxiter=1500, xatol=1e-4, fatol=1e-6))
        v = score(res.x)[0]
        if v > best: best, best_d = v, (res.x + pi) % (2 * pi) - pi
    zb = zeros(t, A, B, E, best_d); sb = gaps(zb)
    # how far zeta's own point is from a local best: climb from d = 0 only
    loc = minimize(lambda d: -score(d, 0.03)[0], np.zeros(len(ps)), method='Nelder-Mead',
                   options=dict(maxiter=1500, xatol=1e-4, fatol=1e-6))
    sl = score(loc.x)[1]
    out[N] = dict(primes=ps, zeros=base, zeta_min_gap=float(s0.min()), zeta_gap_var=float(s0.var()),
                  random_keep_all_zeros=float(okr.mean()),
                  random_min_gap_median_when_kept=float(np.median(rs[okr])) if okr.any() else None,
                  best_min_gap=float(best), best_gap_var=float(sb.var()), best_d=[float(x) for x in best_d],
                  best_distance_from_zeta=float(np.linalg.norm(best_d)),
                  local_from_zeta_min_gap=float(sl.min()), local_from_zeta_distance=float(np.linalg.norm(loc.x)),
                  best_zeros_kept=len(zb))
    keep[f'z0_{N}'] = z0; keep[f'zb_{N}'] = zb; keep[f'zl_{N}'] = zeros(t, A, B, E, loc.x)
    keep[f'rs_{N}'] = rs; keep[f'bd_{N}'] = best_d; keep[f'ld_{N}'] = loc.x
    print(N, {k: v for k, v in out[N].items() if k != 'best_d'}, flush=True)
json.dump(out, open('apart.json', 'w'), indent=1)
np.savez('apart.npz', **keep)
