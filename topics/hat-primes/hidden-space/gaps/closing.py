"""Approach 2: the turn of the prime angles that closes a gap fastest.
Near the peak P between two neighbouring zeros, turning the angles by a small vector d changes the slice by
d . S(t0), with S_p = dZ/dd_p = -2 sum_n n^-1/2 v_p(n) sin(theta - t log n). The fastest closing direction is
-sign(Z) S / |S|, so to first order the pair meets after a turn of
    m_best = P / |S|_2          (all primes at once, steepest direction)
against m_prime = P / max_p |S_p| for the best single prime (section 4). Part A tests this on every gap for
N = 3..14 and checks it against exact bisection along that direction. Part B uses the map t -> m_best as a
detector: on zeta near t = 7005 (main sum plus the first Riemann-Siegel correction) its smallest value
should sit at Lehmer's pair 7005.062 / 7005.100 (Lehmer 1956).
Writes closing.json, closing.npz."""
import json, numpy as np, sympy as sp
from math import pi, log
def theta(t): return t / 2 * np.log(t / (2 * pi)) - t / 2 - pi / 8 + 1 / (48 * t) + 7 / (5760 * t**3)
def setup(N, T0, T1, DT, rs=False):
    t = np.arange(T0, T1, DT); ns = np.arange(1, N + 1)
    ph = theta(t)[None, :] - t[None, :] * np.log(ns)[:, None]
    A = 2 * ns[:, None]**-0.5 * np.cos(ph); B = 2 * ns[:, None]**-0.5 * np.sin(ph)
    ps = list(sp.primerange(2, N + 1))
    E = np.array([[sp.multiplicity(p, int(n)) for p in ps] for n in ns], float)
    R = np.zeros_like(t)
    if rs:  # first Riemann-Siegel correction, held fixed under the turn
        q = np.sqrt(t / (2 * pi)) - N
        R = (-1)**(N - 1) * (t / (2 * pi))**-0.25 * np.cos(2 * pi * (q * q - q - 1 / 16)) / np.cos(2 * pi * q)
    return t, A, B, E, R, ps
def F(c, d):
    t, A, B, E, R, ps = c; ph = E @ d
    return np.cos(ph) @ A - np.sin(ph) @ B + R
def nz(f): return int(np.sum(np.sign(f[1:]) != np.sign(f[:-1])))
def gaps_info(c, DT):
    t, A, B, E, R, ps = c
    f = F(c, np.zeros(len(ps))); S = -(E.T @ B)
    sc = np.nonzero(np.sign(f[1:]) != np.sign(f[:-1]))[0]; rows = []
    for a, b in zip(sc[:-1], sc[1:]):
        k = a + 1 + np.argmax(np.abs(f[a + 1:b + 1])); P = abs(f[k]); Sk = S[:, k]
        g = t[b] - t[a]; s = g * log(t[k] / (2 * pi)) / (2 * pi)
        rows.append(dict(a=int(a), b=int(b), t0=float(t[k]), g=float(g), s=float(s), P=float(P),
                         m_best=float(P / np.linalg.norm(Sk)), m_prime=float(P / np.abs(Sk).max()),
                         dirn=(-np.sign(f[k]) * Sk / np.linalg.norm(Sk)).tolist(),
                         worst_prime=int(ps[np.argmax(np.abs(Sk))])))
    return f, rows
def exact_along(c, row, dirn):
    t, A, B, E, R, ps = c; a, b = row['a'], row['b']; w = b - a
    sl = slice(max(a - w, 0), min(b + w, len(t) - 1))
    cs = (t[sl], A[:, sl], B[:, sl], E, R[sl], ps)
    base = nz(F(cs, np.zeros(len(ps)))); lo, hi = 0.0, None
    for x in np.geomspace(1e-5, 2 * pi, 80):
        if nz(F(cs, x * np.array(dirn))) < base: hi = x; break
        lo = x
    if hi is None: return np.inf
    for _ in range(30):
        mid = (lo + hi) / 2
        if nz(F(cs, mid * np.array(dirn))) < base: hi = mid
        else: lo = mid
    return hi
# Part A
DT = 0.0005; allrows = []; caches = {}
for N in range(3, 15):
    c = setup(N, 2 * pi * N**2 + 0.3, 2 * pi * (N + 1)**2 - 0.3, DT); caches[N] = c
    _, rows = gaps_info(c, DT)
    for r in rows: r['N'] = N
    allrows += rows
mb = np.array([r['m_best'] for r in allrows]); mp_ = np.array([r['m_prime'] for r in allrows])
s = np.array([r['s'] for r in allrows])
rng = np.random.default_rng(11)
sel = np.concatenate([np.argsort(mb)[:30], rng.choice(len(allrows), 30, replace=False)])
ex = np.array([exact_along(caches[allrows[i]['N']], allrows[i], allrows[i]['dirn']) for i in sel])
ok = np.isfinite(ex)
fit = np.polyfit(np.log(s), np.log(mb), 1)
partA = dict(gaps=len(allrows), gain_median=float(np.median(mp_ / mb)), gain_q10=float(np.quantile(mp_ / mb, 0.1)),
             gain_q90=float(np.quantile(mp_ / mb, 0.9)),
             first_order_vs_exact_median=float(np.median(mb[sel][ok] / ex[ok])),
             first_order_vs_exact_q10=float(np.quantile(mb[sel][ok] / ex[ok], 0.1)),
             first_order_vs_exact_q90=float(np.quantile(mb[sel][ok] / ex[ok], 0.9)),
             C_median=float(np.median(mb / s**2)), beta_fit=float(fit[0]),
             smallest=[dict(N=allrows[i]['N'], t0=allrows[i]['t0'], s=allrows[i]['s'], m_best=allrows[i]['m_best'])
                       for i in np.argsort(mb)[:5]])
print('A', partA, flush=True)
# Part B: Lehmer's pair
DTB = 0.0002
c = setup(33, 6990.0, 7020.0, DTB, rs=True)
f, rowsB = gaps_info(c, DTB)
mbB = np.array([r['m_best'] for r in rowsB]); order = np.argsort(mbB)
top = rowsB[order[0]]
ex_top = exact_along(c, top, top['dirn'])
zt = c[0][np.nonzero(np.sign(f[1:]) != np.sign(f[:-1]))[0]]
partB = dict(window=[6990, 7020], N=33, gaps=len(rowsB),
             smallest=[dict(t0=rowsB[i]['t0'], g=rowsB[i]['g'], s=rowsB[i]['s'], P=rowsB[i]['P'], m_best=rowsB[i]['m_best'],
                            m_prime=rowsB[i]['m_prime'], worst_prime=rowsB[i]['worst_prime']) for i in order[:5]],
             top_exact=float(ex_top), median_m_best=float(np.median(mbB)),
             top_is_lehmer=bool(7005.0 < top['t0'] < 7005.15))
print('B', partB, flush=True)
json.dump(dict(A=partA, B=partB), open('closing.json', 'w'), indent=1)
np.savez('closing.npz', s=s, mb=mb, mp=mp_, sel=sel, ex=ex,
         tB=np.array([r['t0'] for r in rowsB]), mbB=mbB, sB=np.array([r['s'] for r in rowsB]),
         fB=f[::10], tfB=c[0][::10])
