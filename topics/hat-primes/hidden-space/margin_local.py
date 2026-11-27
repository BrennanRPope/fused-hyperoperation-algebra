"""Local (per gap) margins. Between each pair of neighbouring zeros of zeta's slice, Z_N has one peak.
Turning prime p's angle by d moves that peak by about d * dZ/dd_p; the pair of zeros is lost when the
peak reaches 0. So the gap's margin is about |peak| / max_p |dZ/dd_p| (first order), and the window's
margin is the smallest over its gaps. Also compares each slice peak with the true |Z| peak there
(mpmath siegelz), to see whether the slice's near-misses are Lehmer-type near-misses of zeta itself.
Writes margin_local.json, margin_local.npz."""
import json, numpy as np, sympy as sp, mpmath as mp
from math import pi, log
def theta(t): return t / 2 * np.log(t / (2 * pi)) - t / 2 - pi / 8 + 1 / (48 * t) + 7 / (5760 * t**3)
DT = 0.0005
rows = []
for N in range(3, 15):
    T0, T1 = 2 * pi * N**2 + 0.3, 2 * pi * (N + 1)**2 - 0.3
    t = np.arange(T0, T1, DT); th = theta(t); ns = np.arange(1, N + 1)
    ph = th[None, :] - t[None, :] * np.log(ns)[:, None]
    f = (ns[:, None]**-0.5 * np.cos(ph)).sum(0)
    ps = list(sp.primerange(2, N + 1))
    D = np.array([(-(ns**-0.5) * np.array([sp.multiplicity(p, int(n)) for n in ns]))[:, None].__mul__(np.sin(ph)).sum(0) for p in ps])
    sc = np.nonzero(np.sign(f[1:]) != np.sign(f[:-1]))[0]
    for a, b in zip(sc[:-1], sc[1:]):
        k = a + 1 + np.argmax(np.abs(f[a + 1:b + 1]))
        peak = abs(f[k]); sens = np.abs(D[:, k]).max()
        tz0, tz1 = t[a], t[b]
        gap = (tz1 - tz0) * log(t[k] / (2 * pi)) / (2 * pi)
        rows.append((N, t[k], peak, sens, peak / sens, gap))
R = np.array(rows)
# true Z at the slice's peaks, for the 40 smallest local margins and 40 random gaps
order = np.argsort(R[:, 4])
rng = np.random.default_rng(0)
pick = np.concatenate([order[:40], rng.choice(len(R), 40, replace=False)])
trueZ = np.array([abs(float(mp.siegelz(R[i, 1]))) for i in pick])
sliceZ = 2 * R[pick, 2]          # Z_N = 2 * sum
out = dict(gaps=len(R),
           corr_logmargin_loggap=float(np.corrcoef(np.log(R[:, 4]), np.log(R[:, 5]))[0, 1]),
           corr_logpeak_loggap=float(np.corrcoef(np.log(R[:, 2]), np.log(R[:, 5]))[0, 1]),
           smallest=[dict(N=int(R[i, 0]), t=float(R[i, 1]), slice_peak=float(2 * R[i, 2]), margin=float(R[i, 4]),
                          gap=float(R[i, 5]), true_peak=float(trueZ[j])) for j, i in enumerate(order[:12])],
           corr_slice_true_peak=float(np.corrcoef(np.log(sliceZ), np.log(trueZ))[0, 1]),
           median_true_peak_smallest40=float(np.median(trueZ[:40])), median_true_peak_random40=float(np.median(trueZ[40:])),
           median_margin_by_N={int(n): float(np.median(R[R[:, 0] == n, 4])) for n in range(3, 15)},
           frac_below_0p05_by_N={int(n): float(np.mean(R[R[:, 0] == n, 4] < 0.05)) for n in range(3, 15)})
json.dump(out, open('margin_local.json', 'w'), indent=1)
np.savez('margin_local.npz', R=R, pick=pick, trueZ=trueZ)
for k, v in out.items(): print(k, v)
