"""Hidden space 2: the infinite torus prod_p S^1 (Bohr). n^{-it} = prod_p (p^{-it})^{v_p(n)}, so a
Dirichlet sum is a function on the torus read along the Kronecker line phi_p = -t log p.
Measures how the zeta zeros sit on finite slices of it, via the Riemann-Siegel main sum
Z_N(t) = 2 sum_{n<=N} n^{-1/2} cos(theta(t) - t log n), N = floor(sqrt(t/2pi)).
Writes multiplicative.json and mult_curves.npz."""
import json, numpy as np, sympy as sp
from math import pi, log, sqrt
Z = np.load('../../primes/zeta_zeros_1000.npy')
def theta(t):
    return t / 2 * np.log(t / (2 * pi)) - t / 2 - pi / 8 + 1 / (48 * t) + 7 / (5760 * t**3)
out = {}
# (a) the torus dimension each height needs: primes up to sqrt(t/2pi)
Ts = [14, 50, 100, 300, 1000, 1419, 10**4, 10**6, 10**12, 10**22]
out['dimension'] = [dict(t=t, N=int(sqrt(t / (2 * pi))), primes=int(sp.primepi(int(sqrt(t / (2 * pi))))))
                    for t in Ts]
# (b) zeros of the main sum vs true zeros, t in [10, 1420]
t = np.arange(10.0, 1420.0, 0.002)
th = theta(t)
Nt = np.floor(np.sqrt(t / (2 * pi))).astype(int)
ZN = np.zeros_like(t)
for n in range(1, Nt.max() + 1):
    ZN += np.where(Nt >= n, 2 * n**-0.5 * np.cos(th - t * log(n)), 0)
sc = np.nonzero(np.sign(ZN[1:]) != np.sign(ZN[:-1]))[0]
zN = (t[sc] + t[sc + 1]) / 2
out['main_sum'] = dict(true_zeros=int(len(Z)), main_sum_zeros=int(len(zN)))
# nearest-match errors and gap-normalised errors
idx = np.clip(np.searchsorted(zN, Z), 1, len(zN) - 1)
near = np.where(np.abs(zN[idx] - Z) < np.abs(zN[idx - 1] - Z), zN[idx], zN[idx - 1])
gap = 2 * pi / np.log(Z / (2 * pi))
err = np.abs(near - Z) / gap
out['main_sum'].update(median_err_gaps=float(np.median(err)), max_err_gaps=float(err.max()),
                       within_quarter_gap=float(np.mean(err < 0.25)))
# count per window between true-zero midpoints: does each true zero get exactly one main-sum zero?
mids = np.concatenate([[10.0], (Z[1:] + Z[:-1]) / 2, [1420.0]])
c = np.histogram(zN, bins=mids)[0]
out['main_sum'].update(windows_exactly_one=float(np.mean(c == 1)), windows_zero=int(np.sum(c == 0)),
                       windows_two_plus=int(np.sum(c >= 2)))
# (c) the 2-torus slice, t in [2pi*4, 2pi*9): N = 2, Z_2 = 2[cos A + 2^{-1/2} cos(A - B)],
# A = theta(t), B = t log 2: the zeros are where the curve (A, B) crosses the zero set on the torus.
t2 = np.linspace(2 * pi * 4, 2 * pi * 9, 20001)
A, B = theta(t2) % (2 * pi), (t2 * log(2)) % (2 * pi)
f2 = np.cos(theta(t2)) + 2**-0.5 * np.cos(theta(t2) - t2 * log(2))
s2 = np.nonzero(np.sign(f2[1:]) != np.sign(f2[:-1]))[0]
z2 = (t2[s2] + t2[s2 + 1]) / 2
true2 = Z[(Z >= t2[0]) & (Z < t2[-1])]
out['torus2'] = dict(zeros_slice=[float(x) for x in z2], zeros_true=[float(x) for x in true2])
# (d) random shifts of the Kronecker curve: put independent phases on each prime and count zeros
# lost (deficit against the count of cos theta's zeros). Zeta is the shift with every phase 0.
rng = np.random.default_rng(1)
shift = []
for N in range(1, 9):
    T0 = 2 * pi * N**2 + 0.5; T1 = 2 * pi * (N + 1)**2 - 0.5
    tt = np.arange(T0, T1, 0.002); thv = theta(tt)
    base = int(np.sum(np.sign(np.cos(thv[1:])) != np.sign(np.cos(thv[:-1]))))
    ps = list(sp.primerange(2, N + 1)); vs = [sp.factorint(n) for n in range(1, N + 1)]
    def count(ph):
        f = np.zeros_like(tt)
        for n in range(1, N + 1):
            c = sum(ph[p] * e for p, e in vs[n - 1].items())
            f += n**-0.5 * np.cos(thv - tt * log(n) + c)
        return int(np.sum(np.sign(f[1:]) != np.sign(f[:-1])))
    z0 = count({p: 0.0 for p in ps})
    trials = [count({p: rng.uniform(0, 2 * pi) for p in ps}) for _ in range(40)] if ps else [z0]
    lost = [base - x for x in trials]
    tz = int(np.sum((Z >= T0) & (Z < T1))) if T1 < Z.max() else None
    shift.append(dict(N=N, primes=len(ps), t_range=[T0, T1], cos_theta_zeros=base, zeta_slice_zeros=z0,
                      true_zeros=tz, random_mean_lost=float(np.mean(lost)),
                      random_frac_none_lost=float(np.mean(np.array(lost) <= 0)),
                      dominant=bool(sum(n**-0.5 for n in range(2, N + 1)) < 1)))
out['shifts'] = shift
json.dump(out, open('multiplicative.json', 'w'), indent=1)
np.savez_compressed('mult_curves.npz', t2=t2, A=A, B=B, z2=z2, true2=true2,
                    t=t[::50], ZN=ZN[::50], zN=zN, err=err, Zt=Z)
print(out['dimension']); print(out['main_sum']); print(out['torus2'])
for s in shift: print(s)
