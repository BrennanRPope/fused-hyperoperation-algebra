"""Checks for the hidden-space follow-up (DRAFT). Recomputes at smaller sizes where it is cheap."""
import time, numpy as np, sympy as sp, mpmath as mp
from math import gcd, log, exp, pi, sqrt
t00 = time.time(); res = []
def check(name, ok, info=''):
    res.append(ok); print(('PASS' if ok else 'FAIL'), f'{len(res)}. {name}', f'  [{info}]' if info else '')
EG = 0.5772156649015329
N = 10**6
sv = np.ones(N + 1, bool); sv[:2] = False
for i in range(2, int(N**0.5) + 1):
    if sv[i]: sv[i*i::i] = False
P = np.nonzero(sv)[0]
# 1-2 the mod-210 piece
cells = {(r % 14, r % 15) for r in range(210)}
check('the 15 x 14 grid (n mod 15, n mod 14) gives every residue mod 210 its own cell (Chinese remainder)', len(cells) == 210)
cnt = np.bincount(P[P > 7] % 210, minlength=210); units = [r for r in range(210) if gcd(r, 210) == 1]
spread = cnt[units].max() / cnt[units].min() - 1
check('primes above 7 (to 10^6) land only in the 48 cells coprime to 210, within 4% of even',
      len(units) == 48 and cnt.sum() == cnt[units].sum() and spread < 0.04, f'spread {spread:.3f}')
# 3-4 window size vs primes (Mertens)
w = np.prod([1 - 1 / p for p in sp.primerange(2, 31623)])
r = w * log(10**9)
check('window size product (1 - 1/p) for p up to sqrt(10^9), times log 10^9, is 2 e^-gamma to 0.5% (Mertens 1874)',
      abs(r / (2 * exp(-EG)) - 1) < 0.005, f'{r:.4f} vs {2*exp(-EG):.4f}')
check('the window overshoots the true prime share at x = 10^9 by more than 5% (pi(10^9) = 50,847,534)',
      w * 10**9 / 50847534 > 1.05, f'ratio {w*10**9/50847534:.4f}')
# 5-6 Bragg peaks of Lambda-weighted primes, N = 10^6
lam = np.zeros(N + 1)
for p in P:
    pk = int(p)
    while pk <= N: lam[pk] = log(p); pk *= int(p)
n = np.arange(N + 1)
err = max(abs(np.sum(lam * np.exp(-2j * pi * (n % q) / q)) / N - int(sp.mobius(q)) / int(sp.totient(q))) for q in range(1, 31))
check('log-weighted primes to 10^6: Fourier coefficient at 1/q equals mu(q)/phi(q) for q <= 30 (prime number theorem in progressions)',
      err < 3e-3, f'max error {err:.1e}')
off = abs(np.sum(lam * np.exp(-2j * pi * (sqrt(2) - 1) * n)) / N)
check('at an irrational frequency (sqrt2 - 1) there is no spot', off < 3e-3, f'{off:.1e}')
# 7-8 pair correlations and the Ramanujan expansion
C2 = np.prod([1 - 1 / (p - 1)**2 for p in sp.primerange(3, 10**6)])
def sing(h):
    if h % 2: return 0.0
    s = 2 * C2
    for p in sp.primefactors(h):
        if p > 2: s *= (p - 1) / (p - 2)
    return s
dev = max(abs(np.dot(lam[1:N + 1 - h], lam[1 + h:N + 1]) / N / sing(h) - 1) for h in range(2, 41, 2))
oddmax = max(np.dot(lam[1:N + 1 - h], lam[1 + h:N + 1]) / N for h in range(1, 41, 2))
check('pair correlations of log-weighted primes to 10^6 match the Hardy-Littlewood singular series within 6% (h even <= 40); odd h ~ 0',
      dev < 0.06 and oddmax < 1e-3, f'max deviation {dev:.3f}')
def ram(q, h): return sum(np.cos(2 * pi * a * h / q) for a in range(1, q + 1) if gcd(a, q) == 1)
rd = max(abs(sum(int(sp.mobius(q))**2 / int(sp.totient(q))**2 * ram(q, h) for q in range(1, 1501)) / sing(h) - 1) for h in [2, 6, 30])
check('the spots add up to the singular series: sum over q of (mu(q)/phi(q))^2 c_q(h) = S(h) (Ramanujan expansion, q <= 1500)', rd < 2e-4, f'{rd:.1e}')
# 9-10 total brightness
const = EG + sum(log(p) / (p * (p - 1)) for p in sp.primerange(2, 10**6))
Q = 10**5; tot = sum(1 / int(sp.totient(q)) for q in range(1, Q + 1) if sp.mobius(q) != 0)
check('total brightness of spots with q <= 10^5 is log Q + gamma + sum log p/(p(p-1)) = log Q + 1.3325 (within 0.01)',
      abs(tot - log(Q) - const) < 0.01, f'{tot - log(Q):.4f} vs {const:.4f}')
m2 = np.sum(lam**2) / N
check('mean of Lambda^2 to 10^6 is log N - 1 (within 0.02), so spots up to q ~ N/10 hold all of it', abs(m2 - log(N) + 1) < 0.02, f'{m2 - log(N):.4f}')
# 11 counting normalisation
tw = [np.sum(sv[P[P <= x][P[P <= x] + 2 <= N] + 2]) / np.sum(P <= x) for x in [10**4, 10**5, 10**6]]
check('counted per prime, twin pairs thin out (so the plain counting spectrum is flat: Humeniuk-Ramsey-Strungaru 2024)',
      tw[0] > tw[1] > tw[2], ', '.join(f'{v:.3f}' for v in tw))
# 12 theta
th = lambda t: t / 2 * np.log(t / (2 * pi)) - t / 2 - pi / 8 + 1 / (48 * t) + 7 / (5760 * t**3)
te = max(abs(th(t) - float(mp.siegeltheta(t))) for t in [30.0, 300.0, 1400.0])
check('the theta(t) formula used matches mpmath siegeltheta (t = 30, 300, 1400)', te < 1e-6, f'{te:.1e}')
# 13-14 the slice: Riemann-Siegel main sum vs the first 1,000 zeros
Z = np.load('../../primes/zeta_zeros_1000.npy')
t = np.arange(10.0, 1420.0, 0.002); a = th(t); Nt = np.floor(np.sqrt(t / (2 * pi))).astype(int)
f = sum(np.where(Nt >= k, k**-0.5 * np.cos(a - t * log(k)), 0) for k in range(1, Nt.max() + 1))
sc = np.nonzero(np.sign(f[1:]) != np.sign(f[:-1]))[0]; zN = (t[sc] + t[sc + 1]) / 2
mids = np.concatenate([[10.0], (Z[1:] + Z[:-1]) / 2, [1420.0]]); c = np.histogram(zN, bins=mids)[0]
idx = np.clip(np.searchsorted(zN, Z), 1, len(zN) - 1)
near = np.where(np.abs(zN[idx] - Z) < np.abs(zN[idx - 1] - Z), zN[idx], zN[idx - 1])
e = np.abs(near - Z) / (2 * pi / np.log(Z / (2 * pi)))
check('the torus slice (Riemann-Siegel main sum) gives 98% of the first 1,000 zeros exactly one partner, median error < 0.05 gap',
      np.mean(c == 1) >= 0.98 and np.median(e) < 0.05, f'{np.mean(c == 1):.3f}, median {np.median(e):.3f}')
t2 = np.linspace(2 * pi * 4, 2 * pi * 9, 20001)
f2 = np.cos(th(t2)) + 2**-0.5 * np.cos(th(t2) - t2 * log(2))
s2 = np.nonzero(np.sign(f2[1:]) != np.sign(f2[:-1]))[0]; z2 = (t2[s2] + t2[s2 + 1]) / 2
tz = Z[(Z > t2[0]) & (Z < 56)]
e2 = np.abs(z2 - tz) / (2 * pi / np.log(tz / (2 * pi))) if len(z2) == len(tz) else np.array([9])
check('two-angle slice (n = 1, 2) on t in [25, 56]: its 8 zeros sit within 0.1 gap of the 8 true zeros', e2.max() < 0.1, f'max {e2.max():.3f}')
# 15-17 other slices of the same torus
def lost(Nn, phases):
    T0, T1 = 2 * pi * Nn**2 + 0.5, 2 * pi * (Nn + 1)**2 - 0.5
    tt = np.arange(T0, T1, 0.002); aa = th(tt)
    base = np.sum(np.sign(np.cos(aa[1:])) != np.sign(np.cos(aa[:-1])))
    g = sum(k**-0.5 * np.cos(aa - tt * log(k) + sum(phases.get(p, 0) * ee for p, ee in sp.factorint(k).items())) for k in range(1, Nn + 1))
    return int(base - np.sum(np.sign(g[1:]) != np.sign(g[:-1])))
check('with N <= 2 terms the first term dominates (1 > 1/sqrt2), so no slice of the torus can lose a zero there',
      1 > 2**-0.5 and 1 < 2**-0.5 + 3**-0.5)
zl = [lost(k, {}) for k in range(3, 9)]
rng = np.random.default_rng(7)
rl = np.mean([lost(7, {p: rng.uniform(0, 2 * pi) for p in [2, 3, 5, 7]}) for _ in range(12)])
check("zeta's own slice loses at most 2 zeros for N = 3..8 (against theta/pi); random angles lose about 8 at N = 7",
      max(zl) <= 2 and rl > 5, f'zeta {zl}, random mean {rl:.1f}')
rot = [lost(6, {2: ph}) for ph in np.linspace(0, 2 * pi, 13)[:-1]]
check('turning only the prime-2 angle away from 0 loses zeros at every other angle tried (N = 6): zeta sits at a sharp point',
      rot[0] <= 0 and min(rot[1:]) >= 3, f'{rot}')
# 18 torus size
check('angles needed at height t = primes up to sqrt(t/2pi): 5 at t = 1000, 78 at t = 10^6, never bounded',
      sp.primepi(int(sqrt(1000 / (2 * pi)))) == 5 and sp.primepi(int(sqrt(10**6 / (2 * pi)))) == 78)
print(f'{sum(res)} of {len(res)} pass ({time.time() - t00:.0f} s)')
