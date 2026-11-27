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
# 19-22 stability margin of zeta's slice (margin.py, margin_local.py)
rows = []
for Nn in range(3, 15):
    T0, T1 = 2 * pi * Nn**2 + 0.3, 2 * pi * (Nn + 1)**2 - 0.3
    tt = np.arange(T0, T1, 0.001); ns = np.arange(1, Nn + 1)
    ph = th(tt)[None, :] - tt[None, :] * np.log(ns)[:, None]
    f = (ns[:, None]**-0.5 * np.cos(ph)).sum(0)
    D = np.array([(-(ns**-0.5) * np.array([sp.multiplicity(p, int(k)) for k in ns]))[:, None] * np.sin(ph) for p in sp.primerange(2, Nn + 1)]).sum(1)
    sc = np.nonzero(np.sign(f[1:]) != np.sign(f[:-1]))[0]
    for a, b in zip(sc[:-1], sc[1:]):
        k = a + 1 + np.argmax(np.abs(f[a + 1:b + 1]))
        rows.append((Nn, tt[k], abs(f[k]), abs(f[k]) / np.abs(D[:, k]).max(), (tt[b] - tt[a]) * log(tt[k] / (2 * pi)) / (2 * pi), abs(f[k]) / np.linalg.norm(D[:, k])))
Rm = np.array(rows)
cg = np.corrcoef(np.log(Rm[:, 3]), np.log(Rm[:, 4]))[0, 1]
check('a gap\'s margin (smallest prime-angle turn that loses its two zeros, first order) follows the gap: correlation of logs > 0.75',
      cg > 0.75, f'{cg:.3f} over {len(Rm)} gaps')
o = np.argsort(Rm[:, 3])[:10]
tp = [abs(float(mp.siegelz(Rm[i, 1]))) for i in o]
rp = [abs(float(mp.siegelz(x))) for x in np.random.default_rng(5).choice(Rm[:, 1], 10, replace=False)]
check('at the 10 most fragile gaps true |Z| peaks are small (median < 0.5); at 10 random gaps they are not (median > 1): Lehmer-type near misses',
      np.median(tp) < 0.5 and np.median(rp) > 1, f'{np.median(tp):.2f} vs {np.median(rp):.2f}')
med = [np.median(Rm[Rm[:, 0] == k, 3]) for k in (3, 14)]
check('the typical gap margin is about one radian and drifts down with height (median at N = 3 above N = 14)', 0.5 < med[1] < med[0] < 2, f'{med[0]:.2f} -> {med[1]:.2f}')
Nn = 9; tt = np.arange(563.9, 564.8, 0.0005); ns = np.arange(1, 10)
e2 = np.array([sp.multiplicity(2, int(k)) for k in ns])
def g(dd): return (ns[:, None]**-0.5 * np.cos(th(tt)[None, :] - tt[None, :] * np.log(ns)[:, None] + dd * e2[:, None])).sum(0)
nz = lambda v: int(np.sum(np.sign(v[1:]) != np.sign(v[:-1])))
zt = Z[(Z > 563.9) & (Z < 564.8)]
check('at t = 564.3 the slice touches zero without crossing where zeta has two zeros; turning prime 2 by -0.0003 makes the pair',
      nz(g(0)) == 0 and nz(g(-3e-4)) == 2 and len(zt) == 2, f'true zeros {np.round(zt, 3).tolist()}')
# 23-29 functions that minimize the gaps (gaps/apart.py, gaps/closing.py, gaps/polys.py)
import json
Ap = json.load(open('gaps/apart.json'))['8']
Nn = 8; tt = np.arange(2 * pi * 64 + 0.3, 2 * pi * 81 - 0.3, 0.001); ns = np.arange(1, 9)
E8 = np.array([[sp.multiplicity(p, int(k)) for p in (2, 3, 5, 7)] for k in ns], float)
ph8 = th(tt)[None, :] - tt[None, :] * np.log(ns)[:, None]
def zs(dd):
    v = (ns[:, None]**-0.5 * np.cos(ph8 + (E8 @ dd)[:, None])).sum(0)
    i = np.nonzero(np.sign(v[1:]) != np.sign(v[:-1]))[0]; return tt[i]
def mingap(z): return float(np.min(np.diff(z) * np.log((z[1:] + z[:-1]) / 2 / (2 * pi)) / (2 * pi)))
z0 = zs(np.zeros(4)); zb = zs(np.array(Ap['best_d']))
check('pushing gaps apart (N = 8): the saved turn keeps all zeros on the line and widens the smallest gap more than 1.5 times',
      len(zb) == len(z0) and mingap(zb) > 1.5 * mingap(z0), f'{len(z0)} zeros, {mingap(z0):.3f} -> {mingap(zb):.3f} mean gaps, turn {np.linalg.norm(Ap["best_d"]):.2f} rad')
rng2 = np.random.default_rng(21)
kept = np.mean([len(zs(rng2.uniform(-pi, pi, 4))) >= len(z0) for _ in range(200)])
check('fewer than 5% of random angle settings keep every zero (N = 8): zeta sits in a small pocket', kept < 0.05, f'{100 * kept:.1f}%')
gain = Rm[:, 3] / Rm[:, 5]; Cm = np.median(Rm[:, 5] / Rm[:, 4]**2)
check('closing a gap: turning all primes at once (steepest direction) beats the best single prime by 1.1 to 1.4 times (median), and the turn is about 0.6 to 1.1 times the gap squared',
      1.1 < np.median(gain) < 1.4 and 0.6 < Cm < 1.1, f'gain {np.median(gain):.2f}, turn / gap^2 {Cm:.2f}')
# Lehmer's pair, from scratch: main sum (N = 33) plus the first Riemann-Siegel correction
tl = np.arange(7004.9, 7005.3, 0.0002); nl = np.arange(1, 34); psl = list(sp.primerange(2, 34))
El = np.array([[sp.multiplicity(p, int(k)) for p in psl] for k in nl], float)
phl = th(tl)[None, :] - tl[None, :] * np.log(nl)[:, None]
q = np.sqrt(tl / (2 * pi)) - 33; Rl = (-1)**32 * (tl / (2 * pi))**-0.25 * np.cos(2 * pi * (q * q - q - 1 / 16)) / np.cos(2 * pi * q)
def Zl(dd): return 2 * (nl[:, None]**-0.5 * np.cos(phl + (El @ dd)[:, None])).sum(0) + Rl
f0 = Zl(np.zeros(len(psl))); k = np.argmin(np.abs(tl - 7005.08))
Sl = -2 * (El.T @ (nl[:, None]**-0.5 * np.sin(phl))); dirn = -np.sign(f0[k]) * Sl[:, k] / np.linalg.norm(Sl[:, k])
zt = [float(x) for x in tl[np.nonzero(np.sign(f0[1:]) != np.sign(f0[:-1]))[0]]]
check("Lehmer's pair: zeta (main sum + first correction) has zeros at 7005.06 and 7005.10; a turn of 0.0010 rad along the steepest direction keeps them, 0.0013 loses them",
      nz(f0) == 2 and abs(zt[0] - 7005.062) < 0.002 and abs(zt[1] - 7005.100) < 0.002 and nz(Zl(0.0010 * dirn)) == 2 and nz(Zl(0.0013 * dirn)) == 0,
      f'zeros {np.round(zt, 3).tolist()}')
def BH(c, PM, twin=False):
    w2 = 2 if c % 2 == 0 else 0; v = (1 - w2 / 2) / (0.5**(2 if twin else 1))
    for p in sp.primerange(3, PM):
        def w(D):
            Dm = D % p; return 1 if Dm == 0 else (2 if pow(Dm, (p - 1) // 2, p) == 1 else 0)
        wf = w(1 - 4 * c) + (w(1 - 4 * (c + 2)) if twin else 0)
        v *= (1 - wf / p) / (1 - 1 / p)**(2 if twin else 1)
    return v
c41 = BH(41, 200000)
check("Euler's n^2 + n + 41: prime for n = 0..39, and its Hardy-Littlewood constant (C/2 in our units) is 3.3198 to 0.2% (OEIS A221712: 3.31977...)",
      all(sp.isprime(n * n + n + 41) for n in range(40)) and not sp.isprime(40 * 40 + 40 + 41) and abs(c41 / 2 / 3.31977 - 1) < 0.002, f'{c41 / 2:.4f}')
M2 = 100000
for cc, twin, tol in ((55661, False, 0.02), (46127, True, 0.06)):
    C_ = BH(cc, 200000, twin)
    fv = [n * n + n + cc for n in range(M2)]
    if not twin:
        cnt = sum(1 for v in fv if sp.isprime(v)); pred = C_ * sum(1 / log(v) for v in fv)
        check(f'n^2 + n + {cc}: prime values for n < 100,000 match the Bateman-Horn prediction within 2%, and its constant beats Euler\'s ({C_:.2f} vs {c41:.2f})',
              abs(cnt / pred - 1) < tol and C_ > c41, f'{cnt} vs {pred:.0f}')
    else:
        cnt = sum(1 for v in fv if sp.isprime(v) and sp.isprime(v + 2)); pred = C_ * sum(1 / log(v)**2 for v in fv)
        check(f'n^2 + n + {cc}: twin prime values (f, f+2) for n < 100,000 match the prediction within 6%; its twin constant {C_:.1f} is 15 times the integers\' 1.32',
              abs(cnt / pred - 1) < tol and C_ > 15 * 1.3203, f'{cnt} vs {pred:.0f}')
print(f'{sum(res)} of {len(res)} pass ({time.time() - t00:.0f} s)')
