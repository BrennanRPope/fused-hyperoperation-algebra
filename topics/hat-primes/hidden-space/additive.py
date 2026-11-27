"""Hidden space 1: the profinite integers Z^ = prod_p Z_p (n -> its residues mod every prime power).
Measures (a) where primes land in a finite piece of it, (b) the prime window's size (Mertens),
(c) the Bragg peaks of the log-weighted primes at a/q and the pair correlations (Hardy-Littlewood).
Writes additive.json."""
import json, numpy as np, sympy as sp
from math import gcd, log, exp
EG = 0.5772156649015329
out = {}
N = 10**8
sieve = np.ones(N + 1, bool); sieve[:2] = False
for i in range(2, int(N**0.5) + 1):
    if sieve[i]: sieve[i*i::i] = False
P = np.nonzero(sieve)[0]
# (a) primes in the mod-210 torus (CRT grid 15 x 14)
M = 210
cnt = np.bincount(P[P > 7] % M, minlength=M)
units = [r for r in range(M) if gcd(r, M) == 1]
out['mod210'] = dict(units=len(units), counts_units=[int(cnt[r]) for r in units],
                     nonunit_hits=int(cnt[[r for r in range(M) if gcd(r, M) > 1]].sum()),
                     spread=float(cnt[units].max() / cnt[units].min() - 1))
# (b) window size prod_{p<=sqrt x}(1-1/p) against the true prime density pi(x)/x
pis = {10**9: 50847534}
rows = []
for e in [3, 4, 5, 6, 7, 8, 9]:
    x = 10**e
    z = int(x**0.5)
    w = np.prod([1 - 1/p for p in sp.primerange(2, z + 1)])
    pix = pis.get(x, int(np.searchsorted(P, x, side='right')))
    rows.append(dict(x=x, window=float(w), density=pix / x, ratio=float(w * x / pix),
                     ratio_loglaw=float(w * log(x))))
out['mertens'] = dict(rows=rows, two_e_minus_gamma=2 * exp(-EG))
# (c) Bragg peaks of Lambda-weighted primes, N1 = 10^7
N1 = 10**7
lam = np.zeros(N1 + 1)
for p in P[P <= N1]:
    pk = int(p)
    while pk <= N1:
        lam[pk] = log(p); pk *= int(p)
n = np.arange(N1 + 1)
peaks = []
for q in range(1, 41):
    for a in range(1, q + 1):
        if gcd(a, q) != 1: continue
        S = np.sum(lam * np.exp(-2j * np.pi * a * (n % q) / q)) / N1
        pred = int(sp.mobius(q)) / int(sp.totient(q))
        peaks.append(dict(a=a, q=q, re=float(S.real), im=float(S.imag), pred=pred))
        break   # one a per q is enough: S(a/q) depends on q only (Ramanujan sum)
out['bragg'] = peaks
# off-rational check: a generic frequency gives no peak
for al in [np.sqrt(2) - 1, (np.sqrt(5) - 1) / 2]:
    S = np.sum(lam * np.exp(-2j * np.pi * al * n)) / N1
    out.setdefault('offpeak', []).append(dict(alpha=float(al), abs=float(abs(S))))
# pair correlations against the singular series S(h) (Hardy-Littlewood 1923, conjecture)
C2 = np.prod([1 - 1 / (p - 1)**2 for p in sp.primerange(3, 10**6)])
def sing(h):
    if h % 2: return 0.0
    s = 2 * C2
    for p in sp.primefactors(h):
        if p > 2: s *= (p - 1) / (p - 2)
    return s
pairs = []
for h in range(1, 61):
    r = float(np.dot(lam[1:N1 + 1 - h], lam[1 + h:N1 + 1]) / N1)
    pairs.append(dict(h=h, measured=r, hl=sing(h)))
out['pairs'] = pairs
out['C2'] = float(C2)
# Ramanujan expansion: S(h) = sum_q mu(q)^2/phi(q)^2 c_q(h) (Bragg peaks -> pair correlations)
def ram(q, h):
    return sum(np.cos(2 * np.pi * a * h / q) for a in range(1, q + 1) if gcd(a, q) == 1)
rexp = []
for h in [2, 4, 6, 30]:
    s = sum(int(sp.mobius(q))**2 / int(sp.totient(q))**2 * ram(q, h) for q in range(1, 3001))
    rexp.append(dict(h=h, ramanujan_sum_q3000=float(s), hl=sing(h)))
out['ramanujan'] = rexp
# total intensity: sum_{q<=Q} mu^2/phi (all peaks with denominator <= Q) vs mean of Lambda^2
tot = []
for Q in [10, 100, 1000, 10000, 100000]:
    s = sum(1 / int(sp.totient(q)) for q in range(1, Q + 1) if sp.mobius(q) != 0)
    tot.append(dict(Q=Q, bragg_total=float(s), logQ=log(Q)))
out['bragg_total'] = tot
out['lambda_sq_mean'] = dict(N=N1, mean=float(np.sum(lam**2) / N1), logN=log(N1))
# counting normalisation (Humeniuk-Ramsey-Strungaru 2024): pairs at h=2 per prime -> 0
cn = []
for x in [10**4, 10**5, 10**6, 10**7, 10**8]:
    Px = P[P <= x]; s = set(Px.tolist()) if x <= 10**6 else None
    tw = int(np.sum(sieve[Px[Px + 2 <= N] + 2]))
    cn.append(dict(x=x, twin_pairs=tw, per_prime=tw / len(Px)))
out['counting'] = cn
json.dump(out, open('additive.json', 'w'), indent=1)
print('mod210 units', len(units), 'spread', out['mod210']['spread'], 'nonunit', out['mod210']['nonunit_hits'])
for r in rows: print('mertens', r)
print('2e^-g', 2 * exp(-EG))
for p in peaks[:12]: print('bragg', p)
print('offpeak', out['offpeak'])
for p in pairs[:8]: print('pair', p)
print('ram', rexp); print('tot', tot, out['lambda_sq_mean']); print('count', cn)
