"""Pass 1, items 1-2: the prime cube.
Corners of the Boolean lattice over the primes are the squarefree numbers n; the corner S has weight prod_{p in S} p^-s.
Folding ∘ₖ over the primes x_p = p^-s gives 1 + k*(fold) = prod_p (1 + k p^-s) = sum over squarefree n of k^omega(n) n^-s =: Z_k(s).
Writes cube.json."""
import json, numpy as np, mpmath as mp
from sympy import primerange
mp.mp.dps = 30
out = {}

def fold(xs, k):
    acc = mp.mpf(0)
    for x in xs: acc = acc + x + k*acc*x
    return acc

# 1a. fold identity and closed forms at s = 2, 3 (P = primes below 2e5)
P = list(primerange(2, 200000))
rows = []
for s in (2, 3):
    for k in (-1, -0.5, 0.5, 1, 2):
        k = mp.mpf(k)
        xs = [mp.mpf(p)**(-s) for p in P]
        f = fold(xs, k)
        prod = mp.mpf(1)
        for x in xs: prod *= (1 + k*x)
        # G_k = prod (1+k p^-s)(1-p^-s)^k, Z_k = zeta^k G_k
        G = mp.mpf(1)
        for x in xs: G *= (1 + k*x)*(1 - x)**k
        rows.append(dict(s=s, k=float(k), fold_lift=float(1 + k*f), product=float(prod),
                         zeta_k_G=float(mp.zeta(s)**k*G)))
out['fold_rows'] = rows
out['closed_k1_s2'] = [float(mp.zeta(2)/mp.zeta(4))]
out['closed_km1_s2'] = [float(1/mp.zeta(2))]

# 1b. local zeros of the cube's factors: 1 + k p^-s = 0  <=>  p^s = -k  <=>  Re s = log|k|/log p
loc = {}
for k in (0.5, 1, 2, 3):
    loc[str(k)] = [(int(p), float(np.log(abs(k))/np.log(p))) for p in primerange(2, 60)]
out['local_zero_lines'] = loc

# 1c. order dial: near a zeta zero rho, Z_k ~ (s-rho)^k.  log|Z_k(rho+h)| / log|h| -> k  (G_k analytic there; check at s = rho + h with h -> 0
# using the factorisation zeta^k G_k, G_k evaluated at Re s = 1/2 + 0.3 where it converges, scaled along a ray from rho)
rho = mp.zetazero(1)
ords = []
for k in (-1, -0.5, 0.5, 1, 2):
    hs = [mp.mpf('1e-3'), mp.mpf('1e-5')]
    vals = [k*mp.log(abs(mp.zeta(rho + h))) for h in hs]
    slope = (vals[1]-vals[0])/(mp.log(hs[1])-mp.log(hs[0]))
    ords.append(dict(k=k, order=float(slope)))
out['order_at_first_zero'] = ords

# 2. Mertens walk: corner signs at k = -1 (Mobius) summed up to x
N = 10**7
mu = np.ones(N+1, dtype=np.int8); mu[0] = 0
isprime = np.ones(N+1, dtype=bool); isprime[:2] = False
for p in range(2, int(N**0.5)+1):
    if isprime[p]: isprime[p*p::p] = False
primes = np.nonzero(isprime)[0]
for p in primes:
    mu[p::p] *= -1
    pp = p*p
    if pp <= N: mu[pp::pp] = 0
M = np.cumsum(mu.astype(np.int64))
x = np.arange(N+1)
ratio = np.abs(M[1:]) / np.sqrt(x[1:])
out['mertens'] = dict(N=N, max_abs_M_over_sqrt_x=float(ratio.max()), argmax=int(ratio.argmax()+1),
                      M_at_N=int(M[N]), checkpoints={str(10**j): int(M[10**j]) for j in range(1, 8)})
# subsample for the figure
idx = np.unique(np.logspace(0, 7, 4000).astype(int))
out['mertens_curve'] = [[int(i), int(M[i])] for i in idx]
# a fair-coin walk on the squarefree corners for comparison (same support, random signs)
rng = np.random.default_rng(1)
sq = (mu != 0)
coin = np.where(sq, rng.choice([-1, 1], size=N+1), 0); coin[0] = 0
C = np.cumsum(coin)
out['coin_curve'] = [[int(i), int(C[i])] for i in idx]
json.dump(out, open('cube.json', 'w'), indent=1)
print(json.dumps({k: out[k] for k in ('fold_rows','order_at_first_zero','mertens')}, indent=1)[:3000])
