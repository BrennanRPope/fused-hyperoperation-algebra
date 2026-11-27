"""Approach 3: functions whose values hit primes with the smallest gaps.
For f(n) = n^2 + n + c the Bateman-Horn (Hardy-Littlewood Conjecture F) prediction is
    #{n <= M : f(n) prime} ~ C(f) * sum_{n<=M} 1/log f(n),   C(f) = prod_p (1 - w(p)/p) / (1 - 1/p),
w(p) = number of roots of f mod p. So the average gap between n's with f(n) prime is about log f(n) / C(f):
a bigger C means smaller gaps. In the hidden space (the profinite integers) C(f) measures how much of the
prime window f's values avoid losing: f misses every residue class 0 mod p for which w(p) = 0.
For twin values f(n), f(n)+2 both prime the constant is prod_p (1 - w2(p)/p) / (1 - 1/p)^2 with w2 the roots of
f(f+2). We rank c <= 100000 by both constants, then count prime values for n < 10^6 with a polynomial sieve
and check the prediction. Writes polys.json, polys.npz."""
import json, numpy as np, sympy as sp
from math import log
CMAX = 100000; cs = np.arange(1, CMAX + 1)
def omega(p, D):     # roots of n^2+n+c mod odd p, from the discriminant D = 1 - 4c
    qr = np.zeros(p, bool); qr[(np.arange(p)**2) % p] = True
    Dm = D % p
    return np.where(Dm == 0, 1, np.where(qr[Dm], 2, 0))
def consts(cs, PMAX):
    logC = np.zeros(len(cs)); logT = np.zeros(len(cs))
    w2 = np.where(cs % 2 == 0, 2, 0)          # p = 2: f = n(n+1)+c has the parity of c
    logC += np.log(np.maximum(1 - w2 / 2, 1e-300)) - np.log(1 - 1 / 2)
    logT += np.log(np.maximum(1 - w2 / 2, 1e-300)) - 2 * np.log(1 - 1 / 2)
    for p in sp.primerange(3, PMAX):
        wf = omega(p, 1 - 4 * cs); wg = omega(p, 1 - 4 * (cs + 2))
        logC += np.log(1 - wf / p) - np.log(1 - 1 / p)
        logT += np.log(np.maximum(1 - (wf + wg) / p, 1e-300)) - 2 * np.log(1 - 1 / p)
    return np.exp(logC), np.exp(logT)
def consts_exact(cl, PMAX):                      # per-element Legendre symbols, for a short list of c
    out = []
    for c in cl:
        w2 = 2 if c % 2 == 0 else 0
        lc = log(max(1 - w2 / 2, 1e-300)) - log(0.5); lt = log(max(1 - w2 / 2, 1e-300)) - 2 * log(0.5)
        for p in PR:
            if p == 2 or p >= PMAX: continue
            def w(D):
                Dm = D % p
                return 1 if Dm == 0 else (2 if pow(Dm, (p - 1) // 2, p) == 1 else 0)
            wf, wg = w(1 - 4 * c), w(1 - 4 * (c + 2))
            lc += log(1 - wf / p) - log(1 - 1 / p)
            lt += log(max(1 - (wf + wg) / p, 1e-300)) - 2 * log(1 - 1 / p)
        out.append((np.exp(lc), np.exp(lt)))
    return np.array([o[0] for o in out]), np.array([o[1] for o in out])
def tonelli(n, p):
    if p % 4 == 3: return pow(n, (p + 1) // 4, p)
    q, s = p - 1, 0
    while q % 2 == 0: q //= 2; s += 1
    z = 2
    while pow(z, (p - 1) // 2, p) != p - 1: z += 1
    m, c, t, r = s, pow(z, q, p), pow(n, q, p), pow(n, (q + 1) // 2, p)
    while t != 1:
        i, t2 = 0, t
        while t2 != 1: t2 = t2 * t2 % p; i += 1
        b = pow(c, 1 << (m - i - 1), p); m, c, t, r = i, b * b % p, t * b * b % p, r * b % p
    return r
PR = list(sp.primerange(2, 1_000_100))
C1, T1 = consts(cs, 2000)                         # quick ranking
topC = cs[np.argsort(-C1)[:8]]; topT = cs[np.argsort(-T1)[:8]]
pick = sorted(set(topC.tolist()) | set(topT.tolist()) | {41, 1, 17, 11})
Cr, Tr = consts_exact(pick, 1_000_000)        # refined constants
# polynomial sieve, n < M
M = 1_000_000; PS = 1_000_100
plist = np.array(PR)
isp_small = np.zeros(PS, bool); isp_small[plist] = True
def prime_values(c):
    n = np.arange(M, dtype=np.int64); f = n * n + n + c
    ok = np.ones(M, bool); D = 1 - 4 * c
    for p in plist.tolist():
        if p == 2:
            if c % 2 == 0: ok[:] = False
            continue
        Dm = D % p
        if Dm == 0: roots = [((-1) * pow(2, -1, p)) % p]
        elif pow(Dm, (p - 1) // 2, p) != 1: continue
        else:
            r = tonelli(Dm, p); inv2 = pow(2, -1, p)
            roots = list({((-1 + r) * inv2) % p, ((-1 - r) * inv2) % p})
        for r0 in roots: ok[r0::p] = False
    small = f < PS                                   # values that are themselves small primes
    ok[small] = isp_small[f[small]]
    return ok, f
res = {}; keep = {}
for c in pick:
    okf, f = prime_values(c); okg, _ = prime_values(c + 2)
    n1 = int(okf.sum()); n2 = int((okf & okg).sum())
    i = pick.index(c)
    pred1 = Cr[i] * np.sum(1 / np.log(np.maximum(f[1:], 3)))
    pred2 = Tr[i] * np.sum(1 / np.log(np.maximum(f[1:], 3))**2)
    first = int(np.argmin(okf)) if not okf.all() else M      # run of primes from n = 0
    idx = np.nonzero(okf)[0]; meangap = float(np.mean(np.diff(idx))) if len(idx) > 1 else None
    res[c] = dict(C=float(Cr[i]), C_twin=float(Tr[i]), primes=n1, predicted=float(pred1), ratio=float(n1 / pred1) if pred1 else None,
                  twins=n2, twins_predicted=float(pred2), twin_ratio=float(n2 / pred2) if pred2 > 0 else None,
                  prime_run_from_0=first, mean_gap_in_n=meangap, discriminant=int(1 - 4 * c))
    keep[f'ok_{c}'] = np.packbits(okf[:2000])          # first 2000 n, for the barcode figure
    print(c, res[c], flush=True)
# integers for comparison: primes below 10^12 have mean gap about log(10^12) = 27.6; twins constant 2*C2 = 1.3203
out = dict(ranking_primes_upto=2000, refined_primes_upto=1_000_000, c_max=CMAX, M=M,
           top_by_C=[int(c) for c in topC], top_by_C_twin=[int(c) for c in topT],
           C_quick_top=[float(x) for x in np.sort(C1)[::-1][:8]], C_twin_quick_top=[float(x) for x in np.sort(T1)[::-1][:8]],
           C_quick_median=float(np.median(C1[cs % 2 == 1])), results={int(k): v for k, v in res.items()})
json.dump(out, open('polys.json', 'w'), indent=1)
np.savez_compressed('polys.npz', C1=C1.astype(np.float32), T1=T1.astype(np.float32), **keep)
