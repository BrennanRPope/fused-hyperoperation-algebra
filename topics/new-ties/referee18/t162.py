"""Theorem 16.2: distinct reduced two-letter words (not NN) have distinct fingerprints."""
import math, time
import mpmath as mp
from common import *

print('=== Theorem 16.2 ===')
# ---- D != 0 : D=0 would force r P 19^{g'} = r' P' 19^{g}; check the exponent map is injective
seen = {}
for (L, M, x) in two_letter_words(300):
    r = rval(x); g = 1 if (x == 'N' and L >= 1) else 0
    key = (vp(r, 2) + 3 * L + 2 * M, (1 if r == 10 else 0) + 2 * L + M, (1 if r == 19 else 0) - g)
    assert key not in seen, (key, seen.get(key), (L, M, x))
    seen[key] = (L, M, x)
print('D=0 impossible: (v2(rP), v5(rP), v19(r)-g) injective on words H<=300 (and by the hand argument in general)')

C4 = matveev_C(4, [1.0, LOG2, LOG5, LOG19])  # times log K
print('Matveev n=4 constant (without the log K factor): %.6e' % C4)

def HS_of(H):  # short word: 2 H_S <= v2(rP) <= 2 + log2|D|, |D| <= H (long word H >= H_S)
    return int((2 + math.log2(H)) // 2)

def logK_of(HS):
    return math.log(20.9) + HS * math.log(200)

def f(H):
    return H * LOG10 + LOG10 - C4 * logK_of(HS_of(H)) * (1 + math.log(3 * H + 1))

lo, hi = 2.0, 1e30
for _ in range(400):
    mid = math.sqrt(lo * hi)
    if f(mid) > 0: hi = mid
    else: lo = mid
# f is not monotone (step in HS); scan upward from hi to be safe
H0 = int(hi)
x = H0
while x < 1e30:
    assert f(x) > 0, x
    x *= 1.01
for j in range(1, 120):   # f drops only where HS_of jumps (powers of 2): check there
    if 2 ** j > H0: assert f(2 ** j) > 0, j
print('a priori: long word H <= %d (%.4e), short word H_S <= %d, log K <= %.1f'
      % (H0, H0, HS_of(H0), logK_of(HS_of(H0))))

Hb = H0
rnd = 0
while True:
    rnd += 1
    vmax = 2 + int(math.floor(math.log2(max(Hb, 1))))   # v2(r_S P_S) <= 2 + v2(D), |D|<=Hb
    shorts = [w for w in two_letter_words(vmax) if vp(rval(w[2]), 2) + 3 * w[0] + 2 * w[1] <= vmax]
    X0 = 3 * Hb + 1
    t = time.time()
    best = 0; Cexps = set(); minlogc = 1e9
    for (L, M, x) in shorts:
        K, g, T = kappa2(L, M, x)
        assert K > 1
        fn = (lambda K=K: [mp.log(K), mp.log(2), mp.log(5), mp.log(19)])
        nb, info = reduce_bound(fn, 4, X0, math.log(0.1), LOG10)
        best = max(best, nb); Cexps.add(info['Cexp']); minlogc = min(minlogc, info['log10c'])
    print('round %d: H<=%d, X0=%d, short words with v2(rP)<=%d: %d lattices, C=10^%s, min log10 c=%.2f -> long H<=%d  (%.1fs)'
          % (rnd, Hb, X0, vmax, len(shorts), sorted(Cexps), minlogc, best, time.time() - t))
    if best >= Hb:
        break
    Hb = best
Hfinal = max(Hb, HS_of(max(Hb, 2)), 2)
fp = {}
for w in two_letter_words(Hfinal):
    k = kappa2(*w)[0]
    assert k not in fp, (w, fp[k])
    fp[k] = w
print('final bound H<=%d ; exact check: all %d two-letter words with H<=%d have distinct fingerprints'
      % (Hb, len(fp), Hfinal))
