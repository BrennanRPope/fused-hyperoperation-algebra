"""sanity tests of the tools: discrete-log bound vs brute force; de Weger lemma
does not exclude a real near-relation."""
import random, math
import mpmath as mp
from common import *
random.seed(7)
# 1. max_v2_K19m against brute force
def brute(K, Mb):
    best = 0
    for m in range(-Mb, Mb + 1):
        if m >= 0: v = vp(K - 19 ** m, 2) if K != 19 ** m else 10 ** 9
        else: v = vp(K * 19 ** (-m) - 1, 2)
        best = max(best, v)
    return best
cnt = 0
for (L, M, x) in two_letter_words(25):
    K = kappa2(L, M, x)[0]
    for Mb in (5, 37, 200):
        b = brute(K, Mb); u = max_v2_K19m(K, Mb)
        assert b <= u, (K, Mb, b, u)
        cnt += 1
print('discrete-log bound >= brute-force max on %d cases' % cnt)
# 2. de Weger with a genuine near relation: K = kappa(V), Lambda = K 19^g/(rP) - 1 ~ 10^-H
for (L, M, x) in [(3, 27, 'N'), (0, 30, 'N'), (12, 0, 'K'), (20, 20, 'N')]:
    K, g, T = kappa2(L, M, x); H = L + M
    fn = lambda: [mp.log(K), mp.log(2), mp.log(5), mp.log(19)]
    nb, info = reduce_bound(fn, 4, 3 * (H + 5) + 1, math.log(0.1), LOG10)
    print('  true near-relation H=%d, lattice bound H<=%d (must be >= %d)' % (H, nb, H))
    assert nb >= H
print('de Weger implementation does not exclude genuine solutions')
