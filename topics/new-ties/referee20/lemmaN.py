"""Lemma N (my version): checks used by the vanishing analysis.
(i)  r P 19^{L'} = rho r' P' 19^{L} has a solution in items only if v2(rho) - 2 v5(rho) = v19(rho);
     for rho = 1 the only solutions are equal items.  (brute force over L, M <= 12; the identity
     (v2 - 2 v5)(rP) = -L - [K] = v19(r) - L - 1 makes it general.)
(ii) NN-type suffix pairs: suffix N N z rest  vs  N (z+m) rest  give rho = 20/19 or 19/20
     (computed from the actual T and R values), which fails (i).
"""
import random, itertools
from fractions import Fraction as F
from engine import Tval, rP, strip, vp
its = [(L, M, y) for L in range(13) for M in range(13) for y in 'NK' if not (y == 'K' and M)]
def val(q, p):
    return vp(q.numerator, p) - vp(q.denominator, p)
for it in its:
    L, M, y = it
    assert vp(rP(it), 2) - 2 * vp(rP(it), 5) == -L - (y == 'K') == (y == 'N') - L - 1
sols = 0
for a in its:
    for b in its:
        rho = F(rP(a) * 19 ** b[0], rP(b) * 19 ** a[0])
        assert val(rho, 2) - 2 * val(rho, 5) == val(rho, 19)
        if rho == 1:
            assert a == b, (a, b)
            sols += 1
print('(i) checked on %d x %d item pairs: consistency identity holds; rho = 1 only for equal items' % (len(its), len(its)))
def plus_m(z):
    L, M, y = z
    return (L, M + 1, 'N') if y == 'N' else (L + 1, 0, 'N')
random.seed(2)
rhos = set()
for _ in range(3000):
    z = random.choice(its[:200])
    rest = [random.choice(its[:200]) for _ in range(random.randint(0, 1))]
    s = [(0, 0, 'N'), z] + rest          # items of suffix N N z rest
    s2 = [plus_m(z)] + rest               # items of N (z+m) rest
    assert strip(Tval(s)) == strip(Tval(s2))
    R = 1
    for it in s: R *= rP(it)
    R2 = 1
    for it in s2: R2 *= rP(it)
    rho = F(Tval(s) * R2, Tval(s2) * R)   # rho = T(suf) R'/(T(suf') R)
    rhos.add(rho)
    assert val(rho, 2) - 2 * val(rho, 5) != val(rho, 19)
print('(ii) NN-type suffix pairs: rho values', sorted(rhos), '-> inconsistent with (i): vanishing impossible')
