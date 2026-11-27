"""Sanity tests of the engine against brute force."""
import random, math
from fractions import Fraction as F
from engine import *
random.seed(1)
# psi via affine maps (independent of Tval)
MAPS = {'l': (F(19, 200), F(0)), 'm': (F(1, 20), F(0)), 'N': (F(1, 19), F(1, 19)), 'K': (F(1, 10), F(1, 10))}
def psi(word):
    s, t = F(1), F(0)
    for ch in word:
        sx, tx = MAPS[ch]; s, t = s * sx, s * tx + t
    return t / s
def wstr(items):
    return 'N' + ''.join('l' * L + 'm' * M + y for (L, M, y) in items)
def rand_item(H=6):
    L = random.randint(0, H); M = random.randint(0, H); y = random.choice('NK')
    if y == 'K': M = 0
    return (L, M, y)
for _ in range(2000):
    n = random.randint(1, 4)
    its = [rand_item() for _ in range(n)]
    A = sum(i[0] for i in its)
    assert psi(wstr(its)) == F(Tval(its), 19 ** A)
    # eps bounds: T = T_a R_a (1+eps), eps < 10^-H_{a+1}/9
    for a in range(n):
        Ta = Tval(its[:a]); R = 1
        for it in its[a:]: R *= rP(it)
        eps = F(Tval(its), Ta * R) - 1
        assert 0 < eps < F(1, 9 * 10 ** gap(its[a]))
    # 2-adic tail congruence: T = 19^{A_b} tau mod 2^{v2(R_b)+g2_b}
    for b in range(1, n + 1):
        tail = its[b:]
        tau = Tval(tail); R = 1
        for it in tail: R *= rP(it)
        G = vp(R, 2) + g2(its[b - 1]) if R % 2 == 0 or g2(its[b - 1]) else 0
        Ab = sum(i[0] for i in its[:b])
        G = sum(g2(it) for it in its[b - 1:])
        assert (Tval(its) - 19 ** Ab * tau) % (1 << G) == 0
print('psi/T, eps and 2-adic tail congruence: ok')
# vmax2 vs brute force
for _ in range(300):
    Q = random.choice([1, 5, 25, 3, 7, 11, 13, 21, 361 * 5, 19 ** 3]) * random.choice([1, -1])
    Q = Q if Q % 2 else Q + 1
    Mb = random.randint(1, 300)
    s, t = st2(Q % (1 << 60) if Q > 0 else Q)
    exc = None
    k = is_19power(F(abs(Q))) if Q > 0 else None
    if k is not None: exc = -k
    vm = vmax2(s, t, Mb, exclude=exc)
    best = 0
    for m in range(-Mb, Mb + 1):
        if m == exc: continue
        x = F(19) ** m * Q - 1
        if x == 0: continue
        best = max(best, vp(x.numerator, 2) - vp(x.denominator, 2))
    assert vm >= best, (Q, Mb, vm, best)
print('vmax2 >= brute-force max: ok')
# LLL sanity: 3-log form small values
X1, e = lll3(100)
best = 99
for b2 in range(-100, 101):
    for b5 in range(-100, 101):
        for b19 in range(-12, 13):
            if b2 == b5 == b19 == 0: continue
            v = abs(b2 * math.log(2) + b5 * math.log(5) + b19 * math.log(19))
            if v > 1e-300:
                h = math.log10(1 / (9 * v))
                assert h <= X1 + 1e-9, (b2, b5, b19, h, X1)
print('lll3(100) = %d, consistent with brute force' % X1)
