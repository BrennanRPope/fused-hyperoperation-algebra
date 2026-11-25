"""Check the four-chart gluing of P^1 x P^1 by computation.

A point of P^1 is a number or INF. Chart False on a factor sees x; chart True
sees 1/x (its x = 0 is the point at infinity). A chart of P^1 x P^1 is a pair
of choices, one per quadrant of the fan with rays +e1, -e1, +e2, -e2.
"""
import itertools, random
from fractions import Fraction
INF = 'inf'

def make(p=None):
    inv = (lambda x: pow(x, p - 2, p)) if p else (lambda x: 1 / Fraction(x))
    def chart1(b, x):
        if not b: return x
        return INF if x == 0 else inv(x)
    def chart2(c, v): return (chart1(c[0], v[0]), chart1(c[1], v[1]))
    def trans(c, d, v):
        return tuple(v[i] if c[i] == d[i] else (0 if v[i] == 0 else inv(v[i])) for i in range(2))
    return chart1, chart2, trans

CH = list(itertools.product([False, True], repeat=2))

def finite_field_checks(p):
    chart1, chart2, trans = make(p)
    F = range(p); pts = set(itertools.product(list(F) + [INF], repeat=2))
    assert len(pts) == (p + 1) ** 2
    hits = {q: 0 for q in pts}
    for c in CH:
        img = [chart2(c, v) for v in itertools.product(F, repeat=2)]
        assert len(set(img)) == p * p, 'chart not injective'
        for q in img: hits[q] += 1
    assert all(h >= 1 for h in hits.values()), 'charts miss a point'
    # a point lies in 4 charts on the torus, 2 on a circle, 1 at a fixed point
    by = {}
    for q, h in hits.items(): by[h] = by.get(h, 0) + 1
    a = p - 1
    assert by == {4: a * a, 2: 4 * a, 1: 4}, by
    # transitions: wherever two charts overlap, trans moves chart c's coords to chart d's
    for c, d in itertools.product(CH, repeat=2):
        for v in itertools.product(F, repeat=2):
            ok = all(v[i] != 0 for i in range(2) if c[i] != d[i])
            if ok: assert chart2(d, trans(c, d, v)) == chart2(c, v)
    # cocycle: going c -> d -> e equals c -> e
    for c, d, e in itertools.product(CH, repeat=3):
        for v in itertools.product(F, repeat=2):
            assert trans(d, e, trans(c, d, v)) == trans(c, e, v)
    return by

def rational_checks(n=2000):
    chart1, chart2, trans = make(None)
    rnd = random.Random(1)
    for _ in range(n):
        u = Fraction(rnd.randint(-50, 50), rnd.randint(1, 50))
        if u == -1: continue
        # gluing x -> 1/x in circle-operation coordinates (x = 1 + u) is the circle inverse
        uinv = -u / (1 + u)
        assert u + uinv + u * uinv == 0
        assert 1 / (1 + u) == 1 + uinv
    return True

if __name__ == '__main__':
    for p in [2, 3, 5, 7, 11, 13]:
        print('F_%d:' % p, 'points by number of charts holding them', finite_field_checks(p), 'total', (p + 1) ** 2)
    print('Circle inverse = gluing map on 2000 random rationals:', rational_checks())
    print('all gluing checks passed')
