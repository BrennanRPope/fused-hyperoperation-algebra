"""Check the three-chart gluing of the projective plane P^2 by computation.

A point [x:y:z] is stored normalized: ('A', a, b) is [a:b:1], ('L', u) is
[1:u:0], ('N',) is [0:1:0]. Charts follow the fan with rays e1, e2, -e1-e2:
c0 sees (x/z, y/z), c1 sees (z/x, y/x), c2 sees (z/y, x/y). Same definitions
as ProjectivePlane.lean.
"""
import itertools
from math import gcd

def field(p):
    inv = lambda x: 0 if x % p == 0 else pow(x, p - 2, p)
    return inv

def make(p):
    inv = field(p)
    def chart(c, v):
        x, y = v
        if c == 0: return ('A', x, y)
        if c == 1: return ('L', y) if x == 0 else ('A', inv(x), y * inv(x) % p)
        if x == 0: return ('N',) if y == 0 else ('L', inv(y))
        return ('A', y * inv(x) % p, inv(x))
    T = {(0, 1): lambda a, b: (inv(a), inv(a) * b % p),
         (0, 2): lambda a, b: (inv(b), a * inv(b) % p),
         (1, 0): lambda v, u: (inv(v), u * inv(v) % p),
         (1, 2): lambda v, u: (v * inv(u) % p, inv(u)),
         (2, 0): lambda t, s: (s * inv(t) % p, inv(t)),
         (2, 1): lambda t, s: (t * inv(s) % p, inv(s))}
    def trans(c, d, v): return v if c == d else T[(c, d)](*v)
    def wind(n, P):
        if P[0] == 'A': return ('A', pow(P[1], n, p), pow(P[2], n, p))
        if P[0] == 'L': return ('L', pow(P[1], n, p))
        return P
    return chart, trans, wind, inv

def checks(p):
    chart, trans, wind, inv = make(p)
    F = range(p); V = list(itertools.product(F, repeat=2)); a = p - 1
    pts = {('A', x, y) for x, y in V} | {('L', u) for u in F} | {('N',)}
    assert len(pts) == p * p + p + 1
    hits = {P: [] for P in pts}
    for c in range(3):
        img = [chart(c, v) for v in V]
        assert len(set(img)) == p * p, 'chart not injective'
        for v, P in zip(V, img): hits[P].append(c)
    assert all(hits.values()), 'charts miss a point'
    by = {}
    for P, h in hits.items(): by[len(h)] = by.get(len(h), 0) + 1
    assert by == {3: a * a, 2: 3 * a, 1: 3}, by
    # transitions: wherever chart c's point is also in chart d, trans moves the coordinates
    for c, d in itertools.product(range(3), repeat=2):
        for v in V:
            if d in hits[chart(c, v)]:
                assert chart(d, trans(c, d, v)) == chart(c, v)
    # cocycle on every triple overlap
    for c, d, e in itertools.product(range(3), repeat=3):
        for v in V:
            h = hits[chart(c, v)]
            if d in h and e in h:
                assert trans(d, e, trans(c, d, v)) == trans(c, e, v)
    # mult-add form of c0 -> c1: (1+u, 1+w) -> (1+u~, 1+(w (+) u~)), u~ = -u/(1+u)
    for u, w in V:
        if (1 + u) % p == 0: continue
        ui = (-u) * inv(1 + u) % p
        ma = (w + ui + w * ui) % p
        assert trans(0, 1, ((1 + u) % p, (1 + w) % p)) == ((1 + ui) % p, (1 + ma) % p)
    # windings commute with every chart change and agree on P^2
    for n in range(1, 2 * p + 2):
        pw = lambda v: (pow(v[0], n, p), pow(v[1], n, p))
        for c, d in itertools.product(range(3), repeat=2):
            for v in V:
                assert trans(c, d, pw(v)) == pw(trans(c, d, v))
            for v in V:
                assert chart(c, pw(v)) == wind(n, chart(c, v))
    # fixed points of the winding by n: g^2 + 3g + 3 with g = gcd(n-1, p-1)
    fixed = {}
    for n in range(2, 2 * p + 2):
        f = sum(1 for P in pts if wind(n, P) == P)
        g = gcd(n - 1, p - 1)
        assert f == g * g + 3 * g + 3, (p, n, f, g)
        fixed[n] = f
    assert fixed[p] == p * p + p + 1   # the winding by p is the Frobenius: it fixes every point
    return by, fixed

if __name__ == '__main__':
    for p in [2, 3, 5, 7, 11, 13]:
        by, fixed = checks(p)
        print('F_%d:' % p, 'points held by 3/2/1 charts', by[3], by.get(2, 0), by[1],
              'total', p * p + p + 1, '| fixed points of winding n=2..%d:' % (p + 1),
              [fixed[n] for n in range(2, p + 2)])
    print('all P^2 gluing checks passed')
