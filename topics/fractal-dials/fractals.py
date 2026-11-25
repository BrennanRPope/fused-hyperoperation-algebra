"""Shared definitions: self-similar sets in [0,1] given as lists of (ratio, offset)."""
import math
from fractions import Fraction as F

# Five sets of dimension 1/2 (sqrt-ratios add to 1), and two of dimension log3/log2 (in the plane).
HALF = {
    "A": dict(name="two pieces at 1/4", maps=[(F(1,4), F(0)), (F(1,4), F(3,4))]),
    "B": dict(name="three at 1/9, evenly spaced", maps=[(F(1,9), F(0)), (F(1,9), F(4,9)), (F(1,9), F(8,9))]),
    "C": dict(name="three at 1/9, middle off-centre", maps=[(F(1,9), F(0)), (F(1,9), F(2,9)), (F(1,9), F(8,9))]),
    "D": dict(name="uneven: 0.09 and 0.49", maps=[(F(9,100), F(0)), (F(49,100), F(51,100))]),
    "E": dict(name="four at 1/16", maps=[(F(1,16), F(0)), (F(1,16), F(5,16)), (F(1,16), F(10,16)), (F(1,16), F(15,16))]),
}

def moran_sum(maps, s):
    return sum(float(r) ** s for r, _ in maps)

def dimension(maps, lo=0.0, hi=1.0):
    for _ in range(200):
        mid = (lo + hi) / 2
        if moran_sum(maps, mid) > 1: lo = mid
        else: hi = mid
    return (lo + hi) / 2

def intervals(maps, level):
    iv = [(F(0), F(1))]
    for _ in range(level):
        iv = [(a + (b - a) * o, a + (b - a) * o + (b - a) * r) for (a, b) in iv for (r, o) in maps]
    return iv

def first_gaps(maps):
    m = sorted(maps, key=lambda t: t[1])
    g = []
    for (r1, o1), (r2, o2) in zip(m, m[1:]):
        g.append(o2 - (o1 + r1))
    return [x for x in g if x > 0]

def gaps_above(maps, thresh):
    """all gap lengths >= thresh (floats), via the scale tree."""
    g0 = [float(x) for x in first_gaps(maps)]
    out = []
    stack = [1.0]
    gmax = max(g0)
    while stack:
        rho = stack.pop()
        if rho * gmax < thresh: continue
        out += [rho * x for x in g0 if rho * x >= thresh]
        stack += [rho * float(r) for r, _ in maps]
    return out

def tube_ratio(maps, d, eps):
    """|K_eps| / eps^(1-d) for K in [0,1] containing 0 and 1, measure zero."""
    big = gaps_above(maps, 2 * eps)
    vol = 2 * eps + (1 - sum(big)) + 2 * eps * len(big)
    return vol / eps ** (1 - d)

def is_lattice(maps, tol=1e-9, qmax=40):
    logs = [math.log(1 / float(r)) for r, _ in maps]
    base = logs[0]
    for x in logs[1:]:
        ratio = x / base
        if not any(abs(ratio - p / q) < tol for q in range(1, qmax) for p in range(1, 4 * qmax)):
            return False
    return True
