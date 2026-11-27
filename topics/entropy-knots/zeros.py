"""Trace monoids as hard-core gases: the Cartier-Foata polynomial mu(t) is the independence
polynomial of the dependence graph at fugacity -t, and the entropy is -log of its nearest zero.
Paths and cycles put the zeros at Beraha numbers 4cos^2(pi/n), the Jones index values."""
import itertools, math, json
import numpy as np, sympy as sp
from fractions import Fraction as Fr

t, d, x = sp.symbols('t d x')

def indep_poly(n, edges):
    """Independence polynomial sum_S x^|S| over independent sets of the graph."""
    E = set(map(frozenset, edges)); p = 0
    for k in range(n + 1):
        cnt = sum(1 for S in itertools.combinations(range(n), k)
                  if all(frozenset(e) not in E for e in itertools.combinations(S, 2)))
        p += cnt * x ** k
    return sp.expand(p)

def mu_of(n, dep_edges):
    return sp.expand(indep_poly(n, dep_edges).subs(x, -t))

def growth(mu):
    roots = [complex(r) for r in sp.Poly(mu, t).nroots(n=30)]
    pos = [r.real for r in roots if abs(r.imag) < 1e-12 and r.real > 0]
    rmin = min(abs(r) for r in roots)
    return 1 / min(pos), rmin, roots

path = lambda m: [(i, i + 1) for i in range(m - 1)]
cycle = lambda m: [(i, (i + 1) % m) for i in range(m)]
def multipartite(sizes):
    lab = [k for k, s in enumerate(sizes) for _ in range(s)]
    n = len(lab)
    return n, [(i, j) for i in range(n) for j in range(i + 1, n) if lab[i] != lab[j]]

def beraha(n): return 4 * math.cos(math.pi / n) ** 2

def nearest_beraha(g):
    ns = range(3, 400); b = min(ns, key=lambda n: abs(beraha(n) - g))
    return b, beraha(b)

def cheb_U(m):  # U_m(d/2) by the Jones-Wenzl / Chebyshev recursion D_m = d D_{m-1} - D_{m-2}
    a, b = sp.Integer(1), d
    if m == 0: return a
    for _ in range(m - 1): a, b = b, sp.expand(d * b - a)
    return b

# ---- strand maps: affine maps g_i acting on coordinates (i, i+1) of Q^(m+1) ----
def strand_maps(m, seed=1):
    rng = np.random.default_rng(seed)
    primes = [7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47]
    maps = []
    for i in range(m):
        M = [[Fr(int(rng.integers(1, 6)), primes[(2 * i) % 12]), Fr(int(rng.integers(1, 6)), primes[(2 * i + 1) % 12])],
             [Fr(int(rng.integers(1, 6)), primes[(2 * i + 3) % 12]), Fr(int(rng.integers(1, 6)), primes[(2 * i + 5) % 12])]]
        v = [Fr(int(rng.integers(1, 9)), 3), Fr(int(rng.integers(1, 9)), 5)]
        maps.append((i, M, v))
    return maps

def apply_map(g, A, b, m):
    """Compose: new = g o (A,b), where (A,b): y = A y0 + b on Q^(m+1)."""
    i, M, v = g
    A = [row[:] for row in A]; b = b[:]
    r0, r1 = A[i], A[i + 1]
    A[i] = [M[0][0] * p + M[0][1] * q for p, q in zip(r0, r1)]
    A[i + 1] = [M[1][0] * p + M[1][1] * q for p, q in zip(r0, r1)]
    b0, b1 = b[i], b[i + 1]
    b[i] = M[0][0] * b0 + M[0][1] * b1 + v[0]
    b[i + 1] = M[1][0] * b0 + M[1][1] * b1 + v[1]
    return A, b

def count_composites(maps, m, L):
    N = m + 1
    I = [[Fr(int(r == c)) for c in range(N)] for r in range(N)]
    level = {(tuple(map(tuple, I)), tuple([Fr(0)] * N)): (I, [Fr(0)] * N)}
    counts = []
    for _ in range(L):
        nxt = {}
        for A, b in level.values():
            for g in maps:
                A2, b2 = apply_map(g, A, b, m)
                key = (tuple(map(tuple, A2)), tuple(b2))
                nxt[key] = (A2, b2)
        level = nxt; counts.append(len(level))
    return counts

def series(mu, L):
    s = sp.series(1 / mu, t, 0, L + 1).removeO()
    return [int(s.coeff(t, k)) for k in range(1, L + 1)]

if __name__ == '__main__':
    out = {}
    print('Paths (dependence graph P_m = Temperley-Lieb generators):')
    for m in range(1, 9):
        mu = mu_of(m, path(m)); g, rmin, roots = growth(mu)
        ident = sp.expand(d ** (m + 1) * mu.subs(t, 1 / d ** 2) - cheb_U(m + 1)) == 0
        print(' m=%d mu=%s growth=%.6f B_%d=%.6f cheb_identity=%s all_real=%s' % (
            m, mu, g, m + 2, beraha(m + 2), ident, all(abs(r.imag) < 1e-12 for r in roots)))
        out['P%d' % m] = dict(mu=str(mu), growth=g, beraha=beraha(m + 2), identity=ident,
                              zeros=sorted(r.real for r in roots))
    print('Cycles:')
    for m in range(3, 9):
        mu = mu_of(m, cycle(m)); g, _, roots = growth(mu)
        print(' C%d growth=%.6f B_%d=%.6f' % (m, g, 2 * m, beraha(2 * m)))
        out['C%d' % m] = dict(mu=str(mu), growth=g, beraha=beraha(2 * m))
    print('Affine centre patterns (complete multipartite dependence graphs), growth < 4:')
    pats = []
    for n in range(2, 9):
        for sizes in sp.utilities.iterables.partitions(n):
            sz = sorted([k for k, v in sizes.items() for _ in range(v)], reverse=True)
            if len(sz) < 2: continue
            N, E = multipartite(sz); mu = mu_of(N, E)
            g, _, roots = growth(mu)
            if g < 4 - 1e-9:
                bn, bv = nearest_beraha(g)
                real = all(abs(r.imag) < 1e-12 for r in roots)
                pats.append(dict(sizes=sz, mu=str(mu), growth=g, beraha_n=bn, is_beraha=abs(bv - g) < 1e-9, real_zeros=real,
                                 zeros=[[r.real, r.imag] for r in roots]))
                print('  classes %s mu=%s growth=%.6f %s real_zeros=%s' % (sz, mu, g,
                      ('= B_%d' % bn) if abs(bv - g) < 1e-9 else '(not Beraha)', real))
    out['multipartite'] = pats
    print('Strand maps g_i on coordinates (i, i+1): distinct composites vs trace count 1/mu:')
    sm = {}
    for m, L in [(2, 8), (3, 8), (4, 7), (5, 6), (6, 6)]:
        maps = strand_maps(m)
        c = count_composites(maps, m, L)
        s = series(mu_of(m, path(m)), L)
        print('  m=%d composites %s trace %s match=%s' % (m, c, s, c == s))
        sm[m] = dict(composites=c, trace=s, match=c == s)
    out['strand_maps'] = sm
    json.dump(out, open('zeros.json', 'w'), indent=1, default=str)
