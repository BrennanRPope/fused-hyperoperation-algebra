"""Positive braid words modulo the braid relations, counted two ways: Burau matrices at a generic
rational t (distinct products), and Deligne's growth series 1/sum_T (-1)^|T| t^len(Delta_T)."""
import itertools, sympy as sp, json, math
from fractions import Fraction as Fr
t = sp.symbols('t')

def burau(i, n, q):
    M = [[Fr(int(r == c)) for c in range(n)] for r in range(n)]
    M[i][i] = 1 - q; M[i][i + 1] = q; M[i + 1][i] = Fr(1); M[i + 1][i + 1] = Fr(0)
    return M

def mul(A, B):
    n = len(A)
    return [[sum(A[r][k] * B[k][c] for k in range(n)) for c in range(n)] for r in range(n)]

def count(n, L, q=Fr(3, 7)):
    gens = [burau(i, n, q) for i in range(n - 1)]
    I = [[Fr(int(r == c)) for c in range(n)] for r in range(n)]
    level = {tuple(map(tuple, I)): I}; out = []
    for _ in range(L):
        nxt = {}
        for A in level.values():
            for G in gens:
                B = mul(A, G); nxt[tuple(map(tuple, B))] = B
        level = nxt; out.append(len(level))
    return out

def deligne_mu(n):
    """Subsets T of the n-1 generators; Delta_T length = sum over connected blocks of k(k+1)/2."""
    mu = 0
    for k in range(n):
        for T in itertools.combinations(range(n - 1), k):
            blocks, run = [], 0
            prev = None
            for s in T:
                if prev is not None and s == prev + 1: run += 1
                else:
                    if prev is not None: blocks.append(run)
                    run = 1
                prev = s
            if prev is not None: blocks.append(run)
            mu += (-1) ** k * t ** sum(b * (b + 1) // 2 for b in blocks)
    return sp.expand(mu)

def series(mu, L):
    s = sp.series(1 / mu, t, 0, L + 1).removeO()
    return [int(s.coeff(t, k)) for k in range(1, L + 1)]

def growth(mu):
    rs = [complex(r) for r in sp.Poly(mu, t).nroots()]
    return 1 / min(r.real for r in rs if abs(r.imag) < 1e-12 and r.real > 0)

if __name__ == '__main__':
    out = {}
    for n, L in [(3, 10), (4, 8), (5, 6)]:
        mu = deligne_mu(n); s = series(mu, L); c = count(n, L)
        g = growth(mu)
        print('B%d+ mu=%s  burau %s  deligne %s match=%s growth=%.6f entropy=%.6f' % (n, mu, c, s, c == s, g, math.log(g)))
        out['B%d' % n] = dict(mu=str(mu), burau=c, deligne=s, match=c == s, growth=g)
    json.dump(out, open('posbraid.json', 'w'), indent=1)
