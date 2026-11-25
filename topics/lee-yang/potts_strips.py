"""Chromatic polynomials of cyclic square-lattice strips (width L, length n, free top/bottom,
periodic along the length), by a transfer matrix on colour partitions. Exact integer polynomials in Q.
The dial reading: k = Q - 1 (Potts edge matrix trace 1 + (Q-1) t)."""
import numpy as np
from functools import lru_cache

def canon(labels):
    m = {}; out = []
    for x in labels:
        if x not in m: m[x] = len(m)
        out.append(m[x])
    return tuple(out)

def polymul_ff(p, b, r):
    # multiply poly p (list of ints, low->high in Q) by (Q-b)(Q-b-1)...(Q-b-r+1)
    for j in range(r):
        c = b + j
        q = [0] * (len(p) + 1)
        for i, a in enumerate(p):
            q[i + 1] += a; q[i] -= c * a
        p = q
    return p

def set_partitions_assign(L, nb, forbid):
    """assign labels to L new vertices in a column (vertical path), each label either an existing
    block 0..nb-1 (not in forbid[i]) or a new block nb, nb+1, ... (canonical order among new).
    adjacent vertices in the column must differ. yields (labels, number_of_new_blocks)."""
    res = []
    def rec(i, cur, nnew):
        if i == L:
            res.append((tuple(cur), nnew)); return
        opts = [x for x in range(nb) if x not in forbid[i]] + list(range(nb, nb + nnew + 1))
        for x in opts:
            if i > 0 and cur[-1] == x: continue
            rec(i + 1, cur + [x], nnew + (1 if x == nb + nnew else 0))
    rec(0, [], 0)
    return res

def cyclic_strip_chromatic(L, n):
    # state: canonical labels of (first column + current column), 2L entries
    states = {}
    for lab, nn in set_partitions_assign(L, 0, [set()] * L):
        st = canon(lab + lab)
        p = polymul_ff([1], 0, nn)
        states[st] = add(states.get(st, [0]), p)
    for _ in range(n - 1):
        new = {}
        for st, p in states.items():
            nb = max(st) + 1
            cur = st[L:]
            forbid = [{cur[i]} for i in range(L)]
            for lab, nn in set_partitions_assign(L, nb, forbid):
                q = polymul_ff(p, nb, nn)
                ns = canon(st[:L] + lab)
                new[ns] = add(new.get(ns, [0]), q)
        states = new
    total = [0]
    for st, p in states.items():
        if all(st[i] != st[L + i] for i in range(L)):
            total = add(total, p)
    return total

def add(a, b):
    n = max(len(a), len(b)); return [(a[i] if i < len(a) else 0) + (b[i] if i < len(b) else 0) for i in range(n)]

def evalp(p, x): return sum(c * x ** i for i, c in enumerate(p))

if __name__ == "__main__":
    # sanity: L=1 cyclic strip = cycle C_n: (Q-1)^n + (-1)^n (Q-1)
    from fractions import Fraction as F
    for n in range(3, 8):
        p = cyclic_strip_chromatic(1, n)
        for Q in range(0, 6):
            assert evalp(p, Q) == (Q - 1) ** n + (-1) ** n * (Q - 1), (n, Q)
    # L=2 n=4: square-lattice 2x4 cyclic = C4 x P2 prism; brute check by counting colourings with Q=3,4
    import itertools
    def brute(L, n, Q):
        V = [(i, j) for i in range(L) for j in range(n)]
        E = [((i, j), (i + 1, j)) for i in range(L - 1) for j in range(n)] + [((i, j), (i, (j + 1) % n)) for i in range(L) for j in range(n)]
        cnt = 0
        for col in itertools.product(range(Q), repeat=len(V)):
            c = dict(zip(V, col))
            if all(c[a] != c[b] for a, b in E): cnt += 1
        return cnt
    for L, n, Q in ((2, 3, 3), (2, 4, 3), (3, 3, 3), (2, 5, 3), (3, 4, 2), (2, 3, 4)):
        assert evalp(cyclic_strip_chromatic(L, n), Q) == brute(L, n, Q), (L, n, Q)
    print("strip transfer matrix matches brute-force colouring counts")
