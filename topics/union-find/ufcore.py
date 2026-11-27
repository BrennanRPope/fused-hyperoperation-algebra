"""Shared code for topics/union-find: Ackermann and the hyperoperation ladder, union-find, and the
dial gas Z(z) = sum_S z^|S| k^e(S) built by Asano merges (the Lee-Yang circle theorem)."""
import random
from fractions import Fraction
import numpy as np


# ---------------------------------------------------------------- Ackermann and the ladder
def H(n, a, b):
    """Goodstein's hyperoperation: H0 = successor, H1 = +, H2 = x, H3 = ^, H4 = tetration, ..."""
    if n == 0:
        return b + 1
    if n == 1:
        return a + b
    if n == 2:
        return a * b
    if n == 3:
        return a ** b
    r = 1
    for _ in range(b):
        r = H(n - 1, a, r)
    return r


def ack(n, m):
    """Ackermann-Peter function, iterative on small arguments."""
    if n == 0:
        return m + 1
    if n == 1:
        return m + 2
    if n == 2:
        return 2 * m + 3
    if n == 3:
        return 2 ** (m + 3) - 3
    r = ack(n - 1, 1)
    for _ in range(m):
        r = ack(n - 1, r)
    return r


def fused(n, a, b):
    """rung n of the fused ladder: H_n read through x -> x + 1 (hyperoperation-ladder page)."""
    return H(n, a + 1, b + 1) - 1


def clog2(n):
    """ceil(log2 n) for integers n >= 1."""
    return (n - 1).bit_length()


def alpha_k(k, n):
    """the inverse-rung hierarchy: alpha_1(n) = ceil(n/2) undoes doubling, alpha_2 = ceil(log2) undoes 2^x,
    alpha_3 = log* undoes towers; alpha_k counts applications of alpha_{k-1} needed to bring n down to 1."""
    if k == 1:
        return -(-n // 2)
    if k == 2:
        return clog2(n)
    c = 0
    while n > 1:
        n = alpha_k(k - 1, n); c += 1
    return c


def alpha(n):
    """inverse Ackermann: the first rung whose inverse brings n down to 3 or less."""
    k = 1
    while alpha_k(k, n) > 3:
        k += 1
    return k


def tower(m):
    r = 1
    for _ in range(m):
        r = 2 ** r
    return r


# ---------------------------------------------------------------- union-find
class UF:
    def __init__(self, n, rank=True, compress=True):
        self.p = list(range(n)); self.r = [0] * n
        self.rank, self.compress = rank, compress
        self.steps = 0

    def find(self, x):
        root = x
        while self.p[root] != root:
            root = self.p[root]; self.steps += 1
        if self.compress:
            while self.p[x] != root:
                self.p[x], x = root, self.p[x]
        return root

    def union(self, x, y):
        a, b = self.find(x), self.find(y)
        if a == b:
            return a
        if self.rank:
            if self.r[a] < self.r[b]:
                a, b = b, a
            self.p[b] = a
            if self.r[a] == self.r[b]:
                self.r[a] += 1
            return a
        self.p[a] = b          # naive: hang the first root under the second
        return b

    def blocks(self):
        g = {}
        for v in range(len(self.p)):
            g.setdefault(self.find(v), []).append(v)
        return sorted(tuple(sorted(b)) for b in g.values())

    def depth(self):
        def d(x):
            c = 0
            while self.p[x] != x:
                x = self.p[x]; c += 1
            return c
        return max(d(v) for v in range(len(self.p)))


def workload(n, rank, compress, seed=0, adversarial=False):
    """average parent-pointer steps per operation."""
    uf = UF(n, rank, compress)
    rnd = random.Random(seed)
    ops = 0
    if adversarial:
        for i in range(1, n):              # chain: old root under new element
            uf.union(i - 1, i); ops += 1
        for _ in range(n):
            uf.find(0); ops += 1
    else:
        for _ in range(2 * n):
            x, y = rnd.randrange(n), rnd.randrange(n)
            if rnd.random() < 0.5:
                uf.union(x, y)
            else:
                uf.find(x)
            ops += 1
    return uf.steps / ops


# ---------------------------------------------------------------- graphs
def cycle(n):
    return n, [(i, (i + 1) % n) for i in range(n)]


def grid(L, M, torus=False):
    idx = lambda i, j: (i % L) * M + (j % M)
    E = set()
    for i in range(L):
        for j in range(M):
            if torus or i + 1 < L:
                E.add(tuple(sorted((idx(i, j), idx(i + 1, j)))))
            if torus or j + 1 < M:
                E.add(tuple(sorted((idx(i, j), idx(i, j + 1)))))
    return L * M, sorted(E)


def petersen():
    E = [(i, (i + 1) % 5) for i in range(5)] + [(5 + i, 5 + (i + 2) % 5) for i in range(5)] + [(i, i + 5) for i in range(5)]
    return 10, E


def complete(n):
    return n, [(i, j) for i in range(n) for j in range(i + 1, n)]


def random_graph(n, p, seed):
    rnd = random.Random(seed)
    return n, [(i, j) for i in range(n) for j in range(i + 1, n) if rnd.random() < p]


def degrees(N, E):
    d = [0] * N
    for u, v in E:
        d[u] += 1; d[v] += 1
    return d


# ---------------------------------------------------------------- the dial gas and its merges
def count_matrix(N, E, inner=False):
    """C[m, b] = number of vertex sets S with |S| = m and b edges leaving S (inner=True: b edges inside S)."""
    S = np.arange(1 << N, dtype=np.int64)
    m = np.zeros(1 << N, dtype=np.int64)
    for v in range(N):
        m += (S >> v) & 1
    b = np.zeros(1 << N, dtype=np.int64)
    for u, v in E:
        bu, bv = (S >> u) & 1, (S >> v) & 1
        b += (bu & bv) if inner else (bu ^ bv)
    C = np.zeros((N + 1, len(E) + 1), dtype=np.int64)
    np.add.at(C, (m, b), 1)
    return C


def components(N, E):
    p = list(range(N))
    def f(x):
        while p[x] != x:
            p[x] = p[p[x]]; x = p[x]
        return x
    for u, v in E:
        p[f(u)] = f(v)
    groups = {}
    for v in range(N):
        groups.setdefault(f(v), []).append(v)
    out = []
    for vs in groups.values():
        lab = {v: i for i, v in enumerate(vs)}
        out.append((len(vs), [(lab[u], lab[v]) for u, v in E if u in lab]))
    return out


def gas_coeffs(N, E, k, dps=60):
    """a_m (mpmath) with Z(w) = sum_m a_m w^m, each vertex at fugacity k^(-d_v/2) w (the scale Asano's proof gives):
    a_m = sum_{|S| = m} a^(edges leaving S), a = k^(-1/2)."""
    import mpmath as mp
    mp.mp.dps = dps
    a = mp.mpf(k) ** mp.mpf(-0.5)
    C = count_matrix(N, E)
    return [sum(int(C[m, b]) * a ** b for b in range(C.shape[1]) if C[m, b]) for m in range(N + 1)]


def gas_coeffs_uniform(N, E, k, dps=60):
    """coefficients of the plain dial gas Z(z) = sum_S z^|S| k^e(S)."""
    import mpmath as mp
    mp.mp.dps = dps
    kk = mp.mpf(k)
    C = count_matrix(N, E, inner=True)
    return [sum(int(C[m, e]) * kk ** e for e in range(C.shape[1]) if C[m, e]) for m in range(N + 1)]


def roots(co, dps=60):
    import mpmath as mp
    mp.mp.dps = dps
    c = list(co)
    while len(c) > 1 and c[-1] == 0:
        c.pop()
    if len(c) <= 1:
        return []
    return [complex(r) for r in mp.polyroots(c[::-1], maxsteps=800, extraprec=600)]


def gas_roots(N, E, k):
    """zeros in w of the scaled gas, component by component (a disjoint union multiplies)."""
    out = []
    for n_, E_ in components(N, E):
        out += roots(gas_coeffs(n_, E_, k))
    return out


def merge_stages(N, E, order_seed=None):
    """Asano's build as a union-find run. Start with every edge on its own two copies of its end vertices
    (a disjoint union of single edges); then union copies of the same vertex one pair at a time.
    Yields (stage graph as (n, edges), number of unions so far)."""
    copies = []                       # copy id -> original vertex
    cE = []
    for u, v in E:
        cE.append((len(copies), len(copies) + 1)); copies += [u, v]
    pairs = []
    by_v = {}
    for c, v in enumerate(copies):
        by_v.setdefault(v, []).append(c)
    for v, cs in by_v.items():
        for c in cs[1:]:
            pairs.append((cs[0], c))
    if order_seed is not None:
        random.Random(order_seed).shuffle(pairs)
    uf = UF(len(copies))

    def stage():
        roots_ = sorted({uf.find(c) for c in range(len(copies))})
        lab = {r: i for i, r in enumerate(roots_)}
        return len(roots_), [(lab[uf.find(x)], lab[uf.find(y)]) for x, y in cE]
    yield stage(), 0
    for i, (x, y) in enumerate(pairs):
        uf.union(x, y)
        yield stage(), i + 1


def gas_coeffs_exact(N, E, a):
    """exact (Fraction) coefficients of sum_S a^(edges leaving S) w^|S|."""
    C = count_matrix(N, E)
    return [sum(int(C[m, b]) * a ** b for b in range(C.shape[1])) for m in range(N + 1)]


def ack_diag_alpha(n):
    """inverse Ackermann on the diagonal: least k with A(k, k) >= n. A(3,3) = 61; A(4,4) > A(4,2) = 2^65536 - 3."""
    for k in range(4):
        if ack(k, k) >= n:
            return k
    if n <= 2 ** 65536 - 3:
        return 4
    return None
