"""Shared helpers for checks2.py (same definitions as checks.py)."""
import itertools, math
import numpy as np

def subsets(parties):
    for r in range(1, len(parties) + 1):
        for T in itertools.combinations(parties, r):
            yield frozenset(T)

def corner_sum(S, parties, k):
    return sum(k ** (len(T) - 1) * S[T] for T in subsets(parties))

def occ(terms, party):
    return tuple(1 if party in t else 0 for t in terms)

def find_contraction(lhs, rhs, parties, count_all=False, limit=10**7):
    """f: {0,1}^|lhs| -> {0,1}^|rhs|, plain Hamming contraction (repeated terms = repeated bits),
    boundary conditions for each party and the purifier (all zeros)."""
    L, R = len(lhs), len(rhs)
    dH = lambda x, y: sum(a != b for a, b in zip(x, y))
    fixed = {tuple([0] * L): tuple([0] * R)}
    for p in parties:
        x, y = occ(lhs, p), occ(rhs, p)
        if x in fixed and fixed[x] != y:
            return None, 0
        fixed[x] = y
    for (x1, y1), (x2, y2) in itertools.combinations(list(fixed.items()), 2):
        if dH(y1, y2) > dH(x1, x2):
            return None, 0
    free = [x for x in itertools.product((0, 1), repeat=L) if x not in fixed]
    free.sort(key=lambda x: min(dH(x, b) for b in fixed))
    codomain = list(itertools.product((0, 1), repeat=R))
    assign = dict(fixed); found = [None]; count = [0]; nodes = [0]
    def rec(i):
        nodes[0] += 1
        if nodes[0] > limit: return True
        if i == len(free):
            count[0] += 1
            if found[0] is None: found[0] = dict(assign)
            return not count_all
        x = free[i]
        for y in codomain:
            if all(dH(y, yy) <= dH(x, xx) for xx, yy in assign.items()):
                assign[x] = y
                if rec(i + 1): return True
                del assign[x]
        return False
    rec(0)
    return found[0], count[0]

def _spectrum(psi, dims, keep):
    n = len(dims); t = psi.reshape(dims)
    keep = sorted(keep); rest = [i for i in range(n) if i not in keep]
    M = np.transpose(t, keep + rest).reshape(int(np.prod([dims[i] for i in keep])), -1)
    s = np.linalg.svd(M, compute_uv=False) ** 2
    return s[s > 1e-14]

def ent(psi, dims, keep):
    s = _spectrum(psi, dims, keep); return float(-(s * np.log2(s)).sum())

def renyi(psi, dims, keep, alpha):
    s = _spectrum(psi, dims, keep)
    if abs(alpha - 1) < 1e-12: return float(-(s * np.log2(s)).sum())
    if alpha > 1e5: return float(-np.log2(s.max()))
    return float(np.log2((s ** alpha).sum()) / (1 - alpha))

def I3_state(psi, dims):
    idx = {"A": 0, "B": 1, "C": 2}
    S = {T: ent(psi, dims, [idx[p] for p in T]) for T in subsets("ABC")}
    return corner_sum(S, "ABC", -1)

def haar(rng, dim):
    v = rng.normal(size=dim) + 1j * rng.normal(size=dim); return v / np.linalg.norm(v)

def sparse_state(rng):
    nterms = int(rng.integers(2, 5))
    v = np.zeros(16, complex); idx = rng.choice(16, nterms, replace=False)
    v[idx] = rng.normal(size=nterms) + 1j * rng.normal(size=nterms)
    return v / np.linalg.norm(v)

def random_graph(rng, max_bulk=7):
    nbk = int(rng.integers(1, max_bulk + 1)); N = 4 + nbk
    W = rng.exponential(1.0, (N, N)) * (rng.random((N, N)) < 0.6)
    W = np.triu(W, 1); W = W + W.T
    return nbk, W

def _cuts(nbk, W, T, boundary=("A", "B", "C", "O")):
    nb = len(boundary)
    side0 = np.array([1 if b in T else 0 for b in boundary] + [0] * nbk)
    out = []
    for m in itertools.product((0, 1), repeat=nbk):
        side = side0.copy(); side[nb:] = m
        out.append(W[np.ix_(side == 1, side == 0)].sum())
    return np.array(out)

def mincut_entropies(nbk, W):
    return {T: float(_cuts(nbk, W, T).min()) for T in subsets("ABC")}

def softmin_entropies(nbk, W, h):
    """S_h(T) = -h log(Z_T/Z_0), Z_T = sum over bulk configurations of exp(-cut/h).
    h -> 0 gives the min cut (Maslov's tropical limit). For h = 1/log D this is the
    large-D random tensor network Renyi-2 entropy in units of log D (Hayden et al. 2016)."""
    def F(T):
        c = _cuts(nbk, W, T); c0 = c.min()
        return c0 - h * math.log(np.exp(-(c - c0) / h).sum())
    F0 = F(frozenset())
    return {T: F(T) - F0 for T in subsets("ABC")}
