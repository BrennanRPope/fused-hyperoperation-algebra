"""Corner-sum threshold k*_n for n named regions (plus purifier O).
f_k(S) = sum over nonempty I subset [n] of k^(|I|-1) S_I.
k*_n(cone) = largest k at which f_k < 0 somewhere on the cone (normalised s = 1).
Cones: vN = SA + elemental SSA with purification symmetry (outer bound for all quantum states);
       HOLO = vN + all MMI instances (outer bound for holographic states; exact for n <= 3).
Inner examples: star graph (perfect tensor / AME), flower graphs, random min-cut graphs, random quantum states.
"""
import itertools, json, sys
import numpy as np
from scipy.optimize import linprog

def setup(n):
    N = n + 1                      # parties 0..n-1 named, n = O
    full = (1 << N) - 1
    # canonical representative: subsets not containing O (bit n)
    reps = [m for m in range(1, 1 << n)]   # nonempty subsets of named = all canonical classes
    idx = {m: i for i, m in enumerate(reps)}
    def var(m):
        m &= full
        if m == 0 or m == full: return None
        if m >> n & 1: m = full ^ m
        return idx[m]
    return N, full, reps, idx, var

def row(var, nvar, terms):
    r = np.zeros(nvar)
    for c, m in terms:
        v = var(m)
        if v is not None: r[v] += c
    return r

def cone(n, holo=False):
    N, full, reps, idx, var = setup(n)
    nv = len(reps); rows = []
    for i, j in itertools.combinations(range(N), 2):
        rest = [x for x in range(N) if x not in (i, j)]
        for r in range(len(rest) + 1):
            for K in itertools.combinations(rest, r):
                k = sum(1 << x for x in K)
                # S_iK + S_jK - S_ijK - S_K >= 0
                rows.append(row(var, nv, [(1, (1<<i)|k), (1, (1<<j)|k), (-1, (1<<i)|(1<<j)|k), (-1, k)]))
    if holo:
        for lab in itertools.product(range(4), repeat=N):
            I = sum(1 << x for x in range(N) if lab[x] == 1)
            J = sum(1 << x for x in range(N) if lab[x] == 2)
            K = sum(1 << x for x in range(N) if lab[x] == 3)
            if not (I and J and K) or not (I < J < K): continue
            rows.append(row(var, nv, [(1, I|J), (1, I|K), (1, J|K), (-1, I), (-1, J), (-1, K), (-1, I|J|K)]))
    A = -np.array(rows)            # A x <= 0
    return A, reps, var

def fvec(n, k, reps):
    return np.array([k ** (bin(m).count('1') - 1) for m in reps])

def lp_min(n, k, A, reps):
    c = fvec(n, k, reps)
    aeq = np.array([[1.0 if bin(m).count('1') == 1 else 0.0 for m in reps]])
    res = linprog(c, A_ub=A, b_ub=np.zeros(len(A)), A_eq=aeq, b_eq=[1.0], bounds=[(0, None)] * len(reps), method='highs')
    return res.fun, res.x

def threshold(fun, lo=-1.0, hi=0.0, tol=1e-10):
    # find largest k in [lo,hi] with fun(k) < 0 ; assumes negative at lo side, >=0 at hi
    assert fun(lo) < -1e-12, 'not negative at lo'
    if fun(hi) < -1e-12: return hi
    for _ in range(60):
        mid = (lo + hi) / 2
        if fun(mid) < -1e-12: lo = mid
        else: hi = mid
    return (lo + hi) / 2

def star_poly(n):
    # perfect tensor / star graph on N=n+1 legs: S_I = min(|I|, N-|I|)
    N = n + 1
    from math import comb
    return [comb(n, j) * min(j, N - j) for j in range(1, n + 1)]   # coefficient of k^(j-1)

def star_root(n):
    co = star_poly(n)
    roots = np.roots(co[::-1])
    real = sorted(r.real for r in roots if abs(r.imag) < 1e-12 and r.real < 0)
    return (real[-1] if real else None), co

if __name__ == '__main__':
    out = {}
    for n in [int(a) for a in sys.argv[1:]] or [3, 4, 5, 6]:
        A_vn, reps, var = cone(n)
        A_h, _, _ = cone(n, holo=True)
        kv = threshold(lambda k: lp_min(n, k, A_vn, reps)[0])
        kh = threshold(lambda k: lp_min(n, k, A_h, reps)[0])
        ks, co = star_root(n)
        # minimiser at the threshold side
        fun, x = lp_min(n, kv - 1e-4, A_vn, reps)
        xs = {format(m, f'0{n}b')[::-1]: round(float(v / max(x)), 4) for m, v in zip(reps, x)}
        out[n] = dict(vN=kv, holo_outer=kh, star=ks, star_poly=co, vN_minimiser=xs)
        print(n, 'vN-LP', kv, 'holo-LP(SA+MMI)', kh, 'star', ks, 'poly', co)
    json.dump(out, open('thresholds.json', 'w'), indent=1)
