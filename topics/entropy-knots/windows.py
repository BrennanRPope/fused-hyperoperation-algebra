"""Superstable windows of the logistic map r x (1 - x), their orbit patterns (the 1D braid type),
the entropy each pattern forces, and the map's entropy from the kneading determinant."""
import mpmath as mp, numpy as np, json, itertools
mp.mp.dps = 60

def F(r, x): return r * x * (1 - x)

def orbit(r, n):
    x = mp.mpf(1) / 2; out = [x]
    for _ in range(n - 1):
        x = F(r, x); out.append(x)
    return out

def g(r, n):
    x = mp.mpf(1) / 2
    for _ in range(n): x = F(r, x)
    return x - mp.mpf(1) / 2

def exact_period(r, n):
    pts = orbit(r, n)
    for m in range(1, n):
        if n % m == 0 and abs(g(r, m)) < mp.mpf(10) ** -30: return False
    return True

def superstable(nmax, rlo=2.0, rhi=4.0, grid=200000):
    found = []
    rs = np.linspace(rlo + 1e-9, rhi - 1e-12, grid)
    for n in range(2, nmax + 1):
        vals = []
        x = np.full_like(rs, 0.5)
        for _ in range(n): x = rs * x * (1 - x)
        v = x - 0.5
        idx = np.where(np.sign(v[:-1]) != np.sign(v[1:]))[0]
        for i in idx:
            try:
                r = mp.findroot(lambda rr: g(rr, n), (mp.mpf(rs[i]), mp.mpf(rs[i + 1])), solver='anderson')
            except Exception:
                continue
            if not (rlo < r < rhi): continue
            if abs(g(r, n)) > mp.mpf(10) ** -40: continue
            if exact_period(r, n) and all(abs(r - f[1]) > mp.mpf(10) ** -20 for f in found if f[0] == n):
                found.append((n, r))
    return sorted(found, key=lambda t: t[1])

def pattern(r, n):
    pts = orbit(r, n)          # pts[0] = 1/2
    order = sorted(range(n), key=lambda i: pts[i])
    rank = {i: k for k, i in enumerate(order)}
    # cyclic permutation on positions 0..n-1: position of point i -> position of point i+1
    perm = [0] * n
    for i in range(n): perm[rank[i]] = rank[(i + 1) % n]
    return perm

def pattern_matrix(perm):
    """Markov matrix of the 'connect-the-dots' (P-linear) map on the n-1 gaps."""
    n = len(perm); A = np.zeros((n - 1, n - 1), dtype=int)
    for j in range(n - 1):
        a, b = sorted((perm[j], perm[j + 1]))
        for k in range(a, b): A[j, k] = 1
    return A

def pattern_entropy(perm):
    A = pattern_matrix(perm)
    if A.shape[0] == 0: return 0.0
    lam = max(abs(np.linalg.eigvals(A.astype(float))))
    return float(np.log(max(lam, 1.0)))

def kneading_entropy(r, N=70):
    """Entropy from the smallest zero of the Milnor-Thurston kneading determinant."""
    c = mp.mpf(1) / 2; x = F(r, c); theta = [1]; s = 1
    coeffs = [1]
    for _ in range(N):
        eps = 1 if x < c else -1          # increasing branch left of c
        s *= eps; coeffs.append(s); x = F(r, x)
    p = np.polynomial.Polynomial([float(cf) for cf in coeffs])
    roots = p.roots()
    real = [z.real for z in roots if abs(z.imag) < 1e-9 and 0 < z.real < 1 + 1e-9]
    if not real: return 0.0, coeffs
    return float(-np.log(min(real))), coeffs

def kneading_exact(r, n, side):
    """Exact kneading determinant at a superstable parameter: the critical orbit has period n;
    the last symbol (the critical point itself) is replaced by side = +1 or -1 (one-sided limits).
    D(t) = (sum_{k<n} theta_k t^k) / (1 - S t^n); returns (entropy, numerator coefficients)."""
    c = mp.mpf(1) / 2; x = F(r, c); eps = []
    for k in range(1, n):
        eps.append(1 if x < c else -1); x = F(r, x)
    eps.append(side)
    theta = [1]
    for e in eps[:-1]: theta.append(theta[-1] * e)
    p = np.polynomial.Polynomial([float(v) for v in theta])
    real = [z.real for z in p.roots() if abs(z.imag) < 1e-9 and 0 < z.real < 1 - 1e-4]
    return (float(-np.log(min(real))) if real else 0.0), theta

if __name__ == '__main__':
    W = superstable(9)
    rows = []
    for n, r in W:
        perm = pattern(r, n)
        hp = pattern_entropy(perm)
        hk_lo, th = kneading_exact(r, n, +1)
        hk_hi, _ = kneading_exact(r, n, -1)
        rows.append(dict(n=n, r=float(r), perm=perm, h_pattern=hp, h_map_below=hk_lo, h_map_above=hk_hi, kneading=th))
    json.dump(rows, open('windows.json', 'w'), indent=1)
    for row in rows:
        print(row['n'], round(row['r'], 6), row['perm'], round(row['h_pattern'], 6), round(row['h_map_below'], 6), round(row['h_map_above'], 6))
    print(len(rows), 'windows')
