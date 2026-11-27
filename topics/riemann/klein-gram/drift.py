"""Lattice count and prime count for the torus, Klein bottle and odd half.

Coefficients (eigenvalue n = m^2 + k^2 counted with multiplicity):
  torus a_T(n) = r2(n);  Klein a_K(n) = r2(n)/2 - eta_n;  odd a_O(n) = r2(n)/2 + eta_n,
  eta_n = (-1)^(j+1) if n = j^2, else 0.
Lattice count N(X) = sum_{n<=X} a(n): poles of F only, so zeros never enter it.
  Exact: N_K(X) = N_T(X)/2 - [floor(sqrt X) odd],  N_O(X) = N_T(X)/2 + [floor(sqrt X) odd].
Prime count psi(X) = sum_{n<=X} c(n), where -F'/F = sum c(n) n^-s (a(1) normalised to 1).
  For the torus c = von Mangoldt of Z[i] (log Norm p on prime-ideal powers).
  Explicit formula: psi(X) = X - sum_rho X^rho/rho + ...   so the zero with largest real part sets the drift.
"""
import json, math, numpy as np, sys

X = int(float(sys.argv[1])) if len(sys.argv) > 1 else 2_000_000

def r2_table(X):
    R = math.isqrt(X)
    m = np.arange(-R, R + 1)
    q = (m[:, None]**2 + m[None, :]**2).ravel()
    q = q[q <= X]
    return np.bincount(q, minlength=X + 1).astype(float)

def coeffs(X):
    r2 = r2_table(X); r2[0] = 0
    et = np.zeros(X + 1); j = np.arange(1, math.isqrt(X) + 1); et[j*j] = (-1.0)**(j + 1)
    return {'torus': r2, 'klein': r2/2 - et, 'odd': r2/2 + et}

def log_derivative_coeffs(a):
    """c(n) with -F'/F = sum c(n) n^-s, F = sum a(n) n^-s (a(1) != 0)."""
    a = a/a[1]; X = len(a) - 1
    logn = np.log(np.maximum(np.arange(X + 1), 1))
    c = np.zeros(X + 1); acc = np.zeros(X + 1)
    for d in range(1, X + 1):
        c[d] = a[d]*logn[d] - acc[d]
        if 2*d <= X:
            k = X//d
            acc[2*d::d] += c[d]*a[2:k + 1]
    return c

if __name__ == '__main__':
    A = coeffs(X)
    grid = np.unique(np.round(np.logspace(1, np.log10(X), 400)).astype(int))
    out = {'X': X, 'grid': grid.tolist()}
    for k, a in A.items():
        N = np.cumsum(a)
        c = log_derivative_coeffs(a)
        psi = np.cumsum(c)
        out[k] = {'N': N[grid].tolist(), 'psi': psi[grid].tolist(),
                  'psi1': np.cumsum(psi - np.arange(X + 1))[grid].tolist(),
                  'c_head': c[:40].tolist()}
        print(k, 'N(X)', N[X], 'psi(X)-X', psi[X] - X, flush=True)
    NT, NK, NO = (np.cumsum(A[k]) for k in ('torus', 'klein', 'odd'))
    odd_root = (np.floor(np.sqrt(np.arange(X + 1) + 0.5)).astype(int) % 2 == 1).astype(float)
    out['lattice_identity_max_err'] = float(max(np.abs(NK - (NT/2 - odd_root)).max(), np.abs(NO - (NT/2 + odd_root)).max()))
    print('lattice identity max error', out['lattice_identity_max_err'])
    json.dump(out, open('drift.json', 'w'))
