"""DRAFT. Conley-Zehnder (Maslov-type) index of each Floquet pair across iterates.

Convention: state x = (q, p), Omega(x, y) = x^T J y with J = [[0, I], [-I, 0]], so x' = J S x for H = x^T S x / 2
and a positive-definite H gives a positive index (the harmonic oscillator over less than one turn has index +1).

For a nondegenerate monodromy M whose unit-circle eigenvalues are simple, Long's iteration formula
(Long 2000, Adv. Math.; book 2002, Thm 8.3.1) reads
    mu(gamma^m) = m (mu(gamma) - r) + 2 sum_j E(m theta_j / 2 pi) - r,      E(a) = least integer >= a,
where r is the number of elliptic pairs and theta_j in (0, 2 pi) is the pair's angle normalised by its Krein sign
(theta = phi for a positive pair at e^{i phi}, 2 pi - phi for a negative one; phi in (0, pi)).
Hyperbolic pairs and complex quartets only enter the linear term.  So, up to a term m * (even integer) that depends
on the trivialisation along the orbit, the iterate indices are
    muhat(m) = sum_elliptic (2 E(m theta_j / 2 pi) - 1) + m * (#inverse-hyperbolic pairs),
and muhat(m) - m * Dhat (Dhat = sum theta_j / pi + #IH, the mean index up to an even integer) is the bounded,
trivialisation-free part: zero for an orbit with no elliptic pair, a sum of one sawtooth per elliptic pair otherwise.
"""
import numpy as np
from scipy.linalg import expm

def Jmat(n):
    Z, I = np.zeros((n, n)), np.eye(n)
    return np.block([[Z, I], [-I, Z]])

def E(a):
    return np.ceil(a - 1e-12)

def krein_pairs(M, Vlift=None, J=None, tol=1e-7):
    """Floquet pairs of a symplectic M (2n x 2n), or of a reduced block M acting on span(Vlift) (columns of
    Vlift represent the reduced space inside a bigger symplectic space with form J).
    Returns a list of dicts: kind in E/H/IH/CU; for E also phi (0..pi), krein sign s, theta (0..2pi)."""
    lam, W = np.linalg.eig(M)
    if Vlift is None:
        Vlift = np.eye(M.shape[0]); J = Jmat(M.shape[0] // 2)
    out, used = [], np.zeros(len(lam), bool)
    for i in np.argsort(-np.abs(lam)):
        if used[i]: continue
        l = lam[i]
        if abs(abs(l) - 1) < tol and abs(l.imag) > tol:                # elliptic: take the e^{i phi}, phi in (0, pi)
            k = i if l.imag > 0 else int(np.argmin(np.abs(lam - np.conj(l)) + used * 9))
            phi = float(np.angle(lam[k]))
            w = Vlift @ W[:, k]
            kappa = float(np.real(-1j * (np.conj(w) @ J @ w)))
            s = 1 if kappa > 0 else -1
            out.append(dict(kind='E', phi=phi, s=s, theta=phi if s > 0 else 2 * np.pi - phi, kappa=kappa))
            used[i] = True; used[int(np.argmin(np.abs(lam - np.conj(l)) + used * 9))] = True
        elif abs(l.imag) > tol:                                         # complex quartet
            out.append(dict(kind='CU')); used[i] = True
            for z in (np.conj(l), 1 / l, 1 / np.conj(l)):
                used[int(np.argmin(np.abs(lam - z) + used * 9))] = True
            out.append(dict(kind='CU2'))
        else:
            out.append(dict(kind='H' if l.real > 0 else 'IH')); used[i] = True
            used[int(np.argmin(np.abs(lam - 1 / l) + used * 9))] = True
    return [p for p in out if p['kind'] != 'CU2']

def muhat(pairs, m):
    m = np.asarray(m)
    v = sum((2 * E(m * p['theta'] / (2 * np.pi)) - 1) for p in pairs if p['kind'] == 'E')
    return v + m * sum(1 for p in pairs if p['kind'] == 'IH')

def mean_hat(pairs):
    return sum(p['theta'] / np.pi for p in pairs if p['kind'] == 'E') + sum(1 for p in pairs if p['kind'] == 'IH')

def deviation(pairs, m):
    return muhat(pairs, m) - np.asarray(m) * mean_hat(pairs)

def decode_r(mu, K=None, grid=200001):
    """How many elliptic pairs do the integers mu(1..K) need?  0 if mu(m) = m c exactly; 1 if
    mu(m) = m c + 2 E(m theta/2pi) - 1 for some theta (c = mu(1) - 1); else 2.  Uses the integers only."""
    mu = np.asarray(mu); K = len(mu); m = np.arange(1, K + 1)
    if np.all(mu == m * mu[0]): return 0
    c = mu[0] - 1
    th = np.linspace(0, 2 * np.pi, grid)[1:-1]
    ok = np.ones(len(th), bool)
    for j in range(K):
        ok &= (m[j] * c + 2 * E(m[j] * th / (2 * np.pi)) - 1) == mu[j]
        if not ok.any(): return 2
    return 1

# ---------- independent check: Robbin-Salamon crossing count on a periodic linear system ----------
def S_of_t(S0, S1, S2, t):
    return S0 + S1 * np.cos(t) + S2 * np.sin(t)

def path(S0, S1, S2, tmax, nstep):
    """Phi(t) for x' = J S(t) x, Phi(0) = I, on a uniform grid (4th-order Magnus-free RK via small expm steps)."""
    n2 = S0.shape[0]; J = Jmat(n2 // 2); h = tmax / nstep
    Ph = [np.eye(n2)]; P = np.eye(n2)
    c = (0.5 - np.sqrt(3) / 6, 0.5 + np.sqrt(3) / 6)                  # 4th-order Magnus with 2 Gauss points
    for k in range(nstep):
        t = k * h
        A1 = J @ S_of_t(S0, S1, S2, t + c[0] * h); A2 = J @ S_of_t(S0, S1, S2, t + c[1] * h)
        Om = h / 2 * (A1 + A2) + np.sqrt(3) / 12 * h * h * (A2 @ A1 - A1 @ A2)
        P = expm(Om) @ P; Ph.append(P)
    return np.array(Ph)

def crossing_index(S0, S1, S2, tmax, nstep):
    """mu_CZ of the path on [0, tmax] = (1/2) sign S(0) + sum over interior crossings of sign(v^T S v), v in ker(Phi - I)."""
    Ph = path(S0, S1, S2, tmax, nstep)
    n2 = S0.shape[0]; ts = np.linspace(0, tmax, nstep + 1)
    mu = 0.5 * np.sum(np.sign(np.linalg.eigvalsh(S_of_t(S0, S1, S2, 0))))
    def sgn(P):                                   # sign det(P - I) from eigenvalues: robust when P is large
        lam = np.linalg.eigvals(P); real = np.abs(lam.imag) < 1e-9 * np.maximum(1, np.abs(lam))
        return int(np.prod(np.sign(lam.real[real] - 1)))
    d = np.array([sgn(P) for P in Ph], float)
    for k in range(1, nstep):
        if np.sign(d[k]) != np.sign(d[k + 1]) and np.sign(d[k]) != 0:
            j = k
            lam, W = np.linalg.eig(Ph[j]); i1 = int(np.argmin(np.abs(lam - 1)))
            v = np.real(W[:, i1]); v /= np.linalg.norm(v)
            mu += np.sign(v @ S_of_t(S0, S1, S2, ts[j]) @ v)
    return mu, Ph[-1]

def pairs_from_traces(M, taus, Vlift, J):
    """Floquet pairs typed by their traces tau = lambda + 1/lambda (robust at catalog precision); each elliptic
    pair gets its Krein sign from the eigenvector of M nearest e^{i phi}, phi = arccos(tau/2)."""
    lam, W = np.linalg.eig(M); out = []
    for t in taus:
        t = complex(t)
        if abs(t.imag) > 1e-9: out.append(dict(kind='CU')); continue
        t = t.real
        if abs(t) >= 2: out.append(dict(kind='H' if t > 2 else 'IH')); continue
        phi = float(np.arccos(t / 2)); k = int(np.argmin(np.abs(lam - np.exp(1j * phi))))
        w = Vlift @ W[:, k]
        kappa = float(np.real(-1j * (np.conj(w) @ J @ w)))
        s = 1 if kappa > 0 else -1
        out.append(dict(kind='E', phi=phi, s=s, theta=phi if s > 0 else 2 * np.pi - phi, kappa=kappa,
                        eig_miss=float(np.abs(lam[k] - np.exp(1j * phi)))))
    if [p['kind'] for p in out] == ['CU', 'CU']: out = out[:1]
    return out
