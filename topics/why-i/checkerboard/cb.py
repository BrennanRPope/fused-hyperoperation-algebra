"""Shared code: Feynman's checkerboard and Kac's telegraph as one weighted-circle chain.
Units c = hbar = 1. Lattice step eps in time and space. Two directions R (moving right), L (moving left).
State vector index: 2*site + d, d = 0 (R), 1 (L). Sites 0..n-1 on a ring, origin at n//2.
One time step = turn (dial factor 1 + k*t, t = -eps*mu*X) then move (R one site right, L one site left).
"""
import numpy as np

def circ(a, b, k):
    """weighted circle operation for (possibly non-commuting) matrices: a o_k b = a + b + k a b"""
    return a + b + k*(a @ b)

def move_matrix(n):
    M = np.zeros((2*n, 2*n), dtype=complex)
    for s in range(n):
        M[2*((s + 1) % n) + 0, 2*s + 0] = 1   # R moves right
        M[2*((s - 1) % n) + 1, 2*s + 1] = 1   # L moves left
    return M

def swap_matrix(n):
    X = np.zeros((2*n, 2*n), dtype=complex)
    for s in range(n):
        X[2*s + 1, 2*s] = 1; X[2*s, 2*s + 1] = 1
    return X

def evolve(N, eps, mu, k, normalized=False, start_dir=0, n=None):
    """Run N steps from one particle at the origin moving in start_dir. Returns array (N+1, n, 2).
    normalized=False: continue weight 1, turn weight -k*eps*mu (pure chain).
    normalized=True : k=-1 -> stochastic step (1 - a eps, a eps); k=-1j -> unitary step (cos, i sin)."""
    n = n or 2*N + 3
    psi = np.zeros((n, 2), dtype=complex); psi[n//2, start_dir] = 1
    w = -k*eps*mu
    if normalized:
        if np.isclose(k, -1): c_, t_ = 1 - eps*mu, eps*mu
        else: c_, t_ = np.cos(eps*mu), 1j*np.sin(eps*mu)
    else: c_, t_ = 1, w
    out = [psi.copy()]
    for _ in range(N):
        turned = np.stack([c_*psi[:, 0] + t_*psi[:, 1], c_*psi[:, 1] + t_*psi[:, 0]], axis=1)
        psi = np.stack([np.roll(turned[:, 0], 1), np.roll(turned[:, 1], -1)], axis=1)
        out.append(psi.copy())
    return np.array(out)

def kernel(x, t, w):
    """Continuum kernel of  d_t G_R = -d_x G_R + w G_L,  d_t G_L = d_x G_L + w G_R, from (R, x=0).
    Smooth parts inside the light cone |x| < t (the R part also has the atom delta(x - t)).
    w real: Kac's telegraph (undamped); w = i m: 1+1 Dirac. Uses I0, I1 of complex argument."""
    from scipy.special import iv
    u = np.sqrt(t**2 - x**2 + 0j)
    GR = (w/2)*np.sqrt((t + x)/(t - x) + 0j)*iv(1, w*u)
    GL = (w/2)*iv(0, w*u)
    return GR, GL
