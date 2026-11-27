"""Reflection positivity for the weighted-circle chain (topics/why-i/reflection-positivity).

A chain of weighted-circle steps a_N o_k ... o_k a_1 has multiplier
    1 + k (a_N o_k ... o_k a_1) = (1 + k a_N) ... (1 + k a_1)        (wAdd_lift)
With steps a_j = eps*H the slice multiplier is T = 1 + k eps H, and a path
x_0 -> x_1 -> ... -> x_N gets weight  u(x_0) * prod_j T[x_{j+1}, x_j] * u(x_N).

The time reflection theta sends x_j -> x_{N-j} and complex-conjugates a functional.
Reflection positivity (RP) asks  <theta F, F> >= 0  for every F that depends only on
the later half of the path.  Site reflection: N even, mirror through x_{N/2}.
Link reflection: N odd, mirror through the middle step.

Everything here is brute force over paths, so the Gram matrices are honest path sums,
not formulas that already assume the answer.
"""
import itertools
import numpy as np


def slice_multiplier(k, eps, H, kind="explicit"):
    """Multiplier of one slice. explicit: 1 + k eps H.  rewind: (1 - k eps H)^-1, the
    weighted-circle inverse of a backward step (Chernoff's implicit step)."""
    n = H.shape[0]
    I = np.eye(n)
    if kind == "explicit":
        return I + k * eps * H
    if kind == "rewind":
        return np.linalg.inv(I - k * eps * H)
    raise ValueError(kind)


def path_weight(T, u, path):
    w = u[path[0]] * u[path[-1]]
    for a, b in zip(path[:-1], path[1:]):
        w = w * T[b, a]
    return w


def gram(T, N, u=None):
    """Gram matrix G[i, j] = <theta F_i, F_j> by summing over every path of length N.

    The F_j are the indicator functionals of every configuration of the later half
    (sites M..N for site reflection, N even; sites M+1..N for link reflection, N odd,
    M = N // 2).  They span every functional of the later half, so RP <=> G is
    Hermitian positive semidefinite."""
    n = T.shape[0]
    if u is None:
        u = np.ones(n)
    M = N // 2
    first = M if N % 2 == 0 else M + 1          # first site of the later half
    half = N - first + 1                          # number of sites in the later half
    configs = list(itertools.product(range(n), repeat=half))
    index = {c: i for i, c in enumerate(configs)}
    G = np.zeros((len(configs), len(configs)), dtype=complex)
    for path in itertools.product(range(n), repeat=N + 1):
        later = path[first:]                      # what F_j looks at
        mirrored = tuple(path[N - s] for s in range(first, N + 1))  # theta F_i looks here
        # theta F_i (x) = conj(F_i(theta x)); with indicators the conj is trivial,
        # and the coefficient vector c enters as conj(c_i) c_j below.
        G[index[mirrored], index[later]] += path_weight(T, u, path)
    return G


def rp_violation(G):
    """0 when G is Hermitian PSD.  Otherwise the worst normalized failure of
    <theta F, F> >= 0: either the pairing is not real (non-Hermitian part) or it
    goes negative (negative eigenvalue of the Hermitian part)."""
    scale = max(np.abs(G).max(), 1e-300)
    herm = (G + G.conj().T) / 2
    anti = (G - G.conj().T) / 2
    neg = max(0.0, -np.linalg.eigvalsh(herm).min())
    return max(neg, np.abs(np.linalg.eigvals(anti)).max()) / scale


def complexify(M):
    """Real 2n x 2n form of a complex n x n matrix, i -> J = [[0,-1],[1,0]] per level."""
    return np.block([[M.real, -M.imag], [M.imag, M.real]])


def J_of(n):
    return complexify(1j * np.eye(n))


def expm(M, terms=30):
    """Matrix exponential by scaling and squaring with a Taylor series (no scipy needed)."""
    M = np.asarray(M)
    norm = np.abs(M).sum(axis=0).max()
    s = max(0, int(np.ceil(np.log2(norm))) + 1) if norm > 0 else 0
    A = M / (2 ** s)
    E = np.eye(M.shape[0], dtype=A.dtype)
    term = np.eye(M.shape[0], dtype=A.dtype)
    for j in range(1, terms):
        term = term @ A / j
        E = E + term
    for _ in range(s):
        E = E @ E
    return E
