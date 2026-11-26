"""Face-power dial: one roundness dial p aimed at any polytope (DRAFT).

Recipe (for a polytope P given by its corners):
  1. Facets u_i . x <= h_i (from the convex hull, in P's own dimension d).
  2. Centre z = analytic centre: the point where sum_i u_i / (h_i - u_i.z) = 0 (it maximises sum log slack).
  3. Scaled normals n_i = u_i / (h_i - u_i.z).  Then facets sit at n_i . y <= 1 (y = x - z) and sum_i n_i = 0.
  4. Whitening W = H^(1/2), H = sum n_i n_i^T (the log-barrier Hessian at z).  In whitened coordinates w = W y the
     normals are N_i = W^-1 n_i and sum N_i N_i^T = I.
  5. Offset a = (R - 1)/2 with R = max over corners v and facets i of -n_i.(v - z)  (a = 0 when P is symmetric).
  6. Shape at dial p:  K_p = { y : M_p( a + n_i . y ) <= a + 1 },  M_p = power mean (1/m sum |.|^p)^(1/p).
     p = 2: exactly the ellipsoid y^T H y <= m(2a+1) (round after whitening).  p -> inf: exactly P.
Only base inputs: p and the shape.  a, z, H are read off the shape (derived), never dials.
"""
import numpy as np
from scipy.spatial import ConvexHull


def affine_basis(V):
    """Orthonormal coordinates on the affine hull of V. Returns (coords, origin, basis)."""
    V = np.asarray(V, float)
    c = V.mean(0)
    U, S, Wt = np.linalg.svd(V - c, full_matrices=False)
    d = int((S > 1e-9 * S[0]).sum())
    B = Wt[:d]
    return (V - c) @ B.T, c, B


def facets_from_vertices(V):
    """Unique facet inequalities u.x <= h (|u| = 1) of conv(V), V full-dimensional."""
    hull = ConvexHull(V)
    eq = hull.equations  # u.x + b <= 0
    U, h = eq[:, :-1], -eq[:, -1]
    nrm = np.linalg.norm(U, axis=1)
    U, h = U / nrm[:, None], h / nrm
    keys = {}
    for u, hh in zip(U, h):
        k = tuple(np.round(np.append(u, hh), 7))
        keys[k] = (u, hh)
    U = np.array([v[0] for v in keys.values()])
    h = np.array([v[1] for v in keys.values()])
    return U, h


def analytic_centre(U, h, z0=None, iters=100):
    z = np.zeros(U.shape[1]) if z0 is None else z0.copy()
    for _ in range(iters):
        s = h - U @ z
        g = (U / s[:, None]).sum(0)          # gradient of -sum log s
        Hs = (U.T / s ** 2) @ U
        dz = np.linalg.solve(Hs, -g)
        t = 1.0
        while np.any(h - U @ (z + t * dz) <= 0):
            t *= 0.5
        z = z + t * dz
        if np.linalg.norm(g) < 1e-14:
            break
    return z


class FaceDial:
    """A polytope as a first-class shape with its face-power dial."""

    def __init__(self, V, name="", whiten=True, a=None):
        V = np.asarray(V, float)
        self.name = name
        X, self.c0, self.B = affine_basis(V)
        self.d = X.shape[1]
        U, h = facets_from_vertices(X)
        z = analytic_centre(U, h, z0=X.mean(0))
        n = U / (h - U @ z)[:, None]
        self.m = len(n)
        self.z = z
        self.n = n
        H = n.T @ n
        self.H = H
        evals, evecs = np.linalg.eigh(H)
        self.tightness = evals.max() / evals.min() - 1  # 0 when H is a multiple of I (round at p=2 without whitening)
        if whiten:
            self.W = evecs @ np.diag(np.sqrt(evals)) @ evecs.T
        else:
            self.W = np.eye(self.d) * np.sqrt(evals.mean())
        self.Winv = np.linalg.inv(self.W)
        self.N = n @ self.Winv            # normals in whitened coordinates: facets at N_i.w <= 1
        Y = (X - z) @ self.W.T            # corners in whitened coordinates
        self.Vw = Y
        self.R = float((-(Y @ self.N.T)).max())
        self.symmetric = abs(self.R - 1) < 1e-7
        self.a = (0.0 if self.symmetric else max(0.0, (self.R - 1) / 2)) if a is None else a
        self.sumN = np.linalg.norm(self.N.sum(0))

    # ---- the dial
    def gauge_value(self, w, p):
        """M_p(a + N.w) for points w (k x d)."""
        t = self.a + np.atleast_2d(w) @ self.N.T
        if np.isinf(p):
            return np.abs(t).max(1)
        A = np.abs(t)
        mx = A.max(1, keepdims=True)
        return mx[:, 0] * (np.mean((A / mx) ** p, axis=1)) ** (1 / p)

    def radial(self, U, p, tol=1e-13):
        """Distance from the centre to the boundary of K_p along unit directions U (k x d), whitened coords."""
        U = np.atleast_2d(U)
        T = U @ self.N.T
        a = self.a
        target = a + 1
        lo = np.zeros(len(U))
        hi = np.ones(len(U))
        f = lambda r: self._mp(a + r[:, None] * T, p) - target
        while np.any(f(hi) < 0):
            hi = np.where(f(hi) < 0, hi * 2, hi)
        for _ in range(200):
            mid = 0.5 * (lo + hi)
            neg = f(mid) < 0
            lo, hi = np.where(neg, mid, lo), np.where(neg, hi, mid)
            if np.max(hi - lo) < tol * np.max(hi):
                break
        return 0.5 * (lo + hi)

    def to_original(self, U):
        """Whitened directions -> unit directions in the shape's original coordinates."""
        D = np.atleast_2d(U) @ self.Winv.T @ self.B
        return D / np.linalg.norm(D, axis=1, keepdims=True)

    def radial_polytope(self, U):
        """Exact radial function of P itself (p = inf)."""
        T = np.atleast_2d(U) @ self.N.T
        return 1.0 / T.max(1)

    @staticmethod
    def _mp(t, p):
        A = np.abs(t)
        if np.isinf(p):
            return A.max(1)
        mx = A.max(1, keepdims=True)
        return mx[:, 0] * (np.mean((A / mx) ** p, axis=1)) ** (1 / p)


def sphere_dirs(d, k, seed=0):
    rng = np.random.default_rng(seed)
    G = rng.standard_normal((k, d))
    return G / np.linalg.norm(G, axis=1, keepdims=True)


def roundness_gap(fd, p, U):
    r = fd.radial(U, p)
    return r.max() / r.min() - 1


# ----------------------------------------------------------------- shapes
def regular_polygon(n):
    t = 2 * np.pi * np.arange(n) / n + np.pi / 2
    return np.column_stack([np.cos(t), np.sin(t)])


def simplex(d):
    return np.eye(d + 1)  # lives in the sum = 1 hyperplane; FaceDial takes the affine hull


def cube(d):
    import itertools
    return np.array(list(itertools.product([-1, 1], repeat=d)), float)


def cross_polytope(d):
    return np.vstack([np.eye(d), -np.eye(d)])


PHI = (1 + 5 ** 0.5) / 2


def icosahedron():
    V = []
    for s1 in (-1, 1):
        for s2 in (-1, 1):
            V += [(0, s1, s2 * PHI), (s1, s2 * PHI, 0), (s2 * PHI, 0, s1)]
    return np.array(V, float)


def dodecahedron():
    import itertools
    V = [c for c in itertools.product([-1, 1], repeat=3)]
    for s1 in (-1, 1):
        for s2 in (-1, 1):
            V += [(0, s1 / PHI, s2 * PHI), (s1 / PHI, s2 * PHI, 0), (s2 * PHI, 0, s1 / PHI)]
    return np.array(V, float)


def permutohedron(n):
    import itertools
    return np.array(list(itertools.permutations(range(1, n + 1))), float)


def loday_associahedron(n):
    """Loday's realisation of the associahedron K_{n+1} (dimension n-1) from planar binary trees with n leaves... built
    from all binary bracketings of n+1 letters: coordinate i = (#leaves left of node i) * (#leaves right of node i)."""
    def trees(lo, hi):  # internal nodes between leaves lo..hi
        if lo == hi:
            return [None]
        out = []
        for k in range(lo, hi):
            for L in trees(lo, k):
                for R in trees(k + 1, hi):
                    out.append((k, L, R, lo, hi))
        return out

    def coords(t, x):
        if t is None:
            return
        k, L, R, lo, hi = t
        x[k] = (k - lo + 1) * (hi - k)
        coords(L, x)
        coords(R, x)

    V = []
    for t in trees(0, n):
        x = np.zeros(n)
        coords(t, x)
        V.append(x)
    return np.array(V)


def zonotope(G):
    import itertools
    G = np.asarray(G, float)
    return np.array([np.array(s) @ G for s in itertools.product([-1, 1], repeat=len(G))])


def cyclic_polytope(npts, d):
    t = np.linspace(-1, 1, npts)
    return np.column_stack([t ** k for k in range(1, d + 1)])


def birkhoff(n):
    import itertools
    V = []
    for perm in itertools.permutations(range(n)):
        M = np.zeros((n, n))
        M[range(n), perm] = 1
        V.append(M.ravel())
    return np.array(V)


def hypersimplex(k, n):
    import itertools
    V = []
    for S in itertools.combinations(range(n), k):
        x = np.zeros(n)
        x[list(S)] = 1
        V.append(x)
    return np.array(V)


# ----------------------------------------------------------------- the temperature (log-sum-exp) route, same normals
def lse_radial(fd, U, beta, tol=1e-13):
    """Boundary of { w : (1/beta) log mean exp(beta N_i.w) <= 1 } along directions U (the free-energy route)."""
    from scipy.special import logsumexp
    T = np.atleast_2d(U) @ fd.N.T
    m = fd.m
    f = lambda r: (logsumexp(beta * r[:, None] * T, axis=1) - np.log(m)) / beta - 1
    lo, hi = np.zeros(len(T)), np.ones(len(T))
    while np.any(f(hi) < 0):
        hi = np.where(f(hi) < 0, hi * 2, hi)
    for _ in range(200):
        mid = 0.5 * (lo + hi)
        neg = f(mid) < 0
        lo, hi = np.where(neg, mid, lo), np.where(neg, hi, mid)
        if np.max(hi - lo) < tol * np.max(hi):
            break
    return 0.5 * (lo + hi)


def circle_dirs(n, offset=0.0):
    t = offset + 2 * np.pi * np.arange(n) / n
    return t, np.column_stack([np.cos(t), np.sin(t)])
