"""Shape dial vs plane dial: canonical (Vinberg) curvature of convex bodies between the disk and the triangle.

For a convex body Omega in the plane (a slice of a convex cone in R^3), Vinberg's characteristic function is
chi(y) = c * area of the polar body of (Omega - y)  =  (c/2) * integral_0^{2pi} (h(phi) - u(phi).y)^(-2) dphi,
with h the support function and u(phi) = (cos phi, sin phi).  The canonical projective metric on Omega is
    g = Hess(log chi) - (1/3) grad(log chi) grad(log chi)^T
(the Hessian metric of log chi on the cone, with the radial direction removed).  It is invariant under every
projective symmetry of Omega.  For the disk it is 3 x the Klein (hyperbolic) metric; for the triangle it is flat.
We report the normalised curvature Khat = 3 K, so the disk gives -1 and the triangle 0.  Leftover turn around a small
loop = -Khat x (area in units where the disk has curvature -1).
"""
import numpy as np
from scipy import special

TWO_PI = 2 * np.pi


class Body:
    """Convex body given by boundary points (counter-clockwise or not); support function sampled on a phi grid."""

    def __init__(self, boundary_xy, nphi=4096, name=""):
        self.name = name
        self.boundary = np.asarray(boundary_xy, float)
        self.phi = np.linspace(0, TWO_PI, nphi, endpoint=False)
        self.u = np.column_stack([np.cos(self.phi), np.sin(self.phi)])
        # support function h(phi) = max over boundary points of u . x
        h = np.empty(nphi)
        B = self.boundary
        for s in range(0, nphi, 512):
            h[s:s + 512] = (self.u[s:s + 512] @ B.T).max(axis=1)
        self.h = h
        self.w = TWO_PI / nphi  # trapezoid weight (periodic)

    def inside(self, y, margin=0.0):
        y = np.atleast_2d(y)
        return np.all(self.h[None, :] - y @ self.u.T > margin, axis=1)

    def metric(self, Y):
        """Canonical metric g at points Y (m x 2). Returns (m x 2 x 2)."""
        Y = np.atleast_2d(Y)
        f = self.h[None, :] - Y @ self.u.T  # m x nphi
        f2, f3, f4 = f ** -2, f ** -3, f ** -4
        I0 = f2.sum(axis=1) * self.w
        I1 = 2 * (f3 @ self.u) * self.w  # m x 2
        uu = np.einsum("ki,kj->kij", self.u, self.u)
        I2 = 6 * np.einsum("mk,kij->mij", f4, uu) * self.w
        gphi = I1 / I0[:, None]
        Hphi = I2 / I0[:, None, None] - np.einsum("mi,mj->mij", gphi, gphi)
        return Hphi - np.einsum("mi,mj->mij", gphi, gphi) / 3.0

    def log_chi(self, Y):
        Y = np.atleast_2d(Y)
        f = self.h[None, :] - Y @ self.u.T
        return np.log((f ** -2).sum(axis=1) * self.w)


def gauss_curvature(metric_fn, Y, step):
    """Gaussian curvature of a 2D metric at points Y via the Brioschi formula, central differences."""
    Y = np.atleast_2d(Y)
    m = len(Y)
    offs = np.array([[i, j] for i in (-1, 0, 1) for j in (-1, 0, 1)], float) * step
    P = (Y[:, None, :] + offs[None, :, :]).reshape(-1, 2)
    G = metric_fn(P).reshape(m, 3, 3, 2, 2)  # index [i+1, j+1]
    E, F, Gg = G[..., 0, 0], G[..., 0, 1], G[..., 1, 1]
    c = (slice(None), 1, 1)

    def du(A):
        return (A[:, 2, 1] - A[:, 0, 1]) / (2 * step)

    def dv(A):
        return (A[:, 1, 2] - A[:, 1, 0]) / (2 * step)

    def duu(A):
        return (A[:, 2, 1] - 2 * A[:, 1, 1] + A[:, 0, 1]) / step ** 2

    def dvv(A):
        return (A[:, 1, 2] - 2 * A[:, 1, 1] + A[:, 1, 0]) / step ** 2

    def duv(A):
        return (A[:, 2, 2] - A[:, 2, 0] - A[:, 0, 2] + A[:, 0, 0]) / (4 * step ** 2)

    E0, F0, G0 = E[c], F[c], Gg[c]
    Eu, Ev, Fu, Fv, Gu, Gv = du(E), dv(E), du(F), dv(F), du(Gg), dv(Gg)
    Evv, Fuv, Guu = dvv(E), duv(F), duu(Gg)
    M1 = np.stack([np.stack([-Evv / 2 + Fuv - Guu / 2, Eu / 2, Fu - Ev / 2], -1),
                   np.stack([Fv - Gu / 2, E0, F0], -1),
                   np.stack([Gv / 2, F0, G0], -1)], -2)
    Z = np.zeros_like(E0)
    M2 = np.stack([np.stack([Z, Ev / 2, Gu / 2], -1),
                   np.stack([Ev / 2, E0, F0], -1),
                   np.stack([Gu / 2, F0, G0], -1)], -2)
    return (np.linalg.det(M1) - np.linalg.det(M2)) / (E0 * G0 - F0 ** 2) ** 2


def khat(body, Y, step=None):
    if step is None:
        n = len(body.h)
        width = (body.h + np.roll(body.h, n // 2)).min()  # smallest width of the body
        step = 1e-3 * width
    return 3.0 * gauss_curvature(body.metric, Y, step)


class PolyBody(Body):
    """Exact polygon: Gauss-Legendre nodes on each arc between edge-normal angles (support function is smooth there)."""

    def __init__(self, V, nodes=200, name="polygon"):
        V = np.asarray(V, float)
        a, b = V[1] - V[0], V[2] - V[1]
        if a[0] * b[1] - a[1] * b[0] < 0:
            V = V[::-1]
        E = np.roll(V, -1, 0) - V
        nrm = np.column_stack([E[:, 1], -E[:, 0]])
        ang = np.sort(np.mod(np.arctan2(nrm[:, 1], nrm[:, 0]), TWO_PI))
        x, w = np.polynomial.legendre.leggauss(nodes)
        phis, ws = [], []
        for a, b in zip(ang, np.append(ang[1:], ang[0] + TWO_PI)):
            phis.append((a + b) / 2 + (b - a) / 2 * x)
            ws.append((b - a) / 2 * w)
        self.phi = np.mod(np.concatenate(phis), TWO_PI)
        self.wv = np.concatenate(ws)
        self.u = np.column_stack([np.cos(self.phi), np.sin(self.phi)])
        self.h = (self.u @ V.T).max(1)
        self.boundary, self.name = V, name

    def metric(self, Y):
        Y = np.atleast_2d(Y)
        f = self.h[None, :] - Y @ self.u.T
        W = self.wv
        I0 = (f ** -2 * W).sum(1)
        I1 = 2 * ((f ** -3 * W) @ self.u)
        uu = np.einsum("ki,kj->kij", self.u, self.u)
        I2 = 6 * np.einsum("mk,kij->mij", f ** -4 * W, uu)
        gp = I1 / I0[:, None]
        Hp = I2 / I0[:, None, None] - np.einsum("mi,mj->mij", gp, gp)
        return Hp - np.einsum("mi,mj->mij", gp, gp) / 3


def width_of(V):
    return min(np.ptp(V @ np.array([np.cos(a), np.sin(a)])) for a in np.linspace(0, np.pi, 180))


def poly_khat(V, Y):
    return 3 * gauss_curvature(PolyBody(V).metric, Y, 1e-3 * width_of(np.asarray(V)))


# ------------------------------------------------------------------ bodies
TRI_ANG = np.radians([90, 210, 330])
TRI_V = np.column_stack([np.cos(TRI_ANG), np.sin(TRI_ANG)])  # circumradius 1, inradius 1/2


def polygon_boundary(V, per_edge=2000):
    pts = []
    for a, b in zip(V, np.roll(V, -1, axis=0)):
        t = np.linspace(0, 1, per_edge, endpoint=False)[:, None]
        pts.append(a + t * (b - a))
    return np.vstack(pts)


def disk(R=1.0, n=6000):
    t = np.linspace(0, TWO_PI, n, endpoint=False)
    return Body(R * np.column_stack([np.cos(t), np.sin(t)]), name="disk")


def triangle():
    return Body(polygon_boundary(TRI_V), name="triangle")


def minkowski(t, n=6000):
    """(1-t) * disk(radius 1/2) + t * triangle(inradius 1/2): boundary via support point sums."""
    phi = np.linspace(0, TWO_PI, n, endpoint=False)
    u = np.column_stack([np.cos(phi), np.sin(phi)])
    disk_pt = 0.5 * u
    tri_pt = TRI_V[np.argmax(u @ TRI_V.T, axis=1)]
    pts = (1 - t) * disk_pt + t * tri_pt
    # add the triangle's edges, scaled, joining consecutive support points (flat sides of the sum)
    extra = []
    for i in range(n):
        a, b = pts[i], pts[(i + 1) % n]
        if np.linalg.norm(b - a) > 1e-3:
            s = np.linspace(0, 1, 400, endpoint=False)[:, None]
            extra.append(a + s * (b - a))
    if extra:
        pts = np.vstack([pts] + extra)
    return Body(pts, name=f"minkowski t={t}")


def cut_triangle(s, n=8000):
    """Triangle (circumradius 1) intersected with the disk of radius s, 1/2 <= s <= 1."""
    t = np.linspace(0, TWO_PI, n, endpoint=False)
    circ = s * np.column_stack([np.cos(t), np.sin(t)])
    tri = polygon_boundary(TRI_V, 3000)
    # keep circle points inside the triangle and triangle points inside the circle
    normals = -TRI_V  # inward-facing? triangle edges have outward normals at -vertex directions, offset 1/2
    def in_tri(P):
        return np.all(P @ (-TRI_V).T <= 0.5 + 1e-12, axis=1)
    pts = np.vstack([circ[in_tri(circ)], tri[np.linalg.norm(tri, axis=1) <= s + 1e-12]])
    return Body(pts, name=f"cut triangle s={s}")


E1 = np.array([1, -1, 0]) / np.sqrt(2)
E2 = np.array([1, 1, -2]) / np.sqrt(6)


def free_energy_body(beta, n=6000):
    """{x in plane sum=0 : (1/beta) log mean exp(-beta x_i) <= 1} in 2D coordinates (vectorised bisection)."""
    ang = np.linspace(0, TWO_PI, n, endpoint=False)
    OM = np.outer(np.cos(ang), E1) + np.outer(np.sin(ang), E2)  # n x 3
    f = lambda r: (special.logsumexp(-beta * r[:, None] * OM, axis=1) - np.log(3)) / beta - 1
    lo, hi = np.zeros(n), np.ones(n)
    while np.any(f(hi) < 0):
        hi = np.where(f(hi) < 0, hi * 2, hi)
    for _ in range(70):
        mid = 0.5 * (lo + hi)
        neg = f(mid) < 0
        lo, hi = np.where(neg, mid, lo), np.where(neg, hi, mid)
    r = 0.5 * (lo + hi)
    return Body(r[:, None] * np.column_stack([np.cos(ang), np.sin(ang)]), name=f"free energy beta={beta}")
