"""Blaschke (equi-affine) structure of a parametrised hypersurface in R^{n+1}.

G_ij = det(x_1..x_n, x_ij),  h = G / |det G|^{1/(n+2)}   (Blaschke metric)
xi   = (1/n) Lap_h x                                     (affine normal)
x_ij = Gamma^k_ij x_k + h_ij xi                          (induced connection)
K    = Gamma - LeviCivita(h)  (difference tensor),  C = -2 h(K.,.)  (cubic form)
"""
import sympy as sp


def structure(pos, coords):
    n = len(coords)
    X = sp.Matrix(pos)
    Xi = [X.diff(c) for c in coords]
    G = sp.zeros(n, n)
    for i in range(n):
        for j in range(n):
            G[i, j] = sp.Matrix.hstack(*Xi, X.diff(coords[i]).diff(coords[j])).det()
    detG = G.det()
    h = G / sp.Abs(detG) ** sp.Rational(1, n + 2)
    return X, Xi, G, h


def affine_normal_at(pos, coords, orient=1):
    """Symbolic X, X_i, h, xi. orient=+1/-1 picks the side that makes h positive definite."""
    n = len(coords)
    X, Xi, G, h = structure(pos, coords)
    h = orient * h
    # sign: make h positive definite where convex (G definite) by orienting with sign of G[0,0]
    hinv = h.inv()
    sq = sp.sqrt(sp.Abs(h.det()))
    lap = sp.zeros(len(pos), 1)
    for a in range(n):
        for b in range(n):
            lap += (sq * hinv[a, b] * Xi[b]).diff(coords[a])
    xi = lap / (sq * n)
    return X, Xi, h, xi


class Surface:
    def __init__(self, pos, coords, ref_point):
        self.pos, self.coords = pos, coords
        _, _, G, _ = structure(pos, coords)
        orient = 1 if sp.N(G[0, 0].subs(ref_point)) > 0 else -1
        self.X, self.Xi, self.h, self.xi = affine_normal_at(pos, coords, orient)

    def at(self, point, prec=30):
        return evaluate_sym(self.X, self.Xi, self.h, self.xi, self.coords, point, prec)


def evaluate_sym(X, Xi, h, xi, coords, point, prec=30):
    n = len(coords)
    Xv = X.subs(point).evalf(prec)
    xiv = xi.subs(point).evalf(prec)
    hv = h.subs(point).evalf(prec)
    # induced connection: x_ij = Gamma^k_ij x_k + h_ij xi  -> solve 3x3 per (i,j)
    M = sp.Matrix.hstack(*[x.subs(point).evalf(prec) for x in Xi], xiv)
    Gam = [[[0] * n for _ in range(n)] for _ in range(n)]  # Gam[k][i][j]
    for i in range(n):
        for j in range(n):
            rhs = X.diff(coords[i]).diff(coords[j]).subs(point).evalf(prec)
            sol = M.LUsolve(rhs)
            for k in range(n):
                Gam[k][i][j] = sol[k]
    # Levi-Civita of h at the point needs dh: evaluate derivatives of h
    dh = [[[sp.N(h[i, j].diff(coords[l]).subs(point), prec) for l in range(n)]
           for j in range(n)] for i in range(n)]
    hinv = hv.inv()
    LC = [[[sum(hinv[k, m] * (dh[m][i][j] + dh[m][j][i] - dh[i][j][m]) for m in range(n)) / 2
            for j in range(n)] for i in range(n)] for k in range(n)]
    K = [[[Gam[k][i][j] - LC[k][i][j] for j in range(n)] for i in range(n)] for k in range(n)]
    K2 = sum(hinv[i, a] * hinv[j, b] * hv[k, c] * K[k][i][j] * K[c][a][b]
             for i in range(n) for j in range(n) for k in range(n)
             for a in range(n) for b in range(n) for c in range(n))
    return Xv, xiv, hv, sp.N(K2, prec)
