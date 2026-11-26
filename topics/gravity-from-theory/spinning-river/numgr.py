"""Exact-derivative, high-precision curvature at a point.

The metric's first and second derivatives are taken symbolically (sympy); everything after that
(inverse, Christoffels, Riemann, Ricci) is done in mpmath at 50 digits at a chosen point. Fast and
accurate to ~1e-40, so "zero" means |value| < 1e-30."""
import itertools

import mpmath as mp
import sympy as sp

mp.mp.dps = 50


def prepare(g, X):
    n = len(X)
    dg = [[[sp.diff(g[a, b], X[c]) for c in range(n)] for b in range(n)] for a in range(n)]
    ddg = [[[[sp.diff(dg[a][b][c], X[d]) for d in range(n)] for c in range(n)] for b in range(n)] for a in range(n)]
    syms = list(X)
    f = sp.lambdify(syms, [g.tolist(), dg, ddg], modules="mpmath")
    return f


def curvature_at(f, point):
    G, dG, ddG = f(*[mp.mpf(p) for p in point])
    n = len(G)
    g = mp.matrix(G)
    gi = g ** -1
    # Gamma^a_bc and its derivative d_d Gamma^a_bc
    low = [[[(dG[d][b][c] + dG[d][c][b] - dG[b][c][d]) / 2 for c in range(n)] for b in range(n)] for d in range(n)]
    dlow = [[[[(ddG[d][b][c][e] + ddG[d][c][b][e] - ddG[b][c][d][e]) / 2 for e in range(n)] for c in range(n)]
             for b in range(n)] for d in range(n)]
    # d_e g^{ad} = -g^{ap} d_e g_pq g^{qd}
    dgi = [[[-sum(gi[a, p] * dG[p][q][e] * gi[q, d] for p in range(n) for q in range(n)) for e in range(n)]
            for d in range(n)] for a in range(n)]
    Gam = [[[sum(gi[a, d] * low[d][b][c] for d in range(n)) for c in range(n)] for b in range(n)] for a in range(n)]
    dGam = [[[[sum(dgi[a][d][e] * low[d][b][c] + gi[a, d] * dlow[d][b][c][e] for d in range(n)) for e in range(n)]
              for c in range(n)] for b in range(n)] for a in range(n)]
    # R^a_bcd = d_c Gam^a_db - d_d Gam^a_cb + Gam^a_ce Gam^e_db - Gam^a_de Gam^e_cb
    Riem = {}
    for a, b, c, d in itertools.product(range(n), repeat=4):
        Riem[a, b, c, d] = (dGam[a][d][b][c] - dGam[a][c][b][d]
                            + sum(Gam[a][c][e] * Gam[e][d][b] - Gam[a][d][e] * Gam[e][c][b] for e in range(n)))
    Ric = [[sum(Riem[a, b, a, d] for a in range(n)) for d in range(n)] for b in range(n)]
    return Riem, Ric, gi


def max_abs(xs):
    return max(abs(x) for x in xs)


def riemann_max(f, point):
    Riem, _, _ = curvature_at(f, point)
    return max_abs(Riem.values())


def ricci_max(f, point):
    _, Ric, _ = curvature_at(f, point)
    return max_abs([x for row in Ric for x in row])
