"""Small general-relativity helper (sympy): Christoffels, Riemann, Ricci, Einstein, divergence."""
import sympy as sp


def christoffel(g, X):
    n = len(X)
    gi = sp.simplify(g.inv())
    G = [[[0] * n for _ in range(n)] for _ in range(n)]
    for a in range(n):
        for b in range(n):
            for c in range(n):
                G[a][b][c] = sp.simplify(sum(gi[a, d] * (sp.diff(g[d, b], X[c]) + sp.diff(g[d, c], X[b])
                                                          - sp.diff(g[b, c], X[d])) for d in range(n)) / 2)
    return G, gi


def riemann(G, X):
    """R^a_{bcd} = d_c G^a_{db} - d_d G^a_{cb} + G^a_{ce} G^e_{db} - G^a_{de} G^e_{cb}."""
    n = len(X)
    R = [[[[0] * n for _ in range(n)] for _ in range(n)] for _ in range(n)]
    for a in range(n):
        for b in range(n):
            for c in range(n):
                for d in range(n):
                    e = sp.diff(G[a][d][b], X[c]) - sp.diff(G[a][c][b], X[d])
                    e += sum(G[a][c][f] * G[f][d][b] - G[a][d][f] * G[f][c][b] for f in range(n))
                    R[a][b][c][d] = sp.simplify(e)
    return R


def ricci(R, n):
    return sp.Matrix(n, n, lambda b, d: sp.simplify(sum(R[a][b][a][d] for a in range(n))))


def einstein_mixed(g, X):
    """Return (G^a_b mixed Einstein tensor, Ricci_ab, scalar R, Christoffels, inverse metric)."""
    n = len(X)
    G, gi = christoffel(g, X)
    R = riemann(G, X)
    Ric = ricci(R, n)
    Rs = sp.simplify(sum(gi[a, b] * Ric[a, b] for a in range(n) for b in range(n)))
    Gmix = sp.Matrix(n, n, lambda a, b: sp.simplify(sum(gi[a, c] * Ric[c, b] for c in range(n))
                                                    - sp.Rational(1, 2) * Rs * (1 if a == b else 0)))
    return Gmix, Ric, Rs, G, gi, R


def divergence_mixed(T, G, X):
    """(nabla_a T^a_b) for a mixed (1,1) tensor T^a_b."""
    n = len(X)
    out = []
    for b in range(n):
        e = sum(sp.diff(T[a, b], X[a]) for a in range(n))
        e += sum(G[a][a][c] * T[c, b] for a in range(n) for c in range(n))
        e -= sum(G[c][a][b] * T[a, c] for a in range(n) for c in range(n))
        out.append(sp.simplify(e))
    return out
