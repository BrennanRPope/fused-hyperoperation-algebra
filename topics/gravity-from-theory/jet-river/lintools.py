"""First-order (linearised) Einstein tensor about a background, without heavy simplification.

g = g0 + eps g1.  Returns G1^a_b = d/d eps of G^a_b at eps = 0.  Much faster than the full computation
for small non-diagonal perturbations (swirl, twist) of the river."""
import sympy as sp


def _chr(gi, g, X):
    n = len(X)
    return [[[sum(gi[a, d] * (sp.diff(g[d, b], X[c]) + sp.diff(g[d, c], X[b]) - sp.diff(g[b, c], X[d]))
                  for d in range(n)) / 2 for c in range(n)] for b in range(n)] for a in range(n)]


def _ricci(G0, G1, X):
    """Ricci to first order: returns (Ric0, Ric1) from background and first-order Christoffels."""
    n = len(X)

    def ric(A, B):  # bilinear part with A, B Christoffel tables
        return sp.Matrix(n, n, lambda b, d: sum(A[a][a][e] * B[e][d][b] - A[a][d][e] * B[e][a][b]
                                                 for a in range(n) for e in range(n)))
    lin = lambda Gm: sp.Matrix(n, n, lambda b, d: sum(sp.diff(Gm[a][d][b], X[a]) - sp.diff(Gm[a][a][b], X[d])
                                                      for a in range(n)))
    Ric0 = lin(G0) + ric(G0, G0)
    Ric1 = lin(G1) + ric(G0, G1) + ric(G1, G0)
    return Ric0, Ric1


def lin_einstein(g0, g1, X, simp=sp.simplify):
    n = len(X)
    gi0 = sp.simplify(g0.inv())
    gi1 = -gi0 * g1 * gi0
    G0 = [[[sp.simplify(e) for e in row] for row in tab] for tab in _chr(gi0, g0, X)]
    # first-order Christoffels: gi1 * dg0 + gi0 * dg1
    A = _chr(gi1, g0, X)
    B = _chr(gi0, g1, X)
    G1 = [[[sp.cancel(A[a][b][c] + B[a][b][c]) for c in range(n)] for b in range(n)] for a in range(n)]
    Ric0, Ric1 = _ricci(G0, G1, X)
    Ric0 = Ric0.applyfunc(sp.simplify)
    Rs0 = sp.simplify(sum(gi0[a, b] * Ric0[a, b] for a in range(n) for b in range(n)))
    Rs1 = sum(gi1[a, b] * Ric0[a, b] + gi0[a, b] * Ric1[a, b] for a in range(n) for b in range(n))
    G1mix = sp.Matrix(n, n, lambda a, b: sum(gi1[a, c] * Ric0[c, b] + gi0[a, c] * Ric1[c, b] for c in range(n))
                      - sp.Rational(1, 2) * Rs1 * (1 if a == b else 0))
    return G1mix.applyfunc(lambda e: sp.factor(simp(e)))
