"""Four-charge black holes of N = 8 supergravity's STU truncation, checked as exact solutions.
e^{-1}L = R - 1/2 (dphi)^2 - 1/4 sum_i e^{a_i.phi} F_i^2, phi in R^3, a_i = tetrahedron (a_i.a_j = 4 delta_ij - 1).
ds^2 = -H^{-1/2} dt^2 + H^{1/2}(dr^2 + r^2 dOmega^2), H = H1 H2 H3 H4, H_i = 1 + k_i/r,
phi = 1/2 sum a_i ln H_i, A^i_t = alpha (1 - 1/H_i)."""
import sympy as sp
from fields import r, th, X

def einstein_fast(g):
    """Einstein tensor of a diagonal static metric, no simplification (evaluated numerically later)."""
    gi = sp.diag(*[1/g[i, i] for i in range(4)]); n = 4
    Gam = [[[sum(gi[a, d]*(sp.diff(g[d, b], X[c]) + sp.diff(g[d, c], X[b]) - sp.diff(g[b, c], X[d])) for d in range(n))/2 for c in range(n)] for b in range(n)] for a in range(n)]
    Ric = sp.Matrix(4, 4, lambda b, c: sum(sp.diff(Gam[a][b][c], X[a]) - sp.diff(Gam[a][b][a], X[c]) + sum(Gam[a][a][d]*Gam[d][b][c] - Gam[a][c][d]*Gam[d][b][a] for d in range(n)) for a in range(n)))
    Rs = sum(gi[i, i]*Ric[i, i] for i in range(n))
    return Ric - Rs*g/2, gi
AV = [sp.Matrix([1, 1, 1]), sp.Matrix([1, -1, -1]), sp.Matrix([-1, 1, -1]), sp.Matrix([-1, -1, 1])]
def background(ks):
    Hs = [1 + k/r for k in ks]; H = sp.prod(Hs)
    g = sp.diag(-H**sp.Rational(-1, 2), H**sp.Rational(1, 2), H**sp.Rational(1, 2)*r**2, H**sp.Rational(1, 2)*r**2*sp.sin(th)**2)
    phi = sum((AV[i]*sp.log(Hs[i]) for i in range(4)), sp.zeros(3, 1))/2
    return g, phi, Hs, H
def check(ks, alpha):
    g, phi, Hs, H = background(ks); G, gi = einstein_fast(g); sq = sp.sqrt(-g.det())
    Fs = []
    for i in range(4):
        At = alpha*(1 - 1/Hs[i]); F = sp.zeros(4); F[0, 1] = -sp.diff(At, r); F[1, 0] = -F[0, 1]; Fs.append(F)
    def F2(F): Fup = gi*F*gi; return sum(F[a, b]*Fup[a, b] for a in range(4) for b in range(4))
    T = sp.zeros(4)
    for c in range(3):
        d = [sp.diff(phi[c], x) for x in X]; d2 = sum(gi[a, b]*d[a]*d[b] for a in range(4) for b in range(4))
        T += sp.Matrix(4, 4, lambda m, n: (d[m]*d[n] - g[m, n]*d2/2)/2)
    for i in range(4):
        e = sp.exp((AV[i].T*phi)[0]); F = Fs[i]; f2 = F2(F)
        T += sp.Matrix(4, 4, lambda m, n: e*(sum(F[m, k]*F[n, l]*gi[k, l] for k in range(4) for l in range(4)) - g[m, n]*f2/4)/2)
    num = lambda ex: max(abs(complex(sp.N(ex.subs({r: rv, th: 0.7}), 30))) for rv in (sp.Rational(37, 7), sp.Rational(91, 5)))
    ein = max(num(x) for x in (G - T))
    sc = 0
    for c in range(3):
        box = sum(sp.diff(sq*gi[a, b]*sp.diff(phi[c], X[b]), X[a]) for a in range(4) for b in range(4))/sq
        src = sum(AV[i][c]*sp.exp((AV[i].T*phi)[0])*F2(Fs[i]) for i in range(4))/4
        sc = max(sc, num(box - src))
    mx = 0
    for i in range(4):
        e = sp.exp((AV[i].T*phi)[0]); Fup = gi*Fs[i]*gi
        mx = max(mx, max(num(sum(sp.diff(sq*e*Fup[m, n], X[m]) for m in range(4))) for n in range(4)))
    return ein, sc, mx
if __name__ == "__main__":
    for al in [1, sp.sqrt(2), 2]:
        print(al, check([sp.Integer(1), sp.Integer(2), sp.Integer(3), sp.Rational(1, 2)], al))
