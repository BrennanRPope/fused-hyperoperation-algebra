"""Turn per lap of a probe to second order (2PN) from the exact orbit equation (du/dphi)^2 = F(u) = W/L^2.
F = c0 + c1 u - c2 u^2 + delta(u);  2 * integral du/sqrt(F) to second order in delta, with finite parts from
derivatives in c0 (u = p + q cos(theta) turns each integral into a polynomial in p and q^2).
Result: turn = pi [ T K^2/L^2 + (T2 + T2e * w) K^4/L^4 ] + O(L^-6), with E^2 = 1 - 2 w K^2/L^2
(w = (1 - e^2)/2 for a Newtonian ellipse of eccentricity e).  K = Newtonian pull."""
import sympy as sp
u, L, E = sp.symbols('u L E', positive=True)
lam, w = sp.symbols('lambda w', positive=True)
p, q2, th = sp.symbols('p q2 theta')
C = sp.symbols('c0:7')
def I(k):
    return sp.integrate(sp.expand((p + sp.sqrt(q2)*sp.cos(th))**k), (th, 0, sp.pi))
def turn_from_coeffs(cs, nmax=4):
    c0, c1, c2 = cs[0], cs[1], cs[2]
    delta = {n: cs[n] for n in range(3, nmax + 1)}
    # I[g] = c2^-1/2 sum g_k I_k(p, q2);  d/dc0 = (1/c2) d/dq2
    def Ig(poly):  # poly: dict power->coef
        return sum(cf*I(k) for k, cf in poly.items())/sp.sqrt(c2)
    d1 = {n: cf for n, cf in delta.items()}
    d2 = {}
    for a, ca in delta.items():
        for b, cb in delta.items():
            d2[a + b] = d2.get(a + b, 0) + ca*cb
    tot = sp.pi/sp.sqrt(c2) - sp.Rational(1, 2)*(-2)*sp.diff(Ig(d1), q2)/c2 + sp.Rational(3, 8)*sp.Rational(4, 3)*sp.diff(Ig(d2), q2, 2)/c2**2
    tot = tot.subs({p: c1/(2*c2), q2: c0/c2 + c1**2/(4*c2**2)})
    return 2*tot - 2*sp.pi
def turn2(W, K):
    """W exact in u, E, L.  Returns (T, T2, T2e)."""
    F = sp.series(W/L**2, u, 0, 6).removeO()
    cs = [sp.expand(F.coeff(u, n)) for n in range(6)]
    cs[2] = -cs[2]
    tr = turn_from_coeffs(cs, nmax=5)
    # orbit scaling: L = Lh/lam, E^2 = 1 - 2 w K^2 lam^2/Lh^2  ; expand in lam
    Lh = sp.symbols('Lh', positive=True)
    tr = tr.subs(E, sp.sqrt(1 - 2*w*K**2*lam**2/Lh**2)).subs(L, Lh/lam)
    ser = sp.series(tr, lam, 0, 5).removeO()
    a2 = sp.simplify(ser.coeff(lam, 2)*Lh**2/(sp.pi*K**2))
    a4 = sp.expand(sp.simplify(ser.coeff(lam, 4)*Lh**4/(sp.pi*K**4)))
    return sp.nsimplify(a2), sp.nsimplify(a4.coeff(w, 0)), sp.nsimplify(a4.coeff(w, 1))
if __name__ == "__main__":
    import sys; sys.path.insert(0, '..')
    from precess import T_numeric
    r = sp.symbols('r', positive=True)
    M = sp.Integer(1)
    A = 1 - 2*M*u; B = 1/(1 - 2*M*u)   # Schwarzschild, Schwarzschild coords
    W = (1/B)*(E**2/A - 1 - L**2*u**2)
    print("Schwarzschild:", turn2(W, M))
