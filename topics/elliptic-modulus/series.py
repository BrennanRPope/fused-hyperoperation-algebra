"""Exact power-series tools for the elliptic modulus page (DRAFT).
Slope-chart generator G_m(u) = sn(u|m) dn(u|m) / cn(u|m)  (the Edwards slope x/y, F3.1 / F3.5)."""
import sympy as sp
u, m, s, x = sp.symbols('u m s x')
def trunc(e, N): return sp.series(sp.expand(e), u, 0, N).removeO()
def jacobi_series(N=14, mval=m):
    sn, cn, dn = u, sp.Integer(1), sp.Integer(1)
    for _ in range(N + 2):
        sn, cn, dn = (trunc(sp.integrate(trunc(cn*dn, N), u), N),
                      trunc(1 - sp.integrate(trunc(sn*dn, N), u), N),
                      trunc(1 - sp.integrate(trunc(mval*sn*cn, N), u), N))
    return sp.expand(sn), sp.expand(cn), sp.expand(dn)
def slope_series(N=14, mval=m):
    sn, cn, dn = jacobi_series(N, mval)
    inv = trunc(sp.series(1/cn, u, 0, N).removeO(), N)
    return sp.expand(trunc(sn*dn*inv, N))
def compose(f, g, N):  # f(g(u)), both series in u
    return sp.expand(trunc(f.subs(u, g), N))
def julia(f, N):
    """Iterative logarithm V (vector field) with exp(V d/du) u = f, as series to order N."""
    cs = sp.symbols('c3:%d' % N)
    V = sum(cs[i-3]*u**i for i in range(3, N))
    # flow map: sum_{n} (V d/du)^n u / n!
    term, total = u, u
    for n in range(1, N):
        term = sp.expand(trunc(V*sp.diff(term, u), N)) / n
        total += term
        if term == 0: break
    eqs = sp.Poly(sp.expand(trunc(total - f, N)), u).all_coeffs()
    sol = sp.solve([sp.expand(trunc(total - f, N)).coeff(u, i) for i in range(3, N)], cs, dict=True)[0]
    return sp.expand(V.subs(sol))
def flow(V, t, N):
    term, total = u, u
    for n in range(1, N):
        term = sp.expand(trunc(t*V*sp.diff(term, u), N)) / n
        total += term
    return sp.expand(trunc(total, N))
