"""Independent check of the scalar-tensor (DEF) 1PN formula in the probe limit, from the exact JNW field."""
import sys; sys.path.insert(0, '..')
import sympy as sp
from fields import jnw, r
from precess import T_series, u, L, E
gam = sp.Rational(3, 5); b = sp.Integer(2)
g, phi, At, M, S = jnw(gam, b)
aA = sp.sqrt(1 - gam**2)/gam           # phi ~ -aA M / r
for aB, bB in [(sp.Rational(1, 2), 0), (sp.Rational(1, 3), sp.Rational(2)), (-sp.Rational(1, 4), sp.Rational(-3, 2))]:
    m = sp.exp(aB*phi + bB*phi**2/2)
    A = -g[0, 0]; B = g[1, 1]; C = g[2, 2]/r**2
    W = sp.simplify(((C**2/B)*(E**2/A - m**2 - L**2/(r**2*C))).subs(r, 1/u))
    T, K = T_series(W)
    s = aA*aB
    def_formula = 6*((1 - s/3)/(1 + s) - bB*aA**2/(6*(1 + s)**2))
    print(aB, bB, "probe T =", sp.nsimplify(T), " DEF =", sp.nsimplify(sp.simplify(def_formula)), " K/M =", sp.simplify(K/M), "vs 1+s =", sp.simplify(1+s))
