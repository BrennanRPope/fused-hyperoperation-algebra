"""Measure H and the affine normal at O on the chord surface for a spread of settings (for figures)."""
import sympy as sp, json, random
from blaschke import Surface
from checks_lib import chord_surface, u, v
random.seed(11)
out = []
ks = [-2, -1, -sp.Rational(1, 2), 0, sp.Rational(1, 2), 1, sp.Rational(3, 2), 2]
for i in range(22):
    k = sp.Rational(random.randint(-20, 20), 10)
    d = sp.Rational(random.randint(-10, 15), 10)
    S = Surface(chord_surface(k, d), [u, v], {u: 0, v: 0})
    X0, xi0, _, _ = S.at({u: 0, v: 0}, prec=20)
    X1, xi1, _, _ = S.at({u: sp.Rational(1, 10), v: sp.Rational(1, 20)}, prec=20)
    H = -(xi1 - xi0).dot(X1 - X0) / (X1 - X0).dot(X1 - X0)
    out.append(dict(k=float(k), d=float(d), H=float(H), xi0=[float(e) for e in xi0]))
    print(out[-1], flush=True)
json.dump(out, open('measured.json', 'w'), indent=1)
