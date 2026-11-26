"""Exact identities and the residue of critical branching laws (part A of parabolic-fingerprints).
Prints PASS/FAIL lines; writes identities.json."""
import json, os
import sympy as sp
HERE = os.path.dirname(os.path.abspath(__file__))
out, ok_all = {}, True
def rep(name, ok, extra=""):
    global ok_all
    ok_all &= bool(ok)
    print(("PASS " if ok else "FAIL ") + name + (("   " + extra) if extra else ""))

s, t, y, z, x, u, Q, r = sp.symbols("s t y z x u Q r")

# 1. diamond hierarchical lattice (map from ../diamond-lattice/dynamics.py): T(t, Q) = lens(Q-2, 1-Q, t^2, t^2)
lens = lambda p, q, a, b: (a + b + p*a*b)/(1 - q*a*b)
T = lens(Q - 2, 1 - Q, t**2, t**2)
T0 = sp.simplify(T.subs(Q, 0))
rep("diamond lattice at Q = 0: T(t) = 2t^2/(1+t^2)", sp.simplify(T0 - 2*t**2/(1 + t**2)) == 0, str(T0))
Ty = sp.simplify(1/T0.subs(t, 1/y))            # in y = 1/t
rep("in y = 1/t the Q = 0 diamond step is y -> (y^2+1)/2, the generating function of critical binary splitting",
    sp.simplify(Ty - (y**2 + 1)/2) == 0, str(Ty))
rep("and y = 2z turns it into z -> z^2 + 1/4, the cusp of the Mandelbrot set",
    sp.simplify(Ty.subs(y, 2*z)/2 - (z**2 + sp.Rational(1, 4))) == 0)
lg = r*x*(1 - x)
rep("logistic map at r = 1 is the same cusp map: x = 1/2 - z gives z -> z^2 + 1/4",
    sp.simplify(sp.Rational(1, 2) - lg.subs({r: 1, x: sp.Rational(1, 2) - z}) - (z**2 + sp.Rational(1, 4))) == 0)

# 2. our class: Poisson branching and tetration at the kiss base eta = e^(1/e)
rep("critical Poisson(1) branching: f(s) = e^(s-1) is the rung step e^x - 1 with x = s - 1",
    sp.simplify(sp.exp((1 + x) - 1) - 1 - (sp.exp(x) - 1)) == 0)
eta_map = sp.exp((sp.E + sp.E*x)/sp.E)          # eta^X with X = e + e x
rep("tetration base eta = e^(1/e) near its kiss point X = e: X = e(1 + x) turns X -> eta^X into x -> e^x - 1 exactly",
    sp.simplify(eta_map/sp.E - 1 - (sp.exp(x) - 1)) == 0)

# 3. geometric branching is Moebius
geo = 1/(2 - s)
rep("critical geometric branching f(s) = 1/(2-s) is Moebius: x = s - 1 gives x -> x/(1 - x)",
    sp.simplify(geo.subs(s, 1 + x) - 1 - x/(1 - x)) == 0)

# 4. iterative residue of a critical branching law: resit = 1 - (2/3) kappa,  kappa = E[Z(Z-1)(Z-2)] / sigma^4
def resit_pgf(f):
    g = sp.series(f.subs(s, 1 + x) - 1, x, 0, 4).removeO()
    a1, a2 = g.coeff(x, 2), g.coeff(x, 3)
    return sp.nsimplify(sp.simplify(1 - a2/a1**2)), sp.simplify(2*a1), sp.simplify(6*a2)  # resit, sigma^2, 3rd factorial moment
laws = {
    "binary splitting (0 or 2, half each)": (s**2 + 1)/2,
    "Poisson(1)": sp.exp(s - 1),
    "geometric, mean 1": 1/(2 - s),
    "negative binomial r = 2, mean 1": (sp.Rational(2, 3)/(1 - sp.Rational(1, 3)*s))**2,
    "negative binomial r = 3, mean 1": (sp.Rational(3, 4)/(1 - sp.Rational(1, 4)*s))**3,
    "binomial(2, 1/2)": ((1 + s)/2)**2,
    "binomial(3, 1/3)": ((2 + s)/3)**3,
    "0 or 3 children (2/3, 1/3)": (2 + s**3)/3,
}
out["branching"] = {}
for name, f in laws.items():
    res, sig2, m3 = resit_pgf(f)
    kappa = m3/sig2**2
    good = sp.simplify(res - (1 - sp.Rational(2, 3)*kappa)) == 0
    out["branching"][name] = {"resit": str(res), "sigma2": str(sig2), "kappa": str(sp.nsimplify(kappa)), "resit_float": float(res)}
    rep(f"resit = 1 - (2/3) kappa for {name}", good, f"sigma^2 = {sig2}, kappa = {sp.nsimplify(kappa)}, resit = {res}")
# negative binomial with mean 1 and shape r: resit = (1 - r)/(3 r)? derive symbolically
nb = (r/(r + 1))**r/(1 - s/(r + 1))**r
g = sp.series(nb.subs(s, 1 + x) - 1, x, 0, 4).removeO()
res_nb = sp.simplify(1 - g.coeff(x, 3)/g.coeff(x, 2)**2)
print("negative binomial, mean 1, shape r: resit =", sp.factor(res_nb), "(r -> infinity:", sp.limit(res_nb, r, sp.oo), "= Poisson; r = 1:", res_nb.subs(r, 1), "= geometric)")
out["negbin_resit"] = str(sp.factor(res_nb))
rep("negative binomial family runs from geometric (resit 0, r = 1) to Poisson (resit 1/3, r -> infinity)",
    res_nb.subs(r, 1) == 0 and sp.limit(res_nb, r, sp.oo) == sp.Rational(1, 3))

json.dump(out, open(os.path.join(HERE, "identities.json"), "w"), indent=1)
print("ALL PASS" if ok_all else "SOME FAILED")
