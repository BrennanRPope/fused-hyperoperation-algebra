"""Part B: iterative residue and horn-map size for named systems that creep into a fixed point.
Uses hornlib.py (copied from ../height-flow/horn_invariants.py). Usage: horns.py NAME  (or 'all'). Writes horns_NAME.json."""
import sys, os, json, math
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from hornlib import Germ, horn, Fr, mp, N
import sympy as sp
mp.mp.dps = 34
X = sp.symbols("X")

def series_co(expr, var=X):
    ser = sp.series(expr, var, 0, N + 1).removeO()
    return [Fr(int(sp.fraction(sp.nsimplify(ser.coeff(var, i)))[0]), int(sp.fraction(sp.nsimplify(ser.coeff(var, i)))[1])) for i in range(N + 1)]

def newton_inv(f, guess_fn):
    def finv(y):
        return mp.findroot(lambda w: f(w) - y, guess_fn(y))
    return finv

germs = {}
# diamond lattice at Q = 0, in its own variable t near the cold point t = 1 (w = t - 1)
germs["diamond_Q0"] = ("diamond lattice Potts step at Q = 0 (t near 1)",
    series_co(2*(1 + X)**2/(1 + (1 + X)**2) - 1),
    lambda w: 2*(1 + w)**2/(1 + (1 + w)**2) - 1,
    lambda v: mp.sqrt((1 + v)/(1 - v)) - 1)
# critical branching laws, x = s - 1
def pgf_germ(name, pgf_sym, pgf_num):
    co = series_co(pgf_sym.subs(X, X) )
    f = lambda w: pgf_num(1 + w) - 1
    return (name, co, f, newton_inv(f, lambda y: y - y**2/2 + y**3/3))
S = sp.symbols("S")
def mk(name, expr_s, num):
    return (name, series_co(expr_s.subs(S, 1 + X) - 1), lambda w: num(1 + w) - 1,
            None)
laws = {
    "poisson": ("critical Poisson(1) branching (= e^x - 1)", sp.exp(S - 1), lambda s_: mp.exp(s_ - 1), lambda y: 1 + mp.log1p(y)),
    "binary": ("critical binary splitting", (S**2 + 1)/2, lambda s_: (s_**2 + 1)/2, lambda y: mp.sqrt(2*y - 1)),
    "negbin2": ("negative binomial r = 2, mean 1", (sp.Rational(2, 3)/(1 - S/3))**2, lambda s_: (mp.mpf(2)/3/(1 - s_/3))**2,
                lambda y: 3*(1 - (mp.mpf(2)/3)/mp.sqrt(y))),
    "binom3": ("binomial(3, 1/3)", ((2 + S)/3)**3, lambda s_: ((2 + s_)/3)**3, lambda y: 3*mp.cbrt(y) - 2),
    "poissonlike": ("law (1/3, 1/2, 0, 1/6): same sigma^2 and residue as Poisson", (2 + 3*S + S**3)/6,
                    lambda s_: (2 + 3*s_ + s_**3)/6, None),
}
for key, (name, e, num, inv_s) in laws.items():
    f = (lambda num_: (lambda w: num_(1 + w) - 1))(num)
    if inv_s is not None:
        finv = (lambda inv_: (lambda v: inv_(1 + v) - 1))(inv_s)
    else:
        finv = newton_inv(f, lambda y: y - y**2/2 + y**3/3)
    germs[key] = (name, series_co(e.subs(S, 1 + X) - 1), f, finv)

# type-I intermittency: logistic map cubed at the period-3 saddle-node, r = 1 + sqrt(8)
def logistic3_germ(which):
    mp.mp.dps = 80
    rr = 1 + mp.sqrt(8)
    lg = lambda v: rr*v*(1 - v)
    f3 = lambda v: lg(lg(lg(v)))
    # tangent points: f3(v) = v and f3'(v) = 1; known window opens with orbit near 0.159, 0.512, 0.956
    guesses = [mp.mpf("0.1599"), mp.mpf("0.5143"), mp.mpf("0.9563")]
    dlg = lambda v: rr*(1 - 2*v)
    df3 = lambda v: dlg(v)*dlg(lg(v))*dlg(lg(lg(v)))
    # f3(v) - v has a double root at the tangency, so solve the simple root of f3'(v) = 1 instead
    v0 = mp.findroot(lambda v: df3(v) - 1, guesses[which])
    tay = mp.taylor(lambda w: f3(v0 + w) - v0, 0, N)
    co = [Fr(0), Fr(1)] + [Fr(mp.nstr(c, 70)) for c in tay[2:]]
    resid_fp, resid_d = abs(f3(v0) - v0), abs(tay[1] - 1)
    mp.mp.dps = 34
    v0s = mp.mpf(v0)
    f = lambda w: f3(v0s + w) - v0s
    finv = newton_inv(f, lambda y: y - co2f * y**2)
    co2f = mp.mpf(co[2].numerator)/co[2].denominator
    return (f"logistic map cubed at the period-3 tangency, orbit point {mp.nstr(v0s, 8)}", co, f, finv), (resid_fp, resid_d)

def run(key):
    extra = {}
    if key.startswith("logistic3_"):
        g_args, (r1, r2) = logistic3_germ(int(key[-1]))
        extra = {"fixed_point_residual": mp.nstr(r1, 3), "derivative_minus_1": mp.nstr(r2, 3)}
    else:
        g_args = germs[key]
    name, co, f, finv = g_args
    g = Germ(name, co, f, finv)
    rho = g.rho
    A, r1_, _ = horn(g, 3.0, +1)
    Ab, rb, _ = horn(g, 2.5, +1)
    Am, rm, _ = horn(g, 3.0, -1)
    prod = A[1]*Am[1]
    size = mp.sqrt(abs(prod))
    spread = abs(A[1] - Ab[1])/max(abs(A[1]), mp.mpf(10)**-30)
    res = {"name": name, "resit": str(rho) if len(str(rho)) < 40 else mp.nstr(mp.mpf(rho.numerator)/rho.denominator, 20),
           "resit_float": float(rho), "a1": float(g.a1),
           "A1": mp.nstr(A[1], 14), "A_-1": mp.nstr(Am[1], 14), "A0+": mp.nstr(A[0], 6), "A0-": mp.nstr(Am[0], 10),
           "size": mp.nstr(size, 12), "size_float": float(size), "strip_spread": mp.nstr(spread, 3),
           "fit_resid": mp.nstr(max(r1_, rm), 3)}
    res.update(extra)
    print(json.dumps(res), flush=True)
    json.dump(res, open(os.path.join(HERE, f"horns_{key}.json"), "w"), indent=1)

if __name__ == "__main__":
    keys = sys.argv[1:]
    if keys == ["all"]:
        keys = list(germs) + ["logistic3_0", "logistic3_1", "logistic3_2"]
    for k in keys:
        run(k)
