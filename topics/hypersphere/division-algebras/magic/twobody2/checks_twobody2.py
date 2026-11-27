"""Checks: two comparable black holes at second order (2PN), from Gupta's Einstein-Maxwell-dilaton Lagrangian.
Writes checks_twobody2_output.txt and results_twobody2.json. Bodies: extremal dilaton black holes with coupling a and
opposite charges; a = sqrt3, 1, 1/sqrt3, 0 are the four-charge family's black holes with n = 1, 2, 3, 4 equal charges."""
import sys, os, json; HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE); sys.path.insert(0, os.path.join(HERE, '..', 'order2')); sys.path.insert(0, os.path.join(HERE, '..'))
import sympy as sp
from pn2 import run
from poly import poly_F
from precess import u, L, E
R = sp.Rational; lines = []; ok_all = True; data = {}
def rep(name, ok, detail=""):
    global ok_all; ok_all &= bool(ok); lines.append(f"{'PASS' if ok else 'FAIL'}  {name}  {detail}"); print(lines[-1], flush=True)
masses = [1, R(1, 2), R(1, 3), R(1, 9)]
def nu_of(m2): return R(m2)/(1 + m2)**2
def fit(case, **kw):
    pts = []
    for m2 in masses:
        T, T2, T2e, rest = run(case, 1, m2, **kw); assert rest == 0
        pts.append((nu_of(m2), sp.nsimplify(T), sp.nsimplify(T2), sp.nsimplify(T2e)))
    n0, n1 = pts[0][0], pts[1][0]
    lin = lambda i: ((pts[0][i]*n1 - pts[1][i]*n0)/(n1 - n0), (pts[0][i] - pts[1][i])/(n0 - n1))
    out = {nm: tuple(sp.nsimplify(sp.simplify(x)) for x in lin(i)) for nm, i in [("T", 1), ("T2", 2), ("T2e", 3)]}
    exact = all(sp.simplify(p[i] - (out[nm][0] + out[nm][1]*p[0])) == 0 for p in pts for nm, i in [("T", 1), ("T2", 2), ("T2e", 3)])
    return out, exact, pts
# W1 Einstein: Damour-Schafer 2PN advance with the mass ratio, T = 6, T2 = 105/2 - 15 nu, T2e = -15 + 6 nu
gr, ex, pts = fit('GR')
ok = ex and gr["T"] == (6, 0) and gr["T2"] == (R(105, 2), -15) and gr["T2e"] == (-15, 6)
rep("W1 the pipeline reproduces Einstein's second-order two-body advance (Damour-Schafer) at four mass ratios", ok, str(gr))
data["GR"] = {k: [str(x) for x in v] for k, v in gr.items()}
# probe values from the exact orbit equation (section 11): n equal charges, visitor opposite on the same n charges
def probe(n):
    one = sp.Integer(1); ks = [one]*n + [0]*(4 - n); ws = [one/n]*n + [0]*(4 - n); ss = [-1]*n + [1]*(4 - n)
    P = sp.expand(poly_F(ks, ws, ss) + L**2*u**2); c = [P.coeff(u, m) for m in range(5)]
    K = c[1].subs(E, 1)/2; T = c[2].subs(E, 1)/K**2
    return sp.nsimplify(T), sp.nsimplify(R(3, 4)*T**2 + 3*c[3].subs(E, 1)/K**3), sp.nsimplify(-sp.diff(c[2], E).subs(E, 1)/K**2)
res = {}; okW2 = True; okW4 = True; okW3 = None
for n, av in [(1, sp.sqrt(3)), (2, 1), (3, 1/sp.sqrt(3)), (4, 0)]:
    out, ex, pts = fit('EMD', a=av); pr = probe(n)
    okW2 &= ex and all(sp.simplify(x - y) == 0 for x, y in zip((out["T"][0], out["T2"][0], out["T2e"][0]), pr))
    T = out["T"][0]
    if n == 1: okW3 = ex and all(v == (0, 0) for v in out.values())
    else: okW4 &= ex and sp.simplify(out["T"][1]) == 0 and sp.simplify(out["T2"][1] - R(3, 2)*T) == 0 and sp.simplify(out["T2e"][1] - T) == 0

    res[n] = dict(a=str(av), probe=[str(x) for x in pr], **{k: [str(x) for x in v] for k, v in out.items()}, points=[[str(x) for x in p] for p in pts])
data["EMD_extremal_opposite"] = res
rep("W2 the paper's Lagrangian, as transcribed, reproduces our exact small-visitor results (section 11) for n = 1, 2, 3, 4 charges, velocity and acceleration terms included", okW2,
    "; ".join(f"n={n}: probe {res[n]['probe']}" for n in res))
rep("W3 two one-charge black holes still close at second order, at every mass ratio tested", okW3, f"n=1: {res[1]['T']}, {res[1]['T2']}, {res[1]['T2e']}")
rep("W4 for n = 2, 3, 4 the mass ratio enters only as (3/2) nu T and nu T: proportional to the first-order turn from pairs", okW4,
    "; ".join(f"n={n}: T2 = {res[n]['T2'][0]} + {res[n]['T2'][1]} nu, T2e = {res[n]['T2e'][0]} + {res[n]['T2e'][1]} nu" for n in (2, 3, 4)))
rep("W5 the rule is special to these extremal pairs: in Einstein's theory the second-order mass-ratio term is -15 nu, not (3/2) x 6 nu", gr["T2"][1] != R(3, 2)*6, f"GR nu term {gr['T2'][1]}")
# W6 simulated: integrate the second-order equal-mass orbits and measure the turn per lap
import numpy as np
from orbit2pn import build, measure
sim = {}; ok6 = True; wv = (1 - 0.3**2)/2
for n, av in [(1, sp.sqrt(3)), (2, 1), (4, 0)]:
    mdl = build('EMD', a=av); T = sp.sympify(res[n]["T"][0]); c2 = sp.sympify(res[n]["T2"][0]) + sp.sympify(res[n]["T2"][1])/4 + (sp.sympify(res[n]["T2e"][0]) + sp.sympify(res[n]["T2e"][1])/4)*sp.nsimplify(wv)
    got = []
    for Jk in [8, 16]:
        t = measure(mdl, Jk, 0.3, 6)[0]; got.append((t - np.pi*float(T)/Jk**2)/(np.pi/Jk**4) if n > 1 else t)
    if n == 1: ok6 &= abs(got[0]) < 1e-5 and abs(got[1]) < abs(got[0])/40
    else: ok6 &= abs(got[1] - float(c2)) < 0.02*abs(float(c2)) and abs(got[1] - float(c2)) < abs(got[0] - float(c2))
    sim[n] = dict(pred=float(c2) if n > 1 else 0.0, J8=got[0], J16=got[1])
data["simulated"] = sim
rep("W6 simulated equal-mass orbits match: one charge turns only at the 3PN truncation level (falls ~64x as J doubles); n = 2 and 4 land on the second-order coefficient", ok6,
    "; ".join(f"n={n}: pred {v['pred']:.4g}, J/k=8: {v['J8']:.4g}, J/k=16: {v['J16']:.4g}" for n, v in sim.items()))
lines.append("ALL PASS" if ok_all else "SOME FAILED"); print(lines[-1])
json.dump(dict(passed=ok_all, data=data), open(os.path.join(HERE, "results_twobody2.json"), "w"), indent=1)
open(os.path.join(HERE, "checks_twobody2_output.txt"), "w").write("\n".join(lines) + "\n")
