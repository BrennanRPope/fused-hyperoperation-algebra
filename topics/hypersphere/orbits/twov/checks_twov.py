"""DRAFT checks for hypersphere/orbits/twov: where gravity's 'velocity disc at 2v' comes from.
Run: python3 checks_twov.py  (sympy, numpy, scipy). Writes results_twov.json for the figure.
Finding: the 4v² is 1·v² + 3·v². The 1 is special relativity, the same as the electric orbit. The 3 is gravity
pulling on energy: a falling particle's local energy grows as 1/√A, so its potential E(1 − 1/√A) is steeper than 1/r.
"""
import json, math, os
import numpy as np, sympy as sp
from scipy import integrate
HERE = os.path.dirname(os.path.abspath(__file__))
results, store = [], {}
def check(name, ok, detail):
    results.append((name, bool(ok))); print(f"[{'PASS' if ok else 'FAIL'}] {name}: {detail}")

r, E, M = sp.symbols('r E M', positive=True)
A = sp.Function('A')(r)
Vg = E * (1 - A**sp.Rational(-1, 2))
check("D1 a geodesic in a clock-warp spacetime IS a special-relativity orbit in the potential E(1 − 1/√A)",
      sp.simplify((E - Vg)**2 - E**2 / A) == 0, "radial equations identical: (E − V)² − 1 − L²/r² = E²/A − 1 − L²/r²")
# D2: in the special-relativity picture, a circular orbit's speed L/(γ r) with γ = E/√A equals the local speed √(rA'/2A)
Mn = 1.0; worst2 = 0
for r0 in (7.0, 12.0, 40.0):
    A0 = 1 - 2 * Mn / r0; dA0 = 2 * Mn / r0**2
    L2 = r0**3 * dA0 / (2 * A0 - r0 * dA0); E0 = math.sqrt(A0 * (1 + L2 / r0**2))
    v_sr = math.sqrt(L2) / ((E0 / math.sqrt(A0)) * r0); v_loc = math.sqrt(r0 * dA0 / (2 * A0))
    worst2 = max(worst2, abs(v_sr - v_loc))
check("D2 that picture's speed is the local speed a static observer measures (its γ = E/√A, the local energy)",
      worst2 < 1e-12, f"circular orbits at r = 7, 12, 40 M; worst difference {worst2:.1e}")
v2 = r * sp.diff(A, r) / (2 * A)
NV = 3 + r * sp.diff(Vg, r, 2) / sp.diff(Vg, r); NA = 3 + r * sp.diff(A, r, 2) / sp.diff(A, r)
check("D3 the energy potential's closure number: N(E(1 − 1/√A)) = N(A) − 3v², any clock rate A",
      sp.simplify(NV - (NA - 3 * v2)) == 0, "symbolic")
As = 1 - 2 * M / r
NVs = sp.simplify(NV.subs(A, As).doit()); v2s = sp.simplify(v2.subs(A, As).doit())
check("D4 Schwarzschild: N = 1 − 3v² (Newton's closed orbits spoiled by 3v²)", sp.simplify(NVs - (1 - 3 * v2s)) == 0,
      f"N = {NVs}, v² = {v2s}")
check("D5 so 1 − 6M/r = (1 − 2M/r) · [(1 − 3v²) − v²]: space × (energy pull 3v² + special relativity v²)",
      sp.simplify((1 - 2 * M / r) * (1 - 3 * v2s - v2s) - (1 - 6 * M / r)) == 0, "identity")
ser = sp.series(1 - As**sp.Rational(-1, 2), M, 0, 3).removeO()
check("D6 the energy potential is steeper than Newton's: 1 − 1/√(1 − 2M/r) = −M/r − (3/2)M²/r² − …",
      sp.simplify(ser - (-M / r - sp.Rational(3, 2) * M**2 / r**2)) == 0, f"series {ser}")

# D7: simulate both pictures for a slightly eccentric orbit and compare the perihelion advance
Mv = 1.0
def Af(x): return 1 - 2 * Mv / x
def dAf(x): return 2 * Mv / x**2
def geodesic(r0, kick=0.02, laps=5):
    A0, dA0 = Af(r0), dAf(r0)
    L2 = r0**3 * dA0 / (2 * A0 - r0 * dA0) * (1 + kick); Lc = math.sqrt(L2); Ec = math.sqrt(A0 * (1 + L2 / r0**2))
    def rhs(t, y):
        x, p, ph = y
        return [p, -0.5 * (Ec**2 * dAf(x) / Af(x)**2 - 2 * L2 / x**3), Lc / x**2]
    T = 2 * math.pi * r0**2 / Lc
    s = integrate.solve_ivp(rhs, (0, laps * T * 1.4), [r0, 0.0, 0], rtol=1e-12, atol=1e-14, max_step=T / 500)
    return s, Ec, Lc
def sr_orbit(r0, Ec, Lc, kick=0.02, laps=5):
    # H = sqrt(1 + p^2) + V(r), V = Ec(1 - 1/sqrt(A)); use proper-time-like parameter: integrate in coordinate time of the SR picture
    def V1(x): return Ec * 0.5 * Af(x)**-1.5 * dAf(x)
    def rhs(t, y):
        x, pr, ph = y
        g = math.sqrt(1 + pr**2 + Lc**2 / x**2)
        return [pr / g, Lc**2 / (g * x**3) - V1(x), Lc / (g * x**2)]
    T = 2 * math.pi * r0**2 / Lc * 3
    s = integrate.solve_ivp(rhs, (0, laps * T * 1.4), [r0, 0.0, 0], rtol=1e-12, atol=1e-14, max_step=T / 1500)
    return s
def advance(s, pidx=1):
    p = s.y[pidx]; ph = s.y[2]
    idx = np.where((p[:-1] < 0) & (p[1:] >= 0))[0]  # periapsis (start is apoapsis: L raised)
    peri = [ph[i] + (ph[i + 1] - ph[i]) * (-p[i] / (p[i + 1] - p[i])) for i in idx]
    return np.diff(peri).mean() - 2 * math.pi
worst = 0; pairs = []
for r0 in (12.0, 20.0, 40.0):
    sg, Ec, Lc = geodesic(r0)
    ss = sr_orbit(r0, Ec, Lc)
    a1, a2 = advance(sg), advance(ss)
    worst = max(worst, abs(a1 - a2) / abs(a1)); pairs.append((r0, a1, a2))
check("D7 simulated: the clock-warp geodesic and the special-relativity orbit in E(1 − 1/√A) precess identically",
      worst < 1e-6, f"eccentric orbits (L raised 2% at r = 12, 20, 40 M); worst relative difference {worst:.1e}")
store["pairs"] = pairs

# D8: first-order budget in units of hodograph area / c^2
x0 = 1e-7; vv = x0 / (1 - 2 * x0)
parts = {"space": 2 * (1 / math.sqrt(1 - 2 * x0) - 1) / vv,
         "special relativity": 2 * (1 / math.sqrt(1 - vv) - 1) / vv,
         "energy pull": 2 * (1 / math.sqrt(1 - 3 * vv) - 1) / vv}
check("D8 first-order budget: orbit 6 = 2 (space) + 1 (special relativity) + 3 (energy pull); gyroscope 3 = 2 + 1; electric 1",
      all(abs(parts[k] - t) < 1e-5 for k, t in zip(parts, (2, 1, 3))), {k: round(v, 5) for k, v in parts.items()})
store["parts"] = parts
json.dump(store, open(os.path.join(HERE, "results_twov.json"), "w"), indent=1, default=float)
print(); print(f"{sum(o for _, o in results)} of {len(results)} checks pass")
