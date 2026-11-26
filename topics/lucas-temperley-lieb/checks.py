"""Lucas weights U_m(k, d) against the Temperley-Lieb algebra built from diagrams (tl.py). Writes checks.json."""
import json, os, sys
import sympy as sp
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
from tl import TL, matchings
d, k, dd = sp.symbols("delta k d")
out, okall = {}, True
def rep(name, ok, extra=""):
    global okall
    okall &= bool(ok); print(("PASS " if ok else "FAIL ") + name + (("   " + extra) if extra else ""), flush=True)

def lucas(m, kk, dv):
    U = [sp.Integer(0), sp.Integer(1)]
    for _ in range(m):
        U.append(sp.expand(kk*U[-1] - dv*U[-2]))
    return U[m]

# 1. Jones-Wenzl projector from its definition alone (an idempotent 1 + ... killed by every e_i), no Chebyshev used
out["jw_traces"] = {}
for m in range(1, 6):
    T = TL(m, d)
    others = [b for b in T.basis if b != T.one]
    cs = sp.symbols(f"c0:{len(others)}")
    P = {T.one: sp.Integer(1)}
    for b, c in zip(others, cs):
        P[b] = c
    eqs = []
    for i in range(m - 1):
        prod = T.mul({T.e(i): 1}, P)
        eqs += list(prod.values())
    sol = sp.solve(eqs, cs, dict=True) if others else [{}]
    assert len(sol) == 1, m
    Pm = {b: sp.simplify(sp.sympify(c).subs(sol[0])) for b, c in P.items()}
    idem = all(sp.simplify(v) == 0 for v in T.add(T.mul(Pm, Pm), Pm, -1).values())
    tr = sp.factor(sp.simplify(T.trace(Pm)))
    U = lucas(m + 1, d, 1)
    out["jw_traces"][m] = str(tr)
    rep(f"Jones-Wenzl p_{m} (built from 'idempotent, killed by every cup-cap'): idempotent, and its closed trace equals the Lucas weight U_{m+1}(k = delta, d = 1)",
        idem and sp.simplify(tr - U) == 0, f"trace = {tr}")

# 2. general d is the same thing rescaled: U_m(k, d) = d^((m-1)/2) U_m(k/sqrt(d), 1)
ok = all(sp.simplify(lucas(m, k, dd) - dd**sp.Rational(m - 1, 2)*lucas(m, k/sp.sqrt(dd), 1)) == 0 for m in range(1, 12))
rep("U_m(k, d) = d^((m-1)/2) * U_m(k / sqrt(d), 1) for m = 1..11: the plane dial (k, d) is loop value delta = k / sqrt(d)", ok)

# 3. Gram (meander) determinant of TL_n: its zeros are exactly the silent layers U_m(delta, 1) = 0, m <= n + 1
def posroots(expr):
    rs = set()
    for fac, _ in sp.factor_list(sp.expand(expr), d)[1]:
        for r in sp.Poly(fac, d).real_roots():
            if float(r) > 1e-12:
                rs.add(round(float(r), 12))
    return sorted(rs)
out["gram"] = {}
xs_ = sp.Symbol("x")
for n in range(1, 7):
    silent = sorted(set(r for m in range(2, n + 2) for r in posroots(lucas(m, d, 1))))
    if n <= 3:   # exact determinant
        G = sp.Matrix(TL(n, d).gram())
        det = sp.factor(G.det(method="berkowitz"))
        roots = posroots(det)
        same = len(roots) == len(silent) and all(abs(a - b) < 1e-9 for a, b in zip(roots, silent))
        out["gram"][n] = {"det": str(det), "positive_roots": roots}
        rep(f"TL_{n} Gram determinant (exact) vanishes exactly at the positive silent layers of U_2..U_{n+1}", same,
            "roots " + ", ".join(f"{r:.6f}" for r in roots))
    else:        # numeric: smallest |eigenvalue| is ~0 at each silent layer and bounded away from 0 elsewhere
        E = np.array([[sp.Poly(g, xs_).degree() for g in row] for row in TL(n, xs_).gram()], dtype=float)
        def smallest(x):
            ev = np.linalg.eigvalsh(x ** E); return np.abs(ev).min() / np.abs(ev).max()
        at_roots = max(smallest(r) for r in silent)
        grid = [x for x in np.linspace(0.05, 2.5, 2451) if min(abs(x - r) for r in silent) > 0.01]
        away = min(smallest(x) for x in grid)
        out["gram"][n] = {"positive_roots": silent, "smallest_at_roots": at_roots, "smallest_away": away}
        rep(f"TL_{n} Gram matrix (numeric, {len(E)} x {len(E)}) is singular exactly at the positive silent layers of U_2..U_{n+1}",
            at_roots < 1e-12 and away > 1e-9, f"smallest relative |eig|: {at_roots:.1e} at the layers, {away:.1e} elsewhere")

# 4. positivity: Gram matrix of TL_n is positive semidefinite exactly for delta >= 2cos(pi/(n+1)) (scan)
out["psd"] = {}
grid = np.linspace(0.02, 2.4, 1191)
for n in range(2, 7):
    T = TL(n, sp.Symbol("x"))
    B = T.basis
    # precompute loop exponents once
    E = np.array([[sp.Poly(g, sp.Symbol("x")).degree() for g in row] for row in T.gram()], dtype=float)
    mins = []
    for x in grid:
        mins.append(np.linalg.eigvalsh(x ** E).min() / x ** E.max())
    mins = np.array(mins)
    thr = 2*np.cos(np.pi/(n + 1))
    # PSD (within rounding) above the threshold, clearly negative somewhere just below it
    above = mins[grid > thr + 1e-3].min() > -1e-10
    below = mins[(grid < thr - 1e-3) & (grid > thr - 0.15)].max() < 0
    out["psd"][n] = {"threshold": float(thr), "grid": grid.tolist(), "min_eig": mins.tolist()}
    rep(f"TL_{n}: Gram matrix positive above delta = 2cos(pi/{n + 1}) = {thr:.6f} and not positive just below it", above and below)

# 5. Lucas weights: U_1..U_n >= 0 exactly when k >= 2cos(pi/(n+1)) (d = 1)
ok = True
for n in range(2, 9):
    thr = 2*np.cos(np.pi/(n + 1))
    for kv in np.linspace(0.01, 2.3, 400):
        allpos = all(float(lucas(m, sp.Float(kv), 1)) >= -1e-12 for m in range(1, n + 2))
        if allpos != (kv >= thr - 1e-9):
            ok = False
rep("Lucas weights U_1..U_{n+1} are all >= 0 exactly when k >= 2cos(pi/(n+1)), the same threshold as the Gram matrix (n = 2..8)", ok)

# 6. four-ended versus six-ended tangles
x, y = sp.symbols("x y")
# TL_2 = span{1, e}: N = x + y e satisfies N^2 = kN - D with k^2 - 4D = (delta y)^2 >= 0
T2 = TL(2, d); N2 = {T2.one: x, T2.e(0): y}
kk = 2*x + y*d; DD = x*(x + y*d)
lhs = T2.mul(N2, N2); rhs = T2.add(T2.scal(kk, N2), {T2.one: -DD})
ok2 = all(sp.simplify(v) == 0 for v in T2.add(lhs, rhs, -1).values()) and sp.simplify(kk**2 - 4*DD - (d*y)**2) == 0
rep("4-ended tangles (TL_2): every element obeys N^2 = kN - D with k^2 - 4D = (delta*y)^2 >= 0, so with real coefficients they never reach inside the parabola (elliptic)", ok2)
T3 = TL(3, d)
e1, e2 = {T3.e(0): 1}, {T3.e(1): 1}
N = T3.add(T3.mul(e1, e2), T3.mul(e2, e1), -1)
N3 = T3.mul(T3.mul(N, N), N)
ok3 = all(sp.simplify(v) == 0 for v in T3.add(N3, T3.scal(1 - d**2, N), -1).values())
rep("6-ended tangles (TL_3): N = e1e2 - e2e1 obeys N^3 = (1 - delta^2) N, so on its 2-dimensional block N^2 = -(delta^2 - 1): k = 0, d = delta^2 - 1, inside the parabola for delta > 1", ok3)

json.dump(out, open(os.path.join(HERE, "checks.json"), "w"))
print("ALL PASS" if okall else "SOME FAILED")
