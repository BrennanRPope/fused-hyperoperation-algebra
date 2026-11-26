"""Polish step for the Richards growth-curve fits in checks.py (section 2).

checks.py fits K, r, t0, nu with scipy's least_squares (default tolerances 1e-8, several starts) and
then refits with nu fixed at 1 (logistic), 0 (Gompertz) and -1/3 (von Bertalanffy) for the F-test
p-values. The polish starts each free fit from the parameters stored in fits.json and each fixed fit
from the same start checks.py uses, and runs Levenberg-Marquardt to xtol 1e-14 (numpy only; scipy is
not installed here, so the F-test tail is computed from the incomplete beta function).
These are real, noisy measurements: digits beyond the second decimal of nu carry no information;
the polish only shows whether the optimizer had converged.
Run: python3 polish.py > polish_output.txt
"""
import json, math, time
import numpy as np, pandas as pd

def lm(fun, x0, lo=None, hi=None, xtol=1e-14, ftol=1e-15, maxit=500):
    """Small Levenberg-Marquardt (numpy only; scipy is not installed here). Central-difference Jacobian,
    box bounds by clipping. Returns (x, residuals, iterations)."""
    x = np.array(x0, float); n = len(x)
    lo = np.full(n, -np.inf) if lo is None else np.array(lo, float)
    hi = np.full(n, np.inf) if hi is None else np.array(hi, float)
    x = np.clip(x, lo, hi); r = fun(x); cost = r @ r; lam = 1e-6
    for it in range(1, maxit + 1):
        J = np.empty((len(r), n))
        for i in range(n):
            h = 1e-6 * max(1.0, abs(x[i])); e = np.zeros(n); e[i] = h
            J[:, i] = (fun(x + e) - fun(x - e)) / (2 * h)
        g = J.T @ r; A = J.T @ J
        free = ~(((x <= lo) & (g > 0)) | ((x >= hi) & (g < 0)))   # coordinates not pinned at a bound
        improved = False
        for _ in range(30):
            dx = np.zeros(n)
            Af = A[np.ix_(free, free)]
            try:
                dx[free] = np.linalg.solve(Af + lam * np.diag(np.diag(Af) + 1e-300), -g[free])
            except np.linalg.LinAlgError:
                lam *= 10; continue
            xn = np.clip(x + dx, lo, hi); rn = fun(xn); cn = rn @ rn
            if np.isfinite(cn) and cn <= cost:
                improved = True; break
            lam *= 10
        if not improved:
            break
        step = np.max(np.abs(xn - x) / (np.abs(x) + 1e-12)); drop = (cost - cn) / max(cost, 1e-300)
        x, r, cost = xn, rn, cn; lam = max(lam / 10, 1e-12)
        if step < xtol or drop < ftol:
            break
    return x, r, it
def _betacf(a, b, x):
    qab, qap, qam = a + b, a + 1, a - 1; c = 1.0; d = 1 - qab * x / qap
    d = 1 / (d if abs(d) > 1e-300 else 1e-300); h = d
    for m in range(1, 400):
        m2 = 2 * m; aa = m * (b - m) * x / ((qam + m2) * (a + m2))
        d = 1 + aa * d; d = 1 / (d if abs(d) > 1e-300 else 1e-300); c = 1 + aa / c if abs(1 + aa / c) > 1e-300 else 1e-300; h *= d * c
        aa = -(a + m) * (qab + m) * x / ((a + m2) * (qap + m2))
        d = 1 + aa * d; d = 1 / (d if abs(d) > 1e-300 else 1e-300); c = 1 + aa / c if abs(1 + aa / c) > 1e-300 else 1e-300
        de = d * c; h *= de
        if abs(de - 1) < 3e-16: break
    return h
def betainc(a, b, x):
    """Regularized incomplete beta I_x(a, b) (Numerical Recipes continued fraction)."""
    if x <= 0: return 0.0
    if x >= 1: return 1.0
    lbt = math.lgamma(a + b) - math.lgamma(a) - math.lgamma(b) + a * math.log(x) + b * math.log1p(-x)
    if x < (a + 1) / (a + b + 2): return math.exp(lbt) * _betacf(a, b, x) / a
    return 1 - math.exp(lbt) * _betacf(b, a, 1 - x) / b
def f_sf(F, d1, d2):
    """Upper tail of the F distribution, = 1 - scipy.stats.f.cdf(F, d1, d2)."""
    if F <= 0: return 1.0
    return betainc(d2 / 2, d1 / 2, d2 / (d2 + d1 * F))

def richards(p, tt):
    K, r, t0, nu = p
    base = 1 + nu*np.exp(-r*(tt-t0))
    if abs(nu) < 1e-8: return K*np.exp(-np.exp(-r*(tt-t0)))
    return K*np.maximum(base, 1e-12)**(-1/nu)

stored = json.load(open("fits.json"))
sets = {'Soybean': ('Time', 'weight'), 'Loblolly': ('age', 'height'), 'Orange': ('age', 'circumference')}
T0 = time.perf_counter(); verdict = {}
for nm, (tc, yc) in sets.items():
    df = pd.read_csv(f'data/{nm}.csv'); tt = df[tc].values.astype(float); yy = df[yc].values.astype(float); n = len(yy)
    t = time.perf_counter()
    p0 = np.array(stored[nm]['p'])
    lo = [yy.max()*0.5, 1e-4, tt.min()-5*(tt.max()-tt.min()), -0.95]; hi = [yy.max()*20, 10, tt.max()*3, 20]
    res = lambda p: richards(p, tt) - yy
    rss0 = float(res(p0) @ res(p0))
    p1, r1, it = lm(res, p0, lo, hi); rss1 = float(r1 @ r1)
    print(f"{nm} (n = {n})")
    print(f"  free fit  before: K={p0[0]:.6g} r={p0[1]:.6g} t0={p0[2]:.6g} nu={p0[3]:.10f}  RSS={rss0:.10g}")
    print(f"            after : K={p1[0]:.6g} r={p1[1]:.6g} t0={p1[2]:.6g} nu={p1[3]:.10f}  RSS={rss1:.10g}  ({it} LM steps)")
    print(f"            nu change {p1[3]-p0[3]:+.2e}; RSS change {rss1-rss0:+.2e}; nu at 2 decimals {p0[3]:.2f} -> {p1[3]:.2f}"
          + ("  (no change)" if f"{p0[3]:.2f}" == f"{p1[3]:.2f}" else "  CHANGED"))
    K0 = yy.max()*1.1; t00 = np.median(tt); r0 = 4/(tt.max()-tt.min())
    for nu, lab in ((1, 'logistic'), (1e-9, 'Gompertz'), (-1/3, 'von Bertalanffy')):
        f = lambda q: richards(list(q)+[nu], tt) - yy
        q, rq, _ = lm(f, [K0, r0, t00], lo[:3], hi[:3]); rssq = float(rq @ rq)
        F = (rssq - rss1)/(rss1/(n-4)); pv = f_sf(F, 1, n-4); pv0 = stored[nm]['pv'][lab]
        verdict[(nm, lab)] = (pv0 > 0.05, pv > 0.05)
        print(f"  p({lab}): before {pv0:.3g}, after {pv:.3g}" + ("" if (pv0 > 0.05) == (pv > 0.05) else "   VERDICT CHANGED"))
    print(f"  time {time.perf_counter()-t:.2f} s")
print(f"total wall time {time.perf_counter()-T0:.2f} s")
print("verdicts (p > 0.05 or not) unchanged" if all(a == b for a, b in verdict.values()) else "SOME VERDICT CHANGED")
