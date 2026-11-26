"""Polish step for the mixture fits in checks.py / fit_mixtures.py (section F) and robustness.py.

What checks.py searched, and how:
  * alpha (Chou-Talalay dial) and the odds dial k: least_squares, default tolerances (1e-8);
  * their 95% ranges: the smallest and largest point of a 0.05-wide grid whose profile F statistic is
    under the 5% critical value. Those ends are limited by the grid spacing, not by the data.
The polish (numpy only; scipy is not installed here):
  1. rebuilds the stored fit: alpha and the two Hill slopes from fits.json held fixed, the other three
     numbers fitted; then frees all six and runs Levenberg-Marquardt to xtol 1e-14;
  2. recomputes p(Loewe), p(HSA), p(Bliss) from the polished fits;
  3. moves each end of each 95% range off the grid: regula falsi on F(alpha) - Fcrit inside the grid
     cell that contains the crossing, to 1e-3 in alpha or 1e-3 in F (a few profile fits per end);
  4. refits robustness.py's version with a lower plateau and reports alpha and its p-values.
The root solve inside the Chou-Talalay model is the same grid-bracket-then-solve as fit_mixtures.py,
vectorised over data points, with bisection to 1e-15 in place of brentq(xtol=1e-14).
These are real, noisy measurements: digits beyond the second decimal of alpha carry no information.
Run: python3 polish.py > polish_output.txt
"""
import json, math, time
import numpy as np, pandas as pd
np.seterr(all="ignore")

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

BASE = np.linspace(60, -60, 2401)
def ct_odds_vec(D1, D2, e1, b1, e2, b2, a):
    """Vectorised copy of fit_mixtures.ct_odds: effect odds o solving v1 + v2 + a v1 v2 = 1."""
    o = np.zeros(len(D1))
    only1 = (D1 > 0) & (D2 <= 0); only2 = (D2 > 0) & (D1 <= 0); both = (D1 > 0) & (D2 > 0)
    o[only1] = (D1[only1]/e1)**b1; o[only2] = (D2[only2]/e2)**b2
    if not both.any(): return o
    d1, d2 = D1[both], D2[both]
    pts = [b1*np.log(d1/e1), b2*np.log(d2/e2)]
    if a != 0 and -1/a > 0: pts += [b1*np.log(d1/(e1*(-1/a))), b2*np.log(d2/(e2*(-1/a)))]
    G = np.sort(np.concatenate([np.broadcast_to(BASE, (len(d1), len(BASE))), np.clip(np.array(pts).T, -60, 60)], axis=1), axis=1)[:, ::-1]
    def g(lo, d1=d1, d2=d2):
        v1 = (d1/e1)[:, None]*np.exp(-lo/b1) if lo.ndim == 2 else d1/e1*np.exp(-lo/b1)
        v2 = (d2/e2)[:, None]*np.exp(-lo/b2) if lo.ndim == 2 else d2/e2*np.exp(-lo/b2)
        return v1 + v2 + a*v1*v2 - 1
    gv = g(G)
    cond = (gv[:, :-1] < 0) & (0 <= gv[:, 1:])
    has = cond.any(axis=1); i = np.argmax(cond, axis=1); r = np.arange(len(d1))
    lo_ = G[r, i+1]; hi_ = G[r, i]                     # g(lo_) >= 0 > g(hi_)
    out = np.full(len(d1), math.exp(-60.0))
    exact = has & (gv[r, i+1] == 0)
    for _ in range(64):                                # bisection, bracket 0.05 wide -> below 1e-15
        mid = 0.5*(lo_ + hi_); gm = g(mid, d1, d2)
        pos = gm >= 0; lo_ = np.where(pos, mid, lo_); hi_ = np.where(pos, hi_, mid)
    out[has] = np.exp(0.5*(lo_ + hi_))[has]; out[exact] = np.exp(G[r, i+1])[exact]
    o[both] = out
    return o

def pred_ct(p, D1, D2):
    d, b1, e1, b2, e2, a = p
    return d/(1 + ct_odds_vec(D1, D2, e1, b1, e2, b2, a))
def pred_odds(p, D1, D2):
    d, b1, e1, b2, e2, k = p
    o1 = np.where(D1 > 0, (np.maximum(D1, 1e-300)/e1)**b1, 0.0); o2 = np.where(D2 > 0, (np.maximum(D2, 1e-300)/e2)**b2, 0.0)
    o = o1 + o2 + k*o1*o2
    return d/(1 + np.maximum(o, -0.999))

def setup(nm):
    df = pd.read_csv(f'data/{nm}.csv'); df.loc[df.pct == 999, 'pct'] = 50
    D = df.dose.values.astype(float); pct = df.pct.values.astype(float); y = df.rgr.values
    e1 = np.median(df[df.pct == 100].dose); e2 = np.median(df[df.pct == 0].dose)
    return D*pct/100, D*(100-pct)/100, y, e1, e2, y[D == 0].mean()

def crit_F(dfd):
    lo, hi = 0.0, 100.0
    for _ in range(200):
        mid = 0.5*(lo + hi)
        if f_sf(mid, 1, dfd) > 0.05: lo = mid
        else: hi = mid
    return 0.5*(lo + hi)

stored = json.load(open("fits.json"))
T0 = time.perf_counter(); report = {}
for nm in ("glymet", "acidiq", "mecter"):
    t = time.perf_counter()
    D1, D2, y, e1, e2, d0 = setup(nm); n = len(y); st = stored[nm]
    out = {}
    for model, pr, kname, cikey, lo5 in (("ct", pred_ct, "alpha", "ci", -3.0), ("odds", pred_odds, "k", "cik", -0.9)):
        lo = [0, 0.2, e1/1e3, 0.2, e2/1e3, lo5]; hi = [10, 20, e1*1e3, 20, e2*1e3, 50]
        res = lambda q: pr(q, D1, D2) - y
        if model == "ct":
            # 1. rebuild the stored answer: alpha and slopes fixed at fits.json values
            b1s, b2s = st["slopes"]; a_s = st["alpha"]
            r3 = lambda q: pr([q[0], b1s, q[1], b2s, q[2], a_s], D1, D2) - y
            q3, _, _ = lm(r3, [d0, e1, e2], [lo[0], lo[2], lo[4]], [hi[0], hi[2], hi[4]], xtol=1e-12)
            pstart = np.array([q3[0], b1s, q3[1], b2s, q3[2], a_s])
        else:
            pstart = np.array([d0, 2, e1, 2, e2, 0.0])     # odds fit: only k was stored, so start where checks.py starts
        rss0 = float(res(pstart) @ res(pstart))
        p, r, it = lm(res, pstart, lo, hi); rss = float(r @ r)
        def rss_fixed(val, start):
            f = lambda q: pr(list(q) + [val], D1, D2) - y
            q, rr, _ = lm(f, start[:5], lo[:5], hi[:5], xtol=1e-9)
            return float(rr @ rr), np.r_[q, val]
        Fstat = lambda r_: (r_ - rss)/(rss/(n-6))
        pvals = {}
        for lab, val in ((("Loewe", 0.0), ("HSA", -1.0)) if model == "ct" else (("Bliss", 1.0),)):
            rr, _ = rss_fixed(val, p); pvals[lab] = f_sf(Fstat(rr), 1, n-6)
        crit = crit_F(n-6)
        ends = []
        for side, edge in ((-1, st[cikey][0]), (+1, st[cikey][1])):
            a_in, a_out = edge, edge + side*0.05          # grid point inside the range, its neighbour outside
            r_in, q_in = rss_fixed(a_in, p); r_out, q_out = rss_fixed(a_out, q_in)
            f_in, f_out = Fstat(r_in) - crit, Fstat(r_out) - crit
            note = "" if (f_in <= 0 < f_out) else f" (grid cell check failed: F-crit {f_in:+.3f} inside, {f_out:+.3f} outside)"
            # regula falsi (Illinois) on F(alpha) - crit between a_in and a_out
            A_, FA, B_, FB, side_kept, qw = a_in, f_in, a_out, f_out, 0, q_in
            for _ in range(12):
                c = B_ - FB*(B_ - A_)/(FB - FA)
                rc, qw = rss_fixed(c, qw); fc = Fstat(rc) - crit
                if fc*FB > 0:
                    B_, FB = c, fc
                    if side_kept == +1: FA /= 2
                    side_kept = +1
                else:
                    A_, FA = c, fc
                    if side_kept == -1: FB /= 2
                    side_kept = -1
                if abs(B_ - A_) < 1e-3 or abs(fc) < 1e-3: break
            ends.append((c, note))
        out[model] = dict(before=st[kname], after=float(p[5]), rss0=rss0, rss=rss, it=it, p=pvals,
                          ci_before=st[cikey], ci_after=[ends[0][0], ends[1][0]], notes=[e[1] for e in ends], slopes=(p[1], p[3]))
    out["time"] = time.perf_counter() - t
    report[nm] = out
    o = out["ct"]; k = out["odds"]
    print(f"{nm} (n = {n})")
    print(f"  alpha  before {o['before']:.6f}  after {o['after']:.6f}  ({o['it']} LM steps; RSS {o['rss0']:.10g} -> {o['rss']:.10g})")
    print(f"         95% range  before {o['ci_before'][0]:.2f} to {o['ci_before'][1]:.2f} (grid)   after {o['ci_after'][0]:.3f} to {o['ci_after'][1]:.3f}{''.join(o['notes'])}")
    print(f"         p(Loewe) before {st['p_loewe']:.3g} after {o['p']['Loewe']:.3g};  p(HSA) before {st['p_hsa']:.3g} after {o['p']['HSA']:.3g}")
    print(f"  odds k before {k['before']:.6f}  after {k['after']:.6f}  ({k['it']} LM steps)")
    print(f"         95% range  before {k['ci_before'][0]:.2f} to {k['ci_before'][1]:.2f} (grid)   after {k['ci_after'][0]:.3f} to {k['ci_after'][1]:.3f}{''.join(k['notes'])}")
    print(f"         p(Bliss) before {st['p_bliss']:.3g} after {k['p']['Bliss']:.3g}")
    print(f"  time {out['time']:.1f} s")
print(f"main fits: wall time {time.perf_counter()-T0:.1f} s")

# ---- robustness.py: the same fit with a lower plateau c (seven numbers) ----
t = time.perf_counter()
printed = {"glymet": (0.025, -0.33, ["8e-06", "0.52", "0.34"]), "acidiq": (0.031, -0.95, ["0.01", "8e-09", "4.4e-16"]),
           "mecter": (-0.029, -0.98, ["0.46", "0.0011", "7.6e-06"])}   # robustness_output.txt
for nm in ("glymet", "acidiq", "mecter"):
    D1, D2, y, e1, e2, d0 = setup(nm); n = len(y)
    def pred7(q):
        d, c, b1, ee1, b2, ee2, a = q
        return c + (d - c)/(1 + ct_odds_vec(D1, D2, ee1, b1, ee2, b2, a))
    lo = [0, -0.5, 0.2, e1/1e3, 0.2, e2/1e3, -3]; hi = [10, 1, 20, e1*1e3, 20, e2*1e3, 50]
    pm = report[nm]["ct"]
    c0, a0, pp0 = printed[nm]
    # start: printed plateau and alpha, the other numbers from the polished main fit (refitted with those two held)
    f5 = lambda q: pred7([q[0], c0, q[1], q[2], q[3], q[4], a0]) - y
    st = report[nm]
    q5, _, _ = lm(f5, [d0, st["ct"]["slopes"][0], e1, st["ct"]["slopes"][1], e2], [lo[0], lo[2], lo[3], lo[4], lo[5]], [hi[0], hi[2], hi[3], hi[4], hi[5]], xtol=1e-8)
    start = np.array([q5[0], c0, q5[1], q5[2], q5[3], q5[4], a0])
    p7, r7, it7 = lm(lambda q: pred7(q) - y, start, lo, hi); rss7 = float(r7 @ r7)
    pv = []
    for a in (-1, -0.5, 0):
        q6, rr, _ = lm(lambda q: pred7(list(q) + [a]) - y, p7[:6], lo[:6], hi[:6], xtol=1e-8)
        F = (float(rr @ rr) - rss7)/(rss7/(n-7)); pv.append(f_sf(F, 1, n-7))
    print(f"{nm} with plateau: c {c0:.3f} -> {p7[1]:.4f}, alpha {a0:.2f} -> {p7[6]:.4f} (main fit {pm['after']:.4f}; shift {p7[6]-pm['after']:+.3f}), "
          f"p(HSA, -0.5, Loewe) {pp0} -> {[f'{x:.2g}' for x in pv]}")
print(f"robustness refits: {time.perf_counter()-t:.1f} s")
print(f"total wall time {time.perf_counter()-T0:.1f} s")
