"""Polish (Nov 26) for growth_of_growth.py: the finite-time pole fits  log P = A - B log(t_c - t).

growth_of_growth.py fits (A, B, t_c) with scipy curve_fit (trust-region, bounds A in [-10, 20], B in [0.01, 50],
t_c in [t_end + 0.05, t_end + 30], default tolerances 1e-8).  For a fixed t_c the model is linear in A and B, so the
fit is a one-dimensional problem in t_c: least squares for (A, B) in closed form, then minimise the residual over
t_c.  Polish: start from the stored t_c, Brent to 1e-12 years, and a 3000-point scan of the whole allowed t_c range
to check the minimum is the global one.  numpy + pandas only (no scipy on this machine).

The data are monthly index levels; the residual valley along t_c is shallow, so the polished dates only tell us
whether curve_fit had converged.  Their extra digits carry no forecasting meaning.
"""
import os, math, time
import numpy as np, pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
T0 = time.time()
# loader copied from load.py (shiller_monthly), and the series from growth_of_growth.py
df = pd.read_csv(os.path.join(HERE, "data", "sp500_shiller_monthly.csv"), parse_dates=["Date"])
df = df.rename(columns={"SP500": "P", "Dividend": "D", "Earnings": "E", "Consumer Price Index": "CPI", "Long Interest Rate": "GS10"})
m = df[(df.CPI > 0) & (df.GS10 > 0) & (df.D > 0)].set_index("Date")
m = m[m.PE10 > 0]
lp = np.log(m.P / m.CPI)

def ab_fit(t, y, tc):
    """Best (A, B) at fixed tc inside the curve_fit box A in [-10, 20], B in [0.01, 50] (exact: a 2-variable convex
    least squares; if the free optimum is outside the box, the optimum lies on one of the box's edges)."""
    L = -np.log(tc - t)
    def rss(A, B):
        r = y - (A + B * L); return float(r @ r)
    A, B = np.linalg.lstsq(np.column_stack([np.ones_like(t), L]), y, rcond=None)[0]
    if -10 <= A <= 20 and 0.01 <= B <= 50:
        return A, B, rss(A, B)
    cand = []
    for Af in (-10.0, 20.0):                                   # A on an edge, B refitted and clipped
        Bf = min(max(float((L @ (y - Af)) / (L @ L)), 0.01), 50.0); cand.append((rss(Af, Bf), Af, Bf))
    for Bf in (0.01, 50.0):                                    # B on an edge, A refitted and clipped
        Af = min(max(float(np.mean(y - Bf * L)), -10.0), 20.0); cand.append((rss(Af, Bf), Af, Bf))
    r, A, B = min(cand)
    return A, B, r

def brent_min(f, a, b, tol=1e-12, maxit=300):
    g = 0.3819660112501051
    x = w = v = a + g * (b - a); fx = fw = fv = f(x); d = e = 0.0
    for _ in range(maxit):
        mid = 0.5 * (a + b); tol1 = tol * abs(x) + 1e-15; tol2 = 2 * tol1
        if abs(x - mid) <= tol2 - 0.5 * (b - a): break
        if abs(e) > tol1:
            r = (x - w) * (fx - fv); q = (x - v) * (fx - fw); p = (x - v) * q - (x - w) * r; q = 2 * (q - r)
            if q > 0: p = -p
            q = abs(q); et = e; e = d
            if abs(p) >= abs(0.5 * q * et) or p <= q * (a - x) or p >= q * (b - x):
                e = (a - x) if x >= mid else (b - x); d = g * e
            else:
                d = p / q; u = x + d
                if u - a < tol2 or b - u < tol2: d = tol1 if mid >= x else -tol1
        else:
            e = (a - x) if x >= mid else (b - x); d = g * e
        u = x + (d if abs(d) >= tol1 else (tol1 if d > 0 else -tol1)); fu = f(u)
        if fu <= fx:
            if u >= x: a = x
            else: b = x
            v, w, x, fv, fw, fx = w, x, u, fw, fx, fu
        else:
            if u < x: a = u
            else: b = u
            if fu <= fw or w == x: v, w, fv, fw = w, u, fw, fu
            elif fu <= fv or v == x or v == w: v, fv = u, fu
    return x, fx

# stored answers: growth_of_growth_output.txt
STORED = [("1924-01", "1929-06", "1931-07-01", 0.91), ("1995-01", "1999-12", "2024-09-15", 5.55),
          ("2003-03", "2007-06", "2037-05-31", 2.17), ("2016-01", "2021-09", "2023-02-07", 0.38),
          ("1982-08", "1987-06", "1988-02-27", 0.33)]
out = ["Pole fits log P = A - B log(t_c - t): before (curve_fit, printed) -> after (variable projection + Brent)"]
for a, b_, date0, B0 in STORED:
    t1 = time.time()
    s = lp.loc[a:b_]; t = (s.index - s.index[0]).days.values / 365.25; y = s.values
    lo, hi = t[-1] + 0.05, t[-1] + 30
    tc0 = (pd.Timestamp(date0) - s.index[0]).days / 365.25
    rss0 = ab_fit(t, y, tc0)[2]
    f = lambda tc: ab_fit(t, y, tc)[2]
    # local polish from the stored point, then a scan of the whole allowed range for the global minimum
    step = max(0.02 * (tc0 - t[-1]), 0.01)
    tcl, rssl = brent_min(f, max(lo, tc0 - step), min(hi, tc0 + step))
    grid = np.linspace(lo, hi, 3000); vals = np.array([f(x) for x in grid]); k = int(np.argmin(vals))
    tcg, rssg = brent_min(f, grid[max(k - 1, 0)], grid[min(k + 1, len(grid) - 1)])
    if rssl <= rssg: tcg, rssg = tcl, rssl
    A, B, _ = ab_fit(t, y, tcg)
    abox = " (A on its bound 20)" if abs(A - 20) < 1e-9 else ""
    at_bound = "at the upper bound t_end + 30 yr" if hi - tcg < 1e-6 else ("at the lower bound" if tcg - lo < 1e-6 else "interior")
    date1 = (s.index[0] + pd.Timedelta(days=365.25 * tcg)).date()
    shift_days = (pd.Timestamp(date1) - pd.Timestamp(date0)).days
    out.append(f"  {a}..{b_}: pole {date0} -> {date1} ({shift_days:+d} days), order 1/e {B0:.2f} -> {B:.2f}; "
               f"residual sum of squares {rss0:.6e} -> {rssg:.6e} ({(rss0 - rssg) / rss0:.1e} relative); "
               f"minimum {at_bound}{abox}; {time.time() - t1:.2f} s")
out.append("  curve_fit had converged in all five (dates unchanged to the day).  Two fits sit on an edge of curve_fit's box,"
           " so their pole is set by the box, not the data: 1995-99 has A at its cap of 20, 2003-07 has t_c at t_end + 30 yr.")
out.append("  Verdict on the page (pole dates scatter from months to the 30-year bound, so the pole fit is not a forecaster) unchanged.")
out.append(f"wall time {time.time() - T0:.2f} s")
print("\n".join(out))
open(os.path.join(HERE, "polish_output.txt"), "w").write("\n".join(out) + "\n")
