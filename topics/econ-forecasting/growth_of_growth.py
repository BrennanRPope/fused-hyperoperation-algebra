"""The next rung: growth whose rate grows (reflexive feedback).  Exponential growth (rung 3)
has a constant rate and no pole.  If the rate rises with the level, dP/dt = k P^(1+e), the path
is hyperbolic and hits a pole at a finite time t_c, P ~ (t_c - t)^(-1/e): the order of the pole
is a dimension again.  Locally this shows up as positive curvature (acceleration) of log price.
Test: does acceleration of the S&P's real log price forecast the next 3 years' real return,
beyond the standard valuation baseline (Shiller CAPE)?"""
import numpy as np, pandas as pd, json, load
out = []
def say(s=""): print(s); out.append(s)
m = load.shiller_monthly()
m = m[m.PE10 > 0]
real = m.P / m.CPI
tr = ((m.P + m.D / 12) / m.P.shift(1)).fillna(1).cumprod() / m.CPI    # real total return index
lp = np.log(real)
H = 36
fut = np.log(tr.shift(-H) / tr)                                          # next 3 years, log real total return
W = 60
tt = np.arange(W) / 12.0
def accel(x):
    c = np.polyfit(tt, x, 2); return 2 * c[0]                           # second derivative of log price, per year^2
acc = lp.rolling(W).apply(accel, raw=True)
trend = lp.rolling(W).apply(lambda x: np.polyfit(tt, x, 1)[0], raw=True)
D = pd.DataFrame({"fut": fut, "lcape": np.log(m.PE10), "acc": acc, "slope": trend}).dropna()
say(f"{D.index[0].date()} to {D.index[-1].date()}, {len(D)} months; target = next {H} months' log real total return")
for y in ["1929-08", "1937-02", "1966-01", "1987-08", "1999-12", "2007-10", "2021-12"]:
    if y in D.index.strftime("%Y-%m"):
        r = D.loc[y].iloc[0]; say(f"  {y}: acceleration {r.acc:+.3f}/yr^2, 5-yr growth rate {r.slope:+.3f}/yr, CAPE {np.exp(r.lcape):.1f}, next 3 yrs {r.fut:+.2f}")
say(f"  share of months with acceleration above +0.05/yr^2: {np.mean(D.acc>0.05):.1%}")

def oos(cols, start="1921-01-01"):
    Y = D.fut.values; X = np.column_stack([np.ones(len(D))] + [D[c].values for c in cols])
    k0 = np.searchsorted(D.index, pd.Timestamp(start)); p = np.full(len(D), np.nan); hm = np.full(len(D), np.nan)
    for t in range(k0, len(D)):
        tr_ = slice(0, t - H)                     # only targets already known at t
        p[t] = X[t] @ np.linalg.lstsq(X[tr_], Y[tr_], rcond=None)[0]; hm[t] = Y[tr_].mean()
    ok = ~np.isnan(p)
    return 1 - np.mean((Y[ok] - p[ok]) ** 2) / np.mean((Y[ok] - hm[ok]) ** 2)
res = {}
for nm, cols in (("CAPE (baseline)", ["lcape"]), ("acceleration only", ["acc"]), ("CAPE + acceleration", ["lcape", "acc"]),
                 ("CAPE + 5-yr growth rate", ["lcape", "slope"]), ("CAPE + rate + acceleration", ["lcape", "slope", "acc"])):
    r2 = oos(cols); res[nm] = r2
    say(f"  {nm:28s} out-of-sample R2 vs historical mean: {r2:+.3f}")
# full-sample fit with block bootstrap se
X = np.column_stack([np.ones(len(D)), D.lcape, D.acc]); Y = D.fut.values
b = np.linalg.lstsq(X, Y, rcond=None)[0]
rng = np.random.default_rng(3); bs = []
for _ in range(500):
    idx = np.concatenate([np.arange(s, s + 60) for s in rng.integers(0, len(D) - 60, len(D) // 60)])
    bs.append(np.linalg.lstsq(X[idx], Y[idx], rcond=None)[0])
se = np.std(bs, axis=0)
say(f"  full sample: coef on log CAPE {b[1]:+.3f} (se {se[1]:.3f}), on acceleration {b[2]:+.3f} (se {se[2]:.3f})")

# finite-time pole fit on the 1990s run-up: log P = A - (1/e) log(t_c - t)
from scipy.optimize import curve_fit
def pole_fit(a, b_):
    s = lp.loc[a:b_]; t = (s.index - s.index[0]).days.values / 365.25
    f = lambda t, A, B, tc: A - B * np.log(tc - t)
    try:
        p, _ = curve_fit(f, t, s.values, p0=(s.values[0], 1.0, t[-1] + 2), bounds=([-10, 0.01, t[-1] + 0.05], [20, 50, t[-1] + 30]), maxfev=20000)
        return s.index[0] + pd.Timedelta(days=365.25 * p[2]), p[1]
    except Exception as e:
        return None, None
for a, b_, peak in (("1924-01", "1929-06", "1929-09"), ("1995-01", "1999-12", "2000-08"), ("2003-03", "2007-06", "2007-10"), ("2016-01", "2021-09", "2022-01"), ("1982-08", "1987-06", "1987-08")):
    tc, B = pole_fit(a, b_)
    say(f"  pole fit on {a}..{b_}: predicted pole {tc.date() if tc is not None else 'none'} (pole order 1/e = {B:.2f}); actual peak {peak}" if tc is not None else f"  pole fit {a}..{b_}: failed")
json.dump(res, open("growth_of_growth.json", "w"), indent=1)
D.to_csv("growth_of_growth.csv")
open("growth_of_growth_output.txt", "w").write("\n".join(out) + "\n")
