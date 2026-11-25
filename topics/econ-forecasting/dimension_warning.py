"""Systemic risk as lost dimension.  For independent risks the portfolio's cube factors:
M_sum - 1 = (M_1 - 1) (+) ... (+) (M_n - 1), and it has n independent axes.  Correlation leans the
cube onto its diagonal.  Effective dimension = participation ratio of the correlation matrix,
d_eff = (sum lambda)^2 / sum lambda^2  (n when independent, 1 when everything moves together).
Test: does low d_eff/n forecast next month's volatility or a drawdown, beyond trailing volatility?"""
import numpy as np, pandas as pd, json, load
out = []
def say(s=""): print(s); out.append(s)

def build(name, min_n, start=None):
    P = load.xts(name)
    if start: P = P.loc[start:]
    R = np.log(P).diff().iloc[1:]
    S = P.pct_change(fill_method=None).iloc[1:]
    ew = S.mean(axis=1, skipna=True)                 # equal-weight portfolio of members with data
    dates = R.index
    W = 250
    rec = []
    for t in range(W, len(dates), 5):                # weekly
        win = R.iloc[t - W:t]
        ok = win.notna().all().values
        if ok.sum() < min_n: continue
        X = win.values[:, ok]
        X = X[:, X.std(axis=0) > 0]
        C = np.corrcoef(X, rowvar=False)
        lam = np.linalg.eigvalsh(C)
        n = C.shape[0]
        pr = lam.sum() ** 2 / (lam ** 2).sum()
        ar = np.sort(lam)[::-1][:max(1, n // 5)].sum() / lam.sum()   # Kritzman absorption ratio (top fifth)
        rec.append((dates[t], n, pr, pr / n, ar))
    D = pd.DataFrame(rec, columns=["date", "n", "d_eff", "d_frac", "absorb"]).set_index("date")
    # realized vol features of the equal-weight portfolio
    lr = np.log1p(ew)
    rv = lambda k: np.sqrt((lr ** 2).rolling(k).mean() * 252)
    feat = pd.DataFrame({"rv1": rv(1), "rv5": rv(5), "rv22": rv(22)})
    fut = np.sqrt((lr ** 2).rolling(21).mean().shift(-21) * 252)
    wealth = np.exp(lr.cumsum())
    fwd_min = wealth[::-1].rolling(63, min_periods=63).min()[::-1].shift(-1)
    dd = (fwd_min / wealth - 1 < -0.10).astype(float).where(fwd_min.notna())
    D = D.join(feat).join(fut.rename("fut_rv")).join(dd.rename("dd10"))
    D["absorb_shift"] = (D.absorb.rolling(3).mean() - D.absorb.rolling(52).mean()) / D.absorb.rolling(52).std()
    return D.dropna()

def ols_oos(D, y, cols, start_frac=0.25):
    """Expanding-window out-of-sample forecasts."""
    X = np.column_stack([np.ones(len(D))] + [D[c].values for c in cols]); Y = D[y].values
    k0 = int(len(D) * start_frac); pred = np.full(len(D), np.nan)
    for t in range(k0, len(D) - 4):
        tr = slice(0, t - 4)    # leave a month gap: targets overlap
        b = np.linalg.lstsq(X[tr], Y[tr], rcond=None)[0]; pred[t] = X[t] @ b
    return pred

def logit_oos(D, y, cols, start_frac=0.25):
    from scipy.optimize import minimize
    X = np.column_stack([np.ones(len(D))] + [(D[c] - D[c].mean()).values / D[c].std() for c in cols]); Y = D[y].values
    k0 = int(len(D) * start_frac); pred = np.full(len(D), np.nan); b = np.zeros(X.shape[1])
    nll = lambda b, X, Y: np.sum(np.logaddexp(0, X @ b) - Y * (X @ b)) + 0.5 * b @ b
    for t in range(k0, len(D) - 13, 13):            # refit each quarter
        tr = slice(0, t - 13)
        b = minimize(nll, b, args=(X[tr], Y[tr]), method="BFGS").x
        pred[t:t + 13] = 1 / (1 + np.exp(-(X[t:t + 13] @ b)))
    return pred

def auc(p, y):
    o = np.argsort(p); r = np.empty(len(p)); r[o] = np.arange(1, len(p) + 1)
    n1 = y.sum(); n0 = len(y) - n1
    return (r[y == 1].sum() - n1 * (n1 + 1) / 2) / (n1 * n0)

results = {}
for name, min_n, start in (("DJ_const", 15, None), ("SP500_const", 150, "1985-01-01")):
    D = build(name, min_n, start)
    for c in ("rv1", "rv5", "rv22", "fut_rv"): D["l" + c] = np.log(D[c].clip(lower=1e-4))
    D["ld"] = np.log(D.d_frac)
    say(f"\n== {name}: {D.index[0].date()} to {D.index[-1].date()}, {len(D)} weeks, stocks per window {D.n.min()}-{D.n.max()} ==")
    say(f" effective dimension: median {D.d_eff.median():.1f} axes ({D.d_frac.median():.0%} of n); lowest {D.d_eff.min():.1f} on {D.d_eff.idxmin().date()}; highest {D.d_eff.max():.1f}")
    say(f" corr(log d_eff/n, log trailing 22-day vol) = {np.corrcoef(D.ld, D.lrv22)[0,1]:.2f}")
    base = ["lrv1", "lrv5", "lrv22"]
    specs = {"HAR (vol only)": base, "HAR + dimension": base + ["ld"], "HAR + absorption shift": base + ["absorb_shift"], "dimension only": ["ld"]}
    m = {}
    Y = D.lfut_rv.values
    preds = {k: ols_oos(D, "lfut_rv", v) for k, v in specs.items()}
    ok = ~np.isnan(preds["HAR (vol only)"])
    bm = np.mean((Y[ok] - np.cumsum(Y)[ok] / np.arange(1, len(Y) + 1)[ok]) ** 2)
    for k, p in preds.items():
        mse = np.mean((Y[ok] - p[ok]) ** 2); m[k] = {"oos_mse": float(mse), "oos_r2": float(1 - mse / bm)}
        say(f"  next-month log vol, {k:24s} out-of-sample MSE {mse:.4f}  R2 {1-mse/bm:.3f}")
    gain = 1 - m["HAR + dimension"]["oos_mse"] / m["HAR (vol only)"]["oos_mse"]
    say(f"  adding dimension to HAR changes MSE by {-100*gain:+.1f}%")
    # in-sample t-stat of dimension with HAC-ish block bootstrap
    X = np.column_stack([np.ones(len(D))] + [D[c].values for c in base + ["ld"]])
    b = np.linalg.lstsq(X, Y, rcond=None)[0]
    rng = np.random.default_rng(1); bs = []
    for _ in range(400):
        idx = np.concatenate([np.arange(s, s + 26) for s in rng.integers(0, len(D) - 26, len(D) // 26)])
        bs.append(np.linalg.lstsq(X[idx], Y[idx], rcond=None)[0][-1])
    say(f"  full-sample coefficient on log(d_eff/n): {b[-1]:+.3f} (block-bootstrap se {np.std(bs):.3f})")
    yd = D.dd10.values
    for k, cols in (("vol only", base), ("vol + dimension", base + ["ld"]), ("vol + absorption shift", base + ["absorb_shift"])):
        p = logit_oos(D, "dd10", cols); ok2 = ~np.isnan(p)
        brier = np.mean((p[ok2] - yd[ok2]) ** 2); a = auc(p[ok2], yd[ok2])
        m["dd_" + k] = {"brier": float(brier), "auc": float(a)}
        say(f"  10% drawdown within 3 months, {k:22s} Brier {brier:.4f}  AUC {a:.3f}  (base rate {yd[ok2].mean():.1%})")
    results[name] = m
    D[["n", "d_eff", "d_frac", "absorb", "rv22", "fut_rv", "dd10"]].to_csv(f"dimension_{name}.csv")
json.dump(results, open("dimension_warning.json", "w"), indent=1)
open("dimension_warning_output.txt", "w").write("\n".join(out) + "\n")
