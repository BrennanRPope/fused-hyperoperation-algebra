"""Long-run discounting: average the discount factors (the cube's corners), not the rates.
Real rate each year = 10-year Treasury yield (-) that year's CPI inflation (Shiller data).
Backtest: at each origin year, using only earlier years, forecast the discount factor that
rolling the real rate forward for T years actually delivered."""
import numpy as np, json, load
rng = np.random.default_rng(7)
out = []
def say(s=""): print(s); out.append(s)
A = load.shiller_annual()
r = ((1 + A.r_nom) / (1 + A.infl) - 1)      # ex-post real rate, mult-add difference
years = r.index.values; rv = r.values
say(f"real rate {years[0]}-{years[-1]}: mean {100*rv.mean():.2f}%, mult-add mean {100*(np.exp(np.log1p(rv).mean())-1):.2f}%, sd {100*rv.std():.2f}%, min {100*rv.min():.1f}%, share of negative years {np.mean(rv<0):.0%}")
ac = np.corrcoef(rv[:-1], rv[1:])[0, 1]; say(f"year-to-year autocorrelation {ac:.2f}")

def boot_paths(hist, T, n, block):
    L = len(hist); paths = np.empty((n, T))
    for i in range(n):
        seq = []
        while len(seq) < T:
            s = rng.integers(0, L - block + 1); seq.extend(hist[s:s + block])
        paths[i] = seq[:T]
    return paths

def forecasts(hist, T, n=4000):
    f = {}
    f["constant mean rate"] = (1 + hist.mean()) ** -T
    f["constant mult-add mean"] = np.exp(-T * np.log1p(hist).mean())
    for b in (1, 10):
        p = boot_paths(hist, T, n, b)
        f[f"averaged factors, block {b}"] = np.mean(np.exp(-np.log1p(p).sum(axis=1)))
    f["averaged, uncertain mean"] = np.mean(np.exp(-np.log1p(uncertain_mean_paths(hist, T, n)).sum(axis=1)))
    return f

def uncertain_mean_paths(hist, T, n, block=10):
    """Block-bootstrap the deviations from the mean, and draw the long-run mean itself from
    its own uncertainty (standard error with overlapping-block correction).  The mean is
    held for the whole path: that persistence is what makes long-horizon rates decline."""
    m = hist.mean(); dev = hist - m
    se = hist.std() / np.sqrt(len(hist) / block)
    mu = m + se * rng.standard_normal(n)
    return mu[:, None] + boot_paths(dev, T, n, block)

H0 = 30  # need 30 years of history before the first origin
res = {}
for T in (10, 20, 30, 40, 50):
    rows = []
    for k in range(H0, len(rv) - T + 1):
        hist = rv[:k]
        realized = np.exp(-np.log1p(rv[k:k + T]).sum())
        rows.append((years[k], realized, forecasts(hist, T)))
    real = np.array([x[1] for x in rows])
    names = list(rows[0][2])
    say(f"\nT = {T} years, {len(rows)} origins {rows[0][0]}-{rows[-1][0]} (overlapping; about {len(rows)//T + 1} independent windows)")
    say(f"  realized discount factor: mean {real.mean():.3f}, median {np.median(real):.3f}")
    res[T] = {}
    for nm in names:
        fc = np.array([x[2][nm] for x in rows])
        ratio = real.mean() / fc.mean()
        mae = np.mean(np.abs(np.log(real / fc)))
        say(f"  {nm:28s} mean forecast {fc.mean():.3f}  realized/forecast {ratio:5.2f}  mean |log error| {mae:.3f}")
        res[T][nm] = dict(mean_forecast=float(fc.mean()), ratio=float(ratio), mae_log=float(mae))
    res[T]["realized_mean"] = float(real.mean())

say("\nCertainty-equivalent rate from the whole record (uncertain mean + block 10): r_CE(T) = DF(T)^(-1/T) - 1")
sched = {}
for T in (1, 10, 25, 50, 75, 100, 150, 200, 300, 400):
    p = uncertain_mean_paths(rv, T, 20000)
    df = np.mean(np.exp(-np.log1p(p).sum(axis=1)))
    sched[T] = float(df ** (-1 / T) - 1)
    say(f"  T={T:3d}: averaged factor {df:.3e}, r_CE {100*sched[T]:.2f}%   (constant mean rate gives {100*rv.mean():.2f}%)")
harm = 1e6
for T in (50, 100, 200):
    pc = harm * (1 + rv.mean()) ** -T; pa = harm * (1 + sched[T]) ** -T
    say(f"  a $1M harm in {T} years: worth ${pc:,.0f} today at a constant mean rate, ${pa:,.0f} with averaged factors ({pa/pc:.1f}x)")
say("\nSame with 10-year blocks only (no uncertainty about the long-run mean):")
sched_b = {}
for T in (1, 10, 50, 100, 200):
    p = boot_paths(rv, T, 20000, 10)
    df = np.mean(np.exp(-np.log1p(p).sum(axis=1))); sched_b[T] = float(df ** (-1 / T) - 1)
    say(f"  T={T:3d}: r_CE {100*sched_b[T]:.2f}%   a $1M harm then is worth ${harm*df:,.0f} today")
json.dump({"backtest": res, "schedule": sched, "schedule_blocks": sched_b}, open("discount_backtest.json", "w"), indent=1)
open("discount_backtest_output.txt", "w").write("\n".join(out) + "\n")
