"""DRAFT. Blow-up hunt analysis helpers. T_trust: last time the grid still resolves the flow (delta*kmax >= 2),
interpolated in log delta. For a real blow-up at time T, T_trust(N) closes in on T (gaps shrink each doubling);
for a smooth flow whose delta shrinks exponentially, T_trust grows by about the same step each doubling."""
import json, glob, numpy as np
L = lambda f: json.load(open(f))
def load_all(d='runs/hunt2'):
    R = {}
    for f in glob.glob(d + '/h_*.json'):
        r = L(f); R[(r['ic'], r['N'], r['nu'], r['theta'])] = r
    return R
def t_trust(r, c=2.0):
    kmax = r['N']/3; rows = [q for q in r['rows'] if np.isfinite(q['delta'])]
    for a, b in zip(rows, rows[1:]):
        if a['delta']*kmax >= c > b['delta']*kmax:
            la, lb, lc = np.log(a['delta']*kmax), np.log(b['delta']*kmax), np.log(c)
            return a['t'] + (b['t'] - a['t'])*(la - lc)/(la - lb)
    return None if rows[-1]['delta']*kmax < c else float('nan')   # nan: still resolved at the end
if __name__ == '__main__':
    R = load_all()
    for k in sorted(R):
        r = R[k]; print(k, 't_trust', None if t_trust(r) is None else round(t_trust(r), 3), 't_end', round(r['t_end'], 3), r['stop'],
                        'max spin', round(max(q['wmax'] for q in r['rows']), 2), 'max grad', round(max(q['gmax'] for q in r['rows']), 2))

def fit_blowup(r, frac=0.4, c=2.0):
    """fit delta(t) on the last `frac` of the resolved stretch (delta*kmax >= c) two ways:
    power law delta = A (T - t)^alpha (a blow-up at T) and exponential delta = A exp(-t/tau) (no blow-up).
    Returns T, alpha, and the rms misfit of each in log delta."""
    kmax = r['N']/3; tt = t_trust(r, c); rows = [q for q in r['rows'] if np.isfinite(q['delta']) and q['delta'] > 0]
    tend = tt if tt is not None and np.isfinite(tt) else rows[-1]['t']
    window = frac*tend; sel = [q for q in rows if tend - window <= q['t'] <= tend]
    t = np.array([q['t'] for q in sel]); ld = np.log([q['delta'] for q in sel])
    if len(t) < 6: return None
    pe = np.polyfit(t, ld, 1); me = float(np.sqrt(np.mean((np.polyval(pe, t) - ld)**2)))
    best = None
    for T in t[-1] + np.geomspace(1e-3, 20, 400):
        x = np.log(T - t); p = np.polyfit(x, ld, 1); m = float(np.sqrt(np.mean((np.polyval(p, x) - ld)**2)))
        if best is None or m < best[2]: best = (float(T), float(p[0]), m)
    return dict(T=best[0], alpha=best[1], misfit_power=best[2], tau=float(-1/pe[0]), misfit_exp=me, t_fit_end=float(t[-1]))

def grid_close(R, ic, theta, nu=0.0, c=2.0, grids=(64, 128, 256)):
    """how far each grid gets (t_trust at threshold c), the ratio of successive steps, and the extrapolated blow-up time.
    ratio < 1: the grids close in on a time T (blow-up); ratio near 1: each finer grid goes the same step further (smooth so far)."""
    t = [t_trust(R[(ic, N, nu, theta)], c) if (ic, N, nu, theta) in R else None for N in grids]
    if any(x is None or not np.isfinite(x) for x in t): return dict(t=t)
    r = (t[2] - t[1])/(t[1] - t[0])
    return dict(t=t, ratio=float(r), T=float(t[2] + (t[2] - t[1])*r/(1 - r)) if r < 1 else float('inf'))
