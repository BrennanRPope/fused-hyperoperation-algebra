"""Polish for stretching_signature.py (Nov 26).
The original fits scan a grid (t* over 400 log-spaced points, the offset c over 300 points) and keep the best.
Here each scan is refined from its best grid point with a golden-section minimiser (relative tolerance 1e-14) on the SAME objective,
using the stored arrival data in stretching_signature.json (no simulations are re-run)."""
import json, math, time
import numpy as np

def slog(x):                     # copied from stretching_signature.py (importing it pulls in scipy, not installed here)
    n=0
    while x>1: x=math.log(x); n+=1
    return n+(x-1)

def golden(obj, lo, hi, tol=1e-14, maxit=200):
    """golden-section minimiser on [lo, hi] (numpy only; scipy is not installed in this environment)"""
    g = (math.sqrt(5) - 1) / 2; a, b = lo, hi; c = b - g*(b - a); d = a + g*(b - a); fc, fd = obj(c), obj(d)
    for _ in range(maxit):
        if abs(b - a) <= tol*max(1.0, abs(a)): break
        if fc < fd: b, d, fd = d, c, fc; c = b - g*(b - a); fc = obj(c)
        else: a, c, fc = c, d, fd; d = a + g*(b - a); fd = obj(d)
    x = c if fc < fd else d
    return x, min(fc, fd)

t0 = time.process_time()
res = json.load(open('stretching_signature.json'))
LAM = {'sabra': 2.0, 'helical': (1 + 5**.5) / 2}

def refine(obj, grid, ok=lambda v: True):
    """best grid point, then golden-section search between its neighbours"""
    vals = np.array([obj(g) for g in grid])
    vals = np.where([ok(g) for g in grid], vals, np.inf)
    i = int(np.argmin(vals)); lo, hi = grid[max(i - 1, 0)], grid[min(i + 1, len(grid) - 1)]
    x, fx = golden(obj, lo, hi)
    if fx < vals[i] and ok(x):
        return vals[i], grid[i], fx, x
    return vals[i], grid[i], vals[i], grid[i]

rows = []
for f in res:
    n = np.array(f['n'], float); T = np.array(f['T'], float); lam = LAM[f['kind']]; span = T[-1] - T[0]
    def geo(ts):
        y = np.log(ts - T); A = np.polyfit(n, y, 1); pred = ts - np.exp(np.polyval(A, n))
        return np.sqrt(np.mean((pred - T)**2))
    def tet(c):
        x = np.array([slog(v*2*math.log(lam) + c) for v in n]); A = np.polyfit(x, T, 1)
        return np.sqrt(np.mean((np.polyval(A, x) - T)**2)) if A[0] > 0 else np.inf
    def r3(c):
        x = np.log(n + c); A = np.polyfit(x, T, 1)
        return np.sqrt(np.mean((np.polyval(A, x) - T)**2)) if A[0] > 0 else np.inf
    g0, ts0, g1, ts1 = refine(geo, T[-1] + np.logspace(-7, 0, 400))
    A1 = np.polyfit(n, np.log(ts1 - T), 1); ratio1 = math.exp(A1[0])
    t_0, _, t_1, _ = refine(tet, np.linspace(-0.9*n[0]*2*math.log(lam), 20, 300))
    q0, _, q1, _ = refine(r3, np.linspace(-0.9*n[0], 50, 300))
    rows.append(dict(kind=f['kind'], eps=f['eps'], seed=f['seed'],
                     geo=(g0/span, g1/span), tet=(t_0/span, t_1/span), r3=(q0/span, q1/span),
                     ratio=(f['ratio'], ratio1), tstar=(ts0, ts1)))

def dg(a, b):
    if a == b: return "no change"
    return f"rel. change {abs(a-b)/abs(a):.1e}"
lines = ["Polish of stretching_signature.py grid scans (t* and c refined by golden section to 1e-14, same objective)", ""]
for r in rows:
    lines.append(f"{r['kind']:8s} eps={r['eps']} seed={r['seed']}:")
    lines.append(f"   t*      {r['tstar'][0]:.10f} -> {r['tstar'][1]:.10f}   ratio r {r['ratio'][0]:.4f} -> {r['ratio'][1]:.4f}")
    for k in ('geo', 'tet', 'r3'):
        a, b = r[k]; lines.append(f"   {k:4s} err {a:.5f} -> {b:.5f}  ({dg(a, b)})")
hel = [r for r in rows if r['kind'] == 'helical']
w0 = sum(r['geo'][0] < r['tet'][0] for r in hel); w1 = sum(r['geo'][1] < r['tet'][1] for r in hel)
lines += ["", f"geometric beats tetration-type: {w0} of {len(hel)} helical runs before, {w1} of {len(hel)} after",
          f"Sabra step ratio r: {rows[0]['ratio'][0]:.4f} -> {rows[0]['ratio'][1]:.4f}",
          "Note: the arrival times are threshold crossings read off a saved time grid, so they carry",
          "sampling noise far above these digits; the polish only shows the grid scans had converged.",
          f"CPU time: {time.process_time() - t0:.2f} s"]
out = "\n".join(lines); print(out); open('polish_output.txt', 'w').write(out + "\n")
