"""Closing rate of the natural constants (linear rule), rungs 4..K.
theta_k = (2 - c_k) * (ln 4)^k  should settle if the gap shrinks like (1/ln 4)^k.
Vectorised float evaluation; c_k by bisection on 'min over x of H_k(b,x) - x <= 0'."""
import numpy as np, json, math, sys
from scipy.optimize import minimize_scalar
np.seterr(over='ignore', invalid='ignore')
def H(k, b, x):
    x = np.asarray(x, dtype=float)
    if k == 3: return np.power(b, x)
    out = np.empty_like(x)
    lo = x <= 0
    out[lo] = 1 + x[lo]
    hi = ~lo & (x <= 60)
    if hi.any(): out[hi] = H(k - 1, b, H(k, b, x[hi] - 1))
    out[x > 60] = np.inf
    return out
def gap(k, b, grid):
    v = H(k, b, grid) - grid
    i = int(np.nanargmin(v))
    a, c = grid[max(i - 2, 0)], grid[min(i + 2, len(grid) - 1)]
    r = minimize_scalar(lambda t: float(H(k, b, np.array([t]))[0] - t), bounds=(a, c), method='bounded', options={'xatol': 1e-13})
    return min(r.fun, v[i]), r.x
grid = np.linspace(1.0005, 4.5, 14000)
mp_ref = {int(k): float(v['c']) for k, v in json.load(open('constants.json'))['linear'].items()}
q = math.log(4)
K = int(sys.argv[1]) if len(sys.argv) > 1 else 60
res = []; prev = None
for k in range(4, K + 1):
    lo, hi = (1.5, 2.0) if prev is None else (2 - 1.6 * (2 - prev), 2 - 0.4 * (2 - prev))
    for _ in range(64):
        mid = (lo + hi) / 2
        if gap(k, mid, grid)[0] <= 0: lo = mid
        else: hi = mid
    c = lo; x = gap(k, c, grid)[1]
    eps = 2 - c; th = eps * q ** k
    line = f"k={k:2d} c_k={c:.15f} 2-c_k={eps:.6e} theta={th:.6f} ratio={'%.6f' % (eps/(2-prev)) if prev else '   -    '} kiss x={x:.4f}"
    if k in mp_ref: line += f"  (40-digit check diff {abs(c-mp_ref[k]):.1e})"
    print(line, flush=True); res.append((k, c, x)); prev = c
json.dump([{"k": k, "c": repr(c), "x": x} for k, c, x in res], open('rate_scan.json', 'w'), indent=0)
print("1/ln4 =", 1 / q)
