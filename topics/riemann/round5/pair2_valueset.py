"""Round 5, pair 2 add-on: the path of P(1.2 + it), 0 < t < 200, against the point 1 and the disc of radius P(1.2)."""
import json, numpy as np, mpmath as mp
mp.mp.dps = 15
ts = np.arange(0.0, 200.0, 0.25)
vals = [complex(mp.primezeta(mp.mpc(1.2, t))) for t in ts]
d = np.abs(np.array(vals) - 1)
out = dict(sigma=1.2, P_sigma=float(mp.primezeta(1.2)), t=ts.tolist(), re=[v.real for v in vals], im=[v.imag for v in vals],
           min_dist_t_gt_1=float(d[ts > 1].min()), t_at_min=float(ts[ts > 1][d[ts > 1].argmin()]))
json.dump(out, open('pair2_valueset.json', 'w'), indent=1)
print({k: v for k, v in out.items() if k not in ('t', 're', 'im')})
