"""Which step sizes t give a genuine pgf F_t (all Taylor coefficients >= 0)? Writes grid_output.json."""
import json, sys, time
import mpmath as mp
sys.path.insert(0, '.')
import flow
from flow import F, taylor, dalpha
mp.mp.dps = 30
flow.NSTEP = 1200
out = {}
ts = [mp.mpf(k) / 20 for k in range(1, 61)] + [mp.mpf(1) / m for m in (3, 4, 6, 8, 10, 30, 100)]
ts = sorted(set(ts))
for t in ts:
    cf = taylor(lambda s: F(t, s), M=128, rho=0.9)
    re = [float(x.real) for x in cf[:56]]
    negs = [k for k, v in enumerate(re) if v < -1e-22]
    out[mp.nstr(t, 6)] = {"coeffs": re, "neg": negs,
                          "min_rel": min(v / max(abs(w) for w in re) for v in re)}
    print(mp.nstr(t, 6), "first negative k:", negs[:6], flush=True)
ad = taylor(dalpha, M=128, rho=0.9)
out["alpha_prime"] = [float(x.real) for x in ad[:56]]
print("alpha' negative k:", [k for k, v in enumerate(out["alpha_prime"]) if v < -1e-22])
json.dump(out, open("grid_output.json", "w"))
