"""Large step sizes: where is the first negative coefficient of F_t? Writes far_output.json."""
import json, sys
import mpmath as mp
sys.path.insert(0, '.')
import flow
from flow import F, taylor
mp.mp.dps = 45
flow.NSTEP = 1500
out = {}
for t in sys.argv[1:]:
    tt = mp.mpf(t)
    cf = taylor(lambda s: F(tt, s), M=1024, rho=0.97)
    re = [float(x.real) for x in cf[:480]]
    scale = [abs(v) for v in re]
    negs = [k for k, v in enumerate(re) if v < 0 and abs(v) > 1e-30]
    out[t] = {"coeffs": re, "neg": negs}
    print(t, "first negative k:", negs[:5], "count<480:", len(negs), flush=True)
json.dump(out, open("far_output_%s.json" % "_".join(sys.argv[1:]), "w"))
