"""How close does each method get to zeta(b-1)/zeta(b) on Re b = 3/2? DRAFT."""
import numpy as np, json, mpmath as mp
from scan import L, ns
ts = np.arange(5.0, 45.01, 0.5); out = dict(t=ts.tolist(), single=[], cusp1=[], cusp2=[])
for t in ts:
    b = 1.5+1j*t; tr = complex(mp.zeta(b-1)/mp.zeta(b)); y = np.array([(c*np.exp(-b*lq)).sum() for lq, c in (L[n] for n in ns)])
    out["single"].append(float(abs(y[-1]-tr)))
    for key, cols in (("cusp1", [lambda n: 1, lambda n: n**(1-b)]), ("cusp2", [lambda n: 1, lambda n: n**(1-b), lambda n: n**(-b)])):
        M = np.array([[f(n) for f in cols] for n in ns], complex); out[key].append(float(abs(np.linalg.lstsq(M, y, rcond=None)[0][0]-tr)))
for k in ("single", "cusp1", "cusp2"): print("%-7s median |error| %.5f  max %.4f" % (k, np.median(out[k]), np.max(out[k])))
json.dump(out, open("compare.json", "w"))
