"""What is the cusp size A(b)? Compare with zeta functions along the line; and test the parabolic-fingerprint
prediction that the Farey cusp (x/(1-x), fingerprint 0,0: sits in a flow) gives no log terms. DRAFT."""
import numpy as np, json, mpmath as mp, sys
from scan import L, ns
d = json.load(open("scan.json")); t = np.array(d["t"]); A = np.array([complex(*v) for v in d["A"]])
z1 = np.array([complex(*v) for v in d["z1"]]); zb = np.array([complex(*v) for v in d["zb"]]); b = 1.5+1j*t
cands = {"1/zeta(b)": 1/zb, "zeta(b-1)": z1, "1/(b-1)": 1/(b-1), "1/((b-1) zeta(b))": 1/((b-1)*zb), "zeta(b)/(b-1)": zb/(b-1), "1/((b-1) zeta(b)^2)": 1/((b-1)*zb**2), "1/zeta(b)^2": 1/zb**2}
out = {}
for k, v in cands.items():
    r = A/v; c = np.median(r.real)+1j*np.median(r.imag); spread = np.median(abs(r-c))/abs(c)
    out[k] = [float(c.real), float(c.imag), float(spread)]; print("%-22s A/cand median %.4f%+.4fi  relative spread %.3f" % (k, c.real, c.imag, spread))
def fit(bb, logs):
    cols = [lambda n: 1, lambda n: n**(1-bb), lambda n: n**(-bb)] + ([lambda n: n**(-bb)*np.log(n), lambda n: n**(1-bb)*np.log(n)] if logs else [])
    M = np.array([[f(n) for f in cols] for n in ns], complex); y = np.array([(c*np.exp(-bb*lq)).sum() for lq, c in (L[n] for n in ns)])
    return np.linalg.lstsq(M, y, rcond=None)[0]
ratios = []
for tt in np.arange(6, 44, 2.0):
    bb = 1.5+1j*tt; p = fit(bb, True); ratios.append([float(tt), float(abs(p[4])/abs(p[1])), float(abs(p[3])/abs(p[2]))])
out["log_terms"] = ratios; print("log-term size relative to its plain partner, median: n^(1-b)log n %.3f, n^(-b)log n %.3f" % (np.median([r[1] for r in ratios]), np.median([r[2] for r in ratios])))
json.dump(out, open("shape.json", "w"), indent=1)
