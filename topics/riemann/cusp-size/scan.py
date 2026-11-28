"""Fit Z_n(b) = F + A n^(1-b) + B n^(-b) over lengths 17..33 along the line Re b = 3/2 and compare
the fitted F with zeta(b-1)/zeta(b), and look at the cusp size A(b). DRAFT."""
import numpy as np, sys, json, mpmath as mp
sys.path.insert(0, "../length-smoothing"); from wide import lev
ns = list(range(17, 34)); L = {n: (np.log(lev(n)[0]), lev(n)[1]) for n in ns}
def fit(b):
    M = np.array([[1, n**(1-b), n**(-b)] for n in ns], complex); y = np.array([(c*np.exp(-b*lq)).sum() for lq, c in (L[n] for n in ns)])
    return np.linalg.lstsq(M, y, rcond=None)[0]
if __name__ == "__main__":
    out = dict(t=[], F=[], A=[], B=[], true=[], z1=[], zb=[])
    for t in np.arange(5.0, 45.01, 0.1):
        b = 1.5+1j*t; F, A, B = fit(b); tr = complex(mp.zeta(b-1)/mp.zeta(b))
        out["t"].append(float(t)); [out[k].append([float(v.real), float(v.imag)]) for k, v in (("F", F), ("A", A), ("B", B), ("true", tr), ("z1", complex(mp.zeta(b-1))), ("zb", complex(mp.zeta(b))))]
    json.dump(out, open("scan.json", "w"))
    F = np.array([complex(*v) for v in out["F"]]); T = np.array([complex(*v) for v in out["true"]]); A = np.array([complex(*v) for v in out["A"]])
    print("median |F - true| / |true|:", np.median(abs(F-T)/abs(T)), " median |true|", np.median(abs(T)))
    for k in range(0, len(F), 40): print("t %.1f |F|=%.3f |true|=%.3f |A|=%.3f" % (out["t"][k], abs(F[k]), abs(T[k]), abs(A[k])))
