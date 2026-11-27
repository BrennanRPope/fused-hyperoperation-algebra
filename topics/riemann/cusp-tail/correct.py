"""Remove the cusp tail: model Z_n(b) = F(b) + A(b) n^(1-b) [+ B(b) n^(-b)] and solve for F over a few lengths.
Zeros of the estimated F against the zeta heights, compared with single chains. DRAFT."""
import numpy as np, sys, json, mpmath as mp
sys.path.insert(0, "../length-smoothing"); from wide import lev
gam = [float(mp.zetazero(k).imag) for k in range(1, 9)]
NMAX = int(sys.argv[1]) if (__name__ == "__main__" and len(sys.argv) > 1) else 33
L = {n: (np.log(lev(n)[0]), lev(n)[1]) for n in range(8, NMAX+1)}
def Z(n, b): lq, c = L[n]; return (c*np.exp(-b*lq)).sum()
def F(ns, b, model):
    cols = {"none": [], "cusp": [lambda n: n**(1-b)], "cusp2": [lambda n: n**(1-b), lambda n: n**(-b)],
            "cusp3": [lambda n: n**(1-b), lambda n: n**(-b), lambda n: n**(-1-b)]}[model]
    M = np.array([[1]+[f(n) for f in cols] for n in ns], complex); y = np.array([Z(n, b) for n in ns])
    return np.linalg.lstsq(M, y, rcond=None)[0][0]
def zero(ns, g, model):
    b = 1.5+1j*g
    for _ in range(40):
        f = F(ns, b, model); h = 1e-6; d = (F(ns, b+h, model)-f)/h; db = f/d; b -= db
        if abs(db) < 1e-11: break
    return b
if __name__ == "__main__":
    out = dict(gam=gam, rows={})
    schemes = [("single", "none", 0), ("cusp pair", "cusp", 1), ("cusp 3", "cusp2", 2), ("cusp 4", "cusp3", 3), ("cusp fit 8", "cusp2", 7)]
    for N in range(16, NMAX+1):
        line = []
        for name, model, k in schemes:
            ns = list(range(N-k, N+1)); z = [zero(ns, g, model) for g in gam]
            out["rows"]["%s N=%d" % (name, N)] = [[v.real-1.5, v.imag-g] for v, g in zip(z, gam)]
            line.append("%s %.5f" % (name, np.sqrt(np.mean([(v.real-1.5)**2 for v in z]))))
        print(N, "  ".join(line), flush=True)
    json.dump(out, open("correct.json", "w"), indent=1)
