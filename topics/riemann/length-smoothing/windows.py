"""Wobble of the chain's zeros next to 3/2 + i gamma, after averaging Z_n over a window of lengths. DRAFT.
Each Z_n is normalised by its growth g^n at beta = 3/2 (g measured), then averaged with flat or triangular weights."""
import numpy as np, json, sys, mpmath as mp
sys.path.insert(0, "../knauf-field"); from knauf import levels
gam = [float(mp.zetazero(k).imag) for k in range(1, 6)]
L = {n: levels(n) for n in range(4, 25)}
Z15 = {n: float((c*q**-1.5).sum()) for n, (q, c) in L.items()}
g = np.exp(np.mean([np.log(Z15[n+1]/Z15[n]) for n in range(14, 24)]))
def newton(terms, s, it=80):
    for _ in range(it):
        f = 0; fp = 0
        for q, c, wt in terms:
            t = wt*c*np.exp(-s*np.log(q)); f += t.sum(); fp += -(t*np.log(q)).sum()
        ds = f/fp; s -= ds
        if abs(ds) < 1e-13: break
    return s
if __name__ == "__main__":
    out = dict(gam=gam, growth=g, rows={})
    print("growth per length at beta=3/2: %.5f (sqrt2 = %.5f)" % (g, 2**0.5), flush=True)
    for kind in ("flat", "tri"):
        for k in range(0, 13, 2):
            for N in range(16, 25):
                ns = list(range(max(4, N-k), N+1)); wts = [1.0 if kind == "flat" else (i+1) for i in range(len(ns))]
                terms = [(*L[n], w*g**(-n)) for n, w in zip(ns, wts)]
                off = [newton(terms, 1.5+1j*gg) for gg in gam]
                out["rows"]["%s k=%d N=%d" % (kind, k, N)] = [[z.real-1.5, z.imag-gg] for z, gg in zip(off, gam)]
            rms = [np.sqrt(np.mean([r[0]**2 for r in out["rows"]["%s k=%d N=%d" % (kind, k, N)]])) for N in range(16, 25)]
            print(kind, "k=%2d" % k, " ".join("%.4f" % v for v in rms), "| mean %.4f" % np.mean(rms), flush=True)
    json.dump(out, open("windows.json", "w"), indent=1)
