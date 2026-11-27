"""Wobble for long chains (to n = NMAX) with fixed and growing windows. DRAFT."""
import numpy as np, json, sys, mpmath as mp
sys.path.insert(0, "../knauf-field"); from knauf import levels
gam = [float(mp.zetazero(k).imag) for k in range(1, 6)]
def lev(n):
    if n <= 20: q, c = levels(n); return q.astype(float), c.astype(float)
    d = np.load("/tmp/knauf_levels/levels_%d.npz" % n); return d["q"].astype(float), d["c"].astype(float)
def newton(terms, s, it=60):
    for _ in range(it):
        f = 0; fp = 0
        for lq, c, wt in terms:
            t = wt*c*np.exp(-s*lq); f += t.sum(); fp -= (t*lq).sum()
        ds = f/fp; s -= ds
        if abs(ds) < 1e-12: break
    return s
if __name__ == "__main__":
    NMAX = int(sys.argv[1]); L = {}
    for n in range(4, NMAX+1): q, c = lev(n); L[n] = (np.log(q), c)
    Z15 = {n: float((c*np.exp(-1.5*lq)).sum()) for n, (lq, c) in L.items()}
    g = np.exp(np.mean([np.log(Z15[n+1]/Z15[n]) for n in range(NMAX-10, NMAX)]))
    print("growth %.5f" % g, flush=True)
    schemes = {"single": lambda N: 0, "k=6": lambda N: 6, "k=8": lambda N: 8, "k=12": lambda N: 12, "k=N/3": lambda N: N//3, "k=N/2": lambda N: N//2}
    out = dict(gam=gam, growth=g, rows={})
    for N in range(16, NMAX+1):
        line = []
        for name, kf in schemes.items():
            k = kf(N); ns = list(range(max(4, N-k), N+1))
            terms = [(*L[n], (i+1)*g**(-(n-N))) for i, n in enumerate(ns)]   # triangle weights
            z = [newton(terms, 1.5+1j*gg) for gg in gam]
            out["rows"]["%s N=%d" % (name, N)] = [[v.real-1.5, v.imag-gg] for v, gg in zip(z, gam)]
            line.append("%s %.4f" % (name, np.sqrt(np.mean([(v.real-1.5)**2 for v in z]))))
        print(N, "  ".join(line), flush=True)
        json.dump(out, open("wide.json", "w"), indent=1)
