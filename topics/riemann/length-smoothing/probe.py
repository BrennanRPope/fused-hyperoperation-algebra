"""Does averaging Knauf chains over length shrink the wobble of their zeros next to 3/2 + i gamma? DRAFT probe."""
import numpy as np, json, sys, mpmath as mp
sys.path.insert(0, "../knauf-field"); from knauf import levels
gam = [float(mp.zetazero(k).imag) for k in range(1, 6)]
L = {n: levels(n) for n in range(6, 25)}
def newton(terms, s, it=80):
    for _ in range(it):
        f = 0; fp = 0
        for q, c, wt in terms:
            t = wt*c*np.exp(-s*np.log(q)); f += t.sum(); fp += -(t*np.log(q)).sum()
        ds = f/fp; s -= ds
        if abs(ds) < 1e-13: break
    return s
def wob(terms): return [newton(terms, 1.5+1j*g).real - 1.5 for g in gam]
out = {}
for N in range(12, 25):
    row = {}
    row["single"] = wob([(*L[N], 1.0)])
    for k in (2, 4, 6):   # window of k+1 lengths, each normalised by 2^(-n/2) (Z_n's growth at beta = 3/2)
        row["window%d" % k] = wob([(*L[n], 2.0**(-n/2)) for n in range(N-k, N+1)])
    row["all_lengths"] = wob([(*L[n], 2.0**(-n/2)) for n in range(6, N+1)])
    out[N] = row
    print(N, "  ".join("%s %.4f" % (k, np.sqrt(np.mean(np.square(v)))) for k, v in row.items()), flush=True)
json.dump(dict(gam=gam, rows=out), open("probe.json", "w"), indent=1)
