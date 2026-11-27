import numpy as np, json, mpmath as mp
from knauf import levels
gam = [float(mp.zetazero(k).imag) for k in range(1,6)]
def newton(q, c, s, it=60):
    lq = np.log(q)
    for _ in range(it):
        t = c*np.exp(-s*lq); f = t.sum(); fp = -(t*lq).sum()
        ds = f/fp; s -= ds
        if abs(ds) < 1e-13: break
    return s, abs(f)
out = {}
for n in range(10, 25):
    q, c = levels(n); row = []
    for g in gam:
        s, r = newton(q, c, 1.5+1j*g)
        row.append([s.real, s.imag, r])
    # first Fisher (pinching) zero: start from previous
    out[n] = row
    print(n, len(q), "  ".join("%.5f%+.5fi" % (a, b-g) for (a, b, _), g in zip(row, gam)))
json.dump(dict(gam=gam, rows=out), open("refined.json","w"))
