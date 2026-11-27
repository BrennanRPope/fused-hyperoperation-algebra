"""Cusp fit over a long window (half the chain), against the triangle window from length-smoothing. DRAFT."""
import numpy as np, json, sys
from correct import zero, gam
NMAX = 33; out = dict(gam=gam, rows={})
w = json.load(open("../length-smoothing/wide.json"))
for N in range(20, NMAX+1):
    line = []
    for name, model in [("cusp half", "cusp"), ("cusp2 half", "cusp2")]:
        ns = list(range(N-N//2, N+1)); z = [zero(ns, g, model) for g in gam]
        out["rows"]["%s N=%d" % (name, N)] = [[v.real-1.5, v.imag-g] for v, g in zip(z, gam)]
        line.append("%s %.5f" % (name, np.sqrt(np.mean([(v.real-1.5)**2 for v in z]))))
    tri = w["rows"]["k=N/2 N=%d" % N]; line.append("tri half(5 zeros) %.5f" % np.sqrt(np.mean([r[0]**2 for r in tri])))
    c5 = out["rows"]["cusp half N=%d" % N][:5]; line.append("cusp half(5) %.5f" % np.sqrt(np.mean([r[0]**2 for r in c5])))
    print(N, "  ".join(line), flush=True)
json.dump(out, open("combo.json", "w"), indent=1)
