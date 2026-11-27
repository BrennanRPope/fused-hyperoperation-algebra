import numpy as np, mpmath as mp, json
from knauf import levels
taus = np.linspace(5, 35, 121)
out = {}
for sg in (1.2, 1.5, 1.8, 2.2):
    tgt = np.array([complex(mp.zeta(sg+1j*t-1)/mp.zeta(sg+1j*t)) for t in taus])
    errs = []
    for n in range(10, 25):
        q, c = levels(n); S = sg + 1j*taus
        Zv = (c[None,:]*np.exp(-np.outer(S, np.log(q)))).sum(1)
        errs.append(float(np.sqrt(np.mean(np.abs(Zv-tgt)**2))))
    out[sg] = errs
    print(sg, " ".join("%.4f" % e for e in errs))
json.dump(out, open("errscan.json","w"))
