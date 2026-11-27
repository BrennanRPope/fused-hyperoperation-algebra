"""Is the self-inversive (Cohn) threshold also where the doubled chain's field zeros actually leave the circle?
Direct test: all roots of Q(w) = P(w) + w^(n+1) P(1/w) on |w| = 1 to 1e-15, scanned on a beta grid."""
import numpy as np, json, mpmath as mp
from doubled import P_coef, two_roots, chain
mp.mp.dps = 50
out = {}
for n in range(6, 15):
    phi = chain(n).astype(float); up = np.array([bin(k).count("1") for k in range(2**n)])
    on = lambda b: max(abs(abs(z)-1) for z in two_roots(n, b, phi, up)) < 1e-15
    grid = [round(0.02*i, 2) for i in range(1, 101)]
    fails = [b for b in grid if not on(b)]
    out[n] = dict(last_fail=max(fails) if fails else None, n_fail=len(fails))
    print(n, "circle fails on grid at beta up to", out[n]["last_fail"], "(fails at", len(fails), "of 100 grid points)", flush=True)
json.dump(out, open("doubled_direct.json", "w"), indent=1)
