import numpy as np, sys, json
sys.path.insert(0, "../knauf-field")
from knauf import chain
import mpmath as mp
mp.mp.dps = 60
def on_circle(n, beta, He, up):
    w = np.exp(-beta*He)
    coef = [mp.mpf(0)]*(n+1)
    for k in range(n+1): coef[k] = mp.fsum([mp.mpf(x) for x in w[up == k]])
    # palindromic: map to z + 1/z = 2 cos theta: roots on circle iff the reduced polynomial has all roots real in [-2,2]
    r = mp.polyroots(coef[::-1], maxsteps=400, extraprec=200)
    return max(abs(abs(z)-1) for z in r) < mp.mpf(10)**-20
out = {}
for n in range(15, 23):
    phi = chain(n).astype(float); H = np.log(phi); He = 0.5*(H + H[::-1]); up = np.array([bin(k).count("1") for k in range(2**n)])
    lo, hi = 1.3, 3.0
    if not on_circle(n, hi, He, up): print(n, "fails even at 3"); continue
    if on_circle(n, lo, He, up): print(n, "holds at", lo); out[n] = lo; continue
    for _ in range(30):
        mid = 0.5*(lo+hi)
        if on_circle(n, mid, He, up): hi = mid
        else: lo = mid
    out[n] = hi; print(n, "Lee-Yang circle holds for beta >= %.6f" % hi, flush=True)
json.dump(out, open("threshold_hi.json", "w"))
