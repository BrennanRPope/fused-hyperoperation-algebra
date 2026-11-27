"""Magnetic Lee-Yang circle test on Knauf's chain: Z(s,w)=sum_states phi^-s w^(number of 1-spins)."""
import numpy as np, json
from knauf import chain
def popcounts(n):
    k = np.arange(2**n); return np.array([bin(x).count("1") for x in k])
out = {}
for n in (8, 10, 12):
    phi = chain(n).astype(float); m = popcounts(n)
    # sanity: does chain index k correspond to spin string k? (Knauf: phi_n(k) for k in binary)
    res = {}
    for s in (0.5, 1.0, 1.5, 2.0, 3.0, 5.0):
        coef = np.bincount(m, weights=phi**(-s), minlength=n+1)
        r = np.roots(coef[::-1])
        res[s] = dict(absmin=float(np.abs(r).min()), absmax=float(np.abs(r).max()), roots=[[float(z.real), float(z.imag)] for z in r])
        print(n, s, "|w| in [%.4f, %.4f]" % (np.abs(r).min(), np.abs(r).max()), " palindromic coef?", np.allclose(coef, coef[::-1]*coef[-1]/coef[0]) )
    out[n] = res
json.dump(out, open("circle.json","w"))
