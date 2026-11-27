"""Is the swing the cusp? Test whether the per-length turn of Z_n(rho) goes like -gamma*log((n+a)/(n+a-1)). DRAFT."""
import numpy as np, sys, mpmath as mp
sys.path.insert(0, "../length-smoothing"); from wide import lev
gam = [float(mp.zetazero(k).imag) for k in range(1, 9)]
L = {n: lev(n) for n in range(8, 34)}
for g in gam[:5]:
    E = np.array([(c*np.exp(-(1.5+1j*g)*np.log(q))).sum() for q, c in (L[n] for n in range(8, 34))])
    ang = np.angle(E[1:]/E[:-1]); ns = np.arange(9, 34)
    print("gamma %.2f" % g)
    print("  n    turn    -g/n   -g*log(n+1/n)  |ratio|")
    for n, a, r in zip(ns, ang, np.abs(E[1:]/E[:-1])):
        if n % 3 == 0: print("  %2d  %+.3f  %+.3f  %+.3f  %.3f" % (n, a, -g/n, -g*np.log((n+1)/n), r))
