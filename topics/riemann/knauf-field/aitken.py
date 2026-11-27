import numpy as np, mpmath as mp
from knauf import levels
L = {n: levels(n) for n in range(8, 25)}
def Zn(n, s):
    q, c = L[n]; return (c*np.exp(-s*np.log(q))).sum()
for s in [3.0, 2.5, 2.2, 1.8, 1.5+10j, 1.5+14.134725j, 1.3+5j, 1.7+20j]:
    seq = np.array([Zn(n, s) for n in range(8, 25)])
    d1 = np.diff(seq); ratio = d1[1:]/d1[:-1]
    ait = seq[2:] - d1[1:]**2/(d1[1:]-d1[:-1])
    target = complex(mp.zeta(s-1)/mp.zeta(s))
    print(s, "target", np.round(target,4), "| Z_24", np.round(seq[-1],3), "| ratio", np.round(ratio[-3:],3), "| Aitken", np.round(ait[-3:],4))
