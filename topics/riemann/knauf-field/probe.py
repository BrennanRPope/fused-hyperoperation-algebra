import numpy as np, mpmath as mp, math, time
def chain(n):
    phi = np.array([1], dtype=np.int64)
    for m in range(n):
        phi = np.concatenate([phi, phi + phi[::-1]])
    return phi
for n in [3, 10, 20]:
    p = chain(n); print(n, p[:8] if n == 3 else (len(p), p.max()))
p = chain(22)
c = np.bincount(p); q = np.nonzero(c)[0]; cq = c[q].astype(float)
print("distinct q", len(q))
for s in [3.0, 2.5]:
    Zn = np.sum(cq * q**-s); ex = float(mp.zeta(s - 1) / mp.zeta(s)); print("s", s, Zn, ex)
