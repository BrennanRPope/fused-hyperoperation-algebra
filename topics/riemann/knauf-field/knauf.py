"""Knauf-field: temperature as our roundness dial p, and finite-chain zeros. DRAFT."""
import numpy as np, mpmath as mp, json, time
def chain(n):
    phi = np.array([1], dtype=np.int64)
    for m in range(n):
        phi = np.concatenate([phi, phi + phi[::-1]])
    return phi
def levels(n):
    p = chain(n); c = np.bincount(p); q = np.nonzero(c)[0]
    return q.astype(float), c[q].astype(float)
def Z(q, c, S):
    """Z(s)=sum c q^-s for array S (complex)."""
    lq = np.log(q); out = np.empty(S.shape, complex); Sf = S.ravel(); o = out.ravel()
    for i in range(0, Sf.size, 2000):
        o[i:i+2000] = (c[None,:]*np.exp(-np.outer(Sf[i:i+2000], lq))).sum(1)
    return out
def zeros_by_winding(F, sig, tau):
    ph = np.angle(F)
    def d(a, b): return np.angle(np.exp(1j*(b-a)))
    w = (d(ph[:-1,:-1], ph[:-1,1:]) + d(ph[:-1,1:], ph[1:,1:]) + d(ph[1:,1:], ph[1:,:-1]) + d(ph[1:,:-1], ph[:-1,:-1]))/(2*np.pi)
    w = np.rint(w).astype(int); I, J = np.nonzero(w)
    return [(0.5*(sig[j]+sig[j+1]), 0.5*(tau[i]+tau[i+1]), int(w[i,j])) for i, j in zip(I, J)]
if __name__ == "__main__":
    t0 = time.time(); res = {}
    sig = np.arange(1.0, 2.61, 0.02); tau = np.arange(0.0, 36.0, 0.03)
    S = sig[None,:] + 1j*tau[:,None]
    for n in (12, 14, 16, 18):
        q, c = levels(n); F = Z(q, c, S)
        zs = zeros_by_winding(F, sig, tau)
        res[n] = dict(levels=len(q), zeros=[z for z in zs if z[2] > 0], poles=[z for z in zs if z[2] < 0])
        print(n, len(q), "zeros", len(res[n]['zeros']), "t", round(time.time()-t0,1))
        for z in sorted(res[n]['zeros'], key=lambda z: z[1])[:40]: print("   %.2f %.2f" % (z[0], z[1]))
    json.dump({str(k): v for k, v in res.items()}, open("finite_zeros.json", "w"))
