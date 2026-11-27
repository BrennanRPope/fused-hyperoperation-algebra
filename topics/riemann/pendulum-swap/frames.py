"""Rewrite Knauf's chain in dual (bond) spins b_i = s_i s_{i+1} (Kramers-Wannier for a line), repeat, and test the
rewind in each frame. DRAFT."""
import numpy as np, sys, mpmath as mp
sys.path.insert(0, "../knauf-field"); from knauf import chain
def walsh(f):
    a = f.astype(float).copy(); h = 1
    while h < len(a):
        for i in range(0, len(a), 2*h):
            x = a[i:i+h].copy(); y = a[i+h:i+2*h].copy(); a[i:i+h] = x+y; a[i+h:i+2*h] = x-y
        h *= 2
    return a/len(a)
def to_bonds(H, n):
    """new bits y_i = x_i xor x_{i+1} (i < n-1), y_{n-1} = x_{n-1}; returns H as a function of the new index."""
    x = np.arange(2**n); y = x ^ (x >> 1)          # Gray-code map: y_i = x_i xor x_{i+1}, top bit kept
    Hb = np.empty_like(H); Hb[y] = H; return Hb
def parts(H):
    return 0.5*(H + H[::-1]), 0.5*(H - H[::-1])     # flip-even, flip-odd (flip = complement all bits)
def odd_share(H):
    He, Ho = parts(H); return np.var(Ho)/(np.var(Ho) + np.var(He))
def doubled(H, n):
    He, Ho = parts(H); M = n - 2*np.array([bin(k).count("1") for k in range(2**n)])
    return np.concatenate([He+Ho, He-Ho]), np.concatenate([M+1, M-1]).astype(float)
def cohn_ok(H, n, b=1.5):
    up = np.array([bin(k).count("1") for k in range(2**n)])
    c = [mp.fsum([mp.mpf(v) for v in np.exp(-b*H[up == k])]) for k in range(n+1)]
    return min(abs(z) for z in mp.polyroots(c[::-1], maxsteps=400, extraprec=300))
def circle_dev(H, n, b=1.5):
    H2, M2 = doubled(H, n); j = np.rint((M2 + n + 1)/2).astype(int)
    c = np.bincount(j, weights=np.exp(-b*H2), minlength=n+2)
    r = mp.polyroots([mp.mpf(v) for v in c[::-1]], maxsteps=400, extraprec=300); return max(abs(abs(z)-1) for z in r)
def first_order(H, n, s0=1.5, T=36.0):
    H2, M2 = doubled(H, n); W = np.exp(-s0*H2); p = W/W.sum()
    a = np.sqrt((p*(H2-(p*H2).sum())**2).sum()/(p*M2**2).sum())
    j = np.rint((M2 + n + 1)/2).astype(int); c = np.bincount(j, weights=W, minlength=n+2)
    r = np.roots(c[::-1]); ph = np.angle(r) % (2*np.pi); zs = []
    for m in range(int(T*a/np.pi)+2):
        for q in ph:
            z = 1j*(q + 2*np.pi*m)/(2*a)
            if 0 < z.imag <= T: zs.append(z)
    z0 = np.array(zs); E = np.exp(np.multiply.outer(z0, a*M2))
    G = (W*H2*E).sum(-1); Fz = (W*a*M2*E).sum(-1); zp = -(z0*G)/Fz
    return z0.imag, zp.real
if __name__ == "__main__":
    import json; out = {}
    for n in (8, 10, 12):
        H = np.log(chain(n).astype(float)); rows = []
        Hk = H.copy()
        for k in range(4):
            y, sp = first_order(Hk, n)
            rows.append(dict(frame=k, odd_share=float(odd_share(Hk)), cohn_minroot=float(cohn_ok(Hk, n)), circle_dev=float(circle_dev(Hk, n)),
                             ferro=bool((-walsh(Hk))[1:].min() > -1e-12), speed_per_height=float(np.median(np.abs(sp)/y))))
            print(n, rows[-1], flush=True)
            Hk = to_bonds(Hk, n)
        out[n] = rows
    json.dump(out, open("frames.json", "w"), indent=1)
