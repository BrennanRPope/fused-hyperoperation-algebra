"""Variant rewind: rotate the magnetic field toward the flip-ODD half of the energy only (Ho), which keeps the
flip pairing z <-> -conj(z) at every angle. Does the line survive? DRAFT."""
import numpy as np, json, sys, time
sys.path.insert(0, "../knauf-field"); from knauf import chain

from rotate import newton
def setup_odd(n, sigma0=1.5):
    phi = chain(n).astype(float); H = np.log(phi); He = 0.5*(H+H[::-1]); Ho = 0.5*(H-H[::-1])
    M = n - 2*np.array([bin(k).count("1") for k in range(2**n)])
    H2 = np.concatenate([H, H[::-1]]); M2 = np.concatenate([M+1, M-1]).astype(float); Ho2 = -np.concatenate([Ho, Ho])   # sign chosen so it reinforces the magnetic field (they are anti-correlated)
    W = np.exp(-sigma0*H2); key = np.round(Ho2, 12) + 1000*M2; u, inv = np.unique(key, return_inverse=True)
    cnt = np.bincount(inv); w = np.bincount(inv, weights=W)
    p = W/W.sum(); a = np.sqrt((p*Ho2**2).sum()/(p*M2**2).sum())
    return dict(w=w, H=np.bincount(inv, weights=Ho2)/cnt, M=np.bincount(inv, weights=M2)/cnt, a=a, n=n)
if __name__ == "__main__":
    import movie
    ns = [int(x) for x in sys.argv[1:]] or [8, 10, 12]; tag = "_".join(map(str, ns))
    out = {}
    for n in ns:
        S = setup_odd(n); out[n] = {}
        for deg in range(0, 91, 10):
            z = movie.zeros_box(S, np.radians(deg)); z = z[(np.abs(z.real) <= 1) & (z.imag > 0.3) & (z.imag <= 36)]
            z = np.unique(np.round(z, 8)); out[n][deg] = [[float(v.real), float(v.imag)] for v in z]
            print("odd path n=%d %2d deg: %d zeros, max |Re z| = %.1e" % (n, deg, z.size, np.abs(z.real).max() if z.size else -1), flush=True)
    json.dump(out, open("oddpath_%s.json" % tag, "w"))
