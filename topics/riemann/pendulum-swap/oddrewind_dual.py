"""Odd-half rewind in the dual (Kramers-Wannier / staggered-flip) frame: X = cos(t) a M' - sin(t) Ho', where ' means
the frame whose flip is the staggered flip. Does every zero stay on the line? Also the full rewind's first-order
speed in both frames. DRAFT."""
import numpy as np, sys, json
sys.path.insert(0, "../rewind-dial"); sys.path.insert(0, "../knauf-field")
from knauf import chain
import movie
from frames import to_bonds, parts
def setup_frame(H, n, sigma0=1.5):
    He, Ho = parts(H); M = n - 2*np.array([bin(k).count("1") for k in range(2**n)])
    H2 = np.concatenate([He+Ho, He-Ho]); M2 = np.concatenate([M+1, M-1]).astype(float); Ho2 = -np.concatenate([Ho, Ho])
    if np.corrcoef(Ho2, M2)[0, 1] < 0: Ho2 = -Ho2
    W = np.exp(-sigma0*H2); key = np.round(Ho2, 12) + 1000*M2; u, inv = np.unique(key, return_inverse=True); cnt = np.bincount(inv)
    p = W/W.sum(); a = np.sqrt((p*Ho2**2).sum()/(p*M2**2).sum())
    return dict(w=np.bincount(inv, weights=W), H=np.bincount(inv, weights=Ho2)/cnt, M=np.bincount(inv, weights=M2)/cnt, a=a, n=n)
if __name__ == "__main__":
    out = {}
    for n in (8, 10):
        H = to_bonds(np.log(chain(n).astype(float)), n); S = setup_frame(H, n); out[n] = {}
        for deg in range(0, 91, 15):
            z = movie.zeros_box(S, np.radians(deg), h=0.02); z = z[(np.abs(z.real) <= 1) & (z.imag > 0.3) & (z.imag <= 36)]
            z = np.unique(np.round(z, 8)); out[n][deg] = [[float(v.real), float(v.imag)] for v in z]
            print("dual-frame odd-half rewind n=%d %2d deg: %d zeros, max |Re z| = %.1e" % (n, deg, z.size, np.abs(z.real).max() if z.size else -1), flush=True)
    json.dump(out, open("oddrewind_dual.json", "w"))
