"""Temperature mirror on the chain's safe-half rewind. DRAFT.
The compass spin makes F_theta(z) = G_theta(z) + G_theta(-z) exactly, where G is the sum over the compass-up copy.
Hermite-Biehler test: if G has no zeros in Re z >= 0 (and is lopsided the right way), then |G(z)| > |G(-z)| for Re z > 0,
so every zero of F sits on Re z = 0, is simple, and interlaces with the zeros of G(z) - G(-z).
We count zeros of G in the right half (the wrong side) (box) and compare with F's off-line zeros. Two splits:
  'compass': G = compass-up copy (the chain's own split);  'sorted': each mirror pair split by the sign of X (most lopsided)."""
import numpy as np, sys, json
sys.path.insert(0, "../knauf-field"); sys.path.insert(0, "../pendulum-swap")
from knauf import chain
from frames import to_bonds, parts
def halves(H, n, sigma0=1.5):
    He, Ho = parts(H); M = n - 2*np.array([bin(k).count("1") for k in range(2**n)])
    M0 = (M + 1).astype(float); O0 = -Ho.copy(); W0 = np.exp(-sigma0*(He + Ho))
    M2 = np.concatenate([M+1, M-1]).astype(float); O2 = -np.concatenate([Ho, Ho])
    if np.corrcoef(O2, M2)[0, 1] < 0: O0, O2 = -O0, -O2
    W2 = np.exp(-sigma0*np.concatenate([He+Ho, He-Ho])); p = W2/W2.sum(); a = np.sqrt((p*O2**2).sum()/(p*M2**2).sum())
    return dict(M0=M0, O0=O0, W0=W0, a=a, n=n)
def frame(name, n):
    H = np.log(chain(n).astype(float)); return halves(H if name == "spin" else to_bonds(H, n), n)
def terms(S, th, split):
    X = np.cos(th)*S["a"]*S["M0"] + np.sin(th)*S["O0"]
    if split == "sorted": X = np.abs(X)           # oriented like the compass split at 0 deg (its zeros sit in Re z < 0)
    key = np.round(X, 11); u, inv = np.unique(key, return_inverse=True)
    return np.bincount(inv, weights=S["W0"]), np.bincount(inv, weights=X)/np.bincount(inv)
def ev(w, X, z):
    out = np.empty(z.size, complex)
    for i in range(0, z.size, 2000): out[i:i+2000] = (w[None, :]*np.exp(np.outer(z[i:i+2000], X))).sum(1)
    return out
def winding(w, X, x0, x1, y0, y1, N=9000):
    e = np.concatenate([x1 + 1j*np.linspace(y0, y1, N), np.linspace(x1, x0, 900) + 1j*y1,
                        x0 + 1j*np.linspace(y1, y0, N), np.linspace(x0, x1, 900) + 1j*y0])
    ph = np.unwrap(np.angle(ev(w, X, e))); return int(np.rint((ph[-1] - ph[0])/(2*np.pi)))
def bad_zeros(S, th, split, R=4.0, y0=0.3, y1=36.0):
    """zeros of G on the wrong side (Re z > 0). Zero of them => Hermite-Biehler holds in the box => F's zeros on the line."""
    w, X = terms(S, th, split); return winding(w, X, 1e-9, R, y0, y1)
def good_zeros(S, th, split, R=4.0, y0=0.3, y1=36.0):
    w, X = terms(S, th, split); return winding(w, X, -R, -1e-9, y0, y1)
if __name__ == "__main__":
    fr, n, split = sys.argv[1], int(sys.argv[2]), sys.argv[3]
    degs = [int(d) for d in sys.argv[4].split(",")] if len(sys.argv) > 4 else list(range(0, 41)) + list(range(45, 91, 5))
    S = frame(fr, n)
    with open("mz_%s_%d_%s.jsonl" % (fr, n, split), "a") as f:
        for d in degs:
            f.write(json.dumps([d, bad_zeros(S, np.radians(d), split)]) + "\n"); f.flush()
