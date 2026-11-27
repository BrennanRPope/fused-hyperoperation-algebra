"""The mirror at the far end (90 deg) on the FULL chain: E(z) = Z_n(3/2 + z) = sum phi^(-3/2 - z).
Mirrored chain F(z) = E(z) + E(-z) = Z(3/2+z) + Z(3/2-z). Its zeros are forced onto Re z = 0 iff E has no zeros on one
closed side, i.e. Knauf's Z_n has no zeros on one side of Re beta = 3/2. Count Z_n's zeros on each side. DRAFT."""
import numpy as np, sys, json
sys.path.insert(0, "../knauf-field"); from knauf import chain
from mirror import winding
def E_terms(n):
    phi = chain(n).astype(float); u, c = np.unique(phi, return_counts=True); L = np.log(u)
    return c*np.exp(-1.5*L), -L
if __name__ == "__main__":
    out = {}
    for n in (6, 8, 10, 12, 14, 16):
        w, X = E_terms(n)
        right = winding(w, X, 1e-9, 1.0, 0.3, 36.0); left = winding(w, X, -1.0, -1e-9, 0.3, 36.0)
        # mirrored chain: zeros on the line vs total
        y = np.arange(0.3, 36.0, 0.002); import mirror
        f = (mirror.ev(w, X, 1j*y) + mirror.ev(w, X, -1j*y)).real; online = int(np.sum(np.sign(f[1:]) != np.sign(f[:-1])))
        wm = np.concatenate([w, w]); Xm = np.concatenate([X, -X]); total = winding(wm, Xm, -1.0, 1.0, 0.3, 36.0)
        out[n] = dict(Z_zeros_right_of_line=right, Z_zeros_left_of_line=left, mirrored_on_line=online, mirrored_total=total)
        print(n, out[n], flush=True)
    json.dump(out, open("farend.json", "w"), indent=1)
