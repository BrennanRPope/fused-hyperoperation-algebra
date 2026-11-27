"""Where are Knauf's zeros that sit right of Re beta = 3/2 (the side the mirror needs empty)? DRAFT."""
import numpy as np, json
from farend import E_terms
import mirror
def zeros(w, X, x0, x1, y0=0.3, y1=36.0, h=0.01):
    xs = np.arange(x0, x1, h); ys = np.arange(y0, y1, h); Z = xs[None, :] + 1j*ys[:, None]
    v = np.abs(mirror.ev(w, X, Z.ravel())).reshape(Z.shape); found = []
    for i in range(1, len(ys)-1):
        for j in range(1, len(xs)-1):
            if v[i, j] <= v[i-1:i+2, j-1:j+2].min():
                z = Z[i, j]
                for _ in range(60):
                    f = mirror.ev(w, X, np.array([z]))[0]; d = mirror.ev(w*X, X, np.array([z]))[0]; z = z - f/d
                if abs(mirror.ev(w, X, np.array([z]))[0]) < 1e-9*w.sum() and x0 < z.real < x1 and y0 < z.imag < y1: found.append(z)
    return np.unique(np.round(np.array(found), 7))
if __name__ == "__main__":
    out = {}
    for n in (8, 10, 12, 14, 16):
        w, X = E_terms(n); zr = zeros(w, X, 0.0, 1.0); zl = zeros(w, X, -1.0, 0.0)
        out[n] = dict(right=[[float(z.real), float(z.imag)] for z in zr], left=[[float(z.real), float(z.imag)] for z in zl])
        print(n, "right of line (Re beta - 3/2, Im):", [(round(z.real, 3), round(z.imag, 2)) for z in zr], "| left count", len(zl), flush=True)
    json.dump(out, open("farend_pos.json", "w"), indent=1)
