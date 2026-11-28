"""Rewind with the next known terms: the Taylor term of h gives c*L, and where a zeta zero sits at the same height
it pairs with it into L log L. Fit F + a L + c L log L over small k and compare with zeta(b-1)/zeta(b). DRAFT."""
import numpy as np, mpmath as mp, json
from rewind import Fk
ks = np.array([0.004, 0.003, 0.002, 0.0015, 0.001, 0.0007, 0.0005, 0.0003, 0.0002])
def rewind(b, model):
    rows = []
    for k in ks:
        f, L = Fk(b, k); rows.append((f, L))
    y = np.array([r[0] for r in rows]); L = np.array([r[1] for r in rows])
    cols = {"L": [np.ones_like(L), L], "L+LlogL": [np.ones_like(L), L, L*np.log(L)], "L+LlogL+L2": [np.ones_like(L), L, L*np.log(L), L**2]}[model]
    return np.linalg.lstsq(np.array(cols, complex).T, y, rcond=None)[0][0], y[-1]
if __name__ == "__main__":
    out = dict(t=[], raw=[], L=[], LlogL=[], LlogL2=[])
    for t in np.arange(5.0, 45.01, 0.5):
        b = 1.5+1j*t; tr = complex(mp.zeta(b-1)/mp.zeta(b)); out["t"].append(float(t))
        e1, raw = rewind(b, "L"); e2, _ = rewind(b, "L+LlogL"); e3, _ = rewind(b, "L+LlogL+L2")
        out["raw"].append(abs(raw-tr)); out["L"].append(abs(e1-tr)); out["LlogL"].append(abs(e2-tr)); out["LlogL2"].append(abs(e3-tr))
    for k in ("raw", "L", "LlogL", "LlogL2"): print("%-7s median error %.2e  max %.2e" % (k, np.median(out[k]), np.max(out[k])))
    json.dump(out, open("extrap.json", "w"))
