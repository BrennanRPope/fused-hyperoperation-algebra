"""Compute the first 5500 zeta zeros with mpmath (Odlyzko's tables are not reachable here)."""
import mpmath, json, sys
from multiprocessing import Pool
mpmath.mp.dps = 20
def z(n):
    return n, float(mpmath.zetazero(n).imag)
if __name__ == "__main__":
    N = 5500
    with Pool(4) as p:
        out = dict(p.imap_unordered(z, range(1, N+1), chunksize=20))
    g = [out[n] for n in range(1, N+1)]
    json.dump({"N": N, "gamma": g, "source": "mpmath.zetazero, dps=20"}, open("zeros.json", "w"))
    print("done", N, g[0], g[-1])
