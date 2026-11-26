"""First N nontrivial zeta zeros (heights gamma_n), mpmath zetazero, 20 digits. -> zeros.json"""
import json, sys
from multiprocessing import Pool
import mpmath as mp
N = int(sys.argv[1]) if len(sys.argv) > 1 else 5000
def g(n):
    mp.mp.dps = 25
    return mp.nstr(mp.zetazero(n).imag, 22)
if __name__ == "__main__":
    with Pool(4) as P:
        gs = P.map(g, range(1, N + 1), chunksize=50)
    json.dump({"N": N, "gamma": gs}, open("zeros.json", "w"))
    print("wrote", N, "zeros; last", gs[-1])
