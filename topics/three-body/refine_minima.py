"""DRAFT: take every local minimum of the landscape (absolute and rotation-free), shoot it to an exact
(relative) periodic orbit, and save minima.json. Run after run_landscape.py."""
import sys, os, json
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import numpy as np
from multiprocessing import Pool
from scipy.ndimage import minimum_filter
from tbpy import refine

L = np.load(os.path.join(HERE, "landscape.npz"))
def minima(F, thr):
    G = np.where(np.isfinite(F), F, 9); m = (G == minimum_filter(G, size=7)) & (G < thr)
    return list(zip(*np.nonzero(m)))
jobs = []
for kind, key, tkey in (("absolute", "dabs", "tabs"), ("rotation-free", "drot", "trot")):
    for (i, j) in minima(L[key], 0.08):
        jobs.append((kind, float(L["p1"][i, j]), float(L["p2"][i, j]), float(L[tkey][i, j]), float(L["arot"][i, j]), float(L[key][i, j])))

def work(job):
    kind, p1, p2, T, a, d0 = job
    rot = kind == "rotation-free"
    try:
        q, r = refine(p1, p2, T, rotating=rot, alpha=a if rot else 0.0)
    except Exception:
        return None
    return dict(kind=kind, start=[p1, p2, T, d0], p1=float(q[0]), p2=float(q[1]), T=float(q[2]),
                alpha=float(q[3]) if rot else 0.0, res=float(r))

if __name__ == "__main__":
    with Pool(4) as pool:
        out = [o for o in pool.map(work, jobs, chunksize=2) if o]
    json.dump(dict(n_jobs=len(jobs), results=out), open(os.path.join(HERE, "minima.json"), "w"))
    print(len(jobs), "minima;", sum(o["res"] < 1e-8 for o in out), "converged")
