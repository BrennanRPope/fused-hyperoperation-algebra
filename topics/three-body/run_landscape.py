"""DRAFT: return-distance landscape over the Suvakov slice (p1, p2), absolute and rotation-free. Writes landscape.npz."""
import sys, time, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from tbpy import *
import numpy as np
N = int(sys.argv[1]) if len(sys.argv) > 1 else 600
Tmax = float(sys.argv[2]) if len(sys.argv) > 2 else 40.0
g = (np.arange(N) + 0.5) / N * 0.8
P1, P2 = np.meshgrid(g, g, indexing="xy")
mask = P1 ** 2 + P2 ** 2 < 0.8333
t = time.time()
outs = landscape(P1[mask], P2[mask], Tmax, 1e-10, 0.5)
full = []
for o in outs:
    a = np.full(P1.shape, np.nan); a[mask] = o; full.append(a)
np.savez_compressed(os.path.join(os.path.dirname(os.path.abspath(__file__)), "landscape.npz"),
                    p1=P1, p2=P2, dabs=full[0], tabs=full[1], drot=full[2], trot=full[3], arot=full[4], Tmax=Tmax)
print("done", N, Tmax, time.time() - t)
