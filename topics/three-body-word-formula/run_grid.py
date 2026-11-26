"""DRAFT. Letter itineraries on the collinear slice: first K syzygy letters for each start (v1, v2).
Output: grid.npz (N x N, K letters as ints +-1, +-2; count of letters found). Takes ~10 min on 4 cores."""
import sys, os, time, numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from tbw import grid
N, K, TMAX = int(sys.argv[1]) if len(sys.argv) > 1 else 560, 32, 70.0
vmax = np.sqrt(5 / 6)
g = (np.arange(N) + 0.5) / N * vmax
V1, V2 = np.meshgrid(g, g, indexing='ij')
m = V1**2 + V2**2 < 5 / 6 - 1e-3
t0 = time.time()
L, c = grid(V1[m], V2[m], TMAX, K)
LL = np.zeros((N, N, K), np.int8); C = np.full((N, N), -1, np.int32)
LL[m] = L; C[m] = c
np.savez_compressed(os.path.join(os.path.dirname(os.path.abspath(__file__)), f'grid{N}.npz'), g=g, L=LL, C=C, K=K, TMAX=TMAX)
print('done', N, m.sum(), time.time() - t0)
