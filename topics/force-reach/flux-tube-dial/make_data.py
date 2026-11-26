"""Heavy runs behind the pictures (about 4 minutes): field maps for s = 2, 1.5, 1.25, 1.1 and vortex tensions n = 1..8.
Run: python3 -I make_data.py  ->  data/"""
import sys, os, json, math
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
from field_solver import solve
from vortex_solver import vortex
out = {}
for s in (1.5, 1.25, 1.1):
    Es, widths = [], []
    for L in (2.0, 3.0, 4.0, 5.0):
        e, Dm, x, res = solve(s, L, N=160, box=20.0, w=0.25, eps=1e-4)
        Es.append(e); mid = Dm[len(x) // 2, :]; widths.append(math.sqrt(np.sum(mid * x ** 2) / np.sum(mid)))
        if L == 4.0: np.save(f'data/Dm_s{s}.npy', Dm)
    out[s] = dict(E=Es, widths=widths)
json.dump({str(k): v for k, v in out.items()}, open('data/field2.json', 'w'))
e, Dm, x, res = solve(2.0, 4.0, N=160, box=20.0, w=0.25, eps=1e-4); np.save('data/Dm_s2.0.npy', Dm)
V = {}
for lam in (0.1, 0.25, 0.5, 0.75, 1.0, 1.5, 2.0, 4.0, 8.0):
    V[str(lam)] = [vortex(n, lam, R=40.0 if lam < 0.3 else 30.0)[0] / math.pi for n in range(1, 9)]
json.dump(V, open('data/vortex_E.json', 'w'))
