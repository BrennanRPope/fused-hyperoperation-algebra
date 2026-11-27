import sys; sys.path.insert(0, '..')
import sympy as sp, itertools
from misalign import W_stu
from precess import T_series
from turn2 import turn2
one = sp.Integer(1)
rows = []
bgs = {"1": [one, 0, 0, 0], "2": [one, one, 0, 0], "2u": [one, 2*one, 0, 0], "3": [one, one, one, 0], "4": [one, one, one, one], "4u": [one, 2*one, 3*one, 5*one]}
probes = {"anti1": ([one, 0, 0, 0], [-1, 1, 1, 1]), "anti12": ([one/2, one/2, 0, 0], [-1, -1, 1, 1]), "anti_all": ([one/4]*4, [-1]*4),
          "anti1_al2": ([one/2, one/2, 0, 0], [-1, 1, 1, 1]), "anti1_x3": ([one/2, 0, one/2, 0], [-1, 1, 1, 1])}
for bn, ks in bgs.items():
    for pn, (ws, ss) in probes.items():
        if not any(ss[j] == -1 and ws[j] != 0 and ks[j] != 0 for j in range(4)): continue
        W = W_stu(ks, ws, ss); T, K = T_series(W); T1, T2, T2e = turn2(W, K)
        print(f"bg {bn:3s} probe {pn:10s} K={str(K):5s} T={str(T1):6s} T2={str(T2):10s} T2e={str(T2e):8s}")
