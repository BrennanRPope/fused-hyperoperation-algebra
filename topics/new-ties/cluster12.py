"""Round 12 (DRAFT): how much of the cluster 0.944-0.963 stays unproved at a given halving depth (same search as
cover12.py). A bigger node budget does not help (2,000 and 8,000 nodes give the same pieces); halving finer does.
Usage: python3 cluster12.py BUDGET LOG2_MIN_WIDTH  -> one line in cluster12_K.txt"""
from fractions import Fraction as F
import sys, time
import cover12 as C
budget = int(sys.argv[1]); k = int(sys.argv[2]); lo, hi = F(944, 1000), F(963, 1000); out = []; t = time.time()
C.cover(lo, hi, budget, F(1, 2 ** k), out)
openw = sum(p[1] - p[0] for p in out if p[2] is None)
line = "k %d budget %d open %.10f of %.3f pieces %d hits %d seconds %.0f" % (
    k, budget, float(openw), float(hi - lo), len(out), sum(len(p[3]) for p in out), time.time() - t)
open("cluster12_%d.txt" % k, "w").write(line + "\n"); print(line, flush=True)
