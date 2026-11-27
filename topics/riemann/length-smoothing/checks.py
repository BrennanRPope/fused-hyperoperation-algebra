"""Checks for length-smoothing. Run: python3 checks.py. DRAFT."""
import numpy as np, json
d = json.load(open("windows.json")); R = d["rows"]; ref = json.load(open("../knauf-field/refined.json")); ok = []
def check(name, cond): ok.append(bool(cond)); print(("PASS " if cond else "FAIL ") + name)
def rms(kind, k, N, j=0): return np.sqrt(np.mean([r[j]**2 for r in R["%s k=%d N=%d" % (kind, k, N)]]))
Ns = range(16, 25); m = lambda kind, k, j=0: np.mean([rms(kind, k, N, j) for N in Ns])
check("single-length zeros match knauf-field's refined zeros (N = 20)", np.allclose([r[0] for r in R["tri k=0 N=20"]], [r[0]-1.5 for r in ref["rows"]["20"]], atol=1e-8))
check("triangle window of 7 or 9 lengths cuts the mean sideways wobble at least 3x", min(m("tri", 6), m("tri", 8)) < m("tri", 0)/3)
check("it also cuts the miss in height at least 3x", m("tri", 8, 1) < m("tri", 0, 1)/3)
check("windows of 7 and 9 lengths beat no averaging at every N = 16..24", all(rms("tri", k, N) < rms("tri", 0, N) for k in (6, 8) for N in Ns))
check("U shape: the widest window (13) is worse than 9 on average", m("tri", 12) > m("tri", 8) and m("flat", 12) > m("flat", 8))
check("the widest window improves steadily as N grows (N = 16 to 22)", all(rms("tri", 12, N) > rms("tri", 12, N+1) for N in range(16, 22)))
check("growth per length of Z_n(3/2) is about 1.153, not sqrt 2", abs(d["growth"] - 1.153) < 0.002)
check("each single-chain zero lands on both sides of the line across n = 10..24", all(len(set(np.sign([ref["rows"][str(n)][j][0]-1.5 for n in range(10, 25)]))) == 2 for j in range(5)))
# --- follow-up: long chains and the swing law ---
import sys; sys.path.insert(0, "../knauf-field"); from knauf import levels
import os
if os.path.exists("/tmp/knauf_levels/levels_22.npz"):
    dd = np.load("/tmp/knauf_levels/levels_22.npz"); q, cc = levels(22)
    check("streamed level counts match the direct chain at length 22", np.array_equal(dd["q"], q.astype(np.int64)) and np.array_equal(dd["c"], cc.astype(np.int64)))
w = json.load(open("wide.json")); s = json.load(open("swing.json"))
rw = lambda name, N: np.sqrt(np.mean([r[0]**2 for r in w["rows"]["%s N=%d" % (name, N)]]))
check("lengths 26..33: growing window (half the length) keeps the wobble below 0.002 at every length", all(rw("k=N/2", N) < 0.002 for N in range(26, 34)))
check("lengths 26..33: single chains stay above 0.005", all(rw("single", N) > 0.005 for N in range(26, 34)))
check("averaging cuts the wobble at least 5x over lengths 26..33", np.mean([rw("k=N/2", N) for N in range(26, 34)]) < np.mean([rw("single", N) for N in range(26, 34)])/5)
check("swing law: turn per length grows with height, straight-line slope about 0.05 per unit of height", -0.06 < s["fit"][0] < -0.045)
check("swing is steady (circular consistency >= 0.75) for the first five zeros except one", sum(v >= 0.75 for v in s["turn_spread"][:5]) >= 4)
print("%d/%d passed" % (sum(ok), len(ok)))
