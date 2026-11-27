"""Round 13 (DRAFT): how close does each bad stretch come to the line-up set R?

For a bad pair (a, b) with stretch rho = lam^a / mu^b, its safe radius is 2^-k for the smallest k (4 <= k <= KMAX = 140)
such that the range search of cover12.py, run WITHOUT the visited-state memo, finishes with no hit on
[c - 2^-k, c + 2^-k], where c is rho rounded to a multiple of 2^-160. Without the memo, a point of R in the range would
give an infinite uncut branch (Theorem 12.2), so a finished search proves the range holds no point of R:
dist(rho, R) >= 2^-k - 2^-161.

Usage: python3 safe13.py AMAX [WORKERS] -> safe13.json   (all bad pairs with 1 <= a <= AMAX)
"""
from fractions import Fraction as F
import sys, json, time, inspect
import mpmath as mp
import cover12 as C

mp.mp.dps = 80
KMAX = 140; BUDGET = 2000; Q = 160
src = inspect.getsource(C.search).replace("        if node in seen: continue\n        seen.add(node); nodes += 1",
                                           "        nodes += 1").replace("def search(", "def search_nomemo(")
assert "seen.add" not in src
exec(src, C.__dict__)
search_nomemo = C.search_nomemo

LL, L20 = mp.log(mp.mpf(19) / 200), mp.log(20)
WIN = [(mp.mpf(9) / 10, mp.mpf(10) / 9), (mp.mpf(9) / 19, mp.mpf(100) / 171)]
def bad_pairs(a):
    out = []
    for w1, w2 in WIN:
        lo = (mp.log(w1) - a * LL) / L20; hi = (mp.log(w2) - a * LL) / L20
        for b in range(int(mp.ceil(lo)), int(mp.floor(hi)) + 1):
            if b >= 1: out.append(b)
    return out

def safe(ab):
    a, b = ab
    r = mp.exp(a * LL + b * L20)
    c = F(int(mp.nint(r * 2 ** Q)), 2 ** Q)
    tot = 0
    for k in range(4, KMAX + 1):
        st, n, hits = search_nomemo(c - F(1, 2 ** k), c + F(1, 2 ** k), maxnodes=BUDGET)
        tot += n
        if st == "ok" and not hits: return (a, b, float(r), k, n, tot)
    return (a, b, float(r), None, None, tot)

if __name__ == "__main__":
    from multiprocessing import Pool
    amax = int(sys.argv[1]); workers = int(sys.argv[2]) if len(sys.argv) > 2 else 4
    pairs = [(a, b) for a in range(1, amax + 1) for b in bad_pairs(a)]
    t0 = time.time()
    with Pool(workers) as P: res = P.map(safe, pairs, chunksize=20)
    none = [p for p in res if p[3] is None]
    out = dict(amax=amax, kmax=KMAX, budget=BUDGET, pairs=len(res), unresolved=len(none),
               deepest=sorted([p for p in res if p[3] is not None], key=lambda p: -p[3])[:20],
               seconds=round(time.time() - t0, 1), per_pair=res)
    json.dump(out, open("safe13.json", "w"))
    print("pairs", len(res), "unresolved", len(none), "max k", max(p[3] for p in res if p[3]), "seconds", out["seconds"])
    for p in out["deepest"][:10]: print("  ", p)
