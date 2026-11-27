"""Round 13 (DRAFT): dist(rho, R) for every bad pair with a <= AMAX, by rdist13.dist_R (two-sided, factor 2).
Usage: python3 dist13.py AMAX [WORKERS] -> dist13.json"""
import sys, json, time
import mpmath as mp
import rdist13 as D
import safe13 as S
mp.mp.dps = 120
LL, L20 = mp.log(mp.mpf(19) / 200), mp.log(20)   # recomputed at 120 digits
def one(ab):
    a, b = ab
    r = mp.exp(a * LL + b * L20)
    lb, ub, xw, yw, n, done = D.dist_R(r)
    return (a, b, float(r), float(mp.log(lb, 2)) if lb > 0 else None, float(mp.log(ub, 2)), n, done, xw, yw)
if __name__ == "__main__":
    from multiprocessing import Pool
    amax = int(sys.argv[1]); w = int(sys.argv[2]) if len(sys.argv) > 2 else 4
    pairs = [(a, b) for a in range(1, amax + 1) for b in S.bad_pairs(a)]
    t0 = time.time()
    with Pool(w) as P: res = P.map(one, pairs, chunksize=50)
    out = dict(amax=amax, pairs=len(res), not_done=sum(1 for p in res if not p[6]),
               seconds=round(time.time() - t0, 1), per_pair=res)
    json.dump(out, open("dist13.json", "w"))
    m = min(res, key=lambda p: p[4])
    print("pairs", len(res), "not_done", out["not_done"], "closest", m[:6], "seconds", out["seconds"])
