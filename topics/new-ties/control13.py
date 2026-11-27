"""Round 13 (DRAFT) control: dist(rho, R) for random stretches, uniform in log over both bad windows
(weighted by log length, as Weyl's theorem weights the bad pairs). Usage: python3 control13.py N [WORKERS] -> control13.json"""
import sys, json, random, time
import mpmath as mp
import rdist13 as D
mp.mp.dps = 120
WIN = [(mp.mpf(9) / 10, mp.mpf(10) / 9), (mp.mpf(9) / 19, mp.mpf(100) / 171)]
def one(seed):
    rng = random.Random(seed)
    L = [float(mp.log(b / a)) for a, b in WIN]
    i = 0 if rng.random() < L[0] / sum(L) else 1
    a, b = WIN[i]
    r = a * mp.exp(mp.mpf(rng.random()) * mp.log(b / a))
    lb, ub, xw, yw, n, done = D.dist_R(r)
    return (seed, "NK"[i], float(r), float(mp.log(lb, 2)) if lb > 0 else None, float(mp.log(ub, 2)), n, done)
if __name__ == "__main__":
    from multiprocessing import Pool
    N = int(sys.argv[1]); w = int(sys.argv[2]) if len(sys.argv) > 2 else 2
    t0 = time.time()
    with Pool(w) as P: res = P.map(one, range(1300000, 1300000 + N), chunksize=50)
    json.dump(dict(n=N, not_done=sum(1 for p in res if not p[6]), seconds=round(time.time() - t0, 1), per_point=res), open("control13.json", "w"))
    print("control", N, "closest", min(res, key=lambda p: p[4])[:5])
