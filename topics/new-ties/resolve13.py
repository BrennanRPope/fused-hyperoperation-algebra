"""Round 13 (DRAFT): finish the safe-radius run for the pairs safe13.py left open at k <= 140.
Same method as safe13.py with finer centring (Q = 340), k up to 320, budget 4000 (run 1) then 8000 for what is still open (run 2: python3 resolve13.py 2 8000). Updates safe13.json in place
(adds key 'resolved_later' listing the pairs finished here)."""
import json, sys
from fractions import Fraction as F
import mpmath as mp
import safe13 as S
mp.mp.dps = 140
S.LL, S.L20 = mp.log(mp.mpf(19) / 200), mp.log(20)
S.KMAX, S.Q = 320, 340
S.BUDGET = int(sys.argv[2]) if len(sys.argv) > 2 else 4000
def safe(ab):
    a, b = ab
    r = mp.exp(a * S.LL + b * S.L20); c = F(int(mp.nint(r * 2 ** S.Q)), 2 ** S.Q); tot = 0
    for k in range(4, S.KMAX + 1):
        st, n, hits = S.search_nomemo(c - F(1, 2 ** k), c + F(1, 2 ** k), maxnodes=S.BUDGET); tot += n
        if st == "ok" and not hits: return [a, b, float(r), k, n, tot]
    return [a, b, float(r), None, None, tot]
if __name__ == "__main__":
    from multiprocessing import Pool
    d = json.load(open("safe13.json"))
    todo = [(p[0], p[1]) for p in d["per_pair"] if p[3] is None]
    with Pool(int(sys.argv[1]) if len(sys.argv) > 1 else 2) as P: res = P.map(safe, todo, chunksize=1)
    rm = {(p[0], p[1]): p for p in res}
    d["per_pair"] = [rm.get((p[0], p[1]), p) for p in d["per_pair"]]
    d.setdefault("resolved_later", []); d["resolved_later"] = [q for q in d["resolved_later"] if q[3] is not None] + [dict(budget=S.BUDGET, pair=q) for q in res]; d["unresolved"] = sum(1 for p in d["per_pair"] if p[3] is None)
    d["deepest"] = sorted([p for p in d["per_pair"] if p[3] is not None], key=lambda p: -p[3])[:20]
    json.dump(d, open("safe13.json", "w")); print(res, "unresolved now", d["unresolved"])
