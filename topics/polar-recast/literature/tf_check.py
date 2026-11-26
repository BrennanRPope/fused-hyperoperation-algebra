"""Test the 'tripartite form' (Hubeny & Liu 2026, arXiv:2609.22504) on all 5-region rules,
in corner-weight pictures. A conditional tripartite information -I3(A:B:C|X) is, in corner
weights, minus the sum of the corners in the upward interval [ABC, ABC+X]. So a rule is in
tripartite form iff its corner group (with multiplicity) splits into upward intervals that each
start at a 3-region corner. Exact search, regions A..E (purifier removed by complement)."""
import itertools, json, re, sys
from collections import Counter
sys.path.insert(0, "../hec5")
src = open("../hec5/hec5_signs.py").read().split("report = {}")[0]
exec(src)
REG = "ABCDE"
def tilings(m, limit=1):
    """m: Counter corner(frozenset)->multiplicity. Return up to `limit` tilings."""
    m = +m
    if not m: return [[]]
    J = min(m, key=lambda s: (len(s), sorted(s)))
    if len(J) != 3: return []
    rest = [r for r in REG if r not in J]
    found = []
    for k in range(len(rest)+1):
        for X in itertools.combinations(rest, k):
            iv = [J | frozenset(Y) for kk in range(k+1) for Y in itertools.combinations(X, kk)]
            if all(m[s] > 0 for s in iv):
                m2 = m.copy()
                for s in iv: m2[s] -= 1
                for t in tilings(m2, limit-len(found)):
                    found.append([("".join(sorted(J)), "".join(X))] + t)
                    if len(found) >= limit: return found
    return found
res = {}; examples = {}
for k, (L, R) in REPS.items():
    if k == 1: continue
    base = parse(L, 1) + Counter(); base.update(parse(R, -1))
    insts = set()
    for p in itertools.permutations(SYM):
        mp = dict(zip(SYM, p))
        insts.add(canon(Counter({frozenset(mp[x] for x in I): c for I, c in base.items()})))
    ok = 0
    for inst in insts:
        cc = corner_coeffs(inst)
        m = Counter({frozenset(J): -v for J, v in cc.items()})
        t = tilings(m)
        if t:
            ok += 1; examples.setdefault(k, {"corners": cc, "tiling": t[0]})
    res[NAME[k]] = (ok, len(insts))
    print(f"{NAME[k]:11s}: {ok}/{len(insts)} split into upward intervals from 3-region corners")
print("total", sum(a for a, b in res.values()), "/", sum(b for a, b in res.values()))
for k, e in examples.items(): print(NAME[k], "e.g.", e["tiling"])
json.dump({"result": res, "examples": {NAME[k]: e for k, e in examples.items()}},
          open("tf_check.json", "w"), indent=1)
