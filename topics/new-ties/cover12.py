"""Round 12 (DRAFT): prove "no second collision" for whole ranges of the stretch at once.

For a range [r1, r2] of the stretch rho, run the Round 11 peel-and-prune search with rho kept as a
parameter: the map is A(rho) = (rho*S, rho*T + U), every peel is exact in (S, T, U), and a branch is
cut only when it is dead for EVERY rho in the range (the tests are on convex/concave/linear functions
of rho, so checking the two ends is exact). If every branch is cut and no branch reaches the identity,
no stretch in the range has a solution, in particular no bad pair with its stretch in the range.

Usage: python3 cover12.py BUDGET LOG2_MIN_WIDTH [WORKERS]   -> cover12.json
Bisects each window until every piece is proved or narrower than 2^-LOG2_MIN_WIDTH.
"""
from fractions import Fraction as F
import sys, json, time, math

LAM, MU = F(19, 200), F(1, 20)
LET = {"l": (LAM, F(0)), "m": (MU, F(0)), "N": (F(1, 19), F(1, 19)), "K": (F(1, 10), F(1, 10))}
REDUCED = {"none": "lmNKE", "lam": "lmNK", "mu": "mN"}
UNREDUCED = {"none": "lmNKE", "lam": "lmNK", "mu": "mNK"}
NINTH = F(1, 9)
WINDOWS = {"N": (F(9, 10), F(10, 9)), "K": (F(9, 19), F(100, 171))}

def nxt(kind, x):
    if x in "NK": return "none"
    return "lam" if x == "l" else "mu"

def search(r1, r2, hforced="N", gforced="NK", maxnodes=20000, start=(F(1), F(0), F(0)), allowed=REDUCED):
    """returns (status, nodes, hits): status 'ok' = every branch cut; 'big' = budget exceeded.
    hits = list of rho in [r1, r2] where both words end with A = identity (a real solution)."""
    rm = (r1 + r2) / 2
    S0, T0, U0 = start
    stack = [(S0, T0, U0, "none", "none", False, False, hforced, gforced)]
    seen = set(); hits = []; nodes = 0
    while stack:
        node = stack.pop()
        if node in seen: continue
        seen.add(node); nodes += 1
        if nodes > maxnodes: return ("big", nodes, hits)
        S, T, U, hs, gs, he, ge, hf, gf = node
        if he and ge:
            r0 = 1 / S
            if r1 <= r0 <= r2 and r0 * T + U == 0: hits.append(r0)
            continue
        f0 = (r1 * T + U, r2 * T + U); f1 = (r1 * (S / 9 + T) + U, r2 * (S / 9 + T) + U)
        hi = (max(f0[0], f1[0]), max(f0[1], f1[1])); lo = (min(f0[0], f1[0]), min(f0[1], f1[1]))
        if he:   # h used up: A(I) must cover I somewhere in the range
            if (lo[0] > 0 and lo[1] > 0) or (hi[0] < NINTH and hi[1] < NINTH): continue
        elif ge: # g used up: A(I) must sit inside I somewhere in the range
            if (f0[0] < 0 and f0[1] < 0) or (f1[0] < 0 and f1[1] < 0) or \
               (f0[0] > NINTH and f0[1] > NINTH) or (f1[0] > NINTH and f1[1] > NINTH): continue
        else:    # A(I) must meet I somewhere in the range
            if (hi[0] < 0 and hi[1] < 0) or (lo[0] > NINTH and lo[1] > NINTH): continue
        side = "g" if (he or (not ge and abs(rm * S) >= 1)) else "h"
        kind = gs if side == "g" else hs; forced = gf if side == "g" else hf
        o = allowed[kind]
        for x in (forced if forced else o):
            if x not in o: continue
            if side == "g":
                if x == "E": stack.append((S, T, U, hs, gs, he, True, hf, None))
                else:
                    sy, cy = LET[x]
                    stack.append((S * sy, S * cy + T, U, hs, nxt(gs, x), he, ge, hf, None))
            else:
                if x == "E": stack.append((S, T, U, hs, gs, True, ge, None, gf))
                else:
                    sx, cx = LET[x]
                    stack.append((S / sx, T / sx, (U - cx) / sx, nxt(hs, x), gs, he, ge, None, gf))
    return ("ok", nodes, hits)

def cover(lo, hi, budget, minw, out):
    st, n, hits = search(lo, hi, maxnodes=budget)
    if st == "ok":
        out.append((lo, hi, n, [str(h) for h in hits])); return
    if hi - lo < minw:
        out.append((lo, hi, None, [])); return
    mid = (lo + hi) / 2
    cover(lo, mid, budget, minw, out); cover(mid, hi, budget, minw, out)

def job(args):
    lo, hi, budget, minw = args
    out = []; cover(lo, hi, budget, minw, out); return out

def logmeas(pieces): return sum(math.log(float(b) / float(a)) for a, b in pieces)

if __name__ == "__main__":
    from multiprocessing import Pool
    budget = int(sys.argv[1]); k = int(sys.argv[2]); workers = int(sys.argv[3]) if len(sys.argv) > 3 else 4
    minw = F(1, 2 ** k); CH = 64
    t0 = time.time(); res = {"budget": budget, "min_width": "2^-%d" % k, "windows": {}}
    tasks = []
    for name, (a, b) in WINDOWS.items():
        for i in range(CH): tasks.append((name, a + (b - a) * F(i, CH), a + (b - a) * F(i + 1, CH)))
    with Pool(workers) as P:
        outs = P.map(job, [(lo, hi, budget, minw) for (_, lo, hi) in tasks], chunksize=1)
    for name, (a, b) in WINDOWS.items():
        pieces = [p for (nm, _, _), o in zip(tasks, outs) if nm == name for p in o]
        ok = [(p[0], p[1]) for p in pieces if p[2] is not None]
        opn = []
        for p in pieces:
            if p[2] is None:
                if opn and opn[-1][1] == p[0]: opn[-1][1] = p[1]
                else: opn.append([p[0], p[1]])
        hits = [h for p in pieces for h in p[3]]
        res["windows"][name] = dict(
            window=[str(a), str(b)], pieces=len(pieces), proved=len(ok), hits=hits,
            proved_fraction=float(sum(q - p for p, q in ok) / (b - a)),
            proved_log_fraction=logmeas(ok) / math.log(float(b) / float(a)),
            max_nodes=max(p[2] for p in pieces if p[2] is not None),
            open=[[str(p), str(q)] for p, q in opn],
            proved_pieces=[[str(p[0]), str(p[1]), p[2]] for p in pieces if p[2] is not None])
        w = res["windows"][name]
        print(name, "proved %.4f (log %.4f)" % (w["proved_fraction"], w["proved_log_fraction"]),
              "pieces", len(pieces), "open intervals", len(opn), "hits", len(hits), flush=True)
    res["seconds"] = round(time.time() - t0, 1)
    json.dump(res, open("cover12.json", "w"), indent=0)
    print("time", res["seconds"])
