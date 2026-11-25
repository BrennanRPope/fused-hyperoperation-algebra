"""orders.py: the best possible order of bodies and of actions, and the exceptions.
For each order we count back edges (uses that run from a later shape to an earlier one).
Exact search: all 720 orders of the six bodies; dynamic programming over subsets for the
eleven lead moves. Null: labels shuffled within each Lean module (2,000 draws), best order
recomputed each time. Output: orders.json, orders_output.txt."""
import json, sys, random, itertools, collections, pathlib
HERE = pathlib.Path(__file__).resolve().parent
G = json.load(open(sys.argv[1])); EX = [json.loads(l) for l in open(sys.argv[2])]
P = {p["name"]: p for p in G["proofs"]}
edges = sorted({(c, x["name"]) for x in EX if x["name"] in P for c in x["consts"] if c in P and c != x["name"]})
lines = []
def say(*a):
    s = " ".join(map(str, a)); print(s); lines.append(s)

def best_order(lab, keys):
    W = collections.Counter((lab[t], lab[u]) for t, u in edges if lab[t] != lab[u])
    n = len(keys); idx = {k: i for i, k in enumerate(keys)}
    INF = 10**9; dp = [INF] * (1 << n); par = [None] * (1 << n); dp[0] = 0
    for S in range(1 << n):
        if dp[S] == INF: continue
        for j in range(n):
            if S >> j & 1: continue
            # placing j after S: back edges = uses from j into members of S
            cost = sum(W[(keys[j], keys[i])] for i in range(n) if S >> i & 1)
            T = S | 1 << j
            if dp[S] + cost < dp[T]: dp[T] = dp[S] + cost; par[T] = (S, j)
    order, S = [], (1 << n) - 1
    while S: S, j = par[S]; order.append(keys[j])
    total = sum(W.values())
    return order[::-1], dp[(1 << n) - 1], total

res = {}
for kind, field in (("body", "body"), ("action", "lead")):
    lab = {n: P[n][field] for n in P}
    keys = sorted(set(lab.values()))
    order, back, total = best_order(lab, keys)
    by_mod = collections.defaultdict(list)
    for n in P: by_mod[P[n]["module"]].append(n)
    random.seed(7); nulls = []
    for _ in range(300):
        l2 = {}
        for m, ns in by_mod.items():
            v = [lab[n] for n in ns]; random.shuffle(v); l2.update(zip(ns, v))
        o2, b2, t2 = best_order(l2, keys); nulls.append(b2 / t2)
    share = back / total
    mu = sum(nulls) / len(nulls); sd = (sum((x - mu) ** 2 for x in nulls) / len(nulls)) ** .5
    say(f"\n{kind}: best order {' -> '.join(order)}")
    say(f"  cross-shape uses {total}, against the order {back} ({100*share:.1f}%);"
        f" null best-order back share mean {100*mu:.1f}%, sd {100*sd:.1f}%, z {(share-mu)/sd:.2f},"
        f" draws at or below observed {sum(1 for x in nulls if x <= share)}/{len(nulls)}")
    W = collections.Counter((lab[t], lab[u]) for t, u in edges if lab[t] != lab[u])
    pairs = []
    for i, a in enumerate(order):
        for b in order[i+1:]:
            pairs.append((a, b, W[(a, b)], W[(b, a)]))
    say("  pair by pair (earlier -> later uses, later -> earlier uses):")
    for a, b, f, r in pairs:
        if f or r: say(f"    {a:>8} -> {b:<8} {f:>4}  back {r:>4}" + ("   <- reversed" if r > f else ""))
    res[kind] = dict(order=order, back=back, total=total, null_mean=mu, null_sd=sd,
                     pairs=[dict(a=a, b=b, fwd=f, back=r) for a, b, f, r in pairs])

# exceptions: cube -> rod (any-size law used by a fixed-size one), hyper -> box
def show(a, b, k=8):
    ex = [(t, u) for t, u in edges if P[t]["body"] == a and P[u]["body"] == b]
    say(f"\nExamples, {a} needed by {b} ({len(ex)}):")
    for t, u in ex[:k]:
        say(f"  {P[u]['short']} [{P[u]['module'].split('.')[-1]}] uses {P[t]['short']}: {P[t]['statement'][:90]}")
    return [dict(needed=t, user=u, stmt=P[t]["statement"], ustmt=P[u]["statement"]) for t, u in ex]
res["ex_cube_rod"] = show("cube", "rod", 12)
res["ex_hyper_box"] = show("hyper", "box", 6)
res["ex_rod_cube"] = show("rod", "cube", 6)
json.dump(res, open(HERE / "orders.json", "w"), indent=1)
open(HERE / "orders_output.txt", "w").write("\n".join(lines) + "\n")
