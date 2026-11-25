"""proof_uses.py: which proof shapes lean on which, read off the compiled library.

Input : gallery.json (every library theorem with its shape family, from the proof
        gallery's own rules, rerun on all 2,915 theorems) and extract.jsonl (Lean's own
        list of the constants each compiled proof uses).
Edge  : T -> U when U's compiled proof uses theorem T ("U needs T"). No hand step.
Output: proof_uses.json and proof_uses_output.txt (family graph, body graph, cycles,
        layers, and a fair null: family labels shuffled within each Lean module).
"""
import json, random, sys, pathlib, collections
HERE = pathlib.Path(__file__).resolve().parent
G = json.load(open(sys.argv[1]))
EX = [json.loads(l) for l in open(sys.argv[2])]
out_lines = []
def say(*a):
    s = " ".join(str(x) for x in a); print(s); out_lines.append(s)

fam = {p["name"]: p["cluster"] for p in G["proofs"]}
bod = {p["name"]: p["body"] for p in G["proofs"]}
mod = {p["name"]: p["module"] for p in G["proofs"]}
names = set(fam)
edges = []
for x in EX:
    u = x["name"]
    if u not in names: continue
    for c in x["consts"]:
        if c in names and c != u:
            edges.append((c, u))
edges = sorted(set(edges))
say("theorems with a shape:", len(names), " uses (T -> U, U's proof uses T):", len(edges))
say("theorems whose proof uses another library theorem:", len({u for _, u in edges}))

BODY_ORDER = ["point", "rod", "rect", "box", "hyper", "cube"]
def rank(b): return BODY_ORDER.index(b) if b in BODY_ORDER else -1

def body_flow(lab):
    up = down = same = 0
    for t, u in edges:
        a, b = rank(lab[t]), rank(lab[u])
        if a < b: up += 1
        elif a > b: down += 1
        else: same += 1
    return up, down, same

up, down, same = body_flow(bod)
say(f"\nBody flow (needed theorem's body -> using theorem's body): up {up}, down {down}, same {same}")
# fair null: shuffle bodies (and families) among theorems of the same module
by_mod = collections.defaultdict(list)
for n in names: by_mod[mod[n]].append(n)
random.seed(1)
nulls = []
for _ in range(2000):
    lab = {}
    for m, ns in by_mod.items():
        vals = [bod[n] for n in ns]; random.shuffle(vals)
        lab.update(zip(ns, vals))
    u2, d2, s2 = body_flow(lab)
    nulls.append(u2 - d2)
mu = sum(nulls) / len(nulls); sd = (sum((v - mu) ** 2 for v in nulls) / len(nulls)) ** 0.5
obs = up - down
p_hi = sum(1 for v in nulls if v >= obs) / len(nulls)
say(f"net upward flow {obs}; null (bodies shuffled within module) mean {mu:.1f}, sd {sd:.1f}, z {(obs-mu)/sd:.2f}, share of nulls at least as high {p_hi:.4f}")

# body x body matrix
BM = collections.Counter((bod[t], bod[u]) for t, u in edges)
say("\nBody matrix (row = needed, column = user):")
say("        " + " ".join(f"{b:>6}" for b in BODY_ORDER))
for a in BODY_ORDER:
    say(f"{a:>7} " + " ".join(f"{BM[(a,b)]:>6}" for b in BODY_ORDER))

# family graph
FM = collections.Counter((fam[t], fam[u]) for t, u in edges)
fams = sorted(set(fam.values()))
size = collections.Counter(fam.values())
# coverage: share of family-Y theorems that use at least one family-X theorem
cover = collections.defaultdict(set)
for t, u in edges: cover[(fam[t], fam[u])].add(u)
users_in = collections.Counter()
for u in {u for _, u in edges}: users_in[fam[u]] += 1

def sccs(nodes, adj):
    idx, low, st, on, res, c = {}, {}, [], set(), [], [0]
    def strong(v):
        idx[v] = low[v] = c[0]; c[0] += 1; st.append(v); on.add(v)
        for w in adj.get(v, ()):
            if w not in idx: strong(w); low[v] = min(low[v], low[w])
            elif w in on: low[v] = min(low[v], idx[w])
        if low[v] == idx[v]:
            comp = []
            while True:
                w = st.pop(); on.discard(w); comp.append(w)
                if w == v: break
            res.append(comp)
    sys.setrecursionlimit(10000)
    for v in nodes:
        if v not in idx: strong(v)
    return res

results = {"edges": len(edges), "body_flow": dict(up=up, down=down, same=same, null_mean=mu, null_sd=sd, p=p_hi),
           "body_matrix": {f"{a}>{b}": BM[(a, b)] for a in BODY_ORDER for b in BODY_ORDER},
           "family_edges": [], "thresholds": []}
for (x, y), w in sorted(FM.items(), key=lambda z: -z[1]):
    results["family_edges"].append(dict(src=x, dst=y, w=w, users=len(cover[(x, y)]), dst_size=size[y]))

say("\nFamily graph at rising thresholds (edge X -> Y kept when at least w uses):")
for th in (1, 2, 3, 5, 8, 12, 20):
    adj = collections.defaultdict(set)
    for (x, y), w in FM.items():
        if w >= th and x != y: adj[x].add(y)
    comps = [c for c in sccs(fams, adj) if len(c) > 1]
    two_way = sorted({tuple(sorted((x, y))) for x in adj for y in adj[x] if x in adj.get(y, ())})
    say(f"  w >= {th:>2}: edges {sum(len(v) for v in adj.values()):>3}, cycles (strong components > 1): {[sorted(c) for c in comps]}")
    say(f"           two-way pairs: {two_way}")
    results["thresholds"].append(dict(w=th, edges=sum(len(v) for v in adj.values()),
                                      components=[sorted(c) for c in comps], two_way=two_way))

say("\nStrongest family uses (needed -> user: uses; users in target / target size):")
for e in results["family_edges"][:40]:
    say(f"  {e['src']:>14} -> {e['dst']:<14} {e['w']:>4}   {e['users']}/{e['dst_size']}")

# does every theorem of family Y that needs anything need some family-X theorem?
say("\nFull coverage claims (every proof in Y that uses any theorem uses one from X), Y with >= 5 users:")
for y in fams:
    if users_in[y] < 5: continue
    for x in fams:
        if x == y: continue
        if len(cover[(x, y)]) == users_in[y]:
            say(f"  every using proof in {y} ({users_in[y]}) uses {x}")

# self-use share
self_w = sum(w for (x, y), w in FM.items() if x == y)
say(f"\nUses inside one family: {self_w} of {len(edges)} ({100*self_w/len(edges):.0f}%)")
# action flow
act = {p["name"]: p["lead"] for p in G["proofs"]}
AM = collections.Counter((act[t], act[u]) for t, u in edges)
acts = sorted(set(act.values()))
say("\nAction matrix (row = needed lead move, column = user):")
say("          " + " ".join(f"{a[:7]:>8}" for a in acts))
for a in acts:
    say(f"{a[:9]:>9} " + " ".join(f"{AM[(a,b)]:>8}" for b in acts))
results["action_matrix"] = {f"{a}>{b}": AM[(a, b)] for a in acts for b in acts}
results["family_sizes"] = dict(size)
results["users_in"] = dict(users_in)
json.dump(results, open(HERE / "proof_uses.json", "w"), indent=1)
open(HERE / "proof_uses_output.txt", "w").write("\n".join(out_lines) + "\n")
