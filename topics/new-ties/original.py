"""Round 11 (DRAFT): the original break, second collisions ruled out run by run (README section 15).

The original break: lam = 19/200, mu = 1/20 (about 0), nu = 1/19 (letter N, about 1), kappa = 1/10 (letter K, about 2).
Round 10 (Theorem P12, Steps 1-2) reduced "the reduced pieces S' never collide" to one shape of equation:
    after cancelling the common beginning,  mu^b h' = lam^a g'   as maps,   a, b >= 1,
where h' and g' are reduced words (the rest of the two words), h' starts with N (the mu-side's run contains mu),
g' starts with N or K. Steps 2-4 dispose of every (a, b) whose ratio rho = lam^a / mu^b misses the two bad windows
[9/10, 10/9] and [9/19, 100/171]. Round 10 had no argument for the "bad pairs" (first one: a = 13, b = 10).

This script settles every bad pair with a <= AMAX exactly, by a finite descent (no floating point anywhere):
  y-model: y = z/18 sends [0, 2] to I = [0, 1/9]; the letters become
      l: y -> lam y,  m: y -> mu y,  N: y -> (y + 1)/19,  K: y -> (y + 1)/10     (so N = (10/19) K as functions).
  To solve h = A o g for reduced words h, g, peel letters off the front: h = x h1 gives h1 = x^-1 o A o g; g = y g1
  gives h = (A o y) o g1. Every solution is a path in this tree. Prune a node when A(I) misses I (then h(I) = A(g(I))
  cannot lie in I), when h is used up and A(I) does not contain I, or when g is used up and A(I) is not inside I.
  The bookkeeping (run state: lam before mu inside a run, a mu-run is never followed by K, a word cannot end inside a
  run) makes the tree walk exactly the reduced words. Nodes are memoised, so the walk ends; "no-collision" means
  every branch was pruned, which proves that this (a, b) has no solution.

  O1  the y-model letters and N = (10/19) K, exactly
  O2  every bad pair (a, b) with a <= AMAX (exact window tests)
  O3  the descent for every pair with 0.3 < rho < 3.5 and a <= AMAX (a superset of the bad pairs)
  O4  sanity: the descent finds collisions when they exist (planted random collisions, and the known first collision
      of the unreduced language)
  O5  the certified bracket: the reduced pieces whose runs have at most AMAX copies of lam have no exact overlaps, so
      Hochman gives a lower bound; the upper bound is s_F (Round 9/10). Interval arithmetic at high precision.
Writes original.json."""
import json, sys, time, random
from fractions import Fraction as F

AMAX = int(sys.argv[1]) if len(sys.argv) > 1 else 10000
LAM, MU, NU, KAP = F(19, 200), F(1, 20), F(1, 19), F(1, 10)
I0, I1 = F(0), F(1, 9)
LET = {"l": (LAM, F(0)), "m": (MU, F(0)), "N": (F(1, 19), F(1, 19)), "K": (F(1, 10), F(1, 10))}
ID = (F(1), F(0))

def comp(f, g): return (f[0] * g[0], f[0] * g[1] + f[1])        # f o g, maps y -> s y + t as (s, t)
def inv(f): return (1 / f[0], -f[1] / f[0])
def img(A):
    a, b = A[0] * I0 + A[1], A[0] * I1 + A[1]
    return (min(a, b), max(a, b))

# run states: "none" (after an N/K letter, or at the start of h', g'), "lam" (inside a run, only lam so far), "mu"
# (inside a run that contains mu). E = the word ends here.
REDUCED = {"none": "lmNKE", "lam": "lmNK", "mu": "mN"}
UNREDUCED = {"none": "lmNKE", "lam": "lmNK", "mu": "mNK"}
def nxt(state, x):
    """state is (kind, run length so far)"""
    if x in "NK": return ("none", 0)
    return ("lam" if x == "l" else "mu", state[1] + 1)

def opts(state, allowed, maxrun):
    o = allowed[state[0]]
    if maxrun is not None and state[1] >= maxrun: o = o.replace("l", "").replace("m", "")
    return o

def descend(A0, hforced=None, gforced=None, allowed=REDUCED, hstate=("none", 0), gstate=("none", 0),
            maxnodes=2_000_000, maxdepth=1000, maxrun=None, depth_prune=False):
    """Search for words h, g in the language `allowed` (runs of at most maxrun class-1 letters, if given) with
    h = A0 o g. Returns (status, nodes, depth, witness). "no-collision" is a proof; "inconclusive" proves nothing.
    depth_prune=True searches only words with (letters of h) + (letters of g) + 2 <= maxdepth, completely (used for the
    planted-collision test, where the tree can also wander down infinite branches: two infinite words with the same
    limit point, e.g. ...N N N against ...K K K when A sends 1/9 to 1/18)."""
    stack = [(A0, hstate, gstate, False, False, hforced, gforced, 0, "", "")]
    seen = set(); nodes = 0; deep = 0
    while stack:
        A, hs, gs, hend, gend, hf, gf, d, hw, gw = stack.pop()
        key = (A, hs, gs, hend, gend, hf, gf, d if depth_prune else 0)
        if key in seen: continue
        seen.add(key); nodes += 1; deep = max(deep, d)
        if depth_prune and d > maxdepth: continue
        if nodes > maxnodes or d > maxdepth: return ("inconclusive", nodes, deep, None)
        if A == ID and not (hend and gend) and hf is None and gf is None and hs[0] == gs[0] == "none":
            return ("collision", nodes, deep, (hw, gw))        # equal tails can be appended to both
        if hend and gend:
            if A == ID: return ("collision", nodes, deep, (hw, gw))
            continue
        lo, hi = img(A)
        if hend:
            if not (lo <= I0 and hi >= I1): continue
        elif gend:
            if not (lo >= I0 and hi <= I1): continue
        elif hi < I0 or lo > I1:
            continue
        side = "g" if (hend or (not gend and abs(A[0]) >= 1)) else "h"
        st = gs if side == "g" else hs
        forced = gf if side == "g" else hf
        o = opts(st, allowed, maxrun)
        for x in (forced if forced else o):
            if x not in o: continue
            if side == "g":
                if x == "E": stack.append((A, hs, gs, hend, True, hf, None, d + 1, hw, gw))
                else: stack.append((comp(A, LET[x]), hs, nxt(gs, x), hend, gend, hf, None, d + 1, hw, gw + x))
            else:
                if x == "E": stack.append((A, hs, gs, True, gend, None, gf, d + 1, hw, gw))
                else: stack.append((comp(inv(LET[x]), A), nxt(hs, x), gs, hend, gend, None, gf, d + 1, hw + x, gw))
    return ("no-collision", nodes, deep, None)

def word_map(w):
    A = ID
    for x in w: A = comp(A, LET[x])
    return A

BAD_N = ((1 - NU) / (1 + NU), (1 + NU) / (1 - NU))          # = [9/10, 10/9]: rho J_N meets J_N
BAD_K = ((1 - NU) / 2, (1 + NU) / (2 - 2 * KAP))             # = [9/19, 100/171]: rho J_K meets J_N
def bad(r): return BAD_N[0] <= r <= BAD_N[1] or BAD_K[0] <= r <= BAD_K[1]

def candidate_b(a):
    # rho = lam^a/mu^b in (0.3, 3.5) forces b within 1 of 0.7858 a; the range below is generous
    return range(max(1, int(a * 0.78) - 6), int(a * 0.79) + 8)

def random_reduced(rng, start, length, maxrun):
    w = start; st = nxt(("none", 0), start)
    while True:
        o = opts(st, REDUCED, maxrun)
        if len(w) >= length and "E" in o: return w
        x = rng.choice([c for c in o if c != "E"]); w += x; st = nxt(st, x)

if __name__ == "__main__":
    OUT = dict(AMAX=AMAX)
    t0 = time.time()
    # O1
    for y in [F(0), F(1, 9), F(3, 7)]:
        assert comp(LET["N"], (y, F(0)))[1] == F(10, 19) * comp(LET["K"], (y, F(0)))[1]
    assert img(LET["K"])[1] == I1 and all(img(LET[x])[0] >= I0 and img(LET[x])[1] <= I1 for x in LET)
    assert BAD_N == (F(9, 10), F(10, 9)) and BAD_K == (F(9, 19), F(100, 171))
    OUT["O1"] = dict(I=[str(I0), str(I1)], N_is_10_19_K=True, bad_N=[str(v) for v in BAD_N], bad_K=[str(v) for v in BAD_K])
    # O2 + O3
    badpairs, checked, worst = [], [], dict(nodes=0, depth=0)
    lp, mp_ = F(1), F(1)
    for a in range(1, AMAX + 1):
        for b in candidate_b(a):
            rho = LAM ** a / MU ** b
            if not (F(3, 10) < rho < F(7, 2)): continue
            st, nodes, depth, wit = descend((rho, F(0)), hforced="N", gforced="NK")
            isbad = bad(rho)
            checked.append((a, b, float(rho), isbad, st, nodes, depth))
            if isbad: badpairs.append((a, b, float(rho), st, nodes, depth))
            if st != "no-collision": print("!!", a, b, float(rho), st, wit, flush=True)
            if nodes > worst["nodes"]: worst = dict(nodes=nodes, depth=depth, a=a, b=b, rho=float(rho))
        if a % 1000 == 0: print("  a =", a, "pairs", len(checked), "bad", len(badpairs), round(time.time() - t0, 1), "s", flush=True)
    # every bad pair is inside the checked range: check that no b outside candidate_b can give rho in (0.3, 3.5)
    for a in range(1, AMAX + 1, 997):
        for b in (min(candidate_b(a)) - 1, max(candidate_b(a)) + 1):
            if b >= 1: assert not (F(3, 10) < LAM ** a / MU ** b < F(7, 2))
    allok = all(c[4] == "no-collision" for c in checked)
    OUT["O2_bad_pairs"] = dict(count=len(badpairs), first=badpairs[:12], last=badpairs[-3:],
                               closest_to_1=sorted(badpairs, key=lambda p: abs(p[2] - 1))[:5])
    OUT["O3_descent"] = dict(pairs_checked=len(checked), all_no_collision=allok, worst=worst,
                             max_depth=max(c[6] for c in checked), seconds=round(time.time() - t0, 1),
                             per_pair=[c for c in checked if c[3]])
    print("O2/O3: bad pairs", len(badpairs), "checked", len(checked), "all no-collision", allok, "worst", worst)
    # O4 sanity
    rng = random.Random(11); found = 0; trials = 300
    for _ in range(trials):
        # planted words with runs of at most 3 letters, searched in that language (the unbounded language has
        # infinite branches inside two runs when A is a pure scaling; the bad-pair searches above never met one)
        h = random_reduced(rng, rng.choice("NK"), rng.randint(2, 12), 3)
        g = random_reduced(rng, rng.choice("NK"), rng.randint(2, 12), 3)
        if h == g: g = g + "N"
        A0 = comp(word_map(h), inv(word_map(g)))
        st = descend(A0, hforced=h[0], gforced=g[0], maxrun=3, maxdepth=len(h) + len(g) + 2, depth_prune=True)[0]
        found += st == "collision"
    st_unred = descend((LAM / MU, F(0)), hforced="NK", gforced="NK", allowed=UNREDUCED)
    assert found == trials and st_unred[0] == "collision"
    OUT["O4_sanity"] = dict(planted=trials, planted_maxrun=3, found=found, unreduced_first_collision=[st_unred[0], st_unred[3]])
    print("O4: planted collisions found", found, "/", trials, "; unreduced language a=b=1:", st_unred[0], st_unred[3])
    json.dump(OUT, open("original_descent.json", "w"), default=str)
    print("total", round(time.time() - t0, 1), "s")
