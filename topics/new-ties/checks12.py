"""Round 12 checks (DRAFT). Reads cover12.json and original_descent.json. Writes checks12.json.
Y1  range search at a single stretch agrees with Round 11's descend on the first 120 bad pairs
Y2  planted solutions: the range search never says "every branch cut" on a range that contains one, and reaches it often
Y3  every touching point built from words of up to 4 letters lies in the unproved set, none in a proved piece
Y4  a second, separately written range search (breadth-first, other peeling rule) re-proves 300 random proved pieces
Y5  the bad pairs with a <= 10000: how many have their stretch in a proved piece, against the log-measure prediction
Y6  the hard spot 18/19: the branch h = N N m^n, g = N K^n is never cut (n <= 200), so no range containing 18/19 can finish
Y8  ratio points x/y of the fractal (x after N, y after N or K), sampled with 30-letter words, never fall in a proved piece
Y9  in a cluster, raising the budget from 2,000 to 8,000 nodes leaves exactly the same pieces unproved
Y7  pieces are disjoint, sorted, and with the open intervals tile each window exactly
"""
from fractions import Fraction as F
import sys, os, json, random, math, time
from collections import deque
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import cover12 as C
_a = sys.argv; sys.argv = [_a[0]]
import original as O
sys.argv = _a

OUT = {}; ok_all = True
def report(key, ok, **kw):
    global ok_all
    ok_all &= bool(ok); OUT[key] = dict(ok=bool(ok), **kw); print(key, "PASS" if ok else "FAIL", kw if len(str(kw)) < 400 else "", flush=True)

cov = json.load(open("cover12.json"))
proved = {w: [(F(a), F(b)) for a, b, _ in cov["windows"][w]["proved_pieces"]] for w in cov["windows"]}
opened = {w: [(F(a), F(b)) for a, b in cov["windows"][w]["open"]] for w in cov["windows"]}
allproved = sorted(p for w in proved for p in proved[w])
def in_proved(r):
    lo, hi = 0, len(allproved)
    while lo < hi:
        mid = (lo + hi) // 2
        if allproved[mid][1] < r: lo = mid + 1
        else: hi = mid
    return lo < len(allproved) and allproved[lo][0] <= r <= allproved[lo][1]

# Y1
od = json.load(open("original_descent.json"))
pairs = [(p[0], p[1]) for p in od["O3_descent"]["per_pair"]][:120]
agree = 0
for a, b in pairs:
    r = O.LAM ** a / O.MU ** b
    s1 = C.search(r, r, maxnodes=200000)[0]; s2 = O.descend((r, F(0)), hforced="N", gforced="NK")[0]
    agree += (s1 == "ok") == (s2 == "no-collision")
report("Y1_agree_with_round11", agree == len(pairs), pairs=len(pairs), agree=agree)

# Y2
rng = random.Random(12); found = 0; inconclusive = 0; false_ok = 0; trials = 120
for _ in range(trials):
    h = O.random_reduced(rng, rng.choice("NK"), rng.randint(2, 9), 3)
    g = O.random_reduced(rng, rng.choice("NK"), rng.randint(2, 9), 3)
    if h == g: g = g + "N"
    sh, th = O.word_map(h); sg, tg = O.word_map(g)
    s = sh / sg; t = th - s * tg            # h = A o g with A(z) = s z + t; range parameter rho multiplies s
    eps = F(1, 10 ** 9)
    st, n, hits = C.search(1 - eps, 1 + eps, hforced=h[0], gforced=g[0], maxnodes=15000, start=(s, F(0), t))
    if F(1) in hits: found += 1          # the planted solution was reached (budget may still run out after it)
    elif st == "ok": false_ok += 1       # would be a bug: "every branch cut" although a solution is in range
    else: inconclusive += 1
st_u, n_u, hits_u = C.search(F(189, 100), F(191, 100), hforced="NK", gforced="NK", maxnodes=200000, allowed=C.UNREDUCED)
unred = F(19, 10) in hits_u
report("Y2_planted", false_ok == 0 and found >= trials // 3 and unred, trials=trials, found=found,
       inconclusive=inconclusive, false_ok=false_ok, unreduced_collision_found=unred, unreduced_status=st_u)

# Y3
def words(start, L):
    out = []
    def rec(w, st):
        if w: out.append((w, st))
        if len(w) == L: return
        for x in (start if not w else C.REDUCED[st]):
            if x == "E" or (w and x not in C.REDUCED[st]): continue
            rec(w + x, C.nxt(st, x))
    rec("", "none"); return out
def ev(w, z):
    for x in reversed(w):
        sx, cx = C.LET[x]; z = sx * z + cx
    return z
touch = set()
for X, sx in words("N", 4):
    for Y, sy in words("NK", 4):
        if sy != "mu": touch.add(ev(X, F(0)) / ev(Y, F(1, 9)))      # left end of X meets right end of rho*Y
        if sx != "mu": touch.add(ev(X, F(1, 9)) / ev(Y, F(0)))
touch = sorted(r for r in touch if O.bad(r))
inside = [r for r in touch if in_proved(r)]
report("Y3_touching_points_unproved", not inside, touching_points=len(touch), in_proved_pieces=len(inside),
       examples=[str(r) for r in touch[:5]])

# Y4
def search_bfs(r1, r2, maxnodes=400000):
    """independent: breadth-first, A = (aS*rho, bT*rho + cU) kept as (S, T, U); uses s > 0 so A(0) < A(1/9);
    peel g when the slope at the LEFT end r1 is >= 1."""
    q = deque([(F(1), F(0), F(0), "none", "none", 0, 0, "N", "NK")]); seen = set(); n = 0
    while q:
        node = q.popleft()
        if node in seen: continue
        seen.add(node); n += 1
        if n > maxnodes: return "big", n
        S, T, U, hs, gs, he, ge, hf, gf = node
        if he and ge:
            if r1 <= 1 / S <= r2 and T / S + U == 0: return "hit", n
            continue
        a0 = lambda r: r * T + U; a1 = lambda r: r * (S / 9 + T) + U
        A0 = (a0(r1), a0(r2)); A1 = (a1(r1), a1(r2))
        if he: okk = min(A0) <= 0 and max(A1) >= O.I1
        elif ge: okk = max(A0) >= 0 and min(A1) <= O.I1
        else: okk = max(A1) >= 0 and min(A0) <= O.I1
        if not okk: continue
        side = "g" if (he or (not ge and r1 * S >= 1)) else "h"
        kind = gs if side == "g" else hs; forced = gf if side == "g" else hf
        for x in (forced or C.REDUCED[kind]):
            if x not in C.REDUCED[kind]: continue
            if x == "E":
                q.append((S, T, U, hs, gs, he, 1, hf, None) if side == "g" else (S, T, U, hs, gs, 1, ge, None, gf)); continue
            s_, c_ = C.LET[x]
            if side == "g": q.append((S * s_, S * c_ + T, U, hs, C.nxt(gs, x), he, ge, hf, None))
            else: q.append((S / s_, T / s_, (U - c_) / s_, C.nxt(hs, x), gs, he, ge, None, gf))
    return "ok", n
rng = random.Random(4); sample = rng.sample(allproved, min(300, len(allproved)))
res4 = [search_bfs(a, b)[0] for a, b in sample]
report("Y4_independent_reproof", all(r == "ok" for r in res4), sampled=len(sample),
       reproved=sum(r == "ok" for r in res4), big=sum(r == "big" for r in res4), hit=sum(r == "hit" for r in res4))

# Y5
bp = [(p[0], p[1]) for p in od["O3_descent"]["per_pair"]]
inp = [ab for ab in bp if in_proved(O.LAM ** ab[0] / O.MU ** ab[1])]
logw = sum(math.log(float(b) / float(a)) for w in cov["windows"] for a, b in [tuple(map(F, cov["windows"][w]["window"]))])
logp = sum(math.log(float(b) / float(a)) for a, b in allproved)
report("Y5_bad_pairs_in_proved_set", True, bad_pairs_a_le_10000=len(bp), in_proved=len(inp),
       fraction=len(inp) / len(bp), log_measure_prediction=logp / logw)

# Y6
r = F(18, 19)
assert ev("NN", F(0)) == r * ev("N", F(1, 9))
A = (r, F(0))   # h = A o g
uncut = True
hw, gw = "NN" + "m" * 200, "N" + "K" * 200
for n in range(1, 201):
    hp = O.word_map(hw[:n + 2]); gp = O.word_map(gw[:n + 1])
    Bm = O.comp(O.inv(hp), O.comp(A, gp))   # remaining map after peeling n+2 letters of h and n+1 of g
    lo, hi = Bm[1], Bm[0] / 9 + Bm[1]
    uncut &= (hi >= 0 and lo <= O.I1)
report("Y6_hard_spot_18_19", uncut, identity="N(N(0)) = (18/19) N(1/9)", branch="h = NN m^n, g = N K^n", n_max=200)

# Y8
def rand_word(rng, start, L):
    w = start; st = C.nxt("none", start)
    while len(w) < L:
        o = [c for c in C.REDUCED[st] if c != "E"]; x = rng.choice(o); w += x; st = C.nxt(st, x)
    return w
rng = random.Random(8); tested = 0; missed = []
for _ in range(3000):
    X = rand_word(rng, "N", 30); Y = rand_word(rng, rng.choice("NK"), 30)
    lo_, hi_ = ev(X, F(0)) / ev(Y, F(1, 9)), ev(X, F(1, 9)) / ev(Y, F(0))   # every x/y with x in X(I), y in Y(I)
    if not (O.bad(lo_) and O.bad(hi_)): continue
    tested += 1
    if in_proved(lo_) and any(a <= lo_ and hi_ <= b for a, b in allproved): missed.append(float(lo_))
report("Y8_attractor_ratios_unproved", not missed, cylinder_pairs_in_windows=tested, inside_proved=len(missed))

# Y9
def opens(budget):
    out = []; C.cover(F(947, 1000), F(9485, 10000), budget, F(1, 2 ** 17), out)
    return [(p[0], p[1]) for p in out if p[2] is None]
o2, o8 = opens(2000), opens(8000)
report("Y9_budget_does_not_help", o2 == o8, region="[0.947, 0.9485]", open_pieces_2000=len(o2), open_pieces_8000=len(o8),
       open_measure=float(sum(q - p for p, q in o2)))

# Y7
tile = True
for w in cov["windows"]:
    a, b = map(F, cov["windows"][w]["window"])
    segs = sorted(proved[w] + opened[w])
    tile &= segs[0][0] == a and segs[-1][1] == b and all(segs[i][1] == segs[i + 1][0] for i in range(len(segs) - 1))
report("Y7_tiling", tile)

OUT["all_pass"] = ok_all
json.dump(OUT, open("checks12.json", "w"), indent=1, default=str)
print("ALL PASS" if ok_all else "SOME FAIL")
