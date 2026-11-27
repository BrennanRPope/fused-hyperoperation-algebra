"""Round 13 checks (DRAFT). Reads safe13.json and dist13.json. Writes checks13.json.
Z1  a real collision is a point of R: N = (10/19) K, and dist_R(10/19) is 0 up to rounding (below 1e-100); so are 40 random finite ratios h(0)/g(0)
Z2  the touching point 18/19 lies in R: the branch-and-bound never proves a positive distance there
Z3  brute force: for 200 random stretches, every ratio u(0)/v(0) with words of up to 7 letters is at least the proved lower bound away
Z4  the closest 12 pairs recomputed in exact rational arithmetic, with the stretch enclosed in an interval of width 2^-400
Z5  two independent methods agree: each safe radius from the range search (Round 12 machinery) is at most the measured distance,
    and on the 30 closest pairs with a <= 10000 the two are reported side by side (gap = k + log2 dist >= 0)
Z6  bookkeeping: 1,406 bad pairs with a <= 10000 (as in Round 12), every pair resolved by both methods
"""
from fractions import Fraction as F
import sys, os, json, random, math, heapq
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import mpmath as mp
import rdist13 as D
mp.mp.dps = 160
OUT = {}; ok_all = True
def report(key, ok, **kw):
    global ok_all
    ok_all &= bool(ok); OUT[key] = dict(ok=bool(ok), **kw); print(key, "PASS" if ok else "FAIL", kw if len(str(kw)) < 500 else "", flush=True)

LETF = {"l": (F(19, 200), F(0)), "m": (F(1, 20), F(0)), "N": (F(1, 19), F(1, 19)), "K": (F(1, 10), F(1, 10))}
def wmap(w):
    a, c = F(1), F(0)
    for x in w: s, t = LETF[x]; a, c = a * s, a * t + c
    return a, c
def rand_word(rng, first, n):
    w = first
    while len(w) < n: w += rng.choice(D.AFTER[w[-1]])
    return w

# Z1
lb, ub, xw, yw, n, done = D.dist_R(mp.mpf(10) / 19, cap=2000)
rng = random.Random(13); worst = 0
for _ in range(40):
    h = rand_word(rng, "N", rng.randint(2, 12)); g = rand_word(rng, rng.choice("NK"), rng.randint(2, 12))
    r = wmap(h)[1] / wmap(g)[1]
    worst = max(worst, float(D.dist_R(mp.mpf(r.numerator) / r.denominator, cap=3000)[1]))
report("Z1_collision_is_in_R", ub < 1e-100 and worst < 1e-100, ten_19_ub=float(ub), worst_ub_random_finite_points=worst)

# Z2
lb, ub, xw, yw, n, done = D.dist_R(mp.mpf(18) / 19, cap=20000)
report("Z2_touching_point_18_19", (not done) and lb == 0, nodes=n, lb=float(lb), ub=float(ub))

# Z3
rng = random.Random(31); bad = 0; pts = []
def allwords(first, L):
    out = []; st = [first]
    while st:
        w = st.pop(); out.append(w)
        if len(w) < L: st.extend(w + x for x in D.AFTER[w[-1]])
    return out
X = sorted({wmap(w)[1] for w in allwords("N", 7)}); Y = sorted({wmap(w)[1] for w in allwords("N", 7) + allwords("K", 7)})
ratios = sorted(float(x / y) for x in X for y in Y)
import bisect
for _ in range(200):
    r = rng.uniform(0.45, 1.12)
    lb = float(D.dist_R(mp.mpf(r))[0]); i = bisect.bisect_left(ratios, r)
    near = min(abs(r - ratios[j]) for j in (i - 1, i) if 0 <= j < len(ratios))
    bad += near < lb * (1 - 1e-12)
report("Z3_brute_force_never_beats_lower_bound", bad == 0, stretches=200, ratio_points=len(ratios), violations=bad)

# Z4 exact
def exact_lb(a, b, cap=100000):
    r = mp.exp(a * mp.log(mp.mpf(19) / 200) + b * mp.log(20)); E = F(1, 2 ** 400)
    c = F(int(mp.nint(r * 2 ** 420)), 2 ** 420); r1, r2 = c - E, c + E
    def lbf(u, v):
        x0, x1 = u[2], u[1] / 9 + u[2]; y0, y1 = v[2], v[1] / 9 + v[2]
        lo, hi = x0 / y1, x1 / y0
        return lo - r2 if r2 < lo else (r1 - hi if r1 > hi else F(0))
    def ub(u, v): q = u[2] / v[2]; return max(abs(q - r1), abs(q - r2))
    heap = []; cnt = 0; best = None
    def push(u, v):
        nonlocal cnt, best
        cnt += 1; heapq.heappush(heap, (lbf(u, v), cnt, u, v)); q = ub(u, v); best = q if best is None or q < best else best
    def ext(w, x): s, t = LETF[x]; return (w[0] + x, w[1] * s, w[1] * t + w[2])
    for v in "NK": push(("N",) + LETF["N"], (v,) + LETF[v])
    nodes = 0
    while heap and nodes < cap:
        l, _, u, v = heapq.heappop(heap)
        if 2 * l >= best: return l, best
        nodes += 1
        if (u[1] / 9) / u[2] >= (v[1] / 9) / v[2]:
            for x in D.AFTER[u[0][-1]]: push(ext(u, x), v)
        else:
            for y in D.AFTER[v[0][-1]]: push(u, ext(v, y))
    return heap[0][0], best
dd = json.load(open("dist13.json"))
close = sorted(dd["per_pair"], key=lambda p: p[4])[:12]; agree = 0; rows = []
for p in close:
    l, u = exact_lb(p[0], p[1])
    le = math.log2(l.numerator) - math.log2(l.denominator) if l > 0 else None
    rows.append((p[0], p[1], round(p[3], 2), round(le, 2) if le else None))
    agree += le is not None and abs(le - p[3]) < 1.01
report("Z4_exact_recheck_closest", agree == len(close), rows=rows)

# Z5
sf = json.load(open("safe13.json")); dmap = {(p[0], p[1]): p for p in dd["per_pair"]}
viol = 0; gaps = []
for p in sf["per_pair"]:
    if p[3] is None: continue
    d = dmap[(p[0], p[1])]
    if 2.0 ** (-p[3]) - 2.0 ** -161 > 2.0 ** d[4] * (1 + 1e-9): viol += 1
for p in sorted([q for q in sf["per_pair"] if q[3]], key=lambda q: -q[3])[:30]:
    gaps.append(p[3] + dmap[(p[0], p[1])][4])
report("Z5_two_methods_agree", viol == 0 and min(gaps) >= -0.01, median_gap=sorted(gaps)[len(gaps)//2], max_gap=max(gaps), pairs=len(sf["per_pair"]), violations=viol,
       gap_k_minus_log2dist_on_30_closest=[round(g, 2) for g in gaps])

# Z6
n10 = sum(1 for p in dd["per_pair"] if p[0] <= 10000)
report("Z6_bookkeeping", n10 == 1406 and sf["pairs"] == 1406 and sf["unresolved"] == 0 and dd["not_done"] == 0,
       pairs_a_le_10000=n10, safe_unresolved=sf["unresolved"], dist_pairs=dd["pairs"], dist_not_done=dd["not_done"])
json.dump(OUT, open("checks13.json", "w"), indent=1, default=str)
print("ALL PASS" if ok_all else "SOME FAIL")
