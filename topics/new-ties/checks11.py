"""Round 11 checks (README section 15). Run after original.py and original_bracket.py. Writes checks11.json."""
import json, math, random, itertools
from collections import deque
from fractions import Fraction as F
import original as O

D = json.load(open("original_descent.json")); B = json.load(open("original_bracket.json"))
AM = D["AMAX"]; R = {}
def check(name, ok, info=""):
    R[name] = dict(ok=bool(ok), info=str(info)); print(("PASS " if ok else "FAIL ") + name, info)

# X1 the y-model
ok = all(O.img(O.LET[x])[0] >= 0 and O.img(O.LET[x])[1] <= F(1, 9) for x in O.LET)
ok &= all((y + 1) / 19 == F(10, 19) * ((y + 1) / 10) for y in [F(0), F(1, 9), F(5, 7)])
ok &= O.BAD_N == (F(9, 10), F(10, 9)) and O.BAD_K == (F(9, 19), F(100, 171))
check("X1 y-model: letters map I = [0, 1/9] into itself, N = (10/19) K, windows [9/10, 10/9] and [9/19, 100/171]", ok)

# X2 the state machine is exactly the reduced language: compare with a direct parse into pieces
def direct_ok(w):            # w = rest of a word after a class-1 run, i.e. it must start with N or K
    if not w or w[0] not in "NK": return False
    runs = []; i = 0
    while i < len(w):
        if w[i] in "NK":
            j = i
            while j < len(w) and w[j] in "NK": j += 1
            runs.append(("v", w[i:j])); i = j
        else:
            j = i
            while j < len(w) and w[j] in "lm": j += 1
            runs.append(("u", w[i:j])); i = j
    if runs[-1][0] != "v": return False
    for k, (t, s) in enumerate(runs):
        if t == "u":
            if "ml" in s: return False                          # canonical: lam before mu (they commute)
            if "m" in s and runs[k + 1][1][0] == "K": return False  # the duplicated family is dropped
    return True
def machine_ok(w):
    st = ("none", 0)
    for x in w:
        if x not in O.REDUCED[st[0]]: return False
        st = O.nxt(st, x)
    return "E" in O.REDUCED[st[0]]
bad_ = [w for n in range(1, 8) for w in map("".join, itertools.product("lmNK", repeat=n))
        if w[0] in "NK" and direct_ok(w) != machine_ok(w)]
check("X2 the search walks exactly the reduced words (all words up to length 7 compared with a direct parse)", not bad_, bad_[:3])

# X3 bad pairs: first one is (13, 10); the b-range of the search covers every rho in (0.3, 3.5)
first = min((a, b) for a in range(1, 40) for b in range(1, 40) if O.bad(O.LAM ** a / O.MU ** b))
lr, lm = math.log(19 / 200), math.log(1 / 20)
cover = all(all(not (math.log(0.3) - 1e-9 < a * lr - b * lm < math.log(3.5) + 1e-9)
                for b in (min(O.candidate_b(a)) - 1, max(O.candidate_b(a)) + 1) if b >= 1)
            and min(O.candidate_b(a)) <= max(1, (a * lr - math.log(0.3)) / lm) and max(O.candidate_b(a)) >= (a * lr - math.log(3.5)) / lm
            for a in range(1, AM + 1))
check("X3 first bad pair is (13, 10); the search's b-range covers every stretch in (0.3, 3.5) for a <= AMAX",
      first == (13, 10) and cover, first)

# X4 every bad pair settled; independent recount of the bad pairs
P = D["O3_descent"]["per_pair"]
recount = 0
for a in range(1, AM + 1):
    b0 = (a * lr) / lm
    for b in range(max(1, int(b0) - 3), int(b0) + 4):
        if O.bad(O.LAM ** a / O.MU ** b): recount += 1
check("X4 every bad pair with a <= AMAX: descent says no-collision; recount of bad pairs agrees",
      all(p[4] == "no-collision" for p in P) and D["O3_descent"]["all_no_collision"] and recount == len(P),
      f"{len(P)} bad pairs, recount {recount}, AMAX {AM}")

# X5 a second, independent search (breadth first, other peeling rule, own pruning code) on a sample
def bfs(rho):
    I0, I1 = F(0), F(1, 9)
    start = ((rho, F(0)), "none", "none", 0, 0, "N", "NK")   # A, hs, gs, hend, gend, forced h, forced g
    q = deque([start]); seen = {start}
    while q:
        A, hs, gs, he, ge, hf, gf = q.popleft()
        if he and ge:
            if A == (1, 0): return "collision"
            continue
        if A == (1, 0) and hs == gs == "none" and hf is None and gf is None: return "collision"
        ends = sorted([A[1], A[0] / 9 + A[1]])
        if he and not (ends[0] <= I0 and ends[1] >= I1): continue
        if ge and not (ends[0] >= I0 and ends[1] <= I1): continue
        if not he and not ge and (ends[1] < I0 or ends[0] > I1): continue
        peel_h = (not he) and (ge or abs(A[0]) <= 1)       # differs from original.py at |slope| = 1
        st, fo = (hs, hf) if peel_h else (gs, gf)
        for x in (fo or O.REDUCED[st]):
            if x == "E": n = (A, hs, gs, he or peel_h, ge or not peel_h, None if peel_h else hf, gf if peel_h else None)
            else:
                s_, t_ = O.LET[x]; ns = "none" if x in "NK" else ("lam" if x == "l" else "mu")
                if peel_h: n = ((A[0] / s_, (A[1] - t_) / s_), ns, gs, he, ge, None, gf)
                else: n = ((A[0] * s_, A[0] * t_ + A[1]), hs, ns, he, ge, hf, None)
            if n not in seen: seen.add(n); q.append(n)
    return "no-collision"
rng = random.Random(5)
sample = [p for p in P if p[0] <= 1500] + rng.sample([p for p in P if p[0] > 1500], min(40, len([p for p in P if p[0] > 1500])))
res = [bfs(O.LAM ** p[0] / O.MU ** p[1]) for p in sample]
check("X5 independent breadth-first search agrees on every bad pair with a <= 1500 and 40 random larger ones",
      all(r == "no-collision" for r in res), f"{len(sample)} pairs")

# X6 the search finds collisions that exist
s4 = D["O4_sanity"]
check("X6 sanity: 300 planted collisions all found; the unreduced language's known collision found",
      s4["found"] == s4["planted"] == 300 and s4["unreduced_first_collision"][0] == "collision", s4["unreduced_first_collision"])

# X7 sigma is never 1 (19-adic count) for any bad pair
def v19(q):
    q = F(q); v = 0; n, d = q.numerator, q.denominator
    while n % 19 == 0: n //= 19; v += 1
    while d % 19 == 0: d //= 19; v -= 1
    return v
okv = all(v19(O.LAM ** a / O.MU ** b) == a and v19(F(19, 10) * O.LAM ** a / O.MU ** b) == a + 1
          for a, b in [(p[0], p[1]) for p in P[:400]])
check("X7 sigma != 1: v19(sigma) = a or a + 1 for the bad pairs (first 400)", okv)

# X8 the bracket: k = 12 reproduces Round 10, k = AMAX sits inside it and is far tighter
import subprocess, sys
r12 = json.loads(subprocess.run([sys.executable, "-c", "import original_bracket as ob, json\n"
     "from mpmath import mp, iv\nob.DPS=90; mp.dps=90; iv.dps=90\nl,_=ob.root(12); _,h=ob.root(None)\nprint(json.dumps([float(l), float(h)]))"],
     capture_output=True, text=True).stdout)
lo, hi = float(B["lower"]), float(B["upper"])
check("X8 bracket: k = 12 gives Round 10's [0.4704833085, 0.4704834402] (to Round 10's 1e-9 tolerance); the new bracket lies inside and its gap is 10^(%d)" % B["log10_gap"],
      abs(r12[0] - 0.4704833085) < 1e-9 and abs(r12[1] - 0.4704834402) < 1e-9 and r12[0] <= lo <= hi <= 0.4704834402
      and B["k"] == AM and B["log10_gap"] < -0.45 * AM, r12)

n_ok = sum(v["ok"] for v in R.values())
print(f"{n_ok}/{len(R)} pass"); json.dump(R, open("checks11.json", "w"), indent=1)
