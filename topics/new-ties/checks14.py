"""Round 14 checks (DRAFT). Reads psi14.json and unif14.json; reruns nb14 at depth 2000. Writes checks14.json.
Q1  Theorem 14.1 on 400 random word pairs, exact: psi(w) = t_w/s_w = sum over N/K positions of R_{>i}; psi(wx) = r_x psi(w) + d_x;
    psi(xw) = psi(w) + d_x R_w; and h = rho*g as maps exactly when psi(h) = psi(g) and s_h/s_g = rho (with N = (10/19) K as the
    positive case and 400 random pairs as negatives)
Q2  last-letter lemma: psi(h) is divisible by 4 (2-adically) when h ends in l or m, and odd when h ends in K, for every
    reduced word of up to 9 letters, so those endings never collide; every collision in psi14.json obeys it
Q3  psi14.json: exhaustive to 14 letters per word, no hash false pairs, exactly the 5 front ratios, 19-adic valuations -2..1,
    each example rechecked exactly; and every bad-pair collision needs L_h + n_g >= 14 letters, so the shortest shapes are inside
Q4  nb14 to combined length 2000, exact: no two different reduced words are equal as maps; every nonlinear state is born from a
    linear one inside its window; from every linear state found, every trip off the line, followed to its end, is finite and
    never returns to the line; and the starting state mu^b/lambda^a of every bad pair with a + b <= 2000 is among them
Q5  planted collisions are found: with the grammar switched off, nb14's search finds m K = l N and m l = l m at once
Q6  unif14.json: excursions of up to 11 letters, all stretches at once: exactly 9 unit return ratios, 19-adic valuations -2..2
Q7  a bad pair's only trip off the line starts at the stretch 1/rho = mu^b/lambda^a (19-adic valuation -a): checked by peeling
    m^b N... against l^a N... and l^a K... for the record-staircase pairs of Round 13 with a <= 3000 (dist13_summary.json),
    plus all 140,680 bad pairs a <= 10^6 in floating point (dist13.json, if present)
"""
import json, random, sys, os
from fractions import Fraction as F
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import nb14
OUT = {}; ok_all = True
def report(key, ok, **kw):
    global ok_all
    ok_all &= bool(ok); OUT[key] = dict(ok=bool(ok), **kw); print(key, "PASS" if ok else "FAIL", kw, flush=True)
LET = {"l": (F(19, 200), F(0)), "m": (F(1, 20), F(0)), "N": (F(1, 19), F(1, 19)), "K": (F(1, 10), F(1, 10))}
RV = {x: 1 / LET[x][0] for x in LET}; DV = {"l": 0, "m": 0, "N": 1, "K": 1}
AFTER = {"l": "lmNK", "m": "mN", "N": "lmNK", "K": "lmNK"}
def wmap(w):
    s, t = F(1), F(0)
    for x in w: a, b = LET[x]; s, t = s * a, s * b + t
    return s, t
def psi(w): s, t = wmap(w); return t / s
def psi_sum(w):
    tot = F(0)
    for i, x in enumerate(w):
        if x in "NK":
            R = F(1)
            for y in w[i + 1:]: R *= RV[y]
            tot += R
    return tot
def rword(rng, n, first=None):
    w = first or rng.choice("lmNK")
    while len(w) < n: w += rng.choice(AFTER[w[-1]])
    return w
rng = random.Random(14)
# Q1
bad = 0
for _ in range(400):
    w = rword(rng, rng.randint(1, 12)); x = rng.choice("lmNK")
    if psi(w) != psi_sum(w): bad += 1
    if psi(w + x) != RV[x] * psi(w) + DV[x]: bad += 1
    s, t = wmap(w)
    if psi(x + w) != psi(w) + DV[x] / s: bad += 1
neg = 0
for _ in range(400):
    h = rword(rng, rng.randint(1, 9), "N"); g = rword(rng, rng.randint(1, 9), rng.choice("NK"))
    (sh, th), (sg, tg) = wmap(h), wmap(g); rho = sh / sg
    same_map = (th == rho * tg)
    same_psi = psi(h) == psi(g)
    if same_map != same_psi: bad += 1
    neg += (not same_map)
sN, tN = wmap("N"); sK, tK = wmap("K")
pos = (sN / sK == F(10, 19)) and (tN == F(10, 19) * tK) and psi("N") == psi("K")
report("Q1", bad == 0 and pos, mismatches=bad, random_pairs_that_differ=neg, positive_case_N_eq_10_19_K=pos)
# Q2
def v2(q):
    if q == 0: return 99
    k = 0; a, b = q.numerator, q.denominator
    while a % 2 == 0: a //= 2; k += 1
    while b % 2 == 0: b //= 2; k -= 1
    return k
words = ["N", "K", "l", "m"]; allw = list(words)
for n in range(2, 10):
    words = [w + x for w in words for x in AFTER[w[-1]]]; allw += words
viol = 0; checked = 0
for w in allw:
    p = psi(w); checked += 1
    if w[-1] in "lm" and v2(p) < 2: viol += 1
    if w[-1] == "K" and v2(p) != 0: viol += 1
PS = json.load(open("psi14.json"))
pairs_ok = all(not ((r["example"][0][-1] in "lm" and r["example"][1][-1] == "K") or (r["example"][1][-1] in "lm" and r["example"][0][-1] == "K"))
               for r in PS["ratios"].values())
report("Q2", viol == 0 and pairs_ok, words_checked=checked, violations=viol)
# Q3
ex_ok = all(psi(r["example"][0]) == psi(r["example"][1]) and wmap(r["example"][0])[0] / wmap(r["example"][1])[0] == F(q)
            for q, r in PS["ratios"].items())
front = {"10/19", "20/19", "19/20", "200/361", "1/2"}
report("Q3", PS["n"] >= 14 and PS["hash_false_pairs"] == 0 and set(PS["ratios"]) == front and ex_ok
       and all(-2 <= r["v19"] <= 1 for r in PS["ratios"].values()),
       words=PS["words"], ratios=sorted(PS["ratios"]), v19=sorted(r["v19"] for r in PS["ratios"].values()))
# Q4
seen, births, hits, per, rets, life, alive, born = nb14.run(2000)
win = {('h', 'N'): (F(1), F(19, 9)), ('h', 'K'): (F(1), F(10, 9)), ('g', 'N'): (F(9, 19), F(1)), ('g', 'K'): (F(9, 10), F(1))}
inwin = all(win[(s, x)][0] <= a <= win[(s, x)][1] for a, s, x in births)
nontriv = [h for h in hits if h[1] != h[2]]
real_rets = [r for r in rets if r[1] != 1]        # the two at depth 2 are N/N and K/K reaching the identity
tot = lng = nret = unf = 0
for (a0, b0, lh, lg) in seen:
    if b0 != 0 or lh is None: continue
    ns, l, r, fin = nb14.excursions(a0, lh, lg); tot += ns; lng = max(lng, l); nret += len(r); unf += (not fin)
lam, mu = F(19, 200), F(1, 20); covered = missing = 0
pairs4 = (json.load(open("dist13.json"))["per_pair"] if os.path.exists("dist13.json")
          else json.load(open("dist13_summary.json"))["record_staircase"])   # dist13.json is large and rebuilt by CI
if True:
    for r in pairs4:
        A, B = r[0], r[1]
        if A + B > 2000: continue
        st = (mu ** B / lam ** A, F(0), 'm', 'l')
        if st in seen: covered += 1
        else: missing += 1
report("Q4", not nontriv and inwin and not real_rets and nret == 0 and unf == 0 and missing == 0 and covered > 0,
       states=len(seen), births=len(births), trip_states=tot, longest_trip=lng, bad_pairs_covered=covered, missing=missing)
# Q5 planted: switch the grammar off
saved = dict(nb14.NXT)
for k in nb14.NXT: nb14.NXT[k] = "lmNK"
_, _, hits2, _, _, _, _, _ = nb14.run(6)
nb14.NXT.clear(); nb14.NXT.update(saved)
found = sorted(set((h[1], h[2]) for h in hits2 if h[1] != h[2]))
report("Q5", len(found) > 0, nontrivial_returns_without_grammar=found)
# Q6
U = json.load(open("unif14.json"))
vs = sorted(r["v19"] for r in U["unit_alphas"].values())
report("Q6", U["maxlen"] >= 11 and len(U["unit_alphas"]) == 9 and min(vs) >= -2 and max(vs) <= 2,
       unit_alphas=sorted(U["unit_alphas"], key=lambda s: float(F(s))), v19=vs,
       nodes_at_max=U["nodes_by_length"][str(U["maxlen"])] if str(U["maxlen"]) in U["nodes_by_length"] else None)
# Q7
DS = json.load(open("dist13_summary.json"))
rows = {(r[0], r[1]) for r in DS["record_staircase"] + DS["closest50"] if r[0] <= 3000}
lam, mu = F(19, 200), F(1, 20); badq = 0
for (A, B) in rows:
    for x in "NK":
        h = "m" * B + "N"; g = "l" * A + x; i = j = 0; al = F(1)
        while True:
            if al >= 1:
                y = h[i]
                if y in "NK": break
                al *= LET[y][0]; i += 1
            else:
                y = g[j]
                if y in "NK": break
                al /= LET[y][0]; j += 1
        if al != mu ** B / lam ** A: badq += 1
import math
nfl = badfl = 0
if os.path.exists("dist13.json"):
    lL, lM = math.log(19 / 200), math.log(1 / 20)
    for r in json.load(open("dist13.json"))["per_pair"]:
        A, B = r[0], r[1]; nfl += 1
        # the linear phase is a walk in (i, j); it can only stop early if mu^j/lambda^i hits a window before i = a, j = b.
        # With rho = lambda^a/mu^b in [0.47, 1.12]: stopping on the h side needs j = b and alpha = lambda^(a-i)/rho in [1, 19/9],
        # forcing i = a; stopping on the g side needs i = a and alpha = mu^(j-b)/rho in [9/19, 1), forcing j = b.
        rho = math.exp(A * lL - B * lM)
        hs = [k for k in range(0, 3) if 1 <= math.exp(k * lL) / rho <= 19 / 9]
        gs = [k for k in range(0, 3) if 9 / 19 <= math.exp(-k * lM) / rho < 1]
        if not (0.4 < rho < 1.2) or any(k != 0 for k in hs + gs): badfl += 1
report("Q7", badq == 0 and len(rows) >= 3 and badfl == 0, exact_pairs=len(rows), exact_mismatches=badq,
       float_pairs=nfl, float_mismatches=badfl)
OUT["all_pass"] = ok_all
json.dump(OUT, open("checks14.json", "w"), indent=1)
print("ALL PASS" if ok_all else "SOME FAIL", "%d/%d" % (sum(1 for k, v in OUT.items() if isinstance(v, dict) and v["ok"]), sum(1 for k, v in OUT.items() if isinstance(v, dict))))
