"""Round 15 checks (DRAFT). Reads kappa15.json, sparse15.json, unif14.json, psi14.json. Writes checks15.json.
Q1  Lemma 15.1 on 400 random pairs (P, Q), exact: T = Q^-1 o P is y -> a y + b with a = s_P/s_Q and b = a psi(P) - psi(Q);
    the real cut (b <= 1/9 and b + a/9 >= 0) is the same as psi(Q)/(psi(P)+1/9) <= a <= (psi(Q)+1/9)/psi(P); and T is the
    identity's multiple (b = 0) exactly when a = psi(Q)/psi(P)
Q2  Proposition 15.2 on 400 random words w, exact: psi(Kw) = psi(Nw); psi(NNw) = psi(KNw) = psi(Nmw) = psi(Kmw);
    mK = lN and ml = lm as maps; psi(lw) = psi(mw) = psi(w); psi(wl) = (200/19) psi(w), psi(wm) = 20 psi(w)
Q3  every reduced word of up to 10 letters ending in N or K has the same fingerprint as its canonical form, and canon is
    idempotent
Q4  kappa15.json: n = 14, 29,246,206 words, no clash, distinct canonical forms = distinct fingerprints at every length, unit
    psi only for N, K, NN, KN; and a second, independent implementation (exact Fractions, all reduced words of up to 9 letters)
    gives the same counts for lengths 1..9
Q5  Theorem 15.5 by brute force: every word with exactly two N/K letters and l, m runs up to 60 between them: psi is a
    {2,5,19}-unit only for NN (first letter N; KN has the same psi); and 1 + 2^x 5^y is never a power of 19 for x, y <= 300
Q6  Corollary 15.4: among reduced words of up to 9 letters starting with N or K, every pair with equal psi has the same
    canonical form before trailing l, m are stripped (so the same trailing l and m counts), and the size ratios s_h/s_g, with
    1 added, are exactly the 9 return stretches of Conjecture 14.4, all with |19-adic valuation| <= 2; those 9 are exactly unif14.json's unit returns and
    contain psi14.json's 5 front ratios
Q7  sparse15.json: 3.2 million words with 2 to 4 N/K letters and long l, m runs: no clash, unit psi only for N and NN
Q8  the clash detector works: with canonical forms switched off, or with a wrong canonical form (second-letter rule left
    out), the same search to 10 letters reports clashes
"""
import json, random, sys, os, ast
from fractions import Fraction as F
from collections import defaultdict
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from kappa15 import kappa_int, canon
OUT = {}; ok_all = True
def report(key, ok, **kw):
    global ok_all
    ok_all &= bool(ok); OUT[key] = dict(ok=bool(ok), **kw); print(key, "PASS" if ok else "FAIL", kw, flush=True)
LET = {"l": (F(19, 200), F(0)), "m": (F(1, 20), F(0)), "N": (F(1, 19), F(1, 19)), "K": (F(1, 10), F(1, 10))}
AFTER = {"l": "lmNK", "m": "mN", "N": "lmNK", "K": "lmNK"}
def wmap(w):
    s, t = F(1), F(0)
    for x in w: a, b = LET[x]; s, t = s * a, s * b + t
    return s, t
def psi(w): s, t = wmap(w); return t / s if w else F(0)
def kap(q):
    q = F(q); n, d = q.numerator, q.denominator
    return kappa_int(n), kappa_int(d)
def v(q, p):
    q = F(q); n, d, k = q.numerator, q.denominator, 0
    while n % p == 0: n //= p; k += 1
    while d % p == 0: d //= p; k -= 1
    return k
def rword(rng, n, first=None):
    w = first or rng.choice("lmNK")
    while len(w) < n: w += rng.choice(AFTER[w[-1]])
    return w
def words_upto(n, starts="lmNK"):
    lev = list(starts)
    for L in range(1, n + 1):
        yield from lev
        lev = [w + x for w in lev for x in AFTER[w[-1]]]
rng = random.Random(15)
# Q1
bad = 0; cut_bad = 0; ret_bad = 0
for _ in range(400):
    P = rword(rng, rng.randint(1, 9), rng.choice("NK")); Q = rword(rng, rng.randint(1, 9), rng.choice("NK"))
    (sP, tP), (sQ, tQ) = wmap(P), wmap(Q)
    a = sP / sQ; b = (tP - tQ) / sQ           # Q^-1(P(y)) = (sP y + tP - tQ)/sQ
    if b != a * psi(P) - psi(Q): bad += 1
    for aa in (a, a * F(rng.randint(1, 50), rng.randint(1, 50))):
        bb = aa * psi(P) - psi(Q)
        if ((bb <= F(1, 9)) and (bb + aa / 9 >= 0)) != (psi(Q) / (psi(P) + F(1, 9)) <= aa <= (psi(Q) + F(1, 9)) / psi(P)): cut_bad += 1
        if (bb == 0) != (aa == psi(Q) / psi(P)): ret_bad += 1
report("Q1_lemma_15_1", bad == 0 and cut_bad == 0 and ret_bad == 0, pairs=400, b_formula_failures=bad, cut_failures=cut_bad, return_failures=ret_bad)
# Q2
bad = 0
for _ in range(400):
    w = rword(rng, rng.randint(1, 10), rng.choice("mN"))
    p = psi(w)
    if psi("K" + w) != psi("N" + w): bad += 1
    if not (psi("NN" + w) == psi("KN" + w) == psi("Nm" + w) == psi("Km" + w)): bad += 1
    if psi("l" + w) != p or psi("m" + w) != p: bad += 1
    if psi(w + "l") != F(200, 19) * p or psi(w + "m") != 20 * p: bad += 1
maps_ok = wmap("mK") == wmap("lN") and wmap("ml") == wmap("lm")
report("Q2_front_relations", bad == 0 and maps_ok, words=400, failures=bad, mK_eq_lN_and_ml_eq_lm=maps_ok)
# Q3
bad = 0; idem = 0; cnt = 0
for w in words_upto(10, "NK"):
    if w[-1] not in "NK": continue
    cnt += 1; c = canon(w)
    if kap(psi(w)) != kap(psi(c)): bad += 1
    if canon(c) != c: idem += 1
report("Q3_canonical_form_keeps_fingerprint", bad == 0 and idem == 0, words=cnt, fingerprint_changes=bad, not_idempotent=idem)
# Q4
K = json.load(open("kappa15.json"))
ok = (K["n"] == 14 and K["words"] == 29246206 and K["clashes_after_canonical_form"] == 0
      and K["distinct_canonical_forms_by_length"] == K["distinct_fingerprints_by_length"]
      and set(K["unit_psi_shapes"]) == {"N", "K", "NN", "KN"})
fp = set(); cf = set(); counts = []
lev = ["N", "K"]
for L in range(1, 10):
    for w in lev:
        if w[-1] in "NK": fp.add(kap(psi(w))); cf.add(canon(w))
    counts.append((len(fp), len(cf)))
    lev = [w + x for w in lev for x in AFTER[w[-1]]]
indep_ok = all(a == b == K["distinct_fingerprints_by_length"][i] for i, (a, b) in enumerate(counts))
report("Q4_kappa15_record", ok and indep_ok, words=K["words"], distinct_fingerprints=K["distinct_fingerprints_by_length"][-1],
       clashes=K["clashes_after_canonical_form"], unit_shapes=sorted(K["unit_psi_shapes"]), independent_counts_1_to_9=[a for a, b in counts])
# Q5
units = []
for x2 in "NK":
    for i in range(61):
        for j in range(61):
            if j and x2 == "K": continue
            c = F(10 if x2 == "K" else 19) * F(200, 19) ** i * F(20) ** j
            if kap(1 + c) == (1, 1): units.append("N" + "l" * i + "m" * j + x2)
mod3 = [(x, y) for x in range(2, 301) for y in range(1, 301) if kappa_int(1 + 2 ** x * 5 ** y) == 1]   # the case v2 >= 2, v5 >= 1
report("Q5_theorem_15_5", units == ["NN"] and mod3 == [], words_checked=2 * 61 * 61 - 60 * 61,
       unit_words_with_two_NK_letters=units, powers_of_19_of_form_1_plus_2x5y=mod3)
# Q6
NINE = {F(1), F(19, 10), F(10, 19), F(20, 19), F(19, 20), F(361, 200), F(200, 361), F(2), F(1, 2)}
by_psi = defaultdict(list)
for w in words_upto(9, "NK"): by_psi[psi(w)].append(w)
def full_canon(w):          # canon() without the final strip of trailing l, m
    w = w.lstrip("lm"); w = "N" + w[1:]
    if len(w) > 1 and w[1] == "N": w = "Nm" + w[2:]
    while True:
        u = w.replace("mK", "lN").replace("ml", "lm")
        if u == w: break
        w = u
    return w
ratios = set(); bad = 0; pairs = 0
for p, ws in by_psi.items():
    for i in range(len(ws)):
        for j in range(i + 1, len(ws)):
            h, g = ws[i], ws[j]; pairs += 1
            if full_canon(h) != full_canon(g): bad += 1
            r = wmap(h)[0] / wmap(g)[0]; ratios.add(r); ratios.add(1 / r)
U = json.load(open("unif14.json")); ua = U["unit_alphas"]
ua = ast.literal_eval(ua) if isinstance(ua, str) else ua
unif_set = {F(k) for k in ua}
P14 = json.load(open("psi14.json"))
p14 = set()
for key in ("ratios", "rho_values", "front_ratios"):
    if key in P14:
        val = P14[key]; p14 = {F(str(k)) for k in (val.keys() if isinstance(val, dict) else val)}; break
ok6 = bad == 0 and ratios | {F(1)} == NINE and NINE == unif_set and all(abs(v(r, 19)) <= 2 for r in NINE) and (not p14 or p14 <= NINE)
report("Q6_corollary_15_4", ok6, pairs_with_equal_psi=pairs, pairs_breaking_canonical_form=bad,
       ratios_seen=sorted(str(r) for r in ratios), unif14_unit_returns=sorted(str(r) for r in unif_set),
       psi14_ratios=sorted(str(r) for r in p14) or "key not found", max_abs_v19=max(abs(v(r, 19)) for r in NINE))
# Q7
S = json.load(open("sparse15.json"))
report("Q7_sparse15_record", S["clashes_after_canonical_form"] == 0 and S["unit_psi_words"] == ["N", "NN"] and S["words"] > 3_000_000,
       words=S["words"], plan=S["plan"], clashes=S["clashes_after_canonical_form"], unit_psi_words=S["unit_psi_words"])
# Q8
def wrong_canon(w):
    w = w.lstrip("lm"); w = "N" + w[1:]
    while True:
        u = w.replace("mK", "lN").replace("ml", "lm")
        if u == w: break
        w = u
    return w.rstrip("lm")
def clashes(cf):
    first = {}; n = 0
    for w in words_upto(10, "NK"):
        if w[-1] not in "NK": continue
        k = kap(psi(w)); c = cf(w)
        if k in first and first[k] != c: n += 1
        else: first.setdefault(k, c)
    return n
c_raw, c_wrong, c_right = clashes(lambda w: w), clashes(wrong_canon), clashes(canon)
report("Q8_planted_clashes_found", c_raw > 0 and c_wrong > 0 and c_right == 0, clashes_without_canonical_form=c_raw,
       clashes_with_wrong_canonical_form=c_wrong, clashes_with_canonical_form=c_right)
OUT["all_pass"] = ok_all
json.dump(OUT, open("checks15.json", "w"), indent=1, default=str)
print("ALL PASS" if ok_all else "SOME FAILED", sum(1 for k, x in OUT.items() if isinstance(x, dict) and x["ok"]), "/", sum(1 for x in OUT.values() if isinstance(x, dict)))
