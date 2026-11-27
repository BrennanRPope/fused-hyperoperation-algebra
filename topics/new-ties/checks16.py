"""Round 16 checks (DRAFT). Reads baker16.json and pairs16.json; reruns small parts independently. Writes checks16.json.
Q1  baker16.json: Matveev bounds (about 5 x 10^12) recomputed from the formula; LLL rounds end at E <= 10 (Case II) and
    H <= 7 (Case I); exact check below them finds no solution; brute force on all 3-letter words with runs up to 12 finds no
    unit psi
Q2  the reduction in Theorem 16.1 on 300 random 3-letter words, exact: psi = 1 + r3 A2 (1 + r2 A1) agrees with the maps;
    19^(L1+L2) psi = 19^(L1+L2) + r3 19^L1 P2 + r2 r3 P1 P2; outside Case I, psi - 1 has positive 2- and 5-adic size;
    in Case I, psi = 20 + 19 r2 A1
Q3  the valuation facts used: v2(19^E - 1) = 1 (E odd) or 2 + v2(E) (E even); v5(19^E - 1) = 0 (E odd) or 1 + v5(E)
    (E even), for E <= 3000; and v19(19^L1 + r2 P1) is 0 or 1 for L1 >= 2 (so L1 <= 2), for L1 <= 40, M1 <= 40
Q4  the LLL step is sound: every reduced basis passes an exact LLL-reducedness test, and repeating Case II's first
    reduction with a 1,000 times larger C again gives a bound below 40
Q5  pairs16.json: 526 small words, every final bound <= 39, and no two of the 860 reduced two-letter words with H <= 39
    share a fingerprint; also no clash among all 2,627 such words with H <= 70 (a margin past the proof's bound)
Q6  the 2-adic step of Theorem 16.2 on 300 random pairs of two-letter words (H <= 60), exact: K = 19^(L-g) mod 2^mu, and
    whenever D != 0, min(mu, mu') <= v2(19^|D| - 1) holds for the congruence 19^D = 1 mod 2^min(mu, mu')
Q7  the exact check catches a planted clash: with non-reduced words (K after m) allowed, the two-letter check finds the
    mK ~ lN pairs at once; and a toy lattice comes back LLL-reduced
"""
import json, random, sys, os, math
from fractions import Fraction as F
from decimal import Decimal
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import baker16 as B
from pairs16 import strip, T_of, words as words2, name
OUT = {}; ok_all = True
def report(key, ok, **kw):
    global ok_all
    ok_all &= bool(ok); OUT[key] = dict(ok=bool(ok), **kw); print(key, "PASS" if ok else "FAIL", kw, flush=True)
def v(n, p):
    n = F(n); a, b, k = n.numerator, n.denominator, 0
    while a % p == 0: a //= p; k += 1
    while b % p == 0: b //= p; k -= 1
    return k
LET = {"l": (F(19, 200), F(0)), "m": (F(1, 20), F(0)), "N": (F(1, 19), F(1, 19)), "K": (F(1, 10), F(1, 10))}
def wmap(w):
    s, t = F(1), F(0)
    for x in w: a, b = LET[x]; s, t = s * a, s * b + t
    return s, t
def psi(w): s, t = wmap(w); return t / s
rng = random.Random(16)
# Q1
J = json.load(open("baker16.json"))
CM = B.matveev_constant()
BE = lambda E: math.ceil(4.25 * (E + 4 + 2 * math.log(max(E, 1), 5)))
E0 = B.largest_ok(lambda E: E * math.log(19) - math.log(2 * (1 + 380 * E * E)) < CM * (1 + math.log(BE(E))), 4, 10 ** 16)
H0 = B.largest_ok(lambda H: H * math.log(10) + math.log(9) < CM * (1 + math.log(3 * H + 5)), 1, 10 ** 16)
ok = (J["matveev_bounds"]["case_II_E"] == E0 and J["matveev_bounds"]["case_I_H"] == H0 and J["case_II_final_bound"] <= 10
      and J["case_I_final_bound"] <= 7 and J["exact_check"]["solutions"] == [] and J["brute_force_runs_up_to_12"]["unit_psi"] == [])
report("Q1_baker16_record", ok, matveev_E=E0, matveev_H=H0, final_E=J["case_II_final_bound"], final_H=J["case_I_final_bound"],
       exact_solutions=J["exact_check"]["solutions"], brute_force_words=J["brute_force_runs_up_to_12"]["words"])
# Q2
bad = 0
for _ in range(300):
    r2, r3 = rng.choice([19, 10]), rng.choice([19, 10]); L1, M1, L2, M2 = (rng.randint(0, 6) for _ in range(4))
    w = "N" + "l" * L1 + "m" * M1 + ("N" if r2 == 19 else "K") + "l" * L2 + "m" * M2 + ("N" if r3 == 19 else "K")
    p = psi(w); A1 = F(200, 19) ** L1 * 20 ** M1; A2 = F(200, 19) ** L2 * 20 ** M2
    if p != 1 + r3 * A2 * (1 + r2 * A1): bad += 1
    P1 = 2 ** (3 * L1 + 2 * M1) * 5 ** (2 * L1 + M1); P2 = 2 ** (3 * L2 + 2 * M2) * 5 ** (2 * L2 + M2)
    if 19 ** (L1 + L2) * p != 19 ** (L1 + L2) + r3 * 19 ** L1 * P2 + r2 * r3 * P1 * P2: bad += 1
    if r3 == 19 and L2 == M2 == 0:
        if p != 20 + 19 * r2 * A1: bad += 1
    elif not (v(p - 1, 2) >= 1 and v(p - 1, 5) >= 1): bad += 1
report("Q2_reduction_16_1", bad == 0, words=300, failures=bad)
# Q3
bad = 0
for E in range(1, 3001):
    D = 19 ** E - 1; e2 = 1 if E % 2 else 2 + v(E, 2); e5 = 0 if E % 2 else 1 + v(E, 5)
    if v(D, 2) != e2 or v(D, 5) != e5: bad += 1
bad19 = 0
for r2 in (19, 10):
    for L1 in range(2, 41):
        for M1 in range(41):
            if v(19 ** L1 + r2 * 2 ** (3 * L1 + 2 * M1) * 5 ** (2 * L1 + M1), 19) > 1: bad19 += 1
report("Q3_valuation_facts", bad == 0 and bad19 == 0, E_checked=3000, failures=bad, v19_failures=bad19)
# Q4
X0 = BE(J["matveev_bounds"]["case_II_E"]); C = 10 ** (3 * len(str(X0)) + 9)
lam, b1sq = B.reduce_once(X0, C)
newE = B.largest_ok(lambda E: Decimal(E) * B.LOG[19] - Decimal(2 * (1 + 380 * E * E)).ln() < -lam.ln(), 1, J["matveev_bounds"]["case_II_E"])
report("Q4_lll_sound", lam is not None and newE < 40, C_digits=len(str(C)) - 1, new_bound=newE, note="reduce_once asserts exact LLL-reducedness")
# Q5
P = json.load(open("pairs16.json"))
seen = {}; cl = 0; n = 0
for w in words2(70):
    n += 1; K = strip(T_of(w))
    if K in seen: cl += 1
    else: seen[K] = w
ok = P["small_words"] == 526 and P["largest_final_bound"] <= 39 and P["exact_check"]["clashes"] == [] and P["exact_check"]["H_up_to"] >= 39 and cl == 0
report("Q5_pairs16_record", ok, small_words=P["small_words"], largest_final_bound=P["largest_final_bound"],
       words_checked_in_proof=P["exact_check"]["words"], extra_words_H_le_70=n, extra_clashes=cl)
# Q6
bad = 0; tested = 0
for _ in range(300):
    ws = []
    for _ in range(2):
        r = rng.choice([19, 10]); L = rng.randint(0, 30); M = 0 if r == 10 else rng.randint(0, 30)
        if (r, L, M) == (19, 0, 0): L = 1
        ws.append((r, L, M))
    vals = []
    for (r, L, M) in ws:
        T = T_of((r, L, M)); g = v(T, 19); K = strip(T); mu = v(r * 2 ** (3 * L + 2 * M) * 5 ** (2 * L + M), 2)
        if (19 ** g * K - 19 ** L) % 2 ** mu: bad += 1          # 19^g K = 19^L mod 2^mu
        vals.append((L - g, mu))
    D = vals[0][0] - vals[1][0]; m = min(vals[0][1], vals[1][1])
    if D != 0:
        tested += 1
        lhs_ok = (pow(19, abs(D), 2 ** m) == 1)
        if lhs_ok and m > v(19 ** abs(D) - 1, 2): bad += 1
report("Q6_two_adic_step", bad == 0, pairs=300, pairs_with_D_nonzero=tested, failures=bad)
# Q7
seen = {}; planted = 0
for H in range(12):
    for L in range(H + 1):
        M = H - L
        for r in (19, 10):
            K = strip(T_of((r, L, M)))
            if K in seen and seen[K] != (r, L, M): planted += 1
            else: seen[K] = (r, L, M)
red = B.lll([[1, 0, 10 ** 12], [0, 1, 3 * 10 ** 12 + 7], [0, 0, 7 * 10 ** 12 + 3]])
report("Q7_planted_clash_found", planted > 0 and B.is_lll_reduced(red), planted_mK_lN_clashes=planted, lll_toy_reduced=B.is_lll_reduced(red))
OUT["all_pass"] = ok_all
json.dump(OUT, open("checks16.json", "w"), indent=1, default=str)
print("ALL PASS" if ok_all else "SOME FAILED", sum(1 for x in OUT.values() if isinstance(x, dict) and x["ok"]), "/", sum(1 for x in OUT.values() if isinstance(x, dict)))
