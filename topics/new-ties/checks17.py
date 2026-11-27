"""Round 17 (DRAFT) checks. Usage: python3 checks17.py (needs pairs17.json from pairs17.py). Writes checks17_output.txt and
checks17.json. Each check re-derives a step of Theorem 17.1 (README section 21) independently where it can."""
import json, math, random, time
from fractions import Fraction as Fr
from decimal import Decimal, getcontext
from pairs17 import (runs, P_of, T2, T3, strip, v, jmax, two_letter_words, three_letter_words, brute, lattice_bound,
                     m2_of, hW_of, one_round, LN, name2, name3)
from baker16 import matveev_constant, largest_ok
getcontext().prec = 700
random.seed(17)
res = []; lines = []
def report(q, ok, msg):
    res.append(dict(check=q, passed=bool(ok), detail=msg)); s = ("PASS " if ok else "FAIL ") + q + ": " + msg
    lines.append(s); print(s, flush=True)
J = json.load(open("pairs17.json"))
Xf = J["final_bound"]

# Q1: F1 and coprimality for two-letter words
bad = [w for w in runs(200) if w != (19, 0, 0) and (v(T2(w), 19) > 1 or T2(w) % 2 == 0 or T2(w) % 5 == 0)]
report("Q1 two-letter words: 19-adic size g' <= 1 and T' prime to 10 (all reduced words but NN, gaps up to 200)", not bad,
       "%d words, %d exceptions" % (sum(1 for _ in runs(200)) - 1, len(bad)))

# Q2: psi formula against the letters as maps; Type A coprimality; Type B 2- and 5-adic sizes
S = {"l": Fr(19, 200), "m": Fr(1, 20), "N": Fr(1, 19), "K": Fr(1, 10)}
Dd = {"l": Fr(0), "m": Fr(0), "N": Fr(1), "K": Fr(1)}
def psi_maps(word):
    ps = Dd[word[0]]
    for c in word[1:]: ps = ps / S[c] + Dd[c]
    return ps
mis = 0
for _ in range(400):
    p = random.choice([w for w in runs(8) if w != (19, 0, 0)]); s = random.choice(list(runs(8)))
    if psi_maps(name3(p, s)) != Fr(T3(p, s), 19 ** (p[1] + s[1])): mis += 1
badA = 0; badB = []; nA = nB = 0
for p, s in three_letter_words(25, 25):
    t = T3(p, s)
    if s == (19, 0, 0):
        nB += 1
        if (v(t, 2) > 2 or v(t, 5) > 1) and name3(p, s) != "NmNN": badB.append(name3(p, s))
    else:
        nA += 1; tp = T2(p)
        if t % 2 == 0 or t % 5 == 0 or tp % 2 == 0 or tp % 5 == 0: badA += 1
report("Q2 psi = T3 / 19^(L1+L2) matches the letters as maps; Type A: T3, T_pre prime to 10; Type B: at most 2 twos and 1 five except NmNN",
       mis == 0 and badA == 0 and not badB, "400 random words, %d mismatches; %d Type A words, %d bad; %d Type B words, exceptions %s"
       % (mis, nA, badA, nB, badB))

# Q3: the small quantities and coefficient sizes used in (L3), (L1), (L2), (LB), (LB1), (2ad), (D)
worst = {"eps3": 0, "eps2": 0, "epsB": 0, "epsV": 0, "g3": 0, "|m|": 0}; okq = True
for p, s in three_letter_words(20, 20):
    (r2, L1, M1), (r3, L2, M2) = p, s; H1, H2 = L1 + M1, L2 + M2; X = max(H1, H2, 1); t = T3(p, s); g3 = v(t, 19)
    P1, P2 = P_of(L1, M1), P_of(L2, M2); tp = T2(p)
    if g3 > 2.02 + 1.8 * (H1 + H2): okq = False
    if abs(g3 - L1 - L2) > 4 * X + 3: okq = False
    worst["|m|"] = max(worst["|m|"], abs(g3 - L1 - L2))
    if s == (19, 0, 0):
        e = Fr(20 * 19 ** L1, 19 * r2 * P1) * 10 ** H1; worst["epsB"] = max(worst["epsB"], float(e)); okq &= e < Fr(106, 1000)
    else:
        e3 = Fr(19 ** L1 * (19 ** L2 + r3 * P2), r2 * r3 * P1 * P2) * 10 ** H1
        e2 = Fr(19 ** (L1 + L2), r3 * P2 * tp) * 10 ** H2
        worst["eps3"] = max(worst["eps3"], float(e3)); worst["eps2"] = max(worst["eps2"], float(e2))
        okq &= e3 <= Fr(11, 100) and e2 < Fr(1, 10)
    worst["g3"] = max(worst["g3"], g3)
for w in two_letter_words(60):
    r, L, M = w; e = Fr(19 ** L, r * P_of(L, M)) * 10 ** (L + M); worst["epsV"] = max(worst["epsV"], float(e)); okq &= e <= Fr(1, 10)
report("Q3 small terms and sizes (all words with gaps up to 20; V up to 60): eps3 <= 0.11 10^-H1 (equality only for NKK; log(1 + eps) < eps keeps every Lambda strictly below), eps2 < 0.1 10^-H2, epsB < 0.106 10^-H1, "
       "eps' <= 0.1 10^-H', g3 <= 2.02 + 1.8 (H1+H2), |g3 - L1 - L2| <= 4X + 3", okq,
       "largest scaled values: " + ", ".join("%s %.4g" % kv for kv in worst.items()))

# Q4: the three non-vanishing steps
su = [w for w in runs(300) if strip(T2(w)) == 1]
eqsol = [(c, cp, L1, Lp, M1, Mp) for c in (0, 1) for cp in (0, 1) for L1 in range(40) for Lp in range(40)
         for M1 in range(40) for Mp in range(40)
         if c + 3 * L1 + 2 * M1 == cp + 3 * Lp + 2 * Mp + 2 and c + 2 * L1 + M1 == cp + 2 * Lp + Mp + 1 and L1 - Lp == 1 - c + cp]
report("Q4 non-vanishing: 19^L + rP is a product of 2s, 5s, 19s only for NN (gaps up to 300; Theorem 15.5 proves it), and the "
       "three equations behind Lambda_B = 0 have no solution", [name2(w) for w in su] == ["NN"] and not eqsol,
       "S-unit T's: %s; solutions of the Type B system with values below 40: %d (the hand proof: Eq1 - 2 Eq2 gives L1 - L' = c' - c, "
       "Eq3 gives L1 - L' = 1 - c + c')" % ([name2(w) for w in su], len(eqsol)))

# Q5: the 2-adic steps
bad5 = 0
for D in range(1, 20001):
    if v(19 ** D - 1, 2) > 2 + math.log2(D) + 1e-12: bad5 += 1
badc = 0
for p, s in three_letter_words(14, 14):
    if s == (19, 0, 0): continue
    t = T3(p, s); g3 = v(t, 19); K = strip(t); nu = v(s[0] * P_of(s[1], s[2]), 2); mm = g3 - p[1] - s[1]
    lhs = K * 19 ** mm - 1 if mm >= 0 else K - 19 ** (-mm)
    if lhs != 0 and v(abs(lhs), 2) < nu: badc += 1
badj = 0; eq = 0
for _ in range(300):
    K = strip(random.randrange(3, 10 ** random.randint(2, 30)) | 1)
    if K == 1: continue
    M = random.randint(0, 400); best = 0
    for m in range(-M, M + 1):
        val = K * 19 ** m - 1 if m >= 0 else K - 19 ** (-m)
        best = max(best, v(abs(val), 2))
    j = jmax(K, M)
    if j < best: badj += 1
    if j == best: eq += 1
report("Q5 2-adic: v2(19^D - 1) <= 2 + log2 D (D <= 20,000); K 19^(g3-L1-L2) = 1 mod 2^v2(r3 P2) (Type A, gaps up to 14); exact "
       "2-adic log bound >= the true maximum (300 random K, M <= 400)", bad5 == 0 and badc == 0 and badj == 0,
       "%d, %d and %d failures; the 2-adic bound was exactly the true maximum in %d of 300" % (bad5, badc, badj, eq))

# Q6: the bound chain: a priori number, then every reduction round except the first redone here; first round's universal step redone
R = J["rounds"]; okr = True; redo = []
m1, _, _ = lattice_bound(6 * R[0]["X"] + 6, [LN[19], LN[2], LN[5]], "0.11")
okr &= min(m1, R[0]["X"]) == R[0]["m1"]
for rec in R[1:]:
    tmp = []; Xn = one_round(rec["X"], tmp); redo.append(tmp[0]["new_X"]); okr &= (Xn == rec["new_X"])
C3 = matveev_constant(n=3, alphas=(19, 2, 5))
report("Q6 bound chain reproduced (a priori X* = %.3e; rounds %s; final %d)" % (J["a_priori"]["X_star"], [r["new_X"] for r in R], Xf),
       okr and Xf == R[-1]["X"] and abs(C3 - J["a_priori"]["C3"]) < 1, "first round's universal lattice step and rounds 2-%d redone: %s"
       % (len(R), redo))

# Q7: the exact check below the final bound
B = brute(Xf)
report("Q7 exact comparison below the final bound H <= %d" % Xf, not B["clashes"] and B["planted_found"] == B["planted_NN_words"]
       and B == J["exact_check"], "%d three-letter words against %d two-letter words: %d clashes; planted N N u2 x3 words found %d of %d"
       % (B["three_letter_words"], B["two_letter_words"], len(B["clashes"]), B["planted_found"], B["planted_NN_words"]))

# Q8: an independent search with fingerprints computed from the letters as maps (Fractions), gaps up to 9 against 2-letter up to 40
def kap(q): return strip(q.numerator), strip(q.denominator)
Vk = {kap(psi_maps(name2(w))) for w in two_letter_words(40)}
hits = []; n8 = 0
for p, s in three_letter_words(9, 9):
    n8 += 1
    if kap(psi_maps(name3(p, s))) in Vk: hits.append(name3(p, s))
report("Q8 independent search from the maps (three-letter gaps up to 9, two-letter gaps up to 40)", not hits, "%d words, %d clashes"
       % (n8, len(hits)))

# Q9: Round 16 under the logarithmic-form correction (Matveev bounds |prod - 1|; |log(1+y)| <= 2|y| costs log 2)
from baker16 import matveev_constant as mc
CM = mc(); l2 = math.log(2)
BE = lambda E: math.ceil(4.25 * (E + 4 + 2 * math.log(max(E, 1), 5)))
E0 = largest_ok(lambda E: E * math.log(19) - math.log(2 * (1 + 380 * E * E)) < CM * (1 + math.log(BE(E))) + l2, 4, 10 ** 16)
H0 = largest_ok(lambda H: H * math.log(10) + math.log(9) < CM * (1 + math.log(3 * H + 5)) + l2, 1, 10 ** 16)
from lattice17 import reduce_fast
def chain(start, X0of, ok):
    cur = start
    for _ in range(8):
        X0 = X0of(cur); C = 10 ** (3 * len(str(X0)) + 6); lam, _ = reduce_fast(X0, C, [LN[19], LN[2], LN[5]])
        if lam is None or lam <= 0: C *= 10 ** 6; lam, _ = reduce_fast(X0, C, [LN[19], LN[2], LN[5]])
        nl = -lam.ln(); new = largest_ok(lambda x: ok(x, nl), 1, cur)
        if new >= cur: break
        cur = new
    return cur
eII = chain(E0, BE, lambda E, nl: Decimal(E) * LN[19] - Decimal(2 * (1 + 380 * E * E)).ln() < nl)
hI = chain(H0, lambda H: 3 * H + 5, lambda H, nl: Decimal(H) * LN[10] + Decimal(9).ln() < nl)
C4 = mc(n=4, alphas=(19, 2, 5)); lkb = lambda H: 3.64 + 5.30 * (1 + math.log2(3 * H + 2) / 2)
H0p = largest_ok(lambda H: H * math.log(10) + math.log(10) < C4 * (1 + math.log(3 * H + 3)) * lkb(H) + l2, 3, 10 ** 20)
Hs = math.floor(1 + math.log2(3 * H0p + 2) / 2); worst16 = 0
from pairs16 import words as w16, T_of as T16, strip as s16
for sw in w16(Hs):
    K = s16(T16(sw))
    if K == 1: continue
    cur = H0p; lK = Decimal(K).ln()
    for _ in range(8):
        X0 = 3 * cur + 3; C = 10 ** (4 * len(str(X0)) + 8); lam, _ = reduce_fast(X0, C, [LN[19], LN[2], LN[5], lK])
        t = 0
        while (lam is None or lam <= 0) and t < 4: C *= 10 ** 6; t += 1; lam, _ = reduce_fast(X0, C, [LN[19], LN[2], LN[5], lK])
        if lam is None or lam <= 0: break
        nl = -lam.ln(); new = largest_ok(lambda H: Decimal(H + 1) * LN[10] < nl, 0, cur)
        if new >= cur: break
        cur = new
    worst16 = max(worst16, cur)
report("Q9 Round 16 survives the log-form correction (the Matveev starts move by at most 1 and the lattice steps land inside the "
       "ranges Round 16 checked exactly)", eII <= 10 and hI <= 7 and Hs <= 30 and worst16 <= 39,
       "corrected starts E <= %d, H <= %d, H <= %.4g; after LLL: E <= %d (Round 16 checked 10), H <= %d (7), short words H <= %d (30), "
       "long word H <= %d (39)" % (E0, H0, H0p, eII, hI, Hs, worst16))

n_ok = sum(r["passed"] for r in res)
lines.append("%d/%d checks pass" % (n_ok, len(res))); print(lines[-1])
open("checks17_output.txt", "w").write("\n".join(lines) + "\n")
json.dump(dict(passed=n_ok, total=len(res), checks=res), open("checks17.json", "w"), indent=1)
