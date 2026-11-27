"""Round 18 (DRAFT) checks. Usage: python3 checks18.py (needs pairs18.json from pairs18.py). Writes checks18_output.txt and
checks18.json. Each check re-derives a step of Theorem 18.1 (README section 22) independently where it can. About 2 minutes,
peak memory well under 1 GB."""
import json, math, time, itertools
from fractions import Fraction as Fr
from decimal import Decimal, getcontext
from pairs17 import runs, P_of, T2, T3, strip, v, name2, name3, LN
from pairs18 import (prefixes, is_typeB, m2_of, size_of, a_priori, one_round, brute, S_MAX, T3_UP, T3_LO)
from lattice17 import reduce_fast
getcontext().prec = 700
res = []; lines = []
def report(q, ok, msg):
    res.append(dict(check=q, passed=bool(ok), detail=msg)); s = ("PASS " if ok else "FAIL ") + q + ": " + msg
    lines.append(s); print(s, flush=True)
J = json.load(open("pairs18.json"))
Xf = J["final_bound"]
vr = lambda r: 1 if r == 19 else 0          # 19-adic size of r
cr = lambda r: 1 if r == 19 else 10         # r without its 19s

# Q1: F2, by the mod-19 casework of README 22.2 and by direct computation
# terms of T3: 19^(L1+L2) | r3 P2 19^L1 | r3 P2 r2 P1, sizes t1 = L1+L2, t2 = v(r3)+L1, t3 = v(r2)+v(r3); P = 10^L mod 19.
case_bad = []; tie_cases = 0
for r2, r3 in itertools.product((10, 19), repeat=2):
    for L1 in range(0, 8):
        for L2 in range(0, 8):
            t = (L1 + L2, vr(r3) + L1, vr(r2) + vr(r3)); g = min(t)
            u = (1, cr(r3) * pow(10, L2, 19), cr(r2) * cr(r3) * pow(10, L1 + L2, 19))
            tied = [i for i in range(3) if t[i] == g]
            if len(tied) == 1: bound = g
            elif tied == [1, 2]: bound = vr(r3) + 1                    # r3 P2 T_pre with v19(T_pre) <= 1 (F1)
            else:
                tie_cases += 1
                if sum(u[i] for i in tied) % 19 == 0: case_bad.append((r2, r3, L1, L2)); continue
                bound = g
            if bound > 2: case_bad.append((r2, r3, L1, L2, bound))
# any tie that involves t1 forces L1 + L2 <= 2 (t1 = t2 < t3 or t1 = t3 or all three), so L1, L2 < 8 covers every tie
direct_bad = 0; nd = 0; hist = {}
for p in prefixes(30):
    for s in runs(30):
        g = v(T3(p, s), 19); nd += 1; hist[g] = hist.get(g, 0) + 1
        if g > 2: direct_bad += 1
report("Q1 F2: v19(T3) <= 2 for every three-letter word (mod-19 casework, then directly for gaps up to 30)",
       not case_bad and direct_bad == 0,
       "%d tie cases, sums of leading units mod 19 never 0; %d words computed, v19 counts %s" %
       (tie_cases, nd, dict(sorted(hist.items()))))

# Q2: Type A coprimality and Type B's delta
badA = 0; nA = 0
for p in prefixes(25):
    for s in runs(25):
        if is_typeB(s): continue
        nA += 1; t = T3(p, s); tp = T2(p)
        if t % 2 == 0 or t % 5 == 0 or tp % 2 == 0 or tp % 5 == 0: badA += 1
deltas = {}; badB = []
for p in prefixes(80):
    t = T3(p, (19, 0, 0)); e2 = v(t, 2); e5 = v(t, 5); d = 2 ** e2 * 5 ** e5
    deltas[d] = deltas.get(d, 0) + 1
    if d > 40 or e2 > 3 or e5 > 1 or (e2 == 3 and name3(p, (19, 0, 0)) != "NmNN"): badB.append(name3(p, (19, 0, 0)))
report("Q2 Type A: T3 and T_pre prime to 10; Type B: delta = 2^e2 5^e5 <= 40, e2 = 3 only for NmNN",
       badA == 0 and not badB, "%d Type A words (gaps <= 25), %d exceptions; Type B prefixes up to 80: deltas %s" %
       (nA, badA, dict(sorted(deltas.items()))))

# Q3: the three error terms, exactly
worstA = Fr(0); worst2 = Fr(0); worstB = Fr(0); badE = 0
for p in prefixes(20):
    r2, L1, M1 = p; H1 = L1 + M1
    for s in runs(20):
        r3, L2, M2 = s; H2 = L2 + M2; t = T3(p, s)
        if is_typeB(s):
            e = Fr(t, 19 * r2 * P_of(L1, M1)) - 1; q = e * 10 ** H1; worstB = max(worstB, q)
            if not (0 < e and q <= Fr(106, 1000)): badE += 1
        else:
            e = Fr(t, r2 * r3 * P_of(L1, M1) * P_of(L2, M2)) - 1; q = e * 10 ** H1; worstA = max(worstA, q)
            e2 = Fr(t, r3 * P_of(L2, M2) * T2(p)) - 1; q2 = e2 * 10 ** H2; worst2 = max(worst2, q2)
            if not (0 < e and q <= Fr(11, 100) and 0 < e2 and q2 <= Fr(1, 10)): badE += 1
report("Q3 error terms: Type A eps <= 0.11 * 10^-H1, eps2 <= 0.1 * 10^-H2; Type B eps <= 0.106 * 10^-H1 (exact, gaps <= 20)",
       badE == 0, "largest eps 10^H1 = %.5f (A), eps2 10^H2 = %.5f, eps 10^H1 = %.5f (B); %d exceptions" %
       (float(worstA), float(worst2), float(worstB), badE))

# Q4: the size bounds used to cap the other word
badS = 0; nS = 0
for p in prefixes(25):
    for s in runs(25):
        H = p[1] + p[2] + s[1] + s[2]; t = T3(p, s); nS += 1
        Sp = t // strip(t)
        if not (T3_LO * 20 ** H <= t <= T3_UP * 200 ** H and Sp <= S_MAX): badS += 1
report("Q4 100 * 20^(H1+H2) <= T3 <= 401 * 200^(H1+H2) and the {2,5,19}-part of T3 is at most 14440", badS == 0,
       "%d words (gaps <= 25), %d exceptions" % (nS, badS))

# Q5: the prefix lemma and the three non-vanishing facts
pairs = {}; coinc = []
for w in runs(40):
    r, L, M = w; key = Fr(19 ** L, r * P_of(L, M))
    if key in pairs: coinc.append((name2(pairs[key]), name2(w)))
    else: pairs[key] = w
nonvan = {}; dup11 = 0; unit_pre = []
for p in prefixes(14):
    if strip(T2(p)) == 1: unit_pre.append(name2(p))
    for s in runs(14):
        t = T3(p, s); r2, L1, M1 = p; r3, L2, M2 = s
        D1 = 19 * r2 * P_of(L1, M1) if is_typeB(s) else r2 * r3 * P_of(L1, M1) * P_of(L2, M2)
        key = Fr(t, D1)                                    # 1 + eps: Lambda11 = 0 exactly when these agree
        if key in nonvan: dup11 += 1
        else: nonvan[key] = (p, s)
report("Q5 prefix lemma (19^L / (rP) determines a reduced gap and letter), Lambda11 != 0 (1 + eps determines the word) and "
       "T_pre is never a product of 2s, 5s and 19s", not coinc and dup11 == 0 and not unit_pre,
       "%d gap-letter pairs (gaps <= 40), %d coincidences; %d words (gaps <= 14), %d equal 1 + eps; %d unit prefixes" %
       (len(pairs), len(coinc), len(nonvan), dup11, len(unit_pre)))

# Q6: the 2-adic and 5-adic facts of the D-step
bad6 = 0
for D in range(-3000, 3001):
    if D == 0: continue
    a = 19 ** abs(D) - 1
    e2 = 1 if D % 2 else 2 + v(abs(D), 2); e5 = 0 if D % 2 else 1 + v(abs(D), 5)
    if v(a, 2) != e2 or v(a, 5) != e5: bad6 += 1
bad6b = 0
for s in runs(60):
    if is_typeB(s): continue
    r3, L2, M2 = s; q = r3 * P_of(L2, M2)
    if v(q, 2) < 2 * (L2 + M2) or v(q, 5) < max(L2 + M2, 1): bad6b += 1
report("Q6 v2(19^D - 1) = 1 or 2 + v2(D), v5(19^D - 1) = 0 or 1 + v5(D) (|D| <= 3000); v2(r3 P2) >= 2 H2, v5(r3 P2) >= max(H2, 1)",
       bad6 == 0 and bad6b == 0, "%d and %d exceptions" % (bad6, bad6b))

# Q7: the D-step congruence K = 19^(L1+L2-g) mod 2^v2(r3 P2) 5^v5(r3 P2), on every Type A word with gaps <= 18
bad7 = 0; n7 = 0
for p in prefixes(18):
    for s in runs(18):
        if is_typeB(s): continue
        r3, L2, M2 = s; q = r3 * P_of(L2, M2); mod = 2 ** v(q, 2) * 5 ** v(q, 5)
        t = T3(p, s); g = v(t, 19); K = t // 19 ** g; e = p[1] + L2 - g
        rhs = pow(19, e, mod) if e >= 0 else pow(pow(19, -1, mod), -e, mod)
        n7 += 1
        if K % mod != rhs: bad7 += 1
report("Q7 D-step congruence on Type A words", bad7 == 0, "%d words (gaps <= 18), %d exceptions" % (n7, bad7))

# Q8: the bound chain, reproduced from scratch and compared with pairs18.json
t0 = time.time(); ap = a_priori(); rr = []; X = ap["X_star"]
while True:
    Xn = one_round(X, rr)
    if Xn >= X: break
    X = Xn
keys = ["X", "m11", "m21", "m2", "short_word_sum_s", "other_word_sum", "new_X"]
same = ap["X_star"] == J["a_priori"]["X_star"] and X == Xf and len(rr) == len(J["rounds"]) and \
       all(a[k] == b[k] for a, b in zip(rr, J["rounds"]) for k in keys)
report("Q8 bound chain reproduced (Matveev start, lattice rounds, final bound)", same,
       "X* = %.4g, then %s, final bound %d (%.0f s)" % (ap["X_star"], " -> ".join(str(r["new_X"]) for r in rr), X, time.time() - t0))

# Q9: the lattice step cannot miss a genuine near-relation (planted test on the same code path)
# a brute minimum of |c log 19 + a log 2 + b log 5| over a large slice of the box |a|, |b|, |c| <= 60 must stay above the bound
X0 = 60; C = 10 ** (3 * len(str(X0)) + 6)
lam, _ = reduce_fast(X0, C, [LN[19], LN[2], LN[5]])
l19, l2, l5 = math.log(19), math.log(2), math.log(5)
best = min(abs(c * l19 + a * l2 + b * l5) for c in range(-X0, X0 + 1) for a in range(-X0, X0 + 1) for b in range(-X0, X0 + 1)
           if (a, b, c) != (0, 0, 0) and abs(c * l19 + a * l2 + b * l5) < 0.5)
small = best
report("Q9 the lattice lower bound is a true lower bound (|x| <= 60: the brute minimum over the whole box is above it)",
       lam is not None and lam <= best and lam <= small,
       "lattice bound %.3e, smallest found %.3e" % (float(lam), float(min(best, small))))

# Q10: the exact comparison, by the T3 formula and independently from the letters as maps
ex = brute(Xf)
S = {"l": Fr(19, 200), "m": Fr(1, 20), "N": Fr(1, 19), "K": Fr(1, 10)}
Dd = {"l": Fr(0), "m": Fr(0), "N": Fr(1, 19), "K": Fr(1, 10)}
def psi_maps(word):
    s, t = Fr(1), Fr(0)
    for c in word: s, t = s * S[c], s * Dd[c] + t
    return t / s
gaps = [""] + ["".join(x) for k in range(1, 6) for x in itertools.product("lm", repeat=k)]
byk = {}; nraw = 0
for g1 in gaps:
    for x2 in "NK":
        if g1 == "" and x2 == "N": continue          # N N ... has only two N/K letters (NN = Nm)
        for g2 in gaps:
            for x3 in "NK":
                w = "N" + g1 + x2 + g2 + x3; ps = psi_maps(w); nraw += 1
                k = Fr(strip(ps.numerator), strip(ps.denominator))
                byk.setdefault(k, set()).add(ps)
raw_bad = sum(1 for k in byk if len(byk[k]) > 1)
report("Q10 every reduced three-letter word with H1 + H2 <= final bound has its own fingerprint; mK words found as their "
       "reduced forms; raw words from the maps: equal fingerprint means equal psi",
       not ex["clashes"] and ex["planted_found"] == ex["planted_mK_words"] and not ex["planted_matched_wrong_word"] and raw_bad == 0,
       "%d words up to H1 + H2 = %d, %d distinct; %d/%d planted found; %d raw words (gaps <= 5 letters), %d fingerprints, %d with two psi" %
       (ex["three_letter_words"], Xf, ex["distinct_fingerprints"], ex["planted_found"], ex["planted_mK_words"], nraw, len(byk), raw_bad))

n_ok = sum(r["passed"] for r in res)
lines.append("%d/%d checks pass" % (n_ok, len(res))); print(lines[-1])
open("checks18_output.txt", "w").write("\n".join(lines) + "\n")
json.dump(dict(passed=n_ok, total=len(res), checks=res), open("checks18.json", "w"), indent=1)
