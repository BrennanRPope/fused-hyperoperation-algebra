"""Round 19 (DRAFT) checks. Usage: python3 checks19.py (needs pairs19.json from pairs19.py). Writes checks19_output.txt and
checks19.json. Each check re-derives a step of Theorems 19.1 to 19.3 (README section 23) independently where it can.
About 5 minutes, peak memory well under 1 GB."""
import json, math, time, itertools, random
from fractions import Fraction as Fr
from decimal import Decimal, getcontext
from pairs17 import strip, v, P_of, LN
from lattice17 import reduce_fast
import pairs19 as p19
from pairs19 import (items_gap, items_upto, T_of, Q_of, K_of, name, hsum, S_MAX, jmax_ratio, sdata, S_bound, flog,
                     size_other, box_words)
getcontext().prec = 700
res = []; lines = []
def report(q, ok, msg):
    res.append(dict(check=q, passed=bool(ok), detail=msg)); s = ("PASS " if ok else "FAIL ") + q + ": " + msg
    lines.append(s); print(s, flush=True)
J = json.load(open("pairs19.json"))
vr = lambda r: 1 if r == 19 else 0
cr = lambda r: 1 if r == 19 else 10
random.seed(19)

def words(n, B, first=True):
    """all words with n items and gap sum <= B"""
    for w in box_words(n, (B,) * n):
        if hsum(w) <= B: yield w

# Q1: F3 by the mod-19 casework of README 23.2, then directly; and the same casework fails for five letters
def casework(n, Lmax):
    """T = sum_j 19^(A_j) prod_(i > j) r_i P_i. Term j has 19-adic size t_j = A_j + #(N after j) and leading unit
    prod_(i > j) c_i 10^(L_i) mod 19 (P = 10^L mod 19). Returns tie patterns and those whose tied units add to 0 mod 19."""
    ties = 0; bad = []
    for rs in itertools.product((19, 10), repeat=n):
        for Ls in itertools.product(range(Lmax + 1), repeat=n):
            A = [0] + list(itertools.accumulate(Ls))
            t = [A[j] + sum(vr(rs[i]) for i in range(j, n)) for j in range(n + 1)]
            u = []
            for j in range(n + 1):
                x = 1
                for i in range(j, n): x = x * cr(rs[i]) * pow(10, Ls[i], 19) % 19
                u.append(x)
            g = min(t); tied = [j for j in range(n + 1) if t[j] == g]
            if len(tied) > 1:
                ties += 1
                if sum(u[j] for j in tied) % 19 == 0: bad.append((rs, Ls))
    return ties, bad
ties4, bad4 = casework(3, 4)          # L = 4 stands for "4 or more": a term with such an L before it has size >= 4 > 3 >= min
ties5, bad5 = casework(4, 5)
hist = {}
for w in words(3, 14):
    g = v(T_of(w), 19); hist[g] = hist.get(g, 0) + 1
fam = {M: v(T_of(((10, 0, 0), (10, 1, 0), (10, 0, 0), (19, 0, M))), 19) for M in (0, 38, 5092)}
report("Q1 F3: v19(T) <= 3 for every four-letter word (casework, then directly); five letters: cancellations exist and "
       "v19 grows", not bad4 and max(hist) <= 3 and len(bad5) > 0 and fam == {0: 3, 38: 4, 5092: 5},
       "four letters: %d tie patterns, tied units never add to 0 mod 19; %d words with gap sum <= 14: v19 counts %s; "
       "five letters: %d of %d tie patterns cancel, e.g. %s; N K l K K m^M N has v19 = %s for M = 0, 38, 5092" %
       (ties4, sum(hist.values()), dict(sorted(hist.items())), len(bad5), ties5,
        name(tuple((r, L, 0) for r, L in zip(*bad5[0]))), [fam[0], fam[38], fam[5092]]))

# Q2: the 2- and 5-parts of T for four letters: 2^e2 5^e5 <= 80, only N K m N N reaches 80, none unless x3 = N N
dmax = 0; at = []; bad2 = 0; n2 = 0
for w in words(3, 16):
    T = T_of(w); e2, e5 = v(T, 2), v(T, 5); d = 2 ** e2 * 5 ** e5; n2 += 1
    if (e2 or e5) and w[2] != (19, 0, 0): bad2 += 1
    if e2 > 4 or e5 > 1: bad2 += 1
    if d > dmax: dmax, at = d, [name(w)]
    elif d == dmax: at.append(name(w))
report("Q2 2^e2 5^e5 <= 80 for four letters (S <= 80 * 19^3 = 548,720)", bad2 == 0 and dmax == 80 and at == ["NKmNN"]
       and S_MAX[4] == 548720, "%d words (gap sum <= 16), %d exceptions; largest %d, at %s" % (n2, bad2, dmax, at))

# Q3: the splitting identity, the size of eta, and the size bounds, on random words with 1 to 4 items
bad3 = 0; n3 = 0; worst = 0.0
for _ in range(4000):
    n = random.randint(1, 3); w = []
    for k in range(n):
        while True:
            it = random.choice(items_gap(random.randint(0, 25)))
            if not (k == 0 and it == (19, 0, 0)): break
        w.append(it)
    w = tuple(w); T = T_of(w); s = hsum(w)
    if not (10 ** n * 20 ** s <= T <= Fr(112, 100) * 19 ** n * 200 ** s): bad3 += 1
    for p in range(n):                                           # 0-based pivot p
        pre, sig = w[:p], w[p + 1:]; Ap = sum(L for _, L, _ in w[:p + 1]); Qp = Q_of(w[p:])
        if T != Qp * T_of(pre) + 19 ** Ap * T_of(sig): bad3 += 1
        eta = Fr(19 ** Ap * T_of(sig), Qp * T_of(pre)); hp = w[p][1] + w[p][2]
        r = float(eta) * 10 ** hp / 0.112 if hp < 300 else 0.0
        worst = max(worst, r); n3 += 1
        if not (0 < eta and math.log1p(float(eta)) < 0.112 * 10.0 ** (-hp) * (1 + 1e-12)): bad3 += 1
report("Q3 T = Q_p T(prefix) + 19^A_p T(suffix); 0 < log(1 + eta) < 0.112 * 10^-h_p; 10^n 20^(sum h) <= T <= 1.12 * 19^n * 200^(sum h)",
       bad3 == 0, "%d splits of 4000 random words, %d exceptions; largest eta * 10^h / 0.112 = %.3f" % (n3, bad3, worst))

# Q4: the suffix-step congruence, word by word: K = 19^(A_p + g_sigma - g) kappa(sigma) mod 2^(v2 Q_p - e) when v2 Q_p > e, v5 Q_p > f
bad4c = 0; n4 = 0; nsk = 0
for nn in (1, 2, 3):
    for w in words(nn, {1: 40, 2: 18, 3: 11}[nn]):
        T = T_of(w); g = v(T, 19); K = strip(T)
        for p in range(nn):
            sig = w[p + 1:]; Ts = T_of(sig); e, f, gs = v(Ts, 2), v(Ts, 5), v(Ts, 19); ks = strip(Ts)
            Qp = Q_of(w[p:]); hp = w[p][1] + w[p][2]; rp = w[p][0] * P_of(w[p][1], w[p][2])
            if v(rp, 2) < 2 * hp or v(rp, 5) < hp: bad4c += 1
            if v(Qp, 2) <= e or v(Qp, 5) <= f: nsk += 1; continue
            N2 = v(Qp, 2) - e; mod = 1 << N2; ex = sum(L for _, L, _ in w[:p + 1]) + gs - g
            rhs = (pow(19, ex, mod) if ex >= 0 else pow(pow(19, -1, mod), -ex, mod)) * ks % mod
            n4 += 1
            if K % mod != rhs or v(T, 2) != e or v(T, 5) != f: bad4c += 1
report("Q4 suffix-step congruence (and e2(T) = e, e5(T) = f) at every split; v2(r_p P_p) >= 2 h_p, v5(r_p P_p) >= h_p",
       bad4c == 0, "%d splits tested (%d skipped as trivial cases), %d exceptions" % (n4, nsk, bad4c))

# Q5: the exact 2-adic step against brute force
bad5j = 0; n5 = 0
kaps = sorted({strip(T_of(w)) for w in words(2, 6)} | {1})
pairs = [(a, b) for a in kaps[:60] for b in kaps[:60] if a != b] + [(random.randrange(1, 10 ** 30, 2), 3) for _ in range(200)]
for a, b in pairs:
    if a % 19 == 0 or b % 19 == 0: continue
    for M in (1, 7, 50, 300):
        mod = 1 << 80; rho = a * pow(b, -1, mod) % mod
        br = max(v((pow(19, m, mod) if m >= 0 else pow(pow(19, -1, mod), -m, mod)) * rho % mod - 1 or mod, 2)
                 for m in range(-M, M + 1))
        n5 += 1
        if jmax_ratio(a, b, M) != br: bad5j += 1
report("Q5 the exact 2-adic discrete logarithm step equals brute force (|m| <= M)", bad5j == 0,
       "%d (ratio, M) cases, %d exceptions" % (n5, bad5j))

# Q6: Lemma N's ingredients: the prefix lemma, and equal-kappa suffixes are equal or differ by one N N (ratio 20/19)
vecs = {}; badp = 0
for it in items_upto(60):
    r, L, M = it; key = (L - vr(r), -(3 * L + 2 * M + (r == 10)), -(2 * L + M + (r == 10)))   # 19^L / (r P) as 19, 2, 5 exponents
    if key in vecs: badp += 1
    vecs[key] = it
lists = [()] + [(t,) for t in items_upto(9)] + [(s, t) for s in items_upto(6) for t in items_upto(6)]
groups = {}
for s in lists: groups.setdefault(K_of(s), []).append(s)
badm = 0; nm = 0; seen = {}
for k, g in groups.items():
    for s1, s2 in itertools.combinations(g, 2):
        mu = Fr(T_of(s1), Q_of(s1)) / Fr(T_of(s2), Q_of(s2)); nm += 1
        seen[str(mu)] = seen.get(str(mu), 0) + 1
        if mu not in (Fr(20, 19), Fr(19, 20)): badm += 1
iota = lambda a2, a5, a19: (a2 - 2 * a5 + a19) % 2
odd = all(iota(k[1], k[2], k[0]) == 1 for k in vecs) and iota(2, 1, -1) == 1
report("Q6 Lemma N: 19^L/(rP) determines the item; different suffixes with one fingerprint differ by 20/19; parity",
       badp == 0 and badm == 0 and odd,
       "%d items, %d collisions; %d suffix lists, %d equal-fingerprint pairs, ratios %s; the parity a2 - 2 a5 + a19 is odd "
       "for every 19^L/(rP) and for 20/19" % (len(vecs), badp, len(lists), nm, seen))

# Q7: the a priori bounds (closed forms in ell = 1 + log(9X + 20)) recomputed, and the dominance step tested
ok7 = True; msg7 = []
for key, fn in (("42", p19.a_priori_42), ("43", p19.a_priori_43), ("44", p19.a_priori_44)):
    A = fn(); ok7 &= A["X_star"] == J["pairs"][key]["a_priori"]["X_star"]
    msg7.append("%s: X* = %.3e (A = %.3g, power %d)" % (A["pair"], A["X_star"], A["A"], A["K"]))
P = p19.Poly([3.0, 0.5, 2.0]); A, K = P.dominate(1 + math.log(29))
ok7 &= all(P.at(l) <= A * l ** K * (1 + 1e-12) for l in (4.37, 5, 10, 100, 1e5))
report("Q7 a priori bounds from Matveev's theorem and the closed-form D-step", ok7, "; ".join(msg7))

# Q8: the bound chains reproduced from scratch (4v2, 4v3 in full; 4v4 from its third round on) and compared with pairs19.json
t0 = time.time(); ok8 = True; msg8 = []
for key, fn, start in (("42", p19.round_42, 0), ("43", p19.round_43, 0), ("44", p19.round_44, 2)):
    R = J["pairs"][key]["rounds"]; X = int(R[start]["X"]); got = []
    while True:
        rec, Xn, boxes = fn(X); got.append(rec)
        if Xn >= X: break
        X = Xn
    keys = [k for k in R[0] if k not in ("seconds",)]
    same = len(got) == len(R) - start and all(str(a[k]) == str(b[k]) for a, b in zip(got, R[start:]) for k in keys) \
        and X == J["pairs"][key]["final_X"]
    ok8 &= same
    msg8.append("%s: %s (final %d)" % (key, " -> ".join(str(r["new_X"]) for r in got), X))
report("Q8 bound chains reproduced", ok8, "; ".join(msg8) + " (%.0f s)" % (time.time() - t0))

# Q9: the lattice step with four logarithms is a true lower bound (brute force over the whole box)
X0 = 9; Tp = T_of(((19, 1, 1),)); lg = Decimal(Tp).ln()
C = 10 ** (4 * len(str(X0)) + 8)
lam, _ = reduce_fast(X0, C, [LN[19], LN[2], LN[5], lg])
fl = [math.log(19), math.log(2), math.log(5), float(lg)]
best = min(abs(x[0] * fl[0] + x[1] * fl[1] + x[2] * fl[2] + x[3] * fl[3])
           for x in itertools.product(range(-X0, X0 + 1), repeat=4) if any(x))
report("Q9 the four-logarithm lattice bound is a true lower bound (prefix N l m N, |x| <= 9, whole box)",
       lam is not None and lam <= best, "lattice bound %.3e, smallest found %.3e" % (float(lam) if lam else -1, best))

# Q10: the exact check: pairs19.json, plus an independent full comparison of all words up to small gap sums, and psi from
# the letters as maps
E = J["exact_check"]
okE = not E["clashes"] and E["planted_found"] == E["planted_disguised_words"] > 0
allK = {}; clash = []; nw = 0
for nn, B in ((0, 0), (1, 30), (2, 16), (3, 9)):
    for w in (words(nn, B) if nn else [()]):
        k = K_of(w); nw += 1
        if k in allK: clash.append((name(w), name(allK[k])))
        else: allK[k] = w
S = {"l": Fr(19, 200), "m": Fr(1, 20), "N": Fr(1, 19), "K": Fr(1, 10)}
Dd = {"l": Fr(0), "m": Fr(0), "N": Fr(1, 19), "K": Fr(1, 10)}
def psi_maps(word):
    s, t = Fr(1), Fr(0)
    for c in word: s, t = s * S[c], s * Dd[c] + t
    return t / s
gaps = [""] + ["".join(x) for k in range(1, 3) for x in itertools.product("lm", repeat=k)]
byk = {}; nraw = 0
for g1, x2, g2, x3, g3, x4 in itertools.product(gaps, "NK", gaps, "NK", gaps, "NK"):
    if g1 == "" and x2 == "N": continue
    w = "N" + g1 + x2 + g2 + x3 + g3 + x4; ps = psi_maps(w); nraw += 1
    byk.setdefault(Fr(strip(ps.numerator), strip(ps.denominator)), set()).add(ps)
raw_bad = sum(1 for k in byk if len(byk[k]) > 1)
report("Q10 the exact comparison: no clash in pairs19.json, planted words found; all words with 1 to 4 letters up to small "
       "gap sums distinct; raw four-letter words from the maps: equal fingerprint means equal psi",
       okE and not clash and raw_bad == 0,
       "pairs19.json: %d short-box words against %d two-, %d three- and %d four-letter words, %d clashes, %d/%d planted found; "
       "independent: %d words, %d clashes; %d raw words (gaps <= 2 letters), %d fingerprints, %d with two psi" %
       (E["short_box_words"], E["two_letter_words"], E["three_letter_words"], E["four_letter_words"], len(E["clashes"]),
        E["planted_found"], E["planted_disguised_words"], nw, len(clash), nraw, len(byk), raw_bad))

n_ok = sum(r["passed"] for r in res)
lines.append("%d/%d checks pass" % (n_ok, len(res))); print(lines[-1])
open("checks19_output.txt", "w").write("\n".join(lines) + "\n")
json.dump(dict(passed=n_ok, total=len(res), checks=res), open("checks19.json", "w"), indent=1)
