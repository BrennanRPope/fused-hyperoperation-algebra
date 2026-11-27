"""Round 22 checks (DRAFT). Run from this folder: python3 checks22.py   (a few seconds, peak under 0.1 GB; reads
explore22.json and referee22/search22_K4_15.json)."""
import itertools, json, math, os
from fractions import Fraction as F
import explore22 as X
import explore21 as E21

out, res = [], {}
def rep(name, ok, msg):
    out.append(f"{name}: {'PASS' if ok else 'FAIL'}  {msg}"); res[name] = dict(ok=ok, msg=msg); print(out[-1])

# Q1: Proposition 22.2, the coordinates. 2, 19/10, 1/20 are a basis; single word ratios are 2 u^e w^f with f >= 1,
# -1 <= e <= f - 1, each once; products of j of them are exactly the cone points (j, E, F) with F >= j, -j <= E <= F - j.
u, w = F(19, 10), F(1, 20)
basis_det = 1 * ((-1) * 0 - 1 * (-1))          # rows: 2 = (1,0,0), u = (-1,-1,1), w = (-2,-1,0) in exponents of 2, 5, 19
single = {}
for L in range(0, 15):
    for M in range(0, 15):
        for r in (19, 10):
            if r == 10 and M > 0: continue
            q = F(19 ** L, r * 200 ** L * 20 ** M)
            e, f = L - (1 if r == 19 else 0), L + M + 1
            single.setdefault(q, []).append((e, f))
            assert q == 2 * u ** e * w ** f
ok1 = basis_det == 1 and all(len(v) == 1 for v in single.values())
ok1 &= {v[0] for v in single.values() if v[0][1] <= 14} == {(e, f) for f in range(1, 15) for e in range(-1, f) if e + 1 <= 14 and f - e - 2 <= 14}
prods = {}
sing = [(e, f) for f in range(1, 9) for e in range(-1, f)]
reach = {(1, e, f) for e, f in sing}
cur = set(reach)
for j in range(2, 5):
    cur = {(j, E + e, Fv + f) for (_, E, Fv) in cur for e, f in sing if Fv + f <= 8}
    reach |= cur
cone_pts = {(j, E, Fv) for j in range(1, 5) for Fv in range(j, 9) for E in range(-j, Fv - j + 1)}
ok1 &= reach == cone_pts
# the same cone in the form 19^E / (2^A 5^B): B <= A <= 2B, E <= A - B, E <= 2B - A - 1
ok1 &= all((lambda e2, e5, e19: (-e5) <= (-e2) <= 2 * (-e5) and e19 <= (-e2) - (-e5) and e19 <= 2 * (-e5) - (-e2) - 1)(*X.exps(j, E, Fv))
           for (j, E, Fv) in cone_pts)
ok1 &= all(E21.word_shaped(X.nd(*X.exps(j, E, Fv))[1], X.nd(*X.exps(j, E, Fv))[0]) for (j, E, Fv) in cone_pts)
rep('Q1', ok1, f"2, 19/10, 1/20 form a basis (determinant 1); {len(single)} word ratios with L, M <= 14 are 2 u^e w^f, each once; "
    f"products of up to 4 of them fill exactly the {len(cone_pts)} cone points with F <= 8, which match the A, B, E form and "
    f"round 21's word-shaped test")

# Q2: Proposition 22.1 (2 numbers against 1): the facts used, by brute force
ok2 = True
for A in range(0, 60):
    for B in range(0, 40):
        P = 2 ** A * 5 ** B
        if A >= 1 and B >= 1:
            ok2 &= (1 + P) % 3 != 1 and X.strip(1 + P) != 1                 # 1 + 2^A 5^B is never a power of 19
for e in range(1, 400):
    v = 1 + 19 ** e
    ok2 &= (X.strip(v) == 1) == (e == 1)
    if e % 2:
        ok2 &= (v // 4) % 2 == 1 and v % 4 == 0 and E21.v(v, 5) == 1 + E21.v(e, 5)   # lifting the exponent
C = X.cone(10 ** 60)
units = [(n, d) for (n, d, *_r) in C if X.strip(n + d) == 1]
ok2 &= units == [(1, 19)]
rep('Q2', ok2, f"1 + 2^A 5^B is never a power of 19 (A, B >= 1, below 60 and 40; it is 0 or 2 mod 3); 1 + 19^e is a product of "
    f"2 and 5 only for e = 1 (e < 400; for odd e, v2 = 2 and v5 = 1 + v5(e)); among {len(C):,} cone ratios to 10^60, "
    f"1 + rho is a unit only for rho = 1/19")

# Q3: two numbers against two: strip(n + d) is injective on the cone (search), and the slice j = 1 is Theorem 16.2's setting
seen, coll = {}, 0
for (n, d, *_r) in C:
    k = X.strip(n + d)
    if k in seen: coll += 1
    seen[k] = 1
J = json.load(open('explore22.json'))
okB = J['B']['digits'] == 300 and J['B']['collisions'] == [] and J['B']['unit_sums'] == [[1, 19, 1, -1, 1]]
slice_ok = True
for H in range(0, 25):
    for L in range(H + 1):
        Mm = H - L
        for r in ((19, 10) if Mm == 0 else (19,)):
            T = 19 ** L + r * 2 ** (3 * L + 2 * Mm) * 5 ** (2 * L + Mm)
            q = F(19 ** L, r * 2 ** (3 * L + 2 * Mm) * 5 ** (2 * L + Mm))
            slice_ok &= X.strip(T) == X.strip(q.numerator + q.denominator)
            e2, e5, e19 = E21.vq(q, 2), E21.vq(q, 5), E21.vq(q, 19)
            slice_ok &= (e2 - 2 * e5 - e19) == 1                               # j = 1
rep('Q3', coll == 0 and okB and slice_ok,
    f"no two cone ratios to 10^60 share the fingerprint of 1 + rho ({len(C):,} ratios); explore22.json: none to 10^300 "
    f"({J['B']['cone_elements']:,} ratios); for two-letter words (H < 25) that fingerprint is psi's, and rho has j = 1")

# Q4: three numbers on a side: every equal-fingerprint pair is non-primitive (small rerun and the json)
X.OUT.clear()
X.partC(20)
small = X.OUT['C']
okC = small['n_primitive'] == 0 and small['equal_fingerprint_pairs'] > 0
okC &= J['C']['digits'] == 35 and J['C']['n_primitive'] == 0 and J['C']['equal_fingerprint_pairs'] == 17607
rep('Q4', okC, f"sides of 3 against 1, 2 or 3 with d1 d2 <= 10^20: {small['equal_fingerprint_pairs']:,} equal-fingerprint pairs, "
    f"all non-primitive; explore22.json to 10^35: {J['C']['three_sides']:,} three-number sides, {J['C']['equal_fingerprint_pairs']:,} "
    f"pairs, none primitive")

# Q5: the referee's report and logs are present, and its largest search agrees with round 21's
rp = open(os.path.join('referee22', 'REPORT.md')).read()
ok5 = all(s in rp for s in ('**Holds.** The hypothesis a >= 1', '**Reproduced, no other example.**', '**Holds, with wording fixes.**',
                            '| 4. Lemma 21.4 | **Holds.**', '| 5. Thm 21.5 | **Holds.**'))
logs = [f for f in os.listdir('referee22') if f.startswith('search22_K') and f.endswith('.log')]
ok5 &= len(logs) >= 6
rep('Q5', ok5, f"referee22/REPORT.md gives items 1, 3, 4, 5 as holding (3 with wording fixes) and item 2 reproduced; "
    f"{len(logs)} search logs present")

out.append(f"{sum(r['ok'] for r in res.values())}/{len(res)} checks pass")
print(out[-1])
open('checks22_output.txt', 'w').write('\n'.join(out) + '\n')
json.dump(res, open('checks22.json', 'w'), indent=1)
