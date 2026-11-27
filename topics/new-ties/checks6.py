"""Round 6 (DRAFT) checks: no accidental collisions for the original examples (README section 10). Expect 13/13."""
import json, math
from fractions import Fraction as F
import nocollide as N

results = []
def check(name, ok, detail=""):
    results.append((name, bool(ok), detail)); print(("PASS " if ok else "FAIL ") + name + (": " + detail if detail else ""))

all_words = N.all_words

# V1-V3 real ping-pong and window
for Q in (2, 3, 12):
    r = N.real_part(Q); check(f"V1 real ping-pong at walls ±{Q} (J1, J2 disjoint, 0 in neither)", r["disjoint_and_miss_0"], f"J1={r['J1']}, J2={r['J2']}")
bad, tail = N.window_misses(N.real_part(3)["window"])
check("V2 walls ±3: (3/4)^n, 5(3/4)^n, (3/4)^n/5 all miss the real window [12/13, 13/12]", not bad and tail)
# V3 2-adic criterion for even Q, fails for odd Q
ev = all(N.two_adic_constants(Q)["criterion"] for Q in range(2, 41, 2)); od = not any(N.two_adic_constants(Q)["criterion"] for Q in range(3, 41, 2))
check("V3 2-adic criterion (v2(kappa) >= 1, |1-kappa|_2 = 1) holds for every even Q up to 40 and no odd Q", ev and od)
# V4 exhaustive: one 2-adic level at even Q
lv = {}
for Q in (2, 4, 6, 12):
    vals = set()
    for n in (1, 2, 3, 4):
        for w in all_words(n):
            vals.add(N.v(N.z1(N.word_apply(w, F(-Q), Q), Q), 2))
    lv[Q] = sorted(vals)
check("V4 exhaustive (runs up to 4, exponents up to 2/4): every word ending in law 2 lands on ONE 2-adic level, Q = 2, 4, 6, 12",
      all(len(x) == 1 for x in lv.values()), str(lv))
vals3 = set()
for n in (1, 2, 3, 4):
    for w in all_words(n): vals3.add(N.v(N.z1(N.word_apply(w, F(-3), 3), 3), 2))
check("V5 at walls ±3 the same words spread over several 2-adic levels (why the 2-adic proof can't work there)", len(vals3) > 2, str(sorted(vals3)))
# V6 5-adic constants
c5 = N.five_adic_constants(3)
check("V6 walls ±3: walls distinct mod 5; an odd c-run sends points near 3 to infinity mod 5 (not 1, not 3)", c5["distinct"] and c5["odd_run_avoids_1_and_Q"] and c5["odd_run_lands_on"] == "infinity")
# V7 exhaustive 5-adic levels in the allowed class
v5 = set()
for n in (1, 2, 3, 4):
    for w in all_words(n, max_e=2, max_k=6, rule=N.p6_rule): v5.add(N.v(N.z1(N.word_apply(w, F(-3), 3), 3), 5))
check("V7 exhaustive: words in the allowed class at walls ±3 land at 5-adic size 1 or 5 only", v5 <= {0, -1}, str(sorted(v5)))
# V8 the excluded pattern really leaves the two levels
lvls = {k: N.v(N.z1(N.word_apply([('A', 0, 1), ('C', k)], F(-3), 3), 3), 5) for k in (2, 4, 6, 8, 12, 14, 16, 18)}
check("V8 an even c-run on a single letter can land deeper (b then c^8: size 25; b then c^18: size 125), so this proof needs the class restriction",
      lvls[8] == -2 and lvls[18] == -3, str(lvls))
# V9 certified lower bound for the free sub-system at ±3
lb = N.sub_lower(3)
check("V9 certified: dimension at walls ±3 is at least 0.4259 (free sub-system, Collatz-Wielandt)", lb["certified"] and F(lb["t_low"]) >= F(4259, 10000), str(lb))
# V10 consistency with the length-14 exact count
cm = json.load(open("collisions_modp.json"))
check("V10 the length-14 exact counts (round 5) match the shared-centre count for all three examples", all(cm[k]["match"] for k in cm))
# V11 Lean
out = open("lean/compile_output.txt").read()
n = out.count("depends on axioms") + out.count("does not depend on any axioms")
check("V11 Lean NewTies.lean builds: 22 theorems, no sorry, standard axioms only", "exit status 0" in out and n == 22 and "sorryAx" not in out, f"{n} theorems")
print(f"{sum(r[1] for r in results)}/{len(results)} pass")
json.dump([dict(name=a, ok=b, detail=c) for a, b, c in results], open("checks6.json", "w"), indent=1)
