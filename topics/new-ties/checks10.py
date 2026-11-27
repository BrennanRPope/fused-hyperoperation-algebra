"""Round 10 checks (README section 14). Run after exact.py. Writes checks10_output.txt and checks10.json."""
import json, random
from fractions import Fraction as F
import exact as X

E = json.load(open("exact.json")); R = []
def check(name, ok, detail=""):
    R.append(dict(name=name, ok=bool(ok), detail=str(detail))); print(("PASS " if ok else "FAIL ") + name, detail)

e1 = E["E1"]; c = e1["conditions"]
check("X1 3-adic break on the surface: (231/3200 about 0) o (3/77 about 1) = (3/80 about 0) o (3/40 about 2), "
      "every multiplier has 3-adic valuation 1", e1["collision_equal"] and c["on_surface"] and set(e1["v3"].values()) == {1},
      e1["collision"])
check("X2 hypotheses of Theorem P12: zones disjoint (I), push (II), shared odd prime 3 (P)",
      c["I_class1_vs_N"] and c["I_N_vs_K"] and c["II"] and 3 in e1["shared_primes"] and c["in_unit"], e1["shared_primes"])
check("X3 certified: exact dimension in [0.4232693, 0.4232694], D* in [0.4531865, 0.4531866], drop > 0.0299",
      0.4232693 < e1["exact_dim"][0] and e1["exact_dim"][1] < 0.4232694 and 0.4531865 < e1["Dstar"][0]
      and e1["Dstar"][1] < 0.4531866 and e1["drop_at_least"] > 0.0299, (e1["exact_dim"], e1["Dstar"]))
check("X4 the reduced-piece weight equals Round 9's family weight f - P1 Q (mu kappa)^s", e1["formulas_agree"], "")
s2 = E["E2_sanity"]
check("X5 sanity (not used in the proof): all reduced words up to length 8 give different maps; counts match the family model",
      s2["all_distinct"] and s2["match"], (s2["reduced_words_upto8"], s2["distinct_counts"]))
# X6: the p-adic step: a word whose outer letter is N or K sends 1 to a 3-adic unit
lam, mu, nu, kap = X.surface(F(3, 80), F(3, 40))
H = {"l": X.hom(0, lam), "m": X.hom(0, mu), "N": X.hom(1, nu), "K": X.hom(2, kap)}
rng = random.Random(10); ok = True
for t in range(400):
    w = rng.choice("NK") + "".join(rng.choice("lmNK") for _ in range(rng.randint(0, 6)))
    f = (F(1), F(0))
    for a in w: f = X.comp(f, H[a])
    ok &= X.vp(f[0] * 1 + f[1], 3) == 0
check("X6 p-adic step: 400 random words with outer letter N or K send 1 to a 3-adic unit (exact)", ok, "")
# X7: the push, exactly, for a = 1..40
ok = all((F(77, 40)**a) * (1 - nu) > 1 + nu and (F(77, 40)**a) * (2 - 2 * kap) > 1 + nu for a in range(1, 41))
check("X7 push: rho = (77/40)^a moves the N-zone and the K-zone clear of the N-zone for a = 1..40 (exact; grows with a)", ok, "")
# X8: sanity, random pairs of different reduced words never give the same map
def random_reduced(rng, npieces):
    out = ""
    for _ in range(npieces):
        al, be = rng.randint(0, 2), rng.randint(0, 2)
        if al + be == 0: al = 1
        v = "".join(rng.choice("NK") for _ in range(rng.randint(1, 3)))
        if be and v[0] == "K": v = "N" + v[1:]
        out += "l" * al + "m" * be + v
    return out
def canon(w):  # class-1 letters commute inside a run
    out, run = [], []
    for a in w:
        if a in "lm": run.append(a)
        else: out += sorted(run); run = []; out.append(a)
    return "".join(out + sorted(run))
ok = True; n = 0
for t in range(3000):
    u, v = random_reduced(rng, rng.randint(1, 4)), random_reduced(rng, rng.randint(1, 4))
    if canon(u) == canon(v): continue
    fu, fv = (F(1), F(0)), (F(1), F(0))
    for a in u: fu = X.comp(fu, H[a])
    for a in v: fv = X.comp(fv, H[a])
    n += 1; ok &= fu != fv
check("X8 sanity: random pairs of different reduced words give different maps (exact)", ok and n > 2500, n)
seq = E["E4_sequence"]; vals = [s["exact_dim"][0] for s in seq]
up = E["E5_original"]["upper"][1]
check("X9 exact points p/(20p+1), p/(10p+1) (p = 3 ... 1009) increase toward the original's family bound",
      all(a < b for a, b in zip(vals, vals[1:])) and 0 < up - vals[-1] < 2e-5, [round(v, 7) for v in vals])
o = E["E5_original"]
check("X10 original break: no shared odd prime; first bad pair (a, b) = (13, 10); new certified lower bound within 2e-7 "
      "of the family bound and far above Round 9's 0.4367",
      o["shared_primes"] == [] and o["first_bad"][:2] == [13, 10] and o["kmax"] == 12 and o["new_lower"][1] < o["upper"][0]
      and o["upper"][1] - o["new_lower"][0] < 2e-7 and o["new_lower"][0] > 0.47, (o["new_lower"], o["upper"]))
okpts = [g for g in E["E3_region"] if g[5] is not None]
check("X11 on the whole grid where (I) and (II) hold, the exact size is strictly below D* (drop > 0.01)",
      len(okpts) > 3000 and all(g[6] - g[5] > 0.01 for g in okpts), len(okpts))
lo = open("lean/compile_output.txt").read()
check("X12 Lean: 53 theorems build with core Lean 4.37, no sorry",
      "53 theorems, exit 0, no sorry" in lo and "sorryAx" not in lo, lo.splitlines()[0])
print(f"{sum(r['ok'] for r in R)}/{len(R)} pass")
json.dump(R, open("checks10.json", "w"), indent=1)
open("checks10_output.txt", "w").write("\n".join(("PASS " if r["ok"] else "FAIL ") + r["name"] + " " + r["detail"] for r in R)
                                       + f"\n{sum(r['ok'] for r in R)}/{len(R)} pass\n")
