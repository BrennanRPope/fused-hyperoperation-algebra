"""Round 9 checks (README section 13). Run after breaks.py. Writes checks9_output.txt and checks9.json."""
import json, random, re
from fractions import Fraction as F
import breaks as Bk
import general as G
import plane as P

B = json.load(open("breaks.json")); R = []
def check(name, ok, detail=""):
    R.append(dict(name=name, ok=bool(ok), detail=str(detail))); print(("PASS " if ok else "FAIL ") + name, detail)

b1 = B["B1_collision"]; n = B["B1_numbers"]
check("Y1 exact collision: (19/200 about 0) o (1/19 about 1) = (1/20 about 0) o (1/10 about 2), inside Lemma A' box",
      b1["equal"] and b1["lemmaA"] and b1["x"] == ["1/200", "9/100"], b1["x"])
check("Y2 certified: D* in [0.50353, 0.50354], family bound below 0.47049, drop at least 0.033",
      0.50353 < n["Dstar"][0] and n["Dstar"][1] < 0.50354 and n["s_family_upper"][1] < 0.47049
      and n["Dstar"][0] - n["s_family_upper"][1] > 0.033, (n["Dstar"], n["s_family_upper"]))
lo_, hi_ = [float(re.findall(r"[\d.]+", x)[0]) for x in n["f_at_Dstar"]]
check("Y3 f = (P1 - 1)(Q - 1) crosses 1 inside the D* bracket (free-product identity)", lo_ >= 1 >= float(re.findall(r"[\d.]+", n["f_at_Dstar"][1])[0]) or (lo_ > 1 and hi_ < 1), n["f_at_Dstar"])
check("Y4 order of bounds: family bound < one-piece bound < D*",
      n["s_family_upper"][1] < n["s_Y_upper"][0] and n["s_Y_upper"][1] < n["Dstar"][0], "")
lb = max(v[1] for v in n["lower_bounds"].values())
check("Y5 lower bound (subsystem without mu) below the family bound", lb < n["s_family_upper"][0], lb)
# Y6 the family: (a lam).(nu w) = (a mu).(kappa w) exactly, random a (class-1 run) and w (word in classes 2, 3)
rng = random.Random(6); ok = True
H = {"l": Bk.hom(F(0), Bk.LAM), "m": Bk.hom(F(0), Bk.MU), "n": Bk.hom(F(1), Bk.NU), "k": Bk.hom(F(2), Bk.KAP)}
def word(ws):
    f = (F(1), F(0))
    for a in ws: f = Bk.comp(f, H[a])
    return f
for t in range(200):
    a = "".join(rng.choice("lm") for _ in range(rng.randint(0, 3)))
    w = "".join(rng.choice("nk") for _ in range(rng.randint(0, 3)))
    ok &= word(a + "l" + "n" + w) == word(a + "m" + "k" + w)
check("Y6 the break spreads to the whole family a.lam.nu.w = a.mu.kappa.w (200 random, exact fractions)", ok, "")
# Y7 one map per class inside the box never breaks: random rational systems, counts = trace counts
rng = random.Random(7); ok = True
for t in range(10):
    cl = [((F(i), F(0)), [(F(rng.randint(1, 10), 100), F(0))]) for i in range(3)]
    ok &= G.distinct_counts(cl, 8) == P.trace_counts([1, 1, 1], 8)
check("Y7 one map per class inside the box: no collisions to length 8 (10 random systems)", ok, "")
r5 = B["B5_random_real"]
check("Y8 random real small parameters: different words never gave equal maps", r5["equal_maps"] == 0 and r5["pairs"] > 2000, r5)
c = B["B2_counts"]
check("Y9 the one-break-plus-family prediction matches the distinct-map counts at every length 1-9",
      c["model_matches"] and c["distinct"][-1] == 76096 and c["trace"][-1] == 151316, c["distinct"])
g = B["B4_surface"]
ok = all(r[2] - r[3] > 0.017 for r in g)
for mu, kap, D, Y in g:
    lam = mu * (2 - kap); nu = kap * mu / lam
    ok &= abs(lam * nu - mu * kap) < 1e-15 and abs(lam * (1 - nu) - 2 * mu * (1 - kap)) < 1e-15 and max(lam, mu, nu, kap) < 1/9
check("Y10 every grid point on the collision surface is a break in the box, with drop above 0.017", ok, f"{len(g)} points")
lean = open("lean/compile_output.txt").read(); src = open("lean/NewTies.lean").read()
check("Y11 Lean builds 41 theorems, no sorry", "exit status 0" in lean and "sorryAx" not in lean
      and src.count("\ntheorem ") == 41 and "example_collision" in src, "")
txt = open("README.md").read()
check("Y12 README section 13 states P10, P11 and the example's numbers", all(s in txt for s in ["## 13.", "Theorem P10", "Corollary P11", "0.47048", "0.50353"]), "")
npass = sum(r["ok"] for r in R)
open("checks9_output.txt", "w").write("\n".join(("PASS " if r["ok"] else "FAIL ") + r["name"] + "  " + r["detail"] for r in R)
                                      + f"\n{npass}/{len(R)} pass\n")
json.dump(R, open("checks9.json", "w"), indent=1)
print(f"{npass}/{len(R)} pass")
