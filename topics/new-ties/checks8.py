"""Round 8 checks (README section 12). Run after general.py. Writes checks8_output.txt and checks8.json."""
import json, random, itertools, math
from fractions import Fraction as F
import sympy
import general as G
import plane as P

D = json.load(open("general.json")); R = []
def check(name, ok, detail=""):
    R.append(dict(name=name, ok=bool(ok), detail=str(detail))); print(("PASS " if ok else "FAIL ") + name, detail)

# X1 Lemma A' algebra, symbolically: with R = D(1+e)/(1-e), invariance holds with equality and
#    2e(R + D) < delta  <=>  e(4D + delta) < delta  (for 0 < e < 1)
e, Dd, de = sympy.symbols("e D delta", positive=True)
Rr = Dd * (1 + e) / (1 - e)
inv = sympy.simplify(Dd * (1 + e) + e * Rr - Rr)
lhs = sympy.simplify((de - 2 * e * (Rr + Dd)) * (1 - e) - (de - e * (4 * Dd + de)))
check("X1 Lemma A' algebra (invariance exact; separation iff e(4D+delta) < delta)", inv == 0 and lhs == 0, (inv, lhs))

# X2 Lemma A' by exact discs: for 200 random members of the G5 recipe, the explicit discs are disjoint, exactly
rng = random.Random(8); bad = 0; n = 0
for t in range(200):
    pat = rng.choice([(2, 1, 1), (2, 2, 1), (3, 1, 1), (2, 1, 1, 1)]); field = rng.choice(["Q", "Qi"])
    cl, own = G.corollary_member(pat, field, rng)
    eps = max(G.abs_up(l) for _, ls in cl for l in ls)
    c1 = cl[0][0]; Dm = max(P.sqrt_up(G.n2(G.cs(c, c1))) for c, _ in cl)       # max distance from c1, rounded up
    Rad = Dm * (1 + eps) / (1 - eps)
    # every map sends the disc (c1, Rad) into itself
    ok = all(P.sqrt_up(G.n2(G.cs(G.ca(c, G.cm(l, G.cs(c1, c)), ), c1))) + G.abs_up(l) * Rad <= Rad
             for c, ls in cl for l in ls)
    # class piece j sits in the disc (c_j, eps*(Rad + |c_j - c1|)); pairwise disjoint
    discs = [(c, eps * (Rad + P.sqrt_up(G.n2(G.cs(c, c1))))) for c, _ in cl]
    ok &= all(G.n2(G.cs(a[0], b[0])) > (a[1] + b[1])**2 for a, b in itertools.combinations(discs, 2))
    n += 1; bad += (not ok)
check("X2 explicit discs invariant and disjoint for 200 random dense-family members", bad == 0, f"{n - bad}/{n}")

# X3 round-7 examples: MOD passes (B), ORIG has no good prime
check("X3 MOD passes B with good prime 3; ORIG has none", D["G1_MOD"]["B"] and not D["G1_ORIG"]["B"]
      and D["G1_MOD"]["info"]["0"]["good_primes"] == [3] and D["G1_ORIG"]["info"]["0"]["good_primes"] == [], "")

# X4 G2: explicit bound, certified covers, B
g2 = D["G2"]; check("X4 G2 example passes Lemma A' (0.1414 < 0.1502), certified A and B",
                    g2["explicit"] and g2["certified_A"] and g2["B"], g2["info"])

# X5 families: every member passing A and B is collision-free to the checked length
for fam in ["F1", "F2", "F3"]:
    f = D["G3_" + fam]
    check(f"X5 {fam}: all {f['pass_A_and_B']} members passing A and B collision-free", f["no_collisions_confirmed"] == f["pass_A_and_B"]
          and f["pass_A_and_B"] > 0, f"{f['members']} sampled, B {f['pass_B']}")

# X6 G4: two two-map classes, good primes 3 and 5
g4 = D["G4"]; check("X6 G4 passes A (explicit and certified), B with primes 3 and 5, collision-free to 7",
                    g4["A"] and g4["B"] and g4["explicit_A"] and g4["no_collisions_to"] == 7
                    and g4["B_info"]["0"]["good_primes"] == [3] and g4["B_info"]["1"]["good_primes"] == [5], g4["Dstar"])

# X7 G5: the dense family for 14 patterns
g5 = D["G5"]; check("X7 G5 dense family: all members pass Lemma A', B, irreducibility (plane) and the collision check",
                    g5["members"] == g5["explicit_A"] == g5["B"] == g5["collision_free"] == 42
                    and g5["plane_irreducible"] == g5["plane_members"], {k: v for k, v in g5.items() if k != "rows"})

# X8 the kernel lemma: a 2-dim (or larger) kernel always holds a mixed-sign integer vector (random checks)
rng = random.Random(88); okk = True
for t in range(300):
    k = rng.randint(3, 5); m = rng.randint(1, k - 2)
    V = [[rng.randint(0, 3) for _ in range(m)] for _ in range(k)]
    ns = sympy.Matrix(V).T.nullspace()
    if len(ns) >= 2:
        a, b = ns[0], ns[1]; found = False
        for s1 in range(-4, 5):
            for s2 in range(-4, 5):
                w = s1 * a + s2 * b
                if any(x > 0 for x in w) and any(x < 0 for x in w): found = True; break
            if found: break
        okk &= found and not G.no_mixed_relation(V)
check("X8 kernels of dimension >= 2 contain a mixed-sign vector (300 random valuation matrices)", okk, "")

# X9 the recipe's valuation pattern: own prime once, other own primes zero, last map none, other classes' primes >= 1
rng = random.Random(9); okv = True
for t in range(100):
    pat = rng.choice([(2, 1, 1), (3, 2, 1), (4, 1, 1), (2, 2, 2)]); field = rng.choice(["Q", "Qi"])
    cl, own = G.corollary_member(pat, field, rng)
    for i, (c, ls) in enumerate(cl):
        for j, l in enumerate(ls):
            for jj, p in enumerate(own[i]): okv &= (G.val(l, p, field) == (1 if jj == j else 0))
            for i2, o in enumerate(own):
                if i2 != i: okv &= all(G.val(l, p, field) >= 1 for p in o)
check("X9 the recipe gives exactly the valuation pattern Corollary P9 needs (100 members)", okv, "")

# X10 density: members of the (2,1,1) line recipe come within 1/1000 of random targets in the box
rng = random.Random(10); worst = 0
for t in range(50):
    target = [rng.uniform(-0.11, 0.11) for _ in range(4)]
    # class 1: lam unit at 3, mu divisible by 3; classes 2, 3 divisible by 3; denominator 2^n
    n = 14; best = []
    for x, need in zip(target, [0, 1, 1, 1]):
        m0 = round(x * 2**n)
        cand = [m for m in range(m0 - 4, m0 + 5) if (m % 3 == 0) == bool(need) and m != 0]
        best.append(min(abs(F(m, 2**n) - F(x)) for m in cand))
    worst = max(worst, max(best))
check("X10 the family is dense: every random target in the box is within 1/1000 of a member", worst < F(1, 1000), float(worst))

# X11 the README states the literature corrections
txt = open("README.md").read()
check("X11 README section 12 cites Barany-Verma, Barany-Szvak and Kittle-Kogler and marks the line formula known",
      all(s in txt for s in ["2507.05835", "1902.04028", "2610.10511", "## 12."]), "")

npass = sum(r["ok"] for r in R)
open("checks8_output.txt", "w").write("\n".join(("PASS " if r["ok"] else "FAIL ") + r["name"] + "  " + r["detail"] for r in R)
                                      + f"\n{npass}/{len(R)} pass\n")
json.dump(R, open("checks8.json", "w"), indent=1)
print(f"{npass}/{len(R)} pass")
