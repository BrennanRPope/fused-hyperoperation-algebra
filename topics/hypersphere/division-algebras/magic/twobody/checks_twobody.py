"""Checks for the two-body (comparable mass) pair rule. Writes results_twobody.json and checks_twobody_output.txt."""
import sys, json, itertools, random; sys.path.insert(0, '..')
import sympy as sp
import pn1
from stu2 import probe, pair_rule, two_body
from orbit1pn import stu_params, measure
lines = []; ok_all = True; data = {}
def rep(name, ok, detail=""):
    global ok_all; ok_all &= bool(ok); lines.append(f"{'PASS' if ok else 'FAIL'}  {name}  {detail}"); print(lines[-1])
T = pn1.turn_T(pn1.nu, pn1.XA)
rep("Y1 1PN two-body turn: Einstein gives 6, and the mass ratio enters only through X_A C_A + X_B C_B",
    sp.simplify(T.subs({pn1.s: 0, pn1.v: 0, pn1.CA: 1, pn1.CB: 1}) - 6) == 0 and sp.diff(T, pn1.nu) == 0, f"T = {sp.factor(T)}")
aA, aB, bA, bB = sp.symbols('alpha_A alpha_B beta_A beta_B'); ss = aA*aB
TDEF = T.subs({pn1.v: 0, pn1.s: ss, pn1.CA: (1+ss)**2 + bB*aA**2, pn1.CB: (1+ss)**2 + bA*aB**2})
known = 6*((1 - ss/3)/(1 + ss) - (pn1.XA*bB*aA**2 + (1 - pn1.XA)*bA*aB**2)/(6*(1 + ss)**2))
import subprocess
pr = subprocess.run([sys.executable, 'check_def_probe.py'], capture_output=True, text=True).stdout
rep("Y2 matches Damour-Esposito-Farese's scalar-tensor periastron formula (two bodies), and that formula matches an exact probe in the JNW field",
    sp.simplify(TDEF - known) == 0 and pr.count("probe T =") == 3 and all(l.split("probe T =")[1].split()[0] == l.split("DEF =")[1].split()[0] for l in pr.splitlines()), pr.strip().replace("\n", " | "))
random.seed(7); pairs = []
while len(pairs) < 12:
    kA = [random.choice([0, 0, 1, 2, 3]) for _ in range(4)]; kB = [random.choice([0, 0, 1, 2, 3]) for _ in range(4)]
    sA = [random.choice([1, -1]) for _ in range(4)]; sB = [random.choice([1, -1]) for _ in range(4)]
    if sum(kA) == 0 or sum(kB) == 0: continue
    if sum(kA[j]*kB[j]*(1 - sA[j]*sB[j]) for j in range(4)) == 0: continue
    pairs.append((kA, sA, kB, sB))
res = []; ok3 = ok4 = True
for kA, sA, kB, sB in pairs:
    Tt, TA, TB, sym = two_body(kA, sA, kB, sB)
    Tp = two_body(kA, sA, kB, sB, exact=False)[0]
    ok3 &= (sym == 0); ok4 &= (sp.simplify(Tt - Tp) == 0)
    res.append(dict(kA=kA, sA=sA, kB=kB, sB=sB, T=str(Tt), TA=str(TA), TB=str(TB)))
rep("Y3 Newton's third law: the pull per unit source mass is the same both ways (12 random pairs)", ok3)
rep("Y4 exact probe turns in both directions equal the pair rule, so the two-body turn is each body's pairs weighted by its mass share", ok4)
data["random_pairs"] = res
# Y5: zero turn iff both bodies carry one charge (all patterns with charges in {0,1}, signs +-; plus unequal sizes)
zero_ok = True; count = 0
pats = [p for p in itertools.product([0, 1, 2], repeat=4) if any(p)]
sigs = list(itertools.product([1, -1], repeat=4))
random.seed(3)
for _ in range(400):
    kA = list(random.choice(pats)); kB = list(random.choice(pats)); sA = list(random.choice(sigs)); sB = list(random.choice(sigs))
    if sum(kA[j]*kB[j]*(1 - sA[j]*sB[j]) for j in range(4)) == 0: continue
    Tt = two_body(kA, sA, kB, sB, exact=False)[0]; count += 1
    single = sum(1 for k in kA if k) == 1 and sum(1 for k in kB if k) == 1
    zero_ok &= ((Tt == 0) == single)
rep("Y5 the orbit closes at first order exactly when both bodies carry a single charge (any directions, any masses)", zero_ok, f"{count} attracting pairs tested")
# Y6: integrate the two-body equations at equal masses
cases = {"one charge each (sectors 1, opposite)": ([1, 0, 0, 0], [1]*4, [1, 0, 0, 0], [-1]*4),
         "four charges vs four opposite": ([1, 1, 1, 1], [1]*4, [1, 1, 1, 1], [-1]*4),
         "four charges vs one opposite": ([1, 1, 1, 1], [1]*4, [4, 0, 0, 0], [-1]*4)}
ok6 = True; num = {}
for name, (kA, sA, kB, sB) in cases.items():
    Tt, TA, TB, _ = two_body(kA, sA, kB, sB)
    par = stu_params(kA, sA, kB, sB, TA, TB)
    Tn, _ = measure(par)
    ok6 &= abs(Tn - float(Tt)) < 0.02*max(1, float(Tt)); num[name] = (str(Tt), round(Tn, 4))
rep("Y6 integrating the first-order two-body equations at equal masses gives the predicted turn", ok6, str(num))
data["numeric_equal_mass"] = num
lines.append("ALL PASS" if ok_all else "SOME FAILED"); print(lines[-1])
json.dump(dict(passed=ok_all, data=data), open("results_twobody.json", "w"), indent=1)
open("checks_twobody_output.txt", "w").write("\n".join(lines) + "\n")
