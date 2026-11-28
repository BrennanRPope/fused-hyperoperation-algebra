"""Checks for cusp-size. Run: python3 checks.py (needs /tmp/knauf_levels to length 33 only for rerunning scripts). DRAFT."""
import numpy as np, json
ok = []
def check(name, cond): ok.append(bool(cond)); print(("PASS " if cond else "FAIL ") + name)
s = json.load(open("scan.json")); c = json.load(open("compare.json")); sh = json.load(open("shape.json"))
F = np.array([complex(*v) for v in s["F"]]); T = np.array([complex(*v) for v in s["true"]]); t = np.array(s["t"])
check("cusp fit recovers zeta(b-1)/zeta(b) on Re b = 3/2 to under 0.2% (median relative error, t = 5..45)", np.median(abs(F-T)/abs(T)) < 0.002)
check("same above t = 13: median absolute error under 0.001", np.median(abs(F-T)[t > 13]) < 0.001)
check("fit beats the single length-33 chain at least 15x (median over the line)", np.median(c["cusp2"]) < np.median(c["single"])/15)
check("even the one-term cusp fit beats the single chain at least 10x", np.median(c["cusp1"]) < np.median(c["single"])/10)
A = np.array([complex(*v) for v in s["A"]])
check("cusp size A is not near zero at the first five zeta zeros", all(abs(A[np.argmin(abs(t-g))]) > 0.1 for g in (14.13, 21.02, 25.01, 30.42, 32.94)))
check("no tried zeta formula matches A (relative spread above 1 for every candidate)", all(v[2] > 1 for k, v in sh.items() if k != "log_terms"))
print("%d/%d passed" % (sum(ok), len(ok)))
