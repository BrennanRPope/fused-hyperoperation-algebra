"""Checks for cusp-tail. Run: python3 checks.py. DRAFT."""
import numpy as np, json, sys, os
sys.path.insert(0, "../knauf-field"); from knauf import levels
from sympy import totient
ok = []
def check(name, cond): ok.append(bool(cond)); print(("PASS " if cond else "FAIL ") + name)
for n in (12, 16, 20):
    q, c = levels(n); d = dict(zip(q.astype(int), c.astype(int)))
    check("length %d: every denominator q <= n+1 appears exactly totient(q) times" % n, all(d.get(m, 0) == int(totient(m)) for m in range(1, n+2)))
    check("length %d: denominators n+2 and n+3 each miss their two edge fractions" % n, all(d.get(m, 0) == int(totient(m)) - 2 for m in (n+2, n+3)))
dm = json.load(open("demod.json")); c = json.load(open("combo.json")); w = json.load(open("../length-smoothing/wide.json")); nl = json.load(open("null.json")); cr = json.load(open("correct.json"))
check("removing the cusp's spin n^(rho-1) makes the error's direction far stiller for all 6 zeros", all(b > a + 0.3 for a, b in zip(dm["raw"], dm["demod"])))
check("for the first four zeros the demodulated direction is steady (>= 0.9)", all(v >= 0.9 for v in dm["demod"][:4]))
rms = lambda rows: np.sqrt(np.mean([r[0]**2 for r in rows]))
check("cusp fit over half the chain: miss below 0.0007 at every length 29..33", all(rms(c["rows"]["cusp half N=%d" % N]) < 0.0007 for N in range(29, 34)))
check("on lengths 29..32 it beats the triangle window by at least 2x (first 5 zeros)", all(rms(c["rows"]["cusp half N=%d" % N][:5]) < rms(w["rows"]["k=N/2 N=%d" % N])/2 for N in range(29, 33)))
check("pair of lengths alone (cusp pair) beats single chains at every length 16..33", all(rms(cr["rows"]["cusp pair N=%d" % N]) < rms(cr["rows"]["single N=%d" % N]) for N in range(16, 34)))
right = np.mean(nl["right n^(1-b)"])
check("null test: every wrong shape except the tiny phase offset does at least 1.3x worse", all(np.mean(v) > 1.3*right for k, v in nl.items() if not k.startswith("right") and "0.5i" not in k))
check("null test: dropping the phase altogether is at least 10x worse", np.mean(nl["no phase n^(-1/2)"]) > 10*right)
print("%d/%d passed" % (sum(ok), len(ok)))
