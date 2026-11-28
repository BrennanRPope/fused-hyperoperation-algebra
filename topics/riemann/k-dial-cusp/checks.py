"""Checks for k-dial-cusp. Run: python3 checks.py. DRAFT."""
import numpy as np, json, mpmath as mp
from identity import tree
from rewind import Fk
ok = []
def check(name, cond): ok.append(bool(cond)); print(("PASS " if cond else "FAIL ") + name)
for k in (0.05, 0.3):
    T = tree(k, 12); q = 1 + k
    check("k = %.2f: every tree sum u o_k v equals the q-integer [p+q]_q, q = 1+k (depth 12)" % k, max(abs(w - (q**m - 1)/k)/((q**m - 1)/k) for w, m, a, b in T) < 1e-12)
T0 = tree(0.0, 12); check("k = 0: tree sums are plain integers p+q, each coprime pair once", all(w == m for w, m, a, b in T0) and len({(a, b) for _, _, a, b in T0}) == len(T0))
b = 1.5+14.134725141734695j; f1, _ = Fk(b, 0.01); f2, _ = Fk(b, 0.01)
L = np.log1p(0.01); m = np.arange(1, 6001); ph = __import__("rewind").PH[1:6001]
direct = complex((ph*np.exp(-b*np.log((np.expm1(L*m))/0.01))).sum())
check("k-dial chain sum = sum phi(m) [m]_q^(-b) matches the Bose-weight form (k/L)^b * sum phi m^(-b) h(Lm) (k = 0.01)", abs(direct - (0.01/L)**b*f1) < 1e-9*abs(direct))
r = json.load(open("rewind.json")); ks = np.array(r["ks"])
for key in ("5.0000", "14.1347", "25.0109"):
    e = np.array([v[1] for v in r["rows"][key]]); sl = np.polyfit(np.log(ks[-5:]), np.log(e[-5:]), 1)[0]
    check("height %s: rewind miss shrinks roughly in step with k (log-log slope %.2f within 0.6..1.3)" % (key, sl), 0.6 < sl < 1.3)
ex = json.load(open("extrap.json")); cs = json.load(open("../cusp-size/compare.json"))
check("at k = 0.0002 the rewind is less accurate than the cusp fit (median over the line)", np.median(ex["raw"]) > np.median(cs["cusp2"]))
check("adding the L and L log L terms does not help (median miss does not drop)", np.median(ex["LlogL"]) >= np.median(ex["raw"]))
print("%d/%d passed" % (sum(ok), len(ok)))
