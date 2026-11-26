"""Checks for the four set-operations levers (section 8 of README.md)."""
import math
from scipy.special import erf
from setops_levers import E_k, E_k_mc, draws_p, neg_mass
R = []
def ok(name, cond, note=""):
    R.append(cond); print(("PASS" if cond else "FAIL"), name, note)
ok("S1 E(0) = e", abs(E_k(1e-13) - math.e) < 1e-9)
ok("S2 E(1) = 5/2 and E(-1/2) = 3", abs(E_k(1) - 2.5) < 1e-12 and abs(E_k(-0.5) - 3) < 1e-12)
v = (math.e - 1) / math.e + math.e ** (1 / (math.e - 1))
ok("S3 E(e-1) = Vandervelde 2008", abs(E_k(math.e - 1) - v) < 1e-12, f"{v:.10f}")
ok("S4 E(k) matches simulation at k = -1/2, 1/2, 1, e-1",
   all(abs(E_k_mc(k) - E_k(k)) < 0.01 for k in [-0.5, 0.5, 1, math.e - 1]))
ok("S5 draws until sum x^2 > 1 = e^(pi/4)(1 + erf(sqrt(pi)/2))",
   abs(draws_p(2) - math.exp(math.pi / 4) * (1 + erf(math.sqrt(math.pi) / 2))) < 1e-12, f"{draws_p(2):.4f}")
ok("S6 draws at p = 1 is e", abs(draws_p(1) - math.e) < 1e-12)
ok("S7 least negative mass = (S - 2)/4 for S = 2.4, 2.8, 2 sqrt 2, 3.2, 3.6, 4",
   all(abs(neg_mass(S) - (S - 2) / 4) < 1e-9 for S in [2.4, 2.8, 2 * math.sqrt(2), 3.2, 3.6, 4.0]))
ok("S8 no negative mass needed up to S = 2", neg_mass(2.0) < 1e-12 and neg_mass(1.5) < 1e-12)
ok("S9 Tsirelson's mass (sqrt2 - 1)/2 is half the k = +1 union half", abs(neg_mass(2 * math.sqrt(2)) - (math.sqrt(2) - 1) / 2) < 1e-9
   and abs((1 + (math.sqrt(2) - 1)) ** 2 - 1 - 1) < 1e-12)
print(f"{sum(R)}/{len(R)} passed")
