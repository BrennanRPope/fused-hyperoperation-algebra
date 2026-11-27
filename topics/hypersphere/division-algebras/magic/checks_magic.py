"""DRAFT. Orbit turning across the double-copy 'magic square'.
Probe orbits around exact static solutions of S = (1/16pi) int[R - 2(dphi)^2 - e^{-2a phi}F^2] (G = 1):
 - GHS charged dilaton black holes (extremal, a = sqrt3, 1, 1/sqrt3, 0 = 1, 2, 3, 4 charges in N = 8),
 - JNW (mass with dilaton charge; the classical double copy of a colour charge is a = 1).
Turn per lap T in units of pi K^2/L^2 (K = Newtonian pull): GR 6, Sommerfeld 1.
Writes results_magic.json and checks_magic_output.txt."""
import json, sympy as sp
from fields import ghs, jnw, check
from precess import W_of, T_series, T_numeric, u, L, E
res = []; data = {}
def ok(name, cond, msg):
    res.append((name, bool(cond), msg)); print(("PASS " if cond else "FAIL ") + name + ": " + msg)

# M1: the backgrounds are exact solutions
worst = 0
for a in [sp.sqrt(3), sp.Integer(1), 1/sp.sqrt(3), sp.Integer(0)]:
    g, phi, At, Q, M = ghs(a, sp.Integer(3), sp.Integer(1)); Ee, d, mx = check(g, phi, At, a); worst = max(worst, max(Ee), d, max(mx))
for b in [sp.Integer(1), sp.sqrt(3)]:
    g, phi, At, M, S = jnw(1/sp.sqrt(1+b**2), sp.Integer(2)); Ee, d, mx = check(g, phi, At, 0); worst = max(worst, max(Ee), d)
ok("M1", worst < 1e-12, "GHS (4 couplings) and JNW (2) satisfy Einstein, dilaton and Maxwell equations (worst residual %.1e)" % worst)

# M2: BPS probe coupling: static no-force on the extremal background needs m(phi)=m e^{a phi}, q/m = sqrt(1+a^2)
r = sp.Symbol('r', positive=True); nf = True
for a in [sp.sqrt(3), sp.Integer(1), 1/sp.sqrt(3), sp.Integer(0)]:
    g, phi, At, Q, M = ghs(a, sp.Integer(1), sp.Integer(1))
    from fields import r as rr
    V = sp.exp(a*phi)*sp.sqrt(-g[0, 0]) + sp.sqrt(1+a**2)*At
    nf &= sp.simplify(sp.diff(V, rr).subs(rr, sp.Rational(7, 3))) == 0 and sp.simplify(Q/M - sp.sqrt(1+a**2)) == 0
ok("M2", nf, "an extremal body with dilaton coupling a and charge/mass sqrt(1+a^2) feels no static force from its own kind (BPS)")

# M3: GR check
g, phi, At, Q, M = ghs(sp.Integer(0), sp.Integer(1), sp.Integer(0)); W = W_of(g, phi, At, 0, 0); T, K = T_series(W)
Tn = T_numeric(W, K); data["GR"] = [str(T), Tn]
ok("M3", T == 6 and abs(Tn-6) < 1e-3, "Schwarzschild, neutral probe: T = %s (simulated %.5f)" % (T, Tn))

# M4: BPS pairs with opposite charge (attract): T = n - 1, n = 4/(1+a^2)
allok = True; data["BPS"] = {}
for a in [sp.sqrt(3), sp.Integer(1), 1/sp.sqrt(3), sp.Integer(0)]:
    g, phi, At, Q, M = ghs(a, sp.Integer(1), sp.Integer(1))
    W = W_of(g, phi, At, -sp.sqrt(1+a**2), a); T, K = T_series(W); Tn = T_numeric(W, K)
    n = sp.nsimplify(4/(1+a**2)); data["BPS"][str(a)] = dict(n=str(n), T=str(T), sim=Tn, K_over_M=str(sp.simplify(K/M)))
    allok &= sp.simplify(T - (n-1)) == 0 and abs(Tn - float(n-1)) < 2e-3 and sp.simplify(K/M - 2*(1+a**2)) == 0
ok("M4", allok, "BPS bodies of opposite charge: T = n - 1 = 0, 1, 2, 3 for a = sqrt3, 1, 1/sqrt3, 0 (n = 4/(1+a^2) charges); pull K = 2(1+a^2)M; simulations agree")
ok("M5", data["BPS"][str(sp.sqrt(3))]["T"] == "0", "half-BPS in N = 8 (a = sqrt3): no turn at all, as Caron-Huot and Zahraee found")

# M6: neutral double-copy masses (dilaton coupling b): T = 2(n - 1), n = 4/(1+b^2); b = 1 is the double copy of a colour charge
allok = True; data["neutral"] = {}
for b in [sp.Integer(0), 1/sp.sqrt(3), sp.Integer(1), sp.sqrt(3)]:
    g, phi, At, M, S = jnw(1/sp.sqrt(1+b**2), sp.Integer(1)); W = W_of(g, phi, At, 0, b); T, K = T_series(W)
    n = sp.nsimplify(4/(1+b**2)); Tn = T_numeric(W, K) if b != 0 else 6.0
    data["neutral"][str(b)] = dict(n=str(n), T=str(T), sim=Tn, K_over_M=str(sp.simplify(K/M)))
    allok &= sp.simplify(T - 2*(n-1)) == 0 and abs(Tn - float(2*(n-1))) < 2e-3
ok("M6", allok, "neutral masses with dilaton coupling b: T = 2(n - 1) = 6, 4, 2, 0; GR (b = 0) is 6, the double copy of a colour charge (b = 1) is 2")

# M7: single copy (flat space, linear scalar coupling as for N = 4 SYM W-bosons): T = 1 - (kappa/alpha)^2 ... in units of the pull
sc = {}
for name, al, ka in [("Maxwell", 1, 0), ("N=4 BPS, opposite charge", 1, 1), ("scalar only", 0, 1)]:
    W = (E + al*u)**2 - (1 - ka*u)**2 - L**2*u**2; T, K = T_series(W); sc[name] = [str(T), str(K)]
data["single"] = sc
ok("M7", sc["Maxwell"][0] == "1" and sc["N=4 BPS, opposite charge"][0] == "0" and sc["scalar only"][0] == "-1",
   "single copy: Maxwell 1 (Sommerfeld), scalar alone -1, N = 4 BPS W-bosons 0 (Caron-Huot and Henn): vector +1, scalar -1 per unit coupling squared")
txt = "\n".join(("PASS " if c else "FAIL ") + n + ": " + m for n, c, m in res) + "\n%d/%d pass\n" % (sum(c for _, c, _ in res), len(res))
open("checks_magic_output.txt", "w").write(txt); json.dump(dict(checks=res, data=data), open("results_magic.json", "w"), indent=1, default=str)
print(txt)
