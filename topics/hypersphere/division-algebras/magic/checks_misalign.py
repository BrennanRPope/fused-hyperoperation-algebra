"""DRAFT. Misaligned charges: bodies whose charges point different ways, orbiting N = 8's four-charge black holes.
Background: STU truncation (stu.py), charges k_i on four U(1)s. Visitor: components j with weights w_j (sum 1) and
signs s_j (+1 same direction as the black hole's charge j, -1 opposite). Units pi K^2/L^2 as before.
Writes results_misalign.json and checks_misalign_output.txt."""
import json, math, sympy as sp
from stu import check as stu_check
from misalign import W_stu
from precess import T_series, T_numeric, u, L, E
res = []; data = {}
def ok(name, cond, msg):
    res.append((name, bool(cond), msg)); print(("PASS " if cond else "FAIL ") + name + ": " + msg, flush=True)
one = sp.Integer(1)
ein, sc, mx = stu_check([one, sp.Integer(2), sp.Integer(3), sp.Rational(1, 2)], 1)
ok("X1", max(ein, sc, mx) < 1e-12, "four-charge black hole (charges 1, 2, 3, 1/2) satisfies Einstein, scalar and Maxwell equations (worst %.1e)" % max(ein, sc, mx))
good = True
for n in range(1, 5):
    ks = [one]*n + [0]*(4-n); ws = [sp.Rational(1, n)]*n + [0]*(4-n); ss = [-1]*4
    T, K = T_series(W_stu(ks, ws, ss)); good &= (T == n-1 and K == 2)
ok("X2", good, "n equal charges, visitor exactly opposite: T = n - 1 = 0, 1, 2, 3 again (matches the dilaton black holes of section 8)")
# X3: the pair rule, symbolic
k = sp.symbols('k1:5', positive=True); w = sp.symbols('w1:5', positive=True); s = sp.symbols('s1:5')
W = W_stu(list(k), list(w), list(s))
F = sp.series((W/L**2).subs(E, 1), u, 0, 4).removeO()
c1 = sp.limit(sp.expand(F.coeff(u, 1)*L**2), L, sp.oo)/2; c2 = -F.coeff(u, 2); c3 = sp.limit(sp.expand(F.coeff(u, 3)), L, sp.oo)
d = sp.limit(sp.expand((c2-1)*L**2), L, sp.oo); sub = {w[3]: 1-w[0]-w[1]-w[2]}
num = sp.expand((-d + 3*c3*c1).subs(sub))
x = [w[j]*(1-s[j]) for j in range(4)]
rule = 2*sum(k[i]*k[j]*(x[i]*(1-w[j]) + x[j]*(1-w[i]) + x[i]*x[j]) for i in range(4) for j in range(i+1, 4))
diff = sp.expand(num - rule.subs(sub))
for si in s: diff = diff.subs(si**2, 1)
kok = sp.expand(c1.subs(sub) - sum(k[j]*x[j] for j in range(4)).subs(sub)) == 0
ok("X3", sp.expand(diff) == 0 and kok, "pair rule, exact for any charges and any visitor: T K^2 = 2 sum_{i<j} k_i k_j [x_i(1-w_j) + x_j(1-w_i) + x_i x_j], K = sum k_j x_j, x_j = w_j(1-s_j)")
# X4: one charge: W exactly quadratic in u with -L^2 u^2: orbits close exactly, all orders
k1 = sp.Symbol('k1', positive=True); allq = True
import itertools
for pat in itertools.product([0, 1, -1], repeat=4):
    if all(p == 0 for p in pat): continue
    nn = sum(1 for p in pat if p); ws = [sp.Rational(abs(p), nn) for p in pat]; ss = [p if p else 1 for p in pat]
    Wq = sp.expand(sp.simplify(W_stu([k1, 0, 0, 0], ws, ss)))
    allq &= Wq.is_polynomial(u) and sp.Poly(Wq, u).degree() <= 2 and sp.expand(Wq.coeff(u, 2) + L**2) == 0
ok("X4", allq, "a one-charge (half-BPS) black hole: for all 80 visitor charge mixes the radial equation is exactly Kepler's (quadratic in 1/r, angular term untouched), so orbits close exactly at every order")
# X5: half-BPS visitor flipped against charge j of a general black hole: T = (sum of the other charges)/k_j
ks = [one, sp.Integer(2), sp.Integer(3), sp.Rational(1, 2)]; good = True; tab = {}
for j in range(4):
    ws = [0]*4; ws[j] = one; ss = [1]*4; ss[j] = -1
    T, K = T_series(W_stu(ks, ws, ss)); pred = (sum(ks)-ks[j])/ks[j]; tab[j+1] = str(T); good &= sp.simplify(T - pred) == 0
ok("X5", good, "one-charge visitor opposite charge j of a black hole with charges (1, 2, 3, 1/2): T = (other charges)/k_j = %s" % ", ".join(tab.values()))
# X6: simulated orbits agree with the series
sims = {}; worst = 0
for ks, ws, ss in [([one]*4, [0, 0, one/2, one/2], [1, 1, 1, -1]), ([one, one, one, 0], [one/3]*3+[0], [-1, -1, 1, 1]), ([one, sp.Integer(2), sp.Integer(3), one/2], [one, 0, 0, 0], [-1, 1, 1, 1])]:
    Wx = W_stu(ks, ws, ss); T, K = T_series(Wx); Tn = T_numeric(Wx, K, Lv=400); sims[str((ks, ws, ss))] = [str(T), Tn]
    worst = max(worst, abs(Tn - float(T))/max(1, abs(float(T))))
ok("X6", worst < 2e-3, "simulated orbits agree with the rule for three mixed cases (worst relative gap %.1e)" % worst)
# direction scan for figure: two equal charges, visitor charge direction psi in the (Q1, Q2) plane
scan = []
for i in range(361):
    psi = math.radians(i); c, sn = math.cos(psi), math.sin(psi); a = abs(c) + abs(sn)
    w1, w2 = abs(c)/a, abs(sn)/a; s1 = 1 if c >= 0 else -1; s2 = 1 if sn >= 0 else -1
    x1, x2 = w1*(1-s1), w2*(1-s2); K = x1 + x2; coef = 2*(x1*(1-w2) + x2*(1-w1) + x1*x2)
    scan.append([i, K, coef/K**2 if K > 1e-12 else None])
data.update(sims=sims, X5=tab, scan=scan)
txt = "\n".join(("PASS " if c else "FAIL ") + n + ": " + m for n, c, m in res) + "\n%d/%d pass\n" % (sum(c for _, c, _ in res), len(res))
open("checks_misalign_output.txt", "w").write(txt); json.dump(dict(checks=res, data=data), open("results_misalign.json", "w"), indent=1, default=str)
print(txt)
