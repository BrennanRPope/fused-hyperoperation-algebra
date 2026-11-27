"""Second order (2PN) and all orders around four-charge black holes. Writes results_order2.json, checks_order2_output.txt.
Exact orbit equation (poly.py): L^2 (du/dphi)^2 = P(u) - L^2 u^2, P = sum_m p_m u^m, p_m built from sets of m charges.
Predictions: T = p2/K^2,  T2 = 3/4 T^2 + 3 p3/K^3,  T2e = -(dp2/dE)/K^2   (all at E = 1, K = p1/2)."""
import sys, json, itertools, random; sys.path.insert(0, '..')
import sympy as sp, mpmath as mp
from misalign import W_stu
from precess import T_series, u, L, E
from turn2 import turn2
from poly import poly_F
lines = []; ok_all = True; data = {}
def rep(name, ok, detail=""):
    global ok_all; ok_all &= bool(ok); lines.append(f"{'PASS' if ok else 'FAIL'}  {name}  {detail}"); print(lines[-1])
one = sp.Integer(1)
# Z1 Schwarzschild against Damour-Schafer test-mass 2PN: turn/pi = 6/j^2 + (105/2 - 15 w)/j^4
M = one; Wsch = (1 - 2*M*u)*(E**2/(1 - 2*M*u) - 1 - L**2*u**2)
rep("Z1 second-order turn formula reproduces Einstein's test-mass result (Damour-Schafer): 6, 105/2, -15", turn2(Wsch, M) == (6, sp.Rational(105, 2), -15), str(turn2(Wsch, M)))
# Z2 exact polynomial orbit equation, all 16 sign patterns, symbolic charges and weights
k = sp.symbols('k1:5', positive=True); w = sp.symbols('w1:5', positive=True); random.seed(1); ok2 = True
for ss in itertools.product([1, -1], repeat=4):
    W = W_stu(list(k), list(w), list(ss)); P = poly_F(list(k), list(w), list(ss))
    vals = {**{kk: sp.Rational(random.randint(1, 9), random.randint(1, 5)) for kk in k}, **{ww: sp.Rational(random.randint(1, 9), 7) for ww in w}, u: sp.Rational(1, 13), E: sp.Rational(9, 10), L: 3}
    ok2 &= abs(sp.N((W - P).subs(vals), 30)) < 1e-20 and sp.degree(P, u) == 4
rep("Z2 the exact orbit equation is a polynomial in 1/r of degree = number of charges (all 16 sign patterns)", ok2)
def pcoef(ks, ws, ss):
    P = sp.expand(poly_F(ks, ws, ss) + L**2*u**2); return [P.coeff(u, m) for m in range(5)]
def exact_turn(ks, ws, ss, Lv, Ev, dps=40):
    mp.mp.dps = dps; P = pcoef(ks, ws, ss)
    cs = [mp.mpf(sp.N(c.subs(E, Ev), dps)) for c in P]; Lv = mp.mpf(sp.N(Lv, dps)); cs[2] -= Lv**2
    poly = [cs[4], cs[3], cs[2], cs[1], cs[0]]
    while poly[0] == 0: poly = poly[1:]
    roots = sorted([mp.re(r) for r in mp.polyroots(poly, maxsteps=200, extraprec=200) if abs(mp.im(r)) < mp.mpf(10)**(-dps//2) and mp.re(r) > 0])
    u1, u2 = roots[0], roots[1]; mid = (u1 + u2)/2; half = (u2 - u1)/2
    f = lambda t: mp.mpf(Lv)*half*mp.sin(t)/mp.sqrt(abs(mp.re(mp.polyval(poly, mid - half*mp.cos(t)))))
    return mp.re(2*mp.quad(f, [0, mp.pi]) - 2*mp.pi)
# Z3 two charges: exact rosette, turn = 2 pi (1/sqrt(1 - p2/L^2) - 1) at every order (strong field)
ok3 = True; z3 = []
for ks, ws, ss in [([one, 2*one, 0, 0], [one, 0, 0, 0], [-1, 1, 1, 1]), ([one, one, 0, 0], [one/3, 2*one/3, 0, 0], [-1, 1, 1, 1]), ([3*one, one, 0, 0], [one/2, 0, one/2, 0], [-1, 1, -1, 1])]:
    Pk = pcoef(ks, ws, ss); Kk = Pk[1].subs(E, 1)/2; Lv = 3*Kk; Ev = sp.sqrt(1 - sp.Rational(3, 5)*Kk**2/Lv**2); p2 = Pk[2].subs(E, Ev)
    pred = 2*mp.pi*(1/mp.sqrt(1 - mp.mpf(sp.N(p2, 40))/mp.mpf(sp.N(Lv, 40))**2) - 1); got = exact_turn(ks, ws, ss, Lv, Ev)
    ok3 &= abs(got - pred) < mp.mpf(10)**-15; z3.append((str(ks), float(got), float(pred)))
rep("Z3 around a two-charge black hole the orbit is an exact precessing ellipse; its turn is closed-form at all orders", ok3, str(z3))
# Z4 one charge: zero turn at every order (strong field)
ks4, ws4, ss4 = [2*one, 0, 0, 0], [one/2, one/4, one/4, 0], [-1, 1, -1, 1]
Pk = pcoef(ks4, ws4, ss4); Kk = Pk[1].subs(E, 1)/2; Lv = 3*Kk; Ev = sp.sqrt(1 - sp.Rational(3, 5)*Kk**2/Lv**2)
got = exact_turn(ks4, ws4, ss4, Lv, Ev)
rep("Z4 around a one-charge black hole the turn is exactly zero in strong field too", abs(got) < mp.mpf(10)**-15, mp.nstr(got, 5))
# Z5 2PN formulas versus the general second-order expansion, many cases
bgs = {"2": [one, one, 0, 0], "2u": [one, 2*one, 0, 0], "3": [one, one, one, 0], "3u": [one, 2*one, 4*one, 0], "4": [one]*4, "4u": [one, 2*one, 3*one, 5*one]}
probes = {"anti1": ([one, 0, 0, 0], [-1, 1, 1, 1]), "anti12": ([one/2, one/2, 0, 0], [-1, -1, 1, 1]), "anti_all": ([one/4]*4, [-1]*4), "anti1_al2": ([one/2, one/2, 0, 0], [-1, 1, 1, 1]), "anti1_x3": ([one/2, 0, one/2, 0], [-1, 1, 1, 1])}
ok5 = True; rows = []
for bn, ks in bgs.items():
    for pn, (ws, ss) in probes.items():
        if not any(ss[j] == -1 and ws[j] != 0 and ks[j] != 0 for j in range(4)): continue
        P = pcoef(ks, ws, ss); K = sp.simplify(P[1].subs(E, 1)/2); p2 = P[2].subs(E, 1); p3 = P[3].subs(E, 1)
        T = sp.simplify(p2/K**2); T2p = sp.simplify(sp.Rational(3, 4)*T**2 + 3*p3/K**3); T2ep = sp.simplify(-sp.diff(P[2], E).subs(E, 1)/K**2)
        W = W_stu(ks, ws, ss); T1, T2, T2e = turn2(W, K)
        good = sp.simplify(T1 - T) == 0 and sp.simplify(T2 - T2p) == 0 and sp.simplify(T2e - T2ep) == 0; ok5 &= good
        rows.append(dict(bg=bn, probe=pn, K=str(K), T=str(T), T2=str(T2), T2e=str(T2e), pairs=str(p2), triples=str(p3)))
rep("Z5 second order: T2 = 3/4 T^2 + 3 p3/K^3 (p3 = triples of charges only) and T2e = -(dp2/dE)/K^2, all cases", ok5, f"{len(rows)} cases")
data["order2_cases"] = rows
# Z6 numeric: exact integration at large L matches the two-term series (3 charges)
ks, ws, ss = [one, 2*one, 4*one, 0], [one, 0, 0, 0], [-1, 1, 1, 1]
P = pcoef(ks, ws, ss); K = P[1].subs(E, 1)/2; T = P[2].subs(E, 1)/K**2
T2 = sp.Rational(3, 4)*T**2 + 3*P[3].subs(E, 1)/K**3; T2e = -sp.diff(P[2], E).subs(E, 1)/K**2
mp.mp.dps = 50; res = []
for Lv in [200, 400]:
    wv = sp.Rational(2, 5); Ev = sp.sqrt(1 - 2*wv*K**2/sp.Integer(Lv)**2)
    got = exact_turn(ks, ws, ss, Lv, Ev, dps=50); x = mp.mpf(sp.N(K, 50))**2/Lv**2
    resid = (got/mp.pi - mp.mpf(sp.N(T, 50))*x)/x**2
    res.append(float(resid))
pred = float(T2 + T2e*sp.Rational(2, 5))
rep("Z6 exact integration of a three-charge orbit at large L matches the second-order prediction", abs(res[-1] - pred) < 0.01*abs(pred) and abs(res[-1]-pred) < abs(res[0]-pred), f"pred {pred:.5f}, L=200: {res[0]:.5f}, L=400: {res[1]:.5f}")
# Z7 third rung: the product of all four charges (p4) first moves the turn at third order, by (p4/K^4)(15/2 - 3w) K^6/L^6
def turn_from_poly(cs, Lv, dps=60):
    mp.mp.dps = dps; cs = list(cs); cs[2] -= Lv**2; poly = [cs[4], cs[3], cs[2], cs[1], cs[0]]
    while poly[0] == 0: poly = poly[1:]
    roots = sorted([mp.re(r) for r in mp.polyroots(poly, maxsteps=300, extraprec=300) if abs(mp.im(r)) < mp.mpf(10)**(-dps//2) and mp.re(r) > 0])
    u1, u2 = roots[0], roots[1]; mid = (u1 + u2)/2; half = (u2 - u1)/2
    f = lambda t: Lv*half*mp.sin(t)/mp.sqrt(abs(mp.re(mp.polyval(poly, mid - half*mp.cos(t)))))
    return mp.re(2*mp.quad(f, [0, mp.pi]) - 2*mp.pi)
ks, ws, ss = [one, 2*one, 3*one, 5*one], [one, 0, 0, 0], [-1, 1, 1, 1]
P = pcoef(ks, ws, ss); K = P[1].subs(E, 1)/2; p4 = P[4].subs(E, 1); wv = sp.Rational(2, 5)
pred7 = float(p4/K**4*(sp.Rational(15, 2) - 3*wv)); res7 = []
for Lv in [100, 200, 400]:
    Ev = sp.sqrt(1 - 2*wv*K**2/sp.Integer(Lv)**2); mp.mp.dps = 60
    cs = [mp.mpf(sp.N(c.subs(E, Ev), 60)) for c in P]; Lm = mp.mpf(Lv)
    full = turn_from_poly(cs, Lm); cut = turn_from_poly(cs[:4] + [mp.mpf(0)], Lm)
    res7.append(float((full - cut)/mp.pi*Lm**6/mp.mpf(sp.N(K, 60))**6))
rep("Z7 the product of all four charges first moves the turn at third order, as the exact equation says (four unequal charges)",
    abs(res7[-1] - pred7) < 0.01*abs(pred7) and abs(res7[-1] - pred7) < abs(res7[0] - pred7),
    f"pred {pred7:.5f}, L=100: {res7[0]:.5f}, L=200: {res7[1]:.5f}, L=400: {res7[2]:.5f}")
data["rung3"] = dict(p4=str(p4), K=str(K), w=str(wv), pred=pred7, measured=res7)
lines.append("ALL PASS" if ok_all else "SOME FAILED"); print(lines[-1])
json.dump(dict(passed=ok_all, data=data), open("results_order2.json", "w"), indent=1)
open("checks_order2_output.txt", "w").write("\n".join(lines) + "\n")
