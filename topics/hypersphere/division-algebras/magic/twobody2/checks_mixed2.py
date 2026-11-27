"""Checks: two comparable four-charge (STU) black holes whose charges point different ways, at second order (2PN),
from the lifted Lagrangian (multi.py: Gupta's halves lifted to four photons and three scalars, one static term corrected).
Writes checks_mixed2_output.txt and results_mixed2.json.  Turns as for a small visitor:
turn/pi = T k^2/J^2 + (T2 + T2e w) k^4/J^4, w = -E J^2/(mu k^2), k = Newtonian pull, X = m_A/(m_A + m_B), nu = X(1 - X)."""
import sys, os, json; HERE = os.path.dirname(os.path.abspath(__file__))
for p in (HERE, os.path.join(HERE, '..', 'order2'), os.path.join(HERE, '..')): sys.path.insert(0, p)
import sympy as sp, numpy as np
import multi, pn2
from poly import poly_F
from precess import u, L, E
R = sp.Rational; lines = []; ok_all = True; data = {}
def rep(name, ok, detail=""):
    global ok_all; ok_all &= bool(ok); lines.append(f"{'PASS' if ok else 'FAIL'}  {name}  {detail}"); print(lines[-1], flush=True)
def probe(A, B):
    """exact small-visitor turns of B round A (section 11 orbit equation)."""
    sA = 1 if sum(A.Q) > 0 else -1
    ks = [abs(x) for x in A.Q]; ss = [(1 if x*sA > 0 else -1) if x != 0 else 1 for x in B.Q]
    P = sp.expand(poly_F(ks, B.w, ss) + L**2*u**2); c = [P.coeff(u, m) for m in range(5)]
    K = c[1].subs(E, 1)/2; T = c[2].subs(E, 1)/K**2
    return sp.nsimplify(T), sp.nsimplify(R(3, 4)*T**2 + 3*c[3].subs(E, 1)/K**3), sp.nsimplify(-sp.diff(c[2], E).subs(E, 1)/K**2)
def scaled(Q, s): return multi.Body([sp.nsimplify(q)*s for q in Q])
Xs = sp.Symbol('X')
pairs = [((1, 1, 1, 0), (0, 0, -1, -1)), ((1, 0, 0, 0), (-1, -1, 0, 0)), ((1, 1, 0, 0), (-1, -1, -1, -1)), ((3, 1, 0, 0), (0, -1, -2, 0)),
         ((1, 2, 3, 0), (-2, 0, -1, -1)), ((1, 1, 0, 0), (-2, -1, 0, 0)), ((1, 1, 0, 0), (0, -1, -1, 0))]
ok1 = ok2 = ok3 = True; det1 = []; det3 = []
for QA, QB in pairs:
    A0 = multi.Body(QA); B0 = multi.Body(QB); pBA = probe(A0, B0); pAB = probe(B0, A0)
    pts = []
    for x in (1, R(1, 2), R(1, 3), R(1, 5)):
        A = A0; B = scaled(QB, x*A0.m/B0.m)
        T, T2, T2e, rest = multi.turns(A, B); assert rest == 0
        pts.append((A.m/(A.m + B.m), sp.nsimplify(T), sp.nsimplify(T2), sp.nsimplify(T2e)))
    fit = {nm: sp.expand(sp.interpolate([(p[0], p[i]) for p in pts], Xs)) for nm, i in [("T", 1), ("T2", 2), ("T2e", 3)]}
    # small-visitor limits: X -> 1 is B round A, X -> 0 is A round B (cubic through four comparable-mass points)
    ends = all(sp.simplify(fit[nm].subs(Xs, 1) - pBA[i]) == 0 and sp.simplify(fit[nm].subs(Xs, 0) - pAB[i]) == 0 for i, nm in enumerate(["T", "T2", "T2e"]))
    ok1 &= ends; nu = Xs*(1 - Xs)
    Tl = Xs*pBA[0] + (1 - Xs)*pAB[0]
    ok2 &= sp.expand(fit["T"] - Tl) == 0 and sp.expand(fit["T2e"] - (Xs*pBA[2] + (1 - Xs)*pAB[2]) - nu*Tl) == 0
    rem = sp.cancel(sp.expand(fit["T2"] - (Xs*pBA[1] + (1 - Xs)*pAB[1]))/nu)
    kap = sp.simplify(rem + R(3, 2)*Tl); ok3 &= not kap.has(Xs)
    det1.append(f"{QA}|{QB}: B round A {pBA}, A round B {pAB}"); det3.append(f"{QA}|{QB}: kappa = {kap}")
    data[f"{QA}|{QB}"] = dict(probe_B_round_A=[str(v) for v in pBA], probe_A_round_B=[str(v) for v in pAB], kappa=str(kap),
                              points=[[str(v) for v in p] for p in pts], T=str(fit["T"]), T2=str(fit["T2"]), T2e=str(fit["T2e"]))
rep("M1 both small-visitor limits come out exactly: the comparable-mass turns, fitted at four mass ratios, end on the exact orbit-equation values", ok1, "; ".join(det1))
rep("M2 first order and the energy part of second order follow mass shares: T = X T(B|A) + (1-X) T(A|B), T2e = X T2e(B|A) + (1-X) T2e(A|B) + nu T", ok2)
rep("M3 energy-free part: T2 = X T2(B|A) + (1-X) T2(A|B) + nu (kappa - 3T/2) with kappa the same at every mass ratio", ok3, "; ".join(det3))
# M4 the same-direction case: kappa = 3 (4 - n)/2 for n equal charges each (from the one-field section, checks_twobody2.py)
ok4 = True; det4 = []
for n in (1, 2, 3, 4):
    QA = tuple([1]*n + [0]*(4 - n)); QB = tuple([-1]*n + [0]*(4 - n)); A = multi.Body(QA); B = multi.Body(QB)
    T, T2, T2e, rest = multi.turns(A, B); pr = probe(A, B)
    kap = sp.simplify((T2 - pr[1])/R(1, 4) + R(3, 2)*T); ok4 &= kap == R(3, 2)*(4 - n) and sp.simplify(T2e - pr[2] - T/4) == 0
    det4.append(f"n={n}: kappa {kap}")
rep("M4 same-direction pairs through the multi-field equations: kappa = 3(4 - n)/2, T2e = probe + nu T (equal masses), as the one-field equations give", ok4, "; ".join(det4))
# M5 simulate one mixed pair at equal masses and measure the turn per lap
from scipy.integrate import solve_ivp
from scipy.optimize import brentq
def build(A, B):
    L1, L2 = multi.full(A, B, 0, 0); k = multi.pull(A, B); qq = sum(x*y for x, y in zip(A.Q, B.Q))
    sub = {pn2.G: 1, pn2.m1: A.m, pn2.m2: B.m, pn2.G12: (k + qq)/(A.m*B.m), pn2.q1: 1, pn2.q2: qq}
    H0, H1, H2, kk = pn2.hamiltonian(L1, L2, sub); r = pn2.r
    pr, J = sp.symbols('pr J', real=True); H = (H0 + H1 + H2).subs({sp.Symbol('px'): pr, sp.Symbol('py'): J/r})
    lf = lambda e: sp.lambdify((r, pr, J), e, 'numpy')
    return lf(H), lf(sp.diff(H, pr)), lf(sp.diff(H, J)), lf(-sp.diff(H, r)), float(kk), float(A.m*B.m/(A.m + B.m))
def measure(model, Jk, ecc=0.3, laps=6):
    f, dr, dphi, dpr, k, mu = model; J = Jk*k; E0 = -(1 - ecc**2)/2*mu*k**2/J**2
    rp0 = J**2/(mu*k*(1 + ecc)); rp = brentq(lambda x: f(x, 0.0, J) - E0, 0.6*rp0, 1.3*rp0)
    rhs = lambda t, y: [dr(y[0], y[1], J), dpr(y[0], y[1], J), dphi(y[0], y[1], J)]
    ev = lambda t, y: y[1]; ev.direction = 1
    Tk = 2*np.pi*(J**2/(mu*k*(1 - ecc**2)))**1.5*np.sqrt(mu/k)
    sol = solve_ivp(rhs, [0, laps*1.2*Tk], [rp, 0.0, 0.0], rtol=1e-12, atol=1e-12, events=ev, dense_output=True, max_step=Tk/400)
    te = [t for t in sol.t_events[0] if t > 0.1*Tk]; ph = [sol.sol(t)[2] for t in te]
    return float(np.mean(np.diff([0.0] + ph) - 2*np.pi))
QA, QB = (1, 1, 0, 0), (0, -1, -1, 0); A = multi.Body(QA); B = multi.Body(QB)
T, T2, T2e, rest = multi.turns(A, B); wv = (1 - 0.3**2)/2; pred = float(T2 + T2e*sp.nsimplify(wv))
mdl = build(A, B); got = [(measure(mdl, Jk) - np.pi*float(T)/Jk**2)/(np.pi/Jk**4) for Jk in (8, 16)]
ok5 = abs(got[1] - pred) < 0.02*abs(pred) and abs(got[1] - pred) < abs(got[0] - pred)
data["simulated"] = dict(pair=f"{QA}|{QB}", pred=pred, J8=got[0], J16=got[1])
rep("M5 simulated equal-mass orbit of a mixed pair lands on the predicted second-order turn", ok5, f"{QA}|{QB}: predicted {pred:.4g}, J/k=8: {got[0]:.4g}, J/k=16: {got[1]:.4g}")
lines.append("ALL PASS" if ok_all else "SOME FAILED"); print(lines[-1])
json.dump(dict(passed=ok_all, data=data), open(os.path.join(HERE, "results_mixed2.json"), "w"), indent=1)
open(os.path.join(HERE, "checks_mixed2_output.txt"), "w").write("\n".join(lines) + "\n")
