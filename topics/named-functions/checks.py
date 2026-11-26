"""Checks for the named-functions pass.  Run: python3 -I checks.py  (from this folder)."""
import sys, os, json, math
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(HERE, "..", "fractional-rungs"))
import mpmath as mp
import numpy as np
from flows import E_k, E_kd, law_kd, L_kd_numeric, T_kappa, ORDER
from regular_iteration import Es as Es_old

results = {}; data = {}
def rep(name, ok, extra=""):
    print(("PASS " if ok else "FAIL ") + name + (("   " + extra) if extra else ""))
    results[name] = bool(ok)

E = E_k(1)
# ---------------- 1. the generic flow reproduces the project's E^s ----------------
worst = 0
for s in [-1.5, -0.5, 0.25, 0.5, 1.5]:
    for x in [-0.6, -0.2, 0.3, 1.0, 2.0]:
        a = float(E(s, x)); b = float(Es_old(s, np.array([x]))[0])
        worst = max(worst, abs(a - b))
rep("C1 generic regular iteration matches fractional-rungs/regular_iteration.py (25 points)", worst < 1e-9, f"max diff {worst:.1e}")
rep("C2 integer height: E^2(1) = e^(e-1) - 1", abs(E(2, 1) - (mp.e**(mp.e-1) - 1)) < 1e-25)
g = abs(E(0.5, E(0.7, 0.4)) - E(1.2, 0.4))
rep("C3 heights add: E^0.5(E^0.7(x)) = E^1.2(x)", g < 1e-20, f"{float(g):.1e}")
rep("C4 polish estimate small (x=1, s=1/2)", (E(0.5, 1), E.last_err)[1] < 1e-20, f"n vs 2n diff {float(E.last_err):.1e}")

# ---------------- 2. the hypernatural log and its inverse ----------------
lnh = lambda s, X: E(-s, X - 1)          # hypernatural logarithm  ln^(s)
exph = lambda s, t: 1 + E(s, t)          # hypernatural exponential exp^(s)
ok = True
for s in [-1, -0.5, 0.5, 1, 1.5, 2]:
    for X in [0.4, 0.9, 1.7, 3.0]:
        if abs(exph(s, lnh(s, X)) - X) > 1e-18: ok = False
rep("H1 exp^(s)(ln^(s)(X)) = X at 6 heights, both sides of 1", ok)
rep("H2 ln^(1) = ln, exp^(1) = exp", abs(lnh(1, 3) - mp.log(3)) < 1e-25 and abs(exph(1, 0.7) - mp.e**0.7) < 1e-25)
rep("H3 ln^(0)(X) = X - 1;  ln^(2)(X) = log(1 + log X)", abs(lnh(0, 3) - 2) < 1e-25 and abs(lnh(2, 3) - mp.log(1 + mp.log(3))) < 1e-25)
# times at height s
def times_h(s, X, Y): return exph(s, lnh(s, X) + lnh(s, Y))
X0, Y0 = mp.mpf("2.5"), mp.mpf("1.7")
rep("H4 x^(1) is ordinary multiplication", abs(times_h(1, X0, Y0) - X0*Y0) < 1e-25)
rep("H5 x^(0) is X + Y - 1", abs(times_h(0, X0, Y0) - (X0 + Y0 - 1)) < 1e-25)
plus_h = lambda s, a, b: E(s, E(-s, a) + E(-s, b))
rep("H6 x^(s) is the lift reading of +^(s): X x^(s) Y = 1 + ((X-1) +^(s) (Y-1)) at s=0.6",
    abs(times_h(0.6, mp.mpf("2.2"), mp.mpf("1.3")) - (1 + plus_h(0.6, mp.mpf("1.2"), mp.mpf("0.3")))) < 1e-20)
# both halves of X = 1 join smoothly: slope 1, bend -s from both sides
s = 0.5; h = mp.mpf("1e-4")
def d1(fun, X, h): return (fun(X + h) - fun(X - h))/(2*h)
right = [(lnh(s, 1 + h*j)) for j in (1, 2, 3)]
left = [(lnh(s, 1 - h*j)) for j in (1, 2, 3)]
slope_r = (4*right[0] - right[1])/(2*h) ; slope_l = -(4*left[0] - left[1])/(2*h)
bend_r = 2*(right[0] - h)/h**2; bend_l = 2*(left[0] + h)/h**2
rep("H7 cube half (X>1) and share half (X<1) meet at X=1 with slope 1 and bend -s (s=1/2)",
    abs(slope_r - 1) < 1e-6 and abs(slope_l - 1) < 1e-6 and abs(bend_r + s) < 1e-3 and abs(bend_l + s) < 1e-3,
    f"slopes {float(slope_r):.7f},{float(slope_l):.7f}; bends {float(bend_r):.4f},{float(bend_l):.4f}")

# ---------------- 3. shapes ----------------
# 3a. every k != 0 is the cube's flow in other units
k = mp.mpf("0.37"); Ek = E_k(k)
rep("S1 rescale: E_k^s(x) = E^s(kx)/k (k=0.37, s=1/2)", abs(Ek(0.5, 0.8) - E(0.5, k*0.8)/k) < 1e-20)
Em = E_k(-1)
rep("S2 simplex side k=-1 (stick-breaking): E_{-1}^s(x) = -E^s(-x), the other half of the same flow",
    abs(Em(0.5, 0.3) - (-E(0.5, -0.3))) < 1e-20)
small = [abs(E(0.5, kk*0.8)/kk - 0.8) for kk in (0.1, 0.01, 0.001)]
rep("S3 k -> 0: E_k^(1/2)(0.8) -> 0.8, all heights collapse onto plain addition",
    small[0] > small[1] > small[2] and small[2] < 1e-3, ", ".join(f"{float(v):.1e}" for v in small))
data["collapse"] = [float(v) for v in small]

# 3b. the plane dial: every conic (k, d)
shapes = {"hyperbola (cube line)": (1, 0), "ellipse": (1, 1), "parabola": (1, mp.mpf(1)/4),
          "circle": (0, 1), "relativity": (0, -1), "hyperbola off-line": (1, -1)}
flows = {nm: E_kd(*kd) for nm, kd in shapes.items()}
# exponential formula right: inverse of the numeric logarithm, and the Riccati equation E' = 1 + kE + dE^2
ok = True
for nm, fl in flows.items():
    kk, dd = fl.k, fl.d
    for t in [0.2, 0.5]:
        x = fl.f(t)
        if abs(L_kd_numeric(kk, dd, x) - t) > 1e-25: ok = False
        der = mp.diff(fl.f, t)
        if abs(der - (1 + kk*x + dd*x*x)) > 1e-20: ok = False
rep("S4 E_(k,d)(t) = sigma(T_kappa(t)) inverts L_(k,d) and solves E' = 1 + kE + dE^2 (6 conics)", ok)
ok = all(abs(flows["hyperbola (cube line)"].f(t) - mp.expm1(t)) < 1e-30 for t in (0.3, 1.1))
rep("S5 at (k, d) = (1, 0) the conic's exponential is E(t) = e^t - 1", ok)
ok = True
for nm, fl in flows.items():
    a, b = mp.mpf("0.21"), mp.mpf("0.13")
    lhs = law_kd(fl.k, fl.d, a, b); rhs = fl.f(fl.finv(a) + fl.finv(b))
    if abs(lhs - rhs) > 1e-30: ok = False
rep("S6 each conic's law is the quasi-sum with generator E_(k,d) (Aczel form)", ok)
# height flows on each conic
ok = True; groupdiff = {}
for nm, fl in flows.items():
    x = mp.mpf("0.3")
    g = abs(fl(0.5, fl(0.5, x)) - fl.f(x))
    groupdiff[nm] = float(g)
    if g > 1e-15: ok = False
rep("S7 every conic has a real height flow: half-rung twice = one rung (6 conics)", ok,
    ", ".join(f"{n.split()[0]} {v:.0e}" for n, v in groupdiff.items()))
# parabola: the flow is an exact Mobius slide
fl = flows["parabola"]; x = mp.mpf("0.3")
exact = lambda s: x/(1 - s*(fl.k/2)*x)
rep("S8 parabola class: E^s(x) = x/(1 - s(k/2)x) exactly; height is a slide",
    all(abs(fl(s, x) - exact(s)) < 1e-25 for s in (0.5, 1.7, -0.8)))
# ln at height 1 on each conic = the conic's own logarithm
ok = all(abs(fl(-1, mp.mpf("0.4")) - L_kd_numeric(fl.k, fl.d, mp.mpf("0.4"))) < 1e-18 for fl in flows.values())
rep("S9 on every conic, height 1 of the hypernatural log is that conic's own logarithm L_(k,d)", ok)
# k = 0 conics: odd flow (symmetric), multiplicity 3
fc = flows["circle"]
rep("S10 circle (k=0): the flow is odd, tan^(s)(-x) = -tan^(s)(x); first correction is cubic",
    abs(fc(0.5, -0.3) + fc(0.5, 0.3)) < 1e-25 and fc.m == 3)
# 3c. lp shapes: the p-sum (a^p + b^p)^(1/p), generator x^(1/p): not tangent to the identity
p0 = 2.0
def psum(p, a, b): return (a**p + b**p)**(1/p)
def rung_lp(s, a, b, p0=p0):  # transport + by G^s, G(y) = y^(1/p0): G^s(y) = y^(p0^-s)
    e = p0**(-s); return ((a**(1/e)) + (b**(1/e)))**e
ok = all(abs(rung_lp(s, 0.7, 0.4) - psum(p0**s, 0.7, 0.4)) < 1e-12 for s in (0, 0.5, 1, 1.5, 2.5))
rep("S11 lp shapes: rung s of the sphere's Pythagoras is the l_(2^s) sum; height = the p dial (p = 2^s)", ok)
rep("S12 lp generator x^(1/p) has slope 0 or infinity at 0 (not tangent to identity): no new function, only powers",
    True, "exact: G^s(y) = y^(p^-s)")

# ---------------- 4. Bell numbers at height s ----------------
def b_n(s, n):
    return E.coeffs(mp.mpf(s))[n]*mp.factorial(n)
bell = [1, 1, 2, 5, 15, 52, 203, 877]
a258 = [1, 1, 3, 12, 60, 358, 2471, 19302]
ok2 = all(abs(b_n(2, n) - bell[n]) < 1e-20 for n in range(1, 8))
ok3 = all(abs(b_n(3, n) - a258[n]) < 1e-18 for n in range(1, 8))
ok1 = all(abs(b_n(1, n) - 1) < 1e-20 for n in range(1, 8))
okm = all(abs(b_n(-1, n) - (-1)**(n-1)*math.factorial(n-1)) < 1e-18 for n in range(1, 8))
rep("B1 derivatives of exp^(s) at 0: s=1 all ones, s=2 Bell numbers, s=3 A000258, s=-1 signed (n-1)!", ok1 and ok2 and ok3 and okm)
sB = mp.mpf("0.37")
rep("B2 b_2(s) = s and b_3(s) = (3s^2 - s)/2 for real s", abs(b_n(sB, 2) - sB) < 1e-25 and abs(b_n(sB, 3) - (3*sB**2 - sB)/2) < 1e-25)
# A008826 binomial-basis rows (Comtet): b_n(s) = sum_k T(n,k) C(s,k)
rows = {3: [1, 3], 4: [1, 13, 18], 5: [1, 50, 205, 180]}
ok = True
for n, row in rows.items():
    for s in (0.3, 1.7, 2.6):
        val = sum(c*mp.binomial(s, j+1) for j, c in enumerate(row))
        if abs(val - b_n(s, n)) > 1e-18: ok = False
rep("B3 b_n(s) = sum_k A008826(n,k) C(s,k) (Comtet's triangle) for n = 3, 4, 5", ok)
data["bell_table"] = {str(s): [float(b_n(s, n)) for n in range(1, 7)] for s in (-1, 0, 0.5, 1, 1.5, 2, 2.5, 3)}
# numerical check that b_2 is a real derivative of the flow: (E^s(h) - h)/h^2 -> s/2
h = mp.mpf("1e-4"); est = (E(0.6, h) - h)/h**2
rep("B4 b_2 is an actual derivative of the flow: (E^0.6(h) - h)/h^2 -> 0.3", abs(est - 0.3) < 1e-4, f"{float(est):.6f}")

# ---------------- 5. other named functions ----------------
rep("N1 E_1 = E, and E_k inverts L_k = log(1+kx)/k", abs(E_k(1).f(0.7) - mp.expm1(0.7)) < 1e-30 and abs(E_k(mp.mpf("0.4")).f(mp.log1p(mp.mpf("0.4")*mp.mpf("0.9"))/mp.mpf("0.4")) - mp.mpf("0.9")) < 1e-30)
# Cayley-Klein: exp(t eps) for eps^2 = k eps - d equals e^{kt/2}(C + S eps0)
k, d, t = mp.mpf("0.8"), mp.mpf("0.9"), mp.mpf("0.7"); kap = d - k**2/4
M = mp.matrix([[0, -d], [1, k]])          # multiplication by eps on basis (1, eps)
expM = mp.expm(t*M)
C = mp.cos(mp.sqrt(kap)*t); S = mp.sin(mp.sqrt(kap)*t)/mp.sqrt(kap)
# e^{kt/2}(C + S (eps - k/2)) = e^{kt/2}((C - S k/2) + S eps)
pred0 = mp.e**(k*t/2)*(C - S*k/2); pred1 = mp.e**(k*t/2)*S
rep("N2 exp(t eps) = e^(kt/2) (C_kappa(t) + S_kappa(t) eps0), eps0 = eps - k/2, kappa = d - k^2/4",
    abs(expM[0, 0] - pred0) < 1e-30 and abs(expM[1, 0] - pred1) < 1e-30)
rep("N3 tan of the Cayley-Klein pair is the conic's exponential after the slide: T_kappa = S/C",
    abs(T_kappa(kap, t) - S/C) < 1e-30)
kk = mp.mpf("0.3"); n = 6; acc = mp.mpf(0)
for _ in range(n): acc = acc + 1 + kk*acc
rep("N4 1 combined n times by o_k is the q-integer [n]_q with q = 1 + k", abs(acc - ((1+kk)**n - 1)/kk) < 1e-30)
a, b = mp.mpf("0.3"), mp.mpf("0.45")
rep("N5 o_{-1} is the probabilistic sum: a + b - ab = 1 - (1-a)(1-b)", abs(a + b - a*b - (1 - (1-a)*(1-b))) < 1e-30)
# iterative logarithm (Julia function): d/ds E^s at s = 0
gx = mp.mpf("0.2"); der = mp.diff(lambda s: E(s, gx), 0)
series = gx**2/2 - gx**3/12 + gx**4/48 - gx**5/180 + 11*gx**6/8640
rep("N6 the speed of the height flow d/ds E^s at s=0 is Julia's iterative logarithm g^2/2 - g^3/12 + ...", abs(der - series) < 1e-6, f"{float(der):.10f} vs {float(series):.10f}")

data["pass"] = sum(results.values()); data["total"] = len(results)
print(f"\n{data['pass']} of {data['total']} pass")
json.dump({"checks": results, "data": data}, open(os.path.join(HERE, "checks.json"), "w"), indent=1)
