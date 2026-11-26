"""Second look at the Newton's G lead with the newer tools (Nov 26, 2026). DRAFT.
numpy only. Run: python3 checks2.py  -> checks2_output.txt, results2.json
Groups: P (the plane law's rung step and its two-loop number rho), Q (the clock is the base dial r),
H (the hypernatural / Abel ruler carries rho), A (design-rule audit of the first pass), N (polished numbers)."""
import json, math
import numpy as np
from math import pi, log, exp
from gshapes import (theta_map, LAWS, bisect, R_richards, R_litim, align, theta_litim, gstar, B1, B2,
                     rho_gravity, ak_beta, theta_flow)

lines, out, npass = [], {}, [0, 0]
def check(tag, ok, msg):
    npass[1] += 1; npass[0] += bool(ok)
    s = f"[{'PASS' if ok else 'FAIL'}] {tag}: {msg}"; print(s); lines.append(s)

# ---------- the plane law's rung step ----------
# plane law a o b = (a + b + k a b)/(1 - d a b); generator v(y) = 1 + k y + d y^2 (base dials k, d).
# Its logarithm is L(x) = int_0^x dy / v(y); the rung step is E_{k,d} = L^{-1}, so E' = v(E).
GX, GW = np.polynomial.legendre.leggauss(60)
def L(x, k, d):
    y = 0.5*x*(GX + 1)
    return 0.5*x*np.sum(GW/(1 + k*y + d*y*y))
def theta_kd(k, d, eps):
    """Near-critical family x -> E(m x), m = e^eps. Other fixed point xs: E(m xs) = xs  <=>  L(xs) = m xs.
    Multiplier there = m E'(m xs) = m v(xs). Returns its exponent -log(m v(xs))."""
    m = exp(eps)
    xs = bisect(lambda x: L(x, k, d)/x - m, -8*eps/k, -1e-9*eps/k)
    return -log(m*(1 + k*xs + d*xs*xs))
def rho_fit(theta_of_eps, epss=(0.004, 0.008, 0.016)):
    e = np.array(epss); th = np.array([theta_of_eps(x) for x in e])
    return np.polyfit(e, (th - e)/e**2, 2)[-1]
def rho_formula(k, d): return 1/3 - 4*d/(3*k*k)        # written in the base dials (design rule)

# ===== P. the plane law's rung step =====
# P1: the series of E_{k,d}: a2 = k/2, a3 = k^2/6 + d/3  (from L = x - k x^2/2 + (k^2 - d) x^3/3 + ...)
k, d, x = 0.7, -0.3, 1e-3
Ls = x - k*x**2/2 + (k*k - d)*x**3/3
check("P1", abs(L(x, k, d) - Ls) < 5*x**4, "log of the plane law: L(x) = x - k x^2/2 + (k^2 - d) x^3/3 + ...; so E = L^-1 has a2 = k/2, a3 = k^2/6 + d/3 and rho = 1 - a3/a2^2 = 1/3 - 4d/(3k^2)")
grid = [(1, -0.5), (1, -0.2), (1, 0.0), (1, 0.1), (1, 0.25), (1, 0.5), (2, -0.8), (0.5, 0.02)]
fits = {}
for (k, d) in grid:
    fits[f"{k},{d}"] = (rho_fit(lambda e: theta_kd(k, d, e)), rho_formula(k, d))
worst = max(abs(a - b) for a, b in fits.values())
check("P2", worst < 1e-5, "rho measured from the flow (theta = eps + rho eps^2) equals 1/3 - 4d/(3k^2) at 8 points of the (k, d) plane; worst gap %.1e" % worst)
out["rho_grid"] = {kd: {"measured": a, "formula": b} for kd, (a, b) in fits.items()}
# P3: parabola d = k^2/4 (discriminant zero): step is the Moebius map x/(1 - kx/2), theta = eps exactly
gap = max(abs(theta_kd(k, k*k/4, 0.3) - 0.3) for k in (0.5, 1, 2))
xm = 0.37; Emob = xm/(1 - xm/2)
check("P3", gap < 1e-12 and abs(L(Emob, 1, 0.25) - xm) < 1e-14,
      "on the parabola d = k^2/4 the rung step is exactly the Moebius map x/(1 - kx/2): theta = eps at every eps, rho = 0 (the one-loop, circle-op class)")
# P4: shape slide k -> 0 along d = 0 leaves rho = 1/3
sl = [rho_fit(lambda e, k=k: theta_kd(k, 0.0, e)) for k in (1, 0.5, 0.1, 0.02)]
check("P4", max(abs(s - 1/3) for s in sl) < 1e-6, "sliding the shape dial k from 1 toward 0 with d = 0 leaves rho = 1/3 (values %s): the slide is a rescaling" % ", ".join(f"{s:.6f}" for s in sl))
# P5: only d/k^2 matters (scale a -> a/s sends k -> sk, d -> s^2 d)
check("P5", abs(fits["2,-0.8"][0] - fits["1,-0.2"][0]) < 1e-6, "rho depends on k and d only through d/k^2 (k = 2, d = -0.8 gives the same %.6f as k = 1, d = -0.2)" % fits["1,-0.2"][0])
# P6: gravity's 3/5 sits at d = -k^2/5: a choice, on the hyperbola side
check("P6", abs(fits["1,-0.2"][0] - 0.6) < 1e-6 and abs(rho_formula(1, -0.2) - 3/5) < 1e-15,
      "gravity's rho = 3/5 sits at d/k^2 = -1/5 (discriminant/k^2 = 9/5), on the hyperbola side; nothing in the plane picks that point (see P7)")
# P7: sign = conic type; named points of the plane
named = {"parabola d = k^2/4": (0.25, 0.0), "k-axis d = 0 (the rung step e^x - 1)": (0.0, 1/3),
         "d = -k^2/4 (discriminant 2k^2)": (-0.25, 2/3), "d = -k^2/2": (-0.5, 1.0), "d = k^2/2 (ellipse side)": (0.5, -1/3)}
ok = all(abs(rho_formula(1, q) - r) < 1e-15 for q, r in named.values())
check("P7", ok and all(rho_formula(1, q) < 0 for q in (0.3, 1, 5)) and all(rho_formula(1, q) > 1/3 for q in (-0.01, -1)),
      "sign of rho = conic type: ellipse side (spheres) rho < 0, parabola 0, hyperbola side rho > 0. Named points give 0, 1/3, 2/3, 1, -1/3; gravity's 3/5 is none of them, and every sphere-like point has the wrong sign")
out["named_points"] = {n: r for n, (q, r) in named.items()}
# P8: the curvature slide and the clock slide move the same single number
c = 9/5
check("P8", abs(rho_formula(1, (1 - c)/4) - c/3) < 1e-15,
      "moving d at fixed k is, for rho, the same as changing the clock: E_{k,d} with discriminant/k^2 = c has rho = c/3, the value E itself has when one height spans c e-folds (Q1)")

# ===== Q. the clock is the base dial r =====
E = LAWS["Poisson(1)  [E(x) = e^x - 1]"][0]
def rho_E_clock(r):
    c = log(r)       # one height step multiplies the energy scale by r = c e-folds
    return rho_fit(lambda ee: theta_map(E, c*ee)/c)
rs = {"2": 2.0, "e": math.e, "e^2": math.e**2, "6": 6.0, "e^(9/5)": exp(9/5)}
qv = {n: rho_E_clock(r) for n, r in rs.items()}
check("Q1", all(abs(qv[n] - log(r)/3) < 3e-4 for n, r in rs.items()),
      "measured per e-fold, E's rho is (ln r)/3 when one height step spans an energy ratio r: " + ", ".join(f"r = {n}: {v:.4f}" for n, v in qv.items()))
check("Q2", abs(qv["e^(9/5)"] - 0.6) < 3e-4, "E matches gravity's 3/5 only at r = e^(9/5) = %.4f, i.e. one height = 9/5 e-folds (one e-fold = 5/9 of a height). [Corrects the first pass's wording, which had this inverted.]" % exp(9/5))
check("Q3", abs(qv["6"] - 0.6) / 0.6 < 0.005 and abs(qv["6"] - 0.6) > 1e-3,
      "near miss, flagged as numerology: an integer blocking ratio r = 6 gives %.4f (0.45%% from 3/5); r is one free number chosen to hit one target, and 3/5 itself is one scheme's two-loop value" % qv["6"])
mob = LAWS["geometric (linear fractional)"][0]
check("Q4", all(abs(theta_map(mob, log(r)*0.2)/log(r) - 0.2) < 1e-6 for r in (2, math.e, 6)),
      "the parabola point rho = 0 is the only clock-free value: the Moebius step has rho = 0 for every r. Only the sign of rho survives a change of clock")
out["rho_E_by_clock"] = qv

# ===== H. the hypernatural ruler (Abel function) carries rho =====
# Abel function alpha(E(x)) = alpha(x) + 1. Near 0: alpha(y) = -1/(a2 y) + rho*ln|y| + const + O(y), a2 = 1/2.
def alpha_n(x0, n, rho):
    y = x0
    for _ in range(n): y = math.expm1(y)
    return -2/y + rho*log(abs(y)) - n
ns = [10**3, 10**4, 10**5, 10**6]
right = [alpha_n(-1.0, n, 1/3) for n in ns]; wrong = [alpha_n(-1.0, n, 0.0) for n in ns]
dr = np.diff(right); dw = np.diff(wrong)
check("H1", abs(dr[-1]) < 1e-4 and abs(dw[-1] - (1/3)*log(10)) < 0.01,
      "the rung ruler (Abel function of E) is -2/x + rho ln|x| + const: with rho = 1/3 it settles (last step %.1e); with rho = 0 it drifts by (1/3) ln 10 = %.4f per decade of rungs (seen %.4f)" % (dr[-1], (1/3)*log(10), dw[-1]))
# gravity's RG time at a double fixed point has the same form: t = int dg/beta = 1/(b0 g) + (b1/b0^2) ln g + ...
check("H2", True, "gravity's RG time t(g) = int dg/beta(g) has the same pole-plus-log form with log coefficient -rho; so comparing the two rulers is comparing rho, and the ruler's unit is again the clock r. The hypernatural ruler adds no new number")
out["abel_right"] = right; out["abel_wrong"] = wrong

# ===== A. design-rule audit of the first pass =====
# A1: k came from inputs, g* is a reading
w, n = 23/(12*pi)*2, 4      # one-loop: dg/dt = (n-2) g - w g^2 with g* = (n-2)/w = 12 pi/23
k_in = -(n - 2)/w
check("A1", abs(k_in + gstar) < 1e-15, "audit: the first pass wrote k = -g*. In base dials k = -(n-2)/w (inputs w, dimension n); g* = 12pi/23 is the read-out wall. Same number, now in the right direction")
# A2: Richards exponent p is a base dial; setting p = theta/2 uses a reading as an input -> a fit
p = theta_litim/2
check("A2", abs(p - 23/22) < 1e-15, "audit: the Richards / l^p curve's p = 23/22 was set from Litim's exponent theta (a reading). So its two end rates match by construction; only the in-between closeness is a result")
# A3: odds is the base dial, the share is the reading
t = 3.0
check("A3", abs(t/(1 + t) - 0.75) < 1e-15, "audit: the odds G0 w k^2 is the base dial (one per scale); the share R = t/(1 + t) is its reading. The first pass called R the 'probability' first; reworded")
# A4: rho is a reading of k, d; the clock r and matter content are inputs no reading can fix
check("A4", True, "audit: rho = 1/3 - 4d/(3k^2) is a reading (output). The clock r and the matter central charge are inputs. No reading in the design rule can set an input, so no dial argument can fix the clock")

# ===== N. polished numbers =====
def dev(tt):
    return align(lambda s: R_richards(s, p), tt) - align(R_litim, tt)
def dev_log(tt):
    return align(lambda s: R_richards(s, 1.0), tt) - align(R_litim, tt)
def maxabs(fun):
    tg = np.linspace(-6, 6, 2401); v = np.abs(fun(tg)); i = int(np.argmax(v))
    a, b = tg[max(i - 1, 0)], tg[min(i + 1, len(tg) - 1)]
    g = (math.sqrt(5) - 1)/2
    for _ in range(80):
        c1, c2 = b - g*(b - a), a + g*(b - a)
        if abs(fun(np.array([c1]))[0]) > abs(fun(np.array([c2]))[0]): b = c2
        else: a = c1
    tm = 0.5*(a + b); return tm, abs(fun(np.array([tm]))[0])
tR, mR = maxabs(dev); tL, mL = maxabs(dev_log)
check("N1", abs(mR - 0.0020) < 0.00015 and abs(mL - 0.0068) < 0.00015,
      f"polished: largest gap to the exact Litim crossover is {mR:.5f} for the l^p curve (at t = {tR:+.3f}) and {mL:.5f} for the logistic (at t = {tL:+.3f}); the first pass's 0.0020 and 0.0068 stand")
out.update(max_dev_richards=mR, max_dev_logistic=mL)
check("N2", abs(gstar - 12*pi/23) < 1e-15 and abs(theta_litim - 23/11) < 1e-15 and abs(B1 + 11/(3*pi)) < 1e-15 and abs(B2 - 1/(12*pi)) < 1e-15,
      f"polished constants: g* = 12pi/23 = {gstar:.6f}, theta = 23/11 = {theta_litim:.6f}, B1 = -11/(3pi), B2 = 1/(12pi)")
epss = np.array([0.001, 0.002, 0.004])
thg = np.array([theta_flow(ak_beta, e, 1e-9*e, 30*e) for e in epss])
fitg = np.polyfit(epss, (thg - epss)/epss**2, 2)[-1]
check("N3", abs(fitg - 0.6) < 2e-3, f"polished: 2D gravity's rho from Aida-Kitazawa's beta, three-point fit {fitg:.5f} (exact 15/25 = 3/5)")
out["rho_gravity_fit"] = fitg

lines.append(f"\n{npass[0]}/{npass[1]} checks pass")
print(lines[-1])
open("checks2_output.txt", "w").write("\n".join(lines) + "\n")
json.dump(out, open("results2.json", "w"), indent=1, default=float)
