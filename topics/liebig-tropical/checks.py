"""Checks for liebig-tropical/README.md (lead R5): is Liebig's law of the minimum a stop on the paired-sum side?

Growth fraction for one nutrient: theta = 1/(1 + u), u = K/S (Monod). Two nutrients, written on the lift ruler
Y = 1/theta = 1 + u:
  multiplicative  Y = Y_A Y_B                      (the old dial at k = 1: u = u_A o_1 u_B)
  additive        Y = Y_A + Y_B - 1                 (the old dial at k = 0)
  Liebig          Y = max(Y_A, Y_B)                 (growth = the smaller theta)
"""
import json
import math

import numpy as np
import sympy as sp

results, store = [], {}


def check(name, ok, detail=""):
    results.append((name, bool(ok)))
    print(("PASS " if ok else "FAIL ") + name + (("   " + detail) if detail else ""))


uA, uB, k, a, p = sp.symbols("u_A u_B k alpha p", real=True)
tA, tB = sp.symbols("t_A t_B", positive=True)

# 1. reproduce the old closure
g = np.logspace(-2, 2, 4000)
ks = np.linspace(0, 1.4, 1401)
miss = [np.max(np.abs(1 / (1 + 2 * g + kk * g * g) - 1 / (1 + g))) for kk in ks]
ib = int(np.argmin(miss))
check("L1 old closure reproduced: on the one dial 1/theta = 1 + u_A + u_B + k u_A u_B, the best k is 0 and it "
      "still misses Liebig by 17% (both nutrients equally scarce)", abs(ks[ib]) < 1e-9 and abs(miss[ib] - 0.1716) < 1e-3,
      f"best k = {ks[ib]:.3f}, miss = {miss[ib]:.4f}")
store["dial_miss"] = {"k": ks[::10].tolist(), "miss": miss[::10]}

# 2. the old dial is the Ali-Mikhail-Haq copula with alpha = 1 - k
amh = tA * tB / (1 - a * (1 - tA) * (1 - tB))
dial = 1 / (1 + (1 / tA - 1) + (1 / tB - 1) + k * (1 / tA - 1) * (1 / tB - 1))
check("L2 the old one-dial growth law, written in theta, is exactly the Ali-Mikhail-Haq copula with alpha = 1 - k",
      sp.simplify(dial.subs(k, 1 - a) - amh) == 0)

# 3. the knob G_p(u) = (1 + u)^p - 1 is a symmetry of the multiplicative law (it commutes with o_1)
x, y = sp.symbols("x y", positive=True)
Gp = lambda v: (1 + v) ** p - 1
o1 = lambda m, n: m + n + m * n
check("L3 the rate knob G_p(u) = (1+u)^p - 1 commutes with o_1: G_p(x o_1 y) = G_p(x) o_1 G_p(y), so turning it "
      "leaves the multiplicative law unchanged", sp.simplify(sp.expand(Gp(o1(x, y)) - o1(Gp(x), Gp(y)))) == 0)

# 4. transporting the additive law by the knob gives Clayton's family
YA, YB = 1 / tA, 1 / tB
Y_knob = (YA ** p + YB ** p - 1) ** (1 / p)            # G_p^-1(G_p(u_A) + G_p(u_B)) on the lift ruler
clayton = (tA ** (-p) + tB ** (-p) - 1) ** (-1 / p)
check("L4 the additive law moved by the knob, u = G_p^-1(G_p(u_A) + G_p(u_B)), is exactly the Clayton copula "
      "theta = (theta_A^-p + theta_B^-p - 1)^(-1/p)", sp.simplify(1 / Y_knob - clayton) == 0)


def clay(ta, tb, pv):
    if pv == 0:
        return ta * tb
    if np.isinf(pv):
        return np.minimum(ta, tb)
    sa, sb = -np.log(ta) * pv, -np.log(tb) * pv          # stable: theta = exp(-(1/p) ln(e^sa + e^sb - 1))
    m = np.maximum(sa, sb)
    return np.exp(-(m + np.log(np.exp(sa - m) + np.exp(sb - m) - np.exp(-m))) / pv)


T = np.linspace(0.01, 1, 200)
TA, TB = np.meshgrid(T, T)
e0 = np.max(np.abs(clay(TA, TB, 1e-7) - TA * TB))
e1 = np.max(np.abs(clay(TA, TB, 1.0) - 1 / (1 / TA + 1 / TB - 1)))
errs = {pv: float(np.max(np.abs(clay(TA, TB, pv) - np.minimum(TA, TB)))) for pv in (1, 3, 10, 30, 100, 300)}
check("L5 one family holds all three classic laws: p -> 0 multiplicative, p = 1 additive, p -> infinity Liebig",
      e0 < 1e-5 and e1 < 1e-12 and errs[300] < 0.01 and all(errs[a_] > errs[b_] for a_, b_ in zip((1, 3, 10, 30, 100), (3, 10, 30, 100, 300))),
      "worst gap to Liebig: " + ", ".join(f"p={k_}: {v:.4f}" for k_, v in errs.items()))
store["clayton_miss"] = errs

# 6. Maslov: on the height ruler s = ln Y the family is a softened max at temperature 1/p, within T ln 2 of the max
worst_ratio = 0
for pv in (2, 5, 20, 80):
    s_true = np.log(1 / clay(TA, TB, pv)); s_max = np.log(np.maximum(1 / TA, 1 / TB))
    gap = np.max(s_true - s_max)
    worst_ratio = max(worst_ratio, gap / (math.log(2) / pv))
check("L6 Maslov: on the height ruler s = ln(1/theta) the family is h*ln(e^(s_A/h) + e^(s_B/h) - 1) with h = 1/p, "
      "which sits above max(s_A, s_B) by at most h ln 2 and cools to it", worst_ratio <= 1 + 1e-9,
      f"largest gap / (h ln 2) = {worst_ratio:.4f}")

# 7. Kendall's tau: AMH (the old dial) tops out at 1/3 at k = 0; Clayton has tau = p/(p+2), reaching 1 at Liebig


def tau_of(Cfun):
    # tau = 1 - 4 * int int dC/du dC/dv du dv (Nelsen 5.1.9), by finite differences on a grid
    n = 600
    uu = (np.arange(n) + 0.5) / n
    U, V = np.meshgrid(uu, uu, indexing="ij")
    h = 1e-6
    Cu = (Cfun(np.minimum(U + h, 1), V) - Cfun(np.maximum(U - h, 0), V)) / (2 * h)
    Cv = (Cfun(U, np.minimum(V + h, 1)) - Cfun(U, np.maximum(V - h, 0))) / (2 * h)
    return 1 - 4 * np.mean(Cu * Cv)


tau_amh = {kk: tau_of(lambda u_, v_, kk=kk: u_ * v_ / (1 - (1 - kk) * (1 - u_) * (1 - v_))) for kk in (2.0, 1.0, 0.5, 0.0)}
tau_cl = {pv: tau_of(lambda u_, v_, pv=pv: clay(u_, v_, pv)) for pv in (1.0, 2.0, 6.0)}
ok_tau = abs(tau_amh[0.0] - 1 / 3) < 3e-3 and abs(tau_amh[1.0]) < 3e-3 and all(abs(tau_cl[pv] - pv / (pv + 2)) < 3e-3 for pv in tau_cl)
check("L7 Kendall's tau (how tied the two limitations are): the old dial runs from -0.18 (k = 2) through 0 "
      "(multiplicative, k = 1) to at most 1/3 (additive, k = 0); the knob family has tau = p/(p+2), meeting the dial "
      "at 1/3 and reaching 1 at Liebig", ok_tau,
      "dial: " + ", ".join(f"k={k_}: {v:.3f}" for k_, v in tau_amh.items()) + "; knob: " + ", ".join(f"p={k_}: {v:.3f}" for k_, v in tau_cl.items()))
store["tau_amh"] = tau_amh; store["tau_clayton"] = tau_cl

# 8. the paired sum taken literally (o-plus of the lift ruler, a power sum with no unit fix) puts "no shortage" wrong
Ylit = lambda ya, yb, pv: (ya ** pv + yb ** pv) ** (1 / pv)
check("L8 the literal paired sum (power sum of 1/theta, no unit fix) gives growth 1/2 at temperature 1 when neither "
      "nutrient is short; its zero sits at the wall, not at 'no shortage'. The knob version has unit 0 = no shortage",
      abs(1 / Ylit(1.0, 1.0, 1.0) - 0.5) < 1e-15 and abs(clay(1.0, 1.0, 3.0) - 1) < 1e-15 and abs(clay(0.4, 1.0, 3.0) - 0.4) < 1e-15)

# 9. Kooijman's synthesizing unit (DEB theory): 1/theta = 1 + u_A + u_B - u_A u_B/(u_A + u_B)
su = lambda ta, tb: 1 / (1 / ta + 1 / tb - 1 - (1 / ta - 1) * (1 / tb - 1) / ((1 / ta - 1) + (1 / tb - 1) + 1e-300))
class _B: pass
best = _B()
_ps = np.linspace(0.5, 5, 901)
_ms = [np.max(np.abs(clay(TA, TB, pv) - su(TA, TB))) for pv in _ps]
best.x, best.fun = float(_ps[int(np.argmin(_ms))]), float(min(_ms))
best_dial = min((np.max(np.abs(1 / (1 / TA + 1 / TB - 1 + kk * (1 / TA - 1) * (1 / TB - 1)) - su(TA, TB))), kk) for kk in np.linspace(0, 1, 101))
tau_su = tau_of(lambda u_, v_: su(np.maximum(u_, 1e-9), np.maximum(v_, 1e-9)))
check("L9 Kooijman's synthesizing unit sits between additive and Liebig: closest knob setting p = %.2f misses by %.3f; "
      "the old dial can only reach it to %.3f; its tau is %.3f" % (best.x, best.fun, best_dial[0], tau_su),
      best.fun < best_dial[0] and 1 / 3 < tau_su < 1)
store["su"] = {"best_p": float(best.x), "miss": float(best.fun), "dial_miss": float(best_dial[0]), "tau": float(tau_su)}

# 10. growth with both nutrients equally scarce: knob family closes the 17% gap
gaps = {pv: float(np.max(np.abs(clay(1 / (1 + g), 1 / (1 + g), pv) - 1 / (1 + g)))) for pv in (1, 3, 10, 100)}
check("L10 on the old test (equal scarcity), the knob family's worst gap to Liebig starts at the old 17% (p = 1) and closes: "
      + ", ".join(f"p={k_}: {v:.3f}" for k_, v in gaps.items()),
      abs(gaps[1] - 0.1716) < 1e-3 and gaps[100] < 0.01 and gaps[3] < gaps[1] and gaps[10] < gaps[3])
store["equal_gaps"] = gaps

print("ALL PASS" if all(o for _, o in results) else "SOME FAIL", f"({sum(o for _, o in results)} of {len(results)})")
json.dump(store, open("checks.json", "w"), indent=1, default=float)
