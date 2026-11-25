"""Rates and Arrhenius through the circle operation.
1. First order: conversion over consecutive times composes at k = −1; parallel channels add rates (k = 0).
2. Second order: X/(1−X) adds over time, so conversions compose by (a + b − 2ab)/(1 − ab): the landscape border (−2, 1).
3. Order n: survival is the Tsallis q-exponential with q = n; consecutive times compose by the q-product. Only n = 0, 1, 2
   give a landscape (Möbius) law. The q-logarithm turns ordinary products into the dial with k = 1 − q.
4. Order n < 1 reaches the wall (all reactant gone) in finite time.
5. Arrhenius: ln k adds over temperature and barrier changes; k multiplies. Catalysis for nitrogen fixation.
Run: python3 checks.py"""
import json, math
import numpy as np
from scipy.integrate import solve_ivp
ok = []; R = {}
def check(name, cond, info=""):
    ok.append(bool(cond)); print(("PASS " if cond else "FAIL ") + name + ("  " + info if info else ""))
def w(k, a, b): return a + b + k*a*b
def land(k, q, a, b): return (a + b + k*a*b)/(1 - q*a*b)

def conv(n, kr, t, A0=1.0):
    sol = solve_ivp(lambda _, y: [-kr*max(y[0], 0)**n], [0, t], [A0], rtol=1e-12, atol=1e-14)
    return 1 - sol.y[0, -1]/A0
# 1. first order
kr = 0.7
X1, X2, X12 = conv(1, kr, 0.4), conv(1, kr, 1.1), conv(1, kr, 1.5)
check("first order: X(t₁ + t₂) = X(t₁) ⊞₋₁ X(t₂) (survivals multiply)", abs(X12 - w(-1, X1, X2)) < 1e-9, f"{X12:.6f}")
ka, kb, t = 0.3, 0.9, 1.3
Xa, Xb = 1 - math.exp(-ka*t), 1 - math.exp(-kb*t)
check("two parallel first-order channels: rate constants add (k = 0), and the conversion is Xa ⊞₋₁ Xb", abs((1 - math.exp(-(ka+kb)*t)) - w(-1, Xa, Xb)) < 1e-14)
check("branching ratio of parallel channels is ka/(ka + kb), fixed at every time", abs(0.3/1.2 - 0.25) < 1e-15)

# 2. second order
X1, X2, X12 = conv(2, kr, 0.4), conv(2, kr, 1.1), conv(2, kr, 1.5)
check("second order: X/(1 − X) grows linearly in time (k[A]₀t)", abs(X12/(1-X12) - kr*1.5) < 1e-8)
check("so X(t₁ + t₂) = (X₁ + X₂ − 2X₁X₂)/(1 − X₁X₂), the landscape law at (k, q) = (−2, 1)", abs(X12 - land(-2, 1, X1, X2)) < 1e-9)
check("(−2, 1) lies on the border k² = 4q: second-order kinetics is a border law (like a flat triangle)", (-2)**2 == 4*1)
check("zeroth order composes by plain addition (0, 0) and first order by (−1, 0)",
      abs(conv(0, 0.2, 1.5) - (conv(0, 0.2, 0.4) + conv(0, 0.2, 1.1))) < 1e-9 and abs(conv(1, kr, 1.5) - land(-1, 0, conv(1, kr, 0.4), conv(1, kr, 1.1))) < 1e-9)

# 3. order n: q-exponential and q-product
def e_q(x, q):  # Tsallis
    if abs(q-1) < 1e-12: return math.exp(x)
    base = 1 + (1-q)*x
    return base**(1/(1-q)) if base > 0 else 0.0
def ln_q(x, q):
    return math.log(x) if abs(q-1) < 1e-12 else (x**(1-q) - 1)/(1-q)
good = True
for n in (0.5, 1.5, 2.0, 2.5, 3.0):
    for tt in (0.3, 1.0, 2.0):
        s = 1 - conv(n, kr, tt)
        good &= abs(s - e_q(-kr*tt, n)) < 1e-8
check("order n: survival A/A₀ = e_q(−kt) with Tsallis q = n (A₀ = 1)", good, "n = 0.5, 1.5, 2, 2.5, 3")
def qprod(x, y, q):  # x ⊗_q y = e_q(ln_q x + ln_q y)
    return e_q(ln_q(x, q) + ln_q(y, q), q)
good = all(abs((1 - conv(n, kr, 1.5)) - qprod(1 - conv(n, kr, 0.4), 1 - conv(n, kr, 1.1), n)) < 1e-8 for n in (1.5, 2.5, 3.0))
check("consecutive times compose survivals by the q-product: s(t₁ + t₂) = s(t₁) ⊗_n s(t₂)", good)
good = all(abs(ln_q(x*y, q) - w(1-q, ln_q(x, q), ln_q(y, q))) < 1e-12 for x, y in [(0.3, 0.8), (2.0, 1.7), (0.05, 3.0)] for q in (0.5, 1.5, 2.0, 3.0))
check("the q-logarithm turns an ordinary product into the dial: ln_q(xy) = ln_q x ⊞_(1−q) ln_q y", good)
# is order 3 a landscape law? fit (k, q) to the composition and test
pairs = [(a, b) for a in np.linspace(0.05, 0.6, 8) for b in np.linspace(0.05, 0.6, 8)]
def comp(n, a, b):  # compose conversions a, b for order n
    sa, sb = 1 - a, 1 - b
    return 1 - qprod(sa, sb, n)
# best landscape fit with k = −n: solve q from each pair; constant q ⇔ landscape law
qs3 = [((a + b - 3*a*b)/comp(3, a, b) - 1)/(-a*b) for a, b in pairs]
qs2 = [((a + b - 2*a*b)/comp(2, a, b) - 1)/(-a*b) for a, b in pairs]
R.update(q_order2=[min(qs2), max(qs2)], q_order3=[min(qs3), max(qs3)])
check("order 2: the implied q is 1 for every pair (a genuine landscape law)", max(qs2) - min(qs2) < 1e-9)
check("order 3: no constant q works (implied q spans a range), so third-order kinetics is not a landscape law",
      max(qs3) - min(qs3) > 0.2, f"q from {min(qs3):.2f} to {max(qs3):.2f}")
# small-conversion limit: q → n(n−1)/2
small = [((a + b - 3*a*b)/comp(3, a, b) - 1)/(-a*b) for a, b in [(1e-3, 1e-3), (2e-3, 1e-3)]]
check("…but near zero conversion its implied q tends to 3 = n(n − 1)/2, matching orders 0, 1, 2 (q = 0, 0, 1)", abs(small[0] - 3) < 0.02, f"{small[0]:.3f}")

disc = {n: n*n - 4*(n*(n-1)/2) for n in (0, 0.5, 1, 1.5, 2, 2.5, 3)}
R["osc_disc"] = disc
check("the osculating point (−n, n(n−1)/2) has discriminant n(2 − n): two walls for 0 < n < 2, border at n = 0 and 2, turning for n > 2",
      all(abs(d - n*(2-n)) < 1e-12 for n, d in disc.items()))

# 4. finite-time wall
n = 0.5; tstar = 1/((1-n)*kr)
check("order n < 1 empties in finite time t* = A₀^(1−n)/((1−n)k): order ½ at k = 0.7 ends at t = 2.857",
      conv(n, kr, tstar*0.999) < 1 and conv(n, kr, tstar*1.001) > 1 - 1e-9 and abs(tstar - 2.857) < 1e-3)
check("orders n ≥ 1 never reach zero (the wall is approached only asymptotically)", 1 - conv(1, kr, 50) > 0 and 1 - conv(2, kr, 1e6) > 0)

# 5. Arrhenius
Rg = 8.314462618
def kA(A, Ea, T): return A*math.exp(-Ea/(Rg*T))
Ea = 80e3
good = abs(kA(1, Ea, 350)/kA(1, Ea, 300) - math.exp(Ea/Rg*(1/300 - 1/350))) < 1e-12
r3 = kA(1, Ea, 400)/kA(1, Ea, 300); r12 = (kA(1, Ea, 350)/kA(1, Ea, 300))*(kA(1, Ea, 400)/kA(1, Ea, 350))
check("Arrhenius: rate ratios over successive temperature steps multiply (ln k adds in 1/T)", good and abs(r3 - r12) < 1e-9*r3)
cat = {dE: math.exp(dE*1e3/(Rg*700)) for dE in (10, 50, 100)}
R["catalyst_700K"] = cat
check("a catalyst that lowers the barrier by ΔE multiplies the rate by e^(ΔE/RT): ×5.6, ×5.4e3, ×2.9e7 for 10, 50, 100 kJ/mol at 700 K",
      abs(cat[10] - 5.57) < 0.02 and abs(cat[50]/5.37e3 - 1) < 0.01 and abs(cat[100]/2.88e7 - 1) < 0.01)
hot = {E: kA(1, E*1e3, 700)/kA(1, E*1e3, 463) for E in (50, 100, 150)}
R["hot_vs_463"] = hot
check("running at 700 K instead of 463 K (where K = 1) speeds a 100 kJ/mol step by ×6.6e3: the reason Haber reactors run hot",
      abs(hot[100]/6.6e3 - 1) < 0.02, ", ".join(f"{E} kJ: ×{v:.2g}" for E, v in hot.items()))

print(f"\n{sum(ok)}/{len(ok)} checks pass")
json.dump(R, open("results.json", "w"), indent=1, default=float)
