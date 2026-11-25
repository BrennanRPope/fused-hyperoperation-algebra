"""Equilibrium and mass action through the circle operation.
1. Binding polynomials: two sites with interaction ω give P − 1 = a ⊞_ω b (the dial of subject 7).
2. Diacids: ω = 4·10^(−ΔpKa); the dial climbs from near-exclusion (oxalic) toward independence as the charges separate.
3. Two sites: positive cooperativity ⇔ ω > 1 ⇔ the binding polynomial's zeros leave the real axis (circle |x| = 1/√ω);
   Hill coefficient 2√ω/(1+√ω).
4. The n-fold dial is the MWC model with exclusive binding and L = k − 1.
5. Haber equilibrium: ammonia fraction against pressure at 700 K, ideal and with vdW fugacities.
Run: python3 checks.py"""
import json, math
import numpy as np
from scipy.optimize import brentq
ok = []; R = {}
def check(name, cond, info=""):
    ok.append(bool(cond)); print(("PASS " if cond else "FAIL ") + name + ("  " + info if info else ""))
def w(k, a, b): return a + b + k*a*b

# 1. binding polynomial
rng = np.random.default_rng(11)
good = all(abs((1 + K1*h + K2*h + om*K1*K2*h*h) - 1 - w(om, K1*h, K2*h)) < 1e-12 for K1, K2, om, h in rng.uniform(0, 3, (40, 4)))
check("two-site binding polynomial 1 + K₁h + K₂h + ωK₁K₂h²: P − 1 = (K₁h) ⊞_ω (K₂h)", good)

# 2. diacids (Wikipedia 'Dicarboxylic acid' table)
di = [("oxalic", 0, 1.27, 4.27), ("malonic", 1, 2.85, 5.05), ("succinic", 2, 4.21, 5.41), ("glutaric", 3, 4.34, 5.41),
      ("adipic", 4, 4.41, 5.41), ("pimelic", 5, 4.50, 5.43), ("suberic", 6, 4.526, 5.498), ("azelaic", 7, 4.550, 5.498), ("sebacic", 8, 4.720, 5.450)]
om = {n: 4*10**(-(p2 - p1)) for n, _, p1, p2 in di}
R["diacid_omega"] = om
# derivation check: macroscopic constants from intrinsic K and ω reproduce the pKa gap
K = 10**5.0
for name, o in om.items():
    beta1, beta2 = 2*K, o*K*K
    pKa2, pKa1 = math.log10(beta1), math.log10(beta2/beta1)
    assert abs((pKa2 - pKa1) - (math.log10(4) - math.log10(o))) < 1e-12
check("symmetric diacid: pKa₂ − pKa₁ = log 4 − log ω, so ω = 4·10^(−ΔpKa) (log 4 = 0.60 is the independent-site gap)", True)
vals = list(om.values())
check("the dial climbs from 0.004 (oxalic) to 0.75 (sebacic): from near-exclusion toward independence as the carboxylates separate",
      vals[0] < 0.01 and vals[-1] > 0.7 and all(o < 1 for o in vals), ", ".join(f"{n} {o:.3f}" for n, o in om.items()))
check("never past 1: like charges only repel, so k < 1 for every diacid", max(vals) < 1)

# 3. two sites: zeros and Hill coefficient
def nH2(o):
    # symmetric: P = 1 + 2x + o x²; Hill slope at half saturation
    Y = lambda x: (x + o*x*x)/(1 + 2*x + o*x*x)
    x0 = brentq(lambda x: Y(x) - 0.5, 1e-9, 1e9)
    h = 1e-6
    lg = lambda x: math.log(Y(x)/(1-Y(x)))
    return (lg(x0*(1+h)) - lg(x0*(1-h)))/(math.log(1+h) - math.log(1-h))
good = True
for o in (0.01, 0.25, 0.5, 1.0, 2.0, 9.0, 100.0):
    good &= abs(nH2(o) - 2*math.sqrt(o)/(1+math.sqrt(o))) < 1e-5
    roots = np.roots([o, 2, 1])
    real = np.all(np.abs(roots.imag) < 1e-12)
    good &= (real == (o <= 1))
    if o > 1: good &= np.allclose(np.abs(roots), 1/math.sqrt(o))
check("two sites: Hill coefficient 2√ω/(1+√ω); zeros real iff ω ≤ 1; for ω > 1 they sit on the circle |h| = 1/√ω", good)
circ = all(abs(abs(z)**2 + z.real) < 1e-12 for o in (1.2, 2, 5, 50) for z in np.roots([o, 2, 1]))
check("all complex zero pairs (ω > 1) lie on one circle |x|² = −Re x, through 0 and −1", circ)
check("so cooperativity (n_H > 1), complex zeros, and the dial past the circle operation (k > 1) are the same condition for two sites", nH2(1.0001) > 1 > nH2(0.9999))
R["nH_diacid"] = {n: nH2(o) for n, o in om.items()}

# 4. the n-fold dial = MWC exclusive binding with L = k − 1
def fold_poly(x, k, n): return 1 + ((1 + k*x)**n - 1)/k
def Y_fold(x, k, n, h=1e-7):
    P = lambda t: fold_poly(t, k, n)
    return x*(P(x*(1+h)) - P(x*(1-h)))/(2*x*h)/(n*P(x))
def Y_mwc(alpha, L, n): return alpha*(1+alpha)**(n-1)/(L + (1+alpha)**n)
good = True
for k in (1.5, 3, 10, 100):
    for n in (2, 3, 4):
        for x in (0.01, 0.1, 1, 5):
            good &= abs(Y_fold(x, k, n) - Y_mwc(k*x, k-1, n)) < 1e-6
check("the dial folded over n sites, 1 + ((1+kx)ⁿ − 1)/k, has exactly the saturation curve of the MWC model with binding only to R and L = k − 1",
      good, "n = 2, 3, 4; k = 1.5 … 100")
# coefficient structure: j ligands carry k^(j−1)
import sympy as sp
xs, ks = sp.symbols('x k')
P4 = sp.expand(1 + ((1 + ks*xs)**4 - 1)/ks)
coeffs = [sp.simplify(P4.coeff(xs, j)) for j in range(5)]
check("fold coefficients: C(n, j)·k^(j−1) for j ≥ 1: only the first ligand pays, every later one gains the factor k",
      coeffs == [1, 4, 6*ks, 4*ks**2, ks**3], str(coeffs))
def nH_fold(k, n):
    Y = lambda x: Y_fold(x, k, n)
    x0 = brentq(lambda x: Y(x) - 0.5, 1e-9, 1e6)
    h = 1e-4
    lg = lambda x: math.log(Y(x)/(1-Y(x)))
    return (lg(x0*(1+h)) - lg(x0*(1-h)))/(math.log(1+h) - math.log(1-h))
hill = {k: nH_fold(k, 4) for k in (1.0001, 3, 10, 100, 1000, 1e4)}
R["hill4"] = hill
check("four sites: Hill coefficient rises from 1 (k = 1) toward 4 as k grows; it passes 2.8 by k ≈ 10³",
      abs(hill[1.0001] - 1) < 1e-3 and hill[1000] > 2.8 and hill[1e4] < 4, ", ".join(f"k={k:g}: {v:.2f}" for k, v in hill.items()))

# 5. Haber equilibrium (K from subject 7, constant ΔH)
Rg = 8.314462618; dH = -91.80e3; dS = -198.11
def xNH3(T, P, Kphi=1.0):
    Kp = math.exp(-(dH - T*dS)/(Rg*T))/Kphi   # in bar^-2, for 2 NH3
    # stoichiometric feed 1 N2 : 3 H2, extent e: N2 1−e, H2 3−3e, NH3 2e, total 4−2e
    f = lambda e: ((2*e)**2*(4-2*e)**2)/((1-e)*(3-3*e)**3*P**2) - Kp
    e = brentq(f, 1e-12, 1-1e-12)
    return 2*e/(4-2*e)
y = {P: xNH3(700, P) for P in (1, 50, 100, 200, 300)}
y_real = xNH3(700, 200, 0.588)
R.update(haber_ideal=y, haber_200_vdw=y_real)
check("Haber at 700 K (ideal gas, constant ΔH): ammonia fraction rises with pressure (Le Chatelier)", all(a < b for a, b in zip(list(y.values()), list(y.values())[1:])),
      ", ".join(f"{P} bar {100*v:.1f}%" for P, v in y.items()))
check("at 200 bar the vdW fugacity factor from subject 6 lifts the fraction", y_real > y[200], f"{100*y[200]:.1f}% → {100*y_real:.1f}%")

print(f"\n{sum(ok)}/{len(ok)} checks pass")
json.dump(R, open("results.json", "w"), indent=1, default=float)
