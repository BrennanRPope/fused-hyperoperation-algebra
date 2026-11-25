"""Acids and bases through the circle operation.
1. Logs: pH turns products of constants into sums; the fraction deprotonated is a logistic in pH; independent groups
   combine at k = −1. The isoelectric point is an average of pKas.
2. Buffers: capacities add (k = 0); mixing two buffers of one acid combines their base/acid ratios by the mediant.
3. Triprotic acids, a real test: the dial folded over three sites (subject 8's MWC form) forces the second pKa gap to be
   log 3 = 0.48, whatever k is. A uniform pairwise interaction instead forces equal gaps. Tabulated pKas for citric,
   phosphoric and arsenic acids have near-equal gaps, far from 0.48: the data choose the pairwise picture, not the fold.
pKa values: one general-chemistry table (OpenStax-style Appendix): citric 3.13, 4.76, 6.40; phosphoric 2.14, 7.20, 12.34;
arsenic 2.21, 6.93, 11.52. Arsenic values differ by about 0.2 between sources.
Run: python3 checks.py"""
import json, math
import numpy as np
from scipy.optimize import brentq
ok = []; R = {}
def check(name, cond, info=""):
    ok.append(bool(cond)); print(("PASS " if cond else "FAIL ") + name + ("  " + info if info else ""))
def w(k, a, b): return a + b + k*a*b
Kw = 1e-14

# 1. logs, logistic, independence, pI
frac = lambda pH, pKa: 1/(1 + 10**(pKa - pH))
check("Henderson–Hasselbalch: pH = pKa + log(base/acid) is the logistic fraction 1/(1 + 10^(pKa − pH))",
      all(abs(pH - (4.76 + math.log10(frac(pH, 4.76)/(1 - frac(pH, 4.76))))) < 1e-12 for pH in (3, 4.76, 6.2)))
check("Ka multiplies when steps are chained (pKa adds): log turns the lift into addition", abs(-math.log10(10**-3.13*10**-4.76) - (3.13 + 4.76)) < 1e-12)
a1, a2 = frac(7.4, 6.0), frac(7.4, 9.0)
check("two independent groups: chance at least one is deprotonated = α₁ ⊞₋₁ α₂", abs(1 - (1 - a1)*(1 - a2) - w(-1, a1, a2)) < 1e-15)
def gly_charge(pH, p1=2.34, p2=9.60):
    return (1 - frac(pH, p2)) - frac(pH, p1)  # NH3+ minus COO−
pI = brentq(gly_charge, 3, 9); R["glycine_pI"] = pI
check("glycine (pKa 2.34, 9.60): net charge zero at pH 5.97, the plain average of the two pKas", abs(pI - 5.97) < 1e-9, f"pI = {pI:.4f}")

# 2. buffers
def pH_solution(Ca, Cb, pKa):
    # weak acid HA (Ca) plus its sodium salt (Cb): charge balance [H] + [Na] = [OH] + [A]
    Ka = 10**-pKa; C = Ca + Cb
    f = lambda pH: 10**-pH + Cb - Kw/10**-pH - C*Ka/(Ka + 10**-pH)
    return brentq(f, 0, 14)
def capacity(Ca, Cb, pKa, d=1e-7):
    # β = dCb/dpH by finite difference on added strong base at fixed total
    p0 = pH_solution(Ca, Cb, pKa); p1 = pH_solution(Ca - d, Cb + d, pKa)
    return d/(p1 - p0)
bA = capacity(0.05, 0.05, 4.76); bB = capacity(0.02, 0.08, 7.20)
pH_mixA = pH_solution(0.05, 0.05, 4.76)
def capacity_mix(d=1e-7):
    def pHmix(extra):
        f = lambda pH: (10**-pH + 0.05 + 0.08 + extra - Kw/10**-pH
                        - 0.10*10**-4.76/(10**-4.76 + 10**-pH) - 0.10*10**-7.20/(10**-7.20 + 10**-pH))
        return brentq(f, 0, 14)
    return pHmix
# capacity is a local slope; two buffers in one solution: β_total = β_A(pH) + β_B(pH) − β_water (water counted once)
def beta_formula(pH, comps):
    h = 10**-pH; b = math.log(10)*(h + Kw/h)
    for C, pKa in comps:
        Ka = 10**-pKa; b += math.log(10)*C*Ka*h/(Ka + h)**2
    return b
pm = capacity_mix(); p0 = pm(0); p1 = pm(1e-7); bmix = 1e-7/(p1 - p0)
check("buffer capacities add (k = 0): the mixture's β equals the sum of each acid's term at the shared pH (water counted once)",
      abs(bmix - beta_formula(p0, [(0.10, 4.76), (0.10, 7.20)]))/bmix < 1e-4, f"β = {bmix:.4f} M")
# mixing two buffers of the same acid: ratios combine by the mediant
c1 = (0.08, 0.02); c2 = (0.01, 0.09)  # (acid, base) moles in equal volumes
med = (c1[1] + c2[1])/(c1[0] + c2[0])
pH_mix = pH_solution((c1[0] + c2[0])/2, (c1[1] + c2[1])/2, 4.76)
check("mixing two acetate buffers: base/acid ratios combine by the mediant (b₁ + b₂)/(a₁ + a₂), not by the circle operation",
      abs(pH_mix - (4.76 + math.log10(med))) < 0.01, f"pH {pH_mix:.3f} vs {4.76 + math.log10(med):.3f}")
check("the mediant lies between the two ratios (0.25 and 9)", c1[1]/c1[0] < med < c2[1]/c2[0])

# 3. triprotic acids: fold vs pairwise
acids = {"citric": (3.13, 4.76, 6.40), "phosphoric": (2.14, 7.20, 12.34), "arsenic": (2.21, 6.93, 11.52)}
gaps = {n: (p[1] - p[0], p[2] - p[1]) for n, p in acids.items()}; R["gaps"] = gaps
def fold_poly(k):  # 1 + ((1 + kx)^3 − 1)/k: coefficients of x, x², x³
    return (3, 3*k, k*k)
def pair_poly(k):  # three equal sites, factor k per occupied pair
    return (3, 3*k, k**3)
def gaps_from(beta):  # protonation β_j (x scaled out); stepwise log K_j; pKa order reversed
    b = (1,) + beta
    logK = [math.log10(b[j]/b[j-1]) for j in (1, 2, 3)]
    return (logK[0] - logK[1], logK[1] - logK[2])
ok_fold = all(abs(gaps_from(fold_poly(k))[1] - math.log10(3)) < 1e-12 for k in (1e-6, 1e-3, 0.1, 1, 7))
check("the dial folded over three sites fixes one pKa gap at log 3 = 0.477, whatever k is", ok_fold)
ok_pair = all(abs(gaps_from(pair_poly(k))[0] - gaps_from(pair_poly(k))[1]) < 1e-12 for k in (1e-6, 1e-3, 0.1, 1, 7))
check("a uniform pairwise interaction gives two equal gaps, each log(3/k)", ok_pair and abs(gaps_from(pair_poly(0.01))[0] - math.log10(300)) < 1e-12)
check("tabulated gaps: citric 1.63 / 1.64, phosphoric 5.06 / 5.14, arsenic 4.72 / 4.59: each pair within 0.15 of equal",
      all(abs(g[0] - g[1]) < 0.15 for g in gaps.values()), "; ".join(f"{n} {g[0]:.2f}/{g[1]:.2f}" for n, g in gaps.items()))
check("…and every gap is at least 1.6, far from the fold's 0.48", all(min(g) > 1.6 for g in gaps.values()))
kpair = {n: 3*10**(-(g[0] + g[1])/2) for n, g in gaps.items()}; R["k_pair"] = kpair
eng = {n: -8.314462618*298.15*math.log(k)/1000 for n, k in kpair.items()}; R["pair_kJ"] = eng
check("fitted pair factors: citric k ≈ 0.07 (6.6 kJ/mol repulsion per pair), phosphoric ≈ 2.4e-5 (26 kJ/mol), arsenic ≈ 6.6e-5 (24 kJ/mol)",
      abs(kpair["citric"] - 0.0698) < 0.002 and abs(kpair["phosphoric"]/2.4e-5 - 1) < 0.06 and abs(kpair["arsenic"]/6.6e-5 - 1) < 0.03,
      "; ".join(f"{n} {k:.2g}, {eng[n]:.1f} kJ/mol" for n, k in kpair.items()))
# predict the third pKa of each acid from the first two under the pairwise picture
pred = {n: p[1] + (p[1] - p[0]) for n, p in acids.items()}; R["pKa3_pred"] = pred
check("pairwise picture predicts pKa₃ from pKa₁ and pKa₂: citric 6.39 (table 6.40), phosphoric 12.26 (12.34), arsenic 11.65 (11.52)",
      all(abs(pred[n] - acids[n][2]) < 0.15 for n in acids), "; ".join(f"{n} {v:.2f}" for n, v in pred.items()))

# titration curve of phosphoric acid for the picture
def ntitr(pH, pKas=(2.14, 7.20, 12.34)):
    h = 10**-pH; Ks = [10**-p for p in pKas]
    terms = [1, Ks[0]/h, Ks[0]*Ks[1]/h**2, Ks[0]*Ks[1]*Ks[2]/h**3]
    Z = sum(terms); return sum(j*t for j, t in enumerate(terms))/Z
R["phos_curve"] = [[ntitr(p), p] for p in np.linspace(0.5, 13.5, 131)]
check("phosphoric titration: half-equivalence points sit at the three pKas (protons removed 0.5, 1.5, 2.5)",
      all(abs(ntitr(p) - (j + 0.5)) < 0.02 for j, p in enumerate((2.14, 7.20, 12.34))))

print(f"\n{sum(ok)}/{len(ok)} checks pass")
json.dump(R, open("results.json", "w"), indent=1, default=float)
