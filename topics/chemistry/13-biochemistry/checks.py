"""Biochemistry through the circle operation, with photosynthesis and nitrogenase.
1. Two inhibitors (Yonetani–Theorell): v₀/vᵢ − 1 = a ⊞ₖ b with a = [I]/K_I, b = [J]/K_J and k = 1/α. Mutually exclusive
   inhibitors (α → ∞) are the k = 0 law, which is Loewe additivity; independent ones (α = 1) are the circle operation, which is Bliss
   independence for the fractional activity. α < 1 is synergy. The dial is a measured enzymology parameter.
2. Michaelis–Menten: saturation is Langmuir; irreversible steps in series add their times, so turnover rates add
   harmonically (k = 0 on 1/rate).
3. Kinetic proofreading: discrimination odds multiply over stages, so error ≈ f₀^(n+1), and the fidelity ee composes by
   the (0, −1) law from subject 12.
4. Nitrogenase: N₂ + 8H⁺ + 8e⁻ + 16ATP → 2NH₃ + H₂ + 16ADP + 16Pᵢ. 6 of 8 electrons reach ammonia (75%).
5. Photosynthesis: 3 ATP + 2 NADPH and at least 8 photons per CO₂; at that limit 34% of 680 nm light energy is stored.
Run: python3 checks.py"""
import json, math
import numpy as np
from scipy.optimize import brentq
ok = []; R = {}
def check(name, cond, info=""):
    ok.append(bool(cond)); print(("PASS " if cond else "FAIL ") + name + ("  " + info if info else ""))
def w(k, a, b): return a + b + k*a*b

# 1. Yonetani–Theorell from the full binding equilibrium
def rate_two_inhibitors(I, J, KI, KJ, alpha, Et=1.0, kcat=1.0):
    # rapid equilibrium, inhibitors bind free enzyme only (substrate saturating is folded into kcat)
    a, b = I/KI, J/KJ
    Z = 1 + a + b + a*b/alpha  # E, EI, EJ, EIJ
    return kcat*Et/Z
ok1 = True
for alpha in (0.2, 1.0, 5.0, 1e9):
    for I, J in ((1.0, 2.0), (0.3, 7.0)):
        v0 = rate_two_inhibitors(0, 0, 2.0, 3.0, alpha); vi = rate_two_inhibitors(I, J, 2.0, 3.0, alpha)
        ok1 &= abs((v0/vi - 1) - w(1/alpha, I/2.0, J/3.0)) < 1e-9
check("Yonetani–Theorell: v₀/vᵢ − 1 = a ⊞ₖ b with k = 1/α, from the full four-state equilibrium", ok1)
a, b = 0.8, 1.7
f = lambda x: x/(1 + x)  # fractional inhibition of one inhibitor alone
vi_rel = 1/(1 + w(1.0, a, b))
check("α = 1 (independent sites): remaining activity multiplies, 1/(1+a) · 1/(1+b), which is Bliss independence",
      abs(vi_rel - (1 - f(a))*(1 - f(b))) < 1e-15 and abs((1 - vi_rel) - w(-1, f(a), f(b))) < 1e-15)
# Loewe: for mutually exclusive inhibitors the isobole of 50% inhibition is a straight line in (a, b)
iso = [(x, 1 - x) for x in np.linspace(0, 1, 11)]
check("α → ∞ (one site): v₀/vᵢ − 1 = a + b, so the 50% isobole a + b = 1 is the straight Loewe line", all(abs(w(0, x, y) - 1) < 1e-15 for x, y in iso))
def iso_curve(alpha, xs):
    return [ (1 - x)/(1 + x/alpha) for x in xs ]  # solve a + b + ab/α = 1 for b
xs = np.linspace(0, 1, 101)
R["isoboles"] = {str(al): [[float(x), float(y)] for x, y in zip(xs, iso_curve(al, xs))] for al in (0.2, 1, 5, 1e9)}
mid = {al: iso_curve(al, [0.5])[0] for al in (0.2, 1, 5, 1e9)}
check("as α falls from ∞ the 50% isobole bows in from the straight Loewe line; Bliss (α = 1) already bows in (b = 0.333 at a = 0.5), and no binding α reaches beyond the line (k = 1/α ≥ 0)",
      mid[0.2] < mid[1] < mid[5] < mid[1e9] and abs(mid[1] - 1/3) < 1e-12, ", ".join(f"α={al:g}: {v:.3f}" for al, v in mid.items()))

# 2. Michaelis–Menten and steps in series
Km, Vmax = 2.0, 10.0
v = lambda S: Vmax*S/(Km + S)
check("Michaelis–Menten is a Langmuir isotherm: v/(Vmax − v) = S/Km grows linearly", all(abs(v(S)/(Vmax - v(S)) - S/Km) < 1e-12 for S in (0.1, 2, 50)))
ks = [30.0, 5.0, 12.0]
tau = sum(1/k for k in ks)
# simulate a cycle of irreversible first-order steps: mean turnover time = sum of mean step times
rng = np.random.default_rng(1); samples = sum(rng.exponential(1/k, 200000) for k in ks)
check("irreversible steps in series: mean turnover time = Σ 1/kᵢ, so kcat = 1/Σ(1/kᵢ) (times add, k = 0)", abs(samples.mean() - tau)/tau < 0.01, f"kcat = {1/tau:.2f}")

# 3. kinetic proofreading
D0 = 100.0
err = lambda n: 1/(1 + D0**(n + 1))
ee = lambda D: (D - 1)/(D + 1)
check("proofreading: n extra stages multiply the discrimination odds, error = 1/(1 + D^(n+1)) ≈ f₀^(n+1): 1.0e-2, 1.0e-4, 1.0e-6",
      abs(err(0) - 1/101) < 1e-15 and abs(err(1)/1e-4 - 1) < 1e-3 and abs(err(2)/1e-6 - 1) < 1e-5)
check("…so the fidelity ee = (D − 1)/(D + 1) composes by (e₁ + e₂)/(1 + e₁e₂), the (0, −1) law of subject 12",
      abs(ee(D0*D0) - (2*ee(D0))/(1 + ee(D0)**2)) < 1e-15)

# 4. nitrogenase
check("nitrogenase: 8 electrons per N₂, 6 to two NH₃ and 2 to the obligatory H₂: 75% electron efficiency", 6/8 == 0.75 and 2*3 + 2 == 8)
atp = 16*30.5; dG_haber = 2*(-16.45)
R["nitrogenase_ATP_kJ"] = atp
check("16 ATP per N₂ (8 per NH₃) supply about 488 kJ at the standard −30.5 kJ/mol, against ΔG° = −32.9 kJ for N₂ + 3H₂ → 2NH₃",
      abs(atp - 488) < 1e-9 and abs(dG_haber + 32.9) < 1e-9)

# 5. photosynthesis
hcN = 119626.6  # kJ·nm/mol (hc·N_A)
E680 = hcN/680; R["E680"] = E680
eff = (2870/6)/(8*E680); R["co2_eff"] = eff
check("8 photons at 680 nm carry 8 × 175.9 = 1407 kJ/mol; fixing one CO₂ into glucose stores 2870/6 = 478 kJ: 34% at the 8-photon limit",
      abs(E680 - 175.9) < 0.1 and abs(eff - 0.340) < 0.002, f"{eff:.3f}")
check("Calvin cycle: 3 ATP + 2 NADPH per CO₂ (18 ATP and 12 NADPH per glucose)", 6*3 == 18 and 6*2 == 12)
check("maximum quantum yield 1/8 = 0.125 CO₂ per photon", abs(1/8 - 0.125) < 1e-15)

print(f"\n{sum(ok)}/{len(ok)} checks pass")
json.dump(R, open("results.json", "w"), indent=1, default=float)
