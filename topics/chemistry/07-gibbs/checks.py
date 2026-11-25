"""Chemical thermodynamics (Gibbs) through the circle operation.
1. Free energies add when reactions add; equilibrium constants multiply (the lift); K − 1 composes by the circle operation.
2. ATP coupling (glutamine synthetase step) as a product of K's.
3. Partition functions: independent parts multiply. Two sites with interaction Boltzmann factor ω:
   Z = 1 + x + y + ω x y, so Z − 1 = x ⊞_ω y. The dial k is e^{−ε/kT}: k = 0 exclusion (addition),
   k = 1 independence (the circle operation), k > 1 attraction. Correlation sign = sign(k − 1).
4. Nitrogen fixation thermodynamics: ΔH, ΔS, ΔG(T) and K(T) for N₂ + 3H₂ → 2NH₃.
Run: python3 checks.py"""
import json, math
import numpy as np
ok = []; R = {}
def check(name, cond, info=""):
    ok.append(bool(cond)); print(("PASS " if cond else "FAIL ") + name + ("  " + info if info else ""))
Rg = 8.314462618; T0 = 298.15
def w(k, a, b): return a + b + k*a*b

# 1. Hess and the lift
rng = np.random.default_rng(7)
good = True
for _ in range(100):
    g1, g2 = rng.normal(0, 20000, 2)
    K1, K2 = math.exp(-g1/(Rg*T0)), math.exp(-g2/(Rg*T0))
    Kt = math.exp(-(g1+g2)/(Rg*T0))
    good &= abs(Kt/(K1*K2) - 1) < 1e-9 and abs(((K1-1) + (K2-1) + (K1-1)*(K2-1)) - (Kt-1)) < 1e-9*max(1, Kt)
check("reactions add ⇒ ΔG adds (k = 0) and K multiplies; K − 1 composes by the circle operation: (K₁−1) ∘₁ (K₂−1) = K₁K₂ − 1", good)

# 2. ATP coupling (Lehninger textbook values, kJ/mol at pH 7)
dG_gln, dG_atp = 14.2e3, -30.5e3
Kg, Ka = math.exp(-dG_gln/(Rg*T0)), math.exp(-dG_atp/(Rg*T0))
R.update(K_gln=Kg, K_atp=Ka, K_coupled=Kg*Ka)
check("glutamine from glutamate + NH₄⁺ (+14.2 kJ/mol, K′ = 3.2e−3) coupled to ATP hydrolysis (−30.5, K′ = 2.2e5): K′ = 7.2e2",
      abs(Kg*Ka - math.exp(16.3e3/(Rg*T0))) < 1e-6*Kg*Ka, f"{Kg:.2e} × {Ka:.2e} = {Kg*Ka:.2e}")

# 3. two interacting sites
def Z2(x, y, om): return 1 + x + y + om*x*y
good = all(abs((Z2(x, y, om) - 1) - w(om, x, y)) < 1e-12 for x, y, om in rng.uniform(0, 3, (50, 3)))
check("two sites with interaction factor ω: Z − 1 = x + y + ωxy = x ⊞_ω y exactly", good)
check("ω = 0 (hard exclusion) gives Z = 1 + x + y: plain addition; ω = 1 (no interaction) gives (1+x)(1+y): the circle operation",
      Z2(0.7, 1.3, 0) == 1 + 0.7 + 1.3 and abs(Z2(0.7, 1.3, 1) - 1.7*2.3) < 1e-12)
def corr(x, y, om):
    Z = Z2(x, y, om); n1 = (x + om*x*y)/Z; n2 = (y + om*x*y)/Z; n12 = om*x*y/Z
    return n12 - n1*n2
good = all(np.sign(round(corr(x, y, om), 14)) == np.sign(round(om - 1, 14)) for x, y in rng.uniform(0.1, 3, (30, 2)) for om in (0, 0.3, 1, 2.5, 10))
check("occupancy correlation ⟨n₁n₂⟩ − ⟨n₁⟩⟨n₂⟩ has the sign of k − 1: the dial reads attraction or repulsion", good)
# exact formula: corr = xy(ω − 1)/Z²
good = all(abs(corr(x, y, om) - x*y*(om-1)/Z2(x, y, om)**2) < 1e-12 for x, y, om in rng.uniform(0, 3, (50, 3)))
check("in fact ⟨n₁n₂⟩ − ⟨n₁⟩⟨n₂⟩ = xy(k − 1)/Z²: the covariance is the dial's distance from the circle operation", good)
# k as Boltzmann factor
eps = {"repulsive 5 kJ/mol": 5e3, "none": 0.0, "attractive 5 kJ/mol": -5e3}
R["k_of_eps"] = {n: math.exp(-e/(Rg*T0)) for n, e in eps.items()}
check("k = e^(−ε/RT): a 5 kJ/mol repulsion sets k = 0.13, an equal attraction k = 7.5 (298 K)",
      abs(R["k_of_eps"]["repulsive 5 kJ/mol"] - 0.133) < 0.001 and abs(R["k_of_eps"]["attractive 5 kJ/mol"] - 7.52) < 0.01)

# 4. Haber thermodynamics (NIST/CRC standard values)
dfH = {"NH3": -45.90e3, "N2": 0.0, "H2": 0.0}
S = {"NH3": 192.77, "N2": 191.61, "H2": 130.68}
dH = 2*dfH["NH3"]; dS = 2*S["NH3"] - S["N2"] - 3*S["H2"]
dG298 = dH - T0*dS
R.update(dH=dH, dS=dS, dG298=dG298)
check("N₂ + 3H₂ → 2NH₃: ΔH° = −91.8 kJ, ΔS° = −198.1 J/K, ΔG°(298) = −32.7 kJ = 2 × CRC ΔfG°(NH₃) = 2 × −16.4",
      abs(dH + 91.8e3) < 1e3 and abs(dS + 198.1) < 0.1 and abs(dG298/2 + 16.4e3) < 0.2e3, f"ΔG° = {dG298/1000:.2f} kJ")
Tstar = dH/dS
R["T_K1"] = Tstar
check("K = 1 near 463 K (ΔH/ΔS, constant-ΔH approximation): above it the equilibrium favours N₂ + H₂", abs(Tstar - 463) < 2, f"{Tstar:.0f} K")
Ks = {T: math.exp(-(dH - T*dS)/(Rg*T)) for T in (298.15, 400, 500, 600, 700, 800)}
R["K"] = Ks
check("K(298) ≈ 5.5e5 but K(700) ≈ 3e−4 (per bar², constant ΔH): exothermic, so heat lowers K (van 't Hoff)",
      Ks[298.15] > 1e5 and Ks[700] < 1e-3, ", ".join(f"{T:.0f} K: {v:.2e}" for T, v in Ks.items()))
# van 't Hoff: lnK linear in 1/T with slope −ΔH/R
Ts = np.array(list(Ks)); lnK = np.log(list(Ks.values()))
slope = np.polyfit(1/Ts, lnK, 1)[0]
check("van 't Hoff: ln K against 1/T has slope −ΔH/R (log side adds; K side multiplies by e^(−ΔH/R·Δ(1/T)))", abs(slope + dH/Rg) < 1e-6*abs(dH/Rg))

print(f"\n{sum(ok)}/{len(ok)} checks pass")
json.dump(R, open("results.json", "w"), indent=1, default=float)
