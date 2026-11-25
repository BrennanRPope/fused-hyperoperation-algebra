"""Electrochemistry and Nernst through the circle operation, with photosynthesis.
1. Nernst and coupling: ΔG = −nFE adds over coupled steps, K multiplies, so K − 1 composes by the circle operation (as in subject 7).
   Potentials themselves combine as an electron-weighted mean (Luther's rule), the slope rule of a Frost diagram.
   A middle state above the chord disproportionates. Data: standard potentials for copper and iron (CRC-style tables).
2. Flux and force: the forward/backward ratio of any step is e^(A/RT); in series the ratios multiply, so (ratio − 1)
   composes by the circle operation. Butler–Volmer: only the fully one-sided case (α = 1) makes the current itself compose by the circle operation in η.
3. Two redox sites: the electron binding polynomial minus one is x ⊞ₖ x. ΔE = E₁ − E₂ = (RT/F) ln(4/k):
   independent sites sit 35.6 mV apart, cooperativity (k > 1) means ΔE < 35.6 mV, inversion (E₂ > E₁) means k > 4.
   This is subject 8's diacid with electrons for protons.
4. Photosynthesis: stage losses compose at k = −1; Z-scheme energetics; Kok's four-flash clock as a circulant.
Run: python3 checks.py"""
import json, math
import numpy as np
ok = []; R = {}
def check(name, cond, info=""):
    ok.append(bool(cond)); print(("PASS " if cond else "FAIL ") + name + ("  " + info if info else ""))
def w(k, a, b): return a + b + k*a*b
F, Rg, T = 96485.33212, 8.314462618, 298.15
f = F/(Rg*T)

# 1. Nernst, Luther, Frost
dG = lambda n, E: -n*F*E
E1, E2 = 0.40, -0.15; n1, n2 = 2, 1
K = lambda n, E: math.exp(n*f*E)
check("coupled steps: ΔG adds, K multiplies, so K − 1 composes by the circle operation", abs((K(n1,E1)*K(n2,E2) - 1) - w(1, K(n1,E1)-1, K(n2,E2)-1)) < 1e-6*K(n1,E1)*K(n2,E2))
luther = lambda n1, E1, n2, E2: (n1*E1 + n2*E2)/(n1 + n2)
cu = luther(1, 0.153, 1, 0.521); fe = luther(1, 0.771, 2, -0.447)
R["cu2_cu"] = cu; R["fe3_fe"] = fe
check("Luther's rule: Cu²⁺/Cu = mean of Cu²⁺/Cu⁺ (0.153) and Cu⁺/Cu (0.521) = 0.337 V; tables give 0.342", abs(cu - 0.337) < 1e-9 and abs(cu - 0.3419) < 0.006)
check("Fe³⁺/Fe = (0.771 + 2·(−0.447))/3 = −0.041 V; tables give −0.037", abs(fe + 0.041) < 1e-3 and abs(fe + 0.037) < 0.005)
# Frost: points (N, nE) relative to the element; slope = potential
frost_cu = [(0, 0.0), (1, 0.521), (2, 0.521 + 0.153)]
mid_above = frost_cu[1][1] > (frost_cu[0][1] + frost_cu[2][1])/2
Kdisp = math.exp(f*(0.521 - 0.153)); R["cu_disp_K"] = Kdisp
check("Frost convexity: Cu⁺ sits above the Cu–Cu²⁺ chord, so 2Cu⁺ → Cu + Cu²⁺ with K = e^(F·0.368/RT) ≈ 1.7e6", mid_above and abs(Kdisp/1.66e6 - 1) < 0.03, f"K = {Kdisp:.3g}")
# mean is not compose by the circle operation: weighted mean of potentials is order-free and idempotent
check("the Luther mean is idempotent (E with E gives E), so potentials do not compose by the circle operation; only K − 1 does", abs(luther(1, 0.3, 2, 0.3) - 0.3) < 1e-15)

# 2. Flux-force and Butler–Volmer
A1, A2 = 3.1, 1.7  # affinities in RT
r1, r2 = math.exp(A1), math.exp(A2)
check("De Donder: in series the forward/backward ratios multiply, so (ratio − 1) composes by the circle operation", abs((r1*r2 - 1) - w(1, r1-1, r2-1)) < 1e-9)
bv = lambda eta, aa: math.exp(aa*f*eta) - math.exp(-(1-aa)*f*eta)
e1, e2 = 0.031, 0.047
check("Butler–Volmer, α = 1 (one-sided): i(η₁+η₂)/i₀ = j₁ ∘₁ j₂", abs(bv(e1+e2, 1) - w(1, bv(e1, 1), bv(e2, 1))) < 1e-12)
gap = bv(e1+e2, 0.5) - w(1, bv(e1, 0.5), bv(e2, 0.5)); R["bv_gap_half"] = gap
check("Butler–Volmer, α = ½: the current does not compose by the circle operation (gap ≠ 0); only the forward/backward ratio e^(Fη/RT) does",
      abs(gap) > 0.05 and all(abs(math.exp(aa*f*0.05)/math.exp(-(1-aa)*f*0.05) - math.exp(f*0.05)) < 1e-9 for aa in (0.2, 0.5, 0.8)), f"gap = {gap:.3f}")
tafel = 1000*math.log(10)/(0.5*f); R["tafel_mV"] = tafel
check("Tafel slope at α = ½, 25 °C: 118 mV per decade", abs(tafel - 118.3) < 0.1, f"{tafel:.1f} mV")

# 3. two redox sites
dE_indep = 1000*math.log(4)/f; R["dE_indep_mV"] = dE_indep
check("independent identical sites (k = 1): E₁ − E₂ = (RT/F) ln 4 = 35.6 mV", abs(dE_indep - 35.6) < 0.05, f"{dE_indep:.2f} mV")
def k_of_dE(dE_mV): return 4*math.exp(-dE_mV/1000*f)
def Z(x, k): return 1 + 2*x + k*x*x
ok3 = True
for dE in (150, 100, 35.6, 0, -50):
    k = k_of_dE(dE); Em = 0.0
    # macroscopic: Z = 1 + e^{f(E1−E)} + e^{f(E1+E2−2E)}; microscopic: x = e^{f(Em−E)}
    E1m = Em + math.log(2)/f; E2m = 2*Em + math.log(k)/f - E1m
    for E in (-0.1, 0.0, 0.07):
        x = math.exp(f*(Em - E))
        ok3 &= abs(Z(x, k) - (1 + math.exp(f*(E1m - E)) + math.exp(f*(E1m + E2m - 2*E)))) < 1e-9*Z(x, k)
    ok3 &= abs(1000*(E1m - E2m) - dE) < 1e-6
check("Z − 1 = x ⊞ₖ x matches the two macroscopic potentials, with E₁ − E₂ = (RT/F) ln(4/k)", ok3)
check("comproportionation constant K_c = e^(F·ΔE/RT) = 4/k (diacid analogue: k = 4·10^(−ΔpKa))", abs(math.exp(0.1*f) - 4/k_of_dE(100)) < 1e-9)
nH = lambda k: 2*math.sqrt(k)/(1 + math.sqrt(k))
def nH_num(k, h=1e-6):
    x0 = 1/math.sqrt(k); th = lambda x: (2*x + 2*k*x*x)/(2*Z(x, k))
    y = lambda x: math.log(th(x)/(1 - th(x)))
    return (y(x0*(1+h)) - y(x0*(1-h)))/(math.log(x0*(1+h)) - math.log(x0*(1-h)))
check("Nernst-curve steepness at the midpoint: n = 2√k/(1 + √k), from 0 (k → 0) to 2 (k → ∞), 1 at 35.6 mV",
      all(abs(nH(k) - nH_num(k)) < 1e-5 for k in (0.08, 1, 4, 30)))
check("cooperative two-electron behaviour (n > 1) ⇔ k > 1 ⇔ ΔE < 35.6 mV ⇔ complex zeros of 1 + 2x + kx²",
      all((nH(k) > 1) == (k > 1) == (4 - 4*k < 0) for k in (0.2, 0.9, 1.1, 5)))
check("potential inversion (E₂ > E₁) ⇔ k > 4", (k_of_dE(-1) > 4) and (k_of_dE(1) < 4))

# 4. photosynthesis
hc = 1239.84193
p680, p700 = hc/680, hc/700
stored = 4*(0.815 + 0.320)  # eV: 4 electrons from water (+0.815 V) to NADP⁺ (−0.320 V)
eff = stored/(4*p680 + 4*p700); R["zscheme_eff"] = eff
check("Z-scheme: 8 photons per O₂ (4 at 680 nm, 4 at 700 nm, 14.4 eV) store 4 × 1.135 eV in NADPH: 31.6%", abs(eff - 0.316) < 0.002, f"{eff:.3f}")
one = 1.135/p680; R["one_photon_eff"] = one
check("thermodynamics alone would allow one 680 nm photon per electron (1.82 eV > 1.135 eV, 62%); the two-photon Z-scheme halves that, so its cost is in the scheme, not the energy",
      p680 > 1.135 and abs(one - 0.623) < 0.002 and abs(eff/one - 0.507) < 0.003)
par = {"C3": 4.6/9.4, "C4": 6.0/12.3}; R["par_fraction"] = par
check("solar-to-biomass ceilings 4.6% (C3) and 6.0% (C4) of full sunlight are 9.4% and 12.3% of PAR: both imply PAR ≈ 49% of sunlight",
      abs(par["C3"] - par["C4"]) < 0.002 and 0.48 < par["C3"] < 0.50)
losses = [0.51, 0.10, 0.20, 0.70, 0.30, 0.25]
tot = 0
for L in losses: tot = w(-1, tot, L)
check("any chain of stage losses composes at k = −1: total loss = L₁ ⊞₋₁ L₂ ⊞₋₁ … = 1 − Π(1 − Lᵢ)", abs(tot - (1 - np.prod([1-L for L in losses]))) < 1e-12)
# Kok clock
a, b = 0.10, 0.05
P = np.roll(np.eye(4), 1, axis=0)
Tm = a*np.eye(4) + (1-a-b)*P + b*P@P
s = np.array([0.25, 0.75, 0, 0]); y = []
for n in range(40):
    # O2 released when the clock wraps past S3: from S3 with ≥1 advance, from S2 with a double hit
    y.append(s[3]*(1 - a) + s[2]*b); s = Tm@s
y = np.array(y); R["kok_yields"] = y[:16].tolist()
check("Kok clock (25% S₀, 75% S₁ dark, 10% misses, 5% double hits): O₂ peaks on flash 3, then every 4 flashes",
      int(np.argmax(y[:6])) == 2 and int(np.argmax(y[4:10])) + 4 == 6)
lam = np.linalg.eigvals(Tm); lam1 = abs(a - b + 1j*(1 - a - b)); R["kok_damp"] = lam1
check("the damping per flash is |α − β + i(1 − α − β)| = 0.851, the circulant eigenvalue at z = i", abs(sorted(abs(lam))[-2] - lam1) < 1e-12 and abs(lam1 - 0.8515) < 1e-3)
check("long-run O₂ per flash = (1 − α + β)/4 (mean clock advance over 4)", abs(y[-1] - (1 - a + b)/4) < 1e-3, f"{y[-1]:.4f}")

print(f"\n{sum(ok)}/{len(ok)} checks pass")
json.dump(R, open("results.json", "w"), indent=1, default=float)
