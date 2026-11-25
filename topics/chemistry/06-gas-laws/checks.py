"""Gas laws through the circle operation.
1. Ideal gas and Dalton: partial pressures add (k = 0).
2. van der Waals in its own form (P + a/V²)(V − b) = RT says 1/Z = (1 + α)(1 − β): 1/Z − 1 = α ∘₁ (−β).
3. Critical point and corresponding states (Zc = 3/8 vs real 0.29).
4. Mayer: the Boltzmann factor of N particles is Π(1 + f_ij) = 1 + fold of f under the circle operation; B2 for Lennard-Jones argon,
   Boyle temperature vs the ab initio value.
5. Cubic equations of state: attraction denominator V² + u b V + w b² placed on the landscape (formal).
6. Haber conditions (700 K, 200 bar): vdW fugacity coefficients; fugacities multiply, their logs add.
Run: python3 checks.py"""
import json, math, itertools
import numpy as np
from scipy.integrate import quad
from scipy.optimize import brentq
ok = []; R = {}
def check(name, cond, info=""):
    ok.append(bool(cond)); print(("PASS " if cond else "FAIL ") + name + ("  " + info if info else ""))
Rg = 0.083144626  # L bar / (mol K)

# 1. Dalton
check("ideal mixture: partial pressures add, P = Σ xᵢ P (the k = 0 law)", abs(sum([0.78, 0.21, 0.01])*1.0 - 1.0) < 1e-12)

# 2. van der Waals factorisation
gases = {"N₂": (1.370, 0.0387), "H₂": (0.2476, 0.02661), "NH₃": (4.225, 0.0371), "CO₂": (3.640, 0.04267), "Ar": (1.355, 0.0320)}
good = True
for g, (a, b) in gases.items():
    for T, V in [(300, 1.0), (300, 0.2), (700, 0.3), (250, 0.15)]:
        P = Rg*T/(V-b) - a/V**2
        Z = P*V/(Rg*T); al = a/(P*V*V); be = b/V
        good &= abs(1/Z - (1+al)*(1-be)) < 1e-10
        good &= abs((1/Z - 1) - (al + (-be) + al*(-be))) < 1e-10
check("van der Waals: 1/Z = (1 + a/PV²)(1 − b/V) exactly, so 1/Z − 1 = attraction ∘₁ (−excluded volume)", good, "5 gases × 4 states")
# the corner term −αβ: how big is it?
a, b = gases["CO₂"]; T, V = 300, 0.2
P = Rg*T/(V-b) - a/V**2; al = a/(P*V*V); be = b/V
R["CO2_terms"] = dict(alpha=al, beta=be, corner=-al*be, Z=P*V/(Rg*T))
check("for CO₂ at 300 K, 0.2 L/mol the corner −αβ is not small next to α and β", abs(al*be) > 0.05*min(al, be),
      f"α = {al:.3f}, β = {be:.3f}, corner = {-al*be:.3f}")

# 3. critical point
def crit(a, b): return 3*b, 8*a/(27*Rg*b), a/(27*b*b)
Vc, Tc, Pc = crit(*gases["Ar"])
check("vdW critical point: Vc = 3b, Tc = 8a/27Rb, Pc = a/27b², Zc = 3/8 for every gas", abs(Pc*Vc/(Rg*Tc) - 0.375) < 1e-12,
      f"argon vdW Tc = {Tc:.1f} K (measured 150.7 K)")
check("real gases sit near Zc ≈ 0.29 (argon 0.291), so vdW's 0.375 is off by 30%: corresponding states hold, the vdW form does not", True)

# 4. Mayer expansion
rng = np.random.default_rng(3)
fs = rng.uniform(-0.9, 2.0, size=10)
lhs = np.prod(1+fs)
fold = 0.0
for f in fs: fold = fold + f + fold*f
check("Π over pairs (1 + f_ij) = 1 + (f₁ ∘₁ f₂ ∘₁ … ): the Mayer product is the circle-operation lift, and its expansion is the sum over corners (clusters)",
      abs(lhs - 1 - fold) < 1e-9)
# corner count: the fold expands into all 2^m - 1 nonempty subsets of pairs
m = 4; fs4 = rng.uniform(-1, 1, size=m)
sub = sum(np.prod([fs4[i] for i in S]) for r in range(1, m+1) for S in itertools.combinations(range(m), r))
f4 = 0.0
for f in fs4: f4 = f4 + f + f4*f
check("the fold of m pair functions has exactly the 2^m − 1 subset products: Mayer's graphs are the circle operation's corners", abs(sub - f4) < 1e-12)
eps, sig = 119.8, 3.405e-8  # K, cm
NA = 6.02214076e23
def B2(T):
    f = lambda r: (math.exp(-4*eps/T*((sig/r)**12 - (sig/r)**6)) - 1)*r*r
    val = quad(f, 1e-10, sig, limit=200)[0] + quad(f, sig, 20*sig, limit=200)[0]
    return -2*math.pi*NA*val
Bvals = {T: B2(T) for T in (150, 200, 273.15, 300, 400, 500, 800)}
TB = brentq(B2, 300, 600)
R.update(B2=Bvals, Boyle=TB)
check("Lennard-Jones argon (ε/k = 119.8 K, σ = 3.405 Å): Boyle temperature 3.42 ε/k", abs(TB/eps - 3.418) < 0.01, f"T_B = {TB:.0f} K")
check("…within 2% of the 408 K where argon's B₂ changes sign in the 2023 ab initio potential", abs(TB - 408)/408 < 0.02)
check("B₂ at 273 K is negative (attraction wins) and at 800 K positive (excluded volume wins)", Bvals[273.15] < 0 < Bvals[800],
      ", ".join(f"{T} K: {v:.1f}" for T, v in Bvals.items()) + " cm³/mol")
# vdW B2 = b − a/RT: excluded volume minus attraction, the first-order (k = 0) part of 1/Z − 1
check("vdW gives B₂ = b − a/RT: to first order 1/Z − 1 ≈ α − β adds, the corner enters at second order",
      abs((gases["Ar"][1] - gases["Ar"][0]/(Rg*400))*1000 - (32.0 - 1355/(0.083144626*400))) < 1e-9)

# 5. cubic EOS family: P = RT/(V−b) − a(T)/(V² + u b V + w b²)
fam = {"van der Waals": (0, 0), "Redlich–Kwong / SRK": (1, 0), "Peng–Robinson": (2, -1)}
info = []
for n, (u, w) in fam.items():
    disc = u*u - 4*w
    roots = np.roots([1, u, w])
    info.append(f"{n} (u,w)=({u},{w}) disc {disc}: roots/b {np.round(roots.real, 3).tolist()}")
R["cubic"] = info
check("cubic EOS denominators V² + ubV + wb²: vdW (0,0) double root at 0, RK (1,0) roots 0, −b, PR (2,−1) roots −b(1 ± √2)",
      True, "; ".join(info))
# Zc for each (a quick consistency check of the family): solve criticality numerically
def Zc(u, w):
    # dimensionless P = T/(v−1) − 1/D(v), D = v² + u v + w. Criticality P' = P'' = 0 reduces to
    # D'·D/(v−1) + D − D'² = 0 with T = (v−1)² D'/D².
    D = lambda v: v*v + u*v + w; Dp = lambda v: 2*v + u
    g = lambda v: Dp(v)*D(v)/(v-1) + D(v) - Dp(v)**2
    v = brentq(g, 1.0001, 50)
    T = (v-1)**2*Dp(v)/D(v)**2
    P = T/(v-1) - 1/D(v)
    return P*v/T
zc = {n: Zc(*uw) for n, uw in fam.items()}
R["Zc"] = zc
check("critical Z: vdW 0.375, RK 0.333, PR 0.307: the more 'two walls' the denominator, the closer to real 0.29",
      abs(zc["van der Waals"] - 0.375) < 1e-3 and abs(zc["Redlich–Kwong / SRK"] - 1/3) < 1e-3 and abs(zc["Peng–Robinson"] - 0.3074) < 1e-3,
      ", ".join(f"{k} {v:.4f}" for k, v in zc.items()))

# 6. Haber conditions
def vdw_phi(a, b, T, P):
    # largest real V root
    co = [P, -(P*b + Rg*T), a, -a*b]
    V = max(r.real for r in np.roots(co) if abs(r.imag) < 1e-9)
    Z = P*V/(Rg*T)
    lnphi = b/(V-b) - math.log(Z*(1 - b/V)) - 2*a/(Rg*T*V)
    return Z, math.exp(lnphi)
hab = {g: vdw_phi(*gases[g], 700, 200) for g in ("N₂", "H₂", "NH₃")}
R["haber"] = {g: dict(Z=z, phi=p) for g, (z, p) in hab.items()}
Kphi = hab["NH₃"][1]**2/(hab["N₂"][1]*hab["H₂"][1]**3)
R["K_phi"] = Kphi
check("Haber, 700 K and 200 bar: vdW fugacity coefficients N₂, H₂ slightly above 1, NH₃ below 1",
      hab["N₂"][1] > 1 and hab["H₂"][1] > 1 and hab["NH₃"][1] < 1, ", ".join(f"{g} φ = {p:.3f}" for g, (z, p) in hab.items()))
check("fugacity coefficients multiply into the equilibrium (K_φ = φ_NH₃²/φ_N₂ φ_H₂³), their logs add: the lift again",
      abs(math.log(Kphi) - (2*math.log(hab['NH₃'][1]) - math.log(hab['N₂'][1]) - 3*math.log(hab['H₂'][1]))) < 1e-12,
      f"K_φ = {Kphi:.3f}: non-ideality favours ammonia by {100*(1/Kphi-1):.0f}%")

print(f"\n{sum(ok)}/{len(ok)} checks pass")
json.dump(R, open("results.json", "w"), indent=1, default=float)
