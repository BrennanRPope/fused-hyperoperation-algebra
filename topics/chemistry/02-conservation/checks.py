"""Conservation of mass and charge through mult-add and the landscape.
Invariant mass is a law of cosines: M^2 = m1^2 + m2^2 + 2 m1 m2 c.
 c = 1: Lavoisier (masses add exactly; the border, a flat triangle)
 c < 1: bound states (turning country, real angle theta)
 c > 1: free pairs in motion (two-walls country, c = cosh of relative rapidity)
Run: python3 checks.py (numpy, mpmath)."""
import numpy as np, mpmath as mp, json
mp.mp.dps = 40
ok = []
def check(name, cond, info=""):
    ok.append(bool(cond)); print(("PASS " if cond else "FAIL ") + name + ("  " + info if info else ""))
res = {}
# CODATA 2018 (MeV) and measured binding energies
me = mp.mpf("0.51099895000"); mp_ = mp.mpf("938.27208816"); mn = mp.mpf("939.56542052")
alpha = mp.mpf("7.2973525693e-3")
mH = mp_ + me - mp.mpf("13.598434005e-6")
def angle(m1, m2, B):
    M = m1 + m2 - B
    c = (M**2 - m1**2 - m2**2)/(2*m1*m2)
    return mp.acos(c), c
# 1. free pairs: random four-vectors, c = (E1E2 - p1.p2)/(m1 m2) = gamma_rel >= 1
rng = np.random.default_rng(3)
good = True
for _ in range(500):
    m1, m2 = rng.uniform(0.1, 5, 2); p1, p2 = rng.normal(size=(2,3))*3
    E1, E2 = np.sqrt(m1**2 + p1@p1), np.sqrt(m2**2 + p2@p2)
    M2 = (E1+E2)**2 - (p1+p2)@(p1+p2)
    c = (E1*E2 - p1@p2)/(m1*m2)
    good &= abs(M2 - (m1**2 + m2**2 + 2*m1*m2*c)) < 1e-9 and c >= 1 - 1e-12
check("free pairs: M² = m1² + m2² + 2m1m2·c with c = γ_rel ≥ 1 (500 random pairs): two-walls country", good)

# 2. hydrogen 1s: the binding angle is the fine-structure constant
th, c = angle(mp_, me, mp.mpf("13.598434005e-6"))
res["hydrogen_theta"] = float(th); res["alpha"] = float(alpha)
check("hydrogen 1s: binding angle θ ≈ α (within 1e-4 relative)", abs(th/alpha - 1) < 1e-4, f"θ = {float(th):.9e}, α = {float(alpha):.9e}, sin θ/α - 1 = {float(mp.sin(th)/alpha-1):.2e}")
# compact two-body Dirac form: cos θ = sqrt(1 - α²) exactly (sin θ = α)
Bpred = mp_ + me - mp.sqrt(mp_**2 + me**2 + 2*mp_*me*mp.sqrt(1-alpha**2))
res["B_pred_eV"] = float(Bpred*1e6)
check("sin θ = α predicts hydrogen's ionisation energy to 3e-6 relative (rest = Lamb shift, nuclear size)",
      abs(Bpred*1e6/mp.mpf("13.598434005") - 1) < 3e-6, f"predicted {float(Bpred*1e6):.9f} eV vs measured 13.598434005 eV")
resid = mp.mpf("13.598434005") - Bpred*1e6   # eV
lamb1s = mp.mpf("8172.874e6")*mp.mpf("4.135667696e-15")  # 1s Lamb shift, 8172.874 MHz
check("what's left over (3.4e-5 eV) is the 1s Lamb shift, a QED effect, to within 2%", abs(-resid/lamb1s - 1) < 0.02, f"left over {float(-resid):.3e} eV, Lamb shift {float(lamb1s):.3e} eV")
# positronium: same angle, very different masses
th_ps, _ = angle(me, me, mp.mpf("6.8029e-6"))
check("positronium: same angle θ ≈ α though the masses are equal", abs(th_ps/alpha - 1) < 2e-4, f"θ/α = {float(th_ps/alpha):.6f}")
# muonic hydrogen-like check with muonium (mu+ e-): binding 13.539 eV
mmu = mp.mpf("105.6583755")
th_mu, _ = angle(mmu, me, mp.mpf("13.5391e-6"))
check("muonium (μ⁺e⁻): θ ≈ α", abs(th_mu/alpha - 1) < 2e-4, f"θ/α = {float(th_mu/alpha):.6f}")
# He+ (Z=2): angle 2α
mHe4 = mp.mpf("3727.3794066")  # alpha particle
th_he, _ = angle(mHe4, me, mp.mpf("54.4177634e-6"))
check("He⁺ (Z = 2): θ ≈ 2α (Coulomb angle is Zα)", abs(th_he/(2*alpha) - 1) < 3e-4, f"θ/(2α) = {float(th_he/(2*alpha)):.6f}")

# 3. chemical bonds and nuclei on the same ruler
table = []
cases = [
 ("H₂ bond (H + H)", mH, mH, mp.mpf("4.4781e-6")),
 ("HCl bond", mH, mp.mpf("35.453")*mp.mpf("931.49410242") - 17*me + 17*me, mp.mpf("4.434e-6")),
 ("Hydrogen atom (p + e)", mp_, me, mp.mpf("13.598434005e-6")),
 ("Deuteron (p + n)", mp_, mn, mp.mpf("2.224566")),
]
for name, m1, m2, B in cases:
    th, c = angle(m1, m2, B)
    mu = m1*m2/(m1+m2)
    approx = mp.sqrt(2*B/mu)
    table.append((name, float(th), float(B/(m1+m2)), float(approx)))
    check(f"{name}: θ = √(2B/μ) to leading order", abs(th/approx - 1) < 0.01, f"θ = {float(th):.3e} rad, mass defect {float(B/(m1+m2)):.2e} of the total")
res["table"] = table
th_bond = table[0][1]
check("chemistry sits ~50× closer to the border than the hydrogen atom, nuclei ~10× farther", table[0][1] < table[2][1]/40 and table[3][1] > 10*table[2][1])

# 4. Lavoisier's balance could never see it: 2H2 + O2 -> 2H2O(l), ΔH = -571.66 kJ
dm = 571.66e3/ (299792458.0**2)  # kg per 2 mol reaction
rel = dm/(2*18.01528e-3)
res["water_mass_defect_rel"] = rel
check("burning hydrogen: mass lost is ~1.8e-10 of the water made (Lavoisier's balances: ~1e-4 at best)", 1.5e-10 < rel < 2e-10, f"{rel:.2e}")

# 5. charge: Q adds, the phase e^{iQφ} multiplies, so its excess mult-adds (rotation)
good = True
for _ in range(300):
    Q1, Q2 = rng.integers(-5, 6, 2); phi = rng.uniform(0, 2*np.pi)
    w1, w2 = np.exp(1j*Q1*phi)-1, np.exp(1j*Q2*phi)-1
    good &= abs((w1 + w2 + w1*w2) - (np.exp(1j*(Q1+Q2)*phi)-1)) < 1e-12
check("charge: excess phase e^{iQφ}−1 mult-adds as charges add (300 cases); |1+w| = 1 is the turning country's circle", good)
# charge and atom balance of a redox equation: MnO4- + 5Fe2+ + 8H+ -> Mn2+ + 5Fe3+ + 4H2O
lhs = dict(Mn=1, O=4, Fe=5, H=8, q=-1+10+8); rhs = dict(Mn=1, Fe=5, O=4, H=8, q=2+15)
check("permanganate–iron(II): atoms and charge both balance (charge +17 each side)", lhs == rhs)
json.dump(res, open("results.json","w"), indent=1)
print(f"{sum(ok)}/{len(ok)} checks pass")
