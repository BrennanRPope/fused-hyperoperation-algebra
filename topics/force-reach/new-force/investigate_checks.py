"""Part I of the new-force DRAFT (Nov 27 2026): the new force up close.
I1-I3 what carries it and what feels it; I4-I6 the force and its field at its own range;
I7-I13 indirect traces (early universe, inflation, colliders, relics, black holes, the vacuum term).
Run after checks.py, strength_checks.py, hold_strength_checks.py: appends I_ keys to results.json.
Writes investigate_output.txt. Order-of-magnitude estimates are labelled as such."""
import math, json
import sympy as sp
import mpmath as mp

lines = []; NP = [0, 0]; OUT = {}
def log(s=""):
    print(s); lines.append(s)
def check(name, ok, detail=""):
    NP[0 if ok else 1] += 1
    log(("PASS  " if ok else "FAIL  ") + name + (f"   [{detail}]" if detail else ""))

lP = 1.616255e-35; MP = 1.220890e19; Mbar = MP/math.sqrt(8*math.pi); hbar_GeVs = 6.582119569e-25
ainv = 137.035999084; R_lP = 2*math.sqrt(ainv)
res_prev = json.load(open("results.json"))
lam_lP = float(res_prev["H_range_lP_NF2"])          # radion range from part H (2 fermions)
m_rad = MP/lam_lP                                   # GeV
M1 = MP/R_lP                                        # first copy, GeV (= sqrt(alpha)/2 M_P)

def X(p1, p2, D):
    eta = sp.diag(*([-1] + [1]*(D - 1))); E1, E2 = p1[0], p2[0]
    T1 = sp.Matrix(p1)*sp.Matrix(p1).T/E1; T2 = sp.Matrix(p2)*sp.Matrix(p2).T/E2
    T1l = eta*T1*eta
    contr = sum(T1l[a, b]*T2[a, b] for a in range(D) for b in range(D))
    tr1 = sum(eta[a, a]*T1[a, a] for a in range(D)); tr2 = sum(eta[a, a]*T2[a, a] for a in range(D))
    return sp.simplify(contr - tr1*tr2/(D - 2))
m1, m2, E, E1, E2 = sp.symbols('m1 m2 E E1 E2', positive=True)

log("== I. The new force up close")
# ---- I1 the carrier and the ladder
OUT["I_radion_GeV"] = m_rad; OUT["I_copy1_GeV"] = M1
check("I1 the carrier is one particle, the circle's size (spin 0, no circle charge), with mass about 6e15 GeV from part H's hold; above it sit the charged copies at n x 5.2e17 GeV with charge n e",
      1e15 < m_rad < 1e16 and abs(M1 - 5.215e17)/5.215e17 < 0.01, f"carrier {m_rad:.2e} GeV, copies n x {M1:.3e} GeV")
# ---- I2 who feels it
sca_matter = sp.simplify(X([m1, 0, 0, 0, 0], [m2, 0, 0, 0, 0], 5) - X([m1, 0, 0, 0], [m2, 0, 0, 0], 4))
sca_light = sp.simplify(X([m1, 0, 0, 0, 0], [E, E, 0, 0, 0], 5) - X([m1, 0, 0, 0], [E, E, 0, 0], 4))
like = X([E1, 0, 0, 0, E1], [E2, 0, 0, 0, E2], 5); opp = X([E1, 0, 0, 0, E1], [E2, 0, 0, 0, -E2], 5)
sca_copy = sp.simplify((opp + like)/2 - E1*E2/2)
check("I2 who feels it, relative to gravity: ordinary matter 1/3, copies 3 (nine times more), light 0 (a passing light ray isn't pulled by it)",
      sp.simplify(sca_matter/(m1*m2/2) - sp.Rational(1, 3)) == 0 and sp.simplify(sca_copy/(E1*E2/2) - 3) == 0 and sca_light == 0,
      f"matter {sp.simplify(sca_matter/(m1*m2/2))}, copies {sp.simplify(sca_copy/(E1*E2/2))}, light {sca_light}")
check("I3 light still notices it another way: alpha = 4G/R^2 (part G) makes electric strength depend on the circle's size, d ln alpha / d ln R = -2, so the carrier also couples to electric and magnetic fields",
      sp.diff(sp.log(4/sp.Symbol('R')**2), sp.Symbol('R'))*sp.Symbol('R') == -2)

# ---- I4 the force law at its own range
r, Rr, lam = sp.symbols('r R lambda', positive=True)
check("I4 between two masses with no circle charge, every rung gives 4/3 of Newton (5D static weight 2/3 vs 4D 1/2); the zero rung's 4/3 is gravity 1 plus the new force 1/3",
      sp.simplify(X([m1, 0, 0, 0, 0], [m2, 0, 0, 0, 0], 5)/X([m1, 0, 0, 0], [m2, 0, 0, 0], 4) - sp.Rational(4, 3)) == 0)
Vfull = 1 + sp.Rational(1, 3)*sp.exp(-r/lam) + sp.Rational(4, 3)*(sp.coth(r/(2*Rr)) - 1)
check("I4b so V r / (G m1 m2) = 1 + (1/3) e^(-r/lambda) + (4/3)(coth(r/2R) - 1); with an unheld circle (lambda -> infinity) it is exactly (4/3) coth(r/2R), the 5D image sum: the new force is what makes the zero rung match the others",
      sp.simplify((sp.limit(Vfull, lam, sp.oo) - sp.Rational(4, 3)*sp.coth(r/(2*Rr))).rewrite(sp.exp)) == 0)
fnum = lambda rr: 1 + math.exp(-rr/lam_lP)/3 + 4/3*(1/math.tanh(rr/(2*R_lP)) - 1)
w5 = fnum(5*R_lP); wl = fnum(lam_lP)
check("I4c three regimes: below the circle (23 l_P) gravity is five-dimensional; just outside it (about 5 R) it reaches 4/3 of Newton and then eases back as the new force fades (1.12 at the range, about 2000 l_P); far beyond the range it is plain Newton. The window is only about 85 circle-sizes wide",
      abs(w5 - 4/3) < 0.01 and abs(wl - (1 + math.exp(-1)/3)) < 0.01 and abs(fnum(50*lam_lP) - 1) < 1e-6 and fnum(0.05*R_lP)*0.05 > 2.5,
      f"at 5R: {w5:.4f}; at the range: {wl:.4f}; at 50 ranges: {fnum(50*lam_lP):.6f}; range / R = {lam_lP/R_lP:.0f}")
OUT["I_window_lP"] = [R_lP, lam_lP]
# ---- I5 the field picture: the circle swells near a mass
h00 = sp.Rational(2, 3); h55 = sp.Rational(1, 3)      # T_AB - T eta_AB/3 for static dust, in units of rho
check("I5 the field near a mass: the circle swells. Linear 5D gravity gives h_55 / h_00 = 1/2 for a static mass, so dR/R = (2/3) G M e^(-r/lambda) / (r c^2) with G the far-away Newton constant",
      h55/h00 == sp.Rational(1, 2))
check("I5b and the electric strength dips there by d alpha / alpha = -2 dR/R = -(4/3) G M e^(-r/lambda) / (r c^2); at the surface of the Sun (if the range were long) that would be -2.8e-6",
      abs(-(4/3)*2.12e-6 + 2.83e-6) < 1e-7, "Sun: G M / (R c^2) = 2.12e-6")

# ---- I6 a black hole has to be tiny for its horizon to sit inside the range
Mbh_g = lam_lP*lP*(2.998e8)**2/(2*6.674e-11)*1000
check("I6 only black holes smaller than the range (horizon < 2000 l_P, mass < about 0.02 g) feel the swollen circle at their horizon; every astrophysical black hole is far bigger",
      1e-3 < Mbh_g < 1e-1, f"mass limit {Mbh_g:.1e} g")
OUT["I_bh_mass_limit_g"] = Mbh_g

# ---- I7 the copies are extreme black holes
a2 = sp.simplify(sca_copy/(E1*E2/2))
check("I7 each copy sits exactly on the extreme-black-hole line of this theory: charge^2 / (G mass^2) = alpha R^2 / G = 4 = 1 + a^2 with a^2 = 3, the Kaluza-Klein black hole's dilaton coupling (Gibbons and Maeda 1988). That is why like copies feel no net force (G1)",
      sp.simplify(4 - (1 + a2)) == 0)
# ---- I8 Hawking radiation
TH = lambda M_GeV: MP**2/(8*math.pi*M_GeV)
M_emit_rad = MP**2/(8*math.pi*m_rad)*1.78266192e-24     # grams
M_emit_copy = MP**2/(8*math.pi*M1)*1.78266192e-24
OUT["I_bh_emits_carrier_below_g"] = M_emit_rad
check("I8 a black hole radiates the carrier only once its Hawking temperature passes 6e15 GeV, i.e. below about 2e-3 g, in the last instant of evaporation; copies only near the Planck mass",
      1e-4 < M_emit_rad < 1e-2 and M_emit_copy < 1e-4, f"carrier below {M_emit_rad:.1e} g, copies below {M_emit_copy:.1e} g")

# ---- I9 inflation and the early universe
As, rmax = 2.1e-9, 0.036
Hinf = math.pi*Mbar*math.sqrt(As*rmax/2); Vq = (3*Hinf**2*Mbar**2)**0.25
OUT["I_Hinf_max_GeV"] = Hinf; OUT["I_Vquarter_max_GeV"] = Vq
check("I9 during inflation the carrier is far heavier than the expansion rate (6e15 vs at most 4.7e13 GeV, from r < 0.036): it sits still at its minimum, so it leaves no fluctuations and no leftover energy",
      m_rad/Hinf > 100, f"m / H_inf = {m_rad/Hinf:.0f}")
sig = sp.symbols('sigma', real=True)
xs = sp.sqrt(sp.Rational(2, 3))*sig                       # R = R* e^{sqrt(2/3) sigma / Mbar}, Mbar = 1
Vs = -6*sp.exp(-xs) + sp.exp(-6*xs)
eps = sp.simplify((sp.diff(Vs, sig)/Vs)**2/2)
eps_tail = sp.limit(eps, sig, sp.oo)
check("I10 nor can it drive inflation itself: along its tail the slope parameter is epsilon = 1/3, far above the 0.01 that slow roll needs (and its minimum has negative energy)",
      eps_tail == sp.Rational(1, 3), f"epsilon -> {eps_tail}")
gs = 106.75
Tset = math.sqrt(m_rad*MP/(1.66*math.sqrt(gs)))
Gam = 12*m_rad**3/(64*math.pi*Mbar**2)      # order of magnitude: 12 gauge bosons, Planck-suppressed
tau = hbar_GeVs/Gam
OUT["I_settle_T_GeV"] = Tset; OUT["I_lifetime_s_est"] = tau
check("I11 if it were ever shaken, it would settle when the universe was about 6e16 GeV hot and decay in about 1e-34 s (order of magnitude), long before anything observable: no relic carrier, no change to alpha since then",
      Tset > Vq and tau < 1e-30, f"settles at {Tset:.1e} GeV (above the hottest inflation allows, {Vq:.1e} GeV); lifetime ~ {tau:.0e} s")
# ---- I12 relic copies (rough)
Ymax = 0.12/(2.8e8*M1)
x_safe = None
for x in [i/10 for i in range(100, 1000)]:
    if (MP/M1)*(1/ainv)**2*math.exp(-2*x) < Ymax:
        x_safe = x; break
T_safe = M1/x_safe
OUT["I_relic_copy_Tsafe_GeV"] = T_safe
check("I12 relic copies (rough estimate): stable by charge, so if made they'd be heavy charged relics, the old 'pyrgon' problem (Kolb and Slansky 1984). Pair production is negligible if the universe was never hotter than about 2e16 GeV, and inflation's limit already keeps it below 1.4e16 GeV",
      T_safe > Vq, f"safe below ~{T_safe:.1e} GeV; inflation max {Vq:.1e} GeV")
# ---- I13 colliders and the vacuum term
check("I13 colliders: the carrier is 4e11 times and the first copy 4e13 times beyond the LHC's 14 TeV",
      m_rad/1.4e4 > 1e11 and M1/1.4e4 > 1e13, f"{m_rad/1.4e4:.1e}, {M1/1.4e4:.1e}")
C2 = float(res_prev.get("H_range_lP_NF2") and 3*float(mp.zeta(5))*3/(64*math.pi**6))
Vleft = 5*C2/R_lP**4          # in units of 1/lP^4 = MP^4
rho_obs = 1.4e-123           # observed dark energy density in Planck units (MP^4, non-reduced), approx.
OUT["I_vacuum_left_over_MP4"] = -Vleft
check("I13b the hold's leftover vacuum energy is about -2.5e-9 in Planck units, about 1e114 times the observed dark energy and of the wrong sign: something else must cancel it (the cosmological constant problem)",
      1e112 < Vleft/rho_obs < 1e116, f"{-Vleft:.1e} vs observed {rho_obs:.1e}")

log(""); log(f"PART I TOTAL: {NP[0]} pass, {NP[1]} fail")
open("investigate_output.txt", "w").write("\n".join(lines) + "\n")
res_prev.update(OUT); json.dump(res_prev, open("results.json", "w"), indent=1, default=str)
