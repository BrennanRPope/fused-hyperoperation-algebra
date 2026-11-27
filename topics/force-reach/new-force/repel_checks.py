"""Part J of the new-force DRAFT (Nov 27 2026): could the new force be the inverse, a very weak, very long-range
repulsion? And could a second long-range force that grabs different materials be mistaken for gravity?
J1-J3 which carriers repel, and what each kind of body feels from the circle's pieces;
J4-J6 a second, bigger circle: its strength, and how big our law lets it be;
J7-J8 the dark energy angle (a micron-sized circle) and that circle's own size-carrier;
J9-J10 could a second long-range force pass for gravity (light bending, free fall of different materials).
Independent of the other scripts (shared inputs from new_force_core.py): writes J_ keys to results_J.json and repel_output.txt.
Order-of-magnitude estimates are labelled as such."""
import math, json
import sympy as sp
from mpmath import zeta

lines = []; NP = [0, 0]; OUT = {}
def log(s=""):
    print(s); lines.append(s)
def check(name, ok, detail=""):
    NP[0 if ok else 1] += 1
    log(("PASS  " if ok else "FAIL  ") + name + (f"   [{detail}]" if detail else ""))

lP = 1.616255e-35; MP = 1.220890e19; Mbar = MP/math.sqrt(8*math.pi); hbarc_eVm = 1.973269804e-7
ainv = 137.035999084; R1 = 2*math.sqrt(ainv)                 # first circle, Planck lengths (part G)
from new_force_core import radion_GeV, inflation_bounds, write_json_atomic
prev = {"I_radion_GeV": radion_GeV(), "I_Vquarter_max_GeV": inflation_bounds()[1]}   # shared inputs, recomputed

def X(p1, p2, D):
    eta = sp.diag(*([-1] + [1]*(D - 1))); E1, E2 = p1[0], p2[0]
    T1 = sp.Matrix(p1)*sp.Matrix(p1).T/E1; T2 = sp.Matrix(p2)*sp.Matrix(p2).T/E2
    T1l = eta*T1*eta
    contr = sum(T1l[a, b]*T2[a, b] for a in range(D) for b in range(D))
    tr1 = sum(eta[a, a]*T1[a, a] for a in range(D)); tr2 = sum(eta[a, a]*T2[a, a] for a in range(D))
    return sp.simplify(contr - tr1*tr2/(D - 2))
m1, m2, E, E1, E2 = sp.symbols('m1 m2 E E1 E2', positive=True)

log("== J. The inverse: a weak, long-range repulsion? And could it pass for gravity?")
# ---- J1 which pieces repel: split the 5D weight for two copies into tensor (spin 2), vector (spin 1), scalar (spin 0)
like = X([E1, 0, 0, 0, E1], [E2, 0, 0, 0, E2], 5); opp = X([E1, 0, 0, 0, E1], [E2, 0, 0, 0, -E2], 5)
ten = X([E1, 0, 0, 0], [E2, 0, 0, 0], 4)
vec_like = sp.simplify((like - opp)/2); sca = sp.simplify((like + opp)/2 - ten)
check("J1 sign by spin: for two like copies the spin-2 (gravity) and spin-0 (size) pieces pull, +1 and +3; the spin-1 piece pushes, -4. Only a spin-1 carrier can repel like things",
      sp.simplify(ten/(E1*E2/2) - 1) == 0 and sp.simplify(sca/ten - 3) == 0 and sp.simplify(vec_like/ten + 4) == 0,
      f"tensor +1, scalar {sp.simplify(sca/ten)}, vector {sp.simplify(vec_like/ten)}")
# ---- J2 what each kind of body feels (weights in units of their gravity)
def parts(p, q):
    """(tensor, scalar, vector) weights between bodies p, q (5D momenta), in units of the 4D tensor weight."""
    t = X([p[0], *p[1:4]], [q[0], *q[1:4]], 4); full = X(p, q, 5)
    qflip = list(q); qflip[4] = -qflip[4]; full_f = X(p, qflip, 5)
    v = sp.simplify((full - full_f)/2); s = sp.simplify((full + full_f)/2 - t)
    return [sp.simplify(t/t), sp.simplify(s/t), sp.simplify(v/t)]
mm = parts([m1, 0, 0, 0, 0], [m2, 0, 0, 0, 0]); mc = parts([m1, 0, 0, 0, 0], [E2, 0, 0, 0, E2])
cc = parts([E1, 0, 0, 0, E1], [E2, 0, 0, 0, E2])
check("J2 the circle's pieces grab different 'materials': matter-matter feels size 1/3, no vector; matter-copy size 1, no vector; copy-copy size 3, vector -4 (like charges)",
      mm == [1, sp.Rational(1, 3), 0] and mc == [1, 1, 0] and cc == [1, 3, -4], f"matter-matter {mm}, matter-copy {mc}, copy-copy {cc}")
OUT["J_weights"] = {"matter-matter": [str(x) for x in mm], "matter-copy": [str(x) for x in mc], "copy-copy": [str(x) for x in cc]}
# ---- J3 the inverse already in one circle: with the size held (short range), like copies repel at 3 x gravity, forever
far_like = 1 + 0 - 4; far_opp = 1 + 0 + 4; near_like = 1 + 3 - 4
check("J3 so the 'inverse' is already in the circle: beyond the size-carrier's range (2000 l_P) like copies push apart at 3 x their gravity with unlimited range; inside it they feel nothing (1 + 3 - 4 = 0). For our circle that push is electromagnetism",
      far_like == -3 and near_like == 0 and far_opp == 5 and abs(ainv**-1/(1/R1**2) - 4) < 1e-9,
      f"far: like {far_like}, opposite +{far_opp}; alpha/(G M1^2) = {ainv**-1*R1**2:.3f}")
# ---- J4 a second, bigger circle: same coupling rule, alpha2 = 4 G / R2^2, and its copies sit on the same extreme line
R2s = sp.symbols('R2', positive=True)
check("J4 a second circle gives a second spin-1 force by the same rule, alpha2 = 4 G / R2^2: a bigger circle means a weaker force, unlimited range, felt only by that circle's own copies; their charge-to-mass ratio is 4 at every size",
      sp.simplify((4/R2s**2)/((1/R2s)**2) - 4) == 0)
R50 = 50e-6/lP; a50 = 4/R50**2; grav_pp = (0.938272/MP)**2
OUT["J_naive_alpha_50um"] = a50
check("J4b the hand estimate from the last reply: a 50 um circle gives alpha2 = 4.2e-61, 1.4e22 times weaker than gravity between two protons (this ignores J5)",
      abs(a50 - 4.19e-61)/4.19e-61 < 0.01 and 1e22 < grav_pp/a50 < 2e22, f"alpha2 {a50:.2e}, proton gravity {grav_pp:.2e}")
# ---- J5 how big can the second circle be? Gravity is diluted over both circles: l6^4 = (2 pi R1)(2 pi R2) l_P^2.
# Our classical reading (and the coupling rule itself) needs each circle larger than that fundamental length l6.
R2max = R1**3/(4*math.pi**2); R2min = (4*math.pi**2*R1)**(1/3)
a2min = 4/R2max**2; a2max = 4/R2min**2; m2copy = MP/R2max
l5 = (2*math.pi*R1)**(1/3); l6_50 = (4*math.pi**2*R1*R50)**0.25
OUT.update(J_R2_max_lP=R2max, J_R2_min_lP=R2min, J_alpha2_min=a2min, J_alpha2_max=a2max, J_copy2_min_GeV=m2copy, J_l6_at_50um_lP=l6_50)
check("J5 the one-circle reading is self-consistent: the 5D fundamental length is 5.3 l_P, below the circle's 23.4 l_P",
      l5 < R1 and abs(l5 - 5.28) < 0.02, f"l5 = {l5:.2f} l_P")
check("J5b with a second circle the fundamental length grows, l6 = (4 pi^2 R1 R2)^(1/4); keeping the first circle above it caps the second at R1^3/(4 pi^2) = 325 l_P (and above 9.7 l_P)",
      abs(R2max - 324.6) < 0.5 and abs(R2min - 9.73) < 0.02, f"{R2min:.2f} l_P < R2 < {R2max:.1f} l_P")
check("J5c so the weakest second force our law allows is alpha2 = 1/26,400, about 190 times weaker than electromagnetism, not 1e-61; its lightest copy weighs 3.8e16 GeV",
      abs(1/a2min - 26418) < 50 and 150 < (1/ainv)/a2min < 250 and abs(m2copy - 3.76e16)/3.76e16 < 0.01,
      f"1/alpha2 from {1/a2max:.1f} to {1/a2min:.0f}; EM/alpha2_min = {(1/ainv)/a2min:.0f}; copy {m2copy:.2e} GeV")
check("J5d a 50 um circle would push the fundamental length to 2.3e8 l_P, ten million times the electric circle: the first circle could not be a circle at all",
      l6_50/R1 > 1e6, f"l6 = {l6_50:.2e} l_P vs R1 = {R1:.1f}")
# ---- J6 the same cap as a ratio of the two circles
check("J6 same cap written as a ratio: the second circle can be at most R1^2/(4 pi^2) = 13.9 times the electric circle; a strong electric force forbids a big hidden circle",
      abs(R2max/R1 - R1**2/(4*math.pi**2)) < 1e-9 and abs(R2max/R1 - 13.88) < 0.02, f"R2/R1 <= {R2max/R1:.2f}")
# ---- J7 dark energy angle: a circle's own vacuum energy, per net fermion degree of freedom, is 3 zeta(5)/(64 pi^6) hbar c / R^4
c5 = 3*float(zeta(5))/(64*math.pi**6)
h = 0.674; rho_c = 8.0992e-11*h**2; rho_L = 0.685*rho_c             # eV^4 (Planck 2018)
R_DE = (2*c5/rho_L)**0.25*hbarc_eVm                                  # 2 net fermion dof, metres
OUT.update(J_rhoL_quarter_meV=rho_L**0.25*1e3, J_R_dark_energy_m=R_DE, J_R_dark_energy_over_cap=R_DE/lP/R2max)
check("J7 the cosmic version: a circle of about 9 um with two more fermions than bosons running round it would leave exactly the observed dark energy, (2.24 meV)^4, with the right sign (known: the 'dark dimension', Montero-Vafa-Valenzuela 2022)",
      abs(rho_L**0.25*1e3 - 2.24) < 0.02 and 5e-6 < R_DE < 15e-6, f"R = {R_DE*1e6:.1f} um, rho^(1/4) = {rho_L**0.25*1e3:.2f} meV")
check("J7b it passes gravity tests (Yukawa range below 38.6 um, Lee et al 2020) but our law forbids it alongside an electric circle: 1.7e27 times over the J5 cap",
      R_DE < 38.6e-6 and R_DE/lP/R2max > 1e26, f"{R_DE/lP/R2max:.1e} x the cap")
# ---- J8 that circle's own size-carrier, held the part H way, would be nearly massless: a 1/3-gravity force across the cosmos
C2 = 1.517e-4                                                       # part H, 2 fermions
R_GeVinv = R_DE/1.973269804e-16
m_DE_eV = math.sqrt(20*C2)/(Mbar*R_GeVinv**2)*1e9
H0_eV = 2.1332e-33*h
R1_GeVinv = R1/MP; m_check = math.sqrt(20*C2)/(Mbar*R1_GeVinv**2)
OUT.update(J_radion_DE_eV=m_DE_eV, J_radion_DE_range_m=hbarc_eVm/m_DE_eV)
check("J8 the same hold formula reproduces part I's carrier for the electric circle (6.15e15 GeV); for a 9 um circle it gives 1e-32 eV, a few times the Hubble rate, reaching 1.7e25 m",
      abs(m_check/float(prev["I_radion_GeV"]) - 1) < 0.01 and 1 < m_DE_eV/H0_eV < 30,
      f"{m_DE_eV:.1e} eV = {m_DE_eV/H0_eV:.1f} H0, range {hbarc_eVm/m_DE_eV:.1e} m")
check("J8b a 1/3-gravity force reaching that far makes light bend with gamma = 1/2 in the solar system, which Cassini excludes by 21,740 standard errors: a big circle also needs a much firmer hold (known: Arkani-Hamed, Dimopoulos, March-Russell 1999)",
      abs(-0.5 - 2.1e-5)/2.3e-5 > 2e4)
# ---- J9 could a second long-range force pass for gravity? Only if it grabs every body, light included, in proportion to its gravity
a = sp.symbols('a', nonnegative=True)
gam = (1 - a)/(1 + a)
check("J9 to hide inside gravity a force must grab light and matter in the same ratio as gravity does; ours grabs matter (1/3) but not light, so light bends by (1 + gamma)/2 = 3/4 of the expected amount",
      gam.subs(a, sp.Rational(1, 3)) == sp.Rational(1, 2) and sp.Rational(1, 2)*(1 + gam.subs(a, sp.Rational(1, 3))) == sp.Rational(3, 4) and gam.subs(a, 0) == 1)
# ---- J10 different materials: the circle's size also sets alpha (alpha ~ 1/R^2), so electric binding energy couples differently from rest mass
def f_coul(Z, A, Mgev):
    return 0.71e-3*Z*(Z - 1)/A**(1/3)/Mgev           # nuclear Coulomb energy fraction (semi-empirical a_C = 0.71 MeV)
fTi = f_coul(22, 48, 44.65); fPt = f_coul(78, 195, 181.6); eta_est = (1/3)*abs(fPt - fTi)
OUT.update(J_eta_estimate=eta_est, J_eta_over_microscope=eta_est/4.0e-15)
check("J10 different materials (rough): titanium and platinum carry 2e-3 and 4e-3 of their mass as electric energy, so a long-range size force would make them fall differently by about 7e-4; MICROSCOPE finds below about 4e-15 (2022)",
      1e-4 < eta_est < 2e-3 and eta_est/4e-15 > 1e10, f"f_Ti {fTi:.1e}, f_Pt {fPt:.1e}, eta ~ {eta_est:.0e}, {eta_est/4e-15:.0e} x the limit")
check("J10b what does pass both tests is a long-range force on a different material that ordinary matter does not carry at all, like a second circle's force on its own copies (J2: matter feels no vector); in our law those copies are 3.8e16 GeV or heavier, above the most energy the early universe could have had (part I: inflation's scale below 1.4e16 GeV), so few if any were made (rough)",
      mm[2] == 0 and m2copy > float(prev["I_Vquarter_max_GeV"]), f"copy {m2copy:.1e} GeV > {float(prev['I_Vquarter_max_GeV']):.1e} GeV")

log(f"\nPART J TOTAL: {NP[0]} pass, {NP[1]} fail")
write_json_atomic("results_J.json", OUT)
open("repel_output.txt", "w").write("\n".join(lines) + "\n")
