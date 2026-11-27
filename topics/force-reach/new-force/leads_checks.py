"""Parts L and M of the new-force DRAFT (Nov 27 2026): chasing the two possible leads from parts J and K.
Lead 1 (part L): the cap on a second circle. What does it predict, and how could it be tested?
  L1 the cap for any number of extra dimensions; L2 the smallest fundamental gravity scale it allows;
  L3 how close that comes to the magnetic weak-gravity ceiling; L4 which experiments could falsify it.
Lead 2 (part M): a cell where our law fixes a strength or a scale.
  M1 the strength of a force from a 2-sphere (non-abelian, SU(2)-like), from the sphere's symmetries;
  M2 the strength window our law allows for a second circle and for a sphere next to the electric circle;
  M3 the edges in closed form; M4 the Standard Model's own strengths at the circle's scale, against the windows;
  M5 the coincidence at the top edge, checked against supersymmetric unification.
Independent of the other scripts (shared inputs from new_force_core.py): writes L_ and M_ keys to results_LM.json and leads_output.txt.
Order-of-magnitude estimates are labelled as such."""
import math
import sympy as sp
from new_force_core import lP, MP, Mbar, ainv, R_lP, write_json_atomic

lines = []; NP = [0, 0]; OUT = {}
def log(s=""):
    print(s); lines.append(s)
def check(name, ok, detail=""):
    NP[0 if ok else 1] += 1
    log(("PASS  " if ok else "FAIL  ") + name + (f"   [{detail}]" if detail else ""))

R1 = R_lP; alpha = 1/ainv
log("== L. Lead 1: the cap on a second circle, what it predicts and how to test it")
# ---- L1 the cap for n extra dimensions of equal radius R next to the electric circle (all lengths in Planck lengths).
# Gravity spreads over the whole hidden space: l_D^(D-2) = l_P^2 * (2 pi R1) * (2 pi R)^n, D = 5 + n.
# The electric circle can only be a circle if it is at least the fundamental length: l_D <= R1.
def Rmax(n, Rl=R1):
    return math.exp(((2 + n)*math.log(Rl) - math.log(2*math.pi))/n)/(2*math.pi)
n_ = sp.symbols('n', positive=True, integer=True); Rs, R1s = sp.symbols('R R1', positive=True)
lD = (R1s*2*sp.pi*(2*sp.pi*Rs)**n_)**(1/(n_ + 3))
cap1 = sp.solve(sp.Eq(lD.subs(n_, 1)**4, R1s**4), Rs)
check("L1 the cap for one extra circle reproduces part J: R2 <= R1^3/(4 pi^2) = 325 l_P",
      len(cap1) == 1 and sp.simplify(cap1[0] - R1s**3/(4*sp.pi**2)) == 0 and abs(Rmax(1) - R1**3/(4*math.pi**2)) < 1e-9,
      f"{Rmax(1):.0f} l_P")
rows = {n: Rmax(n) for n in range(1, 8)}
OUT["L_Rmax_lP"] = {str(n): r for n, r in rows.items()}; OUT["L_Rmax_m"] = {str(n): r*lP for n, r in rows.items()}
check("L1b for any number of extra dimensions the largest allowed size falls, from 325 l_P (one) towards R1/(2 pi) = 3.7 l_P (many): every extra dimension is below 6e-33 m. Our law forbids large extra dimensions of any number",
      all(rows[n] > rows[n + 1] for n in range(1, 7)) and rows[1]*lP < 6e-33 and abs(Rmax(400)/(R1/(2*math.pi)) - 1) < 0.05,
      ", ".join(f"n={n}: {r:.1f} l_P" for n, r in rows.items()))
# ---- L2 the smallest fundamental gravity scale: M_D = 1/l_D >= 1/R1
MD_min = MP/R1
OUT["L_MD_min_GeV"] = MD_min; OUT["L_MD_min_over_CMS"] = MD_min/10.7e3
check("L2 so gravity's true strength scale can be no lower than 1/R1 = 5.2e17 GeV, whatever hides in the extra space. CMS excludes ADD below 10.7 TeV (2 dimensions, 137/fb); our law puts it 5e13 times higher: no extra-dimension signal at any collider",
      abs(MD_min/5.215e17 - 1) < 0.01 and MD_min/10.7e3 > 1e13, f"{MD_min:.3e} GeV, {MD_min/10.7e3:.1e} x the CMS limit")
# ---- L3 the magnetic weak-gravity ceiling, Lambda <= g * reduced Planck mass (Arkani-Hamed, Motl, Nicolis, Vafa 2006)
g = math.sqrt(4*math.pi*alpha); ceil = g*Mbar
check("L3 the magnetic weak-gravity ceiling for electromagnetism is g x reduced Planck mass = 7.4e17 GeV, exactly sqrt(2) above our floor (MP/sqrt(4 x 137) against MP/sqrt(2 x 137)). Floor and ceiling together pin the scale where gravity changes to 5.2e17 to 7.4e17 GeV",
      abs(ceil/MD_min - math.sqrt(2)) < 1e-9, f"ceiling {ceil:.3e} GeV, ratio {ceil/MD_min:.6f}")
OUT["L_MWGC_ceiling_GeV"] = ceil
# ---- L4 falsifiers: what a positive signal would mean
R_DD_low, R_DD_high = 0.1e-6, 10e-6             # dark dimension diameter range, Montero-Vafa-Valenzuela (2022, 2026)
EW_R = 30e-6                                       # Eot-Wash 2020: largest extra dimension, toroidal radius < 30 um
check("L4 the dark dimension (0.1 to 10 um, still allowed: Eot-Wash finds every extra dimension below 30 um) is 1e25 to 1e27 times over our cap. A short-range gravity signal at microns, or any collider sign of extra dimensions, would rule out electromagnetism as gravity's circle",
      R_DD_high < EW_R and R_DD_low/(rows[1]*lP) > 1e25, f"{R_DD_low/(rows[1]*lP):.0e} to {R_DD_high/(rows[1]*lP):.0e} x the cap")
OUT["L_DD_over_cap"] = [R_DD_low/(rows[1]*lP), R_DD_high/(rows[1]*lP)]

log("\n== M. Lead 2: does our law fix a strength or a scale?")
# ---- M1 a force from a 2-sphere: the strength follows from how long the symmetry's arrows are, averaged over the sphere.
# Same rule as the circle, alpha = 4 G / <|K|^2>, with K the arrow of the rotation whose charges are whole numbers.
th, rho = sp.symbols('theta rho', positive=True)
avgK2_S2 = sp.integrate((rho*sp.sin(th))**2*rho**2*sp.sin(th), (th, 0, sp.pi))*2*sp.pi/(4*sp.pi*rho**2)
check("M1 on a sphere of radius rho the rotation arrow has length rho sin(theta); its square averages 2 rho^2/3, so the sphere's SU(2)-like force has alpha = 4G/(2 rho^2/3) = 6G/rho^2 (a circle of radius R gives 4G/R^2)",
      sp.simplify(avgK2_S2 - sp.Rational(2, 3)*rho**2) == 0, f"<|K|^2> = {avgK2_S2}")
# ---- M2 strength windows next to the electric circle (both sizes above the fundamental length)
R2max = R1**3/(4*math.pi**2); R2min = (4*math.pi**2*R1)**(1/3)
rmax = R1**2/math.sqrt(8*math.pi**2); rmin = (8*math.pi**2*R1)**(1/3)
c_lo, c_hi = 4/R2max**2, 4/R2min**2; s_lo, s_hi = 6/rmax**2, 6/rmin**2
OUT.update(M_circle_window=[c_lo, c_hi], M_sphere_window=[s_lo, s_hi], M_sphere_rho_range_lP=[rmin, rmax])
check("M2 our law fixes each extra force's strength to a window: a second circle 1/26,400 to 1/23.7; a 2-sphere (radius 12.3 to 61.7 l_P, gravity spread over 7 dimensions) 1/634 to 1/25.1. Neither can be arbitrarily weak or strong",
      abs(1/c_lo - 26400) < 100 and 23 < 1/c_hi < 24.5 and 600 < 1/s_lo < 660 and 24.5 < 1/s_hi < 25.5,
      f"circle 1/{1/c_lo:.0f} to 1/{1/c_hi:.1f}; sphere 1/{1/s_lo:.0f} to 1/{1/s_hi:.1f}")
# ---- M3 the top edges in closed form
a_ = sp.symbols('alpha', positive=True); R1a = 2/sp.sqrt(a_)
edge_c = sp.simplify(4/((4*sp.pi**2*R1a)**sp.Rational(1, 3))**2); edge_s = sp.simplify(6/((8*sp.pi**2*R1a)**sp.Rational(1, 3))**2)
check("M3 the top edges in closed form: a second circle at most (alpha/pi^4)^(1/3), a sphere at most 6 (alpha/(256 pi^4))^(1/3); both grow only as the cube root of electromagnetism's strength",
      sp.simplify(edge_c - (a_/sp.pi**4)**sp.Rational(1, 3)) == 0 and sp.simplify(edge_s - 6*(a_/(256*sp.pi**4))**sp.Rational(1, 3)) == 0
      and abs(float(edge_c.subs(a_, alpha)) - c_hi) < 1e-12, f"circle {edge_c}, sphere {edge_s}")
# ---- M4 the Standard Model's own strengths, run (one loop) up to the circle's scale 1/R1 = 5.2e17 GeV
MZ = 91.1876; ainvZ = 127.951; s2w = 0.23122; a3Z = 0.1180
inv0 = [0.6*(1 - s2w)*ainvZ, s2w*ainvZ, 1/a3Z]; b_sm = [41/10, -19/6, -7]
L = math.log(MD_min/MZ)/(2*math.pi)
inv_KK = [i - b*L for i, b in zip(inv0, b_sm)]
OUT["M_SM_inv_at_KK"] = inv_KK
inside_s = all(1/s_hi < x < 1/s_lo for x in inv_KK); inside_c = all(1/c_hi < x < 1/c_lo for x in inv_KK)
check("M4 a loose test that could have failed: the Standard Model's three strengths at 5.2e17 GeV (one loop) are 1/35, 1/48 and 1/49, inside both windows. Strong-force-sized (1/3) or very weak (1/10,000) values would have missed the sphere window",
      inside_s and inside_c and 30 < inv_KK[0] < 40 and 45 < inv_KK[1] < 51 and 45 < inv_KK[2] < 52,
      f"1/alpha = {inv_KK[0]:.1f}, {inv_KK[1]:.1f}, {inv_KK[2]:.1f}")
# ---- M5 the coincidence at the top edge: supersymmetric unification (one loop, superpartners at 1 TeV)
b_ms = [33/5, 1, -3]; MS = 1000.0
Ls = math.log(MS/MZ)/(2*math.pi); inv_S = [i - b*Ls for i, b in zip(inv0, b_sm)]
tU = (inv_S[0] - inv_S[1])/(b_ms[0] - b_ms[1]); MU = MS*math.exp(2*math.pi*tU); aU_inv = inv_S[0] - b_ms[0]*tU
a3_at_U = inv_S[2] - b_ms[2]*tU
OUT.update(M_MSSM_MU_GeV=MU, M_MSSM_alphaU_inv=aU_inv, M_MSSM_a3_inv_at_MU=a3_at_U)
check("M5 a coincidence to note, not evidence: with superpartners at 1 TeV the three strengths meet near 1e16 GeV at about 1/26 (one loop; two-loop studies give about 2e16 GeV and 1/24 to 1/25), close to the top edges, 1/23.7 (circle) and 1/25.1 (sphere). The edge is where the hidden sizes equal the fundamental length",
      5e15 < MU < 4e16 and 24 < aU_inv < 27 and abs(a3_at_U - aU_inv) < 1.5 and abs(aU_inv - 1/s_hi) < 1.5,
      f"M_U = {MU:.1e} GeV, 1/alpha_U = {aU_inv:.1f} (strong {a3_at_U:.1f}); edges {1/c_hi:.1f}, {1/s_hi:.1f}")
check("M5b why it stays a coincidence for now: the edge depends on alpha only as a cube root, the windows are wide, and a smooth hidden space cannot give the left-handed particles of the Standard Model (known: Witten 1981, which also shows the smallest such space has 7 dimensions, CP2 x S2 x S1)",
      abs(math.log10(1/s_lo) - math.log10(1/s_hi)) > 1.3)

log(f"\nPARTS L+M TOTAL: {NP[0]} pass, {NP[1]} fail")
write_json_atomic("results_LM.json", OUT)
open("leads_output.txt", "w").write("\n".join(lines) + "\n")
