"""Part N of the new-force DRAFT (Nov 27 2026): does what holds the circle also hold a sphere at the fundamental length?
Section 12 (part M) found that a 2-sphere next to the electric circle gives an SU(2)-like force with strength between 1/634 and 1/25.1,
and that the strong edge (the sphere at the fundamental length) sits near supersymmetric unification. This asks whether a hold puts it there.
  N1 how wide the window is in sphere size; N2 where the Standard Model's strengths sit in it;
  N3 a sphere can be held by its own curvature plus its vacuum energy, with no extra term (unlike the circle, section 7);
  N4 where that hold puts it, and the strength that follows; N5 how many particles that takes (rough);
  N6 the electric circle in the same units; N7 what sets the strength in the end.
Independent of the other scripts (shared inputs from new_force_core.py): writes N_ keys to results_N.json and sphere_hold_output.txt.
Order-of-magnitude estimates are labelled as such."""
import math
import sympy as sp
import mpmath as mp
from new_force_core import ainv, R_lP, MP, write_json_atomic

lines = []; NP = [0, 0]; OUT = {}
def log(s=""):
    print(s); lines.append(s)
def check(name, ok, detail=""):
    NP[0 if ok else 1] += 1
    log(("PASS  " if ok else "FAIL  ") + name + (f"   [{detail}]" if detail else ""))

R1 = R_lP                                   # electric circle, Planck lengths
def l7(rho):                                # 7D fundamental length for circle R1 and sphere rho: l7^5 = l_P^2 (2 pi R1)(4 pi rho^2)
    return (8*math.pi**2*R1*rho**2)**0.2
def s_of_alpha(a):                          # sphere size in fundamental lengths for strength a = 6/rho^2
    rho = math.sqrt(6/a); return rho/l7(rho)
rho_edge = (8*math.pi**2*R1)**(1/3); rho_far = R1**2/math.sqrt(8*math.pi**2)
a_edge = 6/rho_edge**2

log("== N. Does the hold put a sphere at the fundamental length?")
# ---- N1 the whole window is a factor 2.64 in size, and strength goes as the -10/3 power of size in fundamental units
s, A, R1s = sp.symbols('s A R1', positive=True)
rho_s = s**sp.Rational(5, 3)*(8*sp.pi**2*R1s)**sp.Rational(1, 3)        # from rho = s l7 and l7^5 = 8 pi^2 R1 rho^2
alpha_s = sp.simplify(6/rho_s**2)
edge = alpha_s.subs(s, 1)
s_far = rho_far/l7(rho_far)
OUT.update(N_s_window=[1.0, s_far], N_alpha_edge=a_edge)
check("N1 in fundamental lengths the whole sphere window is narrow: from 1 (strength 1/25.1) to 2.64 (1/634, where the fundamental length reaches the electric circle). Strength goes as size^(-10/3), so a 10% change in size moves the strength by 37%",
      sp.simplify(alpha_s/edge - s**sp.Rational(-10, 3)) == 0 and abs(s_far - R1/l7(rho_far)*rho_far/R1) < 1e-12 and 2.6 < s_far < 2.7
      and abs(1.1**(10/3) - 1.374) < 0.01, f"s from 1 to {s_far:.3f}; alpha = alpha_edge * s^(-10/3)")
# ---- N2 the Standard Model's strengths at 5.2e17 GeV (part M, one loop) sit just above the edge
MZ = 91.1876; ainvZ = 127.951; s2w = 0.23122; a3Z = 0.1180
inv0 = [0.6*(1 - s2w)*ainvZ, s2w*ainvZ, 1/a3Z]; b_sm = [41/10, -19/6, -7]
Lg = math.log((MP/R1)/MZ)/(2*math.pi); inv_KK = [i - b*Lg for i, b in zip(inv0, b_sm)]
s_SM = [s_of_alpha(1/x) for x in inv_KK]
OUT["N_s_SM"] = s_SM
check("N2 read as spheres, the Standard Model's strengths at 5.2e17 GeV (1/35, 1/48, 1/49) need spheres only 1.11, 1.21 and 1.22 times the fundamental length: all within a quarter of the edge",
      all(1.05 < x < 1.3 for x in s_SM), ", ".join(f"1/{x:.0f}: s = {y:.2f}" for x, y in zip(inv_KK, s_SM)))
# ---- N3 the hold: in the 4D frame a sphere's size feels its own curvature and its vacuum energy.
# Jordan potential per 4D volume ~ V_int (-R_int + K/rho^7) with R_int = 2/rho^2, V_int ~ rho^2 (units 16 pi G7 = 1);
# the 4D frame multiplies by (V0/V_int)^2, giving U(rho) = -2/rho^4 + K/rho^9.
r, K = sp.symbols('rho K', positive=True)
U = -2/r**4 + K/r**9
crit = sp.solve(sp.diff(U, r), r)
rstar = (sp.Rational(9, 8)*K)**sp.Rational(1, 5)
U2 = sp.simplify(sp.diff(U, r, 2).subs(r, rstar)); Umin = sp.simplify(U.subs(r, rstar))
check("N3 a sphere can be held by its own curvature plus its vacuum energy, with no extra term: U = -2/rho^4 + K/rho^9 has a stable minimum at rho^5 = 9K/8 when K > 0 (more fermions than bosons, as for the circle). The circle has no curvature, which is why section 7 needed a 5D vacuum energy",
      len(crit) == 1 and sp.simplify(crit[0]**5 - sp.Rational(9, 8)*K) == 0 and sp.simplify(U2*rstar**6 - 40) == 0,
      f"rho*^5 = 9K/8, U'' = 40/rho*^6 > 0")
check("N3b it leaves a negative 4D vacuum energy, -(10/9)/rho*^4, which must be cancelled separately, the same open problem as the circle's hold (section 7)",
      sp.simplify(Umin*rstar**4 + sp.Rational(10, 9)) == 0, f"U_min = {sp.simplify(Umin*rstar**4)}/rho*^4")
# ---- N4 where it lands: units with 16 pi G7 = 1 are (16 pi)^(1/5) = 2.19 fundamental lengths, so s* = 2.19 (9K/8)^(1/5)
u = (16*math.pi)**0.2
def s_star(Kv): return u*(9*Kv/8)**0.2
def K_for_s(sv): return (8/9)*(sv/u)**5
K_lo, K_hi = K_for_s(1.0), K_for_s(s_far); K_weak = K_for_s(s_SM[1])
OUT.update(N_unit_over_l7=u, N_K_window=[K_lo, K_hi], N_K_weak=K_weak)
Ks = sp.symbols('K', positive=True)
alpha_K = sp.simplify((u*(sp.Rational(9, 8)*Ks)**sp.Rational(1, 5))**sp.Rational(-10, 3))
check("N4 so the hold sets the strength through K alone, alpha = alpha_edge x s*^(-10/3) with s* = 2.19 (9K/8)^(1/5), which goes as K^(-2/3). The sphere lands in the window for K between 0.018 and 2.3; the weak strength (1/48) needs K = 0.047",
      abs(s_star(K_lo) - 1) < 1e-9 and abs(s_star(K_hi) - s_far) < 1e-9 and abs(float(sp.simplify(sp.diff(sp.log(alpha_K), Ks)*Ks)) + 2/3) < 1e-12
      and 0.016 < K_lo < 0.020 and 2.0 < K_hi < 2.5 and 0.04 < K_weak < 0.055, f"K window {K_lo:.3f} to {K_hi:.2f}; weak needs {K_weak:.3f}")
# ---- N5 how many particles: rough per-particle K from the flat image formula in 7D, 16 pi G7 = 1 units
c7 = float(mp.gamma(3.5)*mp.zeta(7)/mp.pi**3.5)                 # vacuum energy coefficient per degree of freedom, -c/L^7 for period L
K_dof_2pi = c7/(2*math.pi)**7; K_dof_pi = c7/math.pi**7            # period 2 pi rho (around the equator) or pi rho (pole to pole)
N_lo = [K_lo/K_dof_pi, K_lo/K_dof_2pi]; N_weak = [K_weak/K_dof_pi, K_weak/K_dof_2pi]
OUT.update(N_K_per_dof=[K_dof_pi, K_dof_2pi], N_dof_needed_window=N_lo, N_dof_needed_weak=N_weak)
check("N5 rough: one particle contributes K of about 1.6e-7 to 2e-5 (flat-space image formula with period 2 pi rho or pi rho), so even reaching the window takes about 900 to 110,000 more fermion than boson states, and the weak strength about 2,300 to 300,000; with a few particles the sphere sits below the fundamental length, where none of this can be trusted (known pattern: Candelas and Weinberg 1984 needed many fields for a small coupling)",
      700 < N_lo[0] < 1100 and 9e4 < N_lo[1] < 1.4e5 and 2e3 < N_weak[0] < 2.6e3 and 2.5e5 < N_weak[1] < 3.3e5 and s_star(10*K_dof_2pi) < 0.5,
      f"K per state {K_dof_2pi:.1e} to {K_dof_pi:.1e}; window needs {N_lo[0]:.0f} to {N_lo[1]:.0f}; weak needs {N_weak[0]:.0f} to {N_weak[1]:.0f}; 10 states give s* = {s_star(10*K_dof_2pi):.2f}")
# ---- N6 the electric circle in gravity's natural units
l5 = (2*math.pi*R1)**(1/3); u5 = (16*math.pi)**(1/3)
OUT.update(N_R1_over_l5=R1/l5, N_R1_in_16piG5_units=R1/(u5*l5))
check("N6 the electric circle itself sits at 4.4 fundamental lengths (one circle), which is 1.2 in units where 16 pi G5 = 1: also an order-one size, not a big hierarchy",
      4.3 < R1/l5 < 4.5 and 1.1 < R1/(u5*l5) < 1.3, f"R1/l5 = {R1/l5:.2f}, R1/(16 pi G5)^(1/3) = {R1/(u5*l5):.2f}")
# ---- N7 the upshot
check("N7 so the hold does not single out the edge: every strength in the window is a sphere within 2.64 fundamental lengths, and which one the hold picks is set by K, the particle content. Like the conformal window (part K) and section 8, the strength comes down to the 'material'",
      s_far < 3 and abs(float(sp.simplify(sp.diff(sp.log(alpha_K), Ks)*Ks)) + 2/3) < 1e-12)

log(f"\nPART N TOTAL: {NP[0]} pass, {NP[1]} fail")
write_json_atomic("results_N.json", OUT)
open("sphere_hold_output.txt", "w").write("\n".join(lines) + "\n")
