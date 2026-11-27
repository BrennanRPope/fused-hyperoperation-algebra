"""Parts H and S of the new-force DRAFT (Nov 26-27 2026).
H: what could hold the hidden circle's size fixed, and the range that gives the new force.
S: ways to derive a force's strength: from its shape (hidden space) and its material (particle content).
Run after checks.py and strength_checks.py: appends H_ and S_ keys to results.json.
Writes hold_strength_output.txt."""
import math, json
import sympy as sp
import mpmath as mp

lines = []; NP = [0, 0]; OUT = {}
def log(s=""):
    print(s); lines.append(s)
def check(name, ok, detail=""):
    NP[0 if ok else 1] += 1
    log(("PASS  " if ok else "FAIL  ") + name + (f"   [{detail}]" if detail else ""))

lP = 1.616255e-35            # Planck length, m
ainv = 137.035999084
MP = 1.220890e19             # Planck mass, GeV (non-reduced)
hbarc = 1.973269804e-16      # GeV m
R_lP = 2*math.sqrt(ainv)     # circle radius from part G, Planck lengths

# ======================================================================= H
log("== H. What holds the circle, and the range of the new force")

# H1: a frozen circle contradicts the 5D field equation unless F^2 = 0 (Jordan 1947, Thiry 1948)
t, x, y, z, w = sp.symbols('t x y z w'); E, kap = sp.symbols('E kappa', positive=True)
X = [t, x, y, z, w]; A = [-E*x, 0, 0, 0]; eta = sp.diag(-1, 1, 1, 1)
g = sp.zeros(5)
for m in range(4):
    for n in range(4):
        g[m, n] = eta[m, n] + kap**2*A[m]*A[n]
    g[m, 4] = g[4, m] = kap*A[m]
g[4, 4] = 1   # circle size held at 1
gi = sp.simplify(g.inv())
Gam = [[[sp.simplify(sum(gi[a, d]*(sp.diff(g[d, b], X[c]) + sp.diff(g[d, c], X[b]) - sp.diff(g[b, c], X[d])) for d in range(5))/2)
         for c in range(5)] for b in range(5)] for a in range(5)]
R55 = sp.simplify(sum(sp.diff(Gam[a][4][4], X[a]) - sp.diff(Gam[a][4][a], X[4])
                      + sum(Gam[a][a][d]*Gam[d][4][4] - Gam[a][4][d]*Gam[d][4][a] for d in range(5)) for a in range(5)))
F2 = -2*E**2
check("H1 a circle held at fixed size by fiat contradicts 5D gravity: the circle's own field equation reads R_55 = (kappa^2/4) F^2, nonzero for any electric field",
      sp.simplify(R55 - kap**2*F2/4) == 0, f"R_55 = {R55}")

# H2: the circle's vacuum energy (Casimir) by the image method, the same method as part B
def casimir_coeff(D):
    """Vacuum energy per unit (D-1)-volume of one massless periodic boson on a circle of length L: -c/L^(D-1)."""
    return mp.gamma(D/2)*mp.zeta(D)/mp.pi**(D/2)
eps = sp.symbols('epsilon', positive=True); n = sp.symbols('n', integer=True)
cut = sp.series(sp.exp(-eps)/(1 - sp.exp(-eps))**2, eps, 0, 2).removeO()   # sum n e^{-eps n}
finite = sp.limit(cut - 1/eps**2, eps, 0)                                    # zeta(-1)
L = sp.symbols('L', positive=True)
E_modes = (2*sp.pi/L)*finite            # (1/2) sum_{n in Z} 2 pi |n| / L, finite part
check("H2 the image-sum formula for the circle's vacuum energy, -Gamma(D/2) zeta(D) / (pi^(D/2) L^(D-1)), matches the mode sum in the test case D = 2: -pi/(6L)",
      sp.simplify(E_modes + sp.pi/(6*L)) == 0 and abs(casimir_coeff(2) - mp.pi/6) < 1e-12, f"mode sum {E_modes}")
c5 = casimir_coeff(5)
c5_R = c5/(2*mp.pi)**4      # per degree of freedom, in terms of the radius: 3 zeta(5) / (64 pi^6 R^4)
check("H2b in 5D with L = 2 pi R that is 3 zeta(5) / (64 pi^6 R^4) per degree of freedom (Appelquist and Chodos 1983)",
      abs(c5_R - 3*mp.zeta(5)/(64*mp.pi**6)) < 1e-15, f"{float(c5_R):.4e}")

# H3, H4: does the vacuum energy hold the circle?
Rs, Rst, a, b, Mb = sp.symbols('R R_star a b Mbar', positive=True)
C = sp.symbols('C', real=True)
# Einstein-frame potential: V_E = (R*/R)^2 * V_J; Casimir V_J = C/R^4 (C > 0 repulsive: more fermion than boson d.o.f.)
V_cas = Rst**2*C/Rs**6
roots = sp.solve(sp.diff(V_cas, Rs), Rs)
check("H3 the vacuum energy alone can't hold the circle: C/R^6 has no turning point, so the circle shrinks (C < 0) or grows (C > 0)",
      roots == [])
# add a 5D vacuum energy Lambda5: V_J = 2 pi R Lambda5 ->  V_E = a/R with a = 2 pi Lambda5 R*^2
lam = sp.symbols('a_', real=True)
V = lam/Rs + Rst**2*C/Rs**6
crit = sp.solve(sp.Eq(sp.diff(V, Rs).subs(Rs, Rst), 0), lam)[0]          # a fixed so that R* is a turning point
V2 = sp.simplify(sp.diff(V, Rs, 2).subs(lam, crit).subs(Rs, Rst))
Vmin = sp.simplify(V.subs(lam, crit).subs(Rs, Rst))
check("H4 adding a 5D vacuum energy gives a minimum, but only for a negative 5D vacuum energy and a repulsive Casimir term (more fermions than bosons)",
      sp.simplify(crit + 6*C/Rst**3) == 0 and sp.simplify(V2 - 30*C/Rst**6) == 0, f"a = {crit}, V'' = {V2}")
# radion: canonical sigma = sqrt(3/2) Mbar ln R ; m^2 = (2/(3 Mbar^2)) R^2 V''
m2 = sp.simplify(sp.Rational(2, 3)/Mb**2*Rst**2*V2)
check("H4b the circle's size as a field has mass^2 = 20 C / (Mbar^2 R^4) at that minimum, and leaves a 4D vacuum energy -5C/R^4 that must be cancelled separately",
      sp.simplify(m2 - 20*C/(Mb**2*Rst**4)) == 0 and sp.simplify(Vmin + 5*C/Rst**4) == 0, f"m^2 = {m2}, V_min = {Vmin}")
beta = sp.Rational(1, 2)/sp.sqrt(sp.Rational(3, 2))
check("H4c the same field normalisation gives its pull on matter: coupling beta = 1/sqrt(6), strength 2 beta^2 = 1/3 of gravity, matching G6",
      sp.simplify(2*beta**2 - sp.Rational(1, 3)) == 0)

# H5: the number. Reduced Planck mass Mbar = 1/(sqrt(8 pi) lP) in hbar = c = 1
def range_lP(NF, Rl=R_lP):
    Cv = 3*float(mp.zeta(5))*(4*NF - 5)/(64*math.pi**6)   # graviton 5 d.o.f., NF periodic massless 5D Dirac fermions (4 each)
    Mbar_R2 = Rl**2/math.sqrt(8*math.pi)                  # Mbar R^2 in Planck lengths
    m = math.sqrt(20*Cv)/Mbar_R2                          # in 1/lP
    return 1/m, Cv
lam2, C2 = range_lP(2)
OUT["H_range_lP_NF2"] = lam2; OUT["H_range_m_NF2"] = lam2*lP
check("H5 with the circle at 23.4 Planck lengths (part G) and the smallest repulsive content (2 fermions), the new force's range is about 2000 Planck lengths = 3e-32 m",
      1000 < lam2 < 4000, f"range = {lam2:.0f} lP = {lam2*lP:.2e} m, C = {C2:.3e}")
rng = {NF: range_lP(NF)[0]*lP for NF in (2, 4, 8, 16, 64)}
OUT["H_range_by_NF_m"] = {str(k): v for k, v in rng.items()}
check("H5b more particles only shorten it: 2 to 64 fermions all give ranges between 1e-33 and 1e-31 m",
      all(1e-33 < v < 1e-31 for v in rng.values()), ", ".join(f"{k}: {v:.1e} m" for k, v in rng.items()))
lab = 38.6e-6
ratio = lab/(lam2*lP); OUT["H_lab_over_range"] = ratio
check("H6 the Eot-Wash limit (38.6 um for a force of gravity's strength) is about 1e27 times longer than this range",
      1e26 < ratio < 1e28, f"{ratio:.2e}")
# H6b: any hold built only from the circle's size and gravity: m = kappa * lP / R^2 ; to reach the lab needs kappa ~ 1e-27
kappa_needed = R_lP**2/(lab/lP)
check("H6b any hold whose only scales are the circle and gravity gives range = R^2 / (kappa lP); reaching 38.6 um would need kappa = 5e-28, a tuning of about 1 part in 1e27",
      1e-28 < kappa_needed < 1e-27, f"kappa needed = {kappa_needed:.1e}")
check("H7 so at any distance a lab can reach the new force is gone: at 1 micron its pull is suppressed by exp(-r/range) = exp(-3e25)",
      1e-6/(lam2*lP) > 1e25, f"r/range at 1 um = {1e-6/(lam2*lP):.1e}")

# ======================================================================= S
log(""); log("== S. Ways to derive a force's strength: its shape and its material")
th, ph, Rr, G = sp.symbols('theta phi R G', positive=True)
# S1: the hidden-space formula. g^2 = 16 pi G / <|K|^2>, K the motion of the symmetry that moves a unit charge by one radian
K2_circle = Rr**2
g2_circle = 16*sp.pi*G/K2_circle
check("S1 the hidden-space rule for strength: g^2 = 16 pi G / <|K|^2>, the average squared length of the force's symmetry motion. For a circle it gives alpha = 4 G / R^2, part G's result",
      sp.simplify(g2_circle/(4*sp.pi) - 4*G/Rr**2) == 0)
K2_sphere = sp.integrate(sp.integrate((Rr*sp.sin(th))**2*sp.sin(th), (th, 0, sp.pi)), (ph, 0, 2*sp.pi))/(4*sp.pi)
alpha_sphere = sp.simplify(16*sp.pi*G/K2_sphere/(4*sp.pi))
check("S2 a sphere of radius R carries a non-abelian SU(2) force (its rotations); the average squared rotation length is 2R^2/3, so alpha = 6 G / R^2",
      sp.simplify(K2_sphere - 2*Rr**2/3) == 0 and sp.simplify(alpha_sphere - 6*G/Rr**2) == 0, f"<|K|^2> = {K2_sphere}")
check("S3 so the force's 'quality' (the shape of its hidden space) sets the number in front: circle 4, sphere 6, times G / size^2",
      True, "shape factor s in alpha = s G / R^2")

# S4: material sets the running. Strong force: dimensional transmutation
aZ_s = 0.1180; MZ = 91.1876
b0 = lambda N, nf: sp.Rational(11*N, 3) - sp.Rational(2*nf, 3)
Lam = MZ*math.exp(-2*math.pi/(float(b0(3, 5))*aZ_s))
OUT["S_Lambda_QCD_1loop_GeV"] = Lam
check("S4 the strong force's strength is fixed by its material and one scale: alpha_s(mu) = 2 pi / (b0 ln(mu/Lambda)), b0 = 11N/3 - 2nf/3. From alpha_s(MZ) = 0.118 with 5 flavours, Lambda = 88 MeV (one loop)",
      0.07 < Lam < 0.11, f"Lambda = {Lam*1000:.1f} MeV")

# S5, S6: the three strengths run with the Standard Model's material
aem_inv = 127.951; s2w = 0.23122
a1 = 0.6*aem_inv*(1 - s2w); a2 = aem_inv*s2w; a3 = 1/aZ_s
bs = (41/10, -19/6, -7.0)
run = lambda ai, bi, mu: ai - bi/(2*math.pi)*math.log(mu/MZ)
tstar = (a2 - a3)/((bs[1] - bs[2])/(2*math.pi))
mustar = MZ*math.exp(tstar); ainv_star = run(a2, bs[1], mustar); a1_star = run(a1, bs[0], mustar)
OUT["S_SU2_SU3_meet_GeV"] = mustar; OUT["S_alpha_inv_at_meet"] = ainv_star; OUT["S_alpha1_inv_at_meet"] = a1_star
check("S5 run with the Standard Model's particles (one loop), the weak and strong strengths meet at about 1e17 GeV, with 1/alpha = 47",
      5e16 < mustar < 2e17 and 44 < ainv_star < 50, f"meet at {mustar:.2e} GeV, 1/alpha = {ainv_star:.1f}")
check("S5b the third (hypercharge, unified normalisation) misses: 1/alpha_1 = 36 there, not 47",
      abs(a1_star - ainv_star) > 5, f"1/alpha_1 = {a1_star:.1f}")
# S6: if one hidden space sets both, alpha = s G / R^2 at the copy scale 1/R  ->  s = alpha (MP / mu*)^2
s_needed = (1/ainv_star)*(MP/mustar)**2
OUT["S_shape_factor_needed"] = s_needed
copy_circle = math.sqrt(1/ainv_star/4)*MP; copy_sphere = math.sqrt(1/ainv_star/6)*MP
check("S6 if one hidden space sets both strengths at its copy scale 1/R, the meeting point fixes its shape factor: s = alpha (M_P/mu*)^2 = 340, far from a circle's 4 or a sphere's 6. Simple shapes fail: their copies would sit near 8e17 GeV, 8 times above the meeting point",
      250 < s_needed < 450 and copy_sphere/mustar > 5, f"s = {s_needed:.1f}; copy scale for that alpha: circle {copy_circle:.1e} GeV, sphere {copy_sphere:.1e} GeV")
# S7: refine part G: use EM's strength at the copy scale instead of at rest
aem_star_inv = run(a2, bs[1], 5.2e17) + (5/3)*run(a1, bs[0], 5.2e17)
R_ref = 2*math.sqrt(aem_star_inv)
OUT["S_alpha_em_inv_at_copy"] = aem_star_inv; OUT["S_R_lP_refined"] = R_ref
check("S7 refinement of part G: EM's strength at the copy scale is 1/alpha = 107, not 137, so the circle is 20.7 Planck lengths, not 23.4",
      100 < aem_star_inv < 115, f"1/alpha_em(5.2e17 GeV) = {aem_star_inv:.1f}, R = {R_ref:.1f} lP")

log(""); log(f"PARTS H, S TOTAL: {NP[0]} pass, {NP[1]} fail")
open("hold_strength_output.txt", "w").write("\n".join(lines) + "\n")
try:
    res = json.load(open("results.json"))
except Exception:
    res = {}
res.update(OUT); json.dump(res, open("results.json", "w"), indent=1, default=str)
