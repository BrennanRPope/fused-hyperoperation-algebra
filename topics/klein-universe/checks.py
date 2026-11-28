"""DRAFT round: is space itself shaped like a Klein bottle, or twisted at all?
Asked for by Brennan (Nov 27 2026): "Did we determine if reality was really a Klein bottle?" / "This is my hypothesis that I want to check".
Scope widened on his call the same evening to twisted shapes in general and to how a twist of the shape relates to local torsion.

Checks:
 U0 (numbers)  distances in the measured universe (Planck 2018 LCDM): distance to the last-scattering sphere, its diameter,
               the observable universe's radius, and the matched-circle bound in those units.
 U1 (exact)    the Klein space E7 (the simplest right-angled member) and, for contrast, the half-turn space E2: both glue
               flat space with no fixed points; E7 has orientation-reversing loops and its orientable double cover is a
               3-torus of sides a, b, 2c; E2 twists but keeps every loop orientation-preserving. Flat (zero curvature,
               zero torsion) everywhere; the twist lives only in the holonomy of the loops.
 U2 (numbers)  what such a universe shows in the sky: a field built to respect the gluing, sampled on the last-scattering
               sphere, gives matched circle pairs; translation and half-turn pairs match when one circle is read the
               opposite way round, the Klein glide's pairs when both are read the same way round (the flipped pairs).
 U3 (numbers plus exact formula)  the size window for E7: nearest copy of us = min(a, b, 2c, sqrt(4 d^2 + c^2)), d the
               distance to the nearest mirror plane, checked by brute force over the group; with no circles seen above 10
               degrees, the shortest loop of a Klein-shaped universe can still be as short as about half the
               last-scattering diameter, if we sit off its mirror planes.
 U4 (exact plus measured numbers)  handedness: an orientation-reversing loop turns a left-handed spinor into a right-handed
               one (S gamma5 S^-1 = -gamma5 for every reflection, +gamma5 for every rotation). So the gluing needs a
               mirror symmetry of the laws of physics; the Standard Model has none (Jarlskog invariant J != 0).
 U5 (numbers, order of magnitude)  the only loophole, a handedness chosen spontaneously by a field, forces a mirror wall
               somewhere in a Klein-shaped universe; walls from any parity-breaking scale at or above the weak scale are
               too heavy by many orders of magnitude (Zel'dovich-Kobzarev-Okun 1974).
 U7 (exact)    spinors on the twisted spaces: Pin+ and Pin- structures, and the phase picked up going round a twist.
 U6 (numbers, order of magnitude)  a twist of the shape versus local torsion: the shape's twist has zero torsion and zero
               curvature at every point; Einstein-Cartan torsion is a separate, local field tied to spin, with effects only
               above the Cartan density, and lab bounds on background torsion near 1e-31 GeV.
numpy only, plus sympy for the exact gamma-matrix algebra."""
import json, io, time, itertools
from fractions import Fraction as Fr
import numpy as np
import sympy as sp

T0 = time.time(); _log = io.StringIO(); passed = failed = 0; R = {}
def say(s): print(s); _log.write(s + "\n")
def check(name, ok, detail):
    global passed, failed
    passed += bool(ok); failed += (not ok)
    say(f"{'PASS' if ok else 'FAIL'} {name}  | {detail}")

# ------------------------------------------------------------------ U0 distances
say("== U0. Distances in the measured universe (Planck 2018 LCDM)")
c_kms = 299792.458
h, Om = 0.6736, 0.3153                      # Planck 2018 VI, TT,TE,EE+lowE+lensing
Or = 2.469e-5 * (1 + 0.2271 * 3.044) / h**2 # photons + 3.044 massless-neutrino equivalents
OL = 1 - Om - Or
zstar = 1089.92
DH = c_kms / (100 * h) / 1000.0             # Hubble distance in Gpc
def chi(a_lo):
    # chi = D_H * int_{a_lo}^1 da / sqrt(Or + Om a + OL a^4), substituted a = u^2 for the a -> 0 end
    u = np.linspace(np.sqrt(a_lo), 1.0, 400001)
    a = u**2
    f = 2 * u / np.sqrt(Or + Om * a + OL * a**4)
    return DH * np.trapezoid(f, u)
chi_rec_int = chi(1 / (1 + zstar))
chi_hor = chi(0.0)
rstar, theta100 = 144.43e-3, 1.04110        # sound horizon (Gpc) and 100*theta_* from Planck 2018 VI
chi_rec = rstar / (theta100 / 100)          # what Planck's measured angle gives directly
D = 2 * chi_rec
bound = 0.985 * D                           # Vaudrevange, Starkman, Cornish, Spergel 2012; restated by COMPACT Part I (2023)
R["U0"] = dict(chi_rec_from_theta_Gpc=chi_rec, chi_rec_integrated_Gpc=chi_rec_int, lss_diameter_Gpc=D,
               particle_horizon_Gpc=chi_hor, observable_diameter_Gpc=2 * chi_hor, circle_bound_Gpc=bound,
               circle_bound_angle_deg=float(np.degrees(np.arccos(0.985))))
check("U0 (numbers; Planck 2018 parameters)",
      abs(chi_rec_int / chi_rec - 1) < 0.01 and 13.5 < chi_rec < 14.2 and 13.9 < chi_hor < 14.6,
      f"distance to the last-scattering sphere {chi_rec:.3f} Gpc from Planck's measured angle (r*/theta*), {chi_rec_int:.3f} Gpc by "
      f"integrating LCDM; its diameter D = {D:.2f} Gpc; observable universe radius {chi_hor:.2f} Gpc (diameter {2*chi_hor:.2f}). "
      f"No matched circles above 10 degrees means every copy of us is farther than 0.985 D = {bound:.2f} Gpc "
      f"(cos 10 deg = {np.cos(np.radians(10)):.4f})")

# ------------------------------------------------------------------ U1 the spaces themselves
say("== U1. The Klein space E7 and the half-turn space E2, exactly")
a_, b_, c_ = Fr(17, 10), Fr(11, 10), Fr(3, 5)
def aff(L, t): return (tuple(L), tuple(t))
def mul(g, h):   # g after h
    (Lg, tg), (Lh, th) = g, h
    L = tuple(Lg[i] * Lh[i] for i in range(3))            # diagonal linear parts only
    t = tuple(Lg[i] * th[i] + tg[i] for i in range(3))
    return (L, t)
def inv(g):
    L, t = g
    return (L, tuple(-L[i] * t[i] for i in range(3)))
I = aff((1, 1, 1), (0, 0, 0))
def ball(gens, depth):
    gs = gens + [inv(g) for g in gens]; seen = {I}; front = [I]
    for _ in range(depth):
        nxt = []
        for x in front:
            for g in gs:
                y = mul(g, x)
                if y not in seen: seen.add(y); nxt.append(y)
        front = nxt
    return seen
def free(g):
    L, t = g
    return any(L[i] == 1 and t[i] != 0 for i in range(3))   # fixed point exists iff t has no part along Fix(L)
Ta = aff((1, 1, 1), (a_, 0, 0)); Tb = aff((1, 1, 1), (0, b_, 0))
G = aff((-1, 1, 1), (0, 0, c_))            # glide: mirror x -> -x, then slide c along z
H = aff((-1, -1, 1), (0, 0, c_))           # half-turn screw: turn by 180 degrees about z, then slide c
E7 = ball([Ta, Tb, G], 6); E2 = ball([Ta, Tb, H], 6)
det = lambda g: g[0][0] * g[0][1] * g[0][2]
e7_free = all(free(g) for g in E7 if g != I); e2_free = all(free(g) for g in E2 if g != I)
e7_rev = [g for g in E7 if det(g) == -1]; e2_rev = [g for g in E2 if det(g) == -1]
def in_lattice(t): return (t[0] / a_).denominator == 1 and (t[1] / b_).denominator == 1 and (t[2] / (2 * c_)).denominator == 1
e7_pres_ok = all(g[0] == (1, 1, 1) and in_lattice(g[1]) for g in E7 if det(g) == 1)
e7_rev_ok = all((g[1][2] / c_).denominator == 1 and (g[1][2] / c_).numerator % 2 == 1 for g in e7_rev)
R["U1"] = dict(E7_elements=len(E7), E7_reversing=len(e7_rev), E2_elements=len(E2), E2_reversing=len(e2_rev))
check("U1 (exact rational arithmetic over every group element up to word length 6)",
      e7_free and e2_free and len(e7_rev) > 0 and len(e2_rev) == 0 and e7_pres_ok and e7_rev_ok,
      f"E7 (Klein space, generators: slide a, slide b, glide = mirror x->-x then slide c): {len(E7)} elements, none with a fixed "
      f"point, {len(e7_rev)} reverse orientation (each slides an odd multiple of c), and the orientation-keeping ones are exactly "
      f"the slides of the lattice (a, b, 2c), so the orientable double cover is a 3-torus. E2 (half-turn space): {len(E2)} elements, "
      f"none with a fixed point, {len(e2_rev)} reverse orientation. Both spaces are flat: every gluing map is a rigid motion, so "
      f"curvature and torsion are zero at every point; the twist is only the holonomy of a loop (det -1 for the glide, +1 for the half-turn)")

# ------------------------------------------------------------------ U2 the sky: flipped circle pairs
say("== U2. What the sky shows: matched circle pairs and which way round they match")
rng = np.random.default_rng(20261127)
A, B, C, CHI = 1.6, 1.9, 0.7, 1.0          # a deliberately small toy universe, in units of the last-scattering radius
O = np.array([0.15, 0.05, 0.10])
def mk(L, t): return (np.diag(np.array(L, float)), np.array(t, float))
def app(g, x): return g[0] @ x + g[1]
def ginv(g): Li = g[0].T; return (Li, -Li @ g[1])
def field_for(twist):
    ks = []
    for n in range(-10, 11):
        for m in range(-10, 11):
            for l in range(-16, 17):
                k = np.array([2 * np.pi * n / A, 2 * np.pi * m / B, np.pi * l / C])
                kk = np.linalg.norm(k)
                if 0 < kk < 24: ks.append(k)
    ks = np.array(ks); amp = rng.normal(size=len(ks)) / np.linalg.norm(ks, axis=1); ph = rng.uniform(0, 2 * np.pi, len(ks))
    def f(X):
        X = np.atleast_2d(X); Y = (twist[0] @ X.T).T + twist[1]
        return (np.cos(X @ ks.T + ph) + np.cos(Y @ ks.T + ph)) @ amp
    return f, len(ks)
def basis(n):
    e1 = np.cross(n, [0.3, 0.9, 0.31]); e1 /= np.linalg.norm(e1); e2 = np.cross(n, e1); return e1, e2   # (e1, e2, n) right-handed
def circle_pair(g, f, N=720):
    gO = app(g, O); s = np.linalg.norm(gO - O); al = np.arccos(s / (2 * CHI))
    n1 = (gO - O) / s; e1, e2 = basis(n1)
    phi = np.linspace(0, 2 * np.pi, N, endpoint=False)
    P = O + CHI * (np.cos(al) * n1 + np.sin(al) * (np.outer(np.cos(phi), e1) + np.outer(np.sin(phi), e2)))
    gi = ginv(g); Q = (gi[0] @ P.T).T + gi[1]
    on_sphere = np.max(np.abs(np.linalg.norm(Q - O, axis=1) - CHI))
    n2 = (app(gi, O) - O) / s; f1, f2 = basis(n2)
    rel = Q - O; psi = np.unwrap(np.arctan2(rel @ f2, rel @ f1)); sense = np.sign(np.mean(np.diff(psi)))
    T1 = f(P)
    Pc2 = O + CHI * (np.cos(al) * n2 + np.sin(al) * (np.outer(np.cos(phi), f1) + np.outer(np.sin(phi), f2)))
    T2 = f(Pc2)                                                      # the second circle read independently, on its own grid
    def S(sg):
        best = -9
        for d in range(N):
            T2r = T2[(sg * np.arange(N) + d) % N]
            best = max(best, 2 * np.mean(T1 * T2r) / (np.mean(T1**2) + np.mean(T2r**2)))
        return best
    return dict(sep=s, radius_deg=float(np.degrees(al)), on_sphere=on_sphere, sense=int(sense), S_same=S(+1), S_opposite=S(-1),
                n1=n1.tolist(), n2=n2.tolist())
GT, TA, HT = mk((-1, 1, 1), (0, 0, C)), mk((1, 1, 1), (A, 0, 0)), mk((-1, -1, 1), (0, 0, C))
fK, nK = field_for(GT); fH, nH = field_for(HT)
X = rng.normal(size=(200, 3)) * 3
inv_err = max(np.max(np.abs(fK(np.array([app(GT, x) for x in X])) - fK(X))),
              np.max(np.abs(fK(np.array([app(TA, x) for x in X])) - fK(X))),
              np.max(np.abs(fH(np.array([app(HT, x) for x in X])) - fH(X))))
pairs = {"Klein glide": circle_pair(GT, fK), "Klein slide a": circle_pair(TA, fK), "half-turn screw": circle_pair(HT, fH)}
R["U2"] = dict(modes=[nK, nH], invariance_error=float(inv_err), pairs=pairs, observer=O.tolist(), toy=dict(a=A, b=B, c=C, chi=CHI))
ok2 = (inv_err < 1e-9 and all(p["on_sphere"] < 1e-12 for p in pairs.values())
       and pairs["Klein glide"]["sense"] == 1 and pairs["Klein glide"]["S_same"] > 0.999 and pairs["Klein glide"]["S_opposite"] < 0.9
       and pairs["Klein slide a"]["sense"] == -1 and pairs["Klein slide a"]["S_opposite"] > 0.999 and pairs["Klein slide a"]["S_same"] < 0.9
       and pairs["half-turn screw"]["sense"] == -1 and pairs["half-turn screw"]["S_opposite"] > 0.999 and pairs["half-turn screw"]["S_same"] < 0.9)
check("U2 (numbers; matched-circle statistic of Cornish, Spergel and Starkman 1998)", ok2,
      f"a random field built to respect each gluing ({nK} and {nH} modes; gluing error {inv_err:.1e}) sampled on the last-scattering "
      f"sphere of a toy universe (a, b, c = {A}, {B}, {C} in last-scattering radii). Each circle is read round its own centre as seen "
      f"from us. " + "; ".join(f"{k}: copy {v['sep']:.3f} away, circles of radius {v['radius_deg']:.1f} deg, best match read the same way "
      f"round S = {v['S_same']:.3f}, read opposite ways round S = {v['S_opposite']:.3f}" for k, v in pairs.items())
      + ". So slides and the half-turn give pairs that run opposite ways round, and the Klein glide gives pairs that run the same way: "
      "the flipped pairs, which the circle searches test separately (Cornish et al. 2004; Vaudrevange et al. 2012)")

# ------------------------------------------------------------------ U3 the size window
say("== U3. How small a Klein-shaped universe can still be")
def nearest_formula(a, b, c, d): return min(a, b, 2 * c, np.sqrt(4 * d**2 + c**2))
def nearest_brute(a, b, c, O3):
    gens = [mk((1, 1, 1), (a, 0, 0)), mk((1, 1, 1), (0, b, 0)), mk((-1, 1, 1), (0, 0, c))]
    els = {}
    front = [(np.eye(3), np.zeros(3))]
    allg = list(front)
    for _ in range(5):
        nxt = []
        for g in front:
            for q in gens + [ginv(x) for x in gens]:
                y = (q[0] @ g[0], q[0] @ g[1] + q[1]); key = (tuple(np.round(np.diag(y[0]))), tuple(np.round(y[1], 9)))
                if key not in els: els[key] = 1; nxt.append(y); allg.append(y)
        front = nxt
    dmin = min(np.linalg.norm(app(g, O3) - O3) for g in allg if not (np.allclose(g[0], np.eye(3)) and np.allclose(g[1], 0)))
    return dmin
worst = 0
for _ in range(25):
    a, b, c = rng.uniform(0.5, 2.0, 3); x0 = rng.uniform(-a, a)
    d = abs(x0 - np.round(x0 / (a / 2)) * (a / 2))
    O3 = np.array([x0, rng.uniform(0, b), rng.uniform(0, 2 * c)])
    worst = max(worst, abs(nearest_brute(a, b, c, O3) - nearest_formula(a, b, c, d)))
cmin = 0.985 / 2                            # in units of D
amin_at_cmin = 2 * np.sqrt(0.985**2 - cmin**2)
# for the best observer d = a/4 the condition is min(a, b, 2c, sqrt(a^2/4 + c^2)) > 0.985 (units of D)
def allowed_fraction(a, c, bnd=0.985):
    if min(a, 2 * c) <= bnd: return 0.0
    need = (bnd**2 - c**2); dmin = 0.0 if need <= 0 else np.sqrt(need) / 2
    return float(max(0.0, 1 - dmin / (a / 4)))
R["U3"] = dict(brute_vs_formula_worst=float(worst), shortest_loop_min_D=cmin, shortest_loop_min_Gpc=cmin * D,
               a_needed_at_shortest_D=float(amin_at_cmin), a_needed_at_shortest_Gpc=float(amin_at_cmin * D),
               detectable_clone_D=1.2, detectable_clone_Gpc=1.2 * D,
               frac_example=dict(a=2.0, c=0.6, fraction=allowed_fraction(2.0, 0.6)))
check("U3 (exact formula, checked by brute force over the group; bound from Vaudrevange et al. 2012 and COMPACT Part I)",
      worst < 1e-9 and allowed_fraction(amin_at_cmin * (1 - 1e-6), cmin * (1 + 1e-9)) == 0.0 and 0 < allowed_fraction(amin_at_cmin * (1 + 1e-3), cmin * (1 + 1e-9)) < 0.01,
      f"nearest copy of us in E7 = min(a, b, 2c, sqrt(4 d^2 + c^2)) with d our distance from the nearest mirror plane (worst mismatch "
      f"against brute force over 25 random universes and positions: {worst:.1e}). The circle bound needs every term above 0.985 D. "
      f"So the shortest loop of a Klein-shaped universe (the glide, length c, through a mirror plane) can be as short as 0.4925 D = "
      f"{cmin*D:.1f} Gpc, about half the last-scattering diameter, if we sit at least sqrt(0.985^2 - c^2)/2 from every mirror plane; at "
      f"that shortest c this needs a >= {amin_at_cmin:.3f} D = {amin_at_cmin*D:.1f} Gpc. Example: a = 2 D, c = 0.6 D leaves "
      f"{100*allowed_fraction(2.0, 0.6):.0f}% of positions allowed. COMPACT Part IIb: such a universe's correlations could beat cosmic "
      f"variance if the nearest copy is under about 1.2 D = {1.2*D:.1f} Gpc, so the window 0.985 D to 1.2 D ({bound:.1f} to {1.2*D:.1f} Gpc) "
      f"is allowed today and testable")

# ------------------------------------------------------------------ U4 handedness
say("== U4. Handedness across an orientation-reversing loop")
s1, s2, s3 = sp.Matrix([[0, 1], [1, 0]]), sp.Matrix([[0, -sp.I], [sp.I, 0]]), sp.Matrix([[1, 0], [0, -1]])
Z, I2 = sp.zeros(2), sp.eye(2)
def blk(a, b, c, d): return sp.Matrix(sp.BlockMatrix([[a, b], [c, d]]))
g0 = blk(Z, I2, I2, Z); gs = [blk(Z, s, -s, Z) for s in (s1, s2, s3)]
gam = [g0] + gs; g5 = sp.I * g0 * gs[0] * gs[1] * gs[2]
eta = sp.diag(1, -1, -1, -1)
cliff = all(sp.simplify(gam[m] * gam[n] + gam[n] * gam[m] - 2 * eta[m, n] * sp.eye(4)) == sp.zeros(4) for m in range(4) for n in range(4))
refl_ok, flips = True, True
for axis in (1, 2, 3):
    S = gam[axis] * g5; Si = S.inv()
    for mu in range(4):
        want = -gam[mu] if mu == axis else gam[mu]
        refl_ok &= sp.simplify(Si * gam[mu] * S - want) == sp.zeros(4)
    flips &= sp.simplify(S * g5 * Si + g5) == sp.zeros(4)
th = sp.symbols('theta', real=True)
rot_keep = True
for (i, j) in ((1, 2), (2, 3), (3, 1)):
    Srot = sp.cos(th / 2) * sp.eye(4) + sp.sin(th / 2) * gam[i] * gam[j]
    rot_keep &= sp.simplify(Srot * g5 - g5 * Srot) == sp.zeros(4)
# the Standard Model's lack of a mirror symmetry: CP violation, Jarlskog invariant from PDG 2024 Wolfenstein parameters
lam, Aw, rb, eb = 0.22501, 0.826, 0.1591, 0.3523
s12, s23 = lam, Aw * lam**2
z = Aw * lam**3 * (rb + 1j * eb) * np.sqrt(1 - Aw**2 * lam**4) / (np.sqrt(1 - lam**2) * (1 - Aw**2 * lam**4 * (rb + 1j * eb)))
s13, dl = abs(z), np.angle(z); c12, c23, c13 = np.sqrt(1 - s12**2), np.sqrt(1 - s23**2), np.sqrt(1 - s13**2)
V = np.array([[c12 * c13, s12 * c13, s13 * np.exp(-1j * dl)],
              [-s12 * c23 - c12 * s23 * s13 * np.exp(1j * dl), c12 * c23 - s12 * s23 * s13 * np.exp(1j * dl), s23 * c13],
              [s12 * s23 - c12 * c23 * s13 * np.exp(1j * dl), -c12 * s23 - s12 * c23 * s13 * np.exp(1j * dl), c23 * c13]])
J = float(np.imag(V[0, 1] * V[1, 2] * np.conj(V[0, 2]) * np.conj(V[1, 1])))
unit_err = float(np.max(np.abs(V @ V.conj().T - np.eye(3))))
R["U4"] = dict(clifford=bool(cliff), reflection_rep=bool(refl_ok), reflection_flips_chirality=bool(flips), rotation_keeps_chirality=bool(rot_keep),
               jarlskog=J, ckm_unitarity_error=unit_err, epsilon_K=2.228e-3)
check("U4 (exact sympy for the spinor algebra; measured CKM numbers, PDG 2024)", cliff and refl_ok and flips and rot_keep and 2.9e-5 < J < 3.3e-5,
      f"for each mirror x_i -> -x_i the spinor map S = gamma_i gamma5 represents it exactly (S^-1 gamma^mu S = Lambda gamma), and "
      f"S gamma5 S^-1 = -gamma5: a left-handed field carried round the glide loop comes back right-handed. Every rotation (the half-turn "
      f"included) commutes with gamma5 and keeps the hand. So a Klein-shaped space needs laws that are unchanged by a mirror (with the "
      f"CRT theorem, the needed symmetry is CR, equivalent to T). The weak force breaks P outright and the Standard Model breaks CP: "
      f"Jarlskog invariant J = {J:.3e} from the PDG Wolfenstein parameters (unitarity error {unit_err:.0e}); kaon epsilon_K = 2.228e-3 "
      f"(first seen 1964)")

# ------------------------------------------------------------------ U5 the spontaneous loophole: a forced mirror wall
say("== U5. If a field chose the hand, a Klein-shaped universe would need a mirror wall")
Mp = 1.22089e19                              # Planck mass, GeV (G = 1/Mp^2)
H0 = 2.13317e-42 * h                         # GeV
rho_c = 3 * H0**2 * Mp**2 / (8 * np.pi)
sig_energy = rho_c / H0                      # one wall across the horizon carries less than the critical density
sig_cmb = 1e-5 * Mp**2 * H0 / (2 * np.pi)    # delta T / T ~ 2 pi G sigma / H0 below 1e-5 (Zel'dovich-Kobzarev-Okun estimate)
v_ew = 246.22
sig_ew = v_ew**3
R["U5"] = dict(rho_c_GeV4=rho_c, sigma_max_energy_GeV3=sig_energy, sigma_max_cmb_GeV3=sig_cmb,
               sigma_max_energy_MeV=1e3 * sig_energy**(1 / 3), sigma_max_cmb_MeV=1e3 * sig_cmb**(1 / 3),
               sigma_weak_GeV3=sig_ew, ratio_weak_over_cmb=sig_ew / sig_cmb)
check("U5 (numbers, order of magnitude; Zel'dovich, Kobzarev and Okun 1974)", sig_ew / sig_cmb > 1e15 and 0.3 < 1e3 * sig_cmb**(1/3) < 2,
      f"going round the glide must turn the field's choice into its mirror, so it has to flip somewhere: a wall. Critical density "
      f"{rho_c:.2e} GeV^4. A horizon-sized wall stays under it only if its tension is below ({1e3*sig_energy**(1/3):.0f} MeV)^3, and keeps "
      f"the CMB smooth to 1e-5 only below ({1e3*sig_cmb**(1/3):.2f} MeV)^3. A hand chosen at the weak scale or above gives tension about "
      f"v^3 = {sig_ew:.1e} GeV^3, {sig_ew/sig_cmb:.0e} times too heavy. The wall could hide only beyond our horizon, and then the twist "
      f"would be out of sight too")

# ------------------------------------------------------------------ U6 twist of the shape versus local torsion
say("== U6. A twist of the shape is not torsion")
hbar, cc, mn, lP = 1.054571817e-34, 2.99792458e8, 1.67492750e-27, 1.616255e-35
lamC = hbar / (mn * cc); rhoC = mn / (lamC * lP**2)
rho_ns = 1e18
R["U6"] = dict(shape_twist_curvature=0, shape_twist_torsion=0, glide_holonomy_det=-1, halfturn_holonomy_det=1,
               cartan_density_neutron_kg_m3=rhoC, neutron_star_core_kg_m3=rho_ns, ratio=rhoC / rho_ns, lab_torsion_bound_GeV=1e-31)
check("U6 (exact for the shape, order of magnitude for Einstein-Cartan; Hehl et al. 1976, Kostelecky, Russell and Tasson 2008)",
      rhoC / rho_ns > 1e35,
      f"the shape's twist: zero curvature and zero torsion at every point (U1); it shows only in the holonomy of a loop that goes all the "
      f"way round the universe (det -1 for the Klein glide, +1 for the half-turn). Local torsion (Einstein-Cartan) is a field set by spin "
      f"density at each point and zero in empty space; it matters only near the Cartan density m/(lambda_C l_P^2) = {rhoC:.1e} kg/m^3 for "
      f"neutrons, {rhoC/rho_ns:.0e} times a neutron star's core; lab tests bound background torsion near 1e-31 GeV. The repo's "
      f"torsion-dial round found the law's torsion derived (k times the bracket on the chains, forced to zero on the sphere sum), so "
      f"neither the law nor the shape supplies torsion in empty space")

# ------------------------------------------------------------------ U7 spinors and phases on the twisted spaces (for the clock question)
say("== U7. What going round a twist does to a spinor's phase")
def lift_ok(Sa, Sb, SG):
    r1 = sp.simplify(SG * Sa * SG.inv() - Sa.inv()) == sp.zeros(4)      # glide reverses the a-slide
    r2 = sp.simplify(SG * Sb - Sb * SG) == sp.zeros(4)                 # glide and b-slide commute
    r3 = sp.simplify(Sa * Sb - Sb * Sa) == sp.zeros(4)
    return r1 and r2 and r3
Rp = gam[1] * g5; Rm = sp.I * gam[1] * g5                              # two lifts of the same mirror, differing by a phase i
same_mirror = all(sp.simplify(Rp.inv() * gam[mu] * Rp - Rm.inv() * gam[mu] * Rm) == sp.zeros(4) for mu in range(4))
counts, forced, eig = {}, {}, {}
for name, Rl in (("Pin+", Rp), ("Pin-", Rm)):
    n = 0
    for ea, eb, eG in itertools.product((1, -1), repeat=3):
        n += lift_ok(ea * sp.eye(4), eb * sp.eye(4), eG * Rl)
    counts[name] = n
    forced[name] = sp.simplify(Rl * Rl)[0, 0]                         # the slide by 2c = glide twice: its sign is forced
    eig[name] = sorted(set(str(sp.nsimplify(e)) for e in Rl.eigenvals()))
Shalf = gam[1] * gam[2]                                               # spinor lift of the half-turn about z (theta = pi)
half_sq = sp.simplify(Shalf * Shalf)[0, 0]; half_eig = sorted(set(str(e) for e in Shalf.eigenvals()))
R["U7"] = dict(same_mirror=bool(same_mirror), structures=counts, glide_twice_sign={k: str(v) for k, v in forced.items()},
               glide_once_eigenvalues=eig, halfturn_twice_sign=str(half_sq), halfturn_once_eigenvalues=half_eig)
check("U7 (exact sympy; Pin structures as in Kirby and Taylor 1990, Berg, DeWitt-Morette, Gwo and Kriele 2001)",
      same_mirror and counts == {"Pin+": 8, "Pin-": 8} and forced["Pin+"] == 1 and forced["Pin-"] == -1 and half_sq == -1,
      f"on Klein space a spinor can't be a plain (Spin) spinor; it must be a Pin spinor, and there are two kinds, set by whether the mirror's "
      f"spinor map squares to +1 (Pin+) or -1 (Pin-). Both maps represent the same mirror; they differ only by a phase i. Each kind has "
      f"{counts['Pin+']} versions on E7 (one free sign for each of the slide a, slide b and the glide; every one satisfies the gluing rules), "
      f"matching H^1(E7; Z2) = (Z2)^3. Going once round the glide multiplies a spinor by the mirror map: eigenvalues {eig['Pin+']} for Pin+ "
      f"(a sign) and {eig['Pin-']} for Pin- (a quarter turn of phase). Going twice round, which is once round the long slide 2c, gives "
      f"{forced['Pin+']} for Pin+ and {forced['Pin-']} for Pin-, forced, with no choice left: a half turn of phase. For contrast the half-turn "
      f"space's screw lifts to gamma1 gamma2, eigenvalues {half_eig}, and twice round gives {half_sq}, the usual sign of a 360-degree turn. "
      f"None of this touches time: the glide moves only space, so cosmic time stays one single-valued clock on a Klein-shaped universe")

say(f"\n{passed} / {passed + failed} pass  ({time.time() - T0:.1f} s)")
R["summary"] = dict(passed=passed, failed=failed)
open("checks_output.txt", "w").write(_log.getvalue())
json.dump(R, open("results.json", "w"), indent=1, default=float)
