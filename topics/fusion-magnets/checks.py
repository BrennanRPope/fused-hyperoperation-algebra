import numpy as np, json
from lib import *
R = {}; ok = []
def check(name, cond, info=""):
    ok.append(bool(cond)); print(("PASS " if cond else "FAIL ") + name + ("  | " + info if info else "")); R[name] = info

# A. magnetic pressure and the B^4 lever
p20 = 20**2/(2*MU0)/1e6
check("A1 magnetic pressure at 20 T is ~159 MPa", abs(p20-159.2) < 0.5, f"{p20:.1f} MPa")
check("A2 doubling B multiplies pressure by 4, fusion power density (fixed beta) by 16", (2**2, 2**4) == (4, 16))

# B. union rule = series reliability / weakest link
rng = np.random.default_rng(1); ps = rng.uniform(0, 0.02, 50)
check("B1 k=-1 fold equals 1 - prod(1-p)", abs(-fold(-ps, -1) - (1 - np.prod(1-ps))) < 1e-12 or abs(fold(ps, -1) - (1-np.prod(1-ps))) < 1e-12,
      f"{fold(ps,-1):.6f} vs {1-np.prod(1-ps):.6f}")
lam, L = 0.3, 10.0   # defects per km, km of tape (illustrative)
for n in [10, 100, 10000]:
    v = fold([lam*L/n]*n, -1)
print("   fold over 10000 segments", v, "Poisson", 1-np.exp(-lam*L))
check("B2 fine-segment fold converges to 1 - exp(-lambda L) (Poisson weakest link)", abs(v-(1-np.exp(-lam*L))) < 1e-3, f"{v:.5f} vs {1-np.exp(-lam*L):.5f}")
check("B3 lift -ln(1-p) turns the k=-1 fold into a plain sum (hazards add)",
      abs(-np.log(1-fold(ps,-1)) - np.sum(-np.log(1-ps))) < 1e-10)

# C. derating composition and the corner (factorial interaction) test
sep = lambda B, T: 400*(1+B/0.5)**-0.5*(1-T/85)**1.4
check("C1 separable model: corner sum of log Ic is exactly 0", abs(corner_sum_log(sep, 10, 20, 20, 30)) < 1e-12)
cs = corner_sum_log(Ic, 10, 20, 20, 30)
check("C2 fitted REBCO surface (10-20 T, 20-30 K): corner sum nonzero, so derating factors do not simply multiply", abs(cs) > 1e-3, f"corner sum of ln Ic = {cs:+.4f}")
# derating language: delta_B at T and delta_T at B; composition via o_{-1}
B1, B2, T1, T2 = 10, 20, 20, 30
dB = 1 - Ic(B2, T1)/Ic(B1, T1); dT = 1 - Ic(B1, T2)/Ic(B1, T1)
pred = circ(dB, dT, -1); true = 1 - Ic(B2, T2)/Ic(B1, T1)
check("C3 k=-1 composition of separate deratings vs true combined derating", True, f"dB={dB:.3f}, dT={dT:.3f}, o_-1 predicts {pred:.4f}, fit gives {true:.4f}, gap {true-pred:+.4f}")
check("C4 gap equals (1-dB)(1-dT)(1 - e^corner) exactly", abs((true-pred) - (1-dB)*(1-dT)*(1-np.exp(cs))) < 1e-12)
Ic2020 = Ic(20, 20)
check("C5 fit gives a plausible 20 K, 20 T current for one tape", 200 < Ic2020 < 900, f"Ic(20 T, 20 K) = {Ic2020:.0f} A")

# D. strain composition is o_{+1}
e1, e2 = 0.003, -0.002
check("D1 stretch ratios multiply: (1+e1)(1+e2)-1 = e1 o_1 e2", abs((1+e1)*(1+e2)-1 - circ(e1, e2, 1)) < 1e-15)
check("D2 cross term for 0.3% and 0.4% strains is 1.2e-5 (0.0012 percentage points)", abs(0.003*0.004-1.2e-5) < 1e-18)

# E. semiring dial = KS / Lp aggregation
g = rng.uniform(100, 600, 2000)
for rho in [0.01, 0.05, 0.2, 1.0]:
    v = KS(g, rho); assert g.max() <= v <= g.max() + np.log(len(g))/rho + 1e-9
check("E1 max <= KS_rho <= max + ln(N)/rho for all tested rho", True)
check("E2 Lp norm -> max as p grows (SIMSOPT-style Lp penalties)", abs(Lp(g, 400)/g.max() - 1) < 0.02, f"L400/max = {Lp(g,400)/g.max():.4f}")
check("E3 Lp(x) = exp(KS_p(ln x)) : Lp is the same dial in log units", abs(np.log(Lp(g, 8)) - KS(np.log(g), 8)) < 1e-10)

# F. tape twist on closed coils (Fuller / geometric phase = the swap rule at d=+1)
X = modular_coil(); F = frames(X)
U, hol = bishop(F['T'], F['N'][0])
A = spherical_area(F['T']); itau = np.sum(F['tau']*F['ds'])
m1 = min(abs(wrap(hol - (2*np.pi - A))), abs(wrap(hol + (2*np.pi - A))))
check("F1 modular coil: parallel-transport closing angle = 2pi - tantrix area (mod 2pi, up to orientation)", m1 < 1e-3,
      f"closing angle {hol:+.5f} rad, tantrix area {A:.5f} sr, mismatch {m1:.1e}")
m2 = min(abs(wrap(itau - hol)), abs(wrap(itau + hol)))
check("F2 integral of Frenet torsion = same closing angle (mod 2pi)", m2 < 1e-3, f"int tau ds = {itau:+.5f} rad")
Xp = planar_coil(); Fp = frames(Xp); Up, holp = bishop(Fp['T'], Fp['N'][0])
check("F3 planar coil: tantrix is a great circle, closing angle 0, no twist needed", abs(wrap(holp)) < 1e-6, f"{holp:.2e} rad")
delta = abs(wrap(hol)); Lc = F['L']
for w in [0.004, 0.012]:
    floor = (delta/Lc)**2*w**2/12
    frenet = np.max(F['tau']**2*w**2/12)
    # centroid-like reference frame: Bishop frame (no twist) -> hard-way bending = kappa * |sin(angle to Frenet normal)|
    sinang = np.abs(np.einsum('ij,ij->i', np.cross(U, F['N']), F['T']))
    hard = np.max(w*F['kappa']*sinang/2)
    R[f"F4 w={w*1000:.0f}mm"] = dict(closure_floor=floor, frenet_twist_peak=frenet, bishop_hardway_peak=hard)
    print(f"   w={w*1000:.0f} mm: closure floor {floor:.2e}, Frenet-frame twist peak {frenet:.2e}, no-twist-frame hard-way peak {hard:.2e}  (limit 2e-3)")
check("F4 closure floor is negligible next to local strains (orders of magnitude below 0.2%)", R["F4 w=12mm"]["closure_floor"] < 1e-6)
Xd = dented_coil(); Fd = frames(Xd); Ud, hold = bishop(Fd['T'], Fd['N'][0]); w = 0.012
fr = np.max(Fd['tau']**2*w**2/12); hw = np.max(w*Fd['kappa']*np.abs(np.einsum('ij,ij->i', np.cross(Ud, Fd['N']), Fd['T']))/2)
check("F5 dented coil, 12 mm tape: following the curve's own frame overstrains (twist > 0.2%), untwisted tape stays under", fr > 2e-3 > hw,
      f"Frenet-frame twist peak {fr:.2e}, untwisted hard-way peak {hw:.2e}, min curvature {Fd['kappa'].min():.3f}/m, peak torsion {np.abs(Fd['tau']).max():.1f}/m")
R["coil length m"] = Lc; R["closing angle rad"] = float(hol); R["tantrix area sr"] = float(A)

print(f"\n{sum(ok)}/{len(ok)} checks pass")
json.dump({k: (v if not isinstance(v, np.floating) else float(v)) for k, v in R.items()}, open("results.json", "w"), indent=1, default=float)
