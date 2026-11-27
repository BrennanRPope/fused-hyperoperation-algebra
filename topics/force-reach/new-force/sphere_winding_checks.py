"""Part O of the new-force DRAFT (Nov 27 2026): does our law pick one winding on the hidden sphere?
Section 13 (part N) found that a sphere holds itself, but where it rests is set by K, the particle count. The natural way for a
whole number to fix that instead is a winding: the electric circle wraps the sphere n times (a circle bundle over the sphere;
n = 1 is a squashed 3-sphere, n >= 2 a lens space S^3/Z_n). This asks what a winding fixes and whether any n matches the Standard Model.
  O1 the curvature of the wrapped sphere; O2 the winding holds the sphere at a fixed size relative to the circle;
  O3 the sphere's arrows lifted to the wrapped space; O4 the strength ratio that follows; O5 the Standard Model test;
  O6 the fundamental length; O7 no massless fermions; O8 what the circle's hold now has to beat.
Independent of the other scripts (shared inputs from new_force_core.py): writes O_ keys to results_O.json and sphere_winding_output.txt.
Order-of-magnitude estimates are labelled as such."""
import math
import sympy as sp
from new_force_core import ainv, R_lP, MP, write_json_atomic

lines = []; NP = [0, 0]; OUT = {}
def log(s=""):
    print(s); lines.append(s)
def check(name, ok, detail=""):
    NP[0 if ok else 1] += 1
    log(("PASS  " if ok else "FAIL  ") + name + (f"   [{detail}]" if detail else ""))

th, ph, ps = sp.symbols('theta phi psi', real=True)
rho, a, n, c = sp.symbols('rho a n c', positive=True)
X = [th, ph, ps]
A = n/2*sp.cos(th)                                   # monopole of winding n: integral of dA over the sphere = 2 pi n
g = sp.Matrix([[rho**2, 0, 0], [0, rho**2*sp.sin(th)**2 + a**2*A**2, a**2*A], [0, a**2*A, a**2]])   # rho^2 dOmega^2 + a^2 (dpsi + A)^2
gi = sp.simplify(g.inv())

log("== O. Does our law pick one winding on the hidden sphere?")
# ---- O1 curvature
def Gam(i, j, k):
    return sp.simplify(sum(gi[i, l]*(sp.diff(g[l, j], X[k]) + sp.diff(g[l, k], X[j]) - sp.diff(g[j, k], X[l])) for l in range(3))/2)
G = [[[Gam(i, j, k) for k in range(3)] for j in range(3)] for i in range(3)]
def Ric(j, k):
    return sp.simplify(sum(sp.diff(G[i][j][k], X[i]) - sp.diff(G[i][j][i], X[k])
                           + sum(G[i][i][l]*G[l][j][k] - G[i][k][l]*G[l][j][i] for l in range(3)) for i in range(3)))
Rint = sp.simplify(sum(gi[j, k]*Ric(j, k) for j in range(3) for k in range(3)))
r_ = sp.symbols('r', positive=True)
check("O1 wrapping the circle n times around the sphere lowers the hidden space's curvature: R = 2/rho^2 - a^2 n^2/(8 rho^4), where the second term is the circle's own field wrapped on the sphere. Check: n = 1 with rho = a/2 is the round 3-sphere, R = 6/a^2",
      sp.simplify(Rint - (2/rho**2 - a**2*n**2/(8*rho**4))) == 0 and sp.simplify(Rint.subs({n: 1, rho: r_/2, a: r_}) - 6/r_**2) == 0,
      f"R = {sp.factor(Rint)}")
# ---- O2 the hold: in the 4D frame the potential is -R / volume, volume ~ rho^2 a (same normalisation as part N)
V = -Rint/(rho**2*a)
crit = sp.solve(sp.diff(V, rho), rho)
rho_star = n*a*sp.sqrt(sp.Rational(3, 32))
V2 = sp.simplify(sp.diff(V, rho, 2).subs(rho, rho_star))
OUT["O_rho_over_a_per_n"] = float(sp.sqrt(sp.Rational(3, 32)))
check("O2 the winding holds the sphere by itself: with the circle held, V = -2/(a rho^4) + a n^2/(8 rho^6) has one stable minimum at rho = n a sqrt(3/32) = 0.306 n a. This is classical, so the particle count drops out of the sphere's size relative to the circle",
      len(crit) == 1 and sp.simplify(crit[0] - rho_star) == 0 and sp.simplify(V2*rho_star**6*a - 8*2) == 0,
      f"rho*/a = {float(sp.sqrt(sp.Rational(3,32))):.4f} n; V'' = 16/(a rho*^6) > 0")
# ---- O2b how big a count it takes to matter (rough, part N's units and per-state K)
K_dof = [1.5772610071418304e-07, 2.018894089141543e-05]           # part N, N_K_per_dof (period 2 pi rho, pi rho)
R1 = R_lP
def l7(rh, aa=R1): return (8*math.pi**2*aa*rh**2)**0.2
def count_ratio(nv, Nst, Kd):                                       # Casimir term over the wrapped-field term at the rest point
    rh = nv*R1*math.sqrt(3/32); u = (16*math.pi)**0.2*l7(rh)
    return 8*Nst*Kd/((R1/u)**2*nv**2*(rh/u)**3)
cr = {N_: [count_ratio(2, N_, k) for k in K_dof] for N_ in (10, 100, 1000)}
OUT["O_count_ratio_n2"] = cr
check("O2b rough: for n = 2 the vacuum-energy term is under 1% of the winding term at 10 states, about 0.04% to 5% at 100, and 0.4% to 50% at 1,000, so the count only matters if it runs to hundreds or more",
      cr[10][1] < 0.01 and cr[100][1] < 0.06 and 0.2 < cr[1000][1] < 0.8, ", ".join(f"{k} states: {v[0]:.1e} to {v[1]:.1e}" for k, v in cr.items()))
# ---- O3 the sphere's arrows lifted to the wrapped space
Kx = [-sp.sin(ph), -sp.cot(th)*sp.cos(ph), n/2*sp.cos(ph)/sp.sin(th)]
Ky = [sp.cos(ph), -sp.cot(th)*sp.sin(ph), n/2*sp.sin(ph)/sp.sin(th)]
Kz = [0, 1, -c]
def lie(K):
    return sp.Matrix(3, 3, lambda i, j: sp.simplify(sum(K[l]*sp.diff(g[i, j], X[l]) + g[l, j]*sp.diff(K[l], X[i]) + g[i, l]*sp.diff(K[l], X[j]) for l in range(3))))
def comm(U, W): return [sp.simplify(sum(U[l]*sp.diff(W[i], X[l]) - W[l]*sp.diff(U[i], X[l]) for l in range(3))) for i in range(3)]
def avg(K):
    nrm = sp.simplify(sum(g[i, j]*K[i]*K[j] for i in range(3) for j in range(3)))
    return sp.simplify(sp.integrate(sp.integrate(nrm*sp.sin(th), (th, 0, sp.pi)), (ph, 0, 2*sp.pi))/(4*sp.pi))
killing = all(lie(K) == sp.zeros(3) for K in (Kx, Ky, Kz))
xy = comm(Kx, Ky)                                                    # must be -Kz for su(2): fixes c = 0
cf = sp.Symbol('cf')
csol = sp.solve(sp.Eq(xy[2], cf), cf)
Kz0 = [0, 1, 0]
avgs = [avg(Kx), avg(Ky), avg(Kz0)]
target = sp.Rational(2, 3)*rho**2 + a**2*n**2/12
check("O3 the sphere's three rotation arrows still work on the wrapped space once each picks up a push along the circle; keeping their algebra (su(2)) fixes that push, and their average square length grows from 2 rho^2/3 to 2 rho^2/3 + a^2 n^2/12. Together with the circle's own arrow the force group is SU(2) x U(1), as known for a squashed 3-sphere",
      killing and sp.simplify(xy[1] + 1) == 0 and xy[0] == 0 and csol == [0] and all(sp.simplify(x - target) == 0 for x in avgs),
      f"<|K|^2> = {sp.factor(avgs[0])}")
# ---- O4 the strength ratio at the winding's rest point: alpha ~ 1/<|K|^2> (part S), circle arrow |d/dpsi|^2 = a^2
ratio = sp.simplify(a**2/target.subs(rho, rho_star))
rat_n = {k: float(ratio.subs(n, k)) for k in (1, 2, 3, 4)}
OUT.update(O_ratio_formula="48/(7 n^2)", O_ratio_by_n=rat_n)
check("O4 so the winding fixes the ratio of the two strengths with no particle count and no size in it: alpha_SU(2)/alpha_circle = 48/(7 n^2), which is 6.86, 1.71, 0.76, 0.43 for n = 1 to 4",
      sp.simplify(ratio - sp.Rational(48, 7)/n**2) == 0, ", ".join(f"n={k}: {v:.2f}" for k, v in rat_n.items()))
# ---- O5 Standard Model at the circle's copy scale (one loop, as in parts M and N)
MZ = 91.1876; ainvZ = 127.951; s2w = 0.23122; a3Z = 0.1180
inv0 = [0.6*(1 - s2w)*ainvZ, s2w*ainvZ, 1/a3Z]; b_sm = [41/10, -19/6, -7]
Lg = math.log((MP/R1)/MZ)/(2*math.pi); inv_KK = [i - b*Lg for i, b in zip(inv0, b_sm)]
inv1, inv2 = inv_KK[0], inv_KK[1]; invY = inv1*5/3; invEM = invY + inv2
units = {"hypercharge, unit 1/6 (quark doublet)": 36*invY, "hypercharge, unit 1/2 (lepton doublet)": 4*invY,
         "electric, unit e/3 (down quark)": 9*invEM, "electric, unit e (electron, part G)": invEM}
need = {k: math.sqrt(48/(7*(v/inv2))) for k, v in units.items()}      # alpha2/alpha_circle = (1/alpha_circle)/(1/alpha2)
OUT.update(O_SM_inv_Y_2_EM=[invY, inv2, invEM], O_n_needed=need)
best = {k: (round(v) if round(v) > 0 else 1) for k, v in need.items()}
miss = {k: (48/(7*best[k]**2))/(units[k]/inv2) - 1 for k in units}
OUT["O_best_integer_miss"] = {k: [best[k], miss[k]] for k in units}
check("O5 no winding matches the Standard Model: reading the circle's charge unit four simple ways, the measured ratio at 5.2e17 GeV needs n = 0.39, 1.18, 0.58 or 1.75, none of them a whole number. The nearest whole windings miss the measured ratio by 23% to 85% (closest: n = 2 with unit e, 23%; n = 1 with unit 1/2 hypercharge, 39%), far beyond one-loop uncertainty",
      min(abs(m) for m in miss.values()) > 0.2 and abs(need["hypercharge, unit 1/2 (lepton doublet)"] - 1.18) < 0.01 and abs(need["electric, unit e (electron, part G)"] - 1.75) < 0.01,
      "; ".join(f"{k}: n = {need[k]:.2f}, n={best[k]} misses by {100*miss[k]:+.0f}%" for k in units))
# ---- O6 the fundamental length (l7^5 = l_P^2 x volume, volume = 4 pi rho^2 x 2 pi a)
fl = {k: [k*R1*math.sqrt(3/32)/l7(k*R1*math.sqrt(3/32)), R1/l7(k*R1*math.sqrt(3/32))] for k in (1, 2, 3)}
OUT["O_rho_and_a_over_l7"] = fl
check("O6 with the electric circle's size from part G, n = 1 puts the sphere below the fundamental length (0.72), where none of this can be trusted; n = 2 and 3 put it at 1.10 and 1.40",
      fl[1][0] < 1 < fl[2][0] < fl[3][0] and 0.65 < fl[1][0] < 0.8 and 1.05 < fl[2][0] < 1.15, ", ".join(f"n={k}: rho/l7 = {v[0]:.2f}, a/l7 = {v[1]:.2f}" for k, v in fl.items()))
# ---- O7 no massless fermions: curvature stays positive at the rest point, so Lichnerowicz rules out zero modes
R_star = sp.simplify(Rint.subs(rho, rho_star))
check("O7 the winding gives no light matter: at the rest point the hidden space's curvature is still positive, R = 2/(3 rho^2), so by Lichnerowicz's theorem it has no massless fermions at all; the states the wrapping counts (|m n| at each circle level) all sit at the copy scale. This is the known chirality problem for forces made only from the hidden space's shape (Witten 1983)",
      sp.simplify(R_star - sp.Rational(2, 3)/rho_star**2) == 0, f"R* = {sp.simplify(R_star*rho_star**2)}/rho*^2")
# ---- O8 what the circle's hold has to beat
Vstar = sp.simplify(V.subs(rho, rho_star))
expo = sp.simplify(sp.diff(sp.log(-Vstar), a)*a)
check("O8 with the sphere at rest the potential goes as -1/a^5, so it pulls the circle smaller; the circle's own hold (section 7) must now beat that pull, which this check does not do",
      expo == -5 and sp.simplify(Vstar*a**5*n**4 + sp.Rational(2048, 27)) == 0, f"V* = {sp.simplify(Vstar*a**5*n**4)}/(a^5 n^4)")

log(f"\nPART O TOTAL: {NP[0]} pass, {NP[1]} fail")
write_json_atomic("results_O.json", OUT)
open("sphere_winding_output.txt", "w").write("\n".join(lines) + "\n")
