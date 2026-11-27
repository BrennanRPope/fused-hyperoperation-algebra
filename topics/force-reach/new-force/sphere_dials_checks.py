"""Part P of the new-force DRAFT (Nov 27 2026): more than one dial on the wound sphere.
Section 14 (part O) wrapped one circle n times around the hidden sphere: the strength ratio came out 48/(7n^2), and no whole n
matched. Brennan: "Try more than one dial". Here several circles (a torus T^k with any shape h) wrap the sphere, each with its own
winding n_i, so the dials are the windings, how many circles, the angle between them, and which charge one step carries.
  P1 curvature for a 2-torus of any shape; P2 the lifted arrows; P3 the rest point and the strength ratio 48/(7N), N = winding length^2;
  P4 every simple charge reading against whole N; P5 which windings give N = 3; P6 the fundamental length;
  P7 how often a near-hit happens by chance (rough); P8 a continuous dial fits anything.
Independent of the other scripts (shared inputs from new_force_core.py): writes P_ keys to results_P.json and sphere_dials_output.txt."""
import math, itertools
import sympy as sp
from new_force_core import R_lP, MP, write_json_atomic

lines = []; NP = [0, 0]; OUT = {}
def log(s=""):
    print(s); lines.append(s)
def check(name, ok, detail=""):
    NP[0 if ok else 1] += 1
    log(("PASS  " if ok else "FAIL  ") + name + (f"   [{detail}]" if detail else ""))

th, ph, p1, p2 = sp.symbols('theta phi psi1 psi2', real=True)
rho = sp.symbols('rho', positive=True)
n1, n2 = sp.symbols('n1 n2', real=True)
h11, h22 = sp.symbols('h11 h22', positive=True); h12 = sp.symbols('h12', real=True)
X = [th, ph, p1, p2]
A = sp.cos(th)/2                                     # unit monopole: each circle i shifts by n_i A
H = sp.Matrix([[h11, h12], [h12, h22]])               # torus shape: circle lengths and angle
# metric rho^2 dOmega^2 + h_ij (dpsi_i + n_i A dphi)(dpsi_j + n_j A dphi)
E = sp.Matrix([[0, n1*A, 1, 0], [0, n2*A, 0, 1]])     # rows: dpsi_i + n_i A dphi, in (theta, phi, psi1, psi2)
g = sp.diag(rho**2, rho**2*sp.sin(th)**2, 0, 0) + E.T*H*E
g = sp.simplify(g)
gi = sp.simplify(g.inv())
log("== P. More than one dial on the wound sphere")
def Gam(i, j, k):
    return sp.simplify(sum(gi[i, l]*(sp.diff(g[l, j], X[k]) + sp.diff(g[l, k], X[j]) - sp.diff(g[j, k], X[l])) for l in range(4))/2)
G = [[[Gam(i, j, k) for k in range(4)] for j in range(4)] for i in range(4)]
def Ric(j, k):
    return sp.simplify(sum(sp.diff(G[i][j][k], X[i]) - sp.diff(G[i][j][i], X[k])
                           + sum(G[i][i][l]*G[l][j][k] - G[i][k][l]*G[l][j][i] for l in range(4)) for i in range(4)))
Rint = sp.simplify(sum(gi[j, k]*Ric(j, k) for j in range(4) for k in range(4)))
W = (sp.Matrix([[n1, n2]])*H*sp.Matrix([n1, n2]))[0]
import random
random.seed(7)                                       # sympy can't close the trig identity symbolically; test it at 12 random points
PTS = [{rho: random.uniform(0.5, 3), n1: random.uniform(-3, 3), n2: random.uniform(-3, 3), h11: random.uniform(0.5, 3),
        h22: random.uniform(0.5, 3), h12: random.uniform(-0.4, 0.4), th: random.uniform(0.2, 2.9)} for _ in range(12)]
check("P1 with two circles of any lengths and angle, each wound n_i times, the hidden space's curvature is R = 2/rho^2 - W/(8 rho^4), where W = h(n, n) is the squared length of the winding direction in the torus. Only that one length enters, whatever the number of dials",
      all(abs(float((Rint - (2/rho**2 - W/(8*rho**4))).subs(pt))) < 1e-12 for pt in PTS), f"R = 2/rho^2 - ({sp.expand(W)})/(8 rho^4)")
# ---- P2 lifted arrows
Kx = [-sp.sin(ph), -sp.cot(th)*sp.cos(ph), n1/2*sp.cos(ph)/sp.sin(th), n2/2*sp.cos(ph)/sp.sin(th)]
Ky = [sp.cos(ph), -sp.cot(th)*sp.sin(ph), n1/2*sp.sin(ph)/sp.sin(th), n2/2*sp.sin(ph)/sp.sin(th)]
Kz = [0, 1, 0, 0]
def lie(K):
    return sp.Matrix(4, 4, lambda i, j: sp.simplify(sum(K[l]*sp.diff(g[i, j], X[l]) + g[l, j]*sp.diff(K[l], X[i]) + g[i, l]*sp.diff(K[l], X[j]) for l in range(4))))
def comm(U, V2): return [sp.simplify(sum(U[l]*sp.diff(V2[i], X[l]) - V2[l]*sp.diff(U[i], X[l]) for l in range(4))) for i in range(4)]
def avg(K):
    nrm = sp.simplify(sum(g[i, j]*K[i]*K[j] for i in range(4) for j in range(4)))
    return sp.simplify(sp.integrate(sp.integrate(nrm*sp.sin(th), (th, 0, sp.pi)), (ph, 0, 2*sp.pi))/(4*sp.pi))
ok2 = all(lie(K) == sp.zeros(4) for K in (Kx, Ky, Kz)) and [sp.simplify(x) for x in comm(Kx, Ky)] == [0, -1, 0, 0]
avgs = [avg(K) for K in (Kx, Ky, Kz)]
check("P2 the sphere's three arrows lift to the wrapped space and keep their algebra; their average square length is 2 rho^2/3 + W/12, the same form as one circle with a^2 n^2 replaced by W",
      ok2 and all(sp.simplify(x - (sp.Rational(2, 3)*rho**2 + W/12)) == 0 for x in avgs), f"<|K|^2> = 2 rho^2/3 + W/12")
# ---- P3 rest point and ratio
Wv, v2 = sp.symbols('W v2', positive=True)
V = -(2/rho**2 - Wv/(8*rho**4))/rho**2                # times 1/(torus volume), which doesn't depend on rho
crit = sp.solve(sp.diff(V, rho), rho)
rstar = sp.sqrt(3*Wv/32)
ratio = sp.simplify(v2/(sp.Rational(2, 3)*rstar**2 + Wv/12))
check("P3 so the sphere rests at rho^2 = 3W/32 and the sphere force over any circle force (a circle of squared length |v|^2) is 48|v|^2/(7W) = 48/(7N), with N = W/|v|^2 the winding direction's squared length in units of that circle. Several windings only change N: three equal circles wound (1,1,1) act like one circle along the cube's body diagonal, N = 3",
      len(crit) == 1 and sp.simplify(crit[0] - rstar) == 0 and sp.simplify(ratio - sp.Rational(48, 7)*v2/Wv) == 0
      and (sp.Matrix([[1, 1, 1]])*sp.eye(3)*sp.Matrix([1, 1, 1]))[0] == 3, "ratio = 48/(7N)")
# ---- P4 every simple charge reading
MZ = 91.1876; ainvZ = 127.951; s2w = 0.23122; a3Z = 0.1180
inv0 = [0.6*(1 - s2w)*ainvZ, s2w*ainvZ, 1/a3Z]; b_sm = [41/10, -19/6, -7]
R1 = R_lP
Lg = math.log((MP/R1)/MZ)/(2*math.pi); inv_KK = [i - b*Lg for i, b in zip(inv0, b_sm)]
invY = inv_KK[0]*5/3; inv2 = inv_KK[1]; invEM = invY + inv2
readings = {"hypercharge 1/6": 36*invY, "hypercharge 1/3": 9*invY, "hypercharge 1/2": 4*invY, "hypercharge 2/3": 9/4*invY,
            "hypercharge 1": invY, "electric e/3": 9*invEM, "electric e": invEM}
Nneed = {k: 48/(7*(v/inv2)) for k, v in readings.items()}
def best(Nv):
    cands = [m for m in (math.floor(Nv), math.ceil(Nv)) if m >= 1] or [1]
    m = min(cands, key=lambda m: abs(Nv/m - 1)); return m, Nv/m - 1      # ratio miss = (48/7m)/(measured) - 1 = Nneed/m - 1
bests = {k: best(v) for k, v in Nneed.items()}
OUT.update(P_N_needed=Nneed, P_best_whole_N=bests)
near = [k for k, (m, e) in bests.items() if abs(e) < 0.05]
check("P4 against seven simple charge readings, the measured ratio at 5.2e17 GeV needs N = 0.15, 0.62, 1.39, 2.48, 5.57 (hypercharge units 1/6 to 1) and 0.34, 3.07 (electric units e/3, e). N must be a whole number of at least 1 for windings on a square torus, so only the electric-e reading comes close: N = 3, which misses the ratio by 2.5%. Every hypercharge reading misses by 7% or more (closest: unit 1, N = 5.57 against 6)",
      near == ["electric e"] and abs(bests["electric e"][0] - 3) == 0 and 0.02 < abs(bests["electric e"][1]) < 0.03
      and all(abs(e) >= 0.07 for k, (m, e) in bests.items() if k.startswith("hyper")),
      "; ".join(f"{k}: N = {v:.2f}, nearest {bests[k][0]} misses by {100*bests[k][1]:+.1f}%" for k, v in Nneed.items()))
# ---- P5 which windings give N = 3
sq2 = {a*a + b*b for a in range(4) for b in range(4)}
hexs = {a*a + a*b + b*b for a in range(-3, 4) for b in range(-3, 4)}
cube = [(a, b, c) for a, b, c in itertools.product(range(-2, 3), repeat=3) if a*a + b*b + c*c == 3]
check("P5 N = 3 cannot come from one circle or from two square circles (3 is not a sum of two squares). It needs three equal square circles wound once each, along a cube's body diagonal, or two equal circles at 60 degrees (a hexagonal torus) wound along its long diagonal",
      3 not in sq2 and 3 in hexs and (1, 1, 1) in cube, f"two square circles reach N in {sorted(x for x in sq2 if 0 < x < 10)}; hexagonal reaches {sorted(x for x in hexs if 0 < x < 10)}")
# ---- P6 the fundamental length: l_D^(D-2) = l_P^2 x volume, D = 4 + k + 2, volume = (2 pi a)^k sqrt(det h) 4 pi rho^2
a = R1; rho3 = math.sqrt(3*3/32)*a
def lD(k, detfac):
    vol = (2*math.pi*a)**k*detfac*4*math.pi*rho3**2; return vol**(1/(k + 4))
cube_l = lD(3, 1.0); hex_l = lD(2, math.sqrt(3)/2)
fl = {"cube (9D)": [rho3/cube_l, a/cube_l], "hexagonal (8D)": [rho3/hex_l, a/hex_l]}
OUT["P_N3_rho_and_a_over_lD"] = fl
check("P6 but both ways of getting N = 3 put the sphere below the fundamental length: with every circle at the electric circle's size (part G), the sphere rests at 0.53 circle sizes, which is 0.50 fundamental lengths for the cube and 0.68 for the hexagonal torus, where none of this can be trusted",
      fl["cube (9D)"][0] < 0.55 and 0.6 < fl["hexagonal (8D)"][0] < 0.75 and abs(rho3/a - math.sqrt(9/32)) < 1e-12,
      ", ".join(f"{k}: rho/lD = {v[0]:.2f}, a/lD = {v[1]:.2f}" for k, v in fl.items()))
# ---- P7 rough chance of a near-hit: each reading with N >= 1 lands within 2.5% of a whole number with probability ~ 2 x 0.025 x N
pr = [min(1, 0.05*v) for v in Nneed.values() if v >= 1]
p_any = 1 - math.prod(1 - x for x in pr)
OUT["P_chance_of_near_hit"] = p_any
check("P7 rough: with seven readings and whole N to choose from, a 2.5% near-hit somewhere happens by chance about half the time, so N = 3 for the electric reading is not evidence by itself",
      0.35 < p_any < 0.65, f"chance of at least one hit within 2.5%: {p_any:.2f}")
# ---- P8 a continuous dial fits anything
fit = math.sqrt(Nneed["hypercharge 1/2"])
check("P8 a continuous dial uses up the measurement: let the wound circle be longer than the measured one by any factor and the ratio fits exactly (for the unit-1/2 hypercharge reading, a factor 1.18). So only whole-number dials make a test, and they fail it",
      abs(48/(7*fit**2) - readings["hypercharge 1/2"]/inv2) < 1e-12, f"length factor {fit:.3f}")

log(f"\nPART P TOTAL: {NP[0]} pass, {NP[1]} fail")
write_json_atomic("results_P.json", OUT)
open("sphere_dials_output.txt", "w").write("\n".join(lines) + "\n")
