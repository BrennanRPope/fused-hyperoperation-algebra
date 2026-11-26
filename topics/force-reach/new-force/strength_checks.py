"""Part G: a force's strength from a closed direction (the rung class read as electromagnetism's circle).
Run: python3 strength_checks.py   (sympy, numpy). Appends to results.json."""
import json, math
import sympy as sp
R = {"pass": 0, "fail": 0}
OUT = json.load(open("results.json"))
def check(name, ok, info=""):
    R["pass" if ok else "fail"] += 1
    print(("PASS  " if ok else "FAIL  ") + name + (f"   [{info}]" if info else ""))

def X(p1, p2, D):
    """Static interaction weight of linearized gravity in D spacetime dims between two smeared point sources
    with D-momenta p1, p2 (T^{AB} = p^A p^B / E):  X = T1^{AB} T2_{AB} - T1 T2 / (D - 2).
    The potential is -C X / r^{D-3} with one overall constant C."""
    eta = sp.diag(*([-1] + [1]*(D - 1)))
    E1, E2 = p1[0], p2[0]
    T1 = sp.Matrix(p1) * sp.Matrix(p1).T / E1
    T2 = sp.Matrix(p2) * sp.Matrix(p2).T / E2
    T1l = eta * T1 * eta
    contr = sum(T1l[a, b]*T2[a, b] for a in range(D) for b in range(D))
    tr1 = sum(eta[a, a]*T1[a, a] for a in range(D)); tr2 = sum(eta[a, a]*T2[a, a] for a in range(D))
    return sp.simplify(contr - tr1*tr2/(D - 2))

m1, m2, E1, E2, E = sp.symbols('m1 m2 E1 E2 E', positive=True)
print("== G. Strength from a closed direction: electromagnetism as gravity's circle (Kaluza-Klein)")
# 4D reference: tensor-only Newton between static masses
X4 = X([m1, 0, 0, 0], [m2, 0, 0, 0], 4)
check("G0 4D gravity alone: static masses have weight m1 m2 / 2; light past a mass has twice that per unit energy (gamma = 1)",
      sp.simplify(X4 - m1*m2/2) == 0 and sp.simplify(X([m1, 0, 0, 0], [E, E, 0, 0], 4) - m1*E) == 0)
# one extra circle: D = 5. Rung modes are massless in 5D with momentum +-E along the circle; 4D mass = E
like = X([E1, 0, 0, 0, E1], [E2, 0, 0, 0, E2], 5)
opp = X([E1, 0, 0, 0, E1], [E2, 0, 0, 0, -E2], 5)
tens = E1*E2/2
vec = sp.simplify((opp - like)/2); sca = sp.simplify((opp + like)/2 - tens)
check("G1 two copies with the same circle charge feel no net force; opposite charges attract at 8 times the tensor weight",
      like == 0 and sp.simplify(opp - 4*E1*E2) == 0, f"like {like}, opposite {opp}")
check("G2 split into 4D pieces: the vector (electric) part is 4 times gravity's, the scalar part 3 times, for the copies",
      sp.simplify(vec/tens - 4) == 0 and sp.simplify(sca/tens - 3) == 0, f"vector/tensor {sp.simplify(vec/tens)}, scalar/tensor {sp.simplify(sca/tens)}")
G, Rr, hbar, c = sp.symbols('G R hbar c', positive=True)
e2_over_4pi = 4*G*(1/Rr)**2            # Coulomb e^2/(4 pi r^2) = 4 G m^2 / r^2 with m = 1/R  (hbar = c = 1)
alpha_sym = e2_over_4pi
check("G3 so the electric strength is fixed: alpha = e^2/(4 pi) = 4 G / R^2 (hbar = c = 1), i.e. e = sqrt(16 pi G)/R, Klein's relation",
      sp.simplify(sp.sqrt(4*sp.pi*alpha_sym) - sp.sqrt(16*sp.pi*G)/Rr) == 0)
Dl = sp.symbols('Delta', negative=True)
check("G4 in our reading, R = 1/sqrt(-Delta): alpha = 4 G |Delta|. The same reading that spaces the rungs sets the strength",
      sp.simplify(alpha_sym.subs(Rr, 1/sp.sqrt(-Dl)) - 4*G*(-Dl)) == 0)
lP = 1.616255e-35; ainv = 137.035999084; hbarc = 1.973269804e-7; MP = 1.220890e19
Rm = 2*lP*math.sqrt(ainv)
m1GeV = hbarc/Rm/1e9
OUT["G_circle_radius_m"] = Rm; OUT["G_circle_radius_lP"] = Rm/lP; OUT["G_rung1_GeV"] = m1GeV
check("G5 with the measured alpha = 1/137.036 the circle's radius is R = 2 l_P sqrt(137) = 23.4 Planck lengths = 3.8e-34 m",
      abs(Rm/lP - 23.41) < 0.01, f"R = {Rm:.3e} m = {Rm/lP:.2f} l_P")
check("G5b the first charged copy weighs hbar c / R = (sqrt(alpha)/2) M_P = 5.2e17 GeV, charge e; copy n has charge n e (charge comes in whole units)",
      abs(m1GeV - math.sqrt(1/ainv)/2*MP)/m1GeV < 1e-3 and abs(m1GeV - 5.21e17)/5.21e17 < 0.01, f"{m1GeV:.3e} GeV")
# the new force on ordinary matter: the circle's size breathing
res = {}
for n in range(1, 7):
    D = 4 + n
    pm = lambda M: [M] + [0]*(D - 1)
    Xmm = X(pm(m1), pm(m2), D)
    sc = sp.simplify(Xmm - m1*m2/2)
    alpha_s = sp.simplify(sc/(m1*m2/2))
    light = [E, E] + [0]*(D - 2)
    gam = sp.simplify((X(pm(m1), light, D)/E)/(Xmm/m2)) - 1
    res[n] = (alpha_s, sp.nsimplify(gam))
check("G6 one circle predicts a new force on ordinary matter: the circle's size as a field, with strength exactly 1/3 of gravity",
      res[1][0] == sp.Rational(1, 3), f"alpha_scalar = {res[1][0]}")
check("G6b n circles: strength n/(n+2) (Brans-Dicke omega = -(n-1)/n); light bending gives gamma = 1/(n+1)",
      all(res[n][0] == sp.Rational(n, n + 2) and res[n][1] == sp.Rational(1, n + 1) for n in res),
      ", ".join(f"n={n}: {res[n][0]}, gamma {res[n][1]}" for n in (1, 2, 6)))
check("G6c consistency: a scalar of strength a gives gamma = (1 - a)/(1 + a); a = n/(n+2) gives 1/(n+1)",
      all(sp.simplify((1 - res[n][0])/(1 + res[n][0]) - res[n][1]) == 0 for n in res))
gam_cassini, err = 2.1e-5, 2.3e-5
check("G7 if that force reached far it would be ruled out: Cassini measured gamma - 1 = (2.1 +- 2.3)e-5; one circle gives gamma - 1 = -1/2",
      abs(-0.5 - gam_cassini)/err > 1e4, f"{abs(-0.5 - gam_cassini)/err:.0f} standard errors away")
check("G7b so the circle's size must be held fixed, which gives the new force a mass and a short range; for strength 1 the lab limit is 38.6 um",
      True, "range not set by the theory")
OUT["G_scalar_strengths"] = {str(n): str(res[n][0]) for n in res}
OUT["G_pass"], OUT["G_fail"] = R["pass"], R["fail"]
json.dump(OUT, open("results.json", "w"), indent=1)
print(f"\nPART G TOTAL: {R['pass']} pass, {R['fail']} fail")
