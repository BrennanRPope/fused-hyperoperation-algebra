"""Carnot cascade seen as a jet (Brennan, Nov 27: wherever velocity addition appears, consider the jet).
Run: python3 jet_checks.py   (exact fractions where possible, mpmath for tanh)"""
import json, random
from fractions import Fraction as Fr
import mpmath as mp
mp.mp.dps = 40
random.seed(1127)
res = []
def ok(name, c, note=""):
    res.append({"id": name, "pass": bool(c), "note": note}); print(("PASS " if c else "FAIL ") + name)

def eta_heat(etaC, g): return etaC / (2 - g * etaC)          # Schmiedl-Seifert heat-type family
def M(a, b): return ((a, b), (b, a))                          # a I + b J, J = [[0,1],[1,0]], J^2 = I
def mul(X, Y):
    return tuple(tuple(sum(X[i][k] * Y[k][j] for k in range(2)) for j in range(2)) for i in range(2))

# J1: the jet's velocity layer u = eta / (1 + (g-1) eta) is eta_C/(2 - eta_C) = tanh(beta/2) for every g
J1 = True
for _ in range(40):
    etaC, g = Fr(random.randint(1, 99), 100), Fr(random.randint(0, 100), 100)
    e = eta_heat(etaC, g)
    J1 &= e / (1 + (g - 1) * e) == etaC / (2 - etaC)
beta = mp.mpf("0.731"); etaC = 1 - mp.e ** (-beta)
J1 &= abs(etaC / (2 - etaC) - mp.tanh(beta / 2)) < mp.mpf(10) ** -35
ok("J1 every heat-type bound has the same velocity layer tanh(beta/2); gamma sits only in the scale 1 + (gamma-1) eta", J1)

# J2: N_g = (g-1) I + J has trace k = 2g-2 and det d = -g(2-g) (the S1 dial point) and gives the plane law
J2 = True
for _ in range(20):
    g = Fr(random.randint(0, 100), 100); k, d = 2 * g - 2, -g * (2 - g)
    N = M(g - 1, Fr(1))
    J2 &= N[0][0] + N[1][1] == k and N[0][0] * N[1][1] - N[0][1] * N[1][0] == d
    x, y = Fr(random.randint(1, 40), 100), Fr(random.randint(1, 40), 100)
    P = mul(M(1 + (g - 1) * x, x), M(1 + (g - 1) * y, y))
    # read back the scalar: P = (1 + (g-1) z) I + z J up to an overall scale s
    s_plus = P[0][0] - (g - 1) * P[0][1]; z = P[0][1] / s_plus
    J2 &= z == (x + y + k * x * y) / (1 - d * x * y)
ok("J2 the jet 1 + x N_g with N_g = (gamma-1) I + J multiplies to the S1 plane law (k, d) = (2gamma-2, -gamma(2-gamma))", J2)

# J3: stages with different gamma: jets commute (no turn) and the velocities add as rapidities
J3 = True
for _ in range(10):
    n = random.randint(2, 5)
    betas = [mp.mpf(random.randint(5, 120)) / 100 for _ in range(n)]
    gs = [mp.mpf(random.randint(0, 100)) / 100 for _ in range(n)]
    jets = []
    for b, g in zip(betas, gs):
        eC = 1 - mp.e ** (-b); e = eC / (2 - g * eC)
        jets.append(M(1 + (g - 1) * e, e))
    P = jets[0]
    for Jt in jets[1:]: P = mul(P, Jt)
    Q = jets[-1]
    for Jt in reversed(jets[:-1]): Q = mul(Q, Jt)
    J3 &= max(abs(P[i][j] - Q[i][j]) for i in range(2) for j in range(2)) < mp.mpf(10) ** -30
    J3 &= abs(P[0][1] / P[0][0] - mp.tanh(sum(betas) / 2)) < mp.mpf(10) ** -30
ok("J3 mixing stages with different gamma: jets commute (no turn) and the velocity is tanh(sum beta / 2)", J3)

# J4 contrast: realize two dial points by companion matrices instead; they do not commute (a turn appears)
def comp(k, d): return ((Fr(0), -d), (Fr(1), k))
A, B = comp(Fr(-2), Fr(0)), comp(Fr(0), Fr(-1))
AB, BA = mul(A, B), mul(B, A)
ok("J4 contrast: the same two end laws written as companion matrices do not commute, so the common J is a real choice",
   AB != BA)

print(f"\n{sum(r['pass'] for r in res)}/{len(res)} pass")
json.dump(res, open("jet_checks.json", "w"), indent=1)
