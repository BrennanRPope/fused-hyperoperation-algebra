"""Every candidate 'elliptic object' for zeta, placed in the dial plane.  -> dial.json

Dial point of a quadratic factor 1 + k u + d u^2: discriminant Delta = k^2 - 4d,
elliptic (both roots on |u| = d^(-1/2)) iff Delta <= 0.  We also use the normalized
reading k_hat = k / sqrt(d), so that every factor sits on the line d = 1 and the
elliptic segment is |k_hat| <= 2.

A. Local candidates (one factor per prime p, variable u = p^(-s)):
   A1 zeta's Euler factor 1 - u paired with its mirror 1 - p u  -> (1-u)(1-pu)
   A2 a complex Dirichlet character paired with its conjugate   -> 1 - 2Re(chi(p)) u + u^2
   A3 (reference) an elliptic curve E: y^2 = x^3 + x + 1          -> 1 - a_p u + p u^2
B. Archimedean: trivial zeros s = -2m read through the Cayley map z = 1 - 1/s.
C. Global (Cayley) candidate: each nontrivial zero rho, paired with its mirror 1 - rho,
   gives (1 - z w)(1 - w/z) = 1 - (z + 1/z) w + w^2 with z = 1 - 1/rho.
D. Li coefficients from the zeros: lambda_n = sum over pairs 2(1 - T_n(-k/2)), plus a
   density tail; and a detection toy (one double zero split off the line).
"""
import json, math, cmath
import mpmath as mp
from sympy import primerange

out = {}
P = list(primerange(2, 200))

# A1: mirror-paired Euler factor
out["A1"] = [{"p": p, "k": -(p + 1), "d": p, "Delta": (p - 1) ** 2,
              "k_hat": -(p + 1) / math.sqrt(p)} for p in P]

# A2: characters mod 5 (order 4, chi(2) = i) and mod 7 (order 6, chi(3) = e^{i pi/3})
def chi_mod(q, g, order):
    lg = {}
    x = 1
    for e in range(q - 1):
        lg[x] = e; x = x * g % q
    return lambda n: 0 if n % q == 0 else cmath.exp(2j * math.pi * lg[n % q] / order)
A2 = []
for q, g, o in ((5, 2, 4), (7, 3, 6)):
    chi = chi_mod(q, g, o)
    for p in P:
        if p % q == 0: continue
        c = chi(p)
        A2.append({"q": q, "p": p, "k_hat": -2 * c.real, "d": 1, "Delta": 4 * c.real ** 2 - 4})
out["A2"] = A2

# A3: elliptic curve from the first round (weil.json)
w = json.load(open("../weil.json"))
out["A3"] = [{"p": r["p"], "a": r["a"], "k_hat": r["k"] / math.sqrt(r["p"]), "Delta": r["Delta"]}
             for r in w["genus1"]]

# B: trivial zeros through the Cayley map
B = []
for m in range(1, 41):
    z = 1 + 1 / (2 * m)
    B.append({"m": m, "z": z, "k": -(z + 1 / z)})
out["B"] = B

# C: Cayley dial of the nontrivial zeros
Z = json.load(open("zeros.json"))
gam = [float(g) for g in Z["gamma"]]
C = []
for g in gam:
    rho = complex(0.5, g)
    z = 1 - 1 / rho
    k = -(z + 1 / z)
    C.append({"gamma": g, "abs_z": abs(z), "k_re": k.real, "k_im": k.imag,
              "k_plus_2": k.real + 2, "inv_abs_rho2": 1 / abs(rho) ** 2,
              "curv": (1 - (k * k) / 4).real, "curv_formula": g * g / (g * g + 0.25) ** 2})
out["C"] = C

# C-off: what an off-line zero does (toy points, not real zeros)
Coff = []
for g0 in (5.0, 14.134725, 30.0):
    for dl in (0.05, 0.1, 0.2, 0.3, 0.4, 0.45):
        rho = complex(0.5 + dl, g0)
        z = 1 - 1 / rho
        k = -(z + 1 / z)
        gp = complex(g0, -dl)          # complex height: rho = 1/2 + i g'
        c = gp * gp / (gp * gp + 0.25) ** 2
        Coff.append({"gamma0": g0, "delta": dl, "abs_z": abs(z), "k_re": k.real, "k_im": k.imag,
                     "curv_re": c.real, "curv_im": c.imag})
# a real off-line pair (rho real in (0,1), as a toy): hyperbolic
for b in (0.6, 0.7, 0.8, 0.9):
    z = 1 - 1 / b
    Coff.append({"real_rho": b, "z": z, "k_re": -(z + 1 / z), "k_im": 0.0})
out["C_off"] = Coff

# D: Li coefficients from the zeros, pair form 2(1 - T_n(x)), x = -k/2 = cos(theta)
mp.mp.dps = 30
T_last = gam[-1]
thetas = [2 * math.atan(1 / (2 * g)) for g in gam]
def tail(n):
    # zeros above T_last by the Riemann-von Mangoldt density (1/2pi) log(t/2pi)
    f = lambda t: 2 * (1 - mp.cos(n * 2 * mp.atan(1 / (2 * t)))) * mp.log(t / (2 * mp.pi)) / (2 * mp.pi)
    return float(mp.quad(f, [T_last, 4 * T_last, 40 * T_last, mp.inf]))
NS = list(range(1, 1001)) + list(range(1050, 20001, 50))
D = []
for n in NS:
    s = sum(2 * (1 - math.cos(n * th)) for th in thetas)
    D.append({"n": n, "zero_sum": s, "tail": tail(n)})
out["D"] = D
out["T_last"] = T_last
out["N_zeros"] = len(gam)

# D-toy: split a double zero at height g0 off the line by delta (quartet replaces two on-line pairs)
def quartet(n, g0, dl):
    tot = 0
    for b in (0.5 + dl, 0.5 - dl):
        for sg in (1, -1):
            rho = complex(b, sg * g0)
            z = 1 - 1 / rho
            tot += (1 - z ** n).real
    return tot
toy = []
for g0, dl in ((5.0, 0.25), (5.0, 0.45), (14.134725, 0.45), (14.134725, 0.25)):
    th0 = 2 * math.atan(1 / (2 * g0))
    rows = []
    for r in D:
        n = r["n"]
        base = r["zero_sum"] + r["tail"]
        moved = base - 4 * (1 - math.cos(n * th0)) + quartet(n, g0, dl)
        rows.append([n, moved])
    first_neg = next((n for n, v in rows if v < 0), None)
    eps = math.log(abs(1 - 1 / complex(0.5 - dl, g0)))
    toy.append({"gamma0": g0, "delta": dl, "eps": eps, "first_negative_n": first_neg,
                "curve": rows})
out["toy"] = toy
json.dump(out, open("dial.json", "w"))
print("A1 max k_hat", max(r["k_hat"] for r in out["A1"]))
print("C: max |k_re|", max(abs(r["k_re"]) for r in C), "max |k+2 - 1/|rho|^2|",
      max(abs(r["k_plus_2"] - r["inv_abs_rho2"]) for r in C))
for t in toy: print("toy", t["gamma0"], t["delta"], "eps", t["eps"], "first neg n", t["first_negative_n"])
