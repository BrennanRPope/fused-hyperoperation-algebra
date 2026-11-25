"""Organic chemistry through the circle operation.
1. Rings: Hückel ring levels 2cos(2πj/N) and Möbius ring levels 2cos((2j+1)π/N) (units of β); characteristic
   polynomials 2T_N(x/2) ∓ 2. 4n + 2 electrons give a closed shell and positive resonance energy in Hückel rings,
   4n in Möbius rings. (Subject 4 did open chains; this closes them.)
2. Hammett: substituent constants add (k = 0), so rate ratios multiply and ratio − 1 composes by the circle operation.
3. Stereoselectivity: an enantiomer ratio r has ee = (r − 1)/(r + 1) = tanh(ΔΔG‡/2RT). Matched selectivities multiply r,
   so ee composes by (e₁ + e₂)/(1 + e₁e₂): the landscape law at (k, q) = (0, −1), the relativistic velocity-addition law,
   in the two-walls country with walls at ee = ±1. Mismatched pairs subtract. Horeau's duplication is e ⊕ e = 2e/(1 + e²).
4. Multi-step yields: losses fold at k = −1 (as in subject 5); convergent routes keep the longest branch.
Run: python3 checks.py"""
import json, math
import numpy as np
from numpy.polynomial import chebyshev as C
ok = []; R = {}
def check(name, cond, info=""):
    ok.append(bool(cond)); print(("PASS " if cond else "FAIL ") + name + ("  " + info if info else ""))
def w(k, a, b): return a + b + k*a*b
def land(k, q, a, b): return (a + b + k*a*b)/(1 - q*a*b)

# 1. rings
def ring(N, mobius=False):
    A = np.zeros((N, N))
    for i in range(N):
        A[i, (i+1) % N] = A[(i+1) % N, i] = 1
    if mobius: A[0, N-1] = A[N-1, 0] = -1
    return A
ok1 = True
for N in range(3, 11):
    for mob in (False, True):
        ev = np.sort(np.linalg.eigvalsh(ring(N, mob)))
        th = np.sort([2*math.cos(((2*j+1) if mob else 2*j)*math.pi/N) for j in range(N)])
        ok1 &= np.allclose(ev, th, atol=1e-10)
check("Hückel ring levels are 2cos(2πj/N); a Möbius twist shifts them to 2cos((2j+1)π/N), N = 3…10", ok1)
def charpoly_ok(N, mob):
    cp = np.poly(ring(N, mob))[::-1]  # coefficients low→high of det(xI − A)
    TN = C.cheb2poly([0]*N + [1])  # T_N(y), low→high
    # 2 T_N(x/2): scale coefficient of y^m by 2^(1−m)
    p = np.array([2*TN[m]/2**m for m in range(N+1)])
    p[0] += 2 if mob else -2
    return np.allclose(cp, p, atol=1e-8)
check("characteristic polynomial of the ring = 2T_N(x/2) − 2 (Hückel) and 2T_N(x/2) + 2 (Möbius)", all(charpoly_ok(N, m) for N in range(3, 11) for m in (False, True)))
def pi_energy(N, ne, mob):
    ev = sorted(np.linalg.eigvalsh(ring(N, mob)), reverse=True)
    occ = []; left = ne
    for e in ev:
        o = min(2, left); occ.append(o); left -= o
    E = sum(o*e for o, e in zip(occ, ev))
    homo = ev[(ne-1)//2]; lumo = ev[ne//2] if ne % 2 == 0 else ev[(ne-1)//2]
    gap = ev[(ne-1)//2] - ev[ne//2] if ne//2 < len(ev) else 0
    return E, E - ne, gap  # resonance vs ne/2 isolated ethylenes (each 2β per pair)
res = {}
for N in (4, 6, 8, 10, 12):
    res[N] = {"huckel": pi_energy(N, N, False)[:2] + (pi_energy(N, N, False)[2],), "mobius": pi_energy(N, N, True)}
R["rings"] = {N: {k: [float(x) for x in v] for k, v in d.items()} for N, d in res.items()}
check("benzene: π energy 8β, resonance 2β over three ethylenes", abs(res[6]["huckel"][0] - 8) < 1e-9 and abs(res[6]["huckel"][1] - 2) < 1e-9)
check("cyclobutadiene: resonance 0 and a degenerate half-filled pair (gap 0): antiaromatic", abs(res[4]["huckel"][1]) < 1e-9 and abs(res[4]["huckel"][2]) < 1e-9)
check("Möbius cyclobutadiene: levels ±√2 doubly, π energy 4√2 = 5.66β, resonance 1.66β, gap 2√2: aromatic", abs(res[4]["mobius"][0] - 4*math.sqrt(2)) < 1e-9 and abs(res[4]["mobius"][2] - 2*math.sqrt(2)) < 1e-9)
check("4n + 2 rule: Hückel rings with 6 and 10 electrons have a gap; 4, 8, 12 have none. Möbius is the reverse",
      all(res[N]["huckel"][2] > 0.5 for N in (6, 10)) and all(abs(res[N]["huckel"][2]) < 1e-9 for N in (4, 8, 12))
      and all(res[N]["mobius"][2] > 0.5 for N in (4, 8, 12)) and all(abs(res[N]["mobius"][2]) < 1e-9 for N in (6, 10)))

# 2. Hammett
rho = 2.5; s1, s2 = 0.23, 0.37
r1, r2 = 10**(rho*s1), 10**(rho*s2)
check("Hammett: σ adds over substituents (k = 0), so k/k₀ = 10^(ρΣσ) multiplies and k/k₀ − 1 composes by the circle operation", abs((10**(rho*(s1+s2)) - 1) - w(1, r1-1, r2-1)) < 1e-12)

# 3. ee
ee = lambda r: (r - 1)/(r + 1)
check("ee = (r − 1)/(r + 1) = tanh(ln r / 2): ee is a tanh of the free-energy difference", all(abs(ee(r) - math.tanh(math.log(r)/2)) < 1e-14 for r in (1.5, 10, 99)))
ra, rb = 10.0, 5.0
check("matched pair: ratios multiply (10 × 5 = 50) and ee composes by (e₁ + e₂)/(1 + e₁e₂): 0.818 ⊕ 0.667 = 0.961",
      abs(ee(ra*rb) - land(0, -1, ee(ra), ee(rb))) < 1e-14 and abs(ee(ra*rb) - 0.9608) < 1e-4)
check("mismatched pair: ratios divide, ee composes by (e₁ − e₂)/(1 − e₁e₂): 0.818 ⊖ 0.667 = 0.333", abs(ee(ra/rb) - land(0, -1, ee(ra), -ee(rb))) < 1e-14)
check("(0, −1) is in the two-walls country (k² − 4q = 4) with walls at ee = ±1: a perfectly selective step absorbs anything",
      0**2 - 4*(-1) == 4 and abs(land(0, -1, 1.0, 0.37) - 1) < 1e-15 and abs(land(0, -1, -1.0, 0.37) + 1) < 1e-15)
check("this is exactly relativistic velocity addition, with ee as v/c and ln r as twice the rapidity", abs(land(0, -1, 0.6, 0.8) - (0.6+0.8)/(1+0.48)) < 1e-15)
hor = [0.2]
for _ in range(4): hor.append(land(0, -1, hor[-1], hor[-1]))
R["horeau"] = hor
check("Horeau duplication: e → 2e/(1 + e²) is e ⊕ e; 20% ee climbs to 38%, 67%, 92%, 99.7% in four rounds",
      [round(100*x) for x in hor[1:]] == [38, 67, 92, 100] and abs(hor[4] - 0.99696) < 1e-4)
# selectivity in kJ/mol at 298 K
RT = 8.314462618*298.15/1000
need = {e: 2*RT*math.atanh(e) for e in (0.9, 0.99, 0.999)}; R["kJ_for_ee"] = need
check("ee 90%, 99%, 99.9% need ΔΔG‡ = 7.3, 13.1, 18.8 kJ/mol: each extra nine costs close to RT ln 10 = 5.7 kJ/mol",
      abs(need[0.9] - 7.30) < 0.02 and abs(need[0.99] - 13.12) < 0.02 and abs(need[0.999] - 18.83) < 0.02 and abs((need[0.999] - need[0.99]) - RT*math.log(10)) < 0.02)

# 4. yields
ys = [0.9, 0.85, 0.95, 0.8, 0.9, 0.7]
loss = 0
for y in ys: loss = w(-1, loss, 1 - y)
check("six-step linear route: losses fold at k = −1 to 1 − Πyᵢ (here 64% lost)", abs(loss - (1 - np.prod(ys))) < 1e-12 and abs(loss - 0.6338) < 1e-3)
lin = 0.9**8; conv = 0.9**4 * 0.9  # two 4-step branches then one coupling (per unit of the limiting branch)
check("eight 90% steps in a line keep 43%; two parallel four-step branches plus a 90% coupling keep 59% (only the longest branch counts)",
      abs(lin - 0.430) < 1e-3 and abs(conv - 0.590) < 1e-3)

print(f"\n{sum(ok)}/{len(ok)} checks pass")
json.dump(R, open("results.json", "w"), indent=1, default=float)
