"""Information causality and the 2 sqrt 2 ceiling (DRAFT, Nov 28 2026).
Run:  python3 checks.py > checks_output.txt   (writes results.json)

  I1  XOR of independent bits is the law with k = -2; its lift 1 - 2p (the bias) multiplies (piling-up lemma)
  I2  quantum boxes give the two biases (cos phi, sin phi): one turn of the law's rotation lift, the unit circle
  I3  Pawlowski et al.'s nested random access code: Bob's information stays <= 1 bit on and inside the circle,
      and exceeds 1 bit for every box outside it (PR box at once, r = 0.72 after a few levels)
  I4  the largest CHSH allowed by the circle is 2 sqrt 2 (Tsirelson); the local edge is E1 + E2 <= 1 (CHSH 2)
"""
import json
import numpy as np
import sympy as sp
from math import lgamma, log

results, passed, total = {}, 0, 0


def check(name, ok, detail, **extra):
    global passed, total
    total += 1; passed += bool(ok)
    print(f"[{'PASS' if ok else 'FAIL'}] {name}: {detail}")
    results[name] = {"pass": bool(ok), "detail": detail, **extra}


# I1
p, q, k = sp.symbols("p q k")
law = lambda a, b, kk: a + b + kk * a * b
xor_prob = p * (1 - q) + q * (1 - p)
ok1 = sp.expand(xor_prob - law(p, q, -2)) == 0 and sp.expand((1 - 2 * law(p, q, -2)) - (1 - 2 * p) * (1 - 2 * q)) == 0
check("I1", ok1,
      "P(a xor b = 1) = p + q - 2pq = p o_(-2) q exactly, and its fused-shift lift 1 + k p = 1 - 2p, the bias E, "
      "multiplies: E(a xor b) = E(a) E(b) (Matsui's piling-up lemma, 1993). Chaining boxes in Pawlowski's protocol "
      "is the law at k = -2")

# I2
sx = np.array([[0, 1], [1, 0]]); sz = np.diag([1, -1])
phi_plus = np.array([1, 0, 0, 1]) / np.sqrt(2)
corr = lambda A, B: float(phi_plus @ np.kron(A, B) @ phi_plus)
dev = 0
for ph in np.linspace(0, np.pi / 2, 31):
    A0 = np.cos(ph) * sz + np.sin(ph) * sx; A1 = np.cos(ph) * sz - np.sin(ph) * sx
    E1 = (corr(A0, sz) + corr(A1, sz)) / 2                 # Bob asks question 0: want a xor b = 0
    E2 = (corr(A0, sx) - corr(A1, sx)) / 2                 # Bob asks question 1: want a xor b = x
    dev = max(dev, abs(E1 + 1j * E2 - np.exp(1j * ph)))
check("I2", dev < 1e-12,
      f"a maximally entangled pair, Alice's two settings opened by angle phi around Bob's: the biases Bob gets for "
      f"his two questions are exactly (cos phi, sin phi) = (Re, Im) of e^(i phi), to {dev:.0e} over 31 angles. "
      f"Every such quantum box sits on the unit circle E1^2 + E2^2 = 1: one turn of the law's rotation lift "
      f"(imaginary k) hands its length out between the two questions")


# I3
def g(E):  # 1 - h((1+E)/2), bits of information carried by a bias E
    E = abs(E)
    if E >= 1: return 1.0
    if E < 0.5:  # stable series: sum_j E^(2j) / (2j (2j-1) ln 2)
        return sum(E ** (2 * j) / (2 * j * (2 * j - 1)) for j in range(1, 40)) / np.log(2)
    return ((1 + E) * np.log2(1 + E) + (1 - E) * np.log2(1 - E)) / 2


def info(E1, E2, n):
    # Pawlowski et al. 2009: N = 2^n bits, one bit sent; Bob's bit along a path with n-j zeros and j ones is
    # guessed with bias E1^(n-j) E2^j; information <= m = 1 is required
    tot = 0.0
    for j in range(n + 1):
        lb = lgamma(n + 1) - lgamma(j + 1) - lgamma(n - j + 1)
        tot += np.exp(lb) * g(E1 ** (n - j) * E2 ** j)
    return tot


ns = list(range(1, 301))
circle_max = max(max(info(np.cos(ph), np.sin(ph), n) for n in ns[:120]) for ph in np.linspace(0.05, np.pi / 2 - 0.05, 12))
inside_max = max(info(0.5, 0.5, n) for n in ns[:60])
first_viol = {}
for r in [1.0, 0.75, 0.72, 0.71]:
    first_viol[r] = next((n for n in ns if info(r, r, n) > 1), None)
uneq = next((n for n in ns if info(1.03 * np.cos(0.4), 1.03 * np.sin(0.4), n) > 1), None)
curves = {str(r): [info(r, r, n) for n in range(1, 41)] for r in [1.0, 0.75, 0.72, 1 / np.sqrt(2), 0.6]}
ok3 = circle_max <= 1 and inside_max <= 1 and all(v is not None for v in first_viol.values()) and uneq is not None
check("I3", ok3,
      f"nested random access code with one bit sent: on the circle Bob's information never passes "
      f"{circle_max:.4f} bit (12 angles, up to 120 levels); inside it (0.5, 0.5) at most {inside_max:.4f}. Outside "
      f"it, information causality breaks: PR box (1, 1) at level {first_viol[1.0]}, r = 0.75 at level "
      f"{first_viol[0.75]}, r = 0.72 at level {first_viol[0.72]}, r = 0.71 (just past 1/sqrt 2 = 0.7071) at level "
      f"{first_viol[0.71]}; an uneven box 3% outside at level {uneq}. Information causality allows exactly the disc "
      f"(Pawlowski et al., Nature 461, 1101, 2009)", first_violation={str(k): v for k, v in first_viol.items()},
      curves=curves)

# I4
E1s, E2s, lam = sp.symbols("E1 E2 lam", real=True)
sol = sp.solve([sp.diff(2 * (E1s + E2s) - lam * (E1s ** 2 + E2s ** 2 - 1), v) for v in (E1s, E2s)] +
               [E1s ** 2 + E2s ** 2 - 1], [E1s, E2s, lam], dict=True)
chsh_max = max(sp.nsimplify(2 * (s[E1s] + s[E2s])) for s in sol)
check("I4", sp.simplify(chsh_max - 2 * sp.sqrt(2)) == 0,
      f"CHSH = 2 (E1 + E2) for these boxes. Inside the disc its largest value is {chsh_max} (Tsirelson's bound), "
      f"reached at E1 = E2 = 1/sqrt 2, the 45-degree point of the turn; the local boxes fill E1 + E2 <= 1 (CHSH "
      f"2), and the PR box sits at the corner (1, 1) (CHSH 4)")

print(f"\n{passed}/{total} pass")
results["summary"] = f"{passed}/{total} pass"
json.dump(results, open("results.json", "w"), indent=1, default=float)
