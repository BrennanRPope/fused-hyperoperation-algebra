"""Checks for the dial map (DRAFT). Each check backs one concrete claim drawn in figures/.
Run: python3 checks.py  ->  checks_output.txt"""
import sympy as sp

a, b, k, d, x, u, v, s1, s2, t1 = sp.symbols("a b k d x u v s1 s2 t1")
results = []


def check(name, ok):
    results.append((name, bool(ok)))


circ = lambda p, q, kk: p + q + kk * p * q

# 1. a ∘k b is composition of scalings that share one centre, the wall −1/k (fig 3a)
scale = lambda p: (lambda y: (1 + k * p) * y + p)
comp = sp.expand(scale(a)(scale(b)(x)))
check("scalings about one centre compose by a ∘k b", sp.simplify(comp - scale(circ(a, b, k))(x)) == 0)
check("the shared centre is −1/k", sp.simplify(scale(a)(-1 / k) + 1 / k) == 0)
check("maps with one centre commute", sp.simplify(scale(a)(scale(b)(x)) - scale(b)(scale(a)(x))) == 0)

# 2. the shift x ↦ 1 + kx carries ∘k to ordinary multiplication (fig 1, row 2)
check("1 + k(a ∘k b) = (1 + ka)(1 + kb)", sp.expand(1 + k * circ(a, b, k) - (1 + k * a) * (1 + k * b)) == 0)

# 3. different centres: two orders differ, worked example of fig 3b and fig 5 panel 2
raise10 = lambda y: sp.Rational(11, 10) * y
grant = lambda y: y + 500
check("raise then grant on 1000 gives 1600", grant(raise10(1000)) == 1600)
check("grant then raise on 1000 gives 1650", raise10(grant(1000)) == 1650)
A1, C1, A2, C2 = sp.symbols("A1 C1 A2 C2")
f = lambda y: A1 * y + C1
g = lambda y: A2 * y + C2
gap = sp.expand(f(g(x)) - g(f(x)))
check("order gap of two affine maps is the oriented area (A1−1)C2 − (A2−1)C1",
      sp.simplify(gap - ((A1 - 1) * C2 - (A2 - 1) * C1)) == 0)

# 4. Rulecraft: same-target rules fuse at k = −1 (fig 5 panel 3)
rule = lambda s, tg: (lambda y: y + s * (tg - y))
fused = sp.expand(rule(s1, t1)(rule(s2, t1)(x)))
check("two rules with one target fuse to strength s1 + s2 − s1·s2",
      sp.simplify(fused - rule(circ(s1, s2, -1), t1)(x)) == 0)

# 5. the (k, d) plane placements (fig 4); Δ = k² − 4d
Delta = lambda kk, dd: kk ** 2 - 4 * dd
plane = lambda p, q, kk, dd: (p + q + kk * p * q) / (1 - dd * p * q)
check("velocity addition (c = 1) is the point (0, −1)",
      sp.simplify(plane(u, v, 0, -1) - (u + v) / (1 + u * v)) == 0)
check("tan addition is the point (0, 1)",
      sp.simplify(plane(sp.tan(a), sp.tan(b), 0, 1) - sp.tan(a + b)).rewrite(sp.sin).simplify() == 0)
check("Euler characteristic point (2, 1) lies on the parabola Δ = 0", Delta(2, 1) == 0)
check("the χ_y line d = k − 1 passes through (1, 0), (0, −1) and (2, 1)",
      all(dd == kk - 1 for kk, dd in [(1, 0), (0, -1), (2, 1)]))
check("velocity addition has two real walls (Δ > 0)", Delta(0, -1) > 0)
check("tan addition has none (Δ < 0)", Delta(0, 1) < 0)

# 6. the counts drawn in fig 2 and fig 6 match results.py
from results import RESULTS
from collections import Counter
c = Counter(r[2] for r in RESULTS)
check("30 results: 16 use ∘k, 6 the plane law or wider, 8 the picture only",
      len(RESULTS) == 30 and c["op"] == 16 and c["plane"] == 6 and c["pic"] == 8)

# polish: the order gap at 30 digits on random numbers
import random
random.seed(1)
vals = {A1: sp.Float(random.random(), 30), C1: sp.Float(random.random(), 30),
        A2: sp.Float(random.random(), 30), C2: sp.Float(random.random(), 30), x: sp.Float(random.random(), 30)}
check("polish: oriented-area gap holds to 30 digits",
      abs((gap - ((A1 - 1) * C2 - (A2 - 1) * C1)).subs(vals)) < sp.Float("1e-28"))

lines = [f"{'PASS' if ok else 'FAIL'}  {n}" for n, ok in results]
passed = sum(ok for _, ok in results)
lines.append(f"{passed} of {len(results)} pass")
out = "\n".join(lines)
print(out)
with open(__file__.replace("checks.py", "checks_output.txt"), "w") as fh:
    fh.write(out + "\n")
if passed != len(results):
    raise SystemExit(1)
