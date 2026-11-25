"""Patterns across the results: the counting recurrence behind every law.

Every law on the landscape is a ⊕ b = (a + b + k·ab) / (1 − q·ab).
Adding a to itself n times is the n-th power of the 2×2 step
    T(a) = [[1 + k·a, a], [−q·a, 1]],
so it is governed by one Lucas sequence U with
    P = 2 + k·a (trace),  Q = 1 + k·a + q·a² (determinant),
    U0 = 0, U1 = 1, U(n+2) = P·U(n+1) − Q·U(n),
    a^(⊕(n+1)) = a·U(n+1) / (U(n+1) − Q·U(n)).
Fibonacci is the Lucas sequence with P = 1, Q = −1.

All checks use exact fractions unless they say otherwise.
"""
from fractions import Fraction as F
import cmath, itertools, json, math, os, random

HERE = os.path.dirname(os.path.abspath(__file__))
out = {"checks": [], "data": {}}
npass = 0


def check(name, ok, detail=""):
    global npass
    out["checks"].append({"name": name, "ok": bool(ok), "detail": detail})
    print(("PASS " if ok else "FAIL ") + name + (f"  [{detail}]" if detail else ""))
    npass += bool(ok)
    assert ok, name


def law(k, q):
    return lambda a, b: (a + b + k * a * b) / (1 - q * a * b)


def lucas(P, Q, n):
    U = [0, 1]
    for i in range(n):
        U.append(P * U[-1] - Q * U[-2])
    return U


def fold(op, a, n):
    x = a
    for _ in range(n - 1):
        x = op(x, a)
    return x


# 1. The Lucas formula for every law on the landscape map, and random ones
dots = json.load(open(os.path.join(HERE, "..", "landscape", "shapes.json")))["dots"]
bad = 0
tested = 0
rng = random.Random(7)
laws = [(F(d["x"]).limit_denominator(10000), F(d["y"]).limit_denominator(10000), d["t"]) for d in dots]
laws += [(F(rng.randint(-9, 9), rng.randint(1, 5)), F(rng.randint(-9, 9), rng.randint(1, 5)), "random") for _ in range(40)]
for k, q, name in laws:
    op = law(k, q)
    for a in [F(1, 3), F(-2, 7), F(3, 5), F(1, 11)]:
        P, Q = 2 + k * a, 1 + k * a + q * a * a
        U = lucas(P, Q, 12)
        x = a
        for n in range(1, 12):
            den = U[n] - Q * U[n - 1]
            try:
                if den != 0 and a * U[n] / den != x:
                    bad += 1
                tested += 1
                x = op(x, a)
            except ZeroDivisionError:
                break  # ran into a wall
check("n-fold sum = a·U(n)/(U(n) − Q·U(n−1)) for all 16 map laws + 40 random laws",
      bad == 0, f"{tested} exact cases")

# 2. Associativity of the law IS the Lucas addition formula
okadd = True
for P, Q in [(3, 1), (1, -1), (F(5, 2), F(7, 3)), (-2, 5)]:
    U = lucas(P, Q, 30)
    for m in range(0, 12):
        for n in range(0, 12):
            okadd &= U[m + n + 1] == U[m + 1] * U[n + 1] - Q * U[m] * U[n]
check("Lucas addition U(m+n+1) = U(m+1)U(n+1) − Q·U(m)U(n)", okadd)

okcas = True
for P, Q in [(3, 1), (1, -1), (F(5, 2), F(7, 3)), (-2, 5)]:
    U = lucas(P, Q, 30)
    for n in range(0, 25):
        okcas &= U[n + 1] ** 2 - U[n] * U[n + 2] == Q ** n
check("Cassini: U(n+1)² − U(n)U(n+2) = Qⁿ (determinant of the n-step)", okcas)

# 3. The circle operation: the denominator is always 1, U(n) = ((1+a)^n − 1)/a
okm = True
for a in [F(1, 2), F(3), F(-1, 3)]:
    U = lucas(2 + a, 1 + a, 20)
    for n in range(1, 20):
        okm &= U[n] - (1 + a) * U[n - 1] == 1 and a * U[n] == (1 + a) ** n - 1
check("The circle operation: denominator U(n) − (1+a)U(n−1) = 1, so a^(⊕n) = (1+a)ⁿ − 1", okm)

# 4. The Fibonacci law: (k, q) = (−1, −1), adding 1 to itself
fibop = law(F(-1), F(-1))
fib = lucas(1, -1, 40)
okf = all(fold(fibop, F(1), n) == F(fib[n], fib[n + 1]) for n in range(1, 30))
check("law (a + b − ab)/(1 + ab): 1 ⊕ 1 ⊕ … (n times) = F(n)/F(n+1)", okf,
      "1, 1/2, 2/3, 3/5, 5/8, …")
phi = (1 + 5 ** 0.5) / 2
walls = sorted([(-1 + 5 ** 0.5) / 2, (-1 - 5 ** 0.5) / 2])
check("its walls are the roots of 1 − a − a² = 0: 1/φ and −φ",
      all(abs(1 - w - w * w) < 1e-12 for w in walls) and abs(walls[1] - 1 / phi) < 1e-12)
check("the counting marches into the wall 1/φ", abs(float(fold(fibop, F(1), 30)) - 1 / phi) < 1e-12)

# the ladder (simple-laws page): Z -> (1 + 2Z)/(1 + Z) = Z ⊕ 1 for the law at (k, q) = (1, −1)? check mirror
ladder = lambda z: (1 + 2 * z) / (1 + z)
mirror = law(F(1), F(-1))
okl = all(mirror(F(n, d), F(1)) == (F(n, d) + 1 + F(n, d)) / (1 + F(n, d)) for n, d in [(1, 2), (3, 7), (5, 3)])
check("the resistor ladder step is z ⊕ 1 in the mirror law (k, q) = (1, −1)", okl)
okmir = all(-law(F(-1), F(-1))(-a, -b) == law(F(1), F(-1))(a, b)
            for a, b in [(F(1, 3), F(2, 5)), (F(-3, 7), F(1, 9))])
check("mirror: −((−a) ⊕ (−b)) at (k, q) is a ⊕ b at (−k, q)", okmir)

# 5. Smallest irrational walls: integer laws with real walls, disc k² − 4q (in 1 + k a + q a²)
discs = sorted({k * k - 4 * q for k in range(-20, 21) for q in range(-50, 51) if k * k - 4 * q > 0})
nonsq = [d for d in discs if int(math.isqrt(d)) ** 2 != d]
check("smallest non-square wall discriminant of an integer law is 5 (the golden walls)", nonsq[0] == 5,
      f"next: {nonsq[1:5]}")
out["data"]["nonsquare_discs"] = nonsq[:10]

# 6. Turns: the turn number tau = P²/Q. Period n exactly when tau = 4cos²(πm/n) (Beraha numbers)
def period(k, q, a, maxn=40):
    op = law(k, q)
    x = a
    for n in range(1, maxn + 1):
        if abs(x) < 1e-9:
            return n
        try:
            x = op(x, a)
        except ZeroDivisionError:
            return None
    return None

okb = True
for n in range(3, 13):
    t = math.tan(math.pi / n)
    pr = period(0.0, 1.0, t)
    tau = 4 / (1 + t * t)
    okb &= pr == n and abs(tau - 4 * math.cos(math.pi / n) ** 2) < 1e-12
check("slope addition at a = tan(π/n) has period n and turn number 4cos²(π/n)", okb, "n = 3..12")
check("period 5 ⇔ turn number φ² = 2.618 or 1/φ² = 0.382 (golden)",
      abs(4 * math.cos(math.pi / 5) ** 2 - phi ** 2) < 1e-12 and abs(4 * math.cos(2 * math.pi / 5) ** 2 - phi ** -2) < 1e-12)
# U5 factorisation
okU5 = True
for P, Q in [(F(3), F(2)), (F(-1, 2), F(5)), (F(7, 3), F(-2))]:
    okU5 &= lucas(P, Q, 5)[5] == (P * P - Q) ** 2 - Q * P * P
check("U5 = (P² − Q)² − Q·P²: with Q = 1 it is (P² − P − 1)(P² + P − 1), golden roots", okU5)

# 7. Crystallographic restriction: integer (k, q, a) steps can only cycle with period 1, 2, 3, 4 or 6
periods = {}
for k in range(-6, 7):
    for q in range(-6, 7):
        for a in range(-6, 7):
            if a == 0:
                continue
            op = law(F(k), F(q))
            x = F(a)
            for n in range(1, 25):
                if x == 0:
                    periods[n] = periods.get(n, 0) + 1
                    break
                try:
                    x = op(x, F(a))
                except ZeroDivisionError:
                    break
check("integer laws only ever cycle with period 1, 2, 3, 4 or 6 (never 5)",
      set(periods) <= {1, 2, 3, 4, 6}, f"periods seen: {dict(sorted(periods.items()))}")
out["data"]["integer_periods"] = periods

# 8. Fibonacci inside the circle operation itself: x(n+2) = x(n+1) ∘₁ x(n)
ma = lambda a, b: a + b + a * b
xs = [1, 1]
for _ in range(20):
    xs.append(ma(xs[-1], xs[-2]))
fibn = lucas(1, -1, 25)
check("Circle-operation Fibonacci from 1, 1 is 2^F(n) − 1: 1, 1, 3, 7, 31, 255, 8191, …",
      all(xs[n] == 2 ** fibn[n + 1] - 1 for n in range(20)))
okg = all(math.gcd(xs[m], xs[n]) == 2 ** fibn[math.gcd(m + 1, n + 1)] - 1 for m in range(15) for n in range(15))
check("strong divisibility: gcd(x_m, x_n) = x_gcd (inherits Fibonacci's)", okg)
okgen = True
for x0, x1 in [(2, 5), (F(1, 2), F(1, 3))]:
    s = [x0, x1]
    for _ in range(12):
        s.append(ma(s[-1], s[-2]))
    okgen &= all(1 + s[n] == (1 + x0) ** fibn[n - 1] * (1 + x1) ** fibn[n] for n in range(1, 14))
check("general start: 1 + x(n) = (1 + x0)^F(n−1) · (1 + x1)^F(n)", okgen)

# 9. The curve that carries each law (elliptic page): count of the 3 countries over primes
# p − 1 points (two walls, split), p + 1 (turning), p (border). Golden law: split at p ≡ ±1 mod 5
def country_count(k, q, p):
    d = (k * k - 4 * q) % p
    if d == 0:
        return p
    return p - 1 if pow(d, (p - 1) // 2, p) == 1 else p + 1
primes = [p for p in range(7, 400) if all(p % r for r in range(2, int(p ** .5) + 1))]
okq = all((country_count(-1, -1, p) == p - 1) == (p % 5 in (1, 4)) for p in primes)
check("golden law's walls exist mod p exactly when p ≡ ±1 (mod 5): quadratic reciprocity for √5", okq)

out["summary"] = f"{npass} of {len(out['checks'])} checks pass"
print(out["summary"])
json.dump(out, open(os.path.join(HERE, "checks_results.json"), "w"), indent=1, default=str)
