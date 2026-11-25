"""Elliptic counting: the n-fold sum on an elliptic curve is governed by an elliptic
divisibility sequence (EDS, Ward 1948), the elliptic twin of the Lucas sequence.

Curve: y² + a1 xy + a3 y = x³ + a2 x² + a4 x + a6. Point P. Division values ψ_n(P) give
    x(nP) = x(P) − ψ(n−1)ψ(n+1)/ψ(n)²,
and satisfy W(m+n)W(m−n)W(1)² = W(m+1)W(m−1)W(n)² − W(n+1)W(n−1)W(m)².
The landscape laws are the nodal curves (a3 = a4 = a6 = 0); there the EDS is a Lucas
sequence up to a gauge c^(n²−1).
"""
from fractions import Fraction as F
import cmath, json, math, os

HERE = os.path.dirname(os.path.abspath(__file__))
res = {"checks": [], "data": {}}


def check(name, ok, detail=""):
    res["checks"].append({"name": name, "ok": bool(ok), "detail": detail})
    print(("PASS " if ok else "FAIL ") + name + (f"  [{detail}]" if detail else ""))
    assert ok, name


class Curve:
    def __init__(s, a1, a2, a3, a4, a6):
        s.a = [F(a1), F(a2), F(a3), F(a4), F(a6)]
        a1, a2, a3, a4, a6 = s.a
        s.b2 = a1 * a1 + 4 * a2; s.b4 = 2 * a4 + a1 * a3; s.b6 = a3 * a3 + 4 * a6
        s.b8 = a1 * a1 * a6 + 4 * a2 * a6 - a1 * a3 * a4 + a2 * a3 * a3 - a4 * a4
        s.disc = -s.b2 ** 2 * s.b8 - 8 * s.b4 ** 3 - 27 * s.b6 ** 2 + 9 * s.b2 * s.b4 * s.b6

    def on(s, P):
        a1, a2, a3, a4, a6 = s.a; x, y = P
        return y * y + a1 * x * y + a3 * y == x ** 3 + a2 * x * x + a4 * x + a6

    def neg(s, P):
        a1, a2, a3, a4, a6 = s.a
        return None if P is None else (P[0], -P[1] - a1 * P[0] - a3)

    def add(s, P, Q):
        a1, a2, a3, a4, a6 = s.a
        if P is None: return Q
        if Q is None: return P
        x1, y1 = P; x2, y2 = Q
        if x1 == x2:
            if y1 + y2 + a1 * x2 + a3 == 0: return None
            lam = (3 * x1 * x1 + 2 * a2 * x1 + a4 - a1 * y1) / (2 * y1 + a1 * x1 + a3)
        else:
            lam = (y2 - y1) / (x2 - x1)
        nu = y1 - lam * x1
        x3 = lam * lam + a1 * lam - a2 - x1 - x2
        y3 = -(lam + a1) * x3 - nu - a3
        return (x3, y3)

    def psi(s, P, N):
        x, y = P; a1, a2, a3, a4, a6 = s.a
        b2, b4, b6, b8 = s.b2, s.b4, s.b6, s.b8
        p = {0: F(0), 1: F(1), 2: 2 * y + a1 * x + a3}
        p[3] = 3 * x ** 4 + b2 * x ** 3 + 3 * b4 * x ** 2 + 3 * b6 * x + b8
        p[4] = p[2] * (2 * x ** 6 + b2 * x ** 5 + 5 * b4 * x ** 4 + 10 * b6 * x ** 3 + 10 * b8 * x ** 2
                       + (b2 * b8 - b4 * b6) * x + (b4 * b8 - b6 * b6))
        for n in range(5, N + 1):
            m = n // 2
            if n % 2:
                p[n] = p[m + 2] * p[m] ** 3 - p[m - 1] * p[m + 1] ** 3
            else:
                p[n] = p[m] * (p[m + 2] * p[m - 1] ** 2 - p[m - 2] * p[m + 1] ** 2) / p[2]
        return [p[n] for n in range(N + 1)]


def eds_ok(W, N):
    for m in range(2, N):
        for n in range(1, m):
            if m + n > N: break
            if W[m + n] * W[m - n] * W[1] ** 2 != W[m + 1] * W[m - 1] * W[n] ** 2 - W[n + 1] * W[n - 1] * W[m] ** 2:
                return False
    return True


# 1. A smooth curve: y² + y = x³ − x, P = (0, 0) (rank one, the classic EDS)
E = Curve(0, 0, 1, -1, 0); P = (F(0), F(0))
W = E.psi(P, 30)
check("y² + y = x³ − x is smooth (discriminant ≠ 0)", E.disc != 0, f"Δ = {E.disc}")
check("its counting sequence is whole numbers: 1, 1, −1, 1, 2, −1, −3, −5, 7, −4, −23, 29, 59, …",
      all(w.denominator == 1 for w in W) and [int(w) for w in W[1:14]] == [1, 1, -1, 1, 2, -1, -3, -5, 7, -4, -23, 29, 59])
res["data"]["eds_smooth"] = [int(w) for w in W[:26]]
Q = None; ok = True
for n in range(1, 25):
    Q = E.add(Q, P)
    ok &= Q[0] == P[0] - W[n - 1] * W[n + 1] / W[n] ** 2 if n > 1 else True
check("x(nP) = x(P) − W(n−1)W(n+1)/W(n)², n = 2..24 (exact)", ok)
check("the EDS recurrence W(m+n)W(m−n) = W(m+1)W(m−1)W(n)² − W(n+1)W(n−1)W(m)² holds", eds_ok(W, 30))
Wi = [int(w) for w in W]
okg = all(math.gcd(Wi[m], Wi[n]) == abs(Wi[math.gcd(m, n)]) for m in range(1, 26) for n in range(1, 26))
check("strong divisibility gcd(W_m, W_n) = |W_gcd(m,n)|, like Fibonacci", okg)
ok4 = all(Wi[n + 4] * Wi[n] == Wi[n + 3] * Wi[n + 1] * Wi[2] ** 2 - Wi[3] * Wi[1] * Wi[n + 2] ** 2 for n in range(1, 25))
check("it is a Somos-4 sequence: W(n+4)W(n) = W(n+3)W(n+1) + W(n+2)² (since W2 = 1, W1·W3 = −1); its odd terms are the classic Somos-4 numbers 1, 1, 2, 3, 7, 23, 59 up to sign", ok4 and [abs(Wi[2*j-1]) for j in range(1,8)] == [1,1,2,3,7,23,59])

# 2. On the landscape (nodal curves) the EDS is a Lucas sequence up to a gauge c^(n²−1)
def nodal_vs_lucas(k, q, z):
    # curve y² + a1 xy = x³ + a2 x² carries the law at k = −a1, q = −a2 in z = −x/y
    a1, a2 = -k, -q
    E = Curve(a1, a2, 0, 0, 0)
    # the point with z: line y = t x through the node, t = −1/z, second hit x = t² + a1 t − a2
    t = -1 / z
    x = t * t + a1 * t - a2; y = t * x
    Pn = (x, y)
    assert E.on(Pn)
    W = E.psi(Pn, 16)
    a = z
    PP, QQ = 2 + k * a, 1 + k * a + q * a * a
    U = [F(0), F(1)]
    for _ in range(20): U.append(PP * U[-1] - QQ * U[-2])
    r = [W[n] / U[n] if U[n] != 0 else None for n in range(1, 16)]
    # gauge test: r(n) = A^(n²−1) B^(n−1) ⇔ r(n+1) r(n−1) / r(n)² constant
    g = [r[i + 1] * r[i - 1] / r[i] ** 2 for i in range(1, len(r) - 1) if None not in (r[i - 1], r[i], r[i + 1])]
    return len(set(g)) == 1, g[0] if g else None, eds_ok(W, 16)

okn = True; gauges = []
for k, q, z in [(1, 0, F(1, 3)), (0, -1, F(1, 2)), (0, 1, F(2, 5)), (-1, -1, F(1)), (2, 1, F(1, 4)), (F(3, 2), F(-2, 3), F(2, 7))]:
    same, gval, eok = nodal_vs_lucas(F(k), F(q), z)
    okn &= same and eok; gauges.append(str(gval))
check("on every landscape law tested, the curve's EDS = Lucas sequence × gauge A^(n²−1)B^(n−1)", okn,
      "the circle operation, velocity, slope, Fibonacci law, Euler char., a random law")

# 3. Lucas sequences satisfy the Somos-4 relation exactly, with a factor Q
okq = True
for PP, QQ in [(F(3), F(1)), (F(1), F(-1)), (F(5, 2), F(7, 3)), (F(-2), F(5))]:
    U = [F(0), F(1)]
    for _ in range(30): U.append(PP * U[-1] - QQ * U[-2])
    okq &= all(QQ * U[m + 4] * U[m] == U[m + 3] * U[m + 1] * U[2] ** 2 - U[3] * U[1] * U[m + 2] ** 2 for m in range(25))
check("Lucas: Q·U(m+4)U(m) = U(m+3)U(m+1)U(2)² − U(3)U(1)U(m+2)² (the q = 0 edge of Somos-4)", okq)

# 4. Theta tie (gluing number q): W_n = θ(wⁿ) satisfies the EDS recurrence; q = 0 gives Lucas
def theta(z, q, M=80):
    # z is the additive coordinate, w = e^z; w^(1/2) is taken as e^(z/2) so the branch is consistent
    v = cmath.exp(z / 2) - cmath.exp(-z / 2)
    w = cmath.exp(z)
    for m in range(1, M):
        v *= (1 - q ** m * w) * (1 - q ** m / w)
    return v

okt = True; worst = 0
for z, q in [(0.53 + 0.17j, 0.2), (0.08 + 0.98j, 0.35), (0.79, 0.05), (0.31 - 0.3j, 0.5)]:
    Wt = [theta(n * z, q) if n else 0 for n in range(0, 14)]
    for m in range(2, 7):
        for n in range(1, m):
            lhs = Wt[m + n] * Wt[m - n] * Wt[1] ** 2
            rhs = Wt[m + 1] * Wt[m - 1] * Wt[n] ** 2 - Wt[n + 1] * Wt[n - 1] * Wt[m] ** 2
            e = abs(lhs - rhs) / max(1e-300, abs(lhs) + abs(rhs))
            worst = max(worst, e)
okt = worst < 1e-9
check("theta along a line, W(n) = θ(n·z) with gluing q, satisfies the EDS recurrence", okt, f"worst rel. error {worst:.1e}")
z = 0.64 + 0.1j
okl = all(abs(theta(n * z, 0) - 2 * cmath.sinh(n * z / 2)) < 1e-12 * abs(cmath.exp(n * z / 2)) for n in range(1, 10))
check("at q = 0 the theta collapses to wⁿᐟ² − w⁻ⁿᐟ², a Lucas sequence (Q = 1) up to gauge", okl)

# 5. Cycles: torsion. Over the rationals a line law cycles only with period 1,2,3,4,6.
# Rational elliptic curves allow 1..10 and 12 (Mazur 1977): 5 included. Example: y² + y = x³ − x²
E5 = Curve(0, -1, 1, 0, 0); P5 = (F(0), F(0))
W5 = E5.psi(P5, 12)
pts = []; Q = None
for n in range(1, 6):
    Q = E5.add(Q, P5); pts.append(Q)
check("y² + y = x³ − x² has a rational point of order 5: W5 = 0 and 5P = O", W5[5] == 0 and pts[-1] is None,
      f"W1..W6 = {[int(v) for v in W5[1:7]]}")

# 6. The family with a point of order 5 (Tate normal form b = c): y² + (1−b)xy − by = x³ − bx².
#    It is smooth except where b⁵(b² − 11b − 1) = 0, i.e. b = 0 or b = (11 ± 5√5)/2 = φ⁵, −φ⁻⁵.
phi = (1 + 5 ** 0.5) / 2
def disc_tate5(b):
    E = Curve(1 - b, -b, -b, 0, 0)
    return E.disc
okd = all(disc_tate5(F(b)) == F(b) ** 5 * (F(b) ** 2 - 11 * F(b) - 1) for b in [2, -3, F(1, 7), F(5, 2), 9])
check("the 5-torsion family's discriminant is b⁵(b² − 11b − 1)", okd)
check("so it pinches to a node exactly at b = φ⁵ = 11.09 and b = −φ⁻⁵ = −0.090: golden",
      abs(phi ** 5 - (11 + 5 * 5 ** 0.5) / 2) < 1e-12 and abs(-phi ** -5 - (11 - 5 * 5 ** 0.5) / 2) < 1e-12)

# At that node: move the node to the origin, read off the law (k, q) and the point's value z,
# and check the 5-cycle on the line, with turn number φ² or 1/φ² (the golden turns).
def node_law(b):
    a1, a2, a3 = 1 - b, -b, -b
    # f(x,y) = y² + a1 xy + a3 y − x³ − a2 x²; singular point solves f = fx = fy = 0
    # fy = 2y + a1 x + a3 = 0 → y = −(a1 x + a3)/2 ; then fx = a1 y − 3x² − 2 a2 x = 0
    best = None
    for x0 in [i / 50 for i in range(-1000, 1000)]:
        pass
    # solve fx(x, y(x)) = 0: a1(−(a1x+a3)/2) − 3x² − 2a2 x = 0  (quadratic in x)
    A = -3; B = -2 * a2 - a1 * a1 / 2; C = -a1 * a3 / 2
    for x0 in [(-B + s * cmath.sqrt(B * B - 4 * A * C)) / (2 * A) for s in (1, -1)]:
        y0 = -(a1 * x0 + a3) / 2
        f = y0 * y0 + a1 * x0 * y0 + a3 * y0 - x0 ** 3 - a2 * x0 * x0
        if abs(f) < 1e-8:
            best = (x0, y0)
    x0, y0 = best
    # shift: X = x − x0, Y = y − y0. Quadratic part of f: Y² + a1 XY − (3x0 + a2) X²  (cubic −X³)
    # Curve Y² + a1' XY = X³ + a2' X² with a1' = a1, a2' = 3 x0 + a2 → law k = −a1', q = −a2'
    k = -a1; q = -(3 * x0 + a2)
    X, Y = 0 - x0, 0 - y0          # the point P = (0, 0) in shifted coordinates
    z = -X / Y
    return k, q, z

okg5 = True; info = []
for b in [phi ** 5, -phi ** -5]:
    k, q, z = node_law(b)
    PP = 2 + k * z; QQ = 1 + k * z + q * z * z
    tau = (PP * PP / QQ)
    U = [0, 1]
    for _ in range(8): U.append(PP * U[-1] - QQ * U[-2])
    okg5 &= abs(U[5]) < 1e-8 * max(1, abs(U[4])) and min(abs(tau - phi ** 2), abs(tau - phi ** -2)) < 1e-8
    info.append(f"b={b:.4f}: law k={complex(k).real:.4f}, q={complex(q).real:.4f}, turn number {complex(tau).real:.6f}")
check("at both golden nodes the order-5 point becomes a 5-cycle of a line law with turn number φ² or 1/φ²", okg5,
      "; ".join(info))
res["data"]["golden_nodes"] = info

# 7. The 5-fold family as a gluing: b = −R(q)⁵, R = Rogers–Ramanujan continued fraction.
#    The curve with that b is the Tate curve with gluing q⁵ (period 5τ), the point being q itself.
#    Gluing loose (q → 0): b → 0. Gluing tight (q → 1): R → 1/φ, b → −φ⁻⁵, the golden pinch.
def RR(q, N=600):
    v = q ** 0.2
    for n in range(1, N):
        v *= (1 - q ** (5 * n - 1)) * (1 - q ** (5 * n - 4)) / ((1 - q ** (5 * n - 2)) * (1 - q ** (5 * n - 3)))
    return v
def j_tate(Q, N=600):
    E4 = 1 + 240 * sum(sum(d ** 3 for d in range(1, n + 1) if n % d == 0) * Q ** n for n in range(1, 150))
    D = Q
    for n in range(1, N): D *= (1 - Q ** n) ** 24
    return E4 ** 3 / D
def j_b(b):
    a1, a2, a3 = 1 - b, -b, -b
    b2 = a1 * a1 + 4 * a2; b4 = a1 * a3; b6 = a3 * a3; b8 = a2 * a3 * a3
    D = -b2 ** 2 * b8 - 8 * b4 ** 3 - 27 * b6 ** 2 + 9 * b2 * b4 * b6
    return (b2 * b2 - 24 * b4) ** 3 / D
worst = 0; worst2 = 0
for q in [0.03, 0.08, 0.15, 0.25]:
    r = RR(q)
    eta = q ** -1
    for n in range(1, 400): eta *= (1 - q ** n) ** 6 / (1 - q ** (5 * n)) ** 6
    worst = max(worst, abs((1 / r ** 5 - 11 - r ** 5) - eta) / abs(eta))
    worst2 = max(worst2, abs(j_b(-r ** 5) - j_tate(q ** 5)) / abs(j_tate(q ** 5)))
check("Ramanujan: 1/R⁵ − 11 − R⁵ = q⁻¹∏(1−qⁿ)⁶/(1−q⁵ⁿ)⁶, the same b² − 11b − 1 with b = −R⁵", worst < 1e-9,
      f"worst rel. error {worst:.1e}")
check("the 5-fold curve with b = −R(q)⁵ is the Tate curve with gluing q⁵ (same j-invariant)", worst2 < 1e-8,
      f"worst rel. error {worst2:.1e}")
lim = [(q, RR(q)) for q in (0.9, 0.99)]
check("as the gluing tightens (q → 1), R(q) → 1/φ, so b → −φ⁻⁵: the golden pinch", all(abs(v - 1 / phi) < 1e-9 for _, v in lim),
      "; ".join(f"R({q}) = {v:.12f}" for q, v in lim) + f"; 1/φ = {1 / phi:.12f}")

res["summary"] = f"{sum(c['ok'] for c in res['checks'])} of {len(res['checks'])} checks pass"
print(res["summary"])
json.dump(res, open(os.path.join(HERE, "elliptic_counting_results.json"), "w"), indent=1, default=str)
