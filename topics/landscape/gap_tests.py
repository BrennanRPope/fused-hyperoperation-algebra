"""Landscape map: numeric tests of the gaps.

For a growth law Y' = g(Y) started at Y0, measure displacement in units of the starting rate,
a = (Y - Y0)/g(Y0). Two runs of time t1 and t2 then compose by a law F(a, b) with identity 0.
If g is quadratic the law is exactly the Mobius law F_{p,d}, with
  p = g'(Y0),  d = g''(Y0) g(Y0) / 2.
This script integrates real growth laws and checks that claim, then tests the other gap claims.

Run:  python3 gap_tests.py   (writes gap_tests_output.txt and dots.json)
"""
import json
import numpy as np
from scipy.integrate import solve_ivp

out, npass, nall = [], 0, 0


def check(name, cond, detail=""):
    global npass, nall
    nall += 1
    npass += bool(cond)
    out.append(("PASS " if cond else "FAIL ") + name + (f"   {detail}" if detail else ""))


def mobius(a, b, p, q):
    return (a + b + p * a * b) / (1 - q * a * b)


def flow(g, Y0, t):
    sol = solve_ivp(lambda _t, y: [g(y[0])], [0, t], [Y0], rtol=1e-12, atol=1e-14)
    return sol.y[0, -1]


def law_test(name, g, dg, d2g, Y0, t1, t2):
    g0 = g(Y0)
    p, q = dg(Y0), d2g(Y0) * g0 / 2
    a = (flow(g, Y0, t1) - Y0) / g0
    b = (flow(g, Y0, t2) - Y0) / g0
    ab = (flow(g, Y0, t1 + t2) - Y0) / g0          # running t1 then t2 = running t1 + t2
    err = abs(ab - mobius(a, b, p, q))
    disc = p * p - 4 * q
    seq = p * p / disc if disc != 0 else np.inf
    check(f"{name}: two runs compose by F_(p,d), p = {p:.4g}, d = {q:.4g}", err < 1e-8, f"err {err:.1e}, s_eq {seq:.4g}")
    return p, q, seq


out.append("A. Growth laws that live off the straight line (exact Mobius laws)")
r, K = 0.8, 10.0
lg = lambda Y: r * Y * (1 - Y / K)
for y0 in (1.0, 3.0, 5.0, 8.0, 14.0):
    p, q, seq = law_test(f"logistic, start at {y0 / K:.0%} of capacity", lg, lambda Y: r * (1 - 2 * Y / K), lambda Y: -2 * r / K, y0, 0.7, 1.1)
    check(f"   unit-free shape s_eq = (1 - 2 Y0/K)^2 for logistic start {y0 / K:.0%}", abs(seq - (1 - 2 * y0 / K) ** 2) < 1e-9)

gr, cd = 9.81, 0.02            # fall with quadratic drag: v' = g - c v^2
law_test("falling body with quadratic drag, from rest (terminal speed = the wall)", lambda v: gr - cd * v * v,
         lambda v: -2 * cd * v, lambda v: -2 * cd, 0.0, 0.9, 1.6)
law_test("falling body shot down faster than terminal speed (both walls on one side)", lambda v: gr - cd * v * v,
         lambda v: -2 * cd * v, lambda v: -2 * cd, 40.0, 0.3, 0.5)

A_, Q_, R_ = -0.3, 2.0, 0.5    # scalar Kalman-Bucy variance: P' = 2AP + Q - P^2/R
kb = lambda P: 2 * A_ * P + Q_ - P * P / R_
for P0 in (0.0, 0.4, 3.0):
    law_test(f"Kalman-Bucy filter variance from P0 = {P0}", kb, lambda P: 2 * A_ - 2 * P / R_, lambda P: -2 / R_, P0, 0.2, 0.35)

out.append("")
out.append("B. Waves: the log-derivative w = psi'/psi of psi'' = (V - E) psi obeys w' = (V - E) - w^2")
for VmE, w0, lab in ((1.5, 0.3, "forbidden region (E < V): two real walls"),
                     (-1.5, 0.3, "allowed region (E > V): no real wall, it turns"),
                     (0.0, 0.3, "turning point (E = V): one double wall")):
    ricc = lambda w, VmE=VmE: VmE - w * w
    p, q, seq = law_test(f"Schrodinger, {lab}", ricc, lambda w: -2 * w, lambda w: -2.0, w0, 0.15, 0.2)
    disc = p * p - 4 * q
    want = np.sign(VmE)
    check(f"   region sign = sign(V - E): discriminant {disc:+.3g}", np.sign(round(disc, 12)) == want)

out.append("")
out.append("C. The straight line: symmetric noise cannot pass the parity wall")
for L in (2, 3, 4, 5, 8, 50):
    k = -L / (L - 1)
    worst = 1 + k * 1.0           # eigenvalue at error chance 1
    check(f"{L}-letter symmetric channel: weight {k:.4g}, eigenvalue at full error {worst:+.3f} >= -1", worst >= -1 - 1e-12)


def weyl(d):
    X = np.roll(np.eye(d), 1, axis=0)
    Z = np.diag(np.exp(2j * np.pi * np.arange(d) / d))
    return [np.linalg.matrix_power(X, i) @ np.linalg.matrix_power(Z, j) for i in range(d) for j in range(d)]


rng = np.random.default_rng(1)
for d in (2, 3):
    W = weyl(d)

    def dep(rho, pp, W=W, d=d):
        return (1 - pp) * rho + pp / (d * d - 1) * sum(w @ rho @ w.conj().T for w in W[1:])

    M = rng.normal(size=(d, d)) + 1j * rng.normal(size=(d, d))
    rho = M @ M.conj().T
    rho /= np.trace(rho)
    p1, p2 = 0.13, 0.27
    k = -d * d / (d * d - 1)
    err = np.abs(dep(dep(rho, p1), p2) - dep(rho, p1 + p2 + k * p1 * p2)).max()
    check(f"depolarizing qudit d = {d}: weight -d^2/(d^2-1) = {k:.4g}", err < 1e-12, f"err {err:.1e}")
lam = np.exp(2j * np.pi / 3)
check("past the wall the line turns complex: a 3-cycle shift has eigenvalues on the unit circle (a turn), not below -1",
      abs(abs(lam) - 1) < 1e-12 and abs(np.linalg.eigvals(np.roll(np.eye(3), 1, axis=0))).max() - 1 < 1e-12)
kok = 0.056 + 0.888j
check("oxygen clock (Kok cycle, Shinkarev's fit): one flash turns by about a quarter circle",
      80 < np.degrees(np.angle(kok)) < 95, f"{np.degrees(np.angle(kok)):.1f} degrees, period {360 / np.degrees(np.angle(kok)):.2f} flashes")

out.append("")
out.append("D. Fractional rungs are not analytic laws (the structural gap)")
# coefficients of the formal s-th iterate of E(x) = e^x - 1: each coefficient c_k(s) is a polynomial
# in s of degree k - 1, so compute E^n for n = 0..K-1 exactly and interpolate at s = 1/2 and -1/2.
from fractions import Fraction
from math import factorial
K = 26
Ecoef = [Fraction(0)] + [Fraction(1, factorial(n)) for n in range(1, K + 1)]


def mul(A, B):
    C = [Fraction(0)] * (K + 1)
    for i, ai in enumerate(A):
        if ai:
            for j in range(K + 1 - i):
                C[i + j] += ai * B[j]
    return C


def compose(outer, inner):          # outer(inner(x)), inner(0) = 0
    res = [Fraction(0)] * (K + 1)
    pw = [Fraction(1)] + [Fraction(0)] * K
    for n in range(1, K + 1):
        pw = mul(pw, inner)
        if outer[n]:
            res = [r + outer[n] * w for r, w in zip(res, pw)]
    return res


iters = [[Fraction(0), Fraction(1)] + [Fraction(0)] * (K - 1)]
for n in range(1, K):
    iters.append(compose(Ecoef, iters[-1]))


def at(sv, k):                       # Lagrange interpolation of c_k(n), n = 0..k-1
    tot = Fraction(0)
    for i in range(k):
        term = iters[i][k]
        for j in range(k):
            if j != i:
                term *= (sv - j) / Fraction(i - j)
        tot += term
    return tot


for sv in (Fraction(1, 2), Fraction(-1, 2), Fraction(1, 3)):
    roots = [abs(float(at(sv, k))) ** (1 / k) for k in (8, 14, 20, 26)]
    check(f"rung {sv}: n-th roots of its series coefficients keep growing (radius zero, Baker 1958)",
          roots[-1] > roots[0] and roots[-1] > roots[-2], " -> ".join(f"{r:.3f}" for r in roots))
c2 = at(Fraction(1, 2), 2)
check("interpolation sanity: the half rung's x^2 coefficient is s/2 = 1/4", c2 == Fraction(1, 4))
check("interpolation sanity: integer s = 2 gives E(E(x)) exactly", all(at(Fraction(2), k) == iters[2][k] for k in range(2, 12)))

# --- the dots -------------------------------------------------------------------------------
dots = [
    # name, domain, c1, c2, family, source page
    ("Compound growth, returns, interest", "finance", 1, 0, "line", "econ-forecasting"),
    ("Redshifts stacking (1+z multiplies)", "cosmology", 1, 0, "line", "dark-sector"),
    ("K-theory, Todd genus", "topology", 1, 0, "line", "topology"),
    ("Claim-count models (PGF)", "insurance", 1, 0, "line", "actuarial"),
    ("PUE: power and cooling overheads", "engineering", 1, 0, "line", "datacenter-energy"),
    ("Tsallis q_T = 0 (outcomes minus one), Campbell cost", "information", 1, 0, "line", "information-theory"),
    ("Shannon and Renyi entropy, hazards, log-returns", "information", 0, 0, "line", "information-theory"),
    ("Burgers shocks (empty corner)", "fluids", 0, 0, "line", "ns-rungs"),
    ("Union of failures, survival, decrements", "reliability", -1, 0, "line", "information-theory"),
    ("Amplitude damping (T1)", "quantum", -1, 0, "line", "quantum-noise"),
    ("Gini-Simpson (Tsallis q_T = 2), Kraft slack, erasures", "information", -1, 0, "line", "information-theory"),
    ("Slashing, validator churn", "consensus", -1, 0, "line", "consensus"),
    ("Parity, bit flips, binary channel cascade", "codes", -2, 0, "line", "information-theory"),
    ("Bit-flip and phase-flip noise", "quantum", -2, 0, "line", "quantum-noise"),
    ("3-letter channel", "channels", -1.5, 0, "line", "information-theory"),
    ("4-letter channel, depolarizing qubit", "quantum", -4 / 3, 0, "line", "quantum-noise"),
    ("Depolarizing qutrit", "quantum", -9 / 8, 0, "line", "this page"),
    ("Leveraged compounding, 2x", "finance", 2, 0, "line", "this page"),
    ("Partial damage, each hit removes half", "reliability", -0.5, 0, "line", "this page"),
    ("Velocity addition (relativity), stacked reflections, signature", "physics", 0, -1, "axis", "information-theory, topology"),
    ("Slope addition (tangent of angle sum)", "geometry", 0, 1, "axis", "this page"),
    ("Sine addition (third-order point only)", "geometry", 0, -0.5, "axis", "this page"),
    ("Euler characteristic genus (a cusp)", "topology", 2, 1, "parabola", "this page"),
    ("Rung 2 (smooth fused exponent)", "fused tower", 2, 0.5, "rung", "fractional-rungs"),
    ("Rung -1: log(e^a + e^b - 1)", "fused tower", -1, 0.5, "rung", "fractional-rungs"),
    ("Half rung", "fused tower", 0.5, -1 / 16, "rung", "fractional-rungs"),
]
json.dump([dict(name=n, domain=d, c1=a, c2=b, family=f, source=s) for n, d, a, b, f, s in dots],
          open(__file__.replace('gap_tests.py', 'dots.json'), 'w'), indent=1)
out.append("")
out.append(f"{len(dots)} dots written to dots.json")
out.append(f"{npass} of {nall} numeric checks pass")
open(__file__.replace('.py', '_output.txt'), 'w').write("\n".join(out) + "\n")
print("\n".join(out))
