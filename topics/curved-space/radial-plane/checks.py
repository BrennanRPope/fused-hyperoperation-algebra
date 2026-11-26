"""Counting in a black hole's radial plane (DRAFT, Nov 26 2026).

Spherical symmetry turns 4D gravity into a problem on the 2D (t, r) plane. A point of that plane is a whole sphere
of area 4 pi r^2. A "spherical diamond" between two spheres p < q is everything inside the four light shells they
send in and out: in light-cone coordinates (u = t - r*, v = t + r*) it is a rectangle. We count:
  V2 = 2D volume of the rectangle (each point of the plane counted once), and
  V4 = 4D volume = integral of 4 pi r^2 over the rectangle (each sphere counted by its area).
Corners sit on a radial free-fall path with energy E (E = 1: a river parcel, falling from rest at infinity;
E^2 = 1 - 2M/r0: momentarily at rest at r0), at proper time -tau/2 and +tau/2 from the middle point r0.

Run: python3 -W ignore checks.py > checks_output.txt   (a few seconds). Writes results.json.
"""
import json, math
import numpy as np
import sympy as sp
from scipy.integrate import solve_ivp, dblquad
from scipy.special import lambertw

results, PASS, FAIL = {}, 0, 0


def check(name, ok, detail=""):
    global PASS, FAIL
    ok = bool(ok); PASS += ok; FAIL += (not ok)
    results[name] = {"pass": ok, "detail": str(detail)}
    print(("PASS " if ok else "FAIL ") + name + ("  | " + str(detail) if detail != "" else ""))


def rstar(r, M):
    return r if M == 0 else r + 2 * M * math.log(r / (2 * M) - 1)


def r_of(rs, M):
    return rs if M == 0 else 2 * M * (1 + lambertw(math.exp(rs / (2 * M) - 1)).real)


def corners(M, r0, E, tau):
    f = lambda r: 1 - 2 * M / r
    rd0 = -math.sqrt(max(E**2 - f(r0), 0.0))
    rhs = lambda s, y: [y[1], -M / y[0]**2, E / f(y[0])]
    out = []
    for sgn in (-1, 1):
        sol = solve_ivp(rhs, [0, sgn * tau / 2], [r0, rd0, 0.0], rtol=1e-13, atol=1e-14)
        r, _, t = sol.y[:, -1]
        out.append((t - rstar(r, M), t + rstar(r, M)))
    return out


def volumes(M, r0, E, tau):
    (up, vp), (uq, vq) = corners(M, r0, E, tau)
    assert uq > up and vq > vp

    def rf(u, v):
        r = r_of((v - u) / 2, M)
        return r, (1 - 2 * M / r if M else 1.0)
    V2 = dblquad(lambda v, u: rf(u, v)[1] / 2, up, uq, vp, vq, epsabs=1e-15, epsrel=1e-12)[0]
    V4 = dblquad(lambda v, u: 4 * math.pi * rf(u, v)[0]**2 * rf(u, v)[1] / 2, up, uq, vp, vq, epsabs=1e-12, epsrel=1e-12)[0]
    return V2, V4


cases = [(0, 10, 1.0), (0, 10, 1.3), (1, 10, math.sqrt(0.8)), (1, 10, 1.0), (1, 6, 1.0), (1, 6, math.sqrt(2 / 3)),
         (1, 8, 1.2), (1, 25, 1.0)]
small = np.array([0.04, 0.06, 0.08, 0.10, 0.12])


def leading(M, r0, E):
    taus = small * r0
    a2, a4 = [], []
    for t in taus:
        V2, V4 = volumes(M, r0, E, t)
        a2.append(V2 / (t**2 / 2) - 1)
        a4.append(V4 / (2 * math.pi * r0**2 * t**2) - 1)
    c2 = np.polyfit(taus, np.array(a2) / taus**2, 2)[-1]
    c4 = np.polyfit(taus, np.array(a4) / taus**2, 2)[-1]
    return c2, c4


rows = []
w2 = w4 = 0
for (M, r0, E) in cases:
    c2, c4 = leading(M, r0, E)
    K = 2 * M / r0**3
    w2 = max(w2, abs(c2 * 24 + K) * r0**3 / 2)   # error in units of M
    w4 = max(w4, abs(c4 * 24 * r0**2 - (2 * E**2 - 1)))
    rows.append({"M": M, "r0": r0, "E": round(E, 5), "2D coeff x 24": c2 * 24, "-K": -K, "4D coeff x 24 r0^2": c4 * 24 * r0**2,
                 "2E^2-1": 2 * E**2 - 1})
results["leading_rows"] = rows
check("D1 radial-plane count V2 = (tau^2/2)(1 - K tau^2/24 + ...), K = 2M/r0^3, for every free-fall energy (8 cases)",
      w2 < 1e-5, f"max error, in units of M, {w2:.1e}")
check("D2 sphere-weighted count V4 = 2 pi r0^2 tau^2 (1 + (2E^2 - 1) tau^2/(24 r0^2) + ...): the mass drops out "
      "except through the energy E (8 cases, flat and black hole)", w4 < 1e-5, f"max error {w4:.1e}")

# D3 flat space: the formula is exact (sympy), for an observer of rapidity eta crossing r0 at t = 0
r0s, tau, eta, u, v = sp.symbols('r0 tau eta u v', positive=True)
h = tau / 2
tq, rq = h * sp.cosh(eta), r0s - h * sp.sinh(eta)        # moving inward
tp, rp = -h * sp.cosh(eta), r0s + h * sp.sinh(eta)
up_, vp_, uq_, vq_ = tp - rp, tp + rp, tq - rq, tq + rq
V4flat = sp.integrate(sp.integrate(4 * sp.pi * ((v - u) / 2)**2 / 2, (v, vp_, vq_)), (u, up_, uq_))
target = 2 * sp.pi * r0s**2 * tau**2 * (1 + (2 * sp.cosh(eta)**2 - 1) * tau**2 / (24 * r0s**2))
check("D3 flat space, exactly: V4 = 2 pi r0^2 tau^2 (1 + (2 gamma^2 - 1) tau^2/(24 r0^2)), gamma = cosh(eta); "
      "2 gamma^2 - 1 = cosh(2 eta)", sp.simplify(sp.expand(V4flat - target).rewrite(sp.exp)) == 0)

# D4 why: the two mass effects and how they cancel (2D covariant calculus on -f dt^2 + dr^2/f)
t, r, M = sp.symbols('t r M', positive=True)
f = 1 - 2 * M / r
g = sp.diag(-f, 1 / f); X = [t, r]; gi = g.inv()
Gam = [[[sp.simplify(sum(gi[i, l] * (sp.diff(g[l, j], X[k]) + sp.diff(g[l, k], X[j]) - sp.diff(g[j, k], X[l]))
                     for l in range(2)) / 2) for k in range(2)] for j in range(2)] for i in range(2)]
hess = lambda phi: sp.Matrix(2, 2, lambda a, b: sp.simplify(sp.diff(phi, X[a], X[b]) - sum(Gam[c][a][b] * sp.diff(phi, X[c]) for c in range(2))))
Hr = hess(r)
check("D4a the sphere radius bends uniformly: grad grad r = (M/r^2) g", sp.simplify(Hr - (M / r**2) * g) == sp.zeros(2))
gradr2 = sp.simplify(sum(gi[a, b] * sp.diff(r, X[a]) * sp.diff(r, X[b]) for a in range(2) for b in range(2)))
H2 = hess(r**2)
box_r2 = sp.simplify(sum(gi[a, b] * H2[a, b] for a in range(2) for b in range(2)))
check("D4b |grad r|^2 = 1 - 2M/r and box(r^2) = 2, the same as in flat space (a vacuum equation)",
      sp.simplify(gradr2 - f) == 0 and box_r2 == 2, box_r2)
Ev, rdot = sp.symbols('E rdot', real=True)
# unit tangent of a radial path: u^t = E/f, u^r = rdot, with rdot^2 = E^2 - f
Tvec = sp.Matrix([Ev / f, rdot])
TT = sp.simplify((Tvec.T * H2 * Tvec)[0].subs(rdot**2, Ev**2 - f))
check("D4c along the path, T T grad grad (r^2) = 2 rdot^2 - 2M/r", sp.simplify(TT - (2 * (Ev**2 - f) - 2 * M / r)) == 0, TT)
# averaging r^2 over a small 2D diamond: <x^a x^b> = (tau^2/24)(T^a T^b + X^a X^b) = (tau^2/24)(g + 2TT)
sphere_part = sp.Rational(1, 2) * (box_r2 + 2 * TT)
dial_part = -(2 * M / r**3) * r**2                     # 2D area shortfall, -K r^2 with K = 2M/r^3
total = sp.simplify(sphere_part + dial_part)
check("D4d sphere bending 1 + 2 rdot^2 - 2M/r, plus the radial-plane dial -2M/r, totals 2E^2 - 1: M cancels",
      sp.simplify(total - (2 * Ev**2 - 1)) == 0, total)
check("D4e the river parcel (E = 1, rdot^2 = 2M/r) counts exactly like a hovering observer in empty space (coefficient 1)",
      total.subs(Ev, 1) == 1)
ms = sp.simplify((r / 2) * (1 - Ev**2 + (Ev**2 - f)))
check("D5 to read M you also need r and rdot from sphere areas: (r/2)(1 - E^2 + rdot^2) = M (the Misner-Sharp mass)",
      sp.simplify(ms - M) == 0)


# D6 the mass shows up directly one order later
def next_order(M, r0, E):
    taus = np.linspace(0.05, 0.5, 10) * r0
    a2, a4 = [], []
    K = 2 * M / r0**3
    for tt_ in taus:
        V2, V4 = volumes(M, r0, E, tt_)
        a2.append(V2 / (tt_**2 / 2) - 1 + K * tt_**2 / 24)
        a4.append(V4 / (2 * math.pi * r0**2 * tt_**2) - 1 - (2 * E**2 - 1) * tt_**2 / (24 * r0**2))
    p2 = np.polyfit(taus, np.array(a2) / taus**4, 2)[-1] * r0**4
    p4 = np.polyfit(taus, np.array(a4) / taus**4, 2)[-1] * r0**4
    return p2, p4


nxt = {}
for mr in (0.0, 0.02, 0.05, 0.1, 0.2):
    r0 = 10.0
    M = mr * r0
    nxt[mr] = next_order(M, r0, 1.0)
results["next_order_E1"] = {str(k): v for k, v in nxt.items()}
xs = np.array([k for k in nxt if k > 0]); y2 = np.array([nxt[k][0] for k in xs]); y4 = np.array([nxt[k][1] for k in xs])
A = np.vstack([xs, xs**2]).T
c2 = np.linalg.lstsq(A, y2, rcond=None)[0]; c4 = np.linalg.lstsq(A, y4, rcond=None)[0]
results["next_order_fit"] = {"2D": c2.tolist(), "4D": c4.tolist()}
check("D6 at the next order (tau/r0)^4 the mass appears directly; zero in flat space (river parcel, E = 1)",
      abs(nxt[0.0][0]) < 1e-8 and abs(nxt[0.0][1]) < 1e-8 and all(abs(nxt[k][1]) > 1e-6 for k in xs),
      f"fit a(M/r0) + b(M/r0)^2: 2D a = {c2[0]:.5f}, b = {c2[1]:.4f}; 4D a = {c4[0]:.6f}, b = {c4[1]:.4f}")

print(f"\n{PASS} of {PASS + FAIL} checks pass")
results["summary"] = {"pass": PASS, "total": PASS + FAIL}
json.dump(results, open("results.json", "w"), indent=1, default=str)
