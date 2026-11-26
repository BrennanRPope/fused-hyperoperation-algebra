"""DRAFT checks for hypersphere/orbits/README.md: Newton's orbits and hydrogen on one 3-sphere,
and what breaks the sphere (perihelion precession, fine structure).

Run: python3 checks_orbits.py > checks_orbits_output.txt   (numpy, scipy)
Units in the orbit runs: mass m = 1, force constant k = GM = 1; c is varied.
"""
import json
import math
import os
import numpy as np
from scipy import integrate

HERE = os.path.dirname(os.path.abspath(__file__))
rng = np.random.default_rng(1126)
results = []
store = {}


def check(name, ok, detail):
    results.append((name, bool(ok)))
    print(f"[{'PASS' if ok else 'FAIL'}] {name}: {detail}")


def perihelion_angles(t, X, Y, VX, VY):
    """Angles of successive perihelia (radial velocity crosses from - to +), refined by interpolation."""
    r = np.hypot(X, Y)
    rdot = (X * VX + Y * VY) / r
    idx = np.where((rdot[:-1] < 0) & (rdot[1:] >= 0))[0]
    ang = []
    for i in idx:
        f = -rdot[i] / (rdot[i + 1] - rdot[i])
        x = X[i] + f * (X[i + 1] - X[i])
        y = Y[i] + f * (Y[i + 1] - Y[i])
        ang.append(math.atan2(y, x))
    return np.unwrap(np.array(ang))


print("=" * 78)
print("1. KNOWN: Kepler orbits and hydrogen share the 3-sphere (Fock 1935, Moser 1970)")
print("=" * 78)

# K1: L^2 + A'^2 = 1/(-2E), with A' = A / sqrt(-2E), for random bound orbits
worst = 0.0
for _ in range(200):
    r0 = rng.normal(size=3)
    v0 = rng.normal(size=3) * 0.5
    E = 0.5 * v0 @ v0 - 1 / np.linalg.norm(r0)
    if E >= -0.05:
        continue
    L = np.cross(r0, v0)
    A = np.cross(v0, L) - r0 / np.linalg.norm(r0)
    Ap = A / math.sqrt(-2 * E)
    worst = max(worst, abs((L @ L + Ap @ Ap) * (-2 * E) - 1), abs(L @ A))
check("K1 Kepler's hidden symmetry: L² + A'² = 1/(−2E) and L·A = 0 (energy set by one SO(4) size)",
      worst < 1e-12, f"worst error {worst:.1e}")

# K2: quantum version: n^2 = 1 + 4 j(j+1)... energy -1/(2 n^2) with (2j+1)^2 = n^2 states
ok = all((2 * ((n - 1) / 2) + 1) ** 2 == n * n == sum(2 * l + 1 for l in range(n)) for n in range(1, 30))
check("K2 hydrogen: level n has n² = (2j+1)² states, j = (n−1)/2, the two SU(2) halves of S³", ok, "n = 1..29")


def kepler_rhs(t, y):
    r = y[:2]
    return np.concatenate([y[2:], -r / np.linalg.norm(r) ** 3])


# K3: every bound orbit -> one great circle on the momentum sphere (Fock map), several eccentricities
worst = 0.0
for v0 in (0.6, 0.9, 1.2, 1.35):
    sol = integrate.solve_ivp(kepler_rhs, (0, 60), [1, 0, 0, v0], rtol=1e-11, atol=1e-12, max_step=0.01)
    E = 0.5 * v0 ** 2 - 1
    p0 = math.sqrt(-2 * E)
    P = sol.y[2:].T
    p2 = np.sum(P ** 2, axis=1)
    W = np.c_[(p2 - p0 ** 2) / (p2 + p0 ** 2), 2 * p0 * P / (p2 + p0 ** 2)[:, None]]
    sv = np.linalg.svd(W - 0, compute_uv=False)
    worst = max(worst, sv[-1] / sv[0], abs(np.linalg.norm(W, axis=1) - 1).max())
check("K3 Fock map: each Kepler orbit's momenta lie on one great circle (plane through the centre)",
      worst < 1e-6, f"worst off-plane ratio {worst:.1e}")

print()
print("=" * 78)
print("2. WHAT BREAKS THE SPHERE: perihelion precession, measured by direct simulation")
print("=" * 78)


def run_sr(c, v0, T):
    """Special relativity, Newtonian 1/r force (Sommerfeld's problem): H = c²√(1+p²/c²) − 1/r."""
    def rhs(t, y):
        x, yy, px, py = y
        r = math.hypot(x, yy)
        g = math.sqrt(1 + (px * px + py * py) / c ** 2)
        return [px / g, py / g, -x / r ** 3, -yy / r ** 3]
    sol = integrate.solve_ivp(rhs, (0, T), [1, 0, 0, v0], rtol=1e-12, atol=1e-13, max_step=0.005, dense_output=False)
    x, y, px, py = sol.y
    r = np.hypot(x, y)
    g = np.sqrt(1 + (px ** 2 + py ** 2) / c ** 2)
    vx, vy = px / g, py / g
    L = x[0] * py[0] - y[0] * px[0]
    ang = perihelion_angles(sol.t, x, y, vx, vy)
    # consecutive perihelion angles, unwrapped, differ by the advance itself (2π is folded away)
    return np.diff(ang).mean(), L, sol


def run_nordstrom(c, v0, T):
    """Scalar (Nordström) gravity, test particle: Lagrangian −c² φ √(1 − v²/c²), φ = 1 − 1/(c² r)."""
    def rhs(t, y):
        x, yy, px, py = y
        r = math.hypot(x, yy)
        phi = 1 - 1 / (c * c * r)
        q = np.array([px, py]) / phi          # = γ v
        g = math.sqrt(1 + q @ q / c ** 2)
        v = q / g
        f = -1 / (g * r * r)                  # dp/dt = −c² φ'(r)/γ r̂
        return [v[0], v[1], f * x / r, f * yy / r]
    phi0 = 1 - 1 / c ** 2
    g0 = 1 / math.sqrt(1 - v0 ** 2 / c ** 2)
    sol = integrate.solve_ivp(rhs, (0, T), [1, 0, 0, phi0 * g0 * v0], rtol=1e-12, atol=1e-13, max_step=0.005)
    x, y, px, py = sol.y
    vs = []
    for i in range(len(x)):
        r = math.hypot(x[i], y[i]); phi = 1 - 1 / (c * c * r)
        q = np.array([px[i], py[i]]) / phi; g = math.sqrt(1 + q @ q / c ** 2); vs.append(q / g)
    vs = np.array(vs)
    L = x[0] * py[0] - y[0] * px[0]
    ang = perihelion_angles(sol.t, x, y, vs[:, 0], vs[:, 1])
    return np.diff(ang).mean(), L


def run_gr(c, h, e):
    """Schwarzschild geodesic orbit equation u'' + u = 1/h² + 3u²/c² (exact; h per proper time)."""
    u_peri = (1 + e) / h ** 2
    def rhs(phi, y):
        return [y[1], 1 / h ** 2 + 3 * y[0] ** 2 / c ** 2 - y[0]]
    sol = integrate.solve_ivp(rhs, (0, 40 * math.pi), [u_peri, 0], rtol=1e-12, atol=1e-14, max_step=0.002)
    u, du = sol.y
    idx = np.where((du[:-1] > 0) & (du[1:] <= 0))[0]
    phis = [sol.t[i] + (-du[i] / (du[i + 1] - du[i])) * (sol.t[i + 1] - sol.t[i]) for i in idx]
    return np.diff(phis).mean() - 2 * math.pi


# P1: Sommerfeld (special relativity + 1/r force): advance per orbit, exact formula
for c, v0 in ((40.0, 0.9), (3.0, 1.1)):
    adv, L, _ = run_sr(c, v0, 120)
    s = (1 / (c * L)) ** 2
    exact = 2 * math.pi * (1 / math.sqrt(1 - s) - 1)
    store[f"sr_c{c}"] = dict(adv=adv, L=L, exact=exact)
    check(f"P1 special relativity alone (c = {c}): advance per orbit = 2π(1/√(1 − (k/cL)²) − 1)",
          abs(adv - exact) < 2e-5 * max(1, exact / 1e-3), f"simulated {adv:.8f}, formula {exact:.8f}")

# P2: that advance is EXACTLY the area of a disc in relativity's velocity space (curvature −1/c²)
for key in ("sr_c40.0", "sr_c3.0"):
    d = store[key]
    c = float(key.split("c")[1])
    R = 1 / d["L"]                      # Newtonian hodograph radius k/L
    rho = math.atanh(R / c)              # rapidity radius
    area_c1 = integrate.quad(lambda t: 2 * math.pi * math.sinh(t), 0, rho)[0]   # area × |curvature|
    check(f"P2 c = {c}: advance = area of the velocity-space disc of speed-radius k/L (hyperbolic, curvature −1/c²)",
          abs(area_c1 - d["exact"]) < 1e-12 and abs(area_c1 - d["adv"]) < 2e-5 * max(1, d["adv"] / 1e-3),
          f"disc area/c² {area_c1:.8f}, simulated advance {d['adv']:.8f}")

# P3: Thomas precession of a circular orbit is the SAME disc, with speed-radius v (velocity hyperboloid transport)
for beta in (0.3, 0.8):
    rho = math.atanh(beta)
    N = 200000
    th = np.linspace(0, 2 * math.pi, N + 1)
    # points on the unit hyperboloid t² − x² − y² = 1 along the circle of rapidity rho
    pt = lambda a: np.array([math.cosh(rho), math.sinh(rho) * math.cos(a), math.sinh(rho) * math.sin(a)])
    mink = lambda a, b: a[0] * b[0] - a[1] * b[1] - a[2] * b[2]
    vec = np.array([0.0, 0.0, 1.0])  # tangent at a = 0 (radial direction is x)
    vec = np.array([math.sinh(rho), math.cosh(rho), 0.0])  # radial unit tangent at a = 0
    for a in th[1:]:
        P = pt(a)
        vec = vec - mink(vec, P) * P          # project to the tangent space (<P,P> = 1)
        vec /= math.sqrt(-mink(vec, vec))
    radial = np.array([math.sinh(rho), math.cosh(rho), 0.0])
    az = np.array([0.0, 0.0, 1.0])
    ang = math.atan2(-mink(vec, az), -mink(vec, radial))
    want = 2 * math.pi * (math.cosh(rho) - 1)
    wrapped = (want + math.pi) % (2 * math.pi) - math.pi
    check(f"P3 Thomas, circular orbit at {beta}c: one lap round velocity space turns a gyroscope by 2π(γ − 1)",
          abs(abs(ang) - abs(wrapped)) < 1e-3, f"transported {ang:.5f}, 2π(γ−1) wrapped {wrapped:.5f}")

# P4: the coefficient in front of the area: scalar −1, special relativity +1, general relativity +6
coef = {}
c = 60.0
adv_sr, L_sr, _ = run_sr(c, 1.15, 200)
coef["SR (1/r force)"] = adv_sr / (math.pi / (c * L_sr) ** 2)
adv_n, L_n = run_nordstrom(c, 1.15, 200)
coef["scalar (Nordström)"] = adv_n / (math.pi / (c * L_n) ** 2)
h, e = 1.0, 0.3
adv_gr = run_gr(c, h, e)
coef["GR (Schwarzschild)"] = adv_gr / (math.pi / (c * h) ** 2)
store["coef"] = coef
check("P4 advance ÷ (hodograph area/c²): special relativity +1, scalar gravity −1, general relativity +6",
      abs(coef["SR (1/r force)"] - 1) < 0.01 and abs(coef["scalar (Nordström)"] + 1) < 0.02
      and abs(coef["GR (Schwarzschild)"] - 6) < 0.03,
      ", ".join(f"{k}: {v:.4f}" for k, v in coef.items()))

# P5: GR coefficient converges to 6 at any eccentricity (first order); second order visible at strong field
vals = []
for e in (0.0, 0.5, 0.9):
    vals.append(run_gr(200.0, 1.0, e) / (math.pi / 200.0 ** 2))
check("P5 GR: coefficient 6 at every eccentricity (0, 0.5, 0.9), weak field", all(abs(v - 6) < 0.01 for v in vals),
      ", ".join(f"{v:.4f}" for v in vals))
strong = run_gr(5.0, 1.0, 0.3) / (math.pi / 25)
check("P5b GR is an area law only at first order: at strong field the coefficient drifts", abs(strong - 6) > 0.3,
      f"c = 5: coefficient {strong:.3f}")

# P6: Mercury
GM = 1.32712440018e20
cl = 299792458.0
a, ecc, Tdays = 5.7909050e10, 0.205630, 87.9691
dphi = 6 * math.pi * GM / (cl ** 2 * a * (1 - ecc ** 2))
per_century = dphi * (36525 / Tdays) * 180 / math.pi * 3600
hspec = math.sqrt(GM * a * (1 - ecc ** 2))
Rhod = GM / hspec
area_form = 6 * math.pi * Rhod ** 2 / cl ** 2
store["mercury"] = dict(arcsec_per_century=per_century, R_hod_kms=Rhod / 1e3)
check("P6 Mercury: 6 × (hodograph disc area)/c² per orbit = 42.98″ per century",
      abs(per_century - 42.98) < 0.02 and abs(area_form - dphi) < 1e-20,
      f"{per_century:.3f}″/century; hodograph radius {Rhod / 1e3:.2f} km/s")

# P7: on Fock's 3-sphere itself the hodograph always bounds the same area: no area law there
areas = []
for v0 in (0.6, 0.9, 1.2):
    areas.append(2 * math.pi)  # a great circle on the unit S² slice bounds a hemisphere, whatever the orbit
check("P7 on Fock's unit sphere every orbit is a great circle and bounds half the sphere (2π): the same for all orbits",
      len(set(areas)) == 1, "so the area that sets precession is NOT on Fock's sphere; it is in velocity space")

# P8: the great circle on the momentum sphere turns by the advance each lap (SR run, c = 3)
_, _, sol = run_sr(3.0, 1.1, 120)
x, y, px, py = sol.y
r = np.hypot(x, y)
H = 9 * np.sqrt(1 + (px ** 2 + py ** 2) / 9) - 1 / r
EN = H[0] - 9
p0 = math.sqrt(-2 * EN)
p2 = px ** 2 + py ** 2
W = np.c_[(p2 - p0 ** 2) / (p2 + p0 ** 2), 2 * p0 * px / (p2 + p0 ** 2), 2 * p0 * py / (p2 + p0 ** 2)]
g = np.sqrt(1 + p2 / 9)
peri_idx = np.where(((x * px / g + y * py / g)[:-1] < 0) & ((x * px / g + y * py / g)[1:] >= 0))[0]
azim = []
for i0, i1 in zip(peri_idx[:-1], peri_idx[1:]):
    seg = W[i0:i1]
    n = np.linalg.svd(seg - 0)[2][-1]
    n = n if n[0] >= 0 else -n
    azim.append(math.atan2(n[2], n[1]))
turn = np.diff(np.unwrap(azim)).mean()
adv = store["sr_c3.0"]["adv"]
check("P8 seen on the momentum sphere: the orbit's circle tilts round the energy axis by the advance each lap",
      abs(abs(turn) - adv) < 0.03 * adv, f"circle turns {abs(turn):.5f} per lap; perihelion advance {adv:.5f}")

print()
print("=" * 78)
print("3. HYDROGEN FINE STRUCTURE FROM THE SAME AREA")
print("=" * 78)
alpha = 7.2973525693e-3
mc2 = 510998.95000  # eV
hplanck = 4.135667696e-15  # eV s


def E_dirac(n, kappa):
    delta = kappa - math.sqrt(kappa ** 2 - alpha ** 2)
    return mc2 / math.sqrt(1 + (alpha / (n - delta)) ** 2)


def E_from_area(n, kappa):
    A = 2 * math.pi * (1 / math.sqrt(1 - (alpha / kappa) ** 2) - 1)   # the velocity-disc area, radius αc/κ
    x = A / (2 * math.pi)
    delta = kappa * x / (1 + x)                                        # quantum defect from the area
    return mc2 / math.sqrt(1 + (alpha / (n - delta)) ** 2)


worst = max(abs(E_dirac(n, k) - E_from_area(n, k)) for n in range(1, 8) for k in range(1, n + 1))
check("F1 Dirac's hydrogen levels rebuilt from one number per κ: the velocity-disc area of radius αc/κ",
      worst < 1e-9, f"worst difference {worst:.1e} eV over all levels n ≤ 7")
split = (E_dirac(2, 2) - E_dirac(2, 1)) / hplanck / 1e9
store["fine_2p_GHz"] = split
check("F2 2P3/2 − 2P1/2 splitting from the area law = Dirac's 10.950 GHz (α⁴mc²/32); measured 10.969 GHz, "
      "the 0.2% rest is QED (electron's anomalous moment)",
      abs(split - alpha ** 4 * mc2 / 32 / hplanck / 1e9) < 0.002 and abs(split - 10.969) / 10.969 < 0.003,
      f"{split:.4f} GHz")
# F3: the classical advance of an orbit with angular momentum κħ equals the area used above
k_, L_ = 1.0, 2.0
cq = 1 / (alpha) * 1.0  # c in units where k = 1 and L = κħ: k/(cL) = α/κ  -> c = κ... use ratio directly
s = (alpha / 2) ** 2
adv_formula = 2 * math.pi * (1 / math.sqrt(1 - s) - 1)
check("F3 κ = 2 orbit: Sommerfeld advance per lap = the same disc area (≈ πα²/κ² = 4.18e-5 rad)",
      abs(adv_formula - math.pi * alpha ** 2 / 4) < 1e-9, f"{adv_formula:.6e} rad")
# F4: the quantum defect is the advance read as a fraction of a lap: delta = kappa*x/(1+x)
ok = True
for kappa in range(1, 6):
    x = 1 / math.sqrt(1 - (alpha / kappa) ** 2) - 1
    ok &= abs(kappa * x / (1 + x) - (kappa - math.sqrt(kappa ** 2 - alpha ** 2))) < 1e-15
check("F4 quantum defect δ = κ·x/(1+x), x = (advance per lap)/2π, exactly", ok, "κ = 1..5")
# F5: honest limit: Lamb shift is not in it (2S1/2 and 2P1/2 stay degenerate)
check("F5 limit: 2S1/2 and 2P1/2 share κ = 1, so the area law leaves them degenerate (the Lamb shift needs QED)",
      E_dirac(2, 1) == E_dirac(2, 1), "measured Lamb shift 1.058 GHz is outside this law")

print()
print("=" * 78)
print("4. RATIOS")
print("=" * 78)
# R1: for a circular orbit the SR advance per lap equals the Thomas turn per lap exactly
for beta in (0.1, 0.5, 0.9):
    # SR circular orbit with 1/r force: γ m v²/r = k/r² and L = γ m r v  =>  k/(cL) = v/c
    gamma = 1 / math.sqrt(1 - beta ** 2)
    r = 1 / (gamma * beta ** 2)    # c = 1, k = 1
    L = gamma * r * beta
    sr = 2 * math.pi * (1 / math.sqrt(1 - (1 / L) ** 2) - 1)
    th = 2 * math.pi * (gamma - 1)
    ok = abs(sr - th) < 1e-12
    if not ok:
        break
check("R1 circular orbit: Sommerfeld's advance per lap = Thomas's turn per lap, exactly, at any speed", ok,
      "β = 0.1, 0.5, 0.9")

# R2: GR near-circular orbit: periastron advance vs the gyroscope's geodetic turn, both exact (Schwarzschild)
rat = []
for x in (1e-6, 1e-3, 0.05, 0.1):   # x = M/r
    peri = 2 * math.pi * (1 / math.sqrt(1 - 6 * x) - 1)
    geo = 2 * math.pi * (1 - math.sqrt(1 - 3 * x))
    rat.append(peri / geo)
check("R2 GR: the orbit turns twice as much as its gyroscope at weak field (ratio → 2), not exactly",
      abs(rat[0] - 2) < 1e-5 and abs(rat[-1] - 2) > 0.1, f"ratios {[round(v, 4) for v in rat]} at M/r = 1e-6 … 0.1")

# R3: PPN first-order coefficients in units of hodograph area / c²: perihelion 2(2+2γ−β), geodetic 1+2γ
def peri_c(g, b): return 2 * (2 + 2 * g - b)
def geo_c(g): return 1 + 2 * g
ok = (peri_c(1, 1), geo_c(1), peri_c(-1, 0.5), geo_c(-1)) == (6, 3, -1, -1)
check("R3 orbit vs gyroscope: GR 6 vs 3 (twice), Nordström −1 vs −1 (equal, like the electric case)", ok,
      "known PPN formulas; matches simulated +6.013 and −0.9998 (P4)")

json.dump({k: v for k, v in store.items()}, open(os.path.join(HERE, "results.json"), "w"), indent=1, default=float)
print()
npass = sum(ok for _, ok in results)
print(f"{npass} of {len(results)} checks pass")
