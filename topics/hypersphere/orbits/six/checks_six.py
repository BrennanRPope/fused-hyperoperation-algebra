"""DRAFT checks for hypersphere/orbits/six: an exact split of general relativity's 1 - 6M/r.
Run: python3 checks_six.py   (sympy, numpy, scipy). Writes results_six.json for the figures.

T = (Omega_r / Omega_phi)^2 for a near-circular orbit (radial over angular frequency, squared; any clock).
Advance per lap = 2*pi*(1/sqrt(T) - 1). v = local orbital speed measured by a static observer (c = 1).
"""
import json, math, os
import numpy as np, sympy as sp
from scipy import integrate

HERE = os.path.dirname(os.path.abspath(__file__))
results, store = [], {}
def check(name, ok, detail):
    results.append((name, bool(ok)))
    print(f"[{'PASS' if ok else 'FAIL'}] {name}: {detail}")

print("=" * 78); print("1. EXACT ALGEBRA (sympy)"); print("=" * 78)
x = sp.symbols('x', positive=True)
v2 = x / (1 - 2 * x)                     # Schwarzschild circular orbit, local speed^2
check("S1 GR orbit: 1 - 6M/r = (1 - 2M/r)(1 - 4v²) exactly, v = local orbital speed",
      sp.simplify((1 - 2 * x) * (1 - 4 * v2) - (1 - 6 * x)) == 0, "identity in x = M/r")
check("S2 GR gyroscope (Gravity Probe B page): 1 - 3M/r = (1 - 2M/r)(1 - v²) exactly, same space factor",
      sp.simplify((1 - 2 * x) * (1 - v2) - (1 - 3 * x)) == 0, "identity in x = M/r")
check("S3 the orbit's velocity disc reaches light speed (2v = c) exactly at the innermost stable orbit r = 6M",
      sp.simplify(4 * v2.subs(x, sp.Rational(1, 6))) == 1, "v = c/2 there")
check("S4 the gyroscope's velocity disc reaches light speed (v = c) exactly at the photon sphere r = 3M",
      sp.simplify(v2.subs(x, sp.Rational(1, 3))) == 1, "v = c there")

# General static metric -A dt^2 + B dr^2 + r^2 dphi^2
r, E, L = sp.symbols('r E L', positive=True)
A = sp.Function('A')(r); B = sp.Function('B')(r)
G = E**2 / A - 1 - L**2 / r**2
sol = sp.solve([sp.Eq(G, 0), sp.Eq(sp.diff(G, r), 0)], [E**2, L**2], dict=True)[0]
T_gen = sp.simplify((-sp.diff(G, r, 2).subs({E**2: sol[E**2], L**2: sol[L**2]}) / (2 * B)) / (sol[L**2] / r**4))
vloc2 = r * sp.diff(A, r) / (2 * A)
N_A = 3 + r * sp.diff(A, r, 2) / sp.diff(A, r)
check("S5 ANY static round metric: T = (N - 4v²)/B, N = 3 + rA''/A' (Newton's closure number from the clock rate)",
      sp.simplify(T_gen - (N_A - 4 * vloc2) / B) == 0, "symbolic in A(r), B(r)")

print(); print("=" * 78); print("2. SPECIAL RELATIVITY IN ANY POTENTIAL (numeric, exact to rounding)"); print("=" * 78)
rr, pr, Lm = sp.symbols('r p_r L_m', positive=True)
worst = 0; pts_sr = []
for Vexpr in [-1 / rr, -sp.Rational(1, 3) / rr**2, -1 / rr + sp.Rational(3, 10) / rr**2, sp.log(rr) / 5, rr**2 / 10]:
    H = sp.sqrt(1 + pr**2 + Lm**2 / rr**2) + Vexpr
    Hr = sp.diff(H, rr).subs(pr, 0); Hrr = sp.diff(H, rr, 2).subs(pr, 0); Hpp = sp.diff(H, pr, 2).subs(pr, 0); HL = sp.diff(H, Lm).subs(pr, 0)
    for rv in (0.7, 2.0):
        Ls = sp.nsolve(Hr.subs(rr, rv), Lm, 1.0)
        T = float((Hrr * Hpp / HL**2).subs({rr: rv, Lm: Ls}))
        vv = float(((rr * HL)**2).subs({rr: rv, Lm: Ls}))
        N = float((3 + rr * sp.diff(Vexpr, rr, 2) / sp.diff(Vexpr, rr)).subs(rr, rv))
        worst = max(worst, abs(T - (N - vv))); pts_sr.append((N - vv, T))
check("S6 electric-type orbit in ANY potential V: T = N - v², N = 3 + rV''/V' (Newton's closure number)",
      worst < 1e-10, f"5 potentials × 2 radii, worst error {worst:.1e}")
store["sr_points"] = pts_sr

print(); print("=" * 78); print("3. SIMULATED ORBITS (geodesics, small eccentricity)"); print("=" * 78)
def geodesic_advance(Af, dAf, Bf, dBf, r0, kick=1e-3, laps=6):
    """Integrate a geodesic of -A dt^2 + B dr^2 + r^2 dphi^2 in proper time; return mean periapsis advance per lap
    and the local speed of the circular orbit at r0."""
    A0, dA0 = Af(r0), dAf(r0)
    Lc2 = r0**3 * dA0 / (2 * A0 - r0 * dA0); Ec2 = A0 * (1 + Lc2 / r0**2)
    Lc, Ec = math.sqrt(Lc2), math.sqrt(Ec2)
    def rhs(tau, y):
        rr_, pr_, ph = y
        # H = 1/2 ( -E^2/A + p_r^2/B + L^2/r^2 )
        drdt = pr_ / Bf(rr_)
        dprdt = -0.5 * (Ec**2 * dAf(rr_) / Af(rr_)**2 - pr_**2 * dBf(rr_) / Bf(rr_)**2 - 2 * Lc**2 / rr_**3)
        return [drdt, dprdt, Lc / rr_**2]
    Tlap = 2 * math.pi * r0**2 / Lc
    s = integrate.solve_ivp(rhs, (0, laps * Tlap * 1.3), [r0, kick, 0.0], rtol=1e-12, atol=1e-14, max_step=Tlap / 400)
    rs, prs, phs = s.y
    idx = np.where((prs[:-1] < 0) & (prs[1:] >= 0))[0]      # periapsis: p_r crosses 0 upward
    ph_peri = [phs[i] + (phs[i + 1] - phs[i]) * (-prs[i] / (prs[i + 1] - prs[i])) for i in idx]
    adv = np.diff(ph_peri).mean() - 2 * math.pi
    v2loc = r0 * dA0 / (2 * A0)
    return adv, v2loc

M = 1.0
metrics = {
    "Schwarzschild": (lambda r: 1 - 2 * M / r, lambda r: 2 * M / r**2, lambda r: 1 / (1 - 2 * M / r), lambda r: -2 * M / (r - 2 * M)**2),
    "clock warp only (B = 1)": (lambda r: 1 - 2 * M / r, lambda r: 2 * M / r**2, lambda r: 1.0, lambda r: 0.0),
    "Reissner–Nordström Q = 0.6": (lambda r: 1 - 2 / r + 0.36 / r**2, lambda r: 2 / r**2 - 0.72 / r**3,
                                   lambda r: 1 / (1 - 2 / r + 0.36 / r**2), lambda r: -(2 / r**2 - 0.72 / r**3) / (1 - 2 / r + 0.36 / r**2)**2),
}
sim = []
for name, (Af, dAf, Bf, dBf) in metrics.items():
    worst = 0
    for r0 in (10.0, 14.0, 25.0, 60.0):
        adv, vv = geodesic_advance(Af, dAf, Bf, dBf, r0)
        h = 1e-4
        N = 3 + r0 * ((dAf(r0 + h) - dAf(r0 - h)) / (2 * h)) / dAf(r0)
        T = (N - 4 * vv) / Bf(r0)
        pred = 2 * math.pi * (1 / math.sqrt(T) - 1)
        worst = max(worst, abs(adv - pred) / abs(pred)); sim.append((name, r0, vv, adv, pred))
    check(f"S7 {name}: simulated advance per lap = 2π(1/√T − 1), T = (N − 4v²)/B", worst < 2e-3,
          f"r = 10, 14, 25, 60 M; worst relative error {worst:.1e} (kick 1e-3, so small eccentricity)")
store["geodesics"] = sim

# Schwarzschild: the two factors plotted against M/r
xs = np.linspace(0.001, 1 / 3 - 1e-3, 300)
store["curve"] = {"x": xs.tolist()}

print(); print("=" * 78); print("4. FIRST ORDER AND WHAT IT MEANS"); print("=" * 78)
# first-order coefficients in units of the hodograph area / c^2 = pi v^2 (circular): space 2, clock 4, total 6
x0 = 1e-7; vv = x0 / (1 - 2 * x0)
space = 2 * math.pi * (1 / math.sqrt(1 - 2 * x0) - 1) / (math.pi * vv)
clock = 2 * math.pi * (1 / math.sqrt(1 - 4 * vv) - 1) / (math.pi * vv)
total = 2 * math.pi * (1 / math.sqrt(1 - 6 * x0) - 1) / (math.pi * vv)
check("S8 the 6 = 2 (space curvature, as for the gyroscope) + 4 (velocity disc of radius 2v, four times Thomas's area)",
      abs(space - 2) < 1e-5 and abs(clock - 4) < 1e-5 and abs(total - 6) < 1e-5, f"{space:.5f} + {clock:.5f} = {total:.5f}")
x_m = 1.32712440018e20 / (299792458.0**2 * 5.7909e10 * (1 - 0.2056**2))   # M/p for Mercury (p = semi-latus rectum)
adv_m = 2 * math.pi * (1 / math.sqrt((1 - 2 * x_m) * (1 - 4 * x_m / (1 - 2 * x_m))) - 1)
store["mercury_arcsec_century"] = adv_m * 206264.806 * 36525 / 87.969
check("S9 Mercury from the split: (1 − 2M/p)(1 − 4v²) per orbit gives 42.98″ per century (first-order test only)",
      abs(store["mercury_arcsec_century"] - 42.98) < 0.05, f"{store['mercury_arcsec_century']:.2f}″/century")

json.dump(store, open(os.path.join(HERE, "results_six.json"), "w"), indent=1, default=float)
print(); npass = sum(ok for _, ok in results); print(f"{npass} of {len(results)} checks pass")
