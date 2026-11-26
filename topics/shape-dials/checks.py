"""Checks for the shape-dials page (DRAFT). Run: python3 checks.py > checks_output.txt"""
import sys, os, time, json
import numpy as np
from scipy.optimize import minimize_scalar
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(HERE, "..", "simplex"))
from facedial import *
from shapes_catalog import REGULAR, NONREGULAR
from wythoff import CHAINS, wythoff_polytope, roots_for, chamber_corners
import shape_dial as SD

T0 = time.time()
results, fails = [], 0
RES = {}


def check(cid, desc, ok, detail=""):
    global fails
    results.append((cid, ok))
    fails += (not ok)
    print(f"[{'PASS' if ok else 'FAIL'}] {cid}  {desc}" + (f"  ({detail})" if detail else ""))


ALL = {**REGULAR, **NONREGULAR}
DIALS = {nm: FaceDial(V, nm) for nm, V in ALL.items()}
PS = [2, 3, 4, 8, 16, 64, 256]

# ---------------------------------------------------------------- A. round at p = 2, exact at p = inf
worst2, worstinf = 0, 0
for nm, f in DIALS.items():
    U = sphere_dirs(f.d, 2000, seed=1)
    worst2 = max(worst2, roundness_gap(f, 2, U))
    worstinf = max(worstinf, np.max(np.abs(f.radial(U, np.inf) / f.radial_polytope(U) - 1)))
check("A1", f"p = 2 is exactly round for all {len(DIALS)} shapes (after the analytic-centre whitening)", worst2 < 1e-10,
      f"worst radius spread {worst2:.1e}")
check("A2", "p = inf is exactly the polytope for every shape", worstinf < 1e-10, f"worst {worstinf:.1e}")
reg_tight = max(DIALS[nm].tightness for nm in REGULAR)
check("A3", "regular shapes need no whitening: their scaled normals are already a tight frame", reg_tight < 1e-12,
      f"worst {reg_tight:.1e}")
nr_whiten = {nm: DIALS[nm].tightness for nm in NONREGULAR}
need = [nm for nm, v in nr_whiten.items() if v > 1e-9]
check("A4", "permutohedra, Birkhoff B3 and hypersimplices are tight too; associahedra, zonotope, cyclic need whitening",
      set(need) == {"associahedron (3D, Loday)", "associahedron (4D, Loday)", "zonotope (3D, 5 random zones)",
                    "cyclic polytope C(7,3)", "cyclic polytope C(8,4)"}, ", ".join(need))
sumN = max(f.sumN for f in DIALS.values())
check("A5", "scaled normals sum to zero at the analytic centre for every shape", sumN < 1e-9, f"{sumN:.1e}")

# ---------------------------------------------------------------- B. offset a, symmetry, corners
sym_ok = all((f.a == 0) == f.symmetric for f in DIALS.values())
check("B1", "offset a = 0 exactly for the centrally symmetric shapes, a > 0 otherwise", sym_ok)
simp = [DIALS[k].a for k in ["3-gon", "tetrahedron", "5-cell", "5-simplex", "6-simplex"]]
check("B2", "simplex in d dimensions has a = (d - 1)/2", np.allclose(simp, [0.5, 1, 1.5, 2, 2.5]), str(simp))
odd = {k: round(DIALS[k].a, 6) for k in ["5-gon", "7-gon"]}
check("B3", "odd polygons: a = (sec(pi/n) - 1)/2 (pentagon 0.1180, heptagon 0.0549)",
      all(abs(DIALS[f"{n}-gon"].a - (1 / np.cos(np.pi / n) - 1) / 2) < 1e-9 for n in (3, 5, 7)), str(odd))
cfix = 0
for k in ["3-gon", "tetrahedron", "5-cell", "6-simplex"]:
    f = DIALS[k]
    for p in [2, 3, 7.5, 40]:
        cfix = max(cfix, np.max(np.abs(f.gauge_value(f.Vw, p) - (f.a + 1))))
check("B4", "simplex: every corner stays on the boundary at every p (p = 2 is the circumscribed ball)", cfix < 1e-9,
      f"{cfix:.1e}")

# ---------------------------------------------------------------- C. nesting and rate
nest, bound_bad, worst_ratio = 0, 0, {}
for nm, f in DIALS.items():
    U = sphere_dirs(f.d, 600, seed=2)
    rP = f.radial_polytope(U)
    prev = None
    for p in PS:
        r = f.radial(U, p)
        if prev is not None:
            nest = max(nest, np.max(r - prev))
        prev = r
        q = r / rP
        hi = 1 + (f.m ** (1 / p) - 1) * (f.a + 1)
        if q.min() < 1 - 1e-9 or q.max() > hi + 1e-9:
            bound_bad += 1
    worst_ratio[nm] = [float(np.max(f.radial(U, p) / rP) - 1) for p in PS]
RES["sharpness_gap"] = worst_ratio
check("C1", "shapes shrink as p grows (K_q inside K_p for q > p), every shape", nest < 1e-10, f"{nest:.1e}")
check("C2", "rate: r_P <= r_p <= (1 + (m^(1/p) - 1)(a + 1)) r_P in every direction, every shape, p = 2..256",
      bound_bad == 0, f"{bound_bad} violations")

# ---------------------------------------------------------------- D. the cube gives back the existing p dial
dev = 0
for d in (2, 3, 4):
    f = FaceDial(cube(d))
    U = sphere_dirs(d, 800, seed=3)
    for p in [1, 1.25, 1.5, 2, 3, 8, 30]:
        rf = f.radial(U, p)
        Uo = f.to_original(U)
        rl = 1 / np.sum(np.abs(Uo) ** p, axis=1) ** (1 / p)
        q = rf / rl
        dev = max(dev, q.max() / q.min() - 1)
check("D1", "aimed at the cube, the dial is the existing roundness dial p (l^p ball) up to size, for p in [1, 30]",
      dev < 1e-9, f"{dev:.1e}")
f = FaceDial(cube(3))
U = sphere_dirs(3, 4000, seed=4)
r1 = f.radial(U, 1)
Uo = f.to_original(U)
q = r1 / (1 / np.abs(Uo).sum(1))
check("D2", "and its p = 1 end is the cross-polytope (octahedron), as before", q.max() / q.min() - 1 < 1e-9)

# ---------------------------------------------------------------- E. Dikin ellipsoid at p = 2 (unwhitened)
de = 0
for nm in ["associahedron (3D, Loday)", "cyclic polytope C(7,3)", "zonotope (3D, 5 random zones)"]:
    f = DIALS[nm]
    U = sphere_dirs(3, 500, seed=5)
    r = f.radial(U, 2)
    Y = (r[:, None] * U) @ f.Winv.T
    de = max(de, np.max(np.abs(np.einsum("ki,ij,kj->k", Y, f.H, Y) / (f.m * (2 * f.a + 1)) - 1)))
check("E1", "p = 2 is the Dikin ellipsoid y'Hy = m(2a+1) of the log barrier at the analytic centre", de < 1e-9,
      f"{de:.1e}")

# ---------------------------------------------------------------- F. the simplex test
tri = DIALS["3-gon"]
th, U2 = circle_dirs(6000)
O = np.zeros((1, 2))


def body(r):
    return SD.Body(r[:, None] * U2)


def norm(r):
    return r / r.mean()


lim = {}
beta = 3.0
rl = lse_radial(tri, U2, beta)
for a in [2, 8, 32, 128, 512]:
    fa = FaceDial(regular_polygon(3), a=a)
    lim[a] = float(np.max(np.abs(fa.radial(U2, beta * a) / rl - 1)))
RES["limit"] = lim
rates = [lim[a] / lim[4 * a] for a in [2, 8, 32, 128]]
check("F1", "temperature route = face-power with offset a -> inf and p = beta*a (error falls like 1/a)",
      lim[512] < 1e-3 and all(3.3 < x < 4.7 for x in rates), f"errors {[round(v, 5) for v in lim.values()]}")
kh = {b: float(SD.khat(body(lse_radial(tri, U2, b)), O)[0]) for b in [1, 3, 10]}
check("F2", "our temperature route reproduces the simplex page's centre turns (-0.69, -0.37, -0.14 at beta 1, 3, 10)",
      all(abs(kh[b] - v) < 0.01 for b, v in zip([1, 3, 10], [-0.69, -0.37, -0.14])), str({k: round(v, 3) for k, v in kh.items()}))
best = {}
for p in [3, 4, 6, 10, 20, 50]:
    rp = norm(tri.radial(U2, p))
    fun = lambda lb: np.max(np.abs(norm(lse_radial(tri, U2, np.exp(lb))) - rp))
    res = minimize_scalar(fun, bounds=(-6, 8), method="bounded")
    best[p] = (float(np.exp(res.x)), float(res.fun))
RES["fp_vs_temperature"] = best
mx = max(v[1] for v in best.values())
check("F3", "at the triangle's own offset (a = 1/2) face-power stays within 1.2% of the temperature route but is not "
            "equal to it", 1e-3 < mx < 0.012, f"closest-match radius error at most {mx:.4f}")
# blend and trim routes: best match against face-power
def shape_r(bodyobj):
    B = bodyobj.boundary
    ang = np.arctan2(B[:, 1], B[:, 0])
    rr = np.linalg.norm(B, axis=1)
    o = np.argsort(ang)
    return np.interp(th, np.mod(ang[o], 2 * np.pi), rr[o], period=2 * np.pi)
TRI_ROT = None
fp_curve = {p: norm(tri.radial(U2, p)) for p in np.exp(np.linspace(np.log(2.0), np.log(400), 120))}
# align orientation: our triangle's corner directions vs the simplex page's (corner at 90 deg)
corner_ang = np.arctan2(tri.Vw[:, 1], tri.Vw[:, 0])
shift = (np.pi / 2 - corner_ang[0])
def fp_match(rt):
    out = 1e9
    for p, rp in fp_curve.items():
        rr = np.interp(np.mod(th - shift, 2 * np.pi), th, rp, period=2 * np.pi)
        out = min(out, np.max(np.abs(rr - rt)))
    return out
blend = {t: fp_match(norm(shape_r(SD.minkowski(t)))) for t in [0.25, 0.5, 0.75]}
trim = {s: fp_match(norm(shape_r(SD.cut_triangle(s)))) for s in [0.6, 0.75, 0.9]}
RES["blend_trim_distance"] = {"blend": blend, "trim": trim}
check("F4", "the blend and trim routes are separate routes too (closest face-power match off by 1.4% to 5.7%)",
      min(list(blend.values()) + list(trim.values())) > 0.01,
      f"blend {[round(v, 3) for v in blend.values()]}, trim {[round(v, 3) for v in trim.values()]}")
kfp = {}
for p in [2, 2.5, 3, 4, 6, 10, 20, 50, 200]:
    r = tri.radial(U2, p)
    kfp[p] = (float(r.min() / r.max()), float(SD.khat(body(r), O)[0]))
RES["khat_facepower"] = kfp
vals = [v[1] for v in kfp.values()]
check("F5", "centre turn along face-power: exactly -1 at p = 2, rises monotonically toward 0", abs(vals[0] + 1) < 1e-3
      and all(np.diff(vals) > 0) and vals[-1] > -0.03, str([round(v, 3) for v in vals]))

# ---------------------------------------------------------------- G. curved sides: the elliptope (spectral version)
def ell_matrix(v):
    x, y, z = v
    return np.array([[1, x, y], [x, 1, z], [y, z, 1.0]])


def spectral_radial(u, p, a=0.5):
    lo, hi = 0.0, 4.0
    for _ in range(100):
        mid = (lo + hi) / 2
        lam = np.linalg.eigvalsh(ell_matrix(mid * u))
        t = a + (1 - lam)
        val = np.abs(t).max() if np.isinf(p) else np.mean(np.abs(t) ** p) ** (1 / p)
        lo, hi = (mid, hi) if val <= a + 1 else (lo, mid)
    return (lo + hi) / 2


def ell_radial(u):
    lo, hi = 0.0, 4.0
    for _ in range(100):
        mid = (lo + hi) / 2
        lo, hi = (mid, hi) if np.linalg.eigvalsh(ell_matrix(mid * u)).min() >= 0 else (lo, mid)
    return (lo + hi) / 2


Ue = sphere_dirs(3, 300, seed=6)
r2 = np.array([spectral_radial(u, 2) for u in Ue])
check("G1", "elliptope: eigenvalues 1 - lambda play the faces; p = 2 is the ball of radius sqrt 3", np.max(np.abs(r2 - 3 ** 0.5)) < 1e-9)
corners = np.array([[1, 1, 1], [1, -1, -1], [-1, 1, -1], [-1, -1, 1.0]])
cf = max(abs(np.mean(np.abs(0.5 + 1 - np.linalg.eigvalsh(ell_matrix(c))) ** p) ** (1 / p) - 1.5) for c in corners for p in [2, 5, 30])
check("G2", "its four corners (rank-one correlation matrices) stay on the boundary at every p", cf < 1e-9, f"{cf:.1e}")
rinf = np.array([spectral_radial(u, np.inf) for u in Ue])
re_ = np.array([ell_radial(u) for u in Ue])
check("G3", "p = inf is exactly the elliptope (curved sides and four corners)", np.max(np.abs(rinf / re_ - 1)) < 1e-9)
gaps = {p: float(np.max(np.array([spectral_radial(u, p) for u in Ue]) / re_ - 1)) for p in [4, 16, 64]}
RES["elliptope_gap"] = gaps
check("G4", "and it sharpens toward the elliptope as p grows", gaps[4] > gaps[16] > gaps[64] > 0, str({k: round(v, 4) for k, v in gaps.items()}))

# ---------------------------------------------------------------- H. qutrits
rng = np.random.default_rng(8)
qok = 0
for _ in range(200):
    A = rng.standard_normal((3, 3)) + 1j * rng.standard_normal((3, 3))
    Hh = A + A.conj().T
    Hh -= np.trace(Hh).real / 3 * np.eye(3)
    Hh /= np.linalg.norm(Hh)
    lo, hi = 0, 5
    for _ in range(80):
        mid = (lo + hi) / 2
        rho = np.eye(3) / 3 + mid * Hh
        lam = np.linalg.eigvalsh(rho)
        t = 0.5 + 1 - 3 * lam
        lo, hi = (mid, hi) if np.mean(t ** 2) ** 0.5 <= 1.5 else (lo, mid)
    rho = np.eye(3) / 3 + lo * Hh
    qok = max(qok, abs(np.trace(rho @ rho).real - 1))
check("H1", "qutrit density matrices: the p = 2 shape is exactly the pure-state sphere tr(rho^2) = 1", qok < 1e-9, f"{qok:.1e}")

# ---------------------------------------------------------------- W. Wythoff's kaleidoscope: the shape-choice dial
names = {("B3", 0): (8, "cube"), ("B3", 1): (12, "cuboctahedron"), ("B3", 2): (6, "octahedron"),
         ("A3", 0): (4, "tetrahedron"), ("A3", 1): (6, "octahedron"), ("A3", 2): (4, "tetrahedron (dual)"),
         ("H3", 0): (20, "dodecahedron"), ("H3", 1): (30, "icosidodecahedron"), ("H3", 2): (12, "icosahedron"),
         ("B4", 0): (16, "tesseract"), ("B4", 2): (24, "24-cell"), ("B4", 3): (8, "16-cell"),
         ("F4", 0): (24, "24-cell"), ("H4", 0): (600, "120-cell"), ("H4", 3): (120, "600-cell"),
         ("A4", 1): (10, "hypersimplex D(2,5)"), ("A5", 2): (20, "hypersimplex D(3,6)")}
wy_bad = [f"{k}" for k, (nv, _) in names.items() if len(wythoff_polytope(k[0], k[1])) != nv]
check("W1", "depth t walks each family through its regular members: B3 cube/cuboctahedron/octahedron, H3, A3, B4 "
            "(24-cell at t = 2), F4, H4, A (hypersimplices)", not wy_bad, ", ".join(wy_bad))
poly_ok = all(len(wythoff_polytope("I2(n)", 0.5, n)) == 2 * n and len(wythoff_polytope("I2(n)", 0, n)) == n for n in range(3, 13))
check("W2", "in 2D the depth walks n-gon -> regular 2n-gon -> n-gon; the side count n itself stays a whole number", poly_ok)
w2 = 0
for nm in ["A3", "B3", "H3", "A4", "B4", "F4"]:
    d = len(CHAINS[nm]) + 1
    for t in np.linspace(0, d - 1, 9):
        f = FaceDial(wythoff_polytope(nm, t))
        U = sphere_dirs(f.d, 300, seed=9)
        w2 = max(w2, roundness_gap(f, 2, U), f.tightness)
check("W3", "every point on every depth path is round at p = 2 with no whitening (all paths meet at one ball)", w2 < 1e-9,
      f"{w2:.1e}")
perm = FaceDial(wythoff_polytope_generic := __import__("wythoff").orbit(roots_for("A3"),
                (lambda C: (C.sum(0) / np.linalg.norm(C.sum(0))))(chamber_corners(roots_for("A3")))))
check("W4", "the permutohedron is the A-family seed at the chamber's centre (all mirrors ringed)", len(perm.Vw) == 24 and perm.m == 14)
avals = {t: round(FaceDial(wythoff_polytope("A3", t)).a, 4) for t in np.linspace(0, 2, 9)}
check("W5", "offset a moves continuously with depth on the simplex family (1 at tetrahedron, 0 at octahedron)",
      avals[0.0] == 1.0 and avals[1.0] == 0.0, str(avals))

print(f"\n{len(results) - fails} of {len(results)} pass, {time.time() - T0:.1f} s")
with open(os.path.join(HERE, "results.json"), "w") as fh:
    json.dump(RES, fh, indent=1, default=str)
