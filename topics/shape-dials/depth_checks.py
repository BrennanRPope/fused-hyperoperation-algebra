"""Is Wythoff depth a new dial, or derived?  (DRAFT, Nov 26 2026)

Brennan: "make sure it's really a new dial ... in as rigorous a way as you can ... I also don't want things to be
treated as fundamental if they really are derived."

Checks (all numeric; the proofs are written in README section 6):
  M  depth is a Minkowski sum of the family's corner shapes, weighted by chances (odds reading)
  H  the same shape on the face side: corner face sets, scaled by similarity ratios, intersected
  O  the odds reading: complement, rescaling = a constant odds factor, centre of the chamber = permutohedron
  G  the same weighted sums reach non-kaleidoscope shapes we use (Postnikov): associahedron, zonotope, permutohedron
  I  independence of p and depth: what each reaches that the other cannot
  S  the simplex gap: exact only as a -> infinity; at a = 1/2 a real difference (never a circle at finite beta)
Run: python3 depth_checks.py  ->  depth_checks_output.txt, depth_results.json
"""
import sys, os, json, time, itertools
import numpy as np
from scipy.spatial import ConvexHull, HalfspaceIntersection
from scipy.special import iv

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from facedial import FaceDial, lse_radial, circle_dirs, regular_polygon, loday_associahedron, permutohedron, \
    hypersimplex, zonotope, sphere_dirs, cube, cross_polytope, simplex, icosahedron, dodecahedron
from wythoff import CHAINS, roots_for, chamber_corners, orbit, seed_at, wythoff_polytope

T0 = time.time()
OUT, RES = [], {}


def check(cid, desc, ok, detail=""):
    line = f"[{'PASS' if ok else 'FAIL'}] {cid}: {desc}" + (f"  ({detail})" if detail else "")
    print(line); OUT.append(line); RES[cid] = {"ok": bool(ok), "desc": desc, "detail": detail}


def h(V, U):
    """Support function of conv(V) in directions U."""
    return (np.atleast_2d(U) @ np.asarray(V).T).max(1)


def dirs(d, k=3000, seed=1):
    X = np.random.default_rng(seed).normal(size=(k, d))
    return X / np.linalg.norm(X, axis=1, keepdims=True)


FAMS = ["A3", "B3", "H3", "A4", "B4", "F4", "H4"]

# ------------------------------------------------------------------ M: Minkowski sum of corner shapes
worst = 0.0
rng = np.random.default_rng(3)
for fam in FAMS:
    R = roots_for(fam)
    C = chamber_corners(R)                      # unit chamber rays
    corner_orbits = [orbit(R, c) for c in C]    # corner shapes, circumradius 1
    U = dirs(len(C), 2000, seed=len(fam))
    trials = 2 if fam == "H4" else 4
    for _ in range(trials):
        lam = rng.dirichlet(np.ones(len(C)))
        x = lam @ C
        P = orbit(R, x)
        lhs = h(P, U)
        rhs = sum(l * h(V, U) for l, V in zip(lam, corner_orbits))
        worst = max(worst, np.abs(lhs - rhs).max())
check("M1", "For every family, the kaleidoscope shape of any seed is the Minkowski sum of the corner shapes, weighted "
            "by the seed's barycentric coordinates (support functions agree)", worst < 1e-10,
      f"7 families A3 B3 H3 A4 B4 F4 H4, random seeds, 2000 directions each, max error {worst:.1e}")

worst = 0.0
for fam in ["A3", "B3", "H3", "B4", "F4"]:
    R = roots_for(fam); C = chamber_corners(R); U = dirs(len(C), 1500, seed=7)
    corner_orbits = [orbit(R, c) for c in C]
    for depth in np.linspace(0, len(C) - 1, 9):
        k = min(int(np.floor(depth)), len(C) - 2); f = depth - k
        raw = (1 - f) * C[k] + f * C[k + 1]
        P = orbit(R, raw / np.linalg.norm(raw))
        rhs = ((1 - f) * h(corner_orbits[k], U) + f * h(corner_orbits[k + 1], U)) / np.linalg.norm(raw)
        worst = max(worst, np.abs(h(P, U) - rhs).max())
check("M2", "Along the depth path, the shape at depth k + f is (1 - f)·corner_k + f·corner_(k+1), Minkowski sum, up to "
            "size", worst < 1e-10, f"A3 B3 H3 B4 F4, 9 depths each, max error {worst:.1e}")

# M3: the dominant-chamber lemma behind M1: for x, u both in the closed chamber, max_w <w x, u> = <x, u>
worst = 0.0
for fam in FAMS[:6]:
    R = roots_for(fam); C = chamber_corners(R)
    for _ in range(20):
        x = rng.dirichlet(np.ones(len(C))) @ C
        u = rng.dirichlet(np.ones(len(C))) @ C
        worst = max(worst, abs(h(orbit(R, x), u[None])[0] - x @ u))
check("M3", "The lemma behind M1: when the seed and the direction both lie in the chamber, the identity wins the "
            "maximum, so the support function is linear in the seed (Kostant convexity, standard)", worst < 1e-10,
      f"6 families, 20 random pairs each, max error {worst:.1e}")

# ------------------------------------------------------------------ H: face side
def radial_V(V, U):
    hull = ConvexHull(V)
    A, b = hull.equations[:, :-1], -hull.equations[:, -1]
    return 1.0 / (U @ A.T / b).max(1)


worst = 0.0
for fam in ["A3", "B3", "H3", "A4", "B4", "F4"]:
    R = roots_for(fam); C = chamber_corners(R); U = dirs(len(C), 1500, seed=11)
    normal_orbits = [orbit(R, c) for c in C]       # facet normals of every kaleidoscope shape lie in these orbits
    for _ in range(3):
        x = rng.dirichlet(np.ones(len(C))) @ C
        P = orbit(R, x)
        rP = radial_V(P, U)
        rI = np.full(len(U), np.inf)
        for j, Nj in enumerate(normal_orbits):
            s = x @ C[j]                         # offset of the j-th face set = similarity ratio of the polar corner
            T = U @ Nj.T
            rI = np.minimum(rI, s / T.max(1))
        worst = max(worst, np.abs(rP - rI).max())
check("H1", "Face side: the same shape is the intersection of the corner shapes' polars, each scaled by the ratio "
            "<seed, corner> (cube ∩ scaled octahedron = truncated cube)", worst < 1e-9,
      f"A3 B3 H3 A4 B4 F4, random seeds, max radial error {worst:.1e}")

# ------------------------------------------------------------------ O: the odds reading
R = roots_for("B3"); C = chamber_corners(R)
A_, B_ = orbit(R, C[0]), orbit(R, C[1])
U = dirs(3, 1500, seed=5)
f = 0.3
t = f / (1 - f)
hP = (1 - f) * h(A_, U) + f * h(B_, U)
hQ = (1 / (1 + t)) * h(A_, U) + (t / (1 + t)) * h(B_, U)
check("O1", "The depth fraction f is the chance reading t/(1 + t) of an odds t = f/(1 - f) between two corners",
      np.abs(hP - hQ).max() < 1e-12, "B3, f = 0.3, t = 3/7")
alpha, beta_ = 2.0, 0.5         # resize the two corner shapes
t2 = t * alpha / beta_
hR = (1 / (1 + t2)) * alpha * h(A_, U) + (t2 / (1 + t2)) * beta_ * h(B_, U)
ratio = hR / hP
check("O2", "Resizing the corner shapes changes no shape, only multiplies the odds by a fixed factor (a Bayes factor): "
            "depth fractions depend on a size convention, odds ratios do not", np.ptp(ratio) < 1e-12,
      f"corners scaled 2 and 0.5, odds 3/7 -> {t2:.4f}, same shape (ratio spread {np.ptp(ratio):.1e})")
# complement: odds 1/t swaps the corners
hC = (t / (1 + t)) * h(A_, U) + (1 / (1 + t)) * h(B_, U)
hD = (1 - (1 - f)) * h(A_, U) + (1 - f) * h(B_, U)
check("O3", "The odds complement 1/t is the mirror of the depth path (f -> 1 - f)", np.abs(hC - hD).max() < 1e-12)
# centre of the A_(n-1) chamber, equal odds on every corner: the permutohedron
worst = 0.0
for n in [3, 4, 5, 6]:
    Pn = permutohedron(n); U = dirs(n, 1500, seed=n)
    rhs = sum(h(hypersimplex(k, n), U) for k in range(1, n + 1))
    worst = max(worst, np.abs(h(Pn, U) - rhs).max())
check("O4", "Equal odds on every corner of the simplex family give the permutohedron: Perm(n) = Δ(1,n) + Δ(2,n) + … + "
            "Δ(n,n), the hypersimplices summed", worst < 1e-10, f"n = 3..6, max error {worst:.1e}")

# ------------------------------------------------------------------ G: Postnikov's weighted sums reach non-kaleidoscope shapes
worst = 0.0
for n in [3, 4, 5]:
    L = loday_associahedron(n)          # coordinate k: (#leaves left)(#leaves right)
    U = dirs(n, 1500, seed=20 + n)
    # Postnikov: Loday's associahedron = sum over intervals [i, j] of the simplex faces Δ_[i,j]
    rhs = np.zeros(len(U))
    for i in range(n):
        for j in range(i, n):
            E = np.eye(n)[i:j + 1]
            rhs += h(E, U)
    worst = max(worst, np.abs(h(L, U) - rhs).max())
check("G1", "The associahedron (not a kaleidoscope shape) is a Minkowski sum of simplex faces with odds 1 each "
            "(Postnikov 2009; Loday 2004)", worst < 1e-10, f"Loday realisation, n = 3..5, max error {worst:.1e}")
Gm = np.random.default_rng(7).normal(size=(5, 3))
Z = zonotope(Gm); U = dirs(3, 1500, seed=9)
rhs = sum(h(np.array([-g, g]), U) for g in Gm)
check("G2", "Every zonotope is a weighted Minkowski sum of segments", np.abs(h(Z, U) - rhs).max() < 1e-10)

# ------------------------------------------------------------------ I: independence of p and depth
# I1: at 1 < p < inf, face-power shapes are strictly convex (no flat pieces), so never a polytope.
# Hessian of Σ|a + n_i·x|^p = Σ p(p-1)|t_i|^(p-2) n_i n_iᵀ: positive definite wherever every weight is > 0 and the
# normals span. Check both facts on the boundary (weights finite and positive, normals of full rank).
ok_all, min_w, ranks = True, np.inf, []
for name, V in [("cube", cube(3)), ("octahedron", cross_polytope(3)), ("tetrahedron", simplex(3)),
                ("icosahedron", icosahedron()), ("dodecahedron", dodecahedron())]:
    fd = FaceDial(V, name)
    ranks.append(np.linalg.matrix_rank(fd.N))
    U = dirs(3, 400, seed=2)
    for p in [1.5, 3, 10, 40]:
        r = fd.radial(U, p)
        T = np.abs(fd.a + (r[:, None] * U) @ fd.N.T)
        w = T.min()                                  # weights |t|^(p-2) > 0 iff no t vanishes (p < 2) / finite
        min_w = min(min_w, w)
        ok_all &= bool(np.all(np.isfinite(T)) and w > 0)
check("I1", "For 1 < p < ∞ every face-power shape is strictly convex: the Hessian of Σ|a + n_i·x|^p is a positive "
            "sum of n_i n_iᵀ and the normals span, so no flat facet exists and no finite p gives a polytope",
      ok_all and all(rk == 3 for rk in ranks),
      f"5 solids, p = 1.5, 3, 10, 40: normals rank 3 for all, smallest face value on the boundary {min_w:.2e} > 0")


def p1_body(fd):
    """Face-power at p = 1 is a polytope: { w : Σ |a + N_i.w| <= m (a+1) }. Return its vertices (whitened)."""
    m, a, N = fd.m, fd.a, fd.N
    if fd.symmetric:
        # polar of the zonotope (1/m) Σ [-N_i, N_i]; with ± pairs, keep one of each pair (doubling the segment)
        keep, seen = [], []
        for v in N:
            if not any(np.allclose(v, -s) or np.allclose(v, s) for s in seen):
                seen.append(v); keep.append(v)
        G = np.array(keep) * 2 / m
        Zv = zonotope(G)
        hull = ConvexHull(Zv)
        A, b = hull.equations[:, :-1], -hull.equations[:, -1]
        return np.unique(np.round(A / b[:, None], 9), axis=0)
    S = np.array(list(itertools.product([-1, 1], repeat=m)), float)
    A = S @ N
    b = m * (a + 1) - a * S.sum(1)
    hs = HalfspaceIntersection(np.column_stack([A, -b]), np.zeros(N.shape[1]))
    return np.unique(np.round(hs.intersections, 9), axis=0)


def fvec(V):
    hull = ConvexHull(V)
    # merge coplanar triangles into facets
    eq = np.round(hull.equations / np.linalg.norm(hull.equations[:, :-1], axis=1)[:, None], 7)
    facets = {}
    for e, simp in zip(map(tuple, eq), hull.simplices):
        facets.setdefault(e, set()).update(simp)
    sizes = sorted(len(s) for s in facets.values())
    return len(hull.vertices), len(facets), sizes


# known depth points: whole and in-between shapes of A3, B3, H3
targets = {}
for fam in ["A3", "B3", "H3"]:
    for depth in [0, 0.5, 1, 1.5, 2]:
        targets[(fam, depth)] = fvec(wythoff_polytope(fam, depth))
listed = {"tetrahedron": simplex(3), "cube": cube(3), "octahedron": cross_polytope(3),
          "icosahedron": icosahedron(), "dodecahedron": dodecahedron(),
          "cuboctahedron": wythoff_polytope("B3", 1), "icosidodecahedron": wythoff_polytope("H3", 1)}
p1 = {}
for nm, V in listed.items():
    fd = FaceDial(V, nm)
    p1[nm] = fvec(p1_body(fd))
hits = {nm: [k for k, v in targets.items() if v[:2] == p1[nm][:2] and v[2] == p1[nm][2]] for nm in listed}
inbetween = [k for k in targets if k[1] % 1 != 0]
hit_inbetween = [k for nm in hits for k in hits[nm] if k in inbetween]
check("I2", "No setting of p on any listed solid gives an in-between (truncated) depth shape: p in (1, ∞) is strictly "
            "convex (I1), p = ∞ gives the solid itself, and p = 1 gives the polar of a zonotope, whose faces are "
            "centrally symmetric, while every truncation's polar has triangular faces",
      len(hit_inbetween) == 0,
      "p = 1 bodies: " + "; ".join(f"{nm} -> {p1[nm][0]} corners, {p1[nm][1]} faces" for nm in listed))
coinc = {nm: hits[nm] for nm in hits if hits[nm]}
RES["I2"]["p1_matches"] = {nm: [f"{k[0]} depth {k[1]}" for k in v] for nm, v in coinc.items()}
check("I3", "But p = 1 does land on some whole-depth shapes (a crossing, not containment): it gives the polar of a "
            "zonotope, and three of those are kaleidoscope shapes (cube -> octahedron, octahedron -> cuboctahedron, "
            "dodecahedron -> icosidodecahedron)",
      any(coinc.values()), "; ".join(f"{nm} at p = 1 has the combinatorics of {', '.join(f'{a} depth {b}' for a, b in v)}"
                                    for nm, v in coinc.items()))
check("I4", "Depth never reaches a rounded shape: every depth shape is a polytope (the convex hull of a finite orbit), "
            "while p = 2 is a ball for every shape", True, "by construction; checked in W2 of checks.py that p = 2 is round")

# I5: the depth path is a line in the chamber of seeds; the odds reach the whole chamber (omnitruncations included)
R = roots_for("B3"); C = chamber_corners(R)
omni = orbit(R, C.mean(0)); trunc_on_path = orbit(R, seed_at(R, 0.5))
n_omni, n_path = len(omni), len(trunc_on_path)
check("I5", "Odds on all corners reach shapes no depth setting reaches: equal odds on the cube family's three corners "
            "give the truncated cuboctahedron (48 corners), off every depth path (24 corners at most)",
      n_omni == 48 and all(len(wythoff_polytope("B3", d)) <= 24 for d in np.linspace(0, 2, 21)),
      f"omnitruncation {n_omni} corners; depth path at most {max(len(wythoff_polytope('B3', d)) for d in np.linspace(0, 2, 21))}")

# ------------------------------------------------------------------ S: the simplex gap, settled
tri = FaceDial(regular_polygon(3), "triangle")
th, U2 = circle_dirs(2048)
r2 = tri.radial(U2, 2)
check("S1", "Face-power on the triangle at p = 2 is exactly a circle (the circumscribed one)",
      np.ptp(r2) / r2.mean() < 1e-12, f"radius spread {np.ptp(r2) / r2.mean():.1e}")
# Bessel argument: on any circle |x| = rho, F = Σ exp(beta rho cos(θ - θ_i)) = 3[I0 + 2 I3 cos3θ' + 2 I6 cos6θ' + ...]
z = np.logspace(-4, 3, 200)
check("S2", "The temperature route is never exactly a circle at finite β: on any circle its function has a cos 3θ term "
            "with coefficient 6·I₃(βρ) > 0 (Bessel), so a level set cannot be a circle",
      np.all(iv(3, z) > 0), "I₃(z) > 0 for all z > 0 (series with positive terms); checked on 200 values 1e-4..1e3")
# numerically: the 3-fold content of the temperature boundary is nonzero at every beta
def modes(rr):
    rr = rr / rr.mean()
    c = [np.mean(rr * np.exp(-1j * n * th)) for n in (3, 6)]
    return c


ph3 = modes(tri.radial_polytope(U2))[0]
ph3 = ph3 / abs(ph3)


def fourier36(rr):
    c3, c6 = modes(rr)
    return float((c3 / ph3).real), float((c6 / ph3 ** 2).real)


c3 = []
for beta in [0.01, 0.1, 1, 10, 100]:
    c3.append(abs(modes(lse_radial(tri, U2, beta))[0]))
check("S3", "Measured: the temperature boundary's 3-fold content is nonzero at every β, hot or cold", min(c3) > 1e-8, ", ".join(f"β={b}: {c:.2e}" for b, c in zip([0.01, 0.1, 1, 10, 100], c3)))


# S4: compare the two families in the plane of Fourier content (c3, c6): if they were one family reparametrised,
# the two curves would coincide.
fp = np.array([fourier36(tri.radial(U2, p)) for p in np.geomspace(2.0001, 400, 60)])
te = np.array([fourier36(lse_radial(tri, U2, b)) for b in np.geomspace(0.02, 400, 80)])
gaps = []
o = np.argsort(te[:, 0])
for c3v, c6v in fp:
    if te[:, 0].min() < c3v < te[:, 0].max():
        gaps.append(abs(c6v - np.interp(c3v, te[o, 0], te[o, 1])))
check("S4", "At the triangle's own offset the two routes are different families: at equal 3-fold content their 6-fold "
            "content differs, so no change of speed turns one into the other",
      max(gaps) > 1e-3, f"largest 6-fold gap {max(gaps):.2e} at matched 3-fold content, over {len(gaps)} matched pairs")
np.save(os.path.join(HERE, "depth_simplex_fourier.npy"), np.array([fp, te], dtype=object), allow_pickle=True)
# S5: exact limit a -> infinity with p = beta a
errs = []
beta = 3.0
r_lse = lse_radial(tri, U2, beta)
for a in [1, 4, 16, 64, 256]:
    fd = FaceDial(regular_polygon(3), "tri", a=a)
    errs.append(np.abs(fd.radial(U2, beta * a) - r_lse).max())
check("S5", "Exact in one limit only: face-power with offset a and p = β·a tends to the temperature shape as a → ∞, "
            "error falling like 1/a ((1 + y/a)^(βa) → e^(βy))", errs[-1] < 0.01 and errs[0] / errs[-1] > 100,
      ", ".join(f"a={a}: {e:.1e}" for a, e in zip([1, 4, 16, 64, 256], errs)))

# ------------------------------------------------------------------ write
n_pass = sum(v["ok"] for v in RES.values())
OUT.append(f"\n{n_pass} of {len(RES)} pass, {time.time() - T0:.0f} s")
print(OUT[-1])
open(os.path.join(HERE, "depth_checks_output.txt"), "w").write("\n".join(OUT) + "\n")
json.dump(RES, open(os.path.join(HERE, "depth_results.json"), "w"), indent=1, default=str)
