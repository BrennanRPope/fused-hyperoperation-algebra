"""DRAFT round 12 checks: let the law set the glued pieces' sizes, and read the curvature back with Regge's method.

Round 11 glued a 3-sphere from flat tetrahedra whose sizes came from a sphere we picked. Here the sizes come from the
twisted-boost law (../twisted-boosts/) alone. The law u (+) v = (u + v)(1 - q u* v)^-1 is written with quaternions,
so it covers 3D velocities (pure quaternions) as well as 2D ones (complex numbers). Distances are read from the law's
own gap g = (-a) (+) b, midpoints from its "half a kick" g / (1 + sqrt(1 + q |g|^2)), and every tetrahedron's shape
from those lengths alone (no embedding). Regge's deficits then give back curvature K = q.
Run: python3 checks.py  (numpy, sympy).  Writes results.json and checks_output.txt.
"""
import json, io, math, itertools
import numpy as np
import sympy as sp

results, passed, failed = {}, 0, 0
_log = io.StringIO()
def say(s):
    print(s); _log.write(s + "\n")
def check(name, ok, detail=""):
    global passed, failed
    passed, failed = (passed + 1, failed) if ok else (passed, failed + 1)
    say(("PASS " if ok else "FAIL ") + name + (("  | " + detail) if detail else ""))

# ------------------------------------------------------------------ Z1: the law in 3D, exactly
say("== Z. The law with 3D velocities sets its own distances and midpoints")
u = sp.symbols("u1:4", real=True); v = sp.symbols("v1:4", real=True); q = sp.symbols("q", real=True)
U, Vq = sp.Quaternion(0, *u), sp.Quaternion(0, *v)
W = (U + Vq) * (sp.Quaternion(1, 0, 0, 0) - q * U.conjugate() * Vq).inverse()
uv = sum(a * b for a, b in zip(u, v)); uu = sum(a * a for a in u); vv = sum(a * a for a in v)
D = 1 - 2 * q * uv + q**2 * uu * vv
closed = [((1 - 2 * q * uv - q * vv) * u[i] + (1 + q * uu) * v[i]) / D for i in range(3)]
okZ1 = sp.simplify(W.a) == 0 and all(sp.simplify(c - g) == 0 for c, g in zip((W.b, W.c, W.d), closed))
check("Z1 (exact) with 3D velocities (pure quaternions) the law stays 3D: u (+) v = [(1 - 2q u.v - q|v|^2) u + (1 + q|u|^2) v] / (1 - 2q u.v + q^2 |u|^2 |v|^2); at q = -1 this is Ungar's Mobius addition",
      okZ1)

# ------------------------------------------------------------------ numerical law, distance, half-kick
def oplus(a, b, qq):
    """the law on 3-vectors (pure quaternions), in the closed form of Z1; works on arrays of shape (..., 3)"""
    ab = (a * b).sum(-1, keepdims=True); aa = (a * a).sum(-1, keepdims=True); bb = (b * b).sum(-1, keepdims=True)
    return ((1 - 2 * qq * ab - qq * bb) * a + (1 + qq * aa) * b) / (1 - 2 * qq * ab + qq**2 * aa * bb)
def law_dist(a, b, qq):
    g = np.linalg.norm(oplus(-a, b, qq), axis=-1)
    if qq > 0:
        return 2 * np.arctan(math.sqrt(qq) * g) / math.sqrt(qq)
    if qq < 0:
        return 2 * np.arctanh(math.sqrt(-qq) * g) / math.sqrt(-qq)
    return 2 * g
def half_kick(g, qq):
    return g / (1 + np.sqrt(1 + qq * (g * g).sum(-1, keepdims=True)))
def law_mid(a, b, qq):
    return oplus(a, half_kick(oplus(-a, b, qq), qq), qq)

rng = np.random.default_rng(12)
def sample(n, qq):
    x = rng.normal(size=(n, 3)) * 0.35
    if qq < 0:   # keep inside the unit ball, where q = -1 velocities live
        x = x / np.maximum(1, np.linalg.norm(x, axis=1, keepdims=True) / 0.9)
    return x
iso = {}
for qq in (1.0, -1.0):
    a, x, y = sample(500, qq), sample(500, qq), sample(500, qq)
    iso[qq] = float(np.abs(law_dist(oplus(a, x, qq), oplus(a, y, qq), qq) - law_dist(x, y, qq)).max())
check("Z2 the law's own boosts keep the law's own distance d = 2 arctan|(-a)(+)b| (q = +1) or 2 artanh|(-a)(+)b| (q = -1): 500 random triples each",
      max(iso.values()) < 1e-12, f"q = +1: {iso[1.0]:.1e}; q = -1: {iso[-1.0]:.1e}")

# Z3: compare with the textbook geometry, computed without the law
def to_s3(x):        # inverse stereographic: chart -> unit 3-sphere in R^4
    r2 = (x * x).sum(-1, keepdims=True); return np.concatenate([(1 - r2), 2 * x], -1) / (1 + r2)
def to_h3(x):        # Poincare ball -> hyperboloid in R^{1,3}
    r2 = (x * x).sum(-1, keepdims=True); return np.concatenate([(1 + r2), 2 * x], -1) / (1 - r2)
a, b = sample(500, 1.0), sample(500, 1.0)
e_s = float(np.abs(law_dist(a, b, 1.0) - np.arccos(np.clip((to_s3(a) * to_s3(b)).sum(-1), -1, 1))).max())
a, b = sample(500, -1.0), sample(500, -1.0)
A, B = to_h3(a), to_h3(b)
e_h = float(np.abs(law_dist(a, b, -1.0) - np.arccosh(A[:, 0] * B[:, 0] - (A[:, 1:] * B[:, 1:]).sum(-1))).max())
check("Z3 the law's distance is the true distance of the round 3-sphere (q = +1) and of hyperbolic 3-space (q = -1), checked against the embedded spaces without the law",
      e_s < 1e-12 and e_h < 1e-10, f"3-sphere {e_s:.1e}; hyperbolic {e_h:.1e}")
mids = {}
for qq in (1.0, -1.0):
    a, b = sample(500, qq), sample(500, qq)
    m = law_mid(a, b, qq); d = law_dist(a, b, qq)
    mids[qq] = float(max(np.abs(law_dist(a, m, qq) - d / 2).max(), np.abs(law_dist(m, b, qq) - d / 2).max()))
check("Z4 the law's half kick g / (1 + sqrt(1 + q|g|^2)) lands exactly halfway: the midpoint is half the distance from each end",
      max(mids.values()) < 1e-12, f"q = +1: {mids[1.0]:.1e}; q = -1: {mids[-1.0]:.1e}")

# ------------------------------------------------------------------ Regge from lengths alone
PAIRS = [(0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3)]
def tet_shapes(Lsq):
    """Lsq: (n, 4, 4) squared edge lengths -> flat tetrahedra's corner coordinates in R^3 built from the lengths alone"""
    G = (Lsq[:, 0, 1:, None] + Lsq[:, 0, None, 1:] - Lsq[:, 1:, 1:]) / 2
    w, Vec = np.linalg.eigh(G)
    if (w < -1e-14).any():
        raise ValueError("lengths don't make a flat tetrahedron")
    X = Vec * np.sqrt(np.clip(w, 0, None))[:, None, :]             # rows = corners 1..3
    return np.concatenate([np.zeros((len(G), 1, 3)), X], axis=1)
def regge_from_lengths(T, edge_len):
    """T: tetrahedra (corner ids); edge_len(i, j) for id arrays -> edges, lengths, deficits, total volume"""
    n = len(T); Lsq = np.zeros((n, 4, 4))
    for i, j in PAIRS:
        l = edge_len(T[:, i], T[:, j]); Lsq[:, i, j] = Lsq[:, j, i] = l * l
    P = tet_shapes(Lsq)
    E, TH, L = [], [], []
    for i, j in PAIRS:
        k, l_ = [x for x in range(4) if x not in (i, j)]
        e = P[:, j] - P[:, i]; en = e / np.linalg.norm(e, axis=1)[:, None]
        a = P[:, k] - P[:, i]; b = P[:, l_] - P[:, i]
        a = a - (a * en).sum(1)[:, None] * en; b = b - (b * en).sum(1)[:, None] * en
        TH.append(np.arccos(np.clip((a * b).sum(1) / np.linalg.norm(a, axis=1) / np.linalg.norm(b, axis=1), -1, 1)))
        E.append(np.sort(T[:, [i, j]], axis=1)); L.append(np.sqrt(Lsq[:, i, j]))
    E = np.concatenate(E); TH = np.concatenate(TH); L = np.concatenate(L)
    ue, inv = np.unique(E, axis=0, return_inverse=True); inv = inv.ravel()
    Le = np.zeros(len(ue)); Le[inv] = L
    delta = 2 * np.pi - np.bincount(inv, TH, minlength=len(ue))
    vol = np.abs(np.linalg.det(P[:, 1:])) / 6
    return ue, Le, delta, float(vol.sum())

# ------------------------------------------------------------------ the 600-cell's corners, put into the law's chart
def cell600_s3():
    phi = (1 + 5**0.5) / 2
    V = []
    for i in range(4):
        for s in (1, -1):
            x = [0.0] * 4; x[i] = s; V.append(x)
    V += [list(s) for s in itertools.product((0.5, -0.5), repeat=4)]
    base = [phi / 2, 0.5, 1 / (2 * phi), 0.0]
    for p in itertools.permutations(range(4)):
        if sum(1 for i in range(4) for j in range(i + 1, 4) if p[i] > p[j]) % 2:
            continue
        for s in itertools.product((1, -1), repeat=3):
            vals = [base[0] * s[0], base[1] * s[1], base[2] * s[2], 0.0]; x = [0.0] * 4
            for k in range(4):
                x[p[k]] = vals[k]
            V.append(x)
    V = np.unique(np.round(np.array(V), 12), axis=0)
    Dm = np.linalg.norm(V[:, None] - V[None], axis=2); Adj = np.abs(Dm - 1 / phi) < 1e-9
    tets = set()
    for i in range(len(V)):
        nb = np.where(Adj[i])[0]
        for j, k, l in itertools.combinations(nb, 3):
            if Adj[j, k] and Adj[j, l] and Adj[k, l]:
                tets.add(tuple(sorted((i, j, k, l))))
    return V, np.array(sorted(tets))

def refine_with_law(C, T, qq):
    """cut each tetrahedron into 8; new corners are the law's midpoints (half kicks), the octahedron is cut on its
    shortest diagonal measured with the law's distance"""
    C = [c for c in C]; mid = {}
    def m(a, b):
        key = (min(a, b), max(a, b))
        if key not in mid:
            C.append(law_mid(C[key[0]], C[key[1]], qq)); mid[key] = len(C) - 1
        return mid[key]
    out = []
    for a, b, c, d in T:
        ab, ac, ad, bc, bd, cd = m(a, b), m(a, c), m(a, d), m(b, c), m(b, d), m(c, d)
        out += [(a, ab, ac, ad), (b, ab, bc, bd), (c, ac, bc, cd), (d, ad, bd, cd)]
        diags = [(ab, cd), (ac, bd), (ad, bc)]
        i = int(np.argmin([law_dist(C[x], C[y], qq) for x, y in diags]))
        p, pp = diags[i]; (x, xp), (y, yp) = [diags[j] for j in range(3) if j != i]
        out += [(p, pp, x, y), (p, pp, y, xp), (p, pp, xp, yp), (p, pp, yp, x)]
    return np.array(C), np.array(out)

say("\n== Y. A 3-sphere glued from pieces whose sizes come only from the law")
X4, T0 = cell600_s3()
Rm = np.linalg.qr(rng.normal(size=(4, 4)))[0]
X4 = X4 @ Rm.T                                   # random turn, so no corner sits at the chart's far point
chart = X4[:, 1:] / (1 + X4[:, :1])              # stereographic chart: the law's velocities at q = +1
rowsY = {}
for qq in (1.0, 0.25):
    C = chart / math.sqrt(qq)                    # the same corners on the sphere of curvature q (radius 1/sqrt q)
    T = T0; rows = []
    for level in range(4):
        Ca = np.asarray(C)
        ue, Le, de, vol = regge_from_lengths(T, lambda i, j: law_dist(Ca[i], Ca[j], qq))
        rows.append((level, len(T), float((Le * de).sum()), vol, float(de.min()), float(de.max())))
        if level < 3:
            C, T = refine_with_law(Ca, T, qq)
    rowsY[qq] = rows
r1 = rowsY[1.0]; tgt1 = 6 * math.pi**2
e1 = [abs(r[2] / tgt1 - 1) for r in r1]
check("Y1 sizes from the law (q = +1): every edge's length is the law's distance, every new corner a half kick; Regge's action tends to (1/2) x integral of R = 6 pi^2 (R = 6q) at second order",
      e1[-1] < 1e-3 and all(3.3 < e1[i] / e1[i + 1] < 4.6 for i in range(3)),
      "; ".join(f"level {l} ({n} tets): {S:.5f}, off {e * 100:.3f}%" for (l, n, S, *_), e in zip(r1, e1)))
r2 = rowsY[0.25]; tgt2 = 6 * math.pi**2 / math.sqrt(0.25)
e2 = [abs(r[2] / tgt2 - 1) for r in r2]
check("Y2 turning the dial to q = 1/4 with the same corners: the action tends to 6 pi^2 / sqrt(q) = 12 pi^2, a sphere of radius 2, so the law's q is the curvature",
      e2[-1] < 1e-3 and abs(e2[-1] - e1[-1]) < 1e-9,
      "; ".join(f"level {l}: {S:.5f} (vs {tgt2:.5f})" for (l, _, S, *_) in r2))
lam = [r[2] / (3 * r[3]) for r in r1]
check("Y3 with a cosmological constant, the law-built sphere's best size gives Lambda rho^2 -> 1, as Einstein's equations do",
      abs(lam[-1] - 1) < 3e-3 and all(abs(lam[i + 1] - 1) < abs(lam[i] - 1) for i in range(3)),
      "; ".join(f"level {l}: {x:.5f}" for (l, *_), x in zip(r1, lam)))

# ------------------------------------------------------------------ Y4: both ends of the dial, locally
say("\n== X. Both ends of the dial: the deficit's sign follows q")
t = (1 + 5**0.5) / 2
ico = np.array([(-1, t, 0), (1, t, 0), (-1, -t, 0), (1, -t, 0), (0, -1, t), (0, 1, t), (0, -1, -t), (0, 1, -t),
                (t, 0, -1), (t, 0, 1), (-t, 0, -1), (-t, 0, 1)], float)
ico /= np.linalg.norm(ico, axis=1)[:, None]
faces = [(0, 11, 5), (0, 5, 1), (0, 1, 7), (0, 7, 10), (0, 10, 11), (1, 5, 9), (5, 11, 4), (11, 10, 2), (10, 7, 6),
         (7, 1, 8), (3, 9, 4), (3, 4, 2), (3, 2, 6), (3, 6, 8), (3, 8, 9), (4, 9, 5), (2, 4, 11), (6, 2, 10), (8, 6, 7), (9, 8, 1)]
Tstar = np.array([(12, a, b, c) for a, b, c in faces])        # 20 tetrahedra round the centre (corner 12)
def chart_radius(eps, qq):      # a law velocity whose law distance from rest is eps
    if qq > 0:
        return math.tan(math.sqrt(qq) * eps / 2) / math.sqrt(qq)
    if qq < 0:
        return math.tanh(math.sqrt(-qq) * eps / 2) / math.sqrt(-qq)
    return eps / 2
def star_deficit(eps, qq):
    C = np.vstack([ico * chart_radius(eps, qq), np.zeros((1, 3))])
    ue, Le, de, vol = regge_from_lengths(Tstar, lambda i, j: law_dist(C[i], C[j], qq))
    radial = ue[:, 1] == 12
    return float(de[radial].mean()), float(np.ptp(de[radial]))
qs = (1.0, 0.25, 0.0, -0.25, -1.0)
rowsX = {}
for qq in qs:
    rowsX[qq] = [(eps,) + star_deficit(eps, qq) for eps in (0.2, 0.1, 0.05, 0.025)]
flat = max(abs(r[1]) for r in rowsX[0.0])
check("X1 at q = 0 (no curvature) the 20 law-built tetrahedra round a corner close up exactly: zero deficit on all 12 edges from the centre",
      flat < 1e-12, f"worst {flat:.1e}")
cons = {qq: rowsX[qq][-1][1] / (qq * rowsX[qq][-1][0]**2) for qq in qs if qq != 0}
spread = max(cons.values()) - min(cons.values())
check("X2 at every other q the deficit is q x (a fixed number) x size^2: positive (gaps, like a sphere) for q > 0, negative (overlaps, saddle-like) for q < 0, the same number at both ends",
      all(np.sign(rowsX[qq][-1][1]) == np.sign(qq) for qq in cons) and spread < 2e-3 * abs(np.mean(list(cons.values()))),
      "; ".join(f"q = {qq:+.2f}: deficit / (q size^2) = {c:.5f}" for qq, c in cons.items()))
results["sphere_from_law"] = {str(k): v for k, v in rowsY.items()}
results["lambda_rho2"] = lam
results["star"] = {str(k): v for k, v in rowsX.items()}
results["star_constant"] = {str(k): v for k, v in cons.items()}
results["errors"] = {"isometry": iso, "dist_s3": e_s, "dist_h3": e_h, "midpoint": mids}
results["summary"] = {"passed": passed, "failed": failed}
say(f"\n{passed} passed, {failed} failed")
json.dump(results, open("results.json", "w"), indent=1, default=float)
open("checks_output.txt", "w").write(_log.getvalue())
