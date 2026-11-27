"""DRAFT round 11 checks: Einstein's equations on a spacetime glued from flat pieces (Regge's method).

Round 8 got Einstein's equations by varying the loop-read curvature of a smooth metric. Regge (1961) does it with no
smooth metric at all: space is glued from flat simplices, and all curvature sits on the hinges (vertices in 2D,
edges in 3D) as a deficit angle: what a loop round the hinge comes back turned by. The action is
S = sum over hinges (hinge size) x (deficit); varying the edge lengths gives Regge's equations.

2D part (X): a sphere glued from flat triangles, compared with the twisted-boost law's twist (../twisted-boosts/),
which round 9 showed equals the area of each velocity triangle at q = +1.
3D part (Y): a 3-sphere glued from flat tetrahedra (the 600-cell, refined 3 times), and a flat 3-torus.
Run: python3 checks.py  (numpy).  Writes results.json and checks_output.txt.
"""
import json, io, math, itertools
import numpy as np

results, passed, failed = {}, 0, 0
_log = io.StringIO()
def say(s):
    print(s); _log.write(s + "\n")
def check(name, ok, detail=""):
    global passed, failed
    passed, failed = (passed + 1, failed) if ok else (passed, failed + 1)
    say(("PASS " if ok else "FAIL ") + name + (("  | " + detail) if detail else ""))

# ------------------------------------------------------------------ the law (from ../twisted-boosts/, as in round 9)
def oplus(u, v, q):
    return (u + v) / (1 - q * np.conj(u) * v)
def tri_twist(a, b, c, q):
    u = oplus(-a, b, q); w = oplus(-a, c, q); v = oplus(-u, w, q)
    return np.angle((1 - q * u * np.conj(v)) / (1 - q * np.conj(u) * v))
def chart_north(p):
    return complex(p[0], p[1]) / (1 + p[2])
def chart_south(p):
    return complex(p[0], -p[1]) / (1 - p[2])

# ------------------------------------------------------------------ 2D: a sphere glued from flat triangles
def icosphere(level):
    t = (1 + 5**0.5) / 2
    V = [(-1, t, 0), (1, t, 0), (-1, -t, 0), (1, -t, 0), (0, -1, t), (0, 1, t), (0, -1, -t), (0, 1, -t),
         (t, 0, -1), (t, 0, 1), (-t, 0, -1), (-t, 0, 1)]
    V = [np.array(v, float) / np.linalg.norm(v) for v in V]
    F = [(0, 11, 5), (0, 5, 1), (0, 1, 7), (0, 7, 10), (0, 10, 11), (1, 5, 9), (5, 11, 4), (11, 10, 2), (10, 7, 6),
         (7, 1, 8), (3, 9, 4), (3, 4, 2), (3, 2, 6), (3, 6, 8), (3, 8, 9), (4, 9, 5), (2, 4, 11), (6, 2, 10), (8, 6, 7), (9, 8, 1)]
    for _ in range(level):
        cache = {}; F2 = []
        def mid(i, j):
            key = (min(i, j), max(i, j))
            if key not in cache:
                m = V[i] + V[j]; V.append(m / np.linalg.norm(m)); cache[key] = len(V) - 1
            return cache[key]
        for a, b, c in F:
            ab, bc, ca = mid(a, b), mid(b, c), mid(c, a)
            F2 += [(a, ab, ca), (b, bc, ab), (c, ca, bc), (ab, bc, ca)]
        F = F2
    return np.array(V), F

def corner_angle(P, A, B):
    a, b = A - P, B - P
    return math.acos(np.clip(a @ b / np.linalg.norm(a) / np.linalg.norm(b), -1, 1))

def rotate_about(v, axis, ang):
    k = axis / np.linalg.norm(axis)
    return v * math.cos(ang) + np.cross(k, v) * math.sin(ang) + k * (k @ v) * (1 - math.cos(ang))

def holonomy_round_vertex(V, faces_at, v):
    """carry a tangent vector across the flat triangles round vertex v (unfolding over each shared edge);
    return the angle it comes back turned by"""
    def rot(f):                      # the face written as (v, x, y), keeping its orientation
        i = f.index(v); return (f[i], f[(i + 1) % 3], f[(i + 2) % 3])
    ring = [rot(f) for f in faces_at[v]]
    # order the faces cyclically: the next face starts with the edge (v, y) that the current one ends with
    order = [ring[0]]
    while len(order) < len(ring):
        nxt = [f for f in ring if f[1] == order[-1][2]]
        if len(nxt) != 1:
            raise RuntimeError("ring not ordered")
        order.append(nxt[0])
    def normal(f):
        n = np.cross(V[f[1]] - V[f[0]], V[f[2]] - V[f[0]]); return n / np.linalg.norm(n)
    f0 = order[0]; e0 = V[f0[1]] - V[v]
    vec = e0 / np.linalg.norm(e0); start = vec.copy()
    for i in range(len(order)):
        f, g = order[i], order[(i + 1) % len(order)]
        shared = f[2]
        axis = V[shared] - V[v]; n1, n2 = normal(f), normal(g)
        ang = math.atan2(np.cross(n1, n2) @ axis / np.linalg.norm(axis), n1 @ n2)   # fold f's plane onto g's about the edge
        vec = rotate_about(vec, axis, ang)
    n0 = normal(f0)
    return math.atan2(np.cross(start, vec) @ n0, start @ vec)

def orient(V, F):
    out = []
    for a, b, c in F:
        if np.cross(V[b] - V[a], V[c] - V[a]) @ (V[a] + V[b] + V[c]) < 0:
            b, c = c, b
        out.append((a, b, c))
    return out

say("== X. A sphere glued from flat triangles: curvature sits at the corners, and the law's twist carries the same total")
rng = np.random.default_rng(11)
Rm = np.linalg.qr(rng.normal(size=(3, 3)))[0]
rowsX = []
hol_err = None
centers = rng.normal(size=(20, 3)); centers /= np.linalg.norm(centers, axis=1)[:, None]   # 20 caps of radius 0.5 rad
for level in (1, 2, 3, 4, 5):
    V, F = icosphere(level); V = V @ Rm.T; F = orient(V, F)
    nv = len(V)
    angle_sum = np.zeros(nv); twist_share = np.zeros(nv); faces_at = [[] for _ in range(nv)]; valence = np.zeros(nv, int)
    total_twist = 0.0
    for f in F:
        a, b, c = f
        for p, q1, q2 in ((a, b, c), (b, c, a), (c, a, b)):
            angle_sum[p] += corner_angle(V[p], V[q1], V[q2])
            faces_at[p].append(f); valence[p] += 1
        zc = (V[a] + V[b] + V[c])[2]
        ch = chart_north if zc > 0 else chart_south
        w = abs(tri_twist(ch(V[a]), ch(V[b]), ch(V[c]), 1.0))
        total_twist += w
        for p in f:
            twist_share[p] += w / 3
    deficit = 2 * math.pi - angle_sum
    rel = np.abs(deficit / twist_share - 1)
    caps = [abs(deficit[V @ c > math.cos(0.5)].sum() / twist_share[V @ c > math.cos(0.5)].sum() - 1) for c in centers]
    rowsX.append((level, nv, len(F), float(deficit.sum()), float(total_twist), float(max(caps)),
                  float(rel[valence == 5].max()), float(rel[valence == 6].max())))
    if level == 3:
        assert (valence == 5).sum() == 12
        hol = np.array([holonomy_round_vertex(V, faces_at, v) for v in range(nv)])
        hol_err = float(np.abs(hol - deficit).max()); hol_n = nv
        hol_range = (float(deficit.min()), float(deficit.max()))
check("X1 carry an arrow round each corner across the flat triangles: it comes back turned by exactly that corner's deficit angle (all 642 corners, level 3)",
      hol_err < 1e-12, f"worst {hol_err:.1e}; deficits {hol_range[0]:.4f} to {hol_range[1]:.4f} rad")
check("X2 total curvature: deficits add to 4 pi (Descartes), and the law's twist over the same triangles adds to 4 pi, at every refinement",
      all(abs(r[3] - 4 * math.pi) < 1e-9 and abs(r[4] - 4 * math.pi) < 1e-9 for r in rowsX),
      "; ".join(f"level {r[0]}: {r[3] / math.pi:.9f} pi and {r[4] / math.pi:.9f} pi" for r in rowsX))
capsX = [r[5] for r in rowsX]
check("X3 over any patch (20 caps of radius 0.5 rad): the deficits inside match the law's twist over the same patch, with the mismatch falling 4 times per halving",
      capsX[-1] < 2e-4 and all(3.3 < capsX[i] / capsX[i + 1] < 4.5 for i in range(len(capsX) - 1)),
      "; ".join(f"level {l} ({nv} corners): worst cap {c:.1e}" for (l, nv, *_), c in zip(rowsX, capsX)))
check("X4 but not corner by corner: at the 12 five-fold corners the deficit stays about 14.6% off the law's local twist however fine the pieces, and six-fold corners level off near 0.3%",
      all(0.14 < r[6] < 0.15 for r in rowsX[2:]) and abs(rowsX[-1][6] - rowsX[-2][6]) < 2e-3 and 0.002 < rowsX[-1][7] < 0.004,
      "; ".join(f"level {l}: five-fold {r5 * 100:.2f}%, six-fold worst {r6 * 100:.2f}%" for l, _, _, _, _, _, r5, r6 in rowsX))
results["sphere2d"] = rowsX

# ------------------------------------------------------------------ 3D: a 3-sphere glued from flat tetrahedra
def cell600():
    phi = (1 + 5**0.5) / 2
    V = []
    for i in range(4):
        for s in (1, -1):
            v = [0.0] * 4; v[i] = s; V.append(v)
    V += [list(s) for s in itertools.product((0.5, -0.5), repeat=4)]
    base = [phi / 2, 0.5, 1 / (2 * phi), 0.0]
    for p in itertools.permutations(range(4)):
        if sum(1 for i in range(4) for j in range(i + 1, 4) if p[i] > p[j]) % 2:
            continue
        for s in itertools.product((1, -1), repeat=3):
            vals = [base[0] * s[0], base[1] * s[1], base[2] * s[2], 0.0]; v = [0.0] * 4
            for k in range(4):
                v[p[k]] = vals[k]
            V.append(v)
    V = np.unique(np.round(np.array(V), 12), axis=0)
    D = np.linalg.norm(V[:, None] - V[None], axis=2); A = np.abs(D - 1 / phi) < 1e-9
    tets = set()
    for i in range(len(V)):
        nb = np.where(A[i])[0]
        for j, k, l in itertools.combinations(nb, 3):
            if A[j, k] and A[j, l] and A[k, l]:
                tets.add(tuple(sorted((i, j, k, l))))
    return V, np.array(sorted(tets))

def subdivide(V, T):
    """each tetrahedron -> 8 (4 corners + an octahedron cut on its shortest diagonal); new points pushed onto the sphere"""
    V = [np.array(v) for v in V]; mid = {}
    def m(a, b):
        key = (min(a, b), max(a, b))
        if key not in mid:
            p = V[a] + V[b]; V.append(p / np.linalg.norm(p)); mid[key] = len(V) - 1
        return mid[key]
    out = []
    for a, b, c, d in T:
        ab, ac, ad, bc, bd, cd = m(a, b), m(a, c), m(a, d), m(b, c), m(b, d), m(c, d)
        out += [(a, ab, ac, ad), (b, ab, bc, bd), (c, ac, bc, cd), (d, ad, bd, cd)]
        diags = [(ab, cd), (ac, bd), (ad, bc)]
        i = int(np.argmin([np.linalg.norm(V[x] - V[y]) for x, y in diags]))
        p, pp = diags[i]; (x, xp), (y, yp) = [diags[j] for j in range(3) if j != i]
        out += [(p, pp, x, y), (p, pp, y, xp), (p, pp, xp, yp), (p, pp, yp, x)]
    return np.array(V), np.array(out)

PAIRS = [(0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3)]
def regge(V, T):
    """edges, their lengths, their deficit angles (2 pi minus the dihedral angles round them), and total volume.
    Each tetrahedron is the flat simplex on its 4 points; works for points in R^3 or R^4."""
    P = V[T]
    E, TH, L = [], [], []
    for i, j in PAIRS:
        k, l = [x for x in range(4) if x not in (i, j)]
        e = P[:, j] - P[:, i]; en = e / np.linalg.norm(e, axis=1)[:, None]
        a = P[:, k] - P[:, i]; b = P[:, l] - P[:, i]
        a = a - (a * en).sum(1)[:, None] * en; b = b - (b * en).sum(1)[:, None] * en
        th = np.arccos(np.clip((a * b).sum(1) / np.linalg.norm(a, axis=1) / np.linalg.norm(b, axis=1), -1, 1))
        E.append(np.sort(T[:, [i, j]], axis=1)); TH.append(th); L.append(np.linalg.norm(e, axis=1))
    E = np.concatenate(E); TH = np.concatenate(TH); L = np.concatenate(L)
    ue, inv = np.unique(E, axis=0, return_inverse=True)
    inv = inv.ravel()
    Le = np.zeros(len(ue)); Le[inv] = L
    delta = 2 * np.pi - np.bincount(inv, TH, minlength=len(ue))
    Q = P[:, 1:] - P[:, :1]
    vol = np.sqrt(np.abs(np.linalg.det(np.einsum('nij,nkj->nik', Q, Q)))) / 6
    return ue, Le, delta, float(vol.sum())

say("\n== Y. A 3-sphere glued from flat tetrahedra: Regge's action and equations")
V0, T0 = cell600()
ue, Le, de, vol = regge(V0, T0)
exact = 2 * math.pi - 5 * math.acos(1 / 3)
check("Y1 the 600-cell (120 corners, 720 edges, 600 flat tetrahedra): every edge has 5 tetrahedra round it, so every deficit is 2 pi - 5 arccos(1/3)",
      len(V0) == 120 and len(ue) == 720 and len(T0) == 600 and np.abs(de - exact).max() < 1e-12,
      f"deficit {exact:.10f} rad = {math.degrees(exact):.4f} deg on all 720 edges")

# Y2: Schlafli: sum_e L_e d(delta_e) = 0 for any change of shape, so dS/dL_e = delta_e (Regge's equations)
def action_and_grad_check(V, T, seed):
    r = np.random.default_rng(seed)
    dV = r.normal(size=V.shape) * 0.02
    h = 1e-6
    _, Lp, dp, _ = regge(V + h * dV, T); _, Lm, dm, _ = regge(V - h * dV, T); _, L0, d0, _ = regge(V, T)
    dS = ((Lp * dp).sum() - (Lm * dm).sum()) / (2 * h)
    regge_pred = (d0 * (Lp - Lm) / (2 * h)).sum()
    schlafli = (L0 * (dp - dm) / (2 * h)).sum()
    scale = np.abs(d0 * (Lp - Lm) / (2 * h)).sum()
    return dS, regge_pred, schlafli, scale
V1, T1 = subdivide(V0, T0)
sch = [action_and_grad_check(V0, T0, 1), action_and_grad_check(V1, T1, 2)]
check("Y2 Schlafli's identity: for random reshapings, the edges' lengths times the change in their deficits add to zero, so varying the action gives just the deficits (Regge's equations)",
      all(abs(s) < 1e-6 * sc and abs(dS - p) < 1e-6 * sc for dS, p, s, sc in sch),
      "; ".join(f"level {i}: sum L d(delta) = {s:.1e} (scale {sc:.1f}), dS = {dS:.6f} vs sum delta dL = {p:.6f}" for i, (dS, p, s, sc) in enumerate(sch)))

rowsY = []
V, T = V0, T0
for level in range(4):
    ue, Le, de, vol = regge(V, T)
    S = float((Le * de).sum())
    rowsY.append((level, len(T), len(ue), S, vol, float(de.min()), float(de.max())))
    if level < 3:
        V, T = subdivide(V, T)
target = 6 * math.pi**2          # (1/2) x integral of R over the unit 3-sphere: R = 6, volume 2 pi^2
errs = [abs(S / target - 1) for _, _, _, S, _, _, _ in rowsY]
ratios = [errs[i] / errs[i + 1] for i in range(3)]
check("Y3 Regge's action, sum over edges of length x deficit, tends to (1/2) x the integral of R over the 3-sphere = 6 pi^2, with the error falling 4 times per halving (second order)",
      errs[-1] < 1e-3 and all(3.5 < r < 4.5 for r in ratios),
      "; ".join(f"level {l} ({n} tets): {S:.5f} vs {target:.5f}, off {e * 100:.3f}%" for (l, n, _, S, _, _, _), e in zip(rowsY, errs))
      + f"; error ratios {', '.join(f'{r:.2f}' for r in ratios)}")
# Y4: with a cosmological constant. Continuum: S = int (R - 2 Lambda) sqrt g on a 3-sphere of radius rho is
# 12 pi^2 rho - 4 pi^2 Lambda rho^3, stationary at Lambda rho^2 = 1 (the 3D Einstein equation G = -Lambda g).
# Regge: S = 2 rho sum(L delta) - 2 Lambda rho^3 sum(vol) (unit-size lengths and volumes), stationary at
# Lambda rho^2 = sum(L delta) / (3 sum(vol)).
lam = [S / (3 * vol) for _, _, _, S, vol, _, _ in rowsY]
check("Y4 with a cosmological constant: the glued 3-sphere's best size gives Lambda rho^2 -> 1, the answer Einstein's equations give for a smooth 3-sphere",
      abs(lam[-1] - 1) < 3e-3 and all(3.5 < abs(lam[i] - 1) / abs(lam[i + 1] - 1) < 4.5 for i in range(3)),
      "; ".join(f"level {l}: {x:.5f}" for (l, *_), x in zip(rowsY, lam)) + " (off by 4 times less per halving)")

# Y5: flat space solves Regge's equations: a 3-torus glued from cubes, each cut into 6 tetrahedra (Kuhn)
def kuhn_torus(n):
    pts = np.array([(i, j, k) for i in range(n) for j in range(n) for k in range(n)], float)
    idx = lambda i, j, k: ((i % n) * n + (j % n)) * n + (k % n)
    T, shifts = [], []
    for i in range(n):
        for j in range(n):
            for k in range(n):
                for perm in itertools.permutations(range(3)):
                    c = [0, 0, 0]; verts = [(i, j, k)]
                    for ax in perm:
                        c[ax] += 1; verts.append((i + c[0], j + c[1], k + c[2]))
                    T.append([idx(*v) for v in verts]); shifts.append(verts)
    return pts, np.array(T), np.array(shifts, float)
def regge_torus(n):
    pts, T, P = kuhn_torus(n)
    E, TH = [], []
    for i, j in PAIRS:
        k, l = [x for x in range(4) if x not in (i, j)]
        e = P[:, j] - P[:, i]; en = e / np.linalg.norm(e, axis=1)[:, None]
        a = P[:, k] - P[:, i]; b = P[:, l] - P[:, i]
        a = a - (a * en).sum(1)[:, None] * en; b = b - (b * en).sum(1)[:, None] * en
        TH.append(np.arccos(np.clip((a * b).sum(1) / np.linalg.norm(a, axis=1) / np.linalg.norm(b, axis=1), -1, 1)))
        E.append(np.sort(T[:, [i, j]], axis=1))
    E = np.concatenate(E); TH = np.concatenate(TH)
    ue, inv = np.unique(E, axis=0, return_inverse=True)
    return len(T), len(ue), 2 * np.pi - np.bincount(inv.ravel(), TH)
nT, nE, dflat = regge_torus(4)
check("Y5 flat space solves Regge's equations: a 3-torus glued from 64 cubes, each cut into 6 flat tetrahedra, has zero deficit on every edge",
      np.abs(dflat).max() < 1e-12, f"{nT} tetrahedra, {nE} edges, worst deficit {np.abs(dflat).max():.1e}")
results["sphere3d"] = rowsY
results["lambda_rho2"] = lam
results["schlafli"] = [list(map(float, s)) for s in sch]
results["summary"] = {"passed": passed, "failed": failed}

say(f"\n{passed} passed, {failed} failed")
json.dump(results, open("results.json", "w"), indent=1)
open("checks_output.txt", "w").write(_log.getvalue())
