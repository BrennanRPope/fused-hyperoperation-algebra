"""DRAFT round 13 checks: let the law grow the corners too, then glue and read the curvature back with Regge.

Round 12 let the twisted-boost law (../twisted-boosts/) size the glued pieces, but the starting corners (the
600-cell's 120 points) were given. Here the corners are grown from rest by the law's own kicks: start at rest, and
from every point reached apply the same 12 kicks (and their reverses), keeping every new point, until nothing new
appears. Nothing else is put in: 12 kick directions (an icosahedron's) and one kick size.
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

# ------------------------------------------------------------------ the law in 3D (round 12, Z1-Z4)
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

def to_s3(x):        # inverse stereographic: chart -> unit 3-sphere in R^4
    r2 = (x * x).sum(-1, keepdims=True); return np.concatenate([(1 - r2), 2 * x], -1) / (1 + r2)

# ------------------------------------------------------------------ unit quaternions (the 3-sphere), for comparison only
def qm(p, r):
    a1, b1, c1, d1 = np.moveaxis(p, -1, 0); a2, b2, c2, d2 = np.moveaxis(r, -1, 0)
    return np.stack([a1 * a2 - b1 * b2 - c1 * c2 - d1 * d2, a1 * b2 + b1 * a2 + c1 * d2 - d1 * c2,
                     a1 * c2 - b1 * d2 + c1 * a2 + d1 * b2, a1 * d2 + b1 * c2 - c1 * b2 + d1 * a2], -1)
def qsqrt(p):          # the unit quaternion halfway from 1 to p
    s = p + np.array([1.0, 0, 0, 0]); return s / np.linalg.norm(s, axis=-1, keepdims=True)
def to_chart(X):       # unit quaternion -> the law's velocity (stereographic chart, inverse of to_s3)
    return X[..., 1:] / (1 + X[..., :1])

say("== W. Growing the corners from rest with the law's kicks")
rng = np.random.default_rng(13)
a, x = rng.normal(size=(500, 3)) * 0.4, rng.normal(size=(500, 3)) * 0.4
g = qsqrt(to_s3(a))
errW1 = float(np.abs(to_s3(oplus(a, x, 1.0)) - qm(qm(g, to_s3(x)), g)).max())
check("W1 at q = +1 a kick by a acts on the 3-sphere as x -> g x g with g halfway to a: half the kick on each side (500 random pairs)",
      errW1 < 1e-12, f"worst {errW1:.1e}")

t = (1 + 5**0.5) / 2
ico = np.array([(-1, t, 0), (1, t, 0), (-1, -t, 0), (1, -t, 0), (0, -1, t), (0, 1, t), (0, -1, -t), (0, 1, -t),
                (t, 0, -1), (t, 0, 1), (-t, 0, -1), (-t, 0, 1)], float)
ico /= np.linalg.norm(ico, axis=1)[:, None]
def kicks(deg):
    """12 law velocities pointing to an icosahedron's corners, each a law distance of deg degrees from rest"""
    return ico * math.tan(math.radians(deg) / 2)
def grow(K, generations=8, cap=2000):
    """start at rest; from every point reached, apply each kick and its reverse with the law (k (+) p); keep new points.
    The point opposite rest on the 3-sphere lies at infinity in the law's chart: when a kick lands there it is
    recorded as reached, but no kicks are applied from it (its neighbours are reached from other points anyway)."""
    K = np.concatenate([K, -K])
    S = [to_s3(np.zeros(3))]; P = [np.zeros(3)]; front = [np.zeros(3)]; counts = [1]; far = False
    with np.errstate(all="ignore"):
        for _ in range(generations):
            new = []
            for p in front:
                Dn = 1 - 2 * (K @ p) + (K * K).sum(1) * (p @ p)       # the law's denominator (Z1 of round 12)
                for y, dn in zip(oplus(K, np.tile(p, (len(K), 1)), 1.0), Dn):
                    if dn < 1e-9:                                        # the kick lands on the far point
                        far = True; continue
                    sy = to_s3(y)
                    # compare on the sphere, within 1e-7
                    if np.min(np.linalg.norm(np.array(S) - sy, axis=1)) > 1e-7:
                        S.append(sy); P.append(y); new.append(y)
            counts.append(len(P) + far); front = new
            if not new or len(P) > cap:
                break
    return counts, np.array(P), far
def grow_sphere(K, generations=8, cap=20000):
    """the same growth with the kick in its sphere form x -> g x g (W1), where the far point is ordinary"""
    G = qsqrt(to_s3(np.concatenate([K, -K])))
    one = np.array([1.0, 0, 0, 0])
    pts = {tuple(np.round(one, 8)): one}; front = [one]; counts = [1]
    for _ in range(generations):
        new = []
        for x in front:
            for y in qm(qm(G, np.tile(x, (len(G), 1))), G):
                key = tuple(np.round(y, 8))
                if key not in pts:
                    pts[key] = y; new.append(y)
        counts.append(len(pts)); front = new
        if not new or len(pts) > cap:
            break
    return counts, np.array(list(pts.values()))

c72, P, far72 = grow(kicks(72))
c72s, Xs = grow_sphere(kicks(72))
Xc = to_s3(P)
same = all(np.min(np.linalg.norm(Xs - x, axis=1)) < 1e-9 for x in Xc) and len(Xc) + far72 == len(Xs)
# the 600-cell as the binary icosahedral group, in the frame the kicks pick: every grown point is a product of the
# kicks' halves, and the 120 corners sit 36 degrees from their nearest neighbours
ang = np.degrees(np.arccos(np.clip(Xs @ Xs.T, -1, 1))); np.fill_diagonal(ang, 999)
nn = np.sort(ang, axis=1)[:, :12]
check("W2 grow from rest with 12 kicks of 72 degrees toward an icosahedron's corners: it closes after 3 rounds at exactly 120 corners, each with 12 nearest neighbours 36 degrees away: the 600-cell's corners",
      c72[-1] == 120 and c72[-2] == 120 and same and np.abs(nn - 36).max() < 1e-6,
      f"corners after each round {c72}; same 120 in the kick's sphere form; nearest-neighbour angle {nn.min():.6f} to {nn.max():.6f} deg")
results["growth_72"] = c72
scan = []
for deg in range(6, 181, 2):
    cs, _ = grow_sphere(kicks(deg), generations=7, cap=3000)
    scan.append((deg, cs[-1]))
closing = [(d, n) for d, n in scan if n < 3000]
_, X144 = grow_sphere(kicks(144))
same144 = len(X144) == len(Xs) and all(np.min(np.linalg.norm(Xs - x, axis=1)) < 1e-9 for x in X144)
check("W3 of all kick sizes from 6 to 180 degrees (2-degree steps) only 72 and 144 degrees close up, both onto those 120 corners (180 degrees just swaps rest and its opposite)",
      closing == [(72, 120), (144, 120), (180, 2)] and same144,
      f"closing sizes {closing}; 144 degrees gives the same 120 corners: {same144}; every other size passed 3000 points")
results["scan"] = scan
c36, _ = grow_sphere(kicks(36), generations=5, cap=10**6)
check("W4 kicks the size of the nearest-neighbour step (36 degrees) never close: the count grows about 11 times per round",
      c36[-1] > 10**5 and all(c36[i + 1] > 9 * c36[i] for i in range(1, len(c36) - 1)), f"corners after each round {c36}")
results["growth_36"] = c36
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


say("\n== Y. Glue the grown corners and read the curvature back")
# move everything by one small law kick, so the far point comes into the chart: the law's limit c (+) v as |v| -> oo
c = np.array([0.11, -0.07, 0.05])
C = oplus(np.tile(c, (len(P), 1)), P, 1.0)
if far72:
    C = np.vstack([C, -c / (c @ c)])
n = len(C)
with np.errstate(all="ignore"):
    Dl = np.array([law_dist(np.tile(C[i], (n, 1)), C, 1.0) for i in range(n)])
Dl[~np.isfinite(Dl)] = math.pi          # opposite points: the law's gap is infinite (0/0 in the closed form), 2 arctan(oo) = pi
np.fill_diagonal(Dl, np.inf)
step = Dl.min()
Adj = np.abs(Dl - step) < 1e-9
tets = set()
for i in range(n):
    nb = np.where(Adj[i])[0]
    for j, k, l in itertools.combinations(nb, 3):
        if Adj[j, k] and Adj[j, l] and Adj[k, l]:
            tets.add(tuple(sorted((i, j, k, l))))
T0 = np.array(sorted(tets))
ue, Le, de, vol = regge_from_lengths(T0, lambda i, j: law_dist(C[i], C[j], 1.0))
per_edge = np.bincount(np.unique(np.concatenate([np.sort(T0[:, list(p)], axis=1) for p in PAIRS]), axis=0, return_inverse=True)[1].ravel())
check("W5 glue the grown corners wherever the law's distance is smallest (36 degrees): 720 edges and 600 flat tetrahedra, 5 round every edge, every deficit 2 pi - 5 x the dihedral angle",
      Adj.sum() // 2 == 720 and len(T0) == 600 and (per_edge == 5).all() and np.ptp(de) < 1e-12,
      f"step {math.degrees(step):.6f} deg; {Adj.sum() // 2} edges, {len(T0)} tetrahedra; deficit {de.mean():.10f} rad on every edge")
rows = []
Cc, T = C, T0
for level in range(4):
    Ca = np.asarray(Cc)
    ue, Le, de, vol = regge_from_lengths(T, lambda i, j: law_dist(Ca[i], Ca[j], 1.0))
    rows.append((level, len(T), float((Le * de).sum()), vol))
    if level < 3:
        Cc, T = refine_with_law(Ca, T, 1.0)
tgt = 6 * math.pi**2
errs = [abs(r[2] / tgt - 1) for r in rows]
lam = [r[2] / (3 * r[3]) for r in rows]
check("W6 refine with the law's half kicks and read back with Regge: the action tends to Einstein's 6 pi^2 at second order, and with a cosmological constant Lambda rho^2 -> 1",
      errs[-1] < 1e-3 and all(3.3 < errs[i] / errs[i + 1] < 4.6 for i in range(3)) and abs(lam[-1] - 1) < 3e-3,
      "; ".join(f"level {l} ({nt} tets): {S:.5f}, off {e * 100:.3f}%, Lambda rho^2 {x:.5f}" for (l, nt, S, _), e, x in zip(rows, errs, lam)))
results["regge"] = rows
results["lambda_rho2"] = lam
results["summary"] = {"passed": passed, "failed": failed}
say(f"\n{passed} passed, {failed} failed")
json.dump(results, open("results.json", "w"), indent=1, default=float)
open("checks_output.txt", "w").write(_log.getvalue())
