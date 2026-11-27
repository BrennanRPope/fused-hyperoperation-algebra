"""DRAFT round 18 checks: are the law's two 3-spheres one?

Round 17 put light on the 3-sphere the law grows as space (its 3D form, rounds 12-13) and found that round 10's
charge-1/2 twist, a field on the 3-sphere of the law's 2D boosts, is a lowest photon of one hand. The two 3-spheres were
matched by the usual R^4 = C^2, not by the law. Here the law's own operations are compared directly: its 3D sum against
the unit-quaternion product, its 2D boosts against points of the space it grows, its grown corners against a group, and
its own transport against the two Hopf fields.
Run: python3 checks.py  (numpy, scipy).  Writes results.json and checks_output.txt.
"""
import json, io, math, itertools, time
import numpy as np
from scipy.optimize import brentq

results, passed, failed = {}, 0, 0
_log = io.StringIO()
def say(s):
    print(s, flush=True); _log.write(s + "\n")
def check(name, ok, detail=""):
    global passed, failed
    passed, failed = (passed + 1, failed) if ok else (passed, failed + 1)
    say(("PASS " if ok else "FAIL ") + name + (("  | " + detail) if detail else ""))
T_START = time.time()

# ------------------------------------------------------------------ the law in 3D and the corners it grows (rounds 12-13, unchanged)
def oplus(a, b, qq):
    ab = (a * b).sum(-1, keepdims=True); aa = (a * a).sum(-1, keepdims=True); bb = (b * b).sum(-1, keepdims=True)
    return ((1 - 2 * qq * ab - qq * bb) * a + (1 + qq * aa) * b) / (1 - 2 * qq * ab + qq**2 * aa * bb)
def law_dist(a, b, qq):
    g = np.linalg.norm(oplus(-a, b, qq), axis=-1)
    return 2 * np.arctan(math.sqrt(qq) * g) / math.sqrt(qq)
def half_kick(g, qq):
    return g / (1 + np.sqrt(1 + qq * (g * g).sum(-1, keepdims=True)))
def law_mid(a, b, qq):
    return oplus(a, half_kick(oplus(-a, b, qq), qq), qq)
def to_s3(x):        # inverse stereographic: the law's chart -> unit 3-sphere in R^4
    r2 = (x * x).sum(-1, keepdims=True); return np.concatenate([(1 - r2), 2 * x], -1) / (1 + r2)

t_ = (1 + 5**0.5) / 2
ico = np.array([(-1, t_, 0), (1, t_, 0), (-1, -t_, 0), (1, -t_, 0), (0, -1, t_), (0, 1, t_), (0, -1, -t_), (0, 1, -t_),
                (t_, 0, -1), (t_, 0, 1), (-t_, 0, -1), (-t_, 0, 1)], float)
ico /= np.linalg.norm(ico, axis=1)[:, None]
def grow(K, generations=8):
    """round 13: start at rest, apply each kick and its reverse with the law to every point reached, keep new points"""
    K = np.concatenate([K, -K])
    S = [to_s3(np.zeros(3))]; P = [np.zeros(3)]; front = [np.zeros(3)]; far = False
    with np.errstate(all="ignore"):
        for _ in range(generations):
            new = []
            for p in front:
                Dn = 1 - 2 * (K @ p) + (K * K).sum(1) * (p @ p)
                for y, dn in zip(oplus(K, np.tile(p, (len(K), 1)), 1.0), Dn):
                    if dn < 1e-9:
                        far = True; continue
                    sy = to_s3(y)
                    if np.min(np.linalg.norm(np.array(S) - sy, axis=1)) > 1e-7:
                        S.append(sy); P.append(y); new.append(y)
            front = new
            if not new:
                break
    return np.array(P), far
PAIRS = [(0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3)]
def refine_with_law(C, T, qq):
    """round 13: cut each tetrahedron into 8 with the law's half kicks, octahedron cut on its law-shortest diagonal"""
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

P, far = grow(ico * math.tan(math.radians(72) / 2))
c = np.array([0.11, -0.07, 0.05])                      # one small kick brings the far point into the chart (round 13)
C = oplus(np.tile(c, (len(P), 1)), P, 1.0)
if far:
    C = np.vstack([C, -c / (c @ c)])
n = len(C)
with np.errstate(all="ignore"):
    Dl = np.array([law_dist(np.tile(C[i], (n, 1)), C, 1.0) for i in range(n)])
Dl[~np.isfinite(Dl)] = math.pi
np.fill_diagonal(Dl, np.inf)
Adj = np.abs(Dl - Dl.min()) < 1e-9
tets = set()
for i in range(n):
    nb = np.where(Adj[i])[0]
    for j, k, l in itertools.combinations(nb, 3):
        if Adj[j, k] and Adj[j, l] and Adj[k, l]:
            tets.add(tuple(sorted((i, j, k, l))))
LEVELS = [(C, np.array(sorted(tets)))]
for _ in range(2):
    LEVELS.append(refine_with_law(*LEVELS[-1], 1.0))
say(f"corners grown by the law (round 13): {n}; space pieces at levels 0-2: {[len(T) for _, T in LEVELS]}  ({time.time() - T_START:.0f} s)")

import scipy.sparse as sp, scipy.sparse.linalg as sla

# ------------------------------------------------------------------ exterior calculus on the glued pieces (Whitney 1957)
def mesh(L):
    """corners, edges, triangles, pieces of level L, the incidence matrices d0, d1, d2 (integers), pieces oriented outward"""
    C, T = LEVELS[L]; X = to_s3(C); T = np.array(T)
    nv = len(X)
    E = np.unique(np.sort(np.concatenate([T[:, list(p)] for p in PAIRS]), 1), axis=0)
    F = np.unique(np.sort(np.concatenate([T[:, list(f)] for f in [(0,1,2),(0,1,3),(0,2,3),(1,2,3)]]), 1), axis=0)
    eid = {tuple(e): i for i, e in enumerate(E)}; fid = {tuple(f): i for i, f in enumerate(F)}
    ne, nf, nt = len(E), len(F), len(T)
    d0 = sp.csr_matrix((np.r_[-np.ones(ne), np.ones(ne)], (np.r_[np.arange(ne), np.arange(ne)], np.r_[E[:, 0], E[:, 1]])), shape=(ne, nv))
    r, c, v = [], [], []
    for k, (a, b, cc) in enumerate(F):
        for (x, y), s in (((b, cc), 1), ((a, cc), -1), ((a, b), 1)):
            r.append(k); c.append(eid[(x, y)]); v.append(s)
    d1 = sp.csr_matrix((v, (r, c)), shape=(nf, ne))
    # orient each tet by the outward normal of the sphere: det[x0, x1, x2, x3] > 0
    P = X[T]; o = np.sign(np.linalg.det(P)); T = T.copy(); sw = o < 0; T[sw, 0], T[sw, 1] = T[sw, 1].copy(), T[sw, 0].copy()
    r, c, v = [], [], []
    for k, tt in enumerate(T):
        for m in range(4):
            face = [tt[j] for j in range(4) if j != m]
            perm = np.argsort(face); sgn = (-1) ** m * (1 if np.linalg.det(np.eye(3)[perm]) > 0 else -1)
            r.append(k); c.append(fid[tuple(sorted(face))]); v.append(sgn)
    d2 = sp.csr_matrix((v, (r, c)), shape=(nt, nf))
    return dict(X=X, T=T, E=E, F=F, d0=d0, d1=d1, d2=d2, eid=eid)

def whitney(m):
    """Whitney 1-forms on each flat piece, in that piece's own 3D frame (oriented with the sphere's outward normal):
    K = int curl W . curl W (curl-curl), M = int W . W (mass), H = int W . curl W (helicity, the handedness form)"""
    X, T, eid = m["X"], m["T"], m["eid"]; nt = len(T)
    P = X[T]; D = P[:, 1:] - P[:, :1]                         # (nt,3,4)
    Q, _ = np.linalg.qr(np.transpose(D, (0, 2, 1)))           # (nt,4,3) orthonormal frame of each piece's 3-plane
    nrm = P.mean(1)                                          # outward: the centroid points away from the sphere's centre
    s = np.sign(np.linalg.det(np.concatenate([Q, nrm[:, :, None]], 2)))
    Q[:, :, 2] *= s[:, None]                                  # local frame + outward normal positively oriented in R^4
    Dl = np.einsum('nik,nkj->nij', D, Q)                      # edge vectors in local 3D coords
    Gi = np.linalg.inv(Dl)                                     # rows of Dl are edges; grads of lambda_1..3 = columns of inv
    g123 = np.transpose(Gi, (0, 2, 1)); g = np.concatenate([-g123.sum(1, keepdims=True), g123], 1)   # (nt,4,3)
    vol = np.abs(np.linalg.det(Dl)) / 6
    gg = np.einsum('nik,njk->nij', g, g)
    Mll = (np.ones((4, 4)) + np.eye(4)) / 20
    Me = np.zeros((nt, 6, 6)); Ke = np.zeros((nt, 6, 6)); He = np.zeros((nt, 6, 6))
    for a, (i, j) in enumerate(PAIRS):
        for b, (k, l) in enumerate(PAIRS):
            Me[:, a, b] = vol * (Mll[i, k] * gg[:, j, l] - Mll[i, l] * gg[:, j, k] - Mll[j, k] * gg[:, i, l] + Mll[j, l] * gg[:, i, k])
            Ke[:, a, b] = vol * 4 * (gg[:, i, k] * gg[:, j, l] - gg[:, i, l] * gg[:, j, k])
            crs = np.cross(g[:, k], g[:, l])                  # curl W_kl = 2 g_k x g_l ; int W_ij = vol/4 (g_j - g_i)
            He[:, a, b] = vol / 4 * 2 * np.einsum('nk,nk->n', g[:, j] - g[:, i], crs)
    gl = np.array([[eid[(min(t[i], t[j]), max(t[i], t[j]))] for (i, j) in PAIRS] for t in T])
    sg = np.array([[1.0 if t[i] < t[j] else -1.0 for (i, j) in PAIRS] for t in T])
    def asm(Ae):
        Ae = Ae * sg[:, :, None] * sg[:, None, :]
        r = np.repeat(gl, 6, 1).ravel(); c = np.tile(gl, (1, 6)).ravel()
        return sp.csr_matrix((Ae.ravel(), (r, c)), shape=(len(eid), len(eid)))
    M0 = np.bincount(T.ravel(), weights=np.repeat(vol / 4, 4), minlength=len(X))
    return asm(Ke), asm(Me), asm(He), M0, vol.sum()

# ================================================================== G. The law's two 3-spheres
say("\n== G. Are the 3-sphere the law grows and the 3-sphere of its boosts one?")
def qmul(a, b):                                   # unit quaternions as (w, x, y, z): the 3-sphere's own product
    a0, av = a[..., :1], a[..., 1:]; b0, bv = b[..., :1], b[..., 1:]
    return np.concatenate([a0 * b0 - (av * bv).sum(-1, keepdims=True), a0 * bv + b0 * av + np.cross(av, bv)], -1)
def qconj(a):
    return a * np.array([1, -1, -1, -1.0])
pt = lambda x: to_s3(x)                           # the law's chart -> the 3-sphere (round 13's placement)
hk = lambda x: half_kick(x, 1.0)                  # the law's half kick
rg = np.random.default_rng(18)
X1 = rg.normal(size=(2000, 3)) * 0.8; Y1 = rg.normal(size=(2000, 3)) * 0.8
S = pt(oplus(X1, Y1, 1.0)); H_ = pt(hk(X1)); PX, PY = pt(X1), pt(Y1)
sand = np.abs(qmul(qmul(H_, PY), H_) - S).max()
grp = np.abs(qmul(PX, PY) - qmul(qmul(H_, S), qconj(H_))).max()
noneq = np.abs(qmul(PX, PY) - S).max()
check("G1 (exact) the law's sum in 3D at q = +1 is the 3-sphere's own product taken symmetrically: x (+) y lands at h . y . h, where h is the law's half kick of x. The plain product x . y differs from it only by a turn (conjugation by h). So the space the law grows is the unit-quaternion group, and the law is its symmetric product",
      sand < 1e-12 and grp < 1e-12 and noneq > 0.5,
      f"2000 random pairs: h.y.h vs the law {sand:.0e}; x.y vs h (x (+) y) h^-1 {grp:.0e}; x.y vs the law directly, for contrast, {noneq:.2f}")
results["sandwich"] = dict(sandwich=float(sand), group_vs_law=float(grp), plain_vs_law=float(noneq))

emb = lambda u: np.stack([0 * np.real(u), np.real(u), np.imag(u)], -1)
op2 = lambda a, b: (a + b) / (1 - np.conj(a) * b)
U = rg.normal(size=2000) * 0.8 + 1j * rg.normal(size=2000) * 0.8; Vv = rg.normal(size=2000) * 0.8 + 1j * rg.normal(size=2000) * 0.8
e_plane = np.abs(oplus(emb(U), emb(Vv), 1.0) - emb(op2(U, Vv))).max()
Bq = lambda w: pt(hk(emb(w)))                    # a 2D boost, placed in the grown space at its velocity's half kick
Xb = Bq(U)
SU2 = np.stack([np.stack([Xb[:, 0] + 1j * Xb[:, 1], Xb[:, 2] + 1j * Xb[:, 3]], -1),
                np.stack([-(Xb[:, 2] - 1j * Xb[:, 3]), Xb[:, 0] - 1j * Xb[:, 1]], -1)], -2)
Bm = np.stack([np.stack([np.ones_like(U), U], -1), np.stack([-np.conj(U), np.ones_like(U)], -1)], -2) / np.sqrt(1 + abs(U)**2)[:, None, None]
e_pt = np.abs(SU2 - Bm).max()
cc = 1 - np.conj(U) * Vv; al = -np.angle(cc)
fibq = np.stack([np.cos(al), np.sin(al), 0 * al, 0 * al], -1)
e_comp = np.abs(qmul(Bq(U), Bq(Vv)) - qmul(Bq(op2(U, Vv)), fibq)).max()
check("G2 (exact) the 2D law is the 3D law in a plane through rest, and each 2D boost (divided by its size) is the point of the grown space at its velocity's half kick. The law's factorization B(u)B(v) = B(u (+) v) . (fibre turn) is then the 3-sphere's own product, with the fibre turn a right turn about the plane's normal. The two 3-spheres are one: the match round 17 assumed is the law's own, with only the hand left as a choice (G5)",
      e_plane < 1e-9 and e_pt < 1e-12 and e_comp < 1e-12,
      f"2000 random pairs: 2D law vs 3D law in the plane {e_plane:.0e}; boost vs placed point {e_pt:.0e}; factorization as the 3-sphere's product {e_comp:.0e}")
results["one_sphere"] = dict(plane=float(e_plane), point=float(e_pt), composition=float(e_comp))

P13, far13 = grow(ico * math.tan(math.radians(72) / 2))
Qc = np.vstack([pt(P13), [[-1.0, 0, 0, 0]]] if far13 else [pt(P13)])
prod = qmul(np.repeat(Qc, len(Qc), 0), np.tile(Qc, (len(Qc), 1)))
dmin = np.min(np.linalg.norm(prod[:, None, :] - Qc[None], axis=2), axis=1)
inv_ok = np.min(np.linalg.norm(qconj(Qc)[:, None, :] - Qc[None], axis=2), axis=1).max()
rest = np.argmin(np.linalg.norm(Qc - [1, 0, 0, 0], axis=1))
# rotations of 3D space by these group elements (q v q^-1): how many distinct?
Rots = np.round(np.array([[qmul(qmul(q[None], np.r_[0, e][None]), qconj(q)[None])[0][1:] for e in np.eye(3)] for q in Qc]), 9)
n_rot = len(np.unique(Rots.reshape(len(Qc), -1), axis=0))
check("G3 (exact) the corners the law grew (round 13) are a group under the 3-sphere's product: all 14,400 products land on corners, rest is the identity, each corner's reverse kick is its inverse. It is the binary icosahedral group: 120 elements, turning 3D space by the 60 rotations of the icosahedron, two elements per rotation",
      len(Qc) == 120 and dmin.max() < 1e-9 and inv_ok < 1e-9 and np.linalg.norm(Qc[rest] - [1, 0, 0, 0]) < 1e-12 and n_rot == 60,
      f"{len(Qc)} corners; largest distance of a product from the set {dmin.max():.0e}; inverses {inv_ok:.0e}; distinct rotations {n_rot}")
results["group"] = dict(corners=len(Qc), closure=float(dmin.max()), rotations=n_rot, points_wxyz=np.round(Qc, 12).tolist())

def hopf_field(X, hand=1):                       # round 17's two fields: hand +1 = 'Hopf', hand -1 = 'mirror'
    return np.stack([-X[..., 1], X[..., 0], -hand * X[..., 3], hand * X[..., 2]], -1)
loops = []
def law_transport(th, N=3000, order="right", mirror=False):
    """carry a frame round a circle of velocities of angular radius th by composing the law's boosts.
    order "right": g -> g . B(step), round 10's B(u)B(v); order "left": g -> B(step) . g, the same boosts composed the other way.
    mirror: also report the path carried through the parity x -> -x, which on the 3-sphere is q -> conj(q)."""
    r0 = math.tan(th / 2); pts = r0 * np.exp(2j * np.pi * np.arange(N + 1) / N)
    place = lambda w: pt(hk(emb(np.atleast_1d(w))))[0]
    g = place(pts[0]); G = [g]; ph = 0.0
    for k in range(N):
        step = op2(-pts[k], pts[k + 1]); b = place(step)
        g = qmul(g[None], b[None])[0] if order == "right" else qmul(b[None], g[None])[0]; G.append(g)
        c_ = 1 - np.conj(pts[k]) * step; ph += float(np.angle(np.conj(c_) / c_)) / 2
    G = np.array(G)
    if mirror:
        G = qconj(G)
    mid = 0.5 * (G[1:] + G[:-1]); dG = G[1:] - G[:-1]
    return dict(cap=th, order=order, mirrored=mirror, half_area=math.pi * (1 - math.cos(th)), law_phase=ph,
                along_hopf=float((hopf_field(mid, 1) * dG).sum()), along_mirror=float((hopf_field(mid, -1) * dG).sum()),
                worst_step_mirror=float(np.abs((hopf_field(mid, -1) * dG).sum(1)).max()))
for th in (0.3, 0.9, 1.5, 2.2):
    loops.append(law_transport(th))
okG4 = all(abs(r["along_mirror"]) < 1e-12 and abs(r["along_hopf"]) > 0.05 and abs(r["law_phase"] - r["half_area"]) < 1e-5 for r in loops)
check("G4 the law's own transport (carry a frame along a loop of velocities by composing the law's boosts) moves along exactly one of round 17's two Hopf fields: it never moves along the one round 17 called the mirror, so that is the law's connection, and the phase it collects is half the enclosed area (round 10). Round 17 attached the law's twist to the other of the pair; the law's own placement picks the mirror field, curl +2 in round 17's orientation",
      okG4, "; ".join(f"cap {r['cap']}: law's phase {r['law_phase']:.5f} = half area {r['half_area']:.5f}; along the law's transport, mirror field {r['along_mirror']:.0e}, Hopf field {r['along_hopf']:.3f}" for r in loops))
results["transport"] = loops

e_par = np.abs(oplus(-X1, -Y1, 1.0) + oplus(X1, Y1, 1.0)).max()
Rm = np.diag([1.0, 1.0, -1.0])
e_ref = np.abs(oplus(X1 @ Rm.T, Y1 @ Rm.T, 1.0) - oplus(X1, Y1, 1.0) @ Rm.T).max()
e_conj = np.abs(pt(-X1) - qconj(pt(X1))).max()
mirrored = [law_transport(th, mirror=True) for th in (0.9, 2.2)]
lefted = [law_transport(th, order="left") for th in (0.9, 2.2)]
okG5 = (e_par < 1e-12 and e_ref < 1e-12 and e_conj < 1e-12
        and all(abs(r["along_hopf"]) < 1e-12 and abs(r["along_mirror"]) > 0.05 for r in mirrored + lefted))
check("G5 the law has no hand of its own. It is unchanged by a reflection of space and by the parity x -> -x, and the parity is the 3-sphere's conjugation q -> q-bar, which reverses products and so swaps the two Hopf fields: the law's transport seen through the parity runs along the Hopf field instead. Composing the same boosts in the other order does the same. The photon's hand is fixed by the orientation together with the order of composition (round 10's B(u)B(v)), not by the law",
      okG5, f"law under parity {e_par:.0e}, under a reflection {e_ref:.0e}; parity on the 3-sphere vs conjugation {e_conj:.0e}; "
      + "through the parity: along the Hopf field " + ", ".join(f"{r['along_hopf']:.0e}" for r in mirrored)
      + ", along the mirror field " + ", ".join(f"{r['along_mirror']:.3f}" for r in mirrored)
      + "; other order: along the Hopf field " + ", ".join(f"{r['along_hopf']:.0e}" for r in lefted)
      + ", along the mirror field " + ", ".join(f"{r['along_mirror']:.3f}" for r in lefted))
results["mirror"] = dict(parity=float(e_par), reflection=float(e_ref), parity_is_conjugation=float(e_conj), through_parity=mirrored, other_order=lefted)

rows6 = []
for L in range(3):
    m = mesh(L); K, M, H, M0, V = whitney(m); Xm, E = m["X"], m["E"]
    Af = hopf_field(Xm, -1); a = (0.5 * (Af[E[:, 0]] + Af[E[:, 1]]) * (Xm[E[:, 1]] - Xm[E[:, 0]])).sum(1)
    nrm2 = a @ (M @ a); rows6.append(dict(level=L, freq2=float(a @ (K @ a) / nrm2), curl=float(a @ (H @ a) / nrm2)))
check("G6 on the law's pieces, the law's own field is a lowest photon: frequency^2 -> 4 and curl -> +2 with the pieces oriented outward (round 17, E5, found its overlap with the lowest level to be 1.0000). So the identification round 17 assumed is now forced, with the hand corrected",
      abs(rows6[2]["curl"] - 2) < 0.02 and abs(rows6[2]["freq2"] - 4) < 0.05 and rows6[0]["curl"] > rows6[1]["curl"] > rows6[2]["curl"] > 2,
      "; ".join(f"level {r['level']}: frequency^2 {r['freq2']:.4f}, curl {r['curl']:+.4f}" for r in rows6))
results["law_field_on_pieces"] = rows6

results["summary"] = {"passed": passed, "failed": failed, "seconds": round(time.time() - T_START)}
say(f"\n{passed} passed, {failed} failed  ({time.time() - T_START:.0f} s)")
json.dump(results, open("results.json", "w"), indent=1, default=float)
open("checks_output.txt", "w").write(_log.getvalue())
