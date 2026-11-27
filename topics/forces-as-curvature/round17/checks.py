"""DRAFT round 17 checks: electromagnetism as loops on the law's own 3-sphere.

Maxwell's field lives on loops: the potential is a number on each edge (a 1-cochain), the flux through a triangle is
the sum round its edges, and Whitney's forms turn these into a field on the glued flat pieces. Here that is done on the
3-sphere the law grew (round 13) and refined (levels 0-2): the loop bookkeeping is exact, light's spectrum on the closed
space comes out as (n + 1)^2 with 6, 16, 30 modes, each level splits into two hands, and round 10's charge-1/2 twist,
placed on the 3-sphere, is one lowest-frequency photon of one hand (the Hopf field).
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

# ================================================================== E. Light as loops on the law's 3-sphere
say("\n== E. Electromagnetism as loops on the 3-sphere the law grew")
MESH, MAT = {}, {}
for L in range(3):
    MESH[L] = mesh(L); MAT[L] = whitney(MESH[L])

# E1 the loop bookkeeping is exact
rowsE1 = []; okE1 = True
for L in range(3):
    m = MESH[L]; nv, ne, nf, nt = len(m["X"]), len(m["E"]), len(m["F"]), len(m["T"])
    z1 = abs(m["d1"] @ m["d0"]).max(); z2 = abs(m["d2"] @ m["d1"]).max()
    row = dict(level=L, corners=nv, edges=ne, triangles=nf, pieces=nt, euler=nv - ne + nf - nt, d1d0=float(z1), d2d1=float(z2),
               helicity_asym=float(abs(MAT[L][2] - MAT[L][2].T).max()))
    if L <= 1:
        r0 = np.linalg.matrix_rank(m["d0"].toarray()); r1 = np.linalg.matrix_rank(m["d1"].toarray())
        row["b0"] = int(nv - r0); row["b1"] = int(ne - r0 - r1)
    rowsE1.append(row)
    okE1 &= z1 == 0 and z2 == 0 and row["euler"] == 0 and row["helicity_asym"] < 1e-14
okE1 &= all(r["b0"] == 1 and r["b1"] == 0 for r in rowsE1 if "b1" in r)
rng17 = np.random.default_rng(17); m0 = MESH[1]
a_rand = rng17.normal(size=len(m0["E"])); phi = rng17.normal(size=len(m0["X"]))
gauge = float(abs(m0["d1"] @ (a_rand + m0["d0"] @ phi) - m0["d1"] @ a_rand).max())
check("E1 (exact) the loop bookkeeping on the law's pieces: going round a triangle after a gradient gives zero (d1 d0 = 0) and the flux out of a closed piece is zero (d2 d1 = 0), in integers; so a change of gauge leaves every loop's flux unchanged. The glued space has one piece-connected part and no loop that isn't a boundary (b0 = 1, b1 = 0, Euler number 0): the topology of the 3-sphere, so light has no zero modes. The handedness form is symmetric",
      okE1 and gauge < 1e-12,
      "; ".join(f"level {r['level']}: {r['corners']:,} corners, {r['edges']:,} edges, {r['triangles']:,} triangles, {r['pieces']:,} pieces, Euler {r['euler']}" + (f", b0 = {r['b0']}, b1 = {r['b1']}" if "b1" in r else "") for r in rowsE1)
      + f"; gauge change moves loop fluxes by {gauge:.0e}; helicity form symmetric to {max(r['helicity_asym'] for r in rowsE1):.0e}")
results["bookkeeping"] = rowsE1

# E2 light's spectrum on the closed space, and E3 its two hands
spec = {}
CL = [(0, 6, 4), (6, 22, 9), (22, 52, 16)]                     # cluster bounds and the smooth 3-sphere's (n + 1)^2
for L in range(3):
    K, M, H, M0, V = MAT[L]
    vals, vecs = sla.eigsh(K.tocsc(), k=52, M=M.tocsc(), sigma=12.5, which='LM')
    o = np.argsort(vals); vals, vecs = vals[o], vecs[:, o]
    hel = [np.linalg.eigvalsh(vecs[:, a0:a1].T @ H @ vecs[:, a0:a1]) for a0, a1, _ in CL]
    spec[L] = dict(vals=vals, vecs=vecs, hel=hel, mean=[float(vals[a0:a1].mean()) for a0, a1, _ in CL],
                   spread=[float(np.ptp(vals[a0:a1])) for a0, a1, _ in CL], gap=[float(vals[a1] - vals[a1 - 1]) if a1 < 52 else None for a0, a1, _ in CL],
                   vol=float(V))
    say(f"   level {L}: photon levels {', '.join(f'{x:.4f}' for x in spec[L]['mean'])}  ({time.time() - T_START:.0f} s)")
err = [[spec[L]["mean"][i] / CL[i][2] - 1 for i in range(3)] for L in range(3)]
clean = all(spec[L]["spread"][i] < 0.05 * spec[L]["gap"][i] for L in (1, 2) for i in range(2))
check("E2 light on the closed space: the lowest modes come in clusters of 6, 16 and 30, as on the smooth 3-sphere, at frequency^2 -> (n + 1)^2 = 4, 9, 16 (frequency (n + 1)/a for a universe of radius a), the error falling by about 4 per refinement",
      clean and all(abs(err[2][i]) < 0.01 for i in range(3)) and all(abs(err[L][i]) > abs(err[L + 1][i]) for L in (0, 1) for i in range(3))
      and 3.3 < err[1][0] / err[2][0] < 4.7,
      "; ".join(f"level {L}: {', '.join(f'{x:.4f}' for x in spec[L]['mean'])}" for L in range(3))
      + f"; error at level 2 {', '.join(f'{100 * e:.2f}%' for e in err[2])}; clusters of 6, 16, 30 separated by gaps {', '.join(f'{g:.2f}' for g in spec[2]['gap'][:2])}")
hands = []
for L in (1, 2):
    for i, (a0, a1, w) in enumerate(CL):
        h = spec[L]["hel"][i]; npos = int((h > 0).sum()); nneg = int((h < 0).sum())
        hands.append(dict(level=L, cluster=w, plus=npos, minus=nneg, mean_plus=float(h[h > 0].mean()), mean_minus=float(h[h < 0].mean())))
okE3 = all(r["plus"] == r["minus"] for r in hands) and all(abs(r["mean_plus"] - math.sqrt(r["cluster"])) < 0.05 * math.sqrt(r["cluster"]) and abs(r["mean_minus"] + math.sqrt(r["cluster"])) < 0.05 * math.sqrt(r["cluster"]) for r in hands if r["level"] == 2)
check("E3 two hands: within each level the handedness form splits the modes evenly into left and right, with curl = +-(n + 1): 3 + 3 at +-2, 8 + 8 at +-3, 15 + 15 at +-4. Light on the 3-sphere has both hands equally",
      okE3, "; ".join(f"level {r['level']}, (n+1)^2 = {r['cluster']}: {r['plus']} at {r['mean_plus']:+.3f}, {r['minus']} at {r['mean_minus']:+.3f}" for r in hands))
results["spectrum"] = {str(L): dict(levels=spec[L]["mean"], spread=spec[L]["spread"], vol=spec[L]["vol"]) for L in spec}
results["hands"] = hands

# E4 the law's twist is the Hopf fibre: its boosts at q = +1 form the 3-sphere, its twist turns along the fibre
def Bh(u):
    return np.array([[1, u], [-np.conj(u), 1]]) / np.sqrt(1 + abs(u)**2)
def op2(a, b):
    return (a + b) / (1 - np.conj(a) * b)
def sph2(z):
    r2 = abs(z)**2; return np.array([2 * z.real, 2 * z.imag, 1 - r2]) / (1 + r2)
SX, SY, SZ = np.array([[0, 1], [1, 0]]), np.array([[0, -1j], [1j, 0]]), np.diag([1.0, -1.0])
def hopf_pt(g):
    P_ = g @ SZ @ g.conj().T; return np.array([np.trace(P_ @ s).real / 2 for s in (SX, SY, SZ)])
eU = eC = eP = eF = 0.0
for _ in range(500):
    u, v = complex(*rng17.normal(size=2) * 0.8), complex(*rng17.normal(size=2) * 0.8)
    c = 1 - np.conj(u) * v; fib = np.diag([np.conj(c) / abs(c), c / abs(c)])          # diag(T^1/2, T^-1/2), T = c*/c
    g = Bh(u); eU = max(eU, abs(g @ g.conj().T - np.eye(2)).max() + abs(np.linalg.det(g) - 1))
    eC = max(eC, abs(Bh(u) @ Bh(v) - Bh(op2(u, v)) @ fib).max())
    eP = max(eP, abs(hopf_pt(g) - sph2(-np.conj(u))).max())                           # projection: the velocity point (mirrored chart)
    eF = max(eF, abs(hopf_pt(g @ fib) - hopf_pt(g)).max())
caps = []
for theta0 in (0.3, 0.9, 1.5, 2.2):
    N = 3000; r0 = math.tan(theta0 / 2); pts = [r0 * np.exp(2j * math.pi * k / N) for k in range(N + 1)]
    law_phase = 0.0; hol = 0.0
    for k in range(N):
        z1 = pts[k]; z2 = op2(-z1, pts[k + 1]); cc = 1 - np.conj(z1) * z2
        law_phase += float(np.angle(np.conj(cc) / cc)) / 2                                # round 10: the law's common phase, half the twist
        psi1 = Bh(pts[k])[:, 0]; psi2 = Bh(pts[k + 1])[:, 0]
        hol += float(np.angle(np.vdot(psi1, psi2)))                                       # the Hopf connection along the lifted loop
    area = 2 * math.pi * (1 - math.cos(theta0))
    caps.append(dict(cap=theta0, area=area, law=law_phase, hopf=hol))
okcaps = all(abs(abs(r["hopf"]) - r["area"] / 2) < 1e-5 and abs(abs(r["law"]) - r["area"] / 2) < 1e-5 for r in caps)
check("E4 the law's twist is the Hopf fibre: at q = +1 the law's boosts, divided by their size, are the unit 3-sphere SU(2); composing two gives the boost to u (+) v times a turn along the fibre by half the twist; projecting a boost to the velocity sphere gives its velocity, and the fibre turn doesn't move it. Round a loop of velocities, the Hopf connection's phase equals round 10's charge-1/2 phase: half the enclosed area. So round 10's charge 1/2 is the Hopf bundle (one unit of flux)",
      eU < 1e-12 and eC < 1e-12 and eP < 1e-12 and eF < 1e-12 and okcaps,
      f"500 random pairs: unitary {eU:.0e}, composition {eC:.0e}, projection {eP:.0e}, fibre {eF:.0e}; "
      + "; ".join(f"cap {r['cap']}: half area {r['area'] / 2:.5f}, law {abs(r['law']):.5f}, Hopf {abs(r['hopf']):.5f}" for r in caps))
results["hopf_law"] = dict(unitary=eU, composition=eC, projection=eP, fibre=eF, caps=caps)

# E5 placed on the law's 3-sphere, that connection is a lowest photon of one hand
def hopf_field(X, hand=1):                   # A = x1 dx2 - x2 dx1 + hand (x3 dx4 - x4 dx3) = Im(z1* dz1 + hand z2* dz2)
    return np.stack([-X[:, 1], X[:, 0], -hand * X[:, 3], hand * X[:, 2]], 1)
def hopf_map(X):
    z1 = X[:, 0] + 1j * X[:, 1]; z2 = X[:, 2] + 1j * X[:, 3]; w = z1 * np.conj(z2)
    return np.stack([2 * w.real, 2 * w.imag, abs(z1)**2 - abs(z2)**2], 1)
def sph_area(A_, B_, C_):
    return 2 * np.arctan2(np.einsum('ij,ij->i', A_, np.cross(B_, C_)), 1 + (A_ * B_).sum(1) + (B_ * C_).sum(1) + (C_ * A_).sum(1))
hf = []
for L in range(3):
    m = MESH[L]; K, M, H, M0, V = MAT[L]; X, E = m["X"], m["E"]; V0 = spec[L]["vecs"][:, :6]
    row = dict(level=L)
    for hand, name in ((1, "hopf"), (-1, "mirror")):
        Af = hopf_field(X, hand); a = (0.5 * (Af[E[:, 0]] + Af[E[:, 1]]) * (X[E[:, 1]] - X[E[:, 0]])).sum(1)   # exact edge integrals
        Ma = M @ a; nrm2 = a @ Ma; cvec = V0.T @ Ma
        row[name] = dict(freq2=float(a @ (K @ a) / nrm2), curl=float(a @ (H @ a) / nrm2), overlap=float(cvec @ cvec / nrm2))
        if hand == 1:
            flux = m["d1"] @ a; F = m["F"]; Y = hopf_map(X); ar = sph_area(Y[F[:, 0]], Y[F[:, 1]], Y[F[:, 2]])
            sg = np.sign((flux * ar).sum()); row["flux_err"] = float(np.abs(flux - sg * ar / 2).sum() / np.abs(flux).sum()); row["flux_sign"] = float(sg)
    hf.append(row)
okE5 = (abs(hf[2]["hopf"]["overlap"] - 1) < 1e-3 and abs(hf[2]["mirror"]["overlap"] - 1) < 1e-3
        and abs(abs(hf[2]["hopf"]["curl"]) - 2) < 0.02 and hf[2]["hopf"]["curl"] * hf[2]["mirror"]["curl"] < 0
        and abs(hf[2]["hopf"]["freq2"] - 4) < 0.05
        and abs(hf[1]["hopf"]["curl"] + 2) < abs(hf[0]["hopf"]["curl"] + 2) * 1.01)
check("E5 placed on the law's 3-sphere, that connection (the Hopf field) is light: it lies entirely in the lowest level (frequency 2/a) and has one hand, curl -2 with the pieces oriented outward; its mirror image has the other hand, +2. So the law's charge-1/2 twist, carried onto the 3-sphere the law grew, is one of the three lowest photons of one handedness",
      okE5, "; ".join(f"level {r['level']}: overlap with the lowest level {r['hopf']['overlap']:.5f} (mirror {r['mirror']['overlap']:.5f}), frequency^2 {r['hopf']['freq2']:.4f}, curl {r['hopf']['curl']:+.4f} (mirror {r['mirror']['curl']:+.4f})" for r in hf))
# E6 its loops: flux through each triangle is half the area of its image on the velocity sphere
fe = [r["flux_err"] for r in hf]
check("E6 going round the loops of the law's pieces: the flux of the Hopf field through each triangle is half the area that triangle covers on the velocity sphere (round 10's half area, now loop by loop on the 3-sphere), the mismatch falling by about half with each refinement",
      fe[0] > fe[1] > fe[2] and fe[2] < 0.07 and 1.6 < fe[1] / fe[2] < 2.5,
      f"summed mismatch / summed flux at levels 0-2: {', '.join(f'{x:.4f}' for x in fe)} (falls by {fe[0] / fe[1]:.2f}, {fe[1] / fe[2]:.2f}); sign {hf[0]['flux_sign']:+.0f}")
results["hopf_field"] = hf

results["summary"] = {"passed": passed, "failed": failed, "seconds": round(time.time() - T_START)}
say(f"\n{passed} passed, {failed} failed  ({time.time() - T_START:.0f} s)")
json.dump(results, open("results.json", "w"), indent=1, default=float)
open("checks_output.txt", "w").write(_log.getvalue())
