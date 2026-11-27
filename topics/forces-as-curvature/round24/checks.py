"""DRAFT round 24 checks: charge as a source of light on the law's mesh.

Round 10 found the law's twist gives matter charge 1/2; round 17 put light on the 3-sphere the law grows; round 18 found the law's
own field there. Here charges are put on the law's mesh as sources: what Gauss's law allows on the law's closed space, the static
field of charges, light made by a moving charge, and whether round 10's charge 1/2 is itself a source in the law's space.
Run: python3 checks.py  (numpy, scipy).  Writes results.json and checks_output.txt.
"""
import json, io, math, itertools, time
import numpy as np
from scipy.integrate import solve_ivp, quad
from scipy.sparse.csgraph import dijkstra

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
import scipy.sparse.linalg as sla
MESH, MAT = {}, {}
for L in range(3):
    MESH[L] = mesh(L); MAT[L] = whitney(MESH[L])
say(f"meshes and Whitney matrices built  ({time.time() - T_START:.0f} s)")
def pole_pair(X):
    return int(np.argmin(np.linalg.norm(X - [1, 0, 0, 0], axis=1))), int(np.argmin(np.linalg.norm(X - [-1, 0, 0, 0], axis=1)))
def solve_L0(L0, rho):            # L0 phi = rho on the closed mesh: pin corner 0 (constants are the kernel), then remove the mean
    lu = sla.splu(L0[1:, 1:].tocsc()); ph = np.r_[0.0, lu.solve(rho[1:])]; return ph
def chi_of(X):
    return np.arccos(np.clip(X[:, 0], -1, 1))

# ================================================================== Q. Charge as a source of light on the law's mesh
say("\n== Q. Charge on the law's mesh")

# Q1: on the law's closed space, the charges that any electric field can have add up to zero
rowsQ1 = []
for L in range(3):
    m = MESH[L]; d0 = m["d0"]; K, M, H, M0, V = MAT[L]
    L0 = (d0.T @ M @ d0).tocsr()
    colsum = float(np.abs(np.asarray(d0.sum(1))).max())                     # d0 applied to a constant: 0 exactly (integers)
    E = np.random.default_rng(L).normal(size=d0.shape[0]); tot = float(abs((d0.T @ (M @ E)).sum()))   # charge of a random field
    N_, S_ = pole_pair(m["X"]); rho1 = np.zeros(d0.shape[1]); rho1[N_] = 1.0
    sol = sla.lsqr(L0, rho1, atol=1e-14, btol=1e-14, iter_lim=20000)[0]
    resid = float(np.linalg.norm(L0 @ sol - rho1))
    rowsQ1.append(dict(level=L, corners=int(d0.shape[1]), edges=int(d0.shape[0]), grad_of_constant=colsum, total_charge_random_field=tot,
                       single_charge_residual=resid, lower_bound=float(1 / math.sqrt(d0.shape[1]))))
okQ1 = all(r["grad_of_constant"] == 0 and r["total_charge_random_field"] < 1e-9 and r["single_charge_residual"] > 0.9 * r["lower_bound"] for r in rowsQ1)
check("Q1 (exact; textbook) on the law's closed space the charges of any electric field add up to zero: the discrete Gauss law reads charge = -d0^T M E, and d0 of a constant is exactly 0 (integer incidence), so the total of -d0^T M E is 0 for every field. A lone charge has no field at all: L0 phi = (one charge) has no solution, whatever is tried, it misses by at least 1/sqrt(number of corners). Charges on the law's space come in neutral sets (textbook for a closed universe)",
      okQ1, "; ".join(f"level {r['level']} ({r['corners']} corners): d0(constant) = {r['grad_of_constant']:.0f}, total charge of a random field {r['total_charge_random_field']:.0e}, best try at one lone charge misses by {r['single_charge_residual']:.4f} (bound {r['lower_bound']:.4f})" for r in rowsQ1))
results["neutral"] = rowsQ1

# Q2: a + and - pair at the law's two poles: Coulomb on the 3-sphere, phi = cot(chi)/(4 pi)
rowsQ2 = []; prof = {}
for L in range(3):
    m = MESH[L]; d0 = m["d0"]; K, M, H, M0, V = MAT[L]; X = m["X"]
    L0 = (d0.T @ M @ d0).tocsc()
    N_, S_ = pole_pair(X); rho = np.zeros(len(X)); rho[N_], rho[S_] = 1.0, -1.0
    ph = solve_L0(L0, rho); ch = chi_of(X)
    sel = (ch > 0.6) & (ch < math.pi - 0.6)
    exact = 1 / np.tan(ch) / (4 * math.pi)
    off = np.mean(ph[sel] - exact[sel]); err = float(np.abs(ph[sel] - off - exact[sel]).max() / np.abs(exact[sel]).max())
    flux = []
    for c0 in (0.8, 1.4, 2.0):                                        # charge enclosed by a cap of corners chi < c0, read from the field
        inside = (ch < c0).astype(float); flux.append(float(inside @ (L0 @ ph)))
    w = M0[sel]; err2 = float(math.sqrt((w * (ph[sel] - off - exact[sel])**2).sum() / (w * exact[sel]**2).sum()))
    rowsQ2.append(dict(level=L, err=err, err_l2=err2, enclosed=flux, antipodal=float(np.linalg.norm(X[N_] + X[S_]))))
    prof[L] = dict(chi=ch.tolist(), phi=(ph - off).tolist())
okQ2 = all(abs(f - 1) < 1e-9 for r in rowsQ2 for f in r["enclosed"]) and rowsQ2[2]["err"] < rowsQ2[1]["err"] < rowsQ2[0]["err"] and rowsQ2[2]["err_l2"] < rowsQ2[1]["err_l2"] < rowsQ2[0]["err_l2"]
check("Q2 (textbook, discrete) a + charge at the law's rest point and a - charge at its antipode (both corners the law grew) make the 3-sphere's Coulomb field: potential cot(chi)/(4 pi), so the field 1/(4 pi sin^2 chi) carries the same flux, 1, through every sphere around the + charge, which grows to the equator and shrinks again. On the law's mesh the enclosed charge is exactly 1 for every cap, and the potential approaches the smooth one with refinement",
      okQ2, "; ".join(f"level {r['level']}: potential vs cot(chi)/(4 pi) for 0.6 < chi < pi - 0.6, worst {r['err']:.3f}, average {r['err_l2']:.3f} (relative); enclosed charge for caps chi < 0.8, 1.4, 2.0: {', '.join(f'{x:.12f}' for x in r['enclosed'])}" for r in rowsQ2))
results["pair"] = rowsQ2

# Q3: a lone charge is allowed only with an equal and opposite charge spread evenly over the whole space
rowsQ3 = []; prof3 = {}
def E_bg(c):                       # smooth: enclosed charge 1 - vol(cap)/vol(S^3), over the area 4 pi sin^2 chi
    return (1 - (c - math.sin(c) * math.cos(c)) / math.pi) / (4 * math.pi * math.sin(c)**2)
for L in range(3):
    m = MESH[L]; d0 = m["d0"]; K, M, H, M0, V = MAT[L]; X = m["X"]
    L0 = (d0.T @ M @ d0).tocsc()
    N_, S_ = pole_pair(X); rho = -M0 / M0.sum(); rho[N_] += 1.0
    ph = solve_L0(L0, rho); ch = chi_of(X)
    ref = np.array([quad(E_bg, c, math.pi / 2)[0] for c in ch])     # phi(chi) - phi(pi/2) = int_chi^{pi/2} E
    sel = (ch > 0.6)
    off = np.mean(ph[sel] - ref[sel]); err = float(np.abs(ph[sel] - off - ref[sel]).max() / np.abs(ref[sel]).max())
    w = M0[sel]; err2 = float(math.sqrt((w * (ph[sel] - off - ref[sel])**2).sum() / (w * ref[sel]**2).sum()))
    rowsQ3.append(dict(level=L, total=float(rho.sum()), err=err, err_l2=err2, far_field=float(E_bg(math.pi - 1e-6))))
    prof3[L] = dict(chi=ch.tolist(), phi=(ph - off).tolist())
okQ3 = all(abs(r["total"]) < 1e-12 for r in rowsQ3) and rowsQ3[2]["err"] < rowsQ3[1]["err"] < rowsQ3[0]["err"] and rowsQ3[2]["err_l2"] < rowsQ3[1]["err_l2"] < rowsQ3[0]["err_l2"]
check("Q3 (textbook, discrete) a lone charge can sit on the law's space only with an equal and opposite charge spread evenly over all of it; then its field is (1 - (chi - sin chi cos chi)/pi)/(4 pi sin^2 chi), which dies to zero at the far pole instead of focusing there. The law's mesh reproduces it, the error falling with refinement",
      okQ3, "; ".join(f"level {r['level']}: total charge {r['total']:.0e}, potential vs smooth for chi > 0.6 worst {r['err']:.3f}, average {r['err_l2']:.3f}" for r in rowsQ3))
results["background"] = rowsQ3

# Q4: a moving charge makes light, and Gauss's law holds exactly as it moves
L = 1; m = MESH[L]; d0 = m["d0"]; K, M, H, M0, V = MAT[L]; X, Em = m["X"], m["E"]; ne, nv = len(Em), len(X)
Mlu = sla.splu(M.tocsc()); L0 = (d0.T @ M @ d0).tocsc(); L0lu = sla.splu(L0[1:, 1:].tocsc())
def transverse(v):
    rhs = d0.T @ (M @ v); s_ = L0lu.solve(rhs[1:]); return v - d0 @ np.r_[0.0, s_]
N_, S_ = pole_pair(X)
elen = np.linalg.norm(X[Em[:, 1]] - X[Em[:, 0]], axis=1)
G = sp.csr_matrix((np.r_[elen, elen], (np.r_[Em[:, 0], Em[:, 1]], np.r_[Em[:, 1], Em[:, 0]])), shape=(nv, nv))
dist, pred = dijkstra(G, indices=N_, return_predecessors=True)
target = int(np.argmin(np.abs(chi_of(X) - 1.0) + 10 * (np.abs(X[:, 2:]).sum(1) > 0.4 * np.abs(X[:, 1])) ))   # a corner near chi = 1
path = [target]
while path[-1] != N_:
    path.append(int(pred[path[-1]]))
path = path[::-1]
eid = m["eid"]
steps = [(eid[(min(a_, b_), max(a_, b_))], 1.0 if a_ < b_ else -1.0) for a_, b_ in zip(path[:-1], path[1:])]   # edge and the sign of travel along it
def schedule(t, T):                 # the charge leaves the pole at t = 0 and arrives at t = T, one edge after another, smoothly
    n = len(steps); u = np.clip(t / T, 0, 1) * n; k = min(int(u), n - 1); s = u - k
    frac = np.zeros(n); frac[:k] = 1.0; frac[k] = np.sin(math.pi * s / 2)**2 if t < T else 1.0
    dfrac = np.zeros(n)
    if 0 <= t < T:
        dfrac[k] = (math.pi / 2) * math.sin(math.pi * s) * n / T
    return frac, dfrac
def rho_of(frac):
    rho = np.zeros(nv); rho[S_] = -1.0; rho[path[0]] += 1.0
    for (e, sg), f_ in zip(steps, frac):
        a_, b_ = (Em[e, 0], Em[e, 1]) if sg > 0 else (Em[e, 1], Em[e, 0])
        rho[a_] -= f_; rho[b_] += f_
    return rho
def j_of(dfrac):
    j = np.zeros(ne)
    for (e, sg), df in zip(steps, dfrac):
        j[e] += sg * df
    return j
def run_move(T, t_after=3.0):
    rho0 = rho_of(np.zeros(len(steps))); ph0 = np.r_[0.0, L0lu.solve(rho0[1:])]
    y0 = np.concatenate([np.zeros(ne), M @ (d0 @ ph0)])            # temporal gauge: E = -dA/dt; P = M dA/dt; starts in the pair's static field
    def rhs(t, y):
        A, P = y[:ne], y[ne:]; _, df = schedule(t, T)
        return np.concatenate([Mlu.solve(P), -(K @ A) + j_of(df)])
    sol = solve_ivp(rhs, (0, T + t_after), y0, rtol=1e-10, atol=1e-12, method="DOP853", max_step=T / (8 * len(steps)), dense_output=True)
    ts = np.linspace(0, T + t_after, 120); gauss = 0.0; Erad = []
    for t in ts:
        y = sol.sol(t); A, P = y[:ne], y[ne:]; fr, _ = schedule(t, T)
        gauss = max(gauss, float(np.abs(d0.T @ P - rho_of(fr)).max()))
        Ad = Mlu.solve(P); AT, AdT = transverse(A), transverse(Ad)
        Erad.append(0.5 * (AdT @ (M @ AdT) + AT @ (K @ AT)))
    Erad = np.array(Erad); after = ts > T + 0.05
    rho1 = rho_of(np.ones(len(steps))); ph1 = np.r_[0.0, L0lu.solve(rho1[1:])]
    yT = sol.sol(T + t_after); Ad = Mlu.solve(yT[ne:])
    static_err = float(np.linalg.norm(Ad - transverse(Ad) - d0 @ ph1) / np.linalg.norm(d0 @ ph1))
    return dict(T=T, gauss=gauss, E_rad_after=float(Erad[after].mean()), E_rad_drift=float(np.ptp(Erad[after]) / Erad[after].mean()),
                static_after=static_err, t=ts.tolist(), E_rad=Erad.tolist())
rowsQ4 = [run_move(T) for T in (0.5, 1.0, 2.0, 4.0)]
Es = [r["E_rad_after"] for r in rowsQ4]
okQ4 = (all(r["gauss"] < 1e-7 and r["E_rad_drift"] < 1e-6 and r["static_after"] < 1e-6 for r in rowsQ4) and all(a > b for a, b in zip(Es, Es[1:])) and Es[0] > 1e-4)
check("Q4 (textbook, discrete) a moving charge makes light on the law's mesh. A + charge walks from the rest point along the mesh's edges to a corner about 1 radian away while its partner stays at the antipode; the field follows Maxwell's equations with the charge's current as source. Gauss's law holds exactly at every moment; once the charge stops, the gradient part of the field is the new pair's static field and the rest, a gradient-free field of fixed energy, is light that keeps going round the closed space. Quicker moves make more light (Larmor: radiation grows with acceleration)",
      okQ4, f"level {L}, path of {len(steps)} edges to chi = {chi_of(X)[target]:.3f}; " + "; ".join(f"move in time {r['T']}: Gauss {r['gauss']:.0e}, light's energy after {r['E_rad_after']:.2e} (kept to {r['E_rad_drift']:.0e}), static part vs new pair {r['static_after']:.0e}" for r in rowsQ4))
results["moving"] = dict(level=L, path_edges=len(steps), end_chi=float(chi_of(X)[target]), runs=rowsQ4)

# Q5: is round 10's charge 1/2 itself a source in the law's space? The law's own field has no source there
def hopf_field(X, hand=1):
    return np.stack([-X[:, 1], X[:, 0], -hand * X[:, 3], hand * X[:, 2]], 1)
rowsQ5 = []
for L in range(3):
    m = MESH[L]; d0, d1, d2 = m["d0"], m["d1"], m["d2"]; K, M, H, M0, V = MAT[L]; X, Em = m["X"], m["E"]
    Af = hopf_field(X, -1); a = (0.5 * (Af[Em[:, 0]] + Af[Em[:, 1]]) * (X[Em[:, 1]] - X[Em[:, 0]])).sum(1)
    div = d0.T @ (M @ a); scale = np.abs(d0).T @ np.abs(M @ a)
    rowsQ5.append(dict(level=L, electric_source=float(np.abs(div).max() / scale.max()), magnetic_source=float(np.abs(d2 @ (d1 @ a)).max()),
                       flux_scale=float(np.abs(d1 @ a).max())))
al = np.linspace(0, 2 * math.pi, 2001); fib = np.stack([np.cos(al), np.sin(al), 0 * al, 0 * al], 1)
circ = float(np.trapezoid((hopf_field(fib, -1) * np.stack([-np.sin(al), np.cos(al), 0 * al, 0 * al], 1)).sum(1), al))
Jm = np.array([[0, -1, 0, 0], [1, 0, 0, 0], [0, 0, 0, 1], [0, 0, -1, 0]], float)     # the law's field is X -> J X (hand -1)
smooth_div = float(abs(np.trace(Jm))) + float(np.abs(Jm + Jm.T).max())              # antisymmetric J: a rotation, divergence-free on S^3
okQ5 = all(r["magnetic_source"] < 1e-12 and r["electric_source"] < 1e-3 for r in rowsQ5) and smooth_div == 0 and abs(circ - 2 * math.pi) < 1e-9
check("Q5 (exact for the magnetic part; ours as a reading) round 10's charge 1/2 is not a source in the law's space. The law's own field (round 18) has no magnetic source anywhere in the law's 3-space (the flux out of every piece is exactly 0) and no electric source (smoothly it is X -> J X with J antisymmetric, a rotation of the 3-sphere, so divergence-free exactly; on the mesh its divergence is 0 to rounding at levels 0 and 1 and 1.5e-4 at level 2). Its 2 pi of flux, one Dirac unit, appears only on the sphere of directions it projects to (the Hopf base; the circulation round each fibre is 2 pi). So charge 1/2 is how matter couples to the twist, a charge on the sphere of directions, not a charge sitting in space; and Q1 forbids a lone charge on the law's closed space in any case",
      okQ5, "; ".join(f"level {r['level']}: magnetic source {r['magnetic_source']:.0e}, electric source {r['electric_source']:.1e} (relative)" for r in rowsQ5) + f"; smooth: J antisymmetric, trace {np.trace(Jm):.0f}; circulation round a fibre {circ:.10f} = 2 pi")
results["charge_half"] = dict(levels=rowsQ5, fibre_circulation=circ)
results["profiles"] = dict(pair=prof[2], background=prof3[2])

say(f"\n{passed} passed, {failed} failed  ({time.time() - T_START:.0f} s)")
results["summary"] = dict(passed=passed, failed=failed, seconds=round(time.time() - T_START, 1))
json.dump(results, open("results.json", "w"), indent=1)
open("checks_output.txt", "w").write(_log.getvalue())
