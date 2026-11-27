"""DRAFT round 22 checks: light, then matter, in the law's de Sitter spacetime.

Round 20 built de Sitter spacetime from the law (one dial q, a dot product with one minus sign): its slices are 3-spheres of
radius cosh(tau), the law's own expanding space. Round 17 put light on one such 3-sphere, frozen in time; round 18 found the
law's own field there (the 'mirror' Hopf field, curl +2). Here the slices are allowed to grow: Maxwell is read back on the law's
spacetime, light's redshift and horizon are shown, the law's mesh carries the same light, and then free matter is added.
Run: python3 checks.py  (numpy, scipy, sympy).  Writes results.json and checks_output.txt.
"""
import json, io, math, itertools, time
import numpy as np
import sympy as spy
from scipy.integrate import solve_ivp

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
for _ in range(1):
    LEVELS.append(refine_with_law(*LEVELS[-1], 1.0))
say(f"corners grown by the law (round 13): {n}; space pieces at levels 0-1: {[len(T) for _, T in LEVELS]}  ({time.time() - T_START:.0f} s)")

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
# ------------------------------------------------------------------ the law with a Lorentzian dot product (round 20, unchanged)
ETA = np.array([-1.0, 1, 1, 1])
def dot4(a, b):
    return (a * b * ETA).sum(-1, keepdims=True)
def oplus4(a, b, qq=1.0):
    ab = dot4(a, b); aa = dot4(a, a); bb = dot4(b, b)
    return ((1 - 2 * qq * ab - qq * bb) * a + (1 + qq * aa) * b) / (1 - 2 * qq * ab + qq**2 * aa * bb)
def place(x, qq=1.0):             # round 20's placement into 5D (W, T, X, Y, Z), on W^2 - T^2 + X^2 + Y^2 + Z^2 = 1
    e = dot4(x, x); return np.concatenate([(1 - qq * e), 2 * math.sqrt(qq) * x], -1) / (math.sqrt(qq) * (1 + qq * e))
def slice_coords(P):              # global time tau (T = sinh tau) and the angle chi from the pole (1, 0, 0, 0) on that slice's 3-sphere
    tau = np.arcsinh(P[:, 1]); S = P[:, [0, 2, 3, 4]]
    chi = np.arctan2(np.linalg.norm(S[:, 1:], axis=1), S[:, 0]); return tau, chi
gd = lambda t: 2 * np.arctan(np.tanh(t / 2))            # Gudermannian: conformal time eta(t) = int_0^t dt'/cosh t'

# ================================================================== M. Light and matter in the law's de Sitter
say("\n== M. Light, then matter, in the law's de Sitter spacetime")

# M1: Maxwell on the law's spacetime is flat Maxwell in the law's chart (conformal invariance, textbook in 4D)
t_, x_, y_, z_, q_ = spy.symbols('t x y z q', real=True)
X4 = [t_, x_, y_, z_]; eta_m = spy.diag(-1, 1, 1, 1)
ee = -t_**2 + x_**2 + y_**2 + z_**2
Om = 2 / (1 + q_ * ee)
g_m = Om**2 * eta_m; ginv = g_m.inv(); sqrtg = Om**4                     # sqrt(-det g) = Omega^4 in 4D
det_ok = spy.simplify(-g_m.det() - sqrtg**2) == 0
dens = [[spy.simplify(sqrtg * ginv[a, c] * ginv[b, d] - eta_m[a, c] * eta_m[b, d]) for d in range(4) for c in range(4)] for b in range(4) for a in range(4)]
weyl_ok = all(v == 0 for row in dens for v in row)
f = spy.Function('f'); u_ = t_ - x_
A_low = [0, 0, f(u_), 0]                                                  # a light wave of any profile f moving along x, polarised along y
F = spy.Matrix(4, 4, lambda a, b: spy.diff(A_low[b], X4[a]) - spy.diff(A_low[a], X4[b]))
Fup = ginv * F * ginv
div = [spy.simplify(sum(spy.diff(sqrtg * Fup[a, b], X4[a]) for a in range(4)) / sqrtg) for b in range(4)]
bianchi = [spy.simplify(spy.diff(F[1, 2], X4[0]) + spy.diff(F[2, 0], X4[1]) + spy.diff(F[0, 1], X4[2]))]
wave_ok = all(d == 0 for d in div) and all(b == 0 for b in bianchi)
check("M1 (exact, symbolic; textbook) Maxwell's equations on the law's spacetime, metric 4 eta(dx, dx)/(1 + q eta(x, x))^2 (round 20), are exactly the flat ones in the law's chart: sqrt(-g) g^ac g^bd = eta^ac eta^bd for every q, so the curvature drops out of the light equations entirely. A light wave of any profile moving along a light ray of the chart solves them on the curved spacetime. This is Maxwell's conformal invariance in 4D (Bateman 1910); here the conformal factor is the law's own",
      det_ok and weyl_ok and wave_ok,
      f"sqrt(-g) = Omega^4: {det_ok}; sqrt(-g) g^ac g^bd - eta^ac eta^bd = 0 for all 256 entries: {weyl_ok}; wave f(t - x) dy: div F = {div}, dF = {bianchi}")
results["maxwell_chart"] = dict(sqrtg_is_Omega4=bool(det_ok), conformal_identity=bool(weyl_ok), wave_solves=bool(wave_ok))

# M2: on the expanding 3-spheres, each photon mode redshifts as 1/a and the photon count is kept
ts, ks = spy.symbols('t k', positive=True)
a_s = spy.cosh(ts); gd_s = 2 * spy.atan(spy.tanh(ts / 2))
dgd = spy.simplify(spy.diff(gd_s, ts) - 1 / a_s)
A_s = spy.cos(ks * gd_s)
ode = spy.simplify((spy.diff(a_s * spy.diff(A_s, ts), ts) + ks**2 * A_s / a_s).rewrite(spy.exp))
ode_num = max(abs(float(((spy.diff(a_s * spy.diff(A_s, ts), ts) + ks**2 * A_s / a_s)).subs({ts: tv, ks: kv}))) for tv in (0.3, 1.1, 2.7) for kv in (2, 3, 7))
En_s = spy.Rational(1, 2) * (a_s * spy.diff(A_s, ts)**2 + ks**2 * A_s**2 / a_s)
Ea_num = max(abs(float((En_s * a_s).subs({ts: tv, ks: kv})) - kv**2 / 2) for tv in (0.0, 0.3, 1.1, 2.7, 5.0) for kv in (2, 3, 4, 12))
rowsM2 = []
for k in (2, 3, 4, 12, 40):                                       # integrate the mode equation directly, not using the closed form
    rhs = lambda t, y: [y[1] / math.cosh(t), -k**2 * y[0] / math.cosh(t)]          # y = (A, a dA/dt)
    tt = np.linspace(0, 6, 3001)
    sol = solve_ivp(rhs, (0, 6), [1.0, 0.0], t_eval=tt, rtol=1e-11, atol=1e-13, method="DOP853")
    A, P = sol.y; a = np.cosh(tt); E = 0.5 * (P**2 / a + k**2 * A**2 / a)
    iz = np.where(np.sign(A[:-1]) != np.sign(A[1:]))[0]
    zt = [tt[i] - A[i] * (tt[i + 1] - tt[i]) / (A[i + 1] - A[i]) for i in iz]       # zeros of the integrated mode, by interpolation
    per = [(zt[i + 1] - zt[i], math.pi * math.cosh(0.5 * (zt[i] + zt[i + 1])) / k) for i in range(len(zt) - 1)]
    per = [(p_, r_) for p_, r_ in per if r_ < 0.5]                 # while the wave is short against the horizon (half period < 0.5)
    a_short = max([math.cosh(0.5 * (zt[i] + zt[i + 1])) for i in range(len(zt) - 1) if math.pi * math.cosh(0.5 * (zt[i] + zt[i + 1])) / k < 0.5], default=float("nan"))
    rowsM2.append(dict(k=k, closed_form=float(np.abs(A - np.cos(k * gd(tt))).max()), Ea_drift=float(np.abs(E * a / E[0] - 1).max()),
                       E_end_over_E0=float(E[-1] / E[0]), a_end=float(a[-1]), zeros=len(zt), cycles_left_total=k / 4,
                       half_period_vs_pi_a_over_k=float(max([abs(p / r - 1) for p, r in per], default=float("nan"))), short_wave_up_to_a=a_short))
okM2 = (dgd == 0 and ode == 0 and ode_num < 1e-9 and Ea_num < 1e-9 and all(r["closed_form"] < 1e-8 and r["Ea_drift"] < 1e-8 for r in rowsM2)
        and rowsM2[-1]["half_period_vs_pi_a_over_k"] < 0.01)
check("M2 (exact; textbook) on the law's slices (3-spheres of radius a = cosh t, round 20 D6) a photon mode of frequency k on the unit 3-sphere (k = 2 is round 17's lowest) obeys d/dt(a dA/dt) + k^2 A / a = 0, solved by A = cos(k eta) with eta = gd(t) = int dt / a. Its energy is exactly k^2/(2a): it falls as 1/a, the redshift, while its frequency seen locally is k/a, so energy / frequency (the photon count) never changes. Light from the bounce has only k/4 cycles left in the whole future",
      okM2, f"d eta/dt = 1/a: {dgd == 0}; mode equation symbolic {ode == 0}, numeric {ode_num:.0e}; E*a - k^2/2 {Ea_num:.0e}; direct integration to t = 6 (a = {rowsM2[0]['a_end']:.0f}): "
      + "; ".join(f"k {r['k']}: vs cos(k eta) {r['closed_form']:.0e}, E*a drift {r['Ea_drift']:.0e}, E(6)/E(0) = {r['E_end_over_E0']:.5f} = 1/a, zeros {r['zeros']}" for r in rowsM2)
      + f"; k = 40, local half period vs pi a/k while shorter than 0.5 (up to a = {rowsM2[-1]['short_wave_up_to_a']:.1f}): worst {rowsM2[-1]['half_period_vs_pi_a_over_k']:.4f}")
results["redshift_modes"] = rowsM2

# M3: the law's light ray from the bounce crosses exactly a quarter of the 3-sphere: the horizon
n = np.array([1.0, 1.0, 0, 0])                                     # a light direction of the chart (eta(n, n) = 0)
s = np.concatenate([np.linspace(0, 5, 401)[1:], np.geomspace(5, 1e7, 200)])
tau, chi = slice_coords(place(s[:, None] * n))
e_gd = np.abs(chi - gd(tau)).max(); chi_far = chi[-1]
rowsM3 = []
for tau0 in (0.5, 1.0, 2.0):
    p = np.array([[math.tanh(tau0 / 2), 0, 0, 0]])                 # the law's time kick of rapidity tau0 from rest (round 20 D6)
    tp, cp = slice_coords(place(p)); 
    ss = np.geomspace(1e-4, 1e7, 4000)
    with np.errstate(all="ignore"):
        R = place(oplus4(np.tile(p, (len(ss), 1)), ss[:, None] * n))
    ok = np.all(np.isfinite(R), 1); tr, cr = slice_coords(R[ok])
    rowsM3.append(dict(tau0=tau0, start_tau=float(tp[0]), gd_track=float(np.abs(cr - (gd(tr) - gd(tau0))).max()),
                       reach=float(cr[np.argmax(tr)]), predicted=float(math.pi / 2 - gd(tau0)), on_ds=float(np.abs((R[ok]**2 * [1, -1, 1, 1, 1]).sum(1) - 1).max())))
okM3 = e_gd < 1e-12 and abs(chi_far - math.pi / 2) < 1e-6 and all(r["gd_track"] < 1e-9 and abs(r["reach"] - r["predicted"]) < 1e-6 and abs(r["start_tau"] - r["tau0"]) < 1e-12 for r in rowsM3)
check("M3 (exact) the law's own light ray (along a light direction the law is plain addition, round 20 D3) is light's path: the chart line s(1, 1, 0, 0) from rest places at angle chi = gd(tau) on the slice at time tau, so it has travelled exactly the conformal time of M2, and as tau -> infinity it reaches chi = pi/2, a quarter of the way round the 3-sphere and no further. Light sent later, at tau0, by adding the ray to the law's time kick (the law's addition, an isometry), reaches pi/2 - gd(tau0). The far hemisphere is never reached: de Sitter's event horizon, textbook, here traced by the law's own sums",
      okM3, f"from the bounce: chi vs gd(tau) {e_gd:.0e}, chi at s = 1e7 {chi_far:.8f} (pi/2 = {math.pi / 2:.8f}); "
      + "; ".join(f"sent at tau0 {r['tau0']}: tracks gd(tau) - gd(tau0) to {r['gd_track']:.0e}, reaches {r['reach']:.6f} = pi/2 - gd(tau0) {r['predicted']:.6f}" for r in rowsM3))
results["horizon"] = dict(from_bounce_vs_gd=float(e_gd), chi_far=float(chi_far), later=rowsM3)

# M4: the same light on the law's own mesh (round 18's level 1, Whitney forms), with the slices growing
m1 = mesh(1); K, Mm, H, M0, V = whitney(m1); Xm, Em, d0 = m1["X"], m1["E"], m1["d0"]
ne = len(Em)
import scipy.sparse.linalg as sla
Mlu = sla.splu(Mm.tocsc())
def hopf_field(X, hand=1):                       # round 17's fields; round 18: the law's own transport follows hand = -1 ('mirror')
    return np.stack([-X[:, 1], X[:, 0], -hand * X[:, 3], hand * X[:, 2]], 1)
Af = hopf_field(Xm, -1); a_h = (0.5 * (Af[Em[:, 0]] + Af[Em[:, 1]]) * (Xm[Em[:, 1]] - Xm[Em[:, 0]])).sum(1)
lam_h = float(a_h @ (K @ a_h) / (a_h @ (Mm @ a_h))); curl_h = float(a_h @ (H @ a_h) / (a_h @ (Mm @ a_h)))
L0 = (d0.T @ Mm @ d0).tocsc()                      # scalar Laplacian, to split off gradients; constants are its kernel, so pin corner 0
L0lu = sla.splu(L0[1:, 1:].tocsc())
def transverse(v):                                # remove the gradient part (M-orthogonal projection)
    rhs_ = d0.T @ (Mm @ v); sol_ = L0lu.solve(rhs_[1:]); ph = np.concatenate([np.zeros((1,) + sol_.shape[1:]), sol_])
    return v - d0 @ ph
rgm = np.random.default_rng(22)
phi = rgm.normal(size=d0.shape[1])
def run(A0, Adot0, t_end, conformal=False):
    """cosmic time: y = (A, P), P = a M dA/dt; dA/dt = M^-1 P / a, dP/dt = -K A / a  (from the action int a Ad M Ad/2 - A K A/(2a)).
    conformal=True: the static 3-sphere's wave equation M A'' + K A = 0 in eta, with no a anywhere."""
    def rhs(t, y):
        A, P = y[:ne], y[ne:]; a = 1.0 if conformal else math.cosh(t)
        return np.concatenate([Mlu.solve(P) / a, -(K @ A) / a])
    y0 = np.concatenate([A0, Mm @ Adot0])
    return solve_ivp(rhs, (0, t_end), y0, rtol=1e-10, atol=1e-12, method="DOP853", dense_output=True)
T_END = 4.0
tt = np.linspace(0, T_END, 81)
# run 1: the law's field at the bounce, plus a pure-gauge gradient, at rest
A0 = a_h + d0 @ phi
r1 = run(A0, np.zeros(ne), T_END)
Y1 = r1.sol(tt); A1 = Y1[:ne]
proj = (a_h @ (Mm @ A1)) / (a_h @ (Mm @ a_h))
e_proj = float(np.abs(proj - np.cos(math.sqrt(lam_h) * gd(tt))).max())
grad_part = A1 - transverse(A1)
e_gauge = float(np.abs(grad_part - (d0 @ phi)[:, None]).max() / np.abs(d0 @ phi).max())
# run 2: a generic field (random, gradient-free, moving) evolved in cosmic time, against the static 3-sphere in conformal time
B0 = transverse(rgm.normal(size=ne)); B0 /= math.sqrt(B0 @ (Mm @ B0)); Bd0 = transverse(rgm.normal(size=ne)) * 0.3
r2 = run(B0, Bd0, T_END); r2c = run(B0, Bd0, float(gd(T_END)), conformal=True)
Y2 = r2.sol(tt); A2, P2 = Y2[:ne], Y2[ne:]
A2c = r2c.sol(gd(tt))[:ne]
e_conf = float(np.abs(A2 - A2c).max() / np.abs(A2).max())
a_t = np.cosh(tt)
E2 = np.array([0.5 * (P2[:, i] @ Mlu.solve(P2[:, i]) + A2[:, i] @ (K @ A2[:, i])) / a_t[i] for i in range(len(tt))])
e_Ea = float(np.abs(E2 * a_t / (E2[0] * a_t[0]) - 1).max())
gauss = float(np.abs(d0.T @ P2).max() / np.abs(P2).max())
results["mesh_light"] = dict(level=1, edges=int(ne), law_field_freq2=lam_h, law_field_curl=curl_h, mode_vs_cos=e_proj, gauge_kept=e_gauge,
                             cosmic_vs_conformal=e_conf, E_times_a_drift=e_Ea, gauss=gauss, E_end_over_E0=float(E2[-1] / E2[0]), a_end=float(a_t[-1]),
                             proj_t=tt.tolist(), proj=proj.tolist(), E2=E2.tolist())
okM4 = e_proj < 2e-3 and e_gauge < 1e-6 and e_conf < 1e-6 and e_Ea < 1e-7 and gauss < 1e-8 and abs(curl_h - 2) < 0.05
check("M4 the same light on the law's own mesh (round 18's level 1: corners grown and cut by the law, Whitney forms). Discrete Maxwell with the slices growing, M d/dt(a dA/dt) + K A / a = 0, is exactly the frozen 3-sphere's wave equation in conformal time eta = gd(t), so the discrete light inherits M2's redshift. The law's own field (round 18's 'mirror' Hopf field, curl +2) left at rest at the bounce follows cos(sqrt(lambda) gd(t)); a gradient (pure gauge) never moves; a generic field's energy times a stays fixed; Gauss's law (no charge) is kept",
      okM4, f"{ne} edges; law's field: frequency^2 {lam_h:.4f} (smooth: 4), curl {curl_h:+.4f}; its amplitude vs cos(sqrt(lambda) gd(t)) to t = {T_END} (a = {a_t[-1]:.1f}) {e_proj:.0e}; gradient kept to {e_gauge:.0e}; "
      f"generic field: cosmic-time run vs static run in eta {e_conf:.0e}, E*a drift {e_Ea:.0e}, E({T_END})/E(0) = {E2[-1] / E2[0]:.4f} vs 1/a = {1 / a_t[-1]:.4f}, Gauss {gauss:.0e}")

# M5: matter: the law's tilted time lines are free particles, and their momentum falls as 1/a
rowsM5 = []
for beta in (0.5, 1.0, 2.0):
    u = np.concatenate([np.linspace(0, 0.99, 2000), 1 - np.geomspace(1e-2, 1e-7, 400)[1:]])
    v = np.array([math.cosh(beta), math.sinh(beta), 0, 0])         # a time-like direction: the law boosts along it (round 21)
    P = place(u[:, None] * v); tau_, chi_ = slice_coords(P)
    plane_rank = np.linalg.svd(P - 0, compute_uv=False); planarity = float(plane_rank[2] / plane_rank[0])
    # also the law's own sums along the line: (s v) (+) (w v) = ((s + w)/(1 + s w)) v  (law at -q along a time-like direction)
    sv, wv = 0.3, 0.45; lhs = oplus4((sv * v)[None], (wv * v)[None])[0]; e_line = float(np.abs(lhs - (sv + wv) / (1 + sv * wv) * v).max())
    c_, s_ = math.cosh(beta), math.sinh(beta); D = 1 - u**2                       # exact derivatives of the placed line in u
    W, Tt, Xx = (1 + u**2) / D, 2 * u * c_ / D, 2 * u * s_ / D
    dW, dT, dX = 4 * u / D**2, 2 * c_ * (1 + u**2) / D**2, 2 * s_ * (1 + u**2) / D**2
    dchi = ((W * dX - Xx * dW) / (W**2 + Xx**2)) / (dT / np.sqrt(1 + Tt**2))
    vel = np.cosh(tau_) * dchi; gam = 1 / np.sqrt(1 - vel**2)
    pa = gam * vel * np.cosh(tau_)
    sel = (tau_ > 0.05) & (tau_ < 6)
    rowsM5.append(dict(beta=beta, planarity=planarity, line_sum=e_line, p_times_a_dev=float(np.abs(pa[sel] / math.sinh(beta) - 1).max()),
                       chi_end=float(chi_[-1]), gd_beta=float(gd(beta)), speed_at_bounce=float(vel[0]), tanh_beta=math.tanh(beta),
                       tau=tau_[sel][::40].tolist(), p=(gam * vel)[sel][::40].tolist()))
okM5 = all(r["planarity"] < 1e-12 and r["line_sum"] < 1e-12 and r["p_times_a_dev"] < 1e-9 and abs(r["chi_end"] - r["gd_beta"]) < 1e-4 and abs(r["speed_at_bounce"] - r["tanh_beta"]) < 1e-9 for r in rowsM5)
check("M5 matter: a chart line through rest along a time-like direction (where the law boosts, round 21) is a free particle. Its placed points lie in a 2-plane through the 5D origin, so it is a de Sitter geodesic (textbook), and the law's sums along it are rapidity addition. A particle leaving the pole at the bounce with rapidity beta has momentum p = gamma v falling as sinh(beta)/cosh(tau): p a fixed, the same 1/a as light (M2). It slows to rest with the expansion at angle gd(beta) < pi/2; light (beta -> infinity) is the limit pi/2 of M3",
      okM5, "; ".join(f"beta {r['beta']}: off the plane {r['planarity']:.0e}, law's sum along it {r['line_sum']:.0e}, speed at bounce {r['speed_at_bounce']:.4f} = tanh beta {r['tanh_beta']:.4f}, p a vs sinh beta {r['p_times_a_dev']:.0e}, final angle {r['chi_end']:.5f} = gd(beta) {r['gd_beta']:.5f}" for r in rowsM5))
results["matter"] = rowsM5

say(f"\n{passed} passed, {failed} failed  ({time.time() - T_START:.0f} s)")
results["summary"] = dict(passed=passed, failed=failed, seconds=round(time.time() - T_START, 1))
json.dump(results, open("results.json", "w"), indent=1)
open("checks_output.txt", "w").write(_log.getvalue())
