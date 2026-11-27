"""DRAFT round 14 checks: real 4D spacetime, with time inside the glued pieces (Lorentzian Regge calculus).

Rounds 11-13 glued flat 3D pieces into a closed space and read Einstein's equations back with Regge's method,
with the corners grown by the twisted-boost law (../twisted-boosts/) and the pieces refined by its half kicks.
Here each piece also lasts a while: a space piece at size a, the same piece at size a' a moment later, and the
flat 4D slab between them in Minkowski spacetime (a frustum). Gluing the slabs makes a closed expanding
universe built from flat pieces of spacetime, and the curvature sits on two kinds of hinge:
  - timelike hinges (a space edge swept through time), where the deficit is an ordinary angle, and
  - spacelike hinges (a space triangle at one moment), where the deficit is a boost.
This is textbook Lorentzian Regge calculus (Regge 1961, Collins & Williams 1973, Sorkin 1975); what is ours is that
the corners come from the law's kicks, as in round 13. Every angle here is computed from the pieces themselves.
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

# ------------------------------------------------------------------ the law in 3D and the corners it grows (rounds 12-13)
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
for _ in range(3):
    LEVELS.append(refine_with_law(*LEVELS[-1], 1.0))
say(f"corners grown by the law (round 13): {n}; space pieces at levels 0-3: {[len(T) for _, T in LEVELS]}  ({time.time() - T_START:.0f} s)")

# ------------------------------------------------------------------ flat pieces of Minkowski spacetime
ETA = np.array([-1.0, 1, 1, 1, 1])                     # R^{1,4}: time first
def mink(u, v):
    return (u * v * ETA).sum(-1)
def proj_out(w, H):
    """remove from w the part lying in the plane spanned by the rows of H (Minkowski projection)"""
    G = np.einsum('...ik,...jk->...ij', H * ETA, H)
    b = np.einsum('...ik,...k->...i', H * ETA, w)
    return w - np.einsum('...i,...ik->...k', np.linalg.solve(G, b[..., None])[..., 0], H)
def tri_area(P0, P1, P2):
    H = np.stack([P1 - P0, P2 - P0], -2)
    return 0.5 * np.sqrt(np.abs(np.linalg.det(np.einsum('...ik,...jk->...ij', H * ETA, H))))
def simplex_vol(Q):                                    # Q (..., 5, 5): 4-volume of a flat Lorentzian 4-simplex
    H = Q[..., 1:, :] - Q[..., :1, :]
    return np.sqrt(np.abs(np.linalg.det(np.einsum('...ik,...jk->...ij', H * ETA, H)))) / 24
def corners(X, a, ap, tau):
    """the slab over a space piece with corners X (n,4,4): bottom (0, a x), top (dt, a' x), struts of length^2 = -tau^2"""
    dt = math.sqrt(tau**2 + (ap - a)**2)
    B = np.concatenate([np.zeros(X.shape[:-1] + (1,)), a * X], -1)
    U = np.concatenate([np.full(X.shape[:-1] + (1,), dt), ap * X], -1)
    return B, U
def frustum(X, a, ap, tau):
    """everything Regge needs from one slab, computed from its corners alone:
    th (n,6): ordinary dihedral angle at each timelike hinge (edge ij swept in time), At (n,6) its area;
    eB, eT (n,4): the boost by which the slab's wall leans away from a right angle at its bottom and top
    triangles (the triangle opposite corner l), AB, AT (n,4) their areas; V (n,) the slab's 4-volume."""
    B, U = corners(X, a, ap, tau)
    th, At = [], []
    for i, j in PAIRS:
        k, l = [m for m in range(4) if m not in (i, j)]
        H = np.stack([B[:, j] - B[:, i], U[:, i] - B[:, i]], 1)     # the hinge: a space edge and a strut
        u1, u2 = proj_out(B[:, k] - B[:, i], H), proj_out(B[:, l] - B[:, i], H)
        th.append(np.arccos(np.clip(mink(u1, u2) / np.sqrt(mink(u1, u1) * mink(u2, u2)), -1, 1)))
        At.append(tri_area(B[:, i], B[:, j], U[:, j]) + tri_area(B[:, i], U[:, j], U[:, i]))
    eB, AB, eT, AT = [], [], [], []
    for l in range(4):
        f = [m for m in range(4) if m != l]
        for Pp, Qq, oe, oA in ((B, U, eB, AB), (U, B, eT, AT)):
            H = np.stack([Pp[:, f[1]] - Pp[:, f[0]], Pp[:, f[2]] - Pp[:, f[0]]], 1)   # the hinge: a space triangle
            uh = proj_out(Pp[:, l] - Pp[:, f[0]], H); uh = uh / np.sqrt(mink(uh, uh))[:, None]   # into the space piece
            v = proj_out(Qq[:, f[0]] - Pp[:, f[0]], H)                                  # along the wall (timelike)
            beta = mink(v, uh); w = v - beta[:, None] * uh; alpha = np.sqrt(-mink(w, w))
            oe.append(-np.arctanh(beta / alpha)); oA.append(tri_area(Pp[:, f[0]], Pp[:, f[1]], Pp[:, f[2]]))
    Q = [np.stack([B[:, 0], B[:, 1], B[:, 2], B[:, 3], U[:, 3]], 1), np.stack([B[:, 0], B[:, 1], B[:, 2], U[:, 2], U[:, 3]], 1),
         np.stack([B[:, 0], B[:, 1], U[:, 1], U[:, 2], U[:, 3]], 1), np.stack([B[:, 0], U[:, 0], U[:, 1], U[:, 2], U[:, 3]], 1)]
    V = sum(simplex_vol(q) for q in Q)
    return np.array(th).T, np.array(At).T, np.array(eB).T, np.array(AB).T, np.array(eT).T, np.array(AT).T, V

say("\n== L. One slab of spacetime")
X0 = to_s3(LEVELS[1][0])[LEVELS[1][1][:200]]          # 200 real space pieces (level 1)
a, ap, tau = 1.0, 1.08, 0.13
B, U = corners(X0, a, ap, tau)
strut = mink(U - B, U - B)
plan = max(np.linalg.matrix_rank(np.stack([B[k, j] - B[k, i], U[k, i] - B[k, i], U[k, j] - B[k, i]]), tol=1e-10)
           for k in range(len(X0)) for i, j in PAIRS)
# the 4-volume, a second way: the slab is a stack of scaled space pieces, so V = vol3 * height * int (scale)^3
G = X0[:, 1:] - X0[:, :1]
vol3 = np.sqrt(np.abs(np.linalg.det(np.einsum('nik,njk->nij', G, G)))) / 6
nrm = np.linalg.svd(np.concatenate([X0, np.ones((len(X0), 4, 1))], -1))[2][:, -1]   # the piece's hyperplane n.x = h
h = np.abs(nrm[:, 4]) / np.linalg.norm(nrm[:, :4], axis=1)
dt = math.sqrt(tau**2 + (ap - a)**2)
Vclosed = vol3 * np.sqrt(dt**2 - (ap - a)**2 * h**2) * (ap**4 - a**4) / (4 * (ap - a))
V = frustum(X0, a, ap, tau)[6]
errV = float(np.abs(V / Vclosed - 1).max())
check("L1 each slab is a flat piece of Minkowski spacetime: its 4 struts all have length^2 = -tau^2 (they last a proper time tau), its walls are flat, and its 4-volume from 4 glued simplices equals the stacked-pieces formula",
      np.abs(strut + tau**2).max() < 1e-14 and plan == 2 and errV < 1e-12,
      f"200 pieces; strut length^2 {strut.mean():.6f} (= -tau^2 = {-tau**2:.6f}); wall faces have rank {plan}; 4-volume agrees to {errV:.1e}")

rng = np.random.default_rng(14)
Xr = rng.normal(size=(40, 4, 4))                        # 40 slabs over random lopsided pieces, nothing symmetric
Xr = Xr / np.linalg.norm(Xr, axis=-1, keepdims=True) * (1 + 0.15 * rng.normal(size=(40, 4, 1)))
def schlafli(Xr, sgn, h=1e-6):
    p0 = np.array([1.0, 1.3, 0.7]); dp = rng.normal(size=3); dX = rng.normal(size=Xr.shape)
    th0, At0, eB0, AB0, eT0, AT0, _ = frustum(Xr, *p0)
    thp, _, eBp, _, eTp, _, _ = frustum(Xr + h * dX, *(p0 + h * dp))
    thm, _, eBm, _, eTm, _, _ = frustum(Xr - h * dX, *(p0 - h * dp))
    tl = (At0 * (thp - thm)).sum(1); sl = (AB0 * (eBp - eBm)).sum(1) + (AT0 * (eTp - eTm)).sum(1)
    return np.abs(tl + sgn * sl) / (np.abs(tl) + np.abs(sl))
sch_minus, sch_plus = schlafli(Xr, -1).max(), schlafli(Xr, +1).min()
check("L2 Schlafli's identity holds in every slab, with time inside: sum of (area x change of ordinary angle) over the timelike hinges equals sum of (area x change of boost) over the spacelike ones. This fixes how the two kinds of curvature add, with nothing fitted",
      sch_minus < 1e-7 and sch_plus > 0.5,
      f"40 random slabs, random wiggles of every corner: off by at most {sch_minus:.1e} of either side; with the other sign it misses by at least {sch_plus * 100:.0f}%")
results["schlafli_rel_err"] = float(sch_minus)

# ------------------------------------------------------------------ whole slabs of the universe, one for each level of refinement
class Slabs:
    """the slab of universe between two moments, over every space piece of one level. Pieces of the same shape give the
    same numbers, so they are grouped (the grouping is checked against the full sum in L3)."""
    def __init__(self, C, T):
        Xs = to_s3(C)[T]
        Gm = np.einsum('nik,njk->nij', Xs, Xs)
        keys = np.array([np.round(Gm[:, p][:, :, p].reshape(len(T), 16), 9) for p in itertools.permutations(range(4))])
        # canonical key: the lexicographically smallest relabelling
        best = keys[0].copy()
        for k in keys[1:]:
            less = np.zeros(len(T), bool); decided = np.zeros(len(T), bool)
            for col in range(16):
                lt = (k[:, col] < best[:, col]) & ~decided; gt = (k[:, col] > best[:, col]) & ~decided
                less |= lt; decided |= lt | gt
            best[less] = k[less]
        _, first, self.mult = np.unique(best, axis=0, return_index=True, return_counts=True)
        self.X = Xs[first]; self.Xall = Xs; self.T = T
        E = np.concatenate([np.sort(T[:, list(p)], 1) for p in PAIRS])
        self.ue = np.unique(E, axis=0)
        L = np.linalg.norm(to_s3(C)[self.ue[:, 0]] - to_s3(C)[self.ue[:, 1]], axis=1)
        self.elen, self.ecount = np.unique(np.round(L, 10), return_counts=True)
    def parts(self, a, ap, tau, full=False):
        """(timelike part, spacelike part, 4-volume) of the slab's Regge action, with this slab's own share of
        the spacelike hinges at its two ends (the other share belongs to the neighbouring slab)"""
        X = self.Xall if full else self.X; w = 1 if full else self.mult
        th, At, eB, AB, eT, AT, V = frustum(X, a, ap, tau)
        # 2 pi x (area of every timelike hinge), minus every slab's angle there: the deficit, edge by edge
        Xe = np.zeros((len(self.elen), 4, 4)); Xe[:, 0, 0] = 1
        s = np.sqrt(np.clip(1 - (self.elen / 2) ** 2, 0, 1))     # an edge of this chord length on the unit sphere
        Xe[:, 1] = np.stack([1 - self.elen**2 / 2, self.elen * s, 0 * s, 0 * s], -1)
        Bq, Uq = corners(Xe, a, ap, tau)
        Ae = tri_area(Bq[:, 0], Bq[:, 1], Uq[:, 1]) + tri_area(Bq[:, 0], Uq[:, 1], Uq[:, 0])
        tl = 2 * np.pi * (self.ecount * Ae).sum() - (w * (At * th).sum(1)).sum()
        sl = (w * ((AB * eB).sum(1) + (AT * eT).sum(1))).sum()   # Schlafli (L2) fixed this sign
        return tl, sl, (w * V).sum()
    def action(self, a, ap, tau, lam=0.0):
        tl, sl, V = self.parts(a, ap, tau); return tl + sl - lam * V
SL = [Slabs(*lv) for lv in LEVELS]
say(f"space pieces grouped by shape: {[len(s.X) for s in SL]} shapes at levels 0-3  ({time.time() - T_START:.0f} s)")

say("\n== M. Gluing slabs into an expanding universe")
tgt = 6 * math.pi**2
tau = 0.05
full1 = SL[1].parts(1.0, 1.01, tau, full=True); grp1 = SL[1].parts(1.0, 1.01, tau)
static = [SL[L].parts(1.0, 1.0, tau) for L in range(4)]
st_ratio = [s[0] / (tgt * tau) for s in static]
st_err = [abs(r - 1) for r in st_ratio]
r11 = [57.131, 58.643, 59.070, 59.181]                  # round 11's space-only Regge sums for the same kind of slices
check("M1 a slab that doesn't change size: no boosts anywhere, and its action is (proper time) x (round 11's space curvature sum), tending to Einstein's 6 pi^2 a tau at second order",
      max(abs(s[1]) for s in static) < 1e-12 and all(abs(s[0] / tau - r) < 2e-3 for s, r in zip(static, r11))
      and all(3.3 < st_err[i] / st_err[i + 1] < 4.6 for i in range(3)) and abs(grp1[0] + grp1[1] - full1[0] - full1[1]) < 1e-9,
      "; ".join(f"level {L}: {s[0] / tau:.3f} x tau, {r * 100:.3f}% of 6 pi^2" for L, (s, r) in enumerate(zip(static, st_ratio)))
      + f"; grouping by shape matches the full sum to {abs(grp1[0] + grp1[1] - full1[0] - full1[1]):.0e}")
kin = []
for L in range(4):
    row = []
    for d in (0.005, 0.01, 0.02):
        tl, sl, _ = SL[L].parts(1 - d / 2, 1 + d / 2, tau)
        row.append((tl + sl - static[L][0]) / (-tgt * d * d / tau))
    kin.append(row)
kin_err = [abs(r[1] - 1) for r in kin]
tl_alone = SL[3].parts(0.995, 1.005, tau)
check("M2 a growing slab: the spacelike hinges' boosts give the kinetic term -6 pi^2 a (da)^2 / tau, with Einstein's minus sign, tending to it at second order (it is quadratic in the growth rate)",
      all(r[1] > 0 for r in kin) and all(3.3 < kin_err[i] / kin_err[i + 1] < 4.6 for i in range(3)) and kin_err[-1] < 5e-3
      and max(abs(r[0] - r[2]) for r in kin[1:]) < 3e-3,
      "; ".join(f"level {L}: {r[1] * 100:.2f}% of it" for L, r in enumerate(kin))
      + f"; the same to within {max(abs(r[0] - r[2]) for r in kin[1:]) * 100:.2f}% for growth rates da/tau of 0.1, 0.2 and 0.4")
results["static_ratio"] = st_ratio; results["kinetic_ratio"] = kin

say("\n== N. Einstein's equations with a cosmological constant, read back from the glued slabs")
LAM = 3.0                                              # units where the de Sitter throat is at a = 1
def constraint(S, a, ap, tau):                         # vary the struts: d(action)/d(tau) = 0
    h = 1e-5 * tau
    return (S.action(a, ap, tau + h, LAM) - S.action(a, ap, tau - h, LAM)) / (2 * h)
throat = []
for L in range(4):
    tl, _, V = SL[L].parts(1.0, 1.0, tau)
    throat.append(math.sqrt(tl / (LAM * V)))           # a static slab solves the constraint where a tl = Lambda a^3 V
lam_r11 = [1.14085, 1.03649, 1.00917, 1.00230]
th_err = [abs(x - 1) for x in throat]
check("N1 varying the struts (the lapse) gives the Friedmann constraint; a slab at rest solves it at the throat a^2 = 3/Lambda, tending to it at second order, the same number as round 11's Lambda rho^2",
      all(abs(x**2 - r) < 2e-4 for x, r in zip(throat, lam_r11)) and all(3.3 < th_err[i] / th_err[i + 1] < 4.6 for i in range(3)),
      "; ".join(f"level {L}: throat {x:.5f}" for L, x in enumerate(throat)))
results["throat"] = throat

def grow_universe(S, a0, tau, steps, from_rest=False):
    """grow the universe slab by slab: given a slice of size a, the Friedmann constraint (N1) fixes the next"""
    A = [a0]
    for k in range(steps):
        a = A[-1]
        lo = a + (0.3 * tau**2 if (from_rest and k == 0) else 1e-12)   # from rest: skip the trivial root (staying at rest)
        hi = a + 3 * tau * math.sqrt(max(a * a - 1, 0)) + 3 * tau**2
        if constraint(S, a, lo, tau) * constraint(S, a, hi, tau) > 0:
            break                                       # no growing solution: the glued slabs have stopped following de Sitter
        A.append(brentq(lambda x: constraint(S, a, x, tau), lo, hi, xtol=1e-13))
    return np.array(A)
runs = {}
for L, tau_ in ((1, 0.05), (2, 0.05), (3, 0.05), (3, 0.1), (3, 0.025)):
    A = grow_universe(SL[L], math.cosh(0.5), tau_, int(round(1.0 / tau_)))
    runs[(L, tau_)] = (0.5 + tau_ * np.arange(len(A)), A)
gap = {k: float(np.max(np.abs(A / np.cosh(t) - 1))) for k, (t, A) in runs.items()}
g_lv = [gap[(L, 0.05)] for L in (1, 2, 3)]
check("N2 start the glued universe at the size de Sitter has at t = 0.5 and let the Friedmann constraint grow it slab by slab to t = 1.5: it follows de Sitter's a = cosh(t sqrt(Lambda/3)), the gap falling about 4x per level of refinement; at level 3 what is left comes from the space pieces, not the time step",
      all(3.3 < g_lv[i] / g_lv[i + 1] < 5.5 for i in range(2)) and g_lv[-1] < 1e-3 and abs(gap[(3, 0.025)] / gap[(3, 0.05)] - 1) < 0.2,
      "worst gap with tau = 0.05: " + ", ".join(f"level {L} {g * 100:.3f}%" for L, g in zip((1, 2, 3), g_lv))
      + f"; level 3 with tau = 0.1, 0.05, 0.025: {gap[(3, 0.1)] * 100:.3f}%, {gap[(3, 0.05)] * 100:.3f}%, {gap[(3, 0.025)] * 100:.3f}%")
results["de_sitter_gap"] = {f"level {L} tau {tt}": v for (L, tt), v in gap.items()}

def a_residual(S, A, tau):
    """the other Regge equations: vary each slice's size a_n with the struts held. Not used to build the solution."""
    res, scale = [], []
    for k in range(1, len(A) - 1):
        h = 1e-6
        tot = lambda x: S.action(A[k - 1], x, tau, LAM) + S.action(x, A[k + 1], tau, LAM)
        one = lambda x: S.action(x, A[k + 1], tau, LAM)
        res.append((tot(A[k] + h) - tot(A[k] - h)) / (2 * h)); scale.append(abs((one(A[k] + h) - one(A[k] - h)) / (2 * h)))
    return float(np.max(np.abs(res) / np.array(scale)))
ra = {k: a_residual(SL[k[0]], runs[k][1], k[1]) for k in runs}
check("N3 the remaining Regge equations (vary each slice's size with the struts held), never used to grow the solution, hold too: what is left is set by the time step alone, the same at every level and falling about 8x when the step halves",
      max(ra[(L, 0.05)] for L in (1, 2, 3)) / min(ra[(L, 0.05)] for L in (1, 2, 3)) < 1.3
      and 6.5 < ra[(3, 0.1)] / ra[(3, 0.05)] < 9.5 and 6.5 < ra[(3, 0.05)] / ra[(3, 0.025)] < 9.5,
      "; ".join(f"level {L}, tau {tt}: {v * 100:.4f}% of one slab's term" for (L, tt), v in ra.items()))
results["a_equation_residual"] = {f"level {L} tau {tt}": v for (L, tt), v in ra.items()}

rest = {}
for tau_ in (0.1, 0.05):
    A = grow_universe(SL[3], throat[3], tau_, int(round(1.5 / tau_)), from_rest=True)
    t = tau_ * np.arange(len(A))
    rest[tau_] = (float(np.max(np.abs(A / (throat[3] * np.cosh(t / throat[3])) - 1))), (A[1] - A[0]) / tau_**2)
check("N4 caveat: started exactly at rest at the throat, the slabs still follow de Sitter, but the gap only halves when the time step halves, because the first step off the bounce comes out tau^2 where de Sitter has tau^2 / 2",
      1.6 < rest[0.1][0] / rest[0.05][0] < 2.4 and all(abs(r[1] - 1) < 0.1 for r in rest.values()),
      f"level 3, out to t = 1.5: worst gap {rest[0.1][0] * 100:.2f}% with tau = 0.1, {rest[0.05][0] * 100:.2f}% with tau = 0.05; first step {rest[0.1][1]:.3f} and {rest[0.05][1]:.3f} x tau^2")
results["from_rest"] = {str(k): v for k, v in rest.items()}

# pieces against the horizon: de Sitter's horizon is sqrt(3/Lambda) = 1 here, and a piece's size is a x its edge
edge = [float(np.mean(np.linalg.norm(to_s3(Cl)[Tl[:, 0]] - to_s3(Cl)[Tl[:, 1]], axis=1))) for Cl, Tl in LEVELS]
long1 = grow_universe(SL[1], math.cosh(0.5), 0.1, 30)
t1 = 0.5 + 0.1 * np.arange(len(long1)); g1 = np.abs(long1 / np.cosh(t1) - 1); size1 = long1 * edge[1]
small = g1[size1 < 0.5].max(); broke = len(long1) - 1 < 30
check("N5 the glued slabs follow de Sitter only while each piece is small next to the horizon (sqrt(3/Lambda)): at level 1 the gap stays under 1% while pieces are under half the horizon, then grows, and the constraint has no growing solution once pieces are about 2 horizons across",
      small < 0.01 and broke and 1.5 < size1[-1] < 2.5,
      f"level 1 piece edge {edge[1]:.3f} x a; gap {small * 100:.2f}% while pieces < 0.5 horizon; {g1[-1] * 100:.0f}% at the last step (t = {t1[-1]:.1f}), pieces {size1[-1]:.2f} horizons; level 3 pieces are {edge[3]:.3f} x a")
results["horizon_run_level1"] = {"t": list(map(float, t1)), "a": list(map(float, long1)), "piece_over_horizon": list(map(float, size1))}
results["edge_over_a"] = edge
t, A3 = runs[(3, 0.025)]
results["run_level3_tau0.025"] = {"t": list(map(float, t)), "a": list(map(float, A3))}

results["summary"] = {"passed": passed, "failed": failed, "seconds": round(time.time() - T_START)}
say(f"\n{passed} passed, {failed} failed  ({time.time() - T_START:.0f} s)")
json.dump(results, open("results.json", "w"), indent=1, default=float)
open("checks_output.txt", "w").write(_log.getvalue())
