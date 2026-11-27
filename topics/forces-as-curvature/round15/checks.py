"""DRAFT round 15 checks: boost loops on round 14's glued slabs, and time as the fused shift.

Part T. Round 14 found that at a spacelike hinge (a space triangle at one moment) the curvature is a loop of boosts.
Here those loops are gone round, on the real slabs: per hinge, over two hinges, over a whole space piece, and under
refinement, and compared with the twisted-boost law at q = -1 (../twisted-boosts/). The identity it confirms, that
composing boosts gives a boost times a Thomas-Wigner rotation exactly as Moebius addition's gyration predicts, is
known (Ungar); what is checked is that it holds on Regge's pieces, and whether the twist does anything there.
Part S. Time as the fused shift a -> 1 + k a, which turns a o_k b = a + b + k a b into multiplication
(fused_shift in ../riemann/round8/lean/). Each slice's size is read through the shift, to see what it supplies
(the growth, the bounce, de Sitter's length, the breakdown) and what still goes in by hand.
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
for _ in range(3):
    LEVELS.append(refine_with_law(*LEVELS[-1], 1.0))
say(f"corners grown by the law (round 13): {n}; space pieces at levels 0-3: {[len(T) for _, T in LEVELS]}  ({time.time() - T_START:.0f} s)")

# ------------------------------------------------------------------ flat pieces of Minkowski spacetime (round 14, unchanged)
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

# ------------------------------------------------------------------ whole slabs of the universe (round 14, unchanged)
class Slabs:
    """the slab of universe between two moments, over every space piece of one level. Pieces of the same shape give the
    same numbers, so they are grouped (round 14 checked the grouping against the full sum)."""
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
LAM = 3.0                                              # units where the de Sitter throat is at a = 1
def constraint(S, a, ap, tau):                         # vary the struts: d(action)/d(tau) = 0
    h = 1e-5 * tau
    return (S.action(a, ap, tau + h, LAM) - S.action(a, ap, tau - h, LAM)) / (2 * h)
throat = []
for L in range(4):
    tl, _, V = SL[L].parts(1.0, 1.0, 0.05)
    throat.append(math.sqrt(tl / (LAM * V)))           # a static slab solves the constraint where a tl = Lambda a^3 V
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

# ================================================================== T. Boost loops on round 14's glued slabs
say("\n== T. Going round the spacelike hinges: boost loops, composed by the law at q = -1")
def boost(v):
    """Lorentz boost (t, x, y, z) with Einstein velocity v (..., 3), |v| < 1"""
    v = np.atleast_2d(v); s2 = (v * v).sum(-1); g = 1 / np.sqrt(1 - s2)
    B = np.zeros(v.shape[:-1] + (4, 4)); B[..., 0, 0] = g; B[..., 0, 1:] = g[..., None] * v; B[..., 1:, 0] = g[..., None] * v
    vv = np.einsum('...i,...j->...ij', v, v)
    B[..., 1:, 1:] = np.eye(3) + ((g - 1) / np.where(s2 > 0, s2, 1))[..., None, None] * vv
    return B
def rot4(G):
    R = np.zeros(G.shape[:-2] + (4, 4)); R[..., 0, 0] = 1; R[..., 1:, 1:] = G; return R
def gyr(a, b):
    """the law's gyration at q = -1 as a 3x3 matrix: gyr[a,b] c = -(a (+) b) (+) (a (+) (b (+) c)) (Ungar)"""
    ab = oplus(a, b, -1.0); cols = []
    for e in np.eye(3):
        c = np.tile(e, a.shape[:-1] + (1,))
        cols.append(oplus(-ab, oplus(a, oplus(b, c, -1.0), -1.0), -1.0))
    return np.stack(cols, -1)
def half(v):            # Einstein velocity -> the law's half velocity (round 10's half-velocities): |m| = tanh(eta / 2)
    s = np.linalg.norm(v, axis=-1, keepdims=True); eta = np.arctanh(s)
    return np.where(s > 0, np.tanh(eta / 2) / np.where(s > 0, s, 1), 0) * v
def doubled(m):         # half velocity -> Einstein velocity
    return 2 * m / (1 + (m * m).sum(-1, keepdims=True))

def hinge_data(L):
    """every spacelike hinge (space triangle) in the middle slice: its boost deficit, summed from the 4 slab walls
    that meet there (2 slabs below, 2 above), and for every slab above, the hinges on its floor in that slab's own frame"""
    C_, T_ = LEVELS[L]; X = to_s3(C_)[T_]
    TAU_T = 0.2 / 2**L                                           # the time step shrinks with the pieces
    A3 = (math.cosh(0.5 - TAU_T), math.cosh(0.5), math.cosh(0.5 + TAU_T))   # three slices of de Sitter, one step apart
    _, _, _, _, eT, _, _ = frustum(X, A3[0], A3[1], TAU_T)       # walls of the slabs below, at their tops
    _, _, eB, _, _, _, _ = frustum(X, A3[1], A3[2], TAU_T)       # walls of the slabs above, at their floors
    faces = np.stack([np.sort(T_[:, [m for m in range(4) if m != l]], 1) for l in range(4)], 1)   # (nt,4,3)
    uf, inv = np.unique(faces.reshape(-1, 3), axis=0, return_inverse=True); inv = inv.ravel().reshape(-1, 4)
    wedges = np.zeros((len(uf), 4)); fill = np.zeros(len(uf), int)
    for t in range(len(T_)):
        for l in range(4):
            f = inv[t, l]; wedges[f, fill[f]] = eT[t, l]; wedges[f, fill[f] + 1] = eB[t, l]; fill[f] += 2
    eps = wedges.sum(1)
    # each slab above: an orthonormal frame of its floor (the space piece), and the inward direction to each floor triangle
    P = A3[1] * X; D = P[:, 1:] - P[:, :1]
    Q = np.linalg.qr(np.transpose(D, (0, 2, 1)))[0]               # (nt,4,3) orthonormal basis of the piece's 3-plane
    U = np.zeros((len(T_), 4, 3))
    for l in range(4):
        f = [m for m in range(4) if m != l]
        e1, e2 = P[:, f[1]] - P[:, f[0]], P[:, f[2]] - P[:, f[0]]
        w = P[:, l] - P[:, f[0]]
        H = np.stack([e1, e2], 1); G2 = np.einsum('nik,njk->nij', H, H)
        cfs = np.linalg.solve(G2, np.einsum('nik,nk->ni', H, w)[..., None])[..., 0]
        u = w - np.einsum('ni,nik->nk', cfs, H); u /= np.linalg.norm(u, axis=1)[:, None]
        U[:, l] = np.einsum('nkj,nk->nj', Q, u)                    # in the piece's own frame
    return eps, wedges, inv, U
HD = {L: hinge_data(L) for L in range(3)}
eps, wedges, inv, U = HD[1]
m_w = np.tanh(wedges / 2)[..., None] * np.array([1.0, 0, 0])      # the 4 wedge boosts share one plane: one axis
fold = m_w[:, 0]
for w in range(1, 4):
    fold = oplus(fold, m_w[:, w], -1.0)
err_fold = float(np.abs(fold[:, 0] - np.tanh(eps / 2)).max())
g_coll = float(np.abs(gyr(m_w[:, 0], m_w[:, 1]) - np.eye(3)).max())
check("T1 per hinge: the 4 slab walls meeting at a spacelike hinge boost in the one plane at right angles to it, so the law at q = -1 composes them with no twist, and the composite is exactly the hinge's boost deficit",
      err_fold < 1e-14 and g_coll < 1e-14,
      f"{len(eps):,} hinges at level 1; law-composed half velocity vs tanh(deficit / 2): {err_fold:.1e}; twist (gyration - identity): {g_coll:.1e}; deficits {eps.min():.2e} to {eps.max():.2e}")

def law_compose(ms):
    """compose boosts given as the law's half velocities, in order: returns (half velocity, rotation) with
    B(v1) B(v2) ... = B(doubled(m)) R(G)"""
    m, G = ms[0], np.broadcast_to(np.eye(3), ms[0].shape[:-1] + (3, 3)).copy()
    for nxt in ms[1:]:
        gn = np.einsum('...ij,...j->...i', G, nxt)
        G = np.einsum('...ij,...jk->...ik', gyr(m, gn), G); m = oplus(m, gn, -1.0)
    return m, G
rows_T = []
for L in range(3):
    eps_L, _, inv_L, U_L = HD[L]
    V = np.tanh(eps_L[inv_L])[..., None] * U_L                    # (nt,4,3): each floor hinge's boost, in the slab's frame
    Mh = half(V)
    # two hinges: go round one, then the other
    P2 = np.einsum('nij,njk->nik', boost(V[:, 0]), boost(V[:, 1]))
    m2, G2 = law_compose([Mh[:, 0], Mh[:, 1]])
    e2 = float(np.abs(P2 - np.einsum('nij,njk->nik', boost(doubled(m2)), rot4(G2))).max())
    # the whole floor: round all 4 hinges of the piece
    P4 = boost(V[:, 0])
    for w in range(1, 4):
        P4 = np.einsum('nij,njk->nik', P4, boost(V[:, w]))
    m4, G4 = law_compose([Mh[:, w] for w in range(4)])
    e4 = float(np.abs(P4 - np.einsum('nij,njk->nik', boost(doubled(m4)), rot4(G4))).max())
    ang = np.arccos(np.clip((np.trace(G2, axis1=1, axis2=2) - 1) / 2, -1, 1))
    sin12 = np.linalg.norm(np.cross(U_L[:, 0], U_L[:, 1]), axis=1)
    small = 0.5 * np.abs(eps_L[inv_L[:, 0]] * eps_L[inv_L[:, 1]]) * sin12   # small-boost Wigner angle: half e1 e2 sin(angle)
    big = small > 1e-3 * small.max()
    rows_T.append(dict(level=L, hinges=len(eps_L), eps_max=float(np.abs(eps_L).max()), twist_max=float(ang.max()),
                       twist_vs_small=float(np.max(np.abs(ang - small)[big] / small[big])), ratio=float(ang.max() / np.abs(eps_L).max()), err2=e2, err4=e4))
check("T2 two hinges, and a whole space piece's 4 floor hinges, gone round from one slab: the boost loop is a boost times a twist, and both are exactly the law at q = -1 composing the hinges' half velocities. This is Ungar's known identity (Thomas-Wigner rotation = gyration), here on Regge's real pieces",
      max(r["err2"] for r in rows_T) < 1e-12 and max(r["err4"] for r in rows_T) < 1e-12,
      "; ".join(f"level {r['level']}: two hinges {r['err2']:.0e}, four hinges {r['err4']:.0e}" for r in rows_T))
tw = [r["twist_max"] for r in rows_T]; ep = [r["eps_max"] for r in rows_T]
check("T3 under refinement the twist is second order in the deficits (half e1 e2 sin of the angle between the hinges) and shrinks faster than they do, so it carries nothing the smooth limit needs: an expanding universe has no twist to supply",
      all(r["twist_vs_small"] < 0.02 for r in rows_T) and tw[2] / tw[0] < (ep[2] / ep[0]) ** 1.5,
      "; ".join(f"level {r['level']}: largest deficit {r['eps_max']:.2e}, largest twist {r['twist_max']:.2e} rad (small-boost formula to {r['twist_vs_small'] * 100:.2f}%)" for r in rows_T))
results["boost_loops"] = rows_T

# ================================================================== S. Time as the fused shift
say("\n== S. Time as the fused shift a -> 1 + k a")
from fractions import Fraction as Fr
def fuse(a, b, k):                                    # the fused operation a o_k b = a + b + k a b (riemann/round8)
    return a + b + k * a * b
rs = np.random.default_rng(15)
ok_alg = True
for _ in range(200):
    k, a_, b_, tau_ = [Fr(int(rs.integers(-40, 41)), int(rs.integers(1, 30))) for _ in range(4)]
    if k == 0 or 1 + k * tau_ == 0:
        continue
    ok_alg &= (1 + k * a_) * (1 + k * b_) == 1 + k * fuse(a_, b_, k)            # the shift turns o_k into x
    run = Fr(0)
    for n_ in range(1, 7):
        run = fuse(run, tau_, k); ok_alg &= 1 + k * run == (1 + k * tau_) ** n_   # n steps -> (1 + k tau)^n
    inv_ = -tau_ / (1 + k * tau_); ok_alg &= fuse(tau_, inv_, k) == 0              # the step's inverse
    ok_alg &= fuse(a_, -1 / k, k) == -1 / k and 1 + k * (-1 / k) == 0              # the centre absorbs: scale 0
    r_ = 1 + k * tau_
    if k > 0 and tau_ > 0:                                                         # a pure growing run: total conformal time
        ok_alg &= tau_ * r_ / (r_ - 1) == 1 / k + tau_                             # a0 sum tau / a_n = 1/k + tau, exactly
check("S1 the fused shift, exactly (fractions): it turns o_k into multiplication; n time steps tau compose to a size ratio (1 + k tau)^n; a step has an inverse; the centre -1/k is size zero and absorbs; and a pure growing run's horizon (light's total reach, a0 x sum of tau / a_n) is exactly 1/k + tau",
      ok_alg, "200 random rational k, a, b, tau")

def shift_read(t, A, m=None):
    """read a run of slice sizes through the shift: a_n = F (1 + k tau)^n + G (1 + k tau)^-n, the forward run of the step
    plus the run of its inverse. Returns the rate kappa = ln(1 + k tau) / tau, where the two runs balance (the bounce),
    the size there, and how far the sizes are from that form."""
    tau_ = t[1] - t[0]; m = len(A) if m is None else m
    Cn = (A[2:m] + A[:m - 2]) / (2 * A[1:m - 1])                     # the same for every n exactly when a_n has that form
    kap = math.acosh(Cn.mean()) / tau_; R_ = math.exp(kap * tau_); nn = np.arange(m)
    Mx = np.stack([R_**nn, R_**-nn], 1); (F, G), *_ = np.linalg.lstsq(Mx, A[:m], rcond=None)
    return dict(kappa=kap, bounce_t=t[0] + math.log(G / F) / (2 * kap), throat=2 * math.sqrt(F * G),
                F=float(F), G=float(G), t0=float(t[0]), tau=float(tau_), C_spread=float(np.ptp(Cn) / (Cn.mean() - 1)), fit=float(np.abs(Mx @ [F, G] - A[:m]).max() / A[:m].max()))
reads = {}
for L in (1, 2, 3):
    A = grow_universe(SL[L], math.cosh(0.5), 0.025, 10)
    reads[L] = shift_read(0.5 + 0.025 * np.arange(len(A)), A); reads[L]["a"] = list(map(float, A))
r3 = reads[3]
check("S2 read through the fused shift, the glued slabs' sizes (started at t = 0.5, 10 steps of 0.025) are one time step's forward run plus its inverse's run: the rate is sqrt(Lambda/3), the two runs balance at t = 0 (the bounce, never seen by the run), at the throat size 1/k = sqrt(3/Lambda), all converging with refinement",
      r3["fit"] < 1e-5 and abs(r3["kappa"] - 1) < 1e-3 and abs(r3["bounce_t"]) < 0.01 and abs(r3["throat"] - throat[3]) < 2e-4
      and abs(reads[1]["kappa"] - 1) > abs(reads[2]["kappa"] - 1) > abs(r3["kappa"] - 1),
      "; ".join(f"level {L}: rate {r['kappa']:.5f}, bounce at t = {r['bounce_t']:+.4f}, throat {r['throat']:.5f} (static throat {throat[L]:.5f}), off that form by {r['fit']:.0e}" for L, r in reads.items()))
results["shift_read"] = reads

edge = [float(np.mean(np.linalg.norm(to_s3(Cl)[Tl[:, 0]] - to_s3(Cl)[Tl[:, 1]], axis=1))) for Cl, Tl in LEVELS]
drift = {}
for L, tau_, steps in ((1, 0.1, 40), (2, 0.05, 80), (3, 0.05, 80)):
    A = grow_universe(SL[L], math.cosh(0.5), tau_, steps)
    Cn = (A[2:] + A[:-2]) / (2 * A[1:-1])
    drift[L] = dict(piece=list(map(float, A[1:-1] * edge[L])), rate=list(map(float, np.arccosh(Cn) / tau_)), t_end=0.5 + tau_ * (len(A) - 1))
def rate_at(L, p):
    pc, rt = np.array(drift[L]["piece"]), np.array(drift[L]["rate"]); return float(np.interp(p, pc, rt))
p_small = [max([abs(r - 1) for r, p in zip(drift[L]["rate"], drift[L]["piece"]) if p < 0.25]) for L in (2, 3)]
last = [drift[L]["piece"][-1] for L in (1, 2, 3)]; t_end = [drift[L]["t_end"] for L in (1, 2, 3)]
check("S3 the rate read step by step stays within about 1% of sqrt(Lambda/3) while the pieces are under a quarter of 1/k (the size of the absorbing centre: de Sitter's length) and drifts up as they near it; at every level the run ends, with no growing solution, once pieces pass about 1/k, later in time for finer pieces. Round 14's breakdown is the slabs leaving the shift",
      max(p_small) < 0.012 and all(1.0 < p < 2.5 for p in last) and t_end[0] < t_end[1] < t_end[2],
      f"largest drift while pieces < 0.25 x (1/k): level 2 {p_small[0] * 100:.2f}%, level 3 {p_small[1] * 100:.2f}%; rate when pieces = 1/k: level 1 {rate_at(1, 1.0):.3f}, level 2 {rate_at(2, 1.0):.3f}, level 3 {rate_at(3, 1.0):.3f};"
      f" runs end at pieces {last[0]:.2f}, {last[1]:.2f}, {last[2]:.2f} x (1/k), at t = {t_end[0]:.2f}, {t_end[1]:.2f}, {t_end[2]:.2f}")
results["rate_drift"] = drift; results["edge_over_a"] = edge

def two_runs(k, t, tau_):
    z = 1 + k * tau_; nn = round(t / tau_); return 0.5 * (z**nn + z**-nn)
gl = [abs(two_runs(1.0, 1.3, tt) - math.cosh(1.3)) / math.cosh(1.3) for tt in (1e-2, 1e-3, 1e-4)]
ge = [abs(two_runs(1j, 1.3, tt) - math.cos(1.3)) for tt in (1e-2, 1e-3, 1e-4)]
check("S4 what the shift does not settle: with k real the forward and inverse runs give cosh (Lorentzian de Sitter); with k imaginary the same two runs give cos, the Euclidean 4-sphere. So the signature is the choice of k real, still an input; and any time step works, only k sets the rate",
      gl[-1] < 1e-4 and abs(ge[-1]) < 1e-4 and gl[0] > 10 * gl[-1],
      f"k = 1 vs cosh 1.3: {gl[0]:.1e}, {gl[1]:.1e}, {gl[2]:.1e} for tau = 1e-2, 1e-3, 1e-4; k = i vs cos 1.3: {abs(ge[0]):.1e}, {abs(ge[1]):.1e}, {abs(ge[2]):.1e}")

results["summary"] = {"passed": passed, "failed": failed, "seconds": round(time.time() - T_START)}
say(f"\n{passed} passed, {failed} failed  ({time.time() - T_START:.0f} s)")
json.dump(results, open("results.json", "w"), indent=1, default=float)
open("checks_output.txt", "w").write(_log.getvalue())
