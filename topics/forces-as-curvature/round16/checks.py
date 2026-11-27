"""DRAFT round 16 checks: does the law's dial q set the fused shift's k?

Round 15 read the glued slabs' growth through the fused shift a -> 1 + k a and found that at the bounce the space
curvature is q = k^2, so Lambda = 3 k^2. Here the question is whether the law's curvature dial q fixes k (and so
Lambda, de Sitter's length, the bounce size), whether it fixes the time step, and whether it fixes the sign of time.
Part Q only; round 14's slab machinery and round 15's shift reading are reused unchanged.
Run: python3 checks.py  (numpy, scipy, sympy).  Writes results.json and checks_output.txt.
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


# ------------------------------------------------------------------ round 15's reading through the fused shift (unchanged)
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
# ================================================================== Q. Does the law's dial q set the fused shift's k?
from fractions import Fraction as Fr
import sympy as sp
say("\n== Q. Does the law's dial q set the fused shift's k (and so Lambda = 3 k^2)?")

# Q1 the dial is a scale: the law at q is the law at 1 with every velocity scaled by sqrt(q)
def oplus_fr(a, b, qq):
    ab = sum(x * y for x, y in zip(a, b)); aa = sum(x * x for x in a); bb = sum(x * x for x in b)
    den = 1 - 2 * qq * ab + qq**2 * aa * bb
    return [((1 - 2 * qq * ab - qq * bb) * x + (1 + qq * aa) * y) / den for x, y in zip(a, b)]
rq = np.random.default_rng(16); ok_scale = True
for _ in range(200):
    s_ = Fr(int(rq.integers(1, 40)), int(rq.integers(1, 40)))
    a_ = [Fr(int(rq.integers(-30, 31)), int(rq.integers(1, 20))) for _ in range(3)]
    b_ = [Fr(int(rq.integers(-30, 31)), int(rq.integers(1, 20))) for _ in range(3)]
    lhs = oplus_fr(a_, b_, s_ * s_); rhs = [x / s_ for x in oplus_fr([s_ * x for x in a_], [s_ * x for x in b_], Fr(1))]
    ok_scale &= lhs == rhs
def growq(Kk, qq, generations=8, cap=3000):
    """round 13's growth at any dial: apply each kick and its reverse with the law at q, keep new points"""
    Kk = np.concatenate([Kk, -Kk]); front = [np.zeros(3)]; counts = [1]; P = [np.zeros(3)]; far = False
    pos = (lambda y: to_s3(math.sqrt(qq) * y)) if qq > 0 else (lambda y: y)
    S = [pos(np.zeros(3))]
    with np.errstate(all="ignore"):
        for _ in range(generations):
            new = []
            for p in front:
                Dn = 1 - 2 * qq * (Kk @ p) + qq**2 * (Kk * Kk).sum(1) * (p @ p)
                for y, dn in zip(oplus(Kk, np.tile(p, (len(Kk), 1)), qq), Dn):
                    if dn < 1e-9:
                        far = True; continue
                    sy = pos(y)
                    if np.min(np.linalg.norm(np.array(S) - sy, axis=1)) > 1e-7 * (1 + np.linalg.norm(sy)):
                        S.append(sy); P.append(y); new.append(y)
                if len(S) > cap:
                    break
            counts.append(len(S) + (1 if (far and qq > 0) else 0)); front = new
            if not new or len(S) > cap:
                break
    return counts, np.array(P)
kick = ico * math.tan(math.radians(72) / 2)
c1, P1 = growq(kick, 1.0)
dev = {}
for qq in (0.25, 4.0):
    cq, Pq = growq(kick / math.sqrt(qq), qq)
    # same points, in the same order, scaled by 1/sqrt(q)
    dev[qq] = float(np.abs(Pq * math.sqrt(qq) - P1).max()) if len(Pq) == len(P1) else float("inf")
check("Q1 (exact) the dial is a scale: the law at q is the law at q = 1 with every velocity multiplied by sqrt(q). So the corners grown at q (kicks of size tan 36 deg / sqrt q) are round 13's corners shrunk by sqrt(q): a 3-sphere of radius 1/sqrt(q), whose curvature Regge reads back as q (round 12)",
      ok_scale and max(dev.values()) < 1e-12 and c1[-1] == 120,
      f"200 random rational pairs and dials q = s^2, exact; grown corners at q = 1/4 and 4 against round 13's, scaled: {dev[0.25]:.0e}, {dev[4.0]:.0e}; {c1[-1]} corners")
results["dial_scale"] = {"exact_pairs": 200, "corner_dev": dev}

# Q2 with the kick held fixed, only one dial closes space
closure = {}
for qq in (1.0, 0.9, 1.1, 0.25, 4.0, 0.0, -1.0):
    closure[qq] = growq(kick, qq)[0]
closes = {qq: (cnt[-1] == 120) for qq, cnt in closure.items()}
check("Q2 hold the kick fixed and turn the dial: the grown corners close into a finite space (the 600-cell's 120) only at the dial where sqrt(q) x kick = tan 36 deg. Slightly off it (q = 0.9, 1.1), far off it, flat (q = 0) or open (q = -1), they never close. So closure fixes q in units of the kick, and there is no closed space at q <= 0",
      closes[1.0] and not any(closes[qq] for qq in closure if qq != 1.0) and all(closure[qq][-1] > 1000 for qq in closure if qq != 1.0),
      "; ".join(f"q = {qq:g}: {' -> '.join(map(str, cnt[:6]))}{' (closed)' if closes[qq] else ' ...'}" for qq, cnt in closure.items()))
results["closure"] = {str(k_): v for k_, v in closure.items()}

edge0 = [float(np.mean(np.linalg.norm(to_s3(Cl)[Tl[:, 0]] - to_s3(Cl)[Tl[:, 1]], axis=1))) for Cl, Tl in LEVELS]   # mean chord edge, unit sphere
# Q3 the moment of rest: Regge's constraint on the law-sized slice fixes Lambda
QS = (0.25, 1.0, 4.0)
lam_star = {qq: [] for qq in QS}
for qq in QS:
    a0 = 1 / math.sqrt(qq)
    for L in range(4):
        tl, _, V = SL[L].parts(a0, a0, 0.05 * a0)
        lam_star[qq].append(tl / V)                    # a static slab balances where tl = Lambda V
ratios = [[lam_star[qq][L] / qq for L in range(4)] for qq in QS]
spread_q = max(abs(ratios[i][L] - ratios[1][L]) for i in range(3) for L in range(4))
gaps = [ratios[1][L] - 3 for L in range(4)]
check("Q3 a moment of rest: on the slice the law builds at dial q, Regge's constraint for a static slab holds only for one Lambda, and Lambda / q -> 3 under refinement, the same at every dial. This is the textbook constraint for time-symmetric data, 6/a^2 = 2 Lambda, on law-sized pieces: Lambda = 3q, so k = sqrt(q)",
      spread_q < 1e-9 and all(g > 0 for g in gaps) and gaps[0] > gaps[1] > gaps[2] > gaps[3] and gaps[3] < 0.01,
      f"Lambda/q at levels 0-3: {', '.join(f'{r:.4f}' for r in ratios[1])}; the same at q = 1/4, 1, 4 to {spread_q:.0e}; gap to 3 falls by {gaps[1]/gaps[2]:.1f}, {gaps[2]/gaps[3]:.1f} per level")
results["lambda_star"] = {str(k_): v for k_, v in lam_star.items()}

# Q4 grow spacetime from that slice, from rest, and read k through the fused shift
def grow_from(S, a0, tau_, steps, lam, from_rest=True):
    A = [a0]
    for k_ in range(steps):
        a = A[-1]; h3 = lam / 3
        lo = a + (0.3 * tau_**2 * h3 * a if (from_rest and k_ == 0) else 1e-12 * a)
        hi = a + 3 * tau_ * math.sqrt(max(h3 * a * a - 1, 0)) + 3 * tau_**2 * h3 * a
        f = lambda x: (S.action(a, x, tau_ * (1 + 1e-5), lam) - S.action(a, x, tau_ * (1 - 1e-5), lam)) / (2e-5 * tau_)
        if f(lo) * f(hi) > 0:
            break
        A.append(brentq(f, lo, hi, xtol=1e-13 * a))
    return np.array(A)
def read_after_rest(A, tau_, skip=2):
    """read a run that starts at rest through the shift, leaving out its first steps (the step off rest is the one place
    round 14 found a first-order error, N4)"""
    t = tau_ * np.arange(len(A)); return shift_read(t[skip:], A[skip:])
def read_sym(A, tau_):
    """the same run with the first steps kept: mirrored to the past (the inverse run) and read whole"""
    full = np.concatenate([A[:0:-1], A]); t = tau_ * (np.arange(len(full)) - (len(A) - 1))
    return shift_read(t, full)
runs = {}
for qq in QS:
    for L in (1, 2, 3):
        if qq != 1.0 and L != 2:
            continue
        tau_ = 0.025 / math.sqrt(qq)
        A = grow_from(SL[L], 1 / math.sqrt(qq), tau_, 24, lam_star[qq][L])
        r = read_after_rest(A, tau_); r["a"] = list(map(float, A)); r["q"] = qq; r["L"] = L
        r["k_over_sqrtq"] = r["kappa"] / math.sqrt(qq); r["k_over_lam"] = r["kappa"] / math.sqrt(lam_star[qq][L] / 3)
        r["size_sqrtq"] = r["throat"] * math.sqrt(qq); r["bounce_sqrtq"] = r["bounce_t"] * math.sqrt(qq)
        runs[f"q={qq:g},L={L}"] = r
rq1 = [runs[f"q=1,L={L}"] for L in (1, 2, 3)]
same_q = max(abs(runs[f"q={qq:g},L=2"]["k_over_sqrtq"] - runs["q=1,L=2"]["k_over_sqrtq"]) for qq in QS)
check("Q4 grow the glued slabs from that slice, at rest, with that Lambda, and read them through the fused shift: they are one step's forward run plus its inverse run, with k -> sqrt(q), balancing at size 1/sqrt(q) half a step before the law's slice (where the first step off rest is centred), the same at every dial. The dial has set the shift",
      all(r["fit"] < 1e-4 for r in rq1) and abs(rq1[-1]["k_over_sqrtq"] - 1) < 0.003 and abs(rq1[-1]["bounce_sqrtq"]) < 0.6 * 0.025 and abs(rq1[-1]["size_sqrtq"] - 1) < 0.003
      and abs(rq1[0]["k_over_sqrtq"] - 1) > abs(rq1[1]["k_over_sqrtq"] - 1) > abs(rq1[2]["k_over_sqrtq"] - 1) and same_q < 1e-6,
      "; ".join(f"level {r['L']}: k/sqrt(q) {r['k_over_sqrtq']:.5f} (k/sqrt(Lambda/3) {r['k_over_lam']:.5f}), bounce at t sqrt(q) = {r['bounce_sqrtq']:+.4f}, size x sqrt(q) {r['size_sqrtq']:.5f}, off the form by {r['fit']:.0e}" for r in rq1)
      + f"; at level 2 the same for q = 1/4, 1, 4 to {same_q:.0e}")
results["grown_runs"] = runs

# Q5 the time step: free, and nothing depends on it
tau_rows = []
tl2, _, V2 = SL[2].parts(1.0, 1.0, 0.05); lam2 = tl2 / V2; k2 = math.sqrt(lam2 / 3)
for tau_ in (0.05, 0.025, 0.0125):
    A = grow_from(SL[2], 1.0, tau_, int(0.6 / tau_), lam2)
    rr = read_after_rest(A, tau_)
    tau_rows.append(dict(tau=tau_, after_rest=rr["kappa"] / k2, whole=read_sym(A, tau_)["kappa"] / k2, bounce_over_step=rr["bounce_t"] / tau_,
                         tau_over_edge=tau_ / edge0[2]))
spread_tau = np.ptp([r["after_rest"] for r in tau_rows]); whole_gap = [abs(r["whole"] - 1) for r in tau_rows]
check("Q5 the time step is not tied to q: at level 2, steps of 0.05, 0.025 and 0.0125 (from about a quarter down to a fourteenth of a space edge) give the same k after the first step off rest, and read whole (first step included) they close in on it as the step shrinks; the reading balances half a step before the slice at every step. The step is a free choice that the rate does not depend on, so there is nothing for the dial to fix",
      spread_tau < 2e-3 and whole_gap[0] > whole_gap[1] > whole_gap[2] and all(-0.6 < r["bounce_over_step"] < -0.4 for r in tau_rows),
      "; ".join(f"step {r['tau']} ({r['tau_over_edge']:.2f} edge): k after rest x {r['after_rest']:.5f}, whole run x {r['whole']:.5f}, bounce {r['bounce_over_step']:+.3f} steps" for r in tau_rows))
results["time_step"] = tau_rows

# Q6 the sign of time: Friedmann with a moment of rest, both signatures, both signs of space curvature
t_s, Lm, qs = sp.symbols("t Lambda q_s", real=True)
rows6 = []; ok6 = True
for sig, qv in (("Lorentzian", -1), ("Euclidean", 1)):
    for qsv in (1, -1):
        a0 = 1                                                   # moment of rest at a = 1, space curvature q_s = kappa / a^2
        kap = qsv
        # Friedmann: Lorentzian adot^2 = Lambda a^2 / 3 - kappa ; Euclidean adot^2 = kappa - Lambda a^2 / 3
        lam_rest = sp.solve(sp.Eq(0, (Lm * a0**2 / 3 - kap) * (1 if qv == -1 else -1)), Lm)[0]
        acc = (lam_rest / 3) * (1 if qv == -1 else -1)          # a'' = +-(Lambda/3) a
        ksq = sp.nsimplify(acc)                                  # the shift's k^2: a'' = k^2 a
        sol = sp.cosh(sp.sqrt(ksq) * t_s) if ksq > 0 else sp.cos(sp.sqrt(-ksq) * t_s)
        ok6 &= sp.simplify(sp.diff(sol, t_s, 2) - ksq * sol) == 0 and lam_rest == 3 * qsv and ksq == -qv * qsv
        rows6.append(dict(signature=sig, q_v=qv, q_s=qsv, Lambda=int(lam_rest), k2=int(ksq), growth="cosh" if ksq > 0 else "cos",
                          law_closes_space=qsv > 0))
check("Q6 the sign of time is not set by q: at a moment of rest Lambda = 3 q_s in either signature, but the slices then grow as cosh or cos depending on the signature too: k^2 = -q_v q_s, with q_v = -1 for boosts (relativity's end of the dial) and +1 for rotations. The law closes space only at q_s > 0 (Q2), so Lambda > 0; with relativity's boosts that is expansion (de Sitter), with rotations a 4-sphere",
      ok6, "; ".join(f"{r['signature']}, q_s = {r['q_s']:+d}: Lambda = {r['Lambda']:+d}, k^2 = {r['k2']:+d} ({r['growth']}){'' if r['law_closes_space'] else ', no closed law space'}" for r in rows6))
results["signature"] = rows6

# Q7 what stays an input: Lambda in units of the law's own pieces
per_piece = [lam_star[1.0][L] * edge0[L] ** 2 for L in range(4)]
fall = [per_piece[L] / per_piece[L + 1] for L in range(3)]
LAM_OBS_LP2 = 2.9e-122                                             # observed Lambda x Planck length^2 (Planck 2018), approximate
halvings = math.log(per_piece[0] / LAM_OBS_LP2) / math.log(4)
check("Q7 what q does not settle: Lambda in units of the law's own pieces is about 1 on the 600-cell and falls by about 4 with every refinement (each halving of the edges; the first refinement mixes edge lengths), so q sets k only relative to the pieces. The number of refinements, or of pieces, is still an input: if the pieces were Planck-sized, matching today's Lambda would take about 200 halvings",
      per_piece[0] > 0.5 and all(2.5 < f_ < 4.5 for f_ in fall) and all(abs(f_ - 4) < 0.45 for f_ in fall[1:]),
      f"Lambda x edge^2 at levels 0-3: {', '.join(f'{x:.4f}' for x in per_piece)} (falls by {', '.join(f'{f_:.2f}' for f_ in fall)}); halvings to Lambda l_P^2 = 2.9e-122: {halvings:.0f}")
results["per_piece"] = {"lambda_edge2": per_piece, "halvings_to_observed": halvings}

results["summary"] = {"passed": passed, "failed": failed, "seconds": round(time.time() - T_START)}
say(f"\n{passed} passed, {failed} failed  ({time.time() - T_START:.0f} s)")
json.dump(results, open("results.json", "w"), indent=1, default=float)
open("checks_output.txt", "w").write(_log.getvalue())
