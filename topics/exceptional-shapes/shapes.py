"""Shared constructions for the exceptional-shapes pass (24-cell, 600-cell, 120-cell, E8)."""
import itertools
import numpy as np

PHI = (1 + 5 ** 0.5) / 2


def qmul(a, b):
    a = np.asarray(a, float); b = np.asarray(b, float)
    w1, x1, y1, z1 = np.moveaxis(a, -1, 0)
    w2, x2, y2, z2 = np.moveaxis(b, -1, 0)
    return np.stack([w1*w2 - x1*x2 - y1*y2 - z1*z2,
                     w1*x2 + x1*w2 + y1*z2 - z1*y2,
                     w1*y2 - x1*z2 + y1*w2 + z1*x2,
                     w1*z2 + x1*y2 - y1*x2 + z1*w2], axis=-1)


def qconj(a):
    a = np.array(a, float); a[..., 1:] *= -1; return a


def key(v, nd=6):
    return tuple(np.round(np.asarray(v, float) + 0.0, nd) + 0.0)


def keyset(V):
    return {key(v) for v in V}


def uniq(V, nd=6):
    out, seen = [], set()
    for v in V:
        k = key(v, nd)
        if k not in seen:
            seen.add(k); out.append(np.asarray(v, float))
    return np.array(out)


def cross_polytope4():
    V = []
    for i in range(4):
        for s in (1, -1):
            v = np.zeros(4); v[i] = s; V.append(v)
    return np.array(V)


def tesseract_half():
    return np.array([np.array(s) / 2 for s in itertools.product((1, -1), repeat=4)])


def binary_tetrahedral():   # 24 Hurwitz units = 24-cell
    return np.vstack([cross_polytope4(), tesseract_half()])


def dual24():                # 24 points (two of +-1, +-i, +-j, +-k)/sqrt2
    V = []
    for i, j in itertools.combinations(range(4), 2):
        for s, t in itertools.product((1, -1), repeat=2):
            v = np.zeros(4); v[i] = s; v[j] = t; V.append(v / 2 ** 0.5)
    return np.array(V)


def binary_octahedral():     # 48 = Clifford group (mod phase), 24-cell + dual 24-cell
    return np.vstack([binary_tetrahedral(), dual24()])


def even_perms(n=4):
    out = []
    for p in itertools.permutations(range(n)):
        inv = sum(1 for i in range(n) for j in range(i + 1, n) if p[i] > p[j])
        if inv % 2 == 0:
            out.append(p)
    return out


def binary_icosahedral():    # 120 = 600-cell vertices
    V = list(binary_tetrahedral())
    base = np.array([PHI, 1, 1 / PHI, 0]) / 2
    for p in even_perms():
        for s in itertools.product((1, -1), repeat=3):
            v = np.zeros(4)
            vals = base * np.array(list(s) + [1])
            for pos, src in enumerate(p):
                v[pos] = vals[src]
            V.append(v)
    return uniq(V)


def simplex4():
    # regular 4-simplex centred at 0, circumradius 1
    E = np.eye(5) - 1 / 5
    U, S, Vt = np.linalg.svd(E)
    P = E @ Vt[:4].T
    return P / np.linalg.norm(P[0])


def hull_facets(V):
    from scipy.spatial import ConvexHull
    h = ConvexHull(V)
    eq = h.equations  # normal . x + off <= 0
    normals = uniq(np.hstack([eq[:, :4], eq[:, 4:]]), nd=5)
    return h, normals


def onetwenty_cell():
    """Dual of the 600-cell: unit facet normals of the 600-cell's hull."""
    _, F = hull_facets(binary_icosahedral())
    N = F[:, :4]
    return uniq(N / np.linalg.norm(N, axis=1, keepdims=True), nd=5)


def closure_fraction(G, op):
    S = keyset(G)
    hits = 0
    for a in G:
        prods = op(np.broadcast_to(a, G.shape), G)
        hits += sum(key(p) in S for p in prods)
    return hits / len(G) ** 2


def generated(gens, op, cap=5000):
    S = {key(g): np.asarray(g, float) for g in gens}
    frontier = list(S.values())
    while frontier:
        new = []
        for a in frontier:
            for b in list(S.values()):
                for p in (op(a, b), op(b, a)):
                    k = key(p)
                    if k not in S:
                        S[k] = p; new.append(p)
                        if len(S) > cap:
                            return np.array(list(S.values()))
        frontier = new
    return np.array(list(S.values()))


def commutator_subgroup(G):
    comms = [qmul(qmul(a, b), qmul(qconj(a), qconj(b))) for a in G for b in G]
    C = uniq(comms)
    return generated(C, qmul), C


def rot_angle_deg(q):
    w = np.clip(abs(q[0]), 0, 1)
    return np.degrees(2 * np.arccos(w))


# ---------- E8 ----------
def e8_roots():
    R = []
    for i, j in itertools.combinations(range(8), 2):
        for s, t in itertools.product((1, -1), repeat=2):
            v = np.zeros(8); v[i] = s; v[j] = t; R.append(v)
    for s in itertools.product((0.5, -0.5), repeat=8):
        if sum(x < 0 for x in s) % 2 == 0:
            R.append(np.array(s))
    return np.array(R)


def galois(v):
    """Apply sqrt5 -> -sqrt5 to coordinates that lie in Q(sqrt5) with half-integer parts (a + b sqrt5)/2."""
    v = np.asarray(v, float)
    out = np.empty_like(v)
    for idx, x in np.ndenumerate(v):
        # find a, b in (1/4)Z with x = a + b sqrt5
        best = None
        for b4 in range(-12, 13):
            b = b4 / 4
            a = x - b * 5 ** 0.5
            if abs(a * 4 - round(a * 4)) < 1e-7:
                cand = (abs(b4), round(a * 4) / 4, b)
                if best is None or cand < best:
                    best = cand
        out[idx] = best[1] - best[2] * 5 ** 0.5
    return out


def e8_from_icosians():
    """240 icosians of 'modified norm' 1 (2I and 2I/phi), mapped to R^8 as (alpha x, beta conj(x))."""
    G = binary_icosahedral()
    X = np.vstack([G, G / PHI])
    Xb = np.array([galois(x) for x in X])
    al = ((1 + 1 / 5 ** 0.5) / 2) ** 0.5
    be = ((1 - 1 / 5 ** 0.5) / 2) ** 0.5
    return np.hstack([al * X, be * Xb]), X, al, be


def simple_roots(R, rng_seed=1):
    from scipy.optimize import nnls
    rng = np.random.default_rng(rng_seed)
    v = rng.normal(size=R.shape[1])
    P = R[R @ v > 0]
    simple = []
    for i, r in enumerate(P):
        A = np.delete(P, i, axis=0).T
        _, res = nnls(A, r)
        if res > 1e-6:
            simple.append(r)
    return np.array(simple)


def coxeter_plane(R):
    S = simple_roots(R)
    n = R.shape[1]
    C = np.eye(n)
    for s in S:
        C = C @ (np.eye(n) - 2 * np.outer(s, s) / (s @ s))
    w, U = np.linalg.eig(C)
    ang = np.angle(w)
    i = np.argmin(np.where(ang > 1e-9, ang, 99))
    h = 2 * np.pi / ang[i]
    u = U[:, i]
    e1, e2 = np.real(u), np.imag(u)
    e1 /= np.linalg.norm(e1); e2 -= (e2 @ e1) * e1; e2 /= np.linalg.norm(e2)
    return np.vstack([e1, e2]), h, S


# ---------- octonions ----------
LINES = [tuple(sorted(((n) % 7, (n + 1) % 7, (n + 3) % 7))) for n in range(7)]
_ORIENT = [((n) % 7, (n + 1) % 7, (n + 3) % 7) for n in range(7)]


def _otable():
    # index 0 = real unit, 1..7 = e_0..e_6
    T = np.zeros((8, 8), int); Sg = np.zeros((8, 8))
    for a in range(8):
        T[0, a] = a; Sg[0, a] = 1; T[a, 0] = a; Sg[a, 0] = 1
    for a in range(1, 8):
        T[a, a] = 0; Sg[a, a] = -1
    for (i, j, k) in _ORIENT:
        for (p, q, r) in ((i, j, k), (j, k, i), (k, i, j)):
            T[p + 1, q + 1] = r + 1; Sg[p + 1, q + 1] = 1
            T[q + 1, p + 1] = r + 1; Sg[q + 1, p + 1] = -1
    return T, Sg


OT, OS = _otable()


def omul(a, b):
    a = np.asarray(a, float); b = np.asarray(b, float)
    out = np.zeros(np.broadcast(a, b).shape)
    for i in range(8):
        for j in range(8):
            out[..., OT[i, j]] += OS[i, j] * a[..., i] * b[..., j]
    return out
