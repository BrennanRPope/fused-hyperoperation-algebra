"""DRAFT helpers: ctypes wrapper for libtb.so and three-body geometry (Jacobi, Hopf, shape sphere)."""
import ctypes, os
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
_SO = os.path.join(HERE, "libtb.so")
if not os.path.exists(_SO):  # the compiled file isn't kept in the repo; build it from tb.c (README's command)
    import subprocess
    subprocess.run(["gcc", "-O3", "-fopenmp", "-shared", "-fPIC", "-o", _SO, os.path.join(HERE, "tb.c"), "-lm"], check=True)
_lib = ctypes.CDLL(_SO)
_d = ctypes.POINTER(ctypes.c_double)
_lib.integrate.argtypes = [_d, ctypes.c_double, ctypes.c_double]; _lib.integrate.restype = ctypes.c_int
_lib.sample.argtypes = [_d, _d, ctypes.c_int, _d, ctypes.c_double]; _lib.sample.restype = ctypes.c_int
_lib.landscape.argtypes = [_d, _d, ctypes.c_int, ctypes.c_double, ctypes.c_double, ctypes.c_double] + [_d]*5
P = lambda a: a.ctypes.data_as(_d)

def flow(X0, T, tol=1e-13):
    X = np.ascontiguousarray(X0, float).copy()
    if _lib.integrate(P(X), float(T), tol) < 0: raise RuntimeError("integration failed")
    return X

def traj(X0, T, n=2000, tol=1e-13):
    X = np.ascontiguousarray(X0, float).copy()
    ts = np.linspace(0, T, n)[1:]
    out = np.zeros((len(ts), 12))
    if _lib.sample(P(X), P(np.ascontiguousarray(ts)), len(ts), P(out), tol) < 0: raise RuntimeError("fail")
    return np.concatenate([[0.0], ts]), np.vstack([np.asarray(X0, float), out])

def landscape(p1, p2, Tmax, tol=1e-10, leave=0.5):
    p1 = np.ascontiguousarray(p1, float).ravel(); p2 = np.ascontiguousarray(p2, float).ravel()
    n = len(p1); outs = [np.zeros(n) for _ in range(5)]
    _lib.landscape(P(p1), P(p2), n, Tmax, tol, leave, *[P(o) for o in outs])
    return outs  # dabs, tabs, drot, trot, arot

def suvakov(p1, p2):
    return np.array([-1, 0, 1, 0, 0, 0, p1, p2, p1, p2, -2 * p1, -2 * p2], float)

def energy(X):
    r = X[:6].reshape(3, 2); v = X[6:].reshape(3, 2)
    K = 0.5 * (v ** 2).sum()
    U = sum(1 / np.linalg.norm(r[i] - r[j]) for i in range(3) for j in range(i + 1, 3))
    return K - U

def angmom(X):
    r = X[:6].reshape(3, 2); v = X[6:].reshape(3, 2)
    return float((r[:, 0] * v[:, 1] - r[:, 1] * v[:, 0]).sum())

def jacobi(Xs):
    """mass-weighted Jacobi vectors as complex numbers (equal masses): kinetic metric is Euclidean on C^2."""
    Xs = np.atleast_2d(Xs)
    r = Xs[:, :6:2] + 1j * Xs[:, 1:6:2]
    z1 = (r[:, 1] - r[:, 0]) / np.sqrt(2)
    z2 = np.sqrt(2 / 3) * (r[:, 2] - 0.5 * (r[:, 0] + r[:, 1]))
    return z1, z2

def shape(Xs):
    """Hopf map C^2 -> unit shape sphere (Montgomery): w = (|z1|^2-|z2|^2, 2 conj(z1) z2)/|z|^2."""
    z1, z2 = jacobi(Xs)
    R2 = abs(z1) ** 2 + abs(z2) ** 2
    w = np.conj(z1) * z2
    return np.stack([(abs(z1) ** 2 - abs(z2) ** 2) / R2, 2 * w.real / R2, 2 * w.imag / R2], axis=1)

def s3point(Xs):
    """normalized configuration on the unit 3-sphere in C^2 = R^4."""
    z1, z2 = jacobi(Xs)
    R = np.sqrt(abs(z1) ** 2 + abs(z2) ** 2)
    return np.stack([z1.real / R, z1.imag / R, z2.real / R, z2.imag / R], axis=1)

# collision points on the shape sphere (binary collisions) and Euler/Lagrange points
def special_points():
    import itertools
    pts = {}
    def sp(r):
        X = np.zeros(12); X[:6] = np.ravel(r); return shape(X)[0]
    e = 1e-7
    pts["collision 12"] = sp([[0, 0], [e, 0], [1, 0]]) ; pts["collision 23"] = sp([[1, 0], [0, 0], [e, 0]]); pts["collision 13"] = sp([[0, 0], [1, 0], [e, 0]])
    pts["Euler 3 mid"] = sp([[-1, 0], [1, 0], [0, 0]]); pts["Euler 1 mid"] = sp([[0, 0], [1, 0], [-1, 0]]); pts["Euler 2 mid"] = sp([[1, 0], [0, 0], [-1, 0]])
    s = np.sqrt(3) / 2
    pts["Lagrange +"] = sp([[-0.5, 0], [0.5, 0], [0, s]]); pts["Lagrange −"] = sp([[-0.5, 0], [0.5, 0], [0, -s]])
    return pts

def refine(p1, p2, T, rotating=False, alpha=0.0, tol=1e-13):
    """Shoot to an exact periodic orbit in the Suvakov slice: solve X(T) = X0 (or = R_alpha X0)."""
    from scipy.optimize import least_squares
    def rot(X, a):
        c, s = np.cos(a), np.sin(a); Y = X.copy()
        for i in range(6):
            x, y = X[2 * i], X[2 * i + 1]; Y[2 * i] = c * x - s * y; Y[2 * i + 1] = s * x + c * y
        return Y
    def res(q):
        X0 = suvakov(q[0], q[1])
        try:
            X = flow(X0, q[2], tol)
        except RuntimeError:
            return np.full(12, 1e3)
        return X - (rot(X0, q[3]) if rotating else X0)
    x0 = [p1, p2, T] + ([alpha] if rotating else [])
    sol = least_squares(res, x0, method="lm", xtol=1e-15, ftol=1e-15, gtol=1e-15, max_nfev=400)
    return sol.x, float(np.linalg.norm(sol.fun))

# ---------------------------------------------------------------- topology on the shape sphere
def syzygies(t, Xs):
    """Equator crossings of the shape sphere (collinear moments): returns list of (time, middle body 1..3, up?)."""
    w3 = shape(Xs)[:, 2]
    out = []
    for i in np.nonzero(np.sign(w3[1:]) != np.sign(w3[:-1]))[0]:
        if w3[i] == 0 and i > 0: continue
        f = w3[i] / (w3[i] - w3[i + 1]); X = Xs[i] + f * (Xs[i + 1] - Xs[i])
        r = X[:6].reshape(3, 2); d = r[1] - r[0] if np.linalg.norm(r[1] - r[0]) > 1e-12 else r[2] - r[0]
        s = r @ (d / np.linalg.norm(d)); mid = int(np.argsort(s)[1]) + 1
        out.append((t[i] + f * (t[i + 1] - t[i]), mid, bool(w3[i + 1] > w3[i])))
    return out

def free_reduce(w):
    st = []
    for c in w:
        if st and st[-1] == c.swapcase(): st.pop()
        else: st.append(c)
    s = "".join(st)
    while len(s) > 1 and s[0] == s[-1].swapcase(): s = s[1:-1]
    return s

def word_from_syzygies(syz, ignore=3, letters=("a", "b")):
    others = [k for k in (1, 2, 3) if k != ignore]
    m = {others[0]: letters[0], others[1]: letters[1]}
    w = "".join(m[b] if up else m[b].upper() for _, b, up in syz if b != ignore)
    return free_reduce(w)

def canon(w):
    def inv(x): return x[::-1].swapcase()
    def ab(x): return x.translate(str.maketrans("abAB", "baBA"))
    cands = set()
    for v in (w, inv(w), w.swapcase(), inv(w.swapcase())):
        for u in (v, ab(v)):
            for k in range(max(1, len(u))): cands.add(u[k:] + u[:k])
    return min(cands, key=lambda s: (len(s), s)) if cands else ""
