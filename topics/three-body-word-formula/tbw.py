"""DRAFT ctypes wrapper for libtbw.so: syzygy letters from the collinear slice."""
import ctypes, os, numpy as np
H = os.path.dirname(os.path.abspath(__file__))
if not os.path.exists(os.path.join(H, "libtbw.so")):  # the compiled file is not kept in the repo; build it from tbw.c
    import subprocess
    subprocess.run(["gcc", "-O3", "-fopenmp", "-shared", "-fPIC", "-o", os.path.join(H, "libtbw.so"), os.path.join(H, "tbw.c"), "-lm"], check=True)
_l = ctypes.CDLL(os.path.join(H, "libtbw.so"))
_d = ctypes.POINTER(ctypes.c_double); _i = ctypes.POINTER(ctypes.c_int)
_l.letters.argtypes = [ctypes.c_double]*4 + [ctypes.c_int, _i, _d]; _l.letters.restype = ctypes.c_int
_l.letters_grid.argtypes = [_d, _d, ctypes.c_int, ctypes.c_double, ctypes.c_double, ctypes.c_int, _i, _i]
CH = {1: 'a', -1: 'A', 2: 'b', -2: 'B'}
def raw(p1, p2, T, tol=1e-12, maxev=4000):
    L = np.zeros(maxev, np.int32); t = np.zeros(maxev)
    n = _l.letters(p1, p2, T, tol, maxev, L.ctypes.data_as(_i), t.ctypes.data_as(_d))
    ok = n >= 0; n = n if ok else -n - 1
    return ''.join(CH[x] for x in L[:n]), t[:n], ok
def grid(P1, P2, T, maxev, tol=1e-10):
    P1 = np.ascontiguousarray(P1, float).ravel(); P2 = np.ascontiguousarray(P2, float).ravel()
    L = np.zeros((len(P1), maxev), np.int32); c = np.zeros(len(P1), np.int32)
    _l.letters_grid(P1.ctypes.data_as(_d), P2.ctypes.data_as(_d), len(P1), T, tol, maxev, L.ctypes.data_as(_i), c.ctypes.data_as(_i))
    return L, c
def reduce(w):
    st = []
    for ch in w:
        if st and st[-1] == ch.swapcase(): st.pop()
        else: st.append(ch)
    return ''.join(st)
def creduce(w):
    s = reduce(w)
    while len(s) > 1 and s[0] == s[-1].swapcase(): s = s[1:-1]
    return s
_l.integrate.argtypes = [_d, ctypes.c_double, ctypes.c_double]; _l.integrate.restype = ctypes.c_int
def ic(p1, p2): return np.array([-1, 0, 1, 0, 0, 0, p1, p2, p1, p2, -2 * p1, -2 * p2], float)
def flow(X0, T, tol=1e-13):
    X = np.ascontiguousarray(X0, float).copy()
    if T > 0 and _l.integrate(X.ctypes.data_as(_d), float(T), tol) < 0: raise RuntimeError('integration failed')
    return X
def polish(p1, p2, T, tol=1e-13):
    """Light polishing step: Levenberg-Marquardt (Newton-like) shooting on X(T) = X(0) in the slice."""
    from scipy.optimize import least_squares
    def res(q):
        try: return flow(ic(q[0], q[1]), q[2], tol) - ic(q[0], q[1])
        except RuntimeError: return np.full(12, 1e3)
    s = least_squares(res, [p1, p2, T], method='lm', xtol=1e-15, ftol=1e-15, gtol=1e-15, max_nfev=200)
    return s.x, float(np.linalg.norm(s.fun))

_l.retdist_grid.argtypes = [_d, _d, _d, _d, ctypes.c_int, ctypes.c_double, _d, _d]
def retdist(P1, P2, T0, T1, tol=1e-11):
    A = [np.ascontiguousarray(x, float).ravel() for x in (P1, P2, T0, T1)]
    dm = np.zeros(len(A[0])); tm = np.zeros(len(A[0]))
    _l.retdist_grid(*[a.ctypes.data_as(_d) for a in A], len(A[0]), tol, dm.ctypes.data_as(_d), tm.ctypes.data_as(_d))
    return dm, tm
