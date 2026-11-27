"""DRAFT. Core tools for the three-body Conley-index test.

Planar equal-mass three-body problem (G = m = 1). Exact Jacobians from the variational
equations (libtbvar.so), the reduced 4x4 Floquet block, Floquet-pair types, and the
fixed-point index each type forces on every iterate."""
import ctypes, os, subprocess, json
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
_SO = os.path.join(HERE, "libtbvar.so")
if not os.path.exists(_SO):
    subprocess.run(["gcc", "-O3", "-shared", "-fPIC", "-o", _SO, os.path.join(HERE, "tbvar.c"), "-lm"], check=True)
_lib = ctypes.CDLL(_SO)
_d = ctypes.POINTER(ctypes.c_double)
for f in (_lib.integrate, _lib.integrate_var):
    f.argtypes = [_d, ctypes.c_double, ctypes.c_double]; f.restype = ctypes.c_int
_lib.next_hit.argtypes = [_d, ctypes.c_double, ctypes.c_double]; _lib.next_hit.restype = ctypes.c_double
TOL = 1e-13

def flow(X0, T, tol=TOL):
    X = np.ascontiguousarray(X0, float).copy()
    if T <= 0: return X
    if _lib.integrate(X.ctypes.data_as(_d), float(T), tol) < 0: raise RuntimeError("integration failed")
    return X

def flow_var(X0, T, tol=TOL, Phi0=None):
    XP = np.zeros(12 + 144)
    XP[:12] = X0; XP[12:] = (np.eye(12) if Phi0 is None else Phi0).ravel()
    if _lib.integrate_var(XP.ctypes.data_as(_d), float(T), tol) < 0: raise RuntimeError("integration failed")
    return XP[:12].copy(), XP[12:].reshape(12, 12).copy()

def suvakov(p1, p2):
    """Li-Liao / Suvakov start: bodies at (-1,0), (1,0), (0,0); v1 = v2 = (p1,p2), v3 = -2 v1."""
    return np.array([-1, 0, 1, 0, 0, 0, p1, p2, p1, p2, -2 * p1, -2 * p2], float)

def pos(X): return np.asarray(X)[:6].reshape(3, 2)
def vel(X): return np.asarray(X)[6:12].reshape(3, 2)

def energy(X):
    r, v = pos(X), vel(X)
    return 0.5 * (v ** 2).sum() - sum(1 / np.linalg.norm(r[i] - r[j]) for i in range(3) for j in range(i + 1, 3))

def angmom(X):
    r, v = pos(X), vel(X)
    return float((r[:, 0] * v[:, 1] - r[:, 1] * v[:, 0]).sum())

def area(X):
    """twice the oriented area of the triangle (1,2,3); zero exactly at a line-up."""
    r = pos(X); a = r[1] - r[0]; b = r[2] - r[0]
    return a[0] * b[1] - a[1] * b[0]

def field(X):
    """the vector field at X."""
    r = pos(X); acc = np.zeros((3, 2))
    for i in range(3):
        for j in range(i + 1, 3):
            d = r[j] - r[i]; f = d / np.linalg.norm(d) ** 3
            acc[i] += f; acc[j] -= f
    return np.concatenate([X[6:12], acc.ravel()])

def grad_conserved(X):
    """rows: d(COMx), d(COMy), d(Px), d(Py), dH, dL at X (12-vectors)."""
    r, v = pos(X), vel(X)
    G = np.zeros((6, 12))
    G[0, 0:6:2] = 1; G[1, 1:6:2] = 1; G[2, 6::2] = 1; G[3, 7::2] = 1
    f = field(X); G[4, :6] = -f[6:]; G[4, 6:] = X[6:12]          # dH = -accel.dr + v.dv
    G[5, 0:6:2] = v[:, 1]; G[5, 1:6:2] = -v[:, 0]; G[5, 6::2] = -r[:, 1]; G[5, 7::2] = r[:, 0]
    return G

def rot_gen(X):
    r, v = pos(X), vel(X)
    return np.concatenate([np.stack([-r[:, 1], r[:, 0]], 1).ravel(), np.stack([-v[:, 1], v[:, 0]], 1).ravel()])

def reduced_block(X0, Phi):
    """4x4 reduced Floquet block of the monodromy Phi at the periodic point X0:
    Phi on K = {dCOM = dP = dH = dL = 0} modulo span(flow, rotation)."""
    G = grad_conserved(X0)
    _, s, Vt = np.linalg.svd(G); K = Vt[6:].T                    # 12x6 orthonormal basis of K
    f = field(X0); xi = rot_gen(X0)
    A = np.column_stack([f, xi]); A = K @ np.linalg.lstsq(K, A, rcond=None)[0]
    Qa, _ = np.linalg.qr(A)
    Pperp = K @ K.T - Qa @ Qa.T                                  # projector onto K minus span(f, xi)
    w, U = np.linalg.eigh(Pperp); V = U[:, w > 0.5]               # 12x4
    B = np.column_stack([A, V])
    C = np.linalg.lstsq(B, Phi @ V, rcond=None)[0]
    return C[2:, :]

def pair_traces(M4):
    """tau_i = lambda_i + 1/lambda_i for the two Floquet pairs of a 4x4 symplectic block.
    char poly x^4 - a x^3 + b x^2 - a x + 1  =>  tau^2 - a tau + (b - 2) = 0."""
    c = np.poly(M4); a = -c[1]; b = c[2]
    disc = a * a - 4 * (b - 2)
    r = np.sqrt(complex(disc))
    return (a + r) / 2, (a - r) / 2, dict(palindrome_err=float(abs(c[1] - c[3]) + abs(c[4] - 1)), disc=float(disc))

def pair_type(tau):
    if abs(tau.imag) > 1e-9: return 'CU'
    t = tau.real
    return 'E' if abs(t) < 2 else ('H' if t > 2 else 'IH')

def orbit_class(t1, t2):
    if 'CU' in (pair_type(t1), pair_type(t2)): return 'CU'      # complex unstable (Krein quartet)
    return '-'.join(sorted([pair_type(t1), pair_type(t2)], key=['E', 'H', 'IH'].index))

def residue(tau):  # Greene's residue of one pair
    return (2 - tau) / 4

def index_from_traces(t1, t2, n=1):
    """fixed-point index of F^n at a nondegenerate fixed point: sign det(I - DF^n) = prod sign(2 - tau_i^(n)).
    tau^(n) = 2 T_n(tau/2) (Chebyshev); for a complex pair the factor is |2 - tau|^2 > 0."""
    def tn(t):
        lam = np.roots([1, -t, 1])[0]
        return lam ** n + lam ** (-n)
    a, b = tn(t1), tn(t2)
    return int(np.sign(((2 - a) * (2 - b)).real))

def next_hit_c(X0, tol=TOL, tmax=400.0):
    X = np.ascontiguousarray(X0, float).copy()
    t = _lib.next_hit(X.ctypes.data_as(_d), tol, tmax)
    if t < 0: raise RuntimeError("no hit")
    return X, t
