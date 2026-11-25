"""Checks for topics/why-i/dial-map (Nov 25, 2026). Run: python3 checks.py
Which of the framework's dials are deformations in the published (Inonu-Wigner / Faddeev / Mendes) sense?
Rigidity test: H^2(g, g) = 0 (Chevalley-Eilenberg, adjoint coefficients) means no infinitesimal deformation, so the
algebra is stable (Nijenhuis-Richardson / Gerstenhaber); H^2 != 0 means it can be deformed."""
import numpy as np, itertools
from scipy.linalg import expm
results = []
def check(name, ok, detail=""):
    results.append((name, bool(ok))); print(("PASS " if ok else "FAIL ") + name + (("  | " + detail) if detail else ""))

def algebra(n, brackets):
    """structure constants C[i, j, k]: [e_i, e_j] = sum_k C[i,j,k] e_k"""
    C = np.zeros((n, n, n))
    for (i, j), v in brackets.items():
        for kk, c in v.items(): C[i, j, kk] += c; C[j, i, kk] -= c
    return C
def jacobi_ok(C):
    n = C.shape[0]
    for i, j, k in itertools.combinations(range(n), 3):
        t = np.einsum('m,mlr->lr', C[j, k], C[i]) if False else None
    br = lambda x, y: np.einsum('i,j,ijk->k', x, y, C)
    E = np.eye(n)
    return all(np.allclose(br(E[i], br(E[j], E[k])) + br(E[j], br(E[k], E[i])) + br(E[k], br(E[i], E[j])), 0)
               for i, j, k in itertools.combinations(range(n), 3))
def H2_dim(C):
    n = C.shape[0]; E = np.eye(n); br = lambda x, y: np.einsum('i,j,ijk->k', x, y, C)
    pairs = list(itertools.combinations(range(n), 2)); trip = list(itertools.combinations(range(n), 3))
    pidx = {p: a for a, p in enumerate(pairs)}
    # C1 basis: phi_{a,b} maps e_a -> e_b ; C2 basis: psi_{(i,j),b} maps e_i ^ e_j -> e_b
    def psi_eval(psi, x, y):          # psi: (len(pairs), n) array, x, y vectors
        out = np.zeros(n)
        for (i, j), a in pidx.items():
            w = x[i]*y[j] - x[j]*y[i]
            if w: out += w*psi[a]
        return out
    d1 = np.zeros((len(pairs)*n, n*n))
    for col in range(n*n):
        a, b = divmod(col, n); phi = np.zeros((n, n)); phi[a, b] = 1          # phi(e_a) = e_b
        f = lambda x: x @ phi
        for (i, j), r in pidx.items():
            v = br(E[i], f(E[j])) - br(E[j], f(E[i])) - f(br(E[i], E[j]))
            d1[r*n:(r + 1)*n, col] = v
    d2 = np.zeros((len(trip)*n, len(pairs)*n))
    for col in range(len(pairs)*n):
        a, b = divmod(col, n); psi = np.zeros((len(pairs), n)); psi[a, b] = 1
        for t, (i, j, k) in enumerate(trip):
            x, y, z = E[i], E[j], E[k]
            v = (br(x, psi_eval(psi, y, z)) - br(y, psi_eval(psi, x, z)) + br(z, psi_eval(psi, x, y))
                 - psi_eval(psi, br(x, y), z) + psi_eval(psi, br(x, z), y) - psi_eval(psi, br(y, z), x))
            d2[t*n:(t + 1)*n, col] = v
    r1 = np.linalg.matrix_rank(d1, tol=1e-9) if d1.size else 0
    r2 = np.linalg.matrix_rank(d2, tol=1e-9) if d2.size else 0
    return len(pairs)*n - r2 - r1

# ---------- 1+1 kinematics: basis H (time shift), P (space shift), K (boost) ----------
def kin11(k1, k2):   # Cayley-Klein labels: kappa1 (curvature), kappa2 = -1/c^2 (signature)
    H, P, K = 0, 1, 2
    return algebra(3, {(K, H): {P: 1.0}, (K, P): {H: -k2}, (H, P): {K: -k1*k2}})
G11, P11, dS11, AdS11 = kin11(0, 0), kin11(0, -1), kin11(-1, -1), kin11(1, -1)
check("S0 all four 1+1 kinematical algebras (Galilei, Poincare, de Sitter, anti-de Sitter) satisfy the Jacobi identity",
      all(jacobi_ok(C) for C in [G11, P11, dS11, AdS11]))
h = [H2_dim(C) for C in [G11, P11, dS11, AdS11]]
check("S1 1+1: Galilei can be deformed (H^2 != 0), Poincare can still be deformed, de Sitter and anti-de Sitter are rigid (H^2 = 0): "
      "turning on 1/c^2 then the curvature reaches a stable algebra (Bacry-Levy-Leblond 1968; Mendes 1994)",
      h[0] > 0 and h[1] > 0 and h[2] == 0 and h[3] == 0, f"dim H^2 = {h}")

# ---------- 3+1 kinematics, built from the Cayley-Klein matrix realization (axes 0..4) ----------
def ck_struct(kaps):
    """labels kappa_ab = product of kaps[a:b]; generators M_ab = E_ba - kappa_ab E_ab; structure constants solved numerically"""
    N = len(kaps); n = N + 1; Ms = []
    for a in range(n):
        for b in range(a + 1, n):
            M = np.zeros((n, n)); M[b, a] = 1; M[a, b] = -float(np.prod(kaps[a:b])); Ms.append(M)
    B = np.array([M.ravel() for M in Ms]).T; d = len(Ms)
    assert np.linalg.matrix_rank(B) == d
    C = np.zeros((d, d, d))
    for i in range(d):
        for j in range(d):
            Cm = (Ms[i] @ Ms[j] - Ms[j] @ Ms[i]).ravel(); coef, res, *_ = np.linalg.lstsq(B, Cm, rcond=None)
            assert np.allclose(B @ coef, Cm); C[i, j] = coef
    return C
Gal, Poi, dS, AdS = (ck_struct(np.array(k)) for k in [(0, 0, 1, 1), (0, -1, 1, 1), (-1, -1, 1, 1), (1, -1, 1, 1)])
check("S2 the 3+1 algebras (Galilei, Poincare, de Sitter, anti-de Sitter), built as Cayley-Klein matrices with labels (kappa1, kappa2, 1, 1), "
      "are 10-dimensional Lie algebras (Jacobi holds)", all(jacobi_ok(C) for C in [Gal, Poi, dS, AdS]))
h3 = [H2_dim(C) for C in [Gal, Poi, dS, AdS]]
check("S3 3+1: Galilei can be deformed, Poincare can be deformed (toward de Sitter), de Sitter and anti-de Sitter are rigid (H^2 = 0)",
      h3[0] > 0 and h3[1] > 0 and h3[2] == 0 and h3[3] == 0, f"dim H^2 (Gal, Poi, dS, AdS) = {h3}")
Barg11 = algebra(4, {(2, 0): {1: 1.0}, (2, 1): {3: 1.0}})      # H 0, P 1, K 2, M 3: Galilei with the mass as a central charge
check("S3b 1+1 Bargmann (Galilei with mass, the Schrodinger particle's symmetry) is a Lie algebra and can still be deformed (toward Poincare)",
      jacobi_ok(Barg11) and H2_dim(Barg11) > 0, f"dim H^2 = {H2_dim(Barg11)}")

# ---------- the quantum axis, finite-dimensional stand-in ----------
# classical: x, p, 1 commute (abelian); quantum: [x, p] = hbar 1 (Heisenberg). Neither is rigid alone in dim 3;
# Mendes' point is that the relativistic + quantum algebra is still unstable and needs two more constants.
ab3 = algebra(3, {}); heis = algebra(3, {(0, 1): {2: 1.0}})
hq = [H2_dim(ab3), H2_dim(heis)]
check("S4 x, p, 1 commuting (classical) can be deformed into [x, p] = hbar (Heisenberg); the Heisenberg algebra itself is not rigid "
      "either, so the stable quantum statement lives in the full algebra of observables (Moyal: Vey 1975, Bayen et al. 1978), "
      "not in this 3-dimensional piece", hq[0] > 0 and hq[1] > 0, f"dim H^2 (abelian, Heisenberg) = {hq}")

# ---------- our scalar dials k and d are not deformations ----------
one = algebra(1, {})
check("S5 a one-parameter law has a 1-dimensional (abelian) Lie algebra, which has nothing to deform (H^2 = 0 trivially): "
      "so turning k (or d) on from 0 is not a deformation", H2_dim(one) == 0)
rng = np.random.default_rng(5); ok = True
for k in [-1.0, 0.5, 2.0]:
    L = lambda x: np.log1p(k*x)/k
    for _ in range(20):
        a, b = rng.uniform(-.3, .3, 2); ok &= np.isclose(L(a + b + k*a*b), L(a) + L(b))
check("S6 instead, every k is plain addition after the ruler log(1 + k x)/k (over the rationals every formal group law of one "
      "variable is isomorphic to addition): k changes the ruler and the global shape (walls, circle vs line), not the local law", ok)

# ---------- downhill limits are exact (contractions) ----------
ok = True
for v in [0.3, 1.1]:
    for cc in [1e1, 1e2, 1e3]:
        g = 1/np.sqrt(1 - v**2/cc**2); L = np.array([[g, -g*v/cc**2], [-g*v, g]])   # (t, x) Lorentz boost
        Gb = np.array([[1, 0], [-v, 1]])
        ok &= np.abs(L - Gb).max() < 2*(v + v**2)/cc**2 + 1e-12
check("S7 a Lorentz boost goes to a Galilean boost as c -> infinity, with error ~ 1/c^2 (Inonu-Wigner contraction, the cube's 1/c edge)", ok)
ok = True
for k1 in [0.5, -0.5]:
    for R in [1e1, 1e2]:
        C = kin11(k1/R**2, -1); ok &= np.abs(C - P11).max() <= abs(k1)/R**2 + 1e-15
check("S8 de Sitter / anti-de Sitter go to Poincare as the curvature radius grows (kappa1 = +-1/R^2 -> 0): the Cayley-Klein kappa1 edge", ok)

# ---------- the gravity edge: what a scalar chain can and cannot say ----------
import sympy as sp
a_, b_, c_, al, G_, r_ = sp.symbols("a b c alpha G r", real=True)
# G1. Newtonian binding to first order: two masses at distance r weigh mA + mB - G mA mB / (r c^2)  (units c = 1 below)
k_bind = -G_/r_
bind = a_ + b_ - G_*a_*b_/r_
check("G1 Newtonian binding, to first order in G, is exactly the weighted circle operation a + b + k ab with k = -G/(r c^2) "
      "(a ruler choice at fixed separation, ours as a reading)", sp.simplify(bind - (a_ + b_ + k_bind*a_*b_)) == 0)
# G2. the exact self-energy of a shell (Arnowitt-Deser-Misner 1960): bare m0 = m + alpha m^2, alpha = G/(2R).
L = lambda m: m + al*m**2
Linv = lambda y: (-1 + sp.sqrt(1 + 4*al*y))/(2*al)
ok_inv = sp.simplify(L(Linv(a_)) - a_) == 0
F = Linv(L(a_) + L(b_))
ser = sp.series(F.subs({a_: a_*sp.Symbol("t"), b_: b_*sp.Symbol("t")}), sp.Symbol("t"), 0, 4).removeO()
t = sp.Symbol("t")
o2 = sp.expand(ser.coeff(t, 1)*1 + ser.coeff(t, 2)); o3 = sp.expand(ser.coeff(t, 3))
check("G2 ADM dressed-mass law F = L^-1(L(a) + L(b)) with L(m) = m + (G/2R) m^2 agrees with a + b + k ab, k = -G/R, through "
      "second order and differs at third (extra 4 alpha^2 ab(a + b)); L^-1 is an exact closed form",
      ok_inv and sp.simplify(o2 - (a_ + b_ - 2*al*a_*b_)) == 0 and sp.simplify(o3 - 4*al**2*a_*b_*(a_ + b_)) == 0,
      f"third-order term = {sp.factor(o3)}")
# G3. toy self-coupling iteration: 'gravity weighs its own binding energy'. m_{n+1} = m0 - alpha m_n^2 converges to L^-1(m0).
ok = True
for alv in [0.05, 0.2]:
    for m0 in [0.3, 1.0]:
        m = m0
        for _ in range(200): m = m0 - alv*m*m
        ok &= abs(m - float(Linv(sp.Float(m0)).subs(al, alv))) < 1e-10
check("G3 toy bootstrap (ours, a scalar cartoon of Gupta-Feynman-Deser): iterating 'add the binding energy of the current mass' "
      "converges to the ADM closed form; the summed self-coupling is a ruler change, and not the k-ruler", ok)
# G4. kappa-Poincare (1+1, bicrossproduct basis) momentum composition = the order-dial page's affine stretch-and-shift group
kap = 7.3; kk = -1.0
def kp(p, q):            # (p0, p1) + (q0, q1)
    return (p[0] + q[0], p[1] + np.exp(-p[0]/kap)*q[1])
def aff(x, y):           # (a, c) o (b, d) = (a + b + k ab, c + (1 + k a) d)
    return (x[0] + y[0] + kk*x[0]*y[0], x[1] + (1 + kk*x[0])*y[1])
to_aff = lambda p: ((np.exp(-p[0]/kap) - 1)/kk, p[1])
rng = np.random.default_rng(11); ok = True
for _ in range(40):
    p, q, s = rng.normal(0, 2, 2), rng.normal(0, 2, 2), rng.normal(0, 2, 2)
    ok &= np.allclose(to_aff(kp(p, q)), aff(to_aff(p), to_aff(q)))
    ok &= np.allclose(kp(kp(p, q), s), kp(p, kp(q, s)))
p, q = (1.0, 0.5), (-0.7, 2.0); gap = np.subtract(kp(p, q), kp(q, p))
check("G4 kappa-Poincare's momentum composition (p0 + q0, p1 + e^{-p0/kappa} q1) is the affine stretch-and-shift composition "
      "with 1 + k a = e^{-p0/kappa}: associative, order matters (gap only in the shifted slot)", ok and abs(gap[1]) > 1e-3 and abs(gap[0]) < 1e-12,
      f"order gap = {gap}")
ok = True
for K in [1e2, 1e4, 1e6]:
    comp = (p[0] + q[0], p[1] + np.exp(-p[0]/K)*q[1]); ok &= abs(comp[1] - (p[1] + q[1])) < 3/K
check("G5 as kappa (the Planck-scale dial) -> infinity, kappa-Poincare composition -> plain addition (Poincare), error ~ 1/kappa", ok)
print(f"\n{sum(r[1] for r in results)}/{len(results)} pass")
