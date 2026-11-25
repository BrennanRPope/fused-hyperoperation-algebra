"""Checks for topics/why-i/kappa-poincare (Nov 25, 2026). Run: python3 checks.py
Is the weighted circle operation a o_k b = a + b + k ab one of the deformed momentum composition laws of
kappa-Poincare / doubly special relativity (DSR), or a special case of one?"""
import numpy as np, sympy as sp
from scipy.spatial.transform import Rotation as Rot
results = []
def check(name, ok, detail=""):
    results.append((name, bool(ok))); print(("PASS " if ok else "FAIL ") + name + (("  | " + detail) if detail else ""))
circ = lambda a, b, k: a + b + k*a*b
rng = np.random.default_rng(3)

# ---------- K1. massless momenta in kappa-Poincare (Majid-Ruegg bicrossproduct basis) compose by o_k with k = -1/kappa ----------
K = 1.0   # Planck energy kappa (units); coproduct Delta p_i = p_i (x) 1 + e^{-p0/kappa} (x) p_i
casimir = lambda p0, p: (2*K*np.sinh(p0/(2*K)))**2 - np.exp(p0/K)*p**2      # Majid-Ruegg mass Casimir, 1 space direction
pm = lambda p0: K*(1 - np.exp(-p0/K))                                          # massless branch, moving in +x
ok = all(abs(casimir(e, pm(e))) < 1e-9*max(1, np.exp(e)) for e in [0.01, 0.5, 2, 6])
check("K1a massless momentum in kappa-Poincare is |p| = kappa (1 - e^{-p0/kappa}) (zero of the Majid-Ruegg Casimir); "
      "it never exceeds kappa however large the energy", ok and pm(5) < K and pm(40) <= K)
ok = True
for _ in range(50):
    e1, e2 = rng.uniform(0, 5, 2)
    p1, p2 = pm(e1), pm(e2)
    composed = p1 + np.exp(-e1/K)*p2                       # coproduct: p (+) q = p + e^{-p0/kappa} q, energies add
    ok &= np.isclose(composed, pm(e1 + e2)) and np.isclose(composed, circ(p1, p2, -1/K))
check("K1b for two massless particles moving the same way, kappa-Poincare composes their momenta by exactly a + b + k ab "
      "with k = -1/kappa, and the result is again massless (ours as stated; follows in one line from Majid-Ruegg 1994)", ok)
check("K1c the wall of o_k at -1/k = kappa is kappa-Poincare's maximum momentum: kappa o_k anything = kappa",
      all(np.isclose(circ(K, b, -1/K), K) for b in [-3, 0.2, 0.9, 5]))

# ---------- K2. the full momentum group: affine stretch-and-shift = matrix o_k ----------
def kp31(p, q):          # 3+1 bicrossproduct: (p0 + q0, p_vec + e^{-p0/kappa} q_vec)
    return np.r_[p[0] + q[0], p[1:] + np.exp(-p[0]/K)*q[1:]]
def as_matrix(p, k=-1/K):  # element A with lift I + kA = [[e^{-p0/kappa}, k p_vec],[0, I3]]... (5x5 affine form, scaled)
    L = np.eye(4); L[0, 0] = np.exp(-p[0]/K); L[0, 1:] = p[1:]
    return (L - np.eye(4))/k
def from_matrix(A, k=-1/K):
    L = np.eye(4) + k*A; return np.r_[-K*np.log(L[0, 0]), L[0, 1:]]
ok = True
for _ in range(40):
    p, q = rng.normal(0, 1.5, 4), rng.normal(0, 1.5, 4)
    A, B = as_matrix(p), as_matrix(q); kk = -1/K
    ok &= np.allclose(from_matrix(A + B + kk*A@B), kp31(p, q))
check("K2 the whole 3+1 kappa-Poincare momentum composition is a matrix chain A o_k B = A + B + k AB with lift "
      "I + kA = [[e^{-p0/kappa}, p], [0, 1]]: the group AN(3) (Majid-Ruegg 1994; Kowalski-Glikman 2002-03)", ok)

# ---------- K3. Magueijo-Smolin (DSR2) energy composition sits on the repo's (k, d) plane, where the two walls merge ----------
E1, E2, l, w = sp.symbols("E1 E2 l w", real=True)
U = lambda E: E/(1 - l*E); Uinv = lambda y: y/(1 + l*y)
ms = sp.simplify(Uinv(U(E1) + U(E2)))
kd = (E1 + E2 + (-2*l)*E1*E2)/(1 - l**2*E1*E2)
walls = sp.solve(sp.Eq(l**2*w**2 - 2*l*w + 1, 0), w)
check("K3 Magueijo-Smolin's energy composition U^-1(U(E) + U(E')), U(E) = E/(1 - l E) (as in Judes-Visser 2002), is "
      "(E + E' + k EE')/(1 - d EE') with k = -2l, d = l^2: the repo's plane on the line k^2 = 4d (kappa = d - k^2/4 = 0), "
      "where the two walls merge into one double wall at the invariant Planck energy 1/l (ours as a placement)",
      sp.simplify(ms - kd) == 0 and walls == [1/l] and sp.simplify(kd.subs(E2, 1/l) - 1/l) == 0, f"composition = {sp.factor(ms)}")

# ---------- K4. a one-number law can't tell DSR models apart ----------
ok = True
for k in [-1.0, 0.7]:
    r = lambda x: np.log1p(k*x)/k
    for _ in range(20):
        a, b = rng.uniform(-.4, .4, 2); ok &= np.isclose(r(circ(a, b, k)), r(a) + r(b))
lm = 0.3; Un = lambda E: E/(1 - lm*E)
msn = lambda a, b: (a + b - 2*lm*a*b)/(1 - lm**2*a*b)
for _ in range(20):
    a, b = rng.uniform(-.4, .4, 2); ok &= np.isclose(Un(msn(a, b)), Un(a) + Un(b))
check("K4 every one-coordinate composition (o_k, Magueijo-Smolin's energy law, plain addition) is addition after a change "
      "of ruler, so matching a 1D law says nothing model-specific: DSR bases are coordinates on one momentum space "
      "(Kowalski-Glikman 2002; Kowalski-Glikman-Nowak 2002)", ok)

# ---------- K5. relative-locality geometry: torsion and curvature of the composition ----------
def numeric_geometry(comp, n, h=1e-3):
    """Gamma^{ab}_c = -d^2 (p+q)_c / dp_a dq_b at 0; torsion = antisymmetric part; curvature ~ associator at third order"""
    G = np.zeros((n, n, n))
    for a in range(n):
        for b in range(n):
            ea, eb = np.eye(n)[a]*h, np.eye(n)[b]*h
            G[a, b] = -(comp(ea, eb) - comp(ea, 0*eb) - comp(0*ea, eb))/h**2
    T = G - G.transpose(1, 0, 2)
    s = 0.05; assoc = 0
    for _ in range(10):
        p, q, r = rng.normal(0, s, (3, n))
        assoc = max(assoc, np.abs(comp(comp(p, q), r) - comp(p, comp(q, r))).max())
    return G, T, assoc
G1, T1, A1 = numeric_geometry(lambda p, q: np.atleast_1d(circ(p[0], q[0], -1.0)), 1)
G2, T2, A2 = numeric_geometry(lambda p, q: np.r_[p[0] + q[0], p[1] + np.exp(-p[0]/K)*q[1]], 2)
check("K5a scalar o_k: connection Gamma = -k, no torsion, no curvature (associative): trivial momentum-space geometry",
      np.isclose(G1[0, 0, 0], 1.0, atol=1e-2) and np.allclose(T1, 0) and A1 < 1e-12)
check("K5b kappa-Poincare (the affine chain): torsion 1/kappa (it is the order gap, p+q - q+p ~ -T pq) and zero curvature "
      "(associative): a group manifold, flat connection with torsion (Amelino-Camelia, Freidel, Kowalski-Glikman, Smolin 2011)",
      np.isclose(abs(T2[0, 1, 1]), 1/K, atol=1e-2) and A2 < 1e-12, f"T^(01)_1 = {T2[0,1,1]:.4f}")

# ---------- K6. Snyder: coset momentum space, non-associative, so no chain reaches it ----------
beta = 1.0
def boost(p):            # spatial Snyder momenta (2 components) -> rotation taking the pole to the Beltrami point
    v = np.r_[np.sqrt(beta)*p, 0.0]; th = np.arctan(np.linalg.norm(v))
    if np.linalg.norm(v) < 1e-15: return Rot.identity()
    axis = np.cross([0, 0, 1], np.r_[v[:2], 0]); axis /= np.linalg.norm(axis)
    return Rot.from_rotvec(axis*th)
def snyder(p, q):
    x = (boost(p)*boost(q)).apply([0, 0, 1]); return x[:2]/x[2]/np.sqrt(beta)
p, q, r = np.array([.4, .1]), np.array([-.2, .5]), np.array([.3, -.3])
assoc = np.abs(snyder(snyder(p, q), r) - snyder(p, snyder(q, r))).max()
comm = np.abs(snyder(p, q) - snyder(q, p)).max()
G3, T3, A3 = numeric_geometry(lambda a, b: snyder(a, b), 2)
Rp, Rq, Rr = boost(p), boost(q), boost(r)
grp_assoc = np.abs(((Rp*Rq)*Rr).as_matrix() - (Rp*(Rq*Rr)).as_matrix()).max()
check("K6 Snyder's momentum composition (coset of the de Sitter-type space, here its spatial sphere slice) is non-associative "
      "(curvature), with no torsion at leading order; chains of o_k are always associative, so they reach Snyder only by keeping "
      "the leftover rotation, i.e. as the full matrix group, the same move as Thomas-Wigner (Girelli-Livine 2010)",
      assoc > 1e-3 and np.allclose(T3, 0, atol=1e-2) and grp_assoc < 1e-12,
      f"associator = {assoc:.4f}, commutator = {comm:.4f}")

# ---------- K7. the scalar dials compared: spacetime signature vs momentum-space Planck dial ----------
c = 1.0; vadd = lambda u, v: (u + v)/(1 + u*v/c**2)
ok = True
for _ in range(20):
    u, v = rng.uniform(-.9, .9, 2); A = lambda x: 2*x/(1 - x)          # repo velocity coordinate at c = 1; composes by o_k with k = 1
    ok &= np.isclose(A(vadd(u, v)), circ(A(u), A(v), 1.0))
check("K7 the same o_k carries two different physical dials: velocity addition (spacetime signature kappa2 = -1/c^2, repo's "
      "a = 2u/(1 - u), k = 1 at c = 1) and the kappa-Poincare momentum law (Planck energy, k = -1/kappa): DSR keeps the "
      "first and adds the second", ok)
print(f"\n{sum(r[1] for r in results)}/{len(results)} pass")
