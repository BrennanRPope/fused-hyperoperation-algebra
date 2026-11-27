"""DRAFT checks: is torsion a new dial or a derived property of Affine Composition Theory?

Brennan (Nov 27 2026): "I'd like to take a look at torsion" / "That might be a new dial or a derived property".
Round 28 (S4) found the law's 3-sphere connection is the f = 1 member of the Cartan-Schouten family
nabla_X Y = f[X,Y]/2, torsion 2(f - 1). Here we ask where f comes from.

Run: python3 checks.py   (numpy, scipy, sympy).  Writes checks_output.txt and results.json.
"""
import io, json, math
import numpy as np
import sympy as sp
from scipy.linalg import expm, logm
from scipy.integrate import solve_ivp

_log = io.StringIO(); results = {}; passed = failed = 0
def say(s): print(s, flush=True); _log.write(s + "\n")
def check(name, ok, detail=""):
    global passed, failed
    if ok: passed += 1
    else: failed += 1
    say(("PASS " if ok else "FAIL ") + name + (("  | " + detail) if detail else ""))
rng = np.random.default_rng(20261127)

# ---------------------------------------------------------------- the law in 3D (rounds 12-13, 18, 25, unchanged)
def oplus(a, b, qq=1.0):
    ab = (a * b).sum(-1, keepdims=True); aa = (a * a).sum(-1, keepdims=True); bb = (b * b).sum(-1, keepdims=True)
    return ((1 - 2 * qq * ab - qq * bb) * a + (1 + qq * aa) * b) / (1 - 2 * qq * ab + qq**2 * aa * bb)
def to_s3(x):        # inverse stereographic: the law's chart -> unit 3-sphere in R^4 (round 18)
    r2 = (x * x).sum(-1, keepdims=True); return np.concatenate([(1 - r2), 2 * x], -1) / (1 + r2)

# ---------------------------------------------------------------- SU(2) = unit quaternions, as 2x2 matrices
PAULI = [np.array([[0, 1], [1, 0]], complex), np.array([[0, -1j], [1j, 0]]), np.array([[1, 0], [0, -1]], complex)]
E = [-1j * m for m in PAULI]                      # [E_a, E_b] = 2 eps_abc E_c, like the frame of rounds 25, 28
def mat(v): return sum(v[i] * E[i] for i in range(3))
def vec(M): return np.array([np.real(np.trace(M @ np.conj(e).T)) / 2 for e in E])
def br(A, B): return A @ B - B @ A
def quat_to_mat(p):  # p = (w, x, y, z) -> w I + x E1 + y E2 + z E3  (same product as qmul)
    return p[0] * np.eye(2) + mat(p[1:])
def mpow(M, s): return expm(s * logm(M))
def oplus_f(x, y, f):                             # x, y group elements: the f member's composition
    return mpow(x, 1 - f / 2) @ y @ mpow(x, f / 2)

# ================================================================ T1. The Cartan-Schouten family (exact)
fs = sp.symbols('f', real=True)
def eps(i, j, k): return sp.LeviCivita(i, j, k)
def nab(X, Y):   # nabla_X Y for constant-coefficient left-invariant fields: f/2 [X, Y]
    return [fs / 2 * sum(2 * eps(i, j, k) * X[i] * Y[j] for i in range(3) for j in range(3)) for k in range(3)]
def brk(X, Y): return [sum(2 * eps(i, j, k) * X[i] * Y[j] for i in range(3) for j in range(3)) for k in range(3)]
e = [[1 if i == a else 0 for i in range(3)] for a in range(3)]
T = [sp.simplify(nab(e[0], e[1])[k] - nab(e[1], e[0])[k] - brk(e[0], e[1])[k]) for k in range(3)]
def R(X, Y, Z):
    t1 = nab(X, nab(Y, Z)); t2 = nab(Y, nab(X, Z)); t3 = nab(brk(X, Y), Z)
    return [sp.expand(t1[k] - t2[k] - t3[k]) for k in range(3)]
K = sp.factor(sum(R(e[0], e[1], e[1])[k] * e[0][k] for k in range(3)))
geo = all(sp.simplify(c) == 0 for v in e for c in nab(v, v))
ok1 = (sp.simplify(T[2] - 2 * (fs - 1)) == 0 and T[0] == 0 and T[1] == 0 and sp.simplify(K - (2 * fs - fs**2)) == 0 and geo)
check("T1 (exact; textbook, Cartan-Schouten 1926) the family nabla_X Y = f[X,Y]/2 on the 3-sphere: torsion "
      "T(e1,e2) = 2(f - 1) e3, sectional curvature 2f - f^2, and nabla_X X = 0 for every f, so every member has the "
      "same straight lines (great circles); f only changes how a sideways direction is carried", ok1,
      f"T = {T}, K = {K}")
results["T1"] = {"torsion_e3": str(T[2]), "curvature": str(K), "same_geodesics": bool(geo)}

# ================================================================ T2. Each member is a way to compose: x (+)_f y = x^(1 - f/2) y x^(f/2)
def geodesic_loop(Xv, Yv, f):
    X, Y = mat(Xv), mat(Yv)
    def rhs(t, yv):   # parallel transport in the left trivialisation along exp(tX): y' = -(f/2)[X, y]
        return -(f / 2) * vec(br(X, mat(yv)))
    yv1 = solve_ivp(rhs, (0, 1), Yv, rtol=1e-12, atol=1e-14).y[:, -1]
    return expm(X) @ expm(mat(yv1))           # follow the straight line from x in the carried direction
err2 = 0.0
for _ in range(60):
    Xv, Yv = rng.normal(size=3) * 0.8, rng.normal(size=3) * 0.8; f = rng.uniform(-1, 3)
    err2 = max(err2, np.abs(geodesic_loop(Xv, Yv, f) - oplus_f(expm(mat(Xv)), expm(mat(Yv)), f)).max())
check("T2 (numeric; geodesic loops are textbook, Kikkawa 1964, Sabinin 1999; the closed form for every f is a "
      "one-line consequence we didn't find written) carrying y's direction to x with the f member, then going "
      "straight, gives x (+)_f y = x^(1 - f/2) . y . x^(f/2): f = 0 is the product xy, f = 2 is yx, f = 1 is "
      "x^(1/2) y x^(1/2). f says where x is split around y", err2 < 1e-9, f"max error {err2:.1e} over 60 random x, y, f")
results["T2"] = {"max_err": err2}

# ================================================================ T3. The law's 3-sphere sum is the f = 1 member
errs = {0: 0.0, 1: 0.0, 2: 0.0}; errq = 0.0
for _ in range(500):
    a, b = rng.normal(size=3) * 0.9, rng.normal(size=3) * 0.9
    A, B = quat_to_mat(to_s3(a[None])[0]), quat_to_mat(to_s3(b[None])[0])
    L = quat_to_mat(to_s3(oplus(a[None], b[None]))[0])
    for f in errs:
        Mf = oplus_f(A, B, f); errs[f] = max(errs[f], min(np.abs(Mf - L).max(), np.abs(Mf + L).max()))
    qq = rng.uniform(0.1, 4)   # any positive q: the same law at radius 1/sqrt(q)
    errq = max(errq, np.abs(oplus(a, b, qq) - oplus(np.sqrt(qq) * a, np.sqrt(qq) * b) / np.sqrt(qq)).max())
ok3 = errs[1] < 1e-12 and errs[0] > 0.3 and errs[2] > 0.3 and errq < 1e-12
check("T3 (numeric; textbook symmetric-space fact, Loos 1969 / Kikkawa: the geodesic loop of a symmetric space; "
      "ours: that it is the repo's own sum) the law's 3-sphere sum, mapped to unit quaternions, is exactly "
      "x^(1/2) y x^(1/2), the f = 1 member. Plain products xy (f = 0) and yx (f = 2) miss by order 1. Any q > 0 "
      "is the same law rescaled, so q sets the size of the sphere and never moves f", ok3,
      f"f=1: {errs[1]:.1e}; f=0: {errs[0]:.2f}; f=2: {errs[2]:.2f}; q rescaling {errq:.1e}")
results["T3"] = {str(k): v for k, v in errs.items()}; results["T3"]["q_rescale"] = errq

# ================================================================ T4. Torsion = the order gap at second order; the law's is zero
s_, t_ = sp.symbols('s t'); av = sp.symbols('a1:4'); bv = sp.symbols('b1:4'); qs = sp.symbols('q', positive=True)
aS = sp.Matrix(av) * s_; bS = sp.Matrix(bv) * t_
ab = (aS.T * bS)[0]; aa = (aS.T * aS)[0]; bb = (bS.T * bS)[0]
F = ((1 - 2 * qs * ab - qs * bb) * aS + (1 + qs * aa) * bS) / (1 - 2 * qs * ab + qs**2 * aa * bb)
bil = [sp.simplify(sp.diff(Fi, s_, t_).subs({s_: 0, t_: 0})) for Fi in F]
Fsw = F.subs({**{av[i]: bv[i] for i in range(3)}, **{bv[i]: av[i] for i in range(3)}, s_: t_, t_: s_}, simultaneous=True)
gap = [sp.series(sp.expand((F[i] - Fsw[i]).subs({s_: sp.Symbol('h'), t_: sp.Symbol('h')})), sp.Symbol('h'), 0, 3).removeO()
       for i in range(3)]
law_ok = all(c == 0 for c in bil) and all(sp.simplify(g) == 0 for g in gap)
# the f members, numerically: mixed second derivative of log(x (+)_f y) in log coordinates
def Gamma(f, Xv, Yv, h=1e-4):
    def Fl(sv, tv): return vec(logm(oplus_f(expm(mat(sv * Xv)), expm(mat(tv * Yv)), f)))
    return (Fl(h, h) - Fl(h, -h) - Fl(-h, h) + Fl(-h, -h)) / (4 * h * h)
err4 = 0.0; fvals = [-1, 0, 0.5, 1, 1.5, 2, 3]; tors = []
Xv, Yv = rng.normal(size=3), rng.normal(size=3); brXY = vec(br(mat(Xv), mat(Yv)))
for f in fvals:
    Tf = Gamma(f, Xv, Yv) - Gamma(f, Yv, Xv)   # relative-locality torsion (antisymmetric part of the bilinear term)
    err4 = max(err4, np.abs(Tf - (1 - f) * brXY).max()); tors.append(float(np.linalg.norm(Tf) / np.linalg.norm(brXY)))
check("T4 (exact for the law, numeric for the family; the dictionary is Amelino-Camelia, Freidel, Kowalski-Glikman and "
      "Smolin 2011) torsion is the order gap at second order. The law's sum has NO bilinear term (exact, any q), so "
      "a (+) b - b (+) a starts at third order: its torsion is forced to zero by the law itself. The f member's "
      "bilinear order gap is (1 - f)[X, Y]", law_ok and err4 < 1e-6,
      f"law bilinear term {bil}; gap to 2nd order {gap}; family max error {err4:.1e}")
results["T4"] = {"law_bilinear_zero": bool(law_ok), "family_torsion_over_bracket": dict(zip(map(str, fvals), tors))}

# ================================================================ T5. What pins f = 1: the law's identities
def rand_su2(sc=1.0): return expm(mat(rng.normal(size=3) * sc))
fgrid = np.linspace(-1, 3, 81); ai_err, li_err, as_err, cm_err = [], [], [], []
for f in fgrid:
    m_ai = m_li = m_as = m_cm = 0.0
    for _ in range(8):
        x, y, z = rand_su2(), rand_su2(), rand_su2()
        inv = np.linalg.inv
        m_ai = max(m_ai, np.abs(inv(oplus_f(x, y, f)) - oplus_f(inv(x), inv(y), f)).max())     # automorphic inverse
        m_li = max(m_li, np.abs(oplus_f(inv(x), oplus_f(x, y, f), f) - y).max())               # left inverse
        m_as = max(m_as, np.abs(oplus_f(oplus_f(x, y, f), z, f) - oplus_f(x, oplus_f(y, z, f), f)).max())
        xs = mpow(x, rng.uniform(-0.9, 0.9))                                                    # same line as x
        m_cm = max(m_cm, np.abs(oplus_f(x, xs, f) - x @ xs).max())
    ai_err.append(m_ai); li_err.append(m_li); as_err.append(m_as); cm_err.append(m_cm)
ai_err, as_err = np.array(ai_err), np.array(as_err)
i1 = int(np.argmin(np.abs(fgrid - 1))); i0 = int(np.argmin(np.abs(fgrid))); i2 = int(np.argmin(np.abs(fgrid - 2)))
law_ai = 0.0
for _ in range(200):
    a, b = rng.normal(size=3), rng.normal(size=3)
    law_ai = max(law_ai, np.abs(oplus(-a, -b) + oplus(a, b)).max())
ok5 = (ai_err[i1] < 1e-10 and np.delete(ai_err, i1).min() > 1e-3 and max(li_err) < 1e-10 and max(cm_err) < 1e-10
       and as_err[i0] < 1e-10 and as_err[i2] < 1e-10 and np.delete(as_err, [i0, i2]).min() > 1e-3 and law_ai < 1e-12)
check("T5 (numeric; the algebra is textbook loop theory, the conclusion for the law is ours) which members keep the "
      "law's rules. The law's sum has the automorphic inverse, (-a) (+) (-b) = -(a (+) b); among the members only "
      "f = 1 has it. Every member keeps the left inverse. Only f = 0, 2 are associative. And on any single line "
      "(inputs that commute) every member gives the same answer, the scalar law, so f is invisible to a (+)_k b",
      ok5, f"automorphic inverse: {ai_err[i1]:.0e} at f=1, >= {np.delete(ai_err, i1).min():.1e} elsewhere on the grid; "
           f"left inverse {max(li_err):.0e}; associative only at f=0,2 ({as_err[i0]:.0e}, {as_err[i2]:.0e}); "
           f"one line {max(cm_err):.0e}; law's own automorphic inverse {law_ai:.0e}")
results["T5"] = {"f": fgrid.tolist(), "auto_inverse_err": ai_err.tolist(), "assoc_err": as_err.tolist()}

# ================================================================ T6. f is the order-dial page's lambda
def assoc3(f, X, Y, Z, epsl=1e-3):
    x, y, z = expm(mat(epsl * X)), expm(mat(epsl * Y)), expm(mat(epsl * Z))
    return vec(logm(oplus_f(x, oplus_f(y, z, f), f)) - logm(oplus_f(oplus_f(x, y, f), z, f))) / epsl**3
X3, Y3, Z3 = [rng.normal(size=3) for _ in range(3)]
ref = vec(br(br(mat(X3), mat(Y3)), mat(Z3)))
err6 = max(np.abs(assoc3(f, X3, Y3, Z3) - f * (2 - f) / 8 * ref).max() / np.abs(ref).max() for f in fvals)
an, bn, cn = sp.symbols('a b c', commutative=False); lam = sp.symbols('lambda')
mu = lambda x, y: x + y + lam * x * y + (1 - lam) * y * x
A_mut = sp.expand(mu(mu(an, bn), cn) - mu(an, mu(bn, cn)))
mut_ok = sp.expand(A_mut - lam * (1 - lam) * (bn * (an * cn - cn * an) - (an * cn - cn * an) * bn)) == 0
gap_ok = sp.expand(mu(an, bn) - mu(bn, an) - (2 * lam - 1) * (an * bn - bn * an)) == 0
fl = sp.expand((1 - (1 - fs / 2)) * (1 - fs / 2) - (2 * fs - fs**2) / 4) == 0 and sp.expand((2 * (1 - fs / 2) - 1) - (1 - fs)) == 0
check("T6 (exact + numeric; the mutation family is order-dial section 3, textbook quasi-associative algebras) the "
      "torsion parameter is the order-dial page's lambda, with lambda = 1 - f/2. The mutation product "
      "a + b + lambda ab + (1 - lambda) ba has order gap (2 lambda - 1)[a,b] = (1 - f)[a,b] and grouping error "
      "lambda(1 - lambda)[b,[a,c]] = (2f - f^2)/4 [b,[a,c]]; the group members' third-order associator is "
      "(2f - f^2)/8 [[X,Y],Z]. Both are the torsion line and the curvature curve of T1",
      mut_ok and gap_ok and fl and err6 < 5e-3, f"associator vs (2f - f^2)/8 [[X,Y],Z]: rel. error {err6:.1e} (O(eps))")
results["T6"] = {"assoc_rel_err": float(err6)}

# ================================================================ T7. Where the law does carry torsion: its chains, set by k
k_ = sp.symbols('k'); a1, c1, b1, e1 = sp.symbols('a c b e')
aff = sp.Matrix([a1 + b1 + k_ * a1 * b1, c1 + e1 + k_ * a1 * e1])           # (a,c) o_k (b,e), lift I + kA
Mx = lambda u, v: sp.Matrix([[u, v], [0, 0]])
chain_ok = sp.simplify(Mx(a1, c1) + Mx(b1, e1) + k_ * Mx(a1, c1) * Mx(b1, e1) - Mx(aff[0], aff[1])) == sp.zeros(2)
affsw = aff.subs({a1: b1, b1: a1, c1: e1, e1: c1}, simultaneous=True)
aff_gap = sp.simplify(aff - affsw)                                           # = k (0, a e - b c): torsion k on (stretch, shift)
g1 = sp.symbols('g1:4'); h1 = sp.symbols('h1:4')
p3 = sp.Matrix(sp.symbols('p1:3')); r3 = sp.Matrix(sp.symbols('r1:3'))
def affc(P, Q):   # associativity of the affine chain
    return sp.Matrix([P[0] + Q[0] + k_ * P[0] * Q[0], P[1] + Q[1] + k_ * P[0] * Q[1]])
P, Q, S = sp.Matrix([g1[0], h1[0]]), sp.Matrix([g1[1], h1[1]]), sp.Matrix([g1[2], h1[2]])
assoc_ok = sp.simplify(affc(affc(P, Q), S) - affc(P, affc(Q, S))) == sp.zeros(2, 1)
dd = sp.symbols('d'); x_, y_ = sp.symbols('x y')
CM = (x_ + y_ + k_ * x_ * y_) / (1 - dd * x_ * y_)
cm_comm = sp.simplify(CM - CM.subs({x_: y_, y_: x_}, simultaneous=True)) == 0
# quaternion chain a o_k b = a + b + k ab on pure quaternions: order gap 2k (a x b)
Aq, Bq = rng.normal(size=3), rng.normal(size=3); kk = 0.37
qc = lambda A, B: mat(A) + mat(B) + kk * mat(A) @ mat(B)
qgap = qc(Aq, Bq) - qc(Bq, Aq)
q_ok = np.abs(qgap - 2 * kk * mat(np.cross(Aq, Bq)) ).max() < 1e-12
ok7 = chain_ok and assoc_ok and sp.simplify(aff_gap - sp.Matrix([0, k_ * (a1 * e1 - b1 * c1)])) == sp.zeros(2, 1) and cm_comm and q_ok
check("T7 (exact; known, kappa-Poincare's momentum group, Majid-Ruegg 1994, and the repo's why-i/kappa-poincare K5b) "
      "the law's chains A o_k B = A + B + k AB multiply the lifts I + kA, which is the f = 0 member: flat "
      "(associative) with torsion. On the stretch-and-shift chain the order gap is k (ae - bc), so the torsion is "
      "k itself; on quaternions it is 2k (a x b). Any commuting law, the two-dial (a + b + k ab)/(1 - d ab) "
      "included, has none, so d never makes torsion", ok7,
      f"affine gap {list(aff_gap)}; two-dial law commutative: {cm_comm}; quaternion gap = 2k a x b: {q_ok}")
results["T7"] = {"affine_gap": str(list(aff_gap)), "two_dial_commutes": bool(cm_comm)}

# ================================================================ T8. Against Einstein-Cartan and the bounds (order of magnitude)
hbar_c_GeVm = 1.973269804e-16; hbar_c_Jm = 3.16152677e-26; G = 6.67430e-11; c = 2.99792458e8
kappa_E = 8 * math.pi * G / c**4
a0 = 2.51e27                                            # round 27's lower bound on the radius of curvature (m)
flat_member = 2 * hbar_c_GeVm / a0                      # |T| = 2|f - 1|/a at f = 0 or 2
n_nuc = 0.16e45                                         # nucleons per m^3 (0.16 fm^-3), all spins lined up
T_EC = kappa_E * hbar_c_Jm * n_nuc / 2                  # Cartan's equation, torsion ~ kappa * spin density (1/m)
T_EC_GeV = T_EC * hbar_c_GeVm
bound = 1e-31
ok8 = flat_member < bound * 1e-10 and T_EC_GeV < bound * 1e-8
check("T8 (order of magnitude; Einstein-Cartan is textbook, Kibble 1961, Sciama 1962, Hehl et al. 1976; bounds "
      "Kostelecky, Russell and Tasson 2008) Einstein-Cartan ties torsion to spin density by an algebraic equation, "
      "so it is zero in empty space, as the law's f = 1 is. Fully spin-aligned nuclear matter would give about "
      f"{T_EC_GeV:.0e} GeV; the flat members on today's cosmic slices give {flat_member:.1e} GeV; both far below "
      "the ~1e-31 GeV bounds", ok8,
      f"EC: kappa hbar c n/2 = {T_EC:.1e} 1/m = {T_EC_GeV:.1e} GeV; flat members 2 hbar c / a0 = {flat_member:.2e} GeV")
results["T8"] = {"EC_polarized_nuclear_GeV": T_EC_GeV, "flat_members_cosmic_GeV": flat_member, "bound_GeV": bound}

say(f"\n{passed} / {passed + failed} pass")
results["summary"] = {"passed": passed, "failed": failed}
open("checks_output.txt", "w").write(_log.getvalue())
json.dump(results, open("results.json", "w"), indent=1, default=float)
