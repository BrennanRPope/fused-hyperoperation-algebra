"""Why does quantum mechanics need i?  Checks for topics/why-i (Nov 25, 2026).

Conventions follow topics/order-dial: the pair (a, c) is the map x -> (1 + a)x + c,
mult-add's dial-k element a is the map x -> (1 + k a)x + a (centre -1/k).
Run: python3 checks.py  (needs numpy, sympy)
"""
import numpy as np, sympy as sp
rng = np.random.default_rng(20261125)
results = []
def check(name, ok, detail=""):
    results.append((name, bool(ok)))
    print(("PASS " if ok else "FAIL ") + name + (("  | " + detail) if detail else ""))

# ---------------------------------------------------------------- A. the twist
a, b, c, d, x = sp.symbols('a b c d x', real=True)
f = lambda A, C: (lambda X: (1 + A) * X + C)
gap = sp.expand(f(a, c)(f(b, d)(x)) - f(b, d)(f(a, c)(x)))
check("A1 order gap of (a,c),(b,d) is the wedge a*d - b*c", sp.simplify(gap - (a*d - b*c)) == 0, str(gap))
check("A2 the twist is antisymmetric: swapping the order flips its sign",
      sp.simplify(gap + sp.expand(f(b, d)(f(a, c)(x)) - f(a, c)(f(b, d)(x)))) == 0)
k = sp.symbols('k', real=True)
gk = sp.expand(f(k*a, a)(f(k*b, b)(x)) - f(k*b, b)(f(k*a, a)(x)))
check("A3 two mult-add elements with the same dial k have zero twist (shared centre)", gk == 0)
k2 = sp.symbols('k2', real=True)
gk2 = sp.factor(sp.expand(f(k*a, a)(f(k2*b, b)(x)) - f(k2*b, b)(f(k*a, a)(x))))
check("A4 different dials k, k2: twist = (k - k2)*a*b", sp.simplify(gk2 - (k - k2)*a*b) == 0, str(gk2))

# ------------------------------------------- B. area alone does not pick an i
W = np.array([[0., 1.], [-1., 0.]])          # omega(u,v) = u^T W v = u1 v2 - u2 v1
def omega(u, v): return u @ W @ v
def J_tau(t):                                  # compatible complex structures <-> Siegel half-plane point t
    X, Y = t.real, t.imag
    return np.array([[X, -(X*X + Y*Y)], [1., -X]]) / Y
ok_sq = ok_comp = ok_pos = True; Js = []
for _ in range(200):
    t = complex(rng.normal(), abs(rng.normal()) + 0.05)
    J = J_tau(t); Js.append(J)
    ok_sq &= np.allclose(J @ J, -np.eye(2))
    for _ in range(5):
        u, v = rng.normal(size=2), rng.normal(size=2)
        ok_comp &= np.isclose(omega(J @ u, J @ v), omega(u, v))
        ok_pos &= omega(u, J @ u) > 0
check("B1 every J in the Siegel family squares to -1", ok_sq)
check("B2 ... keeps the oriented area: omega(Ju,Jv) = omega(u,v)", ok_comp)
check("B3 ... is tamed by the orientation: omega(u,Ju) > 0", ok_pos)
check("B4 the area form alone allows a whole half-plane of different J's (no unique i)",
      not np.allclose(Js[0], Js[1]))
# a length picks one: g = standard dot product -> J = g^{-1} omega^T... = rotation by 90 degrees
G = np.eye(2); Jg = np.linalg.solve(G, W.T)    # omega(u,v) = g(Ju, v)
ok = all(np.isclose(omega(u, v), (Jg @ u) @ G @ v) for u, v in rng.normal(size=(50, 2, 2)))
check("B5 area + length -> one J: omega(u,v) = g(Ju,v), J = rotation by 90 deg, J^2 = -1",
      ok and np.allclose(Jg @ Jg, -np.eye(2)) and np.allclose(Jg, [[0, -1], [1, 0]]))
Gr = np.array([[2., .7], [.7, 1.]]); Jr = np.linalg.solve(Gr, W.T)
Jr = Jr / np.sqrt(abs(np.linalg.det(Jr)))      # normalise the area scale
check("B6 a different length gives a different J (still J^2 = -1)",
      np.allclose(Jr @ Jr, -np.eye(2)) and not np.allclose(Jr, Jg))
R = np.diag([1., -1.])                          # orientation flip
check("B7 flipping orientation turns the area to -area and J into -J (i into -i)",
      np.allclose(R.T @ W @ R, -W) and np.allclose(R @ Jg @ R, -Jg))
# the twist as an operator squares to minus (area)^2
th = 1.7
check("B8 the twist raised by a length is theta*J; it squares to -(theta)^2", np.allclose((th*Jg) @ (th*Jg), -th**2*np.eye(2)))

# ------------------------------------------ C. the stretch-and-shift group IS the upper half-plane
def comp(g, h):                                 # (lam, s): x -> lam x + s ; g after h
    return (g[0]*h[0], g[0]*h[1] + g[1])
def z_of(g): return g[1] + 1j*g[0]
ok = True
for _ in range(100):
    g = (abs(rng.normal()) + .1, rng.normal()); h = (abs(rng.normal()) + .1, rng.normal())
    ok &= np.isclose(z_of(comp(g, h)), g[0]*z_of(h) + g[1])
check("C1 element (lam,s) <-> point s + i*lam; composing = acting by z -> lam z + s (a holomorphic map)", ok)
# dial lines
ok = True
for kk in [-2., -.5, .3, 1., 4.]:
    A = rng.normal(size=6)
    pts = [z_of((1 + kk*ai, ai)) for ai in A]
    ok &= all(np.isclose((p.imag - 1) - kk*p.real, 0) for p in pts)   # on the line lam = 1 + k s through i
    ok &= np.isclose(1 + kk*(-1/kk), 0)                                 # meets the edge lam = 0 at s = -1/k
    # closure under the group law = mult-add
    for ai, bi in zip(A[:3], A[3:]):
        g = comp((1 + kk*ai, ai), (1 + kk*bi, bi))
        ok &= np.isclose(g[1], ai + bi + kk*ai*bi) and np.isclose(g[0], 1 + kk*g[1])
check("C2 dial-k elements form the straight line from i to the edge point -1/k; the group law on it is a + b + k*ab", ok)
# every one-parameter subgroup is such a line
ok = True
for _ in range(50):
    al, be = rng.normal(size=2); t = rng.normal()
    lam = np.exp(al*t); s = be*(np.exp(al*t) - 1)/al
    ok &= np.isclose(lam, 1 + (al/be)*s)
check("C3 every one-parameter subgroup exp(t(aX+bY)) is a line through i (dial k = alpha/beta)", ok)
# invariant length and area, and J
lam_, s_, L, S = sp.symbols('lam s L S', positive=True)
gmat = sp.Matrix([[1/lam_**2, 0], [0, 1/lam_**2]])  # coords (s, lam): ds^2 = (ds^2 + dlam^2)/lam^2
# left translation (L,S): (lam, s) -> (L lam, L s + S); Jacobian in (s,lam) coords = L*I
Jac = sp.Matrix([[L, 0], [0, L]])
pulled = sp.simplify(Jac.T * gmat.subs(lam_, L*lam_) * Jac)
check("C4 the Poincare length (ds^2 + dlam^2)/lam^2 is left-invariant", sp.simplify(pulled - gmat) == sp.zeros(2))
Wm = sp.Matrix([[0, 1], [-1, 0]]) / lam_**2         # area form ds ^ dlam / lam^2
Jm = sp.simplify(gmat.inv() * Wm.T)
check("C5 invariant area + invariant length give J = rotation by 90 deg = multiplication by i on z = s + i lam",
      Jm == sp.Matrix([[0, -1], [1, 0]]))
# curvature (Brioschi) for an arbitrary left-invariant metric: G/lam^2 with G any constant inner product
def gauss_curv(E, F, Gg, u, v):
    Eu, Ev, Fu, Fv, Gu, Gv = [sp.diff(q, w) for q in (E, F, Gg) for w in (u, v)]
    Evv, Guu, Fuv = sp.diff(E, v, 2), sp.diff(Gg, u, 2), sp.diff(F, u, v)
    M1 = sp.Matrix([[-Evv/2 + Fuv - Guu/2, Eu/2, Fu - Ev/2], [Fv - Gu/2, E, F], [Gv/2, F, Gg]])
    M2 = sp.Matrix([[0, Ev/2, Gu/2], [Ev/2, E, F], [Gu/2, F, Gg]])
    return sp.simplify((M1.det() - M2.det()) / (E*Gg - F**2)**2)
p, q_, r = sp.symbols('p q r', positive=True)
K0 = gauss_curv(1/lam_**2, 0, 1/lam_**2, s_, lam_)
check("C6 Poincare curvature is -1", sp.simplify(K0 + 1) == 0, str(K0))
Kg = gauss_curv(p/lam_**2, r/lam_**2, q_/lam_**2, s_, lam_)
check("C7 ANY left-invariant length (p, r, q constants) gives constant negative curvature (Milnor 1976)",
      Kg.free_symbols <= {p, q_, r} and all(Kg.subs({p: pv, q_: qv, r: rv}) < 0
          for pv, qv, rv in [(1, 1, 0), (2, 1, .5), (.3, 5, -.9), (1, 4, 1.9)]), str(sp.simplify(Kg)))
# dial lines are hypercycles; k=0 is a horocycle; k=inf (pure stretch) is a geodesic
def dist_to_vertical(zz, p0): return np.arcsinh(abs(zz.real - p0)/zz.imag)
ok = True
for kk in [-2., -.5, .3, 1., 4.]:
    pts = [complex(si, 1 + kk*si) for si in np.linspace(-1/kk*0.9 if kk > 0 else -3, 3, 7) if 1 + kk*si > 0]
    ds = [dist_to_vertical(zz, -1/kk) for zz in pts]
    ok &= np.allclose(ds, ds[0])
check("C8 each dial line keeps a constant hyperbolic distance from the vertical geodesic over its centre (a hypercycle)", ok)
print("NOTE C9 (by definition, not a check): k = 0 is the line lam = 1, a horocycle centred at infinity; pure multiplication is the vertical geodesic s = 0")

# ------------------------------------------- D. why the twist must be a turn (a phase)
# a twist that adds (z + z') and acts by size-keeping multiplication must be a unit complex number
zs = rng.normal(size=20)
ok = all(np.isclose(abs(np.exp(1j*0.8*z)), 1) for z in zs) and not all(np.isclose(abs(np.exp(0.8*z)), 1) for z in zs)
check("D1 real characters e^(lz) change size; only e^(ilz) adds twists while keeping size", ok)
def rot(t): return np.array([[np.cos(t), -np.sin(t)], [np.sin(t), np.cos(t)]])
ok = all(np.allclose(rot(u) @ rot(v), rot(u + v)) for u, v in rng.normal(size=(20, 2)))
h = 1e-6; gen = (rot(h) - rot(-h)) / (2*h)
check("D3 written with real numbers, a size-keeping twist is a rotation; its generator J has J^2 = -1",
      ok and np.allclose(gen @ gen, -np.eye(2), atol=1e-6))
# Heisenberg: twist lands in the centre
def H(a_, b_, c_): return np.array([[1, a_, c_], [0, 1, b_], [0, 0, 1.]])
A_, B_ = H(.7, 0, 0), H(0, 1.3, 0)
Cm = A_ @ B_ @ np.linalg.inv(A_) @ np.linalg.inv(B_)
check("D4 Heisenberg commutator of shift a and boost b is the central element with area a*b", np.allclose(Cm, H(0, 0, .7*1.3)))

# ------------------------------------------- E. Stueckelberg: symmetric observables, antisymmetric twist
ok_anti = ok_imag = True
for n in [2, 3, 4, 6]:
    for _ in range(20):
        M1, M2 = rng.normal(size=(n, n)), rng.normal(size=(n, n))
        A1, A2 = M1 + M1.T, M2 + M2.T
        C = A1 @ A2 - A2 @ A1
        ok_anti &= np.allclose(C.T, -C)
        ev = np.linalg.eigvals(C)
        ok_imag &= np.allclose(ev.real, 0, atol=1e-8)
check("E1 the commutator of two real symmetric matrices (real observables) is antisymmetric (a twist)", ok_anti)
check("E2 its eigenvalues are purely imaginary, in +/- pairs: the i appears even in an all-real theory", ok_imag)
C2 = np.array([[0, -1.], [1, 0]]) * 2.5
check("E3 in two dimensions the twist is theta*J, and (theta*J)^2 = -theta^2", np.allclose(C2 @ C2, -6.25*np.eye(2)))
# Frobenius: the operator J commuting with a real irreducible set can be real, complex or quaternionic
qi = np.array([[0, -1, 0, 0], [1, 0, 0, 0], [0, 0, 0, -1], [0, 0, 1, 0.]])
qj = np.array([[0, 0, -1, 0], [0, 0, 0, 1], [1, 0, 0, 0], [0, -1, 0, 0.]])
qk = qi @ qj
check("E4 quaternion units give three different J's (each squares to -1): the twist alone does not choose C over H",
      all(np.allclose(m @ m, -np.eye(4)) for m in (qi, qj, qk)) and not np.allclose(qi @ qj, qj @ qi))

# ------------------------------------------- F. the affine group's own twist stretches, it does not turn
adX = np.array([[0, 0], [0, 1.]])               # basis X (stretch), Y (shift); [X, Y] = Y
check("F1 in the affine group [X,Y] = Y: ad_X has real eigenvalues 0 and 1 (a stretch, not a turn)",
      np.allclose(sorted(np.linalg.eigvals(adX).real), [0, 1]) and np.allclose(np.linalg.eigvals(adX).imag, 0))
print("NOTE F2 (proved in order-dial J5, J6): the affine twist lands on Y, which is not central, unlike Heisenberg")

# ------------------------------------------- G. the twist's sign = the sign of i (Gelfand-Naimark)
N, Lbox = 4096, 400.
xs = np.linspace(-Lbox/2, Lbox/2, N, endpoint=False); dx = xs[1] - xs[0]
xi = 2*np.pi*np.fft.fftfreq(N, d=dx)
def U(lam, s, psi):                            # (U psi)(x) = |lam|^{-1/2} psi((x - s)/lam)
    return np.interp((xs - s)/lam, xs, psi.real, left=0, right=0)/np.sqrt(abs(lam)) + \
        1j*np.interp((xs - s)/lam, xs, psi.imag, left=0, right=0)/np.sqrt(abs(lam))
def positive_part(psi):
    F = np.fft.fft(psi); return np.fft.ifft(np.where(xi > 0, F, 0))
seed = np.exp(-xs**2/30.) * np.exp(1j*1.5*xs)
psi = positive_part(seed)                       # a positive-frequency wave
def neg_weight(phi):
    F = np.abs(np.fft.fft(phi))**2; return F[xi < 0].sum()/F.sum()
w_after = [neg_weight(U(l, s, psi)) for l, s in [(1.7, 3.), (0.6, -5.), (2.5, 10.)]]
check("G1 stretch-and-shift with lam > 0 keeps a positive-frequency wave positive (the irrep pi+)",
      max(w_after) < 1e-3, f"neg-frequency weight after maps: {[f'{w:.1e}' for w in w_after]}")
refl = np.interp(-xs, xs, psi.real) + 1j*np.interp(-xs, xs, psi.imag)
check("G2 an orientation flip (x -> -x, lam < 0) sends it to negative frequency (pi-)", neg_weight(refl) > 0.999)
check("G3 complex conjugation (i -> -i) does the same: the twist's sign and i's sign are one choice",
      neg_weight(np.conj(psi)) > 0.999)
ws_before = [neg_weight(U(1 + kk*aa, aa, psi)) for kk, aa in [(1, -.5), (.5, 1.), (-2, .3)]]   # 1 + k a > 0
ws_past = [neg_weight(U(1 + kk*aa, aa, psi)) for kk, aa in [(1, -2.), (.5, -3.5), (-2, .8)]]    # 1 + k a < 0
check("G4 mult-add elements before the wall (1 + k a > 0) keep the wave's frequency sign; past the wall (1 + k a < 0) they flip it",
      max(ws_before) < 1e-3 and min(ws_past) > 0.999, f"before {[f'{w:.1e}' for w in ws_before]}, past {[f'{w:.4f}' for w in ws_past]}")

# ------------------------------------------- H. holonomy of a twist = e^(i * area)  (Aharonov-Bohm / Berry / Kostant-Souriau)
def holonomy(path, B=1.0):                     # connection A = (B/2)(-y, x) has curvature B
    ph = 0.
    for (x0, y0), (x1, y1) in zip(path, path[1:] + path[:1]):
        xm, ym = (x0 + x1)/2, (y0 + y1)/2
        ph += (B/2)*(-ym*(x1 - x0) + xm*(y1 - y0))
    return np.exp(1j*ph)
sq = [(0, 0), (1.2, 0), (1.2, .8), (0, .8)]
check("H1 going round a loop of oriented area A multiplies by e^(iA)", np.isclose(holonomy(sq), np.exp(1j*1.2*.8)))
check("H2 going round the other way gives e^(-iA): reversing the twist conjugates i", np.isclose(holonomy(sq[::-1]), np.exp(-1j*.96)))


# ------------------------------------------- P. path integrals as a whole-chain mult-add fold
from scipy.linalg import expm
def wstar(a_, b_, kk):                         # weighted circle operation on matrices: "b first, then a"
    return a_ + b_ + kk * (a_ @ b_)
def fold(steps, kk):                           # steps in time order; later steps act on the left
    F = np.zeros_like(steps[0])
    for s_ in steps:
        F = wstar(s_, F, kk)
    return F
n = 3
steps = [rng.normal(size=(n, n)) + 1j*rng.normal(size=(n, n)) for _ in range(5)]
for kk in [1., -1., -1j, 0.37 + 0.2j]:
    prod = np.eye(n, dtype=complex)
    for s_ in steps: prod = (np.eye(n) + kk*s_) @ prod
    if not np.allclose(np.eye(n) + kk*fold(steps, kk), prod): break
else:
    kk = None
check("P1 a chain of weighted circle steps is a product: 1 + k*(a_N o ... o a_1) = (1 + k a_N)...(1 + k a_1), order kept", kk is None)
# Dyson / subset expansion
from itertools import combinations
F = fold(steps[:4], 0.6)
S = np.zeros((n, n), dtype=complex)
for m in range(1, 5):
    for sub in combinations(range(4), m):
        P_ = np.eye(n, dtype=complex)
        for j in sub: P_ = steps[j] @ P_
        S += 0.6**(m - 1) * P_
check("P2 expanded, the chain is a sum over every time-ordered subset of steps (the Dyson series skeleton)", np.allclose(F, S))
# lattice path sum
Lsites, N_, eps = 4, 4, 0.3
Hh = np.zeros((Lsites, Lsites))
for x_ in range(Lsites):
    Hh[x_, (x_+1) % Lsites] = Hh[(x_+1) % Lsites, x_] = -1.0; Hh[x_, x_] = 0.5
kq = -1j
Tm = np.eye(Lsites) + kq*eps*Hh
prodT = np.linalg.matrix_power(Tm, N_)
from itertools import product as iprod
x0, xN = 0, 2; tot = 0
for mid in iprod(range(Lsites), repeat=N_ - 1):
    path = (x0,) + mid + (xN,)
    w = 1
    for u_, v_ in zip(path, path[1:]): w *= Tm[v_, u_]
    tot += w
check("P3 for a particle hopping on a ring, the chain's entry (end <- start) equals the sum over every lattice path of its step weights",
      np.isclose(prodT[xN, x0], tot), f"{prodT[xN, x0]:.6f} vs {tot:.6f}")
Hr = rng.normal(size=(4, 4)); Hr = (Hr + Hr.T)/2; t_ = 1.0
errs_q, errs_h = [], []
for NN in [10, 100, 1000]:
    errs_q.append(np.linalg.norm(np.linalg.matrix_power(np.eye(4) - 1j*Hr*t_/NN, NN) - expm(-1j*Hr*t_)))
    errs_h.append(np.linalg.norm(np.linalg.matrix_power(np.eye(4) - Hr*t_/NN, NN) - expm(-Hr*t_)))
check("P4 dial k = -i (hbar = 1): the chain of H*dt steps converges to the quantum propagator e^(-iHt); error ~ 1/N",
      errs_q[2] < errs_q[1] < errs_q[0] and errs_q[1]/errs_q[2] > 8, f"{[f'{e:.1e}' for e in errs_q]}")
check("P5 dial k = -1 (real): the same chain converges to the heat kernel e^(-Ht) (Feynman-Kac / imaginary time)",
      errs_h[2] < errs_h[1] < errs_h[0] and errs_h[1]/errs_h[2] > 8, f"{[f'{e:.1e}' for e in errs_h]}")
hh, ee = 1.3, 1e-4
first_order = {kk: (abs(1 + kk*hh*ee) - 1)/ee for kk in [-1, -1j, np.exp(-1j*np.pi/3)]}
check("P6 one short step keeps size to first order only if the dial is imaginary: |1 + k h e| = 1 + Re(k) h e + ...",
      all(np.isclose(v_, np.real(k_)*hh, atol=1e-3) for k_, v_ in first_order.items()) and np.isclose(first_order[-1j], 0, atol=1e-3),
      "first-order size change per unit step: " + ", ".join(f"k={complex(k_):.2f}: {float(v_):.3f}" for k_, v_ in first_order.items()))
# Trotter error is set by the twist
A1 = rng.normal(size=(3, 3)); A2 = rng.normal(size=(3, 3)); A1, A2 = (A1 - A1.T)/2, (A2 - A2.T)/2
NN = 2000
err = np.linalg.matrix_power(expm(A1/NN) @ expm(A2/NN), NN) - expm(A1 + A2)
leading = (A1 @ A2 - A2 @ A1)/(2*NN)
check("P7 slicing two non-commuting steps: the error (e^(A/N) e^(B/N))^N - e^(A+B) is led by the twist [A,B]/(2N)",
      np.linalg.norm(err - leading) < 0.15*np.linalg.norm(leading), f"|err| = {np.linalg.norm(err):.2e}, |[A,B]/2N| = {np.linalg.norm(leading):.2e}")
# Chen / Levy area / Heisenberg
vs = rng.normal(size=(7, 2))
def heis(v, cc=0.): return (np.array(v, float), cc)
def hmul(g, h_): return (g[0] + h_[0], g[1] + h_[1] + 0.5*omega(g[0], h_[0]))
G_ = heis([0, 0])
for v in vs: G_ = hmul(heis(v), G_) if False else hmul(G_, heis(v))
pts = np.vstack([[0, 0], np.cumsum(vs, axis=0)])
poly = np.vstack([pts, pts[:1]])
shoelace = 0.5*np.sum(poly[:-1, 0]*poly[1:, 1] - poly[1:, 0]*poly[:-1, 1])
check("P8 multiplying the steps of a path in the Heisenberg group stores the path's signed (Levy) area: the twist of every pair, summed",
      np.isclose(G_[1], shoelace) and np.isclose(G_[1], 0.5*sum(omega(vs[i], vs[j]) for i in range(7) for j in range(i+1, 7))),
      f"central part {G_[1]:.4f}, shoelace area {shoelace:.4f}")
G2 = heis([0, 0])
for v in vs[:3]: G2 = hmul(G2, heis(v))
G3 = heis([0, 0])
for v in vs[3:]: G3 = hmul(G3, heis(v))
G23 = hmul(G2, G3)
check("P9 Chen's identity at level 2: the record of a joined path is the product of the two records", np.allclose(G23[0], G_[0]) and np.isclose(G23[1], G_[1]))
print("NOTE P10 (follows from D1 + P8; order-dial J2, J4): represented with size kept, that area becomes the phase e^(i*area/hbar) of the path")
haz = rng.uniform(0, .2, size=12)
F_km = 0.
for h_ in haz: F_km = h_ + F_km - h_*F_km
check("P11 at k = -1 the chain of hazards is survival analysis (Kaplan-Meier): event probability = 1 - prod(1 - h_j)",
      np.isclose(F_km, 1 - np.prod(1 - haz)))

print(f"\n{sum(o for _, o in results)}/{len(results)} pass")
