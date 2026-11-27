"""Checks for topics/why-i/rotations-force-i (Nov 26, 2026). Run: python3 checks.py
Brennan's question: can we prove quantum mechanics needs complex numbers because its steps are rotations?
Setting: a real state space R^n. Observables are symmetric matrices. A Hamiltonian H generates steps exp(t A(H)).
Assumptions (numbered as in the README):
  A1 rotations:        A(H) is antisymmetric (each step keeps size and can be undone)
  A2 conservation:     [A(H), H] = 0
  A3 Planck-Einstein:  a state of definite energy E turns steadily at rate |E|/hbar, i.e. A(H)^2 = -H^2/hbar^2
  A4 energy zero:      H + c I is the same physics plus a common turn: A(H + cI) = A(H) + c A(I)
  A5 maximal:          every symmetric operator that commutes with the common turn is an observable
  A6 composites:       a constant added to the total energy turns the joint state the same whichever part it's booked to
"""
import numpy as np, sympy as sp, itertools
results = []
def check(name, ok, detail=""):
    results.append((name, bool(ok))); print(("PASS " if ok else "FAIL ") + name + (("  | " + detail) if detail else ""))
rng = np.random.default_rng(26); hbar = 1.0
J2 = np.array([[0., -1.], [1., 0.]])
def Jstd(m): return np.kron(np.eye(m), J2)                       # the standard quarter turn on R^(2m)
def realify(h):                                                  # complex m x m -> real 2m x 2m, i -> J2 in each slot
    m = h.shape[0]; R = np.zeros((2*m, 2*m))
    for a in range(m):
        for b in range(m):
            z = h[a, b]; R[2*a:2*a+2, 2*b:2*b+2] = [[z.real, -z.imag], [z.imag, z.real]]
    return R
def rand_herm(m):
    x = rng.normal(size=(m, m)) + 1j*rng.normal(size=(m, m)); return (x + x.conj().T)/2
def rand_sym(n): x = rng.normal(size=(n, n)); return (x + x.T)/2
def rand_skew(n): x = rng.normal(size=(n, n)); return (x - x.T)/2
from scipy.linalg import expm, null_space, sqrtm

# ---------- 0. complex quantum mechanics satisfies A1-A4 ----------
ok = True
for m in [1, 2, 3]:
    J = Jstd(m); A = lambda H: J @ H / hbar
    for _ in range(10):
        H = realify(rand_herm(m)); c = rng.normal()
        ok &= np.allclose(A(H).T, -A(H)) and np.allclose(A(H) @ H, H @ A(H)) and np.allclose(A(H) @ A(H), -H @ H)
        ok &= np.allclose(A(H + c*np.eye(2*m)), A(H) + c*A(np.eye(2*m)))
        U = expm(A(H)*0.7); ok &= np.allclose(U.T @ U, np.eye(2*m))
check("Q0 ordinary quantum mechanics, written in real numbers with i -> J, satisfies A1-A4 with A(H) = J H / hbar", ok)

# ---------- 1. a rotation needs a plane ----------
ok = all(abs(np.linalg.det(rand_skew(n))) < 1e-12 for n in [1, 3, 5, 7] for _ in range(10))
check("T1 in odd dimension every rotation generator has a frozen direction (det = 0), so no J with J^2 = -1 exists: "
      "a single real energy level (a line) can't turn; turning needs a plane", ok)
ok = True
for n in [2, 4, 6]:
    A = rand_skew(n); w, V = np.linalg.eig(A)
    ok &= np.allclose(w.real, 0) and np.allclose(np.sort(w.imag), -np.sort(w.imag)[::-1])
check("T1b every rotation generator splits space into turning planes, each with a +/- pair of frequencies (Aste 2019; Vaughan 2026)", ok)

# ---------- 2. one dynamics alone: a J always exists (positive frequency), but it depends on H ----------
ok = True
for n in [2, 4, 6]:
    A = rand_skew(n); H = hbar*np.real(sqrtm(-A @ A)); J = A @ np.linalg.inv(H/hbar)
    ok &= np.allclose(J @ J, -np.eye(n)) and np.allclose(J.T, -J) and np.allclose(A @ A, -H @ H/hbar**2) and np.allclose(A @ H, H @ A)
check("T2 for one dynamics, A1-A3 always hold with H = hbar sqrt(-A^2) and J = A |A|^-1 (positive-frequency J: Segal, "
      "Ashtekar-Magnon; Aste 2019; Vaughan 2026): rotations alone give a J, nothing more", ok)
# two dynamics on R^4 with different J's, each fine on its own (real QM)
P = np.eye(4)[[0, 2, 1, 3]]
JA = Jstd(2); JB = P @ JA @ P.T                       # turning planes (1,2),(3,4) vs (1,3),(2,4)
HA = np.diag([1., 1., 2., 2.]); HB = P @ HA @ P.T
okA = np.allclose((JA @ HA) @ (JA @ HA), -HA @ HA) and np.allclose(JA @ HA, HA @ JA)
okB = np.allclose((JB @ HB) @ (JB @ HB), -HB @ HB) and np.allclose(JB @ HB, HB @ JB)
check("C2 counterexample without A4: two Hamiltonians on R^4, each satisfying A1-A3 with its own quarter turn, and the two "
      "quarter turns differ (real quantum theory has no shared i)", okA and okB and not np.allclose(JA, JB))

# ---------- 3. the main theorem: A4 makes the quarter turn universal ----------
# Proof step (a): J := hbar A(I) has J^2 = -1 by A3 applied to H = I.
# Proof step (b): A2 for H + cI and A4 give [J, H] = 0 for every H.
# Proof step (c): A3 for H + cI, linear term in c, gives B J + J B = -2H with B = hbar A(H).
# Proof step (d): B = JH + E with E antisymmetric, [E, H] = 0, E J = -J E; then B^2 = -H^2 forces E^2 = 0, so E = 0.
# Symbolic test of (c)-(d) on R^4: unknown antisymmetric B, generic H commuting with J.
b = sp.symbols("b0:6"); h = sp.symbols("h0:4", real=True)
Bm = sp.Matrix([[0, b[0], b[1], b[2]], [-b[0], 0, b[3], b[4]], [-b[1], -b[3], 0, b[5]], [-b[2], -b[4], -b[5], 0]])
Js = sp.Matrix(Jstd(2).astype(int))
hc = sp.Matrix([[h[0], h[2] + sp.I*h[3]], [h[2] - sp.I*h[3], h[1]]])
Hs = sp.zeros(4, 4)
for a_ in range(2):
    for b_ in range(2):
        z = hc[a_, b_]; re, im = sp.re(z), sp.im(z); Hs[2*a_:2*a_+2, 2*b_:2*b_+2] = sp.Matrix([[re, -im], [im, re]])
assert sp.simplify(Hs*Js - Js*Hs) == sp.zeros(4, 4)
lin = list(Bm*Js + Js*Bm + 2*Hs) + list(Bm*Hs - Hs*Bm)
sol = sp.solve(lin, b, dict=True)
Bgen = Bm.subs(sol[0]); E = sp.simplify(Bgen - Js*Hs)
free = sorted(E.free_symbols - set(h), key=str)
quad = sp.simplify(Bgen*Bgen + Hs*Hs)
# E is what's left after the linear conditions; the quadratic condition must kill it
eqs = [sp.expand(q) for q in quad if q != 0]
sol2 = sp.solve(eqs, free, dict=True) if free else [{}]
only_zero = all(all(sp.simplify(v) == 0 for v in s.values()) for s in sol2) and len(sol2) >= 1
check("T3 symbolic, R^4: once [J, H] = 0, the conditions B J + J B = -2H, [B, H] = 0 and B^2 = -H^2 have the single "
      "solution B = J H (the leftover E anticommuting with J is forced to zero)", only_zero, f"leftover after linear steps: {free}")
# numeric: the key lemma E antisymmetric and E^2 = 0 implies E = 0, since trace(E^T E) = -trace(E^2)
ok = all(np.isclose(np.trace(X.T @ X), -np.trace(X @ X)) for X in [rand_skew(n) for n in range(2, 9)])
check("T3b the last step: for antisymmetric E, |E|^2 = trace(E^T E) = -trace(E^2), so E^2 = 0 forces E = 0", ok)
# degenerate case H = lam I: the linear conditions leave room, and only B^2 = -H^2 removes it
lam = sp.symbols("lam", positive=True)
Hd = lam*sp.eye(4)
lin_d = list(Bm*Js + Js*Bm + 2*Hd) + list(Bm*Hd - Hd*Bm)
sol_d = sp.solve(lin_d, b, dict=True)[0]
Bd = Bm.subs(sol_d); free_d = sorted(Bd.free_symbols - {lam}, key=str)
quad_d = [sp.expand(q) for q in (Bd*Bd + Hd*Hd) if sp.expand(q) != 0]
# over the reals: some entry of B^2 + H^2 is a nonzero multiple of the sum of squares of the leftovers
sumsq = sum(v**2 for v in free_d)
sos_entry = any(sp.simplify(q/sumsq).is_number and sp.simplify(q/sumsq) != 0 for q in quad_d)
real_only_zero = sos_entry and all(sp.simplify(q.subs({v: 0 for v in free_d})) == 0 for q in quad_d)
sol_dq = [sp.factor(q) for q in quad_d][:2]
check("T3a symbolic, R^4, degenerate energy (H = lam I): the linear conditions leave a 2-parameter family B = lam J + E with E "
      "anticommuting with J; the Planck-Einstein condition B^2 = -H^2 then forces E = 0 over the reals",
      len(free_d) == 2 and real_only_zero, f"free after linear steps: {free_d}; entries of B^2 + H^2: {sol_dq}")
# numeric: any perturbation E (antisymmetric, anticommuting with J, commuting with H) breaks Planck-Einstein
ok = True
for m in [2, 3]:
    n = 2*m; J = Jstd(m); lam_ = 1.7; H = lam_*np.eye(n)
    for _ in range(10):
        X = rand_skew(n); E = (X + J @ X @ J)/2            # the part of X that anticommutes with J (still antisymmetric)
        assert np.allclose(E @ J, -J @ E) and np.allclose(E.T, -E)
        Bp = J @ H + 0.1*E
        ok &= np.allclose(Bp @ H, H @ Bp) and np.allclose(Bp.T, -Bp) and not np.allclose(Bp @ Bp, -H @ H)
check("T3c numerically: such a perturbation keeps rotations (A1) and conservation (A2) but breaks Planck-Einstein (A3): "
      "two states of the same energy would turn at different rates", ok)

# ---------- 4. real quantum theory can't meet A4 ----------
def commutant_dim(gens, n):
    M = np.vstack([np.kron(np.eye(n), G) - np.kron(G.T, np.eye(n)) for G in gens]); return n*n - np.linalg.matrix_rank(M)
ok = True
for n in [2, 3, 4]:
    basis = []
    for i in range(n):
        for j in range(i, n):
            S = np.zeros((n, n)); S[i, j] = S[j, i] = 1; basis.append(S)
    ok &= commutant_dim(basis, n) == 1
check("T4 in real quantum theory every symmetric matrix is an observable, and only multiples of the identity commute with "
      "all of them: no quarter turn J can commute with every observable, so A4 (with A3) fails", ok)

# ---------- 5. complex vs quaternionic: A5 decides ----------
# quaternions as 4x4 real: left and right multiplication by i, j, k
def qmat(q):  # left multiplication by quaternion q = (a, b, c, d) on (1, i, j, k)
    a, b_, c, d = q
    return np.array([[a, -b_, -c, -d], [b_, a, -d, c], [c, d, a, -b_], [d, -c, b_, a]])
def rmat(q):
    a, b_, c, d = q
    return np.array([[a, -b_, -c, -d], [b_, a, d, -c], [c, -d, a, b_], [d, c, -b_, a]])
Ri, Rj, Rk = (np.kron(np.eye(2), rmat(e)) for e in [(0, 1, 0, 0), (0, 0, 1, 0), (0, 0, 0, 1)])   # on H^2 = R^8
def sym_commuting(gens, n):
    basis = []
    for i in range(n):
        for j in range(i, n):
            S = np.zeros((n, n)); S[i, j] = S[j, i] = 1; basis.append(S.ravel())
    Bs = np.array(basis).T
    cons = np.vstack([np.vstack([(np.kron(np.eye(n), G) - np.kron(G.T, np.eye(n))) @ Bs]) for G in gens])
    return null_space(cons).shape[1]
dimC = sym_commuting([Ri], 8); dimH = sym_commuting([Ri, Rj, Rk], 8)
# quaternionic theory with the common turn J = Ri still satisfies A1-A4
Hq_basis = null_space(np.vstack([(np.kron(np.eye(8), G) - np.kron(G.T, np.eye(8))) for G in [Ri, Rj, Rk]]))
ok = True
for _ in range(5):
    v = Hq_basis @ rng.normal(size=Hq_basis.shape[1]); H = v.reshape(8, 8); H = (H + H.T)/2
    A = Ri @ H; ok &= np.allclose(A.T, -A) and np.allclose(A @ A, -H @ H) and np.allclose(A @ H, H @ A)
check("C3 without A5: a theory with quaternionic observables and the common turn J = R_i meets A1-A4, but its observables (6 per "
      "quaternionic 2-level system) are fewer than the 16 symmetric operators commuting with J; A5 adds the rest and gives "
      "complex quantum theory", ok and dimC == 16 and dimH == 6, f"commute with J: {dimC}, with all of i, j, k: {dimH}")
# genuine (Adler) quaternionic QM: generators must be quaternion-linear (commute with R_i, R_j, R_k)
Hq_obs = []
for col in range(Hq_basis.shape[1]):
    Hm = Hq_basis[:, col].reshape(8, 8); Hm = (Hm + Hm.T)/2
    if np.linalg.norm(Hm) > 1e-9: Hq_obs.append(Hm)
cm_obs = null_space(np.vstack([(np.kron(np.eye(8), G) - np.kron(G.T, np.eye(8))) for G in Hq_obs]))      # commutant of all quaternionic observables
span_R = np.array([np.eye(8).ravel(), Ri.ravel(), Rj.ravel(), Rk.ravel()]).T
same = cm_obs.shape[1] == 4 and np.linalg.matrix_rank(np.hstack([cm_obs, span_R])) == 4
cm_both = null_space(np.vstack([(np.kron(np.eye(8), G) - np.kron(G.T, np.eye(8))) for G in Hq_obs + [Ri, Rj, Rk]]))
check("C3b genuine quaternionic QM (generators quaternion-linear, as in Adler) fails A4 already: the operators commuting with every "
      "quaternionic observable are spanned by 1, R_i, R_j, R_k, and of these only multiples of 1 are quaternion-linear, so no "
      "quaternion-linear J can implement the energy shift; the C3 theory meets A1-A4 only by giving up quaternion-linearity",
      same and cm_both.shape[1] == 1, f"commutant of observables: {cm_obs.shape[1]}; also quaternion-linear: {cm_both.shape[1]}")
check("T5 with A5, the observables are exactly the complex Hermitian matrices: (n/2)^2 of them on R^n",
      all(sym_commuting([Jstd(m)], 2*m) == m*m for m in [1, 2, 3]))

# ---------- 6. composites: booking the energy zero to either part forces the complex tensor product ----------
ok = True; details = []
for (m1, m2) in [(1, 1), (2, 2), (2, 3)]:
    JA_, JB_ = Jstd(m1), Jstd(m2); nA, nB = 2*m1, 2*m2
    JL = np.kron(JA_, np.eye(nB)); JR = np.kron(np.eye(nA), JB_)
    V = null_space(JL - JR)                                   # where "booked to A" and "booked to B" agree
    d = V.shape[1]; ok &= d == 2*m1*m2                        # real dim of C^m1 (x)_C C^m2
    JJ = V.T @ np.kron(JA_, JB_) @ V; ok &= np.allclose(JJ, -np.eye(d))
    Jv = V.T @ JL @ V
    # local tomography: products of local observables (commuting with the local J) span all joint observables on V
    def herm_basis(m):
        out = []
        for a_ in range(m):
            for b_ in range(a_, m):
                X = np.zeros((m, m), complex); X[a_, b_] = X[b_, a_] = 1; out.append(X)
                if a_ != b_: Y = np.zeros((m, m), complex); Y[a_, b_] = 1j; Y[b_, a_] = -1j; out.append(Y)
        return [realify(x) for x in out]
    prods = [(V.T @ np.kron(X, Y) @ V).ravel() for X in herm_basis(m1) for Y in herm_basis(m2)]
    rank = np.linalg.matrix_rank(np.array(prods)); target = (m1*m2)**2
    ok &= rank == target and sym_commuting([Jv], d) == target; details.append(f"{m1}x{m2}: dim {d}, span {rank}/{target}")
check("T6 composites: requiring the energy-zero turn to be the same whichever part it's booked to (J x 1 = 1 x J) cuts the "
      "real product down to the complex tensor product; there J x J = -1 is fixed, and products of local observables span "
      "every joint observable (local tomography holds); this is the 'balanced tensor product' of Bang-Cho-Baek 2026",
      ok, "; ".join(details))
JL = np.kron(Jstd(1), np.eye(2)); JR = np.kron(np.eye(2), Jstd(1))
check("C4 without A6: on the plain real product R^2 x R^2 the two bookings differ (J x 1 != 1 x J), and J x J is a free "
      "hidden parameter, the extra joint quarter turn of the e-i-pi page", not np.allclose(JL, JR) and
      np.allclose(np.kron(Jstd(1), Jstd(1)) @ np.kron(Jstd(1), Jstd(1)), np.eye(4)))

# ---------- 7. the dial picture: the chain of small steps at the dial J/hbar ----------
ok = True
for m in [1, 2]:
    H = realify(rand_herm(m)); J = Jstd(m); t = 1.3; N = 20000
    step = np.eye(2*m) + (t/N)*J @ H/hbar; U = np.linalg.matrix_power(step, N)
    ok &= np.allclose(U, expm(t*J @ H/hbar), atol=1e-3)
    ok &= abs(np.linalg.norm(step, 2) - 1) < (t/N)**2*np.linalg.norm(H, 2)**2   # each step keeps size to second order
    for _ in range(3):
        Hr = rand_sym(2*m); stepr = np.eye(2*m) + (t/N)*Hr
        ok &= abs(np.linalg.norm(stepr, 2) - 1) > (t/N)*0.01                    # a real (stretching) dial changes size at first order
check("T7 in the framework's terms: the chain of steps 1 + (t/N) J H/hbar keeps size to second order and builds exp(tJH/hbar); "
      "a real dial stretches at first order. The theorem says the dial direction J must be one and the same for every "
      "Hamiltonian and every part of the world", ok)
# ---------- 8. the recent real formulations, tested against A0-A6 ----------
# Conventions follow each paper as stated; dynamics is taken as the image of the complex one (none of the papers defines a Hamiltonian flow).
def gamma_HW(M):                     # Hoffreumon-Woods Eq. II.3: Gamma(H) = 1 (x) Re H + J (x) Im H on R^2 (x) R^d
    return np.kron(np.eye(2), M.real) + np.kron(J2, M.imag)
def J_HW(d): return np.kron(J2, np.eye(d))
ok = True
for d in [2, 3]:
    Jn = J_HW(d)
    ok &= sym_commuting([Jn], 2*d) == d*d                                   # SY_d = symmetric matrices commuting with J (Lemma B.6): dim d^2
    for _ in range(3):
        M = rand_herm(d); G = gamma_HW(M)
        ok &= np.allclose(G, G.T) and np.allclose(G @ Jn, Jn @ G)             # observables commute with J
        ok &= np.allclose(gamma_HW(1j*M), Jn @ G)                             # i -> J
        U = expm(-1j*M*0.7); ok &= np.allclose(gamma_HW(U), expm(-0.7*Jn @ G))  # dynamics: generator -J Gamma(H), so A4 holds
        ok &= np.allclose(gamma_HW(M + 2.5*np.eye(d)), G + 2.5*np.eye(2*d))     # a shift of H is a shift of Gamma(H)
check("R1 Hoffreumon-Woods single system: observables are exactly the symmetric matrices commuting with J = J (x) 1 (A5 holds), "
      "i maps to J, and the image of e^{-iHt} is exp(-t J Gamma(H)), so A1-A4 hold with one shared J", ok)

ok = True
for (d1, d2) in [(2, 2), (2, 3)]:
    J1, J2_ = J_HW(d1), J_HW(d2); K = np.kron(J1, J2_); Pi = (np.eye(4*d1*d2) - K)/2
    ok &= np.allclose(Pi @ Pi, Pi) and round(np.trace(Pi)) == 2*d1*d2       # projector onto the J x J = -1 space, dim 2 d1 d2
    for _ in range(3):
        M1, M2 = rand_herm(d1), rand_herm(d2); G1, G2 = gamma_HW(M1), gamma_HW(M2)
        box = (np.kron(G1, G2) - np.kron(gamma_HW(1j*M1), gamma_HW(1j*M2)))/2   # Eq. B.34
        ok &= np.allclose(box, Pi @ np.kron(G1, G2))                           # = projection onto our W
        W = null_space(np.kron(J1, np.eye(2*d2)) - np.kron(np.eye(2*d1), J2_))
        ev_W = np.sort(np.linalg.eigvalsh(W.T @ np.kron(G1, G2) @ W))
        ev_c = np.sort(np.repeat(np.linalg.eigvalsh(np.kron(M1, M2)), 2))
        ok &= np.allclose(ev_W, ev_c)                                          # on W it is M1 (x)_C M2, each level doubled
check("R2 Hoffreumon-Woods composition (Eq. B.34): Gamma(M1) [box] Gamma(M2) = Pi (Gamma(M1) (x) Gamma(M2)) with Pi = (1 - J (x) J)/2, "
      "the projector onto W = {J (x) 1 = 1 (x) J}; on W it is the complex tensor product, so A6 holds", ok)

ok = True
JF = J2; e = np.eye(2)
for d in [1, 2]:
    # Barrios Hita et al.: R^d (x) R^2_F per system; composite R^d (x) R^2_F (x) R^d (x) R^2_F
    Jloc = np.kron(np.eye(d), JF)
    K = np.kron(Jloc, Jloc)
    flag_kernel = [np.kron(e[0], e[0]) + np.kron(e[1], e[1]), np.kron(e[0], e[1]) - np.kron(e[1], e[0])]
    vecs = []
    for a in range(d):
        for b in range(d):
            for f in flag_kernel:            # |psi'>_AB (x) flag-kernel vector, reordered to A(x)F_A(x)B(x)F_B
                T = np.einsum("a,b,fg->afbg", e[:d][a] if d > 1 else np.ones(1), e[:d][b] if d > 1 else np.ones(1), f.reshape(2, 2))
                vecs.append(T.ravel())
    Kern = np.array(vecs).T
    plus = null_space(K - np.eye(4*d*d))
    ok &= np.linalg.matrix_rank(Kern) == 2*d*d == plus.shape[1]
    ok &= np.linalg.matrix_rank(np.hstack([Kern, plus])) == 2*d*d             # kernel = the J (x) J = +1 space
    # observables T(Pi) = Re (x) I + Im (x) J_F are exactly the symmetric matrices commuting with 1 (x) J_F
    ok &= sym_commuting([Jloc], 2*d) == d*d
check("R3 Barrios Hita et al.: observables Re(Pi) (x) I + Im(Pi) (x) J_F are the symmetric matrices commuting with J = 1 (x) J_F (A5), "
      "and the quotient kernel is exactly the J (x) J = +1 space, so the composite (the quotient) is W, the complex tensor product (A6)", ok)

ok = True
tau = J2
def gam_inv(L): return np.kron(np.eye(2), L.real) + np.kron(tau, L.imag)        # Maioli-Curado-Gazeau Eq. 4/20
for (d1, d2) in [(2, 2), (2, 3)]:
    for _ in range(3):
        A, B = rand_herm(d1), rand_herm(d2); XA, YA, XB, YB = A.real, A.imag, B.real, B.imag
        otK = np.kron(np.eye(2), np.kron(XA, XB) - np.kron(YA, YB)) + np.kron(tau, np.kron(XA, YB) + np.kron(YA, XB))  # Eq. 32
        ok &= np.allclose(otK, gam_inv(np.kron(A, B)))                         # = the realified complex tensor product
        ok &= np.allclose(otK @ np.kron(tau, np.eye(d1*d2)), np.kron(tau, np.eye(d1*d2)) @ otK)
check("R4 Maioli-Curado-Gazeau: their composition rule (Eq. 32) is the real form of the complex tensor product and commutes with "
      "J = tau (x) 1, so their theory too carries one shared J and composes as W", ok)

ok = True
# Why A4 is physical: in real QM with each dynamics carrying its own positive-frequency J, the dynamics depends on the energy zero.
H = np.diag([-1., -1., 2., 2.]); A_true = Jstd(2) @ H                          # a dynamics with one shared J
def polar_J(A): w, Vv = np.linalg.eig(-A @ A); S = np.real(Vv @ np.diag(np.sqrt(w)) @ np.linalg.inv(Vv)); return A @ np.linalg.inv(S), S
Jp, Hp = polar_J(A_true)
ok &= not np.allclose(Jp, Jstd(2))                                              # positive-frequency J flips the negative-energy plane
ok &= np.allclose(Hp, np.diag([1., 1., 2., 2.]))                               # and reads the energy as |E|: the zero is not a convention
check("R5 why A4 matters: reading each dynamics with its own positive-frequency J makes energies |E| and flips the turn of "
      "negative-energy levels, so the physics would depend on the energy zero; A4 rules this out", ok)

ok = True
# A2 is redundant: Planck-Einstein (A3) at H + cI, with the shift rule (A4), gives J^2 = -1, BJ + JB = -2H, B^2 = -H^2.
# The first two alone force [J, H] = 0; then the third forces [B, H] = 0, which is A2. A1 (antisymmetry) is not used.
for n in [2, 4, 6]:
    for _ in range(20):
        Q = rng.normal(size=(n, n)); J = Q @ Jstd(n//2) @ np.linalg.inv(Q)       # any J with J^2 = -1, not even antisymmetric
        B = rng.normal(size=(n, n)); H = -(B @ J + J @ B)/2                      # any B; H as the c^1 coefficient fixes it
        ok &= np.allclose(J @ J, -np.eye(n)) and np.allclose(J @ H, H @ J, atol=1e-8)
        ok &= np.allclose(B @ H - H @ B, -(B @ B @ J - J @ B @ B)/2, atol=1e-8)  # [B,H] = -[B^2, J]/2 for every B
# so if also B^2 = -H^2 (which then commutes with J), [B,H] = 0.
check("T3d energy conservation (A2) follows from A3 and A4: the c^2 and c^1 terms of (B + cJ)^2 = -(H + cI)^2 force [J,H] = 0, "
      "and then the c^0 term forces [B,H] = 0, for any J with J^2 = -1 and any B (Lean: EnergyShift.A2_redundant)", ok)

ok = True
# A4 from composites: booking a constant c to A or to B gives the same joint Hamiltonian, so with additive
# non-interacting generators and faithful local operators, A(H + c) - A(H) cannot depend on H.
nA = nB = 4; IA, IB = np.eye(nA), np.eye(nB)
for _ in range(10):
    HA, HB, c = rand_sym(nA), rand_sym(nB), rng.normal()
    ok &= np.allclose(np.kron(HA + c*IA, IB) + np.kron(IA, HB), np.kron(HA, IB) + np.kron(IA, HB + c*IB))
def own_turn(H):                                    # an H-dependent complex structure: pair H's eigenvectors in order
    w, V = np.linalg.eigh(H); K = np.zeros_like(H)
    for i in range(0, H.shape[0], 2): a, b = V[:, i], V[:, i+1]; K += np.outer(b, a) - np.outer(a, b)
    return K
JB = Jstd(nB//2); mism = []
for _ in range(5):                                  # shift rule A(H + c) = A(H) + c K_H with K_H depending on H
    HA = rand_sym(nA); c = 0.7
    gen_bookA = np.kron(own_turn(HA)*c, IB); gen_bookB = np.kron(IA, JB*c)  # change of joint generator, per booking
    mism.append(np.linalg.norm(gen_bookA - gen_bookB))
ok &= min(mism) > 1e-3                              # the two bookings disagree: the joint generator is not well defined
KAB = np.kron(Jstd(nA//2), Jstd(nB//2)); w_, V_ = np.linalg.eigh(KAB); Wb = V_[:, w_ < 0]   # balanced product W
for _ in range(20):
    X = rng.normal(size=(nA, nA))
    ok &= np.linalg.norm(np.kron(X, IB) @ Wb) > 1e-6 * np.linalg.norm(X)        # X (x) 1 = 0 on W only if X = 0
check("T6b the shift rule A4 follows from composites: booking a constant to A or to B gives the same joint Hamiltonian; "
      "an H-dependent shift turn would give two different joint generators for it; local operators act faithfully on the "
      "balanced product W, so A(H + c) - A(H) is the same for every H (Lean: Booking.shift_independent)", ok)

ok = True
# A1 from stability: without antisymmetry, A3 and A4 still allow B = JH + E with E^2 = 0, E anticommuting with J.
# Shifting the energy zero to the level lam leaves the generator E there, whose flow 1 + tE grows without bound.
J = Jstd(2); E = np.kron(np.array([[0., 1.], [0., 0.]]), np.diag([1., -1.]))   # J-antilinear and nilpotent
lam = 1.3; H = lam*np.eye(4); B = lam*J + E
ok &= np.allclose(E @ J, -J @ E) and np.allclose(E @ E, 0) and not np.allclose(E.T, -E)
ok &= all(np.allclose((B + c*J) @ (B + c*J), -(H + c*np.eye(4)) @ (H + c*np.eye(4))) for c in np.linspace(-3, 3, 13))
ok &= np.allclose(expm(100*(B - lam*J)), np.eye(4) + 100*E) and np.linalg.norm(np.eye(4) + 100*E, 2) > 99
def sym_commuting(Jm):                              # basis of the symmetric matrices S with SJ = JS (the observables, A5)
    n = Jm.shape[0]; basis = []
    for i in range(n):
        for j in range(i, n):
            S = np.zeros((n, n)); S[i, j] = S[j, i] = 1; basis.append(S)
    M = np.array([(S @ Jm - Jm @ S).ravel() for S in basis]).T
    u, sv, vt = np.linalg.svd(M); null = vt[np.sum(sv > 1e-9):]
    return [sum(c*S for c, S in zip(row, basis)) for row in null]
def commutant_dim(ops):                             # dimension of the real matrices commuting with every op
    n = ops[0].shape[0]; I = np.eye(n)
    M = np.vstack([np.kron(O, I) - np.kron(I, O.T) for O in ops])
    return n*n - np.linalg.matrix_rank(M, tol=1e-9)
ok &= commutant_dim(sym_commuting(Jstd(2))) == 2      # orthogonal J: observables pin everything down to a + bJ
P = rng.normal(size=(4, 4)); Jbad = P @ Jstd(2) @ np.linalg.inv(P)             # J^2 = -1 but not orthogonal
ok &= np.allclose(Jbad @ Jbad, -np.eye(4)) and commutant_dim(sym_commuting(Jbad)) > 2   # observables too few
check("T3e antisymmetry (A1) from stability: without A1, A3 and A4 allow B = JH + E with E^2 = 0, and shifting the "
      "energy zero to the level leaves the generator E, whose flow 1 + tE grows without bound; with E = 0, a J that is "
      "not orthogonal leaves too few observables to fix states up to a + bJ, so completeness up to phase forces J^T = -J "
      "(Lean: Stability.shift_to_level, nilpotent_flow, transpose_turn)", ok)

print(f"\n{sum(r[1] for r in results)}/{len(results)} pass")
