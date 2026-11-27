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

# ---------- 3f. Planck-Einstein from "the rate depends only on the energy" ----------
ok = True
Z, Jg = rand_skew(4), rand_skew(4)                  # any A(0) = Z and any A(1) = J
sq = lambda x: (Z + x*Jg) @ (Z + x*Jg)
ok &= np.allclose(sq(3) - 3*sq(2) + 3*sq(1) - sq(0), 0)          # third difference vanishes: (Z + xJ)^2 is quadratic in x
al, be, ga = 2.25, -1.5, 0.25 + 0.0                  # g(E)^2 = al E^2 + be E + ga with ga = al x0^2 (rest term m = 0)
s = np.sqrt(al); x0 = -be/(2*al); ok &= np.isclose(ga, al*x0**2)
Hc = rand_herm(3); Hr = realify(Hc); Jr = Jstd(3); I6 = np.eye(6)
Ag = lambda Hm: s*Jr @ (Hm - x0*np.eye(Hm.shape[0]))                 # the theorem's answer A(H) = J (H - E0) / hbar, hbar = 1/s
g2 = lambda Hm: al*Hm @ Hm + be*Hm + ga*np.eye(Hm.shape[0])
ok &= all(np.allclose(Ag(Hr + c*I6) @ Ag(Hr + c*I6), -g2(Hr + c*I6)) for c in np.linspace(-2, 2, 9))
ok &= np.allclose(Ag(Hr + 0.7*I6), Ag(Hr) + 0.7*s*Jr)                 # shift rule, with A(1) = s J
m = 0.8; E1, E2 = 0.5, -1.2                           # the quaternion-like way out, m != 0, on R^8 = H^2
H8 = np.diag([E1]*4 + [E2]*4); B8 = s*Ri @ H8 + m*Rj; I8 = np.eye(8)
ok &= np.allclose(Rj @ Ri, -Ri @ Rj) and np.allclose(B8.T, -B8) and np.allclose(B8 @ H8, H8 @ B8)
ok &= all(np.allclose((B8 + c*s*Ri) @ (B8 + c*s*Ri), -(s*s*(H8 + c*I8) @ (H8 + c*I8) + m*m*I8)) for c in np.linspace(-2, 2, 9))
fr = np.sort(np.abs(np.linalg.eigvals(B8).imag)); want = np.sort([np.sqrt(s*s*E**2 + m*m) for E in [E1]*4 + [E2]*4])
ok &= np.allclose(fr, want) and not np.allclose(fr, np.sort([s*abs(E) for E in [E1]*4 + [E2]*4]))   # rate sqrt(s^2 E^2 + m^2), not s|E|
def anti_dim(Jm):                                     # dimension of the antisymmetric M with M J + J M = 0
    n = Jm.shape[0]; basis = []
    for i in range(n):
        for j in range(i + 1, n):
            X = np.zeros((n, n)); X[i, j], X[j, i] = 1, -1; basis.append(X)
    Mx = np.array([(X @ Jm + Jm @ X).ravel() for X in basis]).T
    return len(basis) - np.linalg.matrix_rank(Mx, tol=1e-9)
ok &= anti_dim(J2) == 0 and anti_dim(Jstd(2)) == 2    # a level of real dimension 2 has no room for M; dimension 4 does
check("T3f Planck-Einstein from 'the rate depends only on the energy': with the shift rule, A(x 1)^2 is quadratic in x, so "
      "g(E)^2 = alpha (E - E0)^2 + m^2; A(H) = J (H - E0) sqrt(alpha) solves it with hbar = 1/sqrt(alpha); m != 0 needs a second "
      "quarter turn anticommuting with J (R_j beside R_i on R^8, rate sqrt(alpha E^2 + m^2)), which a level of real dimension 2 "
      "cannot hold (Lean: PlanckEinstein.rate_sq_quadratic, M_anticomm_J, M_sq, planck_einstein)", ok)

# ---------- 3g. composites rule out the gap; C5 the energy-only rate rule is still needed ----------
ok = True
def nullsp(Mx, tol=1e-9):
    u, sv, vt = np.linalg.svd(Mx); return vt[np.sum(sv > tol):].T
Ri4, Rj4, I4q = rmat((0, 1, 0, 0)), rmat((0, 0, 1, 0)), np.eye(4)
PA, QB = np.kron(Ri4, I4q), np.kron(I4q, Ri4)                     # J_A x 1 and 1 x J_B on R^4 x R^4
Wb, Wo = nullsp(PA - QB), nullsp(PA + QB)                          # balanced product W and its partner W-
ok &= Wb.shape[1] == 8 and Wo.shape[1] == 8
ok &= np.allclose(Wo @ Wo.T @ np.kron(Rj4, I4q) @ Wb, np.kron(Rj4, I4q) @ Wb)   # a gap term sends W into W-
ok &= np.allclose(Wb @ Wb.T @ PA @ Wb, PA @ Wb)                                  # the J(H - E0) part keeps W
anti = []
for i in range(4):
    for j in range(i + 1, 4):
        X = np.zeros((4, 4)); X[i, j], X[j, i] = 1, -1; anti.append(X)
co = nullsp(np.array([(X @ Ri4 + Ri4 @ X).ravel() for X in anti]).T)
gaps = [sum(c*X for c, X in zip(col, anti)) for col in co.T]        # all antisymmetric M anticommuting with J (R_j, R_k)
Mx = np.array([(np.kron(X, I4q) @ Wb).ravel() for X in gaps] + [(np.kron(I4q, X) @ Wb).ravel() for X in gaps]).T
ok &= len(gaps) == 2 and Mx.shape[1] - np.linalg.matrix_rank(Mx, tol=1e-9) == 0   # only M_A = M_B = 0 keeps the joint generator on W
ok &= not np.allclose(np.kron(Ri4/1.0, I4q) @ Wb, np.kron(I4q, Ri4/2.0) @ Wb)       # booking a constant forces hbar_A = hbar_B
check("T3g composites rule out the gap: on R^4 x R^4 with J = R_i, each gap term M_A x 1, 1 x M_B sends the balanced product W "
      "into W-, so the joint generator stays on W only if M_A = M_B = 0; booking a constant also forces one hbar for every system "
      "(Lean: NoGap.gap_to_odd, gap_sum_zero, gap_zero)", ok)

def herm_eigs(eigs):
    n = len(eigs); X = rng.normal(size=(n, n)) + 1j*rng.normal(size=(n, n)); U, _ = np.linalg.qr(X)
    return U @ np.diag(eigs) @ U.conj().T
def pinch(Hm, Qm, tol=1e-8):                       # keep only the parts of Q inside each energy eigenspace of H
    w, V = np.linalg.eigh(Hm); out = np.zeros_like(Qm); i = 0
    while i < len(w):
        j = i
        while j < len(w) and abs(w[j] - w[i]) < tol: j += 1
        P = V[:, i:j] @ V[:, i:j].conj().T; out += P @ Qm @ P; i = j
    return out
kap, ok = 0.4, True
QA, QB2 = np.diag([0., 1, 2]).astype(complex), np.diag([0., 1]).astype(complex)
clock1 = lambda Hm, Qm: -1j*(Hm + kap*pinch(Hm, Qm))                       # rate = energy + kap x conserved charge
clock2 = lambda Hm: -1j*((1 + kap)*Hm - kap*np.trace(Hm).real/Hm.shape[0]*np.eye(Hm.shape[0]))   # offset by the mean energy
for HA, HB in [(herm_eigs([0., 1, 1]), herm_eigs([0., 1])), (herm_eigs(rng.normal(size=3)), herm_eigs(rng.normal(size=2)))]:
    Ht = np.kron(HA, np.eye(2)) + np.kron(np.eye(3), HB); Qt = np.kron(QA, np.eye(2)) + np.kron(np.eye(3), QB2)
    for G, GA, GB, sh in [(clock1(Ht, Qt), clock1(HA, QA), clock1(HB, QB2), clock1(HA + 0.7*np.eye(3), QA)),
                          (clock2(Ht), clock2(HA), clock2(HB), clock2(HA + 0.7*np.eye(3)))]:
        ok &= np.allclose(G, np.kron(GA, np.eye(2)) + np.kron(np.eye(3), GB))          # composites (C)
        ok &= np.allclose(G.conj().T, -G) and np.allclose(G @ Ht, Ht @ G)              # A1, A2
        ok &= np.allclose(sh, GA - 0.7j*np.eye(3))                                      # A4
r1 = sorted(np.abs(np.linalg.eigvals(clock1(np.diag([1., 1, 0]).astype(complex), QA))))
ra = abs(np.linalg.eigvals(clock2(np.diag([1., 0, 0])))[0]); rb = abs(np.linalg.eigvals(clock2(np.diag([1., 2, 2])))[0])
ok &= np.allclose(r1[1:], [1.0, 1.4]) and not np.isclose(ra, rb)    # two E = 1 states turn at different rates
check("C5 composites do not imply (A3'): a clock that reads energy plus a multiple of a conserved charge, and one offset by the "
      "Hamiltonian's mean energy, both satisfy A1, A2, A4 and composites, yet give states of equal energy different rates", ok)

# ---------- 8. (A0'): the energy-zero shift changes no measured average ----------
ok = True
J6 = Jstd(3); H6 = realify(rand_herm(3)); I6b = np.eye(6)
rot = lambda t: np.cos(t)*I6b + np.sin(t)*J6                      # the shift's common rotation exp(tJ)
A0g = lambda Hm: J6 @ Hm / hbar                                  # complex dynamics, Hamiltonians restricted to those commuting with J
ok &= np.allclose(A0g(H6).T, -A0g(H6)) and np.allclose(A0g(H6) @ H6, H6 @ A0g(H6))
ok &= np.allclose(A0g(H6) @ A0g(H6), -H6 @ H6/hbar**2) and np.allclose(A0g(H6 + 0.3*I6b), A0g(H6) + 0.3*A0g(I6b))   # A1-A4 hold
Obad = rand_sym(6); psi6 = rng.normal(size=6); psi6 /= np.linalg.norm(psi6)    # a real-QM observable not commuting with J
avg = [(rot(t) @ psi6) @ Obad @ (rot(t) @ psi6) for t in np.linspace(0, np.pi, 9)]
ok &= not np.allclose(J6 @ Obad, Obad @ J6) and np.ptp(avg) > 1e-3        # without (A0) the shift is visible
Ogood = realify(rand_herm(3))
ok &= all(np.allclose(rot(t).T @ Ogood @ rot(t), Ogood) for t in np.linspace(0, 2*np.pi, 9))   # commuting observables: unseen at all t
ok &= np.allclose(J6.T @ Ogood @ J6, Ogood) and not np.allclose(J6.T @ Obad @ J6, Obad)       # the quarter turn alone tells them apart
check("T8 (A0') replaces (A0): real QM with all symmetric observables but only complex Hamiltonians meets A1-A4, yet a shift of "
      "the energy zero makes the average of an observable not commuting with J oscillate; invariance of every average under the "
      "shift's rotation, already under its quarter turn J, is equivalent to [J, O] = 0 (Lean: Unseen.comm_of_quarter, invariant_of_comm)", ok)

# ---------- 9. What (A5) excludes: three smaller observable sets ----------
ok = True
Xp = np.array([[0, 1], [1, 0]], complex); Yp = np.array([[0, -1j], [1j, 0]]); Zp = np.diag([1., -1]).astype(complex); Ip = np.eye(2)
k3 = lambda a, b, c: np.kron(np.kron(a, b), c)
gam = [k3(Xp, Ip, Ip), k3(Yp, Ip, Ip), k3(Zp, Xp, Ip), k3(Zp, Yp, Ip), k3(Zp, Zp, Xp), k3(Zp, Zp, Yp)]
G6 = gam[0] @ gam[1] @ gam[2] @ gam[3] @ gam[4] @ gam[5]
gam.append(next(ph*G6 for ph in [1, 1j, -1, -1j] if np.allclose(ph*G6, (ph*G6).conj().T)))      # seventh gamma
ok &= all(np.allclose(gam[a] @ gam[b] + gam[b] @ gam[a], 2*(a == b)*np.eye(8)) for a in range(7) for b in range(7))
V7 = [np.eye(8, dtype=complex)] + gam                            # spin factor: 8 observables on C^8 (A5 asks for 64)
def in_span(Mx, basis):
    Bm = np.array([b.ravel() for b in basis]).T; x = np.linalg.lstsq(Bm, Mx.ravel(), rcond=None)[0]; return np.allclose(Bm @ x, Mx.ravel())
ok &= all(in_span((a @ b + b @ a)/2, V7) for a in V7 for b in V7)                    # Jordan-closed
J16 = Jstd(8); R7 = [realify(o) for o in V7]
ok &= commutant_dim(R7) == 2                                                      # (T) holds: only a 1 + b J
ok &= all(np.allclose((J16 @ H).T, -J16 @ H) and np.allclose((J16 @ H) @ (J16 @ H), -H @ H) for H in R7)   # A(H) = JH: A1, A3
ok &= len(sym_commuting(J16)) == 64                                               # (A5) would need all 64
prods = [np.kron(a, b) for a in V7 for b in V7]
jp = (np.kron(gam[0], gam[0]) @ np.kron(gam[1], gam[1]) + np.kron(gam[1], gam[1]) @ np.kron(gam[0], gam[0]))/2
ok &= np.allclose(jp, np.kron(gam[0] @ gam[1], gam[0] @ gam[1])) and not in_span(jp, prods)   # local tomography + Jordan closure fail
Pss = np.diag([1, 1, 0, 0]).astype(complex)                                       # superselection on C^4: two sectors
def blk(h, sector):                                                               # a 2x2 Hermitian h placed in one sector
    out = np.zeros((4, 4), complex); out[2*sector:2*sector+2, 2*sector:2*sector+2] = h; return out
ssel = [realify(blk(h, s_)) for s_ in (0, 1) for h in [np.diag([1, 0]), np.diag([0, 1]), Xp, Yp]]
ok &= all(np.allclose(realify(Pss) @ S, S @ realify(Pss)) for S in ssel)            # every observable commutes with P
ok &= commutant_dim(ssel) == 4                                                    # (T) fails: P and JP commute too
def qherm(a, d, q):                                                                # quaternionic Hermitian 2x2, left action on H^2 = R^8
    qc = (q[0], -q[1], -q[2], -q[3])
    return np.block([[a*np.eye(4), qmat(q)], [qmat(qc), d*np.eye(4)]])
qobs = [qherm(1, 0, (0, 0, 0, 0)), qherm(0, 1, (0, 0, 0, 0))] + [qherm(0, 0, e) for e in [(1, 0, 0, 0), (0, 1, 0, 0), (0, 0, 1, 0), (0, 0, 0, 1)]]
ok &= all(np.allclose(S, S.T) and np.allclose(S @ Ri, Ri @ S) for S in qobs)       # observables commute with J = R_i
ok &= commutant_dim(qobs) == 4                                                    # (T) fails: 1, R_i, R_j, R_k
check("T9 (A5) is independent: with A(H) = JH, superselection and quaternionic observables fail (T), while the spin factor "
      "span{1, gamma_1..gamma_7} on C^8 is Jordan-closed and passes (T) with 8 observables instead of 64; only locally "
      "tomographic composites with Jordan-closed joint observables exclude it", ok)

# ---------- 10. What (C) asks: independence does the work, faithfulness = no local state excluded ----------
ok = True
nA, nB = 4, 4
JA, JB = Jstd(2), Jstd(2)
K = np.kron(JA, JB); wv, wV = np.linalg.eigh(K); Wb = wV[:, wv < 0]               # balanced product: J_A x J_B = -1
def absH(Hm):                                                                     # |H| through the spectral decomposition
    w_, V_ = np.linalg.eigh(Hm); return V_ @ np.diag(np.abs(w_)) @ V_.T
HA = realify(herm_eigs([-1.0, 2.0])); HB = realify(herm_eigs([-0.5, 0.3]))
polA = lambda Hm: JA @ absH(Hm)/hbar                                              # the rule A(H) = J|H|/hbar
ok &= np.allclose(polA(HA).T, -polA(HA)) and np.allclose(polA(HA) @ polA(HA), -HA @ HA/hbar**2)   # A1, A3
ok &= not np.allclose(polA(HA + 1.5*np.eye(4)), polA(HA) + 1.5*JA/hbar)          # but not A4
Hj = np.kron(HA, np.eye(4)) + np.kron(np.eye(4), HB)
Jw = np.kron(JA, np.eye(4))
Gj = Jw @ absH(Hj)/hbar                                                           # joint rule, a function of the joint H alone
Pw = Wb @ Wb.T
ok &= np.allclose(Pw @ Gj @ Pw, Gj @ Pw)                                          # it keeps W
ok &= np.allclose((Gj @ Gj) @ Wb, -(Hj @ Hj) @ Wb/hbar**2)                       # A3 for the joint system on W
add = np.kron(polA(HA), np.eye(4)) + np.kron(np.eye(4), JB @ absH(HB)/hbar)
ok &= not np.allclose(Gj @ Wb, add @ Wb)                                          # only independence (additivity) fails
# faithfulness: X (x) 1 vanishes on W  <=>  X kills the support of every reduced state of A
def reduced_A(w): M = w.reshape(nA, nB); return M @ M.T
supp = np.hstack([reduced_A(Wb @ rng.normal(size=Wb.shape[1])) for _ in range(6)])
ok &= np.linalg.matrix_rank(supp, tol=1e-9) == nA                                # supports span: no state of A excluded
Wsmall = np.kron(np.eye(4)[:, :2], np.eye(4))                                     # a joint space that excludes half of A
Z = np.diag([0., 0, 1, 1])
ok &= np.allclose(np.kron(Z, np.eye(4)) @ Wsmall, 0)                             # faithfulness fails there,
supp2 = np.hstack([reduced_A(Wsmall @ rng.normal(size=Wsmall.shape[1])) for _ in range(6)])
ok &= np.linalg.matrix_rank(supp2, tol=1e-9) == 2 and np.allclose(Z @ supp2, 0)  # exactly because those states are excluded
check("T6c (C)'s work is done by independence: J|H|/hbar is a function of the joint Hamiltonian, keeps the balanced "
      "product, meets A1 and A3 and faithfulness, yet fails additivity and A4; faithfulness holds exactly when the "
      "supports of A's reduced states span A's state space (composition excludes no state of A)", ok)


# ---------- 11. (A4) from measured averages: (A0') with (T0) and continuity replaces (A4) ----------
ok = True
m3 = 3; J6c = Jstd(m3); I6c = np.eye(2*m3)
def ev(Xg, t): return expm_np(t*Xg)
def expm_np(Mx):                                                                  # matrix exponential by eigendecomposition
    w_, V_ = np.linalg.eig(Mx); return np.real(V_ @ np.diag(np.exp(w_)) @ np.linalg.inv(V_))
def avg_same(Ag, Hm, c, obs, trials=4):                                          # (A0'): shifting the energy zero changes no average
    good = True
    for _ in range(trials):
        psi = rng.normal(size=Hm.shape[0]); psi /= np.linalg.norm(psi); t = rng.normal()
        a = ev(Ag(Hm + c*np.eye(Hm.shape[0])), t) @ psi; b = ev(Ag(Hm), t) @ psi
        good &= all(np.isclose(a @ O @ a, b @ O @ b) for O in obs)
    return good
obs6 = [realify(rand_herm(m3)) for _ in range(4)]
Hc = realify(rand_herm(m3))
ok &= avg_same(lambda Hm: J6c @ Hm/hbar, Hc, 0.8, obs6)                          # complex QM meets (A0')
# the sign step: with B = JH, (B + beta J)^2 = -(H + c)^2 only for beta = c
Bc = J6c @ Hc; c0 = 0.8
ok &= np.allclose((Bc + c0*J6c) @ (Bc + c0*J6c), -(Hc + c0*I6c) @ (Hc + c0*I6c))
ok &= not np.allclose((Bc - c0*J6c) @ (Bc - c0*J6c), -(Hc + c0*I6c) @ (Hc + c0*I6c))
# without continuity: each class of Hamiltonians may run backwards, A(H) = s(H) J H / hbar
K0 = realify(np.diag([1., -1, 0]).astype(complex))                                # traceless, so s is constant on H + c 1
sgn = lambda Hm: 1.0 if np.trace(Hm @ K0) >= 0 else -1.0
Ab = lambda Hm: sgn(Hm)*J6c @ Hm/hbar
Hneg = realify(np.diag([-1., 2, 0.5]).astype(complex)); ok &= sgn(Hneg) < 0
ok &= np.allclose(Ab(Hneg).T, -Ab(Hneg)) and np.allclose(Ab(Hneg) @ Ab(Hneg), -Hneg @ Hneg/hbar**2)   # A1, A3
ok &= avg_same(Ab, Hneg, 0.8, obs6)                                               # (A0') and (T) hold
ok &= not np.allclose(Ab(Hneg + 0.8*I6c), Ab(Hneg) + 0.8*Ab(I6c))                 # (A4) fails; s jumps, so A is not continuous
# without (T): quaternionic observables, with an imaginary unit that turns continuously with H
Ri_, Rj_ = (np.kron(np.eye(2), rmat(e)) for e in [(0, 1, 0, 0), (0, 0, 1, 0)])
def qherm2(a, d, q):
    qc = (q[0], -q[1], -q[2], -q[3]); return np.block([[a*np.eye(4), qmat(q)], [qmat(qc), d*np.eye(4)]])
qobs2 = [qherm2(1, 0, (0, 0, 0, 0)), qherm2(0, 1, (0, 0, 0, 0))] + [qherm2(0, 0, e) for e in [(1, 0, 0, 0), (0, 1, 0, 0), (0, 0, 1, 0), (0, 0, 0, 1)]]
Kq = qherm2(1, -1, (0, 0, 0, 0))                                                  # traceless: theta is constant on H + c 1
Aq = lambda Hm: (np.cos(np.trace(Hm @ Kq)/8)*Ri_ + np.sin(np.trace(Hm @ Kq)/8)*Rj_) @ Hm/hbar
Hq = qherm2(2.0, -1.0, (0.3, 0.5, -0.2, 0.4))
ok &= np.allclose(Aq(Hq).T, -Aq(Hq)) and np.allclose(Aq(Hq) @ Aq(Hq), -Hq @ Hq/hbar**2)   # A1, A3
ok &= avg_same(Aq, Hq, 0.8, qobs2)                                                # (A0'), and A is continuous in H
ok &= not np.allclose(Aq(Hq + 0.8*np.eye(8)), Aq(Hq) + 0.8*Aq(np.eye(8)))          # (A4) fails, because (T0) fails
# (T0): the antisymmetric matrices commuting with every observable; dimension 0 (real), 1 (complex), 3 (quaternionic)
def skew_commutant_dim(obs):
    n = obs[0].shape[0]; basis = []
    for a in range(n):
        for b in range(a + 1, n):
            E = np.zeros((n, n)); E[a, b], E[b, a] = 1, -1; basis.append(E)
    Mlin = np.array([np.concatenate([(X @ O - O @ X).ravel() for O in obs]) for X in basis]).T
    return len(basis) - np.linalg.matrix_rank(Mlin, tol=1e-9)
ok &= skew_commutant_dim([rand_sym(4) for _ in range(4)]) == 0
ok &= skew_commutant_dim(obs6) == 1 and skew_commutant_dim(qobs2) == 3
# real theory, all symmetric observables: a shift that adds nothing to A(H) breaks (A3) at the shifted Hamiltonian
Hr = rand_sym(4); ok &= not np.allclose((Hr + 0.8*np.eye(4)) @ (Hr + 0.8*np.eye(4)), Hr @ Hr)
check("T10 (A4) from measured averages: complex QM meets (A0'); the antisymmetric commutant of the observables has dimension 0, 1, 3 for real, complex, quaternionic; given (A3), a shift that only adds a multiple of J must add cJ; "
      "dropping continuity lets a class of Hamiltonians run backwards, and dropping (T0) lets the imaginary unit turn with H, "
      "and each keeps (A1), (A3) and (A0') while breaking (A4); with no antisymmetric commutant (real theory) a shift must add nothing, which (A3) forbids", ok)

# ---------- 12. (A3') with measured averages: (A3) from E0, and no gap ----------
ok = True
E0 = 0.7
Aoff = lambda Hm: J6c @ (Hm - E0*np.eye(Hm.shape[0]))/hbar                       # energies measured from E0
ok &= np.allclose(Aoff(Hc) @ Aoff(Hc), -(Hc - E0*I6c) @ (Hc - E0*I6c)/hbar**2)   # (A3') with g(E)^2 = (E - E0)^2/hbar^2
ok &= avg_same(Aoff, Hc, 0.8, obs6) and skew_commutant_dim(obs6) == 1            # (A0') and (T0)
ok &= np.allclose(Aoff(Hc + 0.8*I6c), Aoff(Hc) + 0.8*J6c/hbar)                    # (A4_w) with K = J/hbar
ok &= not np.allclose(Aoff(Hc + 0.8*I6c), Aoff(Hc) + 0.8*Aoff(I6c))              # but not (A4), which forces A(0) = 0
# the gapped rule: J = R_i, H quaternion-linear, A(H) = (J(H - E0) + mu R_j)/hbar
mu = 0.5; I8 = np.eye(8)
Ag = lambda Hm: (Ri_ @ (Hm - E0*I8) + mu*Rj_)/hbar
Hg = qherm2(2.0, -1.0, (0.3, 0.5, -0.2, 0.4))
ok &= np.allclose(Ag(Hg).T, -Ag(Hg))
ok &= np.allclose(Ag(Hg) @ Ag(Hg), -((Hg - E0*I8) @ (Hg - E0*I8) + mu**2*I8)/hbar**2)   # (A3') with a gap
ok &= np.allclose(Ag(Hg + 0.8*I8), Ag(Hg) + 0.8*Ri_/hbar)                         # (A4_w) holds, (A4) does not
ok &= not np.allclose(Ag(Hg + 0.8*I8), Ag(Hg) + 0.8*Ag(I8))
def jcomm_obs(Jm, k=4):                                                           # symmetric matrices commuting with Jm
    out = []
    for _ in range(k):
        O = rand_sym(Jm.shape[0]); out.append((O + Jm @ O @ Jm.T)/2)
    return out
obsg = jcomm_obs(Ri_)
ok &= skew_commutant_dim(obsg) == 1                                               # (T0) holds
ok &= not avg_same(Ag, Hg, 0.8, obsg)                                             # but (A0') fails: the shift is visible
# why: the second derivative of V_t = e^{tX} e^{-tY} at 0 is [Y, D] + D^2, and here [Y, D] has a part outside the commutant
Yg = Ag(Hg); Dg = 0.8*Ri_/hbar; V2 = Yg @ Dg - Dg @ Yg + Dg @ Dg
ok &= np.allclose(V2, 0.8*mu*(Rj_ @ Ri_ - Ri_ @ Rj_)/hbar**2 + Dg @ Dg)
ok &= not all(np.allclose(V2 @ O, O @ V2) for O in obsg)
dt = 1e-4; Vt = lambda t: expm_np(t*Ag(Hg + 0.8*I8)) @ expm_np(-t*Yg)
ok &= np.allclose((Vt(dt) - 2*np.eye(8) + Vt(-dt))/dt**2, V2, atol=1e-4)          # finite-difference check of that formula
# without a gap the dynamics commutes with J, as the proof derives
ok &= np.allclose(Aoff(Hc) @ J6c, J6c @ Aoff(Hc))
check("T11 with the weaker (A3'): A(H) = J(H - E0)/hbar meets (A1), (A3'), (A0'), (T0) and (A4_w) but not (A4), which "
      "forces A(0) = 0; the gapped rule meets (A1), (A3'), (T0) and (A4_w) but fails (A0'), because the second derivative "
      "[Y, D] + D^2 of e^{tA(H+c)} e^{-tA(H)} leaves the commutant", ok)

print(f"\n{sum(r[1] for r in results)}/{len(results)} pass")
