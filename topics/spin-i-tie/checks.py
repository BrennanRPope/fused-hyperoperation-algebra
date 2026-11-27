"""Checks for the spin / i tie (DRAFT, Nov 27 2026).

Run:  python3 checks.py > checks_output.txt
Writes results.json beside it. Needs numpy, scipy, sympy.

S1  the shift lemma: a constant shift c*G is invisible iff G commutes
S2  the law's twist: same dial <=> no twist (exact, sympy)
S3  along the law's twist line f = 1 - w: how many quarter turns keep a
    constant potential shift a convention (3, 1, 3, 0)
S4  a rotation of the round twist field is a constant gauge rotation, so
    angular momentum becomes the charge of that gauge rotation
S5  charge + monopole: field angular momentum eg/4pi; Saha-Wilson = Dirac
S6  quantum gravity edge: a constant energy is a phase for quantum
    mechanics and a curvature for gravity (same fused shift, k imaginary
    or real)
S7  a spin following a direction over the sphere carries a monopole
    (Chern number 1) that real amplitudes cannot carry
"""
import json
import numpy as np
import sympy as sp
from scipy.linalg import expm
from scipy import integrate

rng = np.random.default_rng(20261127)
results = {}
passed = 0
total = 0


def check(name, ok, detail):
    global passed, total
    total += 1
    passed += bool(ok)
    print(f"[{'PASS' if ok else 'FAIL'}] {name}: {detail}")
    results[name] = {"pass": bool(ok), "detail": detail}


# ---------------------------------------------------------------- S1
def realify(M):
    """complex n x n -> real 2n x 2n"""
    return np.block([[M.real, -M.imag], [M.imag, M.real]])


n = 3
Hc = rng.normal(size=(n, n)) + 1j * rng.normal(size=(n, n))
Hc = (Hc + Hc.conj().T) / 2
H = realify(-1j * Hc)  # real antisymmetric generator of the dynamics
J = realify(1j * np.eye(n))  # the shared quarter turn, J^2 = -1
Gnc = rng.normal(size=(2 * n, 2 * n))
Gnc = Gnc - Gnc.T  # a random rotation generator that does not commute
dev_central, dev_noncentral = [], []
for c in [0.3, 1.0, 5.0]:
    for t in [0.5, 2.0]:
        dev_central.append(np.abs(expm(t * (H + c * J)) - expm(t * H) @ expm(t * c * J)).max())
        dev_noncentral.append(np.abs(expm(t * (H + c * Gnc)) - expm(t * H) @ expm(t * c * Gnc)).max())
# observables of complex QM (realified Hermitian) : expectation values unchanged by the J shift
O = rng.normal(size=(n, n)) + 1j * rng.normal(size=(n, n))
O = realify((O + O.conj().T) / 2)
psi = rng.normal(size=2 * n)
psi /= np.linalg.norm(psi)
ev = []
for c in [0.0, 1.0, 7.0]:
    p = expm(1.3 * (H + c * J)) @ psi
    ev.append(p @ O @ p)
ev_spread = max(ev) - min(ev)


def commutant_dim(mats):
    d = mats[0].shape[0]
    rows = []
    for M in mats:
        rows.append(np.kron(np.eye(d), M) - np.kron(M.T, np.eye(d)))
    A = np.vstack(rows)
    s = np.linalg.svd(A, compute_uv=False)
    return int(np.sum(s < 1e-9 * s[0])), A


# real QM on R^4: every real symmetric matrix is an observable
basis_real = []
for i in range(4):
    for j in range(i, 4):
        E = np.zeros((4, 4))
        E[i, j] = E[j, i] = 1
        basis_real.append(E)
dim_real, _ = commutant_dim(basis_real)
# complex QM on C^2 realified: Hermitian observables
herm = []
for i in range(2):
    for j in range(i, 2):
        E = np.zeros((2, 2), complex)
        E[i, j] = E[j, i] = 1
        herm.append(realify(E))
        if i != j:
            E = np.zeros((2, 2), complex)
            E[i, j] = 1j
            E[j, i] = -1j
            herm.append(realify(E))
dim_cplx, _ = commutant_dim(herm)
J2 = realify(1j * np.eye(2))
check("S1",
      max(dev_central) < 1e-10 and min(dev_noncentral) > 1e-2 and ev_spread < 1e-12
      and dim_real == 1 and dim_cplx == 2,
      f"shift by the shared J: max |e^(t(H+cJ)) - e^(tH)e^(tcJ)| = {max(dev_central):.1e}, "
      f"expectation spread {ev_spread:.1e}; by a non-commuting G: min deviation {min(dev_noncentral):.2f}; "
      f"commutant of all observables: real QM on R^4 dim {dim_real} (only 1, no quarter turn), "
      f"complex QM on C^2 dim {dim_cplx} (1 and J, J^2 = -1: {np.allclose(J2 @ J2, -np.eye(4))})")
results["S1"].update({"dev_central": float(max(dev_central)), "dev_noncentral_min": float(min(dev_noncentral)),
                      "dim_real": dim_real, "dim_complex": dim_cplx})

# ---------------------------------------------------------------- S2
a, b, k, k1, k2, x = sp.symbols("a b k k1 k2 x")
fa = lambda X, kk: (1 + kk * a) * X + a
fb = lambda X, kk: (1 + kk * b) * X + b
comp_same = sp.expand(fa(fb(x, k), k) - fb(fa(x, k), k))
comp_diff = sp.expand(fa(fb(x, k2), k1) - fb(fa(x, k1), k2))
law = sp.expand(fa(fb(x, k), k).subs(x, 0))
mult = sp.expand((1 + k * a) * (1 + k * b) - (1 + k * law))
check("S2", comp_same == 0 and sp.simplify(comp_diff - (k1 - k2) * a * b) == 0 and law == sp.expand(a + b + k * a * b) and mult == 0,
      f"same dial: twist {comp_same}; dials k1, k2: twist {sp.factor(comp_diff)}; "
      f"shift part of composition = {law} (the law); (1+ka)(1+kb) = 1 + k(a o_k b): {mult == 0}")

# ---------------------------------------------------------------- S3
sig = [np.array([[0, 1], [1, 0]], complex), np.array([[0, -1j], [1j, 0]]), np.array([[1, 0], [0, -1]], complex)]
T = [-0.5j * s for s in sig]  # anti-Hermitian su(2), [T_a, T_b] = eps_abc T_c
eps = np.zeros((3, 3, 3))
for (i, j, kk), s in {(0, 1, 2): 1, (1, 2, 0): 1, (2, 0, 1): 1, (0, 2, 1): -1, (2, 1, 0): -1, (1, 0, 2): -1}.items():
    eps[i, j, kk] = s
assert np.allclose(T[0] @ T[1] - T[1] @ T[0], T[2])


def A_field(X, wfun, sgn):
    r = np.linalg.norm(X)
    w = wfun(r)
    out = []
    for i in range(3):
        Ai = np.zeros((2, 2), complex)
        for aa in range(3):
            for j in range(3):
                Ai += sgn * eps[aa, i, j] * X[j] * (1 - w) / r ** 2 * T[aa]
        out.append(Ai)
    return out


def F_field(X, wfun, sgn, h=1e-5):
    A0 = A_field(X, wfun, sgn)
    dA = []
    for j in range(3):
        e = np.zeros(3)
        e[j] = h
        Ap, Am = A_field(X + e, wfun, sgn), A_field(X - e, wfun, sgn)
        dA.append([(Ap[i] - Am[i]) / (2 * h) for i in range(3)])  # dA[j][i] = d_j A_i
    F = {}
    for i in range(3):
        for j in range(3):
            F[i, j] = dA[i][j] - dA[j][i] + A0[i] @ A0[j] - A0[j] @ A0[i]
    return F


Xp = np.array([0.37, -0.81, 0.52])
# fix the sign convention so that w = -1 (f = 2) is pure gauge
sgn = None
for s_try in (+1, -1):
    Fm = F_field(Xp, lambda r: -1.0, s_try)
    if max(np.abs(M).max() for M in Fm.values()) < 1e-6:
        sgn = s_try
assert sgn is not None


def adj_commutant_dim(mats):
    """dimension of {phi in su(2): [M, phi] = 0 for all M}"""
    rows = []
    for M in mats:
        rows.append(np.array([[np.vdot(T[c], M @ T[bb] - T[bb] @ M).real * 4 for bb in range(3)] for c in range(3)]))
        rows.append(np.array([[np.vdot(T[c], M @ T[bb] - T[bb] @ M).imag * 4 for bb in range(3)] for c in range(3)]))
    s = np.linalg.svd(np.vstack(rows), compute_uv=False)
    return int(np.sum(s < 1e-6 * max(s[0], 1e-30))) if s[0] > 1e-9 else 3


def cov_const_dim(X, wfun, h=1e-5):
    """dimension of covariantly constant adjoint fields allowed at X: the common commutant of
    the curvature F_ij(X) and of its covariant derivatives D_k F_ij (holonomy algebra, first orders)"""
    F = F_field(X, wfun, sgn)
    A0 = A_field(X, wfun, sgn)
    mats = [F[i, j] for i in range(3) for j in range(i + 1, 3)]
    for kk in range(3):
        e = np.zeros(3)
        e[kk] = h
        Fp, Fm = F_field(X + e, wfun, sgn), F_field(X - e, wfun, sgn)
        for i in range(3):
            for j in range(i + 1, 3):
                mats.append((Fp[i, j] - Fm[i, j]) / (2 * h) + A0[kk] @ F[i, j] - F[i, j] @ A0[kk])
    if max(np.abs(M).max() for M in mats) < 1e-6:
        return 3  # flat: every constant shift is a convention (vacuum)
    return adj_commutant_dim(mats)


profiles = {
    "f=0 (w=1, vacuum)": lambda r: 1.0,
    "f=1 (w=0, the law)": lambda r: 0.0,
    "f=2 (w=-1, pure gauge)": lambda r: -1.0,
    "f=0.6 (w=0.4, frozen)": lambda r: 0.4,
    "bubble profile w=(1-r^2)/(1+r^2), at r=0.8": lambda r: (1 - r ** 2) / (1 + r ** 2),
}
Xs = 0.8 * Xp / np.linalg.norm(Xp)
counts = {name: cov_const_dim(Xs, wf) for name, wf in profiles.items()}


# covariant derivative of the hedgehog generator xhat.T, which is the one survivor at f = 1
def D_hedgehog(X, w, h=1e-6):
    def phi(Y):
        yh = Y / np.linalg.norm(Y)
        return sum(yh[aa] * T[aa] for aa in range(3))
    A0 = A_field(X, lambda r: w, sgn)
    tot = 0
    for i in range(3):
        e = np.zeros(3)
        e[i] = h
        Di = (phi(X + e) - phi(X - e)) / (2 * h) + A0[i] @ phi(X) - phi(X) @ A0[i]
        tot += np.abs(Di).max()
    return tot


dh = {w: D_hedgehog(Xs, w) for w in [1.0, 0.5, 0.0, -0.5, -1.0]}
ratio = [dh[w] / abs(w) for w in [1.0, 0.5, -0.5, -1.0]]
xh = Xs / np.linalg.norm(Xs)
Jx = sum(xh[aa] * 2 * T[aa] for aa in range(3))  # quarter turn: (xhat . (-i sigma))
check("S3",
      counts == {"f=0 (w=1, vacuum)": 3, "f=1 (w=0, the law)": 1, "f=2 (w=-1, pure gauge)": 3,
                 "f=0.6 (w=0.4, frozen)": 0, "bubble profile w=(1-r^2)/(1+r^2), at r=0.8": 0}
      and dh[0.0] < 1e-8 and np.ptp(ratio) < 1e-6 and np.allclose(Jx @ Jx, -np.eye(2)),
      f"sign convention A_i = {sgn:+d} eps_aij x_j (1-w)/r^2 T_a (w=-1 flat); "
      f"covariantly constant generators allowed: {counts}; "
      f"|D(xhat.T)| = {dh[0.0]:.1e} at f=1, proportional to |w| elsewhere (ratio spread {np.ptp(ratio):.1e}); "
      f"the survivor 2 xhat.T squares to -1")
results["S3"].update({"counts": counts, "D_hedgehog": {str(kk): float(v) for kk, v in dh.items()}})

# ---------------------------------------------------------------- S4
def lie_xi_A(X, wfun, h=1e-6):
    xi = lambda Y: np.array([-Y[1], Y[0], 0.0])
    dxi = np.array([[0, -1, 0], [1, 0, 0], [0, 0, 0.0]])  # dxi[j][i] = d_i xi^j
    A0 = A_field(X, wfun, sgn)
    out = []
    v = xi(X)
    for i in range(3):
        acc = np.zeros((2, 2), complex)
        for j in range(3):
            e = np.zeros(3)
            e[j] = h
            acc += v[j] * (A_field(X + e, wfun, sgn)[i] - A_field(X - e, wfun, sgn)[i]) / (2 * h)
            acc += A0[j] * dxi[j][i]
        out.append(acc)
    return out


res_plus, res_minus, size = [], [], []
for _ in range(20):
    X = rng.normal(size=3)
    wf = lambda r: np.cos(1.7 * r) * np.exp(-0.2 * r)
    L = lie_xi_A(X, wf)
    A0 = A_field(X, wf, sgn)
    comm = [A0[i] @ T[2] - T[2] @ A0[i] for i in range(3)]
    res_plus.append(max(np.abs(L[i] - comm[i]).max() for i in range(3)))
    res_minus.append(max(np.abs(L[i] + comm[i]).max() for i in range(3)))
    size.append(max(np.abs(L[i]).max() for i in range(3)))
best = min(max(res_plus), max(res_minus))
W_sign = "+" if max(res_plus) < max(res_minus) else "-"
# ad-invariance used to turn J into a surface term: <E,[A,W]> = -<[A,E],W>
ip = lambda X, Y: -2 * np.trace(X @ Y).real
adinv = 0.0
for _ in range(50):
    Am, Em, Wm = [sum(rng.normal() * T[aa] for aa in range(3)) for _ in range(3)]
    adinv = max(adinv, abs(ip(Em, Am @ Wm - Wm @ Am) + ip(Am @ Em - Em @ Am, Wm)))
check("S4", best < 1e-7 and min(size) > 1e-3 and adinv < 1e-12,
      f"rotation about z of the round twist field (any profile w(r)) = constant gauge rotation W = {W_sign}T_3: "
      f"residual {best:.1e} against |L_xi A| >= {min(size):.2f}; ad-invariance {adinv:.1e}. "
      f"So J_z = integral <E, D W> = surface flux <W, E> (Gauss law D.E = 0): angular momentum is the electric "
      f"charge of the compensating U(1), and a static magnetic bubble (E = 0) carries none")

# ---------------------------------------------------------------- S5
def dyon_J(e, g, d):
    # charge e at origin, monopole g at d zhat; Heaviside-Lorentz, c = 1
    def inner(r):
        def f(th):
            R = np.array([r * np.sin(th), 0.0, r * np.cos(th)])
            Rm = R - np.array([0, 0, d])
            B = g * Rm / (4 * np.pi * np.linalg.norm(Rm) ** 3)
            rh = R / r
            return np.sin(th) * (np.cos(th) * (rh @ B) - B[2])
        val, _ = integrate.quad(f, 0, np.pi, limit=400, points=[1e-6] if abs(r - d) < 1e-9 else None)
        return r * val
    tot = 0.0
    for lo, hi in [(0, d * 0.999), (d * 0.999, d * 1.001), (d * 1.001, 50 * d), (50 * d, 4000 * d)]:
        v, _ = integrate.quad(inner, lo, hi, limit=400)
        tot += v
    return (e / 2) * tot


e, g = 1.3, 0.7
Js = {d: dyon_J(e, g, d) for d in [0.5, 1.0, 2.0]}
target = e * g / (4 * np.pi)
rel = max(abs(abs(v) - target) / target for v in Js.values())
# Saha-Wilson: |J| = n hbar/2  vs Dirac: e g / hbar = 2 pi n  (hbar = 1)
saha = sorted({round(4 * np.pi * (nn / 2), 12) for nn in range(1, 6)})
dirac = sorted({round(2 * np.pi * nn, 12) for nn in range(1, 6)})
check("S5", rel < 2e-3 and saha == dirac,
      f"field angular momentum of charge e={e} and monopole g={g}: {', '.join(f'd={d}: {v:+.6f}' for d, v in Js.items())}, "
      f"|J| = eg/4pi = {target:.6f} within {rel:.1e}, independent of separation; allowed eg from J = n/2 (Saha-Wilson): "
      f"{[round(s, 4) for s in saha]}; from a single-valued phase e^(ieg) (Dirac): {[round(s, 4) for s in dirac]}: identical")
results["S5"].update({"J": {str(d): v for d, v in Js.items()}, "target": target})

# ---------------------------------------------------------------- S6
Hq = rng.normal(size=(4, 4)) + 1j * rng.normal(size=(4, 4))
Hq = (Hq + Hq.conj().T) / 2
P = lambda c, t: np.abs(expm(-1j * (Hq + c * np.eye(4)) * t)) ** 2
dP = max(np.abs(P(c, 1.7) - P(0, 1.7)).max() for c in [1.0, 50.0, -300.0])
# gravity: flat FRW with vacuum energy rho (G = 1): H = sqrt(8 pi rho / 3), a = e^{H t}
rhos = [0.01, 0.04, 0.16]
a_t = {r_: float(np.exp(np.sqrt(8 * np.pi * r_ / 3) * 2.0)) for r_ in rhos}
# the same fused shift (1 + k tau)^N, k real (gravity) or imaginary (quantum), t = 2
N = 200000
tt = 2.0
law_real = {r_: float((1 + np.sqrt(8 * np.pi * r_ / 3) * tt / N) ** N) for r_ in rhos}
law_imag = {r_: float(abs((1 - 1j * r_ * tt / N) ** N)) for r_ in rhos}
ok6 = dP < 1e-10 and all(abs(law_real[r_] / a_t[r_] - 1) < 1e-4 for r_ in rhos) \
    and all(abs(law_imag[r_] - 1) < 1e-6 for r_ in rhos) and a_t[0.16] > 2 * a_t[0.01]
check("S6", ok6,
      f"quantum: transition probabilities change by {dP:.1e} when every energy shifts by c (c up to 300); "
      f"gravity: the same vacuum energy rho sets the expansion, a(t=2) = "
      f"{', '.join(f'{v:.3f}' for v in a_t.values())} for rho = {rhos}; the fused shift (1+k tau)^N reproduces both: "
      f"k = sqrt(8 pi rho/3) real gives {', '.join(f'{v:.3f}' for v in law_real.values())}, "
      f"k = -i rho gives modulus {', '.join(f'{v:.7f}' for v in law_imag.values())}")
results["S6"].update({"a_t": a_t, "law_real": law_real, "law_imag": law_imag})

# ---------------------------------------------------------------- S7
def chern_lattice(hfun, Nt=60, Np=120):
    th = np.linspace(1e-4, np.pi - 1e-4, Nt)
    ph = np.linspace(0, 2 * np.pi, Np, endpoint=False)
    U = np.empty((Nt, Np, 2), complex)
    gap = np.inf
    for i, t_ in enumerate(th):
        for j, p_ in enumerate(ph):
            vals, vecs = np.linalg.eigh(hfun(t_, p_))
            gap = min(gap, vals[1] - vals[0])
            U[i, j] = vecs[:, 0]
    tot = 0.0
    for i in range(Nt - 1):
        for j in range(Np):
            j2 = (j + 1) % Np
            l1 = np.vdot(U[i, j], U[i + 1, j]); l2 = np.vdot(U[i + 1, j], U[i + 1, j2])
            l3 = np.vdot(U[i + 1, j2], U[i, j2]); l4 = np.vdot(U[i, j2], U[i, j])
            tot += np.angle(l1 * l2 * l3 * l4)
    return tot / (2 * np.pi), gap


h_cplx = lambda t_, p_: np.sin(t_) * np.cos(p_) * sig[0] + np.sin(t_) * np.sin(p_) * sig[1] + np.cos(t_) * sig[2]
h_real = lambda t_, p_: np.sin(t_) * np.cos(p_) * sig[0] + np.cos(t_) * sig[2]  # real symmetric: no sigma_y
C_c, gap_c = chern_lattice(h_cplx)
# the real family must close its gap: it has only two directions (sigma_x, sigma_z) to point in
gap_r = 2 * np.sqrt(np.sin(np.pi / 2) ** 2 * np.cos(np.pi / 2) ** 2 + np.cos(np.pi / 2) ** 2)
# a real family that stays gapped has Chern number 0 (real eigenvectors carry no Berry curvature)
h_real_gapped = lambda t_, p_: sig[2] + 0.5 * np.sin(t_) * np.cos(p_) * sig[0]
C_r, gap_rg = chern_lattice(h_real_gapped)
# Berry phase of real eigenvectors round any loop is 0 or pi (pi when the loop circles the closing point)
loop_phases = []
for (t0, p0, rad) in [(0.6, 0.0, 0.3), (2.4, 3.0, 0.3), (np.pi / 2, np.pi / 2, 0.2)]:
    ang = np.linspace(0, 2 * np.pi, 600, endpoint=False)
    vs = [np.linalg.eigh(h_real(t0 + rad * np.cos(u), p0 + rad * np.sin(u)))[1][:, 0] for u in ang]
    prod = np.prod([np.vdot(vs[i], vs[(i + 1) % len(vs)]) for i in range(len(vs))])
    loop_phases.append(float(np.angle(prod)))
z2 = all(min(abs(p_), abs(abs(p_) - np.pi)) < 1e-9 for p_ in loop_phases)
check("S7", abs(abs(C_c) - 1) < 1e-6 and gap_c > 1.99 and gap_r < 1e-12 and abs(C_r) < 1e-6 and gap_rg > 0.9 and z2,
      f"a spin-1/2 that follows a direction over the whole sphere: Chern number {C_c:+.6f} (a unit monopole), gap {gap_c:.3f}; "
      f"with real amplitudes only (no sigma_y) a direction-following family closes its gap (gap {gap_r:.1e} at theta = phi = pi/2), "
      f"a real family that stays gapped has Chern number {C_r:+.6f}, and loop phases are only 0 or pi: "
      f"{[round(p_, 6) for p_ in loop_phases]}")
results["S7"].update({"chern_complex": float(C_c), "gap_real": float(gap_r)})

print(f"\n{passed}/{total} pass")
results["summary"] = f"{passed}/{total} pass"
with open("results.json", "w") as fh:
    json.dump(results, fh, indent=1, default=float)
