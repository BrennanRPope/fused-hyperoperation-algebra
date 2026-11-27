"""DRAFT round 25 checks: the other forces' twist and gravity's tidal pull, together.

Brennan asked to merge two open threads because he expects them to depend on each other: (1) the weak and strong forces'
twist (rounds 4-5: a hidden 3-sphere for SU(2), SU(3)'s group space) put on the law's expanding 3-sphere slices
(rounds 20-22), and (2) whether a dial that varies from place to place, plus the twist, can make the tidal pull that round 23
found missing. Here: the law's own transport is the slice's bending, which in the weak force's language is the 'meron'
field; that field on the expanding slices; how its energy feeds the dial; and why only a dial that varies can make a tide.
Run: python3 checks.py  (numpy, scipy, sympy).  Writes results.json and checks_output.txt.
"""
import json, io, math, time
import numpy as np
import sympy as spy
from scipy.integrate import solve_ivp, cumulative_trapezoid

results, passed, failed = {}, 0, 0
_log = io.StringIO()
def say(s):
    print(s, flush=True); _log.write(s + "\n")
def check(name, ok, detail=""):
    global passed, failed
    passed, failed = (passed + 1, failed) if ok else (passed, failed + 1)
    say(("PASS " if ok else "FAIL ") + name + (("  | " + detail) if detail else ""))
T_START = time.time()

# ------------------------------------------------------------------ the law in 3D (rounds 12-13, 18, unchanged)
def oplus(a, b, qq=1.0):
    ab = (a * b).sum(-1, keepdims=True); aa = (a * a).sum(-1, keepdims=True); bb = (b * b).sum(-1, keepdims=True)
    return ((1 - 2 * qq * ab - qq * bb) * a + (1 + qq * aa) * b) / (1 - 2 * qq * ab + qq**2 * aa * bb)
def law_dist(a, b):
    return 2 * np.arctan(np.linalg.norm(oplus(-a, b), axis=-1))
def to_s3(x):        # inverse stereographic: the law's chart -> unit 3-sphere in R^4 (round 18)
    r2 = (x * x).sum(-1, keepdims=True); return np.concatenate([(1 - r2), 2 * x], -1) / (1 + r2)
def qmul(a, b):
    w1, x1, y1, z1 = a; w2, x2, y2, z2 = b
    return np.array([w1*w2 - x1*x2 - y1*y2 - z1*z2, w1*x2 + x1*w2 + y1*z2 - z1*y2,
                     w1*y2 - x1*z2 + y1*w2 + z1*x2, w1*z2 + x1*y2 - y1*x2 + z1*w2])
def gyr(a, b, c):
    return oplus(-oplus(a, b), oplus(a, oplus(b, c)))

# ================================================================== K1. The law's transport is the slice's bending (the meron)
say("== K1. The law's own transport, and the three ways to carry a direction on the 3-sphere")
rg = np.random.default_rng(25)
w_lc = w_left = w_right = 0.0
for _ in range(300):
    x = rg.normal(size=3); x *= rg.uniform(0.05, 3.0) / np.linalg.norm(x)
    v = rg.normal(size=3); e = 1e-6
    law = (to_s3(oplus(x[None], (e * v)[None]))[0] - to_s3(oplus(x[None], (-e * v)[None]))[0]) / (2 * e)
    P = to_s3(x[None])[0]; N0 = np.array([1.0, 0, 0, 0]); v4 = np.concatenate([[0], 2 * v])
    chi = math.acos(np.clip(P @ N0, -1, 1)); u = (P - math.cos(chi) * N0) / math.sin(chi)
    a_ = v4 @ u; rest = v4 - a_ * u
    lc = a_ * (math.cos(chi) * u - math.sin(chi) * N0) + rest      # Levi-Civita transport along the great circle
    s = np.abs(lc).max()
    w_lc = max(w_lc, np.abs(law - lc).max() / s)
    w_left = max(w_left, np.abs(law - qmul(P, v4)).max() / s); w_right = max(w_right, np.abs(law - qmul(v4, P)).max() / s)
w_gyr = 0.0
for _ in range(500):
    a = rg.normal(size=3); a *= rg.uniform(0.02, 0.6) / np.linalg.norm(a)
    b = rg.normal(size=3); b *= rg.uniform(0.02, 0.6) / np.linalg.norm(b)
    Gm = np.array([gyr(a[None], b[None], c[None])[0] for c in np.eye(3)]).T
    ang = math.acos(np.clip((np.trace(Gm) - 1) / 2, -1, 1))
    O = np.zeros(3); B = oplus(a[None], b[None])[0]
    s1, s2, s3 = law_dist(O[None], a[None])[0], law_dist(a[None], B[None])[0], law_dist(O[None], B[None])[0]
    sp_ = (s1 + s2 + s3) / 2
    excess = 4 * math.atan(math.sqrt(max(0.0, math.tan(sp_/2) * math.tan((sp_-s1)/2) * math.tan((sp_-s2)/2) * math.tan((sp_-s3)/2))))
    w_gyr = max(w_gyr, abs(ang - excess))
# exact: the one-parameter family nabla_{X_a} X_b = f eps_abc X_c on the 3-sphere's left-invariant frame ([X_a, X_b] = 2 eps X_c)
fs = spy.symbols('f')
def eps3(i, j, k): return spy.LeviCivita(i, j, k)
def nab(i, vec): return [fs * sum(eps3(i, j, k) * vec[j] for j in range(3)) for k in range(3)]   # nabla_{X_i} of a constant-coefficient field
E_ = [[1 if j == i else 0 for j in range(3)] for i in range(3)]
def brk(i, j): return [2 * eps3(i, j, k) for k in range(3)]
tors = [spy.simplify(nab(0, E_[1])[k] - nab(1, E_[0])[k] - brk(0, 1)[k]) for k in range(3)]
def nab_vec(vecX, vecY):   # nabla_X Y, X, Y constant-coefficient
    return [sum(vecX[i] * nab(i, vecY)[k] for i in range(3)) for k in range(3)]
Rxy_y = [spy.simplify(nab_vec(E_[0], nab_vec(E_[1], E_[1]))[k] - nab_vec(E_[1], nab_vec(E_[0], E_[1]))[k] - nab_vec(brk(0, 1), E_[1])[k]) for k in range(3)]
Ksec = spy.expand(Rxy_y[0])
tors_ok = spy.simplify(tors[2] - (2 * fs - 2)) == 0
K_ok = spy.simplify(Ksec - (2 * fs - fs**2)) == 0
check("K1 (numeric + exact; textbook geometry, our reading) the law's sum carries directions exactly as the 3-sphere's own bending (Levi-Civita) does: the change in to_s3(x (+) e v) matches transport along the great circle to 1e-9 for 300 random x up to 3 (chi up to 143 deg), while plain left or right quaternion multiplication misses by order 1. The law's twist gyr[a, b] is that bending's holonomy: its angle equals the spherical excess of the triangle 0, a, a (+) b to 1e-12 for 500 triangles. Exactly: on the 3-sphere there is a one-parameter family nabla_X Y = f [X, Y]/2 with torsion (2f - 2) and sectional curvature 2f - f^2; f = 0 and f = 2 are the two flat, twisted (torsion +-2) ways (left and right multiplication), f = 1 is torsion-free with curvature 1. The law sits at f = 1, halfway",
      w_lc < 1e-8 and w_left > 0.5 and w_right > 0.5 and w_gyr < 1e-10 and tors_ok and K_ok,
      f"law vs Levi-Civita {w_lc:.1e}; vs left {w_left:.2f}, right {w_right:.2f}; gyration angle vs excess {w_gyr:.1e}; torsion = 2f - 2: {tors_ok}; curvature = 2f - f^2: {K_ok}")
results["K1"] = dict(law_vs_levi_civita=w_lc, law_vs_left=w_left, law_vs_right=w_right, gyration_vs_excess=w_gyr,
                     torsion="2f-2", sectional_curvature=str(Ksec))

# ================================================================== K2. The twist on the law's expanding slices
say("== K2. The weak force's field locked to the law's slices, on expanding space")
t = spy.symbols('t'); fF = spy.Function('f')(t); aF = spy.Function('a')(t)
eta = [-1, 1, 1, 1]
cc = [[[0] * 4 for _ in range(4)] for _ in range(4)]          # [e_i, e_j] = cc[k][i][j] e_k ; e_0 = d/dt, e_a = X_a / a
ad = spy.diff(aF, t)
for A_ in range(1, 4):
    cc[A_][0][A_] = -ad / aF; cc[A_][A_][0] = ad / aF
    for B_ in range(1, 4):
        for C_ in range(1, 4):
            cc[C_][A_][B_] = 2 * eps3(A_ - 1, B_ - 1, C_ - 1) / aF
cl = lambda k, i, j: eta[k] * cc[k][i][j]
Gam = [[[spy.simplify((cl(k, i, j) - cl(i, j, k) + cl(j, k, i)) / 2) for j in range(4)] for i in range(4)] for k in range(4)]
Gup = lambda k, i, j: eta[k] * Gam[k][i][j]
def der(i, ex): return spy.diff(ex, t) if i == 0 else ex * 0
def ym_system(Tgen):
    """A = f sigma^a T_a with T_a matrices obeying [T_a, T_b] = eps T_c; returns F (frame), YM residual, stress tensor."""
    Z = spy.zeros(*Tgen[0].shape)
    Afr = [Z] + [fF / aF * Tgen[b] for b in range(3)]
    com = lambda P, Q: P * Q - Q * P
    F = [[spy.simplify(der(i, Afr[j]) - der(j, Afr[i]) - sum((cc[k][i][j] * Afr[k] for k in range(4)), Z) + com(Afr[i], Afr[j]))
          for j in range(4)] for i in range(4)]
    YM = []
    for j in range(4):
        s = Z
        for i in range(4):
            s = s + eta[i] * (der(i, F[i][j]) - sum((Gup(k, i, i) * F[k][j] for k in range(4)), Z)
                              - sum((Gup(k, i, j) * F[i][k] for k in range(4)), Z) + com(Afr[i], F[i][j]))
        YM.append(spy.simplify(s))
    ip = lambda P, Q: -2 * (P * Q).trace()                       # su(2) normalisation: T_a = -i sigma_a / 2 has ip = delta
    FF = sum(eta[i] * eta[j] * ip(F[i][j], F[i][j]) for i in range(4) for j in range(4))
    Tst = [[spy.simplify(sum(eta[k] * ip(F[i][k], F[j][k]) for k in range(4)) - spy.Rational(1, 4) * eta[i] * (1 if i == j else 0) * FF)
            for j in range(4)] for i in range(4)]
    return F, YM, Tst
sx = spy.Matrix([[0, 1], [1, 0]]); sy = spy.Matrix([[0, -spy.I], [spy.I, 0]]); sz = spy.Matrix([[1, 0], [0, -1]])
T2 = [-spy.I * m / 2 for m in (sx, sy, sz)]
F2, YM2, T2st = ym_system(T2)
psi = spy.Function('psi')
eom_expected = (aF**2 * spy.diff(fF, t, 2) + aF * ad * spy.diff(fF, t) - 2 * fF * (2 - fF) * (fF - 1)) / aF**3
eom_ok = all(spy.simplify(YM2[j] + eom_expected * T2[j - 1]) == spy.zeros(2, 2) for j in range(1, 4))
gauss_ok = YM2[0] == spy.zeros(2, 2)
rho = T2st[0][0]
iso_ok = all(spy.simplify(T2st[i][j] - (rho / 3 if i == j else 0)) == 0 for i in range(1, 4) for j in range(1, 4)) and \
         all(spy.simplify(T2st[0][j]) == 0 for j in range(1, 4))
rho_expected = spy.Rational(3, 2) * (aF**2 * spy.diff(fF, t)**2 + (fF**2 - 2 * fF)**2) / aF**4
rho_ok = spy.simplify(rho - rho_expected) == 0
# SU(3): the same field through its two SU(2) pieces
l = [None] * 9
l[1] = spy.Matrix([[0, 1, 0], [1, 0, 0], [0, 0, 0]]); l[2] = spy.Matrix([[0, -spy.I, 0], [spy.I, 0, 0], [0, 0, 0]])
l[3] = spy.Matrix([[1, 0, 0], [0, -1, 0], [0, 0, 0]]); l[5] = spy.Matrix([[0, 0, -spy.I], [0, 0, 0], [spy.I, 0, 0]])
l[7] = spy.Matrix([[0, 0, 0], [0, 0, -spy.I], [0, spy.I, 0]])
T3a = [-spy.I * l[k] / 2 for k in (1, 2, 3)]
T3b = [-spy.I * l[7], spy.I * l[5], -spy.I * l[2]]
idx = []
for Tg in (T3a, T3b):
    F3, YM3, T3st = ym_system(Tg)
    ok3 = all(spy.simplify(YM3[j] + eom_expected * Tg[j - 1]) == spy.zeros(3, 3) for j in range(1, 4)) and YM3[0] == spy.zeros(3, 3)
    idx.append((ok3, spy.simplify(T3st[0][0] / rho)))
# Chern-Simons number of A = lambda g^-1 dg (f = 2 lambda): tr(A dA + 2/3 A^3) per unit volume = 3f^2 - f^3 (su(2), tr), over 2 pi^2, / 8 pi^2
lam = spy.symbols('lambda')
trT = lambda P: (P).trace()
# A dA: f^2 tr(T_a T_b) sigma^a ^ d sigma^b, d sigma^b = -sum_cd eps_bcd sigma^c ^ sigma^d ; A^3: f^3 eps_abc tr(T_a T_b T_c)
AdA = sum(-trT(T2[a1] * T2[b1]) * eps3(b1, c1, d1) * eps3(a1, c1, d1) for a1 in range(3) for b1 in range(3) for c1 in range(3) for d1 in range(3))
A3 = sum(eps3(a1, b1, c1) * trT(T2[a1] * T2[b1] * T2[c1]) for a1 in range(3) for b1 in range(3) for c1 in range(3))
cs = spy.simplify((fs**2 * AdA + spy.Rational(2, 3) * fs**3 * A3) * 2 * spy.pi**2 / (8 * spy.pi**2))
cs_lam = spy.factor(cs.subs(fs, 2 * lam) / cs.subs(fs, 2))   # normalised so the flat twist lambda = 1 counts 1
cs_ok = spy.simplify(cs_lam - lam**2 * (3 - 2 * lam)) == 0
# numeric: the field on the law's own de Sitter (a = cosh(sqrt q t)/sqrt q, no back-reaction): energy a^4 rho kept, finite conformal time
qn = 1.0
def rhs(tt, y):
    f_, g_ = y; a_ = math.cosh(math.sqrt(qn) * tt) / math.sqrt(qn); a_d = math.sinh(math.sqrt(qn) * tt)
    return [g_, (2 * f_ * (2 - f_) * (f_ - 1) - a_ * a_d * g_) / a_**2]
sol = solve_ivp(rhs, (-12, 12), [1.4, 0.0], rtol=1e-11, atol=1e-13, dense_output=True)
tt = np.linspace(-12, 12, 4001); fv, gv = sol.sol(tt)
av = np.cosh(math.sqrt(qn) * tt) / math.sqrt(qn)
E4 = 1.5 * (av**2 * gv**2 + (fv**2 - 2 * fv)**2)
drift = np.abs(E4 - E4[2000]).max() / E4[2000]
etav = 2 * np.arctan(np.tanh(math.sqrt(qn) * tt / 2))
psiv = fv - 1
crossings = int(np.sum(np.diff(np.sign(psiv)) != 0))
# expected half-periods in conformal time pi for the anharmonic oscillator psi'' = 2 psi (1 - psi^2), amplitude 0.4 (inside the well)
from scipy.integrate import quad
# conformal invariance: the same field integrated directly in conformal time eta, psi'' = 2 psi (1 - psi^2)
sol_eta = solve_ivp(lambda e_, y: [y[1], 2 * y[0] * (1 - y[0]**2)], (etav[0], etav[-1]), [fv[0] - 1, av[0] * gv[0]],
                    rtol=1e-12, atol=1e-14, dense_output=True)
conf_gap = np.abs(sol_eta.sol(etav)[0] - (fv - 1)).max()
Ew = (1 - (fv[0] - 1)**2)**2
pl = math.sqrt(max(0, 1 - math.sqrt(Ew))); ph_ = math.sqrt(1 + math.sqrt(Ew))
Tper = 2 * quad(lambda p: 1 / math.sqrt(max(1e-300, Ew - (1 - p**2)**2)), pl, ph_, limit=400)[0]
eta_span = etav[-1] - etav[0]
exp_turns = eta_span / Tper
ok_ds = drift < 1e-7 and conf_gap < 1e-6 and abs(eta_span - math.pi) < 1e-4
check("K2 (exact + numeric; textbook Yang-Mills on the 3-sphere, Luscher 1977, Cervero-Jacobs 1978) put the weak force's SU(2) field on the law's expanding slices locked to them, A = f(t) sigma^a T_a (colour a = direction a, possible because the slice is SU(2) itself, round 18). The full Yang-Mills equations reduce exactly to one equation, which in conformal time d eta = dt/a is psi'' = 2 psi (1 - psi^2) with f = 1 + psi: a ball rolling in a double well whose two bottoms psi = +-1 (f = 0, 2) are the two flat twists and whose hilltop psi = 0 is the law's f = 1. Gauss's law holds identically. Its stress is exactly isotropic and traceless (pressure = rho/3, like light), rho = (3/2)(a^2 f'^2 + (f^2 - 2f)^2)/a^4. The strong force's SU(3) carries the same field through its SU(2) pieces, with energy 1 or 4 times as much (lambda_1,2,3 or lambda_7,5,2). Chern-Simons number of A = lambda g^-1 dg is lambda^2 (3 - 2 lambda): the law's f = 1 has 1/2, the top of the barrier between vacua 0 and 1. On the law's own de Sitter, a^4 rho is kept to 1e-7 and the field gets only pi of conformal time, so it makes under one full swing (period >= pi near the well bottom) and then freezes",
      eom_ok and gauss_ok and iso_ok and rho_ok and all(o for o, _ in idx) and [int(v) for _, v in idx] == [1, 4] and cs_ok and ok_ds,
      f"YM -> psi'' = 2psi(1-psi^2): {eom_ok}; Gauss: {gauss_ok}; isotropic traceless: {iso_ok}; rho formula: {rho_ok}; SU(3) embeddings ok/energy ratio {[(bool(o), str(v)) for o, v in idx]}; CS = {cs_lam}; de Sitter: a^4 rho drift {drift:.1e}, cosmic-time field vs conformal-time field {conf_gap:.1e}, conformal time available {eta_span:.5f} (pi), swings {exp_turns:.2f} periods")
results["K2"] = dict(eom="psi''=2psi(1-psi^2) in conformal time", isotropic=bool(iso_ok), su3_energy_ratio=[str(v) for _, v in idx],
                     chern_simons=str(cs_lam), de_sitter_energy_drift=float(drift), conformal_gap=float(conf_gap), periods_in_all_time=exp_turns,
                     eta_span=float(eta_span))
SAVE_DS = dict(t=tt[::8].tolist(), eta=etav[::8].tolist(), f=fv[::8].tolist(), E4=E4[::8].tolist())

# ------------------------------------------------------------------ curvature of a diagonal metric (symbolic; round 23's routine, unchanged)
def curvature(gdiag, X):
    n = len(X); g = spy.diag(*gdiag); gi = spy.diag(*[1 / v for v in gdiag])
    Gm = [[[spy.simplify(sum(gi[a, d] * (spy.diff(g[d, b], X[c]) + spy.diff(g[d, c], X[b]) - spy.diff(g[b, c], X[d])) for d in range(n)) / 2)
            for c in range(n)] for b in range(n)] for a in range(n)]
    R = [[[[spy.simplify(spy.diff(Gm[a][b][d], X[c]) - spy.diff(Gm[a][b][c], X[d])
                         + sum(Gm[a][c][e] * Gm[e][b][d] - Gm[a][d][e] * Gm[e][b][c] for e in range(n)))
            for d in range(n)] for c in range(n)] for b in range(n)] for a in range(n)]
    Ric = spy.Matrix(n, n, lambda b, d: spy.simplify(sum(R[a][b][a][d] for a in range(n))))
    Rs = spy.simplify(sum(gi[b, b] * Ric[b, b] for b in range(n)))
    Rl = [[[[g[a, a] * R[a][b][c][d] for d in range(n)] for c in range(n)] for b in range(n)] for a in range(n)]
    return g, gi, Ric, Rs, Rl
def weyl_sq(g, gi, Ric, Rs, Rl, n=4):
    def C(a, b, c, d):
        return (Rl[a][b][c][d] - (g[a, c] * Ric[b, d] - g[a, d] * Ric[b, c] - g[b, c] * Ric[a, d] + g[b, d] * Ric[a, c]) / (n - 2)
                + Rs * (g[a, c] * g[b, d] - g[a, d] * g[b, c]) / ((n - 1) * (n - 2)))
    W2 = 0
    for a in range(n):
        for b in range(n):
            for c in range(n):
                for d in range(n):
                    W2 += gi[a, a] * gi[b, b] * gi[c, c] * gi[d, d] * C(a, b, c, d)**2
    return spy.simplify(W2)

# ================================================================== K3. The twist feeds the dial, evenly: no tide
say("== K3. The twist's energy pulls the law's slices back, evenly")
tS, chS, thS, phS = spy.symbols('t chi theta phi', positive=True)
aS = spy.Function('a')(tS)
gF, giF, RicF, RsF, RlF = curvature([-1, aS**2, aS**2 * spy.sin(chS)**2, aS**2 * spy.sin(chS)**2 * spy.sin(thS)**2], [tS, chS, thS, phS])
W2F = weyl_sq(gF, giF, RicF, RsF, RlF)
# closed universe with the law's Lambda = 3q and the locked twist (kappa = 4 pi G / g^2): H^2 = kappa C / a^4 + q - 1/a^2,
# a''/a = q - kappa C / a^4, f'' from K2, C = a^2 f'^2 + (f^2 - 2f)^2 kept
rows3 = []
for (qq, kap, f0) in ((0.30, 0.05, 1.5), (0.10, 0.20, 1.2), (0.50, 0.02, 1.8)):
    C0 = (f0**2 - 2 * f0)**2                                 # start at rest: f' = 0
    # time symmetry a' = 0: 1/a^2 = q + kappa C / a^4  ->  q a^4 - a^2 + kappa C = 0  (bounce: smaller root of a^2 ... take larger a^2 root of the turning point)
    disc = 1 - 4 * qq * kap * C0
    a0 = math.sqrt((1 + math.sqrt(disc)) / (2 * qq))        # the bounce of the de Sitter-like branch (minimum of a)
    def rhs3(tt_, y):
        a_, ap, f_, fp = y
        return [ap, a_ * (qq - kap * (a_**2 * fp**2 + (f_**2 - 2 * f_)**2) / a_**4),
                fp, (2 * f_ * (2 - f_) * (f_ - 1) - a_ * ap * fp) / a_**2]
    s3 = solve_ivp(rhs3, (0, 6), [a0, 0.0, f0, 0.0], rtol=1e-11, atol=1e-13, dense_output=True)
    tg = np.linspace(0, 6, 1201); a_, ap, f_, fp = s3.sol(tg)
    Cv = a_**2 * fp**2 + (f_**2 - 2 * f_)**2
    ham = ap**2 / a_**2 - (kap * Cv / a_**4 + qq - 1 / a_**2)
    qin = 1 / a0**2; rho0 = 3 * C0 / (2 * a0**4)              # rho in units of 1/g^2; 8 pi G rho / 3 = kappa C / a^4
    rows3.append(dict(q=qq, kappa=kap, f0=f0, a_bounce=a0, C_drift=float(np.abs(Cv - C0).max() / C0),
                      friedmann_residual=float(np.abs(ham).max()), q_in=qin, q_plus_kappaC_over_a4=qq + kap * C0 / a0**4))
ok3 = spy.simplify(W2F) == 0 and all(r_["C_drift"] < 1e-7 and r_["friedmann_residual"] < 1e-8 and abs(r_["q_in"] - r_["q_plus_kappaC_over_a4"]) < 1e-12 for r_ in rows3)
check("K3 (exact + numeric; textbook Friedmann with radiation, our reading) let the locked twist pull back on the law's spacetime: a closed universe with the law's Lambda = 3q and the twist's energy. The slices stay round 3-spheres (Weyl^2 of every such universe is exactly 0, whatever a(t)), so the twist adds no tidal pull. It only changes the dial: at the moment the universe stops shrinking (a' = 0) the slices are the law's 3-sphere at q_in = 1/a^2 = q + 8 pi G rho/3, the same rule round 23 found for dust. Integrated together (Einstein's second equation driving a, Yang-Mills driving f), Einstein's first equation holds to 1e-8 and the twist's conformal energy C is kept to 1e-7, for three settings",
      ok3, f"Weyl^2 of the closed universe: {spy.simplify(W2F)}; " + "; ".join(f"q {r_['q']}, kappa {r_['kappa']}: Friedmann {r_['friedmann_residual']:.1e}, C drift {r_['C_drift']:.1e}, q_in {r_['q_in']:.6f} = q + kappa C/a^4 {r_['q_plus_kappaC_over_a4']:.6f}" for r_ in rows3))
results["K3"] = dict(weyl_sq_closed_universe=str(spy.simplify(W2F)), runs=rows3)

# ================================================================== K4. A dial that varies must make a tide (Schur), and how much
say("== K4. What a varying dial does to the slice")
r_, qs = spy.symbols('r q', positive=True); m_ = spy.Function('m')(r_)
N_ = 1 - 2 * m_ / r_ - qs * r_**2
g3, gi3, Ric3, Rs3, Rl3 = curvature([1 / N_, r_**2, r_**2 * spy.sin(thS)**2], [r_, thS, phS])
K_perp = spy.simplify(Rl3[1][2][1][2] / (g3[1, 1] * g3[2, 2]))       # plane tangent to the sphere of radius r
K_par = spy.simplify(Rl3[0][1][0][1] / (g3[0, 0] * g3[1, 1]))        # plane containing the radial direction
mp = spy.diff(m_, r_)
q_bar = qs + 2 * m_ / r_**3                 # the dial averaged over the ball inside r (2m/r^3 = (8 pi/3) x mean density)
q_loc = qs + spy.Rational(2, 3) * mp / r_**2   # the dial the law would put at r from the local density (8 pi rho/3, rho = m'/(4 pi r^2))
k4a = spy.simplify(K_perp - q_bar) == 0
k4b = spy.simplify(K_par - (qs + mp / r_**2 - m_ / r_**3)) == 0
k4c = spy.simplify(Rs3 - (6 * qs + 4 * mp / r_**2)) == 0                     # = 6q + 16 pi rho: the time-symmetric constraint
tide = spy.simplify(K_perp - K_par)
k4d = spy.simplify(tide - spy.Rational(3, 2) * (q_bar - q_loc)) == 0 and spy.simplify(tide + r_ / 2 * spy.diff(q_bar, r_)) == 0
cS = spy.symbols('c'); k4e = spy.simplify(tide.subs(m_, cS * r_**3).doit()) == 0    # Schur: isotropic everywhere iff m = c r^3, a single dial
Ms = spy.symbols('M', positive=True)
Ric3M = Ric3.subs(m_, Ms).doit(); Rs3M = spy.simplify(Rs3.subs(m_, Ms).doit())
TF = spy.Matrix(3, 3, lambda i, j: spy.simplify(Ric3M[i, j] - Rs3M * g3.subs(m_, Ms)[i, j] / 3))
gi3M = gi3.subs(m_, Ms)
E2 = spy.simplify(sum(gi3M[i, i] * gi3M[j, j] * TF[i, j]**2 for i in range(3) for j in range(3)))
k4f = spy.simplify(8 * E2 - 48 * Ms**2 / r_**6) == 0
check("K4 (exact, symbolic; textbook: Schur 1886, the time-symmetric constraint, and the Weyl tensor's electric part) a moment-of-rest slice of any round lump, dr^2/N + r^2 dOmega^2 with N = 1 - 2m(r)/r - q r^2. Its curvature across the sphere through r is q_bar = q + 2m/r^3, the dial averaged over everything inside r; along the radius it is q + m'/r^2 - m/r^3; their total is 6q + 16 pi rho, Einstein's constraint. The tide (the difference) is exactly 3/2 (q_bar - q_loc), where q_loc = q + 8 pi rho/3 is the dial the law would set from the matter at r, and equally -(r/2) dq_bar/dr: tidal pull is the slope of the averaged dial. It vanishes everywhere only when m = c r^3, one dial for the whole slice (Schur's lemma: a space that bends the same in every direction at each point bends the same everywhere). So a dial that varies cannot stay law-like at any point it varies through. Outside the lump the tide is 3M/r^3 and gives Weyl^2 = 8 E.E = 48 M^2/r^6, round 23's value",
      k4a and k4b and k4c and k4d and k4e and k4f,
      f"K_perp = q_bar: {k4a}; K_par formula: {k4b}; total = 6q + 4m'/r^2: {k4c}; tide = 3/2(q_bar - q_loc) = -(r/2) q_bar': {k4d}; Schur: {k4e}; outside 8 E.E = 48M^2/r^6: {k4f}")
results["K4"] = dict(K_perp=str(K_perp), K_par=str(K_par), scalar=str(Rs3), tide=str(tide), weyl_outside="48*M**2/r**6")

# ================================================================== K5. Where they meet: a bubble of twist
say("== K5. A bubble of the weak force's twist: twist -> dial -> tide")
def ss(x): return np.where(x <= 0, 0.0, np.where(x >= 1, 1.0, x * x * (3 - 2 * x)))
def dss(x): return np.where((x <= 0) | (x >= 1), 0.0, 6 * x * (1 - x))
def make_profile(r1, R):
    def w(r): return np.cos(2 * math.pi * ss((r - r1) / (R - r1)))
    def wp(r): return -np.sin(2 * math.pi * ss((r - r1) / (R - r1))) * 2 * math.pi * dss((r - r1) / (R - r1)) / (R - r1)
    return w, wp
def bubble(w, wp, kap, qq, rmax, n=40001):
    rr = np.linspace(1e-6, rmax, n)
    def mrhs(r, y):
        N = 1 - 2 * y[0] / r - qq * r * r
        return [kap * (N * wp(r)**2 + (1 - w(r)**2)**2 / (2 * r * r))]
    sol = solve_ivp(mrhs, (rr[0], rr[-1]), [0.0], t_eval=rr, rtol=1e-11, atol=1e-14, max_step=rmax / 4000)
    m = sol.y[0]; N = 1 - 2 * m / rr - qq * rr**2
    rho = kap * (N * wp(rr)**2 + (1 - w(rr)**2)**2 / (2 * rr**2)) / (4 * math.pi * rr**2)
    return rr, m, N, rho
qb = 0.04
prof = dict(ball=make_profile(0.0, 1.2), shell=make_profile(0.7, 1.4))
from scipy.optimize import brentq as _brentq
target_M = 0.18
rows5 = {}
for name, (w, wp) in prof.items():
    kap = _brentq(lambda k: bubble(w, wp, k, qb, 1.6, 8001)[1][-1] - target_M, 1e-4, 1.0, xtol=1e-12)
    rr, m, N, rho = bubble(w, wp, kap, qb, 4.0)
    Kperp_num = (1 - N) / rr**2                                    # read from the metric
    Kpar_num = -np.gradient(N, rr) / (2 * rr)
    R3 = 2 * (Kperp_num + 2 * Kpar_num)
    qbar = Kperp_num; qloc = qb + 8 * math.pi * rho / 3
    sel = (rr > 0.05) & (rr < 3.9)
    ham = np.abs(R3 - (6 * qb + 16 * math.pi * rho))[sel].max()
    tide_metric = (Kperp_num - Kpar_num)[sel]
    tide_energy = (1.5 * (qbar - qloc))[sel]
    out = rr > 1.6
    M = m[-1]
    rows5[name] = dict(kappa=kap, M=float(M), ham=float(ham), tide_gap=float(np.abs(tide_metric - tide_energy).max()),
                       outside_vs_3M_r3=float(np.abs((Kperp_num - Kpar_num)[out & (rr < 3.9)] - 3 * M / rr[out & (rr < 3.9)]**3).max()),
                       min_N=float(N[rr < 3.9].min()),
                       inner_tide=float(np.abs((Kperp_num - Kpar_num)[(rr > 0.05) & (rr < 0.6)]).max()) if name == "shell" else None,
                       r=rr[::40].tolist(), qbar=qbar[::40].tolist(), qloc=qloc[::40].tolist(), tide=(Kperp_num - Kpar_num)[::40].tolist(),
                       Kpar=Kpar_num[::40].tolist(), w=w(rr)[::40].tolist())
# the uniform twist (the law's own f = 1 field filling a whole 3-sphere of radius a): w = cos chi, r = a sin chi
kap_u, q_u = 0.05, 0.3
a_u = math.sqrt((1 + math.sqrt(1 - 4 * q_u * kap_u)) / (2 * q_u))       # 1/a^2 = q + kappa/a^4 (rho of the meron: 8 pi rho/3 = kappa/a^4)
wu = lambda r: np.cos(np.arcsin(np.clip(r / a_u, 0, 1)))
wup = lambda r: -(r / a_u) / np.sqrt(np.clip(1 - (r / a_u)**2, 1e-300, None)) / a_u
rr_u, m_u, N_u, rho_u = bubble(wu, wup, kap_u, q_u, 0.95 * a_u)
sel_u = rr_u > 0.02 * a_u
qbar_u = (1 - N_u) / rr_u**2
uni_flat = np.abs(qbar_u[sel_u] - 1 / a_u**2).max(); uni_rho = np.abs(rho_u[sel_u] * 8 * math.pi / 3 - kap_u / a_u**4).max()
okB, okS = rows5["ball"], rows5["shell"]
ok5 = (okB["ham"] < 2e-4 and okS["ham"] < 2e-4 and okB["tide_gap"] < 2e-4 and okS["tide_gap"] < 2e-4 and okB["outside_vs_3M_r3"] < 1e-6
       and okS["outside_vs_3M_r3"] < 1e-6 and abs(okB["M"] - okS["M"]) < 1e-8 and okS["inner_tide"] < 1e-6 and min(okB["min_N"], okS["min_N"]) > 0
       and uni_flat < 1e-7 and uni_rho < 1e-7)
check("K5 (numeric; textbook Einstein-Yang-Mills at a moment of rest, Bartnik-McKinnon 1988 field shape; the chain is our reading) where the two threads meet. The weak force's field, in the standard round shape A ~ (1 - w(r)), costs energy kappa [N w'^2 + (1 - w^2)^2/(2 r^2)] per unit radius (kappa = 4 pi G/g^2). (a) On a whole 3-sphere of radius a, the shape w = cos chi has exactly even energy, 8 pi rho/3 = kappa/a^4, the same as the law's own f = 1 field in K2, the slice is the law's 3-sphere at q_in = q + kappa/a^4 to 1e-7, and there is no tide anywhere. (b) A bubble of the same field, ball-shaped or a hollow shell, in the law's de Sitter (q = 0.04): solving Einstein's constraint, the tide read from the metric equals 3/2 (q_bar - q_loc) from the field's energy to 2e-4 at every radius; it is zero inside the hollow shell (Newton's shell theorem), nonzero in the field, and outside it is 3M/r^3 with M the field's total energy, the same for both shapes when M matches. So the chain runs one way: the twist sets the dial, a varying dial sets the tide, and the tide outside remembers only the total",
      ok5, f"uniform: q_bar - 1/a^2 {uni_flat:.1e}, 8 pi rho/3 - kappa/a^4 {uni_rho:.1e}; " + "; ".join(
          f"{nm}: kappa {v['kappa']:.4f}, M {v['M']:.6f}, constraint {v['ham']:.1e}, tide metric vs energy {v['tide_gap']:.1e}, outside vs 3M/r^3 {v['outside_vs_3M_r3']:.1e}, min N {v['min_N']:.3f}" + (f", tide inside hollow {v['inner_tide']:.1e}" if v['inner_tide'] is not None else "")
          for nm, v in rows5.items()))
results["K5"] = dict(uniform=dict(a=a_u, q=q_u, kappa=kap_u, qbar_gap=float(uni_flat), rho_gap=float(uni_rho)),
                     bubbles={k: {kk: vv for kk, vv in v.items() if kk not in ("r", "qbar", "qloc", "tide", "Kpar", "w")} for k, v in rows5.items()})

# ------------------------------------------------------------------ write out
say(f"\n{passed} passed, {failed} failed  ({time.time() - T_START:.0f} s)")
results["summary"] = dict(passed=passed, failed=failed, seconds=round(time.time() - T_START))
results["figure"] = dict(ds=SAVE_DS, bubbles={k: {kk: v[kk] for kk in ("r", "qbar", "qloc", "tide", "Kpar", "w", "M")} for k, v in rows5.items()}, q=qb)
json.dump(results, open("results.json", "w"), indent=1, default=str)
open("checks_output.txt", "w").write(_log.getvalue())
