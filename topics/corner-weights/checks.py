"""Corner weights: one matrix for size, pair and Lucas weights (Nov 25, 2026).

Corner sum  Sigma_S w(S) prod_{i in S} x_i  over the n-cube.  Claim tested here:
every weighting we have met is  det(I + t M)  read in two independent choices:
  * which matrix M  (its principal minors are the corners): diagonal = no pairs (ours),
    Cauchy = pair weights (Hirota/KdV/KP, Belinski-Zakharov), generic = 3-body terms;
  * which dial t you evaluate at: a number k (our size weight k^(m-1)) or a 2x2 matrix N with
    trace k, det d (the Lucas weights U_m(k,d)).
Run: python3 checks.py > checks_output.txt   (numpy, sympy)
"""
import itertools, math, json
import numpy as np
import sympy as sp

rng = np.random.default_rng(7)
results, store = [], {}
def check(name, ok, info=""):
    results.append((name, bool(ok)))
    print(("PASS " if ok else "FAIL ") + name + (("   " + info) if info else ""))
def info(s): print("   " + s)
det = np.linalg.det
def subsets(n):
    for r in range(n + 1):
        for S in itertools.combinations(range(n), r):
            yield S
def minor(M, S):
    return 1.0 if len(S) == 0 else det(M[np.ix_(S, S)])
def circ(a, b, k): return a + b + k * a * b
def fold(vals, k):
    out = 0
    for v in vals: out = circ(out, v, k)
    return out

print("== A. Corners are principal minors of one matrix ==")
n = 5
M = rng.normal(size=(n, n)) * 0.3
for t in (1.0, 0.7, -1.3):
    s = sum(t ** len(S) * minor(M, S) for S in subsets(n))
    check(f"A1 sum over corners of t^|S| det M_S = det(I + tM)  (t={t})", np.isclose(s, det(np.eye(n) + t * M)))
k = 0.6
lam = np.linalg.eigvals(M)
corner = sum(k ** (len(S) - 1) * minor(M, S) for S in subsets(n) if S)
check("A1b size-weighted corner sum (k^(m-1) on rank-m minors) = our operation folded over the EIGENVALUES",
      np.isclose(corner, fold(lam, k)), f"{corner:.6f} vs {complex(fold(lam, k)).real:.6f}  (known: linear-algebra page)")
x = rng.uniform(0.2, 1.5, size=4)
D = np.diag(x)
check("A2 diagonal M: corners are prod x_i, every pair factor = 1 (our operation, no interaction)",
      all(np.isclose(minor(D, S), np.prod(x[list(S)])) for S in subsets(4)))
# KdV Cauchy matrix
p = np.array([0.6, 1.1, 1.9, 2.8])
Mk = np.array([[2 * p[i] * x[i] / (p[i] + p[j]) for j in range(4)] for i in range(4)])
A = lambda i, j, p=p: ((p[i] - p[j]) / (p[i] + p[j])) ** 2
ok = all(np.isclose(minor(Mk, S), np.prod(x[list(S)]) * math.prod(A(i, j) for i, j in itertools.combinations(S, 2)))
         for S in subsets(4))
check("A3 Cauchy M_ij = 2 p_i x_i/(p_i+p_j): every corner = prod x_i * prod over pairs A_ij (Hirota's KdV weights fall out)", ok)
comb = sum(k ** (len(S) - 1) * np.prod(x[list(S)]) * math.prod(A(i, j) for i, j in itertools.combinations(S, 2))
           for S in subsets(4) if S)
check("A3b the size x pair candidate w(S)=k^(|S|-1) prod A_ij is exactly (det(I+kM)-1)/k = our op on eigenvalues of the Cauchy matrix",
      np.isclose(comb, (det(np.eye(4) + k * Mk) - 1) / k) and np.isclose(comb, fold(np.linalg.eigvals(Mk), k).real))
lamK = np.sort(np.linalg.eigvals(Mk).real)
check("A3c the Cauchy matrix's eigenvalues are real and positive ('dressed waves')", np.all(lamK > 0) and np.allclose(np.linalg.eigvals(Mk).imag, 0),
      "bare " + str(np.round(np.sort(x), 3)) + " dressed " + str(np.round(lamK, 3)))
# A4 pair factor and three-body term for a general matrix
def pairf(M, i, j): return minor(M, (i, j)) / (M[i, i] * M[j, j])
def threebody(M, S=(0, 1, 2)):
    i, j, l = S
    return minor(M, S) / (M[i, i] * M[j, j] * M[l, l]) / (pairf(M, i, j) * pairf(M, i, l) * pairf(M, j, l)) - 1
G = rng.normal(size=(3, 3)); G = G @ G.T + 2 * np.eye(3)
check("A4 any matrix: pair factor A_ij = 1 - m_ij m_ji/(m_ii m_jj); generic matrix has a 3-body term, Cauchy has none",
      abs(threebody(G)) > 1e-3 and abs(threebody(Mk)) < 1e-10, f"generic 3-body {threebody(G):+.4f}, Cauchy {threebody(Mk):+.1e}")
# A5 Sylvester / source rank
pp = np.array([0.5, 1.2, 2.0]); qq = np.array([0.9, 0.4, 1.6])
P, Q = np.diag(pp), np.diag(qq)
r1, s1 = rng.uniform(0.5, 1.5, 3), rng.uniform(0.5, 1.5, 3)
r2, s2 = rng.uniform(0.5, 1.5, 3), rng.uniform(0.5, 1.5, 3)
C = 1 / (pp[:, None] + qq[None, :])
M1 = np.outer(r1, s1) * C
check("A5a Cauchy-like M = r s^T/(p_i+q_j) solves Sylvester P M + M Q = r s^T (one shared channel, rank 1)",
      np.allclose(P @ M1 + M1 @ Q, np.outer(r1, s1)))
check("A5b rank-1 source: corners are pair products for ANY r, s (no 3-body term)", abs(threebody(M1)) < 1e-10, f"{threebody(M1):+.1e}")
tb = []
for eps in (0.0, 0.01, 0.1, 0.3, 1.0):
    M2 = (np.outer(r1, s1) + eps * np.outer(r2, s2)) * C
    tb.append((eps, threebody(M2)))
info("source r s^T + eps r' s'^T (rank 2):  " + ",  ".join(f"eps={e}: 3-body {v:+.2e}" for e, v in tb))
check("A5c a second channel (rank-2 source) creates a 3-body term that grows with its strength", abs(tb[0][1]) < 1e-10 and abs(tb[-1][1]) > 1e-3 and abs(tb[2][1]) > abs(tb[1][1]))
store["source_rank"] = tb
Rd = np.diag(rng.uniform(0.5, 1.5, 3))
Md = Rd * C  # diagonal source -> diagonal M
check("A5d diagonal source (each wave its own channel) gives a diagonal M: our size-only weights",
      np.allclose(Md, np.diag(np.diag(Md))))
# A6 exclusion
pe = np.array([1.0, 1.0 + 1e-7, 2.0]); xe = np.array([0.7, 0.9, 1.1])
Me = np.array([[2 * pe[i] * xe[i] / (pe[i] + pe[j]) for j in range(3)] for i in range(3)])
ev = np.sort(np.abs(np.linalg.eigvals(Me)))
check("A6 two waves at the same speed: pair factor 0 and the matrix loses rank (one dressed eigenvalue -> 0): exclusion",
      A(0, 1, pe) < 1e-13 and ev[0] < 1e-6, f"smallest eigenvalue {ev[0]:.1e}")
# A7 pair factor = velocity subtraction / cross-ratio
u = (p - 1) / (p + 1)
vs = lambda a, b: (a - b) / (1 - a * b)
check("A7a KdV: sqrt(A_ij) = (p_i-p_j)/(p_i+p_j) = relativistic velocity subtraction of u=(p-1)/(p+1) = tanh(theta/2)",
      all(np.isclose((p[i] - p[j]) / (p[i] + p[j]), vs(u[i], u[j])) for i in range(4) for j in range(4)))
cr = lambda i, j: (pp[i] + qq[i]) * (pp[j] + qq[j]) / ((pp[i] + qq[j]) * (pp[j] + qq[i]))
Akp = lambda i, j: (pp[i] - pp[j]) * (qq[i] - qq[j]) / ((pp[i] + qq[j]) * (pp[j] + qq[i]))
check("A7b KP: pair factor = 1 - cross-ratio of the four points (p_i, p_j, -q_i, -q_j); equals the Cauchy minor",
      all(np.isclose(Akp(i, j), 1 - cr(i, j)) and np.isclose(pairf(M1, i, j), Akp(i, j)) for i, j in itertools.combinations(range(3), 2)))

print("\n== B. The size dial k in the soliton world ==")
X_, T_, Y_ = sp.symbols('x t y', real=True)
def tau_from_weights(etas, wfun, N):
    return sum(wfun(S) * sp.exp(sum(etas[i] for i in S)) for S in itertools.combinations(range(N), r) for r in [None]) if False else \
        sum(wfun(S) * sp.exp(sum((etas[i] for i in S), sp.Integer(0))) for r in range(N + 1) for S in itertools.combinations(range(N), r))
def kdv_bilinear(tau):
    # D_x(D_t + D_x^3) tau.tau = 2(tau tau_xt - tau_x tau_t) + 2(tau tau_xxxx - 4 tau_x tau_xxx + 3 tau_xx^2)
    d = lambda *a: sp.diff(tau, *a)
    return 2 * (tau * d(X_, T_) - d(X_) * d(T_)) + 2 * (tau * d(X_, 4) - 4 * d(X_) * d(X_, 3) + 3 * d(X_, 2) ** 2)
pts = [(sp.Rational(3, 10), sp.Rational(1, 10)), (sp.Rational(-6, 5), sp.Rational(2, 5)), (sp.Rational(2, 1), sp.Rational(-3, 10))]
def maxres(expr, tau):
    out = 0
    for xv, tv in pts:
        num = abs(sp.N(expr.subs({X_: xv, T_: tv}), 40)); den = abs(sp.N((tau ** 2).subs({X_: xv, T_: tv}), 40))
        out = max(out, float(num / den))
    return out
ks = [sp.Rational(7, 10), sp.Rational(13, 10), sp.Rational(21, 10)]
eta = [kk * X_ - kk ** 3 * T_ for kk in ks]
Aexact = lambda i, j: ((ks[i] - ks[j]) / (ks[i] + ks[j])) ** 2
pairw = lambda S: sp.prod([Aexact(i, j) for i, j in itertools.combinations(S, 2)]) if len(S) > 1 else 1
kdvres = {}
for kv in (1, sp.Rational(3, 10), 3, sp.Rational(-1, 2)):
    tau = tau_from_weights(eta, lambda S, kv=kv: kv ** len(S) * pairw(S), 3)
    kdvres[str(kv)] = maxres(kdv_bilinear(tau), tau)
info("KdV bilinear residual with corner weights k^|S| * prod A:  " + ", ".join(f"k={a}: {b:.1e}" for a, b in kdvres.items()))
check("B1 the size dial is a gauge for KdV: k^|S| prod A solves KdV for every k (k>0 slides each wave by -log(k)/p_i; k<0 gives singular solitons)",
      max(kdvres.values()) < 1e-25)
xs = rng.uniform(0.3, 1.2, 3); kk_ = 2.3
same = all(np.isclose(kk_ ** len(S) * np.prod(xs[list(S)]), np.prod((kk_ * xs)[list(S)])) for S in subsets(3))
check("B1b corner by corner: k^|S| prod x_i = prod (k x_i), i.e. multiply each wave by k (shift its phase by log k)", same)
tau_g = tau_from_weights(eta, lambda S: {0: 1, 1: 1, 2: 1, 3: 2}[len(S)] * pairw(S), 3)
tau_d = tau_from_weights(eta, lambda S: 1, 3)
rg, rd = maxres(kdv_bilinear(tau_g), tau_g), maxres(kdv_bilinear(tau_d), tau_d)
check("B2 a size weight that is not geometric (g(3) != g(1)g(2)^2/g(1)... here 1,1,2) breaks KdV", rg > 1e-3, f"residual {rg:.2e}")
check("B3 the diagonal matrix (our op, A=1) is not a KdV solution: separate solitons that do not collide", rd > 1e-3, f"residual {rd:.2e}")
store["kdv_gauge"] = kdvres
# B4 interpolation
def tau_M(Msym):
    return (sp.eye(3) + Msym).det()
interp = []
for s in (0, sp.Rational(1, 4), sp.Rational(1, 2), sp.Rational(3, 4), sp.Rational(9, 10), 1):
    Ms = sp.Matrix(3, 3, lambda i, j: (1 if i == j else s) * 2 * ks[i] / (ks[i] + ks[j]) * sp.exp(eta[i]))
    tau = tau_M(Ms)
    interp.append((float(s), maxres(kdv_bilinear(tau), tau)))
info("M(s) = diagonal + s*(Cauchy off-diagonal):  " + ", ".join(f"s={a}: {b:.1e}" for a, b in interp))
check("B4 straight-line blend from our diagonal matrix to the Cauchy matrix solves KdV only at s = 1", all(r > 1e-4 for s_, r in interp[:-1]) and interp[-1][1] < 1e-25)
store["interp_s"] = interp
spread = []
for rr in (1.5, 3, 10, 30, 100, 1000):
    ps = np.array([1.0, rr, rr ** 2]); xx = np.array([0.8, 1.1, 0.6])
    Mc = np.array([[2 * ps[i] * xx[i] / (ps[i] + ps[j]) for j in range(3)] for i in range(3)])
    gmax = max(1 - ((ps[i] - ps[j]) / (ps[i] + ps[j])) ** 2 for i, j in itertools.combinations(range(3), 2))
    spread.append((rr, gmax, float(np.max(np.abs(np.sort(np.linalg.eigvals(Mc).real) - np.sort(xx))))))
info("speed ratio r between neighbours: " + ", ".join(f"r={a}: max coupling {b:.3f}, |dressed-bare| {c:.3f}" for a, b, c in spread))
check("B5 the integrable road to our op: as the speeds spread apart, couplings -> 0, all A -> 1, dressed eigenvalues -> bare waves",
      spread[-1][1] < 0.01 and spread[-1][2] < 0.05 and all(spread[i+1][2] < spread[i][2] for i in range(len(spread)-1)))
store["spread"] = spread

print("\n== C. KP line solitons in different directions ==")
def P_kp(a, b, c): return a ** 4 - 4 * a * c + 3 * b ** 2
def kp_bilinear(tau):
    d = lambda *a: sp.diff(tau, *a)
    xxxx = 2 * (tau * d(X_, 4) - 4 * d(X_) * d(X_, 3) + 3 * d(X_, 2) ** 2)
    xt = 2 * (tau * d(X_, T_) - d(X_) * d(T_))
    yy = 2 * (tau * d(Y_, 2) - d(Y_) ** 2)
    return xxxx - 4 * xt + 3 * yy
pk = [sp.Rational(1, 2), sp.Rational(6, 5), sp.Rational(-1, 5)]
qk = [sp.Rational(3, 10), sp.Rational(-2, 5), sp.Rational(7, 5)]
om = [((pk[i] + qk[i]), (pk[i] ** 2 - qk[i] ** 2), (pk[i] ** 3 + qk[i] ** 3)) for i in range(3)]
check("C1a dispersion: eta = (p+q)x + (p^2-q^2)y + (p^3+q^3)t is on the KP shell", all(sp.simplify(P_kp(*w)) == 0 for w in om))
Ahir = lambda i, j: -P_kp(*[om[i][a] - om[j][a] for a in range(3)]) / P_kp(*[om[i][a] + om[j][a] for a in range(3)])
Acau = lambda i, j: (pk[i] - pk[j]) * (qk[i] - qk[j]) / ((pk[i] + qk[j]) * (pk[j] + qk[i]))
check("C1b Hirota's KP pair factor -P(w_i-w_j)/P(w_i+w_j) = the Cauchy minor (p_i-p_j)(q_i-q_j)/((p_i+q_j)(p_j+q_i))",
      all(sp.simplify(Ahir(i, j) - Acau(i, j)) == 0 for i, j in itertools.combinations(range(3), 2)))
etak = [om[i][0] * X_ + om[i][1] * Y_ + om[i][2] * T_ for i in range(3)]
Mkp = sp.Matrix(3, 3, lambda i, j: (pk[i] + qk[i]) / (pk[i] + qk[j]) * sp.exp(etak[i]))
tauk = (sp.eye(3) + Mkp).det()
ptsy = [(sp.Rational(3, 10), sp.Rational(1, 5), sp.Rational(1, 10)), (sp.Rational(-1, 1), sp.Rational(-3, 5), sp.Rational(2, 5)), (sp.Rational(1, 2), sp.Rational(2, 1), sp.Rational(-1, 5))]
bk = kp_bilinear(tauk)
rk = max(float(abs(sp.N(bk.subs({X_: a, Y_: b, T_: c}), 40)) / abs(sp.N((tauk ** 2).subs({X_: a, Y_: b, T_: c}), 40))) for a, b, c in ptsy)
dirs = [float(om[i][1] / om[i][0]) for i in range(3)]
check("C1c three KP line solitons in three directions: tau = det(I + Cauchy matrix) solves KP exactly", rk < 1e-25,
      f"residual {rk:.1e}, slopes y-direction/x-direction {np.round(dirs, 2)}")
# C2 equal-A plane reachable by KP?
def Akp_num(p_, q_, i, j): return (p_[i] - p_[j]) * (q_[i] - q_[j]) / ((p_[i] + q_[j]) * (p_[j] + q_[i]))
def find_equal(a, starts=1500):
    """Damped Newton from many random starts: three KP waves whose three pair factors all equal a."""
    def F(v):
        p_, q_ = v[:3], v[3:]
        return np.array([Akp_num(p_, q_, i, j) - a for i, j in itertools.combinations(range(3), 2)])
    best = None
    for _ in range(starts):
        v = rng.uniform(-3, 3, 6)
        for _it in range(100):
            f = F(v)
            if not np.all(np.isfinite(f)) or np.max(np.abs(f)) < 1e-13: break
            J = np.zeros((3, 6)); h = 1e-8
            for c in range(6):
                e = np.zeros(6); e[c] = h; J[:, c] = (F(v + e) - f) / h
            step = np.linalg.lstsq(J, f, rcond=None)[0]
            nrm = np.linalg.norm(step)
            if nrm > 0.5: step *= 0.5 / nrm
            v = v - step
        f = F(v); p_, q_ = v[:3], v[3:]
        if (np.all(np.isfinite(f)) and np.max(np.abs(f)) < 1e-10 and min(abs(p_[i] + q_[j]) for i in range(3) for j in range(3)) > 1e-3
                and min(abs(p_[i] - p_[j]) + abs(q_[i] - q_[j]) for i, j in itertools.combinations(range(3), 2)) > 1e-3):
            regular = bool(all(p_ + q_ > 0))
            if best is None or (regular and not best[1]): best = (v, regular)
            if regular: break
    return best
reach = []
for a in (0.1, 0.5, 0.9, -0.3, 1.2, 1.5, 2.5):
    b = find_equal(a)
    reach.append((a, b is not None, bool(b[1]) if b else False))
info("three waves, one shared pair factor A:  " + ", ".join(f"A={a}: {'found' if f else 'none'}{' (regular)' if r else ''}" for a, f, r in reach))
check("C2 numerical finding: in KP, three line solitons can share one pair factor only when A > 1 (none for A < 1 in 1500 starts each;"
      " on a line, KdV, never)", all(f for a, f, _ in reach if a > 1) and not any(f for a, f, _ in reach if a < 1))
store["kp_equalA"] = reach

print("\n== D. Belinski-Zakharov gravitational solitons ==")
rho, z = 0.8, 0.3
ws = np.array([-1.2, 0.4, 1.7])
mu = ws - z + np.sqrt((ws - z) ** 2 + rho ** 2)    # pole trajectories
check("D0 poles solve mu^2 + 2(z - w)mu - rho^2 = 0", np.allclose(mu ** 2 + 2 * (z - ws) * mu - rho ** 2, 0))
g0 = np.diag([-1.0, rho ** 2])
def psi0inv(m):
    return np.diag([-1.0, 1.0 / (rho ** 2 - 2 * z * m - m ** 2)])
def bz(m0s):
    ms = [m0s[k] @ psi0inv(mu[k]) for k in range(len(mu))]
    Gam = np.array([[ms[k] @ g0 @ ms[l] / (rho ** 2 + mu[k] * mu[l]) for l in range(len(mu))] for k in range(len(mu))])
    Nv = [g0 @ ms[k] for k in range(len(mu))]
    Gi = np.linalg.inv(Gam)
    g = g0 - sum(np.outer(Nv[k], Nv[l]) * Gi[k, l] / (mu[k] * mu[l]) for k in range(len(mu)) for l in range(len(mu)))
    return Gam, g
Gam1, g1 = bz([np.array([1.0, 0.0])] * 3)
pbz = lambda k, l: -rho ** 2 * (mu[k] - mu[l]) ** 2 / (rho ** 2 + mu[k] * mu[l]) ** 2
ok = all(np.isclose(minor(Gam1, S), math.prod(Gam1[i, i] for i in S) * math.prod(pbz(i, j) for i, j in itertools.combinations(S, 2)))
         for S in subsets(3) if S)
check("D1 one-channel (static, Weyl-type) solitons: Gamma_kl ~ 1/(rho^2 + mu_k mu_l) is Cauchy; every minor is a pair product,"
      " pair factor -rho^2 (mu_k-mu_l)^2/(rho^2+mu_k mu_l)^2", ok)
prodform = [g0[0, 0] * np.prod(mu ** 2 / rho ** 2), (-1) ** len(mu) * g0[0, 0] * np.prod(rho ** 2 / mu ** 2)]
check("D2 ...while the metric itself is a product of one-body factors (size-only, our op's 'no pairs' case):"
      " g_tt = (-1)^n g0_tt * prod (rho/mu_k)^2 (the sign is BZ's known normalisation)", np.isclose(g1[0, 0], prodform[1]) and np.isclose(g1[0, 1], 0),
      f"g_tt {g1[0,0]:.6f}, one-body product {prodform[1]:.6f}")
Gam2, g2 = bz([np.array([1.0, 0.7]), np.array([1.0, -0.4]), np.array([1.0, 1.3])])
tbz = threebody(Gam2)
check("D3 two-channel (rotating, Kerr-type) solitons: Gamma has a rank-2 source, and three of them carry a genuine 3-body term",
      abs(tbz) > 1e-3, f"3-body {tbz:+.3f}")
check("D3b with two of them the corner is still just a pair (Kerr = 2 solitons sits inside the pair picture)",
      np.isclose(minor(Gam2, (0, 1)), Gam2[0, 0] * Gam2[1, 1] * pairf(Gam2, 0, 1)))
store["bz"] = {"threebody_one_channel": threebody(Gam1), "threebody_two_channel": tbz}

print("\n== E. Lucas weights: evaluate det(I + tM) at a 2x2 matrix t = N ==")
def lucasU(kv, dv, m):
    U = [0, 1]
    for _ in range(m): U.append(kv * U[-1] - dv * U[-2])
    return U
Mr = rng.normal(size=(4, 4)) * 0.35
kv, dv = 0.8, -0.5
Nmat = np.array([[kv, 1.0], [-dv, 0.0]])
U = lucasU(kv, dv, 6)
tot = sum(minor(Mr, S) * np.linalg.matrix_power(Nmat, len(S)) for S in subsets(4))
cay = sum(minor(Mr, S) * (U[len(S)] * Nmat - dv * U[len(S) - 1] * np.eye(2) if len(S) else np.eye(2)) for S in subsets(4))
check("E1 sum of det M_S N^|S| = 'det(I + tM) at t = N' = sum det M_S (U_m N - d U_(m-1) I)  (Cayley-Hamilton)", np.allclose(tot, cay))
plane = lambda a, b: (a + b + kv * a * b) / (1 - dv * a * b)
lamr = np.linalg.eigvals(Mr)
pf = 0
for l in lamr: pf = plane(pf, l)
num = sum(U[len(S)] * minor(Mr, S) for S in subsets(4) if S)
den = 1 - dv * sum(U[len(S) - 1] * minor(Mr, S) for S in subsets(4) if len(S) >= 2)
check("E2 the plane law folded over the eigenvalues = Lucas-weighted corners of principal minors (U_m top, -d U_(m-1) bottom)",
      np.isclose(pf, num / den), f"{complex(pf).real:.6f}")
al, be = np.roots([1, -kv, dv])
ta, tb_ = det(np.eye(4) + al * Mr), det(np.eye(4) + be * Mr)
check("E3 the same value is a ratio of two determinants: (det(I+aM) - det(I+bM))/(a det(I+bM) - b det(I+aM)), a+b=k, ab=d",
      np.isclose(pf, (ta - tb_) / (al * tb_ - be * ta)))
Nres = Mr @ np.linalg.inv(np.eye(4) + be * Mr)
cf = fold(np.linalg.eigvals(Nres), al - be)
check("E4 one matrix for the plane: our op with dial a-b on the resolvent M(I + bM)^-1, mapped back by y -> y/(1 - b y)",
      np.isclose(pf, cf / (1 - be * cf)))
# E5 sine-Gordon at the tangent point (0,1) on the KdV-type Cauchy matrix
psg = [sp.Rational(1, 2), sp.Rational(1, 1), sp.Rational(9, 5)]
etasg = [pq * X_ + T_ / pq for pq in psg]
Asg = lambda i, j: ((psg[i] - psg[j]) / (psg[i] + psg[j])) ** 2
def corners_lucas(Ufun, shift):
    return sum(Ufun(len(S) - shift) * sp.prod([Asg(i, j) for i, j in itertools.combinations(S, 2)]) * sp.exp(sum((etasg[i] for i in S), sp.Integer(0)))
               for r in range(4) for S in itertools.combinations(range(3), r) if len(S) >= shift)
Ut = lambda m: [0, 1, 0, -1, 0, 1][m] if m >= 0 else -1   # (k,d)=(0,1), U_(-1) = -1/d
numT = corners_lucas(Ut, 0); denT = 1 - corners_lucas(Ut, 1).subs({}) + 0
denT = 1 - 1 * sum(Ut(len(S) - 1) * sp.prod([Asg(i, j) for i, j in itertools.combinations(S, 2)]) * sp.exp(sum((etasg[i] for i in S), sp.Integer(0)))
                   for r in range(2, 4) for S in itertools.combinations(range(3), r))
phi = 4 * sp.atan(numT / denT)
res = sp.diff(phi, X_, T_) - sp.sin(phi)
rs = max(float(abs(sp.N(res.subs({X_: a, T_: b}), 40))) for a, b in pts)
check("E5 tangent point (0,1) on the soliton Cauchy matrix: phi = 4 arctan(Lucas-weighted corners) solves sine-Gordon phi_xt = sin phi (3 kinks)",
      rs < 1e-25, f"residual {rs:.1e}")
Msg = np.array([[2 * float(psg[i]) / float(psg[i] + psg[j]) * math.exp(float(etasg[i].subs({X_: 0.3, T_: 0.1}))) for j in range(3)] for i in range(3)])
ev_sg = np.linalg.eigvals(Msg).real
phival = float(phi.subs({X_: sp.Rational(3, 10), T_: sp.Rational(1, 10)}))
check("E5b ...and phi/4 = sum of arctan(dressed eigenvalue): the field is additive in the matrix's eigenvalues",
      np.isclose(math.tan(phival / 4), math.tan(np.sum(np.arctan(ev_sg)))))
Uh = lambda m: [0, 1, 0, 1, 0, 1][m] if m >= 0 else 1    # (0,-1), U_(-1) = -1/d = 1
numH = sum(Uh(len(S)) * sp.prod([Asg(i, j) for i, j in itertools.combinations(S, 2)]) * sp.exp(sum((etasg[i] for i in S), sp.Integer(0)))
           for r in range(1, 4) for S in itertools.combinations(range(3), r))
denH = 1 + sum(Uh(len(S) - 1) * sp.prod([Asg(i, j) for i, j in itertools.combinations(S, 2)]) * sp.exp(sum((etasg[i] for i in S), sp.Integer(0)))
               for r in range(2, 4) for S in itertools.combinations(range(3), r))
phih = 4 * sp.atanh(numH / denH)
resh = sp.diff(phih, X_, T_) - sp.sinh(phih)
pts_h = [(sp.Rational(-3, 1), sp.Rational(-1, 1)), (sp.Rational(-4, 1), sp.Rational(1, 2)), (sp.Rational(-5, 2), sp.Rational(-2, 1))]
vals_h = []
for a, b in pts_h:
    Xv = float((numH / denH).subs({X_: a, T_: b}))
    if abs(Xv) < 0.95:
        vals_h.append(float(abs(sp.N(resh.subs({X_: a, T_: b}), 40))))
rh = max(vals_h) if vals_h else 1.0
check("E6 velocity point (0,-1) on the same matrix: phi = 4 artanh(...) solves sinh-Gordon phi_xt = sinh phi (where it is finite)",
      rh < 1e-25 and len(vals_h) >= 2, f"residual {rh:.1e} at {len(vals_h)} points")
# parabolic
kp_, dp_ = 2.0, 1.0
plane_p = lambda a, b: (a + b + kp_ * a * b) / (1 - dp_ * a * b)
pfp = 0
for l in lamr: pfp = plane_p(pfp, l)
yv = np.trace(Mr @ np.linalg.inv(np.eye(4) + Mr))
check("E7 parabolic point (2,1), U_m = m: the plane law maps to plain addition of l/(1+l); total = tr M(I+M)^-1 = d/dt log det(I+tM) at t=1",
      np.isclose(pfp / (1 + pfp), yv))
Ufib = lucasU(1, -1, 8)[1:9]
check("E8 sample Lucas rows: Fibonacci at (1,-1), U_m = k^(m-1) on d = 0 (our size dial)", Ufib == [1, 1, 2, 3, 5, 8, 13, 21] and
      lucasU(0.7, 0, 5)[1:6] == [0.7 ** (m - 1) for m in range(1, 6)] or np.allclose(lucasU(0.7, 0, 5)[1:6], [0.7 ** (m - 1) for m in range(1, 6)]))

print("\n== F. The size x pair candidate as a law ==")
a_, b_, c_, kq, Aq = sp.symbols('a b c k A')
two = a_ + b_ + kq * Aq * a_ * b_
three_comb = a_ + b_ + c_ + kq * Aq * (a_ * b_ + a_ * c_ + b_ * c_) + kq ** 2 * Aq ** 3 * a_ * b_ * c_
uu = sp.Symbol('u'); iter3 = sp.expand(two.subs({a_: uu, b_: c_}, simultaneous=True).subs(uu, two))
diff = sp.factor(sp.expand(three_comb - iter3))
info(f"two inputs: a + b + kA ab (our op at dial kA).  three inputs minus iterating that: {diff}")
check("F1 with equal A it is NOT an iterated binary law (not a formal group law) unless A = 0 or 1: triple corner k^2 A^3 vs k^2 A^2",
      sp.simplify(diff - kq ** 2 * Aq ** 2 * (Aq - 1) * a_ * b_ * c_) == 0)
check("F2 hard core A = 0 gives plain addition; A = 1 gives our op: the pair dial runs between them, the same two ends as k",
      sp.simplify(three_comb.subs(Aq, 0) - (a_ + b_ + c_)) == 0 and sp.simplify(three_comb.subs(Aq, 1) - sp.expand((( 1 + kq*a_)*(1+kq*b_)*(1+kq*c_) - 1)/kq)) == 0)
# Mayer: log tau minus one-body logs = connected clusters
f = sp.Symbol('f'); xa, xb, xc = sp.symbols('x1 x2 x3'); ee = sp.Symbol('e')
tau3 = 1 + kq * ee * (xa + xb + xc) + kq ** 2 * ee ** 2 * (1 + f) * (xa * xb + xa * xc + xb * xc) + kq ** 3 * ee ** 3 * (1 + f) ** 3 * xa * xb * xc
excess = sp.series(sp.log(tau3) - sp.log(1 + kq * ee * xa) - sp.log(1 + kq * ee * xb) - sp.log(1 + kq * ee * xc), ee, 0, 4).removeO()
e2 = sp.simplify(excess.coeff(ee, 2) - kq ** 2 * f * (xa * xb + xa * xc + xb * xc))
e3 = sp.factor(excess.coeff(ee, 3))
info(f"third-order excess (Mayer clusters on 3 sites) = {e3}")
check("F3 Mayer/cluster reading: A = 1 + f; log tau minus the one-body logs (our flat coordinate) starts with k^2 f sum x_i x_j",
      e2 == 0)

print("\n== G. Bridge to the cube-weights page: reordering factors (bosons) = one determinant (fermions) ==")
zz = rng.uniform(0.5, 3, 4); ww = -rng.uniform(0.2, 2.5, 4)
boson = math.prod((zz[i] - zz[j]) * (ww[i] - ww[j]) / ((zz[i] - ww[j]) * (ww[i] - zz[j])) for i, j in itertools.combinations(range(4), 2)) \
        * math.prod(1 / (zz[i] - ww[i]) for i in range(4))
fermion = det(1 / (zz[:, None] - ww[None, :]))
check("G1 vertex operators (bosons): product of one reordering factor per pair  =  Wick (fermions): det of ONE matrix 1/(z_i - w_j)",
      np.isclose(boson, fermion), f"{boson:.6e} vs {fermion:.6e}  (Cauchy's determinant: the boson-fermion dictionary)")
pairf_b = lambda i, j: (zz[i] - zz[j]) * (ww[i] - ww[j]) / ((zz[i] - ww[j]) * (ww[i] - zz[j]))
pp_, qq_ = zz, -ww
check("G2 with z = p, w = -q the reordering factor is exactly the KP pair factor (and with q = p, KdV's)",
      all(np.isclose(pairf_b(i, j), (pp_[i] - pp_[j]) * (qq_[i] - qq_[j]) / ((pp_[i] + qq_[j]) * (pp_[j] + qq_[i]))) for i, j in itertools.combinations(range(4), 2)))
check("G3 the matrix side needs no infinite operators: a real 4x4 Cauchy matrix gives real pair factors strictly between 0 and 1",
      all(0 < A(i, j) < 1 for i, j in itertools.combinations(range(4), 2)))

print("\n== H. Open question: does 'no three-body term in any triple' force 'pairs only' everywhere? ==")
idx4 = list(itertools.combinations(range(4), 2))
def mat4(c):
    Mx = np.eye(4)
    for v, (i, j) in zip(c, idx4): Mx[i, j] = Mx[j, i] = v
    return Mx
def tri4(c):
    Mx = mat4(c); out = []
    for i, j, l in itertools.combinations(range(4), 3):
        a1, b1, c1 = Mx[i, j], Mx[i, l], Mx[j, l]
        out.append(2*a1*b1*c1 - (a1*a1*b1*b1 + a1*a1*c1*c1 + b1*b1*c1*c1) + a1*a1*b1*b1*c1*c1)
    return np.array(out)
def four4(c):
    Mx = mat4(c); return det(Mx) - math.prod(1 - Mx[i, j]**2 for i, j in idx4)
sols = []
for _ in range(800):
    c = rng.uniform(-0.95, 0.95, 6)
    for _it in range(80):
        fv = tri4(c)
        if np.max(np.abs(fv)) < 1e-14: break
        J = np.zeros((4, 6)); h = 1e-8
        for kk2 in range(6):
            e = np.zeros(6); e[kk2] = h; J[:, kk2] = (tri4(c + e) - fv) / h
        st = np.linalg.lstsq(J, fv, rcond=None)[0]
        if np.linalg.norm(st) > 0.3: st *= 0.3 / np.linalg.norm(st)
        c = c - st
    if np.max(np.abs(tri4(c))) < 1e-12 and np.all(np.abs(c) < 5) and np.min(np.abs(c)) > 1e-3:
        sols.append(abs(four4(c)))
sols = np.array(sols)
nz = int((sols > 1e-8).sum())
info(f"symmetric 4x4, unit diagonal, every triple free of 3-body terms: {len(sols)} solutions, {len(sols)-nz} with no 4-body term, {nz} WITH one (max {sols.max():.3f})")
check("H1 numerical: 'no 3-body term in any triple' does NOT force 'no 4-body term'; pairs-only is a global property (two branches)",
      nz > 0 and len(sols) - nz > 0)
store["h1"] = {"solutions": len(sols), "with_four_body": nz}

print("\n== I. Face weights from vertex operators: each wave = (multiply) x (shift); a pair weight is the cost of one shift passing one multiply ==")
def vo_tau(ps, a, J, kind="kdv", qs=None, beta=0.0):
    """Apply (1 + a_1 X_1)...(1 + a_n X_n) to the function 1, leftmost operator applied last.
    Functions are sums c * exp(b . t) over modes j = 1..J. X = exp(xi . t) * exp(-eps . d/dt [+ beta (eps . d/dt)^2])."""
    js = np.arange(1, J + 1)
    def parts(i):
        if kind == "kdv":
            odd = (js % 2 == 1)
            xi = np.where(odd, 2.0 * ps[i] ** js, 0.0)
            eps = np.where(odd, 2.0 / (js * ps[i] ** js), 0.0)
        else:
            xi = ps[i] ** js - (-qs[i]) ** js
            eps = (ps[i] ** (-js) - (-qs[i]) ** (-js)) / js
        return xi, eps
    terms = {(): (1.0, np.zeros(J))}
    for i in reversed(range(len(ps))):
        xi, eps = parts(i); new = dict(terms)
        for S, (c, b) in terms.items():
            e = eps @ b
            new[(i,) + S] = (a[i] * c * math.exp(-e + beta * e * e), b + xi)
        terms = new
    return {tuple(sorted(S)): c for S, (c, b) in terms.items()}
ps = [2.0, 1.2, 0.5]; av = [0.7, 1.3, 0.9]
Akdv = lambda i, j: ((ps[i] - ps[j]) / (ps[i] + ps[j])) ** 2
conv = []
for J in (5, 21, 81, 321):
    co = vo_tau(ps, av, J)
    err = max(abs(co[S] / (math.prod(av[i] for i in S) * math.prod(Akdv(i, j) for i, j in itertools.combinations(S, 2))) - 1) for S in co if len(S) >= 2)
    conv.append((J, err))
info("KdV vertex operators, truncated to J modes: worst corner error " + ", ".join(f"J={J}: {e:.1e}" for J, e in conv))
check("I1 KdV vertex operators X(p) = exp(sum 2p^j t_j) exp(-sum (2/(j p^j)) d_j) (odd j): every corner of (1+aX)(1+aX)(1+aX)1 is prod a * prod pairs ((p_i-p_j)/(p_i+p_j))^2",
      conv[-1][1] < 1e-9 and conv[0][1] > conv[-1][1])
x_ = 0.6
s_odd = sum(4 * x_ ** j / j for j in range(1, 2000, 2))
check("I1b the one-pair cost is exp(-4 artanh(p'/p)) = ((p-p')/(p+p'))^2: artanh, the flat coordinate of the (0,-1) velocity law",
      np.isclose(math.exp(-s_odd), ((1 - x_) / (1 + x_)) ** 2) and np.isclose(s_odd, 4 * math.atanh(x_)))
ex = []
for gap in (0.3, 0.1, 0.03, 0.01):
    co = vo_tau([1.0, 1.0 - gap], [1, 1], min(4001, int(60 / gap) + 1))  # enough modes, no overflow
    ex.append((gap, co[(0, 1)]))
info("two waves closing in on the same speed, pair cost: " + ", ".join(f"gap {g}: {c:.2e}" for g, c in ex))
check("I1c equal speeds cost everything (pair weight -> 0): X(p) X(p) = 0, exclusion again", ex[-1][1] < 1e-3 and all(ex[i + 1][1] < ex[i][1] for i in range(3)) and all(np.isclose(c, (g / (2 - g)) ** 2, rtol=1e-6) for g, c in ex))
store["vo_conv"] = conv; store["vo_exclusion"] = ex
pk2 = [2.2, 1.6, 0.9]; qk2 = [2.6, 1.3, 0.7]
co = vo_tau(pk2, [1, 1, 1], 400, kind="kp", qs=qk2)
Akp2 = lambda i, j: (pk2[i] - pk2[j]) * (qk2[i] - qk2[j]) / ((pk2[i] + qk2[j]) * (pk2[j] + qk2[i]))
errk = max(abs(co[S] / math.prod(Akp2(i, j) for i, j in itertools.combinations(S, 2)) - 1) for S in co if len(S) >= 2)
check("I2 KP vertex operators X(p,q) = exp(sum (p^j - (-q)^j) t_j) exp(-sum (p^-j - (-q)^-j) d_j / j): corners carry the KP pair factors, three waves",
      errk < 1e-9, f"worst error {errk:.1e}")
ps4 = [3.0, 2.1, 1.4, 0.6]; a4 = [1, 1, 1, 1]
co4 = vo_tau(ps4, a4, 600)
A4 = lambda i, j: ((ps4[i] - ps4[j]) / (ps4[i] + ps4[j])) ** 2
e4 = max(abs(co4[S] / math.prod(A4(i, j) for i, j in itertools.combinations(S, 2)) - 1) for S in co4 if len(S) >= 2)
Mvo = np.array([[2 * ps4[i] / (ps4[i] + ps4[j]) for j in range(4)] for i in range(4)])
same = all(np.isclose(co4[S], minor(Mvo, S)) for S in co4 if S)
check("I2b four waves: vertex-operator corners (16 of them) = principal minors of the Cauchy matrix, one by one (face picture = matrix picture)",
      e4 < 1e-9 and same)
coI = vo_tau(ps, av, 321, beta=0.05)
tbI = coI[(0, 1, 2)] * coI[()] ** 2 * av[0] * av[1] * av[2] / (coI[(0, 1)] * coI[(0, 2)] * coI[(1, 2)]) - 1
check("I3 make the field interacting (a squared term in the shift) and a 3-body term appears (free field gives pairs only; this one interaction breaks it; the converse is not proved)",
      abs(tbI) > 1e-3, f"3-body {tbI:+.3f}")
store["vo_interacting_threebody"] = tbI
sweep = []
for bb in np.linspace(0.0, 0.08, 9):
    cb = vo_tau(ps, av, 321, beta=float(bb))
    sweep.append((float(bb), cb[(0, 1, 2)] * cb[()] ** 2 * av[0] * av[1] * av[2] / (cb[(0, 1)] * cb[(0, 2)] * cb[(1, 2)]) - 1))
store["vo_beta_sweep"] = sweep
# Weyl gravity: psi is a free (harmonic) field; gamma is quadratic in psi -> pairs only
zs, ms_ = [-1.0, 0.4, 1.5], [0.8, 0.5, 1.1]
def grad_gamma(sel, rho, z, h=1e-5):
    psi = lambda r, zz: sum(-ms_[i] / math.hypot(r, zz - zs[i]) for i in sel)
    pr = (psi(rho + h, z) - psi(rho - h, z)) / (2 * h); pz = (psi(rho, z + h) - psi(rho, z - h)) / (2 * h)
    return np.array([rho * (pr * pr - pz * pz), 2 * rho * pr * pz])
r0, z0 = 0.7, 0.2
inter3 = grad_gamma([0, 1, 2], r0, z0) - sum(grad_gamma(list(P), r0, z0) for P in itertools.combinations(range(3), 2)) + sum(grad_gamma([i], r0, z0) for i in range(3))
check("I4 gravity, static (Weyl): the potential psi superposes (free field, our 'no pairs' layer) and the conformal factor gamma is built only from pairs:"
      " its 3-body part is zero", np.max(np.abs(inter3)) < 1e-6, f"3-body part of grad gamma {np.max(np.abs(inter3)):.1e}")

npass = sum(ok for _, ok in results)
print(f"\n{npass} of {len(results)} checks pass")
json.dump({k_: v for k_, v in store.items()}, open("results.json", "w"), indent=1, default=str)
