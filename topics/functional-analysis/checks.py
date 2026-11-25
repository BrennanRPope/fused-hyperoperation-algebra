"""Functional analysis through mult-add: every claim on the page, checked.
Run: python3 checks.py > checks_output.txt
"""
import math, json
import numpy as np
import scipy.linalg as sla
import sympy as sp
from numpy.polynomial.legendre import leggauss
from scipy.stats import unitary_group
from scipy.integrate import quad

rng = np.random.default_rng(11)
results = []
def check(name, ok, detail=""):
    results.append(bool(ok))
    print(("PASS " if ok else "FAIL ") + name + (("   " + detail) if detail else ""))

I = lambda n: np.eye(n)
def ma(A, B): return A + B + A @ B
def qi(A): return np.linalg.inv(I(len(A)) + A) - I(len(A))
def rmat(n, s=1.0): return s * (rng.normal(size=(n, n)) + 1j * rng.normal(size=(n, n))) / math.sqrt(2 * n)
nrm = lambda A: np.linalg.norm(A, 2)

# ---------------- 1. operator mult-add
A, B, C = rmat(6), rmat(6), rmat(6)
check("A * B = A + B + AB = (I + A)(I + B) - I", nrm(I(6) + ma(A, B) - (I(6) + A) @ (I(6) + B)) < 1e-12)
check("associative: (A*B)*C = A*(B*C)", nrm(ma(ma(A, B), C) - ma(A, ma(B, C))) < 1e-12)
check("0 is the identity, -I is the wall: (-I)*B = B*(-I) = -I",
      nrm(ma(A, 0 * A) - A) < 1e-15 and nrm(ma(-I(6), B) + I(6)) < 1e-12 and nrm(ma(B, -I(6)) + I(6)) < 1e-12)
check("order matters: A*B != B*A for random matrices", nrm(ma(A, B) - ma(B, A)) > 0.1)
check("but det and trace do not see the order: det(I + A*B) = det(I + B*A) = det(I+A)det(I+B), tr(A*B) = tr(B*A)",
      abs(np.linalg.det(I(6) + ma(A, B)) - np.linalg.det(I(6) + A) * np.linalg.det(I(6) + B)) < 1e-10
      and abs(np.linalg.det(I(6) + ma(A, B)) - np.linalg.det(I(6) + ma(B, A))) < 1e-10
      and abs(np.trace(ma(A, B)) - np.trace(ma(B, A))) < 1e-12)

# ---------------- 2. quasi-inverses (Perlis 1942: r + s + rs = 0)
Ap = qi(A)
check("quasi-inverse A' = (I+A)^-1 - I satisfies A*A' = A'*A = 0 (Perlis's r + s + rs = 0)",
      nrm(ma(A, Ap)) < 1e-12 and nrm(ma(Ap, A)) < 1e-12)
check("A' = -A(I+A)^-1 = -(I+A)^-1 A",
      nrm(Ap + A @ np.linalg.inv(I(6) + A)) < 1e-12 and nrm(Ap + np.linalg.inv(I(6) + A) @ A) < 1e-12)
check("Jacobson's circle a o b = a + b - ab is the k = -1 dial: a o b = -((-a)*(-b))",
      nrm((A + B - A @ B) + ma(-A, -B)) < 1e-12)
Bp = qi(B)
check("socks and shoes: (A*B)' = B'*A'", nrm(qi(ma(A, B)) - ma(Bp, Ap)) < 1e-10)
check("resolvent identity in mult-add form: A' - B' = (I + A')(B - A)(I + B')",
      nrm(Ap - Bp - (I(6) + Ap) @ (B - A) @ (I(6) + Bp)) < 1e-10)
T = rmat(6, 1.3)
lam, mu = 0.7 + 0.4j, -0.5 + 1.1j
R = lambda z: np.linalg.inv(T - z * I(6))
check("classical Hilbert resolvent identity R(l) - R(m) = (l - m) R(l) R(m), R(z) = (T - z)^-1",
      nrm(R(lam) - R(mu) - (lam - mu) * R(lam) @ R(mu)) < 1e-10)
check("and R(z) = -(1/z)(I + (-T/z)'): the resolvent is a rescaled quasi-inverse",
      nrm(R(lam) + (I(6) + qi(-T / lam)) / lam) < 1e-10)

# ---------------- 3. the wall in the spectrum, Neumann series
Aw = rmat(5, 1.5)
ev = np.linalg.eigvals(Aw)
check("spectrum of I + A is 1 + spectrum of A", np.allclose(np.sort_complex(np.linalg.eigvals(I(5) + Aw)), np.sort_complex(1 + ev)))
S = rmat(5); v = np.linalg.eig(S)
Ah = S @ np.diag([-1, 0.3, 0.5j, 2, -0.4]) @ np.linalg.inv(S)
check("I + A is singular exactly when -1 is an eigenvalue of A (the wall in the spectrum)",
      abs(np.linalg.det(I(5) + Ah)) < 1e-10 and abs(np.linalg.det(I(5) + Ah * 0.9)) > 1e-3)
A8 = rmat(6); A8 = 0.8 * A8 / nrm(A8)
s, errs = np.zeros_like(A8), []
P = I(6)
for n in range(1, 41):
    P = P @ (-A8); s = s + P
    errs.append(nrm(ma(A8, s)))
check("Neumann: for ||A|| = 0.8 the partial sums of sum (-A)^n converge to the quasi-inverse", nrm(s - qi(A8)) < 1e-3)
s2, P2, ok = np.zeros_like(A8), I(6), True
for n in range(1, 12):
    P2 = P2 @ (-A8); s2 = s2 + P2
    ok &= nrm(ma(A8, s2) + np.linalg.matrix_power(-A8, n + 1)) < 1e-12
check("each partial sum misses the wall-free answer by exactly one term: A * s_n = -(-A)^(n+1)", ok)
check("so the miss shrinks at least like 0.8^n", all(errs[n - 1] <= 0.8 ** (n + 1) + 1e-12 for n in range(1, 41)))
Abig = 2.0 * I(4)
check("the ball ||A|| < 1 is safe but not the whole safe region: A = 2I is far out yet quasi-invertible",
      nrm(ma(Abig, qi(Abig))) < 1e-12)
check("the distance from 0 to the wall -I is exactly 1 in operator norm", abs(nrm(-I(6)) - 1) < 1e-15)

# ---------------- 4. spectra: commuting vs not
X = rmat(5); Ac = X @ X + 0.3 * X; Bc = 0.5 * X - 0.2 * X @ X @ X
e = np.linalg.eigvals(X)
pa, pb = e * e + 0.3 * e, 0.5 * e - 0.2 * e ** 3
got = np.sort_complex(np.linalg.eigvals(ma(Ac, Bc)))
check("commuting A, B (polynomials in one X): eigenvalues of A*B are the pairwise mult-adds a_i * b_i",
      np.allclose(got, np.sort_complex(pa + pb + pa * pb), atol=1e-9))
N1 = np.array([[0, 1], [0, 0]]); N2 = np.array([[0, 0], [-1, 0]])
ok = True; pts = []
for t in np.linspace(0.1, 3.9, 20):
    M = ma(N1, t * N2); evs = np.linalg.eigvals(M)
    ok &= np.allclose(abs(1 + evs), 1) and abs(np.linalg.det(np.eye(2) + M) - 1) < 1e-12
check("two operators with spectrum {0} (nilpotent) can mult-add to eigenvalues all round the circle |1+z| = 1", ok)
ev6 = np.linalg.eigvals(ma(N1, N2))
check("at t = 1 those eigenvalues are e^{+-2 pi i/3}: N1*N2 is a sixth mult-add root of 0 (its sixth mult-add power vanishes)",
      np.allclose(np.sort(np.angle(ev6)), [-2 * np.pi / 3, 2 * np.pi / 3])
      and np.allclose(np.linalg.matrix_power(np.eye(2) + ma(N1, N2), 6) - np.eye(2), 0))
ok = True
for _ in range(200):
    Un = np.triu(rng.normal(size=(4, 4)), 1); Ln = np.tril(rng.normal(size=(4, 4)), -1)
    Mn = ma(Un, Ln)
    ok &= abs(np.linalg.det(np.eye(4) + Mn) - 1) < 1e-9 and np.min(abs(np.linalg.eigvals(Mn) + 1)) > 1e-6
check("they never reach the wall: det(I + A*B) = det(I+A)det(I+B) = 1 for nilpotent A, B (200 random 4x4 pairs)", ok)
# characters: circulants
n = 7
def circ(c): return np.array([np.roll(c, k) for k in range(n)]).T
ca, cb = rng.normal(size=n), rng.normal(size=n)
Ca, Cb = circ(ca), circ(cb)
fa, fb = np.fft.fft(ca), np.fft.fft(cb)
fab = np.fft.fft(ma(Ca, Cb)[:, 0])
check("characters carry mult-add to mult-add: Fourier modes of circulants, F(A*B) = F(A) * F(B) mode by mode",
      np.allclose(fab, fa + fb + fa * fb))
# BCH
eps = 1e-2
Xs, Ys = rmat(4) * eps, rmat(4) * eps
Am, Bm = sla.expm(Xs) - I(4), sla.expm(Ys) - I(4)
Lg = sla.logm(I(4) + ma(Am, Bm))
r2 = nrm(Lg - Xs - Ys); r3 = nrm(Lg - Xs - Ys - 0.5 * (Xs @ Ys - Ys @ Xs))
check("log of the lift turns * into + only for commuting operators; otherwise the leftover starts with (1/2)[X, Y] (BCH)",
      r2 > 10 * r3 and r3 < 50 * eps ** 3, f"{r2:.2e} vs {r3:.2e}")
Xc = rmat(4); Yc = 0.3 * Xc @ Xc
check("commuting X, Y: (e^X - I) * (e^Y - I) = e^(X+Y) - I exactly",
      nrm(ma(sla.expm(Xc) - I(4), sla.expm(Yc) - I(4)) - (sla.expm(Xc + Yc) - I(4))) < 1e-10)
k = 0.37
check("dial on operators: det(I + k(A *_k B)) = det(I + kA) det(I + kB)",
      abs(np.linalg.det(I(6) + k * (A + B + k * A @ B)) - np.linalg.det(I(6) + k * A) * np.linalg.det(I(6) + k * B)) < 1e-10)

# ---------------- 5. Volterra: every multiple stays off the wall
m = 400; h = 1.0 / m; xs = (np.arange(m) + 0.5) * h
V = np.tril(np.ones((m, m)), -1) * h + np.eye(m) * h / 2
vn = []
P = np.eye(m)
for n in range(1, 13):
    P = P @ V; vn.append(np.linalg.norm(P, 2) ** (1 / n))
check("Volterra integration V: ||V^n||^(1/n) falls toward 0 (spectral radius 0), under the bound (1/n!)^(1/n) ~ e/n",
      all(vn[i + 1] < vn[i] for i in range(11)) and all(vn[n - 1] <= (1 / math.factorial(n)) ** (1 / n) * 1.01 for n in range(1, 13)), " ".join(f"{x:.3f}" for x in vn[:6]))
check("||V^n|| <= 1/n! (n = 1..8, discretised)",
      all(np.linalg.norm(np.linalg.matrix_power(V, n), 2) <= 1 / math.factorial(n) * 1.01 for n in range(1, 9)))
check("||V|| = 2/pi", abs(np.linalg.norm(V, 2) - 2 / math.pi) < 1e-4)
ok = True
fx = np.cos(3 * xs)
for lamb in (5.0, -5.0, 12.0):
    y = np.linalg.solve(np.eye(m) + lamb * V, fx)
    ex = np.array([math.exp(-lamb * x) + quad(lambda t: math.exp(-lamb * (x - t)) * (-3 * math.sin(3 * t)), 0, x)[0] for x in xs[::40]])
    ok &= np.max(abs(y[::40] - ex)) < 2e-4 * np.max(abs(ex))
check("I + lambda V is invertible for lambda = 5, -5, 12: solution of y + lambda int y = f matches the exact ODE answer", ok)
t, x, s, L = sp.symbols('t x s lambda')
n_ = sp.symbols('n', integer=True, positive=True)
ker = sp.summation((-L) ** n_ * (x - t) ** (n_ - 1) / sp.factorial(n_ - 1), (n_, 1, sp.oo))
check("Neumann series for lambda V converges for every lambda: quasi-inverse kernel -lambda e^{-lambda (x-t)}",
      sp.simplify(ker + L * sp.exp(-L * (x - t))) == 0)
Q = lambda xx, tt: -lamb * np.exp(-lamb * (xx - tt)) * (tt < xx)
lamb = 3.0
Qm = (Q(xs[:, None], xs[None, :]) + np.diag(np.full(m, -lamb / 2))) * h
check("that kernel really is the quasi-inverse: (lambda V) * Q = 0 up to discretisation",
      nrm(ma(lamb * V, Qm)) < 2e-3)

# ---------------- 6. Fredholm determinant: the fold of the eigenvalues
def nystrom(K, mm, a=0.0, b=1.0):
    xg, wg = leggauss(mm); xg = (b - a) / 2 * xg + (a + b) / 2; wg = wg * (b - a) / 2
    sw = np.sqrt(wg); return xg, wg, sw[:, None] * K(xg[:, None], xg[None, :]) * sw[None, :]
G = lambda x, y: np.minimum(x, y) - x * y
xg, wg, Gm = nystrom(G, 400)
eg = np.sort(np.linalg.eigvalsh(Gm))[::-1]
check("string Green kernel min(x,y) - xy: eigenvalues 1/(n^2 pi^2)", np.allclose(eg[:5], 1 / (np.arange(1, 6) * np.pi) ** 2, rtol=3e-4))
ok = True
for z in (0.7, 2.3, 4.0, 7.5):
    ok &= abs(np.linalg.det(np.eye(400) - z * z * Gm) - np.sin(z) / z) < 2e-5
check("det(I - z^2 G) = sin(z)/z (Nystrom, 400 nodes): the sine product of the complex page", ok)
z = 2.3
check("Lidskii: det(I + K) - 1 is the mult-add fold of the eigenvalues of K",
      abs((np.linalg.det(np.eye(400) - z * z * Gm) - 1) - (np.prod(1 - z * z * eg) - 1)) < 1e-10)
check("the determinant vanishes exactly when an eigenvalue of -z^2 G hits the wall: z = n pi",
      all(abs(-(n * np.pi) ** 2 * eg[n - 1] + 1) < 1e-3 for n in (1, 2, 3)) and abs(np.sin(np.pi) / np.pi) < 1e-15)
ok = all(np.sign(np.prod(1 - zz * zz * eg)) == (-1) ** int(zz // np.pi) for zz in (1.0, 4.0, 7.0, 10.0, 13.0))
check("its sign flips each time one more eigenvalue crosses the wall: sign = (-1)^(number past -1)", ok)
tr1 = np.trace(Gm); tr2 = np.trace(Gm @ Gm)
check("cube layers: e1 = tr G = 1/6 and e2 = (tr^2 - tr G^2)/2 = 1/120 match sin z/z = 1 - z^2/6 + z^4/120",
      abs(tr1 - 1 / 6) < 1e-5 and abs((tr1 ** 2 - tr2) / 2 - 1 / 120) < 1e-5)
K1 = lambda x, y: 0.8 * np.exp(-(x - y) ** 2); K2 = lambda x, y: 0.5 * np.cos(x + 2 * y)
_, w1, M1 = nystrom(K1, 60); _, _, M2 = nystrom(K2, 60)
check("Fredholm determinants multiply along mult-add of integral operators: det(I + K*L) = det(I+K) det(I+L)",
      abs(np.linalg.det(np.eye(60) + ma(M1, M2)) - np.linalg.det(np.eye(60) + M1) * np.linalg.det(np.eye(60) + M2)) < 1e-10)
_, _, M1b = nystrom(K1, 120); _, _, M2b = nystrom(K2, 120)
check("and the Nystrom values are converged (60 vs 120 nodes agree to 1e-12)",
      abs(np.linalg.det(np.eye(60) + ma(M1, M2)) - np.linalg.det(np.eye(120) + ma(M1b, M2b))) < 1e-12)

# ---------------- 7. determinantal points: the k = -1 dial
N = 8; a = 0.6
def Kc(t, s):
    d = t - s
    with np.errstate(invalid='ignore', divide='ignore'):
        v = np.sin(N * d / 2) / (2 * np.pi * np.sin(d / 2))
    return np.where(np.abs(d) < 1e-12, N / (2 * np.pi), v)
_, _, Mc = nystrom(Kc, 40, -a, a)
lc = np.sort(np.linalg.eigvalsh(Mc))[::-1]
check("CUE(8) kernel on the arc |theta| < 0.6: eigenvalues lie in [0, 1]", lc.min() > -1e-12 and lc.max() < 1)
def wfold(k, xs_):
    t = 0.0
    for v_ in xs_: t = t + v_ + k * t * v_
    return t
pnone = np.linalg.det(np.eye(40) - Mc)
check("P(at least one eigenvalue in the arc) = 1 - det(I - K) = lambda_1 *_{-1} lambda_2 *_{-1} ... (failure union)",
      abs((1 - pnone) - wfold(-1, lc)) < 1e-12)
pc = np.array([1.0])
for l in lc[:15]: pc = np.convolve(pc, [1 - l, l])
S = 40000
U = unitary_group.rvs(N, size=S, random_state=np.random.default_rng(1))
cnt = np.sum(np.abs(np.angle(np.linalg.eigvals(U))) < a, axis=1)
emp = np.bincount(cnt, minlength=6)[:4] / S
sig = np.sqrt(pc[:4] * (1 - pc[:4]) / S)
check("Monte Carlo, 40000 random unitaries: counts in the arc match independent coin flips with chances lambda_i (within 3 sigma)",
      np.all(abs(emp - pc[:4]) < 3 * sig + 1e-4), " ".join(f"{e_:.4f}/{p_:.4f}" for e_, p_ in zip(emp, pc[:4])))
void = lambda k: np.prod((1 + k * lc) ** (-1 / k))
check("the dial interpolates: void chance prod(1 + k lambda)^(-1/k) is det(I-K) at k = -1",
      abs(void(-1.0) - pnone) < 1e-12)
check("and tends to e^(-sum lambda) (Poisson, no repulsion) as k -> 0", abs(void(1e-7) - math.exp(-lc.sum())) < 1e-6)
check("void chance rises from fermions (k=-1) through Poisson (k=0) to bosons (k=+1, det(I+K)^-1)",
      void(-1) < void(1e-9) < void(1.0) and abs(void(1.0) - 1 / np.linalg.det(np.eye(40) + Mc)) < 1e-12,
      f"{void(-1):.4f} < {void(1e-9):.4f} < {void(1.0):.4f}")
check("lift form: void = (1 + k F_k)^(-1/k) where F_k is the k-dial fold of the eigenvalues",
      all(abs(void(kk) - (1 + kk * wfold(kk, lc)) ** (-1 / kk)) < 1e-12 for kk in (-1.0, -0.5, 0.5, 1.0)))

n_ = len(results); p_ = sum(results)
print(f"\n{p_} of {n_} checks pass")
json.dump({"vn": vn, "neumann_miss": errs[:30], "cue_lambda": lc[:6].tolist(), "cue_pred": pc[:5].tolist(),
           "cue_emp": (np.bincount(cnt, minlength=6)[:5] / S).tolist(), "void": [void(-1), void(1e-9), void(1.0)]},
          open("results.json", "w"), indent=1)
