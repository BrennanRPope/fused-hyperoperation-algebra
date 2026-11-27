"""DRAFT round 8 checks: all of Einstein's equations (and Maxwell's) from varying the loop action.

Round 6 varied the loop-read curvature only for a uniform expanding universe. Here the metric is generic:
all 10 components of a 4D metric vary from place to place. The curvature scalar is read from loops
(chains of short links joined by the law a o b = a + b + ab, later steps on the left) at every point of a
grid, integrated into Einstein's action S = integral sqrt(g) R, and each of the 10 metric components is
nudged by a smooth bump. The change in S is compared with Einstein's tensor (textbook formula, computed
independently by finite differences): delta S = - integral sqrt(g) G^{mu nu} delta g_{mu nu}.
Then the 5D chain of round 3: varying the same loop action in the electromagnetic potential gives Maxwell.
Euclidean signature (the algebra is the same as Lorentzian). Textbook: Hilbert 1915; Kaluza 1921.
Run: python3 checks.py  (numpy, scipy, sympy).  Writes results.json and checks_output.txt.
"""
import json, io, math, itertools
import numpy as np
from scipy.linalg import expm
import sympy as sp

results, passed, failed = {}, 0, 0
_log = io.StringIO()
def say(s):
    print(s); _log.write(s + "\n")
def check(name, ok, detail=""):
    global passed, failed
    passed, failed = (passed + 1, failed) if ok else (passed, failed + 1)
    say(("PASS " if ok else "FAIL ") + name + (("  | " + detail) if detail else ""))

X, Y, D = sp.symbols("x y d", real=True)
W = 0.12
BUMP = sp.exp(-((X - 0.5)**2 + (Y - 0.5)**2) / W**2)

def compile_metric(gsym):
    """metric depending on x, y (and the nudge size d); returns batched g(x,y,d) and its x, y derivatives"""
    n = gsym.shape[0]
    def entrywise(mat):
        fs = [[sp.lambdify((X, Y, D), mat[i, j], "numpy") for j in range(n)] for i in range(n)]
        def ev(x, y, d):
            return np.stack([np.stack([np.broadcast_to(np.asarray(fs[i][j](x, y, d), float), x.shape)
                                       for j in range(n)], -1) for i in range(n)], -2)
        return ev
    return n, entrywise(gsym), entrywise(gsym.diff(X)), entrywise(gsym.diff(Y))

def gamma_batch(M, P, d):
    """Christoffel symbols at points P (shape [m, n]); metric depends on coordinates 1 and 2 (x, y)"""
    n, g, gx, gy = M
    x, y = P[:, 1], P[:, 2]
    G = g(x, y, d); Gi = np.linalg.inv(G)
    dG = np.zeros(G.shape + (n,)); dG[..., 1] = gx(x, y, d); dG[..., 2] = gy(x, y, d)   # dG[m, i, j, l] = d_l g_ij
    lower = 0.5 * (dG + np.transpose(dG, (0, 1, 3, 2)) - np.transpose(dG, (0, 3, 1, 2)))
    return np.einsum("mad,mdbc->mabc", Gi, lower), G, Gi

def holonomy_batch(M, P, mu, nu, eps, d, sub=6):
    n = M[0]; m = P.shape[0]
    out = np.zeros((m, n, n))
    p = P.copy(); p[:, mu] -= eps / 2; p[:, nu] -= eps / 2
    for dirn, sign in ((mu, +1), (nu, +1), (mu, -1), (nu, -1)):
        dd = eps / sub
        for s in range(sub):
            mid = p.copy(); mid[:, dirn] += sign * dd * (s + 0.5)
            Gm, _, _ = gamma_batch(M, mid, d)
            v = expm(-Gm[:, :, dirn, :] * sign * dd) - np.eye(n)
            out = v + out + v @ out                       # the law, later steps on the left
        p[:, dirn] += sign * eps
    return out

def scalar_from_loops(M, P, d, eps=2e-3):
    n = M[0]; m = P.shape[0]
    Rl = np.zeros((m, n, n, n, n))
    for mu, nu in itertools.combinations(range(n), 2):
        H = (8 * holonomy_batch(M, P, mu, nu, eps / 2, d) - holonomy_batch(M, P, mu, nu, eps, d)) / eps**2
        Rl[:, :, :, mu, nu] = -H; Rl[:, :, :, nu, mu] = H
    _, G, Gi = gamma_batch(M, P, d)
    return np.einsum("mbd,mabad->m", Gi, Rl), Rl, G, Gi

def ricci_formula(M, P, d, h=1e-4):
    """textbook Ricci tensor R_bd = d_a Gam^a_db - d_d Gam^a_ab + Gam^a_ae Gam^e_db - Gam^a_de Gam^e_ab, by finite differences"""
    n = M[0]
    Gm, G, Gi = gamma_batch(M, P, d)
    dGam = np.zeros(Gm.shape + (n,))
    for l in (1, 2):
        e = np.zeros(n); e[l] = h
        dGam[..., l] = (gamma_batch(M, P + e, d)[0] - gamma_batch(M, P - e, d)[0]) / (2 * h)
    Ric = (np.einsum("madba->mbd", np.transpose(dGam, (0, 1, 2, 3, 4)).swapaxes(2, 4).swapaxes(3, 4)) if False else
           np.einsum("mabda->mbd", dGam) - np.einsum("maabd->mbd", dGam)
           + np.einsum("maae,medb->mbd", Gm, Gm) - np.einsum("made,meab->mbd", Gm, Gm))
    return Ric, G, Gi

# grid on [0,1]^2 (x, y); other coordinates trivial
NG = 33
xs = np.linspace(0, 1, NG); wq = np.ones(NG); wq[1:-1:2] = 4; wq[2:-1:2] = 2; wq *= (xs[1] - xs[0]) / 3
XX, YY = np.meshgrid(xs, xs, indexing="ij"); WW = np.outer(wq, wq)
def grid_points(n):
    P = np.zeros((NG * NG, n)); P[:, 1] = XX.ravel(); P[:, 2] = YY.ravel(); return P

# ===================================================================== Q. 4D: all of Einstein's equations
say("== Q. Gravity: a generic 4D metric, all 10 components varying; vary the loop action -> Einstein's equations")
def generic4(dsym=None, comp=None):
    s = [[0] * 4 for _ in range(4)]
    k = 0
    for i in range(4):
        for j in range(i, 4):
            k += 1
            val = 0.15 * sp.sin(1.3 * k * X + 0.7 * (j + 1) * Y + 0.4 * k) + 0.1 * sp.cos(0.9 * X * (i + 1) - 1.1 * k * Y)
            s[i][j] = s[j][i] = val
    g = sp.eye(4) + sp.Matrix(s) * 0.8
    if comp is not None:
        i, j = comp
        E = sp.zeros(4); E[i, j] = 1; E[j, i] = 1
        g = g + D * BUMP * E
    return g
P4 = grid_points(4)
M0 = compile_metric(generic4())
x0 = np.array([[0.0, 0.43, 0.61, 0.0]])
Rl0, Rt0, G0, Gi0 = scalar_from_loops(M0, x0, 0.0, eps=1e-3)
Ric0, _, Gi0f = ricci_formula(M0, x0, 0.0)
Rf0 = float(np.einsum("mbd,mbd->m", Gi0f, Ric0)[0])
check("Q1 generic metric (all 10 components vary): curvature scalar from loops = textbook formula",
      abs(Rl0[0] - Rf0) / abs(Rf0) < 2e-5, f"{Rl0[0]:.7f} vs {Rf0:.7f} (relative {abs(Rl0[0]-Rf0)/abs(Rf0):.1e})")
RicL = np.einsum("mabad->mbd", Rt0)
relQ2 = np.abs(RicL - Ric0).max() / np.abs(Ric0).max()
check("Q2 the whole Ricci tensor from loops = textbook formula (10 components)", relQ2 < 2e-5,
      f"max diff {np.abs(RicL - Ric0).max():.1e} of largest {np.abs(Ric0).max():.2f}")
# Einstein tensor (upper) from the formula, on the grid, for the expected variation
RicG, GG, GiG = ricci_formula(M0, P4, 0.0)
RG = np.einsum("mbd,mbd->m", GiG, RicG)
Gup = np.einsum("mab,mcd,mbd->mac", GiG, GiG, RicG) - 0.5 * RG[:, None, None] * GiG
sqrtg = np.sqrt(np.linalg.det(GG))
bump_grid = np.exp(-((XX.ravel() - 0.5)**2 + (YY.ravel() - 0.5)**2) / W**2)
def action4(M, d):
    R, _, G, _ = scalar_from_loops(M, P4, d)
    return float(np.sum(WW.ravel() * np.sqrt(np.linalg.det(G)) * R))
rowsQ = []
dd = 2e-3
for i in range(4):
    for j in range(i, 4):
        Mc = compile_metric(generic4(comp=(i, j)))
        dS = (action4(Mc, dd) - action4(Mc, -dd)) / (2 * dd)
        fac = 1 if i == j else 2
        expect = -float(np.sum(WW.ravel() * sqrtg * fac * Gup[:, i, j] * bump_grid))
        rowsQ.append(((i, j), dS, expect))
        say(f"   g_{i}{j}: loops {dS:+.6f}   Einstein {expect:+.6f}")
scaleQ = max(abs(e) for _, _, e in rowsQ)
worstQ = max(abs(a - b) for _, a, b in rowsQ) / scaleQ
check("Q3 nudging each of the 10 metric components: change in the loop action = -sqrt(g) G^{mu nu} (Einstein's equations)",
      worstQ < 2e-3, f"worst difference {worstQ:.1e} of the largest")
results["einstein"] = {"R_loops": float(Rl0[0]), "R_formula": Rf0, "ricci_maxdiff": float(np.abs(RicL - Ric0).max()),
                       "rows": [(f"{i}{j}", a, b) for (i, j), a, b in rowsQ], "worst_rel": worstQ}

# ===================================================================== S. 5D: Maxwell from the same loop action
say("\n== S. Light: the 5D chain of round 3 (flat 4D + potential A); vary the loop action in A -> Maxwell's equations")
kap = 1.0
def kk5(comp=None):
    A = [0, 0.3 * sp.sin(2.1 * X + 0.8 * Y) + 0.2 * Y**2, 0.25 * sp.cos(1.7 * Y - 0.6 * X) - 0.3 * X * Y, 0.2 * sp.sin(1.1 * X * Y + 0.5)]
    if comp is not None:
        A[comp] = A[comp] + D * BUMP
    g = sp.zeros(5)
    for i in range(4):
        for j in range(4):
            g[i, j] = (1 if i == j else 0) + kap**2 * A[i] * A[j]
        g[i, 4] = g[4, i] = kap * A[i]
    g[4, 4] = 1
    return g, A
P5 = grid_points(5)
_, A0 = kk5()
Afun = [sp.lambdify((X, Y), a, "numpy") for a in A0]
# Maxwell: d_mu F_{mu nu}, with F_{mu nu} = d_mu A_nu - d_nu A_mu (flat Euclidean 4D)
Fs = [[sp.diff(A0[n], [None, X, Y, None][m]) if m in (1, 2) else 0 for n in range(4)] for m in range(4)]
F = [[(Fs[m][n] - Fs[n][m]) for n in range(4)] for m in range(4)]
divF = [sum(sp.diff(F[m][n], [None, X, Y, None][m]) for m in (1, 2)) for n in range(4)]
F2 = sum(F[m][n]**2 for m in range(4) for n in range(4))
F2f = sp.lambdify((X, Y), F2, "numpy")
M5 = compile_metric(kk5()[0])
xs5 = np.array([[0.0, 0.43, 0.61, 0.0, 0.0]])
R5l, _, _, _ = scalar_from_loops(M5, xs5, 0.0, eps=1e-3)
R5pred = -kap**2 * F2f(0.43, 0.61) / 4
check("S1 5D curvature from loops = -(kappa^2/4) F_{mu nu} F_{mu nu} for a varying field (round 3 had a uniform one)",
      abs(R5l[0] - R5pred) < 1e-5 * max(1, abs(R5pred)), f"{R5l[0]:.7f} vs {R5pred:.7f}")
rowsS = []
for comp in (1, 2, 3):
    Mc = compile_metric(kk5(comp=comp)[0])
    def action5(d):
        R, _, G, _ = scalar_from_loops(Mc, P5, d)
        return float(np.sum(WW.ravel() * np.sqrt(np.linalg.det(G)) * R))
    dS = (action5(dd) - action5(-dd)) / (2 * dd)
    divf = sp.lambdify((X, Y), divF[comp], "numpy")
    expect = float(np.sum(WW.ravel() * kap**2 * np.broadcast_to(divf(XX.ravel(), YY.ravel()), XX.ravel().shape) * bump_grid))
    rowsS.append((comp, dS, expect))
    say(f"   A_{comp}: loops {dS:+.6f}   Maxwell {expect:+.6f}")
worstS = max(abs(a - b) for _, a, b in rowsS) / max(abs(b) for _, _, b in rowsS)
check("S2 nudging the potential A: change in the 5D loop action = kappa^2 d_mu F_{mu nu} (Maxwell's equations)",
      worstS < 2e-3, f"worst difference {worstS:.1e} of the largest")
results["maxwell"] = {"R5_loops": float(R5l[0]), "R5_pred": float(R5pred), "rows": rowsS, "worst_rel": worstS}

say(f"\n{passed} passed, {failed} failed")
results["summary"] = {"passed": passed, "failed": failed}
json.dump(results, open("results.json", "w"), indent=1)
open("checks_output.txt", "w").write(_log.getvalue())
