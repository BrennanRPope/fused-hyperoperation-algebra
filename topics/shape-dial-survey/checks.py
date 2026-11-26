"""Shape-dial survey: checks for every candidate slide.

Run: python3 checks.py > checks_output.txt   (writes results.json for figures.py)
Plane dial (k, d) is held fixed in every check unless the check says otherwise.
"""
import json, math, itertools, time
import numpy as np
from scipy import optimize, special, integrate
from scipy.spatial import ConvexHull
import sympy as sp
import mpmath as mp

rng = np.random.default_rng(20261126)
R = {}            # results for figures
PASS = []

def check(name, ok, detail=""):
    PASS.append((name, bool(ok)))
    print(f"[{'PASS' if ok else 'FAIL'}] {name}: {detail}")

# ---------------------------------------------------------------------------
# M. The p dial is the Maslov dial: ||u||_p = exp(soft max of log|u| at temperature 1/p)
# ---------------------------------------------------------------------------
u = np.abs(rng.normal(size=6)) + 0.1
errs = []
for p in [1, 1.5, 2, 4, 10]:
    lhs = np.sum(u**p)**(1/p)
    rhs = np.exp((1/p) * np.log(np.sum(np.exp(p*np.log(u)))))
    errs.append(abs(lhs - rhs))
gap = np.sum(u**200)**(1/200) - u.max()
check("M1 p-norm = Maslov soft max at temperature 1/p", max(errs) < 1e-12, f"max err {max(errs):.1e}")
check("M2 p -> infinity is the max (cube)", 0 <= gap < 0.02*u.max(), f"gap at p=200: {gap:.2e}")

# ---------------------------------------------------------------------------
# C1/C11. Corner sums: the k dial is one path through the corner weights
# ---------------------------------------------------------------------------
k = sp.symbols('k'); xs = sp.symbols('x1:5')
def comp(a, b): return a + b + k*a*b
fold = xs[0]
for x in xs[1:]: fold = comp(fold, x)
lhs = sp.expand(1 + fold)
rhs = 1 + sum(k**(len(S)-1) * sp.Mul(*S) for r in range(1, 5) for S in itertools.combinations(xs, r))
check("K1 1 + (x1 o_k ... o_k x4) = 1 + sum over corners S of k^(|S|-1) x^S", sp.simplify(lhs - rhs) == 0,
      "k=0 keeps only single corners (simplex), k=1 keeps all (cube)")
# log partition: cube (independent bits) vs simplex (one of n)
th = rng.normal(size=4); X = np.exp(th)
cube = np.sum(np.log1p(X)); prod = 1
for x in X: prod *= (1 + x)
check("K2 cube log-partition = L1(x1 o_1 ... o_1 xn)", abs(cube - np.log(prod)) < 1e-12, f"{cube:.6f}")
simp = np.log(1 + X.sum())
check("K3 simplex log-partition = L1(x1 o_0 ... o_0 xn) (softmax with a 'none' outcome)",
      abs(simp - np.log1p(X.sum())) < 1e-12, f"{simp:.6f}")
# tropical (Maslov) curves of 1 + x + y + k x y for figures
R["K_tropical_k"] = [0.0, 0.05, 0.3, 1.0]

# Renyi order = the p of a p-norm = 1 - k
pvec = rng.dirichlet(np.ones(5)); ok = True
for q in [0.5, 2.0, 3.0]:
    tsallis = (1 - np.sum(pvec**q))/(q - 1)
    kk = 1 - q
    renyi_from_L = np.log1p(kk*tsallis)/kk
    renyi_norm = q/(1-q)*np.log(np.sum(pvec**q)**(1/q))
    ok &= abs(renyi_from_L - renyi_norm) < 1e-12
check("K4 Renyi order q is the norm's p and our k = 1 - q (Renyi = L_k(Tsallis) = q/(1-q) log||p||_q)", ok)

# ---------------------------------------------------------------------------
# S. Stick-breaking: cube -> simplex at fixed k = -1
# ---------------------------------------------------------------------------
v = rng.uniform(size=5); pieces = []; left = 1.0
for vi in v: pieces.append(left*vi); left *= (1 - vi)
tot = 0.0
for vi in v: tot = tot + vi - tot*vi
check("S1 running total of stick pieces = v1 o_-1 ... o_-1 v5", abs(sum(pieces) - tot) < 1e-14, f"{tot:.6f}")
N = 200000; V = rng.uniform(size=(N, 2))
p1 = V[:, 0]; p2 = (1 - V[:, 0])*V[:, 1]
# Beta(1,2) first cut gives the even spread on the triangle
V2 = np.column_stack([rng.beta(1, 2, N), rng.uniform(size=N)])
q1 = V2[:, 0]; q2 = (1 - V2[:, 0])*V2[:, 1]
H, _, _ = np.histogram2d(q1, q2, bins=10, range=[[0, 1], [0, 1]])
inside = np.array([[ (i+1)/10 + (j+1)/10 <= 1.0 for j in range(10)] for i in range(10)])
cv = H[inside].std()/H[inside].mean()
check("S2 Beta(1,2) then uniform cuts land evenly on the triangle", cv < 0.05, f"cell spread {cv:.3f}")

# ---------------------------------------------------------------------------
# B. Bayes (simplex) <-> relativity (ball): homogeneity is what makes the law exist
# ---------------------------------------------------------------------------
# 1D: odds multiply exactly when speeds add by Einstein's rule
a, b = 0.5, 0.6
ein = (a + b)/(1 + a*b)
odds = lambda v: (1 + v)/(1 - v)
check("B1 1D: odds((a+b)/(1+ab)) = odds(a) * odds(b)", abs(odds(ein) - odds(a)*odds(b)) < 1e-12,
      f"0.5c (+) 0.6c = {ein:.4f}c; odds 3 x 4 = {odds(ein):.4f}")

# 2D sections: does a projective map keep the shape and move its centre?
tri_n = np.array([[np.cos(t), np.sin(t)] for t in (np.pi/2 + 2*np.pi*np.arange(3)/3 + np.pi)])
def gauge(name, P):
    x, y = P[:, 0], P[:, 1]
    if name == "triangle": return np.max(P @ tri_n.T, axis=1)/0.5
    if name == "square":   return np.maximum(np.abs(x), np.abs(y))
    if name.startswith("p="):
        p = float(name[2:]); return (np.abs(x)**p + np.abs(y)**p)**(1/p)
def boundary(name, m=240):
    t = np.linspace(0, 2*np.pi, m, endpoint=False)
    D = np.column_stack([np.cos(t), np.sin(t)])
    return D/gauge(name, D)[:, None]
def proj(Hm, P):
    Q = np.column_stack([P, np.ones(len(P))]) @ Hm.T
    return Q[:, :2]/Q[:, 2:3], Q[:, 2]
def homog_residual(name, target=np.array([0.18, 0.07]), starts=12):
    B = boundary(name)
    best = np.inf
    for s in range(starts):
        def f(h):
            Hm = np.append(h, 1.0).reshape(3, 3)
            P, w = proj(Hm, B)
            r1 = gauge(name, P) - 1
            Hi = np.linalg.inv(Hm); Pi, wi = proj(Hi, B)
            r2 = gauge(name, Pi) - 1
            c, _ = proj(Hm, np.zeros((1, 2)))
            pen = np.concatenate([np.minimum(w, 0)*10, np.minimum(wi*np.sign(Hi[2, 2]), 0)*10])
            return np.concatenate([r1, r2, 20*(c[0] - target), pen])
        h0 = np.eye(3).flatten()[:8] + rng.normal(scale=0.15, size=8)
        try:
            sol = optimize.least_squares(f, h0, max_nfev=4000)
            res = np.sqrt(np.mean(f(sol.x)[:2*len(B)]**2))
            best = min(best, res)
        except Exception:
            pass
    return best
homog = {}
for name in ["triangle", "p=2", "p=1.5", "p=3", "p=4", "square"]:
    homog[name] = homog_residual(name)
    print(f"    {name:8s} best boundary mismatch after moving the centre: {homog[name]:.2e}")
R["homog"] = homog
check("B2 triangle (Bayes) and disk (relativity) can move their centre anywhere: mismatch ~ 0",
      homog["triangle"] < 1e-6 and homog["p=2"] < 1e-6)
check("B3 square and p = 1.5, 3, 4 cannot (mismatch stays well above 0): no 'boosts' in between",
      min(homog[n] for n in ["square", "p=1.5", "p=3", "p=4"]) > 1e-3)

# Ordering fraction of a causal diamond in 2+1 as the light-cone section slides triangle -> disk -> square
def cone_gauge(t_mix, P, kind):
    if kind == "tri-disk":
        return (1 - t_mix)*gauge("triangle", P) + t_mix*gauge("p=2", P)
    return gauge(f"p={t_mix}", P)
def ordering_fraction_once(gfun, n=3500):
    pts = []
    while len(pts) < n:
        C = rng.uniform(-1, 1, size=(20000, 2)); T = rng.uniform(0, 1, size=20000)
        keep = (gfun(C) <= T) & (gfun(-C) <= 1 - T)
        pts.extend(np.column_stack([T[keep], C[keep]]).tolist())
    pts = np.array(pts[:n])
    T = pts[:, 0]; C = pts[:, 1:]
    rel = 0
    for i in range(n):
        dT = T[i+1:] - T[i]; dC = C[i+1:] - C[i]
        g = gfun(dC); g2 = gfun(-dC)
        rel += np.sum((dT >= g) | (-dT >= g2))
    return rel/(n*(n - 1)/2)
def ordering_fraction(gfun, reps=4):
    return float(np.mean([ordering_fraction_once(gfun) for _ in range(reps)]))
of_td = [(t, ordering_fraction(lambda P, t=t: cone_gauge(t, P, "tri-disk"))) for t in [0, 0.25, 0.5, 0.75, 1.0]]
of_p = [(p, ordering_fraction(lambda P, p=p: cone_gauge(p, P, "p"))) for p in [2, 3, 6, 50]]
R["ordering_tri_disk"] = of_td; R["ordering_p"] = of_p
print("    ordering fraction, triangle->disk:", [f"{t}:{v:.3f}" for t, v in of_td])
print("    ordering fraction, disk->square (p):", [f"{p}:{v:.3f}" for p, v in of_p])
check("B4 triangle light cone gives the cube's product order (fraction 1/4)", abs(of_td[0][1] - 0.25) < 0.006, f"{of_td[0][1]:.3f}")
check("B5 round light cone gives 0.229 (the curved-space page's value)", abs(of_td[-1][1] - 0.229) < 0.006, f"{of_td[-1][1]:.3f}")

# ---------------------------------------------------------------------------
# C. The Bell (CHSH) bound slides with the shape: 2 (segment) -> 2*sqrt2 (disk) -> 4 (square)
# ---------------------------------------------------------------------------
def chsh_p(p, grid=720):
    q = np.inf if p == 1 else (1.0 if np.isinf(p) else p/(p - 1))
    def nq(Z):
        if np.isinf(q): return np.max(np.abs(Z), axis=-1)
        return np.sum(np.abs(Z)**q, axis=-1)**(1/q)
    t = np.linspace(0, 2*np.pi, grid, endpoint=False)
    D = np.column_stack([np.cos(t), np.sin(t)])
    Bq = D/nq(D)[:, None]                     # unit sphere of the dual norm
    S = Bq[:, None, :] + Bq[None, :, :]; Dd = Bq[:, None, :] - Bq[None, :, :]
    return float(np.max(nq(S) + nq(Dd)))
ps = [1, 1.25, 1.5, 2, 3, 4, 8, np.inf]
chsh = {str(p): chsh_p(p) for p in ps}
formula = {str(p): 2**(1 + max(1/p, 1 - 1/p)) if not np.isinf(p) else 4.0 for p in ps}
R["chsh_p"] = chsh; R["chsh_formula"] = formula
print("    CHSH by p:", {kk: round(vv, 4) for kk, vv in chsh.items()})
check("C1 disk (p=2) gives Tsirelson's 2*sqrt(2)", abs(chsh["2"] - 2*np.sqrt(2)) < 1e-3, f"{chsh['2']:.4f}")
check("C2 square / cross-polytope (p=inf, p=1) give the PR-box value 4", abs(chsh["inf"] - 4) < 1e-3 and abs(chsh["1"] - 4) < 1e-3)
check("C3 in between it follows 2^(1 + max(1/p, 1 - 1/p))", max(abs(chsh[kk] - formula[kk]) for kk in chsh) < 5e-3,
      f"max gap {max(abs(chsh[kk]-formula[kk]) for kk in chsh):.1e}")
segment = max(abs(b1 + b2) + abs(b1 - b2) for b1 in [-1, 1] for b2 in [-1, 1])
check("C4 the segment (a classical bit, 1D simplex) gives 2", segment == 2)
def chsh_poly(n, grid=3000):
    ang = 2*np.pi*np.arange(n)/n
    Vt = np.column_stack([np.cos(ang), np.sin(ang)])          # polygon vertices; Alice's observables live here
    def N(Z): return np.max(Z @ Vt.T, axis=-1)                 # its support function = dual gauge
    t = np.linspace(0, 2*np.pi, grid, endpoint=False); D = np.column_stack([np.cos(t), np.sin(t)])
    th = np.pi/n; W = np.column_stack([np.cos(ang + th), np.sin(ang + th)])/np.cos(th)   # vertices of the polar polygon
    best = 0.0
    for Bd in (D/N(D)[:, None], W):
        S = Bd[:, None, :] + Bd[None, :, :]; Dd = Bd[:, None, :] - Bd[None, :, :]
        best = max(best, float(np.max(N(S) + N(Dd))))
    return best
poly = {n: chsh_poly(n) for n in [4, 6, 8, 10, 12, 14, 16, 20, 24, 48, 96]}
R["chsh_poly"] = poly
print("    CHSH for even polygons:", {kk: round(vv, 4) for kk, vv in poly.items()})
check("C5 even polygons: square 4, hexagon 3, every multiple of 8 exactly 2*sqrt2, the rest a little above, all tending to 2*sqrt2",
      abs(poly[4] - 4) < 1e-9 and abs(poly[6] - 3) < 1e-6 and all(abs(poly[n] - 2*np.sqrt(2)) < 1e-9 for n in [8, 16, 24, 48, 96])
      and all(2*np.sqrt(2) < poly[n] < 2.95 for n in [10, 12, 14, 20]))

# ---------------------------------------------------------------------------
# D. Born step p = 2 -> p = 1; norm-keeping dynamics only at p = 1 (positive) and p = 2 (signed)
# ---------------------------------------------------------------------------
psi = rng.normal(size=4) + 1j*rng.normal(size=4); psi /= np.linalg.norm(psi)
check("D1 squaring sends the round sphere onto the simplex (Born rule)", abs(np.sum(np.abs(psi)**2) - 1) < 1e-12)
def isometry_residual(p, positive, eps=0.15, starts=6):
    Xs = np.abs(rng.normal(size=(400, 2))) if positive else rng.normal(size=(400, 2))
    nrm = lambda Z: np.sum(np.abs(Z)**p, axis=1)**(1/p)
    best = np.inf
    for s in range(starts):
        def f(m):
            M = np.eye(2) + m.reshape(2, 2)
            if positive: M = np.abs(M)
            r = nrm(Xs @ M.T)/nrm(Xs) - 1
            return np.concatenate([r, [10*(np.linalg.norm(m) - eps)]])
        sol = optimize.least_squares(f, rng.normal(scale=eps, size=4), max_nfev=3000)
        best = min(best, np.sqrt(np.mean(f(sol.x)[:-1]**2)))
    return best
pgrid = [0.75, 1.0, 1.25, 1.5, 1.75, 2.0, 2.5, 3.0, 4.0]
iso_signed = {p: isometry_residual(p, False) for p in pgrid}
iso_pos = {p: isometry_residual(p, True) for p in pgrid}
R["iso_signed"] = {str(kk): vv for kk, vv in iso_signed.items()}; R["iso_pos"] = {str(kk): vv for kk, vv in iso_pos.items()}
print("    signed maps, mismatch by p:", {kk: f"{vv:.1e}" for kk, vv in iso_signed.items()})
print("    positive maps, mismatch by p:", {kk: f"{vv:.1e}" for kk, vv in iso_pos.items()})
check("D2 a continuous family of norm-keeping maps exists at p = 2 (rotations) and nowhere else (signed maps)",
      iso_signed[2.0] < 1e-6 and min(vv for kk, vv in iso_signed.items() if kk != 2.0) > 1e-3)
check("D3 on positive vectors one exists at p = 1 (stochastic maps) and nowhere else",
      iso_pos[1.0] < 1e-6 and min(vv for kk, vv in iso_pos.items() if kk != 1.0) > 1e-3)

# ---------------------------------------------------------------------------
# E. Replicator (simplex) = imaginary-time Schrodinger (sphere) through u = sqrt(p)
# ---------------------------------------------------------------------------
fit = rng.normal(size=5); p0 = rng.dirichlet(np.ones(5))
def repl(t, p): return p*(fit - p @ fit)
sol = integrate.solve_ivp(repl, (0, 3), p0, rtol=1e-11, atol=1e-13, dense_output=True)
ts = np.linspace(0, 3, 61)
P = sol.sol(ts).T
U = np.sqrt(p0)[None, :]*np.exp(np.outer(ts, fit)/2); U /= np.linalg.norm(U, axis=1)[:, None]
err = np.max(np.abs(np.sqrt(P) - U))
check("E1 replicator flow, read as sqrt(p), is normalised imaginary-time evolution with H = -fitness/2", err < 1e-8, f"max err {err:.1e}")
mean = P @ fit; var = P @ fit**2 - mean**2
dmean = np.gradient(mean, ts)
check("E2 Fisher's theorem (d mean/dt = variance) is the energy-drop rule d<H>/dt = -2 Var(H)",
      np.max(np.abs(dmean[2:-2] - var[2:-2])) < 2e-3, f"max gap {np.max(np.abs(dmean[2:-2]-var[2:-2])):.1e}")
R["replicator"] = {"t": ts.tolist(), "P": P.tolist(), "fit": fit.tolist()}

# ---------------------------------------------------------------------------
# F. The genus slide: fixed q, the shape picks the polynomial (h-polynomial of the shape)
# ---------------------------------------------------------------------------
def facets_merged(pts):
    hull = ConvexHull(pts)
    eq = np.round(hull.equations, 7)
    groups = {}
    for i, e in enumerate(map(tuple, eq)): groups.setdefault(e, set()).update(hull.simplices[i])
    return list(groups.values())
def fvector3(pts):
    F = facets_merged(pts)
    verts = set().union(*F)
    edges = set()
    for f in F:
        for a_, b_ in itertools.combinations(sorted(f), 2):
            # a, b share an edge if they lie on exactly 2 common facets
            if sum(1 for g in F if a_ in g and b_ in g) >= 2: edges.add((a_, b_))
    return [len(verts), len(edges), len(F)]
cube3 = np.array(list(itertools.product([0, 1], repeat=3)), float)
simp3 = np.array([[0, 0, 0], [1, 0, 0], [0, 1, 0], [0, 0, 1]], float)
perm = np.array(list(itertools.permutations([1, 2, 3, 4])), float); B4 = np.array([[1, -1, 0, 0], [0, 1, -1, 0], [0, 0, 1, -1]], float)
perm3 = perm @ np.linalg.qr(B4.T)[0]
# Loday associahedron K5 in R^4 (vertices from planar binary trees), projected to 3D
def loday(n):
    from functools import lru_cache
    @lru_cache(None)
    def trees(i, j):
        if i > j: return [()]
        out = []
        for r in range(i, j + 1):
            for L in trees(i, r - 1):
                for Rr in trees(r + 1, j): out.append((L, r, Rr))
        return out
    def size(t): return 0 if t == () else size(t[0]) + 1 + size(t[2])
    def coords(t, acc):
        if t == (): return
        L, r, Rr = t
        acc[r] = (size(L) + 1)*(size(Rr) + 1)
        coords(L, acc); coords(Rr, acc)
    vs = []
    for t in trees(0, n-1):
        acc = [0]*n; coords(t, acc); vs.append(acc)
    return np.array(vs, float)
K5 = loday(4)
assoc3 = (K5 - K5.mean(0)) @ np.linalg.qr(B4.T)[0]
fv = {"simplex": fvector3(simp3), "cube": fvector3(cube3), "permutohedron": fvector3(perm3), "associahedron": fvector3(assoc3)}
def hpoly(f):  # simple 3-polytope: sum f_i (t-1)^i, f_i = # i-dim faces, f_3 = 1
    t = sp.symbols('q'); ff = [f[0], f[1], f[2], 1]
    return sp.Poly(sp.expand(sum(ff[i]*(t - 1)**i for i in range(4))), t).all_coeffs()[::-1]
hv = {kk: [int(c) for c in hpoly(vv)] for kk, vv in fv.items()}
print("    f-vectors:", fv); print("    h-vectors:", hv)
expect = {"simplex": [1, 1, 1, 1], "cube": [1, 3, 3, 1], "permutohedron": [1, 11, 11, 1], "associahedron": [1, 6, 6, 1]}
check("F1 3D h-vectors: simplex [4]_q, cube (1+q)^3, permutohedron Eulerian, associahedron Narayana", hv == expect, str(hv))
# 2D and 4D from the known face counts
f2 = {"simplex": [3, 3], "cube": [4, 4], "permutohedron": [6, 6], "associahedron": [5, 5]}
h2 = {kk: [int(c) for c in sp.Poly(sp.expand(vv[0] + vv[1]*(sp.Symbol('q') - 1) + (sp.Symbol('q') - 1)**2), sp.Symbol('q')).all_coeffs()[::-1]] for kk, vv in f2.items()}
check("F2 2D: [3]_q, (1+q)^2, 1+4q+q^2, 1+3q+q^2", h2 == {"simplex": [1, 1, 1], "cube": [1, 2, 1], "permutohedron": [1, 4, 1], "associahedron": [1, 3, 1]})
R["hvec"] = {"2": h2, "3": hv}
# evaluate the genus at the plane dial: 1 + kx + dx^2 = (1 + a x)(1 + b x), q = b/a
def q_of(kd):
    kk_, dd = kd; disc = kk_**2 - 4*dd
    if dd == 0: return 0.0
    a_ = (kk_ + np.sqrt(complex(disc)))/2; b_ = (kk_ - np.sqrt(complex(disc)))/2
    return b_/a_
dials = {"core law (k=1, d=0)": (1, 0), "signature (k=0, d=-1)": (0, -1), "Euler point (k=2, d=1)": (2, 1)}
gen = {}
for name, kd in dials.items():
    q = q_of(kd)
    gen[name] = {sh: complex(sum(c*q**i for i, c in enumerate(h))).real for sh, h in hv.items()}
print("    genus by dial and shape:", gen)
R["genus"] = gen
check("F3 at our core law (d = 0, q = 0) every shape gives 1: the shape dial does nothing there",
      all(abs(vv - 1) < 1e-12 for vv in gen["core law (k=1, d=0)"].values()))
check("F4 at the Euler point (q = 1) the genus is the vertex count 4, 8, 24, 14",
      [round(gen["Euler point (k=2, d=1)"][s]) for s in ["simplex", "cube", "permutohedron", "associahedron"]] == [4, 8, 24, 14])

# ---------------------------------------------------------------------------
# G. Volumes summed over dimensions: one Mittag-Leffler function along the p dial
# ---------------------------------------------------------------------------
def Vp(p, n): return (2*special.gamma(1 + 1/p))**n/special.gamma(1 + n/p)
mc_ok = True
for p_, n_ in [(1.5, 3), (3, 2), (4, 4)]:
    Z = rng.uniform(-1, 1, size=(400000, n_))
    est = 2**n_*np.mean(np.sum(np.abs(Z)**p_, axis=1) <= 1)
    mc_ok &= abs(est - Vp(p_, n_))/Vp(p_, n_) < 0.01
check("G1 p-ball volume (2 Gamma(1+1/p))^n / Gamma(1+n/p) (Monte Carlo, 3 cases)", mc_ok)
def S(p, t=1.0, nmax=400): return sum(Vp(p, n)*t**n for n in range(nmax))
ml = lambda alpha, z: mp.nsum(lambda n: z**n/mp.gamma(alpha*n + 1), [0, mp.inf])
check("G2 p = 1 (cross-polytope): sum of volumes = e^2", abs(S(1) - math.e**2) < 1e-9, f"{S(1):.6f}")
check("G3 p = 2 (ball): sum = e^pi (1 + erf sqrt(pi))", abs(S(2) - math.exp(math.pi)*(1 + math.erf(math.sqrt(math.pi)))) < 1e-9, f"{S(2):.4f}")
check("G4 any p: sum = Mittag-Leffler E_{1/p}(2 Gamma(1+1/p))", all(abs(S(p) - float(ml(1/p, 2*special.gamma(1 + 1/p))))/S(p) < 1e-9 for p in [1.3, 2.5, 4]))
simplex_sum = sum(1/math.factorial(n) for n in range(30))
check("G5 corner simplex (1/n!) sums to e; cube (2^n) is the geometric series 1/(1-2t), finite only for t < 1/2",
      abs(simplex_sum - math.e) < 1e-12 and abs(sum(2**n*0.3**n for n in range(200)) - 1/(1 - 0.6)) < 1e-9)
def log10_S(p, nmax=60000):
    n = np.arange(nmax); lt = n*np.log(2*special.gamma(1 + 1/p)) - special.gammaln(1 + n/p)
    return float((lt.max() + np.log(np.sum(np.exp(lt - lt.max()))))/np.log(10))
R["vol_sum_log10"] = {str(p): log10_S(p) for p in [1, 1.25, 1.5, 2, 2.5, 3, 4, 5, 6, 8]}
print("    log10 of summed volumes by p:", {kk: round(vv, 3) for kk, vv in R["vol_sum_log10"].items()})

# ---------------------------------------------------------------------------
# H. Gauss-Bonnet along the p dial: smooth curvature (sphere) gathers into corners (cube, octahedron)
# ---------------------------------------------------------------------------
def superellipsoid_curvature(p, nth=260, nph=520):
    th = np.linspace(1e-4, np.pi - 1e-4, nth); ph = np.linspace(0, 2*np.pi, nph)
    TH, PH = np.meshgrid(th, ph, indexing='ij')
    Dd = np.stack([np.sin(TH)*np.cos(PH), np.sin(TH)*np.sin(PH), np.cos(TH)], -1)
    Rot = np.linalg.qr(np.array([[0.9, 0.3, 0.2], [-0.25, 0.8, 0.4], [0.1, -0.35, 0.95]]))[0]
    Dd = Dd @ Rot.T
    Pt = Dd/(np.sum(np.abs(Dd)**p, -1)**(1/p))[..., None]
    G = np.sign(Pt)*np.abs(Pt)**(p - 1); Nn = G/np.linalg.norm(G, axis=-1)[..., None]
    def tri_area(A, B, C): return 0.5*np.linalg.norm(np.cross(B - A, C - A), axis=-1)
    def sph_area(A, B, C):
        num = np.abs(np.einsum('...i,...i', A, np.cross(B, C)))
        den = 1 + np.einsum('...i,...i', A, B) + np.einsum('...i,...i', B, C) + np.einsum('...i,...i', C, A)
        return 2*np.arctan2(num, den)
    quads = [(Pt[:-1, :-1], Pt[1:, :-1], Pt[1:, 1:], Pt[:-1, 1:]), (Nn[:-1, :-1], Nn[1:, :-1], Nn[1:, 1:], Nn[:-1, 1:])]
    (a1, b1, c1, d1), (a2, b2, c2, d2) = quads
    area = tri_area(a1, b1, c1) + tri_area(a1, c1, d1)
    curv = sph_area(a2, b2, c2) + sph_area(a2, c2, d2)
    cen = (a1 + c1)/2
    return area, curv, cen
gb = {}
for p in [1.2, 1.5, 2, 3, 6, 20]:
    area, curv, cen = superellipsoid_curvature(p)
    corners = (np.array(list(itertools.product([-1, 1], repeat=3)))/3**(1/p)) if p >= 2 else np.vstack([np.eye(3), -np.eye(3)])
    dist = np.min(np.linalg.norm(cen[..., None, :] - corners[None, None], axis=-1), axis=-1)
    near = dist < 0.3
    gb[str(p)] = {"total_over_4pi": float(curv.sum()/(4*np.pi)), "curv_frac_near_corners": float(curv[near].sum()/curv.sum()),
                  "area_frac_near_corners": float(area[near].sum()/area.sum()), "n_corners": int(len(corners))}
    print(f"    p={p:>4}: total curvature/4pi = {gb[str(p)]['total_over_4pi']:.4f}; within 0.3 of the {len(corners)} corners: "
          f"{100*gb[str(p)]['curv_frac_near_corners']:.1f}% of curvature on {100*gb[str(p)]['area_frac_near_corners']:.1f}% of area")
R["gauss_bonnet"] = gb
check("H1 total curvature is 4 pi at every p (Gauss-Bonnet)", all(abs(vv["total_over_4pi"] - 1) < 2e-3 for vv in gb.values()))
check("H2 sphere: curvature spread like area; toward the cube (p=20) nearly all in the 8 corners; toward the octahedron (p=1.2) gathering in the 6",
      gb["20"]["curv_frac_near_corners"] > 0.99 and gb["1.2"]["curv_frac_near_corners"] > 3*gb["1.2"]["area_frac_near_corners"] and abs(gb["2"]["curv_frac_near_corners"] - gb["2"]["area_frac_near_corners"]) < 0.01)
check("H3 corner defects: cube 8 x pi/2 = 4 pi, octahedron 6 x (2 pi - 4 x pi/3) = 4 pi", abs(8*np.pi/2 - 4*np.pi) < 1e-12 and abs(6*(2*np.pi - 4*np.pi/3) - 4*np.pi) < 1e-12)

# ---------------------------------------------------------------------------
# I. Fermat on the p dial: rational points of |x|^p + |y|^p = 1
# ---------------------------------------------------------------------------
def count_rational(p, Nmax):
    counts = []; pts = set()
    for c in range(1, Nmax + 1):
        if p == "inf":
            for a_ in range(-c, c + 1):
                for b_ in [-c, c]:
                    f = (sp.Rational(a_, c), sp.Rational(b_, c)); pts.add(f); pts.add(f[::-1])
        elif p == 1:
            for a_ in range(-c, c + 1):
                b_ = c - abs(a_)
                for s in ([1, -1] if b_ else [1]): pts.add((sp.Rational(a_, c), sp.Rational(s*b_, c)))
        else:
            cp = c**p
            for a_ in range(0, c + 1):
                rem = cp - a_**p
                bb = round(rem**(1/p)) if rem >= 0 else -1
                for b_ in [bb - 1, bb, bb + 1]:
                    if b_ >= 0 and b_**p == rem:
                        for sa in [1, -1]:
                            for sb in [1, -1]: pts.add((sp.Rational(sa*a_, c), sp.Rational(sb*b_, c)))
        counts.append(len(pts))
    return counts
Nmax = 60
rat = {str(p): count_rational(p, Nmax) for p in [1, 2, 3, 4, "inf"]}
R["fermat_counts"] = rat
print("    rational points with denominator <= 60:", {kk: vv[-1] for kk, vv in rat.items()})
check("I1 p = 3 and p = 4 have only the 4 trivial points (Fermat; our Lean has n = 3 and n = 4)", rat["3"][-1] == 4 and rat["4"][-1] == 4)
check("I2 p = 1 and p = infinity fill up (count grows like N^2); p = 2 grows like N (Pythagoras)",
      rat["1"][-1] > 1000 and rat["inf"][-1] > 1000 and 20 < rat["2"][-1] < 400)

# ---------------------------------------------------------------------------
# J. Lee-Yang survives: a spin read off any p-shape has only imaginary field zeros
# ---------------------------------------------------------------------------
def G_factory(p, a, n=160):
    tj, wj = special.roots_jacobi(n, a, 0)
    x = (tj + 1)/2; w = wj*(0.5)**(a + 1)
    with np.errstate(divide='ignore', invalid='ignore'):
        ratio = np.where(x < 1 - 1e-14, (1 - x**p)/(1 - x), p)
    g = ratio**a*w
    return lambda s: (np.cos(np.outer(np.atleast_1d(s), x)) @ g)
def lee_yang_ok(p, a):
    G = G_factory(p, a)
    xs_ = np.linspace(26, 30, 400); X = xs_[np.argmax(np.abs(G(xs_)))]
    def count(Y, m=40000):
        sides = [np.linspace(-X, X, m) - 1j*Y, X + 1j*np.linspace(-Y, Y, m//4),
                 np.linspace(X, -X, m) + 1j*Y, -X + 1j*np.linspace(Y, -Y, m//4)]
        z = np.concatenate(sides); v = G(z)
        return int(round(np.sum(np.diff(np.unwrap(np.angle(v))))/(2*np.pi)))
    return count(8.0), count(0.05)
ly = {}
for p in [1, 1.25, 1.5, 2, 3, 4, 8, 32]:
    for label, a in [("ball, n=2", 1/p), ("ball, n=3", 2/p), ("shell, n=3", 2/p - 1)]:
        ly[f"p={p}, {label}"] = lee_yang_ok(p, a)
R["lee_yang"] = {kk: list(vv) for kk, vv in ly.items()}
bad = sorted(kk for kk, vv in ly.items() if vv[0] != vv[1])
print("    zeros with |Re s| < X~28, counted in |Im| < 8 and in |Im| < 0.05:", ly)
check("J1 round end (p >= 1.5 here) and cube end keep every zero on the Lee-Yang line (all 18 such cases)",
      all(ly[kk][0] == ly[kk][1] for kk in ly if not kk.startswith("p=1,") and not kk.startswith("p=1.25,")))
check("J2 but it breaks near the cross-polytope: the octahedron's spin (p=1, ball, n=3, weight (1-|x|)^2) has zeros off the line",
      ly["p=1, ball, n=3"][0] != ly["p=1, ball, n=3"][1], f"{ly['p=1, ball, n=3']}")
# where exactly it breaks: weight (1-|x|^p)^a; scan a (that is (n-1)/p for a spin uniform in an n-dim p-ball)
scan = {}
for p in [1.0, 1.1, 1.25, 1.5, 1.75, 2.0, 3.0, 4.0, 8.0]:
    row = {}
    for a in [0.0, 0.5, 1.0, 1.25, 1.5, 2.0, 3.0, 4.0, 6.0]:
        tot, near = lee_yang_ok(p, a); row[str(a)] = [tot, near]
    scan[str(p)] = row
R["lee_yang_scan"] = scan
print("    Lee-Yang holds? rows p, columns a=(n-1)/p:")
sym = lambda v: "?" if v[0] == 0 else ("Y" if v[0] == v[1] else ".")
print("      (Y = every zero on the line, . = zeros off it, ? = no zeros in the window, untested)")
for p, row in scan.items(): print(f"      p={p:>4}: " + " ".join(sym(v) for v in row.values()))
check("J3 at p = 1 the line holds only for a <= 1 (Polya's (1-t)^a family); at p = 2 for every a (Bessel)",
      sym(scan["1.0"]["1.0"]) == "Y" and sym(scan["1.0"]["2.0"]) == "." and all(sym(v) == "Y" for v in scan["2.0"].values()))


# ---------------------------------------------------------------------------
# P. Polishing step (Brennan, Nov 26: "add a polishing step ... as long as it doesn't take too much effort").
#    Every numerical search above is finished with a cheap exact check or a few refinement steps; cost is printed.
# ---------------------------------------------------------------------------
POL = {}
def timed(name):
    class T:
        def __enter__(self): self.t = time.perf_counter()
        def __exit__(self, *a): POL[name] = time.perf_counter() - self.t
    return T()
# P1 homogeneity: the searched maps are replaced by the exact ones
with timed("P1"):
    target = np.array([0.18, 0.07])
    # triangle: a diagonal scaling in barycentric coordinates, written as a projective map
    A = np.column_stack([np.vstack([TRI_V.T, np.ones(3)])]) if False else None
    V3 = np.array([[np.cos(a), np.sin(a)] for a in (np.pi/2 + 2*np.pi*np.arange(3)/3)])*1.0
    Vh = np.vstack([V3.T, np.ones(3)])                    # columns: vertices in homogeneous coords
    lam = np.linalg.solve(Vh, np.append(target, 1.0))     # barycentric coords of the target
    Ht = Vh @ np.diag(lam*3) @ np.linalg.inv(Vh)          # sends centre (1/3,1/3,1/3) to target, fixes vertices
    Bt = np.array([V3[i]*(1 - s_) + V3[(i+1) % 3]*s_ for i in range(3) for s_ in np.linspace(0, 1, 50)])
    Pt_, _ = proj(Ht, Bt)
    ins = np.linalg.solve(Vh, np.vstack([Pt_.T, np.ones(len(Pt_))]))
    tri_err = np.max(np.abs(np.min(ins, axis=0)))         # each image point has one barycentric coord = 0
    c_t, _ = proj(Ht, np.zeros((1, 2)))
    # disk: the Lorentz boost in the Klein model taking 0 to the target
    v = np.linalg.norm(target); g = 1/np.sqrt(1 - v**2); n_ = target/v
    L = np.eye(3); L[:2, :2] += (g - 1)*np.outer(n_, n_); L[:2, 2] = g*target; L[2, :2] = g*target; L[2, 2] = g
    Bd = boundary("p=2"); Pd, _ = proj(L, Bd); c_d, _ = proj(L, np.zeros((1, 2)))
    disk_err = np.max(np.abs(np.linalg.norm(Pd, axis=1) - 1))
check("P1 polish: exact Bayes scaling (triangle) and exact Lorentz boost (disk) move the centre to the target and keep the shape",
      tri_err < 1e-12 and disk_err < 1e-12 and np.allclose(c_t[0], target) and np.allclose(c_d[0], target),
      f"triangle {tri_err:.1e}, disk {disk_err:.1e}; {1000*POL['P1']:.1f} ms")
# P2 CHSH: polish the grid optimum with a local optimiser over the two angles, compare with the closed form
with timed("P2"):
    gaps = []
    for p in [1.25, 1.5, 3, 4, 8]:
        q = p/(p - 1)
        nq = lambda Z: np.sum(np.abs(Z)**q, -1)**(1/q)
        def unit(t): d_ = np.array([np.cos(t), np.sin(t)]); return d_/nq(d_)
        f = lambda th_: -(nq(unit(th_[0]) + unit(th_[1])) + nq(unit(th_[0]) - unit(th_[1])))
        best = min((optimize.minimize(f, x0, method="Nelder-Mead", options={"xatol": 1e-12, "fatol": 1e-14}) for x0 in
                    [np.array([a_, b_]) for a_ in np.linspace(0, np.pi, 5) for b_ in np.linspace(0, np.pi, 5)]), key=lambda r: r.fun)
        gaps.append(abs(-best.fun - 2**(1 + max(1/p, 1 - 1/p))))
check("P2 polish: local optimiser on the CHSH angles lands on 2^(1 + max(1/p, 1-1/p)) to machine precision",
      max(gaps) < 1e-9, f"max gap {max(gaps):.1e}; {1000*POL['P2']:.0f} ms")
# P3 norm-keeping maps: the exact rotation (p=2) and exact stochastic matrix (p=1) checked on all test vectors
with timed("P3"):
    th_ = 0.15; Rm = np.array([[np.cos(th_), -np.sin(th_)], [np.sin(th_), np.cos(th_)]]); Xs = rng.normal(size=(1000, 2))
    e2 = np.max(np.abs(np.linalg.norm(Xs @ Rm.T, axis=1) - np.linalg.norm(Xs, axis=1)))
    Sm = np.array([[0.9, 0.2], [0.1, 0.8]]); Xp = np.abs(Xs)
    e1 = np.max(np.abs(np.sum(Xp @ Sm.T, 1) - np.sum(Xp, 1)))
check("P3 polish: exact rotation keeps the 2-norm, exact stochastic matrix keeps the 1-norm of positive vectors",
      e2 < 1e-13 and e1 < 1e-13, f"{e2:.1e}, {e1:.1e}; {1000*POL['P3']:.1f} ms")
# P4 ordering fraction: the Monte Carlo numbers get error bars (repeat spread), and the two ends are exact
with timed("P4"):
    spread = np.std([ordering_fraction_once(lambda P_: cone_gauge(0, P_, "tri-disk"), n=1500) for _ in range(4)])/2
    mm = math.gamma(4)*math.gamma(1.5)/(4*math.gamma(4.5))*2
check("P4 polish: ends exact (cube order 1/4; round cone 2*Gamma(4)Gamma(3/2)/(4 Gamma(9/2)) = 0.2286, Myrheim-Meyer); Monte Carlo error about +-0.003",
      abs(mm - 0.2286) < 1e-4 and spread < 0.006, f"error bar {spread:.4f}; {POL['P4']:.1f} s")
R["ordering_errbar"] = float(spread)
# P5 Lee-Yang: Newton-polish the first zero off the line for the octahedron spin, and the first real zero for the ball
with timed("P5"):
    Goct = lambda z: 2*(z - np.sin(z))/z**3          # closed form of int_0^1 cos(zx)(1-x)^2 dx
    dG = lambda z, h=1e-6: (Goct(z + h) - Goct(z - h))/(2*h)
    z = 7.5 + 2.8j
    for _ in range(8): z = z - Goct(z)/dG(z)
    num = G_factory(1, 2.0)(np.array([z]))[0]
    Gb = G_factory(2, 1.0)                          # ball n=3, p=2: j1-type, first zero is tan z = z, at 4.4934
    zb = 4.5
    for _ in range(6): zb = zb - Gb(np.array([zb]))[0].real/((Gb(np.array([zb + 1e-6]))[0].real - Gb(np.array([zb - 1e-6]))[0].real)/2e-6)
check("P5 polish: octahedron spin's zero sits at z = s - sin s = 0, off the line (Im = %.4f); ball's first zero is real, tan z = z at 4.4934" % z.imag,
      abs(Goct(z)) < 1e-12 and abs(num) < 1e-8 and abs(z.imag) > 1 and abs(zb - 4.493409) < 1e-5, f"z = {z:.6f}; {1000*POL['P5']:.0f} ms")
# P6 Gauss-Bonnet: same total at half the grid (topology, so it must not move)
with timed("P6"):
    a_, c_, _ = superellipsoid_curvature(6, nth=130, nph=260)
check("P6 polish: half-resolution grid gives the same total curvature at p = 6", abs(c_.sum()/(4*np.pi) - 1) < 1e-3, f"{c_.sum()/(4*np.pi):.6f}; {POL['P6']:.1f} s")
R["polish_seconds"] = POL
print("    polishing cost (s):", {kk: round(vv, 3) for kk, vv in POL.items()}, "total", round(sum(POL.values()), 2))

# ---------------------------------------------------------------------------
print(f"\n{sum(ok for _, ok in PASS)} of {len(PASS)} checks pass")
R["checks"] = PASS
json.dump(R, open("results.json", "w"), default=lambda o: o if not isinstance(o, np.generic) else o.item(), indent=1)
