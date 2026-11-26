"""
Decrement paths: checks for the R6 reopen (closed-leads-revisit).

Picture: a policy can end by several causes in a year. Each cause j has an absolute rate q_j.
The point p(t) = (chance each cause has struck by time t) runs through the cube [0, q] from the
empty corner to the far corner. The total chance of ending is the k = -1 fold (the union rule)
v(p) = 1 - prod(1 - p_j). The dependent rate of cause j along a path is
    D_j = integral of  dv/dp_j (p(t)) * p_j'(t) dt.
A timing assumption between integer ages is a choice of path.

Sections
  A. Any number of causes, exactly: uniform timing (UDD) = equal corner split = Shapley value.
  B. The total does not depend on the path; the split does.
  C. Constant force = straight path in log coordinates = split in proportion to the logs.
  D. Every line dial k: the same two facts hold for v_k(S) = (prod(1 + k x) - 1)/k.
  E. Plane dial (d != 0): the straight path no longer gives Shapley. Gap = -d x y (x - y)/6 + O(d^2).
     At k = -(1+theta), d = theta this is competing risks with Ali-Mikhail-Haq dependence.
  F. Simulation of AMH-dependent causes: the path integral is the real dependent rate.

Known: UDD and constant-force formulas (Bowers et al. 1986); Shapley value (1953); equal split of
Harsanyi dividends (Harsanyi 1959); Owen's diagonal formula (Owen 1972); Aumann-Shapley (1974);
AMH copula (Ali, Mikhail, Haq 1978). Ours: the timing assumption read as a path, and section E.
Needs sympy and numpy only.
"""
import itertools, json, math
import numpy as np
import sympy as sp

results = {}
def ok(name, cond):
    results[name] = bool(cond)
    print(f"  [{'PASS' if cond else 'FAIL'}] {name}")
    return cond

t = sp.symbols('t', nonnegative=True)

def shapley_perm(v, n):
    """Shapley value by the average over all orders."""
    phi = [0]*n
    perms = list(itertools.permutations(range(n)))
    for pi in perms:
        S = frozenset()
        for j in pi:
            phi[j] += v(S | {j}) - v(S)
            S = S | {j}
    return [sp.expand(p/len(perms)) for p in phi]

print("== A. Uniform timing = equal corner split = Shapley, n = 2..5, exact ==")
for n in range(2, 6):
    q = sp.symbols(f'q1:{n+1}')
    v = lambda S: 1 - sp.prod([1 - q[i] for i in S]) if S else sp.Integer(0)
    # UDD in each single-decrement table: p_j(t) = t q_j, straight path through the cube
    udd = [sp.expand(q[j]*sp.integrate(sp.expand(sp.prod([1 - t*q[i] for i in range(n) if i != j])), (t, 0, 1)))
           for j in range(n)]
    # equal corner split: each corner S (Harsanyi dividend (-1)^{|S|+1} prod_S q) shared equally
    corner = []
    for j in range(n):
        s = 0
        others = [i for i in range(n) if i != j]
        for r in range(n):
            for T in itertools.combinations(others, r):
                s += (-1)**r * q[j]*sp.prod([q[i] for i in T]) / (r + 1)
        corner.append(sp.expand(s))
    shap = shapley_perm(v, n)
    ok(f"A{n}: n = {n}: UDD = corner split = Shapley (exact polynomials)",
       all(sp.expand(udd[j]-corner[j]) == 0 and sp.expand(udd[j]-shap[j]) == 0 for j in range(n)))
    ok(f"A{n}b: n = {n}: the shares add up to the union 1 - prod(1 - q)",
       sp.expand(sum(udd) - v(frozenset(range(n)))) == 0)

print("\n== B. Path independence of the total ==")
gl_x, gl_w = np.polynomial.legendre.leggauss(200)
tt = 0.5*(gl_x + 1); ww = 0.5*gl_w
def path_split(P, dP, vgrad):
    """D_j = int_0^1 grad_j v(P(t)) * P_j'(t) dt for a path P: t -> point in the cube."""
    pts = P(tt); vel = dP(tt)
    g = vgrad(pts)
    return (g*vel*ww).sum(axis=1)
def union_grad(pts):
    one = 1 - pts
    tot = np.prod(one, axis=0)
    return tot/one
qn = np.array([0.30, 0.50, 0.15, 0.08])
total = 1 - np.prod(1 - qn)
rng = np.random.default_rng(7)
splits = []
for trial in range(25):
    # p_j(t) = q_j t^a (back-loaded) or q_j (1 - (1-t)^a) (front-loaded), a = 1..4: smooth paths
    a = rng.integers(1, 5, size=4).astype(float); f = rng.integers(0, 2, size=4).astype(bool)
    def P(s, a=a, f=f):
        A, S = a[:, None], s[None, :]
        return qn[:, None]*np.where(f[:, None], 1 - (1 - S)**A, S**A)
    def dP(s, a=a, f=f):
        A, S = a[:, None], s[None, :]
        return qn[:, None]*A*np.where(f[:, None], (1 - S)**(A - 1), S**(A - 1))
    splits.append(path_split(P, dP, union_grad))
splits = np.array(splits)
ok("B1: 25 random paths: every total equals 1 - prod(1 - q) to 1e-10", np.all(np.abs(splits.sum(1) - total) < 1e-10))
ok("B2: ... while the split between causes moves (spread > 0.05 for the largest cause)",
   splits[:, 1].max() - splits[:, 1].min() > 0.05)
# symbolic: d/dt [1 - prod(1 - p_j(t))] = sum_j p_j' prod_{k != j}(1 - p_k)  (product rule)
p1, p2, p3 = [sp.Function(f'p{i}')(t) for i in (1, 2, 3)]
lhs = sp.diff(1 - (1-p1)*(1-p2)*(1-p3), t)
rhs = sp.diff(p1, t)*(1-p2)*(1-p3) + sp.diff(p2, t)*(1-p1)*(1-p3) + sp.diff(p3, t)*(1-p1)*(1-p2)
ok("B3: product rule identity behind it (symbolic, any paths)", sp.simplify(lhs - rhs) == 0)

print("\n== C. Constant force = straight path in log coordinates ==")
mu = -np.log(1 - qn)
P = lambda s: 1 - np.exp(-mu[:, None]*s[None, :])
dP = lambda s: mu[:, None]*np.exp(-mu[:, None]*s[None, :])
cf = path_split(P, dP, union_grad)
prop = mu/mu.sum()*total
ok("C1: constant force split = (log(1 - q_j)/sum log(1 - q)) x total", np.allclose(cf, prop, atol=1e-12))
m1, m2, m3 = sp.symbols('mu1:4', positive=True)
M = m1 + m2 + m3
Dj = sp.integrate(m1*sp.exp(-m1*t)*sp.exp(-m2*t)*sp.exp(-m3*t), (t, 0, 1))
ok("C2: symbolic: D_1 = mu_1 (1 - e^{-M})/M", sp.simplify(Dj - m1*(1 - sp.exp(-M))/M) == 0)
# known textbook formula for constant force: q_j^(tau) * log p'_j / log p^(tau)
ok("C3: same as the textbook constant-force formula", np.allclose(cf, total*np.log(1-qn)/np.log(1-total), atol=1e-12))
udd4 = path_split(lambda s: qn[:, None]*s[None, :], lambda s: qn[:, None]*np.ones_like(s)[None, :], union_grad)
print("   UDD split   ", np.round(udd4, 6))
print("   const force ", np.round(cf, 6), " total", round(total, 6))
ok("C4: the two textbook splits differ (largest cause, by more than 1e-3)", abs(udd4[1] - cf[1]) > 1e-3)

print("\n== D. Every line dial k ==")
k = sp.symbols('k')
x = sp.symbols('x1:4')
vk = lambda S: (sp.prod([1 + k*x[i] for i in S]) - 1)/k if S else sp.Integer(0)
shk = shapley_perm(lambda S: sp.cancel(vk(S)), 3)
diag = []
for j in range(3):
    g = sp.diff(sp.cancel(vk(frozenset(range(3)))), x[j]).subs({x[i]: t*x[i] for i in range(3)}, simultaneous=True)
    diag.append(sp.expand(sp.integrate(sp.expand(g*x[j]), (t, 0, 1))))
ok("D1: any k, 3 players: straight path = Shapley (Owen's diagonal formula)",
   all(sp.expand(diag[j] - shk[j]) == 0 for j in range(3)))
ok("D2: any k: Shapley = corners k^{|S|-1} prod x shared equally",
   sp.expand(shk[0] - (x[0] + k*x[0]*x[1]/2 + k*x[0]*x[2]/2 + k**2*x[0]*x[1]*x[2]/3)) == 0)
# straight path in the dial's log coordinate L_k(x) = log(1 + k x)/k -> proportional split
for kv in (-1.0, -0.4, 0.7, 2.0):
    xs = np.array([0.30, 0.50, 0.15])
    L = np.log1p(kv*xs)/kv
    tot = (np.prod(1 + kv*xs) - 1)/kv
    P = lambda s: np.expm1(kv*L[:, None]*s[None, :])/kv
    dP = lambda s: L[:, None]*np.exp(kv*L[:, None]*s[None, :])
    grad = lambda pts: np.prod(1 + kv*pts, axis=0)/(1 + kv*pts)
    sp_ = path_split(P, dP, grad)
    results.setdefault('_D3', []).append(bool(np.allclose(sp_, L/L.sum()*tot, atol=1e-12)))
ok("D3: k = -1, -0.4, 0.7, 2: straight path in log coordinates = split in proportion to L_k", all(results.pop('_D3')))

print("\n== E. Plane dial: the straight path stops being Shapley ==")
X, Y, d = sp.symbols('x y d')
Fpl = (X + Y + k*X*Y)/(1 - d*X*Y)
sh2 = (Fpl + X - Y)/2
g = sp.diff(Fpl, X).subs({X: t*X, Y: t*Y}, simultaneous=True)
diag2 = sp.integrate(sp.series(g*X, d, 0, 3).removeO(), (t, 0, 1))
gap = sp.factor(sp.expand(sp.series(sh2, d, 0, 3).removeO() - diag2))
print("   Shapley - straight path, through d^2:", gap)
ok("E1: gap = -d x y (x - y)/6 + O(d^2), the same for every k",
   sp.expand(gap - (-d*X*Y*(X - Y)*(3*d*X*Y + 5)/30)) == 0)
def plane_grad2(pts, kv, dv):
    a, b = pts
    den = 1 - dv*a*b
    return np.array([(1 + kv*b)/den + dv*b*(a + b + kv*a*b)/den**2,
                     (1 + kv*a)/den + dv*a*(a + b + kv*a*b)/den**2])
def plane_F(a, b, kv, dv): return (a + b + kv*a*b)/(1 - dv*a*b)
xs = np.array([0.30, 0.50])
for th in (0.0, 0.5, 0.9):
    kv, dv = -(1 + th), th
    stl = path_split(lambda s: xs[:, None]*s[None, :], lambda s: xs[:, None]*np.ones_like(s)[None, :],
                     lambda p: plane_grad2(p, kv, dv))
    Ftot = plane_F(*xs, kv, dv)
    shap = np.array([(Ftot + xs[0] - xs[1])/2, (Ftot + xs[1] - xs[0])/2])
    print(f"   AMH theta = {th}: total {Ftot:.6f}  straight path {np.round(stl, 6)}  Shapley {np.round(shap, 6)}")
    results[f'_E_{th}'] = (stl, shap, Ftot)
stl0, sh0, _ = results.pop('_E_0.0')
stl5, sh5, F5 = results.pop('_E_0.5')
stl9, sh9, F9 = results.pop('_E_0.9')
ok("E2: theta = 0 (independent causes): straight path = Shapley", np.allclose(stl0, sh0, atol=1e-12))
ok("E3: theta = 0.5, 0.9: totals still add up, but straight path != Shapley (gap > 1e-3)",
   abs(stl5.sum() - F5) < 1e-12 and abs(stl9.sum() - F9) < 1e-12 and abs(stl5[0] - sh5[0]) > 1e-3 and abs(stl9[0]-sh9[0]) > abs(stl5[0]-sh5[0]))
# the plane law's own log coordinate still gives a proportional split
def plane_L(xv, kv, dv):
    # integral_0^x dz/(1 + k z + d z^2), AMH case: 1 - (1+th) z + th z^2 = (1 - z)(1 - th z)
    th = dv
    return np.log((1 - th*xv)/(1 - xv))/(1 - th)
th = 0.5; kv, dv = -(1 + th), th
L = plane_L(xs, kv, dv)
Linv = lambda s: (np.exp((1 - th)*s) - 1)/(np.exp((1 - th)*s) - th)
ok("E4: the plane law's log coordinate inverts correctly", np.allclose(Linv(L), xs, atol=1e-13)
   and abs(Linv(L.sum()) - plane_F(*xs, kv, dv)) < 1e-13)
h = 1e-6
dLinv = lambda s: (Linv(s + h) - Linv(s - h))/(2*h)
P = lambda s: Linv(L[:, None]*s[None, :])
dP = lambda s: L[:, None]*dLinv(L[:, None]*s[None, :])
spl = path_split(P, dP, lambda p: plane_grad2(p, kv, dv))
ok("E5: AMH theta = 0.5: straight path in the law's log coordinate = proportional split",
   np.allclose(spl, L/L.sum()*plane_F(*xs, kv, dv), atol=1e-8))

print("\n== F. Simulation: AMH-dependent causes, uniform marginal timing ==")
# Latent times T1, T2 with uniform marginals on [0, 1/q_j] (so P(T_j <= t) = q_j t) and an AMH
# survival copula. Dependent rate of cause 1 = P(T1 <= 1, T1 < T2). Sample the copula by the
# conditional method: draw U, then V from C(v | u) by bisection.
def amh_sample(n, th, rng):
    u = rng.uniform(size=n); w = rng.uniform(size=n)
    lo = np.zeros(n); hi = np.ones(n)
    for _ in range(50):
        v = 0.5*(lo + hi)
        # dC/du for C(u, v) = u v/(1 - th (1-u)(1-v))
        den = 1 - th*(1 - u)*(1 - v)
        cond = v/den - u*v*th*(1 - v)/den**2
        hi = np.where(cond > w, v, hi); lo = np.where(cond > w, lo, v)
    return u, 0.5*(lo + hi)
th = 0.5; N = 2_000_000
U, V = amh_sample(N, th, np.random.default_rng(11))
# survival copula on survival chances: S_j(t) = 1 - q_j t, so T_j = (1 - S_j)/q_j
T1 = (1 - U)/xs[0]; T2 = (1 - V)/xs[1]
d1 = np.mean((T1 <= 1) & (T1 < T2)); d2 = np.mean((T2 <= 1) & (T2 < T1))
se = np.sqrt(np.array([d1, d2])/N)
print(f"   simulated dependent rates {d1:.5f} {d2:.5f}  path integral {np.round(stl5, 5)}  Shapley {np.round(sh5, 5)}")
ok("F1: simulation matches the straight-path integral within 3 standard errors", np.all(np.abs(np.array([d1, d2]) - stl5) < 3*se))
ok("F2: ... and rules out the Shapley split (more than 5 standard errors away)", np.all(np.abs(np.array([d1, d2]) - sh5) > 5*se))

n_pass = sum(results.values()); n_all = len(results)
print(f"\n{n_pass}/{n_all} checks pass")
json.dump(results, open(__file__.replace('checks.py', 'checks.json'), 'w'), indent=1)
