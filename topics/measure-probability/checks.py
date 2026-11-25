"""Measure theory and probability through the circle operation: every claim on index.html, checked.
Run: python3 checks.py > checks_output.txt"""
import itertools, math
import numpy as np, sympy as sp
from scipy.optimize import brentq
from scipy.integrate import quad

rng = np.random.default_rng(7)
results = []
def check(name, ok, detail=""):
    results.append(bool(ok))
    print(f"{'PASS' if ok else 'FAIL'}  {name}" + (f"   [{detail}]" if detail else ""))

dial = lambda k, a, b: a + b + k*a*b
def lam_measure(k, gs):           # fold of the dial over the atoms
    m = 0
    for g in gs: m = dial(k, m, g)
    return m
def subsets(n):
    for r in range(n+1):
        for S in itertools.combinations(range(n), r): yield S

def sugeno_lambda(g):
    """Nonzero root of prod(1+l g_i) = 1+l in (-1, inf)."""
    f = lambda l: np.prod(1 + l*np.asarray(g)) - (1 + l)
    s = sum(g)
    if abs(s-1) < 1e-12: return 0.0
    if s > 1: return brentq(f, -1+1e-12, -1e-12)
    hi = 1.0
    while f(hi) < 0: hi *= 2
    return brentq(f, 1e-12, hi)

print("== 1. The Sugeno lambda-measure is the dial with k = lambda")
x, y, z, l = sp.symbols('x y z l')
check("symbolic: disjoint union a ✳λ b is associative and commutative",
      sp.expand(dial(l, dial(l, x, y), z) - dial(l, x, dial(l, y, z))) == 0 and sp.expand(dial(l, x, y) - dial(l, y, x)) == 0)
check("symbolic: lift 1+λ(a✳λb) = (1+λa)(1+λb)", sp.expand(1 + l*dial(l, x, y) - (1+l*x)*(1+l*y)) == 0)
ok = True
for _ in range(300):
    n = rng.integers(2, 7); g = rng.uniform(0.02, 0.8, n); lam = sugeno_lambda(g)
    mu = {S: lam_measure(lam, [g[i] for i in S]) for S in subsets(n)}
    full = tuple(range(n))
    ok &= abs(mu[full] - 1) < 1e-9
    for A in subsets(n):
        for B in subsets(n):
            if set(A) & set(B): continue
            U = tuple(sorted(set(A) | set(B)))
            ok &= abs(mu[U] - dial(lam, mu[A], mu[B])) < 1e-9
check("300 random spaces: μ(A∪B) = μ(A)+μ(B)+λμ(A)μ(B) for every disjoint pair, μ(Ω)=1", ok)

print("== 2. The normalization equation ∏(1+λg_i) = 1+λ")
ok_unique = ok_sign = True
for _ in range(400):
    n = rng.integers(2, 8); g = rng.uniform(0.01, 0.9, n)
    poly = np.poly1d([1.0])
    for gi in g: poly = poly * np.poly1d([gi, 1.0])
    poly = poly - np.poly1d([1.0, 1.0])
    roots = [r.real for r in poly.r if abs(r.imag) < 1e-9 and r.real > -1 + 1e-12 and abs(r.real) > 1e-9]
    ok_unique &= (len(roots) == 1)
    lam = sugeno_lambda(g)
    ok_sign &= ((lam > 0) == (g.sum() < 1))
check("400 random densities: exactly one nonzero root λ > −1", ok_unique)
check("λ > 0 exactly when the densities sum below 1 (λ < 0 when above)", ok_sign)
g1, g2 = sp.symbols('g1 g2', positive=True)
l2 = (1 - g1 - g2)/(g1*g2)
check("two atoms: λ = (1 − g1 − g2)/(g1 g2) solves it", sp.simplify((1+l2*g1)*(1+l2*g2) - (1+l2)) == 0)

print("== 3. Möbius masses are the corners of the circle-operation cube")
ok = ok_bel = ok_dualbel = True
for _ in range(200):
    n = rng.integers(2, 6); g = rng.uniform(0.05, 0.8, n); lam = sugeno_lambda(g)
    m = {S: (lam**(len(S)-1) * np.prod([g[i] for i in S]) if S else 0.0) for S in subsets(n)}
    for A in subsets(n):
        ok &= abs(sum(m[S] for S in subsets(n) if set(S) <= set(A)) - lam_measure(lam, [g[i] for i in A])) < 1e-9
    if lam >= 0: ok_bel &= all(v >= -1e-15 for v in m.values())
    else:
        # dual μ*(A) = 1 − μ(complement); its Möbius masses
        mu = {S: lam_measure(lam, [g[i] for i in S]) for S in subsets(n)}
        full = set(range(n))
        dual = {S: 1 - mu[tuple(sorted(full - set(S)))] for S in subsets(n)}
        md = {A: sum((-1)**(len(A)-len(B))*dual[B] for B in subsets(n) if set(B) <= set(A)) for A in subsets(n)}
        ok_dualbel &= all(v >= -1e-12 for v in md.values())
check("200 spaces: μ(A) = Σ over corners S ⊆ A of λ^(|S|−1)·∏g_S", ok)
check("λ ≥ 0: every Möbius mass ≥ 0, so μ is a belief function", ok_bel)
check("−1 < λ < 0: the dual 1 − μ(Aᶜ) has masses ≥ 0 (μ is a plausibility)", ok_dualbel)

print("== 4. The lift to an ordinary probability; μ is a q-number of P")
ok_add = ok_q = True
for _ in range(200):
    n = rng.integers(2, 6); g = rng.uniform(0.05, 0.8, n); lam = sugeno_lambda(g)
    if abs(lam) < 1e-6: continue
    P = lambda S: math.log(1 + lam*lam_measure(lam, [g[i] for i in S]))/math.log(1 + lam)
    ok_add &= abs(P(tuple(range(n))) - 1) < 1e-9 and abs(sum(P((i,)) for i in range(n)) - 1) < 1e-9
    q = 1 + lam
    for S in subsets(n):
        ok_q &= abs(lam_measure(lam, [g[i] for i in S]) - (q**P(S) - 1)/(q - 1)) < 1e-9
check("P(A) = log(1+λμ(A))/log(1+λ) is additive with P(Ω) = 1", ok_add)
check("μ(A) = (q^P(A) − 1)/(q − 1) with q = 1 + λ = k + 1", ok_q)
ok = True
for m_ in range(1, 9):
    for (lam, g) in [(0.7, 0.2), (-0.4, 0.3), (2.0, 0.05)]:
        qq = 1 + lam*g
        ok &= abs(lam_measure(lam, [g]*m_) - g*sum(qq**j for j in range(m_))) < 1e-12
check("m equal atoms of size g: μ = g·[m]_q with q = 1 + λg (a q-number)", ok)

print("== 5. Inclusion–exclusion becomes multiplication")
ok2 = okn = True
for _ in range(200):
    n = rng.integers(3, 6); g = rng.uniform(0.05, 0.7, n); lam = sugeno_lambda(g)
    mu = lambda S: lam_measure(lam, [g[i] for i in S])
    A = [i for i in range(n) if rng.random() < .6]; B = [i for i in range(n) if rng.random() < .6]
    U = sorted(set(A) | set(B)); I = sorted(set(A) & set(B))
    ok2 &= abs((1+lam*mu(U))*(1+lam*mu(I)) - (1+lam*mu(A))*(1+lam*mu(B))) < 1e-9
    sets = [[i for i in range(n) if rng.random() < .5] for _ in range(3)]
    lhs = 1 + lam*mu(sorted(set().union(*sets)))
    rhs = 1.0
    for r in range(1, 4):
        for c in itertools.combinations(sets, r):
            rhs *= (1 + lam*mu(sorted(set.intersection(*map(set, c)))))**((-1)**(r+1))
    okn &= abs(lhs - rhs) < 1e-9
check("(1+λμ(A∪B))(1+λμ(A∩B)) = (1+λμ(A))(1+λμ(B))", ok2)
check("three sets: 1+λμ(∪) = ∏ over corners of (1+λμ(∩))^(±1)", okn)

print("== 6. Duality: the dual weight is the circle operation's inverse")
ok = True
for _ in range(200):
    n = rng.integers(2, 6); g = rng.uniform(0.05, 0.8, n); lam = sugeno_lambda(g)
    full = set(range(n)); mu = lambda S: lam_measure(lam, [g[i] for i in S])
    dual = lambda S: 1 - mu(sorted(full - set(S)))
    ld = -lam/(1+lam)
    for A in subsets(n):
        for B in subsets(n):
            if set(A) & set(B): continue
            ok &= abs(dual(sorted(set(A)|set(B))) - dial(ld, dual(A), dual(B))) < 1e-9
check("200 spaces: μ*(A) = 1 − μ(Aᶜ) is a λ*-measure with λ* = −λ/(1+λ)", ok)
check("symbolic: λ* = circle inverse of λ, (λ*)* = λ, 1+λ* = 1/(1+λ)",
      sp.simplify(dial(1, l, -l/(1+l))) == 0 and sp.simplify(-(-l/(1+l))/(1-l/(1+l)) - l) == 0 and sp.simplify(1 - l/(1+l) - 1/(1+l)) == 0)

print("== 7. What the weight means: firing atoms and compounding")
ok = True
for _ in range(100):
    n = rng.integers(2, 6); g = rng.uniform(0.3, 0.9, n)
    if g.sum() <= 1: continue
    lam = sugeno_lambda(g); c = -lam*g      # chance each atom fires
    ok &= np.all((c > 0) & (c < 1))
    pany = 1 - np.prod(1 - c)
    for S in subsets(n):
        # exact enumeration over firing patterns
        pS = sum(np.prod([c[i] if f[i] else 1-c[i] for i in range(n)]) for f in itertools.product([0,1], repeat=n) if any(f[i] for i in S))
        ok &= abs(lam_measure(lam, [g[i] for i in S]) - pS/pany) < 1e-9
check("−1<λ<0: μ(A) = P(some atom of A fires)/P(any fires), atoms fire independently w.p. −λg_i", ok)
ok = True
for _ in range(100):
    n = rng.integers(2, 6); g = rng.uniform(0.02, 0.3, n)
    if g.sum() >= 1: continue
    lam = sugeno_lambda(g); grow = lambda S: np.prod([1+lam*g[i] for i in S]) - 1
    ok &= all(abs(lam_measure(lam, [g[i] for i in S]) - grow(S)/grow(range(n))) < 1e-9 for S in subsets(n))
check("λ>0: μ(A) = (compound growth of A's atoms)/(growth of all), atom rates λg_i", ok)

print("== 8. The wall must stay outside [0, 1]")
ok = True
for k in [-1.0, -1.5, -3.0]:
    for _ in range(200):
        n = rng.integers(1, 8); w = -1/k
        g = rng.uniform(0, w*(1-1e-9), n)
        ok &= lam_measure(k, g) < w - 1e-15 and lam_measure(k, g) < 1
check("k ≤ −1: atoms below the wall never reach 1 (no normalized measure)", ok)
from fractions import Fraction as Fr
check("k = −1: 1 − ∏(1 − g_i) = 1 forces some g_i = 1 (exact fractions)",
      lam_measure(-1, [Fr(3,10), Fr(1), Fr(1,2)]) == 1 and lam_measure(-1, [Fr(999,1000)]*50) < 1)
def lens(p, q, a, b): return (a + b + p*a*b)/(1 - q*a*b)
def L1(p, q):   # additive coordinate at 1: ∫0^1 dt/(1+pt+qt²)
    f = lambda t: 1/(1 + p*t + q*t*t)
    ts = np.linspace(0, 1, 2001)
    if np.min(1 + p*ts + q*ts*ts) <= 1e-9: return math.inf
    return quad(f, 0, 1, limit=200)[0]
def normalizable(p, q):
    ts = np.linspace(0, 1, 4001); return np.min(1 + p*ts + q*ts*ts) > 1e-9
ok = True
for _ in range(400):
    p, q = rng.uniform(-3, 3, 2)
    if abs(np.min(1 + p*np.linspace(0,1,4001) + q*np.linspace(0,1,4001)**2)) < 1e-3: continue
    ok &= (math.isfinite(L1(p, q)) == normalizable(p, q))
check("landscape: L(1) = ∫₀¹dt/(1+pt+qt²) finite exactly when 1+pt+qt² > 0 on [0,1]", ok)
th = np.linspace(-0.99, 0.99, 50)
check("the 'or' line k+q = −1 (statistics page) lies on the edge: 1+p+q = 0 there",
      all(abs(1 + (-(1+t)) + t) < 1e-15 for t in th))
def region(p, q): return (q > -1 - p) if p >= -2 else (q > p*p/4)
ok = True
for _ in range(2000):
    p, q = rng.uniform(-6, 3), rng.uniform(-3, 9)
    ts = np.linspace(0, 1, 20001); mn = np.min(1 + p*ts + q*ts*ts)
    if abs(mn) < 1e-3: continue
    ok &= ((mn > 0) == region(p, q))
check("normalizable region = above the line q = −1−p (p ≥ −2) and above the parabola q = p²/4 (p < −2)", ok)
pp = sp.symbols('p')
check("the line meets the border parabola tangentially at the cusp (−2, 1)",
      sp.solve(sp.Eq(-1 - pp, pp**2/4), pp) == [-2] and sp.diff(pp**2/4, pp).subs(pp, -2) == -1)
# a normalizable bent example: tangent addition (0, 1)
ok = True
for _ in range(100):
    n = rng.integers(2, 6); P = rng.dirichlet(np.ones(n))
    mu = lambda S: math.tan(math.pi/4*sum(P[i] for i in S))
    ok &= abs(mu(range(n)) - 1) < 1e-12
    A = [i for i in range(n) if rng.random() < .5]; B = [i for i in range(n) if i not in A and rng.random() < .5]
    ok &= abs(mu(A+B) - lens(0, 1, mu(A), mu(B))) < 1e-12
check("tangent measure μ = tan(πP/4): (μA+μB)/(1−μAμB) on disjoint sets, μ(Ω)=1", ok)
ok = True
for k in [-0.9, -0.5, 0.5, 3.0]:
    gs = rng.uniform(0.01, 0.2, 5); m = lam_measure(k, gs)
    Lk = lambda x: math.log(1 + k*x)/k
    ok &= abs(Lk(m) - sum(Lk(gi) for gi in gs)) < 1e-12
check("dial line: the additive coordinate log(1+kx)/k adds over atoms", ok)

print("== 9. Odds of an 'or' combine by the circle operation")
ok = True
odds = lambda p: p/(1-p); against = lambda p: (1-p)/p
for _ in range(500):
    p, q, r = rng.uniform(0.01, 0.99, 3)
    ok &= abs(odds(dial(-1, p, q)) - dial(1, odds(p), odds(q))) < 1e-8*max(1, odds(dial(-1,p,q)))
    ok &= abs(against(p*q) - dial(1, against(p), against(q))) < 1e-8*max(1, against(p*q))
    ok &= abs(odds(dial(-1, dial(-1, p, q), r)) - dial(1, dial(1, odds(p), odds(q)), odds(r))) < 1e-7*max(1, odds(dial(-1,dial(-1,p,q),r)))
check("independent A, B: odds(A∪B) = odds A ✳ odds B; odds-against(A∩B) = against A ✳ against B; 3 events", ok)
bad = sum(abs(odds(p*q) - dial(1, odds(p), odds(q))) > 1e-6 for p, q in rng.uniform(0.05, 0.95, (200, 2)))
check("odds of an 'and' do NOT the circle operation (dropped claim): fails on 200 of 200 random pairs", bad == 200)
# Monte Carlo
N = 400000; A = rng.random(N) < .3; B = rng.random(N) < .45
pe = (A|B).mean()
check("simulation (400,000 draws): odds of A∪B matches 0.4286 ✳ 0.8182 = 1.5974",
      abs(pe/(1-pe) - dial(1, odds(.3), odds(.45))) < 0.02, f"{pe/(1-pe):.4f}")

print("== 10. Expectation passes through the dial when variables are independent")
ok = okcov = True
for _ in range(300):
    k = rng.uniform(-2, 2)
    xs, ys = rng.normal(size=5), rng.normal(size=4)
    px, py = rng.dirichlet(np.ones(5)), rng.dirichlet(np.ones(4))
    E = sum(px[i]*py[j]*dial(k, xs[i], ys[j]) for i in range(5) for j in range(4))
    ok &= abs(E - dial(k, px@xs, py@ys)) < 1e-10
    J = rng.dirichlet(np.ones(20)).reshape(5, 4)     # dependent joint
    EX = (J.sum(1))@xs; EY = (J.sum(0))@ys; EXY = xs@J@ys
    E2 = sum(J[i, j]*dial(k, xs[i], ys[j]) for i in range(5) for j in range(4))
    okcov &= abs(E2 - dial(k, EX, EY) - k*(EXY - EX*EY)) < 1e-10
check("300 random independent pairs: E[X ✳ₖ Y] = E X ✳ₖ E Y", ok)
check("dependent pairs: the gap is exactly k·Cov(X, Y) (the corner's excess)", okcov)
r = [(-0.5, 0.5), (0.5, 0.5)]           # fair ±50% bet
vals = list(itertools.product(*[[v for v, _ in r]]*6))
EW = sum(np.prod([0.5]*6) * (np.prod([1+v for v in path]) - 1) for path in vals)
med = np.median([np.prod([1+v for v in path]) - 1 for path in vals])
check("martingale: six fair ±50% bets, E[r₁ ✳ … ✳ r₆] = 0 exactly", abs(EW) < 1e-12)
check("yet the typical path loses: median wealth change = 0.75³ − 1 = −0.578", abs(med - (0.75**3 - 1)) < 1e-12)
check("volatility drag: E log(1+r) = ½log(0.75) < 0 = log(1 + E r)", 0.5*math.log(1.5)+0.5*math.log(0.5) < 0)

print("== 11. Yor's formula: stochastic exponentials multiply by the circle operation, step by step")
ok = True
for _ in range(200):
    n = rng.integers(1, 30); dx, dy = rng.normal(0, .1, n), rng.normal(0, .1, n)
    lhs = (np.prod(1+dx) - 1); rhs = (np.prod(1+dy) - 1)
    ok &= abs(dial(1, lhs, rhs) - (np.prod(1 + dial(1, dx, dy)) - 1)) < 1e-10
check("ℰ(X)ℰ(Y) = ℰ(X+Y+[X,Y]): per step Δx ✳ Δy = Δx+Δy+ΔxΔy", ok)
T, sig = 1.0, 0.8; errs = []; qv = []
for n in [100, 1000, 10000, 100000]:
    e1 = e2 = 0.0
    for _ in range(20):
        dW = rng.normal(0, math.sqrt(T/n), n); W = dW.sum()
        e1 += abs(np.prod(1 + sig*dW) - math.exp(sig*W - sig**2*T/2))/20
        e2 += abs((dW**2).sum() - T)/20
    errs.append(e1); qv.append(e2)
check("Brownian (20 paths each): ∏(1+σΔW) → exp(σW − σ²T/2), mean error falls with finer steps", errs[-1] < 0.01 and all(a > b for a, b in zip(errs, errs[1:])), ", ".join(f"{e:.4f}" for e in errs))
check("the summed self-corners ΣΔW² → T (the Itô drag is half of them)", qv[-1] < 0.01 and all(a > b for a, b in zip(qv, qv[1:])), ", ".join(f"{e:.4f}" for e in qv))

print("== 12. Discrete hazards (already on the information page; restated)")
h = rng.uniform(0, .2, 12); F = lam_measure(-1, h)
check("failure by step n = h₁ ✳₋₁ … ✳₋₁ hₙ = 1 − ∏(1 − hᵢ)", abs(F - (1 - np.prod(1-h))) < 1e-14)

print(f"\n{sum(results)}/{len(results)} checks pass")
