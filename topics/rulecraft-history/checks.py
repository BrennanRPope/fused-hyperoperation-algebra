"""Checks for rulecraft-history (DRAFT). Exact arithmetic (sympy/fractions) where possible; numerics marked."""
import itertools, math, random
from fractions import Fraction as Fr
import sympy as sp
import numpy as np
ok = lambda name, cond: print(("PASS " if cond else "FAIL ") + name) or cond
res = []
lam1, lam2, t1, t2, x = sp.symbols('lambda1 lambda2 t1 t2 x')
def nudge(l, t): return lambda y: (1 - l)*y + l*t
# H1 every orientation-preserving affine self-map of [0,1] with 0<A<1 is a nudge: fixed point in [0,1]
A, C = sp.symbols('A C')
p = C/(1 - A)
res.append(ok("H1 affine map Ax+C has centre C/(1-A); nudge form (1-l)x+l t with l=1-A, t=C/(1-A)",
    sp.simplify((1-(1-A))*0 + (1-A)*p - C) == 0))
# H1b maps [0,1] into [0,1] with A in (0,1) => centre in [0,1] (numeric sweep)
good = True
for _ in range(20000):
    a = random.random(); c = random.uniform(-1, 1)
    if 0 <= c <= 1 and 0 <= a + c <= 1:
        good &= (0 <= c/(1-a) <= 1)
res.append(ok("H1b any affine self-map of [0,1] with 0<A<1 has its centre inside [0,1] (20000 samples)", good))
# H2 same target: strengths fuse by the k=-1 law a+b-ab
f = nudge(lam1, t1)(nudge(lam2, t1)(x))
res.append(ok("H2 same target: composite strength = l1 + l2 - l1 l2 (our law at k=-1)",
    sp.simplify(f - ((1 - (lam1 + lam2 - lam1*lam2))*x + (lam1 + lam2 - lam1*lam2)*t1)) == 0))
# H3 swap gap
g1 = nudge(lam2, t2)(nudge(lam1, t1)(x)); g2 = nudge(lam1, t1)(nudge(lam2, t2)(x))
res.append(ok("H3 (rule 1 then 2) - (rule 2 then 1) = l1 l2 (t2 - t1), whatever x", sp.simplify(g1 - g2 - lam1*lam2*(t2 - t1)) == 0))
# H4 n rules: total pull order-free, weights are stick-breaking
n = 5
ls = sp.symbols('l0:%d' % n); ts = sp.symbols('s0:%d' % n)
y = x
for i in range(n): y = nudge(ls[i], ts[i])(y)
y = sp.expand(y)
coeffx = y.coeff(x)
res.append(ok("H4a starting point keeps share prod(1 - l_i), the same for every order", sp.simplify(coeffx - sp.prod([1 - l for l in ls])) == 0))
w_ok = all(sp.simplify(y.coeff(ts[i]) - ls[i]*sp.prod([1 - ls[j] for j in range(i+1, n)])) == 0 for i in range(n))
res.append(ok("H4b weight on target i = l_i * prod over LATER rules (1 - l_j): stick-breaking read backward in time", w_ok))
tot = sp.simplify(sum(y.coeff(ts[i]) for i in range(n)))
fold = ls[0]
for l in ls[1:]: fold = fold + l - fold*l
res.append(ok("H4c weights sum to the k=-1 fold of all strengths (order-free)", sp.simplify(tot - fold) == 0))
# H5 best order = targets ascending, for any strengths; outcome strictly monotone on the weak order
rng = random.Random(7)
def outcome(order, L, T, x0):
    y = x0
    for i in order: y = (1 - L[i])*y + L[i]*T[i]
    return y
mono = True; best = True
for trial in range(300):
    m = 5
    L = [Fr(rng.randint(1, 99), 100) for _ in range(m)]; T = [Fr(rng.randint(0, 1000), 1000) for _ in range(m)]
    x0 = Fr(rng.randint(0, 100), 100)
    outs = {o: outcome(o, L, T, x0) for o in itertools.permutations(range(m))}
    asc = tuple(sorted(range(m), key=lambda i: T[i])); desc = asc[::-1]
    best &= outs[asc] == max(outs.values()) and outs[desc] == min(outs.values())
    for o, v in outs.items():
        for k in range(m - 1):
            if T[o[k]] < T[o[k+1]]:
                o2 = o[:k] + (o[k+1], o[k]) + o[k+2:]
                mono &= outs[o2] < v
res.append(ok("H5a highest outcome = rules sorted by target ascending, lowest = descending (300 random sets of 5, exact)", best))
res.append(ok("H5b putting a higher target earlier always lowers the outcome: strictly monotone on the weak order of S5", mono))
# H6 second order (BCH): order effect = 1/2 sum over pairs l_i l_j (t_later - t_earlier) + O(l^3)
eps = sp.symbols('epsilon')
L3 = [eps*sp.Rational(a) for a in (2, 3, 5)]; T3 = [sp.Rational(1, 7), sp.Rational(4, 5), sp.Rational(1, 3)]
y0 = sp.Rational(1, 2)
def out_sym(order):
    y = y0
    for i in order: y = (1 - L3[i])*y + L3[i]*T3[i]
    return sp.expand(y)
avg = sum(out_sym(o) for o in itertools.permutations(range(3)))/6
o = (0, 1, 2)
pred = sum(sp.Rational(1, 2)*L3[o[i]]*L3[o[j]]*(T3[o[j]] - T3[o[i]]) for i in range(3) for j in range(i+1, 3))
diff = sp.series(out_sym(o) - avg - pred, eps, 0, 3).removeO()
res.append(ok("H6 order effect at second order = 1/2 sum_pairs l_i l_j (t_later - t_earlier) (BCH term)", sp.simplify(diff) == 0))
# H7 gradualism: interleave n small steps; gap shrinks like 1/n, limit = simultaneous
def interleave(order, L, T, x0, nn):
    y = x0
    for _ in range(nn):
        for i in order:
            li = 1 - (1 - L[i])**(1.0/nn)
            y = (1 - li)*y + li*T[i]
    return y
L = [0.6, 0.5]; T = [0.1, 0.9]; x0 = 0.5
gaps = [abs(interleave((0, 1), L, T, x0, nn) - interleave((1, 0), L, T, x0, nn)) for nn in (1, 2, 4, 8, 16, 32, 64)]
ratios = [gaps[i]/gaps[i+1] for i in range(len(gaps)-1)]
mu = [-math.log(1 - l) for l in L]
sim_target = sum(m*t for m, t in zip(mu, T))/sum(mu); sim = sim_target + (x0 - sim_target)*math.exp(-sum(mu))
lim = interleave((0, 1), L, T, x0, 20000)
res.append(ok("H7a interleaving in n steps: gap ratios -> 2 per doubling (gap ~ 1/n) [numeric] %s" % [round(r, 3) for r in ratios], abs(ratios[-1] - 2) < 0.02))
res.append(ok("H7b limit = simultaneous rule: rates add, target = rate-weighted mean [numeric] %.6f vs %.6f" % (lim, sim), abs(lim - sim) < 1e-4))
# H8 pendulum: alternating two nudges -> 2-cycle with swing l1 l2 |t1-t2| / (l1 + l2 - l1 l2)
xs = sp.symbols('xs')
fp1 = sp.solve(sp.Eq(nudge(lam2, t2)(nudge(lam1, t1)(xs)), xs), xs)[0]
fp2 = sp.solve(sp.Eq(nudge(lam1, t1)(nudge(lam2, t2)(xs)), xs), xs)[0]
res.append(ok("H8 alternating rules settle into a 2-cycle; swing = l1 l2 (t2-t1)/(l1+l2-l1 l2) = order gap / total pull",
    sp.simplify(fp1 - fp2 - lam1*lam2*(t2 - t1)/(lam1 + lam2 - lam1*lam2)) == 0))
# H9 coin-flip elections: stationary law is a Bernoulli convolution with ratio 1-l (numeric support check)
def support_gaps(l, depth=14):
    pts = np.array([0.0])
    for _ in range(depth):
        pts = np.unique(np.round(np.concatenate([(1 - l)*pts, (1 - l)*pts + l]), 12))
    return pts
s7 = support_gaps(0.7, 12); gap_mid = (s7 > 0.3 + 1e-9) & (s7 < 0.7 - 1e-9)
res.append(ok("H9a strength 0.7: no reachable position strictly between 0.3 and 0.7 (the forbidden middle)", not gap_mid.any()))
ll = sp.Rational(7, 10)
img0 = (nudge(ll, 0)(0), nudge(ll, 0)(1)); img1 = (nudge(ll, 1)(0), nudge(ll, 1)(1))
res.append(ok("H9b strength l > 1/2: the two rules map [0,1] onto [0,1-l] and [l,1], leaving (1-l, l) empty; at l = 0.7: %s, %s" % (img0, img1),
    img0 == (0, 1 - ll) and img1 == (ll, 1)))
s8 = support_gaps(1 - (math.sqrt(5) - 1)/2, 14)  # strength 1 - 1/phi
res.append(ok("H9c strength 1 - 1/phi (about 0.382): positions overlap (no forbidden band) yet the law is singular (Erdos 1939, cited not checked); %d distinct points at depth 14 vs 2^14 = 16384" % len(s8), len(s8) < 16384))
# H10 repeal window: undoing rule i at the end needs x in [l t, 1 - l(1-t)] (width 1-l)
lam, t, xx = sp.symbols('lam t xx')
inv = (xx - lam*t)/(1 - lam)
res.append(ok("H10 repeal-by-inverse stays in [0,1] iff l t <= x <= 1 - l(1-t); window width 1-l",
    sp.simplify(inv.subs(xx, lam*t)) == 0 and sp.simplify(inv.subs(xx, 1 - lam*(1 - t)) - 1) == 0))
# H11 repeal != never having had it: difference = sum of later order gaps
L = [Fr(3, 10), Fr(1, 2), Fr(2, 5)]; T = [Fr(9, 10), Fr(1, 10), Fr(3, 10)]; x0 = Fr(1, 2)
full = outcome((0, 1, 2), L, T, x0); never = outcome((1, 2), L, T, x0)
repeal = (full - L[0]*T[0])/(1 - L[0])
# moving rule 0 to the end costs the swaps with rules 1 and 2
moved = outcome((1, 2, 0), L, T, x0)
res.append(ok("H11 repeal now = never had it exactly when rule 0 is moved to the end first; mismatch = %s" % (repeal - never),
    (moved - L[0]*T[0])/(1 - L[0]) == never and repeal != never))
# H12 memory dial: l_n = c/(n+n0)^alpha (step-decay exponent, Polyak-Juditsky notation); origin share prod(1-l_n)
def origin_share(alpha, c=0.8, n0=1, N=200000):
    s = 0.0
    for nn in range(N): s += math.log(1 - c/(nn + n0 + 1)**alpha)
    return math.exp(s)
shares = {b: origin_share(b) for b in (0.5, 1.0, 1.5, 2.0)}
res.append(ok("H12 origin forgotten iff sum l_n = inf: share at alpha=0.5,1 -> ~0; alpha=1.5,2 stays >0 [numeric] %s" % {k: f"{v:.2e}" for k, v in shares.items()},
    shares[0.5] < 1e-6 and shares[1.0] < 1e-4 and shares[1.5] > 1e-3 and shares[2.0] > 0.1))
# H13 Polya: l_n = 1/(n+1) gives the plain average of all targets (every rule weighs the same)
T = [Fr(rng.randint(0, 100), 100) for _ in range(12)]
y = Fr(0)
for nn, tt in enumerate(T): y = (1 - Fr(1, nn + 1))*y + Fr(1, nn + 1)*tt
res.append(ok("H13 strength 1/n makes the slider the plain average of every target so far (Polya / Arthur)", y == sum(T)/len(T)))
# H14 hyperbolic: in (dose D, target t) the invariant metric is dD^2 + (t dD + (e^D - 1) dt)^2, curvature -1
D, tt = sp.symbols('D t', positive=True)
E_ = 1 + tt**2; F_ = tt*(sp.exp(D) - 1); G_ = (sp.exp(D) - 1)**2
# Brioschi via sympy: compute Gaussian curvature of metric E dD^2 + 2F dD dt + G dt^2
def gauss(E, F, G, u, v):
    Eu, Ev, Fu, Fv, Gu, Gv = [sp.diff(a, b) for a in (E, F, G) for b in (u, v)]
    Evv = sp.diff(E, v, 2); Guu = sp.diff(G, u, 2); Fuv = sp.diff(F, u, v)
    M1 = sp.Matrix([[-Evv/2 + Fuv - Guu/2, Eu/2, Fu - Ev/2], [Fv - Gu/2, E, F], [Gv/2, F, G]])
    M2 = sp.Matrix([[0, Ev/2, Gu/2], [Ev/2, E, F], [Gu/2, F, G]])
    return sp.simplify((M1.det() - M2.det())/(E*G - F**2)**2)
K = gauss(E_, F_, G_, D, tt)
res.append(ok("H14 reform space in (dose, target) coordinates has curvature %s" % K, sp.simplify(K + 1) == 0))
# H15-H17 rule families by anchors (the (k,d) law): x -> (x + c + k x c)/(1 - d x c)
k_, d_, c_, c2_, X = sp.symbols('k d c c2 X')
law = lambda a, b, k, d: (a + b + k*a*b)/(1 - d*a*b)
roots = sp.solve(1 + k_*X + d_*X**2, X)
res.append(ok("H15 the map x -> x (+) c fixes exactly the roots of 1 + k x + d x^2 (two when Delta = k^2 - 4d > 0, none when < 0)",
    all(sp.simplify(law(r, c_, k_, d_) - r) == 0 for r in roots)))
comm = sp.simplify(law(law(X, c_, k_, d_), c2_, k_, d_) - law(law(X, c2_, k_, d_), c_, k_, d_))
k2, d2 = sp.Rational(1), sp.Rational(-2)
mixed = sp.simplify(law(law(X, c_, -1, 0), c2_, k2, d2) - law(law(X, c2_, k2, d2), c_, -1, 0))
res.append(ok("H16 rules from one family (same anchors) commute; a nudge-to-1 rule and an odds rule with other anchors do not", comm == 0 and mixed != 0))
res.append(ok("H17 x (+) c at k = -1, d = 0 is the nudge toward 1 with strength c", sp.simplify(law(X, c_, -1, 0) - ((1 - c_)*X + c_)) == 0))
print("\n%d of %d checks pass" % (sum(res), len(res)))
