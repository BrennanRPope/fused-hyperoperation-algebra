"""Game theory through mult-add: every claim on the page, checked.
Dial: a ✳k b = a + b + k·ab, lift 1 + k·x."""
import math, random, itertools
from fractions import Fraction as F
import numpy as np
from scipy.optimize import minimize_scalar
from scipy.integrate import quad
OK = [0, 0]
def check(name, cond, extra=""):
    OK[1] += 1; OK[0] += bool(cond)
    print(("PASS " if cond else "FAIL ") + name + (("  " + extra) if extra else ""))
def dial(k, xs):
    acc = 0
    for x in xs: acc = acc + x + k * acc * x
    return acc
random.seed(11)

print("== 1. Multiplicative weights (Hedge): weights live in the lift, the dial is the learning rate")
def hedge(G, eps):
    n = G.shape[1]; w = np.ones(n); r = []
    for g in G:
        p = w / w.sum(); r.append(float(p @ g)); w = w * (1 + eps * g)
    return np.array(r), w
ok1 = ok2 = ok3 = ok4 = ok5 = True; worst_ratio = 0
for trial in range(300):
    n = random.choice([2, 3, 5, 10]); T = random.choice([50, 200, 1000]); eps = random.choice([0.05, 0.1, 0.3, 0.5])
    kind = trial % 3
    if kind == 0: G = np.random.default_rng(trial).uniform(-1, 1, (T, n))
    elif kind == 1:   # one good expert hidden among noise
        G = np.random.default_rng(trial).uniform(-1, 1, (T, n)); G[:, 0] = np.clip(G[:, 0] + 0.3, -1, 1)
    else:             # adversarial-ish: experts alternate who is good, switching every 37 rounds
        G = -np.ones((T, n)) * 0.5
        for t in range(T): G[t, (t // 37) % n] = 1.0
    r, w = hedge(G, eps)
    W = w.sum()
    ok1 &= np.allclose(w, np.prod(1 + eps * G, axis=0))
    ok1 &= np.allclose(1 + eps * np.array([dial(eps, G[:, i]) for i in range(n)]), w)
    ok2 &= math.isclose(W, n * np.prod(1 + eps * r), rel_tol=1e-9)
    algL = np.sum(np.log1p(eps * r)); bestL = np.max(np.sum(np.log1p(eps * G), axis=0))
    ok3 &= algL >= bestL - math.log(n) - 1e-9
    G_alg = r.sum(); G_best = G.sum(axis=0).max()
    regret = G_best - G_alg
    bound = math.log(n) / eps + eps * T
    ok4 &= regret <= bound + 1e-9
    worst_ratio = max(worst_ratio, regret / bound)
    # tuned rate
    epst = min(0.5, math.sqrt(math.log(n) / T)); r2, _ = hedge(G, epst)
    ok5 &= (G_best - r2.sum()) <= 2 * math.sqrt(T * math.log(n)) + 1e-9
check("each expert's weight is prod(1 + eps·g) = 1 + eps·(g1 ✳eps g2 ✳eps ...): the dial value is the learning rate", ok1)
check("the total weight is n·prod(1 + eps·r_t): the algorithm's own gains fold under the same dial", ok2)
check("lifted regret: sum log(1+eps r) >= best expert's sum log(1+eps g) - ln n (300 runs, three kinds of gains)", ok3)
check("plain regret <= ln n / eps + eps·T for gains in [-1, 1], eps <= 1/2", ok4, "worst regret/bound %.2f" % worst_ratio)
check("with eps = sqrt(ln n / T): regret <= 2 sqrt(T ln n)", ok5)
okq = all(e*g - (e*g)**2 <= math.log1p(e*g) <= e*g for e in [0.01, 0.1, 0.5] for g in np.linspace(-1, 1, 41))
check("the corner cost: eps·g - (eps·g)^2 <= log(1 + eps·g) <= eps·g for |eps·g| <= 1/2", okq)

print("== 2. Kelly betting: bets compose by the dial with k = f, ruin is the wall")
okk = True
for _ in range(100):
    f = random.uniform(0.01, 0.9); rs = [random.uniform(-1, 3) for _ in range(6)]
    okk &= math.isclose(np.prod([1 + f * r for r in rs]), 1 + f * dial(f, rs), rel_tol=1e-12)
check("wealth after bets r1..rn at fraction f is 1 + f·(r1 ✳f r2 ✳f ... ✳f rn)", okk)
f = 0.4
check("ruin is the dial's wall: an outcome r = -1/f leaves wealth 0 whatever comes after",
      all(F(1) + F(2, 5) * dial(F(2, 5), [F(-5, 2)] + [F(x) for x in xs]) == 0 for xs in [[1, 2], [-1, 3, 5]]))
okkel = True; okkl = True; okover = []
for p, b in [(0.55, 1.0), (0.6, 1.0), (0.3, 4.0), (0.52, 1.0), (0.7, 0.8)]:
    q = 1 - p; g = lambda f: p * math.log1p(b * f) + q * math.log1p(-f)
    fs = p - q / b
    res = minimize_scalar(lambda f: -g(f), bounds=(0, 0.999), method="bounded", options={"xatol": 1e-12})
    okkel &= abs(res.x - fs) < 1e-6
    pi = 1 / (1 + b)
    okkl &= math.isclose(1 + b * fs, p / pi) and math.isclose(1 - fs, q / (1 - pi))
    kl = p * math.log(p / pi) + q * math.log(q / (1 - pi))
    okkl &= math.isclose(g(fs), kl, rel_tol=1e-12)
check("Kelly's fraction f* = p - q/b maximises the growth rate E log(1 + f·r)", okkel)
check("at f* the wealth multiplies by p/pi on a win and q/(1-pi) on a loss, pi = 1/(1+b) the bookmaker's chance; growth = KL(p || pi)", okkl)
ratios = []
for edge in [0.2, 0.1, 0.05, 0.02, 0.01]:
    p = 0.5 + edge / 2; q = 1 - p; fs = p - q; g = lambda f: p * math.log1p(f) + q * math.log1p(-f)
    ratios.append(abs(g(2 * fs)) / g(fs))
check("betting twice Kelly gives about zero growth, exactly zero to second order in the edge", all(ratios[i+1] < ratios[i] for i in range(4)) and ratios[-1] < 0.001,
      " ".join("%.4f" % r for r in ratios))

okshr = True
for p, b in [(0.55, 1.0), (0.6, 1.0), (0.3, 4.0)]:
    q = 1 - p; g = lambda f: p * math.log1p(b * f) + q * math.log1p(-f); fs = p - q / b
    okshr &= g(2.5 * fs) < 0 and g(3 * fs) < g(2.5 * fs)
check("past about twice Kelly the growth turns negative and keeps falling", okshr)

print("== 3. Elo: rating gaps add, odds multiply, edges add like velocities")
def P(d): return 1 / (1 + 10 ** (-d / 400))
oke = True; okv = True; okp = True; okt = True
for _ in range(200):
    d1, d2 = random.uniform(-600, 600), random.uniform(-600, 600)
    p1, p2, p12 = P(d1), P(d2), P(d1 + d2)
    O = lambda p: p / (1 - p)
    oke &= math.isclose(O(p12), O(p1) * O(p2), rel_tol=1e-9)
    oke &= math.isclose(O(p12) - 1, (O(p1) - 1) + (O(p2) - 1) + (O(p1) - 1) * (O(p2) - 1), rel_tol=1e-9, abs_tol=1e-12)
    s1, s2, s12 = 2 * p1 - 1, 2 * p2 - 1, 2 * p12 - 1
    okv &= math.isclose(s12, (s1 + s2) / (1 + s1 * s2), rel_tol=1e-9, abs_tol=1e-12)
    okt &= math.isclose(s1, math.tanh(d1 * math.log(10) / 800), rel_tol=1e-9, abs_tol=1e-12)
    okp &= math.isclose(p12, p1 * p2 / (p1 * p2 + (1 - p1) * (1 - p2)), rel_tol=1e-9)
    okt &= math.isclose(O(p1), (1 + s1) / (1 - s1), rel_tol=1e-9)
check("Elo odds multiply along a chain, so (odds - 1) mult-adds (Bradley-Terry)", oke)
check("the edge s = 2P - 1 composes as (s1 + s2)/(1 + s1·s2): relativity's velocity addition, landscape point (0, -1)", okv)
check("edge = tanh(gap·ln10/800) and odds = (1 + s)/(1 - s): gap is rapidity, odds is the Doppler lift", okt)
check("win chance composes as p1·p2 / (p1·p2 + (1-p1)(1-p2))", okp)
check("a 200-point gap wins 76%; two 200-point gaps (400) win 91%, not 76% + 26%",
      abs(P(200) - 0.7597) < 1e-4 and abs(P(400) - 0.9091) < 1e-4)
s200 = 2 * P(200) - 1
check("a 200-point edge is 0.52; two compose to (0.52+0.52)/(1+0.52^2) = 0.82 (91%), while adding edges would claim 1.04 (102%)",
      abs(s200 - 0.519) < 1e-3 and abs((2 * s200) / (1 + s200**2) - 0.818) < 1e-3 and abs((1 + 2 * s200) / 2 - 1.02) < 0.005)
check("certainty is a wall: no finite gap reaches P = 1; a 2000-point favourite still loses about 1 game in 100,000", P(2000) < 1 and abs((1 - P(2000)) * 1e5 - 1) < 0.01)

print("== 4. The dial game: v(S) = (prod(1 + k x_i) - 1)/k")
def v(k, xs, S): return dial(k, [xs[i] for i in S])
def shapley_brute(k, xs):
    n = len(xs); phi = [F(0)] * n; perms = list(itertools.permutations(range(n)))
    for perm in perms:
        S = []
        for i in perm:
            phi[i] += v(k, xs, S + [i]) - v(k, xs, S); S.append(i)
    return [p / len(perms) for p in phi]
def esym(vals, m):
    return sum((math.prod(c) if c else 1) for c in itertools.combinations(vals, m)) if m > 0 else 1
def shapley_closed(k, xs):
    n = len(xs); out = []
    for i in range(n):
        others = xs[:i] + xs[i+1:]
        out.append(xs[i] * sum(F(k) ** m * esym(others, m) / (m + 1) for m in range(n)))
    return out
okd = all(v(k, xs, [0, 1, 2]) == v(k, xs, [0]) + v(k, xs, [1, 2]) + k * v(k, xs, [0]) * v(k, xs, [1, 2])
          for k in [F(-1), F(1, 2), F(3)] for xs in [[F(1, 5), F(1, 3), F(2, 7)]])
check("disjoint coalitions combine by v(S ∪ T) = v(S) + v(T) + k·v(S)·v(T) (Sugeno's lambda-rule)", okd)
okm = all(v(k, xs, [0, 2, 1]) - v(k, xs, [0, 2]) == xs[1] * (1 + k * v(k, xs, [0, 2])) for k in [F(-1), F(2)] for xs in [[F(1, 5), F(1, 3), F(2, 7)]])
check("marginal value of joining S is x_i·(1 + k·v(S)) = x_i · prod over S of (1 + k x_j)", okm)
oks = True; okeff = True
for n in [2, 3, 4, 5, 6]:
    for k in [F(-1), F(-1, 2), F(1, 3), F(1), F(5, 2)]:
        xs = [F(random.randint(1, 9), 10) for _ in range(n)]
        if any(1 + k * x <= 0 for x in xs): continue
        sb, sc = shapley_brute(k, xs), shapley_closed(k, xs)
        oks &= (sb == sc); okeff &= (sum(sc) == v(k, xs, list(range(n))))
check("Shapley value in closed form: phi_i = x_i · sum_m k^m e_m(others)/(m+1) (exact, n = 2..6, five dial values)", oks)
check("the closed form shares out exactly v(N) (efficiency)", okeff)
oki = True
for _ in range(20):
    k = random.uniform(-0.9, 3); xs = [random.uniform(0.05, 1) for _ in range(5)]
    for i in range(5):
        others = xs[:i] + xs[i+1:]
        val = xs[i] * quad(lambda t: math.prod(1 + k * t * x for x in others), 0, 1)[0]
        oki &= math.isclose(val, float(shapley_closed(F(k), [F(x) for x in xs])[i]), rel_tol=1e-9)
check("Owen's diagonal: phi_i = x_i · integral_0^1 prod_{j≠i}(1 + k t x_j) dt", oki)
okml = True
for k in [F(-1), F(2, 3), F(2)]:
    xs = [F(1, 4), F(2, 5), F(1, 2), F(1, 3)]; t = F(3, 7); n = 4
    ev = sum(t ** len(S) * (1 - t) ** (n - len(S)) * v(k, xs, list(S)) for m in range(n + 1) for S in itertools.combinations(range(n), m))
    okml &= (ev == dial(k, [t * x for x in xs]))
check("everyone shows up with chance t: the expected worth is the dial fold of t·x_i (Owen's multilinear extension)", okml)
okb = True
for k in [F(-1), F(1, 2), F(3)]:
    xs = [F(1, 4), F(2, 5), F(1, 2), F(1, 3)]; n = 4
    for i in range(n):
        rest = [j for j in range(n) if j != i]
        bz = sum(v(k, xs, list(S) + [i]) - v(k, xs, list(S)) for m in range(n) for S in itertools.combinations(rest, m)) / 2 ** (n - 1)
        okb &= (bz == xs[i] * math.prod(1 + k * xs[j] / 2 for j in rest))
check("Banzhaf value: beta_i = x_i · prod_{j≠i}(1 + k x_j / 2), the dial fold at half participation", okb)
okc = True
for k in [F(-1), F(1, 2), F(2)]:
    xs = [F(1, 4), F(2, 5), F(1, 2)]; n = 3
    div = {S: k ** (len(S) - 1) * math.prod(xs[j] for j in S) for m in range(1, n + 1) for S in itertools.combinations(range(n), m)}
    okc &= (sum(div.values()) == v(k, xs, [0, 1, 2]))
    sh = [sum(d / len(S) for S, d in div.items() if i in S) for i in range(n)]
    okc &= (sh == shapley_closed(k, xs))
check("the game's Harsanyi dividends are the cube's corners k^(|S|-1)·prod x_S, and Shapley splits each corner equally", okc)
xs = [F(1, 2), F(3, 10), F(1, 10)]
r_syn = [p / x for p, x in zip(shapley_closed(F(2), xs), xs)]; r_red = [p / x for p, x in zip(shapley_closed(F(-1), xs), xs)]
check("synergy (k > 0): smaller players get the larger bonus ratio phi_i/x_i; redundancy (k < 0): larger players keep more",
      r_syn[0] < r_syn[1] < r_syn[2] and r_red[0] > r_red[1] > r_red[2],
      "k=2: %s  k=-1: %s" % (" ".join("%.3f" % float(r) for r in r_syn), " ".join("%.3f" % float(r) for r in r_red)))
k = 2.0; xf = [0.5, 0.3, 0.1]; vN = float(v(F(2), [F(x) for x in xf], [0, 1, 2]))
logs = [math.log1p(k * x) for x in xf]; lsplit = [vN * l / sum(logs) for l in logs]
sh = [float(s) for s in shapley_closed(F(2), [F(x) for x in xf])]
check("the log-chart split vN·log(1+k x_i)/sum log(1+k x_j) is also efficient but differs from Shapley",
      math.isclose(sum(lsplit), vN) and max(abs(a - b) for a, b in zip(lsplit, sh)) > 1e-3,
      "Shapley %s  log split %s" % (" ".join("%.4f" % s for s in sh), " ".join("%.4f" % s for s in lsplit)))
check("k = -1 is the 'some cause strikes' game, where Shapley = the actuarial UDD split q_i - corners shared equally",
      shapley_closed(F(-1), [F(3, 10), F(1, 2)]) == [F(3, 10) - F(3, 10) * F(1, 2) / 2, F(1, 2) - F(3, 10) * F(1, 2) / 2])

print("== 5. What did not fit")
b3 = lambda p: 3 * p**2 - 2 * p**3
def bestof(n, p): return sum(math.comb(n, j) * p**j * (1 - p)**(n - j) for j in range(n // 2 + 1, n + 1))
check("winning a series is a binomial tail, not a fold: best-of-3 of best-of-3s is not best-of-9",
      abs(b3(b3(0.6)) - bestof(9, 0.6)) > 0.01, "%.4f vs %.4f" % (b3(b3(0.6)), bestof(9, 0.6)))

print("\n%d/%d checks pass" % tuple(OK))
