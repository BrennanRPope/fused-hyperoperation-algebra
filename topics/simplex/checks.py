"""Checks for the simplex pass (Affine Composition Theory).

Run:  python3 checks.py > checks_output.txt
Needs numpy, scipy.
Core law: a o_k b = a + b + k*a*b.  Ruler L_k(x) = log(1+kx)/k.
"""
import math
import itertools
import numpy as np
from scipy import integrate, special, stats

rng = np.random.default_rng(20261126)
results = []


def check(name, ok, detail):
    results.append(ok)
    print(f"[{'PASS' if ok else 'FAIL'}] {name}: {detail}")


def section(t):
    print("\n" + "=" * 78 + "\n" + t + "\n" + "=" * 78)


def mo(a, b, k):
    return a + b + k * a * b


def L(x, k):
    return np.log1p(k * x) / k


def rand_simplex(n, size=None, alpha=1.0):
    return rng.dirichlet([alpha] * n, size=size)


# ---------------------------------------------------------------------------
section("1. THREE COORDINATE SYSTEMS ON ONE SIMPLEX: barycentric, stick, angle")
# ---------------------------------------------------------------------------
# Stick-breaking: p_m = v_m * prod_{j<m}(1-v_j).  Cumulative sum = v_1 o_{-1} ... o_{-1} v_m.
worst = 0.0
for _ in range(2000):
    n = rng.integers(2, 9)
    v = rng.random(n - 1)
    p = np.empty(n)
    rem = 1.0
    for m in range(n - 1):
        p[m] = v[m] * rem
        rem *= 1 - v[m]
    p[-1] = rem
    P = 0.0
    for m in range(n - 1):
        P = mo(P, v[m], -1)
        worst = max(worst, abs(P - p[: m + 1].sum()))
check("S1 stick-breaking: running total of pieces = v1 o_{-1} v2 o_{-1} ... (the core law at k=-1)",
      worst < 1e-12, f"worst error {worst:.1e} over 2000 random sticks")

# The map is a bijection from the open cube (0,1)^{n-1} to the open simplex; inverse v_m = p_m / (1 - P_{m-1}).
worst = 0.0
for _ in range(1000):
    n = rng.integers(2, 9)
    p = rand_simplex(n)
    v = p[:-1] / (1 - np.concatenate([[0], np.cumsum(p)[:-2]]))
    rem, q = 1.0, np.empty(n)
    for m in range(n - 1):
        q[m] = v[m] * rem
        rem *= 1 - v[m]
    q[-1] = rem
    worst = max(worst, np.abs(q - p).max(), float(np.any((v <= 0) | (v >= 1))))
check("S2 stick map cube -> simplex is one-to-one (inverse v_m = p_m/(1-P_{m-1}) lands in (0,1))",
      worst < 1e-12, f"worst round-trip error {worst:.1e}")

# Angles: v_m = cos^2(theta_m) gives hyperspherical coordinates of sqrt(p) on the positive orthant.
worst = 0.0
for _ in range(1000):
    n = rng.integers(2, 9)
    th = rng.random(n - 1) * np.pi / 2
    x = np.empty(n)
    s = 1.0
    for m in range(n - 1):
        x[m] = s * np.cos(th[m])
        s *= np.sin(th[m])
    x[-1] = s
    v = np.cos(th) ** 2
    rem, p = 1.0, np.empty(n)
    for m in range(n - 1):
        p[m] = v[m] * rem
        rem *= 1 - v[m]
    p[-1] = rem
    worst = max(worst, np.abs(x ** 2 - p).max())
check("S3 stick fraction v = cos^2(angle): stick-breaking is hyperspherical coordinates of sqrt(p)",
      worst < 1e-12, f"worst error {worst:.1e}")

# Uniform on the simplex <-> v_m ~ Beta(1, n-m) independent (cube with a tilt).
n = 4
N = 200000
v = np.column_stack([rng.beta(1, n - 1 - m, N) for m in range(n - 1)])
p = np.empty((N, n))
rem = np.ones(N)
for m in range(n - 1):
    p[:, m] = v[:, m] * rem
    rem = rem * (1 - v[:, m])
p[:, -1] = rem
ks = max(stats.kstest(p[:, i], "beta", args=(1, n - 1)).statistic for i in range(n))
pair = np.mean(p[:, 0] * p[:, 1])  # Dirichlet(1,1,1,1): E[p1 p2] = 1/(n(n+1)) = 1/20
check("S4 uniform simplex = independent Beta(1, n-m) stick fractions (n=4)",
      ks < 0.006 and abs(pair - 0.05) < 0.001, f"largest KS {ks:.4f}; E[p1 p2] = {pair:.4f} (exact 0.0500)")

# ---------------------------------------------------------------------------
section("2. THE SIMPLEX AND THE SPHERE: square roots (Fisher) and forgotten phases (Born)")
# ---------------------------------------------------------------------------
# Fisher length of the 1-simplex: int_0^1 dp / sqrt(p(1-p)) = pi.
val, err = integrate.quad(lambda t: 1 / np.sqrt(t * (1 - t)), 0, 1)
check("F1 Fisher length of the two-outcome simplex is pi", abs(val - np.pi) < 1e-8, f"{val:.10f}")

# Fisher distance = 2 arccos(sum sqrt(p q)) (Bhattacharyya angle); verify by integrating the metric on the
# geodesic, which is the great circle through sqrt(p), sqrt(q).
worst = 0.0
for _ in range(300):
    n = rng.integers(2, 7)
    p, q = rand_simplex(n), rand_simplex(n)
    a, b = np.sqrt(p), np.sqrt(q)
    ang = np.arccos(np.clip(a @ b, -1, 1))
    w = (b - (a @ b) * a)
    w /= np.linalg.norm(w)
    ts = np.linspace(0, ang, 4001)
    path = np.outer(np.cos(ts), a) + np.outer(np.sin(ts), w)
    P = path ** 2
    dP = np.diff(P, axis=0)
    mid = 0.5 * (P[1:] + P[:-1])
    length = np.sum(np.sqrt(np.sum(dP ** 2 / mid, axis=1)))
    worst = max(worst, abs(length - 2 * ang))
check("F2 Fisher distance on the simplex = 2 x the angle between sqrt(p) and sqrt(q) (sphere of radius 2)",
      worst < 1e-5, f"worst gap {worst:.1e} over 300 pairs")

# Fisher volume of Delta_{n-1} = |S^{n-1}| / 2.
for n in (2, 3, 4):
    sphere = 2 * np.pi ** (n / 2) / special.gamma(n / 2)
    # Fisher volume element: prod p_i^{-1/2} dp over the simplex = Dirichlet(1/2) normaliser
    dir_norm = np.prod([special.gamma(0.5)] * n) / special.gamma(n / 2)
    check(f"F3 n={n}: Fisher volume of the simplex = half the unit sphere S^{n-1}",
          abs(dir_norm - sphere / 2) < 1e-12, f"{dir_norm:.6f} vs {sphere/2:.6f}")

# Archimedes / Born: uniform point on the complex sphere in C^n -> |z|^2 uniform on the simplex.
for n in (2, 3):
    z = rng.normal(size=(200000, n)) + 1j * rng.normal(size=(200000, n))
    z /= np.linalg.norm(z, axis=1, keepdims=True)
    P = np.abs(z) ** 2
    ks = stats.kstest(P[:, 0], "beta", args=(1, n - 1)).statistic
    check(f"F4 n={n}: random complex amplitudes give probabilities uniform on the simplex (Archimedes' hat-box)",
          ks < 0.005, f"KS distance to uniform marginal {ks:.4f}")

# Real amplitudes instead: Jeffreys prior Dirichlet(1/2): for n=2 the arcsine law, density 1/(pi sqrt(p(1-p))).
x = rng.normal(size=(200000, 2))
x /= np.linalg.norm(x, axis=1, keepdims=True)
ks = stats.kstest(x[:, 0] ** 2, "arcsine").statistic
check("F5 real amplitudes (no phase) give the arcsine law 1/(pi sqrt(p(1-p))) instead of uniform",
      ks < 0.005, f"KS {ks:.4f}")

# Levy's arcsine law: the fraction of time a long random walk spends positive is arcsine distributed, i.e. its
# simplex angle theta = arcsin(sqrt(frac)) is uniform on [0, pi/2].
steps = 2000
walks = 20000
s = np.cumsum(rng.choice([-1, 1], size=(walks, steps)).astype(np.int16), axis=1)
frac = np.mean((s > 0) | ((s == 0) & (np.roll(s, 1, axis=1) > 0)), axis=1)
theta = np.arcsin(np.sqrt(frac))
ks = stats.kstest(theta, "uniform", args=(0, np.pi / 2)).statistic
check("F6 Levy: time-positive fraction of a random walk has uniform simplex angle (arcsine law)",
      ks < 0.02, f"KS to uniform angle {ks:.4f} ({walks} walks of {steps} steps)")

# Wootters: angle between real amplitude vectors = half the Fisher distance (statistical distance = Hilbert angle).
worst = 0.0
for _ in range(500):
    n = rng.integers(2, 8)
    p, q = rand_simplex(n), rand_simplex(n)
    fs = np.arccos(np.clip(np.sqrt(p) @ np.sqrt(q), -1, 1))
    worst = max(worst, abs(2 * fs - 2 * np.arccos(np.sum(np.sqrt(p * q)))))
check("F7 the quantum angle between sqrt(p) and sqrt(q) is half the Fisher distance (Wootters 1981)",
      worst < 1e-12, f"worst {worst:.1e}")

# ---------------------------------------------------------------------------
section("3. THE SIMPLEX'S OWN RADIUS: log-ratios, Bayes, Gibbs")
# ---------------------------------------------------------------------------
# Aitchison perturbation p (+) q = normalise(p_i q_i).  Ratios r_i = p_i/p_n multiply, so r - 1 composes by o_1.
worst = 0.0
for _ in range(2000):
    n = rng.integers(2, 8)
    p, q = rand_simplex(n), rand_simplex(n)
    pq = p * q / np.sum(p * q)
    rp, rq, rpq = p / p[-1] - 1, q / q[-1] - 1, pq / pq[-1] - 1
    worst = max(worst, np.abs(rpq - mo(rp, rq, 1)).max() / (1 + np.abs(rpq).max()))
check("R1 Bayes / perturbation on the simplex is the core law at k=+1 on (odds - 1)",
      worst < 1e-10, f"worst relative error {worst:.1e}")

worst = 0.0
for _ in range(1000):
    n = rng.integers(2, 8)
    p = rand_simplex(n)
    r = p / p[-1] - 1
    worst = max(worst, np.abs(L(r, 1) - (np.log(p) - np.log(p[-1]))).max())
check("R2 the simplex's log-ratio coordinate (Aitchison) is our ruler L_1 applied to odds - 1",
      worst < 1e-12, f"worst {worst:.1e}")

# Bayes: prior x likelihood normalised = perturbation, and it is a straight line in log-ratio coordinates.
worst = 0.0
for _ in range(500):
    n = 3
    prior = rand_simplex(n)
    lik = rng.random(n)
    post = prior * lik / np.sum(prior * lik)
    # repeated identical evidence: log-ratios move by equal steps
    u = [np.log(prior) - np.log(prior).mean()]
    cur = prior
    for _ in range(5):
        cur = cur * lik / np.sum(cur * lik)
        u.append(np.log(cur) - np.log(cur).mean())
    du = np.diff(np.array(u), axis=0)
    worst = max(worst, np.abs(du - du[0]).max(), np.abs(post - (prior * lik) / (prior @ lik)).max())
check("R3 repeated evidence moves the belief in equal straight steps in log-ratio coordinates",
      worst < 1e-12, f"worst deviation {worst:.1e}")

# Gibbs / softmax: temperature slides the barycentre (hot) to a corner (cold).
E = np.array([0.0, 0.7, 1.9])
ps = [np.exp(-b * E) / np.exp(-b * E).sum() for b in (0.0, 1.0, 50.0)]
check("R4 Gibbs weights: beta=0 at the centre, beta large at the lowest-energy corner",
      np.allclose(ps[0], 1 / 3) and ps[2][0] > 1 - 1e-12,
      f"beta=0 {np.round(ps[0],3)}, beta=1 {np.round(ps[1],3)}, beta=50 {np.round(ps[2],6)}")

# Temperature dial turns the sphere into the simplex: B_beta = {x in plane sum=0 : (1/beta) log mean exp(-beta x_i) <= 1},
# rescaled. Measure roundness = (largest radius)/(smallest radius) of the boundary.
def body_radius(beta, omega):
    # find r with F(r omega) = 1 where F(x) = (1/beta) log mean exp(-beta x)
    lo, hi = 0.0, 1.0
    f = lambda r: (special.logsumexp(-beta * r * omega) - np.log(len(omega))) / beta - 1.0
    while f(hi) < 0:
        hi *= 2
    for _ in range(200):
        mid = 0.5 * (lo + hi)
        if f(mid) < 0:
            lo = mid
        else:
            hi = mid
    return 0.5 * (lo + hi)


e1 = np.array([1, -1, 0]) / np.sqrt(2)
e2 = np.array([1, 1, -2]) / np.sqrt(6)
angles = np.linspace(0, 2 * np.pi, 721)[:-1]
ratios = {}
for beta in (1e-4, 1.0, 10.0, 1000.0):
    rr = [body_radius(beta, np.cos(a) * e1 + np.sin(a) * e2) for a in angles]
    ratios[beta] = max(rr) / min(rr)
check("R5 temperature dial: hot free-energy ball is round, cold one is the simplex (radius ratio 1 -> 2)",
      abs(ratios[1e-4] - 1) < 0.01 and abs(ratios[1000.0] - 2) < 0.01 and ratios[1.0] < ratios[10.0],
      "max/min radius: " + ", ".join(f"beta={b}: {r:.4f}" for b, r in ratios.items()) + " (triangle exact 2)")

# Information: the Hessian of negative entropy on the simplex is the Fisher metric.
worst = 0.0
for _ in range(200):
    n = 4
    p = rand_simplex(n, alpha=4.0)
    # coordinates p_1..p_{n-1}, p_n = 1 - sum
    G = np.diag(1 / p[:-1]) + 1 / p[-1]
    h = 1e-4
    def negH(y):
        q = np.append(y, 1 - y.sum())
        return np.sum(q * np.log(q))
    y0 = p[:-1]
    Hm = np.zeros((n - 1, n - 1))
    for i in range(n - 1):
        for j in range(n - 1):
            ei, ej = np.eye(n - 1)[i] * h, np.eye(n - 1)[j] * h
            Hm[i, j] = (negH(y0 + ei + ej) - negH(y0 + ei - ej) - negH(y0 - ei + ej) + negH(y0 - ei - ej)) / (4 * h * h)
    worst = max(worst, np.abs(Hm - G).max() / np.abs(G).max())
check("R6 the curvature (Hessian) of entropy on the simplex is the Fisher metric", worst < 1e-3,
      f"worst relative gap {worst:.1e}")

# Replicator equation = Fisher (Shahshahani) gradient of mean fitness, for constant fitnesses.
worst = 0.0
for _ in range(500):
    n = 5
    p = rand_simplex(n)
    f = rng.normal(size=n)
    rep = p * (f - p @ f)
    # Shahshahani gradient of fbar = p.f restricted to sum=0: g = P (diag p) f with projection keeping sum 0
    g = p * f - p * (p @ f)
    worst = max(worst, np.abs(rep - g).max())
check("R7 replicator equation is the Fisher-metric gradient of mean fitness (Shahshahani)", worst < 1e-14,
      f"worst {worst:.1e}")

# ---------------------------------------------------------------------------
section("4. e LIVES ON THE SIMPLEX: corner volumes, orderings, the e-page layers")
# ---------------------------------------------------------------------------
for n in (2, 3, 4, 5):
    u = rng.random((400000, n))
    vol = np.mean(u.sum(axis=1) <= 1)
    check(f"E1 n={n}: corner simplex of the unit cube has volume 1/n!", abs(vol - 1 / math.factorial(n)) < 4 * np.sqrt(vol * (1 - vol) / 400000) + 1e-4,
          f"{vol:.5f} vs {1/math.factorial(n):.5f}")

# Cube = n! order simplices (Kuhn): each ordering of coordinates is one simplex of volume 1/n!.
n = 4
u = rng.random((240000, n))
orders = np.argsort(u, axis=1)
codes = (orders * (n ** np.arange(n))).sum(axis=1)
counts = np.unique(codes, return_counts=True)[1]
check("E2 the 4-cube splits into 24 equal simplices, one per ordering of the coordinates (Kuhn)",
      len(counts) == 24 and np.abs(counts / len(u) - 1 / 24).max() < 0.002,
      f"{len(counts)} pieces, shares {counts.min()/len(u):.4f} to {counts.max()/len(u):.4f} (exact {1/24:.4f})")

# Sum of volumes over all dimensions.
S_simp = sum(1 / math.factorial(n) for n in range(60))
S_cross = sum(2 ** n / math.factorial(n) for n in range(80))
S_ball = sum(np.pi ** (n / 2) / special.gamma(n / 2 + 1) for n in range(200))
ball_closed = np.exp(np.pi) * (1 + special.erf(np.sqrt(np.pi)))
check("E3 volumes summed over every dimension: corner simplex -> e, cross-polytope -> e^2, ball -> e^pi(1+erf sqrt pi)",
      abs(S_simp - np.e) < 1e-14 and abs(S_cross - np.e ** 2) < 1e-12 and abs(S_ball - ball_closed) < 1e-10,
      f"simplex {S_simp:.12f} (e {np.e:.12f}); cross {S_cross:.10f} (e^2 {np.e**2:.10f}); ball {S_ball:.10f} vs {ball_closed:.10f}; cube diverges")

# Expected number of uniforms to pass 1 is e (each P(N>n) is a simplex volume).
u = rng.random((300000, 12))
N_ = np.argmax(np.cumsum(u, axis=1) > 1, axis=1) + 1
check("E4 average number of random numbers in [0,1] needed to pass 1 is e (sum of simplex volumes)",
      abs(N_.mean() - np.e) < 0.01, f"{N_.mean():.4f} vs e = {np.e:.4f}")

# The e page: layer k of 1 + (1/n o 1/n o ... o 1/n) (k=1 core law) is C(n,k)/n^k -> 1/k! = volume of the k-simplex.
for kk in (1, 2, 3, 5):
    vals = [math.comb(nn, kk) / nn ** kk for nn in (10, 100, 10000)]
    check(f"E5 e-page layer {kk} -> 1/{kk}! = corner {kk}-simplex volume",
          abs(vals[-1] - 1 / math.factorial(kk)) < 1e-3 * kk, f"n=10,100,10^4: {vals[0]:.5f}, {vals[1]:.5f}, {vals[2]:.6f}; 1/{kk}! = {1/math.factorial(kk):.6f}")

# Volume of t times the corner simplex is t^n/n!, so the generating function over dimensions is e^t (the simplex is
# the shape of the exponential).
t = 1.7
gen = sum(t ** n / math.factorial(n) for n in range(60))
check("E6 sum over dimensions of vol(t * corner simplex) = e^t", abs(gen - np.exp(t)) < 1e-12, f"{gen:.12f} vs {np.exp(t):.12f}")

# ---------------------------------------------------------------------------
section("5. HILBERT'S DISTANCE ON THREE SHAPES: ball = relativity, simplex = Bayes, square")
# ---------------------------------------------------------------------------
def hilbert_dist(p, q, inside):
    """Hilbert distance (1/2 log cross-ratio) in a convex domain given by an inside(x)->bool test."""
    d = q - p
    nd = np.linalg.norm(d)
    if nd < 1e-15:
        return 0.0
    def edge(sign):
        lo, hi = 0.0, 1.0
        while inside(p + sign * hi * d):
            hi *= 2
        for _ in range(200):
            mid = 0.5 * (lo + hi)
            if inside(p + sign * mid * d):
                lo = mid
            else:
                hi = mid
        return 0.5 * (lo + hi)
    tp, tm = edge(1), edge(-1)  # boundary at p + tp d and p - tm d
    # cross ratio with a = p - tm d, b = p + tp d, points at 0 and 1
    return 0.5 * np.log((tp * (1 + tm)) / ((tp - 1) * tm))


# 1D: segment (-1,1): Hilbert distance from 0 to v is artanh v (rapidity); composing translations = Einstein addition.
worst = 0.0
for _ in range(500):
    a, b = rng.uniform(-0.95, 0.95, 2)
    ein = (a + b) / (1 + a * b)  # plane law k=0, d=-1
    worst = max(worst, abs(np.arctanh(ein) - np.arctanh(a) - np.arctanh(b)))
check("H1 segment: Hilbert distance = rapidity, and Einstein addition (our plane law k=0, d=-1) adds it",
      worst < 1e-12, f"worst {worst:.1e}")

# 1D simplex [0,1] is the same segment: p = (1+v)/2, log odds = 2 artanh v.  Bayes on a coin = velocity addition.
worst = 0.0
for _ in range(500):
    v, w = rng.uniform(-0.95, 0.95, 2)
    p, q = (1 + v) / 2, (1 + w) / 2
    O = (p / (1 - p)) * (q / (1 - q))      # Bayes: odds multiply (q as likelihood ratio in odds form)
    vw = (v + w) / (1 + v * w)
    worst = max(worst, abs(O - (1 + vw) / (1 - vw)))
check("H2 a coin's odds update and Einstein velocity addition are one law: odds = (1+v)/(1-v)",
      worst < 1e-9, f"worst {worst:.1e}")

# Simplex: Hilbert distance = 1/2 log max_{ij} p_i q_j / (p_j q_i), invariant under perturbation.
inside_simplex = lambda x: np.all(x > 0)
worst, worst_inv = 0.0, 0.0
for _ in range(300):
    n = rng.integers(2, 6)
    p, q, s = rand_simplex(n), rand_simplex(n), rand_simplex(n)
    dh = hilbert_dist(p, q, inside_simplex)
    lr = np.log(p) - np.log(q)
    formula = 0.5 * (lr.max() - lr.min())
    ps, qs = p * s / (p @ s), q * s / (q @ s)
    worst = max(worst, abs(dh - formula))
    worst_inv = max(worst_inv, abs(hilbert_dist(ps, qs, inside_simplex) - dh))
check("H3 simplex: Hilbert distance = half the spread (max - min) of the log-ratios", worst < 1e-8, f"worst {worst:.1e}")
check("H4 simplex: Bayes updates (perturbations) are distance-preserving shifts", worst_inv < 1e-8, f"worst {worst_inv:.1e}")

# Normed space: distance along 'power' paths scales linearly: d(c, c (+) t q) = t d(c, q) for the centre c.
worst = 0.0
for _ in range(200):
    n = 3
    q = rand_simplex(n)
    c = np.ones(n) / n
    d1 = hilbert_dist(c, q, inside_simplex)
    for t in (0.5, 2.0, 3.0):
        qt = q ** t / np.sum(q ** t)
        worst = max(worst, abs(hilbert_dist(c, qt, inside_simplex) - t * d1))
check("H5 simplex Hilbert geometry is flat (a normed space): tempering a belief by t scales its distance by t",
      worst < 1e-5, f"worst {worst:.1e} (bisection precision near the edge)")

# Ball: Hilbert distance = hyperbolic distance (Klein model): cosh d = (1 - p.q)/sqrt((1-|p|^2)(1-|q|^2)).
inside_ball = lambda x: x @ x < 1
worst = 0.0
for _ in range(300):
    p = rng.normal(size=3); p *= rng.random() ** (1/3) * 0.95 / np.linalg.norm(p)
    q = rng.normal(size=3); q *= rng.random() ** (1/3) * 0.95 / np.linalg.norm(q)
    dh = hilbert_dist(p, q, inside_ball)
    hyp = np.arccosh((1 - p @ q) / np.sqrt((1 - p @ p) * (1 - q @ q)))
    worst = max(worst, abs(dh - hyp))
check("H6 ball: Hilbert distance = hyperbolic distance (Klein model) = relativistic velocity space",
      worst < 1e-8, f"worst {worst:.1e}")

# Einstein addition in 3D does not commute: u (+) v and v (+) u differ by the Thomas-Wigner rotation.
def einstein(u, v):
    gu = 1 / np.sqrt(1 - u @ u)
    return (u + v / gu + (gu / (1 + gu)) * (u @ v) * u) / (1 + u @ v)


worst = 0.0
for _ in range(200):
    a, b = rng.uniform(0.05, 0.95, 2)
    u, v = np.array([a, 0, 0]), np.array([0, b, 0])
    w1, w2 = einstein(u, v), einstein(v, u)
    ang = np.arccos(np.clip(w1 @ w2 / (np.linalg.norm(w1) * np.linalg.norm(w2)), -1, 1))
    gu, gv = 1 / np.sqrt(1 - a * a), 1 / np.sqrt(1 - b * b)
    wigner = np.arccos((gu + gv) / (1 + gu * gv))
    worst = max(worst, abs(ang - wigner), abs(np.linalg.norm(w1) - np.linalg.norm(w2)))
check("H7 ball: velocity addition for perpendicular velocities misses commuting by the Wigner angle arccos((gu+gv)/(1+gu gv))",
      worst < 1e-9, f"worst {worst:.1e}")

# That angle equals the area (defect) of the hyperbolic triangle with corners 0, u, v (Klein coordinates)?
# Test: triangle 0, u, u(+)v, compute its hyperbolic angles from the hyperboloid model.
def to_hyperboloid(x):
    g = 1 / np.sqrt(1 - x @ x)
    return g * np.append(1.0, x)


def mink(a, b):
    return -a[0] * b[0] + a[1:] @ b[1:]


def angle_at(A, B, C):
    # angle at A between geodesics to B and C, hyperboloid model
    tB = B + mink(A, B) * A
    tC = C + mink(A, C) * A
    return np.arccos(np.clip(mink(tB, tC) / np.sqrt(mink(tB, tB) * mink(tC, tC)), -1, 1))


worst = 0.0
for _ in range(200):
    a, b = rng.uniform(0.05, 0.95, 2)
    u, v = np.array([a, 0, 0]), np.array([0, b, 0])
    w = einstein(u, v)
    A, B, C = to_hyperboloid(np.zeros(3)), to_hyperboloid(u), to_hyperboloid(w)
    defect = np.pi - angle_at(A, B, C) - angle_at(B, A, C) - angle_at(C, A, B)
    gu, gv = 1 / np.sqrt(1 - a * a), 1 / np.sqrt(1 - b * b)
    wigner = np.arccos((gu + gv) / (1 + gu * gv))
    worst = max(worst, abs(defect - wigner))
check("H8 ball: that rotation angle is exactly the area of the velocity triangle (swap rule turn = d x area, d=-1)",
      worst < 1e-9, f"worst {worst:.1e}")

# Simplex: perturbation commutes, so no rotation is left over.
worst = 0.0
for _ in range(1000):
    n = rng.integers(2, 7)
    p, q = rand_simplex(n), rand_simplex(n)
    worst = max(worst, np.abs(p * q / (p @ q) - q * p / (q @ p)).max())
check("H9 simplex: Bayes updates commute exactly, so the leftover rotation is zero (no precession)",
      worst == 0.0, f"worst {worst:.1e}")

# Square: Hilbert distance is not a normed space (Foertsch-Karlsson: only simplices are): tempering-type scaling
# along straight rays from the centre fails to be linear in the rapidity-like coordinate, and distance is not
# the max of the two axis rapidities.
inside_sq = lambda x: np.all(np.abs(x) < 1)
gaps = []
for _ in range(400):
    p, q = rng.uniform(-0.9, 0.9, 2), rng.uniform(-0.9, 0.9, 2)
    dh = hilbert_dist(p, q, inside_sq)
    axis = max(abs(np.arctanh(p[0]) - np.arctanh(q[0])), abs(np.arctanh(p[1]) - np.arctanh(q[1])))
    gaps.append(dh - axis)
gaps = np.array(gaps)
check("H10 square: Hilbert distance is NOT the larger of the two axis rapidities (no product rule)",
      np.abs(gaps).max() > 0.05, f"gap from max-of-axes ranges {gaps.min():.3f} to {gaps.max():.3f}")

# Check the square is not a normed space either: ball shapes change with the centre (compare radii in fixed
# directions in the log-odds chart, where the simplex's balls are all the same hexagon).
def ball_shape(center, inside, chart, r=0.4, dirs=48):
    out = []
    for a in np.linspace(0, 2 * np.pi, dirs, endpoint=False):
        dvec = np.array([np.cos(a), np.sin(a)])
        lo, hi = 0.0, 1.0
        # move in Euclidean coords until Hilbert distance r
        while inside(center + hi * dvec) and hilbert_dist(center, center + hi * dvec, inside) < r:
            hi *= 1.5
        for _ in range(80):
            mid = 0.5 * (lo + hi)
            if inside(center + mid * dvec) and hilbert_dist(center, center + mid * dvec, inside) < r:
                lo = mid
            else:
                hi = mid
        out.append(chart(center + lo * dvec) - chart(center))
    return np.array(out)


chart_sq = lambda x: np.arctanh(x)
s0 = ball_shape(np.zeros(2), inside_sq, chart_sq)
s1 = ball_shape(np.array([0.6, 0.0]), inside_sq, chart_sq)
spread = np.abs(np.linalg.norm(s0, axis=1) - np.linalg.norm(s1, axis=1)).max()
check("H11 square: in its own rapidity chart, Hilbert balls change shape as the centre moves (not flat)",
      spread > 0.02, f"largest radius change {spread:.3f}")

# ---------------------------------------------------------------------------
print("\n" + "=" * 78)
print(f"TOTAL: {sum(results)}/{len(results)} checks pass")
