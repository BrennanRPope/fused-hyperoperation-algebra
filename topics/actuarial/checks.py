"""The circle operation (a (+) b = a + b + ab) as actuarial tools, and the infinity test.

Run:  python3 checks.py > checks_output.txt
Needs numpy, scipy, rdata. Data files are in ./data (see data/SOURCES.md).
"""
import itertools, math, os, warnings
import numpy as np
from math import comb, log, exp
warnings.filterwarnings("ignore")
rng = np.random.default_rng(20261125)
HERE = os.path.dirname(os.path.abspath(__file__))

def ma(a, b):            # circle operation
    return a + b + a * b
def ma_all(xs):
    out = 0.0
    for x in xs:
        out = ma(out, x)
    return out
def ma_inv(a):           # circle inverse: a (+) inv = 0
    return -a / (1 + a)
def ma_pow(a, n):        # n-th circle power (any real n)
    return (1 + a) ** n - 1
def cube_sum(ws, sign=1):
    """Sum over nonempty corners S of the cube of prod_{k in S} w_k (times sign^(|S|+1))."""
    tot = 0.0
    n = len(ws)
    for r in range(1, n + 1):
        for S in itertools.combinations(range(n), r):
            tot += sign ** (r + 1) * np.prod([ws[k] for k in S])
    return tot

def load(name, key):
    import rdata
    return rdata.read_rda(os.path.join(HERE, "data", name))[key]

def ok(flag):
    return "True" if flag else "FALSE"

print("== 1. Interest is the circle operation ==")
i, j = 0.03, 0.05
print(f" 3% then 5%: (+) gives {ma(i,j):.6f}; (1.03)(1.05)-1 = {1.03*1.05-1:.6f}")
print(f" force of interest delta = log(1+i) at 5%: {log(1.05):.6f}; logs add: {ok(abs(log(1+ma(i,j))-(log(1+i)+log(1+j)))<1e-15)}")
d = 0.05 / 1.05
print(f" discount rate d = i/(1+i) = {d:.6f}; minus the circle inverse of i = {-ma_inv(0.05):.6f}: {ok(abs(d+ma_inv(0.05))<1e-15)}")
i12 = 12 * ma_pow(0.05, 1 / 12)
print(f" nominal monthly rate i(12) = 12 x (12th circle root of 5%) = {i12:.6f}")
n = 10
s_closed = ma_pow(0.05, n) / 0.05
s_cube = sum(comb(n, k) * 0.05 ** (k - 1) for k in range(1, n + 1))
s_cube_corners = cube_sum([0.05] * n) / 0.05
print(f" 10-year accumulation s_10 at 5%: closed {s_closed:.6f}, layer sum {s_cube:.6f}, corner sum {s_cube_corners:.6f}")
print(f"  layer 1 (simple interest) = {n}; layers 2..10 (interest on interest) = {s_cube-n:.6f}")
layers = [comb(n, k) * 0.05 ** k for k in range(1, n + 1)]
print("  share of (1.05^10 - 1) in each cube layer:", [round(x / ma_pow(0.05, n), 4) for x in layers[:5]], "...")

print("\n== 2. Survival is the circle operation of negated death rates ==")
ilt = load("soa08Act.rda", "soa08Act")
lx = np.array(ilt.lx, float); ages = np.array(ilt.x)
qx = 1 - lx[1:] / lx[:-1]
q = qx[45:48]
three = 1 - lx[48] / lx[45]
print(f" SOA Illustrative Life Table, ages 45-47: q = {np.round(q,6)}")
print(f" 3q45 from the table {three:.8f}; -((-q45)(+)(-q46)(+)(-q47)) = {-ma_all(-q):.8f}; cube (inclusion-exclusion) = {cube_sum(q,-1):.8f}")
mu = -np.log(1 - q)
print(f" force of mortality = -log(1 + (-q)) adds: {ok(abs(-np.log(1-three) - mu.sum())<1e-12)}")

print("\n== 3. Multiple decrements: a timing assumption is a path through the cube ==")
q45 = qx[45]
qp = np.array([q45, 0.060, 0.004])      # absolute rates: death (table), lapse, disability (illustrative)
names = ["death", "lapse", "disability"]
tot = -ma_all(-qp)
print(f" absolute rates q' = {qp}  (death from the SOA table at 45; lapse 6% and disability 0.4% are illustrative)")
print(f" total q(tau) = 1 - prod(1 - q') = {tot:.8f}; plain sum overcounts by {qp.sum()-tot:.8f} (the cube's faces)")
# UDD in each single-decrement table: q_j * integral_0^1 prod_{k!=j}(1 - t q_k) dt
from scipy.integrate import quad
def udd(qp, jj):
    return quad(lambda t: qp[jj] * np.prod([1 - t * qp[k] for k in range(len(qp)) if k != jj]), 0, 1)[0]
udd_rates = np.array([udd(qp, jj) for jj in range(3)])
# corner-splitting: each corner's term shared equally by its members
corner = np.zeros(3)
for r in range(1, 4):
    for S in itertools.combinations(range(3), r):
        term = (-1) ** (r + 1) * np.prod(qp[list(S)])
        for k in S:
            corner[k] += term / r
# Shapley value of the game v(S) = P(some decrement in S happens) = 1 - prod_{S}(1-q')
def v(S):
    return 1 - np.prod([1 - qp[k] for k in S]) if S else 0.0
shap = np.zeros(3)
perms = list(itertools.permutations(range(3)))
for p in perms:
    seen = []
    for k in p:
        shap[k] += (v(seen + [k]) - v(seen)) / len(perms)
        seen.append(k)
print(" UDD dependent rates      ", np.round(udd_rates, 9), "sum", round(udd_rates.sum(), 9))
print(" equal split of each corner", np.round(corner, 9))
print(" Shapley value of the game ", np.round(shap, 9))
print(f"  UDD = equal corner split = Shapley: {ok(np.allclose(udd_rates,corner,atol=1e-13) and np.allclose(udd_rates,shap,atol=1e-13))}")
# constant force: proportional to circle-operation logs
mus = -np.log(1 - qp)
cf = mus / mus.sum() * tot
# Aumann-Shapley along the diagonal in log coordinates: v = 1 - exp(-sum mu_k x_k)
asl = np.array([quad(lambda t: mus[jj] * np.exp(-t * mus.sum()), 0, 1)[0] for jj in range(3)])
print(" constant force rates     ", np.round(cf, 9), "sum", round(cf.sum(), 9))
print(" diagonal path in log coords", np.round(asl, 9))
print(f"  constant force = straight path in circle-operation log coordinates: {ok(np.allclose(cf,asl,atol=1e-13))}")
# a third path: every death happens in the first half of the year, everything else uniform
def path_rates(Q):    # Q(t) -> array of cumulative absolute probabilities, Q(1) = qp
    ts = np.linspace(0, 1, 200001)
    P = np.array([Q(t) for t in ts])
    dP = np.diff(P, axis=0)
    mid = (P[1:] + P[:-1]) / 2
    out = []
    for jj in range(3):
        others = np.prod(np.delete(1 - mid, jj, axis=1), axis=1)
        out.append(np.sum(others * dP[:, jj]))
    return np.array(out)
early = path_rates(lambda t: np.array([qp[0] * min(2 * t, 1), qp[1] * t, qp[2] * t]))
print(" deaths-early path        ", np.round(early, 9), "sum", round(early.sum(), 9))
print(f"  totals agree on all three paths (path-independent): {ok(abs(early.sum()-tot)<1e-9 and abs(cf.sum()-tot)<1e-12 and abs(udd_rates.sum()-tot)<1e-12)}")
# Monte Carlo check of the UDD and constant force splits
N = 4_000_000
def simulate(kind):
    times = np.full((N, 3), np.inf)
    for k in range(3):
        if kind == "udd":
            hit = rng.random(N) < qp[k]
            times[hit, k] = rng.random(hit.sum())
        else:
            tt = rng.exponential(1 / mus[k], N)
            times[tt < 1, k] = tt[tt < 1]
    first = np.argmin(times, axis=1)
    happened = np.isfinite(times.min(axis=1))
    return np.array([np.mean(happened & (first == k)) for k in range(3)])
mc_u, mc_c = simulate("udd"), simulate("cf")
print(" simulation (4 million lives), UDD:", np.round(mc_u, 5), " constant force:", np.round(mc_c, 5))
print(f"  within 3 standard errors: {ok(np.all(np.abs(mc_u-udd_rates) < 3*np.sqrt(udd_rates/N)) and np.all(np.abs(mc_c-cf) < 3*np.sqrt(cf/N)))}")

print("\n== 4. Independent risks: (M - 1) combines by the circle operation, cumulants are the circle-operation log ==")
def M_exp(s, b=1.0): return b / (b - s)
def M_gamma(s, k=2.0, b=2.0): return (1 - s / b) ** (-k)
def M_cpois(s, lam=0.5): return np.exp(lam * (M_exp(s, 3.0) - 1))
s = 0.3
a = [M_exp(s) - 1, M_gamma(s) - 1, M_cpois(s) - 1]
prod = M_exp(s) * M_gamma(s) * M_cpois(s)
print(f" at s = {s}: M_sum - 1 = {prod-1:.10f}; (+) of (M_k - 1) = {ma_all(a):.10f}; cube corners = {cube_sum(a):.10f}")
# Monte Carlo of the sum
NN = 2_000_000
Nn = rng.poisson(0.5, NN)
Xs = rng.exponential(1.0, NN) + rng.gamma(2.0, 1 / 2.0, NN) + np.where(Nn > 0, rng.gamma(np.maximum(Nn, 1), 1 / 3.0), 0.0)
mc = np.mean(np.exp(s * Xs))
print(f" simulated E[exp(s X_sum)] - 1 = {mc-1:.4f} (2 million draws)")
print(f" logs add (cumulant generating functions): {ok(abs(np.log(prod) - sum(np.log(1+x) for x in a))<1e-12)}")

print("\n== 5. Claim counts: one circle power (1 + c(z - 1))^n ==")
def series_coeffs(n, c, K=40):
    """coefficients of P(z) = (1 + c(z-1))^n = (1-c)^n (1 + c z/(1-c))^n."""
    out = []
    for k in range(K):
        gen = 1.0
        for t in range(k):
            gen *= (n - t) / (t + 1)
        out.append((1 - c) ** n * gen * (c / (1 - c)) ** k)
    return np.array(out)
def ab_fit(p):
    k = np.arange(1, len(p))
    m = p[:-1] > 1e-300
    y = p[1:][m] / p[:-1][m]
    A = np.vstack([np.ones(m.sum()), 1 / k[m]]).T
    (aa, bb), *_ = np.linalg.lstsq(A, y, rcond=None)
    return aa, bb
for label, n_, c_ in [("binomial m=5, q=0.3", 5, 0.3), ("neg. binomial r=2.5, beta=0.8", -2.5, -0.8), ("neg. binomial r=0.7, beta=3", -0.7, -3.0)]:
    p = series_coeffs(n_, c_, 30)
    if n_ > 0:
        p = p[: int(n_) + 1]
    aa, bb = ab_fit(p)
    print(f" {label}: sum p = {p.sum():.6f}, fitted a = {aa:.6f}, b = {bb:.6f}, dimension -1 - b/a = {-1-bb/aa:.6f}")
lam = 2.0
pp = np.array([exp(-lam) * lam ** k / math.factorial(k) for k in range(30)])
aa, bb = ab_fit(pp)
print(f" Poisson lambda=2: a = {aa:.2e} (zero), b = {bb:.6f}: dimension -1 - b/a is infinite")
# fractional positive dimension gives negative probabilities
for n_ in [2.5, 3.3, 1.5, 0.5]:
    p = series_coeffs(n_, 0.3, 12)
    neg = [k for k in range(12) if p[k] < 0]
    print(f" dimension {n_}: first negative coefficient at k = {neg[0] if neg else None}  (p_k = {p[neg[0]]:.2e})" if neg else f" dimension {n_}: no negative coefficient up to k=11")
# scan the dimension line for the zero-truncated (a,b,1) family: p_k for k>=1 must all have one sign
def sign_pattern(n_, c_, K=60, start=1):
    x = c_ / (1 - c_)
    signs = set()
    for k in range(start, K):
        gen = 1.0
        for t in range(k):
            gen *= (n_ - t) / (t + 1)
        val = gen * x ** k
        if abs(val) > 1e-300:
            signs.add(int(np.sign(val)))
    return signs
def allowed(n_, truncated):
    # some c (either sign) must give probabilities of one sign; untruncated also needs p_0 = (1-c)^n in [0,1]
    for c_ in [0.3, -0.8]:
        sg = sign_pattern(n_, c_, start=1)
        if len(sg) != 1:
            continue
        if truncated:
            return True
        p0 = (1 - c_) ** n_
        if sg == {1} and p0 <= 1:
            return True
    return False
grid = [-3, -1.2, -0.4, 0.25, 0.5, 0.9, 1.0, 1.4, 2.0, 2.6, 3.0, 4.5]
print(" count models allowed at dimension n:          ", {g: allowed(g, False) for g in grid})
print(" allowed once the zero count is removed:        ", {g: allowed(g, True) for g in grid})
# divisibility: a binomial(m) splits into k iid parts only if k | m
for m in [6, 7]:
    res = []
    for k in [2, 3, 4, 5, 6, 7]:
        p = series_coeffs(m / k, 0.3, 40)
        res.append((k, bool(np.all(p > -1e-15))))
    print(f" binomial m={m}: k-th circle root has no negative coefficient for k =", [k for k, g in res if g])
# real counts: claims per week in the Danish fire data
dd = load("danishuni.rda", "danishuni")
days = dd["Date"].to_numpy()
weeks = ((days - days.min()) // 7).astype(int)
cnt = np.bincount(weeks)
cnt = cnt[:-1]   # drop the last, partial week
vals, freq = np.unique(cnt, return_counts=True)
mean, var = cnt.mean(), cnt.var(ddof=1)
print(f" Danish fire claims per week: {len(cnt)} weeks, mean {mean:.4f}, variance {var:.4f}")
r_hat = mean ** 2 / (var - mean)
beta_hat = (var - mean) / mean
print(f"  overdispersed, so negative binomial: r = {r_hat:.3f}, beta = {beta_hat:.4f}; dimension n = -r = {-r_hat:.3f}; pole of the PGF at z = 1 + 1/beta = {1+1/beta_hat:.3f}")
# chi-square style comparison
from scipy import stats
def nb_pmf(k, r, b): return stats.nbinom.pmf(k, r, 1 / (1 + b))
ks = np.arange(0, cnt.max() + 1)
obs = np.array([np.sum(cnt == k) for k in ks])
e_p = len(cnt) * stats.poisson.pmf(ks, mean)
e_nb = len(cnt) * nb_pmf(ks, r_hat, beta_hat)
ll_p = np.sum(obs * np.log(np.maximum(stats.poisson.pmf(ks, mean), 1e-300)))
ll_nb = np.sum(obs * np.log(np.maximum(nb_pmf(ks, r_hat, beta_hat), 1e-300)))
print(f"  log-likelihood Poisson {ll_p:.1f}, negative binomial {ll_nb:.1f} (moment fit)")
k_tail = 9
print(f"  weeks with {k_tail}+ claims: observed {obs[k_tail:].sum()}, Poisson expects {e_p[k_tail:].sum():.1f}, neg. binomial {e_nb[k_tail:].sum():.1f}")
aa, bb = ab_fit(obs / obs.sum())
yrs = (days[:len(weeks)] - days.min()) / 365.25
print("  claims per year 1980-1990:", np.bincount(((days - days.min()) // 365.25).astype(int))[:11].tolist(), "(a rising trend adds some of the spread)")
print(f"  (a,b) straight-line fit on observed ratios: a = {aa:.3f}, b = {bb:.3f}, dimension -1 - b/a = {-1-bb/aa:.2f} (noisy; uses the sparse tail)")

print("\n== 6. The infinity test: the order where a moment blows up is a dimension ==")
def hill(x, k):
    xs = np.sort(x)[::-1]
    return 1 / np.mean(np.log(xs[:k] / xs[k]))
for dim in [1, 2, 3, 4]:
    nn = 400000
    # uniform point in the unit d-ball: radius R has P(R < r) = r^d
    R = rng.random(nn) ** (1 / dim)
    Xd = 1 / R
    print(f" 1/distance to a random point in a {dim}-ball: Hill tail index (k=4000) = {hill(Xd, 4000):.3f}")
# Lomax moments as a function of the dimension alpha: poles at alpha = 1..k
def lomax_moment(alpha, k, sig=1.0):
    den = np.prod([alpha - t for t in range(1, k + 1)])
    return sig ** k * math.factorial(k) / den
print(" Lomax E[X^3] near its poles:", {al: round(lomax_moment(al, 3), 2) for al in [3.01, 3.5, 2.99, 2.01, 1.01]})
alpha_l, sig = 2.5, 1.0
Y = sig * (rng.random(2_000_000) ** (-1 / alpha_l) - 1)
print(f" Lomax alpha={alpha_l}: survival (1 + x)^-alpha is a circle power with n = -alpha; sample mean {Y.mean():.4f} vs exact {lomax_moment(alpha_l,1):.4f}; E[X^2] exact {lomax_moment(alpha_l,2):.4f}; E[X^3] infinite")
# Danish fire: Hill estimate, and the max/sum ratio for each moment order
x = load("danish.rda", "danish")
print(f" Danish fire, {len(x)} losses (million DKK, 1985 values): mean {x.mean():.3f}, max {x.max():.2f}")
for k in [50, 100, 200, 300, 500]:
    print(f"  Hill tail index with top {k}: {hill(x, k):.3f}")
xs = np.sort(x)[::-1]
for pwr in [0.5, 1, 1.5, 2, 3]:
    ratio = np.max(x ** pwr) / np.sum(x ** pwr)
    print(f"  max/sum of x^{pwr}: {ratio:.4f}  ({'small: moment looks finite' if ratio < 0.05 else 'one claim dominates: moment looks infinite'})")
# ratio as the sample grows (in time order)
for pwr in [1, 2]:
    seq = [np.max(x[:m] ** pwr) / np.sum(x[:m] ** pwr) for m in [250, 500, 1000, 2000, len(x)]]
    print(f"  max/sum of x^{pwr} at n = 250, 500, 1000, 2000, 2167:", np.round(seq, 4))

print("\n== 7. Ruin: Lundberg's root must sit left of the first infinity ==")
from scipy.optimize import brentq
theta = 0.2
def lundberg(M, mean, wall):
    f = lambda R: M(R) - 1 - (1 + theta) * mean * R
    hi = wall * (1 - 1e-12)
    return brentq(f, 1e-12 * wall, hi) if f(hi) > 0 else None
R_exp = lundberg(lambda r: 1 / (1 - r), 1.0, 1.0)
print(f" exponential claims (mean 1): pole of order 1 at s = 1; R = {R_exp:.6f} (theory theta/(1+theta) = {theta/(1+theta):.6f})")
for kk in [2, 3, 5]:
    Rg = lundberg(lambda r, kk=kk: (1 - r / kk) ** (-kk), 1.0, kk)
    print(f" gamma shape {kk} (mean 1): pole of order {kk} at s = {kk}; R = {Rg:.6f}")
# Pollaczek-Khinchine simulation: psi(u) = P(sum of Geometric number of ladder heights > u)
def pk_sim(draw_ladder, u_list, n=400000):
    pgo = 1 / (1 + theta)
    Ncount = rng.geometric(1 - pgo, n) - 1     # P(N = k) = (1-pgo) pgo^k
    tot = np.zeros(n)
    maxk = Ncount.max()
    for kk in range(1, maxk + 1):
        idx = Ncount >= kk
        tot[idx] += draw_ladder(idx.sum())
    return [np.mean(tot > u) for u in u_list]
us = [0, 5, 10, 20, 40]
sim_e = pk_sim(lambda m: rng.exponential(1.0, m), us)
print(" exponential, simulated ruin psi(u):", np.round(sim_e, 5))
print(" exponential, exact e^(-Ru)/(1+theta):", np.round([exp(-R_exp * u) / (1 + theta) for u in us], 5))
# Lomax alpha = 2.5, mean 1/(alpha-1): no Lundberg root; ladder heights are Lomax(alpha-1)
al = 2.5
mu_l = 1 / (al - 1)
print(f" Lomax alpha={al}: M(s) is infinite for every s > 0, so no Lundberg root exists")
us2 = [10, 30, 100, 300]
sim_l = pk_sim(lambda m: (rng.random(m) ** (-1 / (al - 1)) - 1), us2, n=1_000_000)
asym = [(1 / (theta * mu_l)) * (1 + u) ** (-(al - 1)) / (al - 1) for u in us2]
print(" Lomax, simulated psi(u):", np.round(sim_l, 5), " power-law asymptote (1/(theta mu)) int_u S:", np.round(asym, 5))
print(f"  ruin decays like u^-(alpha-1) = u^-{al-1}: one dimension below the tail")
# Danish data: the empirical MGF is always finite, so a finite sample always yields a root
def emp_R(sample):
    m = sample.mean()
    f = lambda R: np.mean(np.exp(R * sample)) - 1 - (1 + theta) * m * R
    hi = 1.0
    while f(hi) < 0:
        hi *= 2
    return brentq(f, 1e-9, hi)
R_all = emp_R(x)
R_drop = emp_R(np.sort(x)[:-1])
R_drop5 = emp_R(np.sort(x)[:-5])
print(f" Danish fire, empirical Lundberg R: all claims {R_all:.5f}; without the largest {R_drop:.5f}; without the top 5 {R_drop5:.5f}")
print("  a finite sample always has a root, and it moves with the largest claims: the infinity is invisible to it")

print("\n== 8. Mortality: hazard = boundary / volume ==")
for kk in [1, 2, 3]:
    w = 100.0
    xx = 60.0
    S = (1 - xx / w) ** kk
    h = 1e-6
    haz = -(np.log((1 - (xx + h) / w) ** kk) - np.log(S)) / h
    print(f" {kk}-cube survival (1 - x/omega)^{kk} at x=60, omega=100: hazard {haz:.6f}; k/(omega - x) = {kk/(w-xx):.6f}")
us07 = load("demoUsa.rda", "demoUsa")
lm = us07["USSS2007M"].to_numpy(float)
agem = us07["age"].to_numpy()
good = np.isfinite(lm) & (lm > 0)
lm, agem = lm[good], agem[good]
mum = -np.log(lm[1:] / lm[:-1])
sel = (agem[:-1] >= 40) & (agem[:-1] <= 95)
slope, icpt = np.polyfit(agem[:-1][sel], np.log(mum[sel]), 1)
print(f" US 2007 male: log hazard is a straight line from 40 to 95 (Gompertz): slope {slope:.4f}/year, doubling every {log(2)/slope:.2f} years")
mu_ilt = -np.log(lx[1:] / lx[:-1])
sel2 = (ages[:-1] >= 30) & (ages[:-1] <= 100)
sl2, _ = np.polyfit(ages[:-1][sel2], np.log(mu_ilt[sel2] - 0.0007), 1)
print(f" SOA Illustrative table: log(hazard - 0.0007) slope {sl2:.5f} = log of {exp(sl2):.5f} per year (Makeham c = 10^0.04 = {10**0.04:.5f})")
print(" a pole at a finite age would make 1/hazard fall linearly to zero; Gompertz makes it decay exponentially and never reach zero")
inv = 1 / mum
for a_ in [80, 90, 100, 105, 108]:
    ii = np.where(agem[:-1] == a_)[0]
    if len(ii):
        print(f"  US 2007 male, 1/hazard at {a_}: {inv[ii[0]]:.3f} years")
