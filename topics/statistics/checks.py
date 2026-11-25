"""Numeric and symbolic checks for the statistics page (math-theory/statistics).
Everything Lean can't do here: real exponentials, limits, the normal law, Poisson,
plus brute-force cross-checks of the Lean identities. Run: python3 checks.py"""
import itertools, math, random
import numpy as np, sympy as sp
from scipy import stats

random.seed(1); np.random.seed(1)
results = []
def check(name, ok, detail=""):
    results.append((name, bool(ok), detail)); print(("PASS " if ok else "FAIL ") + name + ("  " + detail if detail else ""))

def wadd(k, a, b): return a + b + k*a*b

# 1. Piling-up for unequal coins: brute force over all outcomes
for n in range(1, 7):
    ps = [random.random() for _ in range(n)]
    pe = po = pnone = 0.0
    for bits in itertools.product([0,1], repeat=n):
        w = np.prod([p if b else 1-p for p,b in zip(ps,bits)])
        if sum(bits) % 2: po += w
        else: pe += w
        if sum(bits) == 0: pnone += w
    f2 = 0.0; f1 = 0.0
    for p in ps: f2 = wadd(-2, p, f2); f1 = wadd(-1, p, f1)
    check(f"parity, {n} unequal coins: P(even)-P(odd) = 1-2*fold(-2)", abs((pe-po) - (1-2*f2)) < 1e-12)
    check(f"none, {n} unequal coins: P(none) = 1-fold(-1)", abs(pnone - (1-f1)) < 1e-12)

# 2. Inclusion-exclusion with dependent events on a random finite space
for trial in range(5):
    m, n = 40, 5
    w = np.random.rand(m); w /= w.sum()
    A = (np.random.rand(n, m) < np.random.rand(n, 1)).astype(float)
    union = (w * (A.max(axis=0))).sum()
    ie = 0.0
    for r in range(1, n+1):
        for g in itertools.combinations(range(n), r):
            ie += (-1)**(r+1) * (w * A[list(g)].prod(axis=0)).sum()
    check(f"inclusion-exclusion, 5 dependent events, trial {trial}", abs(union-ie) < 1e-12)

# 3. Dependence is the weight: k = -lift, and its range
ok = True; kmin = 0
for _ in range(20000):
    m = 30; w = np.random.rand(m); w /= w.sum()
    a = (np.random.rand(m) < random.random()); b = (np.random.rand(m) < random.random())
    if random.random() < 0.3: b = a | (np.random.rand(m) < 0.2)   # nested-ish, positive dependence
    p, q, pab = (w*a).sum(), (w*b).sum(), (w*(a&b)).sum()
    if p == 0 or q == 0: continue
    lift = pab/(p*q); k = -lift
    ok &= abs((w*(a|b)).sum() - wadd(k, p, q)) < 1e-12
    ok &= (-1/max(p,q) - 1e-12 <= k <= 1e-12)
check("P(A or B) = p (+)_k q with k = -lift, and -1/max(p,q) <= k <= 0", ok)

# 4. The k-moment dial: Poisson, binomial limit, factorial moments
for k in [-2, -1, -0.5, 0.7, 1, 2]:
    for l1, l2 in [(0.3, 1.1), (2.0, 0.5)]:
        ak = lambda lam: (math.exp(k*lam)-1)/k
        check(f"Poisson k-moment adds by (+)_k, k={k}, lambda={l1}+{l2}", abs(ak(l1+l2) - wadd(k, ak(l1), ak(l2))) < 1e-12)
    lam = 1.3
    errs = []
    for n in [10, 100, 1000, 10000]:
        a_bin = ((1 + k*lam/n)**n - 1)/k
        errs.append(abs(a_bin - (math.exp(k*lam)-1)/k))
    check(f"binomial(n, lam/n) k-moment -> Poisson's (e^(k lam)-1)/k, k={k}", errs[-1] < errs[0]/500, f"errors {errs[0]:.2e} -> {errs[-1]:.2e}")
n, p = 7, 0.37
for j in range(0, 8):
    fm = sum(math.comb(x, j) * math.comb(n, x) * p**x * (1-p)**(n-x) for x in range(n+1))
    check(f"binomial E[C(X,{j})] = C(n,{j}) p^{j} (expected lit {j}-corners)", abs(fm - math.comb(n, j)*p**j) < 1e-12)

# 5. Moments -> cumulants is the circle operation's log, log(1+x), on the moment series
s = sp.symbols('s'); m1, m2, m3, m4 = sp.symbols('m1 m2 m3 m4')
M = 1 + m1*s + m2*s**2/2 + m3*s**3/6 + m4*s**4/24
K = sp.series(sp.log(M), s, 0, 5).removeO()
kap = [sp.expand(K.coeff(s, j)*sp.factorial(j)) for j in range(1, 5)]
lean_k = [m1, m2 - m1**2, m3 - 3*m2*m1 + 2*m1**3, m4 - 4*m3*m1 - 3*m2**2 + 12*m2*m1**2 - 6*m1**4]
for j in range(4):
    check(f"cumulant k{j+1} from log(1 + moment series) matches the Lean formula", sp.simplify(kap[j] - lean_k[j]) == 0)
a, b = sp.symbols('a b')
check("MGF multiplies, so (M-1) combines by the circle operation: (MX MY - 1) = (MX-1) (+) (MY-1)",
      sp.simplify((1+a)*(1+b) - 1 - (a + b + a*b)) == 0)

# 6. Normal and CLT
mu, sig = sp.symbols('mu sigma', positive=True)
Mn = sp.exp(mu*s + sig**2*s**2/2)
check("normal MGF - 1 = E(mu s + sigma^2 s^2/2), E(x)=e^x-1 (the circle operation's exponential)", sp.simplify(Mn - 1 - (sp.exp(mu*s + sig**2*s**2/2) - 1)) == 0)
E = lambda x: sp.exp(x) - 1
check("normals add because E(x+y) = E(x) (+) E(y)", sp.simplify(E(a+b) - (E(a) + E(b) + E(a)*E(b))) == 0)
p = 0.3; prev = None; errs = []
for n in [10, 40, 160, 640, 2560]:
    mean, sd = n*p, math.sqrt(n*p*(1-p))
    x = np.arange(n+1); cdf = stats.binom.cdf(x, n, p)
    err = np.max(np.abs(cdf - stats.norm.cdf((x + 0.5 - mean)/sd)))
    errs.append(err)
check("CLT: coin p=0.3, distance to the normal shrinks ~1/sqrt(n)", errs[-1] < errs[0]/6, " -> ".join(f"{e:.4f}" for e in errs))
for n in [5, 50, 500]:
    kurt_n = (1 - 6*p*(1-p))/(n*p*(1-p)); skew2_n = (1-2*p)**2/(n*p*(1-p))
    xs = np.arange(n+1); pm = stats.binom.pmf(xs, n, p); mm = (pm*xs).sum()
    c2 = (pm*(xs-mm)**2).sum(); c3 = (pm*(xs-mm)**3).sum(); c4 = (pm*(xs-mm)**4).sum()
    check(f"binomial n={n}: kurtosis = kurt(coin)/n, skew^2 = skew^2(coin)/n (exact, as in Lean)",
          abs((c4/c2**2 - 3) - kurt_n) < 1e-9 and abs(c3**2/c2**3 - skew2_n) < 1e-9)
for sv in [0.5, 1.0, 2.0]:
    vals = [(1 + sv**2/(2*n))**n for n in [10, 100, 10000]]
    check(f"MGF climb (1 + s^2/2n)^n -> e^(s^2/2) at s={sv}", abs(vals[-1] - math.exp(sv**2/2)) < 1e-3*math.exp(sv**2/2))
p = 0.3; gaps = []
for n in [100, 400, 1600, 6400]:
    sd = math.sqrt(n*p*(1-p)); sv = 1.0
    mgf = ((1-p) + p*math.exp(sv/sd))**n * math.exp(-sv*n*p/sd)
    gaps.append(mgf - math.exp(.5))
    print(f"   standardized coin-sum MGF at s=1, n={n}: {mgf:.5f}  (normal: {math.exp(.5):.5f})")
check("standardized coin-sum MGF at s=1 -> e^(1/2), gap halves each time n is x4 (skewness term, 1/sqrt n)",
      all(1.8 < gaps[i]/gaps[i+1] < 2.2 for i in range(3)), " ".join(f"{g:.5f}" for g in gaps))

# 7. Variance = half the expected squared difference of two copies
for _ in range(5):
    v = np.random.randn(6); w = np.random.rand(6); w /= w.sum()
    var = (w*v**2).sum() - (w*v).sum()**2
    half = 0.5*sum(w[i]*w[j]*(v[i]-v[j])**2 for i in range(6) for j in range(6))
    check("Var(X) = E[(X-X')^2]/2 for a random distribution", abs(var-half) < 1e-12)

# 8. Bayes: odds, log-odds, velocity
logit = lambda p: math.log(p/(1-p))
for p, q in [(0.2, 0.7), (0.9, 0.95), (0.4, 0.1)]:
    pool = p*q/(p*q + (1-p)*(1-q)); r1, r2 = 2*p-1, 2*q-1
    check(f"pool({p},{q}): log-odds add (rapidities)", abs(logit(pool) - logit(p) - logit(q)) < 1e-12)
    check(f"pool({p},{q}): rho = 2p-1 composes by velocity addition", abs((2*pool-1) - (r1+r2)/(1+r1*r2)) < 1e-12)
    check(f"log-odds is twice the rapidity: logit p = 2 artanh(2p-1)", abs(logit(p) - 2*math.atanh(2*p-1)) < 1e-12)

# 9. AMH dependence line
def amh(t, u, v):
    d = 1 - t*(1-u)*(1-v)
    return 0.0 if u*v == 0 else u*v/d
def lens(k, q, a, b): return (a + b + k*a*b)/(1 - q*a*b)
ok = True; okA = True
for _ in range(2000):
    t = random.uniform(-1, 1); p, q, r = random.random(), random.random(), random.random()
    U = lambda x, y: 1 - amh(t, 1-x, 1-y)
    ok &= abs(U(p, q) - lens(-(1+t), t, p, q)) < 1e-12
    okA &= abs(U(U(p, q), r) - U(p, U(q, r))) < 1e-12
check("union under AMH = lens law at (k, d) = (-(1+theta), theta)", ok)
check("AMH union is associative (the lens turns it into multiplication)", okA)
def two_increasing(t, N=60):
    g = np.linspace(0, 1, N); C = np.array([[amh(t, u, v) for v in g] for u in g])
    return (np.diff(np.diff(C, axis=0), axis=1) >= -1e-12).all()
check("AMH is a valid copula on the whole segment theta in [-1, 1]", all(two_increasing(t) for t in np.linspace(-1, 1, 21)))
check("…and stops being one past the cusp (theta = 1.5)", not two_increasing(1.5))
check("theta=-1 gives velocity addition, theta=0 independence, theta=1 the parabola k^2=4d",
      abs(lens(0, -1, .3, .4) - (.7/1.12)) < 1e-12 and abs(lens(-1, 0, .3, .4) - wadd(-1, .3, .4)) < 1e-12 and (-2)**2 == 4*1)

# 10. Entropy ties
var = 1.0
h_norm = 0.5*math.log(2*math.pi*math.e*var)
h_unif = math.log(math.sqrt(12*var)); h_lap = 1 + math.log(2*math.sqrt(var/2))
check("normal has the most entropy at fixed variance (vs uniform, Laplace)", h_norm > h_unif and h_norm > h_lap, f"{h_norm:.4f} > {h_lap:.4f}, {h_unif:.4f}")
n, p = 400, 0.3; pm = stats.binom.pmf(np.arange(n+1), n, p); pm = pm[pm > 0]
H = -(pm*np.log(pm)).sum(); Hn = 0.5*math.log(2*math.pi*math.e*n*p*(1-p))
check("binomial entropy -> normal's 1/2 log(2 pi e n p q)", abs(H - Hn) < 2e-3, f"{H:.5f} vs {Hn:.5f}")
H1 = -(p*math.log(p) + (1-p)*math.log(1-p))
xs = list(itertools.product([0,1], repeat=3)); Hj = 0
for bits in xs:
    w = np.prod([p if b else 1-p for b in bits]); Hj -= w*math.log(w)
check("entropy of 3 independent coins = 3 x one coin (k = 0 addition)", abs(Hj - 3*H1) < 1e-12)

npass = sum(r[1] for r in results)
print(f"\n{npass}/{len(results)} checks pass")
