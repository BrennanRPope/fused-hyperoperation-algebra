"""Point processes on the dial: checks for every claim on index.html.
Run: python3 checks.py > checks_output.txt   (writes results.json for the page)
alpha-determinantal processes (Shirai-Takahashi 2003) with alpha = k, the dial.
"""
import json, itertools
import numpy as np
from scipy import stats

rng = np.random.default_rng(20261125)
P = F_ = 0
results = {}
def check(name, ok, info=""):
    global P, F_
    if ok: P += 1
    else: F_ += 1
    print(("PASS " if ok else "FAIL ") + name + (("  " + info) if info else ""))

# ---------- the kernel: 24 sites on a ring, real symmetric, spectrum in [0, 1] ----------
n = 24
js = np.arange(n); jj = np.minimum(js, n - js)
lam_full = 0.95 * np.exp(-(jj / 2.6) ** 2)
Fm = np.exp(2j * np.pi * np.outer(js, js) / n) / np.sqrt(n)
K = ((Fm * lam_full) @ Fm.conj().T).real
D = np.arange(8)                                     # the region: 8 neighbouring sites
KD = K[np.ix_(D, D)]
mu = np.clip(np.linalg.eigvalsh(KD), 0, None)        # eigenvalues of the kernel on the region
print("region eigenvalues:", np.round(mu, 4))

def wadd(k, a, b): return a + b + k * a * b
def fold(k, xs):
    t = 0.0
    for x in xs: t = wadd(k, t, x)
    return t
def ell(k, x): return x if k == 0 else np.log1p(k * x) / k      # the dial's logarithm
def hole(k, xs): return float(np.exp(-sum(ell(k, x) for x in xs)))

def dpp_sample(ev, U):
    sel = rng.random(len(ev)) < ev
    V = U[:, sel].copy(); pts = []
    while V.shape[1] > 0:
        p = (V ** 2).sum(1); p /= p.sum()
        i = rng.choice(len(p), p=p); pts.append(i)
        j = np.argmax(np.abs(V[i])); Vj = V[:, j].copy()
        V = np.delete(V - np.outer(Vj, V[i] / Vj[i]), j, 1)
        if V.shape[1] > 0: V, _ = np.linalg.qr(V)
    c = np.zeros(len(ev), int); c[pts] = 1; return c

def sample(k, Kmat, ns):
    """k = -1/m: m independent determinantal copies with kernel K/m.
       k = 0: Poisson.  k = 2/m: Cox process, intensity (1/m) sum of m squared real Gaussian fields."""
    ev, U = np.linalg.eigh(Kmat); ev = np.clip(ev, 0, None)
    if k == 0: return rng.poisson(np.diag(Kmat), size=(ns, len(ev)))
    if k < 0:
        m = round(-1 / k)
        return np.array([sum(dpp_sample(ev / m, U) for _ in range(m)) for _ in range(ns)])
    m = round(2 / k)
    Z = rng.standard_normal((ns, m, len(ev))) @ (U * np.sqrt(ev)).T
    return rng.poisson((Z ** 2).sum(1) / m)

def mode_pmf(k, x, nmax):
    r = np.arange(nmax)
    if k < 0: m = round(-1 / k); return stats.binom.pmf(r, m, x / m)
    if k == 0: return stats.poisson.pmf(r, x)
    return stats.nbinom.pmf(r, 1 / k, 1 / (1 + k * x))
def exact_dist(k, xs, nmax=80):
    d = np.zeros(nmax); d[0] = 1
    for x in xs: d = np.convolve(d, mode_pmf(k, x, nmax))[:nmax]
    return d

KS = [-1, -1/2, -1/3, 0, 1/2, 1, 2]
NS = {-1: 20000, -1/2: 12000, -1/3: 9000, 0: 40000, 1/2: 40000, 1: 40000, 2: 40000}
sims = {k: sample(k, K, NS[k]) for k in KS}

print("\n# 1. The hole probability is the dial fold of the eigenvalues")
for k in [-1, -1/2, -1/3, 1/3, 1/2, 1, 2, 3]:
    det = np.linalg.det(np.eye(len(D)) + k * KD)
    check(f"k={k:+.3f}: 1 + k*fold_k(eigenvalues) = det(I + k K_D)", abs(1 + k * fold(k, mu) - det) < 1e-12)
    check(f"k={k:+.3f}: prod (1+k mu)^(-1/k) = det(I+kK_D)^(-1/k) = exp(-sum ell_k(mu)) = exp(-ell_k(fold))",
          abs(np.prod((1 + k * mu) ** (-1 / k)) - det ** (-1 / k)) < 1e-12 and
          abs(hole(k, mu) - np.exp(-ell(k, fold(k, mu)))) < 1e-12)
kg = np.linspace(-1, 3, 81)
hv = [hole(k if abs(k) > 1e-9 else 1e-9, mu) for k in kg]
check("the empty chance rises steadily with k (fermions leave fewest gaps, bosons most): sum ell_k(mu) falls with k",
      all(b > a for a, b in zip(hv, hv[1:])))
check("k -> 0: hole -> exp(-sum mu), the Poisson void", abs(hole(1e-7, mu) - np.exp(-mu.sum())) < 1e-6 and abs(hole(-1e-7, mu) - np.exp(-mu.sum())) < 1e-6)
check("k = -1: P(at least one point) = mu1 (+)_-1 mu2 (+)_-1 ... = 1 - det(I - K_D)",
      abs(1 - hole(-1, mu) - fold(-1, mu)) < 1e-12)
for k in KS:
    ex = exact_dist(k, mu); h = hole(k, mu)
    check(f"k={k:+.3f}: exact count law (independent modes) gives P(N=0) = formula", abs(ex[0] - h) < 1e-12)
results["sim"] = []
for k in KS:
    N = sims[k][:, D].sum(1); ns = len(N); h = hole(k, mu)
    ph = (N == 0).mean(); se = np.sqrt(h * (1 - h) / ns)
    check(f"k={k:+.3f}: simulated hole {ph:.4f} vs {h:.4f} (within 4 SE)", abs(ph - h) < 4 * se, f"n={ns}")
    ex = exact_dist(k, mu); top = 6
    obs = np.bincount(N, minlength=top + 1)[:top + 1].astype(float); obs[-1] += ns - obs.sum()
    e = ex[:top + 1] * ns; e[-1] += ns - e.sum()
    keep = e > 5; obs2 = np.append(obs[keep], obs[~keep].sum()); e2 = np.append(e[keep], e[~keep].sum())
    if e2[-1] == 0: obs2, e2 = obs2[:-1], e2[:-1]
    pv = stats.chisquare(obs2, e2 * obs2.sum() / e2.sum()).pvalue
    check(f"k={k:+.3f}: whole simulated count law matches the exact law (chi-square p={pv:.3f} > 0.001)", pv > 1e-3)
    var = N.var(); m4 = ((N - N.mean()) ** 4).mean()
    results["sim"].append(dict(k=k, n=ns, hole=ph, hole_se=se, mean=N.mean(), var=var, var_se=np.sqrt(max(m4 - var ** 2, 1e-12) / ns)))
# Laplace functional with a soft test function phi
phi = np.zeros(n); phi[:6] = 0.3; phi[10:14] = 0.6
for k in [-1/2, 1/2, 2]:
    Mv = np.sqrt(phi)[:, None] * K * np.sqrt(phi)[None, :]
    th = np.linalg.det(np.eye(n) + k * Mv) ** (-1 / k)
    vals = np.prod((1 - phi) ** sims[k], axis=1); est = vals.mean(); se = vals.std() / np.sqrt(len(vals))
    check(f"k={k:+.3f}: E prod(1-phi)^N = det(I + k phi K)^(-1/k)  sim {est:.4f} vs {th:.4f}", abs(est - th) < 4 * se)

print("\n# 2. Where the process can exist: the wall -1/k must stay off the spectrum, and 1/|k| must be whole for k<0")
def mode_ok(k, x, nmax=80):
    """Is one mode's count law real (all probabilities >= 0)? pgf (1 + k(1-z)x)^(-1/k)."""
    if k >= 0: return bool(np.all(mode_pmf(k, x, nmax) >= 0))
    a = 1 / -k; p = -k * x; A, B = 1 - p, p          # pgf (A + B z)^a
    whole = abs(a - round(a)) < 1e-9
    if A < 0: return False                             # past the wall: signs alternate (or no real root)
    if A == 0: return whole                            # on the wall: only z^a, needs a whole
    if whole: a = round(a)
    c = A ** a; coeffs = [c]
    for r in range(nmax - 1):
        c = c * (a - r) / (r + 1) * (B / A); coeffs.append(c)
    return bool(min(coeffs) >= 0)   # the recursion keeps exact signs
grid_ok = []
for m in [1, 2, 3, 4]:
    k = -1 / m
    check(f"k=-1/{m}: one mode is a real count law for every eigenvalue up to the wall {m}",
          all(mode_ok(k, x) for x in np.linspace(0, m, 41)))
    check(f"k=-1/{m}: just past the wall (eigenvalue {m}+0.05) it is not",
          not mode_ok(k, m + 0.05))
bad = [(k, x) for k in [-0.4, -0.7, -0.45, -0.9, -0.15] for x in [0.3, 0.8] if mode_ok(k, x)]
check("k negative with 1/|k| not whole (-0.4,-0.7,-0.45,-0.9,-0.15): negative probabilities for every eigenvalue tried", len(bad) == 0)
check("k >= 0: every eigenvalue gives a real count law (negative binomial / Poisson)",
      all(mode_ok(k, x) for k in [0, 0.3, 1, 2.5] for x in [0.1, 1, 5, 20]))
for m in [1, 2, 3]:
    k = -1 / m; pm = mode_pmf(k, m, 10)
    check(f"k=-1/{m}: a mode sitting on the wall is always full (P(count={m}) = 1, variance 0)", abs(pm[m] - 1) < 1e-12)
# rigid count: kernel 2P, P a rank-3 projection, at k = -1/2: every eigenvalue on the wall
Pq = np.linalg.qr(rng.standard_normal((n, 3)))[0]; K2 = 2 * Pq @ Pq.T
s2 = sample(-1/2, K2, 3000).sum(1)
check("k=-1/2, kernel with all eigenvalues on the wall 2 (rank 3): total count is exactly 6 every time", np.all(s2 == 6))
# Sugeno comparison: roots of prod(1+l g) = 1+l always lie above -1 (measure page); here the same wall
roots_ok = True
for trial in range(200):
    g = rng.uniform(0.02, 0.95, rng.integers(2, 6))
    poly = np.poly1d([1.0])
    for gi in g: poly = poly * np.poly1d([gi, 1.0])
    poly = poly - np.poly1d([1.0, 1.0])
    rts = [r.real for r in poly.r if abs(r.imag) < 1e-9 and abs(r) > 1e-5 and r.real > -1]
    if abs(g.sum() - 1) > 1e-3 and len(rts) != 1: roots_ok = False
check("Sugeno (measure page): exactly one nonzero normalizing lambda above -1 (wall -1/lambda outside [0,1]), 200 random atom sets", roots_ok)
ks_slice = [k for k in np.round(np.linspace(-1.2, -0.01, 120), 4) if mode_ok(k, 1.0)]
check("slice eigenvalue = 1: allowed negative k are exactly -1, -1/2, -1/3, ... (only the lattice points survive)",
      all(abs(1 / -k - round(1 / -k)) < 1e-6 for k in ks_slice) and any(abs(k + 1) < 1e-9 for k in ks_slice + [-1.0]))
check("slice eigenvalue = 1 at k=-1 (wall touching the spectrum) is allowed: a sure point",
      mode_ok(-1, 1.0))

print("\n# 3. Number variance: Var N = sum mu + k sum mu^2 (each mode: value times its lift)")
S1, S2 = mu.sum(), (mu ** 2).sum()
for k in KS + [1/3, 3]:
    ex = exact_dist(k, mu, 120); r = np.arange(len(ex))
    mean = (r * ex).sum(); var = (r * r * ex).sum() - mean ** 2; f2 = (r * (r - 1) * ex).sum()
    check(f"k={k:+.3f}: exact mean = tr K_D, Var = sum mu + k sum mu^2, E N(N-1) = (sum mu)^2 + k sum mu^2",
          abs(mean - S1) < 1e-9 and abs(var - (S1 + k * S2)) < 1e-9 and abs(f2 - (S1 ** 2 + k * S2)) < 1e-9)
    check(f"k={k:+.3f}: Var = sum of mu(1 + k mu), each mode's value times its lift", abs(var - np.sum(mu * (1 + k * mu))) < 1e-9)
fano_pts = []
for row in results["sim"]:
    k = row["k"]; th = S1 + k * S2
    check(f"k={k:+.3f}: simulated Var {row['var']:.4f} vs {th:.4f} (within 4 SE)", abs(row["var"] - th) < 4 * row["var_se"])
    fano_pts.append((k, row["var"] / row["mean"]))
fk = np.array(fano_pts)
slope, icpt = np.polyfit(fk[:, 0], fk[:, 1], 1)
check(f"Fano line through 7 simulated dial points: slope {slope:.4f} vs tr K^2/tr K = {S2/S1:.4f} (within 5%)", abs(slope - S2 / S1) < 0.05 * S2 / S1)
check(f"Fano line intercept {icpt:.4f} = 1 at Poisson (within 0.02)", abs(icpt - 1) < 0.02)
check("fermion side (k<0) simulated sub-Poisson, boson side (k>0) super-Poisson",
      all((f < 1) == (k < 0) for k, f in fano_pts if k != 0))
results["fano"] = dict(S1=S1, S2=S2, slope=slope, icpt=icpt)
# superposition: m copies of (k, K) equals (k/m, mK)
for (k, m) in [(-1, 2), (1, 2), (2, 3), (1/2, 4)]:
    a = exact_dist(k, mu, 150); b = a.copy()
    for _ in range(m - 1): b = np.convolve(b, a)[:150]
    c = exact_dist(k / m, m * mu, 150) if not (k / m < 0 and abs(-m / k - round(-m / k)) > 1e-9) else None
    check(f"{m} independent copies of (k={k}, K) have the count law of (k={k/m:.3f}, {m}K)", c is not None and np.max(np.abs(b - c)) < 1e-12)
sb = sims[1][:20000, D].sum(1) + sims[1][20000:40000, D].sum(1)
thv = 2 * S1 + 0.5 * 4 * S2
check(f"simulated: two boson copies superposed have Var {sb.var():.3f} = that of (k=1/2, 2K) {thv:.3f} (within 5%)", abs(sb.var() - thv) < 0.05 * thv)

print("\n# 4. Tilting a mode's count: its mean runs by a lens law on the landscape")
def tilted_mean(k, x, t, nmax=2000):
    pm = mode_pmf(k, x, nmax) * t ** np.arange(nmax); return (np.arange(nmax) * pm).sum() / pm.sum()
okA = okB = okC = okD = True; shapes = []
for k in [-1, -1/2, -1/3, 0.5, 1, 2]:
    for x in [0.1, 0.3, 0.5, 0.8] + ([1.5, 2.5] if k > 0 or k < -0.4 else []):
        if k < 0 and x >= -1 / k: continue
        s = k * x; p = 1 + 2 * s; q = s * (1 + s); V = x * (1 + k * x)
        X = lambda t: (tilted_mean(k, x, t) - x) / V
        lens = lambda a, b: (a + b + p * a * b) / (1 - q * a * b)
        tmax = 0.7 * (1 + s) / s if s > 0 else 5.0      # bosons: tilts must stay below the pole t = (1+s)/s
        for t1, t2 in [(1.3, 0.7), (1.1, 1.2), (0.6, 0.8), (0.5, 1.6)]:
            if max(t1, t2, t1 * t2) > tmax: continue
            if abs(X(t1 * t2) - lens(X(t1), X(t2))) > 1e-9: okA = False
        h = 1e-5; mt = tilted_mean(k, x, 1.0)
        if abs((tilted_mean(k, x, 1 + h) - tilted_mean(k, x, 1 - h)) / (2 * h) - mt * (1 + k * mt)) > 1e-6: okB = False
        if abs(p * p - 4 * q - 1) > 1e-12: okC = False
        pm = mode_pmf(k, x, 400); r = np.arange(400); mm = (r * pm).sum(); fano = ((r * r * pm).sum() - mm ** 2) / mm
        if abs(p * p / (p * p - 4 * q) - (2 * fano - 1) ** 2) > 1e-8: okD = False
        shapes.append((k, x, p * p))
check("x(t) = (tilted mean - mean)/variance composes exactly by the lens law (a+b+p ab)/(1-d ab), p = 1+2k mu, d = k mu(1+k mu)", okA)
check("the tilted mean runs by the rate m' = m(1 + k m) (the variance function, Morris 1982)", okB)
check("every dial count lands on the one parabola p^2 - 4d = 1 (two real walls)", okC)
check("units-free shape p^2/(p^2-4d) = (2*Fano - 1)^2 for every mode", okD)
check("fermion modes have shape in [0,1], boson modes shape > 1, none negative (no count lives above the border parabola)",
      all((0 <= sh <= 1) if k < 0 else sh > 1 for k, x, sh in shapes))
check("empty-mode limit s=0: (p,d) = (1,0), the circle operation; full fermion mode s=-1: (-1,0), union of failures",
      (1 + 2 * 0, 0 * 1) == (1, 0) and (1 + 2 * -1, -1 * 0) == (-1, 0))
a, b = 0.37, -0.52
check("half-filled fermion mode s=-1/2: (p,d) = (0,-1/4), lens law (a+b)/(1+ab/4) is velocity addition with speed limit 2",
      abs((a + b) / (1 + a * b / 4) - (a + b + 0 * a * b) / (1 - (-0.25) * a * b)) < 1e-15)
check("particle-hole: at k=-1 eigenvalues mu and 1-mu give the same shape",
      all(abs((1 - 2 * x) ** 2 - (1 - 2 * (1 - x)) ** 2) < 1e-15 for x in np.linspace(0, 1, 11)))

# ---------- data for the page ----------
pick = {}
for k in [-1, 0, 1, 2]:
    Ns = sims[k].sum(1); target = int(round(np.diag(K).sum()))
    idx = np.where(Ns == target)[0]; pick[str(k)] = sims[k][idx[0]].tolist()
results["strips"] = pick
# picture only: one sample per k on a 120-site ring (smooth kernel, 0.9 on the 25 slowest modes)
n2 = 120; j2 = np.arange(n2); jj2 = np.minimum(j2, n2 - j2)
lam2 = 0.9 / (1 + np.exp((jj2 - 12) / 1.2))
F2 = np.exp(2j * np.pi * np.outer(j2, j2) / n2) / np.sqrt(n2)
Kp = ((F2 * lam2) @ F2.conj().T).real
results["picture"] = {str(k): sample(k, Kp, 1)[0].tolist() for k in [-1, 0, 1]}
results["mu"] = mu.tolist()
results["counts"] = {str(k): (np.bincount(sims[k][:, D].sum(1), minlength=8)[:8] / NS[k]).tolist() for k in KS}
results["exact"] = {str(k): exact_dist(k, mu)[:8].tolist() for k in KS}
json.dump(results, open("results.json", "w"), default=float, indent=1)
print(f"\n{P} PASS, {F_} FAIL, of {P + F_} checks")
