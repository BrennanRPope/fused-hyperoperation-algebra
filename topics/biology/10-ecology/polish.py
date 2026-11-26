"""Polish step for the fitted functional-response numbers in checks.py (section 1) and identify.py.

checks.py: three log-scale fits (dial k free, Beddington-DeAngelis k = 0, Crowley-Martin k = 1) by
scipy's least_squares with default tolerances (1e-8). Only k and the AIC values were stored, so the
polish rebuilds the same simulated data (same seed, same random draws, including the Gillespie run
that comes first) and runs Levenberg-Marquardt to xtol 1e-14 from the same starting point.
identify.py: 800 bounded fits on simulated noisy replicates; the 2.5/97.5 percentiles of k are what
the page quotes. The individual fits were not stored, so they are redone the same way (same seed,
same start, same bounds) to 1e-14 and the percentiles compared.
All of this is simulated noisy data: extra digits mean nothing; the polish only confirms convergence.
numpy only (scipy is not installed here).
Run: python3 polish.py > polish_output.txt
"""
import json, math, time
import numpy as np
np.seterr(invalid="ignore", divide="ignore")   # finite-difference probes near k = -2 can step outside the model

def lm(fun, x0, lo=None, hi=None, xtol=1e-14, ftol=1e-15, maxit=500):
    """Small Levenberg-Marquardt (numpy only; scipy is not installed here). Central-difference Jacobian,
    box bounds by clipping. Returns (x, residuals, iterations)."""
    x = np.array(x0, float); n = len(x)
    lo = np.full(n, -np.inf) if lo is None else np.array(lo, float)
    hi = np.full(n, np.inf) if hi is None else np.array(hi, float)
    x = np.clip(x, lo, hi); r = fun(x); cost = r @ r; lam = 1e-6
    for it in range(1, maxit + 1):
        J = np.empty((len(r), n))
        for i in range(n):
            h = 1e-6 * max(1.0, abs(x[i])); e = np.zeros(n); e[i] = h
            J[:, i] = (fun(x + e) - fun(x - e)) / (2 * h)
        g = J.T @ r; A = J.T @ J
        free = ~(((x <= lo) & (g > 0)) | ((x >= hi) & (g < 0)))   # coordinates not pinned at a bound
        improved = False
        for _ in range(30):
            dx = np.zeros(n)
            Af = A[np.ix_(free, free)]
            try:
                dx[free] = np.linalg.solve(Af + lam * np.diag(np.diag(Af) + 1e-300), -g[free])
            except np.linalg.LinAlgError:
                lam *= 10; continue
            xn = np.clip(x + dx, lo, hi); rn = fun(xn); cn = rn @ rn
            if np.isfinite(cn) and cn <= cost:
                improved = True; break
            lam *= 10
        if not improved:
            break
        step = np.max(np.abs(xn - x) / (np.abs(x) + 1e-12)); drop = (cost - cn) / max(cost, 1e-300)
        x, r, cost = xn, rn, cn; lam = max(lam / 10, 1e-12)
        if step < xtol or drop < ftol:
            break
    return x, r, it

wadd = lambda k, a, b: a + b + k*a*b
pars = dict(al=0.8, h=0.3, ga=0.5, tau=0.6)
def feed_exact(N, P, al, h, ga, tau, k):
    x = al*N*h; y = ga*P*tau
    return al*N/(1 + wadd(k, x, y))

# ---- checks.py section 1, reproduced draw for draw ----
T0 = time.perf_counter()
rng = np.random.default_rng(10)
def gillespie(N, P, al, h, ga, tau, k, T=2e5):
    s = 0; t = 0; eaten = 0
    rates = {0: [(1, al*N), (2, ga*P)], 1: [(0, 1/h), (3, k*ga*P)], 2: [(0, 1/tau)], 3: [(1, 1/tau)]}
    while t < T:
        rs = rates[s]; tot = sum(r for _, r in rs); t += rng.exponential(1/tot)
        u = rng.uniform(0, tot); acc = 0
        for nxt, r in rs:
            acc += r
            if u <= acc:
                if s == 0 and nxt == 1: eaten += 1
                s = nxt; break
    return eaten/T
g = gillespie(5, 6, k=0.5, **pars)
t_sim = time.perf_counter() - T0
Ns = np.repeat([0.5, 1, 2, 4, 8, 16], 5); Ps = np.tile([1, 2, 4, 8, 12], 6); ktrue = 0.5
y = feed_exact(Ns, Ps, k=ktrue, **pars)*np.exp(rng.normal(0, 0.04, len(Ns)))
def model(p, k=None):
    al, h, c = np.exp(p[:3]); kk = p[3] if k is None else k
    x = al*Ns*h; z = c*Ps; return al*Ns/(1 + x + z + kk*x*z)
t = time.perf_counter()
fd, rd, itd = lm(lambda p: np.log(model(p)) - np.log(y), [0, 0, 0, 0.3])
fb, rb, _ = lm(lambda p: np.log(model(p, 0)) - np.log(y), [0, 0, 0])
fc, rc, _ = lm(lambda p: np.log(model(p, 1)) - np.log(y), [0, 0, 0])
dt = time.perf_counter() - t
n = len(y); aic = lambda r, kp: n*np.log(np.mean(r**2)) + 2*kp
A = dict(dial=aic(rd, 4), BD=aic(rb, 3), CM=aic(rc, 3))
st = json.load(open("results.json"))["fr"]
print(f"checks.py feeding-experiment fits (data rebuilt in {t_sim:.1f} s; Gillespie rate {g:.4f})")
print(f"  fitted k  before {st['kfit']:.10f}  after {fd[3]:.10f}  change {fd[3]-st['kfit']:+.1e}  "
      + ("(no change at 2 decimals)" if f"{fd[3]:.2f}" == f"{st['kfit']:.2f}" else "CHANGED at 2 decimals"))
for m in ("dial", "BD", "CM"):
    print(f"  AIC {m:4s} before {st['aic'][m]:.6f}  after {A[m]:.6f}  change {A[m]-st['aic'][m]:+.1e}")
print(f"  verdict unchanged: dial AIC lowest = {A['dial'] < min(A['BD'], A['CM'])}; |k - 0.5| < 0.25 = {abs(fd[3]-ktrue) < 0.25}")
print(f"  polish time {dt*1e3:.1f} ms")

# ---- identify.py, reproduced replicate for replicate ----
t = time.perf_counter()
rng = np.random.default_rng(11)
stored = json.load(open("identify.json")); out = {}
for reps in (5, 20):
    Ns = np.repeat([0.5, 1, 2, 4, 8, 16], reps); Ps = np.tile(np.repeat([1, 2, 4, 8, 12], 1), 6*reps//5)
    for sd in (0.1, 0.25):
        ks = []
        for b in range(200):
            y = feed_exact(Ns, Ps, k=0.5, **pars)*np.exp(rng.normal(0, sd, len(Ns)))
            def model(p):
                al, h, c = np.exp(p[:3]); x = al*Ns*h; z = c*Ps; return al*Ns/(1 + x + z + p[3]*x*z)
            f, _, _ = lm(lambda p: np.log(model(p)) - np.log(y), [0, 0, 0, 0.3], [-5, -5, -5, -2], [5, 5, 5, 10])
            ks.append(f[3])
        lo, hi = np.percentile(ks, [2.5, 97.5]); key = f"{len(Ns)} trials, {int(sd*100)}% noise"
        b0 = stored[key]
        same = f"{b0[0]:.2f}-{b0[1]:.2f}" == f"{lo:.2f}-{hi:.2f}"
        print(f"  identify {key:24s} before {b0[0]:.4f} to {b0[1]:.4f}   after {lo:.4f} to {hi:.4f}   "
              + ("(no change at 2 decimals)" if same else "changed at 2 decimals"))
print(f"  identify polish time {time.perf_counter()-t:.1f} s")
print(f"total wall time {time.perf_counter()-T0:.1f} s")
