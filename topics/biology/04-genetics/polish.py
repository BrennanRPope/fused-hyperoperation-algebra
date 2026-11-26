"""Polish step for the one searched number in checks.py: the best lens-law K for the chi-square (m = 4)
interference model (section 3).

checks.py fits K by scipy's bounded Brent search (minimize_scalar, default xatol 1e-5) to recombination
fractions simulated by Monte Carlo. The polish replays the same random draws (seed 7, every earlier
simulation included, so the data are identical; the replay is checked against numbers printed in
checks_output.txt), reproduces the search with an exact port of scipy's bounded Brent routine (scipy is
not installed here), then refines K with Newton steps on a finite-difference derivative.
The data are Monte Carlo estimates with noise ~0.0009 per fraction, so digits of K beyond the second
decimal mean nothing; the polish only confirms the search had converged.
Run: python3 polish.py > polish_output.txt   (about 30 s, almost all of it the replay)
"""
import math, time
import numpy as np

def fminbound(func, x1, x2, xatol=1e-5, maxfun=500):
    """Port of scipy.optimize._minimize_scalar_bounded (bounded Brent)."""
    sqrt_eps = np.sqrt(2.2e-16); golden_mean = 0.5 * (3.0 - np.sqrt(5.0))
    a, b = x1, x2
    fulc = a + golden_mean * (b - a); nfc, xf = fulc, fulc
    rat = e = 0.0; x = xf; fx = func(x); num = 1
    ffulc = fnfc = fx; xm = 0.5 * (a + b)
    tol1 = sqrt_eps * np.abs(xf) + xatol / 3.0; tol2 = 2.0 * tol1
    while np.abs(xf - xm) > (tol2 - 0.5 * (b - a)):
        golden = 1
        if np.abs(e) > tol1:
            golden = 0
            r = (xf - nfc) * (fx - ffulc); q = (xf - fulc) * (fx - fnfc)
            p = (xf - fulc) * q - (xf - nfc) * r; q = 2.0 * (q - r)
            if q > 0.0: p = -p
            q = np.abs(q); r = e; e = rat
            if (np.abs(p) < np.abs(0.5 * q * r)) and (p > q * (a - xf)) and (p < q * (b - xf)):
                rat = (p + 0.0) / q; x = xf + rat
                if ((x - a) < tol2) or ((b - x) < tol2):
                    si = np.sign(xm - xf) + ((xm - xf) == 0); rat = tol1 * si
            else:
                golden = 1
        if golden:
            e = (a - xf) if xf >= xm else (b - xf)
            rat = golden_mean * e
        si = np.sign(rat) + (rat == 0)
        x = xf + si * np.maximum(np.abs(rat), tol1)
        fu = func(x); num += 1
        if fu <= fx:
            if x >= xf: a = xf
            else: b = xf
            fulc, ffulc = nfc, fnfc; nfc, fnfc = xf, fx; xf, fx = x, fu
        else:
            if x < xf: a = x
            else: b = x
            if (fu <= fnfc) or (nfc == xf):
                fulc, ffulc = nfc, fnfc; nfc, fnfc = x, fu
            elif (fu <= ffulc) or (fulc == xf) or (fulc == nfc):
                fulc, ffulc = x, fu
        xm = 0.5 * (a + b); tol1 = sqrt_eps * np.abs(xf) + xatol / 3.0; tol2 = 2.0 * tol1
        if num >= maxfun: break
    return xf

T0 = time.perf_counter()
rng = np.random.default_rng(7)
def simulate_r(lengths, m=0, N=400000):
    """Same random draws, in the same order, as checks.py's simulate_r (two intervals only).
    Counting uses a comparison instead of np.histogram; they differ only if a uniform draw lands
    exactly on an interval edge."""
    L = sum(lengths); e1 = lengths[0]; rate = 2*(m+1); out = np.zeros((N, 2), dtype=np.int64)
    pois, uni, integ, rnd = rng.poisson, rng.uniform, rng.integers, rng.random
    for i in range(N):
        n = pois(rate*L); pts = np.sort(uni(0, L, n)); start = integ(0, m+1); co = pts[start::m+1]
        co = co[rnd(len(co)) < 0.5]
        c1 = int(np.count_nonzero(co < e1)); out[i, 0] = c1; out[i, 1] = len(co) - c1
    return out
odd = lambda x: (x % 2 == 1)
wadd = lambda k, a, b: a + b + k*a*b
lens = lambda k, q, a, b: (a + b + k*a*b)/(1 - q*a*b)
lengths = [0.15, 0.25]
C = simulate_r(lengths, m=0, N=200000)
r1 = odd(C[:, 0]).mean(); r2 = odd(C[:, 1]).mean(); r12 = odd(C[:, 0] + C[:, 1]).mean()
sync = [f"{r12:.4f}"]
for m in (1, 4):
    C = simulate_r(lengths, m=m, N=200000)
    r1 = odd(C[:, 0]).mean(); r2 = odd(C[:, 1]).mean(); dbl = (odd(C[:, 0]) & odd(C[:, 1])).mean(); sync.append(f"{dbl/(r1*r2):.3f}")
for l in (0.1, 0.25, 0.5):
    C = simulate_r([l, l], m=4, N=150000); a = odd(C[:, 0]); b = odd(C[:, 1]); sync.append(f"{(a&b).mean()/(a.mean()*b.mean()):.2f}")
rng.uniform(0, 0.98, 2); rng.uniform(0, 0.45, 3)          # the two draws checks.py makes in section 3
pairs = [(0.1, 0.1), (0.1, 0.3), (0.2, 0.2), (0.3, 0.5), (0.5, 0.5), (0.05, 0.6), (0.15, 0.15), (0.4, 0.4)]
obs = []
for l1, l2 in pairs:
    C = simulate_r([l1, l2], m=4, N=300000)
    obs.append((odd(C[:, 0]).mean(), odd(C[:, 1]).mean(), odd(C[:, 0] + C[:, 1]).mean()))
obs = np.array(obs)
t_replay = time.perf_counter() - T0
print(f"replay of checks.py's simulations: {t_replay:.1f} s")
print(f"  sync check against checks_output.txt: r12 {sync[0]} (printed 0.2749), c(m=1) {sync[1]} (0.758), c(m=4) {sync[2]} (0.463), "
      f"c at 10/25/50 cM {', '.join(sync[3:])} (0.09, 0.63, 0.97)")

f = lambda K: np.sum((lens(-2*K, -4+4*K, obs[:, 0], obs[:, 1]) - obs[:, 2])**2)
miss = lambda K: np.max(np.abs(lens(-2*K, -4+4*K, obs[:, 0], obs[:, 1]) - obs[:, 2]))
K0 = fminbound(f, -4, 1.5)                       # = minimize_scalar(f, bounds=(-4, 1.5), method="bounded")
t = time.perf_counter()
K = K0
for it in range(10):                             # Newton on f'(K) with central differences
    h = 1e-4
    d1 = (f(K+h) - f(K-h))/(2*h); d2 = (f(K+h) - 2*f(K) + f(K-h))/h**2
    step = d1/d2; K -= step
    if abs(step) < 1e-13: break
dt = time.perf_counter() - t
print("best lens-law K for the chi-square m = 4 model")
print(f"  before (bounded Brent, xatol 1e-5): K = {K0:.10f}   sum sq {f(K0):.12e}   worst miss {miss(K0):.6f}")
print(f"  after  (Newton, {it+1} steps)        : K = {K:.10f}   sum sq {f(K):.12e}   worst miss {miss(K):.6f}")
print(f"  change in K {K-K0:+.1e} (search was good to about {-math.log10(abs(K-K0)/abs(K)):.0f} significant digits; Newton to ~{-math.log10(max(abs(step),1e-16)/abs(K)):.0f})")
print(f"  page numbers: best K {K0:.2f} -> {K:.2f}, worst miss {miss(K0):.4f} -> {miss(K):.4f}, Kosambi miss {miss(0):.4f}: "
      + ("no change" if (f"{K0:.2f}", f"{miss(K0):.4f}") == (f"{K:.2f}", f"{miss(K):.4f}") else "CHANGED"))
print(f"  polish time {dt*1e3:.2f} ms; total wall time {time.perf_counter()-T0:.1f} s")
