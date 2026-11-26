"""Lee-Yang at even roundness above a = p (DRAFT). At p = 4 the gap reduces to a two-parameter family:

  (1 - x^4)^c = (1 - x^2)^c (1 + x^2)^theta * (1 + x^2)^m,   m = floor(c), theta = c - m in [0, 1).

Multiplying the weight by (1 + x^2) applies (1 - D^2) to the transform, which keeps real zeros (Laguerre), so
G_c at p = 4 has only real zeros as soon as F(c, theta)(s) = int_0^1 (1 - x^2)^c (1 + x^2)^theta cos(s x) dx does.
At theta = 0 that is the round ball (Bessel), which is the known whole-a case. Here:
  E1 the operator step F(alpha, theta + 1) = (1 - D^2) F(alpha, theta), exactly (series coefficients, 40 digits)
  E2 G_c at p = 4 equals (1 - D^2)^m F(c, theta) at sample (c, s), against direct quadrature
  E3 the series for F (Bessel functions) against direct quadrature, zeros agree
  E4 over a grid of alpha (4.2 to 35) and theta (0 to 1, 21 values): the zero count up to 3 alpha + 40 never
     changes, every zero moves down as theta grows, and no two zeros come closer than 3.1
  E5 at p = 6 the matching factor 1 + x^2 + x^4 gives the operator 1 - D^2 + D^4, which has non-real roots
"""
import json, time
import numpy as np
import mpmath as mp
from scipy.special import jv, gammaln, binom, roots_jacobi
t0 = time.perf_counter(); RES = []; OUT = {}
def check(name, ok, detail=""):
    RES.append(bool(ok)); print(("PASS " if ok else "FAIL ") + name + (f"  [{detail}]" if detail else ""), flush=True)

# --- exact series: F(alpha, theta)(s) = sum_k (-1)^k s^(2k)/(2k)! M_k, M_k = int x^(2k) (1-x^2)^alpha (1+x^2)^theta
mp.mp.dps = 40
def M(alpha, theta, k):
    h = mp.mpf(1)/2; alpha = mp.mpf(alpha); theta = mp.mpf(theta)
    return h*mp.beta(k + h, alpha + 1)*mp.hyp2f1(-theta, k + h, k + alpha + 1 + h, -1)
# (1 - D^2) on a cosine series sum a_k (-1)^k s^2k/(2k)! gives coefficients a_k + a_{k+1}
err1 = 0
for alpha, theta in [("4.3", "0.25"), ("6.7", "0.6"), ("11.1", "0.9")]:
    alpha, theta = mp.mpf(alpha), mp.mpf(theta)
    for k in range(8):
        lhs = M(alpha, theta + 1, k); rhs = M(alpha, theta, k) + M(alpha, theta, k + 1)
        err1 = max(err1, abs(lhs - rhs)/abs(lhs))
check("E1 F(alpha, theta + 1) = (1 - D^2) F(alpha, theta), coefficient by coefficient", err1 < mp.mpf(10)**-35, f"worst relative gap {mp.nstr(err1, 3)}")
OUT["E1"] = float(err1)

# --- E2: G_c(p = 4) = (1 - D^2)^m F(c, theta), checked through the series coefficients and by quadrature
def Gq(p, c, s, N=400):
    x, w = roots_jacobi(N, c, 0.0); x = (x + 1)/2; w = w/2**(c + 1)
    return (w*((1 - x**p)/(1 - x))**c*np.cos(np.outer(s, x))).sum(1)
def F_series_coeffs(alpha, theta, K):
    return [M(alpha, theta, k) for k in range(K)]
err2 = 0
for c in [4.35, 5.8, 6.2]:
    m = int(np.floor(c)); th = c - m; K = 60
    a = F_series_coeffs(c, th, K + m + 2)
    for _ in range(m):                          # apply (1 - D^2) m times: a_k -> a_k + a_{k+1}
        a = [a[k] + a[k + 1] for k in range(len(a) - 1)]
    for s in [2.0, 7.5, 13.0]:
        ser = mp.fsum((-1)**k*mp.mpf(s)**(2*k)/mp.factorial(2*k)*a[k] for k in range(K))
        q = Gq(4, c, np.array([s]))[0]
        err2 = max(err2, abs(float(ser) - q)/max(1e-12, abs(q)))
check("E2 G at p = 4 equals (1 - D^2)^m F(c, theta) (series against quadrature)", err2 < 1e-8, f"worst relative gap {err2:.1e}")
OUT["E2"] = err2

# --- Bessel series for F, normalised by s^(alpha + 1/2): sum_j binom(theta, j)(-1)^j (alpha+1)_j s^-j J_(alpha+j+1/2)(s)
def Fn(alpha, theta, s, J=160):
    tot = np.zeros_like(s)
    for j in range(J):
        cj = binom(theta, j)*(-1)**j
        if cj == 0: continue
        jj = jv(alpha + j + 0.5, s)
        with np.errstate(divide="ignore"):
            tot += cj*np.sign(jj)*np.exp(gammaln(alpha + j + 1) - gammaln(alpha + 1) - j*np.log(s) + np.log(np.abs(jj)))
    return tot
def Fq(alpha, theta, s, N=500):
    x, w = roots_jacobi(N, alpha, alpha)
    return (w*(1 + x**2)**theta*np.cos(np.outer(s, x))).sum(1)
def zeros_of(f, s):
    i = np.where(np.sign(f[1:]) != np.sign(f[:-1]))[0]
    return s[i] - f[i]*(s[i + 1] - s[i])/(f[i + 1] - f[i])
dz = 0
for alpha, theta in [(4.5, 0.3), (9.3, 0.7)]:
    s = np.linspace(0.5, 30, 3000)
    z1 = zeros_of(Fn(alpha, theta, s), s); z2 = zeros_of(Fq(alpha, theta, s), s)
    if len(z1) != len(z2): dz = np.inf; break
    dz = max(dz, float(np.max(np.abs(z1 - z2))))
check("E3 Bessel series for F against direct quadrature: same zeros up to s = 30", dz < 1e-3, f"largest zero shift {dz:.1e}")

# --- E4: the grid
grid = {}
ok4, worst_gap = True, 9
for alpha in [4.2, 5.5, 8, 12, 20, 35]:
    S = 3*alpha + 40
    s = np.linspace(0.5, S, int(S*25))
    th = np.linspace(0, 1, 21)
    Z = [zeros_of(Fn(alpha, t, s), s) for t in th]
    n = sorted(set(len(z) for z in Z))
    if len(n) != 1:
        ok4 = False; grid[str(alpha)] = dict(counts=[len(z) for z in Z]); continue
    Z = np.array(Z)
    mono = bool((np.diff(Z, axis=0) < 0).all()); gap = float(np.min(np.diff(Z, axis=1)))
    ok4 = ok4 and mono and gap > 3.1; worst_gap = min(worst_gap, gap)
    grid[str(alpha)] = dict(count=int(n[0]), monotone=mono, min_gap=gap, first=[float(Z[0, 0]), float(Z[-1, 0])],
                            zeros=Z[:, :12].tolist())
    print(f"   alpha {alpha}: {n[0]} zeros, all move down: {mono}, closest pair {gap:.3f}", flush=True)
check("E4 zeros of F(alpha, theta) slide down as theta goes 0 -> 1, count fixed, never closer than 3.1", ok4, f"closest pair anywhere {worst_gap:.3f}")
OUT["E4"] = grid

# --- E5: p = 6
r = np.roots([1, 0, -1, 0, 1])          # 1 - z^2 + z^4 (from 1 + x^2 + x^4 with x -> -i z)
nonreal = int(np.sum(np.abs(r.imag) > 1e-9))
check("E5 at p = 6 the factor 1 + x^2 + x^4 gives 1 - D^2 + D^4, whose roots are not real", nonreal == 4, f"{nonreal} of 4 roots non-real")
OUT["E5_roots"] = [[float(z.real), float(z.imag)] for z in r]

OUT["elapsed_s"] = time.perf_counter() - t0
json.dump(OUT, open("even_p_check.json", "w"), indent=1)
print(f"{sum(RES)}/{len(RES)} pass, {OUT['elapsed_s']:.0f} s")
