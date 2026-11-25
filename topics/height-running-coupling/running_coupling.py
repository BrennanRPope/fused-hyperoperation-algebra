"""Height flow of E(x) = e^x - 1 read as a running coupling.

Computes (1) the exact beta function v(x) (Julia's iterative logarithm) from the Julia
equation v(E(x)) = E'(x) v(x); (2) its scheme-free invariants b0, rho = Res dx/v;
(3) the two-loop running of x_n = E^n(-1) against one- and two-loop formulas;
(4) the Borel transforms of v and of the formal half-iterate, with ratio tests and Pade
poles, to locate the Borel singularities; (5) the QFT comparison numbers.
Writes running_coupling.json and prints a pass/fail log.
"""
from fractions import Fraction as Fr
import math, json
import mpmath as mp

mp.mp.dps = 260
N = 150          # number of series coefficients
OUT = {}
checks = []
def check(name, ok, detail=""):
    checks.append((name, bool(ok), detail)); print(("PASS " if ok else "FAIL ") + name + ("  " + detail if detail else ""))

# ---------- 1. exact v(x) = sum c_n x^n from the Julia equation ----------
def julia_coeffs(Nmax, num):
    one = num(1)
    E = [num(0)] + [one / math.factorial(i) if num is not Fr else Fr(1, math.factorial(i)) for i in range(1, Nmax + 4)]
    L = Nmax + 2
    # powers of E truncated at degree L
    pows = [[num(0)] * (L + 1) for _ in range(Nmax + 1)]
    pows[0][0] = one
    for k in range(1, Nmax + 1):
        prev = pows[k - 1]
        cur = pows[k]
        for p in range(k - 1, L + 1):
            if prev[p] == 0: continue
            for q in range(1, L + 1 - p):
                cur[p + q] += prev[p] * E[q]
    Ep = [E[i + 1] * (i + 1) for i in range(L)]  # E'(x)
    c = [num(0)] * (Nmax + 1)
    c[2] = one / 2 if num is not Fr else Fr(1, 2)
    for n in range(3, Nmax + 1):
        o = n + 1  # order of x where c_n first appears with nonzero factor (n/2 - 1)
        rest = num(0)
        for k in range(2, n):
            rest += c[k] * pows[k][o]
            rest -= c[k] * Ep[o - k]
        c[n] = -rest / (num(n) / 2 - 1) if num is not Fr else -rest / (Fr(n, 2) - 1)
    return c

cex = julia_coeffs(24, Fr)
OUT["v_exact"] = {n: str(cex[n]) for n in range(2, 25)}
known = {2: Fr(1, 2), 3: Fr(-1, 12), 4: Fr(1, 48), 5: Fr(-1, 180), 6: Fr(11, 8640)}
check("v(x) = x^2/2 - x^3/12 + x^4/48 - x^5/180 + 11x^6/8640 (matches earlier task)", all(cex[k] == v for k, v in known.items()))
c = julia_coeffs(N, mp.mpf)
check("high-precision coefficients agree with exact ones", all(abs(c[n] - mp.mpf(cex[n].numerator) / cex[n].denominator) < mp.mpf(10) ** -200 for n in range(2, 25)))
# verify Julia equation on a numeric point with truncated series (asymptotic, tiny x)
x0 = mp.mpf("-0.01")
vt = lambda x: sum(c[n] * x ** n for n in range(2, 40))
res = vt(mp.e ** x0 - 1) - mp.e ** x0 * vt(x0)
check("Julia equation v(E(x)) = E'(x) v(x) holds at x = -0.01", abs(res) < mp.mpf(10) ** -60, mp.nstr(res, 3))

# ---------- 2. invariants ----------
b0, b1 = c[2], c[3]
rho = -b1 / b0 ** 2          # = Res_{x=0} dx / v(x)
check("rho = Res dx/v = 1/3", abs(rho - mp.mpf(1) / 3) < mp.mpf(10) ** -100)
# invariance under a random 'scheme change' x = y + a y^2 + b y^3: new v~(y) = v(x(y)) / x'(y)
a, b = mp.mpf("0.7"), mp.mpf("-1.3")
def series_mul(p, q, L):
    r = [mp.mpf(0)] * (L + 1)
    for i, pi in enumerate(p):
        if pi == 0: continue
        for j, qj in enumerate(q[:L + 1 - i]): r[i + j] += pi * qj
    return r
L = 8
X = [0, 1, a, b] + [0] * (L - 3); X = [mp.mpf(t) for t in X]
vX = [mp.mpf(0)] * (L + 1); P = [mp.mpf(1)] + [mp.mpf(0)] * L
for n in range(1, L + 1):
    P = series_mul(P, X, L)
    if n >= 2:
        vX = [vX[i] + c[n] * P[i] for i in range(L + 1)]
dX = [X[i + 1] * (i + 1) for i in range(L)] + [mp.mpf(0)]
# divide vX by dX
inv = [mp.mpf(0)] * (L + 1); inv[0] = 1 / dX[0]
for n in range(1, L + 1):
    inv[n] = -sum(dX[k] * inv[n - k] for k in range(1, n + 1)) / dX[0]
vt2 = series_mul(vX, inv, L)
rho2 = -vt2[3] / vt2[2] ** 2
check("b0 and rho unchanged by a scheme change x = y + 0.7y^2 - 1.3y^3; b2 changes",
      abs(vt2[2] - b0) < 1e-50 and abs(rho2 - rho) < 1e-50 and abs(vt2[4] - c[4]) > 1e-3,
      f"b2: {mp.nstr(c[4],6)} -> {mp.nstr(vt2[4],6)}")
OUT["b0"] = 0.5; OUT["b1"] = -1 / 12; OUT["rho"] = 1 / 3

# ---------- 3. running of x_n = E^n(-1): one loop vs two loop ----------
mp.mp.dps = 50
xs = []; x = mp.mpf(-1)
for n in range(0, 200001):
    if n in (1, 2, 3, 5, 10, 20, 50, 100, 200, 500, 1000, 2000, 5000, 10000, 20000, 50000, 100000, 200000): xs.append((n, x))
    x = mp.e ** x - 1
# Fatou coordinate in z = -2/x: Phi(z) = z - rho log z + phi(z); so z_n - rho log z_n = n + C
# fit C from the largest n, then predict
def z_of(xv): return -2 / xv
nbig, xbig = xs[-1]
zb = z_of(xbig)
# include the first asymptotic correction of phi(z): Phi'(z) = 1/V(z), V = 1 + d1/z + d2/z^2 ...
d = [mp.mpf(1)] + [2 * c[k + 2] * (-2) ** k for k in range(1, 30)]
# 1/V series
iv = [mp.mpf(0)] * 30; iv[0] = 1
for n in range(1, 30): iv[n] = -sum(d[k] * iv[n - k] for k in range(1, n + 1))
def Phi(z, K=12):  # asymptotic Fatou coordinate
    s = z - (-iv[1]) * 0  # placeholder
    s = z + iv[1] * mp.log(z)
    for k in range(2, K): s += iv[k] * z ** (1 - k) / (1 - k)
    return s
check("Fatou coordinate log coefficient is -rho = -1/3 (Phi = z - (1/3) log z + ...)", abs(iv[1] + mp.mpf(1) / 3) < 1e-40)
C = Phi(zb) - nbig
run = []
for n, xv in xs:
    one_loop = -2 / mp.mpf(n) if n > 0 else None
    # two-loop closed form: invert z - rho log z = n + C to leading order
    z2 = mp.findroot(lambda z: z - mp.log(z) / 3 - (n + C), n + C + 1)
    two_loop = -2 / z2
    run.append({"n": n, "x": float(xv), "one_loop": float(one_loop), "two_loop_with_constant": float(two_loop)})
OUT["running"] = run
r200 = run[[r["n"] for r in run].index(200)]
check("two-loop formula fits x_200 far better than one-loop -2/n",
      abs(r200["two_loop_with_constant"] - r200["x"]) < 1e-3 * abs(r200["one_loop"] - r200["x"]),
      f"one-loop err {abs(r200['one_loop']-r200['x']):.2e}, two-loop err {abs(r200['two_loop_with_constant']-r200['x']):.2e}")
# the log term: n*(n x_n + 2) / log n -> 2*rho*2? x = -2/z, z = n + (1/3)log n + C' -> n x + 2 = (2/3 log n + 2C')/n
mp.mp.dps = 260

# ---------- 4. Borel transforms ----------
# (a) generator in z: V(z) = 1 + sum_{k>=1} d_k z^-k
dz = [2 * c[k + 2] * (-2) ** k for k in range(0, N - 1)]   # dz[0] = 1
check("V(z) = 1 + (1/3)/z + ...  (map is z -> z + 1 + (1/3)/z + ...)", abs(dz[1] - mp.mpf(1) / 3) < 1e-100)
Bv = [dz[k] / mp.factorial(k - 1) for k in range(1, len(dz))]   # Borel coefficients of zeta^(k-1)
# radius: since the nearest singularities are a conjugate pair +-i*R, take the root of |b_k b_{k+2}|... use limsup |b_k|^{1/k}
def radius_est(B, lo, hi):
    vals = []
    for k in range(lo, hi):
        if B[k] != 0: vals.append((k, mp.mpf(1) / abs(B[k]) ** (mp.mpf(1) / (k + 1))))
    return vals
# better: fit b_k ~ A k^g R^-k cos(pi k/2 + phi). Use pairs of same parity: b_{k+4}/b_k -> R^-4 (k+4/k)^g
def R_from_parity(B, k):
    return (abs(B[k]) / abs(B[k + 4])) ** (mp.mpf(1) / 4)
Rv = [R_from_parity(Bv, k) for k in range(100, 140)]
Rv_last = Rv[-1]
OUT["borel_radius_generator"] = [float(r) for r in Rv[::5]]
check("Borel transform of the generator: nearest singularity at distance 2*pi (ratio test, k~140)",
      abs(Rv_last - 2 * mp.pi) < 0.1, f"R = {mp.nstr(Rv_last, 8)} vs 2pi = {mp.nstr(2*mp.pi, 8)}")
# sign pattern of period 4 means the singularities sit on the imaginary axis (+-2 pi i)
signs = "".join("+" if Bv[k] > 0 else "-" for k in range(0, 24))
OUT["borel_sign_pattern"] = signs
# Fit the large-order shape: if the nearest singularities are (1 - zeta/(+-2 pi i))^(gamma), gamma = -1 + beta' + i sigma,
# then b_j (2 pi)^j = Re[C i^-j j^(beta - i sigma)]. Theory (e^{2 pi i Phi}, Phi = z - rho log z) predicts sigma = 2 pi rho = 2pi/3.
import numpy as np
def fit_shape(B, lo, hi, corr=3):
    """fit u_j = Re[ i^-j j^(beta - i sigma) (C0 + C1/j + C2/j^2 ...) ], complex C's; scan beta, sigma"""
    js = np.arange(lo, hi).astype(float)
    u = np.array([float(B[int(j)] * (2 * mp.pi) ** int(j)) for j in js])
    def resid(beta, sig):
        ph = js * (-np.pi / 2) - sig * np.log(js)
        cols = []
        for m in range(corr):
            w = js ** (beta - m)
            cols += [w * np.cos(ph), -w * np.sin(ph)]
        A = np.vstack(cols).T
        coef, *_ = np.linalg.lstsq(A, u, rcond=None)
        return np.linalg.norm(A @ coef - u) / np.linalg.norm(u), coef
    best = None
    for beta in np.linspace(-3, 3, 121):
        for sig in np.linspace(-4, 4, 401):
            r, coef = resid(beta, sig)
            if best is None or r < best[0]: best = (r, beta, sig, coef)
    # refine
    r0, b0_, s0_, _ = best
    for step in (0.01, 0.002, 0.0004):
        for beta in np.arange(b0_ - 5 * step, b0_ + 5.5 * step, step):
            for sig in np.arange(s0_ - 5 * step, s0_ + 5.5 * step, step):
                r, coef = resid(beta, sig)
                if r < best[0]: best = (r, beta, sig, coef)
        _, b0_, s0_, _ = best
    return best
fr, fbeta, fsig, fC = fit_shape(Bv, 60, len(Bv))
OUT["generator_shape_fit"] = {"rel_residual": fr, "beta": fbeta, "sigma": fsig, "two_pi_rho": 2 * math.pi / 3}
check("large-order oscillation drifts like j^(-i sigma) with sigma = 2 pi rho = 2pi/3 (the log term of the Fatou coordinate)",
      abs(abs(fsig) - 2 * math.pi / 3) < 0.01 and fr < 1e-5, f"sigma = {fsig:.4f} (2pi/3 = {2*math.pi/3:.4f}), beta = {fbeta:.3f}, residual {fr:.1e}")
# Pade poles of the Borel transform
mp.mp.dps = 200
M = 50
coeffs = Bv[:2 * M + 1]
p, q = mp.pade(coeffs, M, M)
roots = mp.polyroots(q[::-1], maxsteps=400, extraprec=600)
poles = [complex(r) for r in roots]
OUT["pade_poles_generator"] = [[pp.real, pp.imag] for pp in poles]
near = sorted(poles, key=abs)[:4]
check("Pade poles nearest the origin sit at +-2 pi i", all(abs(abs(pp.imag) - 2 * math.pi) < 0.05 and abs(pp.real) < 0.05 for pp in near[:2]),
      str([f"{pp.real:+.4f}{pp.imag:+.4f}i" for pp in near[:2]]))
OUT["pade_nearest"] = [[pp.real, pp.imag] for pp in near]

# (b) the formal half-iterate E^{1/2}: E^s = exp(s v d/dx) x, coefficients a_n(s)
mp.mp.dps = 260
s = mp.mpf(1) / 2
Ns = N
def D(f):  # v * f'
    fp = [f[i + 1] * (i + 1) for i in range(Ns)] + [mp.mpf(0)]
    r = [mp.mpf(0)] * (Ns + 1)
    for i in range(2, Ns + 1):
        if c[i] == 0: continue
        ci = c[i]
        for j in range(0, Ns + 1 - i): r[i + j] += ci * fp[j]
    return r
f = [mp.mpf(0)] * (Ns + 1); f[1] = mp.mpf(1)
half = list(f); term = list(f)
for j in range(1, Ns):
    term = D(term); term = [t * s / j for t in term]
    if all(t == 0 for t in term): break
    half = [half[i] + term[i] for i in range(Ns + 1)]
OUT["half_iterate_first"] = [mp.nstr(half[i], 12) for i in range(1, 9)]
# check: half composed with itself = E (to order 12)
def compose(fa, fb, L):
    out = [mp.mpf(0)] * (L + 1); P = [mp.mpf(1)] + [mp.mpf(0)] * L
    for n in range(1, L + 1):
        P = series_mul(P, fb[:L + 1], L)
        out = [out[i] + fa[n] * P[i] for i in range(L + 1)]
    return out
hh = compose(half, half, 14)
check("formal E^(1/2) composed with itself gives e^x - 1 (to x^14)",
      all(abs(hh[i] - mp.mpf(1) / mp.factorial(i)) < mp.mpf(10) ** -150 for i in range(1, 15)))
# Borel radius of the half iterate in x: coefficients grow like n! / (4 pi)^n (x = -2/z, R = 2pi)
# Borel coefficients in z = -2/x: beta_n = a_n (-2)^n / (n-1)!
Bh = [half[n] * (-2) ** n / mp.factorial(n - 1) for n in range(2, Ns + 1)]   # Bh[j] ~ coeff of zeta^(j+1)
Bh2 = [mp.mpf(0)] + Bh  # Bh2[j] = coeff of zeta^j
hr, hbeta, hsig, hC = fit_shape(Bh2, 60, len(Bh2) - 1)
OUT["half_iterate_shape_fit"] = {"rel_residual": hr, "beta": hbeta, "sigma": hsig}
check("formal half-iterate: same Borel lattice (+-2 pi i) and same drift sigma = 2pi/3",
      abs(abs(hsig) - 2 * math.pi / 3) < 0.01 and hr < 1e-5, f"sigma = {hsig:.4f}, beta = {hbeta:.3f}, residual {hr:.1e}")
OUT["growth_generator"] = [[k, float(mp.log10(abs(Bv[k]) * (2 * mp.pi) ** (k + 1)))] for k in range(1, len(Bv)) if Bv[k] != 0]
OUT["abs_c"] = [[n, float(mp.log10(abs(c[n])))] for n in range(2, N + 1) if c[n] != 0]

# ---------- 5. QFT comparison numbers (t = ln mu) ----------
pi = math.pi
qft = {}
qft["QED (alpha)"] = (2 / (3 * pi), 1 / (2 * pi ** 2))
for nf in (3, 5, 6):
    B0, B1 = 11 - 2 * nf / 3, 102 - 38 * nf / 3
    qft[f"QCD nf={nf} (alpha_s)"] = (-B0 / (2 * pi), -B1 / (8 * pi ** 2))
k16 = 16 * pi ** 2
qft["phi^4 single scalar (lambda/4!)"] = (3 / k16, -17 / 3 / k16 ** 2)
qft["E = e^x - 1 (height)"] = (0.5, -1 / 12)
rows = []
for name, (B0_, B1_) in qft.items():
    rows.append({"theory": name, "b0": B0_, "b1": B1_, "same_sign": (B0_ > 0) == (B1_ > 0),
                 "rho_times_sign_b0": -B1_ / B0_ ** 2 * (1 if B0_ > 0 else -1)})
OUT["qft_table"] = rows
for r in rows: print(f"  {r['theory']:34s} b0={r['b0']:+.5f} b1={r['b1']:+.6f}  same sign: {r['same_sign']}")
check("E shares phi^4's two-loop sign pattern (b0, b1 opposite), not QED's or QCD's",
      [r for r in rows if r['theory'].startswith('E')][0]['same_sign'] is False and
      [r for r in rows if r['theory'].startswith('phi')][0]['same_sign'] is False and
      all(r['same_sign'] for r in rows if r['theory'].startswith(('QED', 'QCD'))))

OUT["checks"] = [{"name": n, "pass": ok, "detail": d_} for n, ok, d_ in checks]
json.dump(OUT, open("/mnt/project-files/math-theory/height-running-coupling/running_coupling.json", "w"), indent=1, default=str)
print(f"\n{sum(ok for _, ok, _ in checks)} of {len(checks)} pass")
