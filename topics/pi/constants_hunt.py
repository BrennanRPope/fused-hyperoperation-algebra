"""
Hunting other transcendental constants that the circle operation produces.
Fused Hyperoperation Algebra, Nov 25 2026.

For each candidate: compute to high precision, then run integer-relation searches (PSLQ)
against a basis of known constants, and an algebraic-number search (findpoly).
Writes constants_results.json and constants_output.txt.
"""
import json
import mpmath as mp

DPS = 130
mp.mp.dps = DPS
log = []
OUT = {}


def say(s):
    print(s)
    log.append(s)


# --------------------------------------------------------------------------------------------
# A. The half rung: regular half-iterate of E(x) = e^x - 1 (fused exponential), on x > 0.
#    E has a parabolic fixed point at 0. Abel (Fatou) coordinate: alpha(E(x)) = alpha(x) + 1,
#    alpha(x) = -2/x + (1/3) log x + sum_k c_k x^k  (asymptotic series; Szekeres 1958, Ecalle).
#    E^s(x) = alpha^{-1}(alpha(x) + s). Then B_s(a,b) = E^s(E^-s(a) + E^-s(b)) is the fused
#    operation at height s: B_0 = +, B_1 = the circle operation.
# --------------------------------------------------------------------------------------------
NC = 60    # series terms


def abel_coeffs(nc):
    """Solve alpha(E(x)) - alpha(x) - 1 = 0 order by order for c_1..c_nc (exact rationals via mpmath at high prec)."""
    # work with truncated power series as lists of mpf, length N
    N = nc + 4
    fact = [mp.factorial(k) for k in range(N + 3)]
    # E(x)/x = sum_{k>=0} x^k/(k+1)!
    r = [mp.mpf(1) / fact[k + 1] for k in range(N + 2)]

    def mul(a, b):
        c = [mp.mpf(0)] * N
        for i, ai in enumerate(a[:N]):
            if ai == 0:
                continue
            for j in range(N - i):
                c[i + j] += ai * b[j]
        return c

    def inv(a):
        b = [mp.mpf(0)] * N
        b[0] = 1 / a[0]
        for n in range(1, N):
            s = 0
            for k in range(1, n + 1):
                s += a[k] * b[n - k]
            b[n] = -s / a[0]
        return b

    def logser(a):  # log(a) for a[0] = 1
        # log(1+u) = sum (-1)^{m+1} u^m/m
        u = a[:]
        u[0] = 0
        res = [mp.mpf(0)] * N
        p = [mp.mpf(1)] + [mp.mpf(0)] * (N - 1)
        for m in range(1, N):
            p = mul(p, u)
            for i in range(N):
                res[i] += (-1) ** (m + 1) * p[i] / m
        return res

    rr = r[:N]
    # -2/E + 2/x = (2/x)(1 - 1/r): series in x
    one_over_r = inv(rr)
    t1 = [2 * ((1 if i == 0 else 0) - one_over_r[i]) for i in range(N)]
    # shift: (2/x)*(series) -> coefficient of x^(i-1)
    head = [t1[i + 1] for i in range(N - 1)] + [mp.mpf(0)]    # constant term is t1[1]
    lr = logser(rr)
    base = [head[i] + lr[i] / 3 for i in range(N)]
    base[0] -= 1                                               # the "-1"
    # now sum_k c_k (E^k - x^k) must cancel base. E^k = x^k r^k.
    c = [mp.mpf(0)] * (nc + 1)
    rpow = [None] * (nc + 2)
    rpow[0] = [mp.mpf(1)] + [mp.mpf(0)] * (N - 1)
    for k in range(1, nc + 2):
        rpow[k] = mul(rpow[k - 1], rr)
    # coefficient of x^n: base[n] + sum_{k=1}^{n} c_k [r^k]_{n-k} - c_n (for n with k=n term: r^n_0 = 1 -> cancels)
    # so the x^n equation involves c_1..c_{n-1} and c_n * ([r^n]_0 - 1) = 0 -> c_n enters at order n+1 via [r^n]_1 = n/2
    # equation at order n+1: base[n+1] + sum_{k=1}^{n+1} c_k [r^k]_{n+1-k} - c_{n+1} = 0
    #   = base[n+1] + sum_{k=1}^{n} c_k [r^k]_{n+1-k}  (k = n+1 term cancels) ; c_n appears with [r^n]_1 = n/2
    for n in range(1, nc + 1):
        s = base[n + 1]
        for k in range(1, n):
            s += c[k] * rpow[k][n + 1 - k]
        c[n] = -s / rpow[n][1]
    assert abs(base[0]) < mp.mpf(10) ** (-DPS + 10), base[0]
    return c


C = abel_coeffs(NC)


def alpha_small(x):
    return -2 / x + mp.log(x) / 3 + mp.fsum(C[k] * x ** k for k in range(1, NC + 1))


NIT = 20000


def alpha(x):
    """Abel coordinate for x > 0: pull x toward 0 with log(1+x), NIT times."""
    for _ in range(NIT):
        x = mp.log1p(x)
    return alpha_small(x) + NIT


def alpha_inv(A):
    """x > 0 with alpha(x) = A: solve near 0, then push out with e^x - 1."""
    target = A - NIT
    y = -2 / target                      # first guess from -2/y
    for _ in range(200):
        fy = alpha_small(y) - target
        dfy = 2 / y ** 2 + 1 / (3 * y) + mp.fsum(k * C[k] * y ** (k - 1) for k in range(1, NC + 1))
        step = fy / dfy
        y -= step
        if abs(step) < abs(y) * mp.mpf(10) ** (-DPS + 5):
            break
    for _ in range(NIT):
        y = mp.expm1(y)
    return y


def Es(s, x):
    return alpha_inv(alpha(x) + s)


# sanity: E^1 = e^x - 1, E^{1/2} twice = E
x0 = mp.mpf(1)
chk1 = Es(1, x0) - mp.expm1(1)
h = Es(mp.mpf(1) / 2, x0)
chk2 = Es(mp.mpf(1) / 2, h) - mp.expm1(1)
chk3 = Es(-1, x0) - mp.log(2)
say(f"half-rung machinery: E^1(1) - (e-1) = {mp.nstr(chk1, 3)}, E^1/2(E^1/2(1)) - (e-1) = {mp.nstr(chk2, 3)}, E^-1(1) - ln2 = {mp.nstr(chk3, 3)}")
good_digits = int(-mp.log10(max(abs(chk1), abs(chk2), abs(chk3), mp.mpf(10) ** -DPS))) - 5
say(f"  trusted digits ~ {good_digits}")

half1 = h
Bhalf = Es(mp.mpf(1) / 2, 2 * Es(-mp.mpf(1) / 2, x0))
say(f"E^(1/2)(1)      = {mp.nstr(half1, 70)}")
say(f"B_(1/2)(1,1)    = {mp.nstr(Bhalf, 70)}   (half-way between 1+1 = 2 and 1 (+) 1 = 3)")

# --------------------------------------------------------------------------------------------
# B. other candidates (known constants, used to show what the detector catches)
# --------------------------------------------------------------------------------------------
cands = {
    'E_half_at_1': half1,
    'B_half_1_1': Bhalf,
    'gluing_q_tau_i': mp.exp(-2 * mp.pi),                    # Tate q of y^2 = x^3 - x (j = 1728)
    'gluing_q_tau_rho': -mp.exp(-mp.pi * mp.sqrt(3)),         # j = 0
    'omega_fused_fixed': mp.lambertw(1).real,                 # fixed point of x -> (1/e)^x (fused rung 3, base 1/e)
    'kiss_rung3': mp.exp(1 / mp.e),                           # e^(1/e)
    'turning_period_rung_minus1': 2 * mp.pi / mp.sqrt(1),     # rung s=-1 on the rung curve: 2 pi / sqrt(-s)
}
OUT['values'] = {k: mp.nstr(v, 80) for k, v in cands.items()}

# Feigenbaum delta: on the road from the circle operation's doubling (z -> z^2, c = 0) to the turning law's doubling
# (z -> z^2 - 2, c = -2). Computed here from superstable parameters of z^2 + c (period 2^n).
mp.mp.dps = 40


def superstable(n, guess):
    # c where 0 is periodic with period 2^n: f^{2^n}(0) = 0, Newton in c
    c = guess
    for _ in range(100):
        z, dz = mp.mpf(0), mp.mpf(0)
        for _ in range(2 ** n):
            dz = 2 * z * dz + 1
            z = z * z + c
        step = z / dz
        c -= step
        if abs(step) < mp.mpf(10) ** -35:
            break
    return c


cs = [mp.mpf(0), mp.mpf(-1)]
deltas = []
for n in range(2, 15):
    d = (cs[-2] - cs[-3]) / (cs[-1] - cs[-2]) if len(cs) > 2 else mp.mpf('4.669')
    guess = cs[-1] + (cs[-1] - cs[-2]) / d
    cs.append(superstable(n, guess))
    if len(cs) > 3:
        deltas.append((cs[-2] - cs[-3]) / (cs[-1] - cs[-2]))
say(f"Feigenbaum road (c = 0 circle doubling -> c = -2 turning doubling): superstable c_12 = {mp.nstr(cs[-1], 15)}, delta estimate {mp.nstr(deltas[-1], 12)} (known 4.66920160910)")
OUT['feigenbaum'] = {'c_superstable': [mp.nstr(c, 20) for c in cs], 'delta_estimates': [mp.nstr(d, 15) for d in deltas]}
mp.mp.dps = DPS

# --------------------------------------------------------------------------------------------
# C. integer-relation searches
# --------------------------------------------------------------------------------------------
P = 110   # digits used in PSLQ (well inside the trusted digits)
mp.mp.dps = P + 10
basis_names = ['1', 'pi', 'e', 'ln2', 'ln3', 'euler_gamma', 'zeta3', 'catalan', 'pi^2', 'sqrt2', 'sqrt3', 'gamma(1/4)', 'ln_pi', 'pi*sqrt3']
basis_vals = [mp.mpf(1), +mp.pi, +mp.e, mp.log(2), mp.log(3), +mp.euler, mp.zeta(3), +mp.catalan, mp.pi ** 2,
              mp.sqrt(2), mp.sqrt(3), mp.gamma(mp.mpf(1) / 4), mp.log(mp.pi), mp.pi * mp.sqrt(3)]

results = {}
for name, val in cands.items():
    v = mp.mpf(mp.nstr(val, P + 8)) if mp.im(val) == 0 else val
    res = {}
    for label, Xs, names in [('X', [v], ['X']), ('log|X|', [mp.log(abs(v))], ['log|X|']),
                             ('X and log|X| together', [v, mp.log(abs(v))], ['X', 'log|X|'])]:
        rel = mp.pslq(Xs + basis_vals, maxcoeff=10 ** 4, maxsteps=10 ** 6)
        if rel:
            terms = [f"{c}*{n}" for c, n in zip(rel, names + basis_names) if c]
            res[label] = " + ".join(terms) + " = 0"
        else:
            res[label] = "no relation with coefficients up to 10^4"
    poly = mp.findpoly(v, 12, maxcoeff=10 ** 6)
    res['algebraic_deg<=12'] = str(poly) if poly else "not a root of any integer polynomial of degree <= 12, coefficients <= 10^6"
    results[name] = res
    say(f"\n{name} = {mp.nstr(v, 30)}")
    for k_, v_ in res.items():
        say(f"   {k_}: {v_}")
OUT['pslq'] = results
OUT['pslq_digits'] = P
OUT['basis'] = basis_names
mp.mp.dps = DPS
OUT['abel_coeffs_first'] = [mp.nstr(c, 25) for c in C[1:12]]
OUT['trusted_digits'] = good_digits

json.dump(OUT, open('/mnt/project-files/math-theory/pi/constants_results.json', 'w'), indent=1)
open('/mnt/project-files/math-theory/pi/constants_output.txt', 'w').write("\n".join(log) + "\n")
