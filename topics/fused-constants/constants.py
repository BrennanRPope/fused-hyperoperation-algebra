"""Natural constants up the hyperoperation ladder, ordinary and fused.

Fused rung n (same definition as the ns-rungs task):  a (n) b = H_n(1+a, 1+b) - 1.
Natural constant of rung n (tower test): the largest base whose infinite tower at the
next rung converges = the base where the curve x -> H_n(b, x) just kisses the diagonal.
Fused constant = ordinary constant - 1 (conjugating by x -> 1+x slides every curve and
the diagonal together, so kisses, slopes and fixed points all move by exactly 1).

Real heights need an extension.  Two are used:
  * linear  : H_k(b,x) = 1+x on (-1,0] for every rung k >= 4 (the classic linear approximation)
  * natural : Robbins' natural super-logarithm (agrees with Kneser's analytic tetration), rung 4 only
"""
import mpmath as mp, json
mp.mp.dps = 40

def HD(k, b, x):
    """(value, d/dx) of the rank-k hyperoperation at real height x > -1, linear extension."""
    if x > 200: return mp.inf, mp.inf
    if k == 1: return b + x, mp.mpf(1)
    if k == 2: return b * x, b
    if k == 3:
        v = b ** x; return v, v * mp.log(b)
    if x <= 0: return 1 + x, mp.mpf(1)
    v1, d1 = HD(k, b, x - 1)
    v2, d2 = HD(k - 1, b, v1)
    return v2, d2 * d1

def Hf(k, b, x):  # float version for scanning (b > 1); huge values capped at inf
    if x > 200: return float('inf')
    if k == 3: return b ** x
    if x <= 0: return 1 + x
    return Hf(k - 1, b, Hf(k, b, x - 1))

def kiss(F, b0, x0):
    """Solve F(b,x)=x and F'(b,x)=1 together.  F returns (value, d/dx)."""
    return mp.findroot([lambda b, x: F(b, x)[0] - x, lambda b, x: F(b, x)[1] - 1], (b0, x0))

out = []
def say(s=""): print(s); out.append(s)
data = {}

# ---------- rung 3 ----------
eta = mp.e ** (1 / mp.e)
b3, x3 = kiss(lambda b, x: HD(3, b, x), mp.mpf("1.44"), mp.mpf("2.7"))
say(f"rung 3  kiss base  {mp.nstr(b3, 30)}   at x = {mp.nstr(x3, 30)}")
say(f"        e^(1/e)    {mp.nstr(eta, 30)}   |diff| {mp.nstr(abs(b3 - eta), 3)};  |kiss point - e| = {mp.nstr(abs(x3 - mp.e), 3)}")
d = mp.findroot(lambda t: mp.diff(lambda u: u ** (1 / u), t), 2.7)
say(f"        cube side whose volume can equal its dimension: max of d^(1/d) at d = {mp.nstr(d, 25)}, side {mp.nstr(d ** (1 / d), 25)}")
for u in (mp.mpf("0.37"), mp.mpf("-0.6"), mp.mpf(2)):
    say(f"        stretch check u={mp.nstr(u,3)}: (e^(1/e))^(e(1+u))/e - 1 - (e^u - 1) = {mp.nstr(eta ** (mp.e * (1 + u)) / mp.e - 1 - (mp.e ** u - 1), 3)}")
W = mp.lambertw(1)
say(f"        tower lower edge e^-e = {mp.nstr(mp.e ** -mp.e, 20)};  omega W(1) = {mp.nstr(W, 20)}, (1/e)^omega - omega = {mp.nstr((1 / mp.e) ** W - W, 3)}")
data["rung3"] = {"c": str(b3), "x": str(x3)}

# ---------- rungs 4..: linear extension ----------
say()
say("Linear extension: tower-test constant c_k, kiss point x_k (fused: subtract 1 from both)")
res = {}
prev = None
for k in range(4, 21):
    xs = [1 + i / 500 for i in range(1, 2000)]
    lo, hi = 1.3, 2.0
    for _ in range(45):
        mid = (lo + hi) / 2
        if min(Hf(k, mid, x) - x for x in xs) <= 0: lo = mid
        else: hi = mid
    xg = min(xs, key=lambda x: Hf(k, lo, x) - x)
    b, x = kiss(lambda b, x: HD(k, b, x), mp.mpf(lo), mp.mpf(xg))
    # check: just below b the curve crosses, just above it clears the diagonal
    m_below = min(HD(k, b - mp.mpf("1e-12"), x + mp.mpf(t) / 1000)[0] - (x + mp.mpf(t) / 1000) for t in range(-50, 51))
    m_above = min(HD(k, b + mp.mpf("1e-12"), mp.mpf(xx))[0] - mp.mpf(xx) for xx in xs[::5] if xx <= 4)
    res[k] = (b, x)
    say(f"  k={k:2d}  c_k = {mp.nstr(b, 25)}   x_k = {mp.nstr(x, 18)}   2-c_k = {mp.nstr(2 - b, 8)}"
        + (f"  ratio {mp.nstr((2 - b) / (2 - prev), 6)}" if prev else "")
        + f"   [crosses below: {m_below < 0}, clears above: {m_above > 0}]")
    prev = b
data["linear"] = {k: {"c": str(b), "x": str(x)} for k, (b, x) in res.items()}
say(f"  ln 2 = {mp.nstr(mp.log(2), 8)}  (expected ratio if 2 - c_k shrinks like (ln c)^k near c = 2)")

# ---------- rung 4 again, natural (Robbins = Kneser) extension ----------
say()
def robbins(b, n):
    L = mp.log(b)
    A = mp.matrix(n, n); rhs = mp.matrix(n, 1)
    for j in range(n):
        for kk in range(1, n + 1):
            A[j, kk - 1] = (kk * L) ** j / mp.factorial(j) - (1 if kk == j else 0)
    rhs[0] = 1
    v = mp.lu_solve(A, rhs)
    slog = lambda z: -1 + sum(v[i] * z ** (i + 1) for i in range(n))
    dslog = lambda z: sum((i + 1) * v[i] * z ** i for i in range(n))
    def sexpD(x):
        if x <= 0:
            z = mp.findroot(lambda z: slog(z) - x, 1 + x)
            return z, 1 / dslog(z)
        v1, d1 = sexpD(x - 1)
        v2 = b ** v1
        return v2, v2 * L * d1
    return sexpD

nat = {}
for n in (20, 30, 40):
    bb, xx = kiss(lambda b, x: robbins(b, n)(x), mp.mpf("1.6353"), mp.mpf("3.09"))
    nat[n] = (bb, xx)
    say(f"natural extension, {n} terms: c_4 = {mp.nstr(bb, 18)}  kiss at x = {mp.nstr(xx, 18)}")
say("published (Tetration Forum 'enigma' base, Kneser): 1.63532449671527639934534, upper fixed point 3.08853227180671765")
data["natural4"] = {n: {"c": str(b), "x": str(x)} for n, (b, x) in nat.items()}

# ---------- bound by 2 ----------
say()
say("H_k(2,2), k=1..8: " + ", ".join(mp.nstr(HD(k, mp.mpf(2), mp.mpf(2))[0], 6) for k in range(1, 9)))
say("min over x in (0,4) of H_k(2,x) - x, k=4..12: " +
    ", ".join(f"{min(Hf(k, 2.0, i / 200) - i / 200 for i in range(1, 800)):.4f}" for k in range(4, 13)))

# ---------- collapse onto the step ----------
say()
for k in (3, 5, 8, 12, 20, 30):
    say(f"b=1.8 k={k:2d}: H at x = 0.5, 1.25, 1.5, 1.75, 2, 3: " + ", ".join(f"{Hf(k, 1.8, t):.6f}" for t in (0.5, 1.25, 1.5, 1.75, 2, 3)))

# ---------- slope-one-at-start test on each rung ----------
say()
for k in (4, 5, 6):
    for b in (mp.mpf(2), mp.e):
        say(f"slope test k={k} b={mp.nstr(b, 6)}: slope just left of 0 = 1, just right = {mp.nstr(HD(k, b, mp.mpf('1e-30'))[1], 12)}")
sn = robbins(mp.e, 30)
say(f"natural extension base e: slope of sexp at 0 = {mp.nstr(sn(mp.mpf(0))[1], 12)}  (not 1: the slope test is extension-dependent above rung 3)")

open("constants_output.txt", "w").write("\n".join(out) + "\n")
json.dump(data, open("constants.json", "w"), indent=1)
