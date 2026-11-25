"""Numerical analysis through the circle operation: every claim on the page, checked in real float64.
a ∘₁ b = a + b + ab = (1+a)(1+b) - 1.  Exact references use Fractions (float(Fraction) rounds
correctly, so it is the true correctly rounded answer) or mpmath at 60 digits.
Writes figure data to page_data.json."""
import math, random, json
from fractions import Fraction as F
from mpmath import mp, mpf
mp.dps = 60
EPS = 2.0**-53          # unit roundoff u for float64
OK = [0, 0]
def check(name, cond, extra=""):
    OK[1] += 1; OK[0] += bool(cond)
    print(("PASS " if cond else "FAIL ") + name + (("  " + extra) if extra else ""))
def ex_ma(a, b): return F(a) + F(b) + F(a) * F(b)
def rel(x, ex): return float(abs((F(x) - ex) / ex))
def digits(r): return 16.0 if r == 0 else max(0.0, min(16.0, -math.log10(r)))
random.seed(7)
data = {}

print("== 1. The lift loses digits; the fused form keeps them")
sizes = [10.0**(-k/2) for k in range(1, 33)]
curve = []
for e in sizes:
    L, Fu, Fm = [], [], []
    for _ in range(400):
        a = random.uniform(.5, 1) * e; b = random.uniform(.5, 1) * e * random.choice([-1, 1])
        ex = ex_ma(a, b)
        if ex == 0: continue
        L.append(rel((1 + a) * (1 + b) - 1, ex)); Fu.append(rel(a + b + a * b, ex)); Fm.append(rel(math.fma(a, b, a + b), ex))
    med = lambda v: sorted(v)[len(v)//2]
    curve.append([e, digits(med(L)), digits(med(Fu)), digits(med(Fm)), med(L), max(Fu), max(Fm)])
data["digits"] = [[c[0], round(c[1], 2), round(c[2], 2), round(c[3], 2)] for c in curve]
small = [c for c in curve if c[0] <= 1e-4]
check("for |a|,|b| <= 1e-4, a+b+ab is within 2.5e-16 relative (about one rounding)", all(c[5] < 2.5e-16 for c in small))
check("for |a|,|b| <= 1e-4, fma(a,b,a+b) is within 2.5e-16 relative", all(c[6] < 2.5e-16 for c in small))
ratio = [c[4] * c[0] / EPS for c in curve if 1e-14 < c[0] < 1e-2]
check("the lift (1+a)(1+b)-1 loses about log10(1/|a|) digits: median error x |a| / u stays within [0.05, 5]",
      all(0.05 < r < 5 for r in ratio), "range %.2f..%.2f" % (min(ratio), max(ratio)))
c12 = [c for c in curve if abs(c[0] - 1e-12) / 1e-12 < 1e-6][0]
check("at |a| ~ 1e-12 the lift keeps only about 4 of 16 digits", 2.5 < c12[1] < 5.5, "%.2f digits" % c12[1])

wb = 0
for e in [1e-1, 1e-4, 1e-8, 1e-12]:
    for _ in range(400):
        a = random.uniform(.5, 1) * e * random.choice([-1, 1]); b = random.uniform(.5, 1) * e * random.choice([-1, 1])
        ex = ex_ma(a, b)
        if ex != 0: wb = max(wb, rel(math.fma(b, 1 + a, a), ex))
check("the fma shape with 1 + a inside, fma(b, 1+a, a), rounds a lift first and loses a little (worst between 1e-14 and 1e-11)", 1e-14 < wb < 1e-11, "%.1e" % wb)
check("float spacing: about 2e-28 near 1e-12, 2.2e-16 near 1", 1e-28 < math.ulp(1e-12) < 3e-28 and math.ulp(1.0) == 2.0**-52)
print("== 2. Near the inverse b = -a/(1+a): fma is exact inside the octave 1/2 <= 1+a <= 2")
band = []; inside_ok = 0; inside_n = 0; exact_sum = 0; outside_bad = 0; outside_n = 0; worst_fused_in = 0
for i in range(1200):
    a = random.uniform(-0.95, 6.0)
    b = -a / (1 + a); b = b * (1 + random.uniform(-1e-9, 1e-9))
    ex = ex_ma(a, b)
    if ex == 0: continue
    r_fma = math.fma(a, b, a + b); r_fu = a + b + a * b
    cr = (r_fma == float(ex))
    if -0.5 <= a <= 1.0:
        inside_n += 1; inside_ok += cr; exact_sum += (F(a + b) == F(a) + F(b)); worst_fused_in = max(worst_fused_in, rel(r_fu, ex))
    else:
        outside_n += 1; outside_bad += (not cr)
    band.append([round(a, 4), round(digits(rel(r_fu, ex)), 2), round(digits(rel(r_fma, ex)), 2)])
band.sort(); data["band"] = band[::4]
check("inside the octave, a+b is computed exactly (Sterbenz)", exact_sum == inside_n, "%d/%d" % (exact_sum, inside_n))
check("inside the octave, fma(a,b,a+b) is correctly rounded every time", inside_ok == inside_n, "%d/%d" % (inside_ok, inside_n))
check("inside the octave, plain a+b+ab can lose many digits (worst relative error > 1e-6)", worst_fused_in > 1e-6, "%.1e" % worst_fused_in)
check("outside the octave, fma is not guaranteed (some results not correctly rounded)", outside_bad > 0, "%d/%d" % (outside_bad, outside_n))
# exact algebra behind the band: for 0<=a<=1, -2a <= inv <= -a/2; for -1/2<=a<0, |a|/2 <= inv <= 2|a|
okb = True
for a in [F(k, 64) for k in range(-32, 65) if k != 0]:
    inv = -a / (1 + a)
    okb &= (min(-2*a, -a/2) <= inv <= max(-2*a, -a/2))
check("exact: whenever 1/2 <= 1+a <= 2, the inverse lies between -2a and -a/2 (Sterbenz's window)", okb)

print("== 3. log1p and expm1: the chart that turns the circle operation into adding, tangent to the identity at 0")
okl = True
for _ in range(50):
    a, b = mpf(random.uniform(-.9, 3)), mpf(random.uniform(-.9, 3))
    okl &= abs(mp.log1p(a + b + a*b) - mp.log1p(a) - mp.log1p(b)) < mpf(10)**-50
check("log1p(a ∘₁ b) = log1p(a) + log1p(b) (60 digits, 50 random pairs)", okl)
okt = True
for _ in range(50):
    a, s, t = mpf(random.uniform(-.9, 3)), mpf(random.uniform(-2, 2)), mpf(random.uniform(-2, 2))
    P = lambda t: mp.expm1(t * mp.log1p(a))
    okt &= abs(P(s) + P(t) + P(s)*P(t) - P(s + t)) < mpf(10)**-45
check("fractional circle powers a^(∘₁t) = expm1(t·log1p a) add their exponents", okt)
x = 1e-10
naive_log = math.log(1 + x); good = math.log1p(x); exl = mp.log1p(mpf(x))
check("log(1+x) at x = 1e-10 loses about 7 digits; log1p keeps them",
      abs(naive_log - exl)/exl > 1e-8 and abs(good - exl)/exl < 2*EPS, "%.1e vs %.1e" % (abs(naive_log-exl)/exl, abs(good-exl)/exl))
def kahan_log1p(x):
    w = 1.0 + x
    return x if w == 1.0 else x * math.log(w) / (w - 1.0)
wk = max(float(abs((mpf(kahan_log1p(x)) - mp.log1p(mpf(x))) / mp.log1p(mpf(x)))) for x in [random.uniform(-.5, .5) * 10.0**-random.randint(0, 14) for _ in range(3000)])
check("Kahan's trick log1p(x) = x·log(w)/(w-1), w = fl(1+x), stays within 5u (Goldberg 1991, Thm 4)", wk < 5 * EPS, "worst %.2f u" % (wk / EPS))
okz = all(abs(float(mp.diff(f, 0)) - 1) < 1e-12 and abs(float(f(0))) < 1e-30 for f in [mp.log1p, mp.expm1])
check("log1p and expm1 both fix 0 with slope 1, while the lift x -> 1+x moves 0 to 1", okz)

print("== 4. Compound interest: naive (1+r/n)^n - 1 against expm1(n·log1p(r/n))")
r = 0.05; ci = []
for n in [12, 365, 8760, 525600, 31536000, 10**9]:
    exv = mp.expm1(n * mp.log1p(mpf(r) / n))
    nv = (1 + r/n)**n - 1; gv = math.expm1(n * math.log1p(r/n))
    en = float(abs((nv - exv) / exv)); eg = float(abs((gv - exv) / exv))
    ci.append([n, en, eg])
    print("   n=%-10d naive rel err %.1e   log1p/expm1 rel err %.1e" % (n, en, eg))
data["interest"] = ci
check("log1p/expm1 compound interest is within 4u for every n up to 1e9", all(c[2] < 4 * EPS for c in ci))
check("naive compound interest error grows with n (above 1e-9 by n = 1e9)", ci[-1][1] > 1e-9 and ci[-1][1] > 100 * ci[1][1])

print("== 5. Long products of small returns")
N = 200000
rs = [random.gauss(0, 1e-7) for _ in range(N)]
exP = mpf(1)
for v in rs: exP *= (1 + mpf(v))
exP -= 1
p = 1.0
for v in rs: p *= (1 + v)
naive = p - 1
acc = 0.0
for v in rs: acc = math.fma(acc, v, acc + v)
fold_fma = acc
via_log = math.expm1(math.fsum(math.log1p(v) for v in rs))
er = lambda x: float(abs((mpf(x) - exP) / exP))
print("   naive %.1e   fma fold %.1e   expm1(fsum log1p) %.1e" % (er(naive), er(fold_fma), er(via_log)))
data["longprod"] = [er(naive), er(fold_fma), er(via_log)]
check("200,000 returns of size 1e-7: naive product-minus-1 loses over 6 digits", er(naive) > 1e-10)
check("the fused fold acc ∘₁ r keeps at least 12 digits", er(fold_fma) < 1e-12)
check("expm1(sum of log1p) keeps at least 14 digits", er(via_log) < 1e-14)

# any dial value k: x1 ∘ₖ x2 ∘ₖ ... = expm1(sum log1p(k x)) / k ; k = -1 is "at least one"
M = 100000; pr = [random.uniform(0.5, 1.5) * 1e-12 for _ in range(M)]
exU = 1 - mp.fprod([1 - mpf(q) for q in pr])
naiveU = 1.0
for q in pr: naiveU *= (1 - q)
naiveU = 1 - naiveU
goodU = -math.expm1(math.fsum(math.log1p(-q) for q in pr))
eu = lambda x: float(abs((mpf(x) - exU) / exU))
print("   at least one of 100,000 events of chance 1e-12: naive %.1e   dial fold via log1p %.1e" % (eu(naiveU), eu(goodU)))
data["union"] = [eu(naiveU), eu(goodU)]
check("k = -1 (at least one): naive 1 - prod(1-p) loses over 6 digits for tiny chances", eu(naiveU) > 1e-10)
check("k = -1 via -expm1(sum log1p(-p)) keeps at least 14 digits", eu(goodU) < 1e-14)
okd = True
for k in [-2.0, -1.0, -0.5, 0.7, 3.0]:
    xs = [random.uniform(-0.2, 0.2) for _ in range(30)]
    acc = mpf(0)
    for x in xs: acc = acc + x + k * acc * x
    via = math.expm1(math.fsum(math.log1p(k * x) for x in xs)) / k
    okd &= abs((mpf(via) - acc) / acc) < 1e-14
check("one routine for every dial k: x1 ∘ₖ ... ∘ₖ xn = expm1(sum log1p(k·x))/k (k = -2, -1, -0.5, 0.7, 3)", okd)

print("== 6. Rounding errors are circle-operation numbers (Higham's theta_n)")
okh = True; okexact = True; th = []
for n in [10, 100, 300, 1000]:
    xs = [math.exp(random.gauss(0, 0.05)) for _ in range(n)]
    p = xs[0]; pe = F(xs[0]); deltas = []
    for x in xs[1:]:
        q = p * x; d = F(q) / (F(p) * F(x)) - 1; deltas.append(d); p = q; pe *= F(x)
    theta = F(p) / pe - 1
    fold = F(0)
    for d in deltas: fold = fold + d + fold * d
    okexact &= (fold == theta)
    m = n - 1
    bound = (1 + F(EPS))**m - 1; gam = m * F(EPS) / (1 - m * F(EPS))
    okh &= abs(theta) <= bound <= gam
    th.append([n, float(abs(theta)), float(bound)])
data["theta"] = th
check("the total relative error of a float product equals the exact circle-operation fold of the step errors", okexact)
check("|theta_n| <= (1+u)^n - 1 <= n·u/(1 - n·u) for n = 10..1000", okh)
x, y, d1, d2 = F(3, 7), F(-5, 11), F(1, 1000), F(-3, 2000)
check("product of perturbed values: relative error d1 ∘₁ d2 (exact)", (x*(1+d1))*(y*(1+d2)) == x*y*(1 + d1 + d2 + d1*d2))
check("quotient of perturbed values: relative error (d1 - d2)/(1 + d2), circle-operation subtraction (exact)",
      (x*(1+d1))/(y*(1+d2)) == (x/y)*(1 + (d1 - d2)/(1 + d2)))

print("== 7. Conditioning: two danger points, the identity (way back) and the wall (way in)")
okc = True
for a in [1e-3, 1e-6, 1e-9, 1e-12]:
    w = 1.0 + a; eta = F(w) / (1 + F(a)) - 1
    okc &= ((F(w) - 1 - F(a)) / F(a) == eta * (1 + F(a)) / F(a))
check("exact: the error of fl(1+a)-1 is the lift's rounding times the amplifier (1+a)/a", okc)
okk = True
for u in [-0.5, -0.9, -0.99, -0.999]:
    inv = lambda t: -t / (1 + t); h = 1e-7 * abs(u)
    k_num = abs((inv(u + h) - inv(u)) / h * u / inv(u)); okk &= abs(k_num - 1/abs(1 + u)) / (1/abs(1 + u)) < 1e-4
check("relative condition of the inverse -u/(1+u) is 1/|1+u|: it blows up only at the wall", okk)
kl = [abs(u / ((1 + u) * math.log1p(u))) for u in [-0.9, -0.99, -0.999, -0.9999]]
check("relative condition of log1p near the wall grows without bound", all(kl[i+1] > 5 * kl[i] for i in range(3)), " ".join("%.0f" % k for k in kl))
check("exact: (a+h) ∘₁ b - a ∘₁ b = h(1+b), so a ∘₁ b has relative condition |a(1+b)/(a ∘₁ b)|",
      all((F(a)+h) + F(b) + (F(a)+h)*F(b) - ex_ma(a, b) == h*(1+F(b)) for a, b, h in [(0.3, -0.2, F(1, 99)), (2.5, 7.0, F(-3, 7))]))

print("== 8. Newton's method for the inverse: v <- v ∘₁ (-(u ∘₁ v)); the residual squares")
okn = True
for u, v0 in [(F(3), F(-1, 2)), (F(-1, 3), F(1, 3)), (F(9), F(-4, 5))]:
    v = v0; e = u + v + u*v
    for _ in range(5):
        v = v + (-e) + v*(-e); e2 = u + v + u*v
        okn &= (e2 == -e*e); e = e2
check("exact: each Newton step maps the residual e = u ∘₁ v to -e² (three starts, five steps)", okn)
u = 3.0; v = -0.6; dg = []
for _ in range(7):
    e = math.fma(u, v, u + v); v = math.fma(v, -e, v - e); dg.append(digits(abs((F(v) - F(-3, 4)) / F(-3, 4))))
check("in float64 the correct digits roughly double each step and reach full precision", dg[-1] >= 15.5 and all(dg[i+1] > 1.6 * dg[i] for i in range(2, 5)), " ".join("%.1f" % d for d in dg))
def newton_runs(u, v0, steps=60):
    v = v0
    for _ in range(steps):
        e = u + v + u*v
        if abs(e) > 1e6: return False
        v = v - e - v*e
    return abs(u + v + u*v) < 1e-12
okb2 = all(newton_runs(3.0, (lift / 4.0) - 1) for lift in [0.05, 0.5, 1.5, 1.95]) and \
       not any(newton_runs(3.0, (lift / 4.0) - 1) for lift in [-0.3, 2.05, 3.0])
check("the basin is |u ∘₁ v0| < 1, i.e. 0 < (1+u)(1+v0) < 2: starts inside converge, starts outside fail", okb2)
okg = True
for e in [F(1, 3), F(-2, 5), F(7, 9)]:
    acc = -e
    for k in range(6):
        acc = acc + e**(2**k) + acc * e**(2**k)
        okg &= (acc == -e**(2**(k+1)))
check("exact: Goldschmidt's factors fold to (-e) ∘₁ e ∘₁ e² ∘₁ e⁴ ∘₁ ... = -e^(2^n)", okg)
gw = 0
for _ in range(2000):
    dd = random.uniform(0.5, 1.0); e = 1 - dd; q = 1.0; den = dd; t = e
    for _ in range(6):
        q *= (1 + t); den *= (1 + t); t = t * t
    gw = max(gw, abs(q * dd - 1))
check("Goldschmidt division of 1 by d in [1/2,1] reaches about 1e-15 after 6 factors", gw < 5e-15, "%.1e" % gw)

print("== 9. Horner and the cube")
okh2 = True
for n in [3, 6, 9]:
    bs = [F(random.randint(-9, 9), random.randint(10, 30)) for _ in range(n)]
    acc = F(0)
    for b in bs: acc = acc * (1 + b) + b
    corners = F(0)
    for mask in range(1, 2**n):
        p = F(1)
        for i in range(n):
            if mask >> i & 1: p *= bs[i]
        corners += p
    okh2 &= (acc == corners)
check("a chain of circle operations, each one Horner step acc·(1+b)+b, sums all 2^n - 1 corners of the cube in n steps", okh2)
check("the three fma shapes agree exactly: a·b+(a+b) = a·(1+b)+b = b·(1+a)+a",
      all(F(a)*F(b) + F(a)+F(b) == F(a)*(1+F(b)) + F(b) == F(b)*(1+F(a)) + F(a) for a, b in [(0.1, 0.2), (-0.7, 3.0)]))
okp = all(((1+F(a))**(2*n) - 1) == (lambda m: m + m + m*m)((1+F(a))**n - 1) for a in [0.03, -0.4] for n in [1, 5, 20])
check("Circle powers by squaring: a^(∘₁2n) = a^(∘₁n) ∘₁ a^(∘₁n)", okp)
okbern = all((1+F(a))**n - 1 >= n*F(a) for a in [F(-1), F(-1, 2), F(1, 10), F(3)] for n in range(0, 30))
check("Bernoulli: a^(∘₁n) >= n·a for a >= -1 (compounding never loses to simple interest)", okbern)

json.dump(data, open(__file__.replace("checks.py", "page_data.json"), "w"))
print("\n%d/%d checks pass" % tuple(OK))
