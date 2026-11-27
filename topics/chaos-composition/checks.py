"""Checks for topics/chaos-composition: chaos as repeated self-composition.

Run: python3 checks.py   (sympy, numpy, mpmath). Prints PASS/FAIL per check and a total.
"""
import math, random
from fractions import Fraction
import numpy as np
import sympy as sp
import mpmath as mp

random.seed(1)
results = []


def check(name, ok, detail=""):
    results.append(bool(ok))
    print(("PASS " if ok else "FAIL ") + name + (("  [" + detail + "]") if detail else ""))


x, y, k, c, w = sp.symbols('x y k c w')
circ = lambda a, b: a + b + a * b
circK = lambda kk, a, b: a + b + kk * a * b

# ---------------------------------------------------------------- A. one arrow, iterated
print("A. Iteration is a functor; the dial changes units only")
f = lambda t: circ(t, t)
fk = lambda t: circK(k, t, t)
e1, e2 = x, k * x
for n in range(1, 4):
    e1, e2 = fk(e1), f(e2)
    check(f"A{n} k*(x o_k x)^(n={n}) = (kx) o (kx) iterated", sp.expand(k * e1 - e2) == 0)
lhs = x
for n in range(3):
    lhs = f(lhs)
check("A4 1 + (x o x)^(3) = (1+x)^8", sp.expand(1 + lhs - (1 + x) ** 8) == 0)
check("A5 k = -2 dial square is logistic r = 2", sp.expand(circK(-2, x, x) - 2 * x * (1 - x)) == 0)
check("A6 quadratic family in lift coords: 1 + (y o y + c) = (1+y)^2 + c",
      sp.expand(1 + f(y) + c - ((1 + y) ** 2 + c)) == 0)

# ---------------------------------------------------------------- B. folding
print("B. Logistic r = 4 is folded squaring")
inv = lambda u: -u / (1 + u)
fold = lambda u: u * inv(u) / 4
L4 = lambda t: 4 * t * (1 - t)
u, v = w, fold(w)
for n in range(1, 4):
    u, v = f(u), L4(v)
    check(f"B{n} fold(squaring^{n}) = logistic4^{n}(fold)", sp.simplify(fold(u) - v) == 0)
check("B4 fold(inv w) = fold(w): the fold forgets the circle inverse", sp.simplify(fold(inv(w)) - fold(w)) == 0)
# on |1+x| = 1, i.e. 1 + w = e^{i th}: fold(w) = sin^2(th/2)
th = sp.symbols('theta', real=True)
wth = sp.exp(sp.I * th) - 1
check("B5 on the circle |1+w|=1 the fold is sin^2(theta/2)",
      sp.simplify(sp.expand_complex(fold(wth)) - sp.sin(th / 2) ** 2) == 0)

# ---------------------------------------------------------------- C. Lyapunov
print("C. Lyapunov exponent = growth rate of the composed slope")
mp.mp.dps = 60
xx = mp.mpf('0.1234567')
acc = mp.mpf(0)
N = 4000
for i in range(N):
    acc += mp.log(abs(4 - 8 * xx))
    xx = 4 * xx * (1 - xx)
lam = float(acc / N)
check("C1 logistic r = 4: (1/N) sum log|f'| ~ log 2 (N = 4000, 60 digits)", abs(lam - math.log(2)) < 0.03,
      f"{lam:.4f} vs {math.log(2):.4f}")
# float run, longer
xf, accf, Nf = 0.3141592, 0.0, 2_000_000
arr = np.empty(Nf)
for i in range(Nf):
    arr[i] = xf
    xf = 4 * xf * (1 - xf)
    if xf <= 0 or xf >= 1:
        xf = random.random()
lamf = float(np.mean(np.log(np.abs(4 - 8 * arr))))
check("C2 same, float orbit of 2,000,000 steps", abs(lamf - math.log(2)) < 0.005, f"{lamf:.5f}")
# on the circle the slope of z->z^2 has size exactly 2
zz = complex(math.cos(1.0), math.sin(1.0))
slopes = []
for i in range(50):
    slopes.append(abs(2 * zz)); zz = zz * zz; zz /= abs(zz)
check("C3 squaring on |1+x| = 1: every step's slope has size 2, so Lyapunov = log 2 exactly",
      max(abs(s - 2) for s in slopes) < 1e-12)
# r = 2 (all circle-dial squares): superstable, Lyapunov -> -inf
x2 = 0.3
for i in range(60):
    x2 = 2 * x2 * (1 - x2)
check("C4 logistic r = 2: orbit lands on the superstable centre 1/2 (slope 0 there)", abs(x2 - 0.5) < 1e-12)

# ---------------------------------------------------------------- D. period doubling
print("D. Period doubling: superstable parameters and Feigenbaum's constants")


def superstable(n, guess):
    """r with f_r^(2^n)(1/2) = 1/2, Newton in mpmath."""
    mp.mp.dps = 40
    r = mp.mpf(guess)
    for _ in range(60):
        X, dX = mp.mpf('0.5'), mp.mpf(0)
        for _ in range(2 ** n):
            dX = X * (1 - X) + r * (1 - 2 * X) * dX
            X = r * X * (1 - X)
        step = (X - mp.mpf('0.5')) / dX
        r -= step
        if abs(step) < mp.mpf('1e-30'):
            break
    return r


R = [mp.mpf(2), superstable(1, 3.2)]
for n in range(2, 11):
    guess = R[-1] + (R[-1] - R[-2]) / 4.6692
    R.append(superstable(n, guess))
deltas = [float((R[i] - R[i - 1]) / (R[i + 1] - R[i])) for i in range(1, len(R) - 1)]
check("D1 r_0 = 2: the circle-dial square is the first superstable logistic map", R[0] == 2)
check("D2 r_1 = 1 + sqrt 5", abs(float(R[1]) - (1 + math.sqrt(5))) < 1e-12)
check("D3 superstable gaps shrink by Feigenbaum delta 4.6692", abs(deltas[-1] - 4.669201609) < 1e-4,
      f"{deltas[-1]:.6f}")
# alpha: distance from 1/2 to its nearest orbit neighbour f^(2^(n-1))(1/2)
dists = []
for n in range(1, 10):
    r = R[n]; X = mp.mpf('0.5')
    for _ in range(2 ** (n - 1)):
        X = r * X * (1 - X)
    dists.append(X - mp.mpf('0.5'))
alpha = float(dists[-2] / dists[-1])
check("D4 the rescaling ratio of f o f tends to Feigenbaum alpha -2.5029", abs(alpha + 2.502907875) < 2e-3,
      f"{alpha:.5f}")

# ---------------------------------------------------------------- E. order is invisible to iteration
print("E. f o g and g o f are the same system through g")
fq = lambda t: 1.3 * t * (1 - t) + 0.05
gq = lambda t: 0.9 * t + 0.07 * t * t
ok = True
for _ in range(20):
    x0 = random.random(); a = x0; b = gq(x0)
    for n in range(30):
        a = fq(gq(a)); b = gq(fq(b))
        ok &= abs(gq(a) - b) < 1e-9
check("E1 g((f o g)^n x) = (g o f)^n (g x), random quadratic f, g, 30 steps", ok)
ok = True
for _ in range(200):
    sa, ta, sb, tb = [random.uniform(-2, 2) for _ in range(4)]
    s = sa * sb
    if abs(1 - s) < 1e-3:
        continue
    cAB = (ta + sa * tb) / (1 - s)        # A after B
    cBA = (tb + sb * ta) / (1 - s)
    a_, b_ = sa - 1, sb - 1
    wedge = a_ * tb - b_ * ta
    ok &= abs((cAB - cBA) - wedge / (1 - s)) < 1e-9
    ok &= abs((sb * cAB + tb) - cBA) < 1e-9   # B carries centre(AB) to centre(BA)
check("E2 affine: same stretch both ways; centres differ by wedge/(1 - s_A s_B); B carries one to the other", ok)
# multipliers of periodic orbits agree: Lyapunov is order-free
xa = 0.2; acc1 = acc2 = 0.0
dfq = lambda t: 1.3 * (1 - 2 * t); dgq = lambda t: 0.9 + 0.14 * t
xb = gq(xa)
for n in range(20000):
    acc1 += math.log(abs(dfq(gq(xa)) * dgq(xa))); xa = fq(gq(xa))
    acc2 += math.log(abs(dgq(fq(xb)) * dfq(xb))); xb = gq(fq(xb))
check("E3 Lyapunov exponents of f o g and g o f agree", abs(acc1 - acc2) / 20000 < 1e-3,
      f"{acc1 / 20000:.5f} vs {acc2 / 20000:.5f}")

# ---------------------------------------------------------------- F. order makes chaos
print("F. Branch growth: list corner (log 2) against bag corner (0)")


def composites(maps, n):
    """distinct composites of length n of affine maps x -> s x + t, counted exactly (Fractions)."""
    maps = [(Fraction(s_), Fraction(t_)) for s_, t_ in maps]
    cur = {(Fraction(1), Fraction(0))}
    for _ in range(n):
        cur = {(s_ * S, s_ * T + t_) for (S, T) in cur for (s_, t_) in maps}
    return len(cur)


nF = 12
c_bin = composites([('1/2', 0), ('1/2', '1/2')], nF)
c_tr = composites([(1, 0), (1, 1)], nF)
c_st = composites([('1/2', '1/4'), ('3/10', '7/20')], nF)   # both stretch about 1/2
c_dc = composites([('1/2', 0), ('3/10', '3/5')], nF)        # centres 0 and 6/7
check("F1 halving branches x/2, (x+1)/2 (centres 0, 1): 2^n distinct composites", c_bin == 2 ** nF, str(c_bin))
check("F2 shifts x, x+1 (k = 0, shared centre at infinity): n + 1 composites", c_tr == nF + 1, str(c_tr))
check("F3 two stretches about one centre, slopes 1/2 and 0.3: n + 1 composites", c_st == nF + 1, str(c_st))
check("F4 same two slopes, different centres: 2^n composites", c_dc == 2 ** nF, str(c_dc))


def golden_count(n):
    """exact count of distinct sum_{i<n} e_i s^i, s = 1/phi, in Z[s] with s^2 = 1 - s."""
    A, B = [1], [0]
    for i in range(n):
        A.append(B[-1]); B.append(A[-2] - B[-1])
    sums = np.zeros((1, 2), dtype=np.int64)
    for i in range(n):
        sums = np.unique(np.vstack([sums, sums + np.array([A[i], B[i]])]), axis=0)
    return len(sums)


gc = [golden_count(n) for n in range(10, 23)]
ratios = [gc[i + 1] / gc[i] for i in range(len(gc) - 1)]
phi = (1 + 5 ** 0.5) / 2
check("F5 golden slope 1/phi (a {0,1} coincidence s + s^2 = 1): growth ratio -> phi, not 2",
      abs(ratios[-1] - phi) < 0.01, f"ratio {ratios[-1]:.4f}")
check("F6 golden slope: still exponential, so the semigroup is not a bag", ratios[-1] > 1.5)

# ---------------------------------------------------------------- F'. partial order: shared centres
print("F'. Some branches share a centre: entropy from the centre pattern (Cartier-Foata)")
tsym = sp.symbols('t')


def trace_counts(classes, n):
    """words of length n in the free product of free commutative monoids with these class sizes:
    coefficients of 1/mu(t), mu = 1 + sum((1 - t)^m - 1) (Cartier-Foata Mobius polynomial)."""
    mu = 1 + sum((1 - tsym) ** m - 1 for m in classes)
    ser = sp.series(1 / mu, tsym, 0, n + 1).removeO()
    return [int(ser.coeff(tsym, i)) for i in range(1, n + 1)]


st = lambda s_, m: (Fraction(s_), Fraction(m) * (1 - Fraction(s_)))   # stretch by s about centre m


def composites_list(maps, n):
    cur = {(Fraction(1), Fraction(0))}; out = []
    for _ in range(n):
        cur = {(s_ * S, s_ * T + t_) for (S, T) in cur for (s_, t_) in maps}; out.append(len(cur))
    return out


patterns = {
    (2, 1): [st('1/10', 0), st('1/11', 0), st('1/13', 1)],
    (1, 1, 1): [st('1/10', 0), st('1/11', 1), st('1/13', 3)],
    (2, 2): [st('1/10', 0), st('1/11', 0), st('1/13', 1), st('1/17', 1)],
    (2, 1, 1): [st('1/10', 0), st('1/11', 0), st('1/13', 1), st('1/17', 3)],
    (3,): [st('1/10', 0), st('1/11', 0), st('1/13', 0)],
}
for pat, maps in patterns.items():
    nn = 8
    got, want = composites_list(maps, nn), trace_counts(pat, nn)
    check(f"F'1 centre pattern {pat}: distinct composites = Cartier-Foata trace count", got == want, f"{got[:6]}")
check("F'2 pattern (2,1) counts are the even Fibonacci numbers F(2n+2): 3, 8, 21, 55, ...",
      trace_counts((2, 1), 8) == [int(sp.fibonacci(2 * n + 2)) for n in range(1, 9)])
rates = {}
for pat in [(1, 1), (2,), (2, 1), (1, 1, 1), (2, 2), (2, 1, 1), (3,)]:
    mu = 1 + sum((1 - tsym) ** m - 1 for m in pat)
    pos = [float(r_) for r_ in sp.real_roots(sp.Poly(sp.expand(mu), tsym)) if r_ > 0]
    rates[pat] = 1 / min(pos)
check("F'3 growth rates: (1,1) 2, (2) 1, (2,1) phi^2, (1,1,1) 3, (2,2) 2+sqrt2, (2,1,1) 2+sqrt3, (3) 1",
      all(abs(rates[pt] - v) < 1e-9 for pt, v in [((1, 1), 2), ((2,), 1), ((2, 1), phi ** 2), ((1, 1, 1), 3),
                                                   ((2, 2), 2 + 2 ** .5), ((2, 1, 1), 2 + 3 ** .5), ((3,), 1)]))
big = composites_list([st('1/2', 0), st('1/3', 0), st('1/5', 1)], 8)
check("F'4 the separation assumption matters: slopes 1/2, 1/3, 1/5 meet by accident (986 < 987 at n = 7)",
      big[6] == 986 and trace_counts((2, 1), 8)[6] == 987)
sa_, sb_, ta_, tb_, X_ = sp.symbols('s_a s_b t_a t_b X')
A_ = lambda z: sa_ * z + ta_; B_ = lambda z: sb_ * z + tb_
Ai_ = lambda z: (z - ta_) / sa_; Bi_ = lambda z: (z - tb_) / sb_
comm = sp.simplify(A_(B_(Ai_(Bi_(X_)))) - X_)
check("F'5 the commutator A B A^-1 B^-1 is the shift by the wedge (s_a - 1) t_b - (s_b - 1) t_a",
      sp.expand(comm - ((sa_ - 1) * tb_ - (sb_ - 1) * ta_)) == 0)

# ---------------------------------------------------------------- G. topological entropy from laps
print("G. Topological entropy from laps of f^n (Milnor-Thurston)")


def laps(r, n, M=2 ** 20):
    t = np.linspace(0, 1, M)
    X = t.copy(); D = np.ones(M)
    for _ in range(n):
        D = D * r * (1 - 2 * X); X = r * X * (1 - X)
    s_ = np.sign(D); s_ = s_[s_ != 0]
    return 1 + int(np.sum(s_[1:] != s_[:-1]))


check("G1 r = 4: f^n has 2^n laps (n = 10)", laps(4.0, 10) == 2 ** 10)
l35 = [laps(3.5, n) for n in (10, 11, 12)]
check("G2 r = 3.5 (period-4 part of the cascade): lap ratio falling toward 1, entropy 0",
      l35[2] / l35[1] < l35[1] / l35[0] < 1.35, f"laps {l35}")
l383 = [laps(3.83, n) for n in (11, 12)]
h = math.log(l383[1] / l383[0])
check("G3 r = 3.83 (period-3 window): entropy = log of the golden ratio 0.481", abs(h - math.log(phi)) < 0.03,
      f"{h:.3f}")

# ---------------------------------------------------------------- H. Pollard rho
print("H. Finite sets: every orbit is a rho (Pollard, Floyd)")


def floyd(fun, x0):
    t = h_ = x0
    while True:
        t = fun(t); h_ = fun(fun(h_))
        if t == h_:
            break
    mu, t = 0, x0
    while t != h_:
        t = fun(t); h_ = fun(h_); mu += 1
    lam, h_ = 1, fun(t)
    while t != h_:
        h_ = fun(h_); lam += 1
    return mu, lam


nP = 8051
fx = lambda t: (t * t + 1) % nP
fy = lambda t: (t * t + 2 * t + 1) % nP       # y o y + 1 with y = x - 1: (y o y) = 2y + y^2
t = h2 = 2
d = 1
while d == 1:
    t = fx(t); h2 = fx(fx(h2)); d = math.gcd(abs(t - h2), nP)
check("H1 Pollard rho with x^2 + 1 on 8051 finds the factor 97", d == 97, f"gcd {d}")
check("H2 the circle-op form y o y + 1 (y = x - 1) has the same rho: tail and loop lengths",
      floyd(fx, 2) == floyd(fy, 1), f"{floyd(fx, 2)}")
for p in (97, 83):
    mu, lam = floyd(lambda t: (t * t + 1) % p, 2)
    check(f"H3 mod {p}: rho with tail {mu} and loop {lam} (finite orbits are always rho-shaped)", mu >= 0 and lam >= 1)

# ---------------------------------------------------------------- I. streams
print("I. Doubling is tail on binary streams; folding makes Gray code")
bits = [random.randint(0, 1) for _ in range(64)]
X = sum(Fraction(b, 2 ** (i + 1)) for i, b in enumerate(bits))
ok = True
for n in range(40):
    Xn = (X * 2 ** n) % 1
    digits = [int((Xn * 2 ** (j + 1)) % 2 >= 1) for j in range(20)]
    ok &= digits == bits[n:n + 20]
check("I1 the n-th doubling of 0.b1b2b3... is 0.b(n+1)b(n+2)...: doubling = tail", ok)
Y = X + Fraction(1, 2 ** 30)
gaps = [float(abs(((X - Y) * 2 ** n + Fraction(1, 2)) % 1 - Fraction(1, 2))) for n in range(32)]
check("I2 points 2^-30 apart separate by a factor 2 per step, reaching size 1/2 at step 29",
      all(abs(gaps[n] - 2.0 ** (n - 30)) < 1e-15 for n in range(29)) and abs(gaps[29] - 0.5) < 1e-15)
mp.mp.dps = 80
thv = mp.mpf(int("".join(map(str, bits)), 2)) / mp.mpf(2) ** 64
ok = True
for n in range(40):
    xn = mp.sin(mp.pi * thv * 2 ** n) ** 2
    sym = 1 if xn > mp.mpf('0.5') else 0
    ok &= sym == (bits[n] ^ bits[n + 1])
check("I3 logistic r = 4 itinerary (x > 1/2 or not) is the Gray code b(n) xor b(n+1) of the angle's digits", ok)

# ---------------------------------------------------------------- J. the order dial q
print("J. Turning the order dial q partway: entropy is a switch at q = 1")
import itertools
qs_ = sp.symbols('q')


def W_rs(n, qq):
    """sum over 0/1 words of length n of q^(out-of-order pairs) = sum_k [n,k]_q (Rogers-Szego at 1):
    W_{m+1} = 2 W_m - (1 - q^m) W_{m-1}."""
    a, b = mp.mpf(1), mp.mpf(2)
    if n == 0:
        return a
    for m in range(1, n):
        a, b = b, 2 * b - (1 - mp.mpf(qq) ** m) * a
    return b


def W_brute(n):
    return sp.expand(sum(qs_ ** sum(1 for i in range(n) for j in range(i + 1, n) if wd[i] > wd[j])
                         for wd in itertools.product([0, 1], repeat=n)))


def W_gauss(n):
    tot = 0
    for kk in range(n + 1):
        num = sp.Mul(*[(1 - qs_ ** (n - i)) for i in range(kk)]); den = sp.Mul(*[(1 - qs_ ** (i + 1)) for i in range(kk)])
        tot += sp.cancel(num / den)
    return sp.expand(tot)


check("J1 weighted word count = sum of Gaussian binomials (n = 7, symbolic in q)", sp.expand(W_brute(7) - W_gauss(7)) == 0)
check("J2 the dial's ends: q = 0 gives n + 1 (bag), q = 1 gives 2^n (list), n = 1..10",
      all(W_gauss(n).subs(qs_, 0) == n + 1 and W_gauss(n).subs(qs_, 1) == 2 ** n for n in range(1, 11)))
mp.mp.dps = 30
okb = True
for qq in (0.5, 0.9, 0.99):
    for n in (50, 400, 4000):
        wv = W_rs(n, qq); okb &= (n + 1 <= wv <= (n + 1) / mp.qp(qq))
check("J3 for q < 1: n + 1 <= W_n(q) <= (n + 1)/(q;q)_inf, so growth is linear and entropy is 0 (q = .5, .9, .99)", okb)
r99 = [float(mp.log(W_rs(n, 0.99)) / n) for n in (20, 4000)]
check("J4 q = 0.99 tracks the list early (rate > 0.6 at n = 20) and falls away later (< 0.05 at n = 4000)",
      r99[0] > 0.6 and r99[1] < 0.05, f"{r99[0]:.3f}, {r99[1]:.3f}")
r11 = [float(mp.log(W_rs(n, 1.1)) / n ** 2) for n in (100, 400)]
check("J5 for q > 1 the count grows like q^(n^2/4): log W / n^2 -> log(q)/4 (q = 1.1), faster than any exponential",
      abs(r11[1] - math.log(1.1) / 4) < abs(r11[0] - math.log(1.1) / 4) < 0.002, f"{r11[1]:.5f} vs {math.log(1.1)/4:.5f}")

print()
print(f"{sum(results)} of {len(results)} pass")
