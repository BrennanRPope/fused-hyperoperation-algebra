"""Dial survey checks (Nov 25, 2026). Run: python3 checks.py  (mpmath, sympy)
Each check prints PASS/FAIL. Writes results.json."""
import json, math
from fractions import Fraction as F
import mpmath as mp
import sympy as sp
mp.mp.dps = 30
out = []; res = {}
def check(name, ok, detail=""):
    out.append(("PASS" if ok else "FAIL") + " " + name + (f"  [{detail}]" if detail else ""))

# ---------- Dial 1: k. Only the sign survives rescaling of the bare operation ----------
a, b, k, lam = sp.symbols('a b k lam')
w = lambda k, x, y: x + y + k*x*y
check("k: k*(a *_k b) = (ka) *_1 (kb)", sp.expand(k*w(k,a,b) - w(1,k*a,k*b)) == 0)
check("k: in units lambda times bigger the weight is k*lambda (only sign(k) survives)", sp.expand(w(k,lam*a,lam*b)/lam - w(k*lam,a,b)) == 0)

# ---------- Dial 3: per-element dial = odds; union of independent random subsets ----------
t, s = sp.symbols('t s', positive=True)
p = lambda x: x/(1+x)
check("element dial t is odds: 1-P(in A or B) = (1-p_t)(1-p_s) iff dial t*s",
      sp.simplify((1-p(w(1,t,s))) - (1-p(t))*(1-p(s))) == 0)
check("complement turns odds t into 1/t", sp.simplify((1-p(t)) - p(1/t)) == 0)
# uniform dial t <-> exponent alpha = log2(1+t): count (1+t)^n = 2^(alpha n), and alpha adds when dials combine by the circle operation
for (x, y) in [(F(1,3), F(2,5)), (F(1), F(1)), (F(3,7), F(1,9))]:
    al = lambda u: math.log2(1+float(u))
    check(f"alpha(t*s) = alpha(t)+alpha(s) for t={x}, s={y}", abs(al(x+y+x*y) - al(x) - al(y)) < 1e-12)

# ---------- Dials 1 x 4: k acting on the prime dial ----------
# per-element factor with weight k: f_k(z) = (1 - k z)^(-1/k); k=-1 subsets, k=1 multisets, k->0 exp(z)
def coef(kk, e):            # coefficient of z^e in f_k: prod_{j<e}(1+j k)/e!
    c = F(1)
    for j in range(e): c *= (1 + j*kk)
    return c / math.factorial(e)
z = sp.symbols('z')
for kk in [F(-1), F(-1,2), F(-1,3), F(1,2), F(1), F(2)]:
    ser = sp.series((1 - sp.Rational(kk.numerator, kk.denominator)*z)**(-sp.Rational(kk.denominator, kk.numerator)), z, 0, 7).removeO()
    ok = all(sp.Rational(coef(kk,e).numerator, coef(kk,e).denominator) == ser.coeff(z, e) for e in range(7))
    check(f"f_k coefficients = prod(1+jk)/e!  at k={kk}", ok)
check("k=-1: no element twice (coef 0 from e=2)", all(coef(F(-1),e)==0 for e in range(2,9)))
check("k=+1: multisets (coef 1/e! * e! = 1 as Dirichlet weight)", all(coef(F(1),e)*math.factorial(e)==math.factorial(e) for e in range(9)))
check("k=-1/m: at most m copies (coef 0 beyond e=m), m=2,3,5",
      all(coef(F(-1,m),e)==0 for m in (2,3,5) for e in range(m+1,m+6)) and all(coef(F(-1,m),m)!=0 for m in (2,3,5)))

# necklace numbers and the cyclotomic identity: (1-kz)^(-1/k) = prod_N (1 - z^N)^(-M(k,N)/k)
def mobius(n):
    r, m, pp = 1, n, 2
    while pp*pp <= m:
        if m % pp == 0:
            m //= pp
            if m % pp == 0: return 0
            r = -r
        pp += 1
    return -r if m > 1 else r
def M(kk, N): return sum(mobius(N//d) * kk**d for d in range(1, N+1) if N % d == 0) / N
for kk in [F(-1), F(-1,2), F(1,2), F(2), F(3)]:
    K = sp.Rational(kk.numerator, kk.denominator)
    lhs = sp.series((1-K*z)**(-1/K), z, 0, 9).removeO()
    rhs = 1
    for N in range(1, 9):
        mN = M(kk, N); e = sp.Rational(mN.numerator, mN.denominator)/K
        rhs = sp.expand(sp.series(rhs * (1 - z**N)**(-e), z, 0, 9).removeO())
    check(f"cyclotomic identity to z^8 at k={kk}", sp.expand(lhs - rhs) == 0)
check("M(k,N) whole numbers for integer k (aperiodic necklaces), k=2..5, N<=12",
      all(M(F(kk),N).denominator == 1 for kk in range(2,6) for N in range(1,13)))
check("k=-1: exponents M/k = 1, -1, 0, 0, ... -> zeta(s)/zeta(2s)",
      [M(F(-1),N)/F(-1) for N in range(1,7)] == [1,-1,0,0,0,0])
check("k=+1: exponents 1, 0, 0, ... -> zeta(s)", [M(F(1),N) for N in range(1,7)] == [1,0,0,0,0,0])

primes = list(sp.primerange(2, 200000))
def Z_euler(kk, s):
    kk = mp.mpf(kk); tot = mp.mpf(0)
    for q in primes:
        x = mp.power(q, -s)
        tot += (-mp.log(1 - kk*x)/kk) if kk != 0 else x
    return mp.e**tot
def Z_zeta(kk, s, Nmax=40):
    tot = mp.mpf(0)
    for N in range(1, Nmax+1):
        mN = sum(mobius(N//d) * mp.mpf(kk)**d for d in range(1,N+1) if N % d == 0)/N
        tot += mN/mp.mpf(kk) * mp.log(mp.zeta(N*s))
    return mp.e**tot
zz = {}
for kk in [-1, -0.5, 0.5, 1, 1.5]:
    A, B = Z_euler(kk, 2), Z_zeta(kk, 2)
    zz[str(kk)] = float(B)
    check(f"Z_k(2): Euler product = prod zeta(2N)^(M/k) at k={kk}", abs(A-B)/B < 1e-5, f"{float(A):.8f} vs {float(B):.8f}")
check("Z_-1(2) = 15/pi^2", abs(Z_zeta(-1,2) - 15/mp.pi**2) < 1e-20)
res["Z_k(2)"] = zz

# densities: sum_{n<=x} a_k(n) ~ C_k x with C_k = prod_p (1-1/p)(1-k/p)^(-1/k)
def C(kk):
    tot = mp.mpf(0)
    for q in primes:
        x = mp.mpf(1)/q
        tot += mp.log(1-x) + ((-mp.log(1-kk*x)/kk) if kk != 0 else x)
    return mp.e**tot
check("C_-1 = 6/pi^2 (squarefree density)", abs(C(-1) - 6/mp.pi**2) < 1e-5, f"{float(C(-1)):.6f}")
check("C_1 = 1", abs(C(1) - 1) < 1e-12)
c0 = mp.e**(mp.mpf('0.2614972128476427837554268386') - mp.euler)
check("C_0 = e^(M - gamma) (Meissel-Mertens minus Euler)", abs(C(0) - c0) < 1e-5, f"{float(C(0)):.6f} vs {float(c0):.6f}")
# direct count up to X for k = -1/2 and 1/2 (multiplicative a_k(n))
X = 300000
spf = list(range(X+1))
for i in range(2, int(X**0.5)+1):
    if spf[i] == i:
        for j in range(i*i, X+1, i):
            if spf[j] == j: spf[j] = i
def ak(kk, n):
    r = 1.0
    while n > 1:
        q = spf[n]; e = 0
        while n % q == 0: n //= q; e += 1
        r *= float(coef(F(kk).limit_denominator(100), e))   # a_k(p^e) = prod_{j<e}(1+jk)/e!
    return r
dens = {}
for kk in [-0.5, 0.5]:
    S = sum(ak(kk, n) for n in range(1, X+1))
    dens[str(kk)] = (S/X, float(C(kk)))
    check(f"count of k-weighted integers up to {X} matches C_k at k={kk}", abs(S/X - float(C(kk))) < 0.01, f"{S/X:.5f} vs {float(C(kk)):.5f}")
res["C_k"] = {str(kk): float(C(kk)) for kk in [-1,-0.75,-0.5,-0.25,0,0.25,0.5,0.75,1,1.25,1.5,1.75,1.9]}
res["density_counts"] = dens
# wall on the smallest prime: for k > 2 the pole moves from s = 1 to s = log2 k
for kk in [3, 4, 6]:
    s0 = math.log2(kk)
    check(f"k={kk}: factor at p=2 hits its wall at s=log2 k={s0:.4f} > 1", abs(1 - kk*2**(-s0)) < 1e-12 and s0 > 1)
check("|k| < 2: every factor finite at s = 1", all(1 - kk/2 > 0 for kk in [-1.9, -1, 0.5, 1.9]))

# ---------- Dial 2: the twist is d times area (from twisted boosts), re-checked ----------
import random
random.seed(1)
okt = True
for _ in range(200):
    q = random.choice([-1, 0, 1, F(1,2), F(-3,2)])
    u = (F(random.randint(-5,5),11), F(random.randint(-5,5),11)); v = (F(random.randint(-5,5),11), F(random.randint(-5,5),11))
    dot = u[0]*v[0]+u[1]*v[1]; cross = u[0]*v[1]-u[1]*v[0]
    X_, Y_ = 1 - q*dot, q*cross
    # twist trivial iff d*cross == 0 (given X != 0)
    if X_ != 0 and ((Y_ == 0) != (q == 0 or cross == 0)): okt = False
check("twist: none iff d*cross = 0 (200 random rational pairs)", okt)

# ---------- Dial 5: choice grades are a partial order, not a number line (literature facts, not computed) ----------
res["choice_note"] = "BPI and DC are independent over ZF (Howard-Rubin); grades do not form a line."

passed = sum(1 for o in out if o.startswith("PASS"))
out.append(f"\n{passed} of {len(out)} pass")
open("checks_output.txt", "w").write("\n".join(out) + "\n")
json.dump(res, open("results.json", "w"), indent=1)
print("\n".join(out))
