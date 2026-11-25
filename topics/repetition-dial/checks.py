"""Checks for the repetition dial (FHA 1.0 checklist item 1).

Question: does "how much repetition" (how many times one element may be used)
fall out of the framework, and if so is it k, and does k sweep it smoothly?

Atom with twist k: Z_k(x) = (1 + k x)^(1/k)  (k=1 subsets, k=0 e^x, k=-1 multisets).
"""
import sympy as sp, numpy as np, mpmath as mp
from math import comb, factorial
from fractions import Fraction as Fr

ok = 0; total = 0
def check(name, cond):
    global ok, total
    total += 1; ok += bool(cond)
    print(("PASS " if cond else "FAIL ") + name)

k, x, y, a, b = sp.symbols('k x y a b')

# 1. k-log composition: phi(Z) = (Z^k - 1)/k turns products into a + b + k a b
A = 1 + k*a; B = 1 + k*b
check("1 (1+ka)(1+kb) = 1 + k(a+b+kab)", sp.expand(A*B - (1 + k*(a + b + k*a*b))) == 0)
Zk = (1 + k*x)**(1/k)
check("1b k-log of one atom is its fugacity x", sp.simplify((sp.powdenest(Zk**k, force=True) - 1)/k - x) == 0)

# 2. coefficients of Z_k are polynomials in k: prod_{j<n}(1 - j k)/n!
def coeff_formula(kk, n):
    p = sp.Integer(1)
    for j in range(n): p *= (1 - j*kk)
    return p / sp.factorial(n)
ser = sp.series(Zk, x, 0, 7).removeO()
good = all(sp.simplify(ser.coeff(x, n) - coeff_formula(k, n)) == 0 for n in range(7))
check("2 coefficient of x^n in (1+kx)^(1/k) = prod(1-jk)/n!", good)
check("2b k=0 gives e^x (1/n!)", all(coeff_formula(0, n) == sp.Rational(1, factorial(n)) for n in range(10)))
check("2c k=1 gives 1+x", [coeff_formula(1, n) for n in range(5)] == [1, 1, 0, 0, 0])
check("2d k=-1 gives 1/(1-x) (all 1)", all(coeff_formula(-1, n) == 1 for n in range(10)))

# 3. cells: k = 1/m -> m exclusive cells; k = -1/m -> m bosonic cells
c3 = True
for m in range(1, 7):
    for n in range(0, 12):
        c3 &= coeff_formula(sp.Rational(1, m), n) == sp.Rational(comb(m, n), m**n)
        c3 &= coeff_formula(sp.Rational(-1, m), n) == sp.Rational(comb(m + n - 1, n), m**n)
check("3 k=1/m: C(m,n)/m^n (cap at m copies); k=-1/m: C(m+n-1,n)/m^n (unlimited)", c3)

# 4. sign map: k <= 0 all positive; k > 0 all nonnegative only when 1/k is an integer
def first_negative(kk, nmax=60):
    p = Fr(1)
    for n in range(1, nmax):
        p *= (1 - (n-1)*kk)
        if p < 0: return n
    return None
grid_neg = [Fr(-p, q) for q in range(1, 13) for p in range(0, 2*q+1)]
check("4a every k in [-2,0] (rational grid) has all coefficients > 0", all(first_negative(kk) is None for kk in grid_neg))
grid_pos = sorted({Fr(p, q) for q in range(1, 13) for p in range(1, 2*q+1)})
recip = [kk for kk in grid_pos if kk.numerator == 1]
nonrec = [kk for kk in grid_pos if kk.numerator != 1]
check("4b k = 1/m: no negative coefficient", all(first_negative(kk) is None for kk in recip))
check("4c k > 0 not of form 1/m: some coefficient negative", all(first_negative(kk) is not None for kk in nonrec))
check("4d k=2/3: coefficient of x^3 = -1/54", Fr(1)*(1-Fr(2,3))*(1-2*Fr(2,3))/6 == Fr(-1, 54))

# 5. mean occupation and fluctuation of a k-atom
n_k = sp.simplify(x*sp.diff(sp.log(Zk), x))
check("5a mean copies n = x/(1+kx), i.e. 1/(e^(beta eps) + k)", sp.simplify(n_k - x/(1 + k*x)) == 0)
var_k = sp.simplify(x*sp.diff(n_k, x))
check("5b variance = n(1 - k n): k = (n - var)/n^2 is constant", sp.simplify(var_k - n_k*(1 - k*n_k)) == 0)
check("5c 1/n - 1/x = k exactly", sp.simplify(1/n_k - 1/x - k) == 0)

# 6. Haldane-Wu exclusion statistics g: w^g (1+w)^(1-g) = 1/x, n = 1/(w+g)
mp.mp.dps = 40
def wu_w(g, X):
    # monotone in u = ln w; bisection on a wide bracket
    f = lambda u: g*u + (1-g)*mp.log(1+mp.e**u) + mp.log(X)
    lo, hi = mp.mpf(-400), mp.mpf(400)
    for _ in range(mp.mp.prec + 20):
        mid = (lo + hi)/2
        if f(mid) > 0: hi = mid
        else: lo = mid
    return mp.e**((lo + hi)/2)
def wu_n(g, X):
    w = wu_w(g, X); return 1/(w + g)
def keff(g, X): return 1/wu_n(g, X) - 1/mp.mpf(X)
# Wu's n agrees with x d/dx ln Z for Z = (1+w)/w (check numerically)
def lnZ(g, X): w = wu_w(g, X); return mp.log((1 + w)/w)
for g in [0.25, 0.5, 0.75]:
    X = mp.mpf('0.7')
    d = X*mp.diff(lambda t: lnZ(g, t), X)
    check(f"6a g={g}: n = x dlnZ/dx with Z=(1+w)/w", abs(d - wu_n(g, X)) < 1e-20)
c6 = True
for g in [0.25, 0.5, 0.75]:
    lo = keff(g, mp.mpf('1e-8')); hi = keff(g, mp.mpf('1e8'))
    c6 &= abs(lo - (2*g - 1)) < 1e-6 and abs(hi - g) < 1e-6
check("6b Haldane effective k slides from 2g-1 (dilute) to g (crowded)", c6)
check("6c g->0 gives constant k = -1 (bosons)", all(abs(keff(mp.mpf('1e-30'), X) + 1) < 1e-12 for X in [0.01, 0.3, 0.9]))
check("6d g=1 is exactly the k=1 atom", all(abs(keff(1, X) - 1) < 1e-20 for X in [0.01, 1, 100]))
# algebraic proof for semions (g=1/2): constant k would force w to one value
w_, kk = sp.symbols('w kk')
expr = sp.expand((w_ + sp.Rational(1, 2) - kk)**2 - w_*(1 + w_))
check("6e g=1/2: (w+1/2-k)^2 = w(1+w) is -2k w + (1/2-k)^2 = 0 (one w only)", sp.simplify(expr - (-2*kk*w_ + (sp.Rational(1,2)-kk)**2)) == 0)

# 7. hard cap at m copies (Gentile): Z = 1 + x + ... + x^m
for m in [2, 3]:
    Zg = sum(x**j for j in range(m+1))
    # is it (1+kx)^(1/k) for any k?  c2 would force k=1-2*1 = -1 (m>=2), then c3 must be 1 if m>=3 else 0
    kfit = 1 - 2*1   # from c2 = (1-k)/2 = 1
    c3_pred = coeff_formula(kfit, 3); c3_true = 1 if m >= 3 else 0
    check(f"7a cap m={m} is not a single k (predicted c3={c3_pred}, true {c3_true}; c4 check)",
          (c3_pred != c3_true) or (coeff_formula(kfit, m+1) != 0))
    check(f"7b cap m={m} = boson(x)/boson(x^{m+1}): two k=-1 atoms at sizes 1 and {m+1}",
          sp.simplify(Zg - (1/(1-x))/(1/(1-x**(m+1)))) == 0)
# effective k for the cap drifts too
def keff_cap(m, X):
    X = mp.mpf(X); Z = sum(X**j for j in range(m+1)); n = sum(j*X**j for j in range(m+1))/Z
    return 1/n - 1/X
check("7c cap m=2: effective k drifts (dilute -1 -> crowded 1/2)",
      abs(keff_cap(2, '1e-8') + 1) < 1e-6 and abs(keff_cap(2, '1e8') - 0.5) < 1e-6)

# 8. smooth in k across 0: Z_k(x) and dZ/dk continuous at k=0
Zf = lambda kk, X: mp.e**X if kk == 0 else (1 + kk*X)**(1/kk)
X = mp.mpf('0.6')
left = Zf(-1e-12, X); right = Zf(1e-12, X)
check("8a Z_k continuous at k=0 (-> e^x)", abs(left - mp.e**X) < 1e-9 and abs(right - mp.e**X) < 1e-9)
dk = lambda kk: (Zf(kk + 1e-6, X) - Zf(kk - 1e-6, X))/2e-6
check("8b dZ/dk continuous at k=0 (= -x^2/2 e^x)", abs(dk(1e-3) - dk(-1e-3)) < 1e-2 and abs(dk(2e-4) + X**2/2*mp.e**X) < 1e-2)

# 9. k=-1 is the k=+1 ladder at doubled sizes (size-dial residual link)
prod = 1
for j in range(6): prod *= (1 + x**(2**j))
check("9 1/(1-x) = prod (1+x^(2^j)) through x^63", sp.expand(sp.series(1/(1-x), x, 0, 64).removeO() - sp.expand(prod)) == 0)

print(f"\n{ok}/{total} checks passed")
