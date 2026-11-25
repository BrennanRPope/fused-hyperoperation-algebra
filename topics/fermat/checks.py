"""Checks for 'Fermat's Last Theorem through the circle operation'.
The circle operation: a (+) b = a + b + ab = (1+a)(1+b) - 1.
"""
import cmath, math
from fractions import Fraction
from sympy import symbols, expand, bernoulli, primerange, series, exp, log, factorint, Poly

t, T, x = symbols('t T x')

def madd(a, b):
    return a + b + a*b

# 1. Rotated circle-operation cube: (+)_{k} zeta^k t = t^p for odd p
print("== 1. The circle operation of p rotated copies ==")
for p in [3, 5, 7, 11]:
    z = cmath.exp(2j*math.pi/p)
    for tv in [0.37, 1.9+0.4j]:
        acc = 0
        for k in range(p):
            acc = madd(acc, z**k*tv)
        print(p, tv, abs(acc - tv**p) < 1e-9)
    # layer sums of the cube: elementary symmetric polys of zeta^k
    import itertools
    layers = []
    for r in range(p+1):
        s = sum(math.prod([z**k for k in c]) if c else 1 for c in itertools.combinations(range(p), r))
        layers.append(round(abs(s), 9))
    print("  layer sums |e_r(zeta^k)|:", layers)

# Fermat in circle-operation form: with t=y/x, s=z/x-1: (+)zeta^k t == [p](s) == (1+s)^p - 1
print("Check 3^p..: x^p+y^p=z^p <=> (+)_k zeta^k (y/x) = (1+s)^p - 1, s=z/x-1 (identity check on 3,4,5 squares p=2 analog skipped)")

# 2. Torsion of the circle operation: [p](a) = (1+a)^p - 1; roots a = zeta^k - 1; pi=zeta-1, |prod over k!=0 of (zeta^k-1)| = p
print("\n== 2. p-torsion of the circle operation ==")
for p in [5, 7, 13]:
    z = cmath.exp(2j*math.pi/p)
    prod = 1
    for k in range(1, p):
        prod *= (z**k - 1)
    print(p, "prod_{k=1}^{p-1}(zeta^k-1) =", round(prod.real, 9), round(prod.imag, 9))
    # [p](T) mod p == T^p  (height 1)
    P = Poly(expand((1+T)**p - 1), T)
    print("  [p](T) coeffs mod p:", [c % p for c in P.all_coeffs()])

# 3. Bernoulli numbers from the circle-operation exponential exp_G(x)=e^x-1
print("\n== 3. Bernoulli from circle-operation exponential ==")
ser = series(x/(exp(x)-1), x, 0, 12).removeO()
print([ser.coeff(x, k)*math.factorial(k) for k in range(12)])
# Coleman power series of cyclotomic unit: [a](T)/T. log of it in log-coordinates
a = 2
g = log((exp(a*x)-1)/(a*(exp(x)-1)))
d = series(g.diff(x), x, 0, 9).removeO()
for k in range(1, 9):
    c = d.coeff(x, k-1)*math.factorial(k-1)  # k-th derivative of g at 0
    pred = (a**k - 1)*bernoulli(k)/k if k > 1 else None
    print(k, c, pred)

# irregular primes
print("\n== 4. irregular primes < 200 ==")
irr = []
for p in primerange(5, 200):
    ks = [k for k in range(2, p-2, 2) if Fraction(int(bernoulli(k).p), int(bernoulli(k).q)).numerator % p == 0]
    if ks:
        irr.append((p, ks))
print(irr)

# relative class numbers h^- via h^- = 2p * prod_{chi odd} (-B_{1,chi}/2)
print("\n== 5. relative class numbers h^-(Q(zeta_p)) ==")
def prim_root(p):
    for g in range(2, p):
        if all(pow(g, (p-1)//q, p) != 1 for q in factorint(p-1)):
            return g
def hminus(p):
    g = prim_root(p)
    ind = {}
    v = 1
    for i in range(p-1):
        ind[v] = i; v = v*g % p
    prod = 1
    for j in range(1, p-1, 2):  # odd characters: chi(g)=w^j, j odd
        w = cmath.exp(2j*math.pi*j/(p-1))
        B1 = sum(aa*w**ind[aa] for aa in range(1, p))/p
        prod *= -B1/2
    return 2*p*prod
hm = {}
for p in primerange(3, 72):
    h = hminus(p)
    hm[p] = round(h.real)
    print(p, round(h.real), abs(h.imag) < 1e-3*max(1, abs(h)), "p|h-" if round(h.real) % p == 0 else "")

# 6. Fermat curve point counts over F_q vs Hasse-Weil bound
print("\n== 6. Fermat curve x^n+y^n=z^n points over F_q (projective) ==")
def count(n, q):
    pw = {}
    for v in range(q):
        pw[v] = pow(v, n, q)
    N = 0
    # projective points: z=1 affine + z=0 points
    for xx in range(q):
        for yy in range(q):
            if (pw[xx] + pw[yy] - 1) % q == 0:
                N += 1
    # at infinity: x^n + y^n = 0, (x:y:0) with y=1
    for xx in range(q):
        if (pw[xx] + 1) % q == 0:
            N += 1
    return N
for n in [3, 4, 5]:
    g = (n-1)*(n-2)//2
    for q in [7, 13, 31, 61]:
        N = count(n, q)
        print(n, g, q, N, abs(N-(q+1)), "<=", round(2*g*math.sqrt(q), 2))

# 7. Tate curve torsion structure check (complex model): E = C*/q^Z, p-torsion = zeta^i q^(j/p)
print("\n== 7. Tate curve: p-torsion count ==")
q = 0.3*cmath.exp(0.7j); p = 5
pts = set()
for i in range(p):
    for j in range(p):
        w = cmath.exp(2j*math.pi*i/p)*q**(j/p)
        # check w^p in q^Z
        r = w**p
        k = round(math.log(abs(r))/math.log(abs(q)))
        assert abs(r - q**k) < 1e-9
        pts.add((round(w.real, 8), round(w.imag, 8)))
print("p-torsion points found:", len(pts), "= p^2 =", p*p)

# 8. Genus via Riemann-Hurwitz for Fermat cover of P^1 branched at 0,1,inf
print("\n== 8. Genus of Fermat curves ==")
for n in range(2, 9):
    chi = n*n*(2-3) + 3*n
    print(n, 1 - chi//2, (n-1)*(n-2)//2)
