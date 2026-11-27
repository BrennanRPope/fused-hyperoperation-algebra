"""Is the plane dial's d Hoffman-Ihara's q-deformation?  (DRAFT checks)
Run: python3 checks.py   (sympy, exact)"""
import json, random, sympy as sp
random.seed(5)
res = []
def ok(name, c, note=""):
    res.append({"id": name, "pass": bool(c), "note": note}); print(("PASS " if c else "FAIL ") + name)

q, k, d, x, y, x0, al, be = sp.symbols("q k d x y x0 alpha beta")
m, n = sp.symbols("m n", integer=True, positive=True)
def qint(N, qq): return sum(qq**j for j in range(N))
def line(a, b, kk): return a + b + kk*a*b
def plane(a, b, kk, dd): return (a + b + kk*a*b) / (1 - dd*a*b)

# Q1  q-integers add by our line law with k = q - 1
Q1 = all(sp.expand(qint(M+N, q) - line(qint(M, q), qint(N, q), q - 1)) == 0
         for M in range(0, 7) for N in range(0, 7))
ok("Q1 [m+n]_q = [m]_q o_k [n]_q with k = q - 1  (m, n = 0..6)", Q1)

# Q2  Hoffman-Ihara's q letter product z_i <> z_j = z_{i+j} + (1-q) z_{i+j-1} on the q-MZV letter values
#     f_i(n) = q^{(i-1)n} / [n]_q^i  (Bradley's q-MZV letters); 1 - q = -k
Q2 = True
for N in range(1, 6):
    for i in range(1, 4):
        for j in range(1, 4):
            f = lambda a: q**((a-1)*N) / qint(N, q)**a
            Q2 &= sp.simplify(f(i)*f(j) - f(i+j) - (1-q)*f(i+j-1)) == 0
ok("Q2 the q-MZV letter values obey z_i<>z_j = z_(i+j) + (1-q) z_(i+j-1); the correction is the identity q^n = 1 + k[n]_q, k = q - 1", Q2)
Q2b = all(sp.expand(q**N - (1 + (q-1)*qint(N, q))) == 0 for N in range(0, 8))
ok("Q2b q^n = 1 + (q - 1)[n]_q: the 'cost' term in the letter product is the line dial's k", Q2b)

# Q3  no-go: a quasi-shuffle with scalar letters gives the shadow x + y + k(x<>y) = x + y + c x y (a line law);
#     the plane law with d != 0 is not a polynomial: its x^2 y^2 coefficient is d*(1 + ...) - check the x^2 y coefficient
ser = sp.series(sp.series(plane(x, y, k, d), y, 0, 3).removeO(), x, 0, 3).removeO()
c21 = sp.expand(ser).coeff(x, 2).coeff(y, 1)
ok("Q3 the plane law's x^2 y coefficient is d, so for d != 0 it is not x + y + c xy: no scalar-letter quasi-shuffle produces it", sp.simplify(c21 - d) == 0)

# Q4  the plane law is the addition law of (P,Q)-numbers: N with eigenvalues alpha, beta (k = a+b, d = ab),
#     g = 1 + x0 N, g^n = A_n + B_n N, x_n = B_n / A_n.  Then x_(m+n) = F(x_m, x_n) and
#     B_n = x0 [n]_{P,Q}, P = 1 + x0 alpha, Q = 1 + x0 beta.
def AB(N, a, b, s):
    P, Qb = 1 + s*a, 1 + s*b
    B = (P**N - Qb**N) / (a - b); A = (a*Qb**N - b*P**N) / (a - b)
    return sp.simplify(A), sp.simplify(B)
Q4 = True; Q4b = True
for _ in range(6):
    a, b = sp.Rational(random.randint(-9, 9), random.randint(1, 5)), sp.Rational(random.randint(-9, 9), random.randint(1, 5))
    if a == b: continue
    s = sp.Rational(random.randint(1, 9), 10)
    kk, dd = a + b, a*b
    for M in range(1, 5):
        for N in range(1, 5):
            Am, Bm = AB(M, a, b, s); An, Bn = AB(N, a, b, s); Amn, Bmn = AB(M+N, a, b, s)
            Q4 &= sp.simplify(Bmn/Amn - plane(Bm/Am, Bn/An, kk, dd)) == 0
    # matrix check: (I + s N)^n with N companion of t^2 - k t + d equals A_n I + B_n N
    Nm = sp.Matrix([[0, -dd], [1, kk]])
    for N in range(1, 6):
        A, B = AB(N, a, b, s)
        Q4b &= sp.simplify((sp.eye(2) + s*Nm)**N - (A*sp.eye(2) + B*Nm)) == sp.zeros(2)
        P, Qb = 1 + s*a, 1 + s*b
        Q4b &= sp.simplify(B - s*sp.simplify((P**N - Qb**N)/(P - Qb))) == 0
ok("Q4 the plane law adds the orbit x_n = B_n/A_n of 1 + x0 N exactly: x_(m+n) = F(x_m, x_n)", Q4)
ok("Q4b (1 + x0 N)^n = A_n + B_n N with B_n = x0 [n]_{P,Q}, P, Q = 1 + x0 alpha, 1 + x0 beta (two bases)", Q4b)

# Q5  d = 0 iff one base is 1: then x_n = x0 [n]_q with q = 1 + k x0, the one-base q-deformation
A5 = True
for N in range(1, 7):
    P = 1 + x0*k
    Bn = sp.simplify(((1 + x0*k)**N - 1) / k)
    A5 &= sp.simplify(Bn - x0*qint(N, P)) == 0
ok("Q5 d = 0 (beta = 0, Q = 1): B_n = x0 [n]_q, q = 1 + k x0, A_n = 1: one base, the q-deformation", A5)
# d != 0 means both bases move: d = alpha beta = (P-1)(Q-1)/x0^2
A6 = sp.simplify(al*be - ((1 + x0*al) - 1)*((1 + x0*be) - 1)/x0**2) == 0
ok("Q6 d = (P - 1)(Q - 1)/x0^2: d measures how far the second base Q sits from 1", A6)

# Q7  Mobius picture: x_n = x0 [n]_r / (Q - beta x0 [n]_r) with r = P/Q, so a plane orbit is a q-orbit
#     (q = r) seen through the Mobius map u -> u / (Q - beta u)
Q7 = True
for _ in range(5):
    a, b = sp.Rational(random.randint(1, 9), 3), sp.Rational(-random.randint(1, 9), 4)
    s = sp.Rational(random.randint(1, 5), 10)
    P, Qb = 1 + s*a, 1 + s*b; r = P/Qb
    for N in range(1, 6):
        A, B = AB(N, a, b, s)
        Q7 &= sp.simplify(B/A - s*qint(N, r)/(Qb - b*s*qint(N, r))) == 0
ok("Q7 each plane orbit is a q-orbit with q = P/Q seen through a Mobius map; d sets that map, not q", Q7)

npass = sum(r["pass"] for r in res)
print(f"\n{npass}/{len(res)} pass")
json.dump(res, open("checks.json", "w"), indent=1)
