# Group theory and Lie groups through the circle operation: every claim on the page, checked.
import sympy as sp, math, random, itertools
import numpy as np
from mpmath import mp, quad, mpf, inf, log as mlog, atan as matan, sqrt as msqrt
mp.dps = 30
random.seed(7)
ok = 0; tot = 0
def check(name, cond):
    global ok, tot
    tot += 1; ok += bool(cond); print(('PASS ' if cond else 'FAIL ') + name)
z = lambda e: sp.simplify(sp.expand(e)) == 0
a, b, c, x, y, p, q, k, lam, t = sp.symbols('a b c x y p q k lam t')
N = lambda p, q, a, b: a + b + p*a*b
Dn = lambda p, q, a, b: 1 - q*a*b
F = lambda p, q, a, b: N(p, q, a, b)/Dn(p, q, a, b)
V = lambda p, q, x: 1 + p*x + q*x**2

# --- the landscape laws are groups (recheck of the landscape page) ---
check("F_{p,q} associative, identity 0, inverse -a/(1+pa)",
      z(F(p,q,F(p,q,a,b),c)-F(p,q,a,F(p,q,b,c))) and z(F(p,q,a,0)-a) and z(F(p,q,a,-a/(1+p*a))))

# --- Möbius matrices ---
X = sp.Matrix([[p, 1], [-q, 0]])
I2 = sp.eye(2)
mob = lambda s: I2 + s*X
check("mob(a) = [[1+pa, a],[-qa, 1]] = I + a·X", mob(a) == sp.Matrix([[1+p*a, a], [-q*a, 1]]))
check("mob(a)·mob(b) = (1-qab)·mob(F(a,b))  (the law is matrix multiplication up to scale)",
      sp.simplify(mob(a)*mob(b) - Dn(p,q,a,b)*mob(F(p,q,a,b))) == sp.zeros(2))
check("Cayley-Hamilton: X² - pX + qI = 0", sp.simplify(X*X - p*X + q*I2) == sp.zeros(2))
check("char poly of X is t² - p t + q, discriminant p² - 4q",
      sp.expand(X.charpoly(t).as_expr() - (t**2 - p*t + q)) == 0)
mobfun = lambda M, s: (M[0,0]*s + M[0,1])/(M[1,0]*s + M[1,1])
check("mob(a) acting on b by x ↦ (αx+β)/(γx+δ) gives F(a,b)", z(mobfun(mob(a), b) - F(p,q,a,b)))
# vector field of a Möbius family at the identity
al, be, ga, de = sp.symbols('al be ga de')
Mt = I2 + t*sp.Matrix([[al, be], [ga, de]])
fld = sp.diff(mobfun(Mt, x), t).subs(t, 0)
check("Möbius generator [[α,β],[γ,δ]] moves x at speed β + (α-δ)x - γx²", z(fld - (be + (al-de)*x - ga*x**2)))
check("the landscape generator X moves x at speed 1 + px + qx² (the Riccati rate)",
      z(fld.subs({al: p, be: 1, ga: -q, de: 0}) - V(p,q,x)))

# --- theta algebra: 1 + x·θ with θ² = pθ - q ---
th = sp.symbols('th')
prod = sp.rem(sp.expand((1+a*th)*(1+b*th)), th**2 - p*th + q, th)
check("(1+aθ)(1+bθ) = (1-qab) + (a+b+pab)θ  when θ² = pθ - q", z(prod - (Dn(p,q,a,b) + N(p,q,a,b)*th)))
check("p=1,q=0: θ²=θ is an idempotent; 1+xθ ↦ (1+x, 1) recovers the lift 1+x",
      z(sp.rem(sp.expand((1+a*th)*(1+b*th)), th**2-th, th) - (1 + (a+b+a*b)*th)))
check("p=0,q=-1: θ²=1 (split-complex j), law is (a+b)/(1+ab)", z(F(0,-1,a,b)-(a+b)/(1+a*b)))
check("p=0,q=1: θ²=-1 (i), law is (a+b)/(1-ab), tan addition", z(F(0,1,a,b)-(a+b)/(1-a*b)))
check("p=q=0: θ²=0 (dual numbers), law is a+b", z(F(0,0,a,b)-(a+b)))
nrm = lambda u0, u1: u0**2 + p*u0*u1 + q*u1**2
check("norm x²+pxy+qy² is multiplicative; norm(1+aθ) = 1+pa+qa² (the rate)",
      z(nrm(Dn(p,q,a,b), N(p,q,a,b)) - nrm(1,a)*nrm(1,b)) and z(nrm(1,a)-V(p,q,a)))
check("rate law: V(a)·V(b) = (1-qab)²·V(F(a,b))", z(V(p,q,a)*V(p,q,b) - Dn(p,q,a,b)**2*V(p,q,F(p,q,a,b))))
check("∂F/∂b = V(a)/(1-qab)²", z(sp.diff(F(p,q,a,b), b) - V(p,q,a)/Dn(p,q,a,b)**2))
check("invariant ruler: ∂F/∂b · V(b) = V(F(a,b)), so dx/V(x) is carried by every translation",
      z(sp.diff(F(p,q,a,b), b)*V(p,q,b) - V(p,q,F(p,q,a,b))))

# --- logarithms in closed form ---
def Lclosed(pv, qv, xv):
    pv, qv, xv = mpf(pv), mpf(qv), mpf(xv); D = pv*pv - 4*qv
    if abs(qv) < mpf(10)**-25 and abs(pv) < mpf(10)**-25: return xv
    if D > 0:
        r = msqrt(D); u = (-pv + r)/2; v = (-pv - r)/2   # 1+px+qx² = (1-ux)(1-vx)
        return mlog((1 - v*xv)/(1 - u*xv))/(u - v)
    if D < 0:
        r = msqrt(-D)
        return 2/r*(matan((2*qv*xv + pv)/r) - matan(pv/r))
    return xv/(1 + pv*xv/2)
def Lnum(pv, qv, xv): return quad(lambda s: 1/(1 + pv*s + qv*s*s), [0, xv])
good = True
for (pv, qv) in [(1,0), (0,-1), (0,1), (2,1), (1,-2), (-0.7,0.3), (0.5,0.0625), (1.3,2)]:
    for xv in [0.1, 0.25, -0.2]:
        if abs(Lclosed(pv,qv,xv) - Lnum(pv,qv,xv)) > 1e-20: good = False
check("closed-form logarithm L(x) = ∫₀ˣ ds/(1+ps+qs²) in all three cases (8 points × 3 x)", good)
good = True
for (pv, qv) in [(1,0), (0,-1), (0,1), (2,1), (1,-2), (-0.7,0.3), (0.5,0.0625), (1.3,2)]:
    for _ in range(5):
        av, bv = mpf(random.uniform(-.2,.2)), mpf(random.uniform(-.2,.2))
        Fv = (av+bv+pv*av*bv)/(1-qv*av*bv)
        if abs(Lclosed(pv,qv,Fv) - Lclosed(pv,qv,av) - Lclosed(pv,qv,bv)) > 1e-20: good = False
check("L(F(a,b)) = L(a) + L(b): the logarithm turns every landscape law into addition", good)
check("named rulers: circle-operation log(1+x), tanh law artanh x, tan law arctan x",
      abs(Lclosed(1,0,mpf('0.3'))-mlog(mpf('1.3')))<1e-25 and abs(Lclosed(0,-1,mpf('0.3'))-mp.atanh(mpf('0.3')))<1e-25 and abs(Lclosed(0,1,mpf('0.3'))-matan(mpf('0.3')))<1e-25)
check("border (p²=4q): L = x/(1+px/2), and F is plain addition in L",
      z(sp.diff(x/(1+p*x/2), x) - 1/V(p, p**2/4, x)) and
      z((lambda L: L(F(p,p**2/4,a,b)) - L(a) - L(b))(lambda s: s/(1+p*s/2))))
check("dial: L = log(1+kx)/k → x as k → 0", sp.limit(sp.log(1+k*x)/k, k, 0) == x)
check("elliptic period: ∫ over all reals of dx/(1+px+qx²) = 2π/√(4q-p²)  (3 points)",
      all(abs(quad(lambda s: 1/(1+pv*s+qv*s*s), [-inf, inf]) - 2*mp.pi/msqrt(4*qv-pv*pv)) < 1e-20
          for pv, qv in [(0,1), (1,1), (-2,3)]))
check("real walls exactly when p² - 4q > 0 (roots of 1+px+qx²)",
      all(((pv*pv-4*qv) > 0) == all(abs(sp.im(r)) < 1e-12 for r in sp.Poly(qv*x**2+pv*x+1, x).nroots()) for pv, qv in [(1,0.1),(0,-1),(0,1),(2,1.5),(3,2),(-1,2)]))
check("line case: L runs to ±∞ at the two walls (the circle operation: wall -1)", sp.limit(sp.log(1+x), x, -1, '+') == -sp.oo)
check("circle case: the point ∞ belongs to the group; tan law F(a, ∞) = -1/a",
      sp.limit(F(0,1,a,x), x, sp.oo) == -1/a)

# --- Lie algebra of the line: Lie's sl2 = {1, x, x²}·∂ ---
def br(f, g): return sp.expand(f*sp.diff(g, x) - g*sp.diff(f, x))
A1, B1, C1, A2, B2, C2, A3, B3, C3 = sp.symbols('A1 B1 C1 A2 B2 C2 A3 B3 C3')
P1, P2, P3 = A1+B1*x+C1*x**2, A2+B2*x+C2*x**2, A3+B3*x+C3*x**2
check("bracket of fields α+βx+γx² stays in the span: (αβ'-βα') + 2(αγ'-γα')x + (βγ'-γβ')x²",
      z(br(P1,P2) - ((A1*B2-B1*A2) + 2*(A1*C2-C1*A2)*x + (B1*C2-C1*B2)*x**2)))
check("Jacobi identity for the field bracket", z(br(P1,br(P2,P3)) + br(P2,br(P3,P1)) + br(P3,br(P1,P2))))
# Killing form via ad matrices on basis (1, x, x²)
basis = [sp.Integer(1), x, x**2]
def coords(f): f = sp.Poly(sp.expand(f), x); return [f.coeff_monomial(1), f.coeff_monomial(x), f.coeff_monomial(x**2)]
def ad(f): return sp.Matrix([coords(br(f, e)) for e in basis]).T
Pv = 1 + p*x + q*x**2
check("Killing form of the landscape generator: tr(ad²) = 2(p² - 4q)", z((ad(Pv)*ad(Pv)).trace() - 2*(p**2-4*q)))
Bf = lambda u, v: sp.expand((ad(u)*ad(v)).trace())
check("Killing form is ad-invariant: B([g,f],h) + B(f,[g,h]) = 0", z(Bf(br(P2,P1),P3) + Bf(P1,br(P2,P3))))
X0 = X - (p/2)*I2
check("matrix form: trace-free part X₀ has tr(X₀²) = (p² - 4q)/2", z((X0*X0).trace() - (p**2-4*q)/2))

# --- re-coordinatizing the line: x ↦ x/(1-cx) slides along parabolas of fixed p²-4q ---
pp, qq = p + 2*c, q + p*c + c**2
check("rate in the new chart: (1+cy)²·V(y/(1+cy)) = 1 + (p+2c)y + (q+pc+c²)y²",
      z(sp.simplify((1+c*y)**2*V(p,q,y/(1+c*y))) - V(pp,qq,y)))
check("discriminant unchanged: (p+2c)² - 4(q+pc+c²) = p² - 4q", z(pp**2-4*qq-(p**2-4*q)))
phi = lambda s: s/(1 - c*s)
check("conjugation: F_{p+2c, q+pc+c²}(φa, φb) = φ(F_{p,q}(a,b)), φ(x) = x/(1-cx)",
      z(sp.simplify(F(pp,qq,phi(a),phi(b)) - phi(F(p,q,a,b)))))
check("The circle operation (1,0) with c = -1/2 lands on (0,-1/4); scaled by 2 that is velocity addition (0,-1)",
      (pp.subs({p:1,q:0,c:sp.Rational(-1,2)}), qq.subs({p:1,q:0,c:sp.Rational(-1,2)})) == (0, sp.Rational(-1,4))
      and z(F(0,sp.Rational(-1,4),2*a,2*b) - 2*F(0,-1,a,b)))
check("parity (-2,0) with c = 1 lands exactly on velocity addition (0,-1); same parabola p²-4q = 4",
      (pp.subs({p:-2,q:0,c:1}), qq.subs({p:-2,q:0,c:1})) == (0, -1) and z(F(-2,0,a,b)-(a+b-2*a*b)))
check("units: scaling x ↦ λx sends (p,q) to (λp, λ²q) and p²-4q to λ²(p²-4q)",
      z(lam*F(lam*p,lam**2*q,a,b) - F(p,q,lam*a,lam*b)) and z((lam*p)**2-4*lam**2*q - lam**2*(p**2-4*q)))
pts = [(1,0),(0,-1),(0,1),(2,1),(-3,1),(1,5)]
cls = lambda pv,qv: (pv*pv-4*qv > 0) - (pv*pv-4*qv < 0)
check("conjugation and scaling never change the sign of p²-4q (random c, λ on 6 points)",
      all(cls(pv+2*cv, qv+pv*cv+cv*cv) == cls(pv,qv) == cls(l*pv, l*l*qv) for pv,qv in pts for cv in (-1.5,.3,2) for l in (.5,3)))
check("rung s at (s, s(s-1)/4) has p² - 4q = s (landscape placement)", z((p**2 - 4*(p*(p-1)/4)) - p))

# --- the dial inside the ax+b group ---
T = lambda kk, aa, xx: (1 + kk*aa)*xx + aa
check("T_a(x) = (1+ka)x + a equals a ✳ₖ x", z(T(k,a,x) - (a + x + k*a*x)))
check("T_a ∘ T_b = T_{a ✳ₖ b}", z(T(k,a,T(k,b,x)) - T(k, a+b+k*a*b, x)))
check("T_a fixes -1/k (the wall), for every a", z(T(k,a,-1/k) + 1/k))
check("about the wall it is a pure dilation: T_a(y - 1/k) + 1/k = (1+ka)·y", z(T(k,a,y-1/k) + 1/k - (1+k*a)*y))
check("scaling conjugation: T^k_{λa}(λx) = λ·T^{kλ}_a(x)", z(T(k,lam*a,lam*x) - lam*T(k*lam,a,x)))
check("k = 0 is a translation with no fixed point (a ≠ 0)", sp.solve(sp.Eq(T(0,a,x), x), x) == [] )
check("generators 1 + kx: [x∂, ∂] = -∂, and [1+kx, 1+k'x] = (k'-k)·1, a translation",
      z(br(x, sp.Integer(1)) + 1) and z(br(1+k*x, 1+c*x) - (c-k)))
check("Inönü–Wigner shadow: (a+b)/(1-qab) → a+b as q → 0 (Lorentz → Galilei)", sp.limit(F(0,q,a,b), q, 0) == a+b)

# --- finite fields ---
def primes(n): return [m for m in range(2, n) if all(m % d for d in range(2, int(m**.5)+1))]
def mp_pow(g, n, P):
    r = 0
    for _ in range(n): r = (r + g + r*g) % P
    return r
cyc = True
for P in primes(60):
    if P == 2: continue
    found = any(len({mp_pow(g, j, P) for j in range(P-1)}) == P-1 for g in range(P-1))
    cyc = cyc and found
check("(F_p minus the wall, ✳) is cyclic of order p-1 for every odd prime p < 60", cyc)
def is_primroot(g, P): return len({pow(g, j, P) for j in range(1, P)}) == P-1
check("1 ✳ 1 ✳ ... (n copies) = 2ⁿ - 1 (Mersenne numbers)", all(mp_pow(1, n, 10**9+7) == (2**n-1) % (10**9+7) for n in range(30)))
check("1 generates (F_p minus wall, ✳) exactly when 2 is a primitive root mod p (odd p < 500)",
      all((len({mp_pow(1, j, P) for j in range(P-1)}) == P-1) == is_primroot(2, P) for P in primes(500) if P > 2))
check("small cases: generator 1 for p = 5, 11, 13; for p = 7 use 2 (1 has order 3)",
      all(len({mp_pow(1,j,P) for j in range(P-1)})==P-1 for P in (5,11,13)) and len({mp_pow(2,j,7) for j in range(6)})==6 and mp_pow(1,3,7)==0)
cnt = lambda P: sum(1 for u in range(P) for v in range(P) if (u*u+v*v) % P == 1)
check("circle x²+y²=1 over F_p has p+1 points when p ≡ 3 mod 4, p-1 when p ≡ 1 mod 4 (p < 60)",
      all(cnt(P) == (P+1 if P % 4 == 3 else P-1) for P in primes(60) if P > 2))
def tanlaw_group(P):
    # projective line F_p ∪ {∞}; ∞ coded as P
    INF = P
    def op(u, v):
        if u == INF and v == INF: return 0
        if u == INF: return INF if v == 0 else (-pow(v, P-2, P)) % P
        if v == INF: return INF if u == 0 else (-pow(u, P-2, P)) % P
        num, den = (u+v) % P, (1 - u*v) % P
        return INF if den == 0 else num*pow(den, P-2, P) % P
    els = list(range(P)) + [INF]
    assoc = all(op(op(u,v),w) == op(u,op(v,w)) for u in els for v in els for w in els)
    def order(g):
        r, n = g, 1
        while r != 0: r = op(r, g); n += 1
        return n
    return assoc, max(order(g) for g in els)
check("tan law on F_p ∪ {∞} is a cyclic group of order p+1 for every prime p ≡ 3 mod 4 below 40 (the finite circle)",
      all(tanlaw_group(P) == (True, P+1) for P in primes(40) if P % 4 == 3))
def proj_count(pv, qv, P):
    # points of the projective line minus the walls (roots of 1 + p x + q x^2, with ∞ a wall when q = 0)
    walls = sum(1 for xx in range(P) if (1 + pv*xx + qv*xx*xx) % P == 0) + (1 if qv % P == 0 else 0)
    return P + 1 - walls
def leg(d, P):
    d %= P
    return 0 if d == 0 else (1 if pow(d, (P-1)//2, P) == 1 else -1)
check("finite shadows: the landscape group over F_p has p - ((p²-4q)/p) elements: p-1 line, p border, p+1 circle (all p,q mod 5,7,11,13)",
      all(proj_count(pv,qv,P) == P - leg(pv*pv-4*qv, P) for P in (5,7,11,13) for pv in range(P) for qv in range(P)))

# --- the twisting corner: matrices ---
def ma(A, B): return A + B + A*B
good = True; good2 = True
for _ in range(5):
    A = sp.Matrix(2, 2, [random.randint(-5,5) for _ in range(4)])
    B = sp.Matrix(2, 2, [random.randint(-5,5) for _ in range(4)])
    good = good and (ma(A,B) - ma(B,A) == A*B - B*A) and (I2 + ma(A,B) == (I2+A)*(I2+B))
check("the matrix circle operation: lift (I+A)(I+B) = I + A✳B, and A✳B - B✳A = AB - BA exactly", good)
A = np.array([[0.3,-1.2],[0.7,0.4]]); B = np.array([[-0.5,0.9],[1.1,0.2]]); E = np.eye(2)
errs = []
for tt in [1e-2, 5e-3, 2.5e-3]:
    G = (E+tt*A)@(E+tt*B)@np.linalg.inv(E+tt*A)@np.linalg.inv(E+tt*B) - E
    errs.append(np.abs(G - tt*tt*(A@B-B@A)).max()/tt**3)
check("group commutator (1+tA)(1+tB)(1+tA)⁻¹(1+tB)⁻¹ - 1 = t²[A,B] + O(t³)", max(errs) < 10 and errs[0]/errs[-1] < 1.5)
check("landscape matrices commute (both are I + sX), so the commutator vanishes on the landscape",
      sp.simplify(mob(a)*mob(b) - mob(b)*mob(a)) == sp.zeros(2))
# group algebra: g ↦ g - 1 turns products into the circle operation (permutation matrices of S3)
perm = lambda s: sp.Matrix(3, 3, lambda i, j: 1 if s[j] == i else 0)
gs = [perm(s) for s in itertools.permutations(range(3))]
I3 = sp.eye(3)
check("group algebra: (g-1) ✳ (h-1) = gh - 1 for all g, h in S3 (as permutation matrices)",
      all(ma(g-I3, h-I3) == g*h - I3 for g in gs for h in gs))
check("trace of AB - BA is 0: brackets land in sl2 (5 random pairs)",
      all((lambda A_, B_: (A_*B_-B_*A_).trace() == 0)(sp.Matrix(2,2,[random.randint(-4,4) for _ in range(4)]), sp.Matrix(2,2,[random.randint(-4,4) for _ in range(4)])) for _ in range(5)))

print(f"\n{ok} of {tot} checks pass")
