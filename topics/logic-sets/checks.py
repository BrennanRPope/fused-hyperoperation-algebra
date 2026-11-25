"""Checks for the logic and set theory page. Every claim on index.html is tested here."""
import itertools, random
import sympy as sp
import mpmath as mp

mp.mp.dps = 40
results = []
def check(name, ok):
    results.append((name, bool(ok)))
    print(("PASS" if ok else "FAIL"), "-", name)

a, b, c, k, g, x, y = sp.symbols('a b c k gamma x y')
w = lambda k_, u, v: u + v + k_*u*v          # the dial
ma = lambda u, v: u + v + u*v                # mult-add

# ---------- 1. Boolean logic on the dial ----------
check("dial table on 0,1: 0,1,1 and 1 (+)k 1 = 2 + k",
      [sp.expand(w(k, p, q)) for p, q in [(0,0),(0,1),(1,0),(1,1)]] == [0, 1, 1, 2 + k])
red = lambda e: sp.expand(sp.Poly(sp.expand(e), a, b).as_expr()).subs({a**2: a, b**2: b})
def reduce_idem(e):
    e = sp.expand(e)
    for _ in range(4):
        e = sp.expand(e.subs({a**4: a, b**4: b}).subs({a**3: a, b**3: b}).subs({a**2: a, b**2: b}))
    return e
check("idempotent defect: (a(+)k b)^2 - (a(+)k b) = (k+1)(k+2)ab when a^2=a, b^2=b",
      sp.expand(reduce_idem(w(k,a,b)**2 - w(k,a,b)) - (k+1)*(k+2)*a*b) == 0)
check("roots of the defect are k = -1 (OR) and k = -2 (XOR)", set(sp.solve((k+1)*(k+2), k)) == {-1, -2})
check("on the integers, 1 (+)k 1 in {0,1} only for k in {-2,-1} (k from -50 to 50)",
      [kk for kk in range(-50, 51) if w(kk, 1, 1) in (0, 1)] == [-2, -1])
bits = [0, 1]
check("k=-1 on bits is OR", all(w(-1,p,q) == (p or q) for p in bits for q in bits))
check("k=-2 on bits is XOR", all(w(-2,p,q) == (p ^ q) for p in bits for q in bits))
check("product on bits is AND", all(p*q == (p and q) for p in bits for q in bits))
def fold(kk, n):
    s = 0
    for _ in range(n): s = w(kk, s, 1)
    return s
qnum = lambda q, n: sum(q**i for i in range(n))
check("n trues under (+)k give the q-number [n]_q with q = k+1 (k in -3..3, n<12)",
      all(fold(kk, n) == qnum(kk+1, n) for kk in range(-3, 4) for n in range(12)))
check("OR folds n>=1 trues to 1 ([n]_0 = 1)", all(fold(-1, n) == 1 for n in range(1, 30)))
check("XOR folds n trues to n mod 2 ([n]_-1)", all(fold(-2, n) == n % 2 for n in range(30)))
check("distributivity defect a(b(+)k c) - (ab (+)k ac) = k abc (1-a)",
      sp.expand(a*w(k,b,c) - w(k,a*b,a*c) - k*a*b*c*(1-a)) == 0)
check("absorption defect a (+)-1 ab = a + ab(1-a)", sp.expand(w(-1,a,a*b) - (a + a*b*(1-a))) == 0)
check("OR self: a(+)-1 a - a = a - a^2;  XOR self: a(+)-2 a = 2(a - a^2)",
      sp.expand(w(-1,a,a) - a - (a - a**2)) == 0 and sp.expand(w(-2,a,a) - 2*(a - a**2)) == 0)
check("Boolean ring axioms for (XOR, AND) on bits",
      all(w(-2,p,w(-2,q,r)) == w(-2,w(-2,p,q),r) and p*w(-2,q,r) == w(-2,p*q,p*r) and w(-2,p,p) == 0
          for p in bits for q in bits for r in bits))
check("Boolean lattice laws for (OR, AND) on bits (absorption, distributivity both ways)",
      all(w(-1,p,p*q) == p and p*w(-1,p,q) == p and w(-1,p,q*r) == w(-1,p,q)*w(-1,p,r)
          for p in bits for q in bits for r in bits))
check("OR = (a XOR b) XOR (a AND b) on bits", all(w(-1,p,q) == w(-2,w(-2,p,q),p*q) for p in bits for q in bits))
check("XOR = (a OR b)(1 - ab) on bits", all(w(-2,p,q) == w(-1,p,q)*(1-p*q) for p in bits for q in bits))
check("lattice order on bits: a (+)-1 b = b iff ab = a", all((w(-1,p,q) == q) == (p*q == p) for p in bits for q in bits))
check("XOR lift: 1 - 2(a(+)-2 b) = (1-2a)(1-2b) (bits become signs, XOR becomes x)",
      sp.expand(1 - 2*w(-2,a,b) - (1-2*a)*(1-2*b)) == 0)

# ---------- 2. walls as truth values ----------
check("OR's wall is 1 (true absorbs): 1 (+)-1 b = 1", sp.expand(w(-1,1,b)) == 1)
check("XOR's wall is 1/2: m (+)-2 1/2 = 1/2 for every m", sp.expand(w(-2,a,sp.Rational(1,2))) == sp.Rational(1,2))
# one-time pad: message bit with P(1)=m, independent fair key, cipher = m XOR key
m = sp.Rational(3, 10)
pc = sum(pm*pk for (mb,pm) in [(0,1-m),(1,m)] for (kb,pk) in [(0,sp.S.Half),(1,sp.S.Half)] if mb ^ kb == 1)
check("one-time pad: P(cipher = 1) = 1/2 for message chance 3/10, and it equals m (+)-2 1/2", pc == sp.S.Half == w(-2,m,sp.S.Half))
joint = {(mb, mb ^ kb): 0 for mb in bits for kb in bits}
for mb, pm in [(0,1-m),(1,m)]:
    for kb in bits: joint[(mb, mb ^ kb)] += pm*sp.S.Half
check("one-time pad: cipher independent of message (joint = product of marginals)",
      all(joint[(mb,cb)] == (m if mb else 1-m)*sp.S.Half for mb in bits for cb in bits))

# ---------- 3. indicator functions and sets ----------
random.seed(1)
U = range(40)
ok1 = ok2 = ok3 = True
for _ in range(50):
    A = {u for u in U if random.random() < .5}; B = {u for u in U if random.random() < .4}
    for u in U:
        ia, ib = int(u in A), int(u in B)
        ok1 &= int(u in A | B) == w(-1, ia, ib)
        ok2 &= int(u in A ^ B) == w(-2, ia, ib)
        ok3 &= int(u in A & B) == ia*ib
check("indicators: 1_(A union B) = a (+)-1 b", ok1)
check("indicators: 1_(A sym-diff B) = a (+)-2 b", ok2)
check("indicators: 1_(A intersect B) = ab", ok3)
ok = True
for _ in range(30):
    S = [{u for u in U if random.random() < .5} for _ in range(4)]
    lhs = len(set().union(*S))
    rhs = sum((-1)**(len(T)+1) * len(set.intersection(*[S[i] for i in T]))
              for r in range(1,5) for T in itertools.combinations(range(4), r))
    ok &= lhs == rhs
check("inclusion-exclusion for 4 sets = signed corners of the 4-cube (1,4,6,4 layers)", ok)
ok = True
for na in range(6):
    for nb in range(6):
        A = [('a',i) for i in range(na)]; B = [('b',j) for j in range(nb)]
        D = [('L',p) for p in A] + [('R',q) for q in B] + [('P',(p,q)) for p in A for q in B]
        ok &= len(D) == ma(na, nb)
        prod = [(p,q) for p in A+['*'] for q in B+['*'] if (p,q) != ('*','*')]
        def to_d(pq):
            p,q = pq
            return ('R',q) if p=='*' else ('L',p) if q=='*' else ('P',(p,q))
        ok &= sorted(map(to_d, prod), key=str) == sorted(D, key=str)
check("|A + B + AxB| = a (+) b and (A+*)x(B+*) minus the base point matches it one-to-one (a,b<6)", ok)
check("finite sizes are >= 0: no set has the wall's size -1 or an inverse size -u/(1+u) < 0 (u>0)",
      all(-sp.Rational(u,1+u) < 0 for u in range(1,20)))
# Euler characteristic of polyhedral sets: a set is a list of open cells by dimension
chi = lambda cells: sum((-1)**d * n for d, n in enumerate(cells))
def cprod(X, Y):
    out = [0]*(len(X)+len(Y)-1)
    for i, p in enumerate(X):
        for j, q in enumerate(Y): out[i+j] += p*q
    return out
def cadd(*Zs):
    L = max(map(len, Zs)); return [sum(Z[i] if i < len(Z) else 0 for Z in Zs) for i in range(L)]
I = [0, 1]              # open interval: one open 1-cell
pt = [1]
check("open interval has Euler characteristic -1 (Schanuel's negative set)", chi(I) == -1)
check("half-open interval [0,1) = point + open interval has Euler characteristic 0 (the lift of the wall)", chi(cadd(pt, I)) == 0)
check("wall: chi(I + B + IxB) = -1 for finite B of size 0..9 and for B = I, B = square",
      all(chi(cadd(I, [nb], cprod(I,[nb]))) == -1 for nb in range(10)) and
      chi(cadd(I, I, cprod(I, I))) == -1 and chi(cadd(I, cprod(I,I), cprod(I,cprod(I,I)))) == -1)
X2 = [0, 2]             # two open intervals, chi = -2
check("two open intervals X: chi(X + X + XxX) = 0, so X is its own mult-add inverse at the level of chi",
      chi(X2) == -2 and chi(cadd(X2, X2, cprod(X2, X2))) == 0 and ma(-2, -2) == 0)
check("chi of products multiplies and so chi(A + B + AxB) = chi A (+) chi B (random cell lists)",
      all(chi(cadd(P, Q, cprod(P, Q))) == ma(chi(P), chi(Q))
          for P in ([random.randint(0,3) for _ in range(3)] for _ in range(20))
          for Q in ([random.randint(0,3) for _ in range(3)] for _ in range(5))))

# ---------- 4. fuzzy logic ----------
S = (a + b + (g-2)*a*b)/(1 + (g-1)*a*b)
T = a*b/(g + (1-g)*(a + b - a*b))
lens = lambda p, q, u, v: (u + v + p*u*v)/(1 - q*u*v)
odds = lambda t: t/(1-t)
check("Hamacher t-conorm is the landscape law at p = gamma - 2, q = 1 - gamma", sp.simplify(S - lens(g-2, 1-g, a, b)) == 0)
check("those points lie on the line p + q = -1", sp.expand((g-2) + (1-g)) == -1)
check("gamma = 1 is the probabilistic sum (dial k = -1)", sp.simplify(S.subs(g,1) - w(-1,a,b)) == 0)
check("gamma = 2 is the Einstein sum (a+b)/(1+ab), relativity's point (0,-1)", sp.simplify(S.subs(g,2) - lens(0,-1,a,b)) == 0)
check("gamma = 0 is (a+b-2ab)/(1-ab) at (-2,1), on the border parabola p^2 = 4q",
      sp.simplify(S.subs(g,0) - lens(-2,1,a,b)) == 0 and (-2)**2 == 4*1)
check("miss: 1 - S = (1-a)(1-b)/(1+(gamma-1)ab), independent of gamma in the numerator",
      sp.simplify(1 - S - (1-a)*(1-b)/(1 + (g-1)*a*b)) == 0)
check("Hamacher is the dial on odds: odds(S) = odds a (+)gamma odds b", sp.simplify(odds(S) - w(g, odds(a), odds(b))) == 0)
check("Hamacher t-norm is the De Morgan mirror: T(a,b) = 1 - S(1-a,1-b)",
      sp.simplify(T - (1 - S.subs({a:1-a, b:1-b}, simultaneous=True))) == 0)
check("the t-norm is the dial on odds-against (1-x)/x",
      sp.simplify((1-T)/T - w(g, (1-a)/a, (1-b)/b)) == 0)
check("gamma = 0 adds odds (dial k = 0), gamma = 1 mult-adds odds (k = 1)",
      sp.simplify(odds(S.subs(g,0)) - odds(a) - odds(b)) == 0 and sp.simplify(odds(S.subs(g,1)) - ma(odds(a), odds(b))) == 0)
check("Hamacher's additive generator log(1 + gamma*odds(x)) is the log of the dial's lift on odds",
      all(abs(mp.log(1+gg*odds(mp.mpf(u)+0)) + mp.log(1+gg*odds(mp.mpf(v))) -
              mp.log(1+gg*odds((u+v+(gg-2)*u*v)/(1+(gg-1)*u*v)))) < mp.mpf(10)**-30
          for gg in [mp.mpf('0.5'), 1, 2, 7] for u in [mp.mpf('0.1'), mp.mpf('0.6')] for v in [mp.mpf('0.3'), mp.mpf('0.85')]))
check("odds wall of dial gamma is o = -1/gamma, i.e. x = 1/(1-gamma); gamma = 2 gives x = -1 (light)",
      sp.simplify((-1/g)/(1 + (-1/g)) - 1/(1-g)) == 0 and (1/(1-sp.Integer(2))) == -1)
check("on truth values every Hamacher t-conorm with gamma != 0 is OR",
      all(sp.simplify(S.subs({a:p, b:q}) - (p or q)) == 0 for p in bits for q in bits if (p,q) != (1,1))
      and sp.simplify(S.subs({a:1, b:1})) == 1)
def St(gg, u, v):
    d = 1 + (gg-1)*u*v
    if d == 0 and u == 1 and v == 1 and gg == 0: return Fr(1)   # gamma = 0 corner, value by continuity
    return (u + v + (gg-2)*u*v)/d
from fractions import Fraction as Fr
grid = [Fr(i, 20) for i in range(21)]
def is_tconorm(gg):
    try:
        vals = {(u,v): St(gg,u,v) for u in grid for v in grid}
    except ZeroDivisionError:
        return False
    inrange = all(-1e-12 <= z <= 1+1e-12 for z in vals.values())
    mono = all(vals[(grid[i],v)] <= vals[(grid[i+1],v)] + 1e-12 for i in range(20) for v in grid)
    assoc = all(St(gg,St(gg,u,v),z) == St(gg,u,St(gg,v,z)) for u in grid[::4] for v in grid[::4] for z in grid[::4])
    return inrange and mono and assoc and all(abs(St(gg,u,0)-u) < 1e-12 for u in grid)
check("Hamacher is a t-conorm (in [0,1], monotone, associative, identity 0) for gamma = 0 (corner by continuity), .5, 1, 2, 5, 50",
      all(is_tconorm(gg) for gg in [0, Fr(1,2), 1, 2, 5, 50]))
check("and fails for gamma = -0.5 and -2 (leaves [0,1] or blows up)", not is_tconorm(Fr(-1,2)) and not is_tconorm(-2))
def Tt(gg, u, v):
    d = gg + (1-gg)*(u + v - u*v)
    return u*v/d if d != 0 else 0.0
def two_increasing(gg):
    h = [i/30 for i in range(31)]
    return all(Tt(gg,h[i+1],h[j+1]) - Tt(gg,h[i],h[j+1]) - Tt(gg,h[i+1],h[j]) + Tt(gg,h[i],h[j]) >= -1e-12
               for i in range(30) for j in range(30))
check("the Hamacher t-norm is a copula (2-increasing) for gamma in [0,2] (checked 0, .5, 1, 1.5, 2)",
      all(two_increasing(gg) for gg in [0, .5, 1, 1.5, 2]))
check("but not for gamma = 3 or 10: a t-norm that is not a probability law", not two_increasing(3) and not two_increasing(10))
check("statistics page's AMH theta is 1 - gamma (same line, same points)",
      sp.simplify(T - a*b/(1 - (1-g)*(1-a)*(1-b))) == 0)

# other classical families: off the plane except where they meet Hamacher
mpf = mp.mpf
def off_plane(f):
    pts = [(mpf('0.2'), mpf('0.3')), (mpf('0.5'), mpf('0.6')), (mpf('0.7'), mpf('0.15')), (mpf('0.35'), mpf('0.8'))]
    # solve p,q exactly from the first two points, measure the miss at the others
    M = mp.matrix([[u*v, f(u,v)*u*v] for u,v in pts[:2]]); r = mp.matrix([f(u,v)-u-v for u,v in pts[:2]])
    p, q = mp.lu_solve(M, r)
    return max(abs(f(u,v) - (u+v+p*u*v)/(1-q*u*v)) for u,v in pts[2:])
frank = lambda s: (lambda u,v: 1 - mp.log(1 + (s**(1-u)-1)*(s**(1-v)-1)/(s-1))/mp.log(s))
yager = lambda wv: (lambda u,v: min(mpf(1), (u**wv + v**wv)**(1/mpf(wv))))
ss = lambda pv: (lambda u,v: 1 - ((1-u)**pv + (1-v)**pv - 1)**(1/mpf(pv)))
dombi = lambda l: (lambda u,v: 1 - 1/(1 + ((u/(1-u))**l + (v/(1-v))**l)**(1/mpf(l))))
aa = lambda l: (lambda u,v: 1 - mp.exp(-((-mp.log(1-u))**l + (-mp.log(1-v))**l)**(1/mpf(l))))
check("Dombi (lambda=1) and Schweizer-Sklar (p=-1) t-conorms equal Hamacher gamma=0",
      off_plane(dombi(1)) < mpf(10)**-30 and off_plane(ss(-1)) < mpf(10)**-30 and
      abs(dombi(1)(mpf('0.3'),mpf('0.6')) - St(0, mpf('0.3'), mpf('0.6'))) < mpf(10)**-30)
check("Aczel-Alsina (lambda=1) is the probabilistic sum", abs(aa(1)(mpf('0.3'),mpf('0.6')) - w(-1, mpf('0.3'), mpf('0.6'))) < mpf(10)**-30)
offs = {'Frank s=0.5': frank(mpf('0.5')), 'Frank s=2': frank(mpf(2)), 'Frank s=10': frank(mpf(10)),
        'Yager w=2': yager(2), 'Schweizer-Sklar p=-2': ss(-2), 'Dombi lambda=2': dombi(2), 'Aczel-Alsina lambda=2': aa(2)}
misses = {kname: off_plane(f) for kname, f in offs.items()}
for kname, mv in misses.items(): print("   miss", kname, mp.nstr(mv, 3))
check("Frank, Yager, Schweizer-Sklar, Dombi, Aczel-Alsina (other parameters) leave the landscape (miss > 1e-4)",
      all(mv > mpf('1e-4') for mv in misses.values()))
check("Frank's family hugs the line p + q = -1 without lying on it (fitted p + q within 0.05 of -1 for s = 0.5, 2)",
      all(abs(sum(mp.lu_solve(mp.matrix([[u*v, frank(s)(u,v)*u*v] for u,v in [(mpf('0.2'),mpf('0.3')),(mpf('0.5'),mpf('0.6'))]]),
                              mp.matrix([frank(s)(u,v)-u-v for u,v in [(mpf('0.2'),mpf('0.3')),(mpf('0.5'),mpf('0.6'))]]))) + 1) < .05
          for s in [mpf('0.5'), mpf(2)]))
check("Frank's equation T + S = a + b holds in Hamacher's family only at gamma = 1 (the product)",
      sp.factor(sp.simplify(T + S - a - b)).has(g - 1) and sp.simplify((T + S - a - b).subs(g, 1)) == 0
      and sp.simplify((T + S - a - b).subs({g: 2, a: sp.Rational(1,5), b: sp.Rational(1,3)})) != 0)
check("Frank s -> 1 gives the probabilistic sum (s = 1 + 1e-12 within 1e-9)",
      abs(frank(1 + mpf('1e-12'))(mpf('0.3'), mpf('0.6')) - w(-1, mpf('0.3'), mpf('0.6'))) < mpf('1e-9'))
phi = lambda t: mp.log(t/(1-t)); phinv = lambda z: 1/(1 + mp.exp(-z))
conj = lambda u, v: phinv(w(mpf('0.7'), phi(u), phi(v)))
check("the dial on another coordinate (log-odds) is again associative (conjugation keeps associativity)",
      all(abs(conj(conj(u,v),z) - conj(u,conj(v,z))) < mpf(10)**-30
          for u in [mpf('0.6'), mpf('0.8')] for v in [mpf('0.55'), mpf('0.9')] for z in [mpf('0.7')]))
luk = lambda u, v: min(1, u + v)
check("Lukasiewicz t-conorm min(a+b,1) is the k=0 dial clipped at 1, and not cancellative (0.6+0.6 = 0.6+0.7 = 1)",
      luk(.2,.3) == .2+.3 and luk(.6,.6) == luk(.6,.7) == 1)
check("max is the Frank limit s -> 0 (s = 1e-30 within 1e-3) and is idempotent, so no group law",
      abs(frank(mpf('1e-30'))(mpf('0.3'), mpf('0.7')) - mpf('0.7')) < 1e-3 and max(.4,.4) == .4)

# ---------- 5. logic of the rungs ----------
fexp = lambda u, v: (u+1)**(v+1) - 1
sadd = lambda u, v: u + v + 1
N = range(7)
check("tower law (x * a) * b = x * (a (+) b) for the fused exponent * (x,a,b < 7)",
      all(fexp(fexp(xx,aa_),bb) == fexp(xx, ma(aa_,bb)) for xx in N for aa_ in N for bb in N))
check("sum law (x * a) (+) (x * b) = x * (a + b + 1) (x,a,b < 7)",
      all(ma(fexp(xx,aa_), fexp(xx,bb)) == fexp(xx, sadd(aa_,bb)) for xx in N for aa_ in N for bb in N))
check("definability: b (+) c = d  iff  for all x, (x*b)*c = x*d  (b,c,d < 7, x < 6)",
      all((ma(bb,cc) == dd) == all(fexp(fexp(xx,bb),cc) == fexp(xx,dd) for xx in range(6)) for bb in N for cc in N for dd in range(60)))
check("definability: b + c + 1 = d  iff  for all x, (x*b)(+)(x*c) = x*d",
      all((sadd(bb,cc) == dd) == all(ma(fexp(xx,bb),fexp(xx,cc)) == fexp(xx,dd) for xx in range(6)) for bb in N for cc in N for dd in range(20)))
check("zero is the only right identity of the fused exponent (y < 20)",
      [yy for yy in range(20) if all(fexp(xx,yy) == xx for xx in range(6))] == [0])
check("mult-add alone is x in disguise: (a (+) b) + 1 = (a+1)(b+1) on Nat, lift one-to-one onto positives",
      all(ma(aa_,bb)+1 == (aa_+1)*(bb+1) for aa_ in range(30) for bb in range(30)))

n_pass = sum(ok for _, ok in results)
print(f"\n{n_pass}/{len(results)} checks passed")
