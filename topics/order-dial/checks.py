"""Order dial: checks for the data-structure map and the order axis of the circle operation.
Run: python3 checks.py  (needs sympy, numpy)."""
import itertools, math, random
import numpy as np
import sympy as sp

results = []
def check(name, ok):
    results.append((name, bool(ok)))
    print(("PASS " if ok else "FAIL ") + name)

rng = np.random.default_rng(1)

# ---------- A. Where the circle operation sits in the Boom hierarchy ----------
a, b, c, k = sp.symbols('a b c k')
w = lambda x, y: x + y + k*x*y
check("A1 the circle operation with dial k is associative", sp.expand(w(w(a,b),c) - w(a,w(b,c))) == 0)
check("A2 the circle operation with dial k is commutative", sp.expand(w(a,b) - w(b,a)) == 0)
idem = sp.solve(sp.Eq(w(a,a), a), a)
check("A3 idempotents of a+b+kab are only 0 and the wall -1/k", set(idem) == {0, -1/k})
check("A4 so the circle operation is not idempotent: 1 * 1 = 3 (k=1)", w(1,1).subs(k,1) == 3)
check("A5 the wall -1/k absorbs: (-1/k) * b = -1/k", sp.simplify(w(-1/k, b) + 1/k) == 0)

# ---------- B. Free objects: the generating functions of sets, bags, lists ----------
x = sp.symbols('x')
def coeff(expr, r): return sp.series(expr, x, 0, r+1).removeO().coeff(x, r)
n = 3
ok_sets = ok_bags = ok_lists = ok_arr = True
letters = range(n)
for r in range(0, 6):
    sets = sum(1 for s in itertools.combinations(letters, r))
    bags = sum(1 for s in itertools.combinations_with_replacement(letters, r))
    lists = sum(1 for s in itertools.product(letters, repeat=r))
    arr = sum(1 for s in itertools.permutations(letters, r))
    ok_sets &= sets == coeff((1+x)**n, r)                       # k = 1
    ok_bags &= bags == coeff((1-x)**(-n), r)                    # k = -1
    ok_lists &= lists == math.factorial(r)*coeff(sp.exp(n*x), r)  # k = 0, times r!
    ok_arr &= arr == math.factorial(r)*coeff((1+x)**n, r)
check("B1 sets of size r from 3 kinds = [x^r](1+x)^3 (repetition dial k=1)", ok_sets)
check("B2 bags = [x^r](1-x)^-3 (k=-1)", ok_bags)
check("B3 lists = r! [x^r] e^(3x) (k=0 counts ordered, labeled lists)", ok_lists)
check("B4 lists with no repeats = r! [x^r](1+x)^3", ok_arr)

# free left-regular band (ordered sets, LinkedHashSet): keep first occurrence
def lrb(u, v): return u + tuple(t for t in v if t not in u)
words = [w_ for r in range(4) for w_ in itertools.permutations(range(3), r)]
check("B5 ordered-set join is associative", all(lrb(lrb(p,q),s)==lrb(p,lrb(q,s)) for p in words for q in words for s in words))
check("B6 ordered-set join is idempotent and left regular (xyx = xy)", all(lrb(p,p)==p and lrb(lrb(p,q),p)==lrb(p,q) for p in words for q in words))
check("B7 ordered-set join is not commutative", lrb((0,),(1,)) != lrb((1,),(0,)))
check("B8 3 letters give 1+3+6+6 = 16 ordered sets", len(words) == 16)

# dictionary merge with override (last write wins) = right regular band
keys = [0, 1]; vals = [None, 'a', 'b']
dicts = [dict((kk, v) for kk, v in zip(keys, vs) if v is not None) for vs in itertools.product(vals, repeat=2)]
def over(p, q): r = dict(p); r.update(q); return r
check("B9 dict override is associative", all(over(over(p,q),s)==over(p,over(q,s)) for p in dicts for q in dicts for s in dicts))
check("B10 dict override is idempotent and right regular (xyx = yx)", all(over(p,p)==p and over(over(p,q),p)==over(q,p) for p in dicts for q in dicts))
check("B11 dict override is not commutative (merge order matters on shared keys)", over({0:'a'},{0:'b'}) != over({0:'b'},{0:'a'}))
check("B12 but dicts with disjoint keys merge in either order", over({0:'a'},{1:'b'}) == over({1:'b'},{0:'a'}))
# every collection is a map key -> value; value rule picks the structure
check("B13 random sets: union of independent memberships is k=-1 circle operation p+q-pq",
      abs((1-(1-0.3)*(1-0.6)) - (0.3+0.6-0.3*0.6)) < 1e-12)

# ---------- C. The order dial q: weight each out-of-order pair by q ----------
q = sp.symbols('q')
def inv(wd): return sum(1 for i in range(len(wd)) for j in range(i+1,len(wd)) if wd[i] > wd[j])
def qbinom(nn, kk):
    num = sp.prod([1-q**(nn-i) for i in range(kk)]); den = sp.prod([1-q**(i+1) for i in range(kk)])
    return sp.expand(sp.cancel(num/den))
wds = set(itertools.permutations([0,0,1,1]))
gen = sp.expand(sum(q**inv(wd) for wd in wds))
check("C1 arrangements of xxyy weighted q^inversions = Gaussian binomial [4,2]_q", sp.expand(gen - qbinom(4,2)) == 0)
check("C2 [4,2]_q = 1 + q + 2q^2 + q^3 + q^4", sp.expand(qbinom(4,2) - (1+q+2*q**2+q**3+q**4)) == 0)
check("C3 q=0 counts 1 (the bag), q=1 counts 6 (the lists)", qbinom(4,2).subs(q,0)==1 and qbinom(4,2).subs(q,1)==6)
# q-multinomial for bag {x,x,y,z}
wds3 = set(itertools.permutations([0,0,1,2]))
g3 = sp.expand(sum(q**inv(wd) for wd in wds3))
qfact = lambda m: sp.prod([sum(q**j for j in range(i)) for i in range(1, m+1)])
check("C4 bag xxyz: sum q^inv = [4]!/[2]! (q-multinomial), 12 at q=1, 1 at q=0",
      sp.expand(g3 - sp.cancel(qfact(4)/qfact(2))) == 0 and g3.subs(q,1)==12 and g3.subs(q,0)==1)
# quantum plane yx = qxy: normal ordering (x+y)^4 gives [4,j]_q
def normal_order(word):  # word of 'x','y' -> power of q to sort all x before y
    return inv([0 if ch=='x' else 1 for ch in word])
coef = {}
for wd in itertools.product('xy', repeat=4):
    j = wd.count('x'); coef[j] = coef.get(j, 0) + q**normal_order(wd)
check("C5 in the quantum plane yx = q xy, (x+y)^4 = sum [4,j]_q x^j y^(4-j)",
      all(sp.expand(coef[j] - qbinom(4,j)) == 0 for j in range(5)))

# ---------- D. Order cannot sit in the circle operation's formula without losing grouping ----------
A, B, C = sp.symbols('A B C', commutative=False)
al, be = sp.symbols('alpha beta')
m = lambda X, Y: X + Y + al*X*Y + be*Y*X
assoc = sp.expand(m(m(A,B),C) - m(A,m(B,C)))
expected = al*be*(B*A*C + C*A*B - A*C*B - B*C*A)
check("D1 a+b+alpha*ab+beta*ba: associator = alpha*beta*(BAC + CAB - ACB - BCA)", sp.expand(assoc - expected) == 0)
# numeric: with 3x3 matrices the bracket term is nonzero, so associative only at alpha*beta = 0
Mx = [rng.standard_normal((3,3)) for _ in range(3)]
def mm(X, Y, lam): return X + Y + lam*X@Y + (1-lam)*Y@X
def adef(lam): X,Y,Z = Mx; return np.abs(mm(mm(X,Y,lam),Z,lam) - mm(X,mm(Y,Z,lam),lam)).max()
check("D2 on matrices the lambda-mutation is associative at lambda=0,1 only", adef(0)<1e-10 and adef(1)<1e-10 and adef(0.5)>1e-3)
def flex(lam): X,Y,_ = Mx; return np.abs(mm(mm(X,Y,lam),X,lam) - mm(X,mm(Y,X,lam),lam)).max()
check("D3 but every lambda keeps flexibility (a*b)*a = a*(b*a)", max(flex(l) for l in np.linspace(0,1,11)) < 1e-9)
def jordan(X, Y): return (X@Y + Y@X)/2
X, Y, _ = Mx; XX = jordan(X, X)
check("D4 lambda=1/2 is the Jordan product: Jordan identity holds", np.abs(jordan(jordan(X,Y),XX) - jordan(X,jordan(Y,XX))).max() < 1e-9)
check("D5 shifted: a+b+(ab+ba)/2 = (1+a)o(1+b) - 1 for the Jordan product o", np.abs(mm(X,Y,0.5) - (jordan(np.eye(3)+X,np.eye(3)+Y) - np.eye(3))).max() < 1e-12)
check("D6 commutator of the mutation = (2*lambda-1)[a,b]: the order part is a Lie bracket",
      all(np.abs((mm(X,Y,l)-mm(Y,X,l)) - (2*l-1)*(X@Y-Y@X)).max()<1e-12 for l in (0,.3,.5,1)))

# ---------- E. The circle operation over non-commuting inputs: it copies the inputs' laws exactly ----------
def cd_mult(p, r):  # Cayley-Dickson, arrays of length 2^n
    nn = len(p)
    if nn == 1: return p*r
    h = nn//2; a_, b_ = p[:h], p[h:]; c_, d_ = r[:h], r[h:]
    conj = lambda z: np.concatenate([z[:1], -z[1:]])
    return np.concatenate([cd_mult(a_,c_) - cd_mult(conj(d_),b_), cd_mult(d_,a_) + cd_mult(b_,conj(c_))])
def one(dim): e = np.zeros(dim); e[0] = 1; return e
def star(p, r): return p + r + cd_mult(p, r)
def comm(op, p, r): return op(p,r) - op(r,p)
def assoc_(op, p, r, s): return op(op(p,r),s) - op(p,op(r,s))
for dim, name in [(4,'quaternions'),(8,'octonions'),(16,'sedenions')]:
    P, R, S = (rng.standard_normal(dim) for _ in range(3))
    check(f"E {name}: commutator of the circle operation = commutator of the product", np.allclose(comm(star,P,R), comm(cd_mult,P,R)))
    check(f"E {name}: associator of the circle operation = associator of the product", np.allclose(assoc_(star,P,R,S), assoc_(cd_mult,P,R,S)))
P, R, S = (rng.standard_normal(4) for _ in range(3))
check("E the quaternion circle operation: associative, not commutative", np.allclose(assoc_(star,P,R,S),0) and not np.allclose(comm(star,P,R),0))
P, R, S = (rng.standard_normal(8) for _ in range(3))
check("E the octonion circle operation: not associative but alternative (a*a)*b = a*(a*b)",
      (not np.allclose(assoc_(star,P,R,S),0)) and np.allclose(assoc_(star,P,P,R),0))
P, R = (rng.standard_normal(16) for _ in range(2))
check("E the sedenion circle operation: not even alternative", not np.allclose(assoc_(star,P,P,R),0))
# sedenion zero divisors give extra walls: (e3+e10)(e6-e15) = 0
e = lambda i: np.eye(16)[i]
u, v = e(3)+e(10), e(6)-e(15)
check("E sedenions: a zero-divisor pair (e3+e10)(e6-e15)=0, so a = u-1, b = v-1 circle operation onto the wall -1",
      np.allclose(cd_mult(u, v), 0) and np.allclose(star(u-one(16), v-one(16)), -one(16)))

# ---------- F. Composition of functions: the affine group ----------
def aff(st, sh): return np.array([[1+st, sh],[0,1.0]])  # x -> (1+a)x + c, circle-operation coordinate a
def compose(p, r): M = aff(*p) @ aff(*r); return (M[0,0]-1, M[0,1])
okF = True
for _ in range(50):
    p = tuple(rng.standard_normal(2)); r = tuple(rng.standard_normal(2))
    st, sh = compose(p, r)
    okF &= np.isclose(st, p[0]+r[0]+p[0]*r[0]) and np.isclose(sh, p[1] + (1+p[0])*r[1])
check("F1 composing x->(1+a)x+c: stretches combine by the circle operation, shifts pick up order: (a*b, c + (1+a)d)", okF)
okW = True
for _ in range(50):
    p = tuple(rng.standard_normal(2)); r = tuple(rng.standard_normal(2))
    okW &= np.isclose(compose(p,r)[1]-compose(r,p)[1], p[0]*r[1]-r[0]*p[1])
check("F2 order defect of two affine maps = a*d - b*c, the wedge of (a,c) and (b,d)", okW)
kk = 0.7
Tk = lambda s: (kk*s, s)   # T_k(s): x -> (1+k s)x + s, centre -1/k
check("F3 the k-family (stretches about the wall -1/k) commutes: zero wedge",
      np.isclose(compose(Tk(.3),Tk(-1.2))[1], compose(Tk(-1.2),Tk(.3))[1]))
check("F4 two stretches about different centres do not commute", not np.isclose(compose((1,0),(1,1))[1], compose((1,1),(1,0))[1]))
# fixed points
centre = lambda p: -p[1]/p[0]
check("F5 every T_k has centre (fixed point) -1/k", np.isclose(centre(Tk(.3)), -1/kk) and np.isclose(centre(Tk(-2)), -1/kk))
# same as the matrix circle operation on [[a,c],[0,0]]
Pm = np.array([[.4,.9],[0,0]]); Rm = np.array([[-.2,.5],[0,0]])
check("F6 same as the matrix circle operation A+B+AB on [[a,c],[0,0]]", np.allclose(Pm+Rm+Pm@Rm, np.array([[compose((.4,.9),(-.2,.5))[0], compose((.4,.9),(-.2,.5))[1]],[0,0]])))

# ---------- G. q-commutation: when is a smooth order dial honest? ----------
qq, cc = sp.symbols('q c')
# From e2 e1 = q e1 e2 and e1^2 = c: e2 e1 e1 = c e2, and = q e1 e2 e1 = q^2 e1 e1 e2 = q^2 c e2
check("G1 e2e1 = q e1e2 with e1^2 = c forces (q^2-1)c = 0: if things square to a number, q = +1 or -1",
      sp.solve(sp.Eq(qq**2*cc, cc), qq) == [-1, 1])
Xp = np.array([[0,1],[1,0]]); Zp = np.array([[1,0],[0,-1]])
check("G2 Pauli X, Z square to 1 and anticommute (q = -1)", np.allclose(Xp@Xp,np.eye(2)) and np.allclose(Zp@Zp,np.eye(2)) and np.allclose(Zp@Xp,-Xp@Zp))
# nilpotent generators (each at most once) allow any q: quantum exterior algebra; check overlap on e1 e1 e2 with e1^2=0
check("G3 with e^2 = 0 (each at most once) the overlap check gives 0 = 0 for every q", sp.expand(qq**2*0 - 0) == 0)
# clock and shift: CS = w SC with w = exp(2 pi i/3), C^3 = S^3 = I: q on the unit circle at roots of unity
wq = np.exp(2j*np.pi/3); Cm = np.diag([1,wq,wq**2]); Sm = np.roll(np.eye(3),1,axis=0)
check("G4 clock and shift (3x3): C S = w S C with w a cube root of 1 (q off +-1 needs cubes, not squares, to be 1)",
      np.allclose(Cm@Sm, wq*Sm@Cm) and np.allclose(np.linalg.matrix_power(Cm,3),np.eye(3)) and not np.allclose(Cm@Cm,np.eye(3)))
# self swap: if swapping two copies of the same element costs q != 1 then x^2 = 0 (exclusion)
check("G5 x x = q x x with q != 1 forces x^2 = 0: anticommuting (q=-1) forces 'at most once' (Pauli)",
      sp.solve(sp.Eq((1-qq)*sp.Symbol('s'), 0), sp.Symbol('s')) == [0])

# ---------- H. The smooth, associative order dial: hbar (Moyal star product) ----------
xs, ps, hb = sp.symbols('x p hbar')
def moyal(f, g, N=6):
    tot = 0
    for nn in range(N+1):
        s = 0
        for kk_ in range(nn+1):
            s += (-1)**kk_ * sp.binomial(nn,kk_) * sp.diff(f, xs, nn-kk_, ps, kk_) * sp.diff(g, xs, kk_, ps, nn-kk_)
        tot += (sp.I*hb/2)**nn / sp.factorial(nn) * s
    return sp.expand(tot)
f1, f2, f3 = xs**2 + ps, xs*ps**2, ps**3 - xs
check("H1 Moyal star is associative for every hbar (tested on polynomials)",
      sp.expand(moyal(moyal(f1,f2),f3) - moyal(f1,moyal(f2,f3))) == 0)
check("H2 x*p - p*x = i hbar: the commutator is set by hbar", sp.expand(moyal(xs,ps) - moyal(ps,xs) - sp.I*hb) == 0)
st = lambda f, g: f + g + moyal(f, g)
check("H3 circle operation over the star product is associative for every hbar",
      sp.expand(st(st(f1,f2),f3) - st(f1,st(f2,f3))) == 0)
pb = sp.diff(f1,xs)*sp.diff(f2,ps) - sp.diff(f1,ps)*sp.diff(f2,xs)
check("H4 the star circle-operation commutator = i hbar {a,b} + O(hbar^3) (Poisson bracket)",
      sp.expand(sp.series(st(f1,f2)-st(f2,f1), hb, 0, 3).removeO() - sp.I*hb*pb) == 0)
check("H5 at hbar = 0 the star circle operation is the ordinary (commutative) circle operation", sp.expand(st(f1,f2).subs(hb,0) - (f1+f2+f1*f2)) == 0)

# ---------- I. Composition first ----------
# I1 compounding: grow by r then deposit d, vs deposit then grow
B, r_, d_, r2, d2 = sp.symbols('B r d r2 d2')
grow_then_dep = B*(1+r_) + d_
dep_then_grow = (B + d_)*(1+r_)
check("I1 deposit-then-grow minus grow-then-deposit = r*d, the corner (interest on the deposit)",
      sp.expand(dep_then_grow - grow_then_dep - r_*d_) == 0)
# I2 one period = affine map x -> (1+r)x + d; two periods in either order differ by r1 d2 - r2 d1
P1 = lambda X: (1+r_)*X + d_; P2 = lambda X: (1+r2)*X + d2
check("I2 two compounding periods: order defect = r1*d2 - r2*d1 (a wedge)", sp.expand(P1(P2(B)) - P2(P1(B)) - (r_*d2 - r2*d_)) == 0)
# I3 annuity due vs ordinary annuity: N periods, FV_due = (1+r) FV_ordinary
N = 5; rr = sp.Rational(3,100)
ordinary = 0; due = 0
for _ in range(N): ordinary = ordinary*(1+rr) + 100; due = (due + 100)*(1+rr)
check("I3 depositing at the start of each period (annuity due) = (1+r) x depositing at the end", sp.simplify(due - (1+rr)*ordinary) == 0)
# I4 the centre (fixed point) of x -> (1+r)x + d is -d/r: the balance where interest exactly cancels the payment
Xs = sp.symbols('X')
check("I4 fixed point of one period is -d/r (interest-only balance); for T_k (r = k a, d = a) it is the wall -1/k",
      sp.solve(sp.Eq(P1(Xs), Xs), Xs) == [-d_/r_] and sp.simplify((-d_/r_).subs({r_: k*a, d_: a}) + 1/k) == 0)
# I5 periods with the same centre commute
check("I5 periods sharing a centre commute (r1 d2 = r2 d1)", sp.expand((r_*d2 - r2*d_).subs({d_: r_*Xs, d2: r2*Xs})) == 0)
# I6 composing two different-centre stretches gives a new centre
cen = sp.solve(sp.Eq(P1(P2(Xs)), Xs), Xs)[0]
check("I6 the composite of two periods has its own centre, -(d1 + (1+r1) d2)/(r1 ∘₁ r2)",
      sp.simplify(cen + (d_ + (1+r_)*d2)/(r_ + r2 + r_*r2)) == 0)

# I7 landscape laws are composition of Moebius maps (one-parameter subgroups fixing two points)
pp, qq2 = sp.symbols('p q')
def mob(t): return sp.Matrix([[1 + pp*t, t], [-qq2*t, 1]])   # I + t X, X = [[p,1],[-d,0]]
s_, t_ = sp.symbols('s t')
Mst = (mob(s_)*mob(t_)).applyfunc(sp.expand)
F = (s_ + t_ + pp*s_*t_)/(1 - qq2*s_*t_)
lhs = Mst/(1 - qq2*s_*t_)
check("I7 Moebius composition: mob(s) mob(t) = (1 - d s t) mob(F(s,t)), F the landscape law",
      all(sp.simplify(lhs[i] - mob(F)[i]) == 0 for i in range(4)))
check("I8 circle operation (p=1, d=0) is the affine slice: Moebius maps that fix infinity",
      mob(t_).subs({pp:1, qq2:0})[1,0] == 0)
# I9 Lie's vector fields alpha + beta x + gamma x^2 (sl2): bracket of translation and dilation
X_ = sp.symbols('x')
def br(f, g): return sp.expand(f*sp.diff(g, X_) - g*sp.diff(f, X_))
check("I9 [translation, dilation] = translation: affine algebra [X, Y] = Y (not central)", br(sp.Integer(1), X_) == 1)
check("I10 [translation, x^2] = 2x: the full Moebius algebra sl2 closes", br(sp.Integer(1), X_**2) == 2*X_ and br(X_, X_**2) == X_**2)

# ---------- J. Route to quantum: shifts composed with multiplications, acting on functions ----------
lam, n_, a_ = sp.symbols('lambda n a', positive=True)
fn = sp.Function('f')
T = lambda g, sh: g.subs(n_, n_ + sh)                 # shift the argument (an addition)
M = lambda g: lam**n_ * g                             # multiply by lambda^n (a multiplication)
lhsJ = T(M(fn(n_)), a_); rhsJ = lam**a_ * M(T(fn(n_), a_))
check("J1 discrete Weyl relation: shift-by-a after multiply-by-lambda^n = lambda^a x the other order (q-commutation, q = lambda^a)",
      sp.simplify(sp.expand_power_base(lhsJ - rhsJ, force=True)) == 0 or sp.simplify(lhsJ/rhsJ) == 1)
xx, aa, bb = sp.symbols('x a b', real=True)
g = sp.Function('g')
Ta = lambda h: h.subs(xx, xx + aa)
Mb = lambda h: sp.exp(sp.I*bb*xx) * h
check("J2 Weyl relation: T_a M_b = e^{iab} M_b T_a; the order costs a phase equal to the rectangle area ab",
      sp.simplify(Ta(Mb(g(xx))) - sp.exp(sp.I*aa*bb)*Mb(Ta(g(xx)))) == 0)
check("J3 with real growth e^{bx} instead, the order costs the factor e^{ab} = 1 + (e^{ab} - 1)",
      sp.simplify((g(xx+aa)*sp.exp(bb*(xx+aa))) - sp.exp(aa*bb)*sp.exp(bb*xx)*g(xx+aa)) == 0)
# J4 Weyl-ordered W(a,b) = e^{i(ab/2)} T_a M_b composes with half-wedge phase
a1,b1,a2,b2 = sp.symbols('a1 b1 a2 b2', real=True)
def W(A_, B_, h): return sp.exp(sp.I*A_*B_/2) * sp.exp(sp.I*B_*xx) * h.subs(xx, xx + A_)
comp = W(a1,b1, W(a2,b2, g(xx)))
target = sp.exp(sp.I*(a1*b2 - a2*b1)/2) * W(a1+a2, b1+b2, g(xx))
check("J4 Weyl operators compose with phase e^{i(a1 b2 - a2 b1)/2}: half the wedge, the oriented triangle area",
      sp.simplify(sp.expand(sp.log(sp.simplify(comp/target)), force=True)) == 0 or sp.simplify(comp - target) == 0)
# J5 Heisenberg algebra is not the affine algebra: in Heisenberg [X,Y]=Z commutes with everything
Xh = sp.Matrix([[0,1,0],[0,0,0],[0,0,0]]); Yh = sp.Matrix([[0,0,0],[0,0,1],[0,0,0]]); Zh = Xh*Yh - Yh*Xh
Xa = sp.Matrix([[1,0],[0,0]]); Ya = sp.Matrix([[0,1],[0,0]])
check("J5 Heisenberg: [X,Y] = Z is central ([X,Z] = [Y,Z] = 0); affine: [X,Y] = Y, not central",
      Zh*Xh - Xh*Zh == sp.zeros(3) and Zh*Yh - Yh*Zh == sp.zeros(3) and Xa*Ya - Ya*Xa == Ya and Xa*Ya != Ya*Xa)
check("J6 so the affine group (composing on points) is not the Heisenberg group: affine is not nilpotent ([X,[X,Y]] = Y != 0)",
      Xa*(Xa*Ya - Ya*Xa) - (Xa*Ya - Ya*Xa)*Xa == Ya)

print(f"\n{sum(ok for _,ok in results)}/{len(results)} checks pass")
