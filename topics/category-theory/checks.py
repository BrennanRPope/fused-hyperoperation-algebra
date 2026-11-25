"""Checks for the category theory page. Every claim on index.html is tested here."""
import itertools
from itertools import zip_longest, product
import numpy as np
import sympy as sp

results = []
def check(name, ok):
    results.append((name, bool(ok)))
    print(("PASS" if ok else "FAIL"), "-", name)

ma = lambda a, b: a + b + a*b
NONE = None
# These values as tagged tuples
This = lambda a: ('this', a); That = lambda b: ('that', b); Both = lambda a, b: ('both', a, b)
def these(A, B): return [This(a) for a in A] + [That(b) for b in B] + [Both(a, b) for a in A for b in B]
def opt(A): return [NONE] + [('some', a) for a in A]

# ---------- 1. counting ----------
check("|These A B| = a + b + ab = a (+) b for a, b < 8",
      all(len(these(range(a), range(b))) == ma(a, b) for a in range(8) for b in range(8)))
def to_pair(x):
    if x is NONE: return (NONE, NONE)
    x = x[1]
    if x[0] == 'this': return (('some', x[1]), NONE)
    if x[0] == 'that': return (NONE, ('some', x[1]))
    return (('some', x[1]), ('some', x[2]))
def of_pair(p):
    l, r = p
    if l is NONE and r is NONE: return NONE
    if r is NONE: return ('some', This(l[1]))
    if l is NONE: return ('some', That(r[1]))
    return ('some', Both(l[1], r[1]))
ok = True
for a in range(5):
    for b in range(5):
        A, B = range(a), range(b)
        dom = opt(these(A, B)); cod = list(product(opt(A), opt(B)))
        ok &= sorted(map(str, map(to_pair, dom))) == sorted(map(str, cod))
        ok &= all(of_pair(to_pair(x)) == x for x in dom) and all(to_pair(of_pair(p)) == p for p in cod)
check("Option(These A B) <-> Option A x Option B is a bijection with both round trips (a, b < 5)", ok)
check("(a+1)(b+1) = (a (+) b) + 1 and = 1 + (a+b) + ab: base point, wedge, smash",
      all((a+1)*(b+1) == ma(a,b) + 1 == 1 + (a+b) + a*b for a in range(30) for b in range(30)))
ok = True
for a in range(4):
    for b in range(4):
        pts = list(product(['*'] + list(range(a)), ['*'] + list(range(b))))
        nonbase = [p for p in pts if p != ('*', '*')]
        wedge = [p for p in nonbase if '*' in p]; smash = [p for p in nonbase if '*' not in p]
        ok &= len(nonbase) == ma(a, b) and len(wedge) == a + b and len(smash) == a*b
check("pointed sets: A+ x B+ has a (+) b non-base points = wedge (a+b) + smash (ab) (a, b < 4)", ok)
check("Option of the empty set is the one-point set (unit goes to unit)", len(opt([])) == 1)

# ---------- 2. symmetric monoidal structure ----------
def assoc(t):
    if t[0] == 'this':
        u = t[1]
        if u[0] == 'this': return This(u[1])
        if u[0] == 'that': return That(This(u[1]))
        return Both(u[1], This(u[2]))
    if t[0] == 'that': return That(That(t[1]))
    u, c = t[1], t[2]
    if u[0] == 'this': return Both(u[1], That(c))
    if u[0] == 'that': return That(Both(u[1], c))
    return Both(u[1], Both(u[2], c))
def tmap(f, g, t):
    if t[0] == 'this': return This(f(t[1]))
    if t[0] == 'that': return That(g(t[1]))
    return Both(f(t[1]), g(t[2]))
def swap(t):
    if t[0] == 'this': return That(t[1])
    if t[0] == 'that': return This(t[1])
    return Both(t[2], t[1])
idf = lambda x: x
A, B, C, D = ['a1', 'a2'], ['b1'], ['c1', 'c2'], ['d1']
L = these(these(A, B), C)
R = these(A, these(B, C))
check("associator is a bijection (A(+)B)(+)C -> A(+)(B(+)C)", sorted(map(str, map(assoc, L))) == sorted(map(str, R)) and len(set(map(str, map(assoc, L)))) == len(L))
check("both sides have the cube's 7 corners: a+b+c+ab+ac+bc+abc",
      all(len(these(these(range(a),range(b)),range(c))) == a+b+c+a*b+a*c+b*c+a*b*c for a in range(4) for b in range(4) for c in range(4)))
P = these(these(these(A, B), C), D)
pent = all(assoc(assoc(t)) == tmap(idf, assoc, assoc(tmap(assoc, idf, t))) for t in P)
check("pentagon: the two ways to re-bracket ((A B) C) D agree on every element", pent)
hexa = all(assoc(swap(assoc(t))) == tmap(idf, swap, assoc(tmap(swap, idf, t))) for t in L)
check("hexagon: the braiding (swap) and associator agree both ways round", hexa)
check("swap is its own inverse", all(swap(swap(t)) == t for t in these(A, C)))
unit = these([], B)
check("unit: empty (+) B has exactly the b elements That(b)", unit == [That(b) for b in B])


# tagged dial: 1 + K x (A + B + K x A x B) <-> (1 + K x A) x (1 + K x B)
ok = True
for kk in range(4):
    for a in range(3):
        for b in range(3):
            K, A_, B_ = range(kk), range(a), range(b)
            TK = [('this', x) for x in A_] + [('that', y) for y in B_] + [('both', j, x, y) for j in K for x in A_ for y in B_]
            lhs = [None] + [(i, tt) for i in K for tt in TK]
            def tp(z):
                if z is None: return (None, None)
                i, tt = z
                if tt[0] == 'this': return ((i, tt[1]), None)
                if tt[0] == 'that': return (None, (i, tt[1]))
                return ((i, tt[2]), (tt[1], tt[3]))
            rhs = list(product([None] + [(i, x) for i in K for x in A_], [None] + [(j, y) for j in K for y in B_]))
            img = [tp(z) for z in lhs]
            ok &= len(set(img)) == len(lhs) and sorted(map(str, img)) == sorted(map(str, rhs))
            ok &= len(TK) == a + b + kk*a*b
check("tagged dial: 1 + K(A (+)k B) <-> (1 + KA)(1 + KB) is a bijection for k < 4 tags, a, b < 3 (the lift 1 + kx on sets)", ok)

# ---------- 3. product in Par ----------
def partial_maps(X, Y):
    return [dict(zip(X, vals)) for vals in product([NONE] + list(Y), repeat=len(X))]
def fst(t): return t[1] if t[0] in ('this', 'both') else NONE
def snd(t): return t[1] if t[0] == 'that' else (t[2] if t[0] == 'both' else NONE)
def kc(h, k): return {x: (NONE if h[x] is NONE else k(h[x])) for x in h}
X = ['x1', 'x2']; A2, B1 = ['a1', 'a2'], ['b1']
T = these(A2, B1)
ok = True
for f in partial_maps(X, A2):
    for g in partial_maps(X, B1):
        matches = [h for h in partial_maps(X, T) if kc(h, fst) == f and kc(h, snd) == g]
        ok &= len(matches) == 1
check("in Par, for every pair of partial maps f: X -> A, g: X -> B exactly one h: X -> A(+)B has projections f, g", ok)
check("total number of such h equals |partial X->A| x |partial X->B| = (a+1)^x (b+1)^x",
      len(partial_maps(X, T)) == (len(A2)+1)**2 * (len(B1)+1)**2)
check("the empty set is a zero object in Par: one partial map into it, one out of it",
      len(partial_maps(X, [])) == 1 and len(partial_maps([], A2)) == 1)
check("disjoint union is the coproduct in Par too (counts: partial maps out of A+B = product of partial maps out of each)",
      len(partial_maps(['a1','a2','b1'], ['y1','y2'])) == len(partial_maps(['a1','a2'], ['y1','y2'])) * len(partial_maps(['b1'], ['y1','y2'])))

# ---------- 4. align / zipLongest / outer join ----------
def align(xs, ys):
    out = []
    for i in range(max(len(xs), len(ys))):
        if i < len(xs) and i < len(ys): out.append(Both(xs[i], ys[i]))
        elif i < len(xs): out.append(This(xs[i]))
        else: out.append(That(ys[i]))
    return out
ok = True
for n in range(6):
    for m in range(6):
        xs, ys = list(range(n)), [chr(97+j) for j in range(m)]
        al = align(xs, ys)
        ok &= len(al) == max(n, m)
        ok &= [fst(t) for t in al if fst(t) is not NONE] == xs and [snd(t) for t in al if snd(t) is not NONE] == ys
        ok &= [(t[1], t[2]) for t in al if t[0] == 'both'] == list(zip(xs, ys))
        ok &= [(fst(t) if fst(t) is not NONE else '-', snd(t) if snd(t) is not NONE else '-') for t in al] == list(zip_longest(xs, ys, fillvalue='-'))
check("align: length max(n,m), both lists come back, zip is the both-part, zip_longest is align with defaults", ok)
old = {'alice': 3, 'bob': 5, 'cara': 2}; new = {'bob': 6, 'cara': 2, 'dev': 1}
join = {k: (Both(old[k], new[k]) if k in old and k in new else This(old[k]) if k in old else That(new[k])) for k in sorted(old.keys() | new.keys())}
check("full outer join by key: each row is a These (left only, right only, both); row count = |keys A union keys B|",
      join == {'alice': This(3), 'bob': Both(5, 6), 'cara': Both(2, 2), 'dev': That(1)} and len(join) == len(old.keys() | new.keys()))

# ---------- 5. decategorification ----------
t = sp.symbols('t')
def P(cells): return sum(c*t**i for i, c in enumerate(cells))
def pmul(X, Y):
    out = [0]*(len(X)+len(Y)-1)
    for i, x in enumerate(X):
        for j, y in enumerate(Y): out[i+j] += x*y
    return out
sphere = lambda n: [2] if n == 0 else [1] + [0]*(n-1) + [1]
import random; random.seed(3)
ok = True
for _ in range(40):
    Xc = [1] + [random.randint(0, 3) for _ in range(3)]; Yc = [1] + [random.randint(0, 3) for _ in range(2)]
    ok &= sp.expand(P(pmul(Xc, Yc)) - 1 - ma(P(Xc) - 1, P(Yc) - 1)) == 0
check("cell polynomials multiply under products, so P - 1 combines by the circle operation, as polynomials in t (random complexes)", ok)
chit = lambda c: P(c).subs(t, -1) - 1
check("spheres: reduced Euler characteristic of S^n is (-1)^n (n < 8)", all(chit(sphere(n)) == (-1)**n for n in range(8)))
check("circle is the wall: reduced chi of S^1 x Y is -1 for Y = S^0..S^5 and torus", all(chit(pmul(sphere(1), sphere(n))) == -1 for n in range(6)) and chit(pmul(sphere(1), [1,2,1])) == -1)
check("torus S^1 x S^1 has cells 1, 2, 1; S^2 x S^2 has reduced chi 1 (+) 1 = 3",
      pmul(sphere(1), sphere(1)) == [1, 2, 1] and chit(pmul(sphere(2), sphere(2))) == 3 == ma(1, 1))
check("at t = 1: non-base cells of a product = wedge (a+b) + smash (ab) of non-base cells",
      all(sum(pmul(Xc, Yc)) - 1 == ma(sum(Xc)-1, sum(Yc)-1) for Xc in ([1,a,b] for a in range(3) for b in range(3)) for Yc in ([1,c] for c in range(4))))
# real homology: triangulated torus (3x3 grid) and its suspension
def torus_triangles(n=3):
    v = lambda i, j: (i % n)*n + (j % n)
    tri = []
    for i in range(n):
        for j in range(n):
            tri += [tuple(sorted((v(i,j), v(i+1,j), v(i+1,j+1)))), tuple(sorted((v(i,j), v(i,j+1), v(i+1,j+1))))]
    return tri
def closure(top):
    S = set()
    for s in top:
        for r in range(1, len(s)+1):
            S |= set(itertools.combinations(s, r))
    return S
def betti(simplices):
    by = {}
    for s in simplices: by.setdefault(len(s)-1, []).append(s)
    for d in by: by[d].sort()
    idx = {d: {s: i for i, s in enumerate(by[d])} for d in by}
    rank = {}
    for d in by:
        if d == 0: rank[d] = 0; continue
        M = np.zeros((len(by[d-1]), len(by[d])))
        for j, s in enumerate(by[d]):
            for k in range(len(s)):
                M[idx[d-1][s[:k] + s[k+1:]], j] = (-1)**k
        rank[d] = np.linalg.matrix_rank(M)
    top = max(by)
    return [len(by[d]) - rank[d] - (rank[d+1] if d+1 in rank else 0) for d in range(top+1)]
Tor = closure(torus_triangles())
check("simplicial homology of a triangulated torus: Betti numbers 1, 2, 1 = coefficients of (1+t)^2", betti(Tor) == [1, 2, 1])
top = [s for s in Tor if len(s) == 3]
susp = closure([s + (100,) for s in top] + [s + (101,) for s in top])
b = betti(susp)
check("unreduced suspension of the torus: Betti 1, 0, 2, 1, i.e. reduced homology of S^2 v S^2 v S^3 = Sigma S^1 v Sigma S^1 v Sigma(S^1 ^ S^1)", b == [1, 0, 2, 1])
check("suspension multiplies the reduced polynomial by t: t((1+t)^2 - 1) = 2t^2 + t^3, the homology of S^2 v S^2 v S^3",
      sp.expand(t*((1+t)**2 - 1) - (t*t + t*t + t*(t*t))) == 0)
check("suspension splitting count: t(a (+) b) = ta + tb + t ab", sp.expand(t*ma(sp.Symbol('a'), sp.Symbol('b')) - (t*sp.Symbol('a') + t*sp.Symbol('b') + t*sp.Symbol('a')*sp.Symbol('b'))) == 0)
check("cardinality is a monoidal functor: |A (+) B| = |A| (+) |B| and |Option A| = |A| + 1 is the lift",
      all(len(these(range(a), range(b))) == ma(a, b) and len(opt(range(a))) == a + 1 for a in range(6) for b in range(6)))

wedge = lambda X_, Y_: [X_[0] + Y_[0] - 1] + [ (X_[i] if i < len(X_) else 0) + (Y_[i] if i < len(Y_) else 0) for i in range(1, max(len(X_), len(Y_)))]
check("a point has reduced chi 0 (the identity); S^1 v S^1 has -2; wedge adds reduced chi",
      chit([1]) == 0 and chit(wedge(sphere(1), sphere(1))) == -2 and chit(wedge(sphere(2), [1,2,1])) == chit(sphere(2)) + chit([1,2,1]))
n_pass = sum(ok for _, ok in results)
print(f"\n{n_pass}/{len(results)} checks passed")
