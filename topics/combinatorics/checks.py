"""Combinatorics through mult-add: every claim on index.html, checked.
Run: python3 checks.py > checks_output.txt"""
import itertools, math
from fractions import Fraction as F
import sympy as sp
import networkx as nx

R = []
def check(name, ok):
    R.append(bool(ok)); print(('PASS' if ok else 'FAIL') + '  ' + name)

x, t, s, k, y = sp.symbols('x t s k y')
fact = math.factorial
C = math.comb

# ---------- 1. The cube is the generating function of subsets ----------
for n in range(1, 6):
    xs = sp.symbols('x1:%d' % (n + 1))
    lhs = sp.expand(sp.prod([1 + xi for xi in xs]) - 1)
    rhs = sum(sp.prod(S) for r in range(1, n + 1) for S in itertools.combinations(xs, r))
    check(f'cube: prod(1+x_i)-1 = sum over nonempty subsets, n={n}', sp.expand(lhs - rhs) == 0)
    # fold of mult-add equals the same
    acc = 0
    for xi in xs: acc = acc + xi + acc * xi
    check(f'fold x1 * x2 * ... (mult-add) = cube, n={n}', sp.expand(acc - lhs) == 0)
check('at x_i = 1: 2^n - 1 nonempty subsets, n=1..10',
      all(sum(1 for r in range(1, n + 1) for _ in itertools.combinations(range(n), r)) == 2**n - 1 for n in range(1, 11)))

# ---------- 2. The dial: k^|S| ----------
for n in range(1, 5):
    xs = sp.symbols('x1:%d' % (n + 1))
    lhs = sp.expand(sp.prod([1 + k * xi for xi in xs]))
    rhs = sum(k**len(S) * sp.prod(S) for r in range(0, n + 1) for S in itertools.combinations(xs, r))
    check(f'dial: prod(1+k x_i) = sum k^|S| x^S, n={n}', sp.expand(lhs - rhs) == 0)
# subsets with each chosen element painted one of k colours = (1+k)^n = q^n
ok = True
for n in range(0, 6):
    for kk in range(0, 4):
        cnt = sum(1 for w in itertools.product(range(kk + 1), repeat=n))  # 0 = off, 1..k = on with colour
        cnt2 = sum(C(n, j) * kk**j for j in range(n + 1))
        ok &= (cnt == cnt2 == (kk + 1)**n)
check('painted subsets: sum_j C(n,j) k^j = (1+k)^n = q^n (brute force n<=5, k<=3)', ok)
check('k = -1: alternating subset sum is 0 for n>=1 (n=1..12)', all(sum((-1)**j * C(n, j) for j in range(n + 1)) == 0 for n in range(1, 13)))
check('k = -2: sum (-2)^|S| = (-1)^n (n=0..12)', all(sum((-2)**j * C(n, j) for j in range(n + 1)) == (-1)**n for n in range(13)))
# derangements by k = -1 on fixed points
def derange_bf(n): return sum(1 for p in itertools.permutations(range(n)) if all(p[i] != i for i in range(n)))
Dvals = [sum((-1)**j * C(n, j) * fact(n - j) for j in range(n + 1)) for n in range(9)]
check('derangements = sum_j (-1)^j C(n,j)(n-j)!  (brute force n<=7)', all(Dvals[n] == derange_bf(n) for n in range(8)))
check('derangements 1,0,1,2,9,44,265,1854,14833', Dvals == [1, 0, 1, 2, 9, 44, 265, 1854, 14833])

# ---------- 3. n! as a mult-add fold; layers are Stirling cycle numbers ----------
ma = lambda a, b: a + b + a * b
ok = True
for n in range(2, 13):
    acc = 0
    for i in range(1, n): acc = ma(acc, i)
    ok &= (acc == fact(n) - 1)
check('1 * 2 * ... * (n-1) = n! - 1 under mult-add (n=2..12)', ok)
def cycles(p):
    seen = set(); c = 0
    for i in range(len(p)):
        if i not in seen:
            c += 1; j = i
            while j not in seen: seen.add(j); j = p[j]
    return c
ok = True
for n in range(1, 8):
    cnt = [0] * (n + 1)
    for p in itertools.permutations(range(n)): cnt[cycles(p)] += 1
    for j in range(n):  # layer j = e_j(1..n-1) = permutations with n-j cycles
        ej = sum(math.prod(S) for S in itertools.combinations(range(1, n), j))
        ok &= (ej == cnt[n - j])
check('layer j of the cube with axis weights 1..n-1 = #perms with n-j cycles (n<=7)', ok)
check('row n=4 layers 1,6,11,6', [sum(math.prod(S) for S in itertools.combinations(range(1, 4), j)) for j in range(4)] == [1, 6, 11, 6])
ok = True
for n in range(1, 8):
    for kk in (-1, 2, 3):
        lhs = sum(kk**(n - cycles(p)) for p in itertools.permutations(range(n)))
        ok &= (lhs == math.prod(1 + kk * i for i in range(1, n)))
check('dial: sum over perms k^(n - cycles) = prod_{i<n}(1 + k i), k=-1,2,3, n<=7', ok)
check('k = -1: even and odd permutations balance (sum of signs 0, n=2..7)',
      all(sum((-1)**(n - cycles(p)) for p in itertools.permutations(range(n))) == 0 for n in range(2, 8)))

# ---------- 4. Labelled structures: (1+A)(1+B) - 1 is "an A-part, a B-part, or both" ----------
N = 9
def egf(seq): return sum(sp.Rational(seq(n), fact(n)) * x**n for n in range(1, N + 1))
def coeffs(f): 
    sr = sp.series(f, x, 0, N + 1).removeO()
    return [sp.expand(sr).coeff(x, n) * fact(n) for n in range(N + 1)]
A = egf(lambda n: fact(n))   # nonempty linear orders
B = egf(lambda n: 1)         # nonempty sets
cf = coeffs((1 + A) * (1 + B) - 1)
ok = True
for n in range(1, 7):
    cnt = 0
    for mask in range(2**n):
        S1 = [i for i in range(n) if mask >> i & 1]; S2 = [i for i in range(n) if not mask >> i & 1]
        cnt += fact(len(S1)) * 1  # S1: empty or an order; S2: empty or a set (one way)
    ok &= (cnt == cf[n])
check('EGF mult-add (1+A)(1+B)-1 counts splits into optional A-part and B-part (n<=6)', ok)
check('two-colourings: (e^x-1) * (e^x-1) = e^{2x}-1, counts 2^n', all(c == 2**n for n, c in enumerate(coeffs(sp.exp(2*x)-1)) if n >= 1)
      and sp.simplify((sp.exp(x)) * (sp.exp(x)) - 1 - (sp.exp(2*x) - 1)) == 0)

# ---------- 5. Exponential formula = E applied to connected pieces ----------
Ex = lambda f: sp.exp(f) - 1
Gall = egf(lambda n: 2**C(n, 2))
Ccon = sp.series(sp.log(1 + Gall), x, 0, N + 1).removeO()
con = [sp.expand(Ccon).coeff(x, n) * fact(n) for n in range(N + 1)]
check('connected graphs = log(1+G): 1,1,4,38,728,26704,1866256', con[1:8] == [1, 1, 4, 38, 728, 26704, 1866256])
def con_bf(n):
    E_all = list(itertools.combinations(range(n), 2)); c = 0
    for mask in range(2**len(E_all)):
        g = nx.Graph(); g.add_nodes_from(range(n)); g.add_edges_from(e for i, e in enumerate(E_all) if mask >> i & 1)
        c += nx.is_connected(g)
    return c
check('connected graph counts by brute force n=1..5', all(con_bf(n) == con[n] for n in range(1, 6)))
back = coeffs(Ex(Ccon))
check('E(connected) = all nonempty graphs 2^C(n,2) (n<=9)', all(back[n] == 2**C(n, 2) for n in range(1, N + 1)))
# rooted recurrence (what Lean computes): g_n = sum_j C(n-1,j-1) c_j g_{n-j}, g_0 = 1
g = [1] + [2**C(n, 2) for n in range(1, N + 1)]
check('recurrence g_n = sum C(n-1,j-1) c_j g_(n-j) (n<=9)', all(g[n] == sum(C(n - 1, j - 1) * con[j] * g[n - j] for j in range(1, n + 1)) for n in range(1, N + 1)))

# ---------- 6. Rung flow on x: integer rungs count nested partitions ----------
def bell_bf(n):
    def parts(l):
        if not l: yield []; return
        a, rest = l[0], l[1:]
        for p in parts(rest):
            for i in range(len(p)): yield p[:i] + [[a] + p[i]] + p[i + 1:]
            yield [[a]] + p
    return sum(1 for _ in parts(list(range(n))))
E1 = coeffs(Ex(x)); E2 = coeffs(Ex(Ex(x))); E3 = coeffs(Ex(Ex(Ex(x))))
Bell = [bell_bf(n) for n in range(N + 1)]
check('E^1(x) = e^x - 1: one nonempty set on n points (n<=9)', all(E1[n] == 1 for n in range(1, N + 1)))
check('E^2(x): Bell numbers 1,2,5,15,52,203,877,4140,21147 (brute force)', all(E2[n] == Bell[n] for n in range(1, N + 1)))
def S2(n, j): return sum((-1)**i * C(j, i) * (j - i)**n for i in range(j + 1)) // fact(j)
nested3 = [sum(S2(n, j) * Bell[j] for j in range(0, n + 1)) for n in range(N + 1)]
check('E^3(x): partitions of partitions 1,3,12,60,358,2471,19302,167894 (sum_j S(n,j) B_j)', E3[1:9] == [1, 3, 12, 60, 358, 2471, 19302, 167894] and all(E3[n] == nested3[n] for n in range(1, N + 1)))
Lm1 = coeffs(sp.log(1 + x))
check('E^-1(x) = log(1+x): signed cyclic orders (-1)^(n-1)(n-1)! (n<=9)', all(Lm1[n] == (-1)**(n - 1) * fact(n - 1) for n in range(1, N + 1)))
# Moebius function of the partition lattice: mu(0hat,1hat) = (-1)^(n-1)(n-1)!  (brute force n<=5)
def all_parts(l):
    if not l: yield []; return
    a, rest = l[0], l[1:]
    for p in all_parts(rest):
        for i in range(len(p)): yield p[:i] + [[a] + p[i]] + p[i + 1:]
        yield [[a]] + p
def mu_top(n):
    ps = [frozenset(frozenset(b) for b in p) for p in all_parts(list(range(n)))]
    finer = lambda s, p: all(any(b <= c for c in p) for b in s)
    ps.sort(key=lambda p: -len(p)); mu = {}
    for p in ps: mu[p] = 1 if len(p) == n else -sum(mu[s] for s in mu if s != p and finer(s, p))
    return mu[frozenset([frozenset(range(n))])]
check('mu(0,1) of the partition lattice = (-1)^(n-1)(n-1)! = n![x^n] log(1+x)  (n<=5)', all(mu_top(n) == Lm1[n] for n in range(1, 6)))
# flow polynomials P_n(t) = n! [x^n] E^t(x), via the Carleman matrix
Nf = 8
def mul(a, b):
    c = [F(0)] * (Nf + 1)
    for i, ai in enumerate(a):
        for j, bj in enumerate(b):
            if i + j <= Nf: c[i + j] += ai * bj
    return c
Eser = [F(0)] + [F(1, fact(i)) for i in range(1, Nf + 1)]
M = sp.zeros(Nf, Nf); pw = [F(1)] + [F(0)] * Nf
for j in range(1, Nf + 1):
    pw = mul(pw, Eser)
    for i in range(1, Nf + 1): M[i - 1, j - 1] = sp.Rational(pw[i].numerator, pw[i].denominator)
Nm = M - sp.eye(Nf); L = sp.zeros(Nf, Nf); P = sp.eye(Nf)
for m in range(1, Nf + 1): P = P * Nm; L += ((-1)**(m + 1)) * P / m
v = sp.zeros(Nf, 1); v[0] = 1; cur = v; res = v.copy()
for m in range(1, Nf + 1): cur = L * cur / m; res += cur * t**m
Pn = [None] + [sp.expand(res[i] * fact(i + 1)) for i in range(Nf)]
check('P_n(t) at t = 1,2,3 gives 1, Bell, partitions of partitions (n<=8)', all(Pn[n].subs(t, 1) == 1 and Pn[n].subs(t, 2) == Bell[n] and Pn[n].subs(t, 3) == nested3[n] for n in range(1, Nf + 1)))
check('P_n(t) at t = 0 is x, at t = -1 the signed cyclic orders, at t = -2 log(1+log(1+x))',
      all(Pn[n].subs(t, 0) == (1 if n == 1 else 0) for n in range(1, Nf + 1)) and all(Pn[n].subs(t, -1) == Lm1[n] for n in range(1, Nf + 1))
      and all(Pn[n].subs(t, -2) == c for n, c in enumerate(coeffs(sp.log(1 + sp.log(1 + x)))) if 1 <= n <= Nf))
check('P_2 = t, P_3 = t(3t-1)/2, P_4 = t(2t-1)(3t-1)/2', sp.expand(Pn[2] - t) == 0 and sp.expand(Pn[3] - t*(3*t - 1)/2) == 0 and sp.expand(Pn[4] - t*(2*t - 1)*(3*t - 1)/2) == 0)
check('P_3 has roots 0, 1/3 and P_4 roots 0, 1/3, 1/2 (leading terms positive): negative on (0,1/3) and (1/3,1/2)',
      set(sp.solve(Pn[3], t)) == {0, sp.Rational(1, 3)} and set(sp.solve(Pn[4], t)) == {0, sp.Rational(1, 3), sp.Rational(1, 2)} and Pn[4].subs(t, sp.Rational(5, 12)) < 0)
check('P_3 < 0 on (0, 1/3): P_3(1/4) = -1/32', Pn[3].subs(t, sp.Rational(1, 4)) == sp.Rational(-1, 32))
half = [Pn[n].subs(t, sp.Rational(1, 2)) / fact(n) for n in range(1, 7)]
check('half rung: x + x^2/4 + x^3/48 + 0 x^4 + x^5/3840 - 7x^6/92160', half == [1, sp.Rational(1, 4), sp.Rational(1, 48), 0, sp.Rational(1, 3840), sp.Rational(-7, 92160)])
hs = sum(h * x**(i + 1) for i, h in enumerate([Pn[n].subs(t, sp.Rational(1, 2)) / fact(n) for n in range(1, Nf + 1)]))
comp = sp.series(hs.subs(x, hs), x, 0, Nf + 1).removeO()
check('half rung composed with itself = e^x - 1 to order 8', sp.expand(comp - sp.series(sp.exp(x) - 1, x, 0, Nf + 1).removeO()) == 0)
# group law of the flow, symbolic in s, t
fs = sum(Pn[n].subs(t, s) / fact(n) * x**n for n in range(1, Nf + 1)); ft = sum(Pn[n] / fact(n) * x**n for n in range(1, Nf + 1))
comp = sp.expand(sp.series(fs.subs(x, ft), x, 0, Nf + 1).removeO())
target = sp.expand(sum(Pn[n].subs(t, s + t) / fact(n) * x**n for n in range(1, Nf + 1)))
check('flow group law E^s(E^t(x)) = E^(s+t)(x) to order 8, symbolic s,t', sp.expand(comp - target) == 0)
gen = [sp.expand(sp.diff(Pn[n], t).subs(t, 0) / fact(n)) for n in range(2, 7)]
check('generator (iterative log): x^2/2 - x^3/12 + x^4/48 - x^5/180 + 11x^6/8640', gen == [sp.Rational(1, 2), sp.Rational(-1, 12), sp.Rational(1, 48), sp.Rational(-1, 180), sp.Rational(11, 8640)])

# fractional rungs turn negative: exact signs up to n = 16
Nb = 16
Eb = [F(0)] + [F(1, fact(i)) for i in range(1, Nb + 1)]
def mulb(a, b):
    c = [F(0)] * (Nb + 1)
    for i, ai in enumerate(a):
        if ai:
            for j, bj in enumerate(b):
                if i + j <= Nb: c[i + j] += ai * bj
    return c
Mb = sp.zeros(Nb, Nb); pw = [F(1)] + [F(0)] * Nb
for j in range(1, Nb + 1):
    pw = mulb(pw, Eb)
    for i in range(1, Nb + 1): Mb[i - 1, j - 1] = sp.Rational(pw[i].numerator, pw[i].denominator)
Nm = Mb - sp.eye(Nb); Lb = sp.zeros(Nb, Nb); P = sp.eye(Nb)
for m in range(1, Nb + 1): P = P * Nm; Lb += ((-1)**(m + 1)) * P / m
v = sp.zeros(Nb, 1); v[0] = 1; cur = v; res = v.copy()
for m in range(1, Nb + 1): cur = Lb * cur / m; res += cur * t**m
Pb = [None] + [sp.expand(res[i] * fact(i + 1)) for i in range(Nb)]
def negs(tv): return [n for n in range(1, Nb + 1) if Pb[n].subs(t, tv) < 0]
check('t = 1/2: first negative count at n = 6', negs(sp.Rational(1, 2))[0] == 6)
check('t = 11/10: negative at n = 12', 12 in negs(sp.Rational(11, 10)))
check('t = 5/4: negative at n = 13', 13 in negs(sp.Rational(5, 4)))
check('t = 3/2: positive to n = 15, negative at n = 16', negs(sp.Rational(3, 2)) == [16])
check('integer rungs t = 1..4 stay positive to n = 16', all(not negs(tv) for tv in (1, 2, 3, 4)))
import json
_g = json.load(open('figures.json'))['grid'] if __import__('os').path.exists('figures.json') else None
check('sign map (figures.json from figures.py): no negative cell in the t = 1 and t = 2 columns, n = 1..16', _g is not None and all(r[40] != '-' and r[80] != '-' for r in _g))
check('sign map: t = 0 column is the bare point x (1 then zeros)', _g is not None and _g[0][0] == '+' and all(r[0] == '0' for r in _g[1:]))
import json

# ---------- 7. The dial exponential counts blocks: Stirling numbers ----------
nn = sp.symbols('n', positive=True)
cc = sp.symbols('c')
lim = sp.limit(((1 + k * cc / nn)**nn - 1) / k, nn, sp.oo)
check('n copies of C/n under the k-dial -> (e^{kC} - 1)/k as n -> oo', sp.simplify(lim - (sp.exp(k * cc) - 1) / k) == 0)
lim1 = sp.limit((1 + cc / nn)**nn - 1, nn, sp.oo)
check('k = 1: n copies of C/n mult-added -> E(C) = e^C - 1 (the e staircase)', sp.simplify(lim1 - (sp.exp(cc) - 1)) == 0)
dialE = sp.series((sp.exp(k * (sp.exp(x) - 1)) - 1) / k, x, 0, 8).removeO()
ok = True
for n in range(1, 8):
    ok &= sp.expand(sp.expand(dialE).coeff(x, n) * fact(n) - sum(S2(n, j) * k**(j - 1) for j in range(1, n + 1))) == 0
check('n![x^n] (e^{k(e^x-1)} - 1)/k = sum_j S(n,j) k^(j-1)  (n<=7)', ok)
def blocks_bf(n):
    cnt = {}
    def parts(l):
        if not l: yield []; return
        a, rest = l[0], l[1:]
        for p in parts(rest):
            for i in range(len(p)): yield p[:i] + [[a] + p[i]] + p[i + 1:]
            yield [[a]] + p
    for p in parts(list(range(n))): cnt[len(p)] = cnt.get(len(p), 0) + 1
    return cnt
check('S(n,j) = set partitions with j blocks (brute force n<=7)', all(blocks_bf(n).get(j, 0) == S2(n, j) for n in range(1, 8) for j in range(1, n + 1)))
check('k = -1: sum_j (-1)^(j-1) S(n,j) = 1,0,-1,-1,2,9,9,-50 (complementary Bell, negated)',
      [sum((-1)**(j - 1) * S2(n, j) for j in range(1, n + 1)) for n in range(1, 9)] == [1, 0, -1, -1, 2, 9, 9, -50])

# ---------- 8. Compositions: the cube over parts ----------
ogf = x / (1 - x)
ok = True
for m in range(1, 5):
    f = sp.series((1 + ogf)**m - 1, x, 0, 10).removeO()
    acc = 0
    for _ in range(m): acc = acc + ogf + acc * ogf
    ok &= sp.expand(sp.series(acc, x, 0, 10).removeO() - f) == 0
    for n in range(1, 10):
        bf = sum(1 for w in itertools.product(range(n + 1), repeat=m) if sum(w) == n)
        layers = sum(C(m, j) * C(n - 1, j - 1) for j in range(1, m + 1))
        ok &= (sp.expand(f).coeff(x, n) == bf == layers == C(n + m - 1, m - 1))
check('m parts, each "a positive part" mult-added: weak compositions C(n+m-1,m-1) = sum_j C(m,j)C(n-1,j-1) (m<=4,n<=9)', ok)
check('two parts: pairs (a,b), either or both nonzero, a+b=n: n+1', all(sum(1 for a in range(n + 1) for b in range(n + 1) if a + b == n) == n + 1 for n in range(1, 12)))

print(f'\n{sum(R)} of {len(R)} checks pass')
