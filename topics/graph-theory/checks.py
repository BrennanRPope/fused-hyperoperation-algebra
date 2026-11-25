"""Graph theory through mult-add: every claim on index.html, checked.
Run: python3 checks.py > checks_output.txt"""
import itertools, math, random
import numpy as np
import sympy as sp
import networkx as nx

R = []
def check(name, ok):
    R.append(bool(ok)); print(('PASS' if ok else 'FAIL') + '  ' + name)

Q, v1, v2, v3, t1, t2, p, a, b, J1, J2, x = sp.symbols('Q v1 v2 v3 t1 t2 p a b J1 J2 x')
lens = lambda P, q, a, b: (a + b + P * a * b) / (1 - q * a * b)
ma = lambda a, b: a + b + a * b
fu = lambda a, b: a + b - a * b          # k = -1, failure union / OR

# ---------- multivariate Tutte (Fortuin-Kasteleyn) by brute force on multigraphs ----------
def comps(n, edges):
    par = list(range(n))
    def f(i):
        while par[i] != i: par[i] = par[par[i]]; i = par[i]
        return i
    for (u, w) in edges: par[f(u)] = f(w)
    return len({f(i) for i in range(n)})
def Z(n, E, ws, Qs=Q):
    tot = 0
    for mask in range(2**len(E)):
        A = [E[i] for i in range(len(E)) if mask >> i & 1]
        tot += Qs**comps(n, A) * sp.prod([ws[i] for i in range(len(E)) if mask >> i & 1])
    return sp.expand(tot)

# 1. parallel reduction: two parallel edges = one edge of weight v1 * v2 (mult-add)
E2 = [(0, 1), (0, 1), (1, 2), (2, 0)]
lhs = Z(3, E2, [v1, v2, v3, v3]); rhs = Z(3, [(0, 1), (1, 2), (2, 0)], [ma(v1, v2), v3, v3])
check('parallel edges: Z(G) = Z(G with one edge of weight v1 * v2 = v1+v2+v1v2), triangle', sp.expand(lhs - rhs) == 0)
# 2. series reduction
E4 = [(0, 1), (1, 2), (2, 3), (3, 0)]
vs = v1 * v2 / (Q + v1 + v2)
lhs = Z(4, E4, [v1, v2, v3, v3]); rhs = (Q + v1 + v2) * Z(3, [(0, 2), (2, 1), (1, 0)], [vs, v3, v3])
check('series edges: Z(G) = (Q+v1+v2) * Z(G\') with v = v1 v2/(Q+v1+v2), 4-cycle', sp.simplify(lhs - rhs) == 0)
# 3. deletion-contraction
K4 = list(itertools.combinations(range(4), 2)); ws = sp.symbols('w0:6')
lhs = Z(4, K4, ws)
dele = Z(4, K4[1:], ws[1:])
contr_edges = [tuple(0 if u == 1 else (u if u < 1 else u - 1) for u in e) for e in K4[1:]]
contr = Z(3, contr_edges, ws[1:])
check('deletion-contraction Z(G) = Z(G-e) + v_e Z(G/e), K4 with six weights', sp.expand(lhs - dele - ws[0] * contr) == 0)
# 4. Tutte polynomial relation and chromatic polynomial
X, Y = sp.symbols('X Y')
for G in [nx.complete_graph(4), nx.cycle_graph(5), nx.petersen_graph()]:
    if G.number_of_edges() > 12: continue
    n = G.number_of_nodes(); E = list(G.edges())
    Zs = Z(n, E, [Y - 1] * len(E), Qs=(X - 1) * (Y - 1))
    T = sp.expand(nx.tutte_polynomial(G).subs({sp.Symbol('x'): X, sp.Symbol('y'): Y}))
    check(f'Z = (X-1)^k (Y-1)^n T(X,Y) with Q=(X-1)(Y-1), v=Y-1: {n} vertices, {len(E)} edges', sp.expand(Zs - (X - 1) * (Y - 1)**n * T) == 0)
check('chromatic polynomial = Z at v = -1: K4 gives Q(Q-1)(Q-2)(Q-3)', sp.expand(Z(4, K4, [-1] * 6) - Q*(Q-1)*(Q-2)*(Q-3)) == 0)

# ---------- 5. transmissivity t = v/(Q+v): the parallel law is a landscape law on the "or" line ----------
tt = lambda v: v / (Q + v)
par = sp.simplify(tt(ma(v1, v2)) - lens(Q - 2, 1 - Q, tt(v1), tt(v2)))
check('parallel in t: t(v1 * v2) = (t1+t2+(Q-2)t1t2)/(1+(Q-1)t1t2) = F_{Q-2, 1-Q}', par == 0)
check('series in t: t(v1 v2/(Q+v1+v2)) = t1 t2', sp.simplify(tt(vs) - tt(v1) * tt(v2)) == 0)
check('or line: (Q-2) + (1-Q) = -1 for every Q', sp.expand((Q - 2) + (1 - Q) + 1) == 0)
check('Q = 1 (percolation, reliability): F_{-1,0} = a + b - ab', sp.simplify(lens(-1, 0, a, b) - fu(a, b)) == 0)
check('Q = 2 (Ising): F_{0,-1} = (a+b)/(1+ab), velocity addition', sp.simplify(lens(0, -1, a, b) - (a + b) / (1 + a * b)) == 0)
K1, K2 = 0.37, 0.81
check('Ising: v = e^{2K} - 1, t = v/(2+v) = tanh K, and parallel bonds add K (numeric)',
      abs((math.exp(2*K1)-1)/(1+math.exp(2*K1)) - math.tanh(K1)) < 1e-14 and abs((math.tanh(K1)+math.tanh(K2))/(1+math.tanh(K1)*math.tanh(K2)) - math.tanh(K1+K2)) < 1e-14)
check('edge weight is the flow E of the coupling: E(J1+J2) = E(J1) * E(J2)', sp.simplify(ma(sp.exp(J1) - 1, sp.exp(J2) - 1) - (sp.exp(J1 + J2) - 1)) == 0)
check('lift of the parallel law: 1+v = (1+(Q-1)t)/(1-t), so it multiplies', sp.simplify((1 + (Q - 1) * tt(v1)) / (1 - tt(v1)) - (1 + v1)) == 0)
check('wall 1: a perfect bond t = 1 absorbs: F(a, 1) = 1', sp.simplify(lens(Q - 2, 1 - Q, a, 1) - 1) == 0)
w = -1 / (Q - 1)
check('wall 2: t = -1/(Q-1) absorbs, and it is v = -1 (proper colouring)', sp.simplify(lens(Q - 2, 1 - Q, a, w) - w) == 0 and sp.simplify(tt(-1) - w) == 0)
th = sp.symbols('theta')
amh = lambda th, u, vv: u * vv / (1 - th * (1 - u) * (1 - vv))
check('Potts parallel = union under Ali-Mikhail-Haq dependence, theta = 1 - Q', sp.simplify(1 - amh(1 - Q, 1 - a, 1 - b) - lens(Q - 2, 1 - Q, a, b)) == 0)
check('Q = 0 (spanning forests): F_{-2,1} = (a+b-2ab)/(1-ab), the double-wall cusp k^2 = 4q', sp.simplify(lens(-2, 1, a, b) - (a + b - 2*a*b) / (1 - a*b)) == 0 and (-2)**2 == 4 * 1)
kk_ = sp.symbols('kk')
check('the or line q = -1 - k touches the border parabola q = k^2/4 only at k = -2, q = 1 (Q = 0)', sp.solve(sp.Eq(-1 - kk_, kk_**2 / 4), kk_) == [-2])
check('the parallel law is associative for every Q (symbolic)', sp.simplify(lens(Q-2, 1-Q, lens(Q-2, 1-Q, a, b), x) - lens(Q-2, 1-Q, a, lens(Q-2, 1-Q, b, x))) == 0)

# ---------- 6. two-terminal reliability: series-parallel = (x, *_{-1}) ----------
def rel_bf(n, E, s, tgt, ps):
    tot = 0
    for st in itertools.product([0, 1], repeat=len(E)):
        g = nx.Graph(); g.add_nodes_from(range(n)); g.add_edges_from(E[i] for i in range(len(E)) if st[i])
        pr = math.prod(ps[i] if st[i] else 1 - ps[i] for i in range(len(E)))
        if nx.has_path(g, s, tgt): tot += pr
    return tot
random.seed(7)
def rand_sp(depth):
    # returns (edge list on nodes, s, t, reduction value function) for a random series-parallel network
    if depth == 0:
        return {'E': [(0, 1)], 'n': 2, 'val': lambda ps, i=[0]: None}
    return None
# build SP networks explicitly with expression trees
def build(tree, nodes, s, tg, E, P, out):
    kind = tree[0]
    if kind == 'e':
        E.append((s, tg)); P.append(tree[1]); return tree[1]
    if kind == 's':
        m = nodes[0]; nodes[0] += 1
        return build(tree[1], nodes, s, m, E, P, out) * build(tree[2], nodes, m, tg, E, P, out)
    if kind == 'p':
        return fu(build(tree[1], nodes, s, tg, E, P, out), build(tree[2], nodes, s, tg, E, P, out))
def rtree(d):
    if d == 0 or random.random() < 0.25: return ('e', random.uniform(0.2, 0.95))
    return (random.choice('sp'), rtree(d - 1), rtree(d - 1))
ok = True; tested = 0
while tested < 40:
    tr = rtree(4); E = []; P = []; nodes = [2]
    val = build(tr, nodes, 0, 1, E, P, None)
    if len(E) > 12 or len(E) < 3: continue
    ok &= abs(rel_bf(nodes[0], E, 0, 1, P) - val) < 1e-12; tested += 1
check('40 random series-parallel networks: brute-force reliability = (x, *_{-1}) reduction', ok)
pp = 0.73
check('failure side: series fails = (1-p1) *_{-1} (1-p2), parallel fails = (1-p1)(1-p2)', sp.expand(1 - a*b - fu(1 - a, 1 - b)) == 0 and sp.expand(1 - fu(a, b) - (1 - a)*(1 - b)) == 0)
u = lambda q: 1 / q - 1
check('odds against u = 1/p - 1: series is plain mult-add, u(p1 p2) = u1 * u2', sp.simplify(u(a * b) - ma(u(a), u(b))) == 0)
bridgeE = [(0, 2), (0, 3), (2, 3), (2, 1), (3, 1)]
Rb = sum((p if st[i] else 1 - p) for i in range(0)) if False else None
Rbr = 0
for st in itertools.product([0, 1], repeat=5):
    g = nx.Graph(); g.add_nodes_from(range(4)); g.add_edges_from(bridgeE[i] for i in range(5) if st[i])
    if nx.has_path(g, 0, 1): Rbr += sp.prod([p if s_ else 1 - p for s_ in st])
Rbr = sp.expand(Rbr)
check('Wheatstone bridge reliability = 2p^2 + 2p^3 - 5p^4 + 2p^5 (brute force)', sp.expand(Rbr - (2*p**2 + 2*p**3 - 5*p**4 + 2*p**5)) == 0)
piv = p * fu(p, p)**2 + (1 - p) * fu(p*p, p*p)
check('bridge by pivoting on the middle edge: p*(p*p)^2 [contracted] + (1-p)*(p^2 *_{-1} p^2) [deleted]', sp.expand(piv - Rbr) == 0)
def has_sp_reduction(n, E, s, tg):
    # a series pair (degree-2 non-terminal vertex) or a parallel pair would allow one reduction step
    g = nx.MultiGraph(); g.add_edges_from(E)
    par_pair = any(g.number_of_edges(u_, w_) > 1 for u_, w_ in g.edges())
    ser = any(g.degree(vv) == 2 for vv in g.nodes() if vv not in (s, tg))
    return par_pair or ser
check('the bridge has no series or parallel pair: the dial pair alone cannot reduce it', not has_sp_reduction(4, bridgeE, 0, 1))

# ---------- 7. the diamond lattice: the dial pair iterated ----------
dia = lambda q: fu(q*q, q*q)
check('diamond step p -> (p p) *_{-1} (p p) = 2p^2 - p^4', sp.expand(dia(p) - (2*p**2 - p**4)) == 0)
# level-2 diamond: 16 edges, brute force reliability equals dia(dia(p)) at a few p
def diamond_edges(level):
    E = [(0, 1)]; n = 2
    for _ in range(level):
        newE = []
        for (s_, t_) in E:
            m1, m2 = n, n + 1; n += 2
            newE += [(s_, m1), (m1, t_), (s_, m2), (m2, t_)]
        E = newE
    return n, E
n2, E2d = diamond_edges(2)
ok = all(abs(rel_bf(n2, E2d, 0, 1, [pv] * 16) - float(dia(dia(sp.Float(pv))))) < 1e-12 for pv in (0.3, 0.618, 0.8))
check('level-2 diamond (16 edges), brute force = diamond step applied twice', ok)
g = (sp.sqrt(5) - 1) / 2
check('percolation fixed point (sqrt5 - 1)/2, the golden ratio', sp.simplify(dia(g) - g) == 0)
slope = sp.diff(dia(p), p).subs(p, g)
nu = math.log(2) / math.log(float(slope))
check('slope there 4 - 4p* = %.4f, length exponent ln2/ln(slope) = %.4f' % (float(slope), nu), sp.simplify(slope - (4 - 4*g)) == 0 and abs(nu - 1.6353) < 1e-3)
diaQ = lambda q, Qv: lens(Qv - 2, 1 - Qv, q*q, q*q)
ti = sp.nsolve(diaQ(x, 2) - x, x, 0.55)
check('Ising (Q=2) diamond fixed point solves t^3 + t^2 + t = 1: t* = %.6f' % ti, abs(ti**3 + ti**2 + ti - 1) < 1e-12 and abs(ti - 0.543689) < 1e-6)
# Potts diamond: brute force effective t from Z with terminals (Q symbolic), level 1
def eff_t(n, E, Qs, vv):
    zc = 0; zs = 0
    for mask in range(2**len(E)):
        A = [E[i] for i in range(len(E)) if mask >> i & 1]
        par = list(range(n))
        def f(i):
            while par[i] != i: i = par[i]
            return i
        for (uu, ww) in A: par[f(uu)] = f(ww)
        k_ = len({f(i) for i in range(n)}); wgt = Qs**k_ * vv**len(A)
        if f(0) == f(1): zc += wgt
        else: zs += wgt
    veff = Qs * zc / zs
    return sp.simplify(veff / (Qs + veff))
n1, E1d = diamond_edges(1)
vv = sp.symbols('v')
check('Potts diamond (level 1, any Q): effective t from brute-force Z = F_{Q-2,1-Q}(t^2, t^2)', sp.simplify(eff_t(n1, E1d, Q, vv) - diaQ(tt(vv), Q)) == 0)
def eff_t2(n, E, Qs, vv, s0, s1):
    zc = 0; zs = 0
    for mask in range(2**len(E)):
        A = [E[i] for i in range(len(E)) if mask >> i & 1]
        par = list(range(n))
        def f(i):
            while par[i] != i: i = par[i]
            return i
        for (uu, ww) in A: par[f(uu)] = f(ww)
        wgt = Qs**len({f(i) for i in range(n)}) * vv**len(A)
        if f(s0) == f(s1): zc += wgt
        else: zs += wgt
    veff = Qs * zc / zs
    return sp.simplify(veff / (Qs + veff))
heroE = [(0, 2), (2, 1), (0, 1), (1, 3)]   # top branch 0-2-1 in series, bottom bond 0-1, last bond 1-3
tv_ = tt(vv)
check('hero network, any Q: brute-force Potts effective t (s=0 to end=3) = F(t^2, t) * t', sp.simplify(eff_t2(4, heroE, Q, vv, 0, 3) - lens(Q - 2, 1 - Q, tv_**2, tv_) * tv_) == 0)
check('Q = 0: the diamond step is 2t^2/(1+t^2); its only fixed points are 0 and a double one at 1', sp.simplify(diaQ(x, 0) - 2*x**2/(1 + x**2)) == 0 and sp.simplify(sp.factor(sp.together(sp.simplify(diaQ(x, 0)) - x)) + x*(x - 1)**2/(x**2 + 1)) == 0)
fixed = {}
def tstar(Qv):
    num = sp.numer(sp.together(diaQ(x, sp.nsimplify(Qv)) - x))
    rts = [complex(r) for r in sp.Poly(sp.expand(num), x).nroots()]
    rr = [r.real for r in rts if abs(r.imag) < 1e-12 and 1e-9 < r.real < 1 - 1e-9]
    return rr
for Qv in (0.5, 1, 2, 3, 4):
    rr = tstar(Qv); fixed[Qv] = rr[0] if len(rr) == 1 else None
check('exactly one fixed point t*(Q) strictly inside (0,1) for Q = 0.5,1,2,3,4: ' + ', '.join('%.4f' % fixed[q] for q in fixed), all(fixed[q] is not None for q in fixed))
check('Q = 3: the diamond fixed point is exactly 1/2', sp.simplify(diaQ(sp.Rational(1, 2), 3) - sp.Rational(1, 2)) == 0)
check('t*(Q) decreases as Q grows (0.5 to 4)', all(fixed[a_] > fixed[b_] for a_, b_ in zip(list(fixed)[:-1], list(fixed)[1:])))

# ---------- 8. graph products: the strong product is the Kronecker mult-add ----------
def conormal(G, H):
    P = nx.Graph(); P.add_nodes_from(itertools.product(G.nodes(), H.nodes()))
    for (g1, h1), (g2, h2) in itertools.combinations(P.nodes(), 2):
        if G.has_edge(g1, g2) or H.has_edge(h1, h2): P.add_edge((g1, h1), (g2, h2))
    return P
gs = [nx.path_graph(3), nx.cycle_graph(4), nx.complete_graph(3), nx.star_graph(3), nx.petersen_graph(), nx.path_graph(2)]
ok1 = ok2 = ok3 = ok4 = ok5 = ok6 = True
for G, H in itertools.product(gs, repeat=2):
    n1_, m1 = G.number_of_nodes(), G.number_of_edges(); n2_, m2 = H.number_of_nodes(), H.number_of_edges()
    S = nx.strong_product(G, H); Cp = nx.cartesian_product(G, H); T = nx.tensor_product(G, H)
    ok1 &= S.number_of_edges() == n1_*m2 + n2_*m1 + 2*m1*m2
    ok2 &= (2*S.number_of_edges() + n1_*n2_) == (2*m1 + n1_) * (2*m2 + n2_)
    ok3 &= all(S.degree((g_, h_)) == ma(G.degree(g_), H.degree(h_)) for g_ in G for h_ in H)
    A = nx.to_numpy_array(G, nodelist=sorted(G.nodes())); B = nx.to_numpy_array(H, nodelist=sorted(H.nodes()))
    I1, I2 = np.eye(n1_), np.eye(n2_)
    nl = [(g_, h_) for g_ in sorted(G.nodes()) for h_ in sorted(H.nodes())]
    AS = nx.to_numpy_array(S, nodelist=nl); AC = nx.to_numpy_array(Cp, nodelist=nl); AT = nx.to_numpy_array(T, nodelist=nl)
    ok4 &= np.allclose(AS, np.kron(A, I2) + np.kron(I1, B) + np.kron(A, B)) and np.allclose(AC, np.kron(A, I2) + np.kron(I1, B)) and np.allclose(AT, np.kron(A, B))
    la, lb = np.linalg.eigvalsh(A), np.linalg.eigvalsh(B)
    ok5 &= np.allclose(np.sort(np.linalg.eigvalsh(AS)), np.sort([ma(x_, y_) for x_ in la for y_ in lb]))
    if n1_ * n2_ <= 40:
        O = conormal(G, H); AO = nx.to_numpy_array(O, nodelist=nl)
        J1_, J2_ = np.ones((n1_, n1_)), np.ones((n2_, n2_))
        ok6 &= np.allclose(np.ones_like(AO) - AO, np.kron(J1_ - A, J2_ - B) - (np.eye(n1_*n2_) - np.eye(n1_*n2_)) ) or np.allclose((np.ones_like(AO) - AO), np.kron(J1_ - A, J2_ - B))
check('P3 strong P3 has 20 edges: 12 Cartesian (k = 0 part) + 8 diagonal (the corner)', nx.strong_product(nx.path_graph(3), nx.path_graph(3)).number_of_edges() == 20 and nx.cartesian_product(nx.path_graph(3), nx.path_graph(3)).number_of_edges() == 12 and nx.tensor_product(nx.path_graph(3), nx.path_graph(3)).number_of_edges() == 8)
check('strong product edges = n1 m2 + n2 m1 + 2 m1 m2 (36 pairs of small graphs)', ok1)
check('lifted count 2m + n multiplies: (2m+n)(G strong H) = (2m1+n1)(2m2+n2)', ok2)
check('closed degrees multiply: deg = d1 * d2 = d1 + d2 + d1 d2', ok3)
check('adjacency: strong = A(x)I + I(x)B + A(x)B (Kronecker mult-add); Cartesian = k = 0 part; tensor = the corner', ok4)
check('eigenvalues of the strong product are lambda * mu (lambda+mu+lambda mu)', ok5)
check('OR (co-normal) product: J - A_or = (J - A)(x)(J - B), the k = -1 dial with J in place of I', ok6)

# ---------- 9. independence polynomial: union mult-adds, join adds ----------
def indep(G):
    nl = list(G.nodes()); c = [0] * (len(nl) + 1)
    for r in range(len(nl) + 1):
        for S in itertools.combinations(nl, r):
            if not any(G.has_edge(u_, w_) for u_, w_ in itertools.combinations(S, 2)): c[r] += 1
    return sp.expand(sum(ci * x**i for i, ci in enumerate(c)))
ok_u = ok_j = True
for G, H in itertools.product(gs[:4] + [nx.path_graph(4)], repeat=2):
    IG, IH = indep(G), indep(H)
    U = nx.disjoint_union(G, H); Jn = nx.complement(nx.disjoint_union(nx.complement(G), nx.complement(H)))
    ok_u &= sp.expand(indep(U) - 1 - ma(IG - 1, IH - 1)) == 0
    ok_j &= sp.expand(indep(Jn) - 1 - ((IG - 1) + (IH - 1))) == 0
check('disjoint union: I(G u H) - 1 = (I(G)-1) * (I(H)-1) (mult-add)', ok_u)
check('join: I(G + H) - 1 = (I(G)-1) + (I(H)-1) (plain addition, k = 0)', ok_j)
check('C4 = join of two non-edges: I = 1 + 4x + 2x^2', sp.expand(indep(nx.cycle_graph(4)) - (1 + 4*x + 2*x**2)) == 0)
kz = sp.symbols('k')
def dialgas(G, z, kk):
    nl = list(G.nodes()); tot = 0
    for r in range(len(nl) + 1):
        for S in itertools.combinations(nl, r):
            e = sum(1 for u_, w_ in itertools.combinations(S, 2) if G.has_edge(u_, w_))
            tot += z**r * (kk**e if e else 1)
    return sp.expand(tot)
check('the Lee-Yang dial gas at k = 0 is the independence polynomial (Petersen, C5)',
      all(sp.expand(dialgas(G, x, 0) - indep(G)) == 0 for G in (nx.petersen_graph(), nx.cycle_graph(5))))
check('the dial gas at k = 1 is free: (1+z)^n (C5)', sp.expand(dialgas(nx.cycle_graph(5), x, 1) - (1 + x)**5) == 0)

# ---------- 10. Ihara zeta: Bass's formula, Euler product; no extra structure ----------
uu = sp.symbols('u')
def ihara_checks(G):
    n_ = G.number_of_nodes(); m_ = G.number_of_edges()
    A = sp.Matrix(nx.to_numpy_array(G, nodelist=sorted(G.nodes())).astype(int)); D = sp.diag(*[G.degree(v_) for v_ in sorted(G.nodes())])
    bass = sp.expand((1 - uu**2)**(m_ - n_) * (sp.eye(n_) - uu * A + uu**2 * (D - sp.eye(n_))).det())
    darts = [(a_, b_) for a_, b_ in G.edges()] + [(b_, a_) for a_, b_ in G.edges()]
    idx = {d_: i for i, d_ in enumerate(darts)}
    Bm = sp.zeros(len(darts), len(darts))
    for (a_, b_) in darts:
        for (c_, d_) in darts:
            if b_ == c_ and d_ != a_: Bm[idx[(a_, b_)], idx[(c_, d_)]] = 1
    hash_ = sp.expand((sp.eye(len(darts)) - uu * Bm).det())
    return bass, hash_, Bm
ok = True
for G in (nx.complete_graph(4), nx.complete_bipartite_graph(2, 3), nx.Graph([(0, 1), (1, 2), (2, 0), (2, 3), (3, 4), (4, 2)])):
    bass, hsh, Bm = ihara_checks(G)
    ok &= sp.expand(bass - hsh) == 0
check("Bass: (1-u^2)^(m-n) det(I - uA + u^2(D-I)) = det(I - uB) (Hashimoto), 3 graphs incl. non-regular", ok)
bass, hsh, Bm = ihara_checks(nx.complete_graph(4))
logz = sp.series(-sp.log(hsh), uu, 0, 9).removeO()
Nm = [sp.trace(Bm**mm) for mm in range(1, 9)]
check('Euler product: log zeta = sum N_m u^m / m with N_m = tr B^m (closed non-backtracking walks), K4', all(sp.expand(logz).coeff(uu, mm) * mm == Nm[mm - 1] for mm in range(1, 9)))
check('K4: shortest prime cycles are the triangles: N_1 = N_2 = 0, N_3 = 24 (8 oriented triangles x 3 starts)', Nm[0] == 0 and Nm[1] == 0 and Nm[2] == 24)
q_ = 2
lam = np.linalg.eigvalsh(nx.to_numpy_array(nx.petersen_graph()))
roots_ok = True
for l_ in lam:
    if abs(abs(l_) - 3) < 1e-9: continue
    rts = np.roots([q_, -l_, 1])
    roots_ok &= np.allclose(np.abs(rts), 1 / math.sqrt(q_))
check('Petersen (3-regular, q = 2) is Ramanujan: every nontrivial factor 1 - lambda u + 2u^2 has roots on |u| = 1/sqrt2', roots_ok and max(abs(l_) for l_ in lam if abs(abs(l_) - 3) > 1e-9) <= 2 * math.sqrt(2))
check('cycle C_n: 1/zeta = (1 - u^n)^2 (one prime cycle each way), n = 5', sp.expand(ihara_checks(nx.cycle_graph(5))[1] - (1 - uu**5)**2) == 0)

print(f'\n{sum(R)} of {len(R)} checks pass')
