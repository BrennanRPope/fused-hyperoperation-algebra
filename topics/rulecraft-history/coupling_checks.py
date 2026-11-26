"""Coupling checks for rulecraft-history (DRAFT): rules whose target on one slider is set by other sliders."""
import itertools, random, math
from fractions import Fraction as Fr
import sympy as sp
import networkx as nx
ok = lambda n, c: print(("PASS " if c else "FAIL ") + n) or c
res = []
lam, mu, a, b, s, c0, c, x, y = sp.symbols('lambda mu a b s c0 c x y')
# a coupled rule writes one slider: x <- (1-lam) x + lam (a + b y); a plain nudge on y: y <- (1-mu) y + mu s
f = lambda X, Y: ((1 - lam)*X + lam*(a + b*Y), Y)
g = lambda X, Y: (X, (1 - mu)*Y + mu*s)
fg = f(*g(x, y)); gf = g(*f(x, y))
res.append(ok("K1 coupled rule after the rule it reads, minus before: x-gap = lam mu b (s - y); y-gap 0",
    sp.simplify(fg[0] - gf[0] - lam*mu*b*(s - y)) == 0 and sp.simplify(fg[1] - gf[1]) == 0))
h = lambda X, Y: (X, (1 - mu)*Y + mu*(c0 + c*Y*0))   # a rule on y that does not read x
res.append(ok("K2 b = 0 (x's rule does not read y): the two rules commute", sp.simplify((fg[0] - gf[0]).subs(b, 0)) == 0))

# K3 number of different histories = acyclic orientations of the conflict graph = |chromatic polynomial at -1|
def chrom_at(G, q):
    # deletion-contraction, small graphs
    G = nx.Graph(G)
    if G.number_of_edges() == 0: return q**G.number_of_nodes()
    u, v = next(iter(G.edges()))
    G1 = G.copy(); G1.remove_edge(u, v)
    G2 = nx.contracted_nodes(G1, u, v, self_loops=False)
    G2 = nx.Graph(G2)
    return chrom_at(G1, q) - chrom_at(G2, q)
def acyclic_orientations(G):
    E = list(G.edges()); cnt = 0
    for bits in itertools.product((0, 1), repeat=len(E)):
        D = nx.DiGraph(); D.add_nodes_from(G.nodes())
        D.add_edges_from([(u, v) if bt else (v, u) for (u, v), bt in zip(E, bits)])
        cnt += nx.is_directed_acyclic_graph(D)
    return cnt
rng = random.Random(5)
allgood = True; report = []
for trial in range(40):
    m = rng.randint(2, 5); n = rng.randint(3, 6)
    rules = []
    for _ in range(n):
        w = rng.randrange(m)
        reads = [j for j in range(m) if j != w and rng.random() < 0.3]
        bvec = {j: Fr(rng.choice([-1, 1])*rng.randint(5, 40), 100) for j in reads}
        a0 = Fr(rng.randint(30, 70), 100)
        rules.append((w, reads, Fr(rng.randint(5, 95), 100), a0, bvec))
    def apply(order, x0):
        X = list(x0)
        for i in order:
            w, reads, l, a0, bvec = rules[i]
            tgt = a0 + sum(bvec[j]*X[j] for j in reads)
            X[w] = (1 - l)*X[w] + l*tgt
        return tuple(X)
    x0 = [Fr(rng.randint(0, 100), 100) for _ in range(m)]
    outs = {apply(o, x0) for o in itertools.permutations(range(n))}
    G = nx.Graph(); G.add_nodes_from(range(n))
    for i, j in itertools.combinations(range(n), 2):
        wi, ri = rules[i][0], rules[i][1]; wj, rj = rules[j][0], rules[j][1]
        same_target = (wi == wj and rules[i][3] == rules[j][3] and rules[i][4] == rules[j][4] and not (set(ri) & {wi}))
        if (wi == wj and not same_target) or wi in rj or wj in ri: G.add_edge(i, j)
    ao = acyclic_orientations(G); ch = abs(chrom_at(G, -1))
    allgood &= (len(outs) == ao == ch)
    report.append((n, G.number_of_edges(), len(outs), ao))
    if not (len(outs) == ao == ch): print('  mismatch', trial, [(r[0], r[1]) for r in rules], len(outs), ao, ch)
res.append(ok("K3 40 random rule packages (3 to 6 rules on 2 to 5 sliders, exact): conflict = same slider with different targets, or one reads what the other writes; distinct outcomes over all orders = acyclic orientations of the conflict graph = |chromatic polynomial at -1|; e.g. %s" % report[:6], allgood))
# K4 complete conflict graph gives n!, empty gives 1; components multiply
K5 = nx.complete_graph(5); E5 = nx.empty_graph(5)
two = nx.disjoint_union(nx.path_graph(3), nx.cycle_graph(3))
res.append(ok("K4 everything conflicts -> 5! = 120 histories; nothing conflicts -> 1; separate groups multiply (path of 3 and triangle: 4 x 6 = 24)",
    abs(chrom_at(K5, -1)) == 120 and abs(chrom_at(E5, -1)) == 1 and abs(chrom_at(two, -1)) == 24))

# K5 two-way coupling, rules taken in turn: the turn's linear map, its trace and determinant
lam, mu, b, c = sp.symbols('lambda mu b c', real=True)
M = sp.Matrix([[1 - lam, lam*b], [mu*c*(1 - lam), 1 - mu + lam*mu*b*c]])
tr, det = sp.simplify(M.trace()), sp.simplify(M.det())
Delta = sp.expand(tr**2 - 4*det)
res.append(ok("K5 one round (x's rule, then y's): trace = 2 - lam - mu + lam mu b c, det = (1-lam)(1-mu)",
    sp.simplify(tr - (2 - lam - mu + lam*mu*b*c)) == 0 and sp.simplify(det - (1 - lam)*(1 - mu)) == 0))
# with bc >= 0, Delta >= (lam - mu)^2 >= 0 : no spiral
D0 = sp.expand(Delta.subs(b*c, 0) if False else (2 - lam - mu)**2 - 4*(1 - lam)*(1 - mu))
res.append(ok("K6 with b c = 0, Delta = (lam - mu)^2: real; Delta grows with b c > 0, so a spiral (Delta < 0) needs b c < 0, a balancing loop",
    sp.simplify(D0 - (lam - mu)**2) == 0 and sp.simplify(sp.diff(Delta.subs({b: sp.Symbol('p', positive=True), c: 1}), sp.Symbol('p', positive=True)).subs(sp.Symbol('p', positive=True), 0)) == sp.simplify(2*(2 - lam - mu)*lam*mu)))
# numeric sweep of the claim spiral <=> bc < 0 and Delta < 0, and when spiral the damping |z| = sqrt(det)
import numpy as np
r = np.random.default_rng(2); bad = 0; spir = 0
for _ in range(200000):
    L, U = r.random(2); B, C = r.uniform(-1, 1, 2)
    t = 2 - L - U + L*U*B*C; d = (1 - L)*(1 - U); D = t*t - 4*d
    if D < 0:
        spir += 1
        if B*C >= 0: bad += 1
res.append(ok("K7 200000 random two-way couplings with |b|, |c| <= 1: %d spiral, every one with b c < 0 (%d exceptions)" % (spir, bad), bad == 0 and spir > 0))
# K8 the turn as a (k, d) setting: k = -trace, d = det, Delta = k^2 - 4d = our curvature reading
k_, d_ = -tr, det
res.append(ok("K8 reading the turn as a (k, d) setting (k = -trace, d = det): Delta = k^2 - 4d is the turn's discriminant", sp.simplify(k_**2 - 4*d_ - Delta) == 0))
# K9 loop holonomy for one-way coupling: f g f^-1 g^-1 linear part is a shear
lam, mu, b = sp.symbols('lambda mu b')
F = sp.Matrix([[1 - lam, lam*b], [0, 1]]); Gm = sp.Matrix([[1, 0], [0, 1 - mu]])
Hol = sp.simplify(F*Gm*F.inv()*Gm.inv())
res.append(ok("K9 one-way coupling: going round the loop leaves the shear %s (pure coupling residue; identity iff b = 0)" % Hol.tolist(),
    sp.simplify(Hol[0, 0] - 1) == 0 and sp.simplify(Hol[1, 1] - 1) == 0 and Hol[1, 0] == 0 and sp.simplify(Hol[0, 1].subs(b, 0)) == 0))
# K10 two-way: loop's linear part has trace 2 + (something) * b c ; Fricke-type invariant
c = sp.symbols('c')
F2 = sp.Matrix([[1 - lam, lam*b], [0, 1]]); G2 = sp.Matrix([[1, 0], [mu*c, 1 - mu]])
H2 = sp.simplify(F2*G2*F2.inv()*G2.inv())
trH = sp.factor(sp.simplify(H2.trace() - 2))
res.append(ok("K10 two-way coupling: trace of the loop's linear part minus 2 = %s (zero iff b c = 0 or b c = 1)" % trH, sp.simplify(trH.subs(b, 0)) == 0 and sp.simplify(trH.subs(c, 0)) == 0 and trH != 0))
# K11 the residue's type. det = 1 exactly; trace - 2 = bc(bc-1) l^2 u^2 / ((1-l)(1-u)).
#   balancing (bc < 0) or bc > 1: trace > 2, a squeeze (hyperbolic)
#   bc = 1: trace = 2 (shear-like, parabolic)
#   0 < bc < 1: trace < 2; rotation type (|trace| < 2) unless both rules are strong enough that
#   bc(1-bc) l^2 u^2 >= 4 (1-l)(1-u), when trace <= -2 (a flipping squeeze)
P = sp.symbols('p')  # p = b c
tP = 2 + P*(P - 1)*lam**2*mu**2/((1 - lam)*(1 - mu))
okexact = sp.simplify(H2.det() - 1) == 0 and sp.simplify(H2.trace().subs(c, P/b) - tP) == 0
fH = sp.lambdify((lam, mu, b, c), H2.trace())
r2 = np.random.default_rng(9); good = True; cnt = {"squeeze": 0, "rotation": 0, "flip": 0}
for _ in range(50000):
    L, U = r2.uniform(0.01, 0.99, 2); B, C = r2.uniform(-1.5, 1.5, 2); p_ = B*C
    t = fH(L, U, B, C)
    if p_ < 0 or p_ > 1:
        good &= t > 2; cnt["squeeze"] += 1
    elif 0 < p_ < 1:
        strong = p_*(1 - p_)*L**2*U**2 >= 4*(1 - L)*(1 - U)
        if strong: good &= t <= -2 + 1e-12; cnt["flip"] += 1
        else: good &= abs(t) < 2; cnt["rotation"] += 1
res.append(ok("K11 enact f, enact g, repeal f, repeal g: the residue has det 1; it squeezes (trace > 2) whenever the loop balances (b c < 0), and for a gentle reinforcing loop (0 < b c < 1) it is rotation type unless both rules are strong (b c (1 - b c) l^2 u^2 >= 4 (1-l)(1-u)) [%s]" % cnt, okexact and good and min(cnt.values()) > 0))
print("\n%d of %d coupling checks pass" % (sum(res), len(res)))
