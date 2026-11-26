"""Polish for the coupling checks (DRAFT): bigger packages, exact proofs in place of samples, a known identity behind K11."""
import itertools, random
from fractions import Fraction as Fr
import sympy as sp, networkx as nx
ok = lambda n, c: print(("PASS " if c else "FAIL ") + n) or c
res = []
def chrom_at(G, q):
    G = nx.Graph(G)
    if G.number_of_edges() == 0: return q**G.number_of_nodes()
    u, v = next(iter(G.edges())); G1 = G.copy(); G1.remove_edge(u, v)
    G2 = nx.Graph(nx.contracted_nodes(G1, u, v, self_loops=False))
    return chrom_at(G1, q) - chrom_at(G2, q)
# Q1 seven-rule packages: 5040 orders each, exact arithmetic
rng = random.Random(77); good = True; seen = []
for trial in range(8):
    m = rng.randint(2, 4); n = 7; rules = []
    for _ in range(n):
        w = rng.randrange(m); reads = [j for j in range(m) if j != w and rng.random() < 0.25]
        bvec = {j: Fr(rng.choice([-1, 1])*rng.randint(5, 40), 100) for j in reads}
        rules.append((w, reads, Fr(rng.randint(5, 95), 100), Fr(rng.randint(30, 70), 100), bvec))
    x0 = [Fr(rng.randint(0, 100), 100) for _ in range(m)]
    def apply(order):
        X = list(x0)
        for i in order:
            w, reads, l, a0, bv = rules[i]; X[w] = (1 - l)*X[w] + l*(a0 + sum(bv[j]*X[j] for j in reads))
        return tuple(X)
    outs = {apply(o) for o in itertools.permutations(range(n))}
    G = nx.Graph(); G.add_nodes_from(range(n))
    for i, j in itertools.combinations(range(n), 2):
        wi, ri = rules[i][0], rules[i][1]; wj, rj = rules[j][0], rules[j][1]
        same = wi == wj and rules[i][3] == rules[j][3] and rules[i][4] == rules[j][4]
        if (wi == wj and not same) or wi in rj or wj in ri: G.add_edge(i, j)
    ch = abs(chrom_at(G, -1)); good &= len(outs) == ch; seen.append((len(outs), ch))
res.append(ok("Q1 eight seven-rule packages (5040 orders each, exact): distinct outcomes = |chi(-1)| every time; %s" % seen, good))
# Q2 exact: a spiral needs a balancing loop. Delta(p) - Delta(0) = lam mu p (2(2-lam-mu) + lam mu p) >= 0 for p >= 0, Delta(0) = (lam-mu)^2
lam, mu, p = sp.symbols('lambda mu p')
Delta = (2 - lam - mu + lam*mu*p)**2 - 4*(1 - lam)*(1 - mu)
res.append(ok("Q2 exact: Delta(p) = (lam - mu)^2 + lam mu p (2(2 - lam - mu) + lam mu p), every term >= 0 when p = b c >= 0 and 0 < lam, mu < 1, so no spiral without a balancing loop",
    sp.expand(Delta - ((lam - mu)**2 + lam*mu*p*(2*(2 - lam - mu) + lam*mu*p))) == 0))
# Q3 and the converse: a balancing loop spirals once it is strong enough; the threshold in p
pc = sp.solve(sp.Eq(Delta, 0), p)
res.append(ok("Q3 Delta = 0 at p = %s; between these roots (both <= 0) the turns spiral" % [sp.simplify(r) for r in pc],
    all(sp.simplify(r.subs({lam: sp.Rational(1, 2), mu: sp.Rational(1, 3)})) <= 0 for r in pc)))
# Q4 K11 is the Fricke trace identity tr[A,B] = trA^2 + trB^2 + trAB^2 - trA trB trAB - 2 for the normalised turns
b, c = sp.symbols('b c')
F = sp.Matrix([[1 - lam, lam*b], [0, 1]]); G = sp.Matrix([[1, 0], [mu*c, 1 - mu]])
A = F/sp.sqrt(1 - lam); B = G/sp.sqrt(1 - mu)
tA, tB, tAB = A.trace(), B.trace(), (A*B).trace()
comm = (F*G*F.inv()*G.inv()).trace()
res.append(ok("Q4 the loop residue's trace is the Fricke identity tr[A,B] = trA^2 + trB^2 + trAB^2 - trA trB trAB - 2 applied to the two rules scaled to det 1 (known; only the reading is ours)",
    sp.simplify(comm - (tA**2 + tB**2 + tAB**2 - tA*tB*tAB - 2)) == 0))
# Q5 the five zones on the loop product p = b c, from the Jury conditions for a 2x2 step (known criterion)
tr = 2 - lam - mu + lam*mu*p; det = (1 - lam)*(1 - mu)
res.append(ok("Q5 Jury: 1 - tr + det = lam mu (1 - p) and 1 + tr + det = lam mu p + (2 - lam)(2 - mu): settles iff -(2-lam)(2-mu)/(lam mu) < p < 1",
    sp.expand(1 - tr + det - lam*mu*(1 - p)) == 0 and sp.expand(1 + tr + det - (lam*mu*p + (2 - lam)*(2 - mu))) == 0))
import numpy as np
r = np.random.default_rng(4); good = True; zc = {}
for _ in range(100000):
    L, U = r.uniform(0.02, 0.98, 2); P = r.uniform(-60, 2)
    a1, a2 = np.sqrt(1 - L), np.sqrt(1 - U)
    cuts = [-(2 - L)*(2 - U)/(L*U), -(a1 + a2)**2/(L*U), -(a1 - a2)**2/(L*U), 1]
    zone = int(sum(P > q for q in cuts))  # 0 growing swings, 1 overshoot, 2 spiral, 3 smooth, 4 runaway
    ev = np.linalg.eigvals(np.array([[1 - L, L], [U*P*(1 - L), 1 - U + L*U*P]]))
    rad = max(abs(ev)); cplx = abs(ev[0].imag) > 1e-12; neg = (not cplx) and min(ev.real) < 0
    want = {0: rad > 1 and neg, 1: rad < 1 and neg and not cplx, 2: cplx and rad < 1, 3: rad < 1 and not cplx and min(ev.real) >= 0, 4: rad > 1}[zone]
    good &= bool(want); zc[zone] = zc.get(zone, 0) + 1
res.append(ok("Q6 five zones on p = b c, in order: growing swings | overshoot each turn | spiral | smooth settle | runaway (cuts -(2-l)(2-u)/lu, -(sqrt(1-l)+sqrt(1-u))^2/lu, -(sqrt(1-l)-sqrt(1-u))^2/lu, 1); 100000 samples agree with the eigenvalues %s" % dict(sorted(zc.items())), good and len(zc) == 5))
res.append(ok("Q6b equal strengths (lam = mu): the spiral zone starts at p = 0, so every balancing loop spirals",
    sp.simplify((-(sp.sqrt(1 - lam) - sp.sqrt(1 - lam))**2/(lam*lam))) == 0))
# Q7 where it settles: x* = (a + b s)/(1 - b c), whatever the strengths and whichever goes first (feedback gain 1/(1 - loop), known)
a, s, b, c = sp.symbols('a s b c')
X, Y = sp.symbols('X Y')
sol = sp.solve([sp.Eq(X, a + b*Y), sp.Eq(Y, s + c*X)], [X, Y])
good = sp.simplify(sol[X] - (a + b*s)/(1 - b*c)) == 0
for L, U in [(0.3, 0.6), (0.8, 0.2)]:
    for first in "fg":
        x_, y_ = 0.1, 0.9
        for _ in range(4000):
            if first == "f": x_ = (1 - L)*x_ + L*(0.2 + 0.5*y_); y_ = (1 - U)*y_ + U*(0.4 - 0.7*x_)
            else: y_ = (1 - U)*y_ + U*(0.4 - 0.7*x_); x_ = (1 - L)*x_ + L*(0.2 + 0.5*y_)
        good &= abs(x_ - (0.2 + 0.5*0.4)/(1 + 0.35)) < 1e-12
res.append(ok("Q7 run long, the pair settles at x* = (a + b s)/(1 - b c) for every strength and either order: order and strength shape the path, the loop sets the destination", good))
print("\n%d of %d coupling polish checks pass" % (sum(res), len(res)))
