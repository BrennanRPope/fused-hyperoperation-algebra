"""Checks for topics/union-find: union-find, the inverse of the hyperoperation tower, and Asano's merges
as the union steps of the Lee-Yang circle theorem on the dial gas.

Run: python3 checks.py   (numpy, sympy, mpmath). Prints PASS/FAIL per check and a total."""
import math, random
from fractions import Fraction
import numpy as np
import sympy as sp
from ufcore import (H, ack, fused, alpha_k, ack_diag_alpha, tower, gas_roots, components, UF, workload, cycle, grid, petersen, complete,
                    random_graph, degrees, gas_coeffs, gas_coeffs_uniform, roots, merge_stages, gas_coeffs_exact)

results = []


def check(name, ok, detail=""):
    results.append(bool(ok))
    print(("PASS " if ok else "FAIL ") + name + (("  [" + detail + "]") if detail else ""))


# ---------------------------------------------------------------- A. Ackermann is a row of the ladder
print("A. Union-find's cost and the hyperoperation tower")
okA = all(ack(n, m) == H(n, 2, m + 3) - 3 for n in range(4) for m in range(7)) and \
    all(ack(4, m) == H(4, 2, m + 3) - 3 for m in range(2))
check("A1 Ackermann A(n, m) = 2 [n] (m + 3) - 3 (rungs 0-3, m <= 6; rung 4, m <= 1: A(4,1) = 2^^4 - 3 = 65533)", okA,
      f"A(4,1) = {ack(4, 1)}")
check("A2 so A(n, m) = (fused rung n)(1, m + 2) - 2: Ackermann is the fused ladder at a = 1",
      all(ack(n, m) == fused(n, 1, m + 2) - 2 for n in range(4) for m in range(7)))
check("A3 each alpha_k undoes one rung: alpha_1(2m) = m, alpha_2(2^m) = m, alpha_3(2^^m) = m (m <= 5)",
      all(alpha_k(1, 2 * m) == m and alpha_k(2, 2 ** m) == m for m in range(1, 60)) and
      all(alpha_k(3, tower(m)) == m for m in range(1, 6)))
N80 = 10 ** 80
ks = [alpha_k(k, N80) for k in (2, 3, 4, 5)]
check("A4 for n = 10^80 (about the atoms in the universe): alpha_2 = 266, alpha_3 = log* = 5, alpha_4 = 4, alpha_5 = 3",
      ks == [266, 5, 4, 3], f"{ks}")
check("A5 inverse Ackermann alpha(n) = least k with A(k, k) >= n is 4 for 10^80 and for 2^65535 (A(3,3) = 61, A(4,2) = 2^65536 - 3)",
      ack(3, 3) == 61 and ack_diag_alpha(N80) == 4 and ack_diag_alpha(2 ** 65535) == 4)

# ---------------------------------------------------------------- B. union-find: order-blind result, order-sensitive cost
print("B. Union-find: the partition ignores order, the cost does not")
rnd = random.Random(3)
n = 60
U = [(rnd.randrange(n), rnd.randrange(n)) for _ in range(45)]
parts, depths = set(), []
for s in range(40):
    W = U[:]; random.Random(s).shuffle(W)
    W = [(y, x) if random.Random(1000 + s + i).random() < .5 else (x, y) for i, (x, y) in enumerate(W)]
    uf = UF(n, rank=False, compress=False)
    for x, y in W:
        uf.union(x, y)
    parts.add(tuple(uf.blocks())); depths.append(uf.depth())
check("B1 40 random orders (and directions) of the same 45 unions give one partition", len(parts) == 1)
check("B2 ... but different trees: naive linking depth ranges over several values", len(set(depths)) >= 3,
      f"depths {min(depths)}-{max(depths)}")
chain = UF(64, rank=False, compress=False)
for i in range(1, 64):
    chain.union(i - 1, i)
check("B3 worst order for naive linking: a chain of depth n - 1 (n = 64)", chain.depth() == 63)
rk = UF(1 << 12, rank=True, compress=False)
r2 = random.Random(5)
for _ in range(1 << 13):
    rk.union(r2.randrange(1 << 12), r2.randrange(1 << 12))
check("B4 union by rank alone keeps depth <= log2 n (n = 4096)", rk.depth() <= 12, f"depth {rk.depth()}")
small, big = workload(1 << 10, True, True, 1), workload(1 << 17, True, True, 1)
check("B5 rank + path compression: under 2 pointer steps per operation, flat from n = 2^10 to 2^17",
      big < 2 and big / small < 1.5, f"{small:.3f} -> {big:.3f}")
nv1, nv2 = workload(1 << 7, False, False, adversarial=True), workload(1 << 10, False, False, adversarial=True)
check("B6 naive linking on the chain: cost per operation grows like n (x8 in n gives about x8 in cost)",
      6 < nv2 / nv1 < 10, f"{nv1:.1f} -> {nv2:.1f}")

# ---------------------------------------------------------------- C. Asano's merges and the Lee-Yang circle
print("C. Lee-Yang on the dial gas, built by union steps (Asano)")
x_, y_, a_ = sp.symbols('x y a', real=True)
w_ = x_ + sp.I * y_
gap = sp.expand(sp.Abs(w_ + a_) ** 2 - sp.Abs(1 + a_ * w_) ** 2)
gap = sp.expand((x_ + a_) ** 2 + y_ ** 2 - (1 + a_ * x_) ** 2 - (a_ * y_) ** 2)
check("C1 one edge's geometry: |w + a|^2 - |1 + a w|^2 = (1 - a^2)(|w|^2 - 1), so for a <= 1 the map "
      "w -> -(1 + a w)/(w + a) sends the outside of the circle inside",
      sp.expand(gap - (1 - a_ ** 2) * (x_ ** 2 + y_ ** 2 - 1)) == 0)
k_, z1, z2, w1, w2 = sp.symbols('k z1 z2 w1 w2', positive=True)
edge = 1 + z1 + z2 + k_ * z1 * z2
resc = sp.expand(edge.subs({z1: w1 / sp.sqrt(k_), z2: w2 / sp.sqrt(k_)}))
check("C2 the dial edge 1 + z1 + z2 + k z1 z2 (= 1 + z1 (+)_k z2) at z = w/sqrt(k) is 1 + (w1 + w2)/sqrt(k) + w1 w2:"
      " the Lee-Yang edge with a = k^(-1/2), so a <= 1 exactly when k >= 1",
      sp.simplify(resc - (1 + (w1 + w2) / sp.sqrt(k_) + w1 * w2)) == 0)
# contraction identity on a triangle with a tail (4 vertices, 4 edges)
Ng, Eg = 4, [(0, 1), (1, 2), (0, 2), (2, 3)]
cv = sp.symbols('c0:8'); av = sp.Symbol('a')
P = 1
for i, (u, v) in enumerate(Eg):
    P *= (1 + av * cv[2 * i] + av * cv[2 * i + 1] + cv[2 * i] * cv[2 * i + 1])
P = sp.Poly(sp.expand(P), *cv)
owner = [u for e in Eg for u in e]
W = sp.symbols('W0:4')
con = 0
for mon, coef in zip(P.monoms(), P.coeffs()):
    inS = {}
    okm = True
    for c, e in enumerate(mon):
        if inS.setdefault(owner[c], e) != e:
            okm = False; break
    if okm:
        term = coef
        for v in range(Ng):
            if inS.get(v, 0):
                term *= W[v]
        con += term
target = sum(sp.Mul(*[W[v] for v in range(Ng) if S >> v & 1]) *
             av ** sum(1 for u, v in Eg if (S >> u & 1) != (S >> v & 1)) for S in range(16))
check("C3 Asano contraction (keep a vertex's copies all-in or all-out, then glue them) of the product of edge factors"
      " gives sum_S a^(edges leaving S) w^S, the dial gas at z_v = k^(-d_v/2) w_v", sp.expand(con - target) == 0)
graphs = {'ring C8': cycle(8), 'Petersen': petersen(), 'grid 3x3': grid(3, 3), 'K5': complete(5),
          'random A': random_graph(10, .3, 1), 'random B': random_graph(11, .25, 2), 'random C': random_graph(9, .5, 3)}
worst = 0
for nm, (N, E) in graphs.items():
    for k in (1.5, 3.0, 10.0):
        rs = gas_roots(N, E, k); worst = max(worst, max(abs(abs(r) - 1) for r in rs))
check("C4 Lee-Yang circle: for k = 1.5, 3, 10 on 7 graphs (ring, Petersen, grid, K5, three random), every zero has |w| = 1",
      worst < 1e-9, f"max ||w| - 1| = {worst:.1e}")
Np, Ep = petersen(); Nt, Et = grid(4, 4, torus=True)
rp = roots(gas_coeffs_uniform(Np, Ep, 3.0)); rt = roots(gas_coeffs_uniform(Nt, Et, 2.0))
check("C5 regular graphs give one circle in z: |z| = k^(-d/2) (Petersen d = 3, k = 3; 4x4 torus d = 4, k = 2 gives 1/4,"
      " the lee-yang page's k^-2)", max(abs(abs(r) - 3 ** -1.5) for r in rp) < 1e-9 and max(abs(abs(r) - .25) for r in rt) < 1e-9)
c1 = gas_coeffs_exact(*petersen(), Fraction(1))
check("C6 at k = 1 (the free gas) the gas is exactly (1 + w)^N: every zero at w = -1 (Petersen, N = 10)",
      c1 == [Fraction(math.comb(10, m)) for m in range(11)])
offs = [max(abs(abs(r) - 1) for r in gas_roots(N, E, k)) for (N, E) in (cycle(8), petersen()) for k in (0.5, 0.2)]
check("C7 below k = 1 (repulsive) zeros leave the circle (ring C8 and Petersen at k = 0.5, 0.2)", min(offs) > 1e-2,
      f"min deviation {min(offs):.3f}")
stage_worst, degs_ok = 0, True
last = None
for (n_, E_), i in merge_stages(*petersen(), order_seed=7):
    rs = gas_roots(n_, E_, 3.0); stage_worst = max(stage_worst, max(abs(abs(r) - 1) for r in rs))
    if last is not None and n_ != last - 1:
        degs_ok = False
    last = n_
check("C8 Petersen built by 20 union steps from 15 loose edges (30 copies), k = 3: zeros on |w| = 1 at every stage, "
      "one fewer zero per union", stage_worst < 1e-9 and degs_ok and last == 10, f"max dev {stage_worst:.1e}")
finals = set()
for s in range(10):
    for (n_, E_), i in merge_stages(*grid(3, 3), order_seed=100 + s):
        pass
    finals.add(tuple(gas_coeffs_exact(n_, E_, Fraction(1, 2))))
check("C9 10 random union orders build the same final gas, exactly (grid 3x3, a = 1/2)", len(finals) == 1)

# ---------------------------------------------------------------- D. which moves keep Lee-Yang
print("D. Products and unions keep the circle; sums do not")
P1 = np.poly1d([1, 1]) ** 4; P2 = np.poly1d([1, -1]) ** 4
rsum = np.roots((P1 + P2).coeffs)
check("D1 a sum (mixture) breaks it: (1 + w)^4 + (1 - w)^4 has zeros at |w| = sqrt2 - 1 and sqrt2 + 1, off the circle",
      abs(min(abs(rsum)) - (2 ** .5 - 1)) < 1e-9 and abs(max(abs(rsum)) - (2 ** .5 + 1)) < 1e-9)
rpr = gas_roots(*cycle(5), 2.0) + gas_roots(*complete(4), 4.0)   # a product: its zeros are the two lists together
check("D2 a product (two systems side by side, no union) keeps every zero on the circle", max(abs(abs(r) - 1) for r in rpr) < 1e-9)

print()
print(f"{sum(results)} of {len(results)} pass")
