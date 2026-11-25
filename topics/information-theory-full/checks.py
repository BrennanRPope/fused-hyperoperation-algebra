"""Information theory, full pass, through mult-add: every claim on index.html, checked.
Run: python3 checks.py > checks_output.txt"""
import math, itertools
import numpy as np, sympy as sp, networkx as nx
from scipy.optimize import minimize_scalar

rng = np.random.default_rng(11)
results = []
def check(name, ok, detail=""):
    results.append(bool(ok)); print(f"{'PASS' if ok else 'FAIL'}  {name}" + (f"   [{detail}]" if detail else ""))
dial = lambda k, a, b: a + b + k*a*b

def lnq(x, q):  # Tsallis q-logarithm
    return math.log(x) if abs(q-1) < 1e-14 else (x**(1-q) - 1)/(1-q)
def eq(x, q):
    return math.exp(x) if abs(q-1) < 1e-14 else (1 + (1-q)*x)**(1/(1-q))

print("== 1. The q-logarithm turns products into the dial, k = 1 − q")
x, y, k, q = sp.symbols('x y k q', positive=True)
L = (x**k - 1)/k; Ly = (y**k - 1)/k; Lxy = ((x*y)**k - 1)/k
check("symbolic: ln_q(xy) = ln_q x ✳ₖ ln_q y with k = 1 − q", sp.simplify(sp.expand_power_base(Lxy - dial(k, L, Ly), force=True)) == 0)
ok = True
for _ in range(500):
    qq = rng.uniform(-2, 3); a, b = rng.uniform(.05, 5, 2); kk = 1 - qq
    ok &= abs(lnq(a*b, qq) - dial(kk, lnq(a, qq), lnq(b, qq))) < 1e-9*max(1, abs(lnq(a*b, qq)))
    if 1 + kk*lnq(a, qq) > 0 and 1 + kk*lnq(b, qq) > 0:
        u, v = lnq(a, qq), lnq(b, qq)
        ok &= abs(eq(dial(kk, u, v), qq) - eq(u, qq)*eq(v, qq)) < 1e-8*max(1, a*b)
        ok &= abs(eq(lnq(a, qq), qq) - a) < 1e-9*max(1, a)
check("500 random q: ln_q(ab) = dial; e_q(u ✳ₖ v) = e_q(u)·e_q(v); e_q inverts ln_q", ok)
check("k = −1 (Tsallis q = 2): ln_2 x = 1 − 1/x and (1 − 1/xy) = (1 − 1/x) ✳₋₁ (1 − 1/y)",
      sp.simplify((1 - 1/(x*y)) - dial(-1, 1 - 1/x, 1 - 1/y)) == 0)
check("k = 1 (q = 0): ln_0 x = x − 1, the 'count minus one' that mult-adds: 5·7 − 1 = 4 ✳ 6",
      abs(lnq(7.0, 0) - 6) < 1e-15 and dial(1, 4, 6) == 5*7 - 1)

print("== 2. Tsallis entropy is the expectation of a q-surprise; the measure engine composes it")
ok = okeng = True
for _ in range(300):
    qq = rng.uniform(-1, 3); kk = 1 - qq
    p, r = rng.dirichlet(np.ones(rng.integers(2, 6))), rng.dirichlet(np.ones(rng.integers(2, 6)))
    S = lambda d: sum(pi*lnq(1/pi, qq) for pi in d)
    ok &= abs(S(p) - (np.sum(p**qq) - 1)/(1 - qq)) < 1e-9 if abs(qq-1) > 1e-6 else True
    joint = np.outer(p, r).ravel()
    okeng &= abs(sum(pi*ri*dial(kk, lnq(1/pi, qq), lnq(1/ri, qq)) for pi in p for ri in r) - dial(kk, S(p), S(r))) < 1e-9
    okeng &= abs(S(joint) - dial(kk, S(p), S(r))) < 1e-9
check("S_q = Σ p·ln_q(1/p) = (Σp^q − 1)/(1 − q)", ok)
check("S(A×B) = E[s_A ✳ₖ s_B] = S(A) ✳ₖ S(B)  (Tsallis 1988, via E[X✳Y] = EX✳EY)", okeng)
# dependence breaks it by k·Cov
J = rng.dirichlet(np.ones(12)).reshape(3, 4); pA, pB = J.sum(1), J.sum(0); qq = 0.5; kk = .5
sA = np.array([lnq(1/v, qq) for v in pA]); sB = np.array([lnq(1/v, qq) for v in pB])
lhs = sum(J[i, j]*dial(kk, sA[i], sB[j]) for i in range(3) for j in range(4))
cov = sA@J@sB - (pA@sA)*(pB@sB)
check("correlated sources: E[s_A ✳ s_B] misses S_A ✳ S_B by exactly k·Cov(s_A, s_B)", abs(lhs - dial(kk, pA@sA, pB@sB) - kk*cov) < 1e-12)

print("== 3. Two q's: the project's q is Tsallis' dual index 2 − q")
xs = sp.symbols('xs', positive=True); qs = sp.symbols('qs')
lq = lambda X, Q: (X**(1-Q) - 1)/(1-Q)
check("ln_{2−q}(x) = −ln_q(1/x)  (the known additive duality q ↔ 2 − q)", sp.simplify(lq(xs, 2-qs) + lq(1/xs, qs)) == 0)
check("weights: k(2 − q) = q − 1 = −k(q), so the duality is the dial's sign flip", sp.simplify((1-(2-qs)) + (1-qs)) == 0)
check("project q_P = k + 1 = 2 − q_T for every Tsallis index", all(abs((1 - qt) + 1 - (2 - qt)) < 1e-15 for qt in np.linspace(-2, 3, 11)))

print("== 4. The Sharma–Mittal plane: every entropy there is Rényi pushed through the dial's exponential")
def renyi(p, a): return -np.sum(p*np.log(p)) if abs(a-1) < 1e-12 else math.log(np.sum(p**a))/(1-a)
def sm(p, a, b):
    if abs(b-1) < 1e-12: return renyi(p, a)
    if abs(a-1) < 1e-12: return (math.exp((1-b)*renyi(p, 1)) - 1)/(1-b)
    return (np.sum(p**a)**((1-b)/(1-a)) - 1)/(1-b)
ok = okc = True
for _ in range(300):
    a, b = rng.uniform(.2, 3, 2); kk = 1 - b
    p, r = rng.dirichlet(np.ones(4)), rng.dirichlet(np.ones(3)); j = np.outer(p, r).ravel()
    ok &= abs(sm(p, a, b) - (math.exp(kk*renyi(p, a)) - 1)/kk) < 1e-9
    okc &= abs(sm(j, a, b) - dial(kk, sm(p, a, b), sm(r, a, b))) < 1e-9
check("H_{α,β} = (e^{(1−β)·Rényi_α} − 1)/(1 − β)", ok)
check("H_{α,β}(A×B) = H(A) ✳ₖ H(B) with k = 1 − β, for every α", okc)
p = rng.dirichlet(np.ones(5)); a = 1.7
check("β = α is Tsallis", abs(sm(p, a, a) - (np.sum(p**a) - 1)/(1 - a)) < 1e-12)
check("β = 2 − α is Landsberg–Vedral (1/Σp^α − 1)/(α − 1), weight k = α − 1 (flipped)", abs(sm(p, a, 2-a) - (1/np.sum(p**a) - 1)/(a - 1)) < 1e-12)
check("β → 1 is Rényi; α = β → 1 is Shannon", abs(sm(p, a, 1 + 1e-7) - renyi(p, a)) < 1e-5 and abs(sm(p, 1+1e-7, 1+1e-7) - renyi(p, 1)) < 1e-5)

print("== 5. The q-exponential is compound interest; it is not the rung flow")
ok = True
for n in [1, 2, 5, 12, 100]:
    for xx in [0.3, 1.0, 2.5]:
        mp = 0.0
        for _ in range(n): mp = dial(1, mp, xx/n)
        ok &= abs(eq(xx, 1 - 1/n) - 1 - mp) < 1e-9*(1+mp)
check("k = 1/n: e_q(x) − 1 = n mult-add steps of x/n = (1 + x/n)^n − 1 (rung 2)", ok)
check("as k → 0, e_q(1) → e: 2.0, 2.5937, 2.7048, 2.7181",
      all(abs(eq(1, 1-1/n) - (1+1/n)**n) < 1e-12 for n in [1, 10, 100, 10000]) and abs(eq(1, 1-1e-6) - math.e) < 1e-5)
t = sp.symbols('t'); g = lambda X: (1 + X/2)**2 - 1
ser = sp.series(g(g(t)), t, 0, 4).removeO(); E = sp.series(sp.exp(t) - 1, t, 0, 4).removeO()
check("half-step test: g_{1/2}∘g_{1/2} matches e^x − 1 to x², fails at x³ (1/8 vs 1/6), so the family is not a flow",
      sp.expand(ser - E).coeff(t, 2) == 0 and ser.coeff(t, 3) == sp.Rational(1, 8) and E.coeff(t, 3) == sp.Rational(1, 6))

print("== 6. Trees: grafting composes internal nodes by the dial with k = D − 1")
def rand_full_tree(D, I):          # grow a full D-ary tree with I internal nodes
    G = nx.DiGraph(); G.add_node(0); leaves = [0]; nxt = 1
    for _ in range(I):
        v = leaves.pop(rng.integers(len(leaves)))
        for _ in range(D): G.add_edge(v, nxt); leaves.append(nxt); nxt += 1
    return G
def graft(A, B):
    G = A.copy(); off = max(G.nodes) + 1
    for leaf in [v for v in A.nodes if A.out_degree(v) == 0]:
        mp = {u: (leaf if u == 0 else off + u) for u in B.nodes}
        for (u, w) in B.edges: G.add_edge(mp[u], mp[w])
        off += max(B.nodes) + 1
    return G
inner = lambda G: sum(1 for v in G.nodes if G.out_degree(v) > 0)
leafc = lambda G: sum(1 for v in G.nodes if G.out_degree(v) == 0)
ok = okleaf = okkr = True
for _ in range(60):
    D = int(rng.integers(2, 5)); A = rand_full_tree(D, int(rng.integers(0, 6))); B = rand_full_tree(D, int(rng.integers(0, 6)))
    Gt = graft(A, B)
    ok &= inner(Gt) == dial(D-1, inner(A), inner(B))
    okleaf &= leafc(Gt) == leafc(A)*leafc(B) and all(leafc(T) == 1 + (D-1)*inner(T) for T in (A, B, Gt))
    depth = nx.shortest_path_length(Gt, 0)
    okkr &= abs(sum(D**-depth[v] for v in Gt.nodes if Gt.out_degree(v) == 0) - 1) < 1e-12
check("60 random D-ary trees (D = 2, 3, 4): inner(graft A B) = inner A ✳_{D−1} inner B", ok)
check("leaves multiply and leaves = 1 + (D−1)·inner (the lift)", okleaf)
check("Kraft sum of every full tree is exactly 1, and stays 1 after grafting", okkr)
check("complete depth-n tree has [n]_D = 1 + D + … + D^(n−1) internal nodes (q-number, q = D)",
      all(inner(nx.balanced_tree(D, n, create_using=nx.DiGraph)) == sum(D**j for j in range(n)) for D in [2, 3, 4] for n in range(1, 6)))

print("== 7. Any binary channel: the crossover sum a + b composes by the union rule")
ok = True
for _ in range(500):
    a1, b1, a2, b2 = rng.uniform(0, 1, 4)
    M = lambda a, b: np.array([[1-a, a], [b, 1-b]])
    C = M(a1, b1) @ M(a2, b2)
    ok &= abs((C[0, 1] + C[1, 0]) - dial(-1, a1+b1, a2+b2)) < 1e-12
    ok &= abs(np.linalg.det(C) - (1-a1-b1)*(1-a2-b2)) < 1e-12
check("500 random 2×2 channels: s = a + b of the cascade = s₁ ✳₋₁ s₂, and det = 1 − s multiplies", ok)
C = np.array([[.9,.1],[.2,.8]]) @ np.array([[.85,.15],[.25,.75]])
check("figure example: s₁ = 0.3, s₂ = 0.4 chain to 0.58 = 0.3 ✳₋₁ 0.4", abs(C[0,1] + C[1,0] - 0.58) < 1e-12 and abs(dial(-1, .3, .4) - .58) < 1e-12)
p_, q_ = sp.symbols('p q')
check("BSC (s = 2p): 2p ✳₋₁ 2q = 2(p ✳₋₂ q), the parity law", sp.expand(dial(-1, 2*p_, 2*q_) - 2*dial(-2, p_, q_)) == 0)
check("Z-channel (s = p): p ✳₋₁ q, the union law", all(abs((np.array([[1-a,a],[0,1]])@np.array([[1-c,c],[0,1]]))[0,1] - dial(-1,a,c)) < 1e-12 for a, c in rng.uniform(0,1,(50,2))))
def cap(Mt):
    h = lambda v: -sum(t*math.log2(t) for t in v if t > 0)
    f = lambda u: -(h(np.array([1-u, u])@Mt) - ((1-u)*h(Mt[0]) + u*h(Mt[1])))
    return -minimize_scalar(f, bounds=(1e-9, 1-1e-9), method='bounded').fun
caps = []; Mz = np.array([[1, 0], [0.3, 0.7]]); P = np.eye(2)
for n in range(1, 7): P = P @ Mz; caps.append(cap(P))
check("Z-channel cascade: capacity falls every hop while s climbs by the union law", all(a > b for a, b in zip(caps, caps[1:])), ", ".join(f"{c:.3f}" for c in caps))

print("== 8. Every landscape point gives an entropy (Tempesta's group entropies)")
ok = True
for _ in range(100):
    p, r = rng.dirichlet(np.ones(3)), rng.dirichlet(np.ones(4)); a = rng.uniform(.3, 2.5); c = 0.3
    S = lambda d: math.tan(c*renyi(d, a))          # tangent point (0, 1)
    j = np.outer(p, r).ravel()
    if c*(renyi(j, a)) < 1.5:
        ok &= abs(S(j) - (S(p) + S(r))/(1 - S(p)*S(r))) < 1e-9
check("tan(c·Rényi) composes by the tangent law (a + b)/(1 − ab), landscape point (0, 1)", ok)

print(f"\n{sum(results)}/{len(results)} checks pass")
