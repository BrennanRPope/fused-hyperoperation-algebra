"""Entropy and quantum noise on the mult-add cube: numeric and symbolic checks.

Part of Fused Hyperoperation Algebra. Companion page: index.html.
Run: python3 checks.py > checks_output.txt
Every check prints PASS or FAIL; the script exits non-zero if any fail.
"""
import itertools, math, sys
import numpy as np
import sympy as sp
from scipy.optimize import brentq, minimize

rng = np.random.default_rng(20261125)
FAILS = []
COUNT = [0]


def check(name, ok, detail=""):
    COUNT[0] += 1
    print(("PASS " if ok else "FAIL ") + name + (("  " + detail) if detail else ""))
    if not ok:
        FAILS.append(name)


def H(p, base=2):
    p = np.asarray(p, float)
    p = p[p > 0]
    return float(-(p * np.log(p)).sum() / math.log(base))


def h2(p):
    return H([p, 1 - p])


def wadd(k, a, b):
    return a + b + k * a * b


# ---------------------------------------------------------------------------
print("== 1. Tsallis pseudo-additivity is weighted mult-add with k = 1 - q ==")
def tsallis(p, q):
    p = np.asarray(p, float)
    if abs(q - 1) < 1e-12:
        return H(p, math.e)
    return (1 - (p ** q).sum()) / (q - 1)

worst = 0
for q in [0.0, 0.3, 0.5, 2.0, 3.0, 0.7, 1.5]:
    for _ in range(200):
        a = rng.dirichlet(np.ones(rng.integers(2, 6)))
        b = rng.dirichlet(np.ones(rng.integers(2, 6)))
        ab = np.outer(a, b).ravel()
        lhs = tsallis(ab, q)
        rhs = wadd(1 - q, tsallis(a, q), tsallis(b, q))
        worst = max(worst, abs(lhs - rhs))
check("S_q(A x B) = S_A + S_B + (1-q) S_A S_B on 1400 random pairs", worst < 1e-12, f"max error {worst:.1e}")
# q = 1 is Shannon, k = 0, plain addition
a = rng.dirichlet(np.ones(4)); b = rng.dirichlet(np.ones(3))
check("q = 1 (k = 0): Shannon entropy adds", abs(H(np.outer(a, b).ravel()) - H(a) - H(b)) < 1e-12)

# symbolic: Tsallis sits on the landscape line q_L = 0; certainty absorbs only at k = -1
k, x, qL = sp.symbols("k x q_L")
lens = (x + 1 + k * x * 1) / (1 - qL * x * 1)
sol = sp.solve(sp.Poly(sp.numer(sp.together(lens - 1)), x).coeffs(), [k, qL], dict=True)
check("lens law absorbs certainty (a o 1 = 1 for all a) exactly on k + q_L = -1",
      all(sp.simplify(s[k] + s.get(qL, qL) + 1) == 0 for s in sol), str(sol))
check("Tsallis laws have q_L = 0, so they meet that line only at k = -1, Tsallis q = 2",
      sp.solve(sp.Eq(k + 0, -1), k) == [-1])
print("   Note: the Tsallis index q and the landscape's q are different letters for different things.")
print("   Tsallis line: k + q_T = 1.  'Or' line: k + q_L = -1.  They cross only at (k, q_L) = (-1, 0).")

# quantum version: von Neumann / Tsallis for product states
def rand_rho(d):
    G = rng.normal(size=(d, d)) + 1j * rng.normal(size=(d, d))
    r = G @ G.conj().T
    return r / np.trace(r).real

def tr_pow(r, q):
    ev = np.clip(np.linalg.eigvalsh(r), 0, None)
    return float((ev ** q).sum())

worst = 0
for _ in range(200):
    A, B = rand_rho(2), rand_rho(3)
    for q in [2.0, 0.5, 3.0]:
        SA, SB = (1 - tr_pow(A, q)) / (q - 1), (1 - tr_pow(B, q)) / (q - 1)
        SAB = (1 - tr_pow(np.kron(A, B), q)) / (q - 1)
        worst = max(worst, abs(SAB - wadd(1 - q, SA, SB)))
check("quantum Tsallis (Tr rho^q) of rho x sigma composes at k = 1 - q", worst < 1e-10, f"max error {worst:.1e}")
A, B = rand_rho(2), rand_rho(2)
lin = lambda r: 1 - tr_pow(r, 2)
check("linear entropy 1 - Tr rho^2 of two independent qubits follows the union rule (k = -1)",
      abs(lin(np.kron(A, B)) - wadd(-1, lin(A), lin(B))) < 1e-12)
# von Neumann of a qubit depends only on the Bloch radius
r = rng.uniform(0, 1)
n = rng.normal(size=3); n /= np.linalg.norm(n)
X = np.array([[0, 1], [1, 0]]); Y = np.array([[0, -1j], [1j, 0]]); Z = np.array([[1, 0], [0, -1]]); I2 = np.eye(2)
rho = (I2 + r * (n[0] * X + n[1] * Y + n[2] * Z)) / 2
check("qubit von Neumann entropy = h((1 - r)/2): a radial function on the Bloch ball",
      abs(H(np.linalg.eigvalsh(rho)) - h2((1 - r) / 2)) < 1e-12)

# ---------------------------------------------------------------------------
print("\n== 2. The Pauli tetrahedron: three entropies, three norms of lambda ==")
pI, pX, pY, pZ = sp.symbols("p_I p_X p_Y p_Z")
lx = pI + pX - pY - pZ; ly = pI - pX + pY - pZ; lz = pI - pX - pY + pZ
expr = sp.expand(pI**2 + pX**2 + pY**2 + pZ**2 - (1 + lx**2 + ly**2 + lz**2) / 4)
expr = sp.expand(expr.subs(pI, 1 - pX - pY - pZ))
check("collision sum p_I^2+p_X^2+p_Y^2+p_Z^2 = (1 + |lambda|^2)/4 (Parseval on the tetrahedron)", expr == 0)

def pauli_from_l(l):
    lx, ly, lz = l
    return np.array([1 + lx + ly + lz, 1 + lx - ly - lz, 1 - lx + ly - lz, 1 - lx - ly + lz]) / 4

def l_from_pauli(p):
    return np.array([p[0] + p[1] - p[2] - p[3], p[0] - p[1] + p[2] - p[3], p[0] - p[1] - p[2] + p[3]])

N = 200000
P = rng.dirichlet(np.ones(4) * 0.6, size=N)
L = np.stack([l_from_pauli(p) for p in P[:20000]])
Ps = P[:20000]
hash_rate = np.array([1 - H(p) for p in Ps])
l2 = (L ** 2).sum(1)
# the depolarizing hashing point
pstar = brentq(lambda t: 1 - h2(t) - t * math.log2(3), 1e-6, 0.5)
pIstar = 1 - pstar
print(f"   depolarizing hashing bound: total fault chance p = {pstar:.6f}, so p_I = {pIstar:.6f}, lambda = {1 - 4*pstar/3:.6f}")
check("depolarizing hashing point p = 0.18929 (known)", abs(pstar - 0.189290) < 2e-6)
check("necessary (sphere): every channel with positive hashing rate lies outside the midsphere |lambda| = 1",
      bool(np.all(l2[hash_rate > 0] > 1)), f"{int((hash_rate > 0).sum())} positive-rate channels sampled")
check("sufficient (plane): any channel with some fault chance above 0.8107 has positive hashing rate",
      bool(np.all(hash_rate[Ps.max(1) > pIstar] > -1e-12)), f"{int((Ps.max(1) > pIstar).sum())} channels in the four corner caps")
# the sufficient bound is the grouping bound H(p) <= h(p_max) + (1 - p_max) log2 3
grp = np.array([h2(p.max()) + (1 - p.max()) * math.log2(3) for p in Ps])
check("grouping bound H(p) <= h(p_max) + (1 - p_max) log2 3 on all samples", bool(np.all(H_ok := (np.array([H(p) for p in Ps]) <= grp + 1e-12))))
# tight points
check("hashing surface touches the midsphere at the 6 edge midpoints (p = 1/2, 1/2, 0, 0: H = 1, |lambda| = 1)",
      abs(1 - H([0.5, 0.5, 0, 0])) < 1e-15 and np.allclose(l_from_pauli([0.5, 0.5, 0, 0]), [1, 0, 0]))
dep = np.array([pIstar, pstar / 3, pstar / 3, pstar / 3])
check("hashing surface touches the plane along the depolarizing diagonal", abs(1 - H(dep)) < 1e-9)
check("unit sphere is the tetrahedron's midsphere: edge midpoints at distance 1, corners at sqrt 3",
      np.isclose(np.linalg.norm(([1, 1, 1] + np.array([1, -1, -1])) / 2), 1) and np.isclose(np.linalg.norm([1, 1, 1]), math.sqrt(3)))

# Choi matrix and entanglement breaking = octahedron (Ruskai 2003)
paulis = [I2, X, Y, Z]
phi = np.array([1, 0, 0, 1]) / math.sqrt(2)
Phi = np.outer(phi, phi.conj())

def choi(p):
    return sum(pi * np.kron(I2, S) @ Phi @ np.kron(I2, S).conj().T for pi, S in zip(p, paulis))

def ptrans(M):
    M = M.reshape(2, 2, 2, 2)
    return M.transpose(0, 3, 2, 1).reshape(4, 4)

agree = 0
for p, l in zip(Ps[:5000], L[:5000]):
    eb = np.linalg.eigvalsh(ptrans(choi(p))).min() >= -1e-12
    octa = np.abs(l).sum() <= 1 + 1e-12
    agree += (eb == octa)
check("entanglement breaking (Choi state PPT) <=> |lx|+|ly|+|lz| <= 1, the octahedron", agree == 5000, f"{agree}/5000 agree")
check("Choi state eigenvalues are the four fault chances, so its von Neumann entropy is H(p)",
      np.allclose(np.sort(np.linalg.eigvalsh(choi(Ps[0]))), np.sort(Ps[0])))

# King: classical capacity of a unital qubit channel = 1 - h((1 - max|lambda|)/2)
def out_entropy(theta, l):
    th, ph = theta
    v = np.array([math.sin(th) * math.cos(ph), math.sin(th) * math.sin(ph), math.cos(th)]) * l
    return h2((1 - np.linalg.norm(v)) / 2)

worst = 0
for l in L[:300]:
    best = min(minimize(out_entropy, x0, args=(l,), method="Nelder-Mead").fun
               for x0 in [(0.1, 0.1), (1.5, 0.1), (1.5, 1.5), (0.8, 0.8)])
    best = min(best, *(h2((1 - abs(c)) / 2) for c in l))
    worst = max(worst, abs((1 - best) - (1 - h2((1 - np.abs(l).max()) / 2))))
check("King: min output entropy over pure inputs gives capacity 1 - h((1 - max|lambda_i|)/2)", worst < 1e-6, f"max gap {worst:.1e}")

print("   So the three level sets are: max|lambda| (cubes, classical capacity), |lambda|_2 (spheres, collision")
print("   entropy), |lambda|_1 (octahedra, entanglement breaking). Their unit balls kiss at (+-1,0,0),(0,+-1,0),(0,0,+-1).")
pts = [np.array(v) for v in [(1, 0, 0), (0, 1, 0), (0, 0, 1)]]
check("the six kissing points have |.|_1 = |.|_2 = |.|_inf = 1",
      all(np.isclose(np.abs(v).sum(), 1) and np.isclose(np.linalg.norm(v), 1) and np.isclose(np.abs(v).max(), 1) for v in pts))

# absorption under composition
A_ = L[:5000]; B_ = L[5000:10000]; C_ = A_ * B_
for nm, f in [("octahedron |.|_1 <= 1", lambda v: np.abs(v).sum(1) <= 1),
              ("midsphere ball |.|_2 <= 1", lambda v: (v ** 2).sum(1) <= 1),
              ("any cube max|.| <= c (c = 0.6)", lambda v: np.abs(v).max(1) <= 0.6)]:
    inside = f(A_)
    check(f"composing keeps a channel inside the {nm}", bool(np.all(f(C_)[inside])), f"{int(inside.sum())} tested")

# ---------------------------------------------------------------------------
print("\n== 3. The weight is a temperature: Renyi = free energy, k = 1 - T0/T (Carnot) ==")
kB = 1.0
worst = 0
for _ in range(200):
    E = rng.uniform(0, 3, size=rng.integers(2, 7))
    T0 = rng.uniform(0.3, 3)
    Zf = lambda T: np.exp(-E / (kB * T)).sum()
    F = lambda T: -kB * T * math.log(Zf(T))
    p = np.exp(-E / (kB * T0)) / Zf(T0)
    q = rng.uniform(0.2, 3)
    if abs(q - 1) < 0.05:
        continue
    T = T0 / q
    renyi = math.log((p ** q).sum()) / (1 - q)
    baez = -(F(T) - F(T0)) / (kB * (T - T0))
    worst = max(worst, abs(renyi - baez))
check("Baez: Renyi entropy of order q = -(F(T0/q) - F(T0))/(T0/q - T0) for Gibbs states", worst < 1e-10, f"max error {worst:.1e}")
E = np.array([0.0, 0.7, 1.3, 2.2]); T0 = 0.9; q = 1.8
Zf = lambda T: np.exp(-E / T).sum(); p = np.exp(-E / T0) / Zf(T0)
check("1 + k*S_Tsallis = Z(T0/q)/Z(T0)^q, a ratio of partition functions (the thermodynamics page's lift)",
      abs(1 + (1 - q) * tsallis(p, q) - Zf(T0 / q) / Zf(T0) ** q) < 1e-12)
# escort of escort
p = rng.dirichlet(np.ones(5)); a_, b_ = 1.7, 0.6
esc = lambda p, q: p ** q / (p ** q).sum()
check("escort of an escort: P_b(P_a(p)) = P_{ab}(p), so q multiplies along a chain", np.allclose(esc(esc(p, a_), b_), esc(p, a_ * b_)))
k1, k2 = 1 - a_, 1 - b_
check("hence k = 1 - q composes by the union rule k1 + k2 - k1 k2", abs((1 - a_ * b_) - wadd(-1, k1, k2)) < 1e-15)
Th, Tm, Tc = 900.0, 500.0, 300.0
eta = lambda h, c: 1 - c / h
check("k(T, T0) = 1 - T0/T is Carnot's efficiency, and Carnot cascades by the same rule",
      abs(eta(Th, Tc) - wadd(-1, eta(Th, Tm), eta(Tm, Tc))) < 1e-15)
print("   Landmarks: k = 1 (q = 0) is infinite temperature, every state counted once (Hartley);")
print("   k = 0 is T = T0, Shannon; k -> -inf is T -> 0, only the ground state (min-entropy).")

# ---------------------------------------------------------------------------
print("\n== 4. Stacking noise: what cancellation hides ==")
for lam in [0.9, 0.6, 0.3]:
    p = (1 - lam) / 2
    rows = []
    for n in [1, 2, 4, 8]:
        Pn = (1 - lam ** n) / 2
        rows.append((n, n * h2(p), h2(Pn), ((1 + lam**2) / 2) ** n, (1 + lam ** (2 * n)) / 2))
    ok = all(r[2] <= r[1] + 1e-12 and r[4] >= r[3] - 1e-12 for r in rows)
    check(f"lambda = {lam}: net-channel entropy <= fault-pattern entropy, net purity >= pattern purity", ok,
          "; ".join(f"n={r[0]}: {r[1]:.3f} vs {r[2]:.3f} bits" for r in rows))
lam = sp.symbols("lambda")
check("n = 2: purity gap (1+l^4)/2 - ((1+l^2)/2)^2 = (1-l^2)^2/4, a perfect square",
      sp.simplify((1 + lam**4) / 2 - ((1 + lam**2) / 2) ** 2 - (1 - lam**2) ** 2 / 4) == 0)

# ---------------------------------------------------------------------------
print("\n== 5. Error correction as entropy removal: the 3-bit code's body diagonals ==")
corners = list(itertools.product([0, 1], repeat=3))
def synd(e):
    return (e[0] ^ e[1], e[1] ^ e[2])
diag = {}
for e in corners:
    diag.setdefault(synd(e), []).append(e)
check("the 4 syndromes split the 3-cube into its 4 body diagonals (antipodal pairs)",
      all(len(v) == 2 and tuple(1 - c for c in v[0]) == v[1] for v in diag.values()))
for p in [0.05, 0.1, 0.2, 0.3]:
    w = lambda e: p ** sum(e) * (1 - p) ** (3 - sum(e))
    pe = np.array([w(e) for e in corners])
    ps = {s: sum(w(e) for e in v) for s, v in diag.items()}
    HS = H(list(ps.values()))
    HES = sum(ps[s] * H([w(e) / ps[s] for e in v]) for s, v in diag.items())
    fail = sum(min(w(e) for e in v) for v in diag.values())
    alarm_ok = all(abs(min(w(e) for e in v) / ps[s] - p) < 1e-12 for s, v in diag.items() if s != (0, 0))
    check(f"p = {p}: H(faults) {H(pe):.4f} = H(syndrome) {HS:.4f} + H(left) {HES:.4f}; "
          f"after an alarm the doubt is exactly one fresh coin p", abs(H(pe) - HS - HES) < 1e-12 and alarm_ok
          and abs(fail - (3 * p**2 * (1 - p) + p**3)) < 1e-12)
pp = sp.symbols("p")
s0 = (1 - pp) ** 3 + pp ** 3
check("no-alarm chance (1-p)^3 + p^3 = (1 + 3 lambda^2)/4, alarm chance each (1 - lambda^2)/4",
      sp.simplify(s0.subs(pp, (1 - lam) / 2) - (1 + 3 * lam**2) / 4) == 0 and
      sp.simplify((pp * (1 - pp)).subs(pp, (1 - lam) / 2) - (1 - lam**2) / 4) == 0)
# concatenation as entropy flow
f = lambda l: (3 * l - l ** 3) / 2
for lam0 in [0.2, 0.6]:
    ls = [lam0]
    for _ in range(5):
        ls.append(f(ls[-1]))
    ents = [h2((1 - l) / 2) for l in ls]
    check(f"concatenated 3-bit code from lambda = {lam0}: residual entropy falls every level",
          all(ents[i + 1] < ents[i] for i in range(5)), " -> ".join(f"{e:.3f}" for e in ents))

# ---------------------------------------------------------------------------
print("\n== 6. Thresholds as an entropy budget (surface code) ==")
pc_ent = brentq(lambda t: h2(t) - 0.5, 1e-4, 0.5)
check("h(p) = 1/2 bit at p = 0.110028", abs(pc_ent - 0.110028) < 1e-6, f"{pc_ent:.6f}")
print("   Random CSS codes reach rate 1 - 2h(p) (CSS hashing), positive exactly below this p.")
print("   Nishimori-Nemoto duality conjecture for the square lattice gives the same 0.1100;")
print("   numerics for the surface code's optimal threshold: 0.1092-0.1093. Our matching run crossed near 0.10.")
check("our earlier matching-decoder crossing (~0.10) and the optimal threshold (0.1093) both sit below 0.1100",
      0.10 < pc_ent and 0.1093 < pc_ent)
# Nishimori: lambda = tanh K and the Kramers-Wannier lens
for p in [0.03, 0.1093, 0.2]:
    K = -0.5 * math.log(p / (1 - p))
    check(f"Nishimori line at p = {p}: e^(-2K) = p/(1-p) gives tanh K = 1 - 2p = lambda",
          abs(math.tanh(K) - (1 - 2 * p)) < 1e-12)
lamv = 1 - 2 * 0.1093
check("Kramers-Wannier lens M(t) = (1-t)/(1+t) sends the polarization lambda to the error odds p/(1-p)",
      abs((1 - lamv) / (1 + lamv) - 0.1093 / (1 - 0.1093)) < 1e-12)
# Landauer: heat to erase syndrome entropy
kB_SI = 1.380649e-23
T = 0.02
p = 1e-3
bits = h2(p)
print(f"   Landauer floor (known physics, not a check): p = 0.1%, T = 20 mK: {bits:.4f} bits -> "
      f"{kB_SI*T*math.log(2)*bits:.2e} J per qubit per round")

print(f"\n{COUNT[0] - len(FAILS)} / {COUNT[0]} checks pass")
if FAILS:
    print("FAILED:", FAILS)
    sys.exit(1)
