"""Checks for the exploration: subsystem entropies live on the cube, and holographic
(gravity-dual) states are the ones whose entropy cube passes hypercube contraction tests.

Sections
  A  The k = -1 corner sum of the project's cube is the multipartite information I_n.
     It kills every one-body (additive) piece for n >= 2, so I_n is cutoff free.
  B  Contraction maps between hypercubes (Bao et al. 2015 method): brute-force / backtracking
     search for SA, SSA, MMI (exist) and reversed MMI, I3 >= 0 (none exist); the 5-party
     cyclic inequality (exists).
  C  Random weighted graphs: min-cut entropies always satisfy MMI (I3 <= 0).
  D  Random quantum states: many violate MMI; GHZ has I3 = +1 bit.
  E  AdS3 / CFT2 (Ryu-Takayanagi with geodesics in the hyperbolic disk): three random
     boundary intervals always give I3 <= 0, and I3 does not depend on the cutoff.
Writes results.json for figures.py.
"""
import itertools, json, math, random
import numpy as np

rng = np.random.default_rng(20261125)
random.seed(20261125)
OUT = {}
passes, fails = 0, 0
def check(name, ok, detail=""):
    global passes, fails
    print(("PASS " if ok else "FAIL ") + name + (("  " + detail) if detail else ""))
    if ok: passes += 1
    else: fails += 1

def subsets(parties):
    for r in range(1, len(parties) + 1):
        for T in itertools.combinations(parties, r):
            yield frozenset(T)

# ---------------------------------------------------------------- A
print("== A. k = -1 corner sum = multipartite information")
def corner_sum(S, parties, k):
    """The project's n-fold cube sum with corner weight k^(|T|-1), applied to a set function S."""
    return sum(k ** (len(T) - 1) * S[T] for T in subsets(parties))

# A1: for a product-form function S_T = prod a_i the cube sum is the circle product (sanity, the project's core identity)
a = {"A": 0.3, "B": -0.7, "C": 1.9}
for k in (-1.0, 0.5, 2.0):
    S = {T: math.prod(a[i] for i in T) for T in subsets("ABC")}
    circ = lambda x, y: x + y + k * x * y
    val = circ(circ(a["A"], a["B"]), a["C"])
    check(f"A1 cube sum of product weights = a o_k b o_k c (k={k})", abs(corner_sum(S, "ABC", k) - val) < 1e-12)

# A2: I3 written out equals the k=-1 corner sum
S = {T: rng.random() for T in subsets("ABC")}
I3 = (S[frozenset("A")] + S[frozenset("B")] + S[frozenset("C")] - S[frozenset("AB")] - S[frozenset("AC")]
      - S[frozenset("BC")] + S[frozenset("ABC")])
check("A2 I3 = k=-1 corner sum", abs(I3 - corner_sum(S, "ABC", -1)) < 1e-12)

# A3: k = -1 corner sum annihilates one-body (additive) functions for n >= 2 -> cutoff independence
for n in range(2, 7):
    parties = "ABCDEF"[:n]
    c = {p: rng.normal() for p in parties}
    S1 = {T: sum(c[p] for p in T) for T in subsets(parties)}
    check(f"A3 k=-1 corner sum kills additive pieces, n={n}", abs(corner_sum(S1, parties, -1)) < 1e-10)
# ...and for other k it does not (so k = -1 is singled out)
parties = "ABC"; c = {p: 1.0 for p in parties}
S1 = {T: sum(c[p] for p in T) for T in subsets(parties)}
vals = {k: corner_sum(S1, parties, k) for k in (-2, -1, -0.5, 0, 0.5, 1)}
check("A4 only k=-1 kills the additive piece (n=3)", [k for k, v in vals.items() if abs(v) < 1e-12] == [-1],
      str({k: round(v, 3) for k, v in vals.items()}))
OUT["A_additive_vs_k"] = {str(k): v for k, v in vals.items()}

# ---------------------------------------------------------------- B
print("== B. contraction maps between hypercubes")
def occ(terms, party):
    return tuple(1 if party in t else 0 for t in terms)

def find_contraction(lhs, rhs, parties, wl=None, wr=None, count_all=False, limit=10**7):
    """Search f: {0,1}^|lhs| -> {0,1}^|rhs| with weighted Hamming contraction and boundary conditions.
    Purifier O: all-zeros -> all-zeros.  Returns (found_map or None, number_found)."""
    L, R = len(lhs), len(rhs)
    wl = wl or [1] * L; wr = wr or [1] * R
    dL = lambda x, y: sum(w for w, a, b in zip(wl, x, y) if a != b)
    dR = lambda x, y: sum(w for w, a, b in zip(wr, x, y) if a != b)
    fixed = {tuple([0] * L): tuple([0] * R)}
    for p in parties:
        x, y = occ(lhs, p), occ(rhs, p)
        if x in fixed and fixed[x] != y:
            return None, 0
        fixed[x] = y
    # boundary points must already be mutually consistent
    fx = list(fixed.items())
    for (x1, y1), (x2, y2) in itertools.combinations(fx, 2):
        if dR(y1, y2) > dL(x1, x2):
            return None, 0
    free = [x for x in itertools.product((0, 1), repeat=L) if x not in fixed]
    # order free points: most constrained first (closest to boundary points)
    free.sort(key=lambda x: min(dL(x, b) for b in fixed))
    codomain = list(itertools.product((0, 1), repeat=R))
    assign = dict(fixed); found = [None]; count = [0]; nodes = [0]
    def rec(i):
        nodes[0] += 1
        if nodes[0] > limit: return True
        if i == len(free):
            count[0] += 1
            if found[0] is None: found[0] = dict(assign)
            return not count_all
        x = free[i]
        for y in codomain:
            if all(dR(y, yy) <= dL(x, xx) for xx, yy in assign.items()):
                assign[x] = y
                if rec(i + 1): return True
                del assign[x]
        return False
    rec(0)
    return found[0], count[0]

ineqs = {
    "SA  S_A + S_B >= S_AB": (["A", "B"], ["AB"], "AB"),
    "SSA S_AB + S_BC >= S_B + S_ABC": (["AB", "BC"], ["B", "ABC"], "ABC"),
    "MMI S_AB + S_BC + S_AC >= S_A + S_B + S_C + S_ABC": (["AB", "BC", "AC"], ["A", "B", "C", "ABC"], "ABC"),
    "reversed MMI (false)": (["A", "B", "C", "ABC"], ["AB", "BC", "AC"], "ABC"),
}
B_results = {}
for name, (lhs, rhs, parties) in ineqs.items():
    f, n = find_contraction(lhs, rhs, parties, count_all=True)
    expect = "false" not in name
    check(f"B  contraction map {'exists' if expect else 'does not exist'}: {name}", (f is not None) == expect,
          f"maps found = {n}")
    B_results[name] = {"exists": f is not None, "count": n,
                       "map": {"".join(map(str, k)): "".join(map(str, v)) for k, v in (f or {}).items()}}
# SSA is NOT a contraction-provable *holographic-specific* statement but it holds for all quantum states; fine.
# 5-party cyclic inequality (Bao et al. 2015): sum of 5 cyclic triples >= 5 cyclic pairs + ABCDE
P = "ABCDE"
lhs5 = [P[i] + P[(i + 1) % 5] + P[(i + 2) % 5] for i in range(5)]
rhs5 = [P[i] + P[(i + 1) % 5] for i in range(5)] + ["ABCDE"]
f5, n5 = find_contraction(lhs5, rhs5, P)
check("B  contraction map exists: 5-party cyclic inequality", f5 is not None)
B_results["cyclic5"] = {"lhs": lhs5, "rhs": rhs5, "exists": f5 is not None}
# a false 5-party candidate: drop the ABCDE term from the left of a reversed version
f5b, _ = find_contraction(rhs5, lhs5, P)
check("B  no contraction map: reversed 5-party cyclic (false)", f5b is None)
OUT["B"] = B_results

# ---------------------------------------------------------------- C
print("== C. random graphs: min-cut entropies")
def mincut_entropies(n_bulk, W, boundary=("A", "B", "C", "O")):
    """Vertices: boundary vertices then bulk. Entropy of boundary set T = min cut separating T from the rest."""
    nb = len(boundary); N = nb + n_bulk
    res = {}
    bulk_masks = list(itertools.product((0, 1), repeat=n_bulk))
    for T in subsets(boundary[:-1]):  # parties A,B,C (O is the purifier)
        side0 = np.array([1 if b in T else 0 for b in boundary] + [0] * n_bulk)
        best = math.inf
        for m in bulk_masks:
            side = side0.copy(); side[nb:] = m
            cut = W[np.ix_(side == 1, side == 0)].sum()
            best = min(best, cut)
        res[T] = best
    return res

I3_graphs = []
for trial in range(1500):
    nbk = int(rng.integers(1, 8))
    N = 4 + nbk
    W = rng.exponential(1.0, (N, N)) * (rng.random((N, N)) < 0.6)
    W = np.triu(W, 1); W = W + W.T
    S = mincut_entropies(nbk, W)
    I3_graphs.append(corner_sum(S, "ABC", -1))
I3_graphs = np.array(I3_graphs)
check("C  MMI holds in all 1500 random graphs (I3 <= 0)", I3_graphs.max() <= 1e-9,
      f"max I3 = {I3_graphs.max():.2e}, min = {I3_graphs.min():.2f}")
OUT["C_I3_graphs"] = I3_graphs.tolist()

# ---------------------------------------------------------------- D
print("== D. random quantum states")
def ent(psi, dims, keep):
    n = len(dims)
    t = psi.reshape(dims)
    keep = sorted(keep); rest = [i for i in range(n) if i not in keep]
    M = np.transpose(t, keep + rest).reshape(int(np.prod([dims[i] for i in keep])), -1)
    s = np.linalg.svd(M, compute_uv=False) ** 2
    s = s[s > 1e-14]
    return float(-(s * np.log2(s)).sum())

def I3_state(psi, dims):
    idx = {"A": 0, "B": 1, "C": 2}
    S = {T: ent(psi, dims, [idx[p] for p in T]) for T in subsets("ABC")}
    return corner_sum(S, "ABC", -1)

ghz = np.zeros(16); ghz[0] = ghz[15] = 1 / math.sqrt(2)
w = np.zeros(16)
for i in range(4): w[1 << i] = 0.5
I3_ghz = I3_state(ghz, [2, 2, 2, 2]); I3_w = I3_state(w, [2, 2, 2, 2])
check("D  GHZ (4 qubits) violates MMI: I3 = +1 bit", abs(I3_ghz - 1) < 1e-9, f"I3 = {I3_ghz:.4f}")
check("D  W state (4 qubits) violates MMI", I3_w > 0, f"I3 = {I3_w:.4f}")
# perfect tensor-like state: 4-qubit absolutely-maximally-entangled doesn't exist; use 4 qutrits AME(4,3)
d = 3; ame = np.zeros((3, 3, 3, 3))
for i in range(3):
    for j in range(3):
        ame[i, j, (i + j) % 3, (i + 2 * j) % 3] = 1 / 3
I3_ame = I3_state(ame.reshape(-1), [3, 3, 3, 3]) / math.log2(3)
check("D  4-qutrit perfect tensor (AME) satisfies MMI with I3 = -2 qutrits", abs(I3_ame + 2) < 1e-9, f"I3 = {I3_ame:.4f}")
OUT["D_named"] = {"GHZ": I3_ghz, "W": I3_w, "AME(4,3) [in log2 3 units]": I3_ame}

def haar(dim):
    v = rng.normal(size=dim) + 1j * rng.normal(size=dim)
    return v / np.linalg.norm(v)
I3_haar = np.array([I3_state(haar(16), [2, 2, 2, 2]) for _ in range(4000)])
# sparse random states (few terms) are far from typical and violate more often
def sparse_state(nterms):
    v = np.zeros(16, complex); idx = rng.choice(16, nterms, replace=False)
    v[idx] = rng.normal(size=nterms) + 1j * rng.normal(size=nterms)
    return v / np.linalg.norm(v)
I3_sparse = np.array([I3_state(sparse_state(int(rng.integers(2, 5))), [2, 2, 2, 2]) for _ in range(4000)])
fh, fs = (I3_haar > 1e-9).mean(), (I3_sparse > 1e-9).mean()
check("D  some random quantum states violate MMI", fh > 0 or fs > 0,
      f"Haar 4-qubit: {100*fh:.1f}% violate; sparse 2-4 term: {100*fs:.1f}% violate")
OUT["D_I3_haar"] = I3_haar.tolist(); OUT["D_I3_sparse"] = I3_sparse.tolist()

# ---------------------------------------------------------------- E
print("== E. AdS3 / CFT2: Ryu-Takayanagi with geodesics")
def geod(t1, t2, eps):
    """Regularised length of the hyperbolic-disk geodesic between boundary angles t1, t2 (cutoff eps)."""
    return 2 * math.log((2 / eps) * abs(math.sin((t2 - t1) / 2)))

def rt(endpoints, eps):
    """Min over non-crossing perfect matchings of cyclically ordered endpoints (union of intervals)."""
    pts = sorted(endpoints)
    from functools import lru_cache
    @lru_cache(None)
    def best(i, j):  # match points i..j (inclusive), j-i+1 even
        if i > j: return 0.0, ()
        out = (math.inf, ())
        for m in range(i + 1, j + 1, 2):
            a, pa = best(i + 1, m - 1); b, pb = best(m + 1, j)
            tot = geod(pts[i], pts[m], eps) + a + b
            if tot < out[0]: out = (tot, ((pts[i], pts[m]),) + pa + pb)
        return out
    return best(0, len(pts) - 1)

def region_endpoints(intervals, T):
    e = []
    for p in T: e += list(intervals[p])
    return e

def S_rt(intervals, eps, c_over_6=1.0):
    return {T: c_over_6 * rt(region_endpoints(intervals, T), eps)[0] for T in subsets("ABC")}

def random_intervals():
    cuts = np.sort(rng.random(6)) * 2 * math.pi
    # A = (c0,c1), B = (c2,c3), C = (c4,c5): disjoint, gaps between
    return {"A": (cuts[0], cuts[1]), "B": (cuts[2], cuts[3]), "C": (cuts[4], cuts[5])}

I3_ads, cut_dep = [], []
for _ in range(3000):
    iv = random_intervals()
    S1 = S_rt(iv, 1e-3); S2 = S_rt(iv, 1e-6)
    i1, i2 = corner_sum(S1, "ABC", -1), corner_sum(S2, "ABC", -1)
    I3_ads.append(i1); cut_dep.append(abs(i1 - i2))
I3_ads = np.array(I3_ads)
check("E  MMI holds for 3000 random interval triples in AdS3 (I3 <= 0)", I3_ads.max() <= 1e-9,
      f"max I3 = {I3_ads.max():.2e}")
check("E  I3 does not depend on the cutoff (eps 1e-3 vs 1e-6)", max(cut_dep) < 1e-9, f"max change {max(cut_dep):.1e}")
OUT["E_I3_ads"] = I3_ads.tolist()
# cutoff sweep for one fixed configuration (for the figure)
iv0 = {"A": (0.2, 1.4), "B": (1.6, 3.0), "C": (3.3, 5.6)}
sweep = []
for eps in np.logspace(-1, -8, 29):
    S = S_rt(iv0, eps)
    sweep.append({"eps": eps, "S_A": S[frozenset("A")], "S_ABC": S[frozenset("ABC")],
                  "I2_AB": S[frozenset("A")] + S[frozenset("B")] - S[frozenset("AB")],
                  "I3": corner_sum(S, "ABC", -1)})
OUT["E_sweep"] = sweep
OUT["E_config"] = {k: list(v) for k, v in iv0.items()}
OUT["E_config_chords"] = {"".join(sorted(T)): rt(region_endpoints(iv0, T), 1e-3)[1] for T in subsets("ABC")}

json.dump(OUT, open("results.json", "w"), default=float)
print(f"\n{passes} passed, {fails} failed")
