"""Deeper pass (Nov 25, follow-up). Sections
  F  Full-paper corrections: integer-coefficient rules, plain Hamming distance with repeated terms;
     contraction maps for the published 5-party facets (MMI uplift, cyclic) and none for the false
     "deformed MMI" of Bao-Furuya-Naskar 2025.
  G  Which basis makes the rules shortest? Coefficients in the dial-k corner basis, swept over k.
  H  OURS (reading): in AdS3, gluing adjacent intervals composes cross-ratios by the union rule o_{-1};
     holographic mutual information is (c/3) max(0, log of the flip phi_{-1}(x) = odds); the flip turns
     o_{-1} into o_1 = the height-flow law E(a+b) = E(a) o_1 E(b); and MMI for this family is the
     statement X o_1 Y >= X*Y. Checked against direct Ryu-Takayanagi, vacuum and thermal (BTZ).
  I  Dial sweeps off k = -1: corner sum sign vs k; soft-min (Maslov / random-tensor-network) dial;
     Renyi index; GHZ cube sum = n-fold o_k of 1.
  J  North Star link: entanglement first law delta S = delta <H> from RT in every dimension d (sympy),
     and its second-order failure = relative entropy >= 0 (BTZ, exact).
Writes results2.json for figures2.py.
"""
import itertools, json, math
import numpy as np
import sympy as sp
from checks_lib import subsets, corner_sum, find_contraction, ent

rng = np.random.default_rng(1125)
OUT = {}; P = [0, 0]
def check(name, ok, detail=""):
    P[0 if ok else 1] += 1
    print(("PASS " if ok else "FAIL ") + name + (("  " + detail) if detail else ""))

# ---------------------------------------------------------------- F
print("== F. published facets and the full paper's example")
cases = [
    ("orbit 3, MMI uplift: S_ABC+S_ADE+S_BCDE >= S_A+S_BC+S_DE+S_ABCDE",
     ["ABC", "ADE", "BCDE"], ["A", "BC", "DE", "ABCDE"], "ABCDE", True),
    ("orbit 4, cyclic: S_ABC+S_ABD+S_ACE+S_BDE+S_CDE >= S_AB+S_AC+S_BD+S_CE+S_DE+S_ABCDE",
     ["ABC", "ABD", "ACE", "BDE", "CDE"], ["AB", "AC", "BD", "CE", "DE", "ABCDE"], "ABCDE", True),
    ("BFN 2025 eq 3.2, deformed MMI (false): MMI with an extra S_AB on the right",
     ["AB", "AC", "BC"], ["A", "B", "C", "ABC", "AB"], "ABC", False),
    ("BFN 2025 eq 3.3 (false): S_BC+S_ACD+S_ACE+S_BCDE >= S_AC+S_BCD+S_BCE+S_ACDE",
     ["BC", "ACD", "ACE", "BCDE"], ["AC", "BCD", "BCE", "ACDE"], "ABCDE", False),
]
for name, l, r, parties, expect in cases:
    f, _ = find_contraction(l, r, parties)
    check(f"F  contraction map {'exists' if expect else 'absent'}: {name}", (f is not None) == expect)

# ---------------------------------------------------------------- G
print("== G. rule length in the dial-k corner basis")
def coeffs(lhs, rhs, parties):
    q = {T: 0.0 for T in subsets(parties)}
    for t in lhs: q[frozenset(t)] += 1
    for t in rhs: q[frozenset(t)] -= 1
    return q
def coeffs_k(q, k):
    """Rule sum_T q_T S_T rewritten in J = M_k S, (M_k S)(T) = sum_{U<=T} k^{|T-U|} S(U)."""
    return {U: sum(qT * (-k) ** (len(T) - len(U)) for T, qT in q.items() if U <= T) for U in q}
rules = {"SA": (["A", "B"], ["AB"], "AB"),
         "MMI": (["AB", "BC", "AC"], ["A", "B", "C", "ABC"], "ABC"),
         "MMI uplift": (["ABC", "ADE", "BCDE"], ["A", "BC", "DE", "ABCDE"], "ABCDE"),
         "cyclic (5)": (["ABC", "ABD", "ACE", "BDE", "CDE"], ["AB", "AC", "BD", "CE", "DE", "ABCDE"], "ABCDE")}
ks = np.round(np.linspace(-2, 2, 81), 4)
G = {}
for name, (l, r, p) in rules.items():
    q = coeffs(l, r, p)
    counts = [sum(abs(v) > 1e-9 for v in coeffs_k(q, k).values()) for k in ks]
    G[name] = counts
    kmin = ks[int(np.argmin(counts))]
    at = dict(zip(ks, counts))
    best = [float(k) for k, cnt in zip(ks, counts) if cnt == min(counts)]
    G[name + "_best"] = best
    print(f"INFO G  {name}: {at[-1.0]} terms at k=-1; k=0: {at[0.0]}; k=1: {at[1.0]}; generic: {max(counts)}; shortest at k in {best}")
check("G  SA and MMI are one term at k = -1 (shortest); the MMI uplift is shorter at k = 0 (it is I3 of merged parties)",
      G["SA"][list(ks).index(-1.0)] == 1 and G["MMI"][list(ks).index(-1.0)] == 1 and G["MMI uplift_best"] == [0.0])
OUT["G"] = {"k": ks.tolist(), "counts": G}
q = coeffs(*rules["MMI"]); c = coeffs_k(q, -1)
check("G  MMI at k=-1 is the single top corner (= -I3)", sum(abs(v) > 1e-9 for v in c.values()) == 1
      and abs(c[frozenset("ABC")] + 1) < 1e-12 or abs(c[frozenset("ABC")] - 1) < 1e-12)

# ---------------------------------------------------------------- H
print("== H. cross-ratios, the union rule, the flip and the height flow (ours, reading)")
a1, a2, b1, m, b2 = sp.symbols("a1 a2 b1 m b2", real=True)
def xr(p, q, r, s, dist=lambda u, v: v - u):
    """cross-ratio of A=(p,q), B=(r,s): x = d(p,q)d(r,s) / (d(p,r)d(q,s))"""
    return dist(p, q) * dist(r, s) / (dist(p, r) * dist(q, s))
circ = lambda x, y, k: x + y + k * x * y
x1, x2, x12 = xr(a1, a2, b1, m), xr(a1, a2, m, b2), xr(a1, a2, b1, b2)
check("H1 gluing adjacent B's: x(A, B1uB2) = x(A,B1) o_{-1} x(A,B2)   (sympy, exact)",
      sp.simplify(x12 - circ(x1, x2, -1)) == 0)
flip = lambda x, k: x / (1 + k * x)
X1, X2 = sp.symbols("X1 X2", positive=True)
xa, xb = sp.symbols("xa xb", positive=True)
check("H2 the flip phi_{-1} (odds) carries o_{-1} to o_1   (sympy)",
      sp.simplify(flip(circ(xa, xb, -1), -1) - circ(flip(xa, -1), flip(xb, -1), 1)) == 0)
t1, t2 = sp.symbols("t1 t2", real=True)
E = lambda t: sp.exp(t) - 1
check("H3 height flow: E(t1+t2) = E(t1) o_1 E(t2)   (sympy)", sp.simplify(E(t1 + t2) - circ(E(t1), E(t2), 1)) == 0)
check("H4 odds of x is E(t) with t = -log(1-x): phi_{-1}(x) = E(-log(1-x))   (sympy)",
      sp.simplify(flip(xa, -1) - (sp.exp(-sp.log(1 - xa)) - 1)) == 0)

# numeric: compare with direct RT on the circle (vacuum) and on a line at finite temperature (BTZ)
def geod_len(p, q, dist, eps=1e-4):
    return 2 * math.log(2 * dist(p, q) / eps)
def S_union(pts, dist):
    pts = sorted(pts)
    from functools import lru_cache
    @lru_cache(None)
    def best(i, j):
        if i > j: return 0.0
        return min(geod_len(pts[i], pts[mm], dist) + best(i + 1, mm - 1) + best(mm + 1, j)
                   for mm in range(i + 1, j + 1, 2))
    return best(0, len(pts) - 1) / 2  # units: c/3 per log  (S = (c/6) * length, c/6 = 1/2 * c/3)
def MI(A, B, dist):
    return S_union(list(A), dist) + S_union(list(B), dist) - S_union(list(A) + list(B), dist)
def MI_formula(A, B, dist):
    x = xr(A[0], A[1], B[0], B[1], dist)
    return max(0.0, math.log(x / (1 - x)))   # units of c/3
vac = lambda u, v: abs(math.sin((v - u) / 2))           # boundary circle, vacuum AdS3
beta = 2.0
btz = lambda u, v: abs(math.sinh(math.pi * (v - u) / beta))  # line at temperature 1/beta (BTZ, planar)
errs = {"vacuum": [], "BTZ": []}; i3s = {"vacuum": [], "BTZ": []}; slack = []
for geom, dist, span in [("vacuum", vac, 2 * math.pi), ("BTZ", btz, 3.0)]:
    for _ in range(3000):
        pts = np.sort(rng.random(5)) * span
        A, B1, B2 = (pts[0], pts[1]), (pts[2], pts[3]), (pts[3], pts[4])
        B12 = (pts[2], pts[4])
        for B in (B1, B2, B12):
            errs[geom].append(abs(MI(A, B, dist) - MI_formula(A, B, dist)))
        xa_, xb_ = xr(*A, *B1, dist), xr(*A, *B2, dist)
        errs[geom].append(abs(xr(*A, *B12, dist) - circ(xa_, xb_, -1)))
        I3 = MI(A, B1, dist) + MI(A, B2, dist) - MI(A, B12, dist)
        i3s[geom].append(I3)
        Xa, Xb = flip(xa_, -1), flip(xb_, -1)
        if Xa > 1 and Xb > 1:   # all three connected: I3 = log(Xa Xb / (Xa o_1 Xb))
            slack.append((I3, math.log(Xa * Xb / circ(Xa, Xb, 1))))
    check(f"H5 {geom}: RT mutual information = (c/3) max(0, log odds(x)) and gluing law, 3000 configs",
          max(errs[geom]) < 1e-8, f"max error {max(errs[geom]):.1e}")
    check(f"H6 {geom}: I3(A:B1:B2) <= 0 in every config", max(i3s[geom]) < 1e-9, f"max {max(i3s[geom]):.1e}")
sl = np.array(slack)
check("H7 when all three are connected, I3 = -(c/3) log(1 + 1/X1 + 1/X2): the '+X1+X2' of o_1 is MMI's slack",
      len(sl) > 50 and np.max(np.abs(sl[:, 0] - sl[:, 1])) < 1e-8, f"{len(sl)} configs")
OUT["H"] = {"I3_vac": i3s["vacuum"], "I3_btz": i3s["BTZ"]}

# ---------------------------------------------------------------- I
print("== I. dial sweeps off k = -1")
from checks_lib import mincut_entropies, random_graph, I3_state, haar, sparse_state, renyi
graphsS = []
for _ in range(600):
    nbk, W = random_graph(rng)
    graphsS.append(mincut_entropies(nbk, W))
statesH = [haar(rng, 16) for _ in range(1500)]
statesS = [sparse_state(rng) for _ in range(1500)]
def S_state(psi, alpha=1.0):
    idx = {"A": 0, "B": 1, "C": 2}
    return {T: renyi(psi, [2, 2, 2, 2], [idx[p] for p in T], alpha) for T in subsets("ABC")}
SH = [S_state(p) for p in statesH]; SS = [S_state(p) for p in statesS]
kgrid = np.round(np.linspace(-2, 0.5, 51), 3)
sweep = {"k": kgrid.tolist(), "graph_neg": [], "haar_neg": [], "sparse_neg": [], "graph_pos": []}
for k in kgrid:
    g = np.array([corner_sum(S, "ABC", k) for S in graphsS])
    h = np.array([corner_sum(S, "ABC", k) for S in SH]); s = np.array([corner_sum(S, "ABC", k) for S in SS])
    sweep["graph_neg"].append(float((g <= 1e-9).mean())); sweep["graph_pos"].append(float((g >= -1e-9).mean()))
    sweep["haar_neg"].append(float((h <= 1e-9).mean())); sweep["sparse_neg"].append(float((s <= 1e-9).mean()))
OUT["I_ksweep"] = sweep
kk = dict(zip(kgrid.tolist(), sweep["graph_neg"]))
check("I1 corner sum is <= 0 for every graph only at k = -1 on the grid (and >= 0 for k >= -1/2)",
      [k for k, v in kk.items() if v == 1.0] == [-1.0] and all(
          v == 1.0 for k, v in zip(kgrid, sweep["graph_pos"]) if k >= -0.5),
      f"k with all <= 0: {[k for k, v in kk.items() if v == 1.0]}")

# soft-min dial: S_h(T) = -h log(Z_T / Z_empty), Z = sum over bulk configs exp(-cut/h). h -> 0 is min cut.
from checks_lib import softmin_entropies
hs = [0.0, 0.05, 0.1, 0.2, 0.35, 0.5, 0.75, 1.0, 1.5, 2.0, 3.0]
soft = {"h": hs, "mmi_viol": [], "sa_viol": []}
graphs_soft = [random_graph(rng) for _ in range(400)]
for h in hs:
    vm = vs = 0
    for nbk, W in graphs_soft:
        S = softmin_entropies(nbk, W, h) if h > 0 else mincut_entropies(nbk, W)
        vm += corner_sum(S, "ABC", -1) > 1e-9
        vs += any(S[frozenset(a)] + S[frozenset(b)] < S[frozenset(a + b)] - 1e-9 for a, b in ["AB", "AC", "BC"])
    soft["mmi_viol"].append(vm / len(graphs_soft)); soft["sa_viol"].append(vs / len(graphs_soft))
OUT["I_soft"] = soft
check("I2 soft-min dial (positive weights): no SA or MMI violation at any h tested, h = 0 to 3",
      max(soft["mmi_viol"]) == 0 and max(soft["sa_viol"]) == 0)
viol_signed = 0
for _ in range(300):
    nbk, W = random_graph(rng)
    W = W * rng.choice([1, -1], size=W.shape); W = np.triu(W, 1); W = W + W.T
    viol_signed += corner_sum(softmin_entropies(nbk, W, 1.0), "ABC", -1) > 1e-9
OUT["I_soft_signed_viol"] = viol_signed / 300
check("I2b ...but with mixed-sign weights (not a ferromagnet) MMI fails often", viol_signed > 0,
      f"{100*viol_signed/300:.0f}% of 300 signed graphs at h = 1")
# five-party facets at finite h (6 boundary vertices incl. purifier)
bd6 = ("A", "B", "C", "D", "E", "O")
def F6(nbk, W, T, h):
    nb = len(bd6); side0 = np.array([1 if b in T else 0 for b in bd6] + [0] * nbk); cs = []
    for mbits in itertools.product((0, 1), repeat=nbk):
        sd = side0.copy(); sd[nb:] = mbits; cs.append(W[np.ix_(sd == 1, sd == 0)].sum())
    c = np.array(cs); c0 = c.min()
    return c0 if h == 0 else c0 - h * math.log(np.exp(-(c - c0) / h).sum())
worst5 = -1e9
for _ in range(120):
    nbk = int(rng.integers(1, 6)); N = 6 + nbk
    W = rng.exponential(1, (N, N)) * (rng.random((N, N)) < 0.6); W = np.triu(W, 1); W = W + W.T
    for h in (0.3, 1.0, 5.0):
        for l, r in (rules["MMI uplift"][:2], rules["cyclic (5)"][:2]):
            Sv = {t: F6(nbk, W, frozenset(t), h) for t in set(l + r)}
            worst5 = max(worst5, sum(Sv[t] for t in r) - sum(Sv[t] for t in l))
check("I2c soft-min dial: the five-party facets (uplift, cyclic) also hold at h = 0.3, 1, 5 (120 graphs)",
      worst5 < 1e-9, f"worst RHS-LHS {worst5:.1e}")

alphas = [0.25, 0.5, 0.75, 1.0, 1.5, 2.0, 3.0, 5.0, 1e6]
ren = {"alpha": alphas, "haar_viol": [], "sparse_viol": []}
for al in alphas:
    ren["haar_viol"].append(float(np.mean([corner_sum(S_state(p, al), "ABC", -1) > 1e-9 for p in statesH[:600]])))
    ren["sparse_viol"].append(float(np.mean([corner_sum(S_state(p, al), "ABC", -1) > 1e-9 for p in statesS[:600]])))
OUT["I_renyi"] = ren
check("I3 Renyi index: random states break MMI at every alpha tested (no Renyi index rescues them)",
      min(ren["sparse_viol"]) > 0, str(dict(zip(alphas, [round(v, 3) for v in ren["sparse_viol"]]))))

for n in range(2, 8):
    parties = "ABCDEFG"[:n]
    S1 = {T: 1.0 for T in subsets(parties)}   # GHZ on n parties + purifier: every entropy = 1 bit
    for k in (-1.0, -0.5, 0.7, 2.0):
        val = 1.0
        for _ in range(n - 1): val = circ(val, 1.0, k)
        ok = abs(corner_sum(S1, parties, k) - val) < 1e-9
        if not ok: break
    check(f"I4 GHZ, n={n}: cube sum = 1 o_k 1 o_k ... o_k 1 = ((1+k)^n - 1)/k, equal to 1 at k=-1", ok and
          abs(corner_sum(S1, parties, -1) - 1) < 1e-12)

# ---------------------------------------------------------------- J
print("== J. first law of entanglement from RT (North Star link)")
d, R, r, mm, Gn = sp.symbols("d R r m G", positive=True)
for dv in range(2, 7):
    z = sp.sqrt(R**2 - r**2)
    # area density of the RT hemisphere in ds^2 = (1/z^2)(-f dt^2 + dx^2 + dz^2/f), f = 1 - m z^d, static
    zp = sp.diff(z, r)
    integrand = r**(dv - 2) / z**(dv - 1) * sp.sqrt(1 + zp**2 / (1 - mm * z**dv))
    dA = sp.simplify(sp.diff(integrand, mm).subs(mm, 0))
    dS = sp.integrate(dA, (r, 0, R)) / (4 * Gn)              # per unit angular volume Omega_{d-2}
    eps_E = (dv - 1) * mm / (16 * sp.pi * Gn)                  # energy density of the black brane
    dH = 2 * sp.pi * eps_E * sp.integrate((R**2 - r**2) / (2 * R) * r**(dv - 2), (r, 0, R))
    check(f"J1 d={dv}: delta S (RT area change) = delta <H_modular> (ball, first order)",
          sp.simplify(mm * dS - dH) == 0, f"both = {sp.simplify(mm*dS)} per unit solid angle")
# BTZ second order: Delta H - Delta S = relative entropy >= 0
y = np.linspace(0.01, 4, 200)   # y = pi * l / beta
rel = (y**2 / 6 - np.log(np.sinh(y) / y))  # in units of c/3
check("J2 BTZ: Delta<H> - Delta S = (c/3)[y^2/6 - log(sinh y / y)] >= 0 for all sizes (relative entropy)",
      bool((rel >= -1e-12).all()) and abs(rel[0]) < 1e-6)
OUT["J"] = {"y": y.tolist(), "dH": (y**2 / 6).tolist(), "dS": np.log(np.sinh(y) / y).tolist()}

json.dump(OUT, open("results2.json", "w"), default=float)
print(f"\n{P[0]} passed, {P[1]} failed")
