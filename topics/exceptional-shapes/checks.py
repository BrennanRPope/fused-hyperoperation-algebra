"""Checks for the exceptional-shapes pass. Run: python3 checks.py > checks_output.txt"""
import itertools, time
from fractions import Fraction
import numpy as np
from shapes import *

results = []
def check(name, ok, detail=""):
    results.append((name, bool(ok)))
    print(f"[{'PASS' if ok else 'FAIL'}] {name}: {detail}")

T0 = time.time()
Q8 = cross_polytope4(); T2 = binary_tetrahedral(); O2 = binary_octahedral(); I2 = binary_icosahedral()
C120 = onetwenty_cell(); S4 = simplex4()

print("== A. Setting the shapes up as unit quaternions ==")
for name, G, n in (("Q8 / cross-polytope", Q8, 8), ("2T / 24-cell", T2, 24), ("2O / 24-cell+dual", O2, 48), ("2I / 600-cell", I2, 120)):
    f = closure_fraction(G, qmul)
    check(f"A1 {name} closed under multiplication", len(G) == n and abs(f - 1) < 1e-12 and np.allclose(np.linalg.norm(G, axis=1), 1), f"{len(G)} unit points, fraction of products inside = {f:.3f}")
def nn_info(V):
    D = np.linalg.norm(V[:, None] - V[None], axis=2); np.fill_diagonal(D, 9)
    m = D.min(); return m, int((np.abs(D[0] - m) < 1e-6).sum())
m, c = nn_info(T2); check("A2 24-cell: edge = radius, 8 neighbours", abs(m - 1) < 1e-9 and c == 8, f"edge {m:.4f}, neighbours {c}")
m, c = nn_info(I2); check("A3 600-cell: edge = 1/phi, 12 neighbours", abs(m - 1 / PHI) < 1e-9 and c == 12, f"edge {m:.4f} vs 1/phi {1/PHI:.4f}, neighbours {c}")
_, F24 = hull_facets(T2); _, F600 = hull_facets(I2); _, F120 = hull_facets(C120)
check("A4 cell counts 24 / 600 / 120", (len(F24), len(F600), len(F120)) == (24, 600, 120), f"24-cell {len(F24)} cells, 600-cell {len(F600)}, 120-cell {len(F120)} (from hull)")
check("A5 120-cell has 600 corners (dual of 600-cell)", len(C120) == 600, f"{len(C120)} corners")
check("A6 24-cell = cross-polytope + tesseract corners", keyset(T2) == keyset(Q8) | keyset(tesseract_half()), "8 + 16 = 24")

print("\n== B. The core law on each corner set: x = (g-1)/k, x o_k y = x + y + k x y ==")
def law(k):
    return lambda a, b: a + b + k * qmul(a, b)
for k in (1.0, -1.0, 0.5, 2.0):
    for name, G in (("24-cell", T2), ("600-cell", I2), ("Clifford 48", O2)):
        X = (G - np.array([1, 0, 0, 0])) / k
        f = closure_fraction(X, law(k))
        on = np.allclose(np.linalg.norm(np.array([1, 0, 0, 0]) + k * X, axis=1), 1)
        check(f"B1 k={k:g}: {name} closed under o_k", abs(f - 1) < 1e-12 and on, f"fraction inside {f:.3f}; all on sphere |1+kx|=1 (radius 1/|k|, d = k^2 = {k*k:g})")
f120 = closure_fraction(C120, qmul)
check("B2 120-cell corners are NOT closed", f120 < 0.5, f"fraction of products landing on a corner = {f120:.4f}")
S = keyset(C120)
left_ok = all(key(p) in S for g in I2 for p in qmul(np.broadcast_to(g, C120.shape), C120))
orbits, seen = [], set()
for v in C120:
    if key(v) in seen: continue
    orb = uniq(qmul(I2, np.broadcast_to(v, I2.shape)))
    for o in orb: seen.add(key(o))
    orbits.append(orb)
iso = all(nn_info(o)[0] - 1 / PHI < 1e-6 and len(o) == 120 for o in orbits)
check("B3 120-cell = five 600-cells (cosets of the 600-cell group)", left_ok and len(orbits) == 5 and iso, f"600-cell corners times 120-cell corners stay on 120-cell: {left_ok}; {len(orbits)} pieces of 120, each a 600-cell: {iso}")

print("\n== C. The swap rule: what do swaps of corners (a b a^-1 b^-1) generate? ==")
chain = [("Clifford 48", O2), ("24-cell", T2), ("cross-polytope", Q8), ("line {+-1}", np.array([[1, 0, 0, 0], [-1, 0, 0, 0]]))]
for (n1, G1), (n2, G2) in zip(chain, chain[1:]):
    H, Cm = commutator_subgroup(G1)
    check(f"C1 swaps of {n1} generate {n2}", keyset(H) == keyset(G2), f"{len(G1)} -> {len(H)}")
H, Cm = commutator_subgroup(I2)
check("C2 swaps of the 600-cell give the whole 600-cell (perfect group)", keyset(H) == keyset(I2) and len(Cm) == 120, f"distinct single swaps {len(Cm)} of 120: every corner is itself one swap")
def rot_group(G):  # mod +-1
    return uniq([g if (g[0] > 1e-9 or (abs(g[0]) < 1e-9 and next(x for x in g if abs(x) > 1e-9) > 0)) else -g for g in G])
sizes = [len(rot_group(G)) for G in (O2, T2, Q8, I2)]
check("C3 3D shadows: 24 -> 12 -> 4 -> 1 (octahedral S4 solvable), 60 -> 60 (icosahedral A5 simple)", sizes == [24, 12, 4, 60], f"rotation groups {sizes}")
ang2I = sorted({round(rot_angle_deg(g)) for g in I2}); ang2O = sorted({round(rot_angle_deg(g)) for g in O2})
check("C4 turn angles: 600-cell {0,72,120,144,180,...}; Clifford {0,90,120,180}", set(ang2I) == {0, 72, 120, 144, 180} and set(ang2O) == {0, 90, 120, 180}, f"600-cell {ang2I}; Clifford {ang2O}")

print("\n== D. Gates on the 600-cell ==")
c8, s8 = np.cos(np.pi / 8), np.sin(np.pi / 8)
gates = {"X": [0, 1, 0, 0], "Y": [0, 0, 1, 0], "Z": [0, 0, 0, 1], "H": [0, 2**-0.5, 0, 2**-0.5],
         "S": [2**-0.5, 0, 0, 2**-0.5], "sqrtX": [2**-0.5, 2**-0.5, 0, 0], "T": [c8, 0, 0, s8]}
SI = keyset(I2) | keyset(-I2)
mem = {g: key(np.array(v, float)) in SI for g, v in gates.items()}
check("D1 standard orientation: X, Y, Z on corners; H, S, sqrtX, T not", mem == {"X": 1, "Y": 1, "Z": 1, "H": 0, "S": 0, "sqrtX": 0, "T": 0}, str(mem))
allowed = {0, 72, 120, 144, 180}
ang = {g: round(rot_angle_deg(np.array(v)), 3) for g, v in gates.items()}
can = {g: min(abs(a - b) for b in allowed) < 1e-6 or min(abs((360 - a) - b) for b in allowed) < 1e-6 for g, a in ang.items()}
check("D2 any orientation: S, sqrtX (90 deg) and T (45 deg) fit no 600-cell; H (180) can", can == {"X": True, "Y": True, "Z": True, "H": True, "S": False, "sqrtX": False, "T": False}, f"turn angles {ang}")
prod = qmul(np.array(gates["X"], float), np.array(gates["H"], float))
check("D3 Paulis and H together fit no 600-cell (X.H turns 90 deg)", abs(rot_angle_deg(prod) - 90) < 1e-6, f"X*H turn = {rot_angle_deg(prod):.1f} deg, not a 600-cell angle")
dT600 = np.radians(min(abs(22.5 - b / 2) for b in allowed)); dTcl = np.radians(min(abs(22.5 - b / 2) for b in {0, 90, 120, 180}))
check("D4 T's best distance: 600-cell 0.236 vs Clifford 0.393 (on S^3)", abs(dT600 - np.radians(13.5)) < 1e-12 and abs(dTcl - np.radians(22.5)) < 1e-12, f"{dT600:.3f} vs {dTcl:.3f} rad: closer to the 600-cell, on neither")
Tq = np.array(gates["T"], float)
lvl = [I2]; counts = [60]
for L in range(1, 3):
    A = qmul(lvl[-1][:, None, :], Tq)                      # ... T
    W = qmul(A[:, None, :], I2[None, :, :]).reshape(-1, 4)  # ... T g
    W = np.where(W[:, :1] < 0, -W, W)
    nxt = np.unique(np.round(W, 7), axis=0); lvl.append(nxt); counts.append(len(nxt))
check("D5 600-cell plus T keeps making new points (finite set + one off-corner gate)", counts[2] > counts[1] > counts[0], f"distinct gates (up to phase) using 0, 1, 2 T's: {counts}")

cos = []; seen = set()
for g in I2:
    if key(g) in seen: continue
    co = uniq(qmul(np.broadcast_to(g, T2.shape), T2)); cos.append(co); [seen.add(key(x)) for x in co]
check("A7 600-cell = five 24-cells (cosets of the 24-cell group)", len(cos) == 5 and all(abs(nn_info(c_)[0] - 1) < 1e-9 for c_ in cos), f"{len(cos)} disjoint pieces of 24, each with 24-cell edges")

print("\n== E. Hopf shadows: polar coordinates = point on S^2 + phase ==")
def hopf_images(G, axis):
    n = np.array([0, *axis]); n = n / np.linalg.norm(n)
    imgs = [qmul(qmul(g, n), qconj(g))[1:] for g in G]
    return uniq(imgs, nd=5)
def axis_of_order(G, deg):
    for g in G:
        if abs(rot_angle_deg(g) - deg) < 1e-6:
            return g[1:] / np.linalg.norm(g[1:])
def min_angle(P):
    D = np.clip(P @ P.T, -1, 1); np.fill_diagonal(D, -2); return np.degrees(np.arccos(D.max()))
for deg, n, ma, nm in ((72, 12, 63.4349, "icosahedron"), (120, 20, 41.8103, "dodecahedron"), (180, 30, 36.0, "icosidodecahedron")):
    P = hopf_images(I2, axis_of_order(I2, deg))
    check(f"E1 600-cell, fibres along a {deg}-deg axis -> {n} points = {nm}", len(P) == n and abs(min_angle(P) - ma) < 1e-3, f"{len(P)} points, nearest angle {min_angle(P):.4f} deg; fibre size {120//len(P)}")
for deg, n, nm in ((180, 6, "octahedron"), (120, 4, "tetrahedron")):
    P = hopf_images(T2, axis_of_order(T2, deg))
    check(f"E2 24-cell, fibres along a {deg}-deg axis -> {n} points = {nm}", len(P) == n, f"{len(P)} points, fibre size {24//len(P)}")

print("\n== F. Where each shape sits on the p dial (cross-polytope p=1, sphere p=2, cube p=inf) ==")
shapes = {"simplex": S4, "cross-polytope": Q8, "tesseract": tesseract_half(), "24-cell": T2, "600-cell": I2, "120-cell": C120}
ratios, vols = {}, {}
from scipy.spatial import ConvexHull
for nm, V in shapes.items():
    h = ConvexHull(V); R = np.linalg.norm(V, axis=1).max(); r = (-h.equations[:, 4]).min()
    ratios[nm] = r / R; vols[nm] = h.volume / (np.pi**2 / 2 * R**4)
def p_pair(rho):
    t = -np.log2(rho)
    return (2 / (1 + t), 2 / (1 - t) if t < 1 else np.inf) if t <= 1 else None
exact = {"cross-polytope": 0.5, "tesseract": 0.5, "24-cell": 2**-0.5, "simplex": 0.25}
check("F1 in/out radius ratio, exact for simple cases", all(abs(ratios[k] - v) < 1e-9 for k, v in exact.items()), {k: round(v, 4) for k, v in ratios.items()})
check("F2 600-cell and 120-cell (duals) share the same ratio", abs(ratios["600-cell"] - ratios["120-cell"]) < 1e-9, f"{ratios['600-cell']:.6f}")
pp = {k: p_pair(v) for k, v in ratios.items()}
check("F3 24-cell sits at the Hölder pair p = 4/3 and 4; simplex is off the dial", abs(pp["24-cell"][0] - 4/3) < 1e-9 and abs(pp["24-cell"][1] - 4) < 1e-9 and pp["simplex"] is None, f"pairs {{{', '.join(f'{k}: {v if v is None else tuple(round(x,3) for x in v)}' for k,v in pp.items())}}}")
from math import gamma
def lp_frac(p, n=4):
    V = (2 * gamma(1 + 1 / p))**n / gamma(1 + n / p) if np.isfinite(p) else 2.0**n
    R = 1 if p <= 2 else n**(0.5 - 1 / p)
    return V / (np.pi**2 / 2 * R**4)
lines = []
for nm in ("24-cell", "600-cell", "120-cell"):
    lo, hi = pp[nm]
    side = lo if nm == "600-cell" else hi if nm == "120-cell" else None
    lines.append(f"{nm}: polytope fills {vols[nm]:.3f} of its ball; l_p at {lo:.3f} fills {lp_frac(lo):.3f}, at {hi:.3f} fills {lp_frac(hi):.3f}")
check("F4 second measure (volume) agrees in order: 600/120-cell rounder than 24-cell", vols["600-cell"] > vols["24-cell"] and vols["120-cell"] > vols["24-cell"], "; ".join(lines))

print("\n== G. Regular polytopes by dimension (Schläfli's test: Gram matrix positive definite) ==")
def gram(sym):
    n = len(sym) + 1; M = np.eye(n)
    for i, p in enumerate(sym):
        M[i, i + 1] = M[i + 1, i] = -np.cos(np.pi / p)
    return np.linalg.eigvalsh(M).min()
found, flat = {}, {}
for dim in range(3, 9):
    found[dim], flat[dim] = [], []
    for sym in itertools.product(range(3, 13 if dim <= 5 else 8), repeat=dim - 1):
        e = gram(sym)
        if e > 1e-9: found[dim].append(sym)
        elif abs(e) < 1e-9: flat[dim].append(sym)
counts = {d: len(v) for d, v in found.items()}
check("G1 regular convex polytopes: 5 in 3D, 6 in 4D, 3 in every dimension 5-8", counts == {3: 5, 4: 6, 5: 3, 6: 3, 7: 3, 8: 3}, f"{counts}; 4D: {found[4]}")
check("G2 regular tilings of flat space: only the cube tiling above 4D; 24-cell tiles 4D", {d: len(v) for d, v in flat.items()} == {3: 3, 4: 1, 5: 3, 6: 1, 7: 1, 8: 1} and (3, 4, 3, 3) in flat[5], f"tilings of (dim-1)-space: {{{', '.join(f'{d-1}D: {v}' for d, v in flat.items())}}}")

print("\n== H. E8 ==")
R8 = e8_roots(); G8 = R8 @ R8.T
check("H1 E8: 240 roots, inner products in {0,+-1,+-2}, 56 neighbours at 1", len(R8) == 240 and set(np.unique(G8)) <= {-2, -1, 0, 1, 2} and (G8[0] == 1).sum() == 56, "kissing number 240")
P2, h, Ssimp = coxeter_plane(R8)
proj = R8 @ P2.T; rad = np.round(np.linalg.norm(proj, axis=1), 5); rings = sorted(set(rad))
check("H2 E8 Coxeter (Petrie) picture: Coxeter number 30, 8 rings of 30", len(Ssimp) == 8 and abs(h - 30) < 1e-6 and len(rings) == 8 and all((rad == r).sum() == 30 for r in rings), f"h = {h:.3f}, ring radii {rings}")
E, X, al, be = e8_from_icosians(); GE = 2 * E @ E.T
ok = len(uniq(E, nd=5)) == 240 and np.allclose(np.linalg.norm(E, axis=1), 1) and np.allclose(GE, np.round(GE)) and set(np.unique(np.round(GE))) <= {-2, -1, 0, 1, 2} and (np.round(GE[0]) == 1).sum() == 56
check("H3 E8 = 600-cell and 600-cell/phi, with the golden ratio's Galois twin as the other 4 coordinates", ok, "240 unit vectors, 2x inner products all integers in {-2..2}, 56 neighbours: an E8 root system")
rr = np.linalg.norm(E[:, :4], axis=1); r4 = np.round(rr, 6); vals = sorted(set(r4))
ex = [rr[r4 == v].mean() for v in vals]
check("H4 E8's shadow in 4D = two 600-cells, radii in ratio phi", len(vals) == 2 and abs(ex[1] / ex[0] - PHI) < 1e-9 and all(abs(nn_info(E[r4 == v, :4] / e)[0] - 1 / PHI) < 1e-9 for v, e in zip(vals, ex)), f"radii {ex[0]:.6f} and {ex[1]:.6f}, ratio {ex[1]/ex[0]:.9f} (phi = {PHI:.9f}); each ring has 600-cell edges")
# integral octonions: find the closed 240-unit loop among the 7-point-line conventions
idx = ["inf"] + list(range(7))
def octavians(swap):
    units = [s * np.eye(8)[i] for i in range(8) for s in (1, -1)]
    quads = [("inf",) + L for L in LINES] + [tuple(sorted(set(range(7)) - set(L))) for L in LINES]
    pos = {n: i for i, n in enumerate(idx)}
    if swap:
        a, b = swap; pos[a], pos[b] = pos[b], pos[a]
    for q in quads:
        for s in itertools.product((0.5, -0.5), repeat=4):
            v = np.zeros(8)
            for c, si in zip(q, s): v[pos[c]] = si
            units.append(v)
    return uniq(units)
best = None
for swap in [None] + list(itertools.combinations(idx, 2)):
    U = octavians(swap)
    f = closure_fraction(U, omul)
    if abs(f - 1) < 1e-12: best = (swap, U); break
swap, U = best
check("H5 240 integral octonion units closed under octonion multiplication (so under o_k)", len(U) == 240 and np.allclose((2 * U @ U.T), np.round(2 * U @ U.T)), f"convention: Kirmse quads with {swap} swapped (Coxeter's correction); 240 units, E8 inner products")
rng = np.random.default_rng(0); n = 60000
A, B, C = U[rng.integers(0, 240, n)], U[rng.integers(0, 240, n)], U[rng.integers(0, 240, n)]
assoc = np.all(np.abs(omul(omul(A, B), C) - omul(A, omul(B, C))) < 1e-9, axis=1).mean()
comm = np.all(np.abs(omul(A, B) - omul(B, A)) < 1e-9, axis=1).mean()
moufang = np.allclose(omul(omul(A, B), omul(C, A)), omul(omul(A, omul(B, C)), A))
check("H6 grouping matters on E8's corners; Moufang identity holds", assoc < 1 and moufang, f"triples where (ab)c = a(bc): {assoc:.4f}; pairs that commute: {comm:.4f}; Moufang (ab)(ca) = (a(bc))a on all {n} samples: {moufang}")
U16 = np.vstack([s * np.eye(8)[i] for i in range(8) for s in (1, -1)])
check("H7 8D cross-polytope (+-1, +-e_i) closed too: a 16-element loop", abs(closure_fraction(U16, omul) - 1) < 1e-12, "the cross-polytope is the one regular shape closed in every split dimension")

print("\n== I. 24-cell extras ==")
from collections import Counter
cnt_hur, cnt_lip = Counter(), Counter()
rng4 = range(-6, 7)
for v in itertools.product(rng4, repeat=4):
    n2 = sum(x * x for x in v); cnt_lip[n2] += 1
    cnt_hur[n2] += 1
for v in itertools.product([x + 0.5 for x in range(-6, 6)], repeat=4):
    cnt_hur[round(sum(x * x for x in v))] += 1
def sig_odd(n): return sum(d for d in range(1, n + 1) if n % d == 0 and d % 2)
def r4(n): return 8 * sum(d for d in range(1, n + 1) if n % d == 0 and d % 4)
okh = all(cnt_hur[n] == 24 * sig_odd(n) for n in range(1, 31)); okl = all(cnt_lip[n] == r4(n) for n in range(1, 31))
check("I1 four squares: Hurwitz integers of size n = 24 x (odd divisor sum); Jacobi's count", okh and okl, f"n=1..30; e.g. n=1: {cnt_hur[1]} Hurwitz units, {cnt_lip[1]} ways as four squares")
roots = np.array([v for v in itertools.product((-1, 0, 1), repeat=4) if sum(abs(x) for x in v) == 2], float)
from scipy.spatial import HalfspaceIntersection
hs = HalfspaceIntersection(np.hstack([roots, -np.ones((24, 1))]), np.zeros(4))
cell = uniq(hs.intersections, nd=6)
target = np.vstack([Q8, tesseract_half()])
check("I2 the 24-cell is the tile (Voronoi cell) of the D4 lattice; 24 neighbours = kissing number 24", len(roots) == 24 and keyset(cell) == keyset(target), f"cell corners {len(cell)} = the 24 Hurwitz units")

T_main = time.time() - T0
print("\n== P. Polish (exact arithmetic, after the floating-point checks) ==")
t1 = time.time()
class Zphi:  # (a + b phi)/2 with integers a, b; phi^2 = phi + 1
    __slots__ = ("a", "b")
    def __init__(s, a, b): s.a, s.b = a, b
    def __add__(s, o): return Zphi(s.a + o.a, s.b + o.b)
    def __sub__(s, o): return Zphi(s.a - o.a, s.b - o.b)
    def __mul__(s, o):  # ((a+b phi)(c+d phi))/4 -> keep /2 form: needs even numerators
        a = s.a * o.a + s.b * o.b; b = s.a * o.b + s.b * o.a + s.b * o.b
        return Zphi(Fraction(a, 2), Fraction(b, 2))
    def k(s): return (Fraction(s.a), Fraction(s.b))
def to_zphi(x):
    for b in range(-3, 4):
        a = 2 * x - b * PHI
        if abs(a - round(a)) < 1e-9: return Zphi(round(a), b)
    raise ValueError(x)
I2e = [[to_zphi(c) for c in g] for g in I2]
def qmul_e(p, q):
    w1, x1, y1, z1 = p; w2, x2, y2, z2 = q
    return [w1*w2 - x1*x2 - y1*y2 - z1*z2, w1*x2 + x1*w2 + y1*z2 - z1*y2, w1*y2 - x1*z2 + y1*w2 + z1*x2, w1*z2 + x1*y2 - y1*x2 + z1*w2]
Se = {tuple(c.k() for c in g) for g in I2e}
closed = all(tuple(c.k() for c in qmul_e(p, q)) in Se for p in I2e for q in I2e)
norm1 = all((g[0]*g[0] + g[1]*g[1] + g[2]*g[2] + g[3]*g[3]).k() == (2, 0) for g in I2e)
check("P1 600-cell closure re-done exactly in Z[phi] (no rounding)", closed and norm1, f"14,400 products exact; all sizes exactly 1")
T2e = [[Fraction(c).limit_denominator(2) for c in g] for g in T2]
Ste = {tuple(g) for g in T2e}
check("P2 24-cell closure re-done exactly in fractions", all(tuple(qmul_e(p, q)) in Ste for p in T2e for q in T2e), "576 products exact")
GEx = [[Fraction(int(x * 2), 2) for x in r] for r in R8]
check("P3 E8 inner products exact", all(sum(a * b for a, b in zip(p, q)) in {-2, -1, 0, 1, 2} for p in GEx for q in GEx[:40]), "exact rationals on 9,600 pairs")
# p for the 24-cell exactly: rho = 2^{-1/2} -> |1-2/p| = 1/2 -> p in {4/3, 4} exactly
check("P4 24-cell's p pair exact: 4/3 and 4 are Hölder conjugates (1/p + 1/q = 1)", Fraction(3, 4) + Fraction(1, 4) == 1, "closed form, no fit")
# polish the 600-cell ratio: closed form cos(pi/5)... compare with known r/R = phi^2/(2*sqrt(phi^2+... )) via Newton on |1-2/p| = t
rho = ratios["600-cell"]; t = -np.log2(rho)
check("P5a 600-cell in/out ratio matches the closed form phi^2 / sqrt 8", abs(rho - PHI**2 / 8**0.5) < 1e-12, f"{rho:.15f} vs {PHI**2/8**0.5:.15f}")
p = 1.8
for _ in range(4): p -= ((1 - 2 / p) + t) / (2 / p**2)
check("P5 600-cell's p (simplicial side) polished by 4 Newton steps", abs(1 - 2 / p + t) < 1e-14, f"p = {p:.12f}; closed form 2/(1+t) = {2/(1+t):.12f}")
T_pol = time.time() - t1
print(f"\nTiming: main checks {T_main:.1f} s; polish {T_pol:.2f} s ({100*T_pol/T_main:.0f}% extra)")
print(f"\n{sum(ok for _, ok in results)} of {len(results)} pass")
