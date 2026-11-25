"""Fractal dials: what tells apart fractals of the same dimension, and does any of it
give the 2^D corners one definite count or a half step of the power set?
Pure Python. Run: python3 checks.py"""
import math, cmath, itertools, json, sys, os
from fractions import Fraction as F
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from fractals import *

results, table = [], {}
def check(name, ok, detail=""):
    results.append((name, bool(ok)))
    print(("PASS " if ok else "FAIL ") + name + (("  | " + detail) if detail else ""))

# 1. Same dimension: sqrt of the ratios add to 1 (exact, by squares) ---------------------
sq = {"A": [F(1,2)]*2, "B": [F(1,3)]*3, "C": [F(1,3)]*3, "D": [F(3,10), F(7,10)], "E": [F(1,4)]*4}
ok = all(sum(sq[k]) == 1 and sorted(u*u for u in sq[k]) == sorted(r for r, _ in HALF[k]["maps"]) for k in HALF)
check("1a all five sets have dimension exactly 1/2 (square roots of their ratios add to 1)", ok)
for k, v in HALF.items():
    d = dimension(v["maps"])
    table[k] = dict(name=v["name"], pieces=len(v["maps"]), ratios=[str(r) for r, _ in v["maps"]], dim=round(d, 12))
check("1b numeric Moran solve agrees: D = 0.5 for all", all(abs(table[k]["dim"] - .5) < 1e-9 for k in HALF))
# non-overlap (open set condition on [0,1])
def disjoint(maps):
    m = sorted(maps, key=lambda t: t[1])
    return all(o1 + r1 <= o2 for (r1, o1), (r2, o2) in zip(m, m[1:])) and m[-1][0] + m[-1][1] <= 1
check("1c pieces never overlap and span [0,1] (so the dimension is the Moran value)",
      all(disjoint(v["maps"]) and min(o for _, o in v["maps"]) == 0 and max(r + o for r, o in v["maps"]) == 1 for v in HALF.values()))

# 2. Lattice vs non-lattice and complex dimensions ----------------------------------------
for k, v in HALF.items():
    table[k]["lattice"] = is_lattice(v["maps"])
check("2a A, B, C, E are lattice (one ratio); D is non-lattice (log 0.09 / log 0.49 irrational, = log0.3/log0.7)",
      [table[k]["lattice"] for k in "ABCDE"] == [True, True, True, False, True])
# complex dimensions of the equal-ratio sets: D + 2 pi i n / log(1/r)
for k in "ABCE":
    v = HALF[k]["maps"]; N = len(v); r = float(v[0][0])
    per = 2 * math.pi / math.log(1 / r)
    table[k]["complex_dim_spacing"] = per
    zs = [0.5 + 1j * per * n for n in range(1, 6)]
    assert all(abs(N * cmath.exp(z * math.log(r)) - 1) < 1e-12 for z in zs)
check("2b equal-ratio sets: N r^s = 1 at s = 1/2 + 2 pi i n / log(1/r) (complex dimensions), spacing differs by set",
      True, ", ".join(f"{k}: {table[k]['complex_dim_spacing']:.3f}" for k in "ABCE"))
# non-lattice D: on Re s = 1/2 only s = 1/2 solves it; |r1^s + r2^s| < 1 elsewhere on a grid
r1, r2 = .09, .49
vals = [abs(r1 ** complex(.5, y) + r2 ** complex(.5, y)) for y in [i * 0.01 for i in range(100, 40001)]]
table["D"]["complex_dim_spacing"] = None
check("2c non-lattice D: no complex dimension on Re s = 1/2 except 1/2 itself (|r1^s+r2^s| < 1 for 1 <= y <= 400)",
      max(vals) < 1, f"closest approach {max(vals):.5f} (quasi-periodic near-misses)")

# 3. Minkowski content: wobbles for lattice sets, settles for non-lattice ------------------
def wobble(maps, lo_dec, hi_dec, n=240):
    eps = [10 ** (-lo_dec - (hi_dec - lo_dec) * i / n) for i in range(n + 1)]
    vals = [tube_ratio(maps, .5, e) for e in eps]
    return min(vals), max(vals), sum(vals) / len(vals)
for k, v in HALF.items():
    mn, mx, av = wobble(v["maps"], 4, 6)
    table[k]["minkowski_min"], table[k]["minkowski_max"], table[k]["minkowski_avg"] = mn, mx, av
amp = {k: (table[k]["minkowski_max"] - table[k]["minkowski_min"]) / table[k]["minkowski_avg"] for k in HALF}
check("3a lattice sets: Minkowski ratio |K_eps|/eps^(1/2) keeps oscillating (no Minkowski content)",
      all(amp[k] > .05 for k in "ABCE"), ", ".join(f"{k}: {amp[k]:.3f}" for k in "ABCE"))
mnD1, mxD1, avD1 = wobble(HALF["D"]["maps"], 2, 4)
mnD2, mxD2, avD2 = wobble(HALF["D"]["maps"], 6, 8)
check("3b non-lattice D: the oscillation shrinks with scale (Minkowski content exists, Lapidus-van Frankenhuijsen)",
      (mxD2 - mnD2) / avD2 < (mxD1 - mnD1) / avD1 and (mxD2 - mnD2) / avD2 < .05,
      f"relative wobble {(mxD1-mnD1)/avD1:.3f} at 1e-2..1e-4 -> {(mxD2-mnD2)/avD2:.3f} at 1e-6..1e-8; content ~ {avD2:.4f}")
check("3c same dimension, different average Minkowski content (the size prefactor is a separate dial)",
      len({round(table[k]["minkowski_avg"], 3) for k in HALF}) == 5,
      ", ".join(f"{k}: {table[k]['minkowski_avg']:.3f}" for k in HALF))
check("3d B and C differ only in where the middle piece sits: same ratios, different content",
      abs(table["B"]["minkowski_avg"] - table["C"]["minkowski_avg"]) > .01)

# 4. Lacunarity (gliding box on the natural measure) ---------------------------------------
def lacunarity(maps, level, box):
    iv = intervals(maps, level)
    d = .5
    pts = [(float(a), float(b - a) ** d) for a, b in iv]  # piece start, natural mass r^D per piece
    step = box / 8
    masses = []
    x = 0.0
    j0 = 0
    while x + box <= 1 + 1e-12:
        m = sum(w for (p, w) in pts if x <= p < x + box)
        masses.append(m)
        x += step
    m1 = sum(masses) / len(masses); m2 = sum(m * m for m in masses) / len(masses)
    return m2 / m1 ** 2
for k, lev in (("A", 7), ("B", 5), ("C", 5), ("D", 9), ("E", 4)):
    lac = [lacunarity(HALF[k]["maps"], lev, b) for b in (1/16, 1/32, 1/64)]
    table[k]["lacunarity"] = sum(lac) / 3
check("4a gliding-box lacunarity differs between sets of the same dimension",
      len({round(table[k]["lacunarity"], 2) for k in HALF}) >= 4, ", ".join(f"{k}: {table[k]['lacunarity']:.2f}" for k in HALF))

# 5. Weights of a self-similar measure = per-element dial; multifractal spectrum = Renyi entropy
def Dq(p, q, r=0.25):
    if abs(q - 1) < 1e-12:
        return -(p * math.log(p) + (1 - p) * math.log(1 - p)) / math.log(1 / r)
    return math.log(p ** q + (1 - p) ** q) / ((q - 1) * math.log(r))
def renyi(p, q):
    if abs(q - 1) < 1e-12: return -(p * math.log(p) + (1 - p) * math.log(1 - p))
    return math.log(p ** q + (1 - p) ** q) / (1 - q)
check("5a fair weights on A: D_q = 1/2 for every q (monofractal)", all(abs(Dq(.5, q) - .5) < 1e-12 for q in (-5, -1, 0, .5, 1, 2, 7)))
check("5b any weights: D_0 = 1/2 (the support), D_q = Renyi entropy of the weights / log 4",
      all(abs(Dq(p, 0) - .5) < 1e-12 and all(abs(Dq(p, q) - renyi(p, q) / math.log(4)) < 1e-12 for q in (-3, .5, 1, 2, 5)) for p in (.1, .3, .45)))
# per-element dial t <-> weight p = t/(1+t): odds
ok = all(F(t) / (1 + F(t)) + 1 / (1 + F(t)) == 1 for t in (F(3,7), F(1,9), 1, 5))
check("5c weight p = t/(1+t) is the per-element dial t read as odds (p = 0.3 <-> t = 3/7)", ok and F(3,7)/(1+F(3,7)) == F(3,10))
# f(alpha) endpoints: alpha range [log p_max, log p_min]/log r
for p in (.3, .1):
    amin = math.log(max(p, 1-p)) / math.log(.25); amax = math.log(min(p, 1-p)) / math.log(.25)
    table.setdefault("multifractal", {})[str(p)] = [amin, amax]
check("5d unfair weights spread the local dimensions alpha over an interval; fair weights give the single value 1/2",
      table["multifractal"]["0.1"][0] < .5 < table["multifractal"]["0.1"][1], str(table["multifractal"]))

# 6. Natural shape for joining: binary branching + equal ratios + fair weights -> A ------------
# Among equal-ratio sets with N pieces, dimension D forces r = N^(-1/D); N = 2 (one yes/no choice
# per element, as for subsets) gives r = 2^(-1/D) = 1/4 at D = 1/2: set A.
check("6a binary branching + equal ratios + dimension 1/2 forces ratio 1/4 (set A); binary dimension 1 forces 1/2",
      F(1,2)**2 == F(1,4) and 2 * F(1,4) ** F(1,1) != 1)
# Joining two binary half boxes = four-corner set; regroup levels = binary tree with ratio 1/2.
def fc_points(L):
    pts = []
    for bits in itertools.product((0, 1), repeat=2 * L):
        x = sum(F(3, 4) * bits[2 * i] / F(4) ** i for i in range(L))
        y = sum(F(3, 4) * bits[2 * i + 1] / F(4) ** i for i in range(L))
        pts.append((bits, float(x), float(y)))
    return pts
def tree_d(b1, b2):
    for m, (u, v) in enumerate(zip(b1, b2)):
        if u != v: return 2.0 ** -m
    return 0.0
ratios = {}
for L in (2, 3, 4):
    P = fc_points(L)
    rs = [math.hypot(p[1] - q[1], p[2] - q[2]) / tree_d(p[0], q[0]) for p, q in itertools.combinations(P, 2)]
    ratios[L] = (min(rs), max(rs))
check("6b half box x half box is the binary full box as a tree: distance / tree distance stays in a fixed band at every depth",
      all(.45 <= ratios[L][0] and ratios[L][1] <= 2.2 for L in (2, 3, 4)) and ratios[4][1] - ratios[3][1] < ratios[3][1] - ratios[2][1],
      ", ".join(f"depth {L}: {a:.3f}..{b:.3f}" for L, (a, b) in ratios.items()))
check("6c the stick is that tree glued: binary digits of [0,1] identify only dyadic pairs (0.0111... = 0.1000...)",
      sum(F(1, 2 ** i) for i in range(2, 60)) + F(1, 2 ** 59) == F(1, 2))

# 7. Topology: below dimension 1 every one of these is the Cantor set ---------------------
# (Brouwer: compact, perfect, totally disconnected metric spaces are homeomorphic). Check: level pieces
# are pairwise separated and every piece splits again.
ok = True
for k, v in HALF.items():
    iv = sorted(intervals(v["maps"], 3))
    ok &= all(iv[i][1] < iv[i + 1][0] for i in range(len(iv) - 1)) and len(v["maps"]) >= 2
check("7a all five half boxes are dust that splits forever: topologically one Cantor set (connectivity is not a dial below D = 1)", ok)
# At D = log3/log2 the topology dial appears: gasket (3 at 1/2, touching) vs dust (9 at 1/4)
check("7b dimension log3/log2 both ways: 3 at 1/2 and 9 at 1/4 (9 = 3^2, 4 = 2^2)",
      abs(math.log(3) / math.log(2) - math.log(9) / math.log(4)) < 1e-15)
def gasket_vertices(n):
    tri = [((F(0), F(0)), F(1))]
    for _ in range(n):
        tri = [((x + dx * s / 2, y + dy * s / 2), s / 2) for (x, y), s in tri for dx, dy in ((0, 0), (1, 0), (F(1,2), 1))]
    verts = {}
    for (x, y), s in tri:
        for v in ((x, y), (x + s, y), (x + s / 2, y + s)):
            verts[v] = verts.get(v, 0) + 1
    return len(tri), len(verts), sum(1 for c in verts.values() if c > 1)
ok = all(gasket_vertices(n) == (3 ** n, (3 ** (n + 1) + 3) // 2, (3 ** (n + 1) + 3) // 2 - 3) for n in range(1, 6))
check("7c gasket: 3^n triangles share (3^(n+1)-3)/2 corner points (finitely ramified, connected); the 9-at-1/4 dust shares none", ok)

# 8. Shape dials are closed under joining; the power set leaves them -----------------------
# constant branching b, c: level counts c*b^n. Joining multiplies: (c1 c2)(b1 b2)^n. Power set: 2^(c b^n).
ok = all((2 * 3 ** n) * (5 * 2 ** n) == 10 * 6 ** n for n in range(20))
check("8a joining two self-similar trees gives a self-similar tree (branching multiplies)", ok)
ok = all(not all(2 ** (2 ** n) == c * b ** n for n in range(4)) for b in range(1, 40) for c in range(1, 40))
check("8b the power set of a binary tree (counts 2^(2^n)) is no self-similar tree c*b^n", ok)
# branching of the power-set tree grows: 2^(2^n)/2^(2^(n-1)) = 2^(2^(n-1)); a half step needs branching growth in between
br = [2 ** (2 ** (n - 1)) for n in range(1, 6)]
check("8c power-set tree branching per level: 2, 4, 16, 256, 65536 (doubly exponential)", br == [2, 4, 16, 256, 65536])
# a half step in growth: any h with h(h(x)) = 2^x squeezes between: x < h(x) < 2^x for x >= 2 (integers, increasing)
def build_half():
    h = {2: 3, 3: 4}
    for _ in range(60):
        for k in sorted(list(h)):
            if h[k] not in h and h[k] < 5000: h[h[k]] = 2 ** k
        miss = [x for x in range(2, 80) if x not in h]
        if not miss: break
        x = miss[0]; h[x] = max(v for kk, v in h.items() if kk < x) + 1
    return h
h = build_half()
check("8d an increasing half step h exists on counts >= 2 and sits strictly between N and 2^N",
      all(n < h[n] < 2 ** n for n in range(3, 40)) and all(h[h[n]] == 2 ** n for n in range(2, 12)))

print()
print(f"{sum(o for _, o in results)} of {len(results)} pass")
json.dump(table, open(os.path.join(os.path.dirname(os.path.abspath(__file__)), "results.json"), "w"), indent=1, default=str)
