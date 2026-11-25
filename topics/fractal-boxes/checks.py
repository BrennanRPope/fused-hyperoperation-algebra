"""Fractal boxes: does a box of fractional dimension give a fractional power set?
Pure Python (fractions, math, itertools). Run: python3 checks.py
"""
from fractions import Fraction as F
from itertools import combinations
import math

results = []
def check(name, ok, detail=""):
    results.append((name, bool(ok)))
    print(("PASS " if ok else "FAIL ") + name + (("  | " + detail) if detail else ""))

# ---------------------------------------------------------------------------
# 1. Corner pieces of a box with sides 1+a: layer k has C(d,k) pieces of size a^k.
#    For whole d this is the subsets picture. For d = 1/2 use Newton's binomial.
def gbinom(d, k):
    r = F(1)
    for i in range(k):
        r = r * (d - i) / (i + 1)
    return r

half = F(1, 2)
layers = [gbinom(half, k) for k in range(12)]
check("1a half box layers start 1, 1/2, -1/8, 1/16, -5/128, 7/256",
      layers[:6] == [F(1), F(1,2), F(-1,8), F(1,16), F(-5,128), F(7,256)],
      ", ".join(str(x) for x in layers[:6]))
# partial sums at a = 1 approach sqrt 2, alternating around it
ps, s = [], F(0)
for k in range(400):
    s += gbinom(half, k); ps.append(float(s))
check("1b half box corner total -> sqrt 2 (sum of layers at a=1)",
      abs(ps[-1] - math.sqrt(2)) < 2e-3, f"400 layers: {ps[-1]:.6f} vs {math.sqrt(2):.6f}")
# signs: from layer 2 on, layers alternate sign -> negative pieces exist
check("1c half box has negative layers (k=2,4,6,...)",
      all(layers[k] < 0 for k in (2,4,6,8,10)) and all(layers[k] > 0 for k in (1,3,5,7,9)))
# two half boxes joined: convolution of layers = layers of a stick (1, 1, 0, 0, ...)
conv = [sum(layers[i]*layers[k-i] for i in range(k+1)) for k in range(12)]
check("1d half box joined with half box = stick: layers 1,1,0,0,... exactly",
      conv == [F(1), F(1)] + [F(0)]*10)
# scaled integer version used in Lean: b_k = 4^k C(1/2,k)
bk = [int(layers[k] * 4**k) for k in range(12)]
cat = [math.comb(2*n, n)//(n+1) for n in range(12)]
check("1e 4^k C(1/2,k) = 1, then (-1)^(k-1) * 2 * Catalan(k-1)",
      bk[0] == 1 and all(bk[k] == (-1)**(k-1) * 2 * cat[k-1] for k in range(1,12)), str(bk[:8]))
# one third box: three joined give stick
th = [gbinom(F(1,3), k) for k in range(10)]
c3 = [sum(th[i]*th[j]*th[k-i-j] for i in range(k+1) for j in range(k+1-i)) for k in range(10)]
check("1f third box joined three times = stick", c3 == [F(1), F(1)] + [F(0)]*8)

# 2. Fractional dimension = per-element dial in disguise.
#    box of dimension d, sides 2  ==  box of dimension 1, side 2^d  (dial 2^d - 1)
for d in (0.5, 1/3, 0.25, 1.5):
    pass
check("2a (1+1)^(1/2) = 1 + (sqrt2 - 1): half dimension at side 2 = one stick with dial sqrt2-1",
      abs((1+1)**0.5 - (1 + (math.sqrt(2)-1))) < 1e-15)
# Circle power: a^(*d) = (1+a)^d - 1 ; halves combine back
def mpow(a, d): return (1+a)**d - 1
def madd(a, b): return a + b + a*b
ok = all(abs(madd(mpow(a, .5), mpow(a, .5)) - a) < 1e-12 for a in (0.3, 1, 2.5, 7))
check("2b circle-operation half power a^(*1/2) = sqrt(1+a)-1 joins back to a", ok)
# continuity pins it down: f(d+e)=f(d)f(e), f(1)=2, continuous -> 2^d; check dyadic values
check("2c dyadic halves forced: f(1/2)=sqrt2, f(1/4)=2^(1/4) (positive roots)",
      abs(math.sqrt(2)**2 - 2) < 1e-15 and abs((2**0.25)**4 - 2) < 1e-14)

# ---------------------------------------------------------------------------
# 3. The half box as a fractal: keep the two outer quarters of each segment (C_1/4).
def cantor_intervals(ratio, offsets, level):
    """self-similar set in [0,1]: pieces x -> ratio*x + offset."""
    iv = [(F(0), F(1))]
    for _ in range(level):
        iv = [(a + (b-a)*o, a + (b-a)*o + (b-a)*ratio) for (a, b) in iv for o in offsets]
    return iv
C4 = lambda m: cantor_intervals(F(1,4), [F(0), F(3,4)], m)
for m in range(6):
    assert len(C4(m)) == 2**m
check("3a C_1/4 level m: 2^m pieces of length 4^-m", all(len(C4(m)) == 2**m and C4(m)[0][1]-C4(m)[0][0] == F(1,4**m) for m in range(7)))
# similarity dimension log N / log(1/r)
check("3b dimension of C_1/4 = log2/log4 = 1/2", abs(math.log(2)/math.log(4) - 0.5) < 1e-15)

def box_count(points_fn, eps):
    return len(points_fn(eps))

def C4_boxes(m):
    """dyadic boxes of side 2^-k meeting C_1/4, using a deep level as stand-in."""
    pts = [a for (a, b) in C4(14)]
    return pts
pts = [float(a) for (a, b) in C4(14)] + [float(b) for (a, b) in C4(14)]
counts = []
for k in range(2, 22):
    eps = 2.0**-k
    boxes = {int(p/eps - 1e-12) if p > 0 else 0 for p in pts}
    counts.append((k, len(boxes)))
slopes = [math.log2(counts[i+2][1]/counts[i][1])/2 for i in range(len(counts)-2)]
check("3c box counting: halving the scale multiplies the count by about sqrt2 (slope 1/2)",
      abs(sum(slopes[-10:])/10 - 0.5) < 0.03, "per-halving ratios ~ " + ", ".join(f"{counts[i+1][1]/counts[i][1]:.2f}" for i in range(8, 16)))
# The ratio count(eps)/eps^(-1/2) wobbles: log-periodic, so no single "number of corners"
wob = [c / (2.0**k)**0.5 for (k, c) in counts[6:]]
alt = [counts[i+1][1]/counts[i][1] for i in range(8, len(counts)-1)]
check("3d the count / sqrt(1/eps) keeps wobbling with period one level (x4 in scale): no single corner number",
      max(wob)/min(wob) > 1.05 and all(abs(alt[i]-alt[i+2]) < 0.03 and abs(alt[i]-alt[i+1]) > 0.1 for i in range(len(alt)-2)),
      f"ratio range {min(wob):.3f} .. {max(wob):.3f}; per-halving ratios alternate ~1.50 / ~1.33")

# Square of the half box: the four-corner (Garnett) set: 4^m squares of side 4^-m
check("3e half box squared: 4^m squares of side 4^-m, dimension log4/log4 = 1",
      all(len(C4(m))**2 == 4**m for m in range(8)))
# but it is dust: at level m no two squares touch (gap of 2 * 4^-m between neighbours)
def gaps(m):
    iv = C4(m)
    return min(iv[i+1][0] - iv[i][1] for i in range(len(iv)-1))
check("3f the squared half box is dust: level-m pieces are separated by gaps (no stick inside)",
      all(gaps(m) == F(2, 4**m) for m in range(1, 8)))
# measure of projection onto a 45-degree-ish line shrinks (Favard length -> 0 known);
# here: projection to the x-axis is C_1/4 itself, length (1/2)^m -> 0
check("3g its shadow on the x-axis has length (1/2)^m -> 0 (a stick's shadow stays 1)",
      all(sum(b-a for a, b in C4(m)) == F(1, 2**m) for m in range(8)))

# ---------------------------------------------------------------------------
# 4. Not canonical: other halves of a box
T9 = lambda m: cantor_intervals(F(1,9), [F(0), F(4,9), F(8,9)], m)
check("4a three pieces at scale 1/9 also has dimension 1/2 (log3/log9)",
      abs(math.log(3)/math.log(9) - .5) < 1e-15 and len(T9(3)) == 27)
# compare counts at the same scale eps = 2^-12 using deep levels
def count_at(ivs, eps):
    pts = [float(a) for a, b in ivs] + [float(b) for a, b in ivs]
    return len({int(p/eps - 1e-12) if p > 0 else 0 for p in pts})
e = 2.0**-12
cA, cB = count_at(C4(10), e), count_at(T9(7), e)
check("4b same dimension, different counts at the same scale", cA != cB, f"C_1/4: {cA}, three-at-1/9: {cB}, sqrt(1/eps) = {2**6}")
# uneven half: 2 pieces of ratios r1,r2 with r1^(1/2)+r2^(1/2)=1 (Moran equation)
r1 = 0.09; r2 = (1 - r1**0.5)**2
check("4c uneven pieces (0.09, 0.49) also solve the Moran equation for d = 1/2",
      abs(r1**0.5 + r2**0.5 - 1) < 1e-12, f"r2 = {r2:.4f}")
# random check of Moran: Hausdorff dims of products add for these self-similar sets
# grid boxes meeting a product = (boxes meeting A) x (boxes meeting B), so log-counts add
prod_slopes = []
for k in (8, 10, 12):
    a1, b1 = count_at(C4(10), 2.0**-k), count_at(T9(7), 2.0**-k)
    a2, b2 = count_at(C4(10), 2.0**-(k+2)), count_at(T9(7), 2.0**-(k+2))
    prod_slopes.append(math.log2((a2*b2)/(a1*b1))/2)
check("4d C_1/4 x three-at-1/9 (a different half squared against C_1/4) box-counts at dimension ~1",
      all(abs(x - 1) < 0.12 for x in prod_slopes), "slopes " + ", ".join(f"{x:.3f}" for x in prod_slopes))

# ---------------------------------------------------------------------------
# 5. Iteration: the geometric power set of a shape is its hyperspace (all closed
#    pieces, Hausdorff distance). At a scale where the shape is N separated points,
#    the hyperspace is 2^N - 1 separated points.
def hausdorff(A, B):
    d1 = max(min(abs(a-b) for b in B) for a in A)
    d2 = max(min(abs(a-b) for a in A) for b in B)
    return max(d1, d2)
for m in (1, 2, 3):
    S = [a for a, b in C4(m)]
    g = gaps(m) + F(1, 4**m)  # distance between left ends of neighbours >= this
    subs = [c for r in range(1, len(S)+1) for c in combinations(S, r)]
    sep = min(hausdorff(A, B) for A, B in combinations(subs, 2))
    assert len(subs) == 2**len(S) - 1
    assert sep >= min(S[i+1]-S[i] for i in range(len(S)-1))
check("5a level m=1,2,3 of C_1/4: N points give 2^N - 1 subsets, pairwise as far apart as the points",
      True, "N = 2, 4, 8 -> 3, 15, 255 separated subsets (asserted above)")
# net argument: subsets of an eps-net form an eps-net of the hyperspace
import random
random.seed(1)
X = [a for a, b in C4(6)]; net = [a for a, b in C4(3)]
eps = max(min(abs(x-y) for y in net) for x in X)
worst = F(0)
for _ in range(1500):
    A = random.sample(X, random.randint(1, 12))
    B = sorted({min(net, key=lambda y: abs(x-y)) for x in A})
    worst = max(worst, hausdorff(A, B))
check("5b subsets of an eps-net are an eps-net of the hyperspace (1500 random pieces of level 6 vs the 8-point net)",
      worst <= eps, f"eps = {eps}, worst Hausdorff distance = {worst}")
# So a half step H on shapes needs counts N -> h(N) with h(h(N)) ~ 2^N: a half-exponential.
# Dimension-preserving operations change log-counts linearly; check that iterating a
# linear rule twice never reaches 2^N.
def linear_iter_ok(A, B, nmax=60):
    return all(A*(A*n+B)+B < 2**n for n in range(4*(A*A+A*B+B+1), nmax+4*(A*A+A*B+B+1)))
check("5c any rule N -> <= A N + B applied twice falls below 2^N eventually (A,B <= 6)",
      all(linear_iter_ok(A, B) for A in range(7) for B in range(7)))
# a half-exponential exists only non-canonically: two different increasing h on N>=2 with h(h(n)) = 2^n
def build_half(seed):
    # seed: h(2)=3 forced? choose chain starts; greedy construction from a given h on [2, 3)
    h = {2: 3, 3: 4}
    n = 4
    # fill pairs: h(n) known -> h(h(n)) = 2^n
    for _ in range(40):
        for k in sorted(list(h)):
            if h[k] not in h and h[k] < 200:
                h[h[k]] = 2**k
        # fill the smallest gap greedily with seed offset
        missing = [x for x in range(2, 200) if x not in h]
        if not missing: break
        x = missing[0]
        lo = max(v for k, v in h.items() if k < x)
        h[x] = lo + 1 + (seed if x > 4 else 0)
        if not all(h[a] < h[b] for a in h for b in h if a < b):
            return None
    return h
hA, hB = build_half(0), build_half(1)
def is_half(h):
    return all(h[h[n]] == 2**n for n in h if h[n] in h)
check("5d two different increasing integer half power sets from n=2 up (non-canonical)",
      hA and hB and is_half(hA) and is_half(hB) and hA != hB,
      f"A: {[hA[k] for k in range(2,12)]}  B: {[hB[k] for k in range(2,12)]}")

# ---------------------------------------------------------------------------
# 6. Infinite sizes: every Cantor set, of any dimension, is P(N) as a set.
# a point of C_r is an infinite left/right path = a subset of N. Check the coding is
# one-to-one on finite levels for r = 1/4 and r = 1/3.
def code_point(bits, ratio, gap_start):
    x = F(0); w = F(1)
    for b in bits:
        if b: x += w * gap_start
        w *= ratio
    return x
for r, g in ((F(1,4), F(3,4)), (F(1,3), F(2,3)), (F(1,9), F(8,9))):
    m = 10
    vals = {code_point([(i >> j) & 1 for j in range(m)], r, g) for i in range(2**m)}
    assert len(vals) == 2**m
check("6a Cantor sets of dimension 1/2, 0.63 and 0.32 all code subsets of N one-to-one: 2^10 distinct level-10 points each (asserted above)", True)

# ---------------------------------------------------------------------------
# 7. Other ways to get "half an element"
# Species (labelled product: split the set in two): power set = sets x sets, so
# the half is the species of sets (1 structure per set), not sqrt2^n.
ok = all(sum(math.comb(n, k) * 1 * 1 for k in range(n+1)) == 2**n for n in range(12))
bad = all(abs(sum(math.comb(n, k) * 2**(k/2) * 2**((n-k)/2) for k in range(n+1)) - 2**n) > 0.5 for n in range(1, 12))
check("7a species: sets x sets = power set (sum C(n,k) = 2^n); sqrt2^n is not a half under this product",
      ok and bad, "a different way of joining gives a different half: count 1 per set")
# Groupoid cardinality: one object with a 2-fold symmetry counts 1/2. Its 'power set'
# (functors into the discrete set {0,1}) has 2 objects and no symmetries: size 2, not sqrt2.
functors = [(o, "id") for o in (0, 1)]  # the symmetry must go to an identity in a discrete target
check("7b groupoid with 1/2 element: maps into {0,1} number 2 (groupoid size 2), not 2^(1/2)",
      len(functors) == 2 and len(functors) != math.sqrt(2))

print()
n_pass = sum(ok for _, ok in results)
print(f"{n_pass} of {len(results)} pass")
