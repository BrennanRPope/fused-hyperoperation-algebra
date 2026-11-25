"""Infinities between countable and the continuum, seen through fused mult-add.

Numeric and exact checks behind index.html. Run:  python3 checks.py > checks_output.txt

Mult-add  a (+) b = a + b + ab = (1+a)(1+b) - 1.
Weighted (zoomed) mult-add  a (+)_t b = a + b + t ab ; t = 0 is addition, t = 1 is mult-add.
Fractional rungs B_s from ../fractional-rungs/regular_iteration.py (s = 0 addition, s = 1 mult-add).
"""
import math, sys, itertools, random, json
from fractions import Fraction as Fr
from math import comb, log
import numpy as np

sys.path.insert(0, "../fractional-rungs")
from regular_iteration import Es  # E^s, regular iteration of e^x - 1

RES = {}
def rep(name, ok, extra=""):
    print(("PASS " if ok else "FAIL ") + name + (("   " + extra) if extra else ""))
    RES[name] = bool(ok)

def madd(a, b): return a + b + a*b
def wmadd(t, a, b): return a + b + t*a*b

# ---------------------------------------------------------------- 1. sum of ones vs mult-add of ones
ok = True
for n in range(0, 40):
    s = 0; m = 0
    for _ in range(n):
        s = s + 1
        m = madd(m, 1)
    ok &= (s == n) and (m == 2**n - 1)
rep("n ones added = n ; n ones mult-added = 2^n - 1 (all corners of the n-cube but the origin)", ok)

# corners of the n-cube by grade (number of switches on) = Pascal layers = mult-add cube layers
ok = True
for n in range(1, 16):
    layers = [0]*(n+1)
    for c in itertools.product([0, 1], repeat=n):
        layers[sum(c)] += 1
    ok &= layers == [comb(n, k) for k in range(n+1)]
rep("grade-k corners of the n-cube = C(n,k)  (n<=15, brute force; n=3 gives 1,3,3,1)", ok)

# the grade-1 layer is the 'addition part': exactly n of the 2^n - 1 corners
rep("grade-1 corners = n = the sum-of-ones part; every other corner is a 'corner term'",
    all(comb(n, 1) == n for n in range(1, 50)))

# ---------------------------------------------------------------- 2. Cantor's diagonal, finite shadow
random.seed(7); ok = True
for trial in range(2000):
    n = random.randint(1, 30)
    rows = [[random.randint(0, 1) for _ in range(n)] for _ in range(n)]
    diag = [1 - rows[i][i] for i in range(n)]
    ok &= all(rows[i] != diag for i in range(n))
rep("diagonal flip misses every row of any n x n 0/1 table (2000 random tables)", ok)

# ---------------------------------------------------------------- 3. weighted mult-add: exact identities
ok = True
for t in [Fr(0), Fr(1, 3), Fr(1, 2), Fr(1), Fr(2), Fr(7, 5)]:
    w = Fr(0)
    for n in range(0, 30):
        ok &= (t*w + 1 == (1 + t)**n)
        if t != 0:
            ok &= (w == ((1+t)**n - 1)/t)
        else:
            ok &= (w == n)
        w = wmadd(t, w, 1)
rep("n-fold weighted mult-add of 1:  t*W_n + 1 = (1+t)^n  (t=0 gives W_n = n)", ok)

ok = True
for t in [Fr(1, 3), Fr(1, 2), Fr(2)]:
    for a, b in [(Fr(1), Fr(2)), (Fr(3, 7), Fr(5, 2)), (Fr(0), Fr(4))]:
        ok &= t*wmadd(t, a, b) == madd(t*a, t*b)
rep("zoom: t*(a (+)_t b) = (ta) (+) (tb) -- weighted mult-add is mult-add seen at scale t", ok)

ok = all(sum(comb(n, k)*Fr(t)**k for k in range(n+1)) == (1+Fr(t))**n
         for n in range(0, 25) for t in [Fr(1, 3), Fr(1), Fr(5, 2)])
rep("weighted corner count  sum_k C(n,k) t^k = (1+t)^n", ok)

# ---------------------------------------------------------------- 4. Besicovitch-Eggleston dimension from the cube's layers
def H(p):  # natural-log entropy
    return 0.0 if p in (0, 1) else -p*log(p) - (1-p)*log(1-p)

def D_closed(t):  # dimension of the dominant layer at weight t, closed form
    return (math.log1p(t) - (t*log(t) if t > 0 else 0.0)/(1+t))/log(2)

worst = 0
for t in [0.05, 0.2, 0.5, 1.0, 0.8]:
    p = t/(1+t)
    worst = max(worst, abs(D_closed(t) - H(p)/log(2)))
rep("closed form  D(t) = [ln(1+t) - t ln t/(1+t)] / ln 2  equals  H(t/(1+t))/ln 2", worst < 1e-12, f"worst {worst:.1e}")

# dominant layer of (1+t)^n sits at k/n -> t/(1+t)
ok = True; notes = []
for t in [0.1, 0.25, 0.5, 1.0]:
    n = 4000
    logs = [math.lgamma(n+1) - math.lgamma(k+1) - math.lgamma(n-k+1) + k*log(t) for k in range(n+1)]
    kmax = int(np.argmax(logs))
    notes.append(f"t={t}: k/n={kmax/n:.4f} vs {t/(1+t):.4f}")
    ok &= abs(kmax/n - t/(1+t)) < 1e-3
rep("heaviest layer of the weighted cube sits at density p = t/(1+t) (n=4000)", ok, "; ".join(notes))

# layer growth: (1/n) log2 C(n, pn) -> H(p)/ln2 (the counting fact behind Eggleston's theorem)
rows = []
for p in [0.05, 0.2, 1/3, 0.5]:
    vals = []
    for n in [100, 1000, 10000, 100000]:
        k = round(p*n)
        vals.append((math.lgamma(n+1) - math.lgamma(k+1) - math.lgamma(n-k+1))/(n*log(2)))
    rows.append((p, vals, H(p)/log(2)))
ok = all(abs(v[-1] - target) < 2e-4 for _, v, target in rows)
rep("(1/n) log2 C(n, pn) -> H(p)/ln 2 as n grows", ok,
    "; ".join(f"p={p:.3f}: {v[0]:.4f},{v[1]:.4f},{v[2]:.4f},{v[3]:.5f} -> {tg:.5f}" for p, v, tg in rows))

# smooth and increasing on [0,1]: zoom 0 (addition) -> dimension 0, zoom 1 (mult-add) -> dimension 1
ts = np.linspace(1e-6, 1, 2001)
Ds = np.array([D_closed(t) for t in ts])
rep("D(t) rises smoothly from 0 (t -> 0, addition) to 1 (t = 1, mult-add)",
    np.all(np.diff(Ds) > 0) and Ds[0] < 1e-4 and abs(Ds[-1] - 1) < 1e-12,
    f"D(0.01)={D_closed(0.01):.4f} D(0.1)={D_closed(0.1):.4f} D(0.5)={D_closed(0.5):.4f} D(1)={D_closed(1):.4f}")

# Monte Carlo sanity: random density-p strings mostly land in the dominant layer band
rng = np.random.default_rng(3); n = 20000
for t in [0.25, 1.0]:
    p = t/(1+t)
    x = rng.random((200, n)) < p
    freq = x.mean(axis=1)
    rep(f"Bernoulli({p:.3f}) strings concentrate at density {p:.3f} (t={t})", np.all(np.abs(freq - p) < 0.02),
        f"spread {freq.min():.4f}..{freq.max():.4f}")

# ---------------------------------------------------------------- 5. fractional rungs: counts between n and 2^n
# N_s(n) = B_s-fold of n ones = E^s(n * E^-s(1)).   log(1 + N_s(n)) = E^(s-1)(n * E^-s(1))  (no overflow)
def logN(s, n):
    c = float(Es(-s, np.array([1.0]))[0])
    return float(Es(s - 1, np.array([n*c]))[0])  # natural log of 1 + N_s(n)

ok = True
for n in [1, 2, 5, 10, 40]:
    ok &= abs(logN(0.0, n) - math.log1p(n)) < 1e-7          # s = 0: N = n
    ok &= abs(logN(1.0, n) - n*log(2)) < 1e-7               # s = 1: N = 2^n - 1
rep("rung 0 count = n (countable sum), rung 1 count = 2^n - 1 (continuum corners)", ok)

table = {}
ns = [10, 100, 1000, 10**4, 10**5, 10**6]
for s in [0.0, 0.25, 0.5, 0.75, 1.0]:
    table[s] = [logN(s, n) for n in ns]
print("log(1+N_s(n)) table, n =", ns)
for s, v in table.items():
    print(f"  s={s:4}: " + "  ".join(f"{x:12.4f}" for x in v))

# strictly ordered in s at each n
ok = all(table[a][i] < table[b][i] for i in range(len(ns))
         for a, b in [(0.0, .25), (.25, .5), (.5, .75), (.75, 1.0)])
rep("counts strictly ordered by rung: N_0 < N_.25 < N_.5 < N_.75 < N_1 at every n tested", ok)

# sub-exponential for s < 1: log N / n falls toward 0
ok = True; notes = []
for s in [0.25, 0.5, 0.75]:
    r = [table[s][i]/ns[i] for i in range(len(ns))]
    ok &= all(r[i+1] < r[i] for i in range(len(r)-1)) and r[-1] < 0.05
    notes.append(f"s={s}: " + ",".join(f"{x:.2e}" for x in r))
rep("s<1: log N_s(n) / n keeps falling (entropy dimension 0, slower than any 2^(eps n) in the trend)", ok, "; ".join(notes))

# super-polynomial for s > 0: log N / log n keeps rising
ok = True; notes = []
for s in [0.25, 0.5, 0.75]:
    r = [table[s][i]/log(ns[i]) for i in range(len(ns))]
    ok &= all(r[i+1] > r[i] for i in range(len(r)-1))
    notes.append(f"s={s}: " + ",".join(f"{x:.2f}" for x in r))
rep("s>0: log N_s(n) / log n keeps rising (faster than any fixed power n^k in the trend)", ok, "; ".join(notes))

# gap between neighbouring rungs widens: N_s / N_s' -> 0 for s < s'
ok = True; notes = []
for a, b in [(0.25, 0.5), (0.5, 0.75)]:
    gaps = [table[b][i] - table[a][i] for i in range(len(ns))]
    ok &= all(gaps[i+1] > gaps[i] for i in range(len(gaps)-1))
    notes.append(f"{a}->{b}: " + ",".join(f"{g:.2f}" for g in gaps))
rep("log N_s' - log N_s grows with n (each rung outgrows the one below)", ok, "; ".join(notes))

# half-iterate sandwich used in the prose proof: x <= E^1/2(x) <= e^x - 1 for x > 0
xs = np.array([0.1, 0.5, 1, 2, 5, 10.0])
h = Es(0.5, xs)
rep("x <= E^(1/2)(x) <= e^x - 1 on x > 0 (samples)", np.all(xs <= h + 1e-12) and np.all(h <= np.expm1(xs) + 1e-9),
    ", ".join(f"{v:.4f}" for v in h))

# ---------------------------------------------------------------- 6. same count, different shape: comb vs uniform tree
def comb_tree_levels(n):   # spine keeps splitting off one dead-end branch per level
    return [k + 1 for k in range(n)]
def uniform_split_levels(target):  # every branch splits at once when the target count allows
    lv = []; cur = 1
    for T in target:
        while 2*cur <= T: cur *= 2
        lv.append(cur)
    return lv

n = 64
cl = comb_tree_levels(n); ul = uniform_split_levels(cl)
rep("comb tree and uniform tree track the same count (within a factor 2 at every level)",
    all(c/2 < u <= c for c, u in zip(cl, ul)), f"level 63: comb {cl[-1]}, uniform {ul[-1]}")
# branches surviving to depth d: comb has d+1 (one new dead end per level, countable in the limit);
# uniform tree with k splits has 2^k paths and every path keeps splitting (perfect -> continuum)
# build both trees to depth n as sets of 0/1 prefixes; ask how many depth-d nodes still split somewhere below
def comb_nodes(n):
    lv = [[()]]
    for d in range(n):
        nxt = []
        for node in lv[-1]:
            if all(b == 0 for b in node):      # spine: splits every level
                nxt += [node + (0,), node + (1,)]
            else:                               # dead end: one child forever
                nxt.append(node + (0,))
        lv.append(nxt)
    return lv
def uniform_nodes(n, split_levels):
    lv = [[()]]
    for d in range(n):
        lv.append([node + (b,) for node in lv[-1] for b in ((0, 1) if d in split_levels else (0,))])
    return lv
def still_splitting(levels, d):
    n = len(levels) - 1
    count = 0
    for node in levels[d]:
        kids_by_depth = [[x for x in levels[D] if x[:d] == node] for D in range(d, n + 1)]
        if any(len(kids_by_depth[i+1]) > len(kids_by_depth[i]) for i in range(len(kids_by_depth) - 1)):
            count += 1
    return count
N = 20; d = 6
ct = comb_nodes(N)
spl = {k for k in range(N) if (k + 1) & k == 0}   # split at levels 0,1,3,7,15 (sizes 2,4,8,16,32)
ut = uniform_nodes(N, spl)
cs, us = still_splitting(ct, d), still_splitting(ut, d)
rep("comb: of the depth-6 nodes only the spine splits again (addition, countable); uniform: all do (mult-add, continuum)",
    cs == 1 and us == len(ut[d]), f"comb {cs} of {len(ct[d])}, uniform {us} of {len(ut[d])}; shape, not count, decides (perfect set theorem)")

# gauge sizes of Cantor sets built from rung counts: C_s covered by L_s(n) pieces of size 2^-n.
# h_s(2^-n) = 1/L_s(n);  sum over C_s's own cover in the h_s' gauge = L_s(n)/L_s'(n) -> 0 for s < s'
ok = True; notes = []
for a, b in [(0.0, 0.25), (0.25, 0.5), (0.5, 0.75), (0.75, 1.0)]:
    lr = [table[a][i] - table[b][i] for i in range(len(ns))]   # natural log of the cover sum
    ok &= all(lr[i+1] < lr[i] for i in range(len(lr)-1))
    notes.append(f"C_{a} in gauge h_{b}: log cover sum " + ",".join(f"{r:.1f}" for r in lr))
rep("each rung's Cantor set has zero size in the next rung's gauge (cover sums shrink)", ok, "; ".join(notes))

print()
print(f"{sum(RES.values())} of {len(RES)} checks pass")
json.dump(RES, open("checks_results.json", "w"), indent=1)
